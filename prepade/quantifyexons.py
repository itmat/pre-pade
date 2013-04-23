#!/usr/bin/env python

from __future__ import print_function, division


import logging
import os
import sys
import pysam
import argparse
import re
import pandas as pd

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from collections import defaultdict, namedtuple
from prepade.geneio import parse_rum_index_genes

CIGAR_CHARS = 'MIDNSHP=X'

def load_exon_index(fh, index_filename=None):

    if index_filename is not None and os.path.exists(index_filename):
        logging.info("Loading exons from index at " + index_filename)
        return pd.read_table(index_filename)

    chr_to_exons = defaultdict(list)

    chrs = []
    starts = []
    ends = []

    for gene in parse_rum_index_genes(fh):
        for exon in gene.sub_features:
            chrs.append(gene.ref)
            starts.append(exon.location.start)
            ends.append(exon.location.end)

    df = pd.DataFrame({
        'chromosome' : chrs,
        'start'      : starts,
        'end'        : ends })

    df = df.sort(columns=['chromosome', 'end'])

    if index_filename is not None:
        df.to_csv(index_filename, sep='\t')

    return df



def iterate_over_exons(exons, sam_filename):
    samfile = pysam.Samfile(sam_filename)

    details = open('details', 'w')

    for i in range(len(exons)):
        chr_ = exons.chromosome[i]
        start = exons.start[i]
        end   = exons.end[i]

        read_ids = set()

        unique_read_ids = set()
        multi_read_ids = set()

        for aln in samfile.fetch(str(chr_), start, end):

            num_alns = aln.opt('IH')
            if num_alns > 1:
                multi_read_ids.add(aln.qname)
            else:
                unique_read_ids.add(aln.qname)

        count_u = len(unique_read_ids)
        count_m = len(multi_read_ids)

        if count_u > 0:
            loc = FeatureLocation(start, end)
            
            feat = SeqFeature(ref=str(chr_), location=loc)
#            print('feature:', feat, 'unique:', count_u,
#                  'unique_reads:', unique_read_ids,
#                  'multi_reads:', multi_read_ids,
#                  'multi:', count_m, file=details)
            yield(feat, count_u, count_m)
    details.close()

def read_sam_file(gene_filename, sam_filename, output_fh):

    samfile = pysam.Samfile(sam_filename)

    for aln in samfile.fetch():

        if aln.cigar is not None:
            spans = cigar_to_spans(aln.cigar, aln.pos)
        else:
            spans = None

        print(aln.qname)
        print(samfile.mate(aln))
        
#        for exon in exons:
#            start = exon.location.start
#            end   = exon.location.end
#            count = 0
#            for rec in samfile.fetch(chr_, start, end):
#                process_sam_rec(rec)
#                count += 1#
#
#            if count > 0:
#                print(chr_, start, end, count, sep="\t", file=output_fh)

CigarElem = namedtuple('CigarElem', ['count', 'op'])

def cigar_to_spans(cigar, start):

    spans = []

    for (op, bases) in cigar:
        opname = CIGAR_CHARS[op]
        if opname == 'M':
            end = start + bases - 1
            spans.append(FeatureLocation(start, end))
            start = end

        elif opname in 'DN':
            start = start + bases + 1

        elif opname == 'I':
            start += 1

    feats = []

    for span in spans:
        if len(feats) > 0 and span.start < feats[-1].location.end:
            feats[-1].location.end = span.end
        else:
            feats.append(SeqFeature(location=span))

    return SeqFeature(sub_features=feats)

        

def remove_ds(cigar):
    """Removes D operations from Cigar string, replacing with Ms.

    Replaces all Ds with Ms and then merges adjacent Ms together.

    >>> remove_ds([ (0, 21), (2, 1), (0, 54) ])
    [(0, 76)]

    >>> remove_ds([(4, 8), (0, 4), (2, 1), (0, 63)])
    [(4, 8), (0, 68)]

    >>> remove_ds([(0, 15), (2, 1), (0, 15), (2, 2), (0, 29), (2, 2), (0, 16)])
    [(0, 80)]

    >>> remove_ds([(0, 21), (2, 1), (0, 41), (3, 177), (0, 13)])
    [(0, 63), (3, 177), (0, 13)]

    >>> remove_ds([(0, 41), (3, 354), (0, 20), (2, 1), (0, 14)])
    [(0, 41), (3, 354), (0, 35)]

    >>> remove_ds([(4, 4), (0, 26), (2, 1), (0, 45)])
    [(4, 4), (0, 72)]

    """

    M = 0
    D = 2

    d_to_m = lambda (op, bases): (M, bases) if op == D else (op, bases)
    converted = map(d_to_m, cigar)
    res = []

    res = [converted[0]]

    for (op, bases) in converted[1:]:
        (last_op, last_bases) = res[-1]
        
        if op == last_op:
            res[-1] = (op, bases + last_bases)
        else:
            res.append((op, bases))

    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rum-gene-info', type=file)
    parser.add_argument('--exon-index')
    parser.add_argument('samfile')

    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    args = parser.parse_args()

    output = args.output if args.output is not None else sys.stdout

    exon_index = load_exon_index(args.rum_gene_info, index_filename=args.exon_index)
    print('feature', 'min', 'max', sep='\t', file=output)
    for (exon, count_u, count_m) in iterate_over_exons(exon_index, args.samfile):
        exon_str = '{chr_}:{start}-{end}'.format(
            chr_=exon.ref,
            start=exon.location.start,
            end=exon.location.end)
        print(exon_str, count_u, count_u + count_m, sep='\t', file=output)

#    read_sam_file(args.rum_gene_info, args.samfile, args.output)

