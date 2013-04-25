#!/usr/bin/env python

from __future__ import print_function, division

import logging
import os
import sys
import pysam
import argparse
import re
import pandas as pd
import numpy as np

from itertools import groupby
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from collections import defaultdict, namedtuple
from prepade.geneio import parse_rum_index_genes, read_exons

CIGAR_CHARS = 'MIDNSHP=X'

def load_exon_index(fh):

    chr_to_exons = defaultdict(list)

    chrs = []
    starts = []
    ends = []

    genes = parse_rum_index_genes(fh)
    exons = genes_to_exons(genes)

    for exon in exons:
        chrs.append(gene.ref)
        starts.append(exon.location.start)
        ends.append(exon.location.end)

    df = pd.DataFrame({
        'chromosome' : chrs,
        'start'      : starts,
        'end'        : ends })

    df = df.sort(columns=['chromosome', 'end'])

    return df

def genes_to_exons(genes):
    """Given an iterator over genes, returns an iterator over exons.

    """
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)


def iterate_over_exons(exons, sam_filename):
    samfile = pysam.Samfile(sam_filename)
    details = open('details', 'w')

    seen = set()

    for i, exon in enumerate(exons):

        if i % 1000 == 0:
            logging.info("Done {i} exons".format(i=i))

        key = exon.ref + ':' + str(exon.location.start) + '-' + str(exon.location.end)

        if key in seen:
            continue
        seen.add(key)

        unique_read_ids = set()
        multi_read_ids = set()

        all_spans = []

        try:
            alns = list(samfile.fetch(exon.ref, exon.location.start, exon.location.end))
        except ValueError as e:
            logging.warn('error fetching reads for {ref}:{start}-{end}'.format(
                    ref=exon.ref, 
                    start=exon.location.start,
                    end=exon.location.end), exc_info=e)
            alns = []
        key_fn = lambda x: (x.qname, x.opt('HI'))

        alns = sorted(alns, key=key_fn)

        logging.debug('Got' + str(len(alns)) + 'candidate alignments')

        for (qname, hi), pair in groupby(alns, key=key_fn):

            logging.debug('  candidate %s[%d]', qname, hi)
            pair = list(pair)
            if match(exon, pair):
                logging.debug('    found a match')
                num_alns = pair[0].opt('IH')
                if num_alns > 1:
                    multi_read_ids.add((qname, hi))
                else:
                    unique_read_ids.add((qname, hi))

        count_u = len(unique_read_ids)
        count_m = len(multi_read_ids)

        yield(exon, count_u, count_m)
    details.close()


def read_sam_file(gene_filename, sam_filename, output_fh):

    samfile = pysam.Samfile(sam_filename)

    for aln in samfile.fetch():

        if aln.cigar is not None:
            spans = cigar_to_spans(aln.cigar, aln.pos, strand)
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

def cigar_to_spans(cigar, start, strand):
    spans = []

    if cigar is None:
        return SeqFeature()

    cigar = remove_ds(cigar)

    for (op, bases) in cigar:
        opname = CIGAR_CHARS[op]
        if opname == 'M':
            end = start + bases
            spans.append(FeatureLocation(start, end))
            start = end

        elif opname in 'DN':
            start = start + bases

        elif opname == 'I':
            start += 1

    feats = []

    for span in spans:
        if len(feats) > 0 and feats[-1].location.end + 1>= span.start:
            start = feats[-1].location.start
            end   = span.end
            loc   = FeatureLocation(start, end)
            feats[-1] = SeqFeature(location=loc)
        else:
            feats.append(SeqFeature(location=span))
    start = feats[0].location.start
    end   = feats[-1].location.end

    return SeqFeature(sub_features=feats)

    
def match(exon, alns):
    stack = []

    overlaps = []
    consistents = []

    for aln in alns:
        strand = -1 if aln.is_reverse else 1
        spans = cigar_to_spans(aln.cigar, aln.pos, strand).sub_features
        debug_spans = ', '.join([ str(s.location.start) + '-' + str(s.location.end) for s in spans ])
        logging.debug('    spans are %s', debug_spans)
        overlaps.extend(spans_overlap(exon, spans))
        consistents.extend(spans_are_consistent(exon, spans))

    logging.debug('overlap %s', overlaps)
    logging.debug('consistent %s', consistents)

    overlaps = np.array(overlaps)
    consistents = np.array(consistents)

    return any(overlaps) and all(consistents)

        
def spans_overlap(exon, spans):
    lexon = exon.location.start
    rexon = exon.location.end

    for span in spans:
        lspan = span.location.start
        rspan = span.location.end
        yield not (rspan <= lexon or lspan >= rexon) 


def spans_are_consistent(exon, spans):

    lexon = exon.location.start
    rexon = exon.location.end

    last_span = len(spans) - 1

    for i, span in enumerate(spans):
        lspan = span.location.start
        rspan = span.location.end
        
        is_first = i == 0
        is_last  = i == last_span

        if rspan <= lexon or lspan >= rexon:
            yield True

        else:
            if i == 0:
                l_ok = lspan >= lexon
            else:
                l_ok = lspan == lexon

            if i == last_span:
                r_ok = rspan <= rexon
            else:
                r_ok = rspan == rexon

            yield l_ok and r_ok



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
    parser.add_argument('--exons', type=file)
    parser.add_argument('--exon-index')
    parser.add_argument('samfile')
    parser.add_argument('--log')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    parser.add_argument('--iterate-sam', '-s', default='exons', choices=['exons', 'reads'])

    args = parser.parse_args()

    level = logging.DEBUG if args.debug else logging.INFO 

    logging.basicConfig(level=level,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=args.log,
                        filemode='w')

    console = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)

    output = args.output if args.output is not None else sys.stdout

    if args.rum_gene_info is not None:

        genes = parse_rum_index_genes(args.rum_gene_info)
        exons = genes_to_exons(genes)
    elif 'exons' in args:
        exons = read_exons(args.exons)

    print('feature', 'min', 'max', sep='\t', file=output)
    for (exon, count_u, count_m) in iterate_over_exons(exons, args.samfile):
        exon_str = '{chr_}:{start}-{end}'.format(
            chr_=exon.ref,
            start=exon.location.start+1,
            end=exon.location.end)
        min_count = count_u
        if count_u > 0:
            max_count = count_u + count_m
        else:
            max_count = 0
        print(exon_str, min_count, max_count, sep='\t', file=output)

#    read_sam_file(args.rum_gene_info, args.samfile, args.output)

