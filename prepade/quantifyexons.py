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

def is_consistent(aln, exon):
    exon_start = exon.location.start
    exon_end   = exon.location.end
    strand = -1 if aln.is_reverse else 1
    spans = cigar_to_spans(aln.cigar, aln.pos, strand).sub_features

    for i, span in enumerate(spans):

        span_start = span.location.start 
        span_end   = span.location.end   

        if i == 0:
            consistent_start = span_start >= exon_start
        else:
            consistent_start = span_start == exon_start

        if i == len(spans) - 1:
            consistent_end = span_end <= exon_end
        else:
            consistent_end = span_end == exon_end

        if consistent_start and consistent_end:
            return True

    return False

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

def genes_to_exons(genes):
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)


def iterate_over_exons(exons, sam_filename):
    samfile = pysam.Samfile(sam_filename)
    details = open('details', 'w')

    seen = set()

    for exon in exons:

        key = exon.ref + ':' + str(exon.location.start) + '-' + str(exon.location.end)

        if key in seen:
            continue
        seen.add(key)

        unique_read_ids = set()
        multi_read_ids = set()

        all_spans = []


        for aln in samfile.fetch(exon.ref, exon.location.start, exon.location.end):
            if match(exon, [aln]):
                num_alns = aln.opt('IH')
                if num_alns > 1:
                    multi_read_ids.add(aln.qname)
                else:
                    unique_read_ids.add(aln.qname)

        count_u = len(unique_read_ids)
        count_m = len(multi_read_ids)

        if count_u > 0:
            yield(exon, count_u, count_m)
    details.close()


#            ============  
# o c         aaa   aaa    ?
# o c        aaaaaa        
# o c              aaaaaa  
# o c  aaa   aaa           Ok because non-first read span starts at exon start
# o i  aaa     aaa         Inconsistent because non-first read span starts after exon start
# o i

# n c  ---   ============
# o c        [  (    )  ]
# o c        [(      )  ]
# o c        [  (      )]
# n c   (  ) [          ]
# n c        [          ]  (  )
# o i     (  [    )     ]
# o i        [      (   ]  )

# if it's first span, it must start at or after exon start
# if it's not first span, it must start at exon start
# if it's last span, must end at or before exon end
# if it's not last span, it must end at exon end

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
        overlaps.extend(spans_overlap(exon, spans))
        consistents.extend(spans_are_consistent(exon, spans))

    any_overlap    = any(overlaps)
    all_consistent = all(consistents)

    return any_overlap and all_consistent
        
def spans_overlap(exon, spans):
    lexon = exon.location.start
    rexon = exon.location.end

    for span in spans:
        lspan = span.location.start
        rspan = span.location.end
        yield not (rspan < lexon or lspan > rexon) 


def spans_are_consistent(exon, spans):

    lexon = exon.location.start
    rexon = exon.location.end

    last_span = len(spans) - 1

    for i, span in enumerate(spans):
        lspan = span.location.start
        rspan = span.location.end
        
        is_first = i == 0
        is_last  = i == last_span

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
    parser.add_argument('--exon-index')
    parser.add_argument('samfile')

    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    args = parser.parse_args()

    output = args.output if args.output is not None else sys.stdout

    genes = parse_rum_index_genes(args.rum_gene_info)
    exons = genes_to_exons(genes)

    print('feature', 'min', 'max', sep='\t', file=output)
    for (exon, count_u, count_m) in iterate_over_exons(exons, args.samfile):
        exon_str = '{chr_}:{start}-{end}'.format(
            chr_=exon.ref,
            start=exon.location.start,
            end=exon.location.end)
        print(exon_str, count_u, count_u + count_m, sep='\t', file=output)

#    read_sam_file(args.rum_gene_info, args.samfile, args.output)

