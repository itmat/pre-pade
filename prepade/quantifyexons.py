#!/usr/bin/env python

from __future__ import print_function

import sys
import pysam
import argparse
import re

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from collections import defaultdict, namedtuple
from prepade.geneio import parse_rum_index_genes

def read_chr_to_genes(fh):
    res = defaultdict(list)
    for (chr_, gene) in parse_rum_index_genes(fh):
        res[chr_].append(gene)
    return res

def read_sam_file(gene_filename, sam_filename, output_fh):

    samfile = pysam.Samfile(sam_filename)

    with open(gene_filename) as fh:
        for (chr_, gene) in parse_rum_index_genes(fh):
            exons = gene.sub_features
            for exon in exons:
                start = exon.location.start
                end   = exon.location.end
                count = 0
                for rec in samfile.fetch(chr_, start, end):
                    count += 1

                if count > 0:
                    print(chr_, start, end, count, sep="\t", file=output_fh)

CigarElem = namedtuple('CigarElem', ['count', 'op'])

def parse_cigar(cigar):
    pat = re.compile("\d+[A-Z]")

    for op in pat.findall(cigar):
        yield CigarElem(int(op[:-1]), op[-1])

def cigar_to_spans(cigar, start):

    spans = []

    for x in parse_cigar(cigar):
        if x.op == 'M':
            end = start + x.count - 1
            spans.append(FeatureLocation(start, end))
            start = end

        elif x.op in 'DN':
            start = start + x.count + 1

        elif x.op == 'I':
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

    >>> remove_ds('21M1D54M')
    '76M'

    >>> remove_ds('8S4M1D63M')
    '8S68M'

    >>> remove_ds('15M1D15M2D29M2D16M')
    '80M'

    >>> remove_ds('21M1D41M177N13M')
    '63M177N13M'

    >>> remove_ds('41M354N20M1D14M')
    '41M354N35M'

    >>> remove_ds('4S26M1D45M')
    '4S72M'

    """
    d_to_m = lambda x: CigarElem(x.count, 'M') if x.op == 'D' else x
    parsed = parse_cigar(cigar)
    converted = map(d_to_m, parsed)

    res = []

    for x in converted:
        if x.op == 'M' and len(res) > 0 and res[-1].op == 'M':
            res[-1] = CigarElem(res[-1].count + x.count, 'M')
        else:
            res.append(x)

    return ''.join([str(x.count) + x.op for x in res])
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rum-gene-info')
    parser.add_argument('samfile')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    args = parser.parse_args()

    read_sam_file(args.rum_gene_info, args.samfile, args.output)
    
