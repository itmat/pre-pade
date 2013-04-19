#!/usr/bin/env python

from __future__ import print_function

import sys
import pysam
import argparse

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from collections import defaultdict
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rum-gene-info')
    parser.add_argument('samfile')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    args = parser.parse_args()

    read_sam_file(args.rum_gene_info, args.samfile, args.output)
    
