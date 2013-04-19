#!/usr/bin/env python

from __future__ import print_function

import argparse

from prepade.geneio import parse_rum_index_genes

def extract_exons(in_fh, out_fh):

    for (chr_, gene) in parse_rum_index_genes(in_fh):
        for exon in gene.sub_features:
            loc = exon.location
            print("{0}:{1}-{2}".format(chr_, loc.start, loc.end), 
                  file=out_fh)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
Reads the gene annotation file from a RUM index and produces a file
listing exon locations.""")

    parser.add_argument('--rum-gene-info', type=file, required=True,
                        help="Gene annotation file from RUM index")
    parser.add_argument('--output', '-o', type=argparse.FileType('w'),
                        help="Output file, defaults to STDOUT")
    args = parser.parse_args()

    extract_exons(args.rum_gene_info, args.output)
