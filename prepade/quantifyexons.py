#!/usr/bin/env python

from __future__ import print_function, division

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

def process_sam_rec(rec):
    print("Cigar: " + str(rec.cigar))
#    cigar = remove_ds(rec.cigarstring)
#    spans = cigar_to_spans(cigar, rec.pos)



def read_sam_file(gene_filename, sam_filename, output_fh):

    samfile = pysam.Samfile(sam_filename)

    with open(gene_filename) as fh:
        
        for exon in exons:
            start = exon.location.start
            end   = exon.location.end
            count = 0
            for rec in samfile.fetch(chr_, start, end):
                process_sam_rec(rec)
                count += 1

            if count > 0:
                print(chr_, start, end, count, sep="\t", file=output_fh)

CigarElem = namedtuple('CigarElem', ['count', 'op'])

class FeatureIndex(object):
    
    def __init__(self, features):
        self.features = sorted(features)

        
def bsearch_boundary(items, key, f=cmp):
    
    """Returns the first i for which f(items[i - 1]) < key and f(items[i])
    >= key. An alternate comparison function can be given as f.

    >>> bsearch_boundary('abcdefg', 'c')
    2

    >>> bsearch_boundary('cdefg', 'a')
    0
    
    >>> bsearch_boundary('abcde', 'g')
    5

    >>> bsearch_boundary('abfg', 'e')
    2

    """

    p = 0
    q = len(items)


    while True:

        r = p + ((q - p) // 2)

        # Find the smallest i where item i-1's end is less than my end
        # and item i's end is greater than or equal to than my end.

        if r == 0 or f(items[r - 1], key) < 0:

            #            key
            # [ ... prev this ...]
            if f(items[r], key) >= 0:
                return r

            #                   key
            # [ ... prev this ] 
            elif r == len(items) - 1:
                return len(items)

            #                  key
            # [ ... prev this ...  ]
            else:
                p = r

        # key
        #    last
        else:
            q = r

def bsearch_range(items, key, f_left, f_right):
    p = bsearch(items, key, f_left)
    q = bsearch(items, key, f_right)
    return items[p:q]

    

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
    parser.add_argument('--rum-gene-info', type=file)
    parser.add_argument('samfile')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    args = parser.parse_args()

    chr_to_genes = read_chr_to_genes(args.rum_gene_info)



    read_sam_file(args.rum_gene_info, args.samfile, args.output)
    
