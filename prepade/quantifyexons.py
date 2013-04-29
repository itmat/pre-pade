#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import logging
import numpy as np
import pandas as pd
import pysam
import sys

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from collections import defaultdict
from itertools import groupby, ifilter
from prepade.geneio import parse_rum_index_genes, read_exons, ExonIndex

class CigarOp:
    """Constants representing the CIGAR operations."""

    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    EQUAL = 7
    X = 8


def genes_to_exons(genes):
    """Given an iterator over genes, returns an iterator over exons.

    """
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)


def qname_and_hi(aln):
    """Returns a tuple of the qname and 'HI' tag for the given AlignedRead."""
    return (aln.qname, aln.opt('HI'))


def sam_iter(samfile):
    while True:
        rec = samfile.next()
        yield rec

def iterate_over_sam(exons, sam_filename):
    """Return exon quantifications by iterating over SAM file.

    Reads all the exons into an index in memory. Then iterates over
    all the aligments in order in the given SAM file, accumulating
    counts for each exon.

    :param exons: 
      Sequence of SeqFeature objects representing exons.

    :param sam_filename:
      Identifies the SAM file. Reads that are part of the same pair
      must be next to each other in the file, and the HI and IH tags
      must be used.

    :return:
      An iterator over (exon, unique_count, non_unique_count) tuples,
      where exon is a SeqFeature object, unique count is the number of
      unique mappers that overlap it consistently, and
      non_unique_count is the number of non-unique mappers that
      overlap it consistently
    
    """
    logging.info('Building exon index')
    idx = ExonIndex(exons)
    samfile = pysam.Samfile(sam_filename)
    unique_counts = defaultdict(lambda: 0)
    multi_counts = defaultdict(lambda: 0)

    mapped = ifilter(lambda x: not x.is_unmapped, sam_iter(samfile))

    for key, pair in groupby(mapped, key=qname_and_hi):

        pair = list(pair)
        ref = samfile.getrname(pair[0].tid)
        num_alns = pair[0].opt('IH')
        spans = []
        for aln in pair:
            spans.extend(cigar_to_spans(aln.cigar, aln.pos))

        overlapping = set()

        for span in spans:
            for exon in idx.get_exons(ref, span.start, span.end):
                key = (exon.ref, 
                       exon.location.start,
                       exon.location.end)
                    
                if matches(exon, pair):
                    overlapping.add(key)

        for key in overlapping:
            if num_alns > 1:
                multi_counts[key] += 1
            else:
                unique_counts[key] += 1

    for key in unique_counts:
        (ref, start, end) = key
        exon = SeqFeature(ref=ref, location=FeatureLocation(start, end))
        yield(exon, unique_counts[key], multi_counts[key])

def iterate_over_exons(exons, bam_filename):
    """Return exon quantifications by iterating over the exons file.

    Iterates over all the exons in the given sequence, and for each
    one, fetches all overlapping reads from the specified BAM file.

    :param exons: 
      Sequence of SeqFeature objects representing exons.

    :param sam_filename:
      Identifies the BAM file. It must be an indexed BAM file, and the
      HI and IH tags must be used.

    :return:
      An iterator over (exon, unique_count, non_unique_count) tuples,
      where exon is a SeqFeature object, unique count is the number of
      unique mappers that overlap it consistently, and
      non_unique_count is the number of non-unique mappers that
      overlap it consistently

    """
    samfile = pysam.Samfile(bam_filename)

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

        alns = sorted(alns, key=qname_and_hi)

        if __debug__:
            logging.debug('Got %d candidate alignments', len(alns))

        for (qname, hi), pair in groupby(alns, key=qname_and_hi):
            if __debug__:
                logging.debug(' ')
                logging.debug('  candidate %s[%d]', qname, hi)
            pair = list(pair)
            if matches(exon, pair):
                if __debug__:
                    logging.debug('    found a match!')
                num_alns = pair[0].opt('IH')
                if num_alns > 1:
                    multi_read_ids.add((qname, hi))
                else:
                    unique_read_ids.add((qname, hi))

        count_u = len(unique_read_ids)
        count_m = len(multi_read_ids)

        yield(exon, count_u, count_m)


def cigar_to_spans(cigar, start):
    """Converts a CIGAR data structure and position to list of FeatureLocations.

    :param cigar:
      List of tuples like what is returned by the 'cigar' property
      from a Pysam alignment.

    :param start:
       Start location, using 0-based indexing.

    :return:

      List of FeatureLocation objects representing the spans of the
      reference that this segment matches.
      
    """
    spans = []
    if __debug__:
        logging.debug('cigar is %s, start is %d', cigar, start)

    if cigar is None:
        return []

    cigar = remove_ds(cigar)

    for (op, bases) in cigar:

        if op == CigarOp.M:
            end = start + bases
            spans.append(FeatureLocation(start, end))
            start = end

        elif op == CigarOp.D or op == CigarOp.N:
            start = start + bases

    res = []

    for span in spans:
        if len(res) > 0 and res[-1].end + 1>= span.start:
            start = res[-1].start
            end   = span.end
            res[-1] = FeatureLocation(start, end)
        else:
            res.append(span)

    return res

def format_spans(spans):
    """Format a list of FeatureLocations representing spans as a string."""
    return ', '.join([ '%s-%s' % (s.start, s.end) for s in spans ])
    
def matches(exon, alns):
    """Returns True if the given alignments match the given exon.

    At least one segment from the alignments must overlap the exon,
    and all of the segments must be 'consistent' with the exon,
    meaning that the read must not cross a junction, and a gap in the
    read cannot occur in the middle of the exon.

    :param exon:
       A SeqFeature representing the exon.

    :param alns:
      List of alignments to check. There should typically be either
       one or two alns; one in the case of single-end reads or a read
       where only one of the pair aligned, and two where both reads
       aligned.

    :return:
       Boolean indicating whether the fragment matches.

    """
    if __debug__:
        logging.debug('    exon is %d-%d', exon.location.start, exon.location.end)

    overlap    = []
    consistent = []

    span_groups = [ cigar_to_spans(aln.cigar, aln.pos) for aln in alns ]

    for spans in span_groups:
        my_overlap = spans_overlap(exon.location, spans)
        overlap.extend(my_overlap)

        if __debug__:
            logging.debug('    spans are %s, overlap is %s', 
                          format_spans(spans), my_overlap)

    if not any(overlap):
        return False

    for spans in span_groups:
        consistent.extend(spans_are_consistent(exon.location, spans))

    if __debug__:
        logging.debug('  consistent is %s', consistent)

    return all(consistent)

        
def spans_overlap(exon, spans):
    """Return bools indicating which spans overlap the exon.

    :param exon:
      A FeatureLocation representing an exon.

    :param spans:
      A list of FeatureLocations, each representing a read segment.

    :return:

      Array of booleans, one for each read segment, indicating which
      segments overlap the exon.

    """
    for span in spans:
        yield not (span.end <= exon.start or span.start >= exon.end) 


def spans_are_consistent(exon, spans):

    """Return bools indicating which read segments are consistent with exon.

    :param exon:
      A FeatureLocation representing an exon.

    :param spans:
      A list of FeatureLocations, each representing a read segment.

    :return:

      Array of booleans, one for each read segment. Any segment that
      spans a junction at either edge of the exon is considered
      inconsistent. Any segment other than the first that starts in
      the middle of the exon is inconsistent, and likewise any segment
      other than the last that ends in the middle of the exon is
      considered inconsistent. All other exons are consistent.
    
    """
    last_span = len(spans) - 1

    for i, span in enumerate(spans):
        
        if span.end <= exon.start or span.start >= exon.end:
            yield True

        else:
            if i == 0:
                l_ok = span.start >= exon.start
            else:
                l_ok = span.start == exon.start

            if i == last_span:
                r_ok = span.end <= exon.end
            else:
                r_ok = span.end == exon.end

            yield l_ok and r_ok


def remove_ds(cigar):
    """Removes D operations from CigarOp string, replacing with Ms.

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

    d_to_m = lambda (op, bases): (CigarOp.M, bases) if op == CigarOp.D else (op, bases)
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

def main():

    parser = argparse.ArgumentParser(
        description="""
        
        This program expects two types of input: a set of alignments
        which can be in SAM or BAM format, and a list of exons, which
        can be in a variety of formats.

        Alignments
        ----------

        The alignments can either be given as a SAM file, a BAM file,
        or an indeed BAM file. If the input is an indexed BAM file, we
        will take advantage of the index by iterating over the list of
        exons and fetching any overlapping reads from the BAM file. If
        it is a non-indexed BAM file or a SAM file, or if you give the
        --ignore-index option, we will instead load all the exons into
        memory, and then iterate over all the alignments to accumulate
        the counts.

        The HI and IH tags must be set properly in the file. Assuming
        those tags are set, any indexed BAM file should be
        acceptable. If the input is given as a SAM file or a BAM file
        without an index, then the forward and reverse reads from the
        same fragment must be adjacent to each other in the file.

        Note that using an indexed BAM file may actually be slower
        than iterating over all the reads in the SAM file. 

        """)
    input_grp = parser.add_mutually_exclusive_group(required=True)

    input_grp.add_argument('--rum-gene-info', 
                       type=file,
                       help="Input as a gene file from a RUM index")

    input_grp.add_argument('--exons',
                       help="Input as a list of exons, one on each line, like 'chr1:234-678'",
                       type=file)

    parser.add_argument('samfile')

    parser.add_argument('--debug', '-d', action='store_true',
                        help="Turn on debug-level logging")
    parser.add_argument('--output', '-o', 
                        type=argparse.FileType('w'),
                        help="Location of output file")
    parser.add_argument('--ignore-index', action='store_true',
                        help="""By default, we expect the alignments to be given as an indexed BAM file. We iterate over the exons and select the overlapping reads from the BAM file using the index. Giving the --iterate-sam flag causes us to instead read in all the exons into memory and then iterate over a SAM file.""")


    args = parser.parse_args()

    level = logging.DEBUG if args.debug else logging.INFO 

    logging.basicConfig(level=level,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
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

    if args.iterate_sam:
        counts = iterate_over_sam(exons, args.samfile)

    else:
        counts = iterate_over_exons(exons, args.samfile)
            
    print('feature', 'min', 'max', sep='\t', file=output)
    for (exon, count_u, count_m) in counts:
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


if __name__ == '__main__':
    main()
