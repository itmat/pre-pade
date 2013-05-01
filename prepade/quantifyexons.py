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
from prepade.clutils import UsageException, setup_logging
from prepade.geneio import (
    parse_rum_index_genes, read_exons, ExonIndex, genes_to_exons)
from prepade.samutils import (
    AlignmentFileType, input_file_ordering, has_hi_and_ih_tags, qname_and_hi, 
    sam_iter)


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

    logging.info("Iterating over alignments to accumulate counts")
    for pair in alns_by_qname_and_hi(mapped):

        pair = list(pair)
        ref = samfile.getrname(pair[0].tid)
        num_alns = pair[0].opt('IH')
        spans = []
        for aln in pair:
            spans.extend(cigar_to_spans(aln.cigar, aln.pos))

        overlapping = set()

        # Go through all the spans and add to 'overlapping' a key for
        # each exon that a span overlaps.
        for span in spans:
            for exon in idx.get_exons(ref, span.start, span.end):
                key = (exon.ref, exon.location.start, exon.location.end)
                    
                if matches(exon, pair):
                    overlapping.add(key)

        # Overlapping now has a key for each exon that any span from
        # this pair of alignments overlaps.
        for key in overlapping:
            if num_alns > 1:
                multi_counts[key] += 1
            else:
                unique_counts[key] += 1

    logging.info("Done accumulating counts")

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
    
        for (qname, hi), pair in groupby(alns, key=qname_and_hi):

            pair = list(pair)
            if matches(exon, pair):
                num_alns = pair[0].opt('IH')
                if num_alns > 1:
                    multi_read_ids.add((qname, hi))
                else:
                    unique_read_ids.add((qname, hi))

        count_u = len(unique_read_ids)
        count_m = len(multi_read_ids)

        yield(exon, count_u, count_m)


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

    overlap    = []
    consistent = []

    span_groups = [ cigar_to_spans(aln.cigar, aln.pos) for aln in alns ]

    for spans in span_groups:
        my_overlap = spans_overlap(exon.location, spans)
        overlap.extend(my_overlap)

    if not any(overlap):
        return False

    for spans in span_groups:
        consistent.extend(spans_are_consistent(exon.location, spans))

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



def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

Generates quantifications, given a list of exons and a set of
alignments. The alignments can either be given as (1) a SAM or BAM
file where all alignments for the same read are on consecutive lines,
or as (2) an indexed BAM file (which would have to be sorted by
chromosomal position in order to be indexed). The HI and IH tags must
be set properly in the file, regardless of the input format. You can
give the list of exons in one of two ways: (1) a file with one exon on
each line, in the format CHR:START-END, or (2) the 'gene info' file
from a RUM index.""")

    parser.add_argument('--rum-gene-info', 
                       action='store_true',
                       help="Indicate that exons is a 'gene info' file from a RUM index")

    parser.add_argument('-d', '--debug', action='store_true',
                        help="Turn on debug-level logging")

    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        help="Location of output file")

    parser.add_argument('-l', '--log',
                        help="Write log messages here; defaults so sys.stderr")

    parser.add_argument('exons')
    parser.add_argument('alignments')

    args = parser.parse_args()
    setup_logging(args)
    output = get_output_fh(args)

    logging.info("Checking alignment file")

    ordering = require_sam_ordering_and_hi_tags(args.alignments)

    if args.rum_gene_info:
        with open(args.exons) as infile:
            genes = parse_rum_index_genes(infile)
            exons = list(genes_to_exons(genes))

    else:
        with open(args.exons) as infile:
            exons = list(read_exons(infile))
                
    if ordering == AlignmentFileType.ORDERED_BY_READ_NAME:
        counts = iterate_over_sam(exons, args.alignments)
    else:
        counts = iterate_over_exons(exons, args.alignments)

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
    try:
        main()
    except UsageException as e:
        print(e, file=sys.stderr)
