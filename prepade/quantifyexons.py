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
from itertools import groupby, ifilter, chain
from prepade.clutils import (
    UsageException, setup_logging, get_output_fh, 
    require_sam_ordering_and_hi_tags)
from prepade.geneio import (
    parse_rum_index_genes, read_exons, ExonIndex, genes_to_exons)
from prepade.samutils import (
    AlignmentFileType, input_file_ordering, has_hi_and_ih_tags, qname_and_hi, 
    sam_iter, alns_by_qname_and_hi, cigar_to_spans, spans_for_aln)


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

    logging.info("Iterating over alignments to accumulate counts")
    for pair in alns_by_qname_and_hi(sam_iter(samfile)):

        pair = list(pair)
        ref = samfile.getrname(pair[0].tid)
        num_alns = pair[0].opt('IH')
        spans = []
        for aln in pair:
            spans.extend(cigar_to_spans(aln.cigar, aln.pos))

        seen = set()

        # Go through all the spans and find all the exons overlapped
        # by any span. For each of those exons, check to see if it
        # 'matches'. If it does, increment one of the counts (unique
        # or multi depending on how many alignments there are for this
        # read). Since we might come across the same exon more than
        # once, keep a set of 'seen' exons so we don't double-count.
        for span in spans:
            for exon in idx.get_exons(ref, span.start, span.end):
                key = (exon.ref, exon.location.start, exon.location.end)

                if key not in seen:
                    if matches_exon(exon, pair):
                        if num_alns > 1:
                            multi_counts[key] += 1
                        else:
                            unique_counts[key] += 1                        
                        seen.add(key)

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
            if matches_exon(exon, pair):
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

def introns_for_exons(exon_starts, exon_ends):
    return (exon_starts[1:], exon_ends[:-1])

def spans_intersect(start, end,
                    span_starts, span_ends):

    ends_before = span_ends <= start
    starts_after = span_starts >= end

    return np.logical_not(
        np.logical_or(starts_before, ends_after))


def compare_aln_to_transcript(transcript, aln):

    spans = cigar_to_spans(aln.cigar, aln.pos)

    exons = transcript.sub_features

    n = len(exons)

    overlap    = np.zeros(len(exons), bool)
    consistent = np.ones(len(exons), bool)

    span_starts = np.array([s.start for s in spans], int)
    span_ends   = np.array([s.end for s in spans], int)

    exon_starts = np.array([e.location.start for e in transcript.sub_features], int)
    exon_ends   = np.array([e.location.end   for e in transcript.sub_features], int)

    intron_starts = exon_ends[:-1]
    intron_ends   = exon_starts[1:]

    gap_starts = span_ends[:-1]
    gap_ends   = span_starts[1:]

    # Find the first exon that any of my spans intersect. If none of
    # them intersect the exon, then we can't call this a hit.
    first_transcribed_exon = None
    for i in range(n):

        hit = np.any(spans_intersect(exon_starts[i], exon_ends[i], span_starts, span_ends))
        
        exon_hits[i] = hit
        if hit:
            last_transcribed_exon = i
            if first_transcribed_exon is None:
                first_transcribed_exon = i

    # Starting at the first exon we hit, any gaps in the read must align exactly
    # with the introns.

    spanned_intron_starts = intron_starts[first_transcribed_exon : last_transcribed_exon]
    spanned_intron_ends = intron_ends[first_transcribed_exon : last_transcribed_exon]

    return (first_transcribed_exon is not None and 
            gap_starts == spanned_intron_starts and
            gap_ends   == spanned_intron_ends)


def matches_exon(exon, alns):
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

    overlap = []
    consistent = []

    # Get a list of lists of spans; one list of spans for each aligned
    # read for this fragment.
    span_groups = map(spans_for_aln, alns)
        
    # Make sure at least one of the spans overlaps the exon
    for spans in span_groups:
        overlap.extend(spans_overlap(exon.location, spans))
    if not any(overlap):
        return False

    # Now make sure that all of the spans are "consistent" with the
    # exon. We check each group separately, because internal read
    # segments are treated differently from the first and last
    # segments of each read.

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

        # If the span doesn't overlap the exon at all, it certainly
        # isn't inconsistent with it.
        if span.end <= exon.start or span.start >= exon.end:
            yield True

        else:
            # If it's the first span for the read, it can start either
            # at or after the exon start. If it's not the first span,
            # then it must start at the exon start, indicating that an
            # intron was skipped.
            if i == 0:
                l_ok = span.start >= exon.start
            else:
                l_ok = span.start == exon.start

            # And the opposite is true for right edge of the span
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
