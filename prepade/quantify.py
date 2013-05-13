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
from collections import defaultdict, namedtuple
from itertools import groupby, ifilter, chain
from prepade.clutils import (
    UsageException, setup_logging, get_output_fh, 
    require_sam_ordering_and_hi_tags)
from prepade.geneio import (
    parse_rum_index_genes, read_exons, ExonIndex, TranscriptIndex, 
    genes_to_exons, parse_gtf_to_genes)
from prepade.samutils import (
    AlignmentFileType, input_file_ordering, has_hi_and_ih_tags, qname_and_hi, 
    sam_iter, alns_by_qname_and_hi, cigar_to_spans, spans_for_aln, aln_to_span_ndarray)

def guess_input_file_type(args):
    
    if args.model.endswith('.gtf') or args.model.endswith('.GTF'):
        return 'gtf'

    elif args.rum_gene_info:
        return 'rum_gene_info'
    
    else:
        return 'exon_list'

class FeatureReadCounter(object):

    def count_towards_min(self, aln):
        try:
            num_alns = aln.opt('IH')
        except KeyError:
            num_alns = aln.opt('NH')
            
        return num_alns == 1


class ExonReadCounter(FeatureReadCounter):

    def __init__(self, exon_index, min_overlap):
        self.exon_index = exon_index
        self.unique_counts = defaultdict(lambda: 0)
        self.multi_counts = defaultdict(lambda: 0)
        self.min_overlap = min_overlap

    def add(self, ref, alns):
        pair = alns

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
            for exon in self.exon_index.get_exons(ref, span.start, span.end):
                key = (exon.id, exon.ref, exon.location.start, exon.location.end)

                if key not in seen:
                    if matches_exon(exon, pair, self.min_overlap):
                        if self.count_towards_min(alns[0]):
                            self.unique_counts[key] += 1                        
                        else:
                            self.multi_counts[key] += 1

                        seen.add(key)
        

    def __iter__(self):
        for key in self.unique_counts:
            (exon_id, ref, start, end) = key
            exon = SeqFeature(id=exon_id, ref=ref, location=FeatureLocation(start, end))
            yield(exon, self.unique_counts[key], self.multi_counts[key])        

class TranscriptReadCounter(FeatureReadCounter):


    def __init__(self, index):
        self.index = index
        self.unique_counts = defaultdict(lambda: 0)
        self.multi_counts = defaultdict(lambda: 0)        
        self.key_to_transcript = {}
        
    def add(self, ref, alns):
        pair = alns
        
        span_groups = [ aln_to_span_ndarray(aln) for aln in alns ]

        seen = set()

        for spans in span_groups:

            for span in spans:
                (start, end) = span
                for t in self.index.get_transcripts(ref, start, end):
                    key = t.id
                    self.key_to_transcript[key] = t
                    if key not in seen:
                        seen.add(key)
                        (hit, details) = compare_alns_to_transcript(t, span_groups)

                        if hit:
                            if self.count_towards_min(alns[0]):
                                self.unique_counts[key] += 1                        
                            else:
                                self.multi_counts[key] += 1

        
    def __iter__(self):
        for key in self.unique_counts:
            yield(self.key_to_transcript[key], self.unique_counts[key], self.multi_counts[key])

def iterate_over_sam(sam_filename, counters):

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

    samfile = pysam.Samfile(sam_filename)

    logging.info("Iterating over alignments to accumulate counts")


    for i, pair in enumerate(alns_by_qname_and_hi(sam_iter(samfile))):

        pair = list(pair)
        ref = samfile.getrname(pair[0].tid)
        for counter in counters:
            counter.add(ref, pair)
        if ((i + 1) % 10000) == 0:
            logging.info("Done " + str(i + 1) + " alignments")

    logging.info("Done accumulating counts")

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
                try:
                    num_alns = pair[0].opt('IH')
                except KeyError:
                    num_alns = pair[0].opt('NH')
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

def first_intersection(introns, spans):
    m = len(introns)
    n = len(spans)
    i = 0
    j = 0
    while i < m and j < n:
        intron_start = introns[i, 0]
        intron_end   = introns[i, 1]
        span_start   = spans[j, 0]
        span_end     = spans[j, 1]

        # If the span ends before the intron starts, advance to next
        # span
        if spans[j, 1] <= introns[i, 0]:
            j += 1
            
        # If the intron ends before the span starts, advance to next
        # intron
        elif introns[i, 1] <= spans[j, 0]:
            i += 1

        else:
            return i

    return None


def exons_intersect(target, spans):
    start = target[0]
    end   = target[1]

    for s, e in spans:
        if not (e <= start or s >= end):
            return True

    return False


TranscriptMatch = namedtuple(
    'TranscriptMatch',
    ['decision',
     'first_exon_hit',
     'spans',
     'gaps',
     'exons', 
     'introns',
     'intron_hits'])

def spans_to_gaps(spans):
    res = np.zeros((len(spans) - 1, 2))
    res[:, 0] = spans[:-1, 1]
    res[:, 1] = spans[1:, 0]
    return res

def compare_alns_to_transcript(transcript, alns):
    matches = [compare_aln_to_transcript(transcript, aln)
               for aln in alns]

    decisions = [ m.decision for m in matches ]
    decision = any(decisions) and not any(d == False for d in decisions)
    return (decision, matches)

def compare_aln_to_transcript(transcript, spans):

    if isinstance(transcript, SeqFeature):
        exons = np.zeros((len(transcript.sub_features), 2), int)
        exons[:, 0] = [e.location.start for e in transcript.sub_features]
        exons[:, 1]  = [e.location.end   for e in transcript.sub_features]        

        return compare_aln_to_transcript(exons, spans)

    elif not isinstance(transcript, np.ndarray):
        return compare_aln_to_transcript(np.array(transcript, int), spans)
    
    if isinstance(spans[0], FeatureLocation):
        tmp = spans
        spans = np.zeros((len(tmp), 2), int)
        spans[:, 0] = [s.start for s in tmp]
        spans[:, 1] = [s.end   for s in tmp]        
        return compare_aln_to_transcript(transcript, spans)

    elif not isinstance(spans, np.ndarray):
        return compare_aln_to_transcript(transcript, np.array(spans, int))

    exons = transcript

    introns = spans_to_gaps(exons)
    gaps    = spans_to_gaps(spans)

    # If we hit any introns we have to call it a miss
    first_intron_hit = first_intersection(introns, spans)

    # If the read doesn't overlap any exons, we can't call it a
    # confirmation. However we won't call it a rejection either if
    # there are no gaps in the read (which would suggest an intron not
    # present in this transcript) and if the read does not fall in an
    # intron region.

    if first_intron_hit is not None:
        decision = False
        first_exon_hit = None

    else:
        first_exon_hit = first_intersection(exons, spans)

        if first_exon_hit is None:
            if len(gaps) > 0:
                decision = False
            else:
                decision = None
        else:
            num_gaps = len(gaps)
            decision = (first_exon_hit + num_gaps <= len(introns) and 
                        np.all(gaps == introns[first_exon_hit : first_exon_hit + len(gaps)]))
        
    return TranscriptMatch(decision=decision, 
                           first_exon_hit=first_exon_hit,
                           spans=spans, gaps=gaps, exons=exons, introns=introns, intron_hits=(first_intron_hit is not None))
 
def matches_exon(exon, alns, min_overlap=1):
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
        overlap.extend(spans_overlap(exon.location, spans, min_overlap))
    if not any(overlap):
        return False

    # Now make sure that all of the spans are "consistent" with the
    # exon. We check each group separately, because internal read
    # segments are treated differently from the first and last
    # segments of each read.

    for spans in span_groups:
        consistent.extend(spans_are_consistent(exon.location, spans))

    return all(consistent)

        
def spans_overlap(exon, spans, min_overlap=1):
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
        yield not (span.end < exon.start + min_overlap or span.start > exon.end - min_overlap) 


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

class BedFileWriter(object):

    def __init__(self, fh):

        self.fh = fh
        print('featureType', 'chrom', 'chromStart', 'chromEnd', 'name', 'min_count', 'max_count',
              'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
              'blockSizes', 'blockStarts', sep='\t', file=self.fh)

    def write(self, feature_type, chrom, chromStart, chromEnd, name, min_count, max_count, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts):
        print(feature_type, chrom, chromStart, chromEnd, name, min_count, max_count, strand, 
              thickStart, thickEnd, itemRgb, blockCount, blockSizes, 
              blockStarts, sep='\t', file=self.fh)        

    def write_exon(self, exon, count_u, count_m):
        chrom = exon.ref
        chromStart = exon.location.start
        chromEnd   = exon.location.end
        name = exon.id

        min_count = count_u
        if count_u > 0:
            max_count = count_u + count_m
        else:
            max_count = 0

        strand = exon.strand if exon.strand is not None else ''
        thickStart = ''
        thickEnd = ''
        itemRgb = '0,0,0'
        blockCount= 0
        blockSizes = ''
        blockStarts = ''

        self.write('exon', chrom, chromStart, chromEnd, name, min_count, max_count, strand, 
              thickStart, thickEnd, itemRgb, blockCount, blockSizes, 
              blockStarts)

    def write_transcript(self, transcript, count_u, count_m):
        chrom = transcript.ref
        chromStart = transcript.location.start
        chromEnd   = transcript.location.end
        name = transcript.id

        min_count = count_u
        if count_u > 0:
            max_count = count_u + count_m
        else:
            max_count = 0

        strand = transcript.strand if transcript.strand is not None else ''
        thickStart = ''
        thickEnd = ''
        itemRgb = '0,0,0'
        exon_locs = [ sf.location for sf in transcript.sub_features ]
        blockCount = len(exon_locs)
        blockSizes  = [ el.end - el.start for el in exon_locs ]
        blockStarts = [ el.start - chromStart for el in exon_locs ]

        blockSizes = ','.join(map(str, blockSizes))
        blockStarts = ','.join(map(str, blockStarts))

        self.write('transcript', chrom, chromStart, chromEnd, name, min_count, max_count, strand, 
              thickStart, thickEnd, itemRgb, blockCount, blockSizes, 
              blockStarts)


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

Generates quantifications, given a set of alignments and either a
transcript model or a list of exons. The alignments can either be
given as (1) a SAM or BAM file where all alignments for the same read
are on consecutive lines, or as (2) an indexed BAM file (which would
have to be sorted by chromosomal position in order to be indexed). The
HI and either IH or NH tags must be set properly in the file,
regardless of the input format. You can supply a transcript model in
GTF format, or a "gene info" file from a RUM index. If you want exons
and not transcripts, you can just give a file that has one exon on
each line, in the format CHR:START-END.

The output is given in a format similar to BED, but with a couple
extra columns that contain some additional information. It should be
easy to transform the output into BED format with a few simple UNIX
commands.""")

    parser.add_argument('--rum-gene-info', 
                       action='store_true',
                       help="Indicate that exons is a 'gene info' file from a RUM index")

    parser.add_argument('-d', '--debug', action='store_true',
                        help="Turn on debug-level logging")

    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        help="Location of the output file; defaults to sys.stdout")

    parser.add_argument('-l', '--log',
                        help="Write log messages here; defaults so sys.stderr")

    parser.add_argument('model',
                        help="The transcript model or list of exons")
    parser.add_argument('alignments',
                        help="SAM or BAM file containing alignments")

    parser.add_argument('--no-transcripts',
                        action='store_true',
                        help="Don't do transcripts, just exons")

    parser.add_argument('--min-overlap',
                        type=int,
                        help="Minimum overlap required to count an exon")

    args = parser.parse_args()
    setup_logging(args)

    logging.info("Checking alignment file")

    ordering = require_sam_ordering_and_hi_tags(args.alignments)

    input_file_type = guess_input_file_type(args)

    # Parse the exons and optionally transcripts from the model file.
    with open(args.model) as infile:

        if input_file_type == 'rum_gene_info':
            logging.info("Loading transcript model from RUM index at " + args.model)
            genes = list(parse_rum_index_genes(infile))
            exons = list(genes_to_exons(genes))
    
        elif input_file_type == 'gtf':
            logging.info("Loading transcript model from GTF file at " + args.model)
            genes = list(parse_gtf_to_genes(infile))
            exons = list(genes_to_exons(genes))
        
        else:
            logging.info("Loading exons from " + args.model)
            genes = None
            exons = list(read_exons(infile))

    if ordering == AlignmentFileType.ORDERED_BY_READ_NAME:

        if genes is None:
            logging.info('Building exon index')
            idx = ExonIndex(exons)
        else:
            logging.info('Building transcript index')
            idx = TranscriptIndex(genes)

        logging.info("For exon quantification, only counting reads that overlap the exon by " + str(args.min_overlap) + " bases or more")

        exon_counter = ExonReadCounter(idx, args.min_overlap)
        transcript_counter = None
        counters = [ exon_counter ]
        if genes is not None:
            if args.no_transcripts:
                logging.info("Ignoring transcripts, just doing exons")
            else:
                transcript_counter = TranscriptReadCounter(idx)
                counters.append(transcript_counter)
            
        iterate_over_sam(args.alignments, counters)

    else:
        exon_counter = iterate_over_exons(exons, args.alignments)
        if genes is not None:
            logging.warn(
                "I currently can't do transcript quantification for indexed " +
                "BAM files; only for SAM or BAM files sorted by read name. " +
                "I will just output exon quantifications for this input file.")
        transcript_counter = None

    output = get_output_fh(args)
    writer = BedFileWriter(output)
    
    for (exon, count_u, count_m) in exon_counter:
        writer.write_exon(exon, count_u, count_m)

    if transcript_counter is not None:
        for (transcript, count_u, count_m) in transcript_counter:
            writer.write_transcript(transcript, count_u, count_m)        

if __name__ == '__main__':
    try:
        main()
    except UsageException as e:
        print(e, file=sys.stderr)

