from Bio.SeqFeature import FeatureLocation, SeqFeature
from collections import namedtuple, defaultdict
from itertools import groupby

import logging
import numpy as np
import pandas as pd
import re

def parse_gtf_attributes(attr_string):
    attrs = {}
    attr_pat = re.compile('(gene_id|transcript_id|exon_number)[= ]+(.*)')

    attrs = {}

    for field in attr_string.split(';'):
        field = field.strip()
        m = attr_pat.match(field)
        if m is not None:
            (name, value) = m.groups()
            if value[0] == '"' and value[-1] == '"':
                value = value[1:-1]
            attrs[name] = value
    return attrs
                
def parse_gtf_to_genes(fh):

    names=['seqname', 'source', 'feature', 'start', 
           'end', 'score', 'strand', 'frame', 'attribute']

    seqnames       = []
    starts         = []
    ends           = []
    strands        = []
    gene_ids       = []
    transcript_ids = []
    exon_numbers   = []

    for line in fh:
        line = line.rstrip()

        fields = line.split('\t', 9)

        (seqname, source, feature, start, 
         end, score, strand, frame, attr_string) = fields

        attrs = parse_gtf_attributes(attr_string)
        if strand == '+':
            strand = 1
        else:
            strand = -1

        if 'exon' in feature:
            seqnames.append(intern(seqname))
            starts.append(int(start) - 1)
            ends.append(int(end))
            strands.append(strand)
            gene_ids.append(intern(attrs['gene_id']))
            transcript_ids.append(intern(attrs['transcript_id']))
            exon_numbers.append(intern(attrs['exon_number']))

    df = pd.DataFrame(
        { 'seqname' : seqnames,
          'start'   : starts,
          'end'     : ends,
          'strand'  : strands,
          'gene_id' : gene_ids,
          'transcript_id' : transcript_ids,
          'exon_number' : exon_numbers })

    df = df.sort(['gene_id', 'transcript_id', 'exon_number'])

    for (gene_id, transcript_id), grp in df.groupby(['gene_id', 'transcript_id']):
        exons = []
        ref = None

        strand = grp.strand[grp.index[0]]

        for i in grp.index:
            ref = grp.seqname[i]
            exons.append(SeqFeature(
                id='{0}:{1}-{2}'.format(ref, grp.start[i] + 1, grp.end[i]),
                ref=ref,
                strand=strand,
                location=FeatureLocation(grp.start[i], grp.end[i])))
        yield SeqFeature(
            ref=ref,
            strand=strand,
            id=transcript_id,
            location=FeatureLocation(min(grp.start), max(grp.end)),
            sub_features=exons)
        

    # for line in fh:
    #     line = line.rstrip()
    #     (seqname, source, feature, start, end, score, strand, frame, attribute) = line.split('\t')

    #     seqname = intern(seqname)
    #     source  = intern(seqname)
    #     feature = intern(feature)
    #     start   = int(start) - 1
    #     end     = int(end)
    #     score = float(score)
    #     frame = int(frame)
    #     attributes = attribute.split(';\s*')
        
    #     data[seqname][gene_id][transcript_id][exon_number] = start, end
        
    #     if feature != 'exon':
    #         raise Exception("I only recognize files that have only exon features")

def parse_rum_index_genes(fh):
    for line in fh:

        (chr_, strand, start, end, num_exons, exon_starts, exon_ends, gene_name) = line.rstrip().split("\t")

        strand = 1 if strand == '+' else -1 if strand == '-' else None

        exon_starts = map(int, exon_starts.rstrip(",").split(","))
        exon_ends   = map(int, exon_ends.rstrip(",").split(","))

        exon_starts = [ x - 1 for x in exon_starts ]

        exon_locs = [FeatureLocation(*x, strand=strand) for x in zip(exon_starts, exon_ends)]

        exons = [ SeqFeature(id="{0}:{1}-{2}".format(chr_, loc.start, loc.end), 
                             ref=chr_, location=loc, strand=strand, type='exon')
                  for loc in exon_locs ]
        gene_loc = FeatureLocation(int(start) - 1, int(end), strand=strand)
        gene = SeqFeature(id=gene_name,
                          ref=chr_, location=gene_loc, strand=strand, type='gene', sub_features=exons)

        yield gene

def read_exons(fh):
    pat = re.compile('(.*):(\d+)-(\d+)')
    for line in fh:
        line = line.rstrip()
        m = pat.match(line)
        if m is None:
            raise Exception("Can't parse exon from '{line}'".format(line=line))
        (chr_, start, end) = m.groups()

        yield SeqFeature(ref=chr_, 
                         id=line,
                         location=FeatureLocation(int(start) -1, int(end)))

class TranscriptIndex(object):

    def __init__(self, transcripts):

        self.exon_to_transcripts = defaultdict(list)
        exons = []
        for t in transcripts:
            for e in t.sub_features:
                key = (e.ref, e.strand, e.location.start, e.location.end)
                self.exon_to_transcripts[key].append(t)
                exons.append(e)

        self.exon_index = ExonIndex(exons)

    def get_transcripts(self, ref, strand, start, end):
        
        seen = set()

        exons = self.exon_index.get_exons(ref, strand, start, end)
        for e in exons:
            key = (e.ref, e.strand, e.location.start, e.location.end)
            for t in self.exon_to_transcripts[key]:
                if t.id not in seen:
                    seen.add(t.id)
                    yield(t)

    def get_exons(self, ref, strand, start, end):
        return self.exon_index.get_exons(ref, strand, start, end)

class ExonIndex(object):

    def __init__(self, exons):

        # From the list of n exons, build a list of 2n 'events'. Each
        # exon gets n event representing its start, and on
        # representing its end.
        events = []
        Event = namedtuple('Event', ['etype', 'location', 'exon'])
        logging.debug("Converting exon locations into start and end events")
        for exon in exons:
            start = exon.location.start
            end   = exon.location.end

            events.append(Event('start', start, exon))
            events.append(Event('end',   end,   exon))

        # Now sort the events by location, and iterate over them. As we
        # go, keep a set of the exons that we're currently
        # overlapping. Build a list of sorted, disjoint spans, and
        # associate with each span the list of exons that are present in
        # all bases of that span.

        key_fn = lambda x: (x.exon.ref, x.exon.strand, x.location)

        logging.debug('Sorting events')
        events = sorted(events, key=key_fn)

        overlap = defaultdict(set)

        starts = defaultdict(list)
        keys   = defaultdict(list)

        old = set()

        logging.debug("Overlap list based on sorted events")
        for (ref, strand, location), values in groupby(events, key_fn):

            for event in values:
                start = event.exon.location.start
                end   = event.exon.location.end
                key = (event.exon.id, ref, strand, start, end)

                if event.etype == 'start':
                    overlap[ref].add(key)
                else:
                    try:
                        overlap[ref].remove(key)
                    except KeyError as e:
                        pass
            starts[(ref, strand)].append(int(location))
            keys[(ref, strand)].append(set(overlap[ref]))

        for k in starts:
            starts[k] = np.array(starts[k], int)

        self.start = starts
        self.keys  = keys
        
    def get_exons(self, ref, strand, start, end):

        key = (ref, strand)
        starts = self.start[(ref, strand)]
        exon_keys = self.keys[(ref, strand)]

        n = len(starts)
        
        rs = np.searchsorted(starts, [start], side='right')
        r = max(rs[0] - 1, 0)

        seen = set()

        while r < n and starts[r] < end:
            for key in exon_keys[r]:
                if key not in seen:
                    (exon_id, exon_ref, exon_strand, exon_start, exon_end) = key
                    seen.add(key)
                    yield(SeqFeature(id=exon_id, ref=exon_ref, strand=exon_strand, location=FeatureLocation(exon_start, exon_end)))
            r += 1

def genes_to_exons(genes):
    """Given an iterator over genes, returns an iterator over exons.

    """
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)
