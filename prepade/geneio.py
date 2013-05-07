from Bio.SeqFeature import FeatureLocation, SeqFeature
from collections import namedtuple, defaultdict
from itertools import groupby

import logging
import re
import numpy as np

def parse_rum_index_genes(fh):
    for line in fh:
        (chr_, strand, start, end, num_exons, exon_starts, exon_ends, gene_name) = line.split("\t")

        strand = 1 if strand == '+' else -1 if strand == '-' else None

        exon_starts = map(int, exon_starts.rstrip(",").split(","))
        exon_ends   = map(int, exon_ends.rstrip(",").split(","))

        exon_starts = [ x - 1 for x in exon_starts ]

        exon_locs = [FeatureLocation(*x, strand=strand) for x in zip(exon_starts, exon_ends)]

        exons = [ SeqFeature(id="{0}:{1}-{2}".format(chr_, loc.start, loc.end), 
                             ref=chr_, location=loc, strand=strand, type='exon')
                  for loc in exon_locs ]
        gene_loc = FeatureLocation(int(start), int(end), strand=strand)
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
                key = (e.ref, e.location.start, e.location.end)
                self.exon_to_transcripts[key].append(t)
                exons.append(e)

        self.exon_index = ExonIndex(exons)

    def get_transcripts(self, ref, start, end):
        
        seen = set()

        exons = self.exon_index.get_exons(ref, start, end)
        for e in exons:
            key = (e.ref, e.location.start, e.location.end)
            for t in self.exon_to_transcripts[key]:
                if t.id not in seen:
                    seen.add(t.id)
                    yield(t)

    def get_exons(self, ref, start, end):
        return self.exon_index.get_exons(ref, start, end)

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

        key_fn = lambda x: (x.exon.ref, x.location)

        logging.debug('Sorting events')
        events = sorted(events, key=key_fn)

        overlap = defaultdict(set)

        starts = defaultdict(list)
        keys   = defaultdict(list)


        old = set()

        logging.debug("Overlap list based on sorted events")
        for (ref, location), values in groupby(events, key_fn):

            for event in values:
                start = event.exon.location.start
                end   = event.exon.location.end
                key = (event.exon.id, ref, start, end)

                if event.etype == 'start':
                    overlap[ref].add(key)
                else:
                    try:
                        overlap[ref].remove(key)
                    except KeyError as e:
                        pass

            starts[ref].append(location)
            keys[ref].append(set(overlap[ref]))

        self.start = starts
        self.keys  = keys
        
    def get_exons(self, ref, start, end):

        n = len(self.start[ref])
        p = 0
        q = len(self.start[ref])
        
        while (p < q):
            r = p + (q - p) // 2

            if r + 1 < n and self.start[ref][r + 1] <= start:
                p = r

            elif self.start[ref][r] > start:
                q = r

            else:
                break

        seen = set()

        while r < n and self.start[ref][r] < end:
            for key in self. keys[ref][r]:
                if key not in seen:
                    (exon_id, exon_ref, exon_start, exon_end) = key
                    seen.add(key)
                    yield(SeqFeature(id=exon_id, ref=exon_ref, location=FeatureLocation(exon_start, exon_end)))
            r += 1

def genes_to_exons(genes):
    """Given an iterator over genes, returns an iterator over exons.

    """
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)
