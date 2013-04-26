from Bio.SeqFeature import FeatureLocation, SeqFeature
from collections import namedtuple
from itertools import groupby
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

        exons = [ SeqFeature(ref=chr_, location=loc, strand=strand, type='exon')
                  for loc in exon_locs ]
        gene_loc = FeatureLocation(int(start), int(end), strand=strand)
        gene = SeqFeature(ref=chr_, location=gene_loc, strand=strand, type='gene', sub_features=exons)

        yield gene

def read_exons(fh):
    pat = re.compile('(.*):(\d+)-(\d+)')
    for line in fh:
        m = pat.match(line)
        (chr_, start, end) = m.groups()

        yield SeqFeature(ref=chr_, 
                         location=FeatureLocation(int(start) -1, int(end)))

class ExonIndex(object):

    def __init__(self, exons):

        # From the list of n exons, build a list of 2n 'events'. Each
        # exon gets n event representing its start, and on
        # representing its end.
        events = []
        Event = namedtuple('Event', ['etype', 'location', 'key'])
        for exon in exons:
            ref   = exon.ref
            start = exon.location.start
            end   = exon.location.end
            key = '%s:%d-%d' % (ref, start + 1, end)
            events.append(Event('start', start, key))
            events.append(Event('end',   end,   key))

        # Now sort the events by location, and iterate over them. As we
        # go, keep a set of the exons that we're currently
        # overlapping. Build a list of sorted, disjoint spans, and
        # associate with each span the list of exons that are present in
        # all bases of that span.
        events = sorted(events, key=lambda x: x.location)

        overlap = set()

        starts = []
        keys   = []

        for location, values in groupby(events, lambda x: x.location):
        
            for event in values:
                if event.etype == 'start':
                    overlap.add(event.key)
                else:
                    overlap.remove(event.key)

            starts.append(location)
            keys.append(set(overlap))



        self.start = np.array(starts, int)

        self.keys  = keys
        
    def get_exons(self, start, end):
        n = len(self.start)
        p = 0
        q = len(self.start)
        
        while (p < q):
            r = p + (q - p) // 2

            if r + 1 < n and self.start[r + 1] <= start:
                p = r

            elif self.start[r] > start:
                q = r

            else:
                break

        keys = set()

        while r < n and self.start[r] < end:
            keys.update(self.keys[r])
            r += 1

        return read_exons(keys)
