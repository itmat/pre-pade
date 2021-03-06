from Bio.SeqFeature import FeatureLocation, SeqFeature
from collections import namedtuple, defaultdict
from itertools import groupby, chain

import logging
import numpy as np
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

    data = {}

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
            gene_id       = attrs['gene_id']
            transcript_id = attrs['transcript_id']
            exon_number   = int(attrs['exon_number'])

            key = (gene_id, transcript_id, seqname, strand)
            if key not in data:
                data[key] = {}
            data[key][int(exon_number)] = [int(start) - 1, int(end)]

    for key in sorted(data):
        gene_id, transcript_id, chrom, strand = key
        exons = []
        ref = None

        gene_start = np.inf
        gene_end   = - np.inf

        for exon_number in sorted(data[key]):
            start, end = data[key][exon_number]
            gene_start = min(gene_start, start)
            gene_end   = max(gene_end, end)
            exons.append(SeqFeature(
                id='{0}:{1}-{2}'.format(chrom, start + 1, end),
                ref=chrom,
                strand=strand,
                location=FeatureLocation(start, end)))

        yield SeqFeature(
            ref=chrom,
            strand=strand,
            id=transcript_id,
            location=FeatureLocation(gene_start, gene_end),
            sub_features=exons)
        

def parse_rum_index_genes(fh):
    for line in fh:

        (chr_, strand, start, end, num_exons, exon_starts, exon_ends, gene_name) = line.rstrip().split("\t")

        strand = 1 if strand == '+' else -1 if strand == '-' else None

        exon_starts = map(int, exon_starts.rstrip(",").split(","))
        exon_ends   = map(int, exon_ends.rstrip(",").split(","))

        # RUM index represents exons with 1-based indexing, but we use
        # 0-based internally.
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
    """Parse an exon-only input file, returning the exons as an iterator
    of SeqFeature objects.

    """
    pat = re.compile('(.*):(\d+)-(\d+)')
    for line in fh:
        line = line.rstrip()
        m = pat.match(line)
        if m is None:
            raise Exception("Can't parse exon from '{line}'".format(line=line))
        (chr_, start, end) = m.groups()

        yield SeqFeature(ref=chr_, 
                         id=line,
                         strand=0,
                         location=FeatureLocation(int(start) -1, int(end)))

class TranscriptIndex(object):

    def __init__(self, transcripts):

        self.exon_to_transcripts = defaultdict(list)
        exons = []
        for t in transcripts:
            for e in t.sub_features:
                
                key = Exon(id=e.id, ref=e.ref, strand=e.strand, start=e.location.start, end=e.location.end)
                self.exon_to_transcripts[key].append(t)
                exons.append(e)

        self.exon_index = ExonIndex(exons)

    def get_transcripts(self, ref, strand, start, end):
        
        seen = set()

        for exon in self.exon_index.get_exons(ref, strand, start, end):
            
            for t in self.exon_to_transcripts[exon]:
                if t.id not in seen:
                    seen.add(t.id)
                    yield(t)

    def get_exons(self, ref, strand, start, end):
        """Get a list of exons on the given strand if the given chromosome
        that overlap the range specified by start and end.

        """
        return self.exon_index.get_exons(ref, strand, start, end)


Exon = namedtuple('Exon', ['id', 'ref', 'strand', 'start', 'end'])

class ExonIndex(object):

    valid_strands = set([-1, 0, 1, None])

    def check_strand(self, strand):
        if strand not in self.valid_strands:
            raise Exception(
                ("Unknown strand {strand}, should be one of "+
                 "{valid_strands}").format(
                     strand=strand,
                     valid_strands=", ".join(map(str, self.valid_strands))))
            
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

        # Map from (ref, strand) to running set of exons that overlap
        # our current position.
        overlap = defaultdict(set)

        # Map from (ref, strand) to list of span start positions.
        starts = defaultdict(list)

        # Map from (ref, strand) to list of exon keys that overlap the
        # corresponding position
        keys   = defaultdict(list)

        # Map from (ref strand) to list of exon keys that overlap the
        # corresponding position and don't overlap the previous
        # position
        new_keys = defaultdict(list)

        logging.debug("Overlap list based on sorted events")
        for (ref, strand, location), values in groupby(events, key_fn):
            self.check_strand(strand)
            for event in values:
                start = event.exon.location.start
                end   = event.exon.location.end
                key = Exon(id=event.exon.id, start=int(start), end=int(end), strand=strand, ref=ref)

                if event.etype == 'start':
                    overlap[ref].add(key)
                else:
                    try:
                        overlap[ref].remove(key)
                    except KeyError as e:
                        pass

            starts[(ref, strand)].append(int(location))

            these_keys = set(overlap[ref])
            last_keys = new_keys[(ref, strand)][-1] if len(new_keys[(ref, strand)]) > 0 else set()
            keys[(ref, strand)].append(these_keys)

            new_keys[(ref, strand)].append(these_keys.difference(last_keys))

        for k in starts:
            starts[k] = np.array(starts[k], int)

        self.start = starts
        self.keys  = keys
        self.new_keys = new_keys
        
    def get_exons(self, ref, strand, start, end):
        self.check_strand(strand)

        if strand is None or strand == 0:
            for exon in chain(
                self.get_exons(ref,  1, start, end),
                self.get_exons(ref, -1, start, end)):
                yield exon
                
        key       = (ref, strand)
        starts    = self.start[key]
        exon_keys = self.keys[key]
        new_keys  = self.new_keys[key]

        n = len(starts)
        
        rs = np.searchsorted(starts, [start], side='right')
        r = max(rs[0] - 1, 0)

        first = True
        while r < n and starts[r] < end:
            if first:
                keys = exon_keys[r]
                first = False
            else:
                keys = new_keys[r]
            r += 1

            for key in keys:
                yield key


def genes_to_exons(genes):
    """Given an iterator over genes, returns an iterator over exons.

    """
    for gene in genes:
        for exon in gene.sub_features:
            yield(exon)
