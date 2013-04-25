from Bio.SeqFeature import FeatureLocation, SeqFeature
import re

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
