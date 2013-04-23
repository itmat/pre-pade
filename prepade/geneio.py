from Bio.SeqFeature import FeatureLocation, SeqFeature

def parse_rum_index_genes(fh):
    for line in fh:
        (_chr, strand, start, end, num_exons, exon_starts, exon_ends, gene_name) = line.split("\t")

        strand = 1 if strand == '+' else -1 if strand == '-' else None

        exon_starts = map(int, exon_starts.rstrip(",").split(","))
        exon_ends   = map(int, exon_ends.rstrip(",").split(","))

        exon_locs = [FeatureLocation(*x, strand=strand) for x in zip(exon_starts, exon_ends)]

        exons = [ SeqFeature(location=loc, strand=strand, type='exon')
                  for loc in exon_locs ]
        gene_loc = FeatureLocation(int(start), int(end), strand=strand)
        gene = SeqFeature(location=gene_loc, strand=strand, type='gene', sub_features=exons)

        yield (_chr, gene)
