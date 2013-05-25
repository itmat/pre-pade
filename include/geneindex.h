enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

struct Exon {
  char *gene_id;
  char *transcript_id;
  int exon_number;
  char *source;
  char *feature;
  char *chrom;
  enum Strand strand;
  int start;
  int end;

  // Extra, not part of GTF file, used for indexing
  int min_start;
};

struct ExonDB {
  struct Exon *exons;
  int exons_len;
  int exons_cap;
};
