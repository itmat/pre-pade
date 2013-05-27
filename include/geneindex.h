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

struct ExonCursor {

  // The database we searched in
  struct ExonDB *exondb;
  
  // Search criteria
  char *chrom;
  int start;
  int end;
  int allow;

  // Next match to be returned
  struct Exon *next;
};


struct Exon * search_exons(struct ExonDB *exondb, char *chrom, int start, int end);
struct Exon *next_exon(struct ExonCursor *cursor, int *flags);
