#ifndef GENEINDEX_H
#define GENEINDEX_H

#define CROSS_EXON_START  1
#define CROSS_EXON_END    2
#define START_IN_EXON     4
#define END_IN_EXON       8
#define WRONG_CHROMOSOME 16
#define START_AFTER_EXON 32
#define END_BEFORE_EXON  64

enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

struct Exon {
  char *gene_id;
  char *transcript_id;
  int   exon_number;
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

  struct ExonIndexEntry *index;
  int index_len;
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

struct ExonIndexEntry {
  char *chrom;
  int start;
  int end;
  struct Exon *exon;
};


int cmp_index_entry(struct ExonIndexEntry *key,
                    struct ExonIndexEntry *entry);

struct Exon * search_exons(struct ExonDB *exondb, char *chrom, int start, int end);
struct Exon *next_exon(struct ExonCursor *cursor, int *flags);

#endif
