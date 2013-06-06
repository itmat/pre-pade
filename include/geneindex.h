#ifndef GENEINDEX_H
#define GENEINDEX_H

#define CROSS_EXON_START  1
#define CROSS_EXON_END    2
#define START_IN_EXON     4
#define END_IN_EXON       8
#define WRONG_CHROMOSOME 16
#define START_AFTER_EXON 32
#define END_BEFORE_EXON  64

#define ALLOW_ALL START_IN_EXON | END_IN_EXON | CROSS_EXON_START | CROSS_EXON_END

enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

typedef struct Exon Exon;
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

  int min_count;
  int max_count;
};

typedef struct ExonList ExonList;
struct ExonList {
  struct Exon *items;
  int len;
  int cap;
};

typedef struct ExonDB ExonDB;
struct ExonDB {
  ExonList exons;

  struct ExonIndexEntry *index;
  int index_len;
};

typedef struct ExonCursor ExonCursor;

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

typedef struct ExonIndexEntry ExonIndexEntry;
struct ExonIndexEntry {
  char *chrom;
  int start;
  int end;
  struct Exon *exon;
};


typedef struct ExonMatch ExonMatch;

struct ExonMatch {
  Exon *exon;
  int overlap;
  int conflict;
};

typedef struct ExonMatches ExonMatches;

struct ExonMatches {
  ExonMatch *items;
  int len;
  int cap;
};



int cmp_index_entry(struct ExonIndexEntry *key,
                    struct ExonIndexEntry *entry);

int search_exons(struct ExonCursor *cursor,
                 struct ExonDB *exondb, char *chrom, int start, int end, int allow);

struct Exon *next_exon(struct ExonCursor *cursor, int *flags);

int parse_gtf_attr_str(char *str, char *name, char **dest);
int parse_gtf_attr_int(char *str, char *name, int *value);

void consolidate_exon_matches(ExonMatches *matches);
int cmp_match_by_exon(ExonMatch *a, ExonMatch *b);

#endif
