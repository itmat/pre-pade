#ifndef GENEINDEX_H
#define GENEINDEX_H

#include <stdio.h>

#define CROSS_EXON_START  1
#define CROSS_EXON_END    2
#define START_IN_EXON     4
#define END_IN_EXON       8
#define WRONG_CHROMOSOME 16
#define START_AFTER_EXON 32
#define END_BEFORE_EXON  64

#define ALLOW_ALL START_IN_EXON | END_IN_EXON | CROSS_EXON_START | CROSS_EXON_END

#define LOG_MESSAGE(fmt, level, args...) fprintf(stderr, "%s %s (%s:%d): "fmt, level, __FUNCTION__, __FILE__, __LINE__, args)

#if defined( LOG_LEVEL_TRACE )
#  define LOG_LEVEL_DEBUG
#  define LOG_TRACE(fmt, args...) LOG_MESSAGE(fmt, "TRACE", args)
#else
#  define LOG_TRACE(fmt, args...)
#endif

#if defined(LOG_LEVEL_DEBUG)
#  define LOG_LEVEL_INFO
#  define LOG_DEBUG(fmt, args...) LOG_MESSAGE(fmt, "DEBUG", args)
#else
#  define LOG_DEBUG(fmt, args...)
#endif

#if defined( LOG_LEVEL_INFO )
#  define LOG_INFO(fmt, args...) LOG_MESSAGE(fmt, "INFO", args)
#else
#  define LOG_INFO(fmt, args...)
#endif



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


struct QuantCandidate {
  
  char *quant_type;
  
  // The feature
  char *gene_id;
  char *transcript_id;
  int exon_number;
  char *chrom;
  int start;
  int end;
  
  char *read_id;
  int alignment_number;
  int num_alignments;
  int consistent;
  
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
