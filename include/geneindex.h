#ifndef GENEINDEX_H
#define GENEINDEX_H

#include <stdio.h>
#include "sam.h"

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


#define STRAND_PLUS '+'
#define STRAND_MINUS '-'


enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

typedef struct Span Span;
struct Span {
  int start;
  int end;
};



typedef struct Exon Exon;
typedef struct ExonList ExonList;
typedef struct ExonDB ExonDB;
typedef struct Transcript Transcript;
typedef struct ExonCursor ExonCursor;
typedef struct Quant Quant;

struct Quant {
  int min, max;
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

  Transcript *transcript;

  // Extra, not part of GTF file, used for indexing
  int min_start;

  Quant exon_quant;
  Quant junction_quant;
};



struct ExonList {
  struct Exon *items;
  int len;
  int cap;
};


struct ExonDB {
  ExonList exons;

  struct ExonIndexEntry *index;
  int index_len;

  Transcript *transcripts;
  int num_transcripts;
};



struct Transcript {
  Exon **exons;
  int exons_len;
  int exons_cap;
  char *id;
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
void parse_gtf_file(ExonDB *exondb, char *filename);
void index_exons(ExonDB *exondb);
void init_exon_matches(ExonMatches *matches);
void find_candidates(ExonMatches *matches, ExonDB *db, char *ref,
                     Span *spans, int num_fwd_spans, int num_rev_spans);

#define MAX_SPANS_PER_READ 1000

typedef struct CigarCursor CigarCursor;
struct CigarCursor {
  bam1_t *read;
  int i;
  int start;
  int end;
  int pos;
};

int next_fragment(bam1_t **reads, samfile_t *samfile, int n);
void init_cigar_cursor(struct CigarCursor *c, bam1_t *read);
int next_span(struct CigarCursor *c);
int extract_spans(Span *spans, bam1_t *read, int n);

int cmp_exon(Exon *e, char *chrom, int start, int end);
void add_match(ExonMatches *matches, Exon *exon, int overlap, int conflict);
Exon *next_exon_in_transcript(Exon *e);
void add_transcripts(ExonDB *db);

void incr_quant(Quant *q, int unique);

#endif
