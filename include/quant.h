#ifndef GENEINDEX_H
#define GENEINDEX_H

#include <stdio.h>
#include "sam.h"

/* These flags are used to characterize the overlap of a read span and
   a span in the model.

   CROSS_EXON_START and CROSS_EXON_END indicate that the read span
   crosses a splice point. This is generally not allowed, although we
   may want to consider allowing it for the first and last exon in the
   model, since apparently the start of the first exon and end of the
   last exon are hard to pin down.

   START_IN_EXON and END_IN_EXON indicate that the read span starts or
   ends inside the exon. Generally the first span is allowed to start
   in an exon and the last is allowed to end in an exon. 

   WRONG_CHROMOSOME means it's on the wrong chromosome, and
   START_AFTER_EXON and END_BEFORE_EXON indicate that the read span is
   on the right chromosome but doesn't overlap at all.
*/
#define CROSS_EXON_START  1
#define CROSS_EXON_END    2
#define START_IN_EXON     4
#define END_IN_EXON       8
#define WRONG_CHROMOSOME 16
#define START_AFTER_EXON 32
#define END_BEFORE_EXON  64

// Allow any overlap at all, even if it might contradict the junction
// model.
#define ALLOW_ALL START_IN_EXON | END_IN_EXON | CROSS_EXON_START | CROSS_EXON_END


/* These are some logging utilities. Compile with -DLOG_LEVEL_TRACE or
   -DLOG_LEVEL_DEBUG to turn on. Warning: TRACE may result in a
   tremendous amount of output, so you probably only want to use it
   when you're debugging a specific case. I would recommend tailoring
   your input files to include only the cases you're interested in if
   you're doing trace-level debugging.*/
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

/* Just so we don't have to type "struct" all over the place. I don't
   do any typedefs to pointers. */
typedef struct Span Span;
typedef struct Region Region;
typedef struct RegionList RegionList;
typedef struct GeneModel GeneModel;
typedef struct Transcript Transcript;
typedef struct RegionCursor RegionCursor;
typedef struct Quant Quant;
typedef struct IndexEntry IndexEntry;
typedef struct RegionMatch RegionMatch;
typedef struct RegionMatches RegionMatches;
typedef struct CigarCursor CigarCursor;

enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

// A read span
struct Span {
  int start, end;
};

// Stores quantification counts. At each level (e.g. exon, junction),
// we store the minimum and maximum possible counts. Minimum includes
// only unique mappers, and maximum includes unique and non-unique.
struct Quant {
  int min, max;
};

// Represents an exon, read from the GTF file, with some additional
// information that allows us to index it for searching, and some
// fields for storing quantification counts.
struct Region {

  char *gene_id;
  char *transcript_id;
  int   exon_number;
  char *source;
  char *feature;
  char *chrom;
  enum Strand strand;
  int start;
  int end;

  // Pointer to the transcript that this exon is a part of. Each exon
  // belongs to exactly one transcript. We do not attempt to
  // consolidate exons with identical coordinates together.
  Transcript *transcript;

  // Extra, not part of GTF file, used for indexing. The exons are
  // sorted by end coordinate, and min_start indicates the lowest
  // starting coordinate for any exon that ends at or after the end
  // point of this exon.
  int min_start;

  Quant exon_quant;
  Quant junction_quant;
};


// List of exons. 
struct RegionList {
  struct Region *items;
  int len; // Current length of items
  int cap; // Current capacity (number of allocated items). len <= cap.

  struct IndexEntry *index;
  int index_len;


};

/* This is the main data structure used to store our transcript
   model. It has a list of exons, sorted by end coordinate, which
   allows searching based on coordinate.  */
struct GeneModel {
  RegionList exons;

  Transcript *transcripts;
  int num_transcripts;
};


/* A transcript is just a list of pointers to the exons that it is
   made of, along with a transcript id. */
struct Transcript {
  Region **exons;
  int exons_len;
  int exons_cap;
  char *id;
};


/* Used to store the state of a query in an GeneModel. */
struct RegionCursor {

  // The database we searched in
  struct RegionList *list;
  
  // The search coordinates
  char *chrom;
  int start, end;

 // Flags indicating what kinds of overlap are allowed
  int allow;

 // Next match to be returned
  struct Region *next;
};

/* The index is made up of a list of IndexEntry structs, which
   divide the genome into non-overlapping segments. This allows us to
   do a binary search on the index given some coordinates.  */
struct IndexEntry {

  // The coordinates of this range
  char *chrom;
  int start, end;

  // The first exon that overlaps the range specified above. To find
  // all matching exons given a span, we search for the IndexEntry
  // that contains the start coordinate. That entry points to the
  // first matching exon. Then we advance through the exon list (just
  // doing exon++ since they're stored in sorted order) until we get
  // to the first one where exon->min_start is greater than the end
  // point of the query.
  struct Region *exon;
};

// Represents a match (or mismatch) of an exon for a particular read span. 
struct RegionMatch {

  // The exon we found
  Region *exon;

  // Does the exon overlap the read span?
  int overlap;

  // Does the exon conflict with the read span?
  int conflict;
};

// A list of RegionMatch structs.
struct RegionMatches {
  RegionMatch *items;
  int len;
  int cap;
};

// Cursor used for converting a CIGAR string into a list of Span structs.
struct CigarCursor {
  bam1_t *read;
  int i;
  int start;
  int end;
  int pos;
};


/* These functions should all be documented where they're defined. */


int cmp_index_entry(struct IndexEntry *key,
                    struct IndexEntry *entry);
int search_exons(struct RegionCursor *cursor,
                 struct RegionList *list, char *chrom, int start, int end, int allow);
struct Region *next_exon(struct RegionCursor *cursor, int *flags);
int parse_gtf_attr_str(char *str, char *name, char **dest);
int parse_gtf_attr_int(char *str, char *name, int *value);
void consolidate_exon_matches(RegionMatches *matches);
int cmp_match_by_exon(RegionMatch *a, RegionMatch *b);
void parse_gtf_file(GeneModel *gm, char *filename);
void index_exons(RegionList *exons);
void init_exon_matches(RegionMatches *matches);
void find_candidates(RegionMatches *matches, GeneModel *gm, char *ref,
                     Span *spans, int num_fwd_spans, int num_rev_spans);
int next_fragment(bam1_t **reads, samfile_t *samfile, int n);
void init_cigar_cursor(struct CigarCursor *c, bam1_t *read);
int next_span(struct CigarCursor *c);
int extract_spans(Span *spans, bam1_t *read, int n);
int cmp_exon(Region *e, char *chrom, int start, int end);
void add_match(RegionMatches *matches, Region *exon, int overlap, int conflict);
Region *next_exon_in_transcript(Region *e);
void add_transcripts(GeneModel *gm);
void incr_quant(Quant *q, int unique);
int matches_junction(Region *left, Span *spans, int num_fwd_spans, int num_rev_spans, int min_overlap);
int load_model(GeneModel *gm, char *filename);

#endif

