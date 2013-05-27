#include <stdio.h>
#include "geneindex.h"

struct Assertion {
  char *name;
  int passed;
};

struct {
  struct Assertion assertions[100];
  int assertions_len;
} RESULTS;

void add_assertion(int passed, char *name) {
  struct Assertion *assertion = RESULTS.assertions + RESULTS.assertions_len++;
  assertion->name = name;
  assertion->passed = passed;
}

void assert_equals(int a, int b, char *name) {
  add_assertion(a == b, name);
}

void assert_str_equals(char *a, char *b, char *name) {
  add_assertion(!strcmp(a, b), name);
}

void assert_exon_ptr_equals(struct Exon *a, struct Exon *b, char *name) {
  add_assertion(a == b, name);
}

int check_results() {
  int failures = 0;

  int i;
  struct Assertion *a = RESULTS.assertions;
  for (i = 0; i < RESULTS.assertions_len; i++) {
    
    if (!a->passed) {
      printf("Failed: %s\n", a->name);
      failures++;
    }
    a++;
  }

  printf("Ran %d tests\n", RESULTS.assertions_len);
  if (failures) {
    printf("%d failures!\n", failures);
    return 1;
  }
  return 0;
}

void test_create_index() {
  struct ExonDB exondb;
  parse_gtf_file(&exondb, "testdata/arabidopsis.gtf");
  index_exons(&exondb);  

  assert_equals(13214, exondb.exons_len, "Number of exons loaded");

  struct Exon *exon = exondb.exons;
  assert_str_equals("1", exon->chrom, "First exon chromosome");
  assert_equals(28692193, exon->start, "First exon start");
  assert_equals(28692362, exon->end, "First exon end");
  assert_str_equals("protein_coding", exon->source, "First exon source");
  assert_str_equals("exon", exon->feature, "First exon feature");

  exon = search_exons(&exondb, "1", 0, 0);
  assert_exon_ptr_equals(exon, exondb.exons, "Search for first exon");
}


int main (int argc, char **argv) {
  test_create_index();
  return check_results();
}

const int CROSS_EXON_START =  1;
const int CROSS_EXON_END   =  2;
const int START_IN_EXON    =  4;
const int END_IN_EXON      =  8;
const int WRONG_CHROMOSOME = 16;
const int START_AFTER_EXON = 32;
const int END_BEFORE_EXON  = 64;


/* Compare the given exon to the specified range and return flags
 * indicating if and how it overlaps. 0 means it matches exactly.
 */
int cmp_exon(struct Exon *e, char *chrom, int start, int end) {

  int result = 0;

  // ...eeee
  //         rrrr...
  if ( e->end < start )
    result |= START_AFTER_EXON;

  //         eeee...
  // ...rrrr
  if ( end < e->start )
    result |= END_BEFORE_EXON;

  // ...eeee
  //   ...rrrr...
  if ( start < e->end && e->end < end )
    result |= CROSS_EXON_END;

  //      eeee...
  // ...rrrr...
  if ( start < e->start && e->start < end ) 
    result |= CROSS_EXON_START;

  //  ...eeee...
  //       rrrr...
  if ( e->start < start && start < e->end ) 
    result |= START_IN_EXON;

  //   ...eeee...
  // ...rrrr
  if ( e->start < end && end < e->end ) 
    result |= END_IN_EXON;

  return result;
}

struct Exon *finish_cursor(struct ExonCursor *cursor) {
  cursor->next = NULL;
  cursor->allow = 0;
  return NULL;
}

struct Exon *next_exon(struct ExonCursor *cursor, int *flags) {

  struct Exon   *e;
  struct ExonDB *db = cursor->exondb;
  struct Exon *last_exon = db->exons + db->exons_len;

  int allow = cursor->allow;
  int disallow = ~allow;

  while (cursor->next <= last_exon) {

    struct Exon *exon = cursor->next++;
    int cmp = cmp_exon(exon, cursor->chrom, cursor->start, cursor->end);

    if ( cmp & WRONG_CHROMOSOME ) {
      // If it's on the wrong chromosome, we're definitely done.
      return finish_cursor(cursor);
    }
   
    else if (e->min_start > cursor->end) {
      // At this point all subsequent exons will start after me, so
      // we're done.
      return finish_cursor(cursor);
    }

    else if (cmp & disallow) {
      // If this match would be disallowed according to the flags the
      // user passed in, just skip it.
      continue;
    }

    else {
      // Otherwise it's a match
      *flags = cmp;
      return exon;
    }
  }

  // If we get to this point, we've passed the last exon in the
  // database, so we're done.
  return finish_cursor(cursor);
}


/*	
	gene_id "AT1G76480"; 
        transcript_id "AT1G76480.1"; 
        exon_number "6"; 
        gene_name "F14G6.8"; 
        transcript_name "F14G6.8-202"; 
        seqedit "false";

*/
