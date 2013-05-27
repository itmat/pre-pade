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
  char *msg;
  asprintf(&msg, "%s: expected %d, got %d", name, a, b);
  add_assertion(a == b, msg);
}

void assert_str_equals(char *a, char *b, char *name) {
  char *msg;
  asprintf(&msg, "%s: expected %s, got %s", name, a, b);
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

  exon = search_exons(&exondb, "foo", 0, 0);
  assert_exon_ptr_equals(NULL, exon, "Search for first exon");

  exon = search_exons(&exondb, "0", 0, 0);
  assert_exon_ptr_equals(exondb.exons, exon, "Search for first exon");

  
}

int test_compare_index_entry() {
  struct ExonIndexEntry entry = { "chr2", 100, 200 };

  struct ExonIndexEntry left_chrom  = { "chr1", 100 };
  struct ExonIndexEntry left        = { "chr2", 50 };
  struct ExonIndexEntry inside      = { "chr2", 150 };
  struct ExonIndexEntry right       = { "chr2", 250 };
  struct ExonIndexEntry right_chrom = { "chr3", 100 };
  
  assert_equals(-1, cmp_index_entry(&left_chrom, &entry), "Left chrom");
  assert_equals(-1, cmp_index_entry(&left, &entry), "Left");
  assert_equals(0, cmp_index_entry(&inside, &entry), "Inside");
  assert_equals(1, cmp_index_entry(&right, &entry), "Right");
  assert_equals(1, cmp_index_entry(&right_chrom, &entry), "Right chrom");
}

int test_compare_exon() {
  struct Exon e;
  e.chrom = "chr1";
  e.start = 100;
  e.end   = 200;

  struct ExonCompTestCase {
    int expected;
    char *chrom;
    int start;
    int end;
    char *name;
  } cases[] = {
    { 0,                "chr1", 100, 200, "Exact match" },
    { WRONG_CHROMOSOME, "chr2", 100, 200, "Wrong chromosome" },
    { START_AFTER_EXON, "chr1", 250, 300, "Start after exon" },
    { END_BEFORE_EXON,  "chr1",   0,  75, "End before exon" },
    { CROSS_EXON_END | 
      START_IN_EXON,    "chr1", 150, 250, "Start in exon, cross exon end" },
    { CROSS_EXON_START | 
      END_IN_EXON,      "chr1",  50, 150, "Cross exon start, end in exon" },
    { START_IN_EXON |
      END_IN_EXON,      "chr1", 125, 175, "Start and end in exon" },
    { CROSS_EXON_START |
      CROSS_EXON_END,   "chr1", 50, 250, "Cross exon start and end" }
  };

  int n = sizeof(cases) / sizeof(struct ExonCompTestCase);
  int i;

  for (i = 0; i < n; i++) {
    struct ExonCompTestCase *tc = cases + i;
    int got = cmp_exon(&e, tc->chrom, tc->start, tc->end);
    assert_equals(tc->expected, got, cases[i].name);
  }
}

int main (int argc, char **argv) {
  test_create_index();
  test_compare_exon();
  test_compare_index_entry();
  return check_results();
}

/*	
	gene_id "AT1G76480"; 
        transcript_id "AT1G76480.1"; 
        exon_number "6"; 
        gene_name "F14G6.8"; 
        transcript_name "F14G6.8-202"; 
        seqedit "false";

*/
