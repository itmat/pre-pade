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
}


int main (int argc, char **argv) {
  test_create_index();

  return check_results();
  
}


/*
	protein_coding	exon	28692193	28692362
	.	-	.	 gene_id "AT1G76480"; transcript_id "AT1G76480.1"; exon_number "6"; gene_name "F14G6.8"; transcript_name "F14G6.8-202"; seqedit "false";

*/
