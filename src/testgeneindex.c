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

void test_create_index() {
  struct ExonDB exondb;
  parse_gtf_file(&exondb, "testdata/arabidopsis.gtf");
  index_exons(&exondb);  
  assert_equals(30000, exondb.exons_len, "Number of exons loaded");
  
}


int main (int argc, char **argv) {
  test_create_index();

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
