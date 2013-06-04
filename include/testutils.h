#include "geneindex.h"

struct Assertion {
  char *name;
  int passed;
};

struct {
  struct Assertion assertions[100];
  int assertions_len;
} RESULTS;


void add_assertion(int passed, char *name);
void assert_equals(int a, int b, char *name);
void assert_str_equals(char *a, char *b, char *name);
void assert_exon_ptr_equals(struct Exon *a, struct Exon *b, char *name);
int check_results();
