#include <string.h>
#include "testutils.h"
#include "quant.h"

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
  add_assertion(!strcmp(a, b), msg);
}

void assert_not_null(void *x, char *name) {
  char *msg;
  asprintf(&msg, "%s: expected non null, was null", name);
  add_assertion(x != NULL, msg);
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

