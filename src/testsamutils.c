#include <stdio.h>
#include "samutils.h"
#include "geneindex.h"
#include "sam.h"
#include "testutils.h"

void test_next_fragment_paired() {
  char *sam_filename = "testdata/test_next_fragment_paired.sam";
  samfile_t *samfile = samopen(sam_filename, "r", NULL);  
  bam1_t *reads[] = { bam_init1(), bam_init1() };

  int len;

  len = next_fragment(reads, samfile, 2);
  assert_equals(2, len, "Num reads");
  assert_str_equals("seq.1", bam1_qname(reads[0]), "qname");
  assert_str_equals("seq.1", bam1_qname(reads[1]), "qname");

  len = next_fragment(reads, samfile, 2);
  assert_equals(2, len, "Num reads");
  assert_str_equals("seq.2", bam1_qname(reads[0]), "qname");
  assert_str_equals("seq.2", bam1_qname(reads[1]), "qname");

  len = next_fragment(reads, samfile, 2);
  assert_equals(0, len, "Num reads");

  samclose(samfile);
}


void test_next_fragment_single() {
  char *sam_filename = "testdata/test_next_fragment_single.sam";
  samfile_t *samfile = samopen(sam_filename, "r", NULL);  
  bam1_t *reads[] = { bam_init1(), bam_init1() };

  int len;

  len = next_fragment(reads, samfile, 2);
  assert_equals(1, len, "Num reads");
  assert_str_equals("seq.1", bam1_qname(reads[0]), "qname");

  len = next_fragment(reads, samfile, 2);
  assert_equals(1, len, "Num reads");
  assert_str_equals("seq.2", bam1_qname(reads[0]), "qname");

  len = next_fragment(reads, samfile, 2);
  assert_equals(0, len, "Num reads");

  samclose(samfile);
}


void test_cigar_to_spans() {
  char *sam_filename = "work/RUM.sam";
  samfile_t *samfile = samopen(sam_filename, "r", NULL);  
  bam1_t *rec = bam_init1();

  struct SpanAssertion {
    int read_num;
    int num_spans;
    struct Span spans[10];
  };

  struct SpanAssertion cases[] = {
    { 102, 1, { { 12465667, 12465724 } } },
    { 104, 1, { { 2095233, 2095289 } } },
    { 128, 1, { { 152316, 152373 } } },
    { 162, 1, { { 14232813, 14232886 } } },
    { 172, 2, { { 3619619, 3619627 },
                { 3619984, 3620048  } } },
    { 642, 1, { { 15291546, 15291622 } } },
    { 670, 2, { { 3950665, 3950724 },
                { 3951436, 3951453 } } }
  };

  int num_cases = sizeof(cases) / sizeof(struct SpanAssertion);
  int read_num = 0;
  int case_num = 0;
  CigarCursor curs;
  while (case_num < num_cases &&
         samread(samfile, rec) > 0) {

    if (cases[case_num].read_num == read_num) {

      int num_spans = cases[case_num].num_spans;
      Span *span;

      init_cigar_cursor(&curs, rec);

      for (span = cases[case_num].spans; span < cases[case_num].spans + num_spans; span++) {

        assert_equals(1, next_span(&curs), "Should have found a span");
        assert_equals(span->start, curs.start, "Start");
        assert_equals(span->end, curs.end, "End");
      }

      assert_equals(0, next_span(&curs), "No more spans");

      case_num++;
    }
    read_num++;
  }

}

void test_cigar_to_spans2() {

  char *sam_filename = "testdata/cigar_bug.sam";
  samfile_t *samfile = samopen(sam_filename, "r", NULL);  
  bam1_t *rec = bam_init1();

  struct SpanAssertion {
    int read_num;
    int num_spans;
    struct Span spans[10];
  };

  struct SpanAssertion cases[] = {
    { 0, 2, { { 46383142, 46383163 },
              { 46384677, 46384749 } } }
  };

  int num_cases = sizeof(cases) / sizeof(struct SpanAssertion);
  int read_num = 0;
  int case_num = 0;
  CigarCursor curs;
  while (case_num < num_cases &&
         samread(samfile, rec) > 0) {

    if (cases[case_num].read_num == read_num) {

      int num_spans = cases[case_num].num_spans;
      Span *span;

      init_cigar_cursor(&curs, rec);

      for (span = cases[case_num].spans; span < cases[case_num].spans + num_spans; span++) {

        if (next_span(&curs)) {
          assert_equals(span->start, curs.start, "Start");
          assert_equals(span->end, curs.end, "End");
        }
        
      }

      assert_equals(0, next_span(&curs), "No more spans");

      case_num++;
    }
    read_num++;
  }

}

int main(int argc, char **argv) {

  test_cigar_to_spans();
  test_cigar_to_spans2();
  test_next_fragment_paired();
  test_next_fragment_single();
  check_results();
  return 0;
}
