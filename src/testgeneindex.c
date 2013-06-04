#include <stdio.h>
#include "geneindex.h"
#include "testutils.h"

void test_create_index() {
  struct ExonDB db;
  parse_gtf_file(&db, "testdata/arabidopsis.gtf");
  index_exons(&db);  

  assert_equals(13214, db.exons_len, "Number of exons loaded");

  struct Exon *exon = db.exons;
  assert_str_equals("1", exon->chrom, "First exon chromosome");
  assert_equals(28692193, exon->start, "First exon start");
  assert_equals(28692362, exon->end, "First exon end");
  assert_str_equals("protein_coding", exon->source, "First exon source");
  assert_str_equals("exon", exon->feature, "First exon feature");

  int i = 1;
  int chroms_decrease = 0;
  int start_gte_end = 0;
  int entries_decrease = 0;
  
  struct ExonIndexEntry *entry = db.index;
  printf("Index len is %d\n", db.index_len);

  for (entry = db.index + 1; entry < db.index + db.index_len; entry++) {
    if (entry->start >= entry->end)
      start_gte_end++;

    int chrom_cmp = strcmp((entry - 1)->chrom, entry->chrom);

    if (chrom_cmp == 0 && (entry - 1)->end > entry->start)
      entries_decrease++;

    else if (chrom_cmp > 0) 
      chroms_decrease++;
  }

  assert_equals(0, chroms_decrease, "Chromosomes decrease");
  assert_equals(0, start_gte_end, "Entries with start greater than end");
  assert_equals(0, entries_decrease, "Entries decrease");

  struct ExonCursor cursor;
  int flags;
  //  search_exons(&cursor, &db, "chrfoobar", 0, 0, 0);
  // assert_exon_ptr_equals(NULL, next_exon(&cursor, &flags), "Unknown chrom");

  //search_exons(&cursor, &db, "1", 0, 10, 0);
  //assert_exon_ptr_equals(NULL, next_exon(&cursor, &flags), 
  //                         "Not found");

  search_exons(&cursor, &db, "1", 28692193, 28692362, 0);
  exon = next_exon(&cursor, &flags);

  assert_str_equals("1", exon->chrom, "Chromosome");
  print_exon(db.exons);
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
    { CROSS_EXON_END   | 
      START_IN_EXON,    "chr1", 150, 250, "Start in exon, cross exon end" },
    { CROSS_EXON_START | 
      END_IN_EXON,      "chr1",  50, 150, "Cross exon start, end in exon" },
    { START_IN_EXON    |
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

void test_parse_gtf_attr() {
  
  char *in = "gene_id \"ATMG00130\"; transcript_id \"ATMG00130.1\"; exon_number \"17\"; gene_name \"ORF121A\"; transcript_name \"ORF121A-201\"; seqedit \"false\";";

  char *value = parse_gtf_attr_str(in, "gene_id");
  assert_str_equals("ATMG00130", value, "Gene id");
  if (value) free(value);

  in = "gene_id=\"ATMG00130\";";
  value = parse_gtf_attr_str(in, "gene_id");
  assert_str_equals("ATMG00130", value, "Gene id");
  if (value) free(value);

  in = "gene_id = \"ATMG00130\";";
  value = parse_gtf_attr_str(in, "gene_id");
  assert_str_equals("ATMG00130", value, "Gene id");
  if (value) free(value);

  in = "gene_id ATMG00130;";
  value = parse_gtf_attr_str(in, "gene_id");
  assert_str_equals("ATMG00130", value, "Gene id");
  if (value) free(value);

  in = "gene_id ATMG00130";
  value = parse_gtf_attr_str(in, "gene_id");
  assert_equals(NULL, value, "Gene id");
  if (value) free(value);

  in = "exon_number \"17\";";
  int exon_number;
  int status = parse_gtf_attr_int(in, "exon_number", &exon_number);
  assert_equals(1, status, "Status");
  assert_equals(17, exon_number, "Exon number");
}



int main (int argc, char **argv) {
  test_create_index();
  test_compare_exon();
  test_compare_index_entry();
  test_parse_gtf_attr();
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
