#include <stdio.h>
#include "geneindex.h"
#include "testutils.h"

void test_create_index() {
  struct ExonDB db;
  parse_gtf_file(&db, "testdata/arabidopsis.gtf");
  index_exons(&db);  
  add_transcripts(&db);

  assert_equals(13214, db.exons.len, "Number of exons loaded");

  struct Exon *exon = db.exons.items;
  assert_str_equals("1", exon->chrom, "First exon chromosome");
  assert_equals(28692193, exon->start, "First exon start");
  assert_equals(28692362, exon->end, "First exon end");
  assert_str_equals("protein_coding", exon->source, "First exon source");
  assert_str_equals("exon", exon->feature, "First exon feature");

  int chroms_decrease = 0;
  int start_gte_end = 0;
  int entries_decrease = 0;
  
  struct ExonIndexEntry *entry = db.index;

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

  assert_str_equals("AT1G76510.1", db.transcripts[3].id, "Transcript id");
  assert_equals(14, db.transcripts[3].exons_len, "Num exons");
  assert_equals(3, db.transcripts[3].exons[11]->exon_number, "Exon number");
  assert_equals((size_t)(db.transcripts + 3), (size_t)(db.transcripts[3].exons[11]->transcript), "Exon's pointer to transcript");

}

void test_compare_index_entry() {
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

void test_compare_exon() {
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

  char *value;
  parse_gtf_attr_str(in, "gene_id", &value);
  assert_str_equals("ATMG00130", value, "Gene id");
  free(value);

  in = "gene_id=\"ATMG00130\";";
  parse_gtf_attr_str(in, "gene_id", &value);
  assert_str_equals("ATMG00130", value, "Gene id");
  free(value);

  in = "gene_id = \"ATMG00130\";";
  parse_gtf_attr_str(in, "gene_id", &value);
  assert_str_equals("ATMG00130", value, "Gene id");
  free(value);

  in = "gene_id ATMG00130;";
  parse_gtf_attr_str(in, "gene_id", &value);
  assert_str_equals("ATMG00130", value, "Gene id");
  free(value);

  in = "gene_id ATMG00130";
  parse_gtf_attr_str(in, "gene_id", &value);
  assert_str_equals("", value, "Gene id");
  free(value);

  in = "exon_number \"17\";";
  int exon_number;
  int status = parse_gtf_attr_int(in, "exon_number", &exon_number);
  assert_equals(1, status, "Status");
  assert_equals(17, exon_number, "Exon number");
}

void test_exon_matches() {
  ExonMatches ms;
  Exon a;
  Exon b;
  Exon c;
  a.chrom = "a";
  b.chrom = "b";
  c.chrom = "c";
  init_exon_matches(&ms);

  add_match(&ms, &a, 1, 0);
  add_match(&ms, &b, 1, 0);
  add_match(&ms, &c, 1, 1);
  add_match(&ms, &a, 1, 0);
  add_match(&ms, &b, 1, 1);
  add_match(&ms, &c, 1, 1);
  add_match(&ms, &a, 1, 0);
  add_match(&ms, &b, 1, 0);
  add_match(&ms, &c, 1, 1);

  assert_equals(16, ms.cap, "Capacity");

  assert_str_equals("a", ms.items[0].exon->chrom, "Chromosome");
  assert_str_equals("b", ms.items[1].exon->chrom, "Chromosome");
  assert_str_equals("c", ms.items[2].exon->chrom, "Chromosome");

  consolidate_exon_matches(&ms);
  assert_equals(3, ms.len, "Len");

  ExonMatch *ma = NULL, *mb = NULL, *mc = NULL;
  int i;

  for (i = 0; i < 3; i++) {
    ExonMatch *m = ms.items + i;
    char *chrom = m->exon->chrom;
    
    if      (!strcmp(chrom, "a")) ma = m;
    else if (!strcmp(chrom, "b")) mb = m;
    else if (!strcmp(chrom, "c")) mc = m;
  }

  assert_not_null(ma, "ma");
  assert_not_null(mb, "mb");
  assert_not_null(mc, "mc");

  assert_str_equals("a", ma->exon->chrom, "chrom a");
  assert_str_equals("b", mb->exon->chrom, "chrom b");
  assert_str_equals("c", mc->exon->chrom, "chrom c");

  assert_equals(0, ma->conflict, "conflict a");
  assert_equals(1, mb->conflict, "conflict b");
  assert_equals(3, mc->conflict, "conflict c");
}

void test_matches_junction() {
  ExonDB db;
  parse_gtf_file(&db, "testdata/arabidopsis.gtf");
  index_exons(&db);  
  add_transcripts(&db);

  printf("%d: %d-%d\n", db.exons.items[0].exon_number, db.exons.items[0].start, db.exons.items[0].end);
  
}

int main (int argc, char **argv) {
  test_create_index();
  test_compare_exon();
  test_compare_index_entry();
  test_parse_gtf_attr();
  test_exon_matches();
  test_matches_junction();
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

