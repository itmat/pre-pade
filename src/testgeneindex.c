#include <stdio.h>
#include "quant.h"
#include "testutils.h"

void test_create_index() {
  struct GeneModel gm;

  load_model(&gm, "testdata/arabidopsis.gtf");

  assert_equals(13214, gm.exons.len, "Number of exons loaded");
  
  struct Region *exon = gm.exons.items;
  assert_str_equals("1", exon->chrom, "First exon chromosome");
  assert_equals(28692192, exon->start, "First exon start");
  assert_equals(28692362, exon->end, "First exon end");
  assert_str_equals("protein_coding", exon->source, "First exon source");
  assert_str_equals("exon", exon->feature, "First exon feature");

  int chroms_decrease = 0;
  int start_gte_end = 0;
  int entries_decrease = 0;
  
  struct IndexEntry *entry = gm.exons.index;

  for (entry = gm.exons.index + 1; entry < gm.exons.index + gm.exons.index_len; entry++) {
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

  struct RegionCursor cursor;
  int flags;
  //  search_exons(&cursor, &gm, "chrfoobar", 0, 0, 0);
  // assert_exon_ptr_equals(NULL, next_exon(&cursor, &flags), "Unknown chrom");

  //search_exons(&cursor, &gm, "1", 0, 10, 0);
  //assert_exon_ptr_equals(NULL, next_exon(&cursor, &flags), 
  //                         "Not found");

  search_exons(&cursor, &gm.exons, "1", 28692192, 28692362, 0);

  exon = next_exon(&cursor, &flags);

  assert_str_equals("1", exon->chrom, "Chromosome");

  assert_str_equals("AT1G76510.1", gm.transcripts[3].id, "Transcript id");
  assert_equals(14, gm.transcripts[3].exons_len, "Num exons");
  assert_equals(3, gm.transcripts[3].exons[11]->exon_number, "Region number");
  assert_equals((size_t)(gm.transcripts + 3), (size_t)(gm.transcripts[3].exons[11]->transcript), "Region's pointer to transcript");

}

void test_compare_index_entry() {
  struct IndexEntry entry = { "chr2", 100, 200 };

  struct IndexEntry left_chrom  = { "chr1", 100 };
  struct IndexEntry left        = { "chr2", 50 };
  struct IndexEntry inside      = { "chr2", 150 };
  struct IndexEntry right       = { "chr2", 250 };
  struct IndexEntry right_chrom = { "chr3", 100 };
  
  assert_equals(-1, cmp_index_entry(&left_chrom, &entry), "Left chrom");
  assert_equals(-1, cmp_index_entry(&left, &entry), "Left");
  assert_equals(0, cmp_index_entry(&inside, &entry), "Inside");
  assert_equals(1, cmp_index_entry(&right, &entry), "Right");
  assert_equals(1, cmp_index_entry(&right_chrom, &entry), "Right chrom");

}

void test_compare_exon() {
  struct Region e;
  e.chrom = "chr1";
  e.start = 100;
  e.end   = 200;

  struct RegionCompTestCase {
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

  int n = sizeof(cases) / sizeof(struct RegionCompTestCase);
  int i;

  for (i = 0; i < n; i++) {
    struct RegionCompTestCase *tc = cases + i;
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
  assert_equals(17, exon_number, "Region number");
}

void test_exon_matches() {
  RegionMatches ms;
  Region a;
  Region b;
  Region c;
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

  assert_str_equals("a", ms.items[0].region->chrom, "Chromosome");
  assert_str_equals("b", ms.items[1].region->chrom, "Chromosome");
  assert_str_equals("c", ms.items[2].region->chrom, "Chromosome");

  consolidate_exon_matches(&ms);
  assert_equals(3, ms.len, "Len");

  RegionMatch *ma = NULL, *mb = NULL, *mc = NULL;
  int i;

  for (i = 0; i < 3; i++) {
    RegionMatch *m = ms.items + i;
    char *chrom = m->region->chrom;
    
    if      (!strcmp(chrom, "a")) ma = m;
    else if (!strcmp(chrom, "b")) mb = m;
    else if (!strcmp(chrom, "c")) mc = m;
  }

  assert_not_null(ma, "ma");
  assert_not_null(mb, "mb");
  assert_not_null(mc, "mc");

  assert_str_equals("a", ma->region->chrom, "chrom a");
  assert_str_equals("b", mb->region->chrom, "chrom b");
  assert_str_equals("c", mc->region->chrom, "chrom c");

  assert_equals(0, ma->conflict, "conflict a");
  assert_equals(1, mb->conflict, "conflict b");
  assert_equals(3, mc->conflict, "conflict c");
}

void test_matches_junction() {
  GeneModel gm;
  load_model(&gm, "testdata/arabidopsis.gtf");

  const int j_start = 28696939;
  const int j_end   = 28697164;

  assert_str_equals("AT1G76490.1", gm.transcripts[1].id, "Transcript id\n");
  assert_equals(j_start, gm.transcripts[1].exons[0]->end, "Intron start\n");
  assert_equals(j_end, gm.transcripts[1].exons[1]->start, "Intron end\n");
  Region *exon = gm.transcripts[1].exons[0];

  Region *e2 = next_exon_in_transcript(exon);
  assert_equals(2, e2->exon_number, "Region number 2");
  Region *e3 = next_exon_in_transcript(e2);
  assert_equals(3, e3->exon_number, "Region number 3");
  Region *e4 = next_exon_in_transcript(e3);
  assert_equals(4, e4->exon_number, "Region number 4");
  Region *e5 = next_exon_in_transcript(e4);
  assert_equals(0, e5, "No exon 5");

  Span min_match[] = { { j_start - 10, j_start },
                       { j_end, j_end + 10 } };

  assert_equals(1, matches_junction(exon, min_match, 2, 0, 8), "Matches fwd only");
  assert_equals(1, matches_junction(exon, min_match, 2, 0, 10), "Barely matches fwd only");
  assert_equals(0, matches_junction(exon, min_match, 2, 0, 11), "Not enough overlap");

  assert_equals(1, matches_junction(exon, min_match, 0, 2, 8), "Matches rev only");
  assert_equals(1, matches_junction(exon, min_match, 0, 2, 10), "Barely matches rev only");
  assert_equals(0, matches_junction(exon, min_match, 0, 2, 11), "Not enough overlap");

  Span match_at_tail[] = { 
    { 0, 10 },
    { j_start - 10, j_start },
    { j_end, j_end + 10 } };

  assert_equals(1, matches_junction(exon, match_at_tail, 3, 0, 8), "Matches fwd only");
  assert_equals(1, matches_junction(exon, match_at_tail, 0, 3, 8), "Matches rev only");

  assert_equals(1, matches_junction(exon, match_at_tail, 1, 2, 8), "1 fwd, 2 rev");
  assert_equals(0, matches_junction(exon, match_at_tail, 2, 1, 8), "1 fwd, 2 rev, junction in unread portion of fragment");

  Span match_at_head[] = { 
    { j_start - 10, j_start },
    { j_end, j_end + 10 },
    { j_end + 1000, j_end + 2000 } };
  assert_equals(1, matches_junction(exon, match_at_head, 3, 0, 8), "Matches fwd only");
  assert_equals(1, matches_junction(exon, match_at_head, 0, 3, 8), "Matches rev only");
  assert_equals(0, matches_junction(exon, match_at_head, 1, 2, 8), "2 fwd, 1 rev, junction in unread portion");
  assert_equals(1, matches_junction(exon, match_at_head, 2, 1, 8), "2 fwd, 1 rev");

  Span no_spans[] = {};
  assert_equals(0, matches_junction(exon, no_spans, 0, 0, 0), "No spans");

  Span crosses_junction[] = { 
    { 0, 100 },
    { j_start - 10, j_end + 10 },
    { j_end + 1000, j_end + 2000 } };

  assert_equals(0, matches_junction(exon, no_spans, 0, 0, 0), "Crosses junction");

  
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

