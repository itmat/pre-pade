#include <stdio.h>

#include "samutils.h"
#include "geneindex.h"
#include "sam.h"

// Exon level counts
// Transcript level counts
// Junction level counts

int main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "Usage: %s GTF_FILE SAM_FILE\n", argv[0]);
    return 1;
  }

  char *gtf_filename = argv[1];
  char *sam_filename = argv[2];  

  struct ExonDB db;
  parse_gtf_file(&db, gtf_filename);
  index_exons(&db);

  samfile_t *samfile = samopen(sam_filename, "r", NULL);
  
  struct Span *read_spans = calloc(MAX_SPANS_PER_READ, sizeof(struct Span));

  bam1_t *reads[] = { bam_init1(), bam_init1() };
  int num_reads;

  ExonMatches matches;
  init_exon_matches(&matches);  

  while (num_reads = next_fragment(reads, samfile, 2)) {

    int num_fwd_spans = extract_spans(read_spans, reads[0], MAX_SPANS_PER_READ);
    int num_rev_spans = 0;

    if (num_reads == 2) {
      num_rev_spans = extract_spans(read_spans + num_fwd_spans, reads[1], MAX_SPANS_PER_READ - num_fwd_spans);
    }

    char *qname = bam1_qname(reads[0]);
    LOG_TRACE("On read %s\n", qname);
    int pos = reads[0]->core.pos;
    int n_cigar = reads[0]->core.n_cigar;
    uint32_t *cigar = bam1_cigar(reads[0]);
    int hi = bam_aux2i(bam_aux_get(reads[0], "HI"));
    int ih = bam_aux2i(bam_aux_get(reads[0], "IH"));
    int i;
    int tid = reads[0]->core.tid;
    char *ref = samfile->header->target_name[tid];
    ref = ref ? ref : "";

    find_candidates(&matches, &db, ref, read_spans, 
                    num_fwd_spans, num_rev_spans);

    /*
    struct ExonCursor exon_curs;
    int span_num;
    for (span_num = 0; span_num < num_fwd_spans + num_rev_spans; span_num++) {
      Span *span = read_spans + span_num;

      struct Exon *exon;
      search_exons(&exon_curs, &db, ref, span->start, span->end, ALLOW_ALL);
      int flags = 0;
      while (exon = next_exon(&exon_curs, &flags)) {
        printf("%s\t", qname);
        printf("%d\t", ih);
        printf("%d\t", hi);
        printf("%s\t", ref);
        printf("%d\t", span->start);
        printf("%d\t", span->end);
        printf("%s\t", exon->gene_id);
        printf("%s\t", exon->transcript_id);
        printf("%d\t", exon->exon_number);
        printf("%d\t", exon->start);
        printf("%d\t", exon->end);
        printf("%d\t%d\t%d\t%d\n",
               (flags & CROSS_EXON_START) > 0,
               (flags & CROSS_EXON_END) > 0,
               (flags & START_IN_EXON) > 0,
               (flags & END_IN_EXON) > 0
               );
      }
      
      }
    */
  }
  LOG_INFO("Cleaning up %s\n", "");

  bam_destroy1(reads[0]);
  bam_destroy1(reads[1]);
  //samclose(samfile);

  return 0;
}


