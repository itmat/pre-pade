#include <stdio.h>

#include "samutils.h"
#include "geneindex.h"
#include "sam.h"

int main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "Usage: %s GTF_FILE SAM_FILE\n", argv[0]);
    return 1;
  }

  struct ExonDB db;
  parse_gtf_file(&db, "testdata/arabidopsis.gtf");

  char *sam_filename = argv[2];

  samfile_t *samfile = samopen(sam_filename, "r", NULL);
  
  printf("Initializing bam\n");
  bam1_t *rec = bam_init1();

  int count = 0;
  
  while (samread(samfile, rec) > 0) {

    char *qname = bam1_qname(rec);
    int pos = rec->core.pos;
    int n_cigar = rec->core.n_cigar;
    uint32_t *cigar = bam1_cigar(rec);
    int hi = bam_aux2i(bam_aux_get(rec, "HI"));
    int ih = bam_aux2i(bam_aux_get(rec, "IH"));
    int i;
    int tid = rec->core.tid;
    char *ref = samfile->header->target_name[tid];
    ref = ref ? ref : "";

    struct CigarCursor span;
    init_cigar_cursor(&span, rec);

    while (next_span(&span)) {
      printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
             qname,
             ih, hi,
             ref,
             span.order,
             span.start,
             span.end);

      //search_exons(struct ExonCursor *cursor,
      // &db, char *chrom, int start, int end, int allow);
      /*

       * fragment name
       * num alns
       * aln num

       * chrom
       * span num
       * span start
       * span end

       
        gene id
        transcript id
        num exons
        exon number
        exon start
        exon end
        
        crosses exon start
        crosses exon end
        starts in exon
        ends in exon
        
       */

    }
  }

  bam_destroy1(rec);
  samclose(samfile);

  return 0;
}


