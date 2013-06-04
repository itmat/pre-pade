#include <stdio.h>

#include "samutils.h"
#include "geneindex.h"
#include "sam.h"

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
  
  fprintf(stderr, "Initializing bam\n");
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

    struct ExonCursor exon_curs;
    while (next_span(&span)) {

      struct Exon *exon;
      search_exons(&exon_curs, &db, ref, span.start, span.end, ALLOW_ALL);
      int flags = 0;
      while (exon = next_exon(&exon_curs, &flags)) {
        printf("%s\t", qname);
        printf("%d\t", ih);
        printf("%d\t", hi);
        printf("%s\t", ref);
        printf("%d\t", span.order);
        printf("%d\t", span.start);
        printf("%d\t", span.end);
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
       * exon start
       * exon end
        
       * crosses exon start
       * crosses exon end
       * starts in exon
       * ends in exon
        
       */

    }
  }

  bam_destroy1(rec);
  samclose(samfile);

  return 0;
}


