#include <stdio.h>
#include "sam.h"
#include "samutils.h"

int init_cigar_cursor(struct CigarCursor *c, bam1_t *read) {
  c->read = read;
  c->start = read->core.pos;
  c->i = 0;
  c->end = 0;
  printf("Pos is %d\n", read->core.pos);
}

int next_span(struct CigarCursor *c) {

  bam1_t *read = c->read;

  int n = read->core.n_cigar;
  uint32_t *cigar = bam1_cigar(read);

  int i;
  int found_match = 0;

  for (i = c->i; i < n; i++) {
    int op    = bam_cigar_op(cigar[i]);
    int oplen = bam_cigar_oplen(cigar[i]);  

    switch (op) {
    
    case BAM_CMATCH:
      c->end = c->start + oplen;
      found_match = 1;
      break;
      
    case BAM_CINS:     
    case BAM_CDEL:
    case BAM_CREF_SKIP:
      break;

    case BAM_CSOFT_CLIP:
      break; 
      
    case BAM_CHARD_CLIP:
    case BAM_CPAD: 
    case BAM_CEQUAL:
    case BAM_CDIFF:
    case BAM_CBACK:
      break;
    }
  }
  fprintf(stderr, "Here I am\n");
  c->i = i;
  return found_match;
}
