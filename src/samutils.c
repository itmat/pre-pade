#include <stdio.h>
#include "sam.h"
#include "samutils.h"

int next_span(struct CigarCursor *c) {

  bam1_t *read = c->read;

  if (c->i == 0) {
    c->start = read->core.pos;
  }

  int n = read->core.n_cigar;
  uint32_t *cigar = bam1_cigar(read);

  if (c->i > n) 
    return 0;

  int op    = bam_cigar_op(cigar[c->i]);
  int oplen = bam_cigar_oplen(cigar[c->i]);  

  c->end = oplen;

}
