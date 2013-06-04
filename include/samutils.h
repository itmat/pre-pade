#ifndef GENEINDEX_H
#define GENEINDEX_H

#include "sam.h"

struct CigarCursor {
  bam1_t *read;
  int i;
  int start;
  int end;
};

int init_cigar_cursor(struct CigarCursor *c, bam1_t *read);
int next_span(struct CigarCursor *c);

#endif
