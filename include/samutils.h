#ifndef SAMUTILS_H
#define SAMUTILS_H

#include "sam.h"

#define MAX_SPANS_PER_READ 1000

struct CigarCursor {
  bam1_t *read;
  int i;
  int start;
  int end;
  int order;
};

struct Span {
  int start;
  int end;
};

int next_fragment(bam1_t **reads, samfile_t *samfile, int n);
int init_cigar_cursor(struct CigarCursor *c, bam1_t *read);
int next_span(struct CigarCursor *c);

#endif
