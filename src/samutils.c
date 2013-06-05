#include <stdio.h>
#include "sam.h"
#include "samutils.h"

int init_cigar_cursor(struct CigarCursor *c, bam1_t *read) {
  c->read = read;
  c->start = c->end = read->core.pos;
  c->i = 0;
  c->order = -1;
}

int next_fragment(bam1_t **reads, samfile_t *samfile, int n) {

  int num_reads = 0;
  int i;
  int is_last = 0;

  for (i = 0; i < n && samread(samfile, *(reads + i)) > 0; i++) {
    num_reads++;
    if (reads[i]->core.flag & BAM_FREAD2) 
      return num_reads;
  }

  fprintf(stderr, "Error: I found a fragment with more than two segments\n");
  exit(1);
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
    int end_span = 0;

    switch (op) {
    
    case BAM_CMATCH:
      c->end += oplen;
      found_match = 1;
      break;
      
    case BAM_CINS:     
      break;

    case BAM_CDEL:
      c->end += oplen;
      break;

    case BAM_CREF_SKIP:
      end_span = 1;
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

    if (end_span) 
      break;
  }
  c->i = i;
  c->order++;
  return found_match;
}
