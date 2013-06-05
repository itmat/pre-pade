#include <stdio.h>
#include "sam.h"
#include "samutils.h"

int extract_spans(Span *spans, bam1_t *read, int n) {
  CigarCursor c;
  init_cigar_cursor(&c, read);
  int i = 0;
  while (next_span(&c)) {
    spans[i].start = c.start;
    spans[i].end   = c.end;
    i++;
  }
  return i;
}

int init_cigar_cursor(CigarCursor *c, bam1_t *read) {
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
    int is_first = reads[i]->core.flag & BAM_FREAD1;
    int is_last  = reads[i]->core.flag & BAM_FREAD2;
    
    num_reads++;

    // BAM_FREAD1 is supposed to indicate that it's the first read in
    // a template, and BAM_FREAD2 is supposed to indicate that it's
    // the last read in a template. So if BAM_FREAD1 is set, that
    // should mean that this is the second (and last) read. If neither
    // BAM_FREAD1 nor BAM_FREAD2 are set, we will assume that the
    // reads are not paired. This is consistent with SAM specification
    // v1.4-r985, section 1.4.2, bullet 2. In practice, this is how
    // RUM sets the flags. I'm not sure if other tools set it similarly.
    if (is_last || (!is_first || is_last))
      return num_reads;
  }

  if (!i) 
    return 0;

  fprintf(stderr, "Error: I found a fragment (%s) with %d reads\n", bam1_qname(reads[i]), i);
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
