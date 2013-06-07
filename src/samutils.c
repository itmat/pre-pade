#include <stdio.h>
#include <string.h>
#include "geneindex.h"

int extract_spans(Span *spans, bam1_t *read, int n) {
  CigarCursor c;
  init_cigar_cursor(&c, read);
  int i = 0;
  printf("in extract spans\n");
  while (next_span(&c)) {
    spans[i].start = c.start;
    spans[i].end   = c.end;
    i++;
  }
  return i;
}

void init_cigar_cursor(CigarCursor *c, bam1_t *read) {
  c->read = read;
  c->pos = read->core.pos;
  c->start = c->end = 0;
  c->i = 0;
}

int next_fragment(bam1_t **reads, samfile_t *samfile, int n) {

  int num_reads = 0;
  int i;

  if (samread(samfile, *reads) > 0) {
    num_reads++;
    
    if ( ! (reads[i]->core.flag & BAM_FPAIRED) )
      return 1;

    else if (samread(samfile, *(reads + 1)) > 0) {
      char *qname[] = { bam1_qname(reads[0]),  bam1_qname(reads[1]) };
      int      hi[] = { bam_aux2i(bam_aux_get(reads[0], "HI")),
                        bam_aux2i(bam_aux_get(reads[1], "HI")) };

      if (strcmp(qname[0], qname[1]) || hi[0] != hi[1]) {
        fprintf(stderr, "Error: paired flag is set, but next read is different. %s(%d) vs %s(%d) \n",
                qname[0], hi[0], qname[1], hi[1]);

        exit(1);
      }
      else {
        return 2;
      }
    }
    else {
      fprintf(stderr, "Error: paired flag is set, but there's no mate\n");
      exit(1);
    }
  }
  return 0;
}

int next_span(struct CigarCursor *c) {

  bam1_t *read = c->read;

  int n = read->core.n_cigar;
  uint32_t *cigar = bam1_cigar(read);

  int i;
  int found_match = 0;

  c->start = c->end = c->pos;

  for (i = c->i; i < n; i++) {
    int op    = bam_cigar_op(cigar[i]);
    int oplen = bam_cigar_oplen(cigar[i]);  
    int end_span = 0;

    switch (op) {
    
    case BAM_CMATCH:
      LOG_TRACE("  matching %d\n", oplen);
      c->end += oplen;
      c->pos += oplen;
      found_match = 1;
      break;
      
    case BAM_CINS:     
      LOG_TRACE("  inserting %d\n", oplen);
      break;

    case BAM_CDEL:
      LOG_TRACE("  deleting %d\n", oplen);
      c->end += oplen;
      c->pos += oplen;
      break;

    case BAM_CREF_SKIP:
      LOG_TRACE("  skipping %d\n", oplen);
      end_span = 1;
      c->pos += oplen;
      i++;
      break;

    case BAM_CSOFT_CLIP:
      LOG_TRACE("  soft-clipping %d\n", oplen);
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

  return found_match;
}
