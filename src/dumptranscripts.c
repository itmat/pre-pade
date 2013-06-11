#include <stdio.h>

#include "quant.h"

int main(int argc, char **argv) {
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s GTF_FILE\n", argv[0]);
    return 1;
  }
  GeneModel gm;
  load_model(&gm, argv[1]);

  int i, j;

  printf("num_transcripts: %d\n", gm.num_transcripts);
  printf("transcripts:\n");

  for (i = 0; i < gm.num_transcripts; i++) {
    Transcript *t = gm.transcripts + i;
    printf("- id: '%s'\n", t->id);
    printf("  num_exons: %d\n", t->exons_len);
    printf("  exons:\n");
    for (j = 0; j < t->exons_len; j++) {
      Region *e = t->exons[j];
      printf("  - exon_number: %d\n", e->exon_number);
      printf("    chrom: %s\n", e->chrom);
      printf("    strand: %c\n", e->strand);
      printf("    start: %d\n", e->start);
      printf("    end: %d\n", e->end);
    }
    printf("\n");
  
  }

  return 0;
}
