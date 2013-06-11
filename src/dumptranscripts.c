#include <stdio.h>

#include "quant.h"

int main(int argc, char **argv) {
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s GTF_FILE\n", argv[0]);
    return 1;
  }
  ExonDB db;
  parse_gtf_file(&db, argv[1]);
  index_exons(&db);
  add_transcripts(&db);

  int i, j;

  printf("num_transcripts: %d\n", db.num_transcripts);
  printf("transcripts:\n");

  for (i = 0; i < db.num_transcripts; i++) {
    Transcript *t = db.transcripts + i;
    printf("- id: '%s'\n", t->id);
    printf("  num_exons: %d\n", t->exons_len);
    printf("  exons:\n");
    for (j = 0; j < t->exons_len; j++) {
      Exon *e = t->exons[j];
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
