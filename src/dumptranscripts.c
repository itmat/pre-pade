#include <stdio.h>
#include "quant.h"

void dump_transcripts(FILE *file, GeneModel *gm) {
  fprintf(file, "num_transcripts: %d\n", gm->num_transcripts);
  fprintf(file, "transcripts:\n");
  int i, j;

  for (i = 0; i < gm->num_transcripts; i++) {
    Transcript *t = gm->transcripts + i;
    fprintf(file, "- id: '%s'\n", t->id);
    fprintf(file, "  num_exons: %d\n", t->exons_len);
    fprintf(file, "  exons:\n");
    for (j = 0; j < t->exons_len; j++) {
      Region *e = t->exons[j];
      fprintf(file, "  - exon_number: %d\n", e->exon_number);
      fprintf(file, "    chrom: %s\n", e->chrom);
      fprintf(file, "    strand: %c\n", e->strand);
      fprintf(file, "    start: %d\n", e->start);
      fprintf(file, "    end: %d\n", e->end);
    }
    fprintf(file, "\n");
  
  }
}

void dump_regions(FILE *file, RegionList *list, char *type) {
  int i;
  fprintf(file, "num_%s: %d\n", type, list->len);
  fprintf(file, "%s:\n", type);
  
  for (i = 0; i < list->len; i++) {
    Region *r = list->items + i;
    fprintf(file, "- chrom: %s\n", r->chrom);
    fprintf(file, "  start: %d\n", r->start);
    fprintf(file, "  end: %d\n",   r->end);
  }

  fprintf(file, "index_len: %d\n", list->index_len);
  fprintf(file, "index:\n");
  for (i = 0; i < list->index_len; i++) {
    IndexEntry *e = list->index + i;
    fprintf(file, "- chrom: %s\n", e->chrom);
    fprintf(file, "  start: %d\n", e->start);
    fprintf(file, "  end: %d\n", e->end);
    Region *r = e->region;
    fprintf(file, "  %s:\n", type);
    fprintf(file, "    chrom: %s\n", r->chrom);
    fprintf(file, "    start: %d\n", r->start);
    fprintf(file, "    end: %d\n", r->end);
  }
}

void dump_exons(FILE *file, GeneModel *gm) {
  dump_regions(file, &gm->exons, "exon");
}

void dump_introns(FILE *file, GeneModel *gm) {
  dump_regions(file, &gm->introns, "exon");
}


void dump(char *filename, GeneModel *gm, void (*dumper_fn)(FILE *, GeneModel *), char *type) {

  fprintf(stderr, "Dumping %s to %s\n", type, filename);

  FILE *file = fopen(filename, "w");
  if (!file) {
    perror(filename);
    exit(1);
  }

  dumper_fn(file, gm);
  
  fclose(file);
}

int main(int argc, char **argv) {
  
  if (argc != 3) {
    fprintf(stderr, "Usage: %s GTF_FILE PREFIX\n", argv[0]);
    return 1;
  }
  GeneModel gm;

  load_model(&gm, argv[1]);
  dump("transcripts.yaml", &gm, dump_transcripts, "transcripts");
  dump("exons.yaml",       &gm, dump_exons,       "exons");
  dump("introns.yaml",     &gm, dump_introns,     "introns");

  return 0;
}
