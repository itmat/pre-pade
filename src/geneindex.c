#include <stdio.h>
#include <strings.h>
#include <stdlib.h>

#include "geneindex.h"

int print_exon(struct Exon *exon) {
  printf("%s:%d-%d\n", exon->chrom, exon->start, exon->end);
}

int parse_gtf_file(struct ExonDB *exondb, char *filename) {

  printf("Loading GTF file %s\n", filename);

  FILE *file = fopen(filename, "r");
  if (!file) {
    perror(filename);
  }

  int exons_len = 0;
  int exons_cap = 1000;
  struct Exon *exons = calloc(exons_cap, sizeof(struct Exon));

  while (1) {
    char *line = NULL;
    size_t linecap = 0;
    // TODO: Free line
    int linelen = getline(&line, &linecap, file);
    if (linelen < 0) {
      break;
    }

    if (exons_len == exons_cap) {
      struct Exon *old_exons = exons;
      int old_exons_cap = exons_cap;
      printf("Growing exons to %d\n", exons_cap);

      exons_cap *= 2;

      exons = calloc(exons_cap, sizeof(struct Exon));
      memcpy(exons, old_exons, old_exons_cap * sizeof(struct Exon));
      free(old_exons);
    }

    const int num_fields = 9;
    char *fields[num_fields];
    char *tok = line;
    
    int i;
    for (i = 0; i < num_fields; i++) {
      fields[i] = tok;
      if (i < num_fields - 1) {
        char *end = index(tok, '\t');
        *end = 0;
        tok = end + 1;
      }
    }
  
    struct Exon *exon = exons + exons_len;

    exon->chrom   = fields[0];
    exon->source  = fields[1];
    exon->feature = fields[2];
    exon->start   = atoi(fields[3]);
    exon->end     = atoi(fields[4]);
    
    exons_len++;
  }

  exondb->exons_len = exons_len;
  exondb->exons_cap = exons_cap;
  exondb->exons = exons;

}

int cmp_exons_by_end(struct Exon *a, struct Exon *b) {
  int str = strcmp(a->chrom, b->chrom);
  if (str) {
    return str;
  }
  return a->end - b->end;
}

int index_exons(struct ExonDB *exondb) {
  printf("Indexing %d exons\n", exondb->exons_len);
  int n = exondb->exons_len;
  int i;

  printf("Sorting exons by start pos\n");
  struct Exon *exons = exondb->exons;
  qsort(exons, n, sizeof(struct Exon), cmp_exons_by_end);

  int min_start = 0;
  char *chrom = NULL;

  for (i = n - 1; i >= 0; i--) {

    if (!chrom || strcmp(chrom, exons[i].chrom)) {
      chrom = exons[i].chrom;
      min_start = exons[i].start;
    }

    else if (exons[i].start < min_start) {
      min_start = exons[i].start;
    }
    exons[i].min_start = min_start;
  }

}

/*
 * Returns a pointer to the first exon whose end is greater than my
 * start.
 */
struct Exon * search_exons(struct ExonDB *exondb, char *chrom, int start, int end) {

  struct Exon *p = exondb->exons;
  struct Exon *q = p + exondb->exons_len - 1;
  struct Exon *e;

  while (p < q) {
    e = p + (q - p) / 2;
    
    // If this exon ends before I start, eliminate it and all exons to
    // the left of it
    if (e->end <= start) 
      p = e + 1;

    // If this exon ends after I start and either it's the first exon
    // or the one to the left ends before I start, return this exon.
    else if (e == exondb->exons ||
             (e - 1)->end <= start)
      return e;

    // Otherwise eliminate this exon and all those to the right of it
    else
      q = e - 1;
  }

  printf(stderr, "I SHOULD NEVER GET TO HERE");
  return NULL;
}

