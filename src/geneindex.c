#include <stdio.h>
#include <strings.h>
#include <stdlib.h>

enum Strand {
  NONE,
  UNKNOWN,
  FORWARD,
  REVERSE
};

struct Exon {
  char *gene_id;
  char *transcript_id;
  int exon_number;
  char *source;
  char *feature;
  char *chrom;
  enum Strand strand;
  int start;
  int end;
};

int print_exon(struct Exon *exon) {
  printf("%s:%d-%d\n", exon->chrom, exon->start, exon->end);
}

int parse_gtf_file(char *filename) {

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
    size_t linelen = getline(&line, &linecap, file);
    printf("Len is %d\n", linelen);
    if (linelen < 0) {
      printf("Done\n");
      return 0;
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
    print_exon(exon);
  }

}

int main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "Bad stuff\n");
    return -1;
  }

  parse_gtf_file(argv[1]);
  return 0;

}

