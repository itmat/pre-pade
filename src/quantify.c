#include <stdio.h>
#include "geneindex.h"

int main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "Bad stuff\n");
    return -1;
  }

  struct ExonDB exondb;

  parse_gtf_file(&exondb, argv[1]);
  index_exons(&exondb);
  return 0;

}

