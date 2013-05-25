#include <stdio.h>
#include "geneindex.h"

int main (int argc, char **argv) {
  struct ExonDB exondb;

  parse_gtf_file(&exondb, "testdata/arabidopsis.gtf");
  index_exons(&exondb);  

}
