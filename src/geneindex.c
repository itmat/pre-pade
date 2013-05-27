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

  printf("Finding minimum start positions\n");

  struct Exon *exon;

  for (exon = exons + n - 1; exon >= exons; exon--) {

    if (!chrom || strcmp(chrom, exon->chrom)) {
      chrom = exon->chrom;
      min_start = exon->start;
    }

    else if (exon->start < min_start) {
      min_start = exon->start;
    }

    exon->min_start = min_start;
  }

  printf("Building index\n");
  struct ExonIndexEntry *index = calloc(n, sizeof(struct ExonIndexEntry));
  struct ExonIndexEntry *entry = index;

  for (exon = exons; exon < exons + n; exon++) {

    // If it's the first exon for this chromosome, the start
    // position is the start of this exon.
    if (exon == exons || strcmp(exon->chrom, (exon - 1)->chrom)) {
      entry->chrom = exon->chrom;
      entry->start = exon->start;
      entry->end   = exon->end;
      entry++;
    }

    // Otherwise it's the end of the last exon.
    else if (exon->end != (exon - 1)->end) {
      int last_end = (entry - 1)->end;
      entry->chrom = exon->chrom;
      entry->start = last_end < exon->start ? exon->start : last_end;
      entry->end   = exon->end;
      entry++;
    }
    
  }

  exondb->index_len = entry - index;
  exondb->index = index;
}

int cmp_exon_end(struct Exon *exon, char *chrom, int pos) {
  int str = strcmp(exon->chrom, chrom);
  if (str) {
    return str;
  }
  return exon->end - pos;
}

/*
 * Returns a pointer to the first exon whose end is greater than my
 * start.
 */
int search_exons(struct ExonCursor *cursor,
                 struct ExonDB *exondb, char *chrom, int start, int end, int allow) {

  struct ExonIndexEntry key;
  key.chrom = chrom;
  key.start = key.end = start;

  struct ExonIndexEntry *entry = 
    bsearch(&key, exondb->index, exondb->index_len, 
            sizeof(struct ExonIndexEntry), cmp_index_entry);

  cursor->exondb = exondb;
  cursor->chrom = chrom;
  cursor->start = start;
  cursor->end = end;
  cursor->allow = allow;
  cursor->next = entry;
  
  return 0;
}

/* Compare the given exon to the specified range and return flags
 * indicating if and how it overlaps. 0 means it matches exactly.
 */
int cmp_exon(struct Exon *e, char *chrom, int start, int end) {

  int result = 0;

  if (strcmp(chrom, e->chrom))
    result |= WRONG_CHROMOSOME;

  // ...eeee
  //         rrrr...
  if ( e->end < start )
    result |= START_AFTER_EXON;

  //         eeee...
  // ...rrrr
  if ( end < e->start )
    result |= END_BEFORE_EXON;

  // ...eeee
  //   ...rrrr...
  if ( start < e->end && e->end < end )
    result |= CROSS_EXON_END;

  //      eeee...
  // ...rrrr...
  if ( start < e->start && e->start < end ) 
    result |= CROSS_EXON_START;

  //  ...eeee...
  //       rrrr...
  if ( e->start < start && start < e->end ) 
    result |= START_IN_EXON;

  //   ...eeee...
  // ...rrrr
  if ( e->start < end && end < e->end ) 
    result |= END_IN_EXON;

  return result;
}

struct Exon *finish_cursor(struct ExonCursor *cursor) {
  cursor->next = NULL;
  cursor->allow = 0;
  return NULL;
}

struct Exon *next_exon(struct ExonCursor *cursor, int *flags) {

  struct Exon   *e;
  struct ExonDB *db = cursor->exondb;
  struct Exon *last_exon = db->exons + db->exons_len;

  int allow = cursor->allow;
  int disallow = ~allow;

  while (cursor->next != NULL && cursor->next <= last_exon) {

    struct Exon *exon = cursor->next++;

    int cmp = cmp_exon(exon, cursor->chrom, cursor->start, cursor->end);

    if ( cmp & WRONG_CHROMOSOME ) {
      // If it's on the wrong chromosome, we're definitely done.
      return finish_cursor(cursor);
    }
   
    else if (e->min_start > cursor->end) {
      // At this point all subsequent exons will start after me, so
      // we're done.
      return finish_cursor(cursor);
    }

    else if (cmp & disallow) {
      // If this match would be disallowed according to the flags the
      // user passed in, just skip it.
      continue;
    }

    else {
      // Otherwise it's a match
      *flags = cmp;
      return exon;
    }
  }

  // If we get to this point, we've passed the last exon in the
  // database, so we're done.
  return finish_cursor(cursor);
}

int cmp_index_entry(struct ExonIndexEntry *key,
                    struct ExonIndexEntry *entry) {
  int cmp = strcmp(key->chrom, entry->chrom);
  int pos = key->start;

  if (cmp)
    return cmp;

  else if (pos < entry->start)
    return -1;

  else if (pos > entry->end)
    return 1;
  
  else
    return 0;
}
