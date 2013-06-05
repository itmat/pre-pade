#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <regex.h>

#include "geneindex.h"

int print_exon(Exon *exon) {
  printf("%s:%d-%d\n", exon->chrom, exon->start, exon->end);
}

int parse_gtf_attr(char *str, char *field_name, char **start, int *len) {
  char *p = strstr(str, field_name);

  if (!p)
    return 0;

  // Advance past the spaces or equal sign separating key from value
  p += strlen(field_name);
 
  while (*p == ' ') p++;
  if (*p == '=') p++;
  while (*p == ' ') p++;

  char *q = index(p, ';');

  if (!q)
    return 0;

  q--;

  while (p < q && *q == ' ') q--;

  if (*p == '"' && *q == '"') {
    p++;
    q--;
  }

  *start = p;
  *len = 1 + q - p;

  return *len;
}

int parse_gtf_attr_str(char *str, char *name, char **dest) {
  char *start;
  int len;
  if (parse_gtf_attr(str, name, &start, &len)) {
    *dest = strndup(start, len);
    return 1;
  }
  else {
    *dest = strdup("");
    return 0;
  }
}

int parse_gtf_attr_int(char *str, char *name, int *value) {
  char *start;
  int len;
  if (parse_gtf_attr(str, name, &start, &len)) {
    *value = atoi(start);
    return 1;
  }
  else {
    return 0;
  }
}

int parse_gtf_file(ExonDB *exondb, char *filename) {

  fprintf(stderr, "Loading GTF file %s\n", filename);

  FILE *file = fopen(filename, "r");
  if (!file) {
    perror(filename);
  }

  int exons_len = 0;
  int exons_cap = 1000;
  Exon *exons = calloc(exons_cap, sizeof(Exon));

  char *line = NULL;
  size_t linecap = 0;

  while (1) {

    int linelen = getline(&line, &linecap, file);
    if (linelen < 0) {
      break;
    }

    if (exons_len == exons_cap) {
      Exon *old_exons = exons;
      int old_exons_cap = exons_cap;
      fprintf(stderr, "Growing exons to %d\n", exons_cap);

      exons_cap *= 2;

      exons = calloc(exons_cap, sizeof(Exon));
      memcpy(exons, old_exons, old_exons_cap * sizeof(Exon));
      free(old_exons);
    }

    const int num_fields = 9;
    char *fields[num_fields];
    char *tok = line;    
    int i;
    Exon *exon = exons + exons_len;

    for (i = 0; i < num_fields; i++) {
      char *end = index(tok, i == num_fields - 1 ? 0 : '\t');

      switch (i) {
      case 0:
        exon->chrom = strndup(tok, end - tok);
        break;

      case 1:
        exon->source = strndup(tok, end - tok);
        break;

      case 2:
        exon->feature = strndup(tok, end - tok);
        break;

      case 3:
        exon->start = atoi(tok);
        break;

      case 4:
        exon->end = atoi(tok);
        break;

      case 8:

        if ( !parse_gtf_attr_str(tok, "gene_id", &exon->gene_id) ) {
          fprintf(stderr, "Warning: can't find gene_id for %s\n", line);
        }
        if ( !parse_gtf_attr_str(tok, "transcript_id", &exon->transcript_id) ) {
          fprintf(stderr, "Warning: can't find transcript_id for %s\n", line);
        }
        if ( !parse_gtf_attr_int(tok, "exon_number", &exon->exon_number) ) {
          fprintf(stderr, "Warning: can't find exon_number for %s\n", line);
        }

      }

      tok = end + 1;
    }
  
    exons_len++;
  }

  free(line);

  exondb->exons.len = exons_len;
  exondb->exons.cap = exons_cap;
  exondb->exons.items = exons;
}

int cmp_exons_by_end(Exon *a, Exon *b) {
  int str = strcmp(a->chrom, b->chrom);
  if (str) {
    return str;
  }
  return a->end - b->end;
}

int index_exons(ExonDB *exondb) {
  fprintf(stderr, "Indexing %d exons\n", exondb->exons.len);
  int n = exondb->exons.len;
  int i;

  fprintf(stderr, "Sorting exons by start pos\n");
  Exon *exons = exondb->exons.items;
  qsort(exons, n, sizeof(Exon), ( int (*)(const void *, const void*) ) cmp_exons_by_end);

  int min_start = 0;
  char *chrom = NULL;

  fprintf(stderr, "Finding minimum start positions\n");

  Exon *exon;

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

  fprintf(stderr, "Building index\n");
  struct ExonIndexEntry *index = calloc(n, sizeof(ExonIndexEntry));
  ExonIndexEntry *entry = index;

  for (exon = exons; exon < exons + n; exon++) {

    // If it's the first exon for this chromosome, the start
    // position is the start of this exon.
    if (exon == exons || strcmp(exon->chrom, (exon - 1)->chrom)) {
      entry->chrom = exon->chrom;
      entry->start = exon->start;
      entry->end   = exon->end;
      entry->exon = exon;
      entry++;
    }

    // Otherwise it's the end of the last exon.
    else if (exon->end != (exon - 1)->end) {
      int last_end = (entry - 1)->end;
      entry->chrom = exon->chrom;
      entry->start = last_end < exon->start ? exon->start : last_end;
      entry->end   = exon->end;
      entry->exon  = exon;
      entry++;
    }
    
  }

  exondb->index_len = entry - index;
  exondb->index = index;
}

int cmp_exon_end(Exon *exon, char *chrom, int pos) {
  int str = strcmp(exon->chrom, chrom);
  if (str) {
    return str;
  }
  return exon->end - pos;
}



/*
int find_candidates() {

  int num_spans = num_fwd_spans + num_rev_spans;

  Span *first_fwd_span = num_fwd_spans ? spans                     : NULL;
  Span *last_fwd_span  = num_fwd_spans ? spans + num_fwd_spans - 1 : NULL;
  Span *first_rev_span = num_rev_spans ? spans + num_fwd_spans     : NULL;
  Span *last_rev_span  = num_rev_spans ? spans + num_spans - 1     : NULL;

  for (span_num = 0; span_num < num_spans; span_num++) {
      Span *span = read_spans + span_num;

      struct Exon *exon;
      search_exons(&exon_curs, &db, ref, span->start, span->end, ALLOW_ALL);

      int flags = 0;
      
      Exon *allow;
      Exon *deny;

      while (exon = next_exon(&exon_curs, &flags)) {      

        int match = 1;

        // If it crosses either the start or end of the exon, we can't
        // count it.
        if (flags & (CROSS_EXON_START | CROSS_EXON_END))
          match = 0;
        
        // If the span starts in the exon, then it can only be a match
        // if it's the first span of either the forward or reverse read
        if ( ( flags & START_IN_EXON ) && 
             ! ( span == first_fwd_span || span == first_rev_span) )
          match = 0;

        // If the span ends in the exon, then it can only be a match
        // if it's the last span of either the forward or reverse read
        else if ( (flags & END_IN_EXON) &&
             ! ( span == last_fwd_span || span == last_rev_span ) )
          match = 0;

        if (match) {
          add_candidate();
        }

      }
  }
}

*/

/*
 * Returns a pointer to the first exon whose end is greater than my
 * start.
 */
int search_exons(ExonCursor *cursor,
                 ExonDB *exondb, char *chrom, int start, int end, int allow) {
  ExonIndexEntry key;
  key.chrom = chrom;
  key.start = key.end = start;
  ExonIndexEntry *entry = 
    bsearch(&key, exondb->index, exondb->index_len, 
            sizeof(ExonIndexEntry), 
            (int (*) (const void *, const void *))cmp_index_entry);
  cursor->exondb = exondb;
  cursor->chrom = chrom;
  cursor->start = start;
  cursor->end = end;
  cursor->allow = allow;

  cursor->next = entry ? entry->exon : NULL;
  
  return 0;
}

/* Compare the given exon to the specified range and return flags
 * indicating if and how it overlaps. 0 means it matches exactly.
 */
int cmp_exon(Exon *e, char *chrom, int start, int end) {

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

Exon *finish_cursor(ExonCursor *cursor) {
  cursor->next = NULL;
  cursor->allow = 0;
  return NULL;
}

Exon *next_exon(ExonCursor *cursor, int *flags) {

  ExonDB *db = cursor->exondb;
  Exon *last_exon = db->exons.items + db->exons.len;

  int allow = cursor->allow;
  int disallow = ~allow;

  while (cursor->next != NULL && cursor->next <= last_exon) {
    Exon *exon = cursor->next++;
    int cmp = cmp_exon(exon, cursor->chrom, cursor->start, cursor->end);

    if ( cmp & WRONG_CHROMOSOME ) {
      // fprintf(stderr, "  wrong chrom\n");   
      // If it's on the wrong chromosome, we're definitely done.
      return finish_cursor(cursor);
    }

    else if (exon->min_start > cursor->end) {
      // At this point all subsequent exons will start after me, so
      // we're done.
      //      fprintf(stderr, "  past\n");   
      return finish_cursor(cursor);
    }

    else if (cmp & disallow) {
      // fprintf(stderr, "  miss\n");   
      // If this match would be disallowed according to the flags the
      // user passed in, just skip it.
      continue;
    }

    else {
      //      fprintf(stderr, "Match\n");
      // Otherwise it's a match
      *flags = cmp;
      return exon;
    }

  }

  // If we get to this point, we've passed the last exon in the
  // database, so we're done.
  return finish_cursor(cursor);
}

int cmp_index_entry(ExonIndexEntry *key,
                    ExonIndexEntry *entry) {
  //    fprintf(stderr, "Comparing %s:%d-%d and %s:%d-%d\n",
  //           key->chrom, key->start, key->end,
  
  //            entry->chrom, entry->start, entry->end);
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
