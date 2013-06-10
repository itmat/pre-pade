#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <regex.h>
#include "samutils.h"
#include "geneindex.h"


void print_exon(Exon *exon) {
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

void parse_gtf_file(ExonDB *exondb, char *filename) {

  LOG_INFO("Loading GTF file %s\n", filename);

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
      LOG_DEBUG("Growing exons to %d\n", exons_cap);

      exons_cap *= 2;

      exons = calloc(exons_cap, sizeof(Exon));
      memcpy(exons, old_exons, old_exons_cap * sizeof(Exon));
      free(old_exons);
    }

    const int num_fields = 9;

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

      case 6:
        switch (strlen(tok) > 0 ? tok[0] : 0) {
        case '+': exon->strand = '+'; break;
        case '-': exon->strand = '-'; break;
        default: 
          fprintf(stderr, "Warning: Couldn't parse strand from %s\n", tok);
        }
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

  if (line) 
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

int cmp_transcript(Transcript *a, Transcript *b) {
  return strcmp(a->id, b->id);
}

int cmp_id_to_transcript(char *id, Transcript *b) {
  return strcmp(id, b->id);
}

int cmp_exon_ptr_end(Exon **a, Exon **b) {
  int str = strcmp((*a)->chrom, (*b)->chrom);
  if (str) {
    return str;
  }
  return (*a)->end - (*b)->end;
}

void add_transcripts(ExonDB *db) {
  LOG_DEBUG("Adding transcripts%s\n", "");

  Exon *exons = db->exons.items;
  int i, n = db->exons.len;
  Transcript *p, *q;

  // First make a list of all the transcript IDs from our exons  
  Transcript *transcripts = calloc(n, sizeof(Transcript));
  for (i = 0; i < n; i++) {
    transcripts[i].id = exons[i].transcript_id;
    transcripts[i].exons_cap = 1;
  }

  // Now sort the transcripts and go through them removing all duplicates and accumulating the exon counts

  qsort(transcripts, n, sizeof(Transcript), ( int (*)(const void *, const void*) ) cmp_transcript);

  for (p = transcripts, q = p + 1; q < transcripts + n; q++) {
    // If p and q have the same transcript id, just increment the
    // exons_cap counter for p and let q advance.
    if ( !strcmp(p->id, q->id) ) {
      p->exons_cap++;
    }

    // Otherwise we've found a new transcript. Advance p, and set its values to q's.
    else {
      memcpy(++p, q, sizeof(Transcript));
    }
  }

  // p is now the end of the real list of transcripts, so set the count appropriately
  int num_transcripts = p + 1 - transcripts;
  LOG_DEBUG("Found %d transcripts\n", num_transcripts);
  // Allocate space for the exon pointers in each transcript
  for (i = 0; i < num_transcripts; i++) {
    transcripts[i].exons = calloc(transcripts[i].exons_cap, sizeof(Exon*));
  }

  // Now for each exon, find its transcript, set its pointer to that
  // transcript, and add the exon to that transcript's list.
  for (i = 0; i < n; i++) {

    Exon *e = exons + i;

    Transcript *t = bsearch(e->transcript_id, transcripts, num_transcripts, 
                            sizeof(Transcript), 
                            (int (*) (const void *, const void *))cmp_id_to_transcript);
    e->transcript = t;
    t->exons[t->exons_len++] = e;
  }

  for (i = 0; i < num_transcripts; i++) {
    qsort(transcripts[i].exons, transcripts[i].exons_cap, sizeof(Exon*), ( int (*)(const void *, const void*) ) cmp_exon_ptr_end);
  }

  db->num_transcripts = num_transcripts;
  db->transcripts = transcripts;
}


void index_exons(ExonDB *exondb) {
  LOG_INFO("Indexing %d exons\n", exondb->exons.len);
  int n = exondb->exons.len;

  Exon *exons = exondb->exons.items;
  qsort(exons, n, sizeof(Exon), ( int (*)(const void *, const void*) ) cmp_exons_by_end);

  int min_start = 0;
  char *chrom = NULL;

  LOG_INFO("Finding minimum start positions for %d exons\n", exondb->exons.len);

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

  struct ExonIndexEntry *index = calloc(n, sizeof(ExonIndexEntry));
  ExonIndexEntry *entry = index;

  for (exon = exons; exon < exons + n; exon++) {

    // If it's the first exon for this chromosome, the start
    // position is the start of the chromosome.
    if (exon == exons || strcmp(exon->chrom, (exon - 1)->chrom)) {
      entry->chrom = exon->chrom;
      entry->start = 0;
      entry->end   = exon->end;
      entry->exon = exon;
      entry++;
    }

    // Otherwise it's the end of the last exon.
    else if (exon->end != (exon - 1)->end) {
      entry->chrom = exon->chrom;
      entry->start = (entry - 1)->end;
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

void init_exon_matches(ExonMatches *matches) {
  matches->cap = 1;
  matches->len = 0;
  matches->items = calloc(matches->cap, sizeof(ExonMatch));
}

void add_match(ExonMatches *matches, Exon *exon, int overlap, int conflict) {

  if (matches->len == matches->cap) {

    ExonMatch *old = matches->items;
    int old_cap = matches->cap;
    
    matches->cap *= 2;
    
    matches->items = calloc(matches->cap, sizeof(ExonMatch));
    memcpy(matches->items, old, old_cap * sizeof(ExonMatch));
    free(old);
  }

  ExonMatch *m = matches->items + matches->len++;
  m->exon = exon;
  m->overlap = overlap;
  m->conflict = conflict;

};

Exon *next_exon_in_transcript(Exon *e) {

  Transcript *t = e->transcript;

  int i;
  for (i = 0; i + 1 < t->exons_len; i++) {
    if (t->exons[i] == e) {
      return t->exons[i + 1];
    }
  }

  return NULL;
}



void find_candidates(ExonMatches *matches, ExonDB *db, char *ref,
                    Span *spans, int num_fwd_spans, int num_rev_spans) {

  matches->len = 0;
  int num_spans = num_fwd_spans + num_rev_spans;

  Span *first_fwd_span = num_fwd_spans ? spans                     : NULL;
  Span *last_fwd_span  = num_fwd_spans ? spans + num_fwd_spans - 1 : NULL;
  Span *first_rev_span = num_rev_spans ? spans + num_fwd_spans     : NULL;
  Span *last_rev_span  = num_rev_spans ? spans + num_spans - 1     : NULL;
  Span *span;

  for (span = spans; span < spans + num_spans; span++) {

    LOG_TRACE("  Looking at span %d-%d\n", span->start, span->end);

    struct Exon *exon;
    ExonCursor exon_curs;
    search_exons(&exon_curs, db, ref, span->start, span->end, ALLOW_ALL);

    int flags = 0;

    while ((exon = next_exon(&exon_curs, &flags))) {      

      LOG_TRACE("    Looking at exon %d-%d, flags are %d\n", exon->start, exon->end, flags);      
      int conflict = 0;

      // If it crosses either the start or end of the exon, we can't
      // count it.
      if (flags & (CROSS_EXON_START | CROSS_EXON_END))
        conflict = 1;

      // If the span starts in the exon, then it can only be a match
      // if it's the first span of either the forward or reverse read
      if ( ( flags & START_IN_EXON ) && 
           ! ( span == first_fwd_span || span == first_rev_span) )
        conflict = 1;

      // If the span ends in the exon, then it can only be a match
      // if it's the last span of either the forward or reverse read
      else if ( (flags & END_IN_EXON) &&
                ! ( span == last_fwd_span || span == last_rev_span ) )
        conflict = 1;

      add_match(matches, exon, 1, conflict);

    }

  }

  consolidate_exon_matches(matches);
}

int cmp_match_by_exon(ExonMatch *a, ExonMatch *b) {
  return b->exon - a->exon;
}

void consolidate_exon_matches(ExonMatches *matches) {
  qsort(matches->items, matches->len, sizeof(ExonMatch), ( int (*)(const void *, const void*) ) cmp_match_by_exon);  

  ExonMatch *p   = matches->items;
  ExonMatch *q   = matches->items + 1;
  ExonMatch *end = matches->items + matches->len;

  while (q < end) {
    if (cmp_match_by_exon(p, q)) {
      p++;
      p->exon = q->exon;
      p->overlap = q->overlap;
      p->conflict = q->conflict;
    }
    else {
      matches->len--;
      p->overlap += q->overlap;
      p->conflict += q->conflict;
    }
    q++;
  }

}

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
  LOG_TRACE("  Bsearch came back with %s:%d-%d\n",
            entry ? entry->exon->chrom : "",
            entry ? entry->exon->start : 0,
            entry ? entry->exon->end : 0);
            
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


int matches_junction(Exon *left, Span *spans, int num_fwd_spans, int num_rev_spans, int min_overlap) {

  Exon *right = next_exon_in_transcript(left);

  Span *last_fwd_span = spans + num_fwd_spans - 1;
  Span *last_rev_span = spans + num_fwd_spans + num_rev_spans - 1;

  // If there's no exon to the right of this one in the transcript,
  // then there's no junction to match to.
  if (!right)
    return 0;

  // Find the first span for which the end of the span is greater than or
  // equal to the end of the left exon
  Span *span = spans;

  for (span = spans; span < last_rev_span; span++) {

    // If the span's end is still to the left of the left exon's end,
    // then continue to the next span.
    if (span->end < left->end)
      continue;

    // If it's to the right of the left exon's end, then we can't
    // confirm the junction with a gap in this read.
    else if (span->end > left->end) 
      return 0;

    // Otherwise the end of this span matches the end of the exon. If
    // it's the last span in the forward read, than we can't use the
    // gap after it to confirm the junction, because it may simply be
    // that the end of the forward read coincidentally matches the end
    // of the exon.
    if (span == last_fwd_span)
      return 0;
    
    // Now we know we have a span with a gap to the right of it, where
    // the right edge of the span matches the right edge of the
    // exon. If the left edge of the next span matches the left edge
    // of the next exon, and both this span and the next are at least
    // as long as the min_overlap, then we have a hit. Otherwise we don't.
    return ((span+1)->start == right->start &&
            span->end - span->start >= min_overlap &&
            (span+1)->end - (span+1)->start >= min_overlap);
  }

  return 0;
}


int load_model(ExonDB *db, char *filename) {
  parse_gtf_file(db, filename);
  index_exons(db);
  add_transcripts(db);
  return 0;
}
