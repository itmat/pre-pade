#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <regex.h>
#include "quant.h"

void print_exon(Region *exon) {
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

void parse_gtf_file(GeneModel *gm, char *filename) {

  LOG_INFO("Loading GTF file %s\n", filename);

  FILE *file = fopen(filename, "r");
  if (!file) {
    perror(filename);
  }

  int exons_len = 0;
  int exons_cap = 1000;
  Region *exons = calloc(exons_cap, sizeof(Region));

  char *line = NULL;
  size_t linecap = 0;

  while (1) {

    int linelen = getline(&line, &linecap, file);
    if (linelen < 0) {
      break;
    }

    if (exons_len == exons_cap) {
      Region *old_exons = exons;
      int old_exons_cap = exons_cap;
      LOG_DEBUG("Growing exons to %d\n", exons_cap);

      exons_cap *= 2;

      exons = calloc(exons_cap, sizeof(Region));
      memcpy(exons, old_exons, old_exons_cap * sizeof(Region));
      free(old_exons);
    }

    const int num_fields = 9;

    char *tok = line;    
    int i;
    Region *exon = exons + exons_len;

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
        exon->start = atoi(tok) - 1;
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

  gm->exons.len = exons_len;
  gm->exons.cap = exons_cap;
  gm->exons.items = exons;
}

int cmp_exons_by_end(Region *a, Region *b) {
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

int cmp_exon_ptr_end(Region **a, Region **b) {
  int str = strcmp((*a)->chrom, (*b)->chrom);
  if (str) {
    return str;
  }
  return (*a)->end - (*b)->end;
}


void add_transcripts(GeneModel *gm) {
  LOG_DEBUG("Adding transcripts%s\n", "");

  Region *exons = gm->exons.items;
  int i, n = gm->exons.len;
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
    transcripts[i].exons = calloc(transcripts[i].exons_cap, sizeof(Region*));
  }

  // Now for each exon, find its transcript, set its pointer to that
  // transcript, and add the exon to that transcript's list.
  for (i = 0; i < n; i++) {

    Region *e = exons + i;

    Transcript *t = bsearch(e->transcript_id, transcripts, num_transcripts, 
                            sizeof(Transcript), 
                            (int (*) (const void *, const void *))cmp_id_to_transcript);
    e->transcript = t;
    t->exons[t->exons_len++] = e;
  }

  for (i = 0; i < num_transcripts; i++) {
    qsort(transcripts[i].exons, transcripts[i].exons_cap, sizeof(Region*), ( int (*)(const void *, const void*) ) cmp_exon_ptr_end);
  }

  gm->num_transcripts = num_transcripts;
  gm->transcripts = transcripts;
}

void add_introns(GeneModel *gm) {
  Transcript *t;
  int num_introns = 0;
  fprintf(stderr, "Adding introns\n");
  for (t = gm->transcripts; t < gm->transcripts + gm->num_transcripts; t++) {
    int introns_in_this_transcript = t->exons_len - 1;
    if (introns_in_this_transcript > 0) {
      num_introns += introns_in_this_transcript;
    }
  }
  LOG_DEBUG("Allocating %d introns\n", num_introns);

  Region *introns = calloc(num_introns, sizeof(Region));

  int j = 0;
  LOG_DEBUG("Populating introns%s\n", "");
  for (t = gm->transcripts; t < gm->transcripts + gm->num_transcripts; t++) {
    int i;
    for (i = 1; i < t->exons_len; i++) {
      introns[j].start   = t->exons[i-1]->end;
      introns[j].end     = t->exons[i]->start;
      introns[j].chrom   = t->exons[i]->chrom;
      introns[j].gene_id = t->exons[i]->gene_id;
      introns[j].transcript_id = t->exons[i]->transcript_id;
      introns[j].exon_number = t->exons[i-1]->exon_number;
      j++;
    }
  }
  fprintf(stderr, "J is %d\n", j);
  LOG_DEBUG("Sorting %d introns\n", num_introns);
  qsort(introns, num_introns, sizeof(Region), ( int (*)(const void *, const void*) ) cmp_exons_by_end);

  gm->introns.cap = num_introns;
  gm->introns.len = num_introns;
  gm->introns.items = introns;

  LOG_DEBUG("Indexing introns%s\n", "");
  index_regions(&gm->introns);

}

void index_regions(RegionList *exons) {
  LOG_INFO("Indexing %d exons\n", exons->len);
  int n = exons->len;

  qsort(exons->items, n, sizeof(Region), ( int (*)(const void *, const void*) ) cmp_exons_by_end);

  int min_start = 0;
  char *chrom = NULL;

  LOG_INFO("Finding minimum start positions for %d exons\n", exons->len);

  Region *exon;

  for (exon = exons->items + n - 1; exon >= exons->items; exon--) {

    if (!chrom || strcmp(chrom, exon->chrom)) {
      chrom = exon->chrom;
      min_start = exon->start;
    }

    else if (exon->start < min_start) {
      min_start = exon->start;
    }

    exon->min_start = min_start;
  }

  IndexEntry *index = calloc(n, sizeof(IndexEntry));
  IndexEntry *entry = index;

  for (exon = exons->items; exon < exons->items + n; exon++) {

    // If it's the first exon for this chromosome, the start
    // position is the start of the chromosome.
    if (exon == exons->items || strcmp(exon->chrom, (exon - 1)->chrom)) {
      entry->chrom = exon->chrom;
      entry->start = 0;
      entry->end   = exon->end;
      entry->region = exon;
      entry++;
    }

    // Otherwise it's the end of the last exon.
    else if (exon->end != (exon - 1)->end) {
      entry->chrom = exon->chrom;
      entry->start = (entry - 1)->end;
      entry->end   = exon->end;
      entry->region  = exon;
      entry++;
    }
    
  }

  exons->index_len = entry - index;
  exons->index = index;
}

int cmp_exon_end(Region *exon, char *chrom, int pos) {
  int str = strcmp(exon->chrom, chrom);
  if (str) {
    return str;
  }
  return exon->end - pos;
}

void init_exon_matches(RegionMatches *matches) {
  matches->cap = 1;
  matches->len = 0;
  matches->items = calloc(matches->cap, sizeof(RegionMatch));
}

void add_transcript(Transcript **transcripts, int *cap, int *len, 
                   Transcript *transcript) {

  if (*len == *cap) {
    Transcript **old = transcripts;
    int old_cap = *cap;
    
    (*cap) *= 2;
    transcripts = calloc(*cap, sizeof(Transcript*));
    memcpy(old, transcripts, old_cap * sizeof(Transcript*));
    free(old);
  }


  transcripts[(*len)++] = transcript;
}


void add_match(RegionMatches *matches, Region *exon, int overlap, int conflict) {
  LOG_TRACE("Adding match %s\n", "");
  if (matches->len == matches->cap) {

    RegionMatch *old = matches->items;
    int old_cap = matches->cap;
    
    matches->cap *= 2;
    
    matches->items = calloc(matches->cap, sizeof(RegionMatch));
    memcpy(matches->items, old, old_cap * sizeof(RegionMatch));
    free(old);
  }

  RegionMatch *m = matches->items + matches->len++;
  m->region = exon;
  m->overlap = overlap;
  m->conflict = conflict;
  LOG_TRACE("Done adding match %s\n", "");
};

Region *next_exon_in_transcript(Region *e) {

  Transcript *t = e->transcript;

  int i;
  for (i = 0; i + 1 < t->exons_len; i++) {
    if (t->exons[i] == e) {
      return t->exons[i + 1];
    }
  }

  return NULL;
}



void find_candidates(RegionMatches *matches, RegionList *list, char *ref,
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

    struct Region *exon;
    RegionCursor exon_curs;
    search_exons(&exon_curs, list, ref, span->start, span->end, ALLOW_ALL);

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

int cmp_match_by_exon(RegionMatch *a, RegionMatch *b) {
  return b->region - a->region;
}

void consolidate_exon_matches(RegionMatches *matches) {
  qsort(matches->items, matches->len, sizeof(RegionMatch), ( int (*)(const void *, const void*) ) cmp_match_by_exon);  

  RegionMatch *p   = matches->items;
  RegionMatch *q   = matches->items + 1;
  RegionMatch *end = matches->items + matches->len;

  while (q < end) {
    if (cmp_match_by_exon(p, q)) {
      p++;
      p->region = q->region;
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
int search_exons(RegionCursor *cursor,
                 RegionList *list, char *chrom, int start, int end, int allow) {
  IndexEntry key;
  key.chrom = chrom;
  key.start = key.end = start;
  IndexEntry *entry = 
    bsearch(&key, list->index, list->index_len, 
            sizeof(IndexEntry), 
            (int (*) (const void *, const void *))cmp_index_entry);
  LOG_TRACE("  Bsearch came back with %s:%d-%d\n",
            entry ? entry->region->chrom : "",
            entry ? entry->region->start : 0,
            entry ? entry->region->end : 0);
            
  cursor->list = list;
  cursor->chrom = chrom;
  cursor->start = start;
  cursor->end = end;
  cursor->allow = allow;

  cursor->next = entry ? entry->region : NULL;
  
  return 0;
}

/* Compare the given exon to the specified range and return flags
 * indicating if and how it overlaps. 0 means it matches exactly.
 */
int cmp_exon(Region *e, char *chrom, int start, int end) {

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

Region *finish_cursor(RegionCursor *cursor) {
  cursor->next = NULL;
  cursor->allow = 0;
  return NULL;
}

Region *next_exon(RegionCursor *cursor, int *flags) {

  RegionList *list = cursor->list;
  Region *last_exon = list->items + list->len;

  int allow = cursor->allow;
  int disallow = ~allow;

  while (cursor->next != NULL && cursor->next < last_exon) {
    Region *r = cursor->next++;
    int cmp = cmp_exon(r, cursor->chrom, cursor->start, cursor->end);

    if ( cmp & WRONG_CHROMOSOME ) {
      // fprintf(stderr, "  wrong chrom\n");   
      // If it's on the wrong chromosome, we're definitely done.
      return finish_cursor(cursor);
    }

    else if (r->min_start > cursor->end) {
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
      return r;
    }

  }

  // If we get to this point, we've passed the last exon in the
  // database, so we're done.
  return finish_cursor(cursor);
}

int cmp_index_entry(IndexEntry *key,
                    IndexEntry *entry) {

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


int matches_junction(Region *left, Span *spans, int num_fwd_spans, int num_rev_spans, int min_overlap) {
  
  Region *right = next_exon_in_transcript(left);

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
    else if (span->end > left->end)  {
      return 0;
    }

    // Otherwise the end of this span matches the end of the exon. If
    // it's the last span in the forward read, than we can't use the
    // gap after it to confirm the junction, because it may simply be
    // that the end of the forward read coincidentally matches the end
    // of the exon.
    if (span == last_fwd_span) {
      return 0;
    }
    
    // Now we know we have a span with a gap to the right of it, where
    // the right edge of the span matches the right edge of the
    // exon. If the left edge of the next span matches the left edge
    // of the next exon, and both this span and the next are at least
    // as long as the min_overlap, then we have a hit. Otherwise we don't.
    const int matches_right = (span+1)->start == right->start;
    const int left_overlap  = span->end       - span->start;
    const int right_overlap = (span + 1)->end - (span + 1)->start;

    return matches_right && 
      left_overlap  >= min_overlap &&
      right_overlap >= min_overlap;
  }

  return 0;
}


int load_model(GeneModel *gm, char *filename) {
  parse_gtf_file(gm, filename);
  index_regions(&gm->exons);
  add_transcripts(gm);
  add_introns(gm);
  return 0;
}

int extract_spans(Span *spans, bam1_t *read, int n) {
  CigarCursor c;
  init_cigar_cursor(&c, read);
  int i = 0;
  while (next_span(&c)) {
    spans[i].start = c.start;
    spans[i].end   = c.end;
    i++;
  }
  return i;
}

void init_cigar_cursor(CigarCursor *c, bam1_t *read) {
  c->read = read;
  c->pos = read->core.pos;
  c->start = c->end = 0;
  c->i = 0;
}

int next_fragment(bam1_t **reads, samfile_t *samfile, int n) {

  int num_reads = 0;

  if (samread(samfile, *reads) > 0) {

    num_reads++;
    
    if ( ! (reads[0]->core.flag & BAM_FPAIRED) )
      return 1;

    else if (samread(samfile, *(reads + 1)) > 0) {

      char *qname[] = { bam1_qname(reads[0]),  bam1_qname(reads[1]) };
      int      hi[] = { bam_aux2i(bam_aux_get(reads[0], "HI")),
                        bam_aux2i(bam_aux_get(reads[1], "HI")) };

      if (strcmp(qname[0], qname[1]) || hi[0] != hi[1]) {
        fprintf(stderr, "Error: paired flag is set, but next read is different. %s(%d) vs %s(%d) \n",
                qname[0], hi[0], qname[1], hi[1]);

        exit(1);
      }
      else {
        return 2;
      }
    }
    else {
      fprintf(stderr, "Error: paired flag is set, but there's no mate\n");
      exit(1);
    }
  }
  return 0;
}

int next_span(struct CigarCursor *c) {

  bam1_t *read = c->read;

  int n = read->core.n_cigar;
  uint32_t *cigar = bam1_cigar(read);

  int i;
  int found_match = 0;

  c->start = c->end = c->pos;

  for (i = c->i; i < n; i++) {
    int op    = bam_cigar_op(cigar[i]);
    int oplen = bam_cigar_oplen(cigar[i]);  
    int end_span = 0;

    switch (op) {
    
    case BAM_CMATCH:
      LOG_TRACE("  matching %d\n", oplen);
      c->end += oplen;
      c->pos += oplen;
      found_match = 1;
      break;
      
    case BAM_CINS:     
      LOG_TRACE("  inserting %d\n", oplen);
      break;

    case BAM_CDEL:
      LOG_TRACE("  deleting %d\n", oplen);
      c->end += oplen;
      c->pos += oplen;
      break;

    case BAM_CREF_SKIP:
      LOG_TRACE("  skipping %d\n", oplen);
      end_span = 1;
      c->pos += oplen;
      i++;
      break;

    case BAM_CSOFT_CLIP:
      LOG_TRACE("  soft-clipping %d\n", oplen);
      break; 
      
    case BAM_CHARD_CLIP:
    case BAM_CPAD: 
    case BAM_CEQUAL:
    case BAM_CDIFF:
    case BAM_CBACK:
      break;
    }

    if (end_span) 
      break;
  }
  c->i = i;

  return found_match;
}


int find_contiguous_exons(int *match_start,
                          Transcript *t, int start_exon,
                          Span *spans, int num_spans) {
  
  Span *s = spans;
  int flags;
  int i = start_exon;
  int n = t->exons_len;

  // Advance to the first exon that overlaps the first span
  while (i < n) {
    flags = cmp_exon(t->exons[i], t->exons[i]->chrom, s->start, s->end);
    if (flags & START_AFTER_EXON) 
      i++;
    else
      break;
  }

  // If no such exon exists, then we can't call a match for the
  // transcript.
  if ( i == n || (flags & END_BEFORE_EXON ) )
    return 0;

  *match_start = i;

  // If the span crosses the left edge of the exon and it isn't the
  // first exon in the transcript, then we can't call a match.
  if (i > 0 && (flags & CROSS_EXON_START ) )
    return 0;

  // Now we know the first span overlaps an exon in the transcript and
  // doesn't conflict with its left edge.

  // Now the sequence of spans should match the sequence of exons. For
  // each one except the last, the end should line up with the exons
  // end, and for each one except the first, the start should line up
  // with the exon's start.
  while ( s < spans + num_spans ) {
    if (s->end != t->exons[i]->end) 
      return 0;
    
    // Advance to the next span and next exon. If there's no next
    // exon, we have a mismatch
    s++;
    i++;
    
    if (i == n)
      return 0;

    // Make sure the start of all exons except the first line up.
    if (s->start != t->exons[i]->start)
      return 0;
  }

  flags = cmp_exon(t->exons[i], t->exons[i]->chrom, s->start, s->end);

  int is_last_exon = i + 1 == n;

  if ( !is_last_exon && (flags & CROSS_EXON_END) ) {
    return 0;
  }

  return i - (*match_start) + 1;
}

int matches_transcript(Transcript *transcript, 
                       Span *spans, 
                       int num_fwd_spans, 
                       int num_rev_spans) {

  int start_exon = 0;

  if (num_fwd_spans) {
    int match_start;
    int len = find_contiguous_exons(&match_start, 
                                    transcript, start_exon,
                                    spans, num_fwd_spans);
    if (!len)
      return 0;

    start_exon = match_start + len - 1;
  }
  
  if (num_rev_spans) {
    int match_start;
    int len = find_contiguous_exons(&match_start, 
                                    transcript, start_exon,
                                    spans + num_fwd_spans, num_rev_spans);
    
    if (!len)
      return 0;
  }
  
  return 1;
}
