#include <stdio.h>
#include <getopt.h>
#include "sam.h"
#include "quant.h"

#define DEFAULT_MIN_OVERLAP 8

// Maximum number of spans that are allowed to appear in 
#define MAX_SPANS_PER_READ 1000

struct Args {

  // File names for GTF and SAM files. Should always be provided.
  char *gtf_filename;
  char *sam_filename;

  // Optional name of output file. May be NULL, in which case we use stdout.
  char *out_filename;

  // Optional name for file to write debugging output. If not given,
  // we won't write debugging output.
  char *details_filename;

  unsigned int  types_specified;
  unsigned int  do_types;
  int min_overlap;
};


/*
 * Some static structures to represent the types of quantification we
 * can do. We use this in parse_args to determine the valid values for
 * the -t/--type argument. 
 */

#define QUANT_TYPE_EXON       0
#define QUANT_TYPE_INTRON     1
#define QUANT_TYPE_JUNCTION   2
#define QUANT_TYPE_TRANSCRIPT 3

const char *QUANT_TYPE_NAMES[] = { "exon", "intron", "junction", "transcript" };
const int NUM_QUANT_TYPES = sizeof(QUANT_TYPE_NAMES) / sizeof(char*);

/*
 * Increment the max counter in the given Quant, and increment the min
 * counter if unique is true.
 */
void incr_quant(Quant *q, int unique) {
  q->max++;
  if (unique)
    q->min++;
}

void print_row(FILE *file, int type, Region **rs, int n, Quant *quant) {

  if (!quant->min)
    return;

  int i;

  Region *r = rs[0];
  fprintf(file, "%s\t%s\t%s\t%s\t", 
          QUANT_TYPE_NAMES[type], r->gene_id, r->transcript_id, r->chrom);

  // Print region numbers, starts, and ends. Exons and introns should
  // have 1 region, start, and end. Junctions should have two regions,
  // starts, and ends (the values for the exons on either side of the
  // junction). Right now transcripts have 0, but we might want to add
  // the coordinates for all the exons in the transcript.

  switch (type) {
  case QUANT_TYPE_EXON:
  case QUANT_TYPE_INTRON:
    fprintf(file, "%d\t%d\t%d", r->exon_number, r->start, r->end);
    break;

  case QUANT_TYPE_JUNCTION: 
    fprintf(file, "%d\t%d\t%d",
            rs[0]->exon_number, 
            rs[0]->end - 1,
            rs[1]->start);
    break;

  case QUANT_TYPE_TRANSCRIPT:
    fprintf(file, "\t");
    for (i = 0; i < n; i++) 
      fprintf(file, i == 0 ? "%d" : ",%d", rs[i]->start);
    fprintf(file, "\t");
    for (i = 0; i < n; i++) 
      fprintf(file, i == 0 ? "%d" : ",%d", rs[i]->end);
  }
  
  fprintf(file, "%d\t%d\n", quant->min, quant->max);

}

void print_quants(FILE *file, GeneModel *gm) {

  // Header row
  fprintf(file, "%s\t", "type");
  fprintf(file, "%s\t", "gene_id");
  fprintf(file, "%s\t", "transcript_id");
  fprintf(file, "%s\t", "exon_number");
  fprintf(file, "%s\t", "chrom");
  fprintf(file, "%s\t", "start");
  fprintf(file, "%s\t", "end");
  fprintf(file, "%s\t", "min_count");
  fprintf(file, "%s\n", "max_count");

  int i;

  // Exon counts
  LOG_INFO("Printing exons%s\n", "");
  for (i = 0; i < gm->exons.len; i++) {
     Region *exon = gm->exons.items + i;
     print_row(file, QUANT_TYPE_EXON, &exon, 1, &exon->exon_quant);
  }

  // Intron quants
  LOG_INFO("Printing introns%s\n", "");
  for (i = 0; i < gm->introns.len; i++) {
    Region *intron = gm->introns.items + i;
    print_row(file, QUANT_TYPE_INTRON, &intron, 1, &intron->exon_quant);
  }

  // Junction quants
  LOG_INFO("Printing junctions%s\n", "");
  for (i = 0; i < gm->num_transcripts; i++) {
    Transcript *t = gm->transcripts + i;
    int j;
    for (j = 0; j + 1 < t->exons_len; j++) {
      print_row(file, QUANT_TYPE_JUNCTION, t->exons + j, 2,
                &t->exons[j]->junction_quant);
    }
  }

  LOG_INFO("Printing transcripts%s\n", "");
  for (i = 0; i < gm->num_transcripts; i++) {
    Transcript *t = gm->transcripts + i;
    print_row(file, QUANT_TYPE_TRANSCRIPT, t->exons, t->exons_len, &t->quant);
  }


}


struct option longopts[] = {
  { "details",     required_argument, NULL, 'd' },
  { "output",      required_argument, NULL, 'o' },
  { "type",        no_argument,       NULL, 't' },
  { "min-overlap", required_argument, NULL, 'm' },
  { NULL,          0,                 NULL,  0  }
};


void usage(char *prog, int retval) {
  fprintf(stderr, "\nUsage: %s [OPTIONS] GTF_FILE SAM_FILE\n", prog);

  struct option *opt = longopts;
  const int size = 1000;
  char buf[size];
  
  printf("\nOptions:\n");

  for (opt = longopts; opt->name; opt++) {

    char *help = "";
    char *arg = "";

    switch (opt->val) {
    case 'd': 
      help = "Print each read/exon candidate pairing to FILE.";
      arg  = "FILE";
      break;

    case 'o':
      help = "Write output to FILE rather than stdout.";
      arg  = "FILE";
      break;

    case 'm':
      help = buf;
      snprintf(help, sizeof(buf), "Minimum overlap for junction hit, default %d", DEFAULT_MIN_OVERLAP);
      arg  = "BASES";
      break;

    case 't': {
      // TODO: This is really goofy, but I'm trying to avoid a possible buffer overflow.
      int cap = size;
      help = buf;
      help[0] = 0;
      cap -= snprintf(help, cap, "Quantify the given type of structure. Will quantify all if no type is\n    specified. Can be specified multiple times with different types. The\n    following are valid types:\n");
      int i;
      for (i = 0; i < NUM_QUANT_TYPES; i++) {
        char *tmp = strdup(help);
        cap -= snprintf(help, cap, "%s      %s\n", tmp, QUANT_TYPE_NAMES[i]);
        free(tmp);
      }
      arg = "TYPE";
      break;
    }
    default: break;
    }

    if (strcmp(arg, "")) {
      printf("  -%c, --%s %s\n    %s\n\n", opt->val, opt->name, arg, help);
    }
    else {
      printf("  -%c, --%s\n    %s\n\n", opt->val, opt->name, help);
    }
  }
  exit(retval);
}

void parse_args(struct Args *args, int argc, char **argv) {
  
  args->sam_filename = NULL;
  args->gtf_filename = NULL;
  args->out_filename = NULL;
  args->details_filename = NULL;
  args->min_overlap = DEFAULT_MIN_OVERLAP;

  int ch;
  args->types_specified = 0;
  int i;
  
  while ((ch = getopt_long(argc, argv, "d:o:t:m:", longopts, NULL)) != -1) {
    switch(ch) {
    case 'd':
      args->details_filename = optarg;
      break;

    case 'o':
      args->out_filename = optarg;
      break;

    case 'm':
      args->min_overlap = atoi(optarg);
      break;

    case 't':
      for (i = 0; i < sizeof(QUANT_TYPE_NAMES) / sizeof(char*); i++)
        if (!strcmp(QUANT_TYPE_NAMES[i], optarg)) 
          break;
      if (i == NUM_QUANT_TYPES) {
        fprintf(stderr, "Invalid value \"%s\" for --type or -t option. ", optarg);
        fprintf(stderr, "Valid options are:\n");
        for (i = 1; i < NUM_QUANT_TYPES; i++) {
          fprintf(stderr, "  %s\n", QUANT_TYPE_NAMES[i]);
        }
        exit(1);
      }
      args->types_specified |= (1 << (unsigned int)i);
      
    default:
      printf("Bad usage\n");
    }
  }

  if (args->types_specified) {
    fprintf(stderr, "Changed: %d\n", args->types_specified);
    args->do_types = args->types_specified;
  }
  else {
    args->do_types = ~0;
    fprintf(stderr, "Default\n");
  }


  if (argc - optind != 2) {
    usage(argv[0], 1);
  }
  
  args->gtf_filename = argv[optind];
  args->sam_filename = argv[optind + 1];
  
  fprintf(stderr, "GTF input file: %s\n", args->gtf_filename);
  fprintf(stderr, "SAM input file: %s\n", args->sam_filename);
  fprintf(stderr, "Output file: %s\n", args->out_filename ? args->out_filename : "(stdout)");
  fprintf(stderr, "Min overlap: %d\n", args->min_overlap);
}

FILE *open_details_file(char *filename) {
  
  FILE *details_file = fopen(filename, "w");
  if (details_file) {
    fprintf(details_file, "%s\t", "gene_id");
    fprintf(details_file, "%s\t", "transcript_id");
    fprintf(details_file, "%s\t", "exon_number");
    fprintf(details_file, "%s\t", "chrom");
    fprintf(details_file, "%s\t", "start");
    fprintf(details_file, "%s\t", "end");
    fprintf(details_file, "%s\t", "read_name");
    fprintf(details_file, "%s\t", "alignment_number");
    fprintf(details_file, "%s\n", "consistent");
    return details_file;
  }

  perror(filename);
  exit(1);
}


void print_match_details(FILE *file, bam1_t **reads, int num_reads, 
                         RegionMatches *matches) {
  int i;
  for (i = 0; i < matches->len; i++) {

    int hi = bam_aux2i(bam_aux_get(reads[0], "HI"));
    char *qname = bam1_qname(reads[0]);
    Region *exon = matches->items[i].region;

    fprintf(file, "%s\t", exon->gene_id);
    fprintf(file, "%s\t", exon->transcript_id);
    fprintf(file, "%d\t", exon->exon_number);
    fprintf(file, "%s\t", exon->chrom);
    fprintf(file, "%d\t", exon->start);
    fprintf(file, "%d\t", exon->end);
    fprintf(file, "%s\t", qname);
    fprintf(file, "%d\t", hi);
    fprintf(file, "%d\n", !matches->items[i].conflict);
  }
}

void print_types_quantified(unsigned int types) {
  int i;

  fprintf(stderr, "Types quantified:\n");
  for (i = 0; i < NUM_QUANT_TYPES; i++) {
    fprintf(stderr, "  %10s: %s\n", 
            QUANT_TYPE_NAMES[i],
            (types & (1 << i)) ? "yes" : "no");
  }
}

int cmp_transcript_ptr(Transcript **a, Transcript **b) {
  return strcmp((*a)->id, (*b)->id);
}

void get_distinct_transcripts(TranscriptMatches *transcript_matches,
                              RegionMatches *matches) {
  transcript_matches->len = 0;

  int i, j;
  for (i = 0; i < matches->len; i++) {
    RegionMatch *m = matches->items + i;
    Region *exon = m->region;
    Transcript *t = exon->transcript;
    add_transcript(transcript_matches, t);
  }  

  qsort(transcript_matches->transcripts, 
        transcript_matches->len, 
        sizeof(Transcript*), 
        ( int (*)(const void *, const void*) ) cmp_transcript_ptr);

  int n = transcript_matches->len;

  Transcript **ts = transcript_matches->transcripts;

  for (i = 0, j = 1; j < n; j++) {
    if ( strcmp(ts[i]->id, ts[j]->id) ) {
      ts[++i] = ts[j];
    }
    else {
      transcript_matches->len--;
    }
  }
} 

void accumulate_counts(GeneModel *gm, samfile_t *samfile, FILE *details_file, 
                       unsigned int types, int min_overlap) {
  struct Span read_spans[MAX_SPANS_PER_READ];

  bam1_t *reads[] = { bam_init1(), bam_init1() };
  int num_reads;

  const int do_exons       = types & (1 << QUANT_TYPE_EXON);
  const int do_introns     = types & (1 << QUANT_TYPE_INTRON);
  const int do_junctions   = types & (1 << QUANT_TYPE_JUNCTION);
  const int do_transcripts = types & (1 << QUANT_TYPE_TRANSCRIPT);
  print_types_quantified(types);

  RegionMatches matches;
  TranscriptMatches transcript_matches;
  init_exon_matches(&matches);
  init_transcript_matches(&transcript_matches);


  while ((num_reads = next_fragment(reads, samfile, 2))) {

    int num_fwd_spans = extract_spans(read_spans, reads[0], MAX_SPANS_PER_READ);
    int num_rev_spans = 0;

    if (num_reads == 2) {
      num_rev_spans = extract_spans(read_spans + num_fwd_spans, reads[1], MAX_SPANS_PER_READ - num_fwd_spans);
    }

    int i;
    int tid = reads[0]->core.tid;
    char *ref = samfile->header->target_name[tid];
    int num_alns = bam_aux2i(bam_aux_get(reads[0], "IH"));
    if (!num_alns) {
      num_alns = bam_aux2i(bam_aux_get(reads[0], "NH"));
    }

    ref = ref ? ref : "";

    LOG_TRACE("Finding exons%s\n", "");
    find_candidates(&matches, &gm->exons, ref, read_spans, num_fwd_spans, num_rev_spans);

    if (details_file) {
      print_match_details(details_file, reads, num_reads, &matches);
    }

    for (i = 0; i < matches.len; i++) {

      int consistent = !matches.items[i].conflict;
      Region *exon = matches.items[i].region;

      if (do_exons && consistent) {
        incr_quant(&exon->exon_quant, num_alns == 1);
      }

      if (do_junctions &&
          matches_junction(exon, read_spans, num_fwd_spans, num_rev_spans, min_overlap)) {
        incr_quant(&exon->junction_quant, num_alns == 1);
      }
    }


    if (do_transcripts) {

      get_distinct_transcripts(&transcript_matches, &matches);
      Transcript **ts = transcript_matches.transcripts;
      int n = transcript_matches.len;
      for (i = 0; i < n; i++) {
        Transcript *t = ts[i];
        if (matches_transcript(t, read_spans, num_fwd_spans, num_rev_spans)) {
          incr_quant(&t->quant, num_alns == 1);
        }
      }
    }

    if (do_introns) {
      LOG_TRACE("Finding introns%s\n", "");      
      find_candidates(&matches, &gm->introns, ref, read_spans, 
                      num_fwd_spans, num_rev_spans);

      for (i = 0; i < matches.len; i++) {
        Region *intron = matches.items[i].region;
        incr_quant(&intron->exon_quant, num_alns == 1);
      }

    }

  }
  bam_destroy1(reads[0]);
  bam_destroy1(reads[1]);

}

int main(int argc, char **argv) {
  
  struct Args args;
  parse_args(&args, argc, argv);

  FILE *out;
  if (args.out_filename) {
    out = fopen(args.out_filename, "w");
    if (!out) {
      perror(args.out_filename);
      exit(1);
    }
  }
  else 
    out = stdout;

  struct GeneModel gm;
  load_model(&gm, args.gtf_filename);

  FILE *details_file;
  if (args.details_filename) {
    fprintf(stderr, "Writing debugging output to %s\n", args.details_filename);
    details_file = open_details_file(args.details_filename);
  }
  else {
    details_file = NULL;
    fprintf(stderr, "Not producing debugging output\n");
  }


  samfile_t *samfile = samopen(args.sam_filename, "r", NULL);
  if (!samfile) {
    perror(args.sam_filename);
    exit(1);
  }

  accumulate_counts(&gm, samfile, details_file, args.do_types, args.min_overlap);
 
  print_quants(out, &gm);
  LOG_INFO("Cleaning up %s\n", "");

  samclose(samfile);

  return 0;
}



