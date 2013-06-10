#include <stdio.h>
#include <getopt.h>

#include "samutils.h"
#include "geneindex.h"
#include "sam.h"

#define MIN_OVERLAP 8
// Exon level counts
// Transcript level counts
// Junction level counts

struct Args {
  char *gtf_filename;
  char *sam_filename;
  char *out_filename;
  char *details_filename;
  char *index_filename;
  int   exons;
  int   junctions;
};

void incr_quant(Quant *q, int unique) {
  q->max++;
  if (unique)
    q->min++;
}

void print_exon_quants(FILE *file, ExonDB *db) {

  int i;
  for (i = 0; i < db->exons.len; i++) {
    
    Exon *exon = db->exons.items + i;
    
    if (exon->exon_quant.min) {
      fprintf(file, "%s\t", "exon");
      fprintf(file, "%s\t", exon->gene_id);
      fprintf(file, "%s\t", exon->transcript_id);
      fprintf(file, "%d\t", exon->exon_number);
      fprintf(file, "%s\t", exon->chrom);
      fprintf(file, "%d\t", exon->start);
      fprintf(file, "%d\t", exon->end);
      fprintf(file, "%d\t", exon->exon_quant.min);
      fprintf(file, "%d\n", exon->exon_quant.max);
    }
  }

  for (i = 0; i < db->num_transcripts; i++) {
    
    Transcript *t = db->transcripts + i;
    
    Exon *left, *right;
    for (left = t->exons[0]; (right = next_exon_in_transcript(left)); left = right) {
      if (left->junction_quant.min) {
        fprintf(file, "%s\t",    "junction");
        fprintf(file, "%s\t",    left->gene_id);
        fprintf(file, "%s\t",    left->transcript_id);
        fprintf(file, "%d,%d\t", left->exon_number, right->exon_number);
        fprintf(file, "%s\t",    left->chrom);
        fprintf(file, "%d\t",    left->end - 1);
        fprintf(file, "%d\t",    right->start);
        fprintf(file, "%d\t",    left->junction_quant.min);
        fprintf(file, "%d\n",    left->junction_quant.max);
      }      
    }
  }
}



void parse_args(struct Args *args, int argc, char **argv) {
  args->sam_filename = NULL;
  args->gtf_filename = NULL;
  args->out_filename = NULL;
  args->details_filename = NULL;
  args->index_filename = NULL;

  static struct option longopts[] = {
    { "details", required_argument, NULL, 'd' },
    { "output",  required_argument, NULL, 'o' },
    { "index",   required_argument, NULL, 'x' },
    { "exons",         no_argument, NULL, 0   },
    { "junctions",     no_argument, NULL, 0   },
    { NULL,      0,                 NULL, 0   }
  };

  int ch;
  while ((ch = getopt_long(argc, argv, "d:o:x:", longopts, NULL)) != -1) {
    switch(ch) {
    case 'd':
      args->details_filename = optarg;
      break;

    case 'o':
      args->out_filename = optarg;
      break;

    case 'x':
      args->index_filename = optarg;
      break;

    default:
      printf("Bad usage\n");
    }

  }

  if (argc - optind != 2) {
    fprintf(stderr, "Usage: %s GTF_FILE SAM_FILE\n", argv[0]);
    exit(1);
  }
  
  args->gtf_filename = argv[optind];
  args->sam_filename = argv[optind + 1];
  
  fprintf(stderr, "GTF input file: %s\n", args->gtf_filename);
  fprintf(stderr, "SAM input file: %s\n", args->sam_filename);
  fprintf(stderr, "Output file: %s\n", args->out_filename ? args->out_filename : "(stdout)");
}

FILE *open_details_file(char *filename) {
  
  FILE *details_file = fopen(filename, "w");
  if (details_file)
    return details_file;
  
  perror(filename);
  exit(1);
}

void dump_index(char *filename, ExonDB *db) {
  
  fprintf(stderr, "Dumping index to %s\n", filename);
  
  FILE *index_file = fopen(filename, "w");
  if (!index_file) {
    perror(filename);
    exit(1);
  }
  ExonIndexEntry *e;
  for (e = db->index; e < db->index + db->index_len; e++) {
    fprintf(index_file,  "%s:%d-%d\t%s:%d-%d\n",
            e->chrom, e->start, e->end, 
            e->exon->chrom, e->exon->start, e->exon->end);
  }
  fclose(index_file);

}

int main(int argc, char **argv) {
  
  struct Args args;
  parse_args(&args, argc, argv);


  struct ExonDB db;
  load_model(&db, args.gtf_filename);

  FILE *details_file = args.details_filename ? 
    open_details_file(args.details_filename) : NULL;

  if (args.index_filename)
    dump_index(args.index_filename, &db);
 
  samfile_t *samfile = samopen(args.sam_filename, "r", NULL);
  
  struct Span *read_spans = calloc(MAX_SPANS_PER_READ, sizeof(struct Span));

  bam1_t *reads[] = { bam_init1(), bam_init1() };
  int num_reads;

  ExonMatches matches;
  init_exon_matches(&matches);  

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
  }

  while ((num_reads = next_fragment(reads, samfile, 2))) {

    int num_fwd_spans = extract_spans(read_spans, reads[0], MAX_SPANS_PER_READ);
    int num_rev_spans = 0;

    if (num_reads == 2) {
      num_rev_spans = extract_spans(read_spans + num_fwd_spans, reads[1], MAX_SPANS_PER_READ - num_fwd_spans);
    }

    char *qname = bam1_qname(reads[0]);
    LOG_TRACE("On read %s\n", qname);
    int hi = bam_aux2i(bam_aux_get(reads[0], "HI"));
    int num_alns = bam_aux2i(bam_aux_get(reads[0], "IH"));
    if (!num_alns) {
      num_alns = bam_aux2i(bam_aux_get(reads[0], "NH"));
    }

    int i;
    int tid = reads[0]->core.tid;
    char *ref = samfile->header->target_name[tid];
    ref = ref ? ref : "";

    find_candidates(&matches, &db, ref, read_spans, 
                    num_fwd_spans, num_rev_spans);

    for (i = 0; i < matches.len; i++) {

      int consistent = !matches.items[i].conflict;
      Exon *exon = matches.items[i].exon;

      if (details_file) {
        fprintf(details_file, "%s\t", exon->gene_id);
        fprintf(details_file, "%s\t", exon->transcript_id);
        fprintf(details_file, "%d\t", exon->exon_number);
        fprintf(details_file, "%s\t", exon->chrom);
        fprintf(details_file, "%d\t", exon->start);
        fprintf(details_file, "%d\t", exon->end);
        fprintf(details_file, "%s\t", qname);
        fprintf(details_file, "%d\t", hi);
        fprintf(details_file, "%d\n", consistent);
      }
      if (consistent) {
        incr_quant(&exon->exon_quant, num_alns == 1);
      }

      if (matches_junction(exon, read_spans, num_fwd_spans, num_rev_spans, MIN_OVERLAP)) {
        fprintf(stderr, "Found a junction!\n");
        incr_quant(&exon->junction_quant, num_alns == 1);
      }
    }
  }

  printf("%s\t", "type");
  printf("%s\t", "gene_id");
  printf("%s\t", "transcript_id");
  printf("%s\t", "exon_number");
  printf("%s\t", "chrom");
  printf("%s\t", "start");
  printf("%s\t", "end");
  printf("%s\t", "min_count");
  printf("%s\n", "max_count");

  print_exon_quants(stdout, &db);
  LOG_INFO("Cleaning up %s\n", "");

  bam_destroy1(reads[0]);
  bam_destroy1(reads[1]);
  //samclose(samfile);



  return 0;
}



