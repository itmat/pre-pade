#include <stdio.h>
#include <getopt.h>

#include "samutils.h"
#include "geneindex.h"
#include "sam.h"

// Exon level counts
// Transcript level counts
// Junction level counts

struct Args {
  char *gtf_filename;
  char *sam_filename;
  char *out_filename;
  char *details_filename;
  char *index_filename;
};

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
    { NULL,      0,                 NULL, 0   }
  };

  int bflag, ch, fd;
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

int main(int argc, char **argv) {
  
  struct Args args;
  parse_args(&args, argc, argv);

  FILE *details_file = NULL;
  if (args.details_filename) {
    details_file = fopen(args.details_filename, "w");
    if (!details_file) {
      perror(args.details_filename);
      exit(1);
    }
  }


  struct ExonDB db;
  parse_gtf_file(&db, args.gtf_filename);
  index_exons(&db);

  if (args.index_filename) {
    fprintf(stderr, "Dumping index to %s\n", args.index_filename);
    
    FILE *index_file = fopen(args.index_filename, "w");
    if (!index_file) {
      perror(args.index_filename);
      exit(1);
    }
    ExonIndexEntry *e;
    for (e = db.index; e < db.index + db.index_len; e++) {
      fprintf(index_file,  "%s:%d-%d\t%s:%d-%d\n",
              e->chrom, e->start, e->end, 
              e->exon->chrom, e->exon->start, e->exon->end);
    }
  }

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


  while (num_reads = next_fragment(reads, samfile, 2)) {

    int num_fwd_spans = extract_spans(read_spans, reads[0], MAX_SPANS_PER_READ);
    int num_rev_spans = 0;

    if (num_reads == 2) {
      num_rev_spans = extract_spans(read_spans + num_fwd_spans, reads[1], MAX_SPANS_PER_READ - num_fwd_spans);
    }

    char *qname = bam1_qname(reads[0]);
    LOG_TRACE("On read %s\n", qname);
    int pos = reads[0]->core.pos;
    int n_cigar = reads[0]->core.n_cigar;
    uint32_t *cigar = bam1_cigar(reads[0]);
    int hi = bam_aux2i(bam_aux_get(reads[0], "HI"));
    int ih = bam_aux2i(bam_aux_get(reads[0], "IH"));
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
        exon->max_count++;
        if (ih == 1)
          exon->min_count++;
      }
    }
    /*
    struct ExonCursor exon_curs;
    int span_num;
    for (span_num = 0; span_num < num_fwd_spans + num_rev_spans; span_num++) {
      Span *span = read_spans + span_num;

      struct Exon *exon;
      search_exons(&exon_curs, &db, ref, span->start, span->end, ALLOW_ALL);
      int flags = 0;
      while (exon = next_exon(&exon_curs, &flags)) {
        printf("%s\t", qname);
        printf("%d\t", ih);
        printf("%d\t", hi);
        printf("%s\t", ref);
        printf("%d\t", span->start);
        printf("%d\t", span->end);
        printf("%s\t", exon->gene_id);
        printf("%s\t", exon->transcript_id);
        printf("%d\t", exon->exon_number);
        printf("%d\t", exon->start);
        printf("%d\t", exon->end);
        printf("%d\t%d\t%d\t%d\n",
               (flags & CROSS_EXON_START) > 0,
               (flags & CROSS_EXON_END) > 0,
               (flags & START_IN_EXON) > 0,
               (flags & END_IN_EXON) > 0
               );
      }
      
      }
    */
  }
  int i;

  printf("%s\t", "gene_id");
  printf("%s\t", "transcript_id");
  printf("%s\t", "exon_number");
  printf("%s\t", "chrom");
  printf("%s\t", "start");
  printf("%s\t", "end");
  printf("%s\t", "min_count");
  printf("%s\n", "max_count");

  for (i = 0; i < db.exons.len; i++) {
    Exon *exon = db.exons.items + i;
    if (exon->min_count) {
      printf("%s\t", exon->gene_id);
      printf("%s\t", exon->transcript_id);
      printf("%d\t", exon->exon_number);
      printf("%s\t", exon->chrom);
      printf("%d\t", exon->start);
      printf("%d\t", exon->end);
      printf("%d\t", exon->min_count);
      printf("%d\n", exon->max_count);
    }
  }

  LOG_INFO("Cleaning up %s\n", "");

  bam_destroy1(reads[0]);
  bam_destroy1(reads[1]);
  //samclose(samfile);

  return 0;
}


