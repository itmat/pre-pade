Pre-PADE
========

Pre-PADE is a set of tools that are used to create input data for the
PADE (Pallarel Analysis for Differential Effects) algorithm.

Specifically it focuses on preparing sets of quantitative
high-throughput genomics data, such as BED files or SAM/BAM files.

Installing
==========

You will need to compile the main pre-pade program, which is called
"quantify". It depends on the samtools library and header files. If
you don't already have a copy of samtools, please download it from
here: http://samtools.sourceforge.net/. Pre-pade has been tested with
samtools 0.1.19. Once you download samtools, please compile it. You
should be able to do that simply by typing "make" in the samtools
directory.

Once you've compiled samtools, you're ready to compile prepade. In the
root directory of the pre-pade distribution, type::

  make SAM_DIR=DIR

where DIR is the path to the samtools directory.

This should create a program called "quantify" in the bin
directory. You can thin run it by doing::

  bin/quantify

Running
=======

quantify requires two input files: the gene model as a GTF file, and a
SAM file containing the aligned reads.

The GTF file should have the "gene_id", "transcript_id", and
"exon_number" tags set in the last field. We use these values to
determine the structure of the gene model.

The SAM file must be sorted by read id, not by chromosomal
coordinates.

To get help on running quantify, just run "bin/quantify" without any
options.

Output File
===========

quantify produces a tab-delimited output file, with one row for each
feature that was covered by at least one read.

* type - The type of the feature. Either "exon", "intron", "junction",
  or "transcript".

* gene_id - The gene id, taken from the gene_id attribute of the GTF
  file.

* transcript_id - The transcript id, taken from the transcript_id
  attribute of the GTF file.

* chrom - The chromosome.

* exon_number - Taken from the exon_number attribute in the GTF
  file. The value of this field is interpreted differently depending
  on the value in the "type" field. When type is "exon", exon_number
  is a single number that identifies the exon. When type is "intron"
  or "junction", it is two numbers separated by a comma, identifying
  the exons on the left and right side of the intron. When type is
  "transcript", exon_number is blank.

* start, end - The start and end coordinates of the feature. Like
  exon_number, these fields are also interpreted differently depending
  on the value in the "type" field. When type is "exon", these fields
  simply contain the start and end coordinates of the exon. When type
  is "intron", they are the start and end coordinates of the
  intron. When type is "junction", start is the coordinate of the
  rightmost base in the left exon, and end is the leftmost base in the
  right exon. When type is "transcript", start and end contain
  comma-separated lists of the start and end points of the exons.

* min - The number of uniquely aligned reads that hit this feature.

* max - The total number of reads that hit this feature, including
  reads that have a unique alignment as well as reads that have
  multiple alignments.
