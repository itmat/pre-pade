Pre-PADE
========

Pre-PADE is a set of tools that are used to create input data for the
PADE (Pallarel Analysis for Differential Effects) algorithm.

Specifically it focuses on preparing sets of quantitative
high-throughput genomics data, such as BED files or SAM/BAM files.

Installing
==========

Pre-PADE is built using CMake. To download and install it, please
follow these steps:

1. Clone the repository

        git clone git@github.com:itmat/pre-pade.git

2. Acquire samtools

        git submodule init
        git submodule update

3. Install cmake
4. mkdir build; cd build
5. cmake ..
6. make

This should create the main program, build/bin/quantify.

Pre-pade has been tested with samtools 0.1.19. Once you download
samtools, please compile it.

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

* gene_id - From the gene_id attribute of the GTF file.

* transcript_id - From the transcript_id attribute of the GTF file.

* chrom - The chromosome.

* exon_number - Taken from the exon_number attribute in the GTF
  file. The value of this field is interpreted differently depending
  on the value in the "type" field. When type is "exon", exon_number
  identifies the exon. When type is "intron" or "junction", it
  identifies the exon to the left of the intron or junction. When type
  is "transcript", exon_number is blank.

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

