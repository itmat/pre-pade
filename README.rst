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
