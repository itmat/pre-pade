Pre-PADE
========

Pre-PADE is a set of tools that are used to create input data for the
PADE (Pallarel Analysis for Differential Effects) algorithm.

Specifically it focuses on preparing sets of quantitative
high-throughput genomics data, such as BED files or SAM/BAM files.

Installing
==========

Pre-PADE is built using CMake version 2.8.x.

To download and install it, please follow these steps:

1. Clone the repository

```bash
git clone git@github.com:itmat/pre-pade.git
```

2. Acquire samtools. It's linked as a git submodule, and you can pull it into the project by running:

```bash
git submodule init
git submodule update
```

3. Install cmake. You can download it here: http://www.cmake.org/cmake/resources/software.html.

4. Build the software

```bash
make
```

__Note:__ The above actually executes on a Makefile, which is a wrapper on the CMake files and build process, which itself creates it's own Makefiles for the actual build. A little wacky, we know, but it works for us.

5. (optional) Run the tests:

```bash
make test
```

6. (optional) Install the `quantify` executable to `/usr/local/bin`, run:

```bash
make install
```

Running
=======

`quantify` requires two input files: the gene model as a GTF file, and a
SAM file containing the aligned reads.

The GTF file should have the "gene_id", "transcript_id", and
"exon_number" tags set in the last field. We use these values to
determine the structure of the gene model.

The SAM file must be sorted by read id, not by chromosomal
coordinates.

To get help on running `quantify`, just run `./build/bin/quantify` without any
options.

Output File
===========

`quantify` produces a tab-delimited output file, with one row for each
feature that was covered by at least one read.

* __type__ - The type of the feature. Either "exon", "intron", "junction",
  or "transcript".

* __gene_id__ - From the gene_id attribute of the GTF file.

* __transcript_id__ - From the transcript_id attribute of the GTF file.

* __chrom__ - The chromosome.

* __exon_number__ - Taken from the exon_number attribute in the GTF
  file. The value of this field is interpreted differently depending
  on the value in the "type" field. When type is "exon", exon_number
  identifies the exon. When type is "intron" or "junction", it
  identifies the exon to the left of the intron or junction. When type
  is "transcript", exon_number is blank.

* __start, end__ - The start and end coordinates of the feature. Like
  exon_number, these fields are also interpreted differently depending
  on the value in the "type" field. When type is "exon", these fields
  simply contain the start and end coordinates of the exon. When type
  is "intron", they are the start and end coordinates of the
  intron. When type is "junction", start is the coordinate of the
  rightmost base in the left exon, and end is the leftmost base in the
  right exon. When type is "transcript", start and end contain
  comma-separated lists of the start and end points of the exons.

* __min__ - The number of uniquely aligned reads that hit this feature.

* __max__ - The total number of reads that hit this feature, including
  reads that have a unique alignment as well as reads that have
  multiple alignments.

## Scripts

To run scripts on the cluster:

    virtualenv ~/my_python2.6 --system-site-packages
    export PYTHONPATH=~/my_python2.6/lib/python2.6/site-packages/
    easy_install --install-dir ~/my_python2.6/lib/python2.6/site-packages argparse

    