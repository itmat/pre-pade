Pre-PADE
========

Pre-PADE is a set of tools that are used to create input data for the
PADE (Pallarel Analysis for Differential Effects) algorithm.

Specifically it focuses on preparing sets of quantitative
high-throughput genomics data, such as BED files or SAM/BAM files.

Installing
==========

If you want to install pre-pade tools into a development environment
(a local directory, not a system install) do:

  python setup.py develop

Running
=======

Currently the prepade tools are available as runnable python modules,
not as executable scripts. We might make executable wrappers for them
at some point.

Please run:

  python -m prepade.quantifyexons
