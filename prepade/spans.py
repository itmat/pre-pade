from __future__ import print_function
from prepade.samutils import sam_iter, spans_for_aln
import pysam

samfile = pysam.Samfile("work/RUM.sam")

i = 0
for rec in sam_iter(samfile):
    spans = spans_for_aln(rec)
    spans_str = ", ".join([ "{0}-{1}".format(s.start, s.end) for s in spans ])
    print(i, rec.cigar, spans_str, sep="\t")
    i += 1
