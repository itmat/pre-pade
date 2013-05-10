import pysam
from Bio.SeqFeature import FeatureLocation, SeqFeature
from itertools import groupby, ifilter, islice
import numpy as np

class AlignmentFileType:
    INDEXED               = 0
    ORDERED_BY_READ_NAME  = 1

class CigarOp:
    """Constants representing the CIGAR operations."""

    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    EQUAL = 7
    X = 8


def qname_and_hi(aln):
    """Returns a tuple of the qname and 'HI' tag for the given AlignedRead."""
    return (aln.qname, aln.opt('HI'))


def sam_iter(samfile, skip_unmapped=True):
    while True:
        rec = samfile.next()
        if skip_unmapped and rec.is_unmapped:
            continue
        else:
            yield rec

def has_hi_and_ih_tags(filename):
    
    samfile = pysam.Samfile(filename)
    
    try:
        mapped = sam_iter(samfile)
    
        count = 0

        for aln in islice(mapped, 1000):
            aln.opt('HI')
            try:
                aln.opt('IH')
            except KeyError:
                try:
                    aln.opt('NH')
                except KeyError:
                    raise Exception("Input " + filename + " does not seem to have IH or NH tags")
                
                    

        return True
    finally:
        samfile.close()

def input_file_ordering(filename):
    
    """Return AlignmentFileType value indicating ordering, or raise if unordered.

    :param filename:

      The input file to check.

    :return:

      * AlignmentFileType.INDEXED if the file is an indexed BAM file.
    
      * AlignmentFileType.ORDERED_BY_READ_NAME if the file is a SAM or
        BAM file that appears to be ordered by read name.

    :raises UsageException:
      when the file does not appear to be BAM or SAM file that is
      either indexed or ordered by read name.

    """

    if input_file_is_indexed(filename):
        return AlignmentFileType.INDEXED

    elif input_file_is_ordered_by_read_name(filename):
        return AlignmentFileType.ORDERED_BY_READ_NAME

    else:
        return None


def input_file_is_indexed(filename):
    """Returns true of the filename is an indexed BAM file."""

    count = 0    
    samfile = pysam.Samfile(filename)         

    mapped = list(islice(sam_iter(samfile), 100))
    rec = mapped[0]
    
    try:
        samfile.fetch(samfile.getrname(rec.tid), rec.pos)
        return True

    except ValueError:
        return False

    finally:
        samfile.close()

def filter_mapped(alns):
    return ifilter(lambda x: not x.is_unmapped, alns)

def grouped_by_qname(alns):
    return groupby(alns, lambda a: a.qname)

def input_file_is_ordered_by_read_name(filename):
    """Checks the first 10,000 alignments and returns true if they are
    ordered by read name.

    """
    read_names = []
    samfile = pysam.Samfile(filename)
    try:
        
        # Open the sam file, turn it into an iterator, filter to
        # include only mapped alns, group them by read name, and
        # select the first 10000 of those groups.
        groups = islice(
            grouped_by_qname(
                    sam_iter(samfile)),
            10000)

        # Each of the alignments in each group should have an IH tag
        # which indicates how many alignments there are for that read,
        # and an HI tag which assigns a number to that alignment. If
        # the alignments are all grouped by qname, then if we look at
        # all the alignments in each group, the HI tags should number
        # 1 through IH.
        for (qname, grp) in groups:
            grp = list(grp)
            try:
                num_alns = grp[0].opt('IH')
            except KeyError:
                num_alns = grp[0].opt('NH')
            his = set([ a.opt('HI') for a in grp ])
            if len(his) != num_alns:
                return False
        return True

    finally:
        samfile.close()

def alns_by_qname_and_hi(alns):
    """Group alignments by qname and HI tag

    :param alns:
      An iterator over alignments. All alignments with the same qname
      must be contiguous. Within a group of alignments with the same
      qname, the ordering of the reads by HI tag does not matter.

    :return:
      An iterator over lists of alignments, grouped by qname and the
      value of the HI tags.

    """
    for (qname, hi), grp in groupby(alns, key=lambda x: (x.qname, x.opt('HI'))):
        yield grp

def spans_for_aln(aln):
    return cigar_to_spans(aln.cigar, aln.pos)

def cigar_to_spans(cigar, start):
    """Converts a CIGAR data structure and position to list of FeatureLocations.

    :param cigar:
      List of tuples like what is returned by the 'cigar' property
      from a Pysam alignment.

    :param start:
       Start location, using 0-based indexing.

    :return:

      List of FeatureLocation objects representing the spans of the
      reference that this segment matches.
      
    """
    spans = []

    if cigar is None:
        return []

    cigar = remove_ds(cigar)

    for (op, bases) in cigar:

        if op == CigarOp.M:
            end = start + bases
            spans.append(FeatureLocation(start, end))
            start = end

        elif op == CigarOp.D or op == CigarOp.N:
            start = start + bases

    res = []

    for span in spans:
        if len(res) > 0 and res[-1].end + 1>= span.start:
            start = res[-1].start
            end   = span.end
            res[-1] = FeatureLocation(start, end)
        else:
            res.append(span)

    return res

def aln_to_span_ndarray(aln):
    locs = cigar_to_spans(aln.cigar, aln.pos)
    res = np.zeros((len(locs), 2), int)
    for i, loc in enumerate(locs):
        res[i, 0] = loc.start
        res[i, 1] = loc.end
    return res

def remove_ds(cigar):
    """Removes D operations from CigarOp string, replacing with Ms.

    Replaces all Ds with Ms and then merges adjacent Ms together.

    >>> remove_ds([ (0, 21), (2, 1), (0, 54) ])
    [(0, 76)]

    >>> remove_ds([(4, 8), (0, 4), (2, 1), (0, 63)])
    [(4, 8), (0, 68)]

    >>> remove_ds([(0, 15), (2, 1), (0, 15), (2, 2), (0, 29), (2, 2), (0, 16)])
    [(0, 80)]

    >>> remove_ds([(0, 21), (2, 1), (0, 41), (3, 177), (0, 13)])
    [(0, 63), (3, 177), (0, 13)]

    >>> remove_ds([(0, 41), (3, 354), (0, 20), (2, 1), (0, 14)])
    [(0, 41), (3, 354), (0, 35)]

    >>> remove_ds([(4, 4), (0, 26), (2, 1), (0, 45)])
    [(4, 4), (0, 72)]

    """

    d_to_m = lambda (op, bases): (CigarOp.M, bases) if op == CigarOp.D else (op, bases)
    converted = map(d_to_m, cigar)
    res = []

    res = [converted[0]]

    for (op, bases) in converted[1:]:
        (last_op, last_bases) = res[-1]
        
        if op == last_op:
            res[-1] = (op, bases + last_bases)
        else:
            res.append((op, bases))

    return res


