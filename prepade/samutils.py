import pysam
from itertools import groupby, ifilter, islice

class AlignmentFileType:
    INDEXED               = 0
    ORDERED_BY_READ_NAME  = 1

def qname_and_hi(aln):
    """Returns a tuple of the qname and 'HI' tag for the given AlignedRead."""
    return (aln.qname, aln.opt('HI'))


def sam_iter(samfile):
    while True:
        rec = samfile.next()
        yield rec

def has_hi_and_ih_tags(filename):
    
    samfile = pysam.Samfile(filename)
    
    try:
        mapped = ifilter(lambda x: not x.is_unmapped, sam_iter(samfile))    
    
        count = 0

        for aln in islice(mapped, 1000):
            aln.opt('HI')
            aln.opt('IH')

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

    mapped = list(islice(ifilter(lambda x: not x.is_unmapped, sam_iter(samfile)), 100))
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
                filter_mapped(
                    sam_iter(samfile))),
            10000)

        # Each of the alignments in each group should have an IH tag
        # which indicates how many alignments there are for that read,
        # and an HI tag which assigns a number to that alignment. If
        # the alignments are all grouped by qname, then if we look at
        # all the alignments in each group, the HI tags should number
        # 1 through IH.
        for (qname, grp) in groups:
            grp = list(grp)
            num_alns = grp[0].opt('IH')
            his = set([ a.opt('HI') for a in grp ])
            if len(his) != num_alns:
                return False
        return True

    finally:
        samfile.close()
