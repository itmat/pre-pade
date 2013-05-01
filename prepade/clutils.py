class UsageException(Exception):
    pass

def setup_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO 

    logging.basicConfig(level=level,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=args.log,
                        filemode='w')

    console = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)

    
def get_output_fh(args):
    return args.output if args.output is not None else sys.stdout    

def require_sam_ordering_and_hi_tags(filename):
    ordering = input_file_ordering(filename)
    has_tags = has_hi_and_ih_tags(filename)

    # Check the alignments file to make sure it's ordered properly and
    # that it has IH and HI tags.
    if ordering is None and not has_tags:
        raise UsageException("""
The alignments must be either an indexed BAM file, or a BAM or SAM
file where all alignments for a single fragment are on consecutive
lines (e.g. sorted by read name). {alignments} does not seem to meet
either of thoseformats. Also, the HI and IH tags must be set for all
aligned reads,and that does not seem to be the case either. You will
probably need to preprocess the input file to add those tags.
""".format(args))

    elif ordering is None:
        raise UsageException("""

The alignments must be either an indexed BAM file, or a BAM or SAM
file where all alignments for a single fragment are on consecutive
lines (e.g. ordered by read name). {alignments} does not seem to meet
either of those formats. You should be able to sort it by read name by
doing:
        
  sort -n {alignments} OUTPUT_FILE
        """.format(**(args.__dict__)))

    elif not has_tags:
        raise UsageException("""
The HI and IH tags must be set for all aligned reads, and that does
not seem to be the case in {filename}. You will probably need to
preprocess the input file to add those tags.""")

    if ordering == AlignmentFileType.INDEXED:
        logging.info("The alignment file is an indexed BAM file")
    elif ordering == AlignmentFileType.ORDERED_BY_READ_NAME:
        logging.info("The alignment file appears to be ordered by read name")

    return ordering
