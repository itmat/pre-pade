#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import pysam
import logging
import os
import sys 
import re
from math import ceil
import numpy as np

"""Samples from a RUM SAM file where the sequences only come from
the canonical chromosomes. E.g. not the contigs or mitochondial
chromosomes. Likely only to work for human and mouse.
"""

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser(description="Downsamples SAM alignment files to by the given probability that a read should be output.")

    args.add_argument(
        '-i','--input',
        type=str,
        help="Input SAM/BAM file for sampling. Default: STDIN"
    )
    args.add_argument(
        '-o','--output',
        type=str,
        help="Output SAM/BAM filename. Default: STDOUT"
    )
    args.add_argument(
        '-b','--input_bam',
        action='store_true',
        help="Whether the input is a BAM file. Normally this is automatically determined by the input file name, but is handy for using STDIN as input."
    )
    args.add_argument(
        '-B','--output_bam',
        action='store_true',
        help="Whether the input is a BAM file. Normally this is automatically determined by the input file name, but is handy for using STDOUT as output."
    )
    args.add_argument(
        '-t','--total_reads',
        type=int,
        required=True,
        help="The number of read pairs that are in the SAM file"
    )
    args.add_argument(
        '-l','--output_limit',
        type=int,
        required=True,
        help="The maximum number of reads (pairs) that should be output. Defaults to value of --total_reads"
    )
    args.add_argument(
        '-n','--bin_num',
        type=int,
        default=1e7,
        help="The maximum number of reads (pairs) that should be output. Defaults to value of --total_reads"
    )
    args.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Print verbose information"
    )
    args.add_argument(
        '--debug', '-d',
        action='store_true',
        help="Print debugging information"
    )
    return args.parse_args()

def setup_logging(args):
    '''Sets up normal and verbose logging'''
    logging.basicConfig(level=logging.ERROR,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def get_output_indexes(num_fragments,output_ratio):
    """Returns an array of dtype=bool of size=size. Each index
    represents whether the current read should be output.

    >>> x = get_output_indexes(10,10)
    >>> len(x)
    10
    >>> x[0] == True
    True
    """
    result = np.random.random_sample(num_fragments)
    if logging.getLogger().level == logging.DEBUG:
        result[ result > output_ratio ] = 0
        result[result != 0] = 1
        nz = np.count_nonzero(result)
        logging.debug("Number of target indexes = %d" % nz)
        result = result == 1
    else:
        result = result <= output_ratio
    return result


def get_next_alignments(samfile,first_entry):
    """Returns a set of mapping entries for a read from a SAM file that map to one of the standard chromosomes.
    """
    entries = []
    last_entry = first_entry
    try:
        if first_entry == None:
            first_entry = samfile.next()
            last_entry = first_entry
        entries.append(first_entry)
        while True:
            tmp_entry = samfile.next()
            # go on to the next entry if not on a standard CHR
            if last_entry.qname == tmp_entry.qname:
                entries.append(tmp_entry)
                last_entry = tmp_entry
            else:
                last_entry = tmp_entry
                break
    except StopIteration, e:
        return (entries,last_entry)
    return (entries,last_entry)

def main():
    args = get_arguments()
    setup_logging(args)

    logging.debug(args)

    # handle STDIN/STDOUT streams
    if not args.input:
        args.input = '-'
    logging.debug("I=%s ; O=%s" % (args.input, args.output))

    # When input filename given, handle missing output filename
    if not args.output:
        if args.input == '-':
            args.output = '-'
        else:
            d,f = os.path.split(args.input)
            f = re.sub(r"\.([sb]am)$",r'.sampled.\1',f)
            args.output = os.path.join(d,f)

    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if args.input_bam or re.search(r'\.bam$',args.input):
        read_mode = 'rb'
    if args.output_bam or re.search(r'\.bam$',args.output):
        write_mode = 'wb'

    # Input/output/rejected SAM/BAM files
    src =  pysam.Samfile(args.input,read_mode)
    target = pysam.Samfile(args.output,write_mode,template=src)


    # some handy vars for bookkeeping and workflow
    total_reads = args.total_reads
    output_limit = args.output_limit
    output_ratio = output_limit /  total_reads 

    # set up the counts for a entry bin size == 1e7
    bin_size = args.bin_num
    if total_reads < bin_size:
        bin_size = total_reads
    num_bins = int(ceil(total_reads / bin_size))
    logging.debug("Number of bins: %d" % num_bins)

    # make an array of the number of reads each bin should output
    bin_sizes = [bin_size] * num_bins
    # adjust the number of reads output from the last bin to the size of the bin
    last_bin_read_num = total_reads % bin_size
    if last_bin_read_num > 0:
        logging.debug("Size of last bin: %d" % last_bin_read_num)
        bin_sizes[-1] = last_bin_read_num

    # grab the first read entry in the SAM file
    try:
        last_entry = src.next()
    except StopIteration, e:
        logging.error("No sequencese to output in this SAM/BAM file!")
        logging.error(e)
        exit()

    # Total number of reads output
    total_output_tally = 0
    last_tally = 0
    add_last_tally = False
    if output_ratio == 1.0:
        total_output_tally =  output_all_valid_entries(src,target,valid_chrs)
    else:
        for bin_index, bin_size in enumerate(bin_sizes):
            target_tally = int(ceil(bin_size * output_ratio))
            if add_last_tally:
                target_tally += last_tally
            last_tally = 0
            current_tally = 0
            output_indexes = get_output_indexes(bin_size,(target_tally / bin_size))
            # now output the reads
            last_entry_modifier = False
            logging.debug("Getting %d from %d fragments" % (target_tally, bin_size))

            for output_flag in output_indexes:
                entries, tmp_last_entry = get_next_alignments(src,last_entry)
                if output_flag or last_entry_modifier:
                    # output the reads to target
                    output_written = False
                    for e in entries:
                        target.write(e)
                        output_written = True
                    if output_written:
                        total_output_tally += 1
                        current_tally += 1
                        last_entry_modifier = False
                    else:
                        # designate that the next entry should be output
                        last_entry_modifier = True
                last_entry = tmp_last_entry
            if current_tally <= target_tally:
                add_last_tally = True
                last_tally = target_tally - current_tally
            logging.debug("Output %d" % current_tally)
    # output the number of 
    src.close()
    target.close()
    sys.stderr.write("Number of entries actually output %d\n" % total_output_tally)
if __name__ == '__main__':
    main()
