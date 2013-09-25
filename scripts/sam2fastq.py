#!/usr/bin/env python
from __future__ import print_function
import argparse
import pysam
import logging
import re

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "sam_file",
        type=str,
        help="The sam file that needs to be turned into fastq format"
    )
    #args.add_argument(
    #    '-m','--mapping_stats',
    #    required=True,
    #    type=str,
    #    help="mapping_stat_report.csv from pull_mapping_stats.py"
    #)
    args.add_argument(
        '-d','--debug',
        action='store_true',
        help="Print debugging information"
    )
    args.add_argument(
        '-v', '--verbose',
        action='store_true',
        help="Print verbose information. Basically this will print the mapping stats table to your screen."
    )
    return args.parse_args()

def setup_logging(args):
    '''Sets up normal and verbose logging'''
    logging.basicConfig(level=logging.ERROR,
                        format='%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def main():
    args = get_arguments()
    setup_logging(args)

    logging.debug(args)


    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if re.search(r'\.bam$',args.sam_file):
        read_mode = 'rb'
    src =  pysam.Samfile(args.sam_file,read_mode)
    while True:
      try:
          entry = src.next()
          if is_read1(entry):
            logging.debug("I am read1 :" + entry)
          else:
            logging.debug("I am read2 :" + entry)
      except StopIteration, e:
          break


if __name__ == '__main__':
    main()
