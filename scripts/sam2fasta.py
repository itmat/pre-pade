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
    args.add_argument(
        '-r','--read_length',
        required=True,
        type=int,
        help="Read length in original fasta file"
    )
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
    
    fwd_all_qnames = set([])
    rev_all_qnames = set([])
    fwd_f = open(args.sam_file + "_fwd.fa", 'w' )
    rev_f = open(args.sam_file + "_rev.fa", 'w' )

    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if re.search(r'\.bam$',args.sam_file):
        read_mode = 'rb'
    src =  pysam.Samfile(args.sam_file,read_mode)
    
    while True:
      try:
          entry = src.next()
          logging.debug(entry.is_read1)
          
          if entry.is_read1 and not (entry.qname in fwd_all_qnames):
              logging.debug("I am read1 : " + str(entry))
              fwd_f.write(">" + str(entry.qname) + "\n")
              query = entry.query
              while len(query) < args.read_length:
                  query = query + "N"
              fwd_f.write(str(query)+ "\n")
              fwd_all_qnames.add(entry.qname)
            #fwd_f.write("+\n")
            #fwd_f.write(str(entry.qual) + "\n")
          elif not (entry.qname in rev_all_qnames):
            logging.debug("I am read2 :" + str(entry) + "\n")
            rev_f.write(">" + entry.qname + "\n")
            query = entry.query
            while len(query) < args.read_length:
                query = query + "N"
            rev_f.write(query + "\n")
            rev_all_qnames.add(entry.qname)
            #rev_f.write("+\n")
            #rev_f.write(entry.qual + "\n")
      except StopIteration, e:
          break
    fwd_f.close
    rev_f.close
      

if __name__ == '__main__':
    main()
