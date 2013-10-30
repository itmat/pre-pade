#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import logging
import os
import sys
import re
from math import ceil
import numpy as np
import pysam

"""Filters a SAM/BAM file, sorted in read ID order, to sequences that only come from the canonical chromosomes. E.g. not the contigs or mitochondial chromosomes. Likely only to work for human and mouse.
"""

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser(description="Filters a SAM/BAM alignment file (sorted by read IDs) to canonical chromosomal alignments. Optionally outputs rejected reads to another file.")
    args.add_argument(
        "input",
        type=str,
        help="Input SAM/BAM file for filtering."
    )
    args.add_argument(
        'output_prefix',
        type=str,
        help="Output SAM/BAM filename prefix. Two files will be produced. <prefix>.uniq.sam and <prefix>.nuniq.sam"
    )
    args.add_argument(
        '-B', '--output-bam',
        action='store_true',
        help="Output files are BAM format, instead of SAM. E.b. <prefix>.uniq.bam and <prefix>.nuniq.bam"
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

def get_next_alignments(samfile,last_entry):
    """Returns a set of mapping entries for a read from a SAM file that map to one of the standard chromosomes.
    """
    entries = []
    try:
        if last_entry == None:
            last_entry = samfile.next()
        entries.append(last_entry)
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
        logging.info("Reached end of SAM file " + samfile.filename)
        entries.append(last_entry)
        last_entry = None
    return (entries,last_entry)

def determine_entries_status2(entry):
    status = 'unmapped'
    tags = entry.tags
    #logging.debug("NAME: " + valid_entry.qname)
    for tag in old_tags:
        if (tag[0] == 'HI' or tag[0] == "NH"):
            if tag[1] == 1:
                status = 'uniq'
            else:
                status = 'nuniq'
    return status

def determine_entries_status(entries):
    count = 0
    status = 'unmapped'
    if not entries[0].is_unmapped:
        status = 'uniq'
    # check to see if this has more that one entry
    if len(entries) > 2 and entries[0].is_paired:
        status = 'nuniq'
    elif len(entries) > 1 and not entries[0].is_paired:
        status = 'nuniq'
    e = entries[0]
    logging.debug("E {e}: length={l} paired={p} secondary={s} mapped={m} status={t}".format(
        e=e.qname,
        l=len(entries),
        p=e.is_paired,
        s=e.is_secondary,
        m=e.is_unmapped,
        t=status))
    return status


def write_entries(entries,samfile):
    for entry in entries:
        samfile.write(entry)

def main():
    args = get_arguments()
    setup_logging(args)
    logging.debug(args)

    # define output filenames
    uniq_fn = nuniq_fn = ''
    if args.output_bam:
        uniq_fn = args.output_prefix + ".uniq.bam"
        nuniq_fn = args.output_prefix + ".nuniq.bam"
    else:
        uniq_fn = args.output_prefix + ".uniq.sam"
        nuniq_fn = args.output_prefix + ".nuniq.sam"

    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if re.search(r'\.bam$',args.input):
        read_mode = 'rb'
    if args.output_bam or re.search(r'\.bam$',uniq_fn):
        write_mode = 'wb'

    # Input/output/rejected SAM/BAM files
    src =  pysam.Samfile(args.input,read_mode)
    uniq = pysam.Samfile(uniq_fn,write_mode,template=src)
    nuniq = pysam.Samfile(nuniq_fn,write_mode,template=src)

    # Total number of reads output
    total_uniq = total_unmapped = total_nu =  0

    # Classic Python iteration using PySAM's fetch() function is
    # **much** slower than using the loop below.
    last_entry = None
    try: 
        while True:
            #entries, last_entry = get_next_alignments(src,last_entry)
            #if not last_entry:
            #    break
            entry = src.next()
            uniq_written = False
            nu_written = False
            status = determine_entries_status2(entry)
            writer = None
            if status == 'uniq':
                write_entries(entry, uniq)
                total_uniq += 1
            elif status == 'nuniq':
                write_entries(entry, nuniq)
                total_nu += 1
            elif 'unmapped':
                total_unmapped += 1

    except StopIteration,e:
        pass

    src.close()
    uniq.close()
    nuniq.close()

    sys.stderr.write("Number of unique entries: %d\n" % total_uniq)
    sys.stderr.write("Number of non-unique entries: %d\n" % total_nu)
    sys.stderr.write("Number of unmapped entries: %d\n" % total_unmapped)

if __name__ == '__main__':
    main()
