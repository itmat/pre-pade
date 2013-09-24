#!/usr/bin/env python
from __future__ import print_function, division
import glob
import os
import re
import argparse
import logging
import csv
import locale

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "GLOB_PATTERN",
        type=str,
        help="A glob pattern like you would give to 'ls' to find a set of RUM mapping_stats.txt files. Example: '**/Sample_*/mapping_stats.txt'"
    )
    args.add_argument(
        '--output','-o',
        default="mapping_stat_report.csv",
        type=str,
        help="An output CSV file name. Default: 'mapping_stat_report.csv'"
    )
    #args.add_argument(
    #    '--species','-s',
    #    default="hg19",
    #    type=str,
    #    help="Please specify the species. Default: 'hg19'"
    #)
    args.add_argument(
        '--debug', '-d',
        action='store_true',
        help="Print debugging information"
    )
    args.add_argument(
        '--verbose', '-v',
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

def format_percentage(num):
    '''Formats a percentage for printing'''
    logging.debug("Formatting as %: " + str(num))
    return ("%0.2f" % (num * 100))  + "%"

def format_integer(num):
    '''Formats an integer for printing'''
    logging.debug("Formatting as Integer: " + str(num))
    return locale.format('%d', num, grouping=True)

def get_sample_name(bdir):
    logging.debug("Grabbing sample name from basedir: " + bdir)
    sn  = None
    if os.path.exists(os.path.join(bdir,'rum_job_report.txt')):
        fn = os.path.join(bdir,'rum_job_report.txt')
        logging.debug("Grabbing sample name from file: " + fn )
        data = open(fn,'rU').readlines()
        sn= re.findall(r'Job name : (\S+)', data[9])[0]
    elif os.path.exists(os.path.join(bdir,'rum.log_master')):
        fn = os.path.join(bdir,'rum.log_master')
        logging.debug("Grabbing sample name from file: " + fn )
        data = open(fn,'rU').readlines()
        sn =  re.findall(r'name: (\S+)', data[6])[0]
    else:
        raise Exception("Log file not found. Can't determine sample name")
    logging.debug("Grabbed sample_name: " + sn)
    return sn

def main():
    '''Run routine to pull mapping statistics from RUM mapping_stats.txt files.
    '''
    args = get_arguments()
    locale.setlocale(locale.LC_ALL,'en_US')
    setup_logging(args)
    report = csv.writer(open(args.output,'wb'))

    min_reads_num = float('inf')
    max_reads_num = 0
    min_aligned_num = float('inf')
    max_aligned_num = 0
    min_aligned_percentage = 1
    max_aligned_percentage = 0
    min_chrm_percentage = 1
    max_chrm_percentage = 0
    min_contig_percentage = 1
    max_contig_percentage = 0

    min_non_unique = float('inf')
    max_non_unique = 0
    min_non_unique_per = 1
    max_non_unique_per = 0
    min_unique = float('inf')
    max_unique = 0
    min_unique_per = 1
    max_unique_per = 0

    mapping_stats_files = []

    # grab the list of mapping stats files from the glob pattern
    for fn in glob.glob(args.GLOB_PATTERN):
        mapping_stats_files.append(fn)
    headers = ['Sample Name',
        "Total read pairs",
        "Total aligned fragments",
        "% aligned (total)",
        "Total unique",
        "% unique",
        "Total non unique",
        "% non unique",
        "% aligned (chrM)",
        "% aligned (contigs)" ]
    report.writerow(headers)
    logging.info("\t".join(headers))
    logging.info("-" * 180)

    for i,ms_fn in enumerate(mapping_stats_files):
        logging.debug("Grabbing mapping stats from: " + ms_fn)
        # get the sample name
        base_dir = os.path.dirname(ms_fn)
        sample_name = get_sample_name(base_dir)

        # now read in the mapping_stats and pull the information
        data = open(ms_fn,'rU').readlines()

        # grab the total reads and percentage mapped for this data set
        total_reads_num = re.findall(r'Number of read pairs\: (\S+)',data[0])[0]
        total_reads_num = int(total_reads_num.replace(',',''))

        total_aligned_num, total_aligned_percentage  = re.findall(r'At least one of forward or reverse mapped: (\S+) (\S+)',data[26])[0]
        total_aligned_num= int(total_aligned_num.replace(',',''))
        total_aligned_percentage = total_aligned_num / total_reads_num
        logging.debug(sample_name + ": total reads = " + str(total_reads_num))
        logging.debug(sample_name + ": total aligned num = " + str(total_aligned_num))

        logging.debug(sample_name + ": total aligned % = " + str(total_aligned_percentage))
        
        total_unique, total_unique_percentage  = re.findall(r'At least one of forward or reverse mapped: (\S+) (\S+)',data[11])[0]
        total_unique = int(total_unique.replace(',',''))
        total_unique_percentage = total_unique / total_reads_num

        logging.debug(sample_name + ": total unique = " + str(total_unique))
        logging.debug(sample_name + ": total unique % = " + str(total_unique_percentage))

        total_non_unique, total_non_unique_percentage  = re.findall(r'Total number consistent ambiguous: (\S+) (\S+)',data[18])[0]
        total_non_unique = int(total_non_unique.replace(',',''))
        total_non_unique_percentage = total_non_unique / total_reads_num
        
        logging.debug(sample_name + ": total non unique = " + str(total_non_unique))
        logging.debug(sample_name + ": total non unique % = " + str(total_non_unique_percentage))

        # grab the number of aligned reads for chrM and contigs
        chrm_reads_num = 0
        contig_reads_num = 0
        re_chr = re.compile(r'chr(\w+)\s+(\d+)')
        re_chrreg = re.compile(r'[^1-9XYM]{1,2}')
        # logging.debug(data[33:-1])
        for l in data[33:-1]:
            if "---------" in l:
                break
            m = re_chr.match(l)
            if not m:
                #no match to a Chr line
                continue
            elif m.group(1) == "M": 
                # This is the Mitochondrial
                chrm_reads_num += int(m.group(2))
            else:
                mm = re_chrreg.match(m.group(1))
                if mm:
                    logging.debug(mm.groups())
                    logging.debug("Matching chr num: " + l) 
                    contig_reads_num += int(m.group(2))
        logging.debug(sample_name + ": chrM aligned num = " + str(chrm_reads_num))
        logging.debug(sample_name + ": contigs aligned num = " + str(contig_reads_num))

        chrm_reads_percentage =  chrm_reads_num / total_reads_num
        contig_reads_percentage = contig_reads_num / total_reads_num

        # Adjust for Min & Max numbers
        if min_reads_num > total_reads_num:
            min_reads_num = total_reads_num
        if max_reads_num <  total_reads_num:
            max_reads_num =  total_reads_num

        if min_aligned_num > total_aligned_num:
            min_aligned_num =  total_aligned_num
        if max_aligned_num < total_aligned_num:
            max_aligned_num =  total_aligned_num

        if min_aligned_percentage > total_aligned_percentage:
            min_aligned_percentage = total_aligned_percentage
        if max_aligned_percentage <  total_aligned_percentage:
            max_aligned_percentage =  total_aligned_percentage

        if min_chrm_percentage > chrm_reads_percentage:
            min_chrm_percentage = chrm_reads_percentage
        if max_chrm_percentage <  chrm_reads_percentage:
            max_chrm_percentage =  chrm_reads_percentage

        if min_contig_percentage > contig_reads_percentage:
            min_contig_percentage = contig_reads_percentage
        if max_contig_percentage <  contig_reads_percentage:
            max_contig_percentage =  contig_reads_percentage

        if min_unique > total_unique:
            min_unique = total_unique
        if max_unique < total_unique:
            max_unique = total_unique

        if min_unique_per > total_unique_percentage:
            min_unique_per = total_unique_percentage
        if max_unique_per < total_unique_percentage:
            max_unique_per = total_unique_percentage

        if min_non_unique > total_non_unique:
            min_non_unique = total_non_unique
        if max_non_unique < total_non_unique:
            max_non_unique = total_non_unique

        if min_non_unique_per > total_non_unique_percentage:
            min_non_unique_per = total_non_unique_percentage
        if max_non_unique_per < total_non_unique_percentage:
            max_non_unique_per = total_non_unique_percentage


        # now print out the information
        report.writerow([sample_name,
            total_reads_num,
            total_aligned_num,
            total_aligned_percentage,
            total_unique,
            total_unique_percentage,
            total_non_unique,
            total_non_unique_percentage,
            chrm_reads_percentage,
            contig_reads_percentage])
        logging.info("\t".join([sample_name,
            format_integer(total_reads_num),
            format_integer(total_aligned_num),
            format_percentage(total_aligned_percentage),
            format_integer(total_unique),
            format_percentage(total_unique_percentage),
            format_integer(total_non_unique),
            format_percentage(total_non_unique_percentage),
            format_percentage(chrm_reads_percentage),
            format_percentage(contig_reads_percentage)]))

    report.writerow(['Minimums',
        min_reads_num,
        min_aligned_percentage])
    report.writerow(['Maximums',
        max_reads_num,
        max_aligned_percentage])
    logging.info("-" * 180 )
    logging.info("\t".join(['Minimums',
        format_integer(min_reads_num),
        format_integer(min_aligned_num),
        format_percentage(min_aligned_percentage),
        format_integer(min_unique),
        format_percentage(min_unique_per),
        format_integer(min_non_unique),
        format_percentage(min_non_unique_per),
        format_percentage(min_chrm_percentage),
        format_percentage(min_contig_percentage)]))
    logging.info("\t".join(['Maximums',
        format_integer(max_reads_num),
        format_integer(max_aligned_num),
        format_percentage(max_aligned_percentage),
        format_integer(max_unique),
        format_percentage(max_unique_per),                    
        format_integer(max_non_unique),
        format_percentage(max_non_unique_per),
        format_percentage(max_chrm_percentage),
        format_percentage(max_contig_percentage),
        ]))

if __name__ == '__main__':
    main()
