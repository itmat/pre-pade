#!/usr/bin/env python
from __future__ import print_function
import argparse
import logging
import glob
import os
import re

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "GLOB_PATTERN",
        type=str,
        help="A glob pattern like you would give to 'ls' to find a set of RUM.sam files. Example: \
'**/Sample_*/RUM.sam'"
    )
    args.add_argument(
        '--output','-o',
        default="run_job.sh",
        type=str,
        help="An output file_name. Default: 'run_job.sh'"
    )
    args.add_argument(
        '--sungrid', '-g',
        action='store_true',
        help="Sun grid engine cluster?"
    )
    args.add_argument(
        '--species', '-s',
        type=str,
        choices= ["hg19","mm9","drosophila"],
        help="Species used in experiment"
    )
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
    '''Run routine to set up pipeline for post mapping.
    '''
    args = get_arguments()
    #locale.setlocale(locale.LC_ALL,'en_US')
    setup_logging(args)

    rum_sam_files = []
    for fn in glob.glob(args.GLOB_PATTERN):
        rum_sam_files.append(fn)
    for i,ms_fn in enumerate(rum_sam_files):
        logging.debug("Working on: " + ms_fn)
        base_dir = os.path.dirname(ms_fn)
        sample_name = get_sample_name(base_dir)
        job_name = sample_name + "_pipeline"
        LSF_header=("#!/bin/sh\n\n"
                    "#BSUB -J " + job_name  + "\n"
                    "#BSUB -oo outfile.%J\n"
                    "#BSUB -eo errorfile.%J\n"
                    "#BSUB -q max_mem30\n")
        if args.sungrid:
            LSF_header=("#!/bin/sh\n\n"
                    "#$ -N " + job_name  + "\n"
                    "#$ -V\n"
                    "#$ -cwd\n"
                    "#$ -j y\n"
                    "#$ -l h_vmem=15G\n")
        logging.debug(LSF_header)
        out_file = base_dir + "/" + job_name + "_jobfile"
        logging.debug("Outfile: " + out_file)
        settings=("module load python-2.7.5\n"
                  "export PYTHONPATH=~/my_python2.7/lib/python2.7/site-packages/\n")
        if args.sungrid:
            settings=""
        f = open(out_file,'w')
        f.write(LSF_header)
        f.write(settings)
        species = ""
        if args.species:
            species = "-s "
        f.write("sam_filter.py -i RUM.sam -r RUM.rejected.sam -v " + species + "\n")
        f.write("sam_filter_uniq.py RUM.filtered.sam RUM.separated\n")
        f.write("sam2mappingstats.pl RUM.separated.nuniq.sam > RUM.separated.nuniq_mapping_stats\n")
        f.write("sam2mappingstats.pl RUM.separated.uniq.sam > RUM.separated.uniq_mapping_stats\n")
        f.flush()
        os.fsync(f.fileno())
        f.close
        old_dir = os.getcwd()
        os.chdir(base_dir)
        logging.debug("Current work dir: " + os.getcwd())
        if args.sungrid:
            commando = "qsub " + job_name + "_jobfile"
        else:
            commando = "bsub < " + job_name + "_jobfile"
        logging.debug(commando)
        #logging.debug(os.system(commando))
        os.chdir(old_dir)
        
if __name__ == '__main__':
    main()
