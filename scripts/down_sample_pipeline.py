#!/usr/bin/env python
import argparse
import logging
import csv

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "GLOB_PATTERN",
        type=str,
        help="A glob pattern like you would give to 'ls' to find a set of sam files. Example: \
'**/Sample_*/RUM.separated.uniq.sam'"
    )
    args.add_argument(
        '-m','--mapping_stats',
        required=True,
        type=argparse.FileType('r'),
        help="mapping_stat_report.csv from pull_mapping_stats.py"
    )
    args.add_argument(
        '-n','--non_unique',
        action='store_true',
        help="Run in non_unique mode"
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
    args = get_arguments()
    setup_logging(args)
    logging.debug(args)

    # Take care of run mode
    if args.non_unique:
        field = "Total non unique"
    else:
        field = "Total unique"
    my_dict = csv.DictReader(args.mapping_stats)
    minimum = 0
    for row in my_dict:
        #logging.debug(row)
        if row['Sample Name'] == 'Minimums':
            minimum = row[field]
    logging.debug("The minimum is: " + minimum)

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
                    "#BSUB -eo errorfile.%J\n")
        settings=("module load python-2.7.5\n"
                  "export PYTHONPATH=~/my_python2.7/lib/python2.7/site-packages/\n")
        read_number = 0
        for row in my_dict:
        #logging.debug(row)
        if row['Sample Name'] == sample_name:
            read_number = row[field]

        f = open(out_file,'w')
        f.write(LSF_header)
        f.write(settings)
        f.write("sam_filter.py -i " + mn_fn + " -t " + read_number + " -l " + 
                minimum + " > " + mn_fn.sampled.sam + "\n")
        f.write(

if __name__ == '__main__':
    main()
