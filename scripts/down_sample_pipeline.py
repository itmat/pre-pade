#!/usr/bin/env python
import argparse
import logging
import csv
import glob
import os
import re

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
        '-m','--mapping_stats_unique',
        required=True,
        type=argparse.FileType('r'),
        help="mapping_stat_report.csv from pull_mapping_stats.py"
    )
    args.add_argument(
        '-n','--mapping_stats_non_unique',
        required=True,
        type=argparse.FileType('r'),
        help="mapping_stat_report.csv from pull_mapping_stats.py"
    )
    args.add_argument(
        '--sungrid', '-g',
        action='store_true',
        help="Sun grid engine cluster?"
    )
    args.add_argument(
        '-a','--gtf_file',
        required=True,
        type=argparse.FileType('r'),
        help="gtf file with gene annotations."
    )
    args.add_argument(
        '-f','--fai_file',
        required=True,
        type=argparse.FileType('r'),
        help="fai file with chromosome annotation (obtained from `samtools faidx *.fa`)."
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
    
    current_dir = os.getcwd()

    
    field = "Total unique"
    my_dict_unique = csv.DictReader(args.mapping_stats_unique)
    minimum_unique = 0
    for row in my_dict_unique:
        #logging.debug(row)
        if row['Sample Name'] == 'Minimums':
            minimum_unique = row[field]
    logging.debug("The minimum (unique) is: " + minimum_unique)
    
    field = "Total non unique"
    my_dict_non_unique = csv.DictReader(args.mapping_stats_non_unique)
    minimum_non_unique = 0
    for row in my_dict_non_unique:
        #logging.debug(row)
        if row['Sample Name'] == 'Minimums':
            minimum_non_unique = row[field]
    logging.debug("The minimum (non_unique) is: " + minimum_non_unique)
    

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
                    "#BSUB -oo " + base_dir + "/outfile.%J\n"
                    "#BSUB -eo " + base_dir + "/errorfile.%J\n"
                    "#BSUB -q max_mem30\n")
        if args.sungrid:
            LSF_header=("#!/bin/sh\n\n"
                    "#$ -N " + job_name  + "\n"
                    "#$ -V\n"
                    "#$ -cwd\n"
                    "#$ -j y\n"
                    "#$ -l h_vmem=15G\n")
        settings=("module load python-2.7.5\n"
                  "export PYTHONPATH=~/my_python2.7/lib/python2.7/site-packages/\n")
        if args.sungrid:
            settings = ""
        read_number = 0
        args.mapping_stats_unique.seek(0)
        field = "Total unique"
        for row in my_dict_unique:
            if row['Sample Name'] == sample_name:
                read_number = row['Total read pairs']
                break
        out_file = base_dir + "/down_sample_pipeline.sh"
        f = open(out_file,'w')
        f.write(LSF_header)
        f.write(settings)
        uniq_name = ms_fn
        
        #f.write("sam_sampler.py -i " + ms_fn + " -t " + str(read_number) + 
        #        " -l " + str(minimum_unique) + " -o " + ms_fn + ".sampled.sam" + "\n")
        
        read_number = 0
        args.mapping_stats_non_unique.seek(0)
        field = "Total non unique"
        for row in my_dict_non_unique:
            if row['Sample Name'] == sample_name:
                read_number = row['Total read pairs']
                break
        nonuniq_name = re.sub('uniq', 'nuniq',ms_fn)
        logging.debug(nonuniq_name)
        f.write("sam_sampler.py -i " + nonuniq_name + " -t " + str(read_number) +
                " -l " + str(minimum_non_unique) + " -o " + nonuniq_name + ".sampled.sam" + "\n")
        f.write("grep -v ^@ " + nonuniq_name + ".sampled.sam > " + nonuniq_name + ".sampled.sam_no_header\n")
        f.write("cat " + ms_fn + ".sampled.sam " + nonuniq_name + ".sampled.sam_no_header > " + ms_fn + "_combined.sam\n")
        
        #f.write("sam2fasta.py -r 51 " + ms_fn + "_combined.sam\n")
        #f.write("rum_runner align --name " + sample_name + " -i /home/apps/RUM/indexes_2.x/drosophila/ --chunks 10 --platform LSF -o " + base_dir + "/rum_merged " + ms_fn + "_combined.sam_fwd.fa " + ms_fn + "_combined.sam_rev.fa\n")
        
        # COVERAGE FILES
        f.write("/home/hayer/tools/sam2cov/sam2cov -r 1 -p " + ms_fn + "_ " + str(args.fai_file.name) + " "  + ms_fn + "_combined.sam\n" )
        # QUANTIFY
        f.write("quantify -t exon -t junction -t transcript " + str(args.gtf_file.name) + " " +  ms_fn + "_combined.sam -o " + ms_fn + "_quantified\n")
        #f.write("htseq-count --stranded=no " + ms_fn + "_combined.sam " + str(args.gtf_file.name) + " > " + ms_fn + "_quantified_htseq\n")
        # MAPPING STATS
        f.write("sam2mappingstats.pl " + ms_fn + "_combined.sam > " + ms_fn + "_mapping_stats.txt\n")
        
        f.flush()
        os.fsync(f.fileno())
        f.close
        #os.chdir(base_dir)
        #logging.debug("Current work dir: " + os.getcwd())
        commando = "bsub < " + out_file
        if args.sungrid:
            commando = "qsub " + out_file
        logging.debug(commando)
        logging.debug(os.system(commando))
        #os.chdir(current_dir)

if __name__ == '__main__':
    main()
