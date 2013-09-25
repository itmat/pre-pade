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
        help="A glob pattern like you would give to 'ls' to find a set of mapping_stats files. Example: \
'**/Sample_*/RUM.separated.uniq_mapping_stats'"
    )
    args.add_argument(
        '-m', '--mapping_stats',
        required=True,
        type=str,
        help="mapping_stat_report.csv from pull_mapping_stats.py"
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


if __name__ == '__main__':
    main()