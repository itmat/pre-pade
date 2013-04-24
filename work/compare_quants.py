from __future__ import print_function

import pandas as pd
import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Filters a CSV file for columns that meet some criteria""")

    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    parser.add_argument('old_quants', type=file)
    parser.add_argument('new_quants', type=file)

    args = parser.parse_args()

    df_a = pd.read_table(args.old_quants, index_col=args.join)
    df_b = pd.read_table(args.new_quants, index_col=args.join)

    joined = pd.merge(df_a, df_b, how='outer', suffixes=['_left', '_right'], left_index=True, right_index=True)
    output = args.output if args.output is not None else sys.stdout
    joined.to_csv(args.output, sep='\t')
