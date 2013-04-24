from __future__ import print_function

import pandas as pd
import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Filters a CSV file for columns that meet some criteria""")

    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    parser.add_argument('--diffs', type=argparse.FileType('w'))
    parser.add_argument('old_quants', type=file)
    parser.add_argument('new_quants', type=file)

    args = parser.parse_args()

    df_a = pd.read_table(args.old_quants, index_col='feature')
    df_b = pd.read_table(args.new_quants, index_col='feature')

    joined = pd.merge(df_a, df_b, how='outer', suffixes=['_old', '_new'], left_index=True, right_index=True)
    output = args.output if args.output is not None else sys.stdout
    joined['min_old'][joined['min_old'].isnull()] = 0
    
    joined['diff'] = joined['min_new'] - joined['min_old']
    joined['diff_squared'] = joined['diff'] ** 2

    var = joined['diff_squared'].sum()

    avg_var = var / len(joined)
    print(var, avg_var)

    joined.to_csv(args.output, sep='\t')

    joined[joined['diff'] != 0].to_csv(args.diffs, sep='\t')
