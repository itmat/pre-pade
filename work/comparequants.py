from __future__ import print_function, division

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
    df_b = pd.read_table(args.new_quants, index_col='name')

    joined = pd.merge(df_a, df_b, how='outer', suffixes=['_old', '_new'], left_index=True, right_index=True)
    output = args.output if args.output is not None else sys.stdout

    print(joined)

    for key in ['min', 'min_count', 'max', 'max_count']:
        joined[key][joined[key].isnull()] = 0
    
    for metric in ['min', 'max']:

        joined[metric + '_diff'] = joined[metric + '_count'] - joined[metric]
        joined[metric+ '_diff_squared'] = joined[metric + '_diff'] ** 2
        var = joined[metric + '_diff_squared'].sum()

        avg_var = var / len(joined)

        print(metric, ':', var, avg_var)

    min_is_diff = joined['min_diff'] != 0
    max_is_diff = joined['max_diff'] != 0
    
    joined[min_is_diff | max_is_diff].to_csv(args.diffs, sep='\t')

    print('Diffs with min: ' + str(sum(min_is_diff)) + ' (' + str(sum(min_is_diff) / len(joined)) + ')')

    joined.to_csv(args.output, sep='\t')


