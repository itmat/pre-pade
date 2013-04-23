from __future__ import print_function

import pandas as pd
import argparse
import sys

def int_transform(fieldnames):
    def transform(x):
        for fieldname in fieldnames:
            x[fieldname] = int(x[fieldname])
        return x
    return transform    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Filters a CSV file for columns that meet some criteria""")

    parser.add_argument('--join', '-j')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    parser.add_argument('file_a', type=file)
    parser.add_argument('file_b', type=file)

    kwargs = { 
        'delimiter' : '\t'
    }

    args = parser.parse_args()

    df_a = pd.read_table(args.file_a, index_col=args.join)
    df_b = pd.read_table(args.file_b, index_col=args.join)

    joined = pd.merge(df_a, df_b, how='outer', suffixes=['_left', '_right'], left_index=True, right_index=True)
    output = args.output if args.output is not None else sys.stdout
    joined.to_csv(args.output, sep='\t')
