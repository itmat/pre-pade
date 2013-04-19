from __future__ import print_function

import csv
import argparse

def int_transform(fieldnames):
    def transform(x):
        for fieldname in fieldnames:
            x[fieldname] = int(x[fieldname])
        return x
    return transform    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Filters a CSV file for columns that meet some criteria""")

    parser.add_argument('--where', '-w')
    parser.add_argument('input', type=file)
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    parser.add_argument('--int', action='append',
                        help='Specify fields that should be treated as ints')

    kwargs = { 
        'delimiter' : '\t'
    }

    args = parser.parse_args()
    print(args.int)
    transform = int_transform(args.int)

    reader = csv.DictReader(args.input, **kwargs)

    writer = csv.DictWriter(args.output, reader.fieldnames, **kwargs)

    writer.writeheader()

    f = eval("lambda x: " + args.where)

    for x in reader:
        x = transform(x)
        if f(x):
            writer.writerow(x)



    
