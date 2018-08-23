import argparse
import pandas as pd
from gtfparse import read_gtf

def main(gtf, tsv):
    gtf = args.gtf
    tsv = args.tsv
    df = read_gtf(gtf)
    df.to_csv(tsv, sep = '\t')

def parse_arguments():
    parser = argparse.ArgumentParser(description = 'Convert a GTF file to a TSV file')
    parser.add_argument('gtf', type=int, help='Input GTF file')
    parser.add_argument('tsv', type=int, help='Output TSV file')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(**args)
