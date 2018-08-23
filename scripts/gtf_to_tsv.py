import argparse
import pandas as pd
from gtfparse import read_gtf

def main(gtf, tsv):
    df = read_gtf(gtf)
    df.to_csv(tsv, sep = '\t')

def parse_arguments():
    parser = argparse.ArgumentParser(description = 'Convert a GTF file to a TSV file')
    parser.add_argument('gtf', type=str, help='Input GTF file')
    parser.add_argument('tsv', type=str, help='Output TSV file')
    args = parser.parse_args()
    return vars(args)

if __name__ == "__main__":
    args = parse_arguments()
    main(**args)
