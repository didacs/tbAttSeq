#!/usr/bin/env python
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert CSV to FASTA adding left and righ flank sequences '
                                                 'to re-create full-length oligos')
    parser.add_argument('--attb_oligos', type=str, help='Path to the input CSV file')
    parser.add_argument('--attb_flank_left', type=str, help='Left flank sequence')
    parser.add_argument('--attb_flank_right', type=str, help='Right flank sequence')
    parser.add_argument('--output_fasta', type=str, help='Output fasta file')
    return parser.parse_args()


def write_fasta(df, output_fasta):
    with open(output_fasta, 'w') as fa:
        for _, row in df.iterrows():
            sequence = row['full_seq'].upper()
            name = row['name']
            fa.write(f'>{name}\n{sequence}\n')

if __name__ == '__main__':
    args = parse_arguments()
    df = pd.read_csv(args.attb_oligos)
    df["full_seq"] = args.attb_flank_left + df["sequence"] + args.attb_flank_right
    write_fasta(df, args.output_fasta)
