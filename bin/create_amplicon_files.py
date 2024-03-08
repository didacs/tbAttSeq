#!/usr/bin/env python
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert CSV to FASTA for amplicons')
    parser.add_argument('--attb_list', type=str, help='Path to the input CSV file')
    parser.add_argument('--attp', type=str, help='Path to the attP oligo fasta file.')
    parser.add_argument('--output', type=str, help='Path to the output FASTA file')
    return parser.parse_args()

def replace_sequence(original, replacement, target):
    position = target.find(original)
    if position == -1:
        return target  # Original sequence not found
    return target[:position] + replacement + target[position + len(original):]

def generate_fasta(df, fasta_file, attp_oligo_fn):

    attL_dict = {}
    attR_dict = {}
    attB_dict = {}

    for i, row in df.iterrows():
        seq = row['sequence'].upper()
        name = row['name']
        seq_5 = seq[0:24]
        seq_3 = seq[24:]

        attB = seq
        attL = seq_5 + "CTCAGTGGTGTACGGTACAAACCCA"
        attR = "GTGGTTTGTCTGGTCAACCACCGCGGT" + seq_3

        attB_dict.setdefault(attB, []).append(name)
        attL_dict.setdefault(attL, []).append(name)
        attR_dict.setdefault(attR, []).append(name)

        print(attB_dict)

    with open(fasta_file, 'w') as outfile:
        for attB, indices in attB_dict.items():
            name = '_'.join(map(str, indices))
            outfile.write(f"attB_{name}\t{attB}\n")

        for attL, indices in attL_dict.items():
            name = '_'.join(map(str, indices))
            outfile.write(f"attL_{name}\t{attL}\n")

        for attR, indices in attR_dict.items():
            name = '_'.join(map(str, indices))
            outfile.write(f"attR_{name}\t{attR}\n")

if __name__ == "__main__":
    args = parse_arguments()
    input_df = pd.read_csv(args.attb_list)
    input_df['sequence'] = input_df['sequence'].str.upper()
    input_df['name'] = input_df['name'].str.upper()
    input_df = input_df.drop_duplicates()
    generate_fasta(input_df, args.output, args.attp)

