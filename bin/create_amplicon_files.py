#!/usr/bin/env python
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert CSV to FASTA for amplicons')
    parser.add_argument('--attb_list', type=str, help='Path to the input CSV file')
    parser.add_argument('--attp', type=str, help='Path to the attP oligo fasta file.')
    parser.add_argument('--output', type=str, help='Path to the output FASTA file')
    parser.add_argument('--method', type=str, help='method (cs2 or direct)')
    return parser.parse_args()

def replace_sequence(original, replacement, target):
    position = target.find(original)
    if position == -1:
        return target  # Original sequence not found
    return target[:position] + replacement + target[position + len(original):]

def generate_fasta(df, fasta_file, attp_oligo_fn, method):

    attL_dict = {}
    attR_dict = {}
    attB_dict = {}

    if method == 'cs2':
        with open(attp_oligo_fn,'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                attp_oligo = line.strip()
                break
        for i, row in df.iterrows():
            seq = row['sequence'].upper()
            oligo = row['oligo'].upper()
            seq_5 = seq[0:24]
            seq_3 = seq[24:]

            attB = oligo
            attL = seq_5 + "CTCAGTGGTGTACGGTACAAACCCA"
            attR = "GTGGTTTGTCTGGTCAACCACCGCGGT" + seq_3

            modified_oligo_L = replace_sequence(seq, attL, oligo)
            modified_oligo_R = replace_sequence(seq, attR, oligo)

            attB_dict.setdefault(attB, []).append(i)
            attL_dict.setdefault(modified_oligo_L, []).append(i)
            attR_dict.setdefault(modified_oligo_R, []).append(i)

        with open(fasta_file, 'w') as outfile:
            outfile.write(f"attP\t{attp_oligo}\n")
            for attB, indices in attB_dict.items():
                name = '_'.join(map(str, indices))
                outfile.write(f"attB_{name}\t{attB}\n")

            for attL, indices in attL_dict.items():
                name = '_'.join(map(str, indices))
                outfile.write(f"attL_{name}\t{attL}\n")

            for attR, indices in attR_dict.items():
                name = '_'.join(map(str, indices))
                outfile.write(f"attR_{name}\t{attR}\n")

    elif method == 'direct':
        for i, row in df.iterrows():
            seq = row['sequence'].upper()
            seq_5 = seq[0:24]
            seq_3 = seq[24:]

            attB = seq
            attL = seq_5 + "CTCAGTGGTGTACGGTACAAACCCA"
            attR = "GTGGTTTGTCTGGTCAACCACCGCGGT" + seq_3

            attB_dict.setdefault(attB, []).append(i)
            attL_dict.setdefault(attL, []).append(i)
            attR_dict.setdefault(attR, []).append(i)

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
    input_df['oligo'] = input_df['oligo'].str.upper()
    input_df = input_df.drop_duplicates()
    generate_fasta(input_df, args.output, args.attp, args.method)
