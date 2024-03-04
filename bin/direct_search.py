#!/usr/bin/env python

import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import tqdm

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_sequence(sequence, flank_left, flank_right):
    left_index = sequence.find(flank_left)
    right_index = sequence.find(flank_right, left_index + 1)
    if left_index != -1 and right_index != -1:
        return sequence[left_index + len(flank_left):right_index]
    return None

def process_fastq(file, attb_flank_left, attb_flank_right,attp_flank_left, attp_flank_right):
    attb_flank_left = attb_flank_left.upper()
    attb_flank_right = attb_flank_right.upper()
    attb_flank_left_rc = reverse_complement(attb_flank_left)
    attb_flank_right_rc = reverse_complement(attb_flank_right)

    attp_flank_left = attp_flank_left.upper()
    attp_flank_right = attp_flank_right.upper()
    attp_flank_left_rc = reverse_complement(attp_flank_left)
    attp_flank_right_rc = reverse_complement(attp_flank_right)

    sequences = {}
    with gzip.open(file, "rt") as handle:
        for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
            sequence = str(record.seq)
            target_seq = find_sequence(sequence, attb_flank_left, attb_flank_right)
            if not target_seq:
                target_seq = find_sequence(sequence, attb_flank_right_rc, attb_flank_left_rc)
                if target_seq:
                    target_seq = reverse_complement(target_seq)
            if not target_seq:
                target_seq = find_sequence(sequence, attb_flank_left, attp_flank_right)
                if not target_seq:
                    target_seq = find_sequence(sequence, attp_flank_right_rc, attb_flank_left_rc)
                    if target_seq:
                        target_seq = reverse_complement(target_seq)
            if not target_seq:
                target_seq = find_sequence(sequence, attp_flank_left, attb_flank_right)
                if not target_seq:
                    target_seq = find_sequence(sequence, attb_flank_right_rc, attp_flank_left_rc)
                    if target_seq:
                        target_seq = reverse_complement(target_seq)
            if target_seq:
                sequences[target_seq] = sequences.get(target_seq, 0) + 1
    return sequences


def cross_reference_amplicons(sequences, amplicons_file):
    amplicons = pd.read_csv(amplicons_file, sep="\t", names=["name", "sequence"])
    amplicons["count"] = amplicons["sequence"].map(sequences).fillna(0)
    print(amplicons)
    return amplicons

def calculate_recombination_percentage(amplicons, sample_name):
    recombination_data = []
    for index in amplicons['name'].str.extract(r'_(\d+)')[0].unique():
        subset = amplicons[amplicons['name'].apply(lambda x: x.count('_') <= 1) & amplicons['name'].str.contains(f'_{index}$')]
        print(subset)
        attL_count = subset[subset['name'].str.contains('attL')]['count'].sum()
        attR_count = subset[subset['name'].str.contains('attR')]['count'].sum()
        attB_count = subset[subset['name'].str.contains('attB')]['count'].sum()
        print(subset[subset['name'].str.contains('attB')]['sequence'])
        attB_sequence = subset[subset['name'].str.contains('attB')]['sequence'].reset_index(drop=True).iloc[0]
        max_l_r = max(attL_count, attR_count)
        recombination_percentage = 100 * max_l_r / (max_l_r + attB_count) if (max_l_r + attB_count) != 0 else 0
        recombination_data.append({
            "index": index, 
            "attB_sequence": attB_sequence, 
            "attL_count": attL_count, 
            "attR_count": attR_count, 
            "attB_count": attB_count, 
            "recombination %": recombination_percentage
        })
    
    recombination_data_df = pd.DataFrame(recombination_data)

    recombination_data_df.to_csv(f'{sample_name}_recombination_data.csv',index=False)

    return pd.DataFrame(recombination_data)

def main():
    parser = argparse.ArgumentParser(description="Process FASTQ file and calculate recombination percentages.")
    parser.add_argument("--fastq_file", help="Path to the .fastq.gz file")
    parser.add_argument("--attb_flank_left", help="Left flanking sequence (12 bp)")
    parser.add_argument("--attb_flank_right", help="Right flanking sequence (12 bp)")
    parser.add_argument("--attp_flank_left", help="Left flanking sequence (12 bp)")
    parser.add_argument("--attp_flank_right", help="Right flanking sequence (12 bp)")
    parser.add_argument("--amplicons_file", help="Path to the amplicons.txt file")
    parser.add_argument("--sample_name", help="sample name")
    args = parser.parse_args()

    sequences = process_fastq(args.fastq_file, args.attb_flank_left, args.attb_flank_right,args.attp_flank_left, args.attp_flank_right)
    amplicons = cross_reference_amplicons(sequences, args.amplicons_file)
    recombination = calculate_recombination_percentage(amplicons, args.sample_name)
    print(recombination)

if __name__ == "__main__":
    main()
