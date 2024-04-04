#!/usr/bin/env python
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--dsb_counts', help='dsb_counts', required=True, nargs = '+')
    parser.add_argument('--output', help='output CSV file', required=True)
    return parser.parse_args()

def collate_dsb_counts(dsb_counts):
    # Read each CSV file into separate DataFrames
    dfs = [pd.read_csv(csv) for csv in dsb_counts]

    # Concatenate the DataFrames into a single DataFrame
    df = pd.concat(dfs, ignore_index=True)

    # pivot df
    results = (df.pivot_table(
        index=['index','quant_window'],
        columns='sample_name',
        values=['support_dsb_total', 'support_dsb_%', 'total_reads'])
                .reset_index())

    # Rename columns
    base_columns = ['index', 'quant_window']
    results_columns = ['{}-{}'.format(col[1], col[0])
                       for col in results.columns
                       if not col[0] in base_columns]
    results.columns = base_columns + results_columns
    return results

if __name__ == "__main__":
    args = parse_arguments()
    results = collate_dsb_counts(args.dsb_counts)
    results.to_csv(args.output, index=False)
