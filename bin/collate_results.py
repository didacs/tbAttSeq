#!/usr/bin/env python

import argparse
import pandas as pd

def merge_csv_files_and_rearrange_columns(csv_files, sample_names, output_excel):
    if len(csv_files) != len(sample_names):
        raise ValueError("The number of CSV files and sample names must match.")

    # Placeholder for merged DataFrame
    merged_df = None

    # Iterate through CSV files and their corresponding sample names
    for csv_file, sample_name in zip(csv_files, sample_names):
        # Load CSV into DataFrame
        df = pd.read_csv(csv_file)

        # Update column names to include sample name, except for 'index' and 'attB_sequence'
        df.rename(columns=lambda x: x if x in ['index', 'attB_sequence'] else f"{sample_name}_{x}", inplace=True)

        # Merge DataFrames
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=['index', 'attB_sequence'], how='outer')

    # Generate new column order: 'index', 'attB_sequence' followed by sorted remaining columns
    base_columns = ['index', 'attB_sequence']
    data_columns = sorted([col for col in merged_df.columns if col not in base_columns],
                          key=lambda x: x.split('_')[-2])  # Sort by suffix
    new_column_order = base_columns + data_columns

    # Reorder the DataFrame columns
    merged_df = merged_df[new_column_order]

    # Write the merged and rearranged DataFrame to an Excel file
    merged_df.to_csv(output_excel, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple CSV files and rearrange columns based on suffix.")
    parser.add_argument("-f", "--files", nargs="+", required=True, help="The path to the CSV files to be merged.")
    parser.add_argument("-n", "--names", nargs="+", required=True, help="The corresponding sample names for each CSV file.")
    parser.add_argument("-o", "--output", required=True, help="The path to the output Excel file.")

    args = parser.parse_args()

    merge_csv_files_and_rearrange_columns(args.files, args.names, args.output)