#!/usr/bin/env python
import argparse
import collections
import os
import sys

import HTSeq
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sample_name', help='sample name', required=True)
    parser.add_argument('--bam', help='input bam file', required=True)
    parser.add_argument('--min_mapq', help='minimum read mapq value', default=30)
    parser.add_argument('--dinucleotide_position',
                        help='position of dinucleotide, 1-based (assumes the same position for all oligos)', required=True, type=int)
    return parser.parse_args()


def parse_bam(bam, min_mapq, dinucleotide_position):
    """
    Parses the bam file, filters unclipped reads (no soft-clipped bases at either end),
    and stores start and end positions.

    It then quantifies number of reads that start or end within a 6 bp quantification window
    centered around the dinucleotide. Reads that end within the quantification window support a left fragment,
    and reads that start within the quantification window support right fragment.

    Args:
        bam: alignment file of reads mapped against unrecombined substrate sequences
        min_mapq: minimum read mapq value [defafult: 30]
        dinucleotide_position: 1-based position of dinucleotide. Assumes the same position for all oligos
    """

    unclipped_start = HTSeq.GenomicArray("auto", stranded=False)
    unclipped_end = HTSeq.GenomicArray("auto", stranded=False)
    unclipped_bounds = HTSeq.GenomicArray("auto", stranded=False)
    coverage = HTSeq.GenomicArray("auto", stranded=False)

    total_reads = collections.Counter()
    total_unclipped = collections.Counter()
    support_dsb_left = collections.Counter()
    support_dsb_right = collections.Counter()

    # loop through each read
    for read in HTSeq.BAM_Reader(bam):
        # make sure read is aligned, is primary alignment, and mapq is above threshold
        if (not read.aligned or read.not_primary_alignment or read.aQual < min_mapq):
            continue

        total_reads[read.iv.chrom] += 1  # chromoseome is oligo id

        # make sure there are no soft-clipped bases at either end of the read
        if not (read.cigar[0].type == 'M' and read.cigar[-1].type == 'M'):
            # TODO: compute lenght of clipped sequence metrics
            continue

        # unclipped reads passed all filters
        total_unclipped[read.iv.chrom] += 1

        # positions are 0-based
        # keep start and end separately used for DSB table
        unclipped_start[HTSeq.GenomicPosition(read.iv.chrom, read.iv.start)] += 1
        unclipped_end[HTSeq.GenomicPosition(read.iv.chrom, read.iv.end)] += 1

        # start and end for all unclipped reads for read_boundaries table
        unclipped_bounds[HTSeq.GenomicPosition(read.iv.chrom, read.iv.start)] += 1
        unclipped_bounds[HTSeq.GenomicPosition(read.iv.chrom, read.iv.end)] += 1

        # used to compute the read coverage at each position
        coverage[read.iv] += 1

    # Define boundaries for a 6 bp quantification window around the dinucleotide
    di = dinucleotide_position - 1  # 0-based
    left_bound = di - 2
    right_bound = di + 4
    # count reads that end within the window (left fragment), reads that start within the window (right fragment)
    for chrom in unclipped_start.chrom_vectors:
        iv = HTSeq.GenomicInterval(chrom, left_bound, right_bound)
        for _, counts in unclipped_end[iv].steps():
            support_dsb_left[chrom] += int(counts)
        for _, counts in unclipped_start[iv].steps():
            support_dsb_right[chrom] += int(counts)

    # {sample_name}.dsb_counts.csv
    # name, total_reads, unclipped, quant_window, support_dsb_left, support_dsb_right, support_dsb_%
    total_reads = pd.DataFrame.from_dict(total_reads, orient='index', columns=['total_reads']).reset_index()
    total_unclipped = pd.DataFrame.from_dict(total_unclipped, orient='index', columns=['unclipped']).reset_index()
    support_dsb_left = pd.DataFrame.from_dict(support_dsb_left, orient='index',
                                              columns=['support_dsb_left']).reset_index()
    support_dsb_right = pd.DataFrame.from_dict(support_dsb_right, orient='index',
                                               columns=['support_dsb_right']).reset_index()
    df_dsb = pd.merge(total_reads, total_unclipped)
    df_dsb['quant_window'] = f"{left_bound+1}-{right_bound}" # 1-based closed-ended
    df_dsb = pd.merge(df_dsb, support_dsb_left)
    df_dsb = pd.merge(df_dsb, support_dsb_right)
    df_dsb['support_dsb_total'] = df_dsb['support_dsb_left'] + df_dsb['support_dsb_right']
    df_dsb['support_dsb_%'] = df_dsb.apply(lambda x: x.support_dsb_total / x.total_reads * 100, axis=1)

    # {sample_name}.read_boundaries.csv:
    # Number of unclipped reads starting or ending at each position along oligos
    l = []
    for iv, total_coverage in coverage.steps():
        for iiv, unclipped in unclipped_bounds[iv].steps():
            if total_coverage > 0:
                l.append({
                    'sample_name': args.sample_name,
                    'chrom': iiv.chrom,
                    'start': iiv.start + 1,  # make it 1-based
                    'fraction_unclipped': unclipped / total_coverage,
                    'unclipped': int(unclipped),
                    'total_coverage': int(total_coverage),
                })
    df_bounds = pd.DataFrame(l)
    
    return df_dsb, df_bounds


if __name__ == "__main__":
    args = parse_arguments()
    df_dsb, df_bounds = parse_bam(args.bam, args.min_mapq, args.dinucleotide_position)
    # save DSB table
    dsb_outfile = f"{args.sample_name}.dsb_counts.csv"
    df_dsb.to_csv(dsb_outfile, index=0)
    # save read boundaries table
    bounds_outfile = f"{args.sample_name}.read_boundaries.csv"
    df_bounds.to_csv(bounds_outfile, index=0)
