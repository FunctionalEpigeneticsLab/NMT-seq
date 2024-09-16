#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import argparse

def get_merged_summary(path_trim,path_map,path_allc,path_cov):
    df_trim = pd.read_csv(path_trim, sep='\t')
    df_map = pd.read_csv(path_map, sep='\t')
    df_allc = pd.read_csv(path_allc, sep='\t')
    df_cov=pd.read_csv(path_cov,sep='\t',header=None)
    df_cov.columns=['barcode','coverage_unfiltered','depth_unfiltered','coverage_filtered','depth_filtered']

    # merge dfs on barcode
    df_all = df_trim.merge(df_map, on='barcode', how='right')
    df_all['uniq_map_rate'] = df_all['mapped_uniq_reads']/df_all['trimmed_reads']*100
    df_all['dup_rate'] = (df_all['mapped_uniq_reads']-df_all['dedup_uniq_reads'])/df_all['mapped_uniq_reads']*100
    df_all['filtered_perc'] = (df_all['dedup_uniq_reads']-df_all['filtered_uniq_reads'])/df_all['dedup_uniq_reads']*100
    df_all['filtered_uniq_to_trimmed'] = df_all['filtered_uniq_reads']/df_all['trimmed_reads']*100
    df_all = df_all.merge(df_cov, on='barcode', how='left').merge(df_allc, on='barcode', how='right')
    return df_all

def main():
    parser = argparse.ArgumentParser(description="Merge summary stat files")
    parser.add_argument('--input_dir', default='./summary', help='Directory containing the summary files')
    parser.add_argument('--file_trim', default='trimmed_reads_summary.all_barcodes.tsv')
    parser.add_argument('--file_map', default='uniq_mapped_reads_summary.all_barcodes.tsv')
    parser.add_argument('--file_allc', default='allc_summary.tsv')
    parser.add_argument('--file_cov', default='combined.coverage.genome.wide.tsv')
    parser.add_argument('--out', default='./summary/metadata.allCells.csv', help='Output file path')
    args = parser.parse_args()

    summary_dir = args.input_dir
    file_trim = args.file_trim
    file_map = args.file_map
    file_allc = args.file_allc
    file_cov = args.file_cov
    outfile = args.out

    path_trim = os.path.join(summary_dir, file_trim)
    path_map = os.path.join(summary_dir, file_map)
    path_allc = os.path.join(summary_dir, file_allc)
    path_cov = os.path.join(summary_dir, file_cov)

    df_all = get_merged_summary(path_trim,path_map,path_allc,path_cov)
    df_all.round(4).to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
