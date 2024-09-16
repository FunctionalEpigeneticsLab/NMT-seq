#!/usr/bin/env python3

import gzip
import os
import argparse
from collections import defaultdict
from multiprocessing import Pool
import glob

def count_barcodes(fq1_path):
    pattern = "_val_"
    prefix = fq1_path.split(pattern)[0].split('/')[-1]

    barcode_counts = defaultdict(int)

    with gzip.open(fq1_path, 'rt') as fq1:
        while True:
            try:
                r1_lines = [fq1.readline().strip() for _ in range(4)]

                if len(r1_lines) < 4:
                    break

                read_name = r1_lines[0]
                if read_name == '':
                    break

                parts = read_name.split(':')

                if len(parts) < 2:
                    print(f"Warning: read name does not contain a barcode: {read_name}")
                    continue

                barcode = parts[0].split('@')[1].replace('+', '_')

                barcode_counts[barcode] += 1

            except StopIteration:
                break

    prefixed_counts = [(f"{prefix}.{barcode}", count) for barcode, count in barcode_counts.items()]

    return prefixed_counts

def write_counts_to_tsv(results, output_file):
    all_counts = [item for sublist in results for item in sublist]

    all_counts.sort(key=lambda x: x[1], reverse=True)

    with open(output_file, 'w') as out_tsv:
        out_tsv.write("barcode\ttrimmed_reads\n")
        for barcode, count in all_counts:
            out_tsv.write(f"{barcode}\t{count*2}\n") # times 2 to have the count of reads instead of read pairs

def process_file(file):
    return count_barcodes(file)

def main():
    parser = argparse.ArgumentParser(description="Count occurrences of barcodes in multiple fastq files and output a concatenated TSV")
    parser.add_argument('--input_dir', required=True, help='Directory containing the input fastq files')
    parser.add_argument('--out', default='.', help='Output the concatenated barcode counts TSV file')
    args = parser.parse_args()

    output_file = args.out 
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fq1_files = glob.glob(os.path.join(args.input_dir, '*/*_val_1.fq.gz'))

    with Pool() as pool:
        results = pool.map(process_file, fq1_files)

    write_counts_to_tsv(results, output_file)

if __name__ == "__main__":
    main()

