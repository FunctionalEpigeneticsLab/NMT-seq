#!/usr/bin/env python3

import gzip
import sys
import re
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import argparse

def count_contexts(input_file):
    context_patterns = {
        'WCGN': re.compile(r'[AT]CG.'),
        'GCHN': re.compile(r'GC[ATC].'),
        'HCHN': re.compile(r'[ATC]C[ATC].'),
        'CCGN': re.compile(r'CCG.'),
        'GCGN': re.compile(r'GCG.')
    }

    context_counts = {key: {'cov': 0, 'mc': 0} for key in context_patterns}
    general_counts = {'NCGN': {'cov': 0, 'mc': 0}, 'NCHN': {'cov': 0, 'mc': 0}}

    ncgn_pattern = re.compile(r'.CG.')
    nchn_pattern = re.compile(r'.C[ATC].')

    with gzip.open(input_file, 'rt') as f:
        for line in f:
            columns = line.strip().split('\t')
            context = columns[3]
            mc = int(columns[4])
            cov = int(columns[5])

            for name, pattern in context_patterns.items():
                if pattern.match(context):
                    context_counts[name]['cov'] += cov
                    context_counts[name]['mc'] += mc

            if ncgn_pattern.match(context):
                general_counts['NCGN']['cov'] += cov
                general_counts['NCGN']['mc'] += mc

            if nchn_pattern.match(context):
                general_counts['NCHN']['cov'] += cov
                general_counts['NCHN']['mc'] += mc

    return context_counts, general_counts

def compute_methylation_frequency(mc, cov):
    return round(mc / cov, 4) if cov > 0 else 0.0

def process_file(input_file):
    context_counts, general_counts = count_contexts(input_file)

    data = [
        os.path.basename(input_file).split('.allc')[0],  # Barcode
        general_counts['NCGN']['cov'],
        general_counts['NCGN']['mc'],
        compute_methylation_frequency(general_counts['NCGN']['mc'], general_counts['NCGN']['cov']),
        general_counts['NCHN']['cov'],
        general_counts['NCHN']['mc'],
        compute_methylation_frequency(general_counts['NCHN']['mc'], general_counts['NCHN']['cov']),
        context_counts['WCGN']['cov'],
        context_counts['WCGN']['mc'],
        compute_methylation_frequency(context_counts['WCGN']['mc'], context_counts['WCGN']['cov']),
        context_counts['GCHN']['cov'],
        context_counts['GCHN']['mc'],
        compute_methylation_frequency(context_counts['GCHN']['mc'], context_counts['GCHN']['cov']),
        context_counts['HCHN']['cov'],
        context_counts['HCHN']['mc'],
        compute_methylation_frequency(context_counts['HCHN']['mc'], context_counts['HCHN']['cov']),
        context_counts['CCGN']['cov'],
        context_counts['CCGN']['mc'],
        compute_methylation_frequency(context_counts['CCGN']['mc'], context_counts['CCGN']['cov']),
        context_counts['GCGN']['cov'],
        context_counts['GCGN']['mc'],
        compute_methylation_frequency(context_counts['GCGN']['mc'], context_counts['GCGN']['cov']),
    ]
    return data


def main(input_dir, output_file):
    header = [
        'barcode', 'NCGN_cov', 'NCGN_mc', 'NCGN_perc', 'NCHN_cov', 'NCHN_mc', 'NCHN_perc',
        'WCGN_cov', 'WCGN_mc', 'WCGN_perc', 'GCHN_cov', 'GCHN_mc', 'GCHN_perc',
        'HCHN_cov', 'HCHN_mc', 'HCHN_perc', 'CCGN_cov', 'CCGN_mc', 'CCGN_perc',
        'GCGN_cov', 'GCGN_mc', 'GCGN_perc'
    ]

    all_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.allc.tsv.gz')]

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_file, all_files))

    df = pd.DataFrame(results, columns=header)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mC summary stats from allc files")
    parser.add_argument('--input_dir', required=True, help="Directory containing .allc.tsv.gz files")
    parser.add_argument('--out', required=True, help="Output file path for the summary")

    args = parser.parse_args()

    input_dir = args.input_dir
    output_file = args.out

    main(input_dir, output_file)
