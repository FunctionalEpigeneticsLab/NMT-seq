#!/usr/bin/env python3

import gzip
import sys
import re

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
    return round(mc / cov, 2) if cov > 0 else 0.0

def main(input_file, output_file):
    context_counts, general_counts = count_contexts(input_file)

    with open(output_file, 'w') as out_f:
        out_f.write(f"{input_file}")
        for key in ['NCGN', 'NCHN']:
            cov = general_counts[key]['cov']
            mc = general_counts[key]['mc']
            freq = compute_methylation_frequency(mc, cov)
            out_f.write(f"\t{cov}\t{mc}\t{freq:.4f}")
        
        for key in ['WCGN', 'GCHN', 'HCHN', 'CCGN', 'GCGN']:
            cov = context_counts[key]['cov']
            mc = context_counts[key]['mc']
            freq = compute_methylation_frequency(mc, cov)
            out_f.write(f"\t{cov}\t{mc}\t{freq:.4f}")
        
        out_f.write("\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: extract_mC_stat_from_allc.py input_allc_file.tsv.gz output_file.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)

