#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a Venn diagram from QC data.')
parser.add_argument('qc_file', type=str, help='Path to the QC CSV file')
parser.add_argument('output_dir', type=str, help='Directory to save the Venn diagram image')

# Parse arguments
args = parser.parse_args()

# Ensure the output directory exists
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Read the QC file
df = pd.read_csv(args.qc_file)
df.rename(columns={'Unnamed: 0': 'Sample'}, inplace=True)

# Extract categories
Pass = df[df['QC'].str.contains('Pass')]['Sample'].tolist()
Low_nCount = df[df['QC'].str.contains('Low_Uniq_map_readpairs')]['Sample'].tolist()
High_ERCC = df[df['QC'].str.contains('High_ERCC')]['Sample'].tolist()
Low_nFeature = df[df['QC'].str.contains('Low_nFeature')]['Sample'].tolist()
High_MT = df[df['QC'].str.contains('High_MT')]['Sample'].tolist()

# Venn diagram
subsets = [set(High_ERCC), set(Low_nFeature), set(High_MT)]
venn = venn3(subsets=subsets,
             set_labels=('High fraction of\nERCC reads', #'Low number of\nmapped reads',
                         'Low number of\ndetected genes',
                         'High fraction of\nmitochondrial reads'))
for text in venn.set_labels:
    text.set_fontsize(8)
venn3_circles(subsets=subsets, linestyle="dashed", linewidth=2)

plt.tight_layout()

# Save the Venn diagram
output_png = os.path.join(args.output_dir, 'venn_qc.1.png')
output_pdf = os.path.join(args.output_dir, 'venn_qc.1.pdf')
plt.savefig(output_png, dpi=200)
plt.savefig(output_pdf)

# Another Venn diagram
plt.clf()
subsets = [set(Low_nCount), set(Low_nFeature), set(High_MT)]
venn = venn3(subsets=subsets,
             set_labels=('Low number of\nmapped reads',
                         'Low number of\ndetected genes',
                         'High fraction of\nmitochondrial reads'))
for text in venn.set_labels:
    text.set_fontsize(8)
venn3_circles(subsets=subsets, linestyle="dashed", linewidth=2)

plt.tight_layout()

# Save the Venn diagram
output_png = os.path.join(args.output_dir, 'venn_qc.2.png')
output_pdf = os.path.join(args.output_dir, 'venn_qc.2.pdf')
plt.savefig(output_png, dpi=200)
plt.savefig(output_pdf)

