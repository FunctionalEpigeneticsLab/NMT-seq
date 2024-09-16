import pandas as pd
import os
import glob

# Pattern to match the RSEM result files
pattern = '*/RSEM/*.Aligned.toTranscriptome.out.rsem.genes.results'

# Initialize an empty DataFrame for merged results
merged_df = pd.DataFrame()
# List to store sample names
sample_names = []

# Iterate over each file that matches the pattern
for file in glob.glob(pattern):
    # Extract sample name from the file name
    sample_name = os.path.basename(file).split('.Aligned.toTranscriptome.out.rsem.genes.results')[0]
    sample_name = sample_name.split('_S')[0]
    print(sample_name)
    sample_names.append(sample_name)

    # Read the file
    df = pd.read_csv(file, sep='\t', usecols=['gene_id', 'TPM'])

    # Rename the 'TPM' column to the sample name
    df.rename(columns={'TPM': sample_name}, inplace=True)

    # Merge with the main DataFrame
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')

# Sort the sample names
sample_names.sort()

# Reorder the columns in the merged DataFrame
merged_df = merged_df[['gene_id'] + sample_names]

# Save the merged DataFrame to a new CSV file
#merged_df.to_csv('merged_TPM_values.csv', index=False)
# Save the merged DataFrame to a new TSV file
merged_df.to_csv('RSEM.TPM.merged.all.tsv', sep='\t', index=False)
