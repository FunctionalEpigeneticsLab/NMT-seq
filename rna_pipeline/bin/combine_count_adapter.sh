#!/bin/bash

# Define the combined file
combined_file="adapter_counts_combined.tsv"

# Write column names line to the combined file
echo -e "Sample\tBatch\traw_read_count\tadapter_f_middle" > "$combined_file"

# Loop through all tsv files
for tsv_file in */COUNT_ADAPTER/*count_adapter.tsv; do
    # Extract sample name from the first line
    sample_name=$(awk 'NR==1 {print $1}' "$tsv_file")

    # Extract batch from the sample name
    batch=$(echo "$sample_name" | cut -d '_' -f 1)

    # Extract column values from the second line onwards and join them with tabs
    values=$(tail -n +2 "$tsv_file" | cut -f 2- | tr '\n' '\t')

    # Remove the last tab from the values
    values=$(echo "$values" | sed 's/\t$//')

    # Write combined line to the combined file
    echo -e "${sample_name}\t${batch}\t${values}" >> "$combined_file"
done

# Sort combined data by sample name
sort -k1,1 "$combined_file" -o "$combined_file"

echo "Combined file generated: $combined_file"

