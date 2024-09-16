#!/bin/bash

# Define the combined file
combined_file="sequence_counts_combined.tsv"

# Write column names line to the combined file
#echo -e "Sample\tBatch\tread_count\tadapter_f\tadapter_r\tprimer_f\tprimer_r\tadapter_sub1_f\tadapter_sub1_r\tadapter_sub2_f\tadapter_sub2_r\tprimer_sub1_f\tprimer_sub1_r\tprimer_sub2_f\tprimer_sub2_r\tprimer_sub3_f\tprimer_sub3_r\tpolyG10\tpolyA10\tpolyT10\tpolyG50\tpolyA50\tpolyT50" > "$combined_file"
echo -e "Sample\tBatch\tread_count\tadapter_r\tprimer_f\tprimer_r\tadapter_sub1_f\tadapter_sub1_r\tadapter_sub2_f\tadapter_sub2_r\tprimer_sub1_f\tprimer_sub1_r\tprimer_sub2_f\tprimer_sub2_r\tprimer_sub3_f\tprimer_sub3_r\tpolyG10\tpolyA10\tpolyT10\tpolyG50\tpolyA50\tpolyT50" > "$combined_file"

# Loop through all tsv files
for tsv_file in */COUNT_SEQ/*count_seq.tsv; do
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

