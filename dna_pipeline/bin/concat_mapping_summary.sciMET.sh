#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_file>"
    exit 1
fi

input_folder="$1"
output_file="$2"

mkdir -p "$(dirname "$output_file")"

header="barcode\tmapped_uniq_reads\tdedup_uniq_reads\tfiltered_uniq_reads"

temp_file="${output_file}.tmp"
for file in "$input_folder"/*.mapping_summary.tsv; do
    if [ -f "$file" ]; then
        cat "$file" >> "$temp_file"
    fi
done

sort -k4nr "$temp_file" > "$output_file"

sed -i "1i$header" "$output_file"

rm "$temp_file"

