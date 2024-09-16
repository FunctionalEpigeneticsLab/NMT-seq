#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_file>"
    exit 1
fi

input_dir="$1"
output_file="$2"

mkdir -p "$(dirname "$output_file")"

header="barcode\tNCGN_cov\tNCGN_mc\tNCGN_perc\tNCHN_cov\tNCHN_mc\tNCHN_perc\tWCGN_cov\tWCGN_mc\tWCGN_perc\tGCHN_cov\tGCHN_mc\tGCHN_perc\tHCHN_cov\tHCHN_mc\tHCHN_perc\tCCGN_cov\tCCGN_mc\tCCGN_perc\tGCGN_cov\tGCGN_mc\tGCGN_perc"
echo -e "$header" > "$output_file"

for file in "$input_dir"/*.allc.mC_stat.tsv; do
    if [ -f "$file" ]; then
        awk -F'\t' -v OFS='\t' '{print substr($1, 1, index($1, ".allc")-1), $2, $3, $4, $5, $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' "$file" >> "$output_file"
    fi
done

