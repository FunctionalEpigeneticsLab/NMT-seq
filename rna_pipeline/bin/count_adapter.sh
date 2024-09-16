#!/bin/bash

if [ -z "$1" ]; then
    echo -e "Usage: $0 <filename>"
    exit 1
fi

fq="$1"

sample=$(basename "$fq" | sed 's/_001.*$//')

zcat "$fq" > ${sample}.temp.fq

# Extract substrings

total=$(( $(cat "${sample}.temp.fq" | wc -l) / 4 ))

a_f_count=$(grep -P '^(.{0,122})CTGTCTCTTATACACATCT' ${sample}.temp.fq | wc -l) #{0,131}

# Write to tsv

echo -e "\t${sample}" > "${sample}.count_adapter.tsv"

echo -e "read_count\t$total" >> "${sample}.count_adapter.tsv"

echo -e "adapter_f_middle\t$a_f_count" >> "${sample}.count_adapter.tsv"

rm ${sample}.temp.fq
echo "Counts written to: ${sample}.count_adapter.tsv"
