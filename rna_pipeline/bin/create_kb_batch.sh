#!/bin/bash

echo -n > batch.txt

for file in *_R1_001.fastq.gz; do
    base=$(echo "$file" | sed 's/_R1_001.fastq.gz//')

    r1_file="${base}_R1_001.fastq.gz"
    r2_file="${base}_R2_001.fastq.gz"

    if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
        sample_id=$(echo "$base" | sed 's/_R[12]_001.fastq.gz//')

        echo -e "$sample_id\t$r1_file\t$r2_file" >> batch.txt
    else
        echo "Warning: paired file for $file not found, skipping"
    fi
done

