#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <filename> <sample>"
    exit 1
fi

filename="$1"
sample="$2"

outfile="$filename.genome.wide"

awk -v sample="$sample" '
BEGIN {
    totalLength = 0;
    totalWeightedCoverage = 0;
    totalWeightedDepth = 0;
}
NR > 1 { # Skip the header line
    chrLength = $3 - $2 + 1;
    weightedCoverage = chrLength * $6;
    weightedDepth = chrLength * $7;
    totalLength += chrLength;
    totalWeightedCoverage += weightedCoverage;
    totalWeightedDepth += weightedDepth;
}
END {
    genomeWideCoverage = totalWeightedCoverage / totalLength;
    genomeWideDepth = totalWeightedDepth / totalLength;
    printf "%s\t%.2f\t%.6f\n", sample, genomeWideCoverage, genomeWideDepth > "'$outfile'";
}
' "$filename"

