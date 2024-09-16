#!/bin/bash

inputFile="$1"
outputFile="$2"

awk -F'\t' '{
    if($2 in counts) {
        ids[$2] = ids[$2] "," $1;
    } else {
        ids[$2] = $1;
    }
    counts[$2]++;
}
END {
    for (key in counts) {
        printf "%s\t%d\t%s\n", key, counts[key], ids[key];
    }
}' "$inputFile" | sort -t$'\t' -k2,2nr > "$outputFile"
