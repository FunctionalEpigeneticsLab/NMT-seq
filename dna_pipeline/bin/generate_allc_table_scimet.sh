#!/bin/bash
rm -f allc_table.txt
for i in $(ls *_allc.tsv.gz); do
    basename=$(basename "$i")
    extracted_part=$(echo "$basename" | sed -n 's/.*nsrt\.\([^\.]*\)\.sc_allc\.tsv\.gz/\1/p')
    echo -e "${extracted_part}\t$(pwd)/${i}" >> allc_table.txt
done

