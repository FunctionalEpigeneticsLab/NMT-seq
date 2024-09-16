#!/bin/bash

for i in $(ls *_allc.tsv.gz); do echo -e "$(echo "$i" | cut -f1 -d '.')\t$(pwd)/${i}" >> allc_table.txt; done
