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

#a_f_count=$(grep 'CTGTCTCTTATACACATCT' ${sample}.temp.fq | wc -l)
a_r_count=$(grep 'AGATGTGTATAAGAGACAG' ${sample}.temp.fq | wc -l)
p_f_count=$(grep 'AAGCAGTGGTATCAACGCAGAGTAC' ${sample}.temp.fq | wc -l)
p_r_count=$(grep 'GTACTCTGCGTTGATACCACTGCTT' ${sample}.temp.fq | wc -l)

asub1_f_count=$(grep 'CTGTCTCTTATA' ${sample}.temp.fq | wc -l)
asub1_r_count=$(grep 'TATAAGAGACAG' ${sample}.temp.fq | wc -l)
asub2_f_count=$(grep 'CACATCT' ${sample}.temp.fq | wc -l)
asub2_r_count=$(grep 'AGATGTG' ${sample}.temp.fq | wc -l)

psub1_f_count=$(grep 'GTACTCTGCGTTGATAC' ${sample}.temp.fq | wc -l)
psub1_r_count=$(grep 'GTACTCTGCGTTGATAC' ${sample}.temp.fq | wc -l)
psub2_f_count=$(grep 'CCCAT' ${sample}.temp.fq | wc -l)
psub2_r_count=$(grep 'ATGGG' ${sample}.temp.fq | wc -l)
psub3_f_count=$(grep 'AAGCAGTG' ${sample}.temp.fq | wc -l)
psub3_r_count=$(grep 'CACTGCTT' ${sample}.temp.fq | wc -l)

g10_count=$(grep 'GGGGGGGGGG' ${sample}.temp.fq | wc -l)
a10_count=$(grep 'AAAAAAAAAA' ${sample}.temp.fq | wc -l)
t10_count=$(grep 'TTTTTTTTTT' ${sample}.temp.fq | wc -l)
g50_count=$(grep 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG' ${sample}.temp.fq | wc -l)
a50_count=$(grep 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' ${sample}.temp.fq | wc -l)
t50_count=$(grep 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT' ${sample}.temp.fq | wc -l)

# Write to tsv

echo -e "\t${sample}" > "${sample}.count_seq.tsv"

echo -e "read_count\t$total" >> "${sample}.count_seq.tsv"

#echo -e "adapter_f\t$a_f_count" >> "${sample}.count_seq.tsv"
echo -e "adapter_r\t$a_r_count" >> "${sample}.count_seq.tsv"
echo -e "primer_f\t$p_f_count" >> "${sample}.count_seq.tsv"
echo -e "primer_r\t$p_r_count" >> "${sample}.count_seq.tsv"

echo -e "adapter_sub1_f\t$asub1_f_count" >> "${sample}.count_seq.tsv"
echo -e "adapter_sub1_r\t$asub1_r_count" >> "${sample}.count_seq.tsv"
echo -e "adapter_sub2_f\t$asub2_f_count" >> "${sample}.count_seq.tsv"
echo -e "adapter_sub2_r\t$asub2_r_count" >> "${sample}.count_seq.tsv"

echo -e "primer_sub1_f\t$psub1_f_count" >> "${sample}.count_seq.tsv"
echo -e "primer_sub1_r\t$psub1_r_count" >> "${sample}.count_seq.tsv"
echo -e "primer_sub2_f\t$psub2_f_count" >> "${sample}.count_seq.tsv"
echo -e "primer_sub2_r\t$psub2_r_count" >> "${sample}.count_seq.tsv"
echo -e "primer_sub3_f\t$psub3_f_count" >> "${sample}.count_seq.tsv"
echo -e "primer_sub3_r\t$psub3_r_count" >> "${sample}.count_seq.tsv"

echo -e "polyG10\t$g10_count" >> "${sample}.count_seq.tsv"
echo -e "polyA10\t$a10_count" >> "${sample}.count_seq.tsv"
echo -e "polyT10\t$t10_count" >> "${sample}.count_seq.tsv"
echo -e "polyG50\t$g50_count" >> "${sample}.count_seq.tsv"
echo -e "polyA50\t$a50_count" >> "${sample}.count_seq.tsv"
echo -e "polyT50\t$t50_count" >> "${sample}.count_seq.tsv"

rm ${sample}.temp.fq
echo "Counts written to: ${sample}.count_seq.tsv"

