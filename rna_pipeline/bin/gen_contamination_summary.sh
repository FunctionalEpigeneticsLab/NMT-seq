#!/bin/bash

# Define the list of species
species_list=("Homo sapiens" "Cyprinus" "Takifugu" "Oncorhynchus" "Salmo " "Coregonus" "Sprattus" "Danio" "Delftia" "Clostridium" "Aeromonas" "Comamonas" "Escherichia" "Acidovorax" "Bacillus" "Pseudomonas" "Pseudomonadaceae" "Carboxydocella" "Schlegelella" "Stutzerimonas" "Naegleria" "Sphingopyxis" "Azospira" "Cutibacterium" "Dechloromonas" "Bacterium" "Shigella" "Niallia" "Digitaria" "Malassezia" "Rhodoferax" "Daphnia" "Exidia" "Pelosinus" "Rhodococcus" "Corynebacterium" "Janthinobacterium" "Diaphorobacter" "Halopseudomonas" "Cladosporium")

# Define the output file
output_file="blast_per_plate_species_summary.tsv"

# Prepare the header for the output file
echo -n "Plate" > "$output_file"
for species in "${species_list[@]}"; do
    echo -n -e "\t$species" >> "$output_file"
done
# Add the header for the total hits column
echo -e "\twith_hits_total" >> "$output_file"

# Find all relevant files and process them
find . -type f -name "unmapped.blastn_out_count.summary" | while read file; do

    # Extract the substring from the file path
    extracted_string=$(echo "$file" | sed -n 's/\(GC[^_]*\).*/\1/p')
    echo -n "$extracted_string" >> "$output_file"
    
    # Initialize the variable to store the sum of all hits
    all_hits_total=0
    
    # For each species, calculate the total and append to the line
    for species in "${species_list[@]}"; do
        total=$(grep "$species" "$file" | cut -f2 | awk '{Total=Total+$1} END {print Total}')
        # If total is empty, print 0
        total=${total:-0}
        echo -n -e "\t$total" >> "$output_file"
        # Add to the all_hits_total
        all_hits_total=$((all_hits_total + total))
    done
    
    # Append the total hits for the file at the end of the line
    echo -e "\t$all_hits_total" >> "$output_file"
done
