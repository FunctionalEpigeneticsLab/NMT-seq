cd $vc_workdir
touch "$cat_output"

# Iterate over all files matching the pattern, excluding the output file
for file in $vc_workdir/*_VariantCalling.tsv; do
    # Exclude the output file from processing
    if [ "$file" != "$cat_output" ]; then
        # Extract the common identifier from the filename
        identifier=$(basename "$file" | sed 's/_VariantCalling\.tsv//')
        echo "Processing file: $file, Identifier: $identifier"
        
        # Add the identifier as a new column and store in a temporary file
        awk -v id="$identifier" '{print id, $0}' "$file" > "$file.tmp"
        
        # Append the modified file to the concatenated output file
        cat "$file.tmp" >> "$cat_output"
        
        # Remove the temporary file
        rm "$file.tmp"
    fi
done
