#!/bin/bash

passing_cells_csv=""
input_dir=""
output_dir=""

while getopts "C:I:O:" opt; do
  case $opt in
    C)
      passing_cells_csv=$OPTARG
      ;;
    I)
      input_dir=$OPTARG
      ;;
    O)
      output_dir=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ -z "$passing_cells_csv" ] || [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 -C <passing_cells_csv> -I <input_dir> -O <output_dir>"
  exit 1
fi

mkdir -p "$output_dir"

awk -F',' 'NR > 1 && $4 > 20000 { gsub(/\+/, "_", $1); print $1 }' "$passing_cells_csv" | while read -r entry; do
  find "$input_dir" -type f -name "*$entry*" | while read -r file_path; do
    filename=$(basename "$file_path")
    ln -s "$file_path" "$output_dir/$filename"
  done
done

echo "Done."
