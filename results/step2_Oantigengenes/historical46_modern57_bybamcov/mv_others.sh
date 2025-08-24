#!/bin/bash
set -euo pipefail

# File with sample names (without extensions)
keep_list="m55_h40.txt"
dest_dir="others"

# Make directory if not exists
mkdir -p "$dest_dir"

# Read sample names to keep into array
mapfile -t keep_samples < "$keep_list"

# Make an associative array for fast lookup
declare -A keep_map
for s in "${keep_samples[@]}"; do
    keep_map["${s}_gene_coverage.tsv"]=1
done

# Loop through all *_gene_coverage.tsv files
for file in *_gene_coverage.tsv; do
    [[ -e "$file" ]] || continue  # skip if no matching files
    if [[ -z "${keep_map[$file]+_}" ]]; then
        echo "Moving $file to $dest_dir/"
        mv "$file" "$dest_dir/"
    fi
done

echo "✅ Done. All unmatched files moved to $dest_dir/"
