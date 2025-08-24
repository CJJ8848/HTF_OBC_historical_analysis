#!/bin/bash

# Input file
input="allOTU5candidates40_withdatesandlocs_uniq.txt"

# Columns to calculate means for
columns=("At_percent" "At_covered" "At_depth" "Ps_percent" "Ps_covered" "Ps_depth")

# Print header
echo -e "Column\tMean"

# Loop through each column
for col in "${columns[@]}"; do
    # Use awk to extract the column index based on the header
    col_index=$(head -1 "$input" | tr '\t' '\n' | awk -v target="$col" '{if ($0 == target) print NR}')
    
    if [ -z "$col_index" ]; then
        echo -e "$col\tColumn not found"
        continue
    fi

    # Calculate mean using awk (skip header)
    mean=$(awk -v col="$col_index" 'NR > 1 {sum += $col; count++} END {if (count > 0) print sum / count}' "$input")

    echo -e "$col\t$mean"
done