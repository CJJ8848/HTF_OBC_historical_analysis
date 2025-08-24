#!/bin/bash
set -euo pipefail
#
	#1.	Concatenates the depth and covprop files.
	#2.	Adds suffixes _depth and _covprop to gene rows.
	#3.	Filters out unwanted genes (espE_2, tagG_2, tagH_2).
	#4.	Creates a binary matrix where a gene-sample pair is 1 if:
	#•	depth >= 0.75
	#•	covprop >= 50
#Otherwise, the value is 0.
# Define directories and input/output files
RESULTS_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical46_modern57_bybamcov/coverage_h40_m57"
DEPTH_FILE="${RESULTS_DIR}/relative_depth_matrix.tsv"
COV_FILE="${RESULTS_DIR}/covprop_matrix.tsv"
CAT_FILE="${RESULTS_DIR}/../final_cov_depth_wide_cat.tsv"
OUT_FILE="${RESULTS_DIR}/../final_sixgenes_binary_matrix.tsv"

# Step 1: Create concatenated file with renamed rows
echo " Creating combined file: $CAT_FILE"
head -n1 "$DEPTH_FILE" > "$CAT_FILE"
tail -n +2 "$DEPTH_FILE" | awk -F'\t' -v OFS='\t' '{ $1 = $1 "_depth"; print }' >> "$CAT_FILE"
tail -n +2 "$COV_FILE" | awk -F'\t' -v OFS='\t' '{ $1 = $1 "_covprop"; print }' >> "$CAT_FILE"

# Step 2: Prepare binary output
echo " Creating binary matrix based on depth >= 0.75 and covprop >= 50"
{
    # Read all lines into an array
    mapfile -t all_lines < "$CAT_FILE"

    # Extract header
    header="${all_lines[0]}"
    echo "$header" > "$OUT_FILE"

    # Extract gene list (excluding unwanted ones)
    declare -A depth_lines cov_lines
    for line in "${all_lines[@]:1}"; do
        gene=$(echo "$line" | cut -f1)
        case "$gene" in
            *tagG_2*|*tagH_2*) continue ;; # skip unwanted
        esac
        if [[ "$gene" == *_depth ]]; then
            depth_lines["$gene"]="$line"
        elif [[ "$gene" == *_covprop ]]; then
            cov_lines["$gene"]="$line"
        fi
    done

    # Process binary matrix
    for dkey in "${!depth_lines[@]}"; do
        gbase="${dkey%_depth}"
        ckey="${gbase}_covprop"
        if [[ -z "${cov_lines[$ckey]+set}" ]]; then
            echo "⚠️ Missing covprop for $gbase, skipping" >&2
            continue
        fi

        dvals=(${depth_lines[$dkey]})
        cvals=(${cov_lines[$ckey]})

        # Start new binary line
        binary_line=("${gbase}")

        for ((i = 1; i < ${#dvals[@]}; i++)); do
            d="${dvals[$i]}"
            c="${cvals[$i]}"

            # Check if both values are numeric
            if [[ "$d" =~ ^[0-9.]+$ && "$c" =~ ^[0-9.]+$ ]]; then
                if [[ $(echo "$d >= 0.75" | bc -l) -eq 1 && $(echo "$c >= 50" | bc -l) -eq 1 ]]; then
                    binary_line+=("1")
                else
                    binary_line+=("0")
                fi
            else
                binary_line+=("0")
            fi
        done

        # Write binary line
        echo -e "${binary_line[*]}" | tr ' ' '\t' >> "$OUT_FILE"
    done
}

echo "inary matrix written to: $OUT_FILE"
