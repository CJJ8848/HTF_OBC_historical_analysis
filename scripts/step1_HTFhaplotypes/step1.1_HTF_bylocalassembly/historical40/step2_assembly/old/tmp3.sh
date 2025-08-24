  reference_genome=$tailocin_fasta
  paf_file="${mapping_dir}/${samplename}_mapped.paf"
  fasta_out="${msa_dir}/${samplename}_tailocin_region.fasta"
  fasta_out1="${msa_dir}/${samplename}_tailocin_region_allconcatenated.fasta"
  $tools/minimap2/minimap2 -cx asm5 "$contig_file" "$reference_genome" > "$paf_file"


  if [[ $(less "$paf_file" | wc -l) -eq 0 ]]; then
    echo "$samplename" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi
  # Calculate cumulative start positions
  declare -A cumulative_starts
  cumulative_starts["tailocin"]=0
  for i in "${!segment_names[@]}"; do
    if [[ $i -gt 0 ]]; then
      prev_segment="${segment_names[$((i-1))]}"
      cumulative_starts["${segment_names[$i]}"]=$((cumulative_starts["$prev_segment"] + segment_lengths["$prev_segment"]))
    fi
  done

  total_length=0
  for length in "${segment_lengths[@]}"; do
    total_length=$((total_length + length))
  done
  final_sequence=$(printf 'N%.0s' $(seq 1 $total_length))

  replace_sequence() {
    local start=$1
    local seq=$2
    final_sequence="${final_sequence:0:start}${seq}${final_sequence:$(($start + ${#seq}))}"
  }

  tailocin_count=0
  htf_count=0
  taf_count=0
  declare -a htf_list
  declare -a tfa_list
  declare -A nonN_counts

  while read -r line; do
    ref_name=$(echo "$line" | awk '{print $1}')
    ref_start=$(echo "$line" | awk '{print $3}')
    ref_end=$(echo "$line" | awk '{print $4}')
    strand=$(echo "$line" | awk '{print $5}')
    contig_name=$(echo "$line" | awk '{print $6}')
    contig_start=$(echo "$line" | awk '{print $8}')
    contig_end=$(echo "$line" | awk '{print $9}')

    if [[ "$strand" == "-" ]]; then
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | seqtk seq -r - | tail -n +2 | tr -d '\n' )
    else
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | tail -n +2 | tr -d '\n')
    fi

    segment_start=$((ref_start + cumulative_starts["$ref_name"]))
    replace_sequence $segment_start "$contig_seq"

    if [[ "$ref_name" == "tailocin" ]]; then
        tailocin_count=$((tailocin_count + 1))
    elif [[ "$ref_name" == HTF* ]]; then
        htf_count=$((htf_count + 1))
        htfname=${ref_name%%:*}
        htf_list+=("$htfname")
    elif [[ "$ref_name" == TFA* ]]; then
        tfa_count=$((tfa_count + 1))
        tfaname=${ref_name%%:*}
        tfa_list+=("$tfaname")
    fi
  done < "$paf_file"

  # Write the final sequence to fasta_out with individual segment headers
  echo ">${samplename}" > "$fasta_out1"
  echo "$final_sequence" >> "$fasta_out1"
  
  ## Write the final sequence to fasta_out with individual segment headers
  {
    for segment in "${segment_names[@]}"; do
      start=${cumulative_starts[$segment]}
      length=${segment_lengths[$segment]}
      segment_sequence=${final_sequence:$start:$length}
      echo ">$segment"
      echo "$segment_sequence"
    done
  } > "$fasta_out"

 # if [[ $tailocin_count -eq 0 && $htf_count -eq 0 && $tfa_count -eq 0 ]]; then
 #   rm $fasta_out $fasta_out1
 #   echo "$samplename no mapped contig"
 # fi

  # Append to the combined MSA file before formatting individual segments
  cat "$fasta_out1" >> "$msa_dir/all_modern76fa_samples_tailocin.fasta"

  # Output the mapping summary and lists to the contigmapping file
  {
    echo "$samplename: tailocin mapped contigs = $tailocin_count, TFA mapped contigs = $tfa_count, HTF mapped contigs = $htf_count"
    
    echo "HTF List:"
    for htf_segment in "${htf_list[@]}"; do
        echo "$htf_segment"
    done

    echo "TFA List:"
    for tfa_segment in "${tfa_list[@]}"; do
        echo "$tfa_segment"
    done
  } >> "$contigmapping_file"

  # Process segments to calculate non-N counts and proportions
  {
    for segment in "${tfa_list[@]}"; do
      sequence=$(grep -A1 ">${segment}" "$fasta_out" | tail -n1)
      nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
      proportion=$(echo "scale=2; $nonN_count/${segment_lengths[$segment]}" | bc)
      nonN_counts["$segment"]=$proportion
      echo "$samplename >${segment}:$proportion"
    done

    for segment in "${htf_list[@]}"; do
      sequence=$(grep -A1 ">${segment}" "$fasta_out" | tail -n1)
      nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
      proportion=$(echo "scale=2; $nonN_count/${segment_lengths[$segment]}" | bc)
      nonN_counts["$segment"]=$proportion
      echo "$samplename >${segment}:$proportion"
    done
  } >> "$nonN_file"

  longest_tfa=$(printf "%s\n" "${!nonN_counts[@]}" | grep "^TFA" | while read segment; do echo "$segment ${nonN_counts[$segment]}"; done | sort -k2,2nr | head -n1 | awk '{print $1}')
  longest_htf=$(printf "%s\n" "${!nonN_counts[@]}" | grep "^HTF" | while read segment; do echo "$segment ${nonN_counts[$segment]}"; done | sort -k2,2nr | head -n1 | awk '{print $1}')
  # Filter out segments that are still all Ns
  if [[ -n "$longest_tfa" ]]; then
    longest_tfa_seq=$(grep -A1 ">${longest_tfa%%:*}" "$fasta_out" | tail -n1)
    longest_tfa_nonN=$(echo "$longest_tfa_seq" | tr -cd 'ATCGatcg' | wc -c)
    if [[ "$longest_tfa_nonN" -eq 0 ]]; then
      longest_tfa=""
    fi
  fi

  if [[ -n "$longest_htf" ]]; then
    longest_htf_seq=$(grep -A1 ">${longest_htf%%:*}" "$fasta_out" | tail -n1)
    longest_htf_nonN=$(echo "$longest_htf_seq" | tr -cd 'ATCGatcg' | wc -c)
    if [[ "$longest_htf_nonN" -eq 0 ]]; then
      longest_htf=""
    fi
  fi

  # Create the final FASTA file with the selected segments
  final_fasta="${haplotype_dir}/${samplename}.final.fasta"
  {
    echo ">tailocin"
    grep -A1 ">tailocin" "$fasta_out" | tail -n1

    if [[ -n "$longest_tfa" ]]; then
      echo ">${longest_tfa%%:*}"
      grep -A1 ">${longest_tfa%%:*}" "$fasta_out" | tail -n1
    fi

    if [[ -n "$longest_htf" ]]; then
      echo ">${longest_htf%%:*}"
      grep -A1 ">${longest_htf%%:*}" "$fasta_out" | tail -n1
    fi
  } > "$final_fasta"

 # if [[ $tailocin_count -eq 0 && $htf_count -eq 0 && $tfa_count -eq 0 ]]; then
 #   rm $final_fasta 
 #   echo "$samplename no mapped contig"
 # fi

  echo "The final sequence has been saved to $final_fasta"
done

#
haplotype_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85/haplotype_selected"

# Create output files for HTF and TFA
htf_fasta="${haplotype_dir}/all_HTF_samples.fasta"
tfa_fasta="${haplotype_dir}/all_TFA_samples.fasta"
> "$htf_fasta"
> "$tfa_fasta"

# Iterate over each final fasta file and extract HTF and TFA sequences
for final_fa in "${haplotype_dir}"/*.final.fasta; do
  samplename=$(basename "$final_fa" .final.fasta)
  
  # Extract HTF sequence
  grep -A1 ">HTF" "$final_fa" | sed "s/^>/>${samplename}|/" >> "$htf_fasta"
  
  # Extract TFA sequence
  grep -A1 ">TFA" "$final_fa" | sed "s/^>/>${samplename}|/" >> "$tfa_fasta"
done

echo "HTF and TFA multi-sample FASTA files have been generated in $haplotype_dir."
#try use directly the msa by cat merge
# Step 5: Build a phylogenetic tree
iqtree -s "$haplotype_dir/selected_modern_all_samples_tailocin.fasta" -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/tailocin_tree_selected"

