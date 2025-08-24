#!/bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=24:30:0
#$ -S /bin/bash
#$ -N vcftofastaandclustalotree
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs


# Remember first PL0042 cant use --meta, check readme
# and PL0066 and PL0222 have zero contig, need to be excluded in MSA

# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract"
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"

trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq/"
assembly_dir="${output_dir}/assemblies"
contig_stats_dir="${output_dir}/contig_stats"
mapping_dir="${output_dir}/mappings"
msa_dir="${output_dir}/msa"
new_tree_dir="${output_dir}/tree"
#reference_genome="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
tailocin_region="${output_dir}/regions.txt"
#22:15502-33558
#001.p11.B8.fa
#001.p12.G4.fa
#001.p12.D2.fa
#001.p21.B5.fa
#001.p12.C5.fa

tailocin_fasta="${output_dir}/tailocin_region.fa"
nocontigs_file="${msa_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"
#  Extract the Tailocin region from the reference genome
#samtools faidx "$reference_genome"  -r "$tailocin_region" > "$tailocin_fasta"
#need to rename from 22:15502-33558 to tailocin
#sed -i '1s/.*/>tailocin/' "$tailocin_fasta"
# Create output directories
# mkdir -p $output_dir $tmp_dir $readids_dir $trimmed_fastq_dir 
# rm -r msa/ mappings/ assemblies/ contig_stats/ tree/
#rm -r $assembly_dir $contig_stats_dir $mapping_dir $msa_dir $tree_dir
rm -r $mapping_dir $msa_dir $new_tree_dir 

mkdir -p $mapping_dir $msa_dir $new_tree_dir
# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
> "$contigmapping_file" 
# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/tools
#SPAdes
for bam in "$bam_dir"/*.bam; do
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.bam)
  subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
  # for 143 samples_cat_merged use --12 
  # --12 <file_name> File with interlaced forward and reverse paired-end reads.
  # --merged <file_name> File with merged paired reads. If the properties of the library permit, overlapping paired-end reads can be merged using special software.
  # Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) must be provided for the same library for SPAdes to correctly detect the original read length.
  # but no additional files, --12 is ok
  # --isolate - isolate (standard) bacterial data;
  # --meta The --meta mode is designed for metagenomic data, which typically involves a heterogeneous mixture of reads. also good for dealing with a mixture of read qualities
  # --only-error-correction Performs read error correction only.
  # --only-assembler Runs assembly module only. if you have high quality reads, but here we have a mix of qualities, better run error correction
  # checked --meta gave more info in low quality reads like in120, it was nothing but here 30 contigs 
  # Determine the correct raw FASTQ file name pattern
  # if [[ "$samplename" == *.* ]]; then
  #   spades.py --12 "$subset_fastq" --meta  -o "$assembly_dir/$samplename"
  # else
  # single end has not meta just use multicell/isolate as default
  #   spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"
  # fi
  # checked, if the fastq only contains 1:N, it also could be processed with --12, in which case no even coverage is required and --meta could be used and final we doubled the fastq, so everything ok
  #but the PL0042 process is frozen even using -t 2 and require 50Gb, so remove --meta for it, and it works
  
  # Step 2: Check the number of contigs in each assembly
  contig_file="$assembly_dir/$samplename/contigs.fasta"
#  contig_count=$(grep -c "^>" "$contig_file")
#  echo "$samplename: $contig_count contigs" >> "$contig_stats_dir/contig_counts.txt"
  # Check if the contig file exists
  if [[ ! -f "$contig_file" ]]; then
    echo "$samplename" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi

  reference_genome="${output_dir}/tailocin_region.fa"


  # Define your output files
  paf_file="${mapping_dir}/${samplename}_mapped.paf"
  fasta_out="${msa_dir}/${samplename}_tailocin_region.fasta"

  #  Run minimap2 to get the mapping
  $tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"

  # Define the order and length of each segment
  declare -A segment_lengths
  segment_lengths["tailocin"]=18057
  segment_lengths["001.p11.B8.fa"]=498
  segment_lengths["001.p12.G4.fa"]=513
  segment_lengths["001.p12.D2.fa"]=513
  segment_lengths["001.p21.B5.fa"]=513
  segment_lengths["001.p12.C5.fa"]=513

  segment_names=("tailocin" "001.p11.B8.fa" "001.p12.G4.fa" "001.p12.D2.fa" "001.p21.B5.fa" "001.p12.C5.fa")

  # Calculate cumulative start positions
  declare -A cumulative_starts
  cumulative_starts["tailocin"]=0
  for i in "${!segment_names[@]}"; do
    if [[ $i -gt 0 ]]; then
        prev_segment="${segment_names[$((i-1))]}"
        cumulative_starts["${segment_names[$i]}"]=$((cumulative_starts["$prev_segment"] + segment_lengths["$prev_segment"]))
    fi
  done

  # Initialize the final sequence with 'N's for the total length
  total_length=0
  for length in "${segment_lengths[@]}"; do
    total_length=$((total_length + length))
  done
  final_sequence=$(printf 'N%.0s' $(seq 1 $total_length))

  # Function to replace part of the sequence with the contig sequence
  replace_sequence() {
    local start=$1
    local seq=$2
    final_sequence="${final_sequence:0:start}${seq}${final_sequence:$(($start + ${#seq}))}"
  }
  #count:
  tailocin_count=0
  haplotype_count=0
  hname=''
  # Read the PAF file and fill in the sequence for each segment
  while read -r line; do
    contig_name=$(echo "$line" | awk '{print $1}')
    contig_start=$(echo "$line" | awk '{print $3}')
    contig_end=$(echo "$line" | awk '{print $4}')
    ref_name=$(echo "$line" | awk '{print $6}')
    ref_start=$(echo "$line" | awk '{print $8}')
    ref_end=$(echo "$line" | awk '{print $9}')

    # Extract the contig sequence from the contig file
    contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | tail -n +2 | tr -d '\n')
    
    # Determine the start position in the final sequence
    segment_start=$((ref_start + cumulative_starts["$ref_name"]))
    
    # Replace the corresponding part of the final sequence with the contig sequence
    replace_sequence $segment_start "$contig_seq"
    # Update the counters
    #count:
    if [[ "$ref_name" == "tailocin" ]]; then
            tailocin_count=$((tailocin_count + 1))
        else
	    hname=$ref_name:$ref_start-$ref_end
            haplotype_count=$((haplotype_count + 1))
    fi
  done < "$paf_file"
 
  # Save the final sequence to a FASTA file
  echo ">${samplename}" > "$fasta_out"
  echo "$final_sequence" >> "$fasta_out"
  
  # Output the mapping summary
  echo "$samplename: tailocin mapped contigs = $tailocin_count, $hname mapped contigs = $haplotype_count" >> "$contigmapping_file"
  if [[ $tailocin_count -eq 0 && $haplotype_count -eq 0 ]]; then
    rm $fasta_out
    echo "$samplename no mapped contig"
  fi
  echo "The final sequence has been saved to $fasta_out"

done


# Concatenate all sequences for MSA
cat "$msa_dir"/*.fasta > "$msa_dir/all_samples_tailocin.fasta"

#try use directly the msa by cat merge
# Step 5: Build a phylogenetic tree
iqtree -s "$msa_dir/all_samples_tailocin.fasta" -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/tailocin_tree_bycatmergefasta"


