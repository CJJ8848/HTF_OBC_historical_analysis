#!/bin/bash
#$ -l tmem=50G
#$ -l h_vmem=50G
#$ -l h_rt=24:00:0
#$ -S /bin/bash
#$ -N spadeh40
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/logs
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/logs
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/logs



# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40"
#"/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40"
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"

output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical"

tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"

trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq/"
assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"
nocontigs_file="${mapping_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"

#rm -r $assembly_dir $mapping_dir $msa_dir $new_tree_dir $haplotype_dir
#rm -r $mapping_dir $msa_dir $new_tree_dir $vcf_dir

mkdir -p $assembly_dir  $mapping_dir 


# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
> "$contigmapping_file" 
reference_genome="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/espE2_p25c2.fasta"
# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools
#SPAdes
for bam in "$bam_dir"/*.bam; do
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.bam)
  subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
  # single end has not meta just use multicell/isolate as default
  spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  --careful -k 21,33
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



  # Define your output files
  paf_file="${mapping_dir}/${samplename}_mapped.paf"
  #  Run minimap2 to get the mapping
  #sergio's idea: map the reference short chunk of hyplotypes to the long contig, in principle works for both historical and modern samples.
  #$tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"
  $tools/minimap2/minimap2 -cx asm5 "$contig_file" "$reference_genome" > "$paf_file"


  if [[ $(less "$paf_file" | wc -l) -eq 0 ]]; then
    echo "$samplename nomappedcontig" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi
done  



#then summary like in modern57
