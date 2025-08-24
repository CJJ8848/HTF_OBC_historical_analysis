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

# ================= USER CONFIG =================
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"        # Output folder
reference_genome="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"

trimmed_fastq_dir="$output_dir/fastq_espE2only/"

tmp_dir="${output_dir}/tmp"

assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"
nocontigs_file="${mapping_dir}/nocontigs.txt"

#rm -r $assembly_dir $mapping_dir $msa_dir $new_tree_dir $haplotype_dir
#rm -r $mapping_dir $msa_dir $new_tree_dir $vcf_dir

mkdir -p $assembly_dir  $mapping_dir $tmp_dir

# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools
#SPAdes
for subset_fastq in "$trimmed_fastq_dir"/*_espE2mapped.fastq; do
  samplename=$(basename "$subset_fastq" _espE2mapped.fastq)
  # single end has not meta just use multicell/isolate as default
#  spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  --careful -k 21,33
  spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  -k 17,21
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
  $tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"
  #$tools/minimap2/minimap2 -cx asm5 "$contig_file" "$reference_genome" > "$paf_file"


  if [[ $(less "$paf_file" | wc -l) -eq 0 ]]; then
    echo "$samplename nomappedcontig" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi
done  



#then summary like in modern57
