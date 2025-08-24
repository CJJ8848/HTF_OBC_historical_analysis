#!/bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=24:30:0
#$ -S /bin/bash
#$ -N vcftofastaandclustalotree
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs



# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract"
mkdir -p $output_dir
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"

trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq/"
assembly_dir="${output_dir}/assemblies"
contig_stats_dir="${output_dir}/contig_stats"
mapping_dir="${output_dir}/mappings"
msa_dir="${output_dir}/msa"
new_tree_dir="${output_dir}/tree"
nocontigs_file="${msa_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"
haplotype_dir="${output_dir}/haplotype_selected"
nonN_file="${haplotype_dir}/nonN_TFAandHTF.txt"

# rm -r msa/ mappings/ assemblies/ contig_stats/ tree/
rm -r $assembly_dir $contig_stats_dir $mapping_dir $msa_dir $new_tree_dir $haplotype_dir
#rm -r $mapping_dir $msa_dir $new_tree_dir $vcf_dir

mkdir -p $assembly_dir $contig_stats_dir $mapping_dir $msa_dir $new_tree_dir $haplotype_dir

tailocin_fasta="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tailocin_region.fa"
reference_genome="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
tailocin_region="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/regions.txt"

#22:15502-33559 but all the 40 genes
#TFA_p23.B8:1-513
#TFA_p26.D6:1-513
#TFA_p21.F9:1-498
#TFA_p25.C2:1-513
#TFA_p5.D5:1-513
#TFA_p25.A12:1-546
#TFA_p7.G11:1-546
#HTF_p23.B8:1-1803
#HTF_p26.D6:1-1803
#HTF_p21.F9:1-1245
#HTF_p25.C2:1-1803
#HTF_p5.D5:1-1803
#HTF_p25.A12:1-1383
#HTF_p7.G11:1-1830
#segment_lengths["TFA_22_fragment_1_21_sp|P03740|TFA_LAMBD_Tail_fiber_assembly_protein"]=$((18230 - 17718 + 1))
#segment_lengths["HTF_22_fragment_1_22_hypothetical_protein"]=$((20043 - 18241 + 1))
while read -r region name; do
  samtools faidx "$reference_genome" "$region" | sed "1s/.*/>$name/" >> "$tailocin_fasta"
done < "$tailocin_region"

# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
> "$contigmapping_file" 
> "$nonN_file"

# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/tools




# Define regions and sequences
declare -A segment_lengths
segment_lengths["tailocin"]=18057
segment_lengths["TFA_p23.B8"]=513
segment_lengths["TFA_p26.D6"]=513
segment_lengths["TFA_p21.F9"]=498
segment_lengths["TFA_p25.C2"]=513
segment_lengths["TFA_p5.D5"]=513
segment_lengths["TFA_p25.A12"]=546
segment_lengths["TFA_p7.G11"]=546
segment_lengths["HTF_p23.B8"]=1803
segment_lengths["HTF_p26.D6"]=1803
segment_lengths["HTF_p21.F9"]=1245
segment_lengths["HTF_p25.C2"]=1803
segment_lengths["HTF_p5.D5"]=1803
segment_lengths["HTF_p25.A12"]=1383
segment_lengths["HTF_p7.G11"]=1830
segment_lengths["TFA_22_fragment_1_21_sp|P03740|TFA_LAMBD_Tail_fiber_assembly_protein"]=$((18230 - 17718 + 1))
segment_lengths["HTF_22_fragment_1_22_hypothetical_protein"]=$((20043 - 18241 + 1))
segment_names=(
    "tailocin"
    "TFA_p23.B8"
    "TFA_p26.D6"
    "TFA_p21.F9"
    "TFA_p25.C2"
    "TFA_p5.D5"
    "TFA_p25.A12"
    "TFA_p7.G11"
    "HTF_p23.B8"
    "HTF_p26.D6"
    "HTF_p21.F9"
    "HTF_p25.C2"
    "HTF_p5.D5"
    "HTF_p25.A12"
    "HTF_p7.G11"
    "TFA_22_fragment_1_21_sp|P03740|TFA_LAMBD_Tail_fiber_assembly_protein"
    "HTF_22_fragment_1_22_hypothetical_protein"
)



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
  #but the PL0042 process is frozen even using -t 2 and require 50Gb, so remove --meta for it, and it works
  if [[ "$samplename" == *.* ]]; then
    subset_fastq1="${trimmed_fastq_dir}/${samplename}_subsetR1.fastq.gz"
    subset_fastq2="${trimmed_fastq_dir}/${samplename}_subsetR2.fastq.gz"
    
    spades.py --merge "$subset_fastq" -1 $subset_fastq1 -2 $subset_fastq2  --careful -o "$assembly_dir/$samplename" -k 21,33
  else
  # single end has not meta just use multicell/isolate as default
    spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  --careful -k 21,33
  fi
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
  #sergio's idea: map the reference short chunk of hyplotypes to the long contig, in principle works for both historical and modern samples.
  #$tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"
  $tools/minimap2/minimap2 -cx asm5 "$contig_file" "$reference_genome" > "$paf_file"
  # Define the order and length of each segment

  #!/bin/bash





#try use directly the msa by cat merge
# Step 5: Build a phylogenetic tree
iqtree -s "$msa_dir/all_samples_tailocin.fasta" -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/tailocin_tree_bycatmergefasta"


