#!/bin/bash

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract"
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"
trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq"

# Create output directories

#but since SPAdes --isolate requires uniform coverage, the length of regions should be as similar as possible so we use tailocin genes:
formatted_regions="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/others/other_tailocininfo/formatted_regions.txt"
regions=()
while IFS= read -r line; do
  regions+=("$line")
done < "$formatted_regions"



for bam in "$bam_dir"/*.bam; do
  # Extract sample name from BAM file name
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.bam)
  
  # Initialize read IDs file for the sample
  read_ids="${tmp_dir}/${samplename}_read_ids.txt"
  > "$read_ids"
  
  # Check references in the BAM file
  if [[ "$samplename" == *.* ]]; then
      
    raw_fastq1="${fastq_dir}/${samplename}.R1.fastq.gz"
    raw_fastq2="${fastq_dir}/${samplename}.R2.fastq.gz"

    # Find the raw FASTQ file using the determined pattern
    subset_fastq1="${trimmed_fastq_dir}/${samplename}_subsetR1.fastq.gz"
    subset_fastq2="${trimmed_fastq_dir}/${samplename}_subsetR2.fastq.gz"
  
    # Run seqtk subseq and handle errors
    /SAN/ugi/plant_genom/jiajucui/tools/seqtk/seqtk subseq "$raw_fastq1" "${readids_dir}/${samplename}_all_read_ids.txt" | gzip > "$subset_fastq1"
    /SAN/ugi/plant_genom/jiajucui/tools/seqtk/seqtk subseq "$raw_fastq2" "${readids_dir}/${samplename}_all_read_ids.txt" | gzip > "$subset_fastq2"
    echo "finish $samplename" 
    else 
    echo 'allgood'
    fi
#  if [[ "$samplename" != *.* ]]; then
#  zcat $subset_fastq | sed 's/1:N/2:N/g'| gzip > ${trimmed_fastq_dir}/${samplename}_subset.fastq.R2.gz
#  zcat $subset_fastq ${trimmed_fastq_dir}/${samplename}_subset.fastq.R2.gz | gzip > ${trimmed_fastq_dir}/double_except143/${samplename}_subset.fastq.gz

#  #mkdir -p double_except143
#  #zless double_except143/HB0863_subset.fastq.gz | grep ':N' | wc -l
#  else
#  cp $subset_fastq  ${trimmed_fastq_dir}/double_except143/
#  echo 'allgood'
#  fi

  # Remove intermediate files
done

echo "Processing complete. Check the directories for results."
