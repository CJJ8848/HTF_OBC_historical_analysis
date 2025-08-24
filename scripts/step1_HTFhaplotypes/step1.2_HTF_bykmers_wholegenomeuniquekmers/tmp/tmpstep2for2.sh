#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=5:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
#$ -N keykmerstep2
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail
# --- Configuration ---fff
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
cd $wd

KMER=31
JFSIZE=300000000

h40_Pv_fastq=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/Pv_only_fastq
m57_raw_fastq=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57
#this is the fastq that coverted from bam files mapped to OTU5 ref with HTF and TFA... for h40
#for m57 we use the pv raw fastq
mkdir -p  step2_isolate_kmers step2_isolate_kmers/isolate_jf step2_isolate_kmers/isolate_jf/dump step3_query step3_query/match_results
dump_out=$wd/step2_isolate_kmers/isolate_jf/dump
# Define your target samples
for iso_name in p1.G2 p23.B2; do
    r1="$m57_raw_fastq/${iso_name}.pair1.truncated.gz"
    r2="$m57_raw_fastq/${iso_name}.pair2.truncated.gz"

    # Check if both R1 and R2 exist
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Missing R1 or R2 for $iso_name, skipping..."
        continue
    fi

    echo "Processing $iso_name..."
    #zcat "$r1" "$r2" | jellyfish count \
    #-m $KMER -s $JFSIZE -C \
    #-o $wd/step2_isolate_kmers/isolate_jf/${iso_name}.jf /dev/fd/0
    jellyfish dump -c "$wd/step2_isolate_kmers/isolate_jf/${iso_name}.jf" > "$dump_out/${iso_name}_dump.txt"


done

