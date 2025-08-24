#!/bin/bash
set -euo pipefail

# ================= USER CONFIG =================
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2025_h40_afterrmAT/"
outbase="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"        # Output folder
threads=4
mkdir -p $outbase
# ================= SETUP =================
espE2_ref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"
mkdir -p "$outbase"/{fastq_from_bam,mapped_bam_espE2,fastq_espE2only,tmp}

# ================= Step 1: Build Reference =================
bwa index "$espE2_ref"

# ================= Step 2: Process BAMs =================
for bam in "$bam_dir"/*.bam; do
    sample=$(basename "$bam" _after_removal_mappedAt.bam)
    echo "🔍 Processing sample: $sample"

    # Step 2.1: BAM → FASTQ
    fastq_raw="$outbase/fastq_from_bam/${sample}.fastq"
    samtools fastq "$bam" > "$fastq_raw"

    # Step 2.2: bwa aln
    sai_tmp="$outbase/tmp/${sample}.sai"
    bwa aln -t $threads "$espE2_ref" "$fastq_raw" > "$sai_tmp"

    # Step 2.3: bwa samse + samtools sort
    bam_out="$outbase/mapped_bam_espE2/${sample}.espE2.bam"
    bwa samse "$espE2_ref" "$sai_tmp" "$fastq_raw" | samtools sort -@ $threads -o "$bam_out"
    samtools index "$bam_out"

    # Step 2.4: extract mapped reads only
    fastq_out="$outbase/fastq_espE2only/${sample}_espE2mapped.fastq"
    samtools fastq -F 4 "$bam_out" > "$fastq_out"

    echo "✅ Extracted mapped reads to: $fastq_out"
done

# ================= Final Cleanup =================
echo "🧹 Cleaning up intermediate files..."
rm -rf "$outbase/tmp"
rm -rf "$outbase/fastq_from_bam"

echo "🎉 Done! All espE2-mapped reads are in: $outbase/fastq_espE2only/"


  


