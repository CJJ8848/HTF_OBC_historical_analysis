#!/bin/bash
set -euo pipefail

# ==================== PATH CONFIG ====================
wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"
#copy from /Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step3_querynewref_distribution/wgfinal_combinedepthandprop_avg/s4_epse2info/HTF_bykmer_final_filtered_corrected64.tsv
#HTF_FILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/HTF_haplotypes/form57_h36_HTF_bykmer_tmp_addthe2modernbymapping.txt"
HTF_FILE=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/HTF_haplotypes/form57_h35_HTF_bykmer.txt
tm1=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/h40sample_bestHTF.txt
tm2=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/m57sample_bestHTF.txt
{
  # Print header (from the first file)
  head -n 1 "$tm1"

  # Concatenate both files, remove all header lines and blanks
  tail -n +2 "$tm1"
  tail -n +2 "$tm2"
} | grep -v '^$' > "$HTF_FILE"
BINARY_FILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h36_epsE2rescued_binary_matrix.tsv"
ESP_FASTA_FILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"
PAF_DIR="${wd}/mappings"
CURATED_LENGTH_FILE1="${PAF_DIR}/h15_espE2_length.txt"
CURATED_LENGTH_FILE="${PAF_DIR}/h15_espE2_length_dotname.txt"
OUTFILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/combined_HTF_oantigen_m57_h35.txt"

# ==================== STEP 1: Generate h15 length file ====================
> "$CURATED_LENGTH_FILE1"
for paf in "${PAF_DIR}"/*_mapped.paf; do
    sample=$(basename "$paf" _mapped.paf)

    if [[ "$sample" == "PL0235" ]]; then
        length="NA"
    else
        length=$(awk 'NR==2 {print $7}' "$paf")
    fi

    echo -e "${sample}\t${length}" >> "$CURATED_LENGTH_FILE1"
done

> "$CURATED_LENGTH_FILE"
while read -r line; do
    converted=$(echo "$line" | sed 's/_/./')
    echo "$converted" >> "$CURATED_LENGTH_FILE"
done < "$CURATED_LENGTH_FILE1"

# ==================== STEP 2: Combine Gene Table ====================
genes=(WfgD rmlC_1 tagG_1 tagH_1 spsA m57_h36_espE2_rescued)
mkdir -p tmp_gene_matrix

# Split matrix by gene
for gene in "${genes[@]}"; do
    awk -v g="$gene" -F"\t" 'NR==1{for(i=2;i<=NF;i++) header[i]=$i} $1==g {for(i=2;i<=NF;i++) print header[i]"\t"$i}' "$BINARY_FILE" > "tmp_gene_matrix/${gene}.tsv"
done

# Generate espE2 info from fasta
awk '
    /^>/ {
        if (seq != "") {
            print id "\t1\t" length(seq);
        }
        id = substr($1, 2);
        seq = "";
        next;
    }
    {
        seq = seq $0;
    }
    END {
        if (seq != "") {
            print id "\t1\t" length(seq);
        }
    }
' "$ESP_FASTA_FILE" > espE2_info.tsv

# Load curated lengths
declare -A curated_lengths
while read -r s len; do
    curated_lengths["$s"]=$len
done < "$CURATED_LENGTH_FILE"

# Write header
echo -e "sample\tHTF_haplotype\tHTFgroup_Oantigen_PA\twfgD_PA\trmlC1_PA\ttagG1_PA\ttagH1_PA\tspsA_PA\tespE2_PA\tespE2length\tnote" > "$OUTFILE"

# Combine table
tail -n +2 "$HTF_FILE" | while IFS=$'\t' read -r sample htf _; do
    obc="NA"
    [[ "$htf" == *"1830"* || "$htf" == *"1383"* ]] && obc="OBC_absent"
    [[ "$htf" == *"1803"* || "$htf" == *"1245"* ]] && obc="OBC_present"
    row="$sample\t$htf\t$obc"

    for gene in "${genes[@]}"; do
        val=$(awk -v s="$sample" '$1 == s {print $2}' "tmp_gene_matrix/${gene}.tsv")
        row+="\t${val:-0}"
    done

    # Determine espE2 presence
    # Determine espE2 presence
    espE2_PA_val=$(echo -e "$row" | awk -F'\t' '{print $(NF)}')

    # Default to 0
    curated_len="0"

    # Try to get curated length from mapping-based results
    if [[ ${curated_lengths[$sample]+_} ]]; then
        curated_len="${curated_lengths[$sample]}"
    else
        # Try softmatch as fallback
        esp_line=$(awk -v s="$sample" '$1 ~ s' espE2_info.tsv | head -n1)
        if [[ -n "$esp_line" ]]; then
            curated_len=$(echo "$esp_line" | cut -f3)
        fi
    fi
    [[ "$espE2_PA_val" == "NA"  ]] && curated_len="NA"

    row+="\t${curated_len:-NA}\t"

    echo -e "$row"
done >> "$OUTFILE"
sed -i 's/^64\.GBR_1933b_S36\s\+HTF_p21\.F9\s*(1245)/64.GBR_1933b_S36\tHTF_p7.G11 (1830)/' "$OUTFILE"
#change the 64.GBR_1933b_S36	HTF_p21.F9 (1245)	OBC_absent	0	0	0	0	0	0	0	
#to  HTF_p7.G11 (1830)

echo "✅ Final combined table with espE2 lengths written to: $OUTFILE"
