---
title: "HTF–O–Antigen–Gene Presence Summary and espE2 Integration"
output: html_document
---

## Step 1: Input Files

### (a) 5-Gene Binary Matrix

**File path:**

/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_binary_matrix.tsv

- Matrix format:
  - Rows: genes
  - Columns: sample names
  - Missing samples → NA

### (b) HTF Haplotypes Table

**File path:**

/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/form57_h37_HTF_bykmer_tmp_addthe2modernbymapping.txt

**Source (copied from):**

/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step3_querynewref_distribution/wgfinal_combinedepthandprop_avg/s4_epse2info/HTF_bykmer_final_filtered_corrected64.tsv

---

## Step 2: Derive O-Antigen Classification from HTF

HTF Haplotypes → OBC Classification  
*1830, *1383 → OBC_absent  
*1803, *1245 → OBC_present

Column name: HTFgroup_Oantigen_PA

---

## Step 3: espE2 Presence and Length

**File:**

/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt

- If a sequence is present → espE2_PA = 1  
- If absent → espE2_PA = 0  
- Length is calculated from the FASTA entry  
- Column note is added manually

---

## Step 4: Output Summary Table

Combined columns:  
sample, HTFgroup_Oantigen_PA, wfgD_PA, rmlC1_PA, tagG1_PA, tagH1_PA, spsA_PA, espE2_PA, espE2length, note

**Output file:**

/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/combined_HTF_oantigen_m57_h37_needOantigenfivegenePAin2modern.txt

---

## Step 5: Append espE2 to Binary Matrix

Start from:  
final_binary_matrix.tsv

Add row espE2, where:  
- 1 = present in the 39 from completeopenframe.18empty_39fasta.txt  
- 0 = all others

**Output path:**

/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/final_sixgenes_binary_matrix.tsv

---

## Reference Files

Directory:  
/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/

Includes:  
- 39fasta_emboss_transeq-I20250724-161752-0560-416209-p1m.out.txt  
- completeopenframe.18empty_39fasta.txt  
- summarytable  
- readme

---

## Working Directory

```r
WORKDIR <- ""
