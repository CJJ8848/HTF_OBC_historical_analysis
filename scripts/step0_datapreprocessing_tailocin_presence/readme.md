# Data Preprocessing: Historical *Pseudomonas viridiflava* ATUE5 Genomes

This pipeline describes the steps used to extract, authenticate, and analyze historical *Pseudomonas viridiflava* ATUE5 genomes from *Arabidopsis thaliana* herbarium-derived metagenomes collected between 1817 and 2015.

---

## [Step 1: Extract ATUE5-Mapped Reads](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats)

A total of 413 globally distributed *A. thaliana* herbarium specimens were screened for the presence of *P. viridiflava* ATUE5.

1. **Host depletion**:
   - Raw metagenomic reads were mapped to the *A. thaliana* TAIR10 reference genome using BWA aln.
   - Unmapped reads were retained.

2. **ATUE5 mapping**:
   - Host-depleted reads were mapped to the *P. viridiflava* ATUE5 reference genome (p25.C2).
   - Samples with ≥65% genome breadth at sufficient depth were retained.
   - **Output**: 46 high-quality historical *P. viridiflava* genomes.

---

## [Step 2: Ancient DNA Authentication](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/scripts/step0_datapreprocessing_tailocin_presence/step2_mapdamage)

Ancient DNA authenticity was assessed using mapDamage:

- Quantified 5′ C→T and 3′ G→A substitution frequencies.
- Verified fragment length distributions.
- All 46 samples showed characteristic ancient DNA patterns. [See DNA damage results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot).  

---

## [Step 3: Historical ATUE5 phlogeny reconstruction](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction)

To confirm strain identity:

- A maximum likelihood phylogeny was built using 296 biallelic SNPs from the 46 genomes, alongside 55 modern ATUE5 and 30 non-ATUE5 reference genomes.
- Tree inference was performed with IQ-TREE using the TPM3+ASC+R2 substitution model.
- Samples clustering within the ATUE5 clade were classified as ATUE5. [See phylogeny results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/results/step0_datapreprocessing_tailocin_presence/step3_h46_m85_tree).
- **Output**: [40 authenticated historical ATUE5 genomes](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/h40_metadata.txt).  

---

## [Step 4: Tailocin Region Detection](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/scripts/step0_datapreprocessing_tailocin_presence/step4_tailocin_present_stats.sh)

We assessed the presence of the tailocin gene cluster:

- All 40 authenticated genomes showed strong coverage of the tailocin region. [See all details](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1eefb4421e942ba509d83aa258371c6b0b1ab0f7/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence/tailocin_coverage_summary.tsv).
  - **Average covered proportion**: 0.81  
  - **Average depth**: 28.57×  

