# Data Preprocessing: Historical *Pseudomonas viridiflava* ATUE5 Genomes

This pipeline describes the steps used to extract, authenticate, and analyze historical *Pseudomonas viridiflava* ATUE5 genomes from *Arabidopsis thaliana* herbarium-derived metagenomes collected between 1817 and 2015.

---

## Step 1: Extract ATUE5-Mapped Reads

A total of 413 globally distributed *A. thaliana* herbarium specimens were screened for the presence of *P. viridiflava* ATUE5.

1. **Host depletion**:
   - Raw metagenomic reads were mapped to the *A. thaliana* TAIR10 reference genome using BWA.
   - Unmapped reads were retained.

2. **ATUE5 mapping**:
   - Host-depleted reads were mapped to the *P. viridiflava* ATUE5 reference genome (p25.C2).
   - Samples with ≥65% genome breadth at sufficient depth were retained.
   - **Output**: 46 high-quality historical *P. viridiflava* genomes.

---

## Step 2: Ancient DNA Authentication

Ancient DNA authenticity was assessed using mapDamage:

- Quantified 5′ C→T and 3′ G→A substitution frequencies.
- Verified fragment length distributions.
- All 46 samples showed characteristic ancient DNA patterns.  

---

## Step 3: Historical ATUE5 phlogeny reconstruction

To confirm strain identity:

- A maximum likelihood phylogeny was built using 296 biallelic SNPs from the 46 genomes, alongside 55 modern ATUE5 and 30 non-ATUE5 reference genomes.
- Tree inference was performed with IQ-TREE using the TPM3+ASC+R2 substitution model.
- Samples clustering within the ATUE5 clade were classified as ATUE5.
- **Output**: 40 authenticated historical ATUE5 genomes.  

---

## Step 4: Tailocin Region Detection

We assessed the presence of the tailocin gene cluster:

- All 40 authenticated genomes showed strong coverage of the tailocin region.
  - **Average covered proportion**: 0.81  
  - **Average depth**: 28.57×  

