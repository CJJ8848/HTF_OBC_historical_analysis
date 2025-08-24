"HTF–O–Antigen–Gene Presence Summary and espE2 Integration"

## Step 1: Input Files

### (a) OBC-Gene Binary Matrix (including espE2 and espE2 length)

**File path:**

./results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h36_epsE2rescued_binary_matrix.tsv

notes for espE2 length calculation:
./results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt

- If a sequence is present → espE2_PA = 1  
- If absent → espE2_PA = 0  
- Length is calculated from the FASTA entry  
- Column note is added manually

### (b) HTF Haplotypes Table

**File path:**

./results/step3_combine/HTF_haplotypes/form57_h35_HTF_bykmer.txt

---

## Step 2: Derive O-Antigen Classification from HTF

HTF Haplotypes → OBC Classification  
*1830, *1383 → OBC_absent  
*1803, *1245 → OBC_present

Column name: HTFgroup_Oantigen_PA

---

## Step 4: Output Summary Table

Combined columns:  
sample, HTFgroup_Oantigen_PA, wfgD_PA, rmlC1_PA, tagG1_PA, tagH1_PA, spsA_PA, espE2_PA, espE2length, note

**Output file:**

./results/step3_combine/combined_HTF_oantigen_m57_h35.txt

---
