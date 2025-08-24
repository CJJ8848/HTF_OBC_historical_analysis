## historical46_modern57_bybamcov/

From **step1 to step7**, we assessed the presence/absence of six tailocin-associated genes across 46 historical and 57 modern isolates based on BAM coverage.

- In **step1** and **step2**, we identified the coordinates of the six genes within the *Pseudomonas viridiflava* p25.C2 reference genome (which includes multiple contigs).
- In **steps 3–6**, we calculated:
  - **Covered proportion** for each gene
  - **Relative read depth** of each isolate’s mapped library over the gene region

> The BAM files were previously generated using `bwa aln`. Example commands are provided in `step0/step1.sh`.

- In **step7**, we filtered results using the following thresholds:
  - (i) **Coverage ≥ 50%**
  - (ii) **Mean depth ≥ 75% of the isolate’s genome-wide average**

Genes not meeting both thresholds were scored as **absent**. The final result is a **binary gene presence/absence matrix** across isolates.

---

## checkespE2_contigs/

The **espE2** locus required additional manual curation due to its size (~4.5 kb) and high sequence divergence.

- For **modern isolates**:
  - Reads mapping to the p25.C2 *espE2* reference were extracted.
  - If a **fully mapped variant** was found, it was retained.
  - For **partially mapped** hits, the region was extended from the mapped anchor to recover full-length gene sequences (~4,527 ± 200 bp) within the FASTA contigs.
  - This yielded a diverse set of modern *espE2* haplotypes.

- For **historical isolates** (which lack assemblies):
  - The curated **modern haplotype set** was used as an expanded reference.
  - Reads mapping to any haplotype were extracted and **locally assembled** with **SPAdes v3.13.0**.
  - Contigs were then aligned back to the haplotype reference with **Minimap2 v2.1** to identify the best-matching allele.

- For **both modern and historical** datasets:
  - *espE2* sequences were **translated into amino acids** using Expasy.
  - **ORFs were identified** via start/stop codon detection.
  - Sequences were **aligned** with **Clustal Omega v1.2.4**.
  - Final multiple sequence alignments and allele length distributions are provided in the **Supplementary Materials**.
