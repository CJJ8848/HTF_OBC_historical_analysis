## HTF Haplotype Recovery by Local Assembly

### Historical Samples

For historical samples, we extracted reads mapped to the reference tailocin region containing seven HTF haplotypes. Assemblies were performed using these reads, and contigs were aligned to the reference using Minimap2.

- **Step 1:** Successfully assembled tailocin regions for **30 out of 40 historical samples** using locally mapped reads. 

- **Step 2:** Identified HTF haplotypes in **24 historical assemblies**.

- **Step 3:** For MSA, after applying a **65% coverage threshold** to the best-matching reference haplotype, **10 historical HTFs** were retained.

- **Step 4:** One of the nucleotide sequences (`PL0046`) was found to be incomplete after translation to amino acids and was excluded.  
  **Final count:  24 historical HTFs retrieved (9/24 high-confidence)**

Later, HTF identity was further confirmed via k-mer presence analysis, which recovered **35 historical HTFs** (including some not retained here due to quality).

---

### Modern Samples (n = 57)

For modern isolates, we aligned local assemblies to the HTF reference using Minimap2.

- For each sample, the **top-matching HTF and TFA reference** was selected based on **coverage proportion**.
- Most samples showed **>98% coverage** to their top HTF reference, along with presence of the matching TFA region from the same haplotype.
- Occasionally, the best-matching TFA did not correspond to the same HTF as the top hit.
- Final haplotype calls were validated using the **HTF k-mer confidence matrix**.
