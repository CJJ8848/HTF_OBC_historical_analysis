## HTF k-mer identification and additional analysis pipeline

---

### `step0_wholegenome_7refs/`

- `step1_generatepaf.sh`  
- `step2_summarypaf_corrected_byhand.sh`  
- `step3_searchanddelete.sh`  

> In **step0**, I run `minimap2` to find the coordinates for the HTF region in each reference genome, then delete them to generate 7 whole-genome backgrounds that exclude the self HTF region for downstream k-mer calling.

---

### `step1.1_refkmers_all_to_wgunique.sh`

> In **step1.1**, I run the k-mer calling to generate 7 sets of kmers that are unique from each other and also exclusive from the 7 representative OTU5 whole genome backgrounds.

---

### `step1.2_refkmers_hamming_filtering/`

> In **step1.2**, I compute the hamming distance matrix and visualize it, then use R to iteratively filter out any k-mer pairs with distance ≤1.  
> This yields the final `kmers_unique_hamming2` set.

---

### `step2.1_generateisolate_jf_kmersm57.sh`  
### `step2.2_generateisolate_jf_kmersh40.sh`

> Here I generate `.jf` files for isolate k-mers.  
> For H40 samples, I remove duplicates and use the `rmdup` FASTQ to generate `.jf` files.

---

### `step3_wg_background_distribution/`

- `readme`  
- `step1_dump_wgdistribution_Rvisual.R`  
- `step2_calculate_mean_wg_kmercount_asthreshold.sh`  

> In **step3**, I analyze whole-genome (WG) k-mer distributions.  
> First, I visualize the distribution in R, then compute the mean WG k-mer count as a threshold for filtering.

---

### `step4_query_and_filterby_wgthreshold/`

> In **step4**, I query each isolate’s k-mers and filter them using the average WG k-mer threshold.  
> Then I:
> - Summarize the raw k-mer counts of each HTF–isolate pair
> - Use R to calculate the proportion of HTF-specific k-mers to the total matched HTF k-mers for that isolate
> - Use this index to rank the best HTF matches  
> - Generate mixture histograms of HTF and WG k-mer distributions

> To do this, I run `step1.sh` to obtain `wg_p25.c2_depthsummary` (long-format k-mer counts), then run `mix.R` to plot the distribution.

---

### `step5/`

> In **step5**, I:
> - Run coinfection detection
> - Break ambiguous isolates
> - Remove `64.GBR` as a uncertain calls  
> - Confirm `PL0240` as a confident break
> - Calculate HTF length group frequency differences
> - Plot boxplots and raw data versus time and geography

---

### Final Combine

#### `combine/step1/`

> In `combine`, I merge the HTF results with O-antigen data.  
> This is the **final main figure**.

#### `../step2/`

> In `step2`, I combine **local assembly** and **k-mer** results to produce a **confidence matrix of HTF calls**.

---

**Done**
