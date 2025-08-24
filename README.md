
HTF and O-antigen Profiling of ***Pseudomonas viridiflava*** ATUE5 in historical ***Arabidopsis thaliana*** Metagenomes
==================================================================================

This repository contains scripts, data, and results for profiling the hypothetical tail fiber (HTF) haplotypes and O-antigen biosynthesis gene content of ***P. viridiflava*** isolates from historical herbarium and modern ***Arabidopsis thaliana*** samples.

Analysis pipeline
--------------
HTF Haplotype Assignment:

- [Local Assembly Approach](https://github.com/CJJ8848/HTF_OBC_historical_analysis/blob/db46ec374808d77ba3394c4c6c6bd663b6b60972/scripts/step1_HTFhaplotypes/Readme_step1.1_HTF_bylocalassembly/step1.1_HTF_bylocalassembly.md):
    - Extract reads mapping to HTF/TFA regions
    - Assemble with SPAdes (skip for modern samples, since we have modern assemblies)
    - Assign best haplotype based on covered proportion (minimap2)
- [K-mer Based Approach](https://github.com/CJJ8848/HTF_OBC_historical_analysis/blob/b163377c1c2f099000a1be1d23b7e0a641bad992/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/Readme_step1.2_HTF_bykmers_wholegenomeuniquekmers.md):
    - Build whole-genome-exclusive HTF-unique kmers
    - Apply iterative Hamming ≥ 2 filtering across haplotypes
    - Query isolate .jf files (Jellyfish)
    - Assign dominant haplotype by max (HTF matched kmers/total matched kmers) proportion
    - Detect coinfections if multiple haplotypes exceed threshold
    - HTF length group frequency distribution

[O-antigen Gene P/A Detection](https://github.com/CJJ8848/HTF_OBC_historical_analysis/blob/b163377c1c2f099000a1be1d23b7e0a641bad992/scripts/Readme_step2_Oantigengenes/step2_Oantigengenes.md):

- A gene is considered present if:
    (i) coverage ≥ 50%
    (ii) mean depth ≥ 75% of genome-wide average
- espE2 handled separately via extended mapping and contig rescue

[Combined Analysis](https://github.com/CJJ8848/HTF_OBC_historical_analysis/blob/b163377c1c2f099000a1be1d23b7e0a641bad992/scripts/step3_combine/Readme_step3_combine_generateasummarytable_HTF_Oantigen.md):

- Merge HTF and OBC profiles
- Output metadata tables and combined heatmaps

Large data including raw fastq, fasta, reference and so on are stored on Zenodo. Link: 

Directory Structure
-------------------
data/

    - modern57.txt, 40h_OTU5.txt: lists of modern/historical samples.
    - readme.txt: description of input formats.

scripts/

    - step0_h40_fastq/: extract ATUE5-mapped reads from historical plant herbarium metagenomes and identify tailocin regions in ATUE5 genomes
    - step1_HTFhaplotypes/: HTF haplotype detection
        - step1.1_HTF_bylocalassembly/: local assembly and mapping-based assignment
        - step1.2_HTF_bykmers_wholegenomeuniquekmers/: kmer filtering, querying, and assignment
    - step2_Oantigengenes/: detection of six O-antigen genes P/A including espE2 rescue
    - step3_combine/: integration of HTF and O-antigen data, generation of summary and plots

results/

    - figures_tables/: Contains all figures and tables used in the manuscript and all analyses
    - step0_h40_fastq/: summary table of tailocin presence in ATUE5 genomes
    - step1_HTFhaplotypes/: HTF haplotype results (assembly/kmer-based)
    - step2_Oantigengenes/: binary P/A gene matrix and espE2 analysis
    - step3_combine/: combined tables and plots (for manuscript figures)

