
HTF and O-antigen Profiling of Pseudomonas viridiflava in Arabidopsis Metagenomes
==================================================================================

This repository contains scripts, data, and results for profiling the hypothetical tail fiber (HTF) haplotypes and O-antigen biosynthesis gene content of P. viridiflava isolates from historical herbarium and modern Arabidopsis samples.

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
    - step0_h40_fastq/: summary table of tailocin presence in ATUE5 genomes
    - step1_HTFhaplotypes/: HTF haplotype results (assembly/kmer-based)
    - step2_Oantigengenes/: binary P/A gene matrix and espE2 analysis
    - step3_combine/: combined tables and plots (for manuscript figures)

Analysis pipeline
--------------
HTF Haplotype Assignment:
1. Local Assembly Approach:
    - Extract reads mapping to HTF/TFA regions
    - Assemble with SPAdes (skip for modern samples, since we have modern assemblies)
    - Assign best haplotype based on covered proportion (minimap2)
2. K-mer Based Approach:
    - Build whole-genome-exclusive HTF-unique kmers
    - Apply iterative Hamming ≥ 2 filtering across haplotypes
    - Query isolate .jf files (Jellyfish)
    - Assign dominant haplotype by max (HTF matched kmers/total matched kmers) proportion
    - Detect coinfections if multiple haplotypes exceed threshold
    - HTF length group frequency distribution
O-antigen Gene P/A Detection:
- A gene is considered present if:
    (i) coverage ≥ 50%
    (ii) mean depth ≥ 75% of genome-wide average
- espE2 handled separately via extended mapping and contig rescue

Combined Analysis:
- Merge HTF and OBC profiles
- Output metadata tables and combined heatmaps

Large data including raw fastq, fasta, reference and so on are stored on Zenodo. Link: 

