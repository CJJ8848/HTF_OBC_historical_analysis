# Load libraries
library(Biostrings)
library(ggplot2)
library(dplyr)

# Set paths
fasta_file <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_refkmers_hamming_filtering/step1_cal_hamming_pairwise_matrix/Rvisualization/runlocationofkmers_before_after_hamming/7HTF_dna.txt"
kmer_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2/intermediate_kmersfiltering/kmers_unique/"
filtered_kmer_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2"
out_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/Rvisualization/runlocationofkmers"
# Load sequences
htf_seqs <- readDNAStringSet(fasta_file)
names(htf_seqs) <- gsub("\\s.*", "", names(htf_seqs))

# Function to load kmers and find positions
find_kmer_positions <- function(htf_name, kmer_file, status_label) {
  kmers <- readLines(kmer_file)
  seq <- htf_seqs[[htf_name]]
  positions <- sapply(kmers, function(k) {
    pos <- matchPattern(k, seq)
    if (length(pos) > 0) start(pos)[1] else NA
  })
  data.frame(
    kmer = kmers,
    position = positions,
    HTF = htf_name,
    Status = status_label
  ) |> filter(!is.na(position))
}

# Get HTF names
htf_names <- names(htf_seqs)

# Original kmers
df_all <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(kmer_dir, paste0(htf, "_unique.txt"))
  if (file.exists(kmer_file)) {
    find_kmer_positions(htf, kmer_file, status_label = "beforehamming")
  }
}))

# Filtered kmers
df_filtered <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(filtered_kmer_dir, paste0(htf, "_unique_filtered_iterative2.txt"))
  if (file.exists(kmer_file)) {
    find_kmer_positions(htf, kmer_file, status_label = "afterhamming")
  }
}))

# Combine with new labels for y-axis
df_all$HTFstatus <- paste0(df_all$HTF, "_beforehamming")
df_filtered$HTFstatus <- paste0(df_filtered$HTF, "_afterhamming")

df_combined <- bind_rows(df_all, df_filtered)

# Plot combined
p_combined <- ggplot(df_combined, aes(x = position, y = HTFstatus, color = Status)) +
  geom_point(alpha = 0.7, size = 1.1) +
  scale_color_manual(values = c("beforehamming" = "darkblue", "afterhamming" = "darkred")) +
  theme_minimal() +
  labs(
    title = "K-mer Position in HTF References (Before vs After Filtering)",
    x = "Position",
    y = "HTF (Status)",
    color = "K-mer Set"
  )

# Save
ggsave(file.path(out_dir, "HTF_kmer_positions_combined.pdf"), p_combined, width = 10, height = 6)

cat("✅ Combined plot saved to:", out_dir, "\n")
