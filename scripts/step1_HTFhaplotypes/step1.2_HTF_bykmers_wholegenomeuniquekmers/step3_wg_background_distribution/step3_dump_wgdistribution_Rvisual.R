library(data.table)
library(ggplot2)
library(scales)

# === Paths ===
# === Paths ===
#dump_dir <- "/Users/jiajuncui/Desktop/others/tmphernan/2025_June/step9_coinfection/dump/tmp1"
dump_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
out_dir <- file.path(dirname(dump_dir), "distribution")
out_dir_filtered <- file.path(out_dir, "filtered")
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_filtered, showWarnings = FALSE)

# === List all *_dump.txt files ===
dump_files <- list.files(dump_dir, pattern = "dump.txt$", full.names = TRUE)

# Assuming dump_files is a list of full paths
for (file in dump_files) {
  sample <- sub("_dump.txt$", "", basename(file))
  df <- fread(file, header = FALSE, col.names = c("kmer", "count"))
  
  if (nrow(df) == 0) next
  
  # Frequency table of depth
  freq_df <- as.data.table(table(df$count))
  setnames(freq_df, c("count", "frequency"))
  freq_df[, count := as.numeric(count)]
  freq_df <- freq_df[count > 0]
  
  if (nrow(freq_df) == 0) next

  # === Compute mean and SD ===
  depth_values <- rep(freq_df$count, freq_df$frequency)
  depth_mean <- mean(depth_values)
  depth_sd <- sd(depth_values)
  threshold <- depth_mean + 2 * depth_sd
  
  # === Filtered distribution: count >1 and ≤ mean+2SD ===
  freq_df_filt <- freq_df[count > 1 & count <= threshold]
  if (nrow(freq_df_filt) == 0) next
  
  p_filt <- ggplot(freq_df_filt, aes(x = count, y = frequency)) +
    geom_col(fill = "#33a02c") +
    labs(
      title = paste("Filtered (count >1, ≤ mean+2sd):", sample),
      x = "K-mer Depth",
      y = "Number of K-mers"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = file.path(out_dir_filtered, paste0(sample, "_filtered_kmer_depth_dist.pdf")),
         plot = p_filt, width = 7, height = 5)
}
cat("✅ Done: Plots saved to:\n  -", out_dir, "\n  -", out_dir_filtered, "\n") 
