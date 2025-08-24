# === Libraries ===
library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
setwd('/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/suppfig1_confidence_matrix_HTF')
# === 1️⃣ Load data ===
df_full <- fread("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/suppfig1_confidence_matrix_HTF/modern_57_h35_htf_supp_mapping_kmer_with_lengths.tsv")
colnames(df_full)[1] <- "Strain"

# === 2️⃣ Prepare matrix using existing length columns ===
df_full <- df_full %>%
  select(Strain,
         `From bigdataset` = HTF_length_from_bigdataset,
         `From mapping` = mapping_length,
         `From k-mer` = kmer_length)

# === 3️⃣ Define color mapping ===
base_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

# === 🔁 Heatmap plot function ===
plot_htf_heatmap <- function(df_mod, output_pdf) {
  # Reshape to long then wide
  df_long <- df_mod %>%
    pivot_longer(cols = -Strain, names_to = "method", values_to = "length")
  
  df_wide <- df_long %>%
    pivot_wider(names_from = Strain, values_from = length)
  
  # Get observed lengths
  all_lens <- unique(unlist(df_wide[, -1]))
  all_lens <- all_lens[!is.na(all_lens)]
  
  # Assign grey for other lengths
  other_lens <- setdiff(all_lens, names(base_colors))
  if (length(other_lens) > 0) {
    other_names <- paste0("other_", other_lens)
    names(other_names) <- other_lens
    other_colors <- setNames(colorRampPalette(c("#bbbbbb", "#444444"))(length(other_lens)), other_names)
  } else {
    other_colors <- character(0)
  }
  
  all_colors <- c(base_colors, other_colors, "NA" = "white")
  
  # Recode display matrix
  df_display <- df_wide
  for (i in 2:ncol(df_display)) {
    df_display[[i]] <- sapply(df_display[[i]], function(val) {
      if (is.na(val)) {
        return("NA")
      } else if (val %in% names(base_colors)) {
        return(val)
      } else {
        return(paste0("other_", val))
      }
    })
  }
  
  # === Sort strain order by kmer length group + modern first ===
  df_kmer <- df_mod %>%
    mutate(
      is_modern = grepl("^p", Strain),
      kmer_group = ifelse(`From k-mer` %in% names(base_colors),
                          `From k-mer`,
                          paste0("other_", `From k-mer`))
    ) %>%
    arrange(
      factor(kmer_group, levels = c(names(base_colors), sort(setdiff(unique(kmer_group), names(base_colors))))),
      desc(is_modern),
      Strain
    ) %>%
    filter(Strain %in% colnames(df_display))
  
  strain_order <- df_kmer$Strain
  
  # Define color of column names
  label_colors <- ifelse(grepl("^p", strain_order), "black", "#804111")  # black for modern, blue for historical
  names(label_colors) <- strain_order
  
  # Build matrix and plot
  mat <- as.matrix(df_display[, strain_order])
  rownames(mat) <- c("From bigdataset", "From mapping", "From k-mer")
  
  ht <- Heatmap(
    mat,
    name = "HTF length",
    col = all_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = gpar(col = label_colors),
    heatmap_legend_param = list(title = "HTF length"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey"))
    }
  )
  
  pdf(output_pdf, width = length(strain_order) * 0.25 + 2, height = 4)
  draw(ht)
  dev.off()
  
  cat("✅ Saved:", output_pdf, "\n")
}
# === 🟢 1. All samples, modern first inside HTF groups ===
plot_htf_heatmap(df_full, "modern57_h35_htf_heatmap_by_length_all.pdf")

# === 🟢 2. Modern53 only (from file), excluding 64.* ===
modern53_list <- fread("modern53-infig.txt", header = FALSE)[[1]]
df_m53 <- df_full %>% filter(Strain %in% modern53_list)

plot_htf_heatmap(df_m53, "modern53_htf_heatmap_by_length.pdf")

# === 🟢 3. Only historical (not starting with p, not 64.*) ===
df_hist_only <- df_full %>%
  filter(!grepl("^p", Strain) & !grepl("^64\\.", Strain))

plot_htf_heatmap(df_hist_only, "hist34_htf_heatmap_by_length.pdf")
