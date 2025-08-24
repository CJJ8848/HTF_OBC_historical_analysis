# Load libraries
library(dplyr)
library(ggplot2)
library(binom)

# ================================
# 1️⃣ Load modern + historical prop tables
# ================================
# Modern
modern_df <- read.table("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummarymodern57_tailocin_kmer_propnorm.tsv" ,header = TRUE, sep = "\t"
)

# Historical
historical_df <- read.table(
"/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummaryhistorical40_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)

#rm historical_df$Isolate=="64.GBR_1933b_S36" & historical_df$Reference=="HTF_p21.F9"
#for now, we believe 64 is unknown
historical_df <- historical_df[!(historical_df$Isolate == "64.GBR_1933b_S36"), ]

# ================================
# 2️⃣ Add IsolateType & merge
# ================================

modern_df$IsolateType <- "Modern"
historical_df$IsolateType <- "Historical"

df_all <- bind_rows(modern_df, historical_df)

# ================================
# 🔥 Apply exclusion filter
# ================================
base_exclude <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258", "PL0027","PL0053","30.ESP_1983b","PL0026")
excluded_samples <- c(base_exclude, paste0(base_exclude, "|NA"))

df_all <- df_all %>% filter(!(Isolate %in% excluded_samples))

# ================================
# 3️⃣ Annotate length category
# ================================

df_all <- df_all %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5", "HTF_p23.B8", "HTF_p26.D6", "HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245",
    TRUE ~ "Unknown"
  ))

df_all <- df_all %>% filter(LengthGroup != "Unknown")

# ================================
# 4️⃣ Get dominant haplotype per isolate
# ================================

dominant_df <- df_all %>%
  group_by(Isolate) %>%
  filter(PropNorm == max(PropNorm, na.rm = TRUE)) %>%
  ungroup()

# ================================
# 5️⃣ Calculate frequency + CI per LengthGroup + IsolateType
# ================================

# Get total N per IsolateType
total_counts <- dominant_df %>%
  distinct(Isolate, IsolateType) %>%
  count(IsolateType) %>%
  rename(TotalN = n)

# Count each length group per Modern/Historical
freq_table_length <- dominant_df %>%
  group_by(LengthGroup, IsolateType) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  left_join(total_counts, by = "IsolateType") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower,
    Upper = ci$upper
  ) %>%
  select(LengthGroup, IsolateType, Frequency, TotalN, Proportion, Lower, Upper)

print(freq_table_length)

# ================================
# 6️⃣ Define colors
# ================================

length_colors <- c(
  "1830" = "pink",
  "1383" = "#2166AC",
  "1803" = "#D6604D",
  "1245" = "#FFCC00"
)



# ✅ Keep your data processing steps same

# ✅ Here’s the updated plotting part:
p_length <- ggplot(freq_table_length, aes(
  x = LengthGroup,
  y = Proportion,
  color = LengthGroup,
  alpha = IsolateType
)) +
  geom_point(
    size = 6,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.1,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(
    aes(label = paste0(Frequency, "/", TotalN)),
    position = position_dodge(width = 0.5),
    vjust = -1,
    size = 5
  ) +
  scale_color_manual(values = length_colors) +
  scale_alpha_manual(values = c("Modern" = 1, "Historical" = 0.6)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "HTF Length Group Frequencies (Modern vs Historical) (95% CI)",
    x = "HTF Length Group",
    y = "Proportion of Isolates",
    color = "HTF Length Group",
    alpha = "Isolate Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12)
  )

print(p_length)

# Save
ggsave(
  "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step1_htf_4lengthgroup_frequencies_pointplot.pdf",
  plot = p_length,
  width = 8, height = 6
)
