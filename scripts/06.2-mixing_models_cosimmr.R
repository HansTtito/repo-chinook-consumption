# =============================================================================
# STABLE ISOTOPE MIXING MODELS - COSIMMR
# Diet estimation of Chinook salmon using delta-13C, delta-15N, delta-34S
# Tissue: muscle
# Sources: 5 (Anchovy, Sardine, Zooplankton, Fish, Invertebrates)
# Covariate: body mass (Mass_g)
# Selected model: Model 8 (muscle, 5 sources, with covariate)
# =============================================================================

library(tidyverse)
library(cosimmr)
library(ggplot2)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("==============================================================\n")
cat("LOADING DATA\n")
cat("==============================================================\n\n")

raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^Sar") ~ "Sardine",
  str_detect(raw_data$Sample, "^Zoo") ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito") ~ "Phytoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchovy",
  TRUE ~ "Unknown"
)

raw_data$pct_sum <- raw_data$pctN + raw_data$pctC + raw_data$pctS

cat("Dataset dimensions:", nrow(raw_data), "rows x", ncol(raw_data), "cols\n")
print(table(raw_data$group))

sources_freshwater <- read.csv("data/data_raw/diet-isotopes/datos_rios.csv")
cat("\nFreshwater sources loaded:\n")
print(sources_freshwater)

# =============================================================================
# 2. QUALITY CONTROL
# =============================================================================

cat("\n==============================================================\n")
cat("QUALITY CONTROL\n")
cat("==============================================================\n\n")

qc_flags <- raw_data %>%
  mutate(
    flag_pct_sum      = ifelse(pct_sum > 100, "pct_sum>100", ""),
    flag_CN_extreme   = case_when(
      group %in% c("Salmon_muscle") & CtoN > 10 ~ "C/N_extreme",
      group %in% c("Sardine", "Anchovy") & CtoN > 8 ~ "C/N_extreme",
      TRUE ~ ""
    ),
    flag_d34S_missing  = ifelse(is.na(d34S), "d34S_missing", ""),
    flag_impossible    = case_when(
      pctN < 0 | pctN > 100 ~ "N_impossible",
      pctC < 0 | pctC > 100 ~ "C_impossible",
      pctS < 0 | pctS > 100 ~ "S_impossible",
      TRUE ~ ""
    ),
    all_flags = paste(flag_pct_sum, flag_CN_extreme, flag_d34S_missing, flag_impossible, sep = ";"),
    all_flags = str_remove_all(all_flags, "^;+|;+$|;{2,}"),
    all_flags = ifelse(all_flags == "", "OK", all_flags)
  )

problems <- qc_flags %>%
  filter(all_flags != "OK") %>%
  select(ID, Sample, group, pct_sum, CtoN, d34S, all_flags)

if (nrow(problems) > 0) {
  cat("Samples with issues:\n")
  print(problems)
} else {
  cat("No QC issues found\n")
}

exclude_samples <- qc_flags %>%
  filter(str_detect(all_flags, "pct_sum>100")) %>%
  pull(ID)

if (length(exclude_samples) > 0) {
  cat("\nExcluded samples:\n")
  print(qc_flags[qc_flags$ID %in% exclude_samples, c("ID", "Sample", "group", "pct_sum", "all_flags")])
}

clean_data <- qc_flags %>% filter(!ID %in% exclude_samples)

cat("\nDataset after QC:", nrow(clean_data), "samples\n")
print(table(clean_data$group))

# =============================================================================
# 3. LIPID CORRECTION FOR delta-13C
# =============================================================================

cat("\n==============================================================\n")
cat("LIPID CORRECTION\n")
cat("==============================================================\n\n")
cat("  - Salmon muscle: optimized McConnaughey for Chinook (D=6.31, I=0.103; Lerner 2024)\n")
cat("  - Sardine/Anchovy/Zooplankton/Phytoplankton: no correction applied\n\n")

lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term       <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent        <- 93 / (1 + (1 / lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  final_correction <- ifelse(CN_ratio > 1, correction_factor, 0)
  return(d13C + final_correction)
}

clean_data <- clean_data %>%
  mutate(
    d13C_corrected = case_when(
      group == "Salmon_muscle" ~ lipid_correction_optimized_chinook(d13C, CtoN),
      TRUE ~ d13C
    )
  )

correction_summary <- clean_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    d13C_original_mean  = round(mean(d13C, na.rm = TRUE), 2),
    d13C_corrected_mean = round(mean(d13C_corrected, na.rm = TRUE), 2),
    mean_correction     = round(mean(d13C_corrected - d13C, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Lipid correction impact by group:\n")
print(correction_summary)

# =============================================================================
# 4. PREPARE MUSCLE CONSUMERS
# =============================================================================

cat("\n==============================================================\n")
cat("PREPARING MUSCLE CONSUMERS\n")
cat("==============================================================\n\n")

consumers_muscle <- clean_data %>%
  filter(group == "Salmon_muscle") %>%
  select(ID, Sample, age, fork_length, weight, d13C_corrected, d15N, d34S, CtoN) %>%
  filter(!is.na(d13C_corrected) & !is.na(d15N) & !is.na(d34S)) %>%
  rename(Mass_g = weight)

cat("Muscle consumers:", nrow(consumers_muscle), "individuals\n")
cat("Mass range:", min(consumers_muscle$Mass_g, na.rm = TRUE), "-",
    max(consumers_muscle$Mass_g, na.rm = TRUE), "g\n\n")

# =============================================================================
# 5. CALCULATE MARINE SOURCES
# =============================================================================

cat("=== CALCULATING MARINE SOURCES ===\n\n")

sources_marine <- clean_data %>%
  filter(group %in% c("Sardine", "Anchovy", "Zooplankton")) %>%
  group_by(group) %>%
  summarise(
    n          = n(),
    d13C_mean  = mean(d13C_corrected, na.rm = TRUE),
    d13C_sd    = sd(d13C_corrected, na.rm = TRUE),
    d15N_mean  = mean(d15N, na.rm = TRUE),
    d15N_sd    = sd(d15N, na.rm = TRUE),
    d34S_mean  = mean(d34S, na.rm = TRUE),
    d34S_sd    = sd(d34S, na.rm = TRUE),
    habitat    = "marine",
    .groups    = "drop"
  )

cat("Marine sources:\n")
print(sources_marine)

# =============================================================================
# 6. PREPARE FRESHWATER SOURCES
# =============================================================================

cat("\n=== FRESHWATER SOURCES ===\n\n")

sources_freshwater <- sources_freshwater %>%
  rename(group = Source) %>%
  mutate(habitat = "freshwater", n = NA) %>%
  select(group, n, d13C_mean, d13C_sd, d15N_mean, d15N_sd, d34S_mean, d34S_sd, habitat)

cat("Freshwater sources:\n")
print(sources_freshwater)

# =============================================================================
# 7. COMBINE SOURCES (5 dietary sources: Anchovy, Sardine, Zooplankton, Fish, Invertebrates)
# =============================================================================

cat("\n=== COMBINING SOURCES ===\n\n")
cat("Note: Phytoplankton and Primary_producers excluded — not direct prey of Chinook salmon\n\n")

sources_5 <- bind_rows(sources_marine, sources_freshwater) %>%
  filter(!group %in% c("Phytoplankton", "Primary_producers"))

cat("Final source list (5 sources):\n")
print(sources_5)

source_names_5 <- sources_5$group

source_means_5 <- as.matrix(sources_5[, c("d13C_mean", "d15N_mean", "d34S_mean")])
rownames(source_means_5) <- source_names_5

source_sds_5 <- as.matrix(sources_5[, c("d13C_sd", "d15N_sd", "d34S_sd")])
rownames(source_sds_5) <- source_names_5

# =============================================================================
# 8. TROPHIC ENRICHMENT FACTORS (TEFs)
# =============================================================================

cat("\n==============================================================\n")
cat("TROPHIC ENRICHMENT FACTORS (TEFs)\n")
cat("==============================================================\n\n")

TEF_means <- c(d13C = 1.3, d15N = 3.5, d34S = 1.3)
TEF_sds   <- c(d13C = 0.5, d15N = 0.5, d34S = 1.3)

cat("TEFs applied:\n")
cat("  delta-13C:", TEF_means["d13C"], "±", TEF_sds["d13C"], "‰  (Lerner et al. 2021)\n")
cat("  delta-15N:", TEF_means["d15N"], "±", TEF_sds["d15N"], "‰  (Lerner et al. 2021)\n")
cat("  delta-34S:", TEF_means["d34S"], "±", TEF_sds["d34S"], "‰  (McCutchan et al. 2003)\n\n")

correction_means_5 <- matrix(rep(TEF_means, nrow(sources_5)),
                              nrow = nrow(sources_5), ncol = 3, byrow = TRUE)
rownames(correction_means_5) <- source_names_5
colnames(correction_means_5) <- c("d13C", "d15N", "d34S")

correction_sds_5 <- matrix(rep(TEF_sds, nrow(sources_5)),
                             nrow = nrow(sources_5), ncol = 3, byrow = TRUE)
rownames(correction_sds_5) <- source_names_5
colnames(correction_sds_5) <- c("d13C", "d15N", "d34S")

# =============================================================================
# 9. MODEL 7: MUSCLE, 5 SOURCES, WITHOUT COVARIATE (for comparison)
# =============================================================================

cat("\n==============================================================\n")
cat("MODEL 7: MUSCLE, 5 SOURCES, WITHOUT BODY MASS COVARIATE\n")
cat("==============================================================\n\n")

cosimmr_7 <- with(
  consumers_muscle,
  cosimmr_load(
    formula           = cbind(d13C_corrected, d15N, d34S) ~ 1,
    source_names      = source_names_5,
    source_means      = source_means_5,
    source_sds        = source_sds_5,
    correction_means  = correction_means_5,
    correction_sds    = correction_sds_5,
    concentration_means = NULL
  )
)

set.seed(26051993)
cosimmr_7_out <- cosimmr_ffvb(cosimmr_7)
saveRDS(cosimmr_7_out, "output/cosimmr_model_7_muscle_5sources_nocov.rds")

print(summary(cosimmr_7_out, type = "statistics"))
plot(cosimmr_7_out, type = "prop_histogram", obs = 1)

# =============================================================================
# 10. MODEL 8: MUSCLE, 5 SOURCES, WITH BODY MASS COVARIATE  ← SELECTED MODEL
# =============================================================================

cat("\n==============================================================\n")
cat("MODEL 8: MUSCLE, 5 SOURCES, WITH BODY MASS COVARIATE (selected)\n")
cat("sigma = 0.762; 95% CI width = 0.383\n")
cat("==============================================================\n\n")

cosimmr_8 <- with(
  consumers_muscle,
  cosimmr_load(
    formula           = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
    source_names      = source_names_5,
    source_means      = source_means_5,
    source_sds        = source_sds_5,
    correction_means  = correction_means_5,
    correction_sds    = correction_sds_5,
    concentration_means = NULL
  )
)

set.seed(26051993)
cosimmr_8_out <- cosimmr_ffvb(cosimmr_8)
saveRDS(cosimmr_8_out, "output/cosimmr_model_8_muscle_5sources_cov.rds")

print(summary(cosimmr_8_out, type = "statistics"))
plot(cosimmr_8, colour_by_cov = TRUE, cov_name = "Mass_g")
plot(cosimmr_8_out, type = "covariates_plot", cov_name = "Mass_g", one_plot = TRUE, n_pred = 100)

# =============================================================================
# 11. EXPORT RESULTS
# =============================================================================

cat("\n==============================================================\n")
cat("EXPORTING RESULTS\n")
cat("==============================================================\n\n")

write.csv(sources_5, "output/sources_5.csv", row.names = FALSE)

sink("output/mixing_model_results.txt")
cat("COSIMMR MIXING MODEL RESULTS\n")
cat("=============================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("--- SOURCES (5) ---\n")
print(sources_5)

cat("\n--- TEFs ---\n")
cat("delta-13C:", TEF_means["d13C"], "±", TEF_sds["d13C"], "‰\n")
cat("delta-15N:", TEF_means["d15N"], "±", TEF_sds["d15N"], "‰\n")
cat("delta-34S:", TEF_means["d34S"], "±", TEF_sds["d34S"], "‰\n\n")

cat("--- MODEL 7: Muscle, 5 sources, no covariate ---\n")
print(summary(cosimmr_7_out, type = "statistics"))

cat("\n--- MODEL 8: Muscle, 5 sources, with body mass covariate (SELECTED) ---\n")
print(summary(cosimmr_8_out, type = "statistics"))

sink()

cat("✅ Results exported to output/mixing_model_results.txt\n\n")
cat("SUMMARY:\n")
cat("  - Selected model: Model 8 (muscle, 5 sources, body mass covariate)\n")
cat("  - Muscle consumers:", nrow(consumers_muscle), "\n")
cat("  - Sources:", paste(source_names_5, collapse = ", "), "\n")
