# =============================================================================
# ANÁLISIS EXPLORATORIO COMPLETO: CORRELACIONES Y ESTADÍSTICAS
# Sistema pelágico: Salmón Chinook - Isótopos y Morfometría
# =============================================================================

library(tidyverse)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(patchwork)
library(vegan)

# =============================================================================
# 1. CARGAR Y PREPARAR DATOS
# =============================================================================

cat("==============================================================\n")
cat("CARGANDO DATOS\n")
cat("==============================================================\n\n")

# Cargar datos
raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

# Crear columna de grupo/especie
raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^Sar") ~ "Sardine",
  str_detect(raw_data$Sample, "^Zoo") ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito") ~ "Phytoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchovy",
  TRUE ~ "Unknown"
)

# Calcular suma de porcentajes
raw_data$pct_sum <- raw_data$pctN + raw_data$pctC + raw_data$pctS

cat("Dataset creado exitosamente\n")
cat("Dimensiones:", nrow(raw_data), "filas x", ncol(raw_data), "columnas\n")
print(table(raw_data$group))

# =============================================================================
# 2. CONTROL DE CALIDAD (QC)
# =============================================================================

cat("\n==============================================================\n")
cat("CONTROL DE CALIDAD\n")
cat("==============================================================\n\n")

# Crear tabla de flags QC
qc_flags <- raw_data %>%
  mutate(
    flag_pct_sum = ifelse(pct_sum > 100, "suma_%>100", ""),
    flag_CN_extreme = case_when(
      group == "Salmon_muscle" & CtoN > 10 ~ "C/N_extremo",
      group %in% c("Sardine", "Anchovy") & CtoN > 8 ~ "C/N_extremo", 
      TRUE ~ ""
    ),
    flag_d34S_missing = ifelse(is.na(d34S), "δ34S_faltante", ""),
    all_flags = paste(flag_pct_sum, flag_CN_extreme, flag_d34S_missing, sep = ";"),
    all_flags = str_remove_all(all_flags, "^;+|;+$|;{2,}"),
    all_flags = ifelse(all_flags == "", "OK", all_flags)
  )

# Mostrar problemas identificados
problems <- qc_flags %>%
  filter(all_flags != "OK") %>%
  select(ID, Sample, group, pct_sum, CtoN, d34S, all_flags)

if(nrow(problems) > 0) {
  cat("Muestras con problemas identificados:\n")
  print(problems)
} else {
  cat("No se encontraron problemas en QC\n")
}

# Decisiones de exclusión
exclude_samples <- qc_flags %>%
  filter(str_detect(all_flags, "suma_%>100")) %>%
  pull(ID)

if(length(exclude_samples) > 0) {
  cat("\nMuestras excluidas por QC:\n")
  print(qc_flags[qc_flags$ID %in% exclude_samples, c("ID", "Sample", "group", "pct_sum", "all_flags")])
} else {
  cat("\nNinguna muestra excluida por QC\n")
}

# Dataset limpio
clean_data <- qc_flags %>%
  filter(!ID %in% exclude_samples)

cat("\nDataset después de QC:", nrow(clean_data), "muestras\n")
print(table(clean_data$group))

# =============================================================================
# 3. CORRECCIÓN LIPÍDICA
# =============================================================================

cat("\n==============================================================\n")
cat("CORRECCIÓN LIPÍDICA\n")
cat("==============================================================\n\n")

# McConnaughey optimizado para Chinook salmon
lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent <- 93 / (1 + (1/lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  final_correction <- ifelse(CN_ratio > 3.5, correction_factor, 0)
  return(d13C + final_correction)
}

# Aplicar correcciones
clean_data <- clean_data %>%
  mutate(
    d13C_corrected = case_when(
      group == "Salmon_muscle" ~ lipid_correction_optimized_chinook(d13C, CtoN),
      TRUE ~ d13C
    )
  )

# Resumen del impacto
correction_summary <- clean_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    d13C_original_mean = round(mean(d13C, na.rm = TRUE), 2),
    d13C_corrected_mean = round(mean(d13C_corrected, na.rm = TRUE), 2),
    mean_correction = round(mean(d13C_corrected - d13C, na.rm = TRUE), 2),
    .groups = 'drop'
  )

cat("Impacto de la corrección lipídica:\n")
print(correction_summary)

# =============================================================================
# 4. PREPARAR DATOS DE SALMONES
# =============================================================================

cat("\n==============================================================\n")
cat("PREPARANDO DATOS DE SALMONES PARA ANÁLISIS\n")
cat("==============================================================\n\n")

# Datos específicos de salmones con edad (solo músculo)
salmon_data <- clean_data %>%
  filter(group == "Salmon_muscle") %>%
  filter(!is.na(age)) %>%
  mutate(
    age_factor = factor(age, levels = c(1, 2, 3))
  )

cat("Salmones con edad:", nrow(salmon_data), "\n")
print(table(salmon_data$age_factor))

# Dataset completo para correlaciones
salmon_complete <- salmon_data %>%
  select(age, fork_length, total_length, weight, 
         d13C_corrected, d15N, d34S, CtoN, pctN, pctC) %>%
  filter(complete.cases(.))

cat("\nDatos completos para correlaciones:", nrow(salmon_complete), "individuos\n\n")

# =============================================================================
# 5. ESTADÍSTICAS DESCRIPTIVAS
# =============================================================================

cat("==============================================================\n")
cat("ESTADÍSTICAS DESCRIPTIVAS POR EDAD\n")
cat("==============================================================\n\n")

age_stats <- salmon_data %>%
  group_by(age_factor) %>%
  summarise(
    n = n(),
    # Morfometría
    mean_length = round(mean(fork_length, na.rm = TRUE), 1),
    sd_length = round(sd(fork_length, na.rm = TRUE), 1),
    mean_weight = round(mean(weight, na.rm = TRUE), 0),
    sd_weight = round(sd(weight, na.rm = TRUE), 0),
    # Isótopos
    mean_d13C = round(mean(d13C_corrected, na.rm = TRUE), 2),
    sd_d13C = round(sd(d13C_corrected, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
    mean_d34S = round(mean(d34S, na.rm = TRUE), 2),
    sd_d34S = round(sd(d34S, na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(age_stats)

# Guardar estadísticas
write.csv(age_stats, "output/statistics_by_age.csv", row.names = FALSE)
cat("\n✅ Guardado: output/statistics_by_age.csv\n")

# =============================================================================
# 6. ANÁLISIS DE CORRELACIONES - EDAD vs ISÓTOPOS
# =============================================================================

cat("\n==============================================================\n")
cat("CORRELACIONES EDAD vs ISÓTOPOS\n")
cat("==============================================================\n\n")

# δ13C vs edad
cor_d13C_age <- cor.test(salmon_complete$age, salmon_complete$d13C_corrected, 
                         method = "spearman")
cat("δ13C vs Edad:\n")
cat(sprintf("  Spearman r = %.3f, p = %.4f", 
            as.numeric(cor_d13C_age$estimate), 
            cor_d13C_age$p.value))
if(cor_d13C_age$p.value < 0.001) {
  cat(" ***\n")
} else if(cor_d13C_age$p.value < 0.01) {
  cat(" **\n")
} else if(cor_d13C_age$p.value < 0.05) {
  cat(" *\n")
} else {
  cat("\n")
}

# δ15N vs edad
cor_d15N_age <- cor.test(salmon_complete$age, salmon_complete$d15N, 
                         method = "spearman")
cat("\nδ15N vs Edad:\n")
cat(sprintf("  Spearman r = %.3f, p = %.4f", 
            as.numeric(cor_d15N_age$estimate), 
            cor_d15N_age$p.value))
if(cor_d15N_age$p.value < 0.001) {
  cat(" ***\n")
} else if(cor_d15N_age$p.value < 0.01) {
  cat(" **\n")
} else if(cor_d15N_age$p.value < 0.05) {
  cat(" *\n")
} else {
  cat("\n")
}

# δ34S vs edad (si hay datos suficientes)
if(sum(!is.na(salmon_complete$d34S)) > 5) {
  cor_d34S_age <- cor.test(salmon_complete$age, salmon_complete$d34S, 
                           method = "spearman")
  cat("\nδ34S vs Edad:\n")
  cat(sprintf("  Spearman r = %.3f, p = %.4f", 
              as.numeric(cor_d34S_age$estimate), 
              cor_d34S_age$p.value))
  if(cor_d34S_age$p.value < 0.001) {
    cat(" ***\n")
  } else if(cor_d34S_age$p.value < 0.01) {
    cat(" **\n")
  } else if(cor_d34S_age$p.value < 0.05) {
    cat(" *\n")
  } else {
    cat("\n")
  }
}

# δ34S vs edad (si hay datos suficientes)
if(sum(!is.na(salmon_complete$d34S)) > 5) {
  cor_d34S_age <- cor.test(salmon_complete$age, salmon_complete$d34S, 
                           method = "spearman")
  cat("\nδ34S vs Edad:\n")
  cat(sprintf("  Spearman r = %.3f, p = %.4f", 
              cor_d34S_age$estimate, cor_d34S_age$p.value))
  if(cor_d34S_age$p.value < 0.001) cat(" ***\n")
  else if(cor_d34S_age$p.value < 0.01) cat(" **\n")
  else if(cor_d34S_age$p.value < 0.05) cat(" *\n")
  else cat("\n")
}

# =============================================================================
# 7. MATRIZ DE CORRELACIONES COMPLETA
# =============================================================================

cat("\n==============================================================\n")
cat("MATRIZ DE CORRELACIONES COMPLETA\n")
cat("==============================================================\n\n")

# Calcular correlaciones con Spearman
salmon_complete$total_length <- as.numeric(salmon_complete$total_length)  # Asegurar tipo numérico


cor_matrix <- cor(salmon_complete, method = "spearman", use = "complete.obs")

cat("Matriz de correlaciones (Spearman):\n")
print(round(cor_matrix, 3))

# Calcular p-values
cor_test_results <- rcorr(as.matrix(salmon_complete), type = "spearman")
cor_pvalues <- cor_test_results$P

cat("\n\nP-values:\n")
print(round(cor_pvalues, 4))

# =============================================================================
# 8. TABLA RESUMEN PARA MANUSCRITO
# =============================================================================

cat("\n==============================================================\n")
cat("TABLA RESUMEN PARA MANUSCRITO\n")
cat("==============================================================\n\n")

# Variables isotópicas
isotopes <- c("d13C_corrected", "d15N", "d34S")

# Tabla de correlaciones con edad
cor_summary <- data.frame()

for(iso in isotopes) {
  if(sum(!is.na(salmon_complete[[iso]])) > 5) {
    test <- cor.test(salmon_complete$age, salmon_complete[[iso]], 
                     method = "spearman")
    
    cor_summary <- rbind(cor_summary, data.frame(
      Variable = iso,
      Spearman_r = round(test$estimate, 3),
      p_value = round(test$p.value, 4),
      p_formatted = formatC(test$p.value, format = "f", digits = 4),
      Significance = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01 ~ "**",
        test$p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      stringsAsFactors = FALSE
    ))
  }
}

print(cor_summary)

# Guardar tabla
write.csv(cor_summary, "output/correlations_age_isotopes.csv", row.names = FALSE)
cat("\n✅ Guardado: output/correlations_age_isotopes.csv\n")

# =============================================================================
# 9. TESTS ESTADÍSTICOS ADICIONALES
# =============================================================================

cat("\n==============================================================\n")
cat("TESTS ESTADÍSTICOS ADICIONALES\n")
cat("==============================================================\n\n")

# Kruskal-Wallis por edad
cat("DIFERENCIAS POR EDAD (Kruskal-Wallis):\n\n")

variables_test <- c("d13C_corrected", "d15N", "d34S", "fork_length", "weight")

for(var in variables_test) {
  if(sum(!is.na(salmon_data[[var]])) > 5) {
    kw_test <- kruskal.test(salmon_data[[var]] ~ salmon_data$age_factor)
    cat(sprintf("%s: χ² = %.2f, df = %d, p = %.4f", 
                var, kw_test$statistic, kw_test$parameter, kw_test$p.value))
    if(kw_test$p.value < 0.001) cat(" ***\n")
    else if(kw_test$p.value < 0.01) cat(" **\n")
    else if(kw_test$p.value < 0.05) cat(" *\n")
    else cat("\n")
  }
}

# PERMANOVA multivariado (si hay suficientes datos)
if(nrow(salmon_data) > 10) {
  cat("\n\nPERMANOVA - Composición isotópica ~ Edad:\n")
  
  iso_matrix <- salmon_data %>%
    select(d13C_corrected, d15N, d34S) %>%
    filter(complete.cases(.))
  
  age_vector <- salmon_data$age_factor[complete.cases(salmon_data[c("d13C_corrected", "d15N", "d34S")])]
  
  if(nrow(iso_matrix) > 5 && length(unique(age_vector)) > 1) {
    perm_age <- adonis2(iso_matrix ~ age_vector, permutations = 999, method = "euclidean")
    print(perm_age)
    
    # Guardar resultado
    capture.output(perm_age, file = "output/permanova_results.txt")
    cat("\n✅ Guardado: output/permanova_results.txt\n")
  } else {
    cat("Insuficientes datos para PERMANOVA\n")
  }
}

# =============================================================================
# 10. VISUALIZACIONES
# =============================================================================

cat("\n==============================================================\n")
cat("GENERANDO VISUALIZACIONES\n")
cat("==============================================================\n\n")

# Crear directorio si no existe
if(!dir.exists("output")) dir.create("output")

# 10.1 Matriz de correlaciones (corrplot)
cat("Generando matriz de correlaciones...\n")

png("output/correlation_matrix_complete.png", width = 12, height = 12, 
    units = "in", res = 300)

corrplot(cor_matrix, method = "color", type = "upper", 
         order = "original", 
         tl.col = "black", tl.srt = 45,
         addCoef.col = "black", 
         number.cex = 0.7,
         col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
         title = "Correlation Matrix: Morphometry and Isotopes (Spearman)",
         mar = c(0,0,2,0),
         diag = FALSE)

dev.off()
cat("✅ Guardado: output/correlation_matrix_complete.png\n")

# 10.2 Correlaciones edad vs isótopos (individuales)
cat("\nGenerando gráficos individuales edad vs isótopos...\n")

# δ13C vs edad
p_d13C <- ggplot(salmon_data, aes(x = age, y = d13C_corrected)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue") +
  labs(
    title = "δ13C vs Age",
    subtitle = sprintf("Spearman r = %.3f, p = %.4f", 
                      cor_d13C_age$estimate, cor_d13C_age$p.value),
    x = "Age (years)",
    y = "δ13C corrected (‰)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11)
  )

# δ15N vs edad
p_d15N <- ggplot(salmon_data, aes(x = age, y = d15N)) +
  geom_point(size = 3, alpha = 0.7, color = "coral") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    title = "δ15N vs Age",
    subtitle = sprintf("Spearman r = %.3f, p = %.4f", 
                      cor_d15N_age$estimate, cor_d15N_age$p.value),
    x = "Age (years)",
    y = "δ15N (‰)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11)
  )

# Combinar gráficos
p_combined <- p_d13C | p_d15N

ggsave("output/correlations_age_isotopes.png", p_combined, 
       width = 12, height = 5, dpi = 300)

cat("✅ Guardado: output/correlations_age_isotopes.png\n")

# 10.3 Boxplots por edad
cat("\nGenerando boxplots por edad...\n")

p_boxplots <- salmon_data %>%
  select(age_factor, d13C_corrected, d15N, d34S) %>%
  pivot_longer(cols = c(d13C_corrected, d15N, d34S), 
               names_to = "isotope", values_to = "value") %>%
  mutate(isotope_label = case_when(
    isotope == "d13C_corrected" ~ "δ13C (‰)",
    isotope == "d15N" ~ "δ15N (‰)",
    isotope == "d34S" ~ "δ34S (‰)"
  )) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = age_factor, y = value, fill = age_factor)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~isotope_label, scales = "free_y", ncol = 3) +
  scale_fill_viridis_d(name = "Age") +
  labs(
    title = "Isotopic Composition by Age",
    x = "Age (years)",
    y = "Isotopic Value"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

ggsave("output/boxplots_isotopes_by_age.png", p_boxplots, 
       width = 12, height = 4, dpi = 300)

cat("✅ Guardado: output/boxplots_isotopes_by_age.png\n")

# =============================================================================
# 11. RESUMEN EJECUTIVO
# =============================================================================

cat("\n==============================================================\n")
cat("GENERANDO RESUMEN EJECUTIVO\n")
cat("==============================================================\n\n")

sink("output/resumen_ejecutivo_correlaciones.txt")

cat("RESUMEN EJECUTIVO: ANÁLISIS DE CORRELACIONES\n")
cat("=============================================\n\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("MUESTRAS:\n")
cat("  - Total salmones con edad:", nrow(salmon_data), "\n")
cat("  - Edad 1:", sum(salmon_data$age == 1), "\n")
cat("  - Edad 2:", sum(salmon_data$age == 2), "\n")
cat("  - Edad 3:", sum(salmon_data$age == 3), "\n\n")

cat("CORRELACIONES EDAD vs ISÓTOPOS (Spearman):\n\n")
print(cor_summary)

cat("\n\nESTADÍSTICAS DESCRIPTIVAS POR EDAD:\n\n")
print(age_stats)

cat("\n\nMATRIZ DE CORRELACIONES COMPLETA:\n\n")
print(round(cor_matrix, 3))

cat("\n\nARCHIVOS GENERADOS:\n")
cat("  - correlation_matrix_complete.png\n")
cat("  - correlations_age_isotopes.png\n")
cat("  - boxplots_isotopes_by_age.png\n")
cat("  - correlations_age_isotopes.csv\n")
cat("  - statistics_by_age.csv\n")

sink()

cat("✅ Guardado: output/resumen_ejecutivo_correlaciones.txt\n")

# =============================================================================
# 12. FINALIZACIÓN
# =============================================================================

cat("\n==============================================================\n")
cat("ANÁLISIS COMPLETADO\n")
cat("==============================================================\n\n")

cat("RESULTADOS CLAVE PARA MANUSCRITO:\n\n")

cat("δ13C vs Edad:\n")
cat(sprintf("  Spearman r = %.2f, p = %.3f", 
            cor_d13C_age$estimate, cor_d13C_age$p.value))
if(cor_d13C_age$p.value < 0.01) cat(" **\n") else cat("\n")

cat("\nδ15N vs Edad:\n")
cat(sprintf("  Spearman r = %.2f, p = %.3f", 
            cor_d15N_age$estimate, cor_d15N_age$p.value))
if(cor_d15N_age$p.value < 0.05) cat(" *\n") else cat("\n")

cat("\nTodos los archivos guardados en: output/\n")
cat("\n¡Análisis completado exitosamente!\n\n")