  # =============================================================================
  # ANÁLISIS ISOTÓPICO COMPLETO: DIETA DE SALMÓN EN ZONA ESTUARINA
  # Mixing Models Bayesianos con δ13C, δ15N, δ34S usando COSIMMR
  # Tejidos: MÚSCULO + HÍGADO
  # Fuentes: Marinas + Agua Dulce
  # Covariable: Masa del pez
  # =============================================================================

  # Cargar librerías necesarias
  library(tidyverse)
  library(cosimmr)
  library(ggplot2)

  # =============================================================================
  # 1. CARGAR DATASETS
  # =============================================================================

  cat("==============================================================\n")
  cat("CARGANDO DATOS\n")
  cat("==============================================================\n\n")

  # Cargar datos individuales (salmón + fuentes marinas)
  raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
  names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

  # Crear columna de grupo/especie
  raw_data$group <- case_when(
    str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
    str_detect(raw_data$Sample, "^h-") ~ "Salmon_liver", 
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

  # Cargar fuentes de agua dulce (ya agregadas)
  cat("\n=== CARGANDO FUENTES DE AGUA DULCE ===\n")
  sources_freshwater <- read.csv("data/data_raw/diet-isotopes/datos_rios.csv")
  cat("Fuentes de agua dulce cargadas:\n")
  print(sources_freshwater)

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
        group %in% c("Salmon_muscle", "Salmon_liver") & CtoN > 10 ~ "C/N_extremo",
        group %in% c("Sardine", "Anchovy") & CtoN > 8 ~ "C/N_extremo", 
        TRUE ~ ""
      ),
      flag_d34S_missing = ifelse(is.na(d34S), "δ34S_faltante", ""),
      flag_impossible = case_when(
        pctN < 0 | pctN > 100 ~ "N_imposible",
        pctC < 0 | pctC > 100 ~ "C_imposible", 
        pctS < 0 | pctS > 100 ~ "S_imposible",
        TRUE ~ ""
      ),
      all_flags = paste(flag_pct_sum, flag_CN_extreme, flag_d34S_missing, flag_impossible, sep = ";"),
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

  # Dataset limpio: MANTENER MÚSCULO + HÍGADO + fuentes
  clean_data <- qc_flags %>%
    filter(!ID %in% exclude_samples)

  cat("\nDataset después de QC (músculo + hígado + fuentes):", nrow(clean_data), "muestras\n")
  print(table(clean_data$group))

  # =============================================================================
  # 3. CORRECCIÓN LIPÍDICA PARA δ13C
  # =============================================================================

  cat("\n==============================================================\n")
  cat("CORRECCIÓN LIPÍDICA\n")
  cat("==============================================================\n\n")
  cat("MÉTODOS:\n")
  cat("  - Salmón (músculo + hígado): McConnaughey optimizado para Chinook\n")
  cat("  - Sardina/Anchoveta: Ecuación Clupeiformes (Sardenne et al. 2023)\n")
  cat("  - Zooplancton/Fitoplancton: SIN corrección\n\n")

  # McConnaughey optimizado para Chinook salmon
  lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
    lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
    L_percent <- 93 / (1 + (1/lipid_term))
    correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
    final_correction <- ifelse(CN_ratio > 1, correction_factor, 0)
    return(d13C + final_correction)
  }

  # Corrección lipídica para Clupeiformes (sardina/anchoveta)
  lipid_correction_clupeiformes <- function(d13C_bulk, CN_bulk) {
    d13C_corrected <- d13C_bulk + 0.894 * CN_bulk - 2.377
    return(d13C_corrected)
  }

  # Aplicar correcciones
  clean_data <- clean_data %>%
    mutate(
      needs_correction = group %in% c("Salmon_muscle", "Salmon_liver", "Sardine", "Anchovy"),
      d13C_corrected = case_when(
        # Salmón (músculo y hígado): parámetros optimizados para Chinook
        group %in% c("Salmon_muscle", "Salmon_liver") ~ lipid_correction_optimized_chinook(d13C, CtoN),
        
        # Sardina y Anchoveta: ecuación Clupeiformes
        # group %in% c("Sardine", "Anchovy") ~ lipid_correction_clupeiformes(d13C, CtoN),
        
        # Zooplancton y Fitoplancton: SIN corrección
        TRUE ~ d13C
      )
    )

  # Resumen del impacto
  correction_summary <- clean_data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      needs_correction = first(needs_correction),
      d13C_original_mean = round(mean(d13C, na.rm = T), 2),
      d13C_corrected_mean = round(mean(d13C_corrected, na.rm = T), 2),
      mean_correction = round(mean(d13C_corrected - d13C, na.rm = T), 2),
      max_correction = round(max(d13C_corrected - d13C, na.rm = T), 2),
      .groups = 'drop'
    )

  cat("Impacto de la corrección lipídica por grupo:\n")
  print(correction_summary)

  # =============================================================================
  # 4. PREPARAR CONSUMIDORES (MÚSCULO Y HÍGADO)
  # =============================================================================

  cat("\n==============================================================\n")
  cat("PREPARACIÓN DE CONSUMIDORES\n")
  cat("==============================================================\n\n")

  # Consumidores: músculo
  consumers_muscle <- clean_data %>%
    filter(group == "Salmon_muscle") %>%
    select(ID, Sample, age, fork_length, weight, d13C_corrected, d15N, d34S, CtoN) %>%
    filter(!is.na(d13C_corrected) & !is.na(d15N) & !is.na(d34S)) %>%
    rename(Mass_g = weight)

  cat("Consumidores MÚSCULO:", nrow(consumers_muscle), "individuos\n")
  cat("Rango de masa:", min(consumers_muscle$Mass_g, na.rm=T), "-", 
      max(consumers_muscle$Mass_g, na.rm=T), "g\n\n")

  # Consumidores: hígado
  consumers_liver <- clean_data %>%
    filter(group == "Salmon_liver") %>%
    select(ID, Sample, age, fork_length, weight, d13C_corrected, d15N, d34S, CtoN) %>%
    filter(!is.na(d13C_corrected) & !is.na(d15N) & !is.na(d34S)) %>%
    rename(Mass_g = weight)

  cat("Consumidores HÍGADO:", nrow(consumers_liver), "individuos\n")
  cat("Rango de masa:", min(consumers_liver$Mass_g, na.rm=T), "-", 
      max(consumers_liver$Mass_g, na.rm=T), "g\n\n")

  # =============================================================================
  # 5. CALCULAR FUENTES MARINAS
  # =============================================================================

  cat("=== CALCULANDO FUENTES MARINAS ===\n\n")

  sources_marine <- clean_data %>%
    filter(group %in% c("Sardine", "Anchovy", "Zooplankton", "Phytoplankton")) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      d13C_mean = mean(d13C_corrected, na.rm = TRUE),
      d13C_sd = sd(d13C_corrected, na.rm = TRUE),
      d15N_mean = mean(d15N, na.rm = TRUE),
      d15N_sd = sd(d15N, na.rm = TRUE),
      d34S_mean = mean(d34S, na.rm = TRUE),
      d34S_sd = sd(d34S, na.rm = TRUE),
      habitat = "marine",
      .groups = 'drop'
    )

  cat("Fuentes marinas calculadas:\n")
  print(sources_marine)

  # =============================================================================
  # 6. PREPARAR FUENTES DE AGUA DULCE
  # =============================================================================

  cat("\n=== PREPARANDO FUENTES DE AGUA DULCE ===\n\n")

  sources_freshwater <- sources_freshwater %>%
    rename(group = Source) %>%
    mutate(
      habitat = "freshwater",
      n = NA
    ) %>%
    select(group, n, d13C_mean, d13C_sd, d15N_mean, d15N_sd, d34S_mean, d34S_sd, habitat)

  cat("Fuentes de agua dulce:\n")
  print(sources_freshwater)

  # =============================================================================
  # 7. COMBINAR TODAS LAS FUENTES
  # =============================================================================

  cat("\n=== COMBINANDO TODAS LAS FUENTES ===\n\n")

  sources_all <- bind_rows(sources_marine, sources_freshwater)

  cat("Total de fuentes disponibles:", nrow(sources_all), "\n")
  cat("Fuentes completas:\n")
  print(sources_all)

  # =============================================================================
  # 8. DEFINIR TEFs
  # =============================================================================

  cat("\n==============================================================\n")
  cat("FACTORES DE ENRIQUECIMIENTO TRÓFICO (TEFs)\n")
  cat("==============================================================\n\n")

  TEF_means <- c(d13C = 1.3, d15N = 3.5, d34S = 1.3)
  TEF_sds <- c(d13C = 0.5, d15N = 0.5, d34S = 1.3)

  cat("TEFs aplicados:\n")
  cat("  δ13C:", TEF_means["d13C"], "±", TEF_sds["d13C"], "‰\n")
  cat("  δ15N:", TEF_means["d15N"], "±", TEF_sds["d15N"], "‰\n")
  cat("  δ34S:", TEF_means["d34S"], "±", TEF_sds["d34S"], "‰\n\n")

  # =============================================================================
  # 9. PREPARAR MATRICES PARA COSIMMR
  # =============================================================================

  cat("==============================================================\n")
  cat("PREPARANDO DATOS PARA COSIMMR\n")
  cat("==============================================================\n\n")

  source_names <- sources_all$group

  source_means <- as.matrix(sources_all[, c("d13C_mean", "d15N_mean", "d34S_mean")])
  rownames(source_means) <- source_names

  source_sds <- as.matrix(sources_all[, c("d13C_sd", "d15N_sd", "d34S_sd")])
  rownames(source_sds) <- source_names

  correction_means <- matrix(rep(TEF_means, nrow(sources_all)), 
                            nrow = nrow(sources_all), ncol = 3, byrow = TRUE)
  rownames(correction_means) <- source_names
  colnames(correction_means) <- c("d13C", "d15N", "d34S")

  correction_sds <- matrix(rep(TEF_sds, nrow(sources_all)), 
                          nrow = nrow(sources_all), ncol = 3, byrow = TRUE)
  rownames(correction_sds) <- source_names
  colnames(correction_sds) <- c("d13C", "d15N", "d34S")

  # Fuentes simplificadas
  sources_simplified <- sources_all %>%
    group_by(habitat) %>%
    summarise(
      n_sources = n(),
      d13C_mean = mean(d13C_mean, na.rm = TRUE),
      d13C_sd = mean(d13C_sd, na.rm = TRUE),
      d15N_mean = mean(d15N_mean, na.rm = TRUE),
      d15N_sd = mean(d15N_sd, na.rm = TRUE),
      d34S_mean = mean(d34S_mean, na.rm = TRUE),
      d34S_sd = mean(d34S_sd, na.rm = TRUE),
      .groups = 'drop'
    )

  source_names_simple <- sources_simplified$habitat
  source_means_simple <- as.matrix(sources_simplified[, c("d13C_mean", "d15N_mean", "d34S_mean")])
  rownames(source_means_simple) <- source_names_simple
  source_sds_simple <- as.matrix(sources_simplified[, c("d13C_sd", "d15N_sd", "d34S_sd")])
  rownames(source_sds_simple) <- source_names_simple

  correction_means_simple <- matrix(rep(TEF_means, 2), nrow = 2, ncol = 3, byrow = TRUE)
  rownames(correction_means_simple) <- source_names_simple
  colnames(correction_means_simple) <- c("d13C", "d15N", "d34S")

  correction_sds_simple <- matrix(rep(TEF_sds, 2), nrow = 2, ncol = 3, byrow = TRUE)
  rownames(correction_sds_simple) <- source_names_simple
  colnames(correction_sds_simple) <- c("d13C", "d15N", "d34S")

  # =============================================================================
  # PARTE A: MODELOS PARA MÚSCULO
  # =============================================================================

  cat("\n==============================================================\n")
  cat("PARTE A: ANÁLISIS DE MÚSCULO\n")
  cat("==============================================================\n\n")

  # =============================================================================
  # MODELO 1: MÚSCULO SIN COVARIABLE
  # =============================================================================

  cat("--- MODELO 1: MÚSCULO SIN COVARIABLE ---\n\n")

  cosimmr_1 <- with(
    consumers_muscle, 
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ 1,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 1...\n")
  saveRDS(cosimmr_1_out,"output/cosimmr_model_1_muscle.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_1_out, type = "statistics"))
  plot(cosimmr_1_out, type = "prop_histogram", obs = 1)

  # =============================================================================
  # MODELO 2: MÚSCULO CON COVARIABLE (DETALLADO)
  # =============================================================================

  cat("\n--- MODELO 2: MÚSCULO CON COVARIABLE (detallado) ---\n\n")

  cosimmr_2 <- with(
    consumers_muscle,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_2_out <- cosimmr_ffvb(cosimmr_2)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 2...\n")
  saveRDS(cosimmr_2_out,"output/cosimmr_model_2_muscle.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_2_out, type = "statistics"))
  plot(cosimmr_2, colour_by_cov = TRUE, cov_name = 'Mass_g')
  plot(cosimmr_2_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # MODELO 3: MÚSCULO SIMPLIFICADO (MARINO VS DULCE)
  # =============================================================================

  cat("\n--- MODELO 3: MÚSCULO SIMPLIFICADO (Marino vs Dulce) ---\n\n")

  cosimmr_3 <- with(
    consumers_muscle,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_simple,
      source_means = source_means_simple,
      source_sds = source_sds_simple,
      correction_means = correction_means_simple,
      correction_sds = correction_sds_simple,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_3_out <- cosimmr_ffvb(cosimmr_3)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 3...\n")
  saveRDS(cosimmr_3_out,"output/cosimmr_model_3_muscle.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_3_out, type = "statistics"))
  plot(cosimmr_3_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # PARTE B: MODELOS PARA HÍGADO
  # =============================================================================

  cat("\n==============================================================\n")
  cat("PARTE B: ANÁLISIS DE HÍGADO\n")
  cat("==============================================================\n\n")

  # =============================================================================
  # MODELO 4: HÍGADO SIN COVARIABLE
  # =============================================================================

  cat("--- MODELO 4: HÍGADO SIN COVARIABLE ---\n\n")

  cosimmr_4 <- with(
    consumers_liver, 
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ 1,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_4_out <- cosimmr_ffvb(cosimmr_4)
  cat("✅ Completado\n\n")


  cat("Guardando Modelo FFVB 4...\n")
  saveRDS(cosimmr_4_out,"output/cosimmr_model_4_liver.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_4_out, type = "statistics"))
  plot(cosimmr_4_out, type = "prop_histogram", obs = 1)

  # =============================================================================
  # MODELO 5: HÍGADO CON COVARIABLE (DETALLADO)
  # =============================================================================

  cat("\n--- MODELO 5: HÍGADO CON COVARIABLE (detallado) ---\n\n")

  cosimmr_5 <- with(
    consumers_liver,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_5_out <- cosimmr_ffvb(cosimmr_5)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 5...\n")
  saveRDS(cosimmr_5_out,"output/cosimmr_moelel_5_liver.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_5_out, type = "statistics"))
  plot(cosimmr_5, colour_by_cov = TRUE, cov_name = 'Mass_g')
  plot(cosimmr_5_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # MODELO 6: HÍGADO SIMPLIFICADO (MARINO VS DULCE)
  # =============================================================================

  cat("\n--- MODELO 6: HÍGADO SIMPLIFICADO (Marino vs Dulce) ---\n\n")

  cosimmr_6 <- with(
    consumers_liver,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_simple,
      source_means = source_means_simple,
      source_sds = source_sds_simple,
      correction_means = correction_means_simple,
      correction_sds = correction_sds_simple,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_6_out <- cosimmr_ffvb(cosimmr_6)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 6...\n")
  saveRDS(cosimmr_6_out,"output/cosimmr_moelel_6_liver.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_6_out, type = "statistics"))
  plot(cosimmr_6_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # COMPARACIÓN MÚSCULO VS HÍGADO
  # =============================================================================

  cat("\n==============================================================\n")
  cat("COMPARACIÓN MÚSCULO VS HÍGADO\n")
  cat("==============================================================\n\n")

  cat("INTERPRETACIÓN:\n")
  cat("  - MÚSCULO: Integra dieta de ~2-6 meses\n")
  cat("  - HÍGADO: Refleja dieta de ~2-4 semanas\n")
  cat("  - Diferencias = Cambio dietario reciente\n\n")

  # Extraer proporciones promedio
  cat("--- PROPORCIONES PROMEDIO (sin covariable) ---\n\n")

  # Músculo
  posterior_muscle <- cosimmr_1_out$output$BUGSoutput$sims.list$p_mean
  props_muscle <- colMeans(posterior_muscle)
  names(props_muscle) <- source_names

  # Hígado
  posterior_liver <- cosimmr_4_out$output$BUGSoutput$sims.list$p_mean
  props_liver <- colMeans(posterior_liver)
  names(props_liver) <- source_names

  # Comparación
  comparison_df <- data.frame(
    Source = source_names,
    Muscle_mean = round(props_muscle, 3),
    Liver_mean = round(props_liver, 3),
    Difference = round(props_liver - props_muscle, 3)
  )

  print(comparison_df)

  # Gráfico comparativo
  comparison_long <- comparison_df %>%
    select(Source, Muscle_mean, Liver_mean) %>%
    pivot_longer(cols = c(Muscle_mean, Liver_mean), 
                names_to = "Tissue", values_to = "Proportion")

  p_comparison <- ggplot(comparison_long, aes(x = Source, y = Proportion, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Muscle_mean" = "darkred", "Liver_mean" = "darkorange"),
                      labels = c("Músculo", "Hígado")) +
    labs(
      title = "Comparación Músculo vs Hígado",
      subtitle = "Proporciones dietarias promedio",
      x = "Fuente",
      y = "Proporción",
      fill = "Tejido"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p_comparison)

  # =============================================================================
  # EXPORTAR RESULTADOS
  # =============================================================================

  cat("\n==============================================================\n")
  cat("EXPORTANDO RESULTADOS\n")
  cat("==============================================================\n\n")

  # if(!dir.exists("output")) dir.create("output")

  # Resumen completo
  sink("output/resultados_completo_musculo_higado.txt")
  cat("RESULTADOS COMPLETOS: MÚSCULO + HÍGADO\n")
  cat("======================================\n\n")
  cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("\n--- FUENTES ---\n")
  print(sources_all)

  cat("\n\n=== MÚSCULO ===\n")
  cat("\nModelo 1 (sin covariable):\n")
  print(summary(cosimmr_1_out, type = "statistics"))

  cat("\n\nModelo 2 (con covariable Mass_g):\n")
  print(summary(cosimmr_2_out, type = "statistics"))

  cat("\n\nModelo 3 (simplificado):\n")
  print(summary(cosimmr_3_out, type = "statistics"))

  cat("\n\n=== HÍGADO ===\n")
  cat("\nModelo 4 (sin covariable):\n")
  print(summary(cosimmr_4_out, type = "statistics"))

  cat("\n\nModelo 5 (con covariable Mass_g):\n")
  print(summary(cosimmr_5_out, type = "statistics"))

  cat("\n\nModelo 6 (simplificado):\n")
  print(summary(cosimmr_6_out, type = "statistics"))

  cat("\n\n=== COMPARACIÓN MÚSCULO VS HÍGADO ===\n")
  print(comparison_df)

  sink()

  # Guardar datos
  # write.csv(consumers_muscle, "output/consumers_muscle.csv", row.names = FALSE)
  # write.csv(consumers_liver, "output/consumers_liver.csv", row.names = FALSE)
  # write.csv(sources_all, "output/sources_all.csv", row.names = FALSE)
  # write.csv(comparison_df, "output/comparison_muscle_liver.csv", row.names = FALSE)

  # # Guardar gráfico
  # ggsave("output/comparison_muscle_liver.png", plot = p_comparison, 
  #        width = 10, height = 6, dpi = 300)


  # =============================================================================
  # PARTE C: MODELOS CON 5 FUENTES (SIN PHYTOPLANKTON NI PRIMARY_PRODUCERS)
  # =============================================================================

  cat("\n==============================================================\n")
  cat("PARTE C: ANÁLISIS CON 5 FUENTES (versión original)\n")
  cat("==============================================================\n\n")

  cat("Excluyendo: Phytoplankton y Primary_producers\n")
  cat("Fuentes incluidas: Anchovy, Sardine, Zooplankton, Fish, Invertebrates\n\n")

  # Filtrar solo las 5 fuentes originales
  sources_5 <- sources_all %>%
    filter(!group %in% c("Phytoplankton", "Primary_producers"))

  cat("Fuentes (5):\n")
  print(sources_5)

  # Preparar matrices para 5 fuentes
  source_names_5 <- sources_5$group

  source_means_5 <- as.matrix(sources_5[, c("d13C_mean", "d15N_mean", "d34S_mean")])
  rownames(source_means_5) <- source_names_5

  source_sds_5 <- as.matrix(sources_5[, c("d13C_sd", "d15N_sd", "d34S_sd")])
  rownames(source_sds_5) <- source_names_5

  correction_means_5 <- matrix(rep(TEF_means, nrow(sources_5)), 
                              nrow = nrow(sources_5), ncol = 3, byrow = TRUE)
  rownames(correction_means_5) <- source_names_5
  colnames(correction_means_5) <- c("d13C", "d15N", "d34S")

  correction_sds_5 <- matrix(rep(TEF_sds, nrow(sources_5)), 
                            nrow = nrow(sources_5), ncol = 3, byrow = TRUE)
  rownames(correction_sds_5) <- source_names_5
  colnames(correction_sds_5) <- c("d13C", "d15N", "d34S")

  # Fuentes simplificadas (5 → 2)
  sources_5_simplified <- sources_5 %>%
    group_by(habitat) %>%
    summarise(
      n_sources = n(),
      d13C_mean = mean(d13C_mean, na.rm = TRUE),
      d13C_sd = mean(d13C_sd, na.rm = TRUE),
      d15N_mean = mean(d15N_mean, na.rm = TRUE),
      d15N_sd = mean(d15N_sd, na.rm = TRUE),
      d34S_mean = mean(d34S_mean, na.rm = TRUE),
      d34S_sd = mean(d34S_sd, na.rm = TRUE),
      .groups = 'drop'
    )

  source_names_5_simple <- sources_5_simplified$habitat
  source_means_5_simple <- as.matrix(sources_5_simplified[, c("d13C_mean", "d15N_mean", "d34S_mean")])
  rownames(source_means_5_simple) <- source_names_5_simple
  source_sds_5_simple <- as.matrix(sources_5_simplified[, c("d13C_sd", "d15N_sd", "d34S_sd")])
  rownames(source_sds_5_simple) <- source_names_5_simple

  correction_means_5_simple <- matrix(rep(TEF_means, 2), nrow = 2, ncol = 3, byrow = TRUE)
  rownames(correction_means_5_simple) <- source_names_5_simple
  colnames(correction_means_5_simple) <- c("d13C", "d15N", "d34S")

  correction_sds_5_simple <- matrix(rep(TEF_sds, 2), nrow = 2, ncol = 3, byrow = TRUE)
  rownames(correction_sds_5_simple) <- source_names_5_simple
  colnames(correction_sds_5_simple) <- c("d13C", "d15N", "d34S")

  # =============================================================================
  # MODELO 7: MÚSCULO SIN COVARIABLE (5 fuentes)
  # =============================================================================

  cat("\n--- MODELO 7: MÚSCULO SIN COVARIABLE (5 fuentes) ---\n\n")

  cosimmr_7 <- with(
    consumers_muscle, 
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ 1,
      source_names = source_names_5,
      source_means = source_means_5,
      source_sds = source_sds_5,
      correction_means = correction_means_5,
      correction_sds = correction_sds_5,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_7_out <- cosimmr_ffvb(cosimmr_7)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 7...\n")
  saveRDS(cosimmr_7_out,"output/cosimmr_moel_7_muscle_5sources.rds")
  cat("✅ Guardado con éxito\n\n")


  print(summary(cosimmr_7_out, type = "statistics"))
  plot(cosimmr_7_out, type = "prop_histogram", obs = 1)

  # =============================================================================
  # MODELO 8: MÚSCULO CON COVARIABLE (5 fuentes, detallado)
  # =============================================================================

  cat("\n--- MODELO 8: MÚSCULO CON COVARIABLE (5 fuentes) ---\n\n")

  cosimmr_8 <- with(
    consumers_muscle,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_5,
      source_means = source_means_5,
      source_sds = source_sds_5,
      correction_means = correction_means_5,
      correction_sds = correction_sds_5,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_8_out <- cosimmr_ffvb(cosimmr_8)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 8...\n")
  saveRDS(cosimmr_8_out,"output/cosimmr_moel_8_muscle_5sources.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_8_out, type = "statistics"))
  plot(cosimmr_8, colour_by_cov = TRUE, cov_name = 'Mass_g')
  plot(cosimmr_8_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # MODELO 9: MÚSCULO SIMPLIFICADO (5 fuentes → 2)
  # =============================================================================

  cat("\n--- MODELO 9: MÚSCULO SIMPLIFICADO (5 fuentes → Marino vs Dulce) ---\n\n")

  cosimmr_9 <- with(
    consumers_muscle,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_5_simple,
      source_means = source_means_5_simple,
      source_sds = source_sds_5_simple,
      correction_means = correction_means_5_simple,
      correction_sds = correction_sds_5_simple,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_9_out <- cosimmr_ffvb(cosimmr_9)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 9...\n")
  saveRDS(cosimmr_9_out,"output/cosimmr_moel_9_muscle_5sources.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_9_out, type = "statistics"))
  plot(cosimmr_9_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # MODELO 10: HÍGADO SIN COVARIABLE (5 fuentes)
  # =============================================================================

  cat("\n--- MODELO 10: HÍGADO SIN COVARIABLE (5 fuentes) ---\n\n")

  cosimmr_10 <- with(
    consumers_liver, 
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ 1,
      source_names = source_names_5,
      source_means = source_means_5,
      source_sds = source_sds_5,
      correction_means = correction_means_5,
      correction_sds = correction_sds_5,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_10_out <- cosimmr_ffvb(cosimmr_10)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 10...\n")
  saveRDS(cosimmr_10_out,"output/cosimmr_moel_10_liver_5sources.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_10_out, type = "statistics"))
  plot(cosimmr_10_out, type = "prop_histogram", obs = 1)

  # =============================================================================
  # MODELO 11: HÍGADO CON COVARIABLE (5 fuentes, detallado)
  # =============================================================================

  cat("\n--- MODELO 11: HÍGADO CON COVARIABLE (5 fuentes) ---\n\n")

  cosimmr_11 <- with(
    consumers_liver,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_5,
      source_means = source_means_5,
      source_sds = source_sds_5,
      correction_means = correction_means_5,
      correction_sds = correction_sds_5,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_11_out <- cosimmr_ffvb(cosimmr_11)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 11...\n")
  saveRDS(cosimmr_11_out,"output/cosimmr_moel_11_liver_5sources.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_11_out, type = "statistics"))
  plot(cosimmr_11, colour_by_cov = TRUE, cov_name = 'Mass_g')
  plot(cosimmr_11_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
  # MODELO 12: HÍGADO SIMPLIFICADO (5 fuentes → 2)
  # =============================================================================

  cat("\n--- MODELO 12: HÍGADO SIMPLIFICADO (5 fuentes → Marino vs Dulce) ---\n\n")

  cosimmr_12 <- with(
    consumers_liver,
    cosimmr_load(
      formula = cbind(d13C_corrected, d15N, d34S) ~ Mass_g,
      source_names = source_names_5_simple,
      source_means = source_means_5_simple,
      source_sds = source_sds_5_simple,
      correction_means = correction_means_5_simple,
      correction_sds = correction_sds_5_simple,
      concentration_means = NULL
    )
  )

  cat("Corriendo FFVB...\n")
  set.seed(26051993)
  cosimmr_12_out <- cosimmr_ffvb(cosimmr_12)
  cat("✅ Completado\n\n")

  cat("Guardando Modelo FFVB 12...\n")
  saveRDS(cosimmr_12_out,"output/cosimmr_moel_12_liver_5sources.rds")
  cat("✅ Guardado con éxito\n\n")

  print(summary(cosimmr_12_out, type = "statistics"))
  plot(cosimmr_12_out, type = 'covariates_plot', cov_name = 'Mass_g', one_plot = TRUE, n_pred = 100)

  # =============================================================================
# COMPARACIÓN: 7 FUENTES VS 5 FUENTES
# =============================================================================

cat("\n==============================================================\n")
cat("COMPARACIÓN: 7 FUENTES VS 5 FUENTES\n")
cat("==============================================================\n\n")

cat("MÚSCULO (sin covariable):\n\n")

# Proporciones con 7 fuentes
posterior_7f <- cosimmr_1_out$output$BUGSoutput$sims.list$p_mean
props_7f <- colMeans(posterior_7f)
names(props_7f) <- source_names

# Proporciones con 5 fuentes
posterior_5f <- cosimmr_7_out$output$BUGSoutput$sims.list$p_mean
props_5f <- colMeans(posterior_5f)
names(props_5f) <- source_names_5

comparison_sources <- data.frame(
  Version = c(rep("7_fuentes", 7), rep("5_fuentes", 5)),
  Source = c(source_names, source_names_5),
  Proportion = c(props_7f, props_5f)
)

print(comparison_sources)

# Gráfico comparativo
p_sources_comp <- ggplot(comparison_sources, aes(x = Source, y = Proportion, fill = Version)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("7_fuentes" = "steelblue", "5_fuentes" = "coral")) +
  labs(
    title = "Comparación: 7 fuentes vs 5 fuentes",
    subtitle = "Músculo, sin covariable",
    x = "Fuente",
    y = "Proporción",
    fill = "Versión"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sources_comp)

# =============================================================================
# ACTUALIZAR EXPORTACIÓN DE RESULTADOS
# =============================================================================

cat("\n==============================================================\n")
cat("ACTUALIZANDO RESULTADOS EXPORTADOS\n")
cat("==============================================================\n\n")

# Actualizar resumen completo
sink("output/resultados_completo_musculo_higado.txt")
cat("RESULTADOS COMPLETOS: MÚSCULO + HÍGADO (12 MODELOS)\n")
cat("====================================================\n\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("\n--- FUENTES (7) ---\n")
print(sources_all)

cat("\n--- FUENTES (5) ---\n")
print(sources_5)

cat("\n\n=== PARTE A: MÚSCULO (7 fuentes) ===\n")
cat("\nModelo 1 (sin covariable):\n")
print(summary(cosimmr_1_out, type = "statistics"))
cat("\nModelo 2 (con covariable):\n")
print(summary(cosimmr_2_out, type = "statistics"))
cat("\nModelo 3 (simplificado):\n")
print(summary(cosimmr_3_out, type = "statistics"))

cat("\n\n=== PARTE B: HÍGADO (7 fuentes) ===\n")
cat("\nModelo 4 (sin covariable):\n")
print(summary(cosimmr_4_out, type = "statistics"))
cat("\nModelo 5 (con covariable):\n")
print(summary(cosimmr_5_out, type = "statistics"))
cat("\nModelo 6 (simplificado):\n")
print(summary(cosimmr_6_out, type = "statistics"))

cat("\n\n=== PARTE C: MÚSCULO (5 fuentes) ===\n")
cat("\nModelo 7 (sin covariable):\n")
print(summary(cosimmr_7_out, type = "statistics"))
cat("\nModelo 8 (con covariable):\n")
print(summary(cosimmr_8_out, type = "statistics"))
cat("\nModelo 9 (simplificado):\n")
print(summary(cosimmr_9_out, type = "statistics"))

cat("\n\n=== PARTE C: HÍGADO (5 fuentes) ===\n")
cat("\nModelo 10 (sin covariable):\n")
print(summary(cosimmr_10_out, type = "statistics"))
cat("\nModelo 11 (con covariable):\n")
print(summary(cosimmr_11_out, type = "statistics"))
cat("\nModelo 12 (simplificado):\n")
print(summary(cosimmr_12_out, type = "statistics"))

cat("\n\n=== COMPARACIÓN MÚSCULO VS HÍGADO (7 fuentes) ===\n")
print(comparison_df)

sink()

# Guardar datos adicionales
write.csv(sources_5, "output/sources_5.csv", row.names = FALSE)
write.csv(comparison_sources, "output/comparison_7vs5_sources.csv", row.names = FALSE)

# Guardar gráfico adicional
ggsave("output/comparison_7vs5_sources.png", plot = p_sources_comp, 
       width = 12, height = 6, dpi = 300)

cat("\n✅ Archivos actualizados:\n")
cat("   - resultados_completo_musculo_higado.txt (ACTUALIZADO con 12 modelos)\n")
cat("   - sources_5.csv (NUEVO)\n")
cat("   - comparison_7vs5_sources.csv (NUEVO)\n")
cat("   - comparison_7vs5_sources.png (NUEVO)\n\n")

cat("==============================================================\n")
cat("ANÁLISIS COMPLETO FINALIZADO\n")
cat("==============================================================\n\n")

cat("RESUMEN FINAL:\n")
cat("  - 12 modelos ejecutados:\n")
cat("    * 6 con 7 fuentes (Parte A+B)\n")
cat("    * 6 con 5 fuentes (Parte C)\n")
cat("  - Tejidos analizados: Músculo + Hígado\n")
cat("  - Fuentes (7): Anchovy, Sardine, Zooplankton, Phytoplankton,\n")
cat("                 Fish, Invertebrates, Primary_producers\n")
cat("  - Fuentes (5): Anchovy, Sardine, Zooplankton, Fish, Invertebrates\n")
cat("  - Consumidores músculo:", nrow(consumers_muscle), "\n")
cat("  - Consumidores hígado:", nrow(consumers_liver), "\n")


