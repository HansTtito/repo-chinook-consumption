# =============================================================================
# COMPARACIÓN COMPLETA DE MODELOS DIETARIOS
# Modelos con covariable Mass_g: 7 fuentes vs 5 fuentes
# Tejidos: Músculo + Hígado
# Incluye: Gráficos, Bondad de Ajuste, DIC, Decisión Final
# =============================================================================

library(tidyverse)
library(cosimmr)
library(ggplot2)
library(patchwork)

# =============================================================================
# 1. CARGAR MODELOS CON COVARIABLE
# =============================================================================

cat("==============================================================\n")
cat("CARGANDO MODELOS CON COVARIABLE (Mass_g)\n")
cat("==============================================================\n\n")

# MÚSCULO
cat("Cargando modelos de MÚSCULO...\n")
model_2_muscle_7f <- readRDS("output/cosimmr_model_2_muscle.rds")           # 7 fuentes
model_8_muscle_5f <- readRDS("output/cosimmr_moel_8_muscle_5sources.rds")  # 5 fuentes

# HÍGADO
cat("Cargando modelos de HÍGADO...\n")
model_5_liver_7f <- readRDS("output/cosimmr_moelel_5_liver.rds")           # 7 fuentes
model_11_liver_5f <- readRDS("output/cosimmr_moel_11_liver_5sources.rds")  # 5 fuentes

cat("✅ Modelos cargados exitosamente\n\n")

# =============================================================================
# 2. GENERAR GRÁFICOS DE COVARIABLE
# =============================================================================

cat("==============================================================\n")
cat("GENERANDO GRÁFICOS DE CAMBIO DIETARIO POR PESO\n")
cat("==============================================================\n\n")

create_diet_plot_custom <- function(model, title, n_pred = 100, 
                                   show_ribbon = TRUE, 
                                   show_means_only = FALSE) {
  
  # Crear gráfico base
  p <- plot(model, 
            type = 'covariates_plot', 
            cov_name = 'Mass_g', 
            one_plot = TRUE, 
            n_pred = n_pred)
  
  # Si solo queremos medias, reconstruir el gráfico
  if(show_means_only) {
    
    # Extraer datos del gráfico original
    plot_data <- p$data
    
    # Crear gráfico solo con líneas (sin ribbons)
    p <- ggplot(plot_data, aes(x = cov, y = mean, color = Source)) +
      geom_line(size = 1.2) +
      scale_color_discrete() +
      labs(
        x = "Mass_g",
        y = "Proportion",
        color = "Source"
      ) +
      theme_bw()
    
  }
  
  # Aplicar estilo final
  p <- p +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# CON RIBBONS (gráficos completos)
p_muscle_7f <- create_diet_plot_custom(model_2_muscle_7f, 
                                       "MÚSCULO: 7 fuentes", 
                                       show_means_only = FALSE)

p_muscle_5f <- create_diet_plot_custom(model_8_muscle_5f, 
                                       "MÚSCULO: 5 fuentes", 
                                       show_means_only = FALSE)

p_liver_7f <- create_diet_plot_custom(model_5_liver_7f, 
                                      "HÍGADO: 7 fuentes", 
                                      show_means_only = FALSE)

p_liver_5f <- create_diet_plot_custom(model_11_liver_5f, 
                                      "HÍGADO: 5 fuentes", 
                                      show_means_only = FALSE)

# SIN RIBBONS (solo líneas)
p_muscle_7f_clean <- create_diet_plot_custom(model_2_muscle_7f, 
                                             "MÚSCULO: 7 fuentes - Solo medias", 
                                             show_means_only = TRUE)

p_muscle_5f_clean <- create_diet_plot_custom(model_8_muscle_5f, 
                                             "MÚSCULO: 5 fuentes - Solo medias", 
                                             show_means_only = TRUE)

p_liver_7f_clean <- create_diet_plot_custom(model_5_liver_7f, 
                                            "HÍGADO: 7 fuentes - Solo medias", 
                                            show_means_only = TRUE)

p_liver_5f_clean <- create_diet_plot_custom(model_11_liver_5f, 
                                            "HÍGADO: 5 fuentes - Solo medias", 
                                            show_means_only = TRUE)

# =============================================================================
# 3. MOSTRAR GRÁFICOS INDIVIDUALES
# =============================================================================

cat("==============================================================\n")
cat("VISUALIZACIÓN DE GRÁFICOS\n")
cat("==============================================================\n\n")

cat("--- MÚSCULO (7 fuentes) ---\n")
print(p_muscle_7f)

cat("\n--- MÚSCULO (7 fuentes) - Solo medias ---\n")
print(p_muscle_7f_clean)

cat("\n--- MÚSCULO (5 fuentes) ---\n")
print(p_muscle_5f)

cat("\n--- MÚSCULO (5 fuentes) - Solo medias ---\n")
print(p_muscle_5f_clean)

cat("\n--- HÍGADO (7 fuentes) ---\n")
print(p_liver_7f)

cat("\n--- HÍGADO (7 fuentes) - Solo medias ---\n")
print(p_liver_7f_clean)

cat("\n--- HÍGADO (5 fuentes) ---\n")
print(p_liver_5f)

cat("\n--- HÍGADO (5 fuentes) - Solo medias ---\n")
print(p_liver_5f_clean)

# =============================================================================
# 4. PANELES COMPARATIVOS
# =============================================================================

cat("\n\nCreando paneles comparativos 2x2...\n")

# Panel con ribbons
panel_all <- (p_muscle_7f | p_muscle_5f) / (p_liver_7f | p_liver_5f) +
  plot_annotation(
    title = 'Cambios Dietarios por Peso: Comparación Completa',
    subtitle = 'Músculo vs Hígado | 7 fuentes vs 5 fuentes',
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(panel_all)

# Panel solo con medias
panel_all_clean <- (p_muscle_7f_clean | p_muscle_5f_clean) / 
                   (p_liver_7f_clean | p_liver_5f_clean) +
  plot_annotation(
    title = 'Cambios Dietarios por Peso: Comparación Completa (Solo Medias)',
    subtitle = 'Músculo vs Hígado | 7 fuentes vs 5 fuentes',
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(panel_all_clean)

# =============================================================================
# 5. BONDAD DE AJUSTE: MÉTRICAS BÁSICAS
# =============================================================================

cat("\n==============================================================\n")
cat("ANÁLISIS DE BONDAD DE AJUSTE\n")
cat("==============================================================\n\n")

cat("Extrayendo métricas de los modelos...\n\n")

# Función para extraer información del modelo FFVB
extract_model_info <- function(model, model_name, tissue, n_sources) {
  
  # Extraer información básica
  n_obs <- model$input$n_obs
  
  # Para FFVB, calculamos métricas basadas en la posterior
  posterior_samples <- model$output$BUGSoutput$sims.list
  
  # Calcular métricas de las proporciones
  p_mean <- posterior_samples$p_mean
  
  # Ancho promedio de intervalos de credibilidad (95%)
  credible_widths <- apply(p_mean, 2, function(x) {
    quantile(x, 0.975) - quantile(x, 0.025)
  })
  mean_credible_width <- mean(credible_widths)
  
  # Desviación estándar promedio de las proporciones
  mean_posterior_sd <- mean(apply(p_mean, 2, sd))
  
  # Métricas de sigma (residuales)
  sigma <- posterior_samples$sigma
  mean_sigma <- mean(apply(sigma, 2, mean))
  
  # Número efectivo de muestras
  n_posterior_samples <- nrow(p_mean)
  
  # Calcular variación entre observaciones individuales
  p_individual <- posterior_samples$p
  
  # Variabilidad individual promedio
  individual_variability <- numeric(n_obs)
  for(i in 1:n_obs) {
    p_obs <- p_individual[i, , ]
    individual_variability[i] <- mean(apply(p_obs, 2, sd))
  }
  mean_individual_var <- mean(individual_variability)
  
  # Crear data frame con resultados
  data.frame(
    Model = model_name,
    Tissue = tissue,
    N_sources = n_sources,
    N_obs = n_obs,
    N_posterior = n_posterior_samples,
    Mean_credible_width = round(mean_credible_width, 3),
    Mean_posterior_SD = round(mean_posterior_sd, 3),
    Mean_sigma = round(mean_sigma, 3),
    Mean_individual_var = round(mean_individual_var, 3),
    stringsAsFactors = FALSE
  )
}

# Extraer información de los 4 modelos
results <- bind_rows(
  extract_model_info(model_2_muscle_7f, "Modelo 2", "Músculo", 7),
  extract_model_info(model_8_muscle_5f, "Modelo 8", "Músculo", 5),
  extract_model_info(model_5_liver_7f, "Modelo 5", "Hígado", 7),
  extract_model_info(model_11_liver_5f, "Modelo 11", "Hígado", 5)
)

cat("==============================================================\n")
cat("TABLA COMPARATIVA DE BONDAD DE AJUSTE\n")
cat("==============================================================\n\n")
print(results)

# =============================================================================
# 6. INTERPRETACIÓN DE MÉTRICAS BÁSICAS
# =============================================================================

cat("\n\n==============================================================\n")
cat("INTERPRETACIÓN DE MÉTRICAS\n")
cat("==============================================================\n\n")

cat("N_posterior:\n")
cat("  - Número de muestras de la distribución posterior\n")
cat("  - FFVB genera ~3600 muestras\n\n")

cat("Mean_credible_width:\n")
cat("  - Ancho promedio de intervalos de credibilidad (95%)\n")
cat("  - Menor = estimaciones más precisas\n")
cat("  - Refleja incertidumbre en las proporciones\n\n")

cat("Mean_posterior_SD:\n")
cat("  - Desviación estándar promedio de las proporciones\n")
cat("  - Menor = estimaciones más consistentes\n\n")

cat("Mean_sigma:\n")
cat("  - Desviación estándar de los residuales\n")
cat("  - Menor = mejor ajuste del modelo a los datos\n")
cat("  - Mide error de predicción\n\n")

cat("Mean_individual_var:\n")
cat("  - Variabilidad promedio entre individuos\n")
cat("  - Mayor = más heterogeneidad en la dieta individual\n\n")

# =============================================================================
# 7. COMPARACIÓN POR TEJIDO (MÉTRICAS BÁSICAS)
# =============================================================================

cat("==============================================================\n")
cat("COMPARACIÓN 7 vs 5 FUENTES (por tejido)\n")
cat("==============================================================\n\n")

# MÚSCULO
cat("--- MÚSCULO ---\n")
muscle_comp <- results %>% filter(Tissue == "Músculo")

cat("Ancho IC (95%) con 7 fuentes:", muscle_comp$Mean_credible_width[1], "\n")
cat("Ancho IC (95%) con 5 fuentes:", muscle_comp$Mean_credible_width[2], "\n")
cat("Sigma con 7 fuentes:", muscle_comp$Mean_sigma[1], "\n")
cat("Sigma con 5 fuentes:", muscle_comp$Mean_sigma[2], "\n")

if(muscle_comp$Mean_sigma[1] < muscle_comp$Mean_sigma[2]) {
  cat("→ Modelo con 7 fuentes tiene mejor ajuste (menor sigma)\n")
} else {
  cat("→ Modelo con 5 fuentes tiene mejor ajuste (menor sigma)\n")
}

# HÍGADO
cat("\n--- HÍGADO ---\n")
liver_comp <- results %>% filter(Tissue == "Hígado")

cat("Ancho IC (95%) con 7 fuentes:", liver_comp$Mean_credible_width[1], "\n")
cat("Ancho IC (95%) con 5 fuentes:", liver_comp$Mean_credible_width[2], "\n")
cat("Sigma con 7 fuentes:", liver_comp$Mean_sigma[1], "\n")
cat("Sigma con 5 fuentes:", liver_comp$Mean_sigma[2], "\n")

if(liver_comp$Mean_sigma[1] < liver_comp$Mean_sigma[2]) {
  cat("→ Modelo con 7 fuentes tiene mejor ajuste (menor sigma)\n")
} else {
  cat("→ Modelo con 5 fuentes tiene mejor ajuste (menor sigma)\n")
}

# =============================================================================
# 8. VISUALIZACIÓN DE MÉTRICAS BÁSICAS
# =============================================================================

cat("\n\n==============================================================\n")
cat("GENERANDO GRÁFICOS DE MÉTRICAS BÁSICAS\n")
cat("==============================================================\n\n")

# Gráfico de Mean_sigma (bondad de ajuste)
p_sigma <- ggplot(results, aes(x = Model, y = Mean_sigma, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Mean_sigma, 3)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Músculo" = "darkred", "Hígado" = "darkorange")) +
  facet_wrap(~N_sources, labeller = labeller(N_sources = c("5" = "5 fuentes", "7" = "7 fuentes"))) +
  labs(
    title = "Comparación Sigma: Bondad de Ajuste",
    subtitle = "Menor sigma = mejor ajuste a los datos",
    x = "Modelo",
    y = "Mean Sigma (residuales)",
    fill = "Tejido"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0)
  )

print(p_sigma)

# Gráfico de incertidumbre (ancho de intervalos)
p_uncertainty <- ggplot(results, aes(x = Model, y = Mean_credible_width, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Mean_credible_width, 3)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Músculo" = "darkred", "Hígado" = "darkorange")) +
  labs(
    title = "Incertidumbre en las Estimaciones",
    subtitle = "Ancho promedio de intervalos de credibilidad (95%)",
    x = "Modelo",
    y = "Ancho promedio IC 95%",
    fill = "Tejido"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_uncertainty)

# Panel combinado de diagnósticos básicos
p_diagnostics_basic <- p_sigma / p_uncertainty +
  plot_annotation(
    title = "Diagnósticos Básicos de Bondad de Ajuste",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(p_diagnostics_basic)

# =============================================================================
# 9. CRITERIOS DE INFORMACIÓN (DIC) - ¡LO MÁS IMPORTANTE!
# =============================================================================

cat("\n==============================================================\n")
cat("CRITERIOS DE INFORMACIÓN: DIC\n")
cat("==============================================================\n\n")

cat("EXTRAYENDO DIC (Deviance Information Criterion)...\n\n")

# Función para extraer DIC
extract_DIC <- function(model, model_name, tissue, n_sources) {
  DIC <- model$output$BUGSoutput$DIC
  pD <- model$output$BUGSoutput$pD  # Effective number of parameters
  
  data.frame(
    Model = model_name,
    Tissue = tissue,
    N_sources = n_sources,
    DIC = round(DIC, 2),
    pD = round(pD, 2),
    stringsAsFactors = FALSE
  )
}

# Calcular DIC para todos los modelos
DIC_results <- bind_rows(
  extract_DIC(model_2_muscle_7f, "Modelo 2", "Músculo", 7),
  extract_DIC(model_8_muscle_5f, "Modelo 8", "Músculo", 5),
  extract_DIC(model_5_liver_7f, "Modelo 5", "Hígado", 7),
  extract_DIC(model_11_liver_5f, "Modelo 11", "Hígado", 5)
)

# Calcular Delta DIC (diferencia con el mejor modelo)
DIC_results <- DIC_results %>%
  mutate(Delta_DIC = round(DIC - min(DIC), 2))

cat("==============================================================\n")
cat("TABLA DE DIC\n")
cat("==============================================================\n\n")
print(DIC_results)

cat("\n\nINTERPRETACIÓN DIC:\n")
cat("  - DIC = Deviance Information Criterion\n")
cat("  - DIC más bajo = MEJOR modelo (balancea ajuste + complejidad)\n")
cat("  - pD = número efectivo de parámetros\n")
cat("  - Delta_DIC es la diferencia con el mejor modelo:\n")
cat("      * Delta_DIC < 2: modelos EQUIVALENTES → usar más SIMPLE\n")
cat("      * Delta_DIC 2-7: evidencia DÉBIL\n")
cat("      * Delta_DIC > 10: evidencia FUERTE para el mejor\n\n")

# Identificar mejor modelo global por DIC
best_by_DIC <- DIC_results %>%
  arrange(DIC) %>%
  slice(1)

cat("MEJOR MODELO GLOBAL SEGÚN DIC:\n")
cat("  →", best_by_DIC$Model, "(", best_by_DIC$Tissue, ",", 
    best_by_DIC$N_sources, "fuentes)\n")
cat("  DIC:", best_by_DIC$DIC, "\n")
cat("  pD:", best_by_DIC$pD, "parámetros efectivos\n\n")

# =============================================================================
# 10. COMPARACIÓN DIC POR TEJIDO
# =============================================================================

cat("==============================================================\n")
cat("COMPARACIÓN DIC POR TEJIDO\n")
cat("==============================================================\n\n")

# --- MÚSCULO ---
cat("--- MÚSCULO ---\n")
muscle_DIC <- DIC_results %>% filter(Tissue == "Músculo") %>% arrange(N_sources)
print(muscle_DIC[, c("Model", "N_sources", "DIC", "pD", "Delta_DIC")])

diff_muscle <- abs(muscle_DIC$DIC[1] - muscle_DIC$DIC[2])
cat("\nDiferencia DIC entre modelos:", round(diff_muscle, 2), "\n")

if(diff_muscle < 2) {
  cat("→ DECISIÓN: Modelos EQUIVALENTES\n")
  cat("   Recomendación: Usar 5 fuentes (principio de parsimonia)\n")
} else if(diff_muscle < 7) {
  winner_muscle <- muscle_DIC %>% arrange(DIC) %>% slice(1)
  cat("→ DECISIÓN: Evidencia DÉBIL para", winner_muscle$N_sources, "fuentes\n")
  cat("   Considerar contexto biológico antes de decidir\n")
} else {
  winner_muscle <- muscle_DIC %>% arrange(DIC) %>% slice(1)
  cat("→ DECISIÓN: Evidencia FUERTE para", winner_muscle$N_sources, "fuentes\n")
  cat("   DIC significativamente mejor\n")
}

# --- HÍGADO ---
cat("\n--- HÍGADO ---\n")
liver_DIC <- DIC_results %>% filter(Tissue == "Hígado") %>% arrange(N_sources)
print(liver_DIC[, c("Model", "N_sources", "DIC", "pD", "Delta_DIC")])

diff_liver <- abs(liver_DIC$DIC[1] - liver_DIC$DIC[2])
cat("\nDiferencia DIC entre modelos:", round(diff_liver, 2), "\n")

if(diff_liver < 2) {
  cat("→ DECISIÓN: Modelos EQUIVALENTES\n")
  cat("   Recomendación: Usar 5 fuentes (principio de parsimonia)\n")
} else if(diff_liver < 7) {
  winner_liver <- liver_DIC %>% arrange(DIC) %>% slice(1)
  cat("→ DECISIÓN: Evidencia DÉBIL para", winner_liver$N_sources, "fuentes\n")
  cat("   Considerar contexto biológico antes de decidir\n")
} else {
  winner_liver <- liver_DIC %>% arrange(DIC) %>% slice(1)
  cat("→ DECISIÓN: Evidencia FUERTE para", winner_liver$N_sources, "fuentes\n")
  cat("   DIC significativamente mejor\n")
}

# =============================================================================
# 11. VISUALIZACIÓN DE DIC
# =============================================================================

cat("\n\n==============================================================\n")
cat("GENERANDO GRÁFICOS DE DIC\n")
cat("==============================================================\n\n")

# Gráfico de DIC por modelo
p_DIC <- ggplot(DIC_results, aes(x = Model, y = DIC, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(DIC, 1)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Músculo" = "darkred", "Hígado" = "darkorange")) +
  geom_hline(yintercept = min(DIC_results$DIC), 
             linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 3.5, y = min(DIC_results$DIC) + 5, 
           label = "Mejor DIC", color = "blue", fontface = "bold") +
  labs(
    title = "Comparación DIC: Selección de Modelo",
    subtitle = "Menor DIC = mejor modelo (penaliza complejidad)",
    x = "Modelo",
    y = "DIC (Deviance Information Criterion)",
    fill = "Tejido"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(p_DIC)

# Gráfico de Delta DIC
p_delta_DIC <- ggplot(DIC_results, aes(x = Model, y = Delta_DIC, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Delta_DIC, 2)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Músculo" = "darkred", "Hígado" = "darkorange")) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "orange", size = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = 3.5, y = 2.5, label = "Equivalente (<2)", 
           color = "orange", fontface = "italic") +
  annotate("text", x = 3.5, y = 10.5, label = "Fuerte (>10)", 
           color = "red", fontface = "italic") +
  labs(
    title = "Delta DIC: Diferencia Respecto al Mejor Modelo",
    subtitle = "Umbral <2 = equivalentes, >10 = evidencia fuerte",
    x = "Modelo",
    y = "Delta DIC",
    fill = "Tejido"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(p_delta_DIC)

# Panel comparativo de DIC
p_DIC_panel <- p_DIC / p_delta_DIC +
  plot_annotation(
    title = "Análisis Completo de DIC",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(p_DIC_panel)

# =============================================================================
# 12. DECISIÓN FINAL INTEGRADA
# =============================================================================

cat("\n==============================================================\n")
cat("DECISIÓN FINAL: INTEGRACIÓN DE TODOS LOS CRITERIOS\n")
cat("==============================================================\n\n")

# Combinar todas las métricas
final_comparison <- results %>%
  left_join(DIC_results, by = c("Model", "Tissue", "N_sources")) %>%
  arrange(Tissue, N_sources)

cat("TABLA COMPLETA DE COMPARACIÓN:\n")
print(final_comparison)

cat("\n\nCRITERIOS DE DECISIÓN (en orden de importancia):\n")
cat("1. DIC: Balancea ajuste y complejidad (CRITERIO PRINCIPAL)\n")
cat("2. Delta_DIC: Magnitud de diferencia entre modelos\n")
cat("3. Sigma: Bondad de ajuste puro (sin penalización)\n")
cat("4. Credible Width: Incertidumbre en estimaciones\n")
cat("5. Interpretabilidad biológica\n")
cat("6. Principio de parsimonia (Navaja de Occam)\n\n")

# Crear resumen de decisión por tejido
decision_summary <- final_comparison %>%
  group_by(Tissue) %>%
  summarise(
    Best_by_DIC = N_sources[which.min(DIC)],
    Best_by_Sigma = N_sources[which.min(Mean_sigma)],
    DIC_difference = abs(diff(DIC)),
    Sigma_difference = abs(diff(Mean_sigma)),
    Recommendation = case_when(
      DIC_difference < 2 ~ "5 fuentes (equivalentes por DIC, usar parsimonia)",
      DIC_difference < 7 & Best_by_DIC == 5 ~ "5 fuentes (DIC ligeramente mejor)",
      DIC_difference < 7 & Best_by_DIC == 7 ~ "7 fuentes (DIC ligeramente mejor, evaluar contexto)",
      Best_by_DIC == 5 ~ "5 fuentes (DIC significativamente mejor)",
      TRUE ~ "7 fuentes (DIC significativamente mejor)"
    ),
    Confidence = case_when(
      DIC_difference < 2 ~ "Alta (modelos equivalentes)",
      DIC_difference < 7 ~ "Media (diferencia moderada)",
      TRUE ~ "Alta (diferencia sustancial)"
    ),
    .groups = 'drop'
  )

cat("==============================================================\n")
cat("RESUMEN DE DECISIÓN POR TEJIDO\n")
cat("==============================================================\n\n")
print(decision_summary)

cat("\n\n==============================================================\n")
cat("RECOMENDACIÓN FINAL\n")
cat("==============================================================\n\n")

for(i in 1:nrow(decision_summary)) {
  tissue <- decision_summary$Tissue[i]
  rec <- decision_summary$Recommendation[i]
  conf <- decision_summary$Confidence[i]
  dic_diff <- round(decision_summary$DIC_difference[i], 2)
  
  cat("--- ", toupper(tissue), " ---\n", sep = "")
  cat("Recomendación:", rec, "\n")
  cat("Confianza:", conf, "\n")
  cat("Diferencia DIC:", dic_diff, "\n")
  
  # Añadir contexto adicional
  tissue_data <- final_comparison %>% filter(Tissue == tissue)
  
  cat("\nDatos de soporte:\n")
  cat("  Modelo 7 fuentes: DIC =", tissue_data$DIC[tissue_data$N_sources == 7], 
      ", Sigma =", tissue_data$Mean_sigma[tissue_data$N_sources == 7], "\n")
  cat("  Modelo 5 fuentes: DIC =", tissue_data$DIC[tissue_data$N_sources == 5], 
      ", Sigma =", tissue_data$Mean_sigma[tissue_data$N_sources == 5], "\n")
  
  cat("\n")
}

cat("\n==============================================================\n")
cat("CONSIDERACIONES ADICIONALES\n")
cat("==============================================================\n\n")

cat("VENTAJAS DE 5 FUENTES:\n")
cat("  ✓ Más parsimonioso (Navaja de Occam)\n")
cat("  ✓ Más fácil de interpretar biológicamente\n")
cat("  ✓ Menor riesgo de sobreajuste\n")
cat("  ✓ Más robusto con muestras pequeñas\n")
cat("  ✓ Mejor para comunicación científica\n\n")

cat("VENTAJAS DE 7 FUENTES:\n")
cat("  ✓ Mayor detalle en resolución de dieta\n")
cat("  ✓ Puede capturar señales sutiles\n")
cat("  ✓ Mejor si fuentes son isotópicamente distinguibles\n")
cat("  ✓ Más realista si todas las fuentes son ecológicamente relevantes\n\n")

cat("DECISIÓN FINAL DEBE CONSIDERAR:\n")
cat("  1. ¿Phytoplankton y Primary_producers son isotópicamente distintos?\n")
cat("  2. ¿Salmón realmente consume estas fuentes directamente?\n")
cat("  3. ¿Las 7 fuentes aportan información única o redundante?\n")
cat("  4. ¿El objetivo es precisión máxima o interpretabilidad?\n")
cat("  5. Contexto del estudio y audiencia del reporte\n\n")

# =============================================================================
# 13. GUARDAR TODOS LOS RESULTADOS
# =============================================================================

cat("==============================================================\n")
cat("GUARDANDO TODOS LOS RESULTADOS\n")
cat("==============================================================\n\n")

# Crear carpeta de salida si no existe
if(!dir.exists("output")) dir.create("output")

# Guardar tablas
write.csv(results, "output/model_comparison_basic_metrics.csv", row.names = FALSE)
write.csv(DIC_results, "output/DIC_comparison.csv", row.names = FALSE)
write.csv(final_comparison, "output/final_model_comparison_complete.csv", row.names = FALSE)
write.csv(decision_summary, "output/decision_summary.csv", row.names = FALSE)

# Guardar gráficos de dieta
ggsave("output/diet_by_weight_muscle_7sources.png", 
       plot = p_muscle_7f, width = 10, height = 6, dpi = 300)
ggsave("output/diet_by_weight_muscle_5sources.png", 
       plot = p_muscle_5f, width = 10, height = 6, dpi = 300)
ggsave("output/diet_by_weight_liver_7sources.png", 
       plot = p_liver_7f, width = 10, height = 6, dpi = 300)
ggsave("output/diet_by_weight_liver_5sources.png", 
       plot = p_liver_5f, width = 10, height = 6, dpi = 300)

# Guardar paneles comparativos
ggsave("output/diet_by_weight_ALL_COMPARISON.png", 
       plot = panel_all, width = 16, height = 12, dpi = 300)
ggsave("output/diet_by_weight_ALL_COMPARISON_clean.png", 
       plot = panel_all_clean, width = 16, height = 12, dpi = 300)

# Guardar gráficos de métricas básicas
ggsave("output/comparison_sigma.png", 
       plot = p_sigma, width = 10, height = 6, dpi = 300)
ggsave("output/comparison_uncertainty.png", 
       plot = p_uncertainty, width = 10, height = 6, dpi = 300)
ggsave("output/comparison_diagnostics_basic.png", 
       plot = p_diagnostics_basic, width = 10, height = 10, dpi = 300)

# Guardar gráficos de DIC
ggsave("output/comparison_DIC.png", 
       plot = p_DIC, width = 10, height = 6, dpi = 300)
ggsave("output/comparison_Delta_DIC.png", 
       plot = p_delta_DIC, width = 10, height = 6, dpi = 300)
ggsave("output/comparison_DIC_panel.png", 
       plot = p_DIC_panel, width = 10, height = 12, dpi = 300)

cat("✅ Archivos guardados exitosamente:\n\n")

cat("TABLAS CSV:\n")
cat("  - model_comparison_basic_metrics.csv\n")
cat("  - DIC_comparison.csv\n")
cat("  - final_model_comparison_complete.csv\n")
cat("  - decision_summary.csv\n\n")

cat("GRÁFICOS DE DIETA:\n")
cat("  - diet_by_weight_muscle_7sources.png\n")
cat("  - diet_by_weight_muscle_5sources.png\n")
cat("  - diet_by_weight_liver_7sources.png\n")
cat("  - diet_by_weight_liver_5sources.png\n")
cat("  - diet_by_weight_ALL_COMPARISON.png\n")
cat("  - diet_by_weight_ALL_COMPARISON_clean.png\n\n")

cat("GRÁFICOS DE MÉTRICAS BÁSICAS:\n")
cat("  - comparison_sigma.png\n")
cat("  - comparison_uncertainty.png\n")
cat("  - comparison_diagnostics_basic.png\n\n")

cat("GRÁFICOS DE DIC:\n")
cat("  - comparison_DIC.png\n")
cat("  - comparison_Delta_DIC.png\n")
cat("  - comparison_DIC_panel.png\n\n")

# =============================================================================
# 14. REPORTE FINAL EN TEXTO
# =============================================================================

cat("==============================================================\n")
cat("GENERANDO REPORTE FINAL\n")
cat("==============================================================\n\n")

sink("output/REPORTE_FINAL_COMPARACION_MODELOS.txt")

cat("=============================================================================\n")
cat("REPORTE FINAL: COMPARACIÓN DE MODELOS DIETARIOS\n")
cat("Modelos Bayesianos con Covariable (Mass_g)\n")
cat("7 fuentes vs 5 fuentes | Músculo vs Hígado\n")
cat("=============================================================================\n\n")

cat("Fecha de análisis:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=============================================================================\n")
cat("1. MODELOS ANALIZADOS\n")
cat("=============================================================================\n\n")

cat("Modelo 2: Músculo, 7 fuentes, covariable Mass_g\n")
cat("Modelo 8: Músculo, 5 fuentes, covariable Mass_g\n")
cat("Modelo 5: Hígado, 7 fuentes, covariable Mass_g\n")
cat("Modelo 11: Hígado, 5 fuentes, covariable Mass_g\n\n")

cat("7 fuentes: Anchovy, Sardine, Zooplankton, Phytoplankton,\n")
cat("           Fish, Invertebrates, Primary_producers\n")
cat("5 fuentes: Anchovy, Sardine, Zooplankton, Fish, Invertebrates\n\n")

cat("=============================================================================\n")
cat("2. TABLA COMPLETA DE COMPARACIÓN\n")
cat("=============================================================================\n\n")

print(final_comparison)

cat("\n=============================================================================\n")
cat("3. ANÁLISIS DIC (CRITERIO PRINCIPAL)\n")
cat("=============================================================================\n\n")

print(DIC_results)

cat("\n--- INTERPRETACIÓN DIC ---\n\n")
cat("DIC más bajo = mejor modelo (balancea ajuste y complejidad)\n")
cat("Delta_DIC < 2: modelos equivalentes → usar más simple\n")
cat("Delta_DIC 2-7: evidencia débil\n")
cat("Delta_DIC > 10: evidencia fuerte\n\n")

cat("=============================================================================\n")
cat("4. DECISIÓN POR TEJIDO\n")
cat("=============================================================================\n\n")

print(decision_summary)

cat("\n=============================================================================\n")
cat("5. RECOMENDACIÓN FINAL\n")
cat("=============================================================================\n\n")

for(i in 1:nrow(decision_summary)) {
  tissue <- decision_summary$Tissue[i]
  rec <- decision_summary$Recommendation[i]
  conf <- decision_summary$Confidence[i]
  dic_diff <- round(decision_summary$DIC_difference[i], 2)
  
  cat("--- ", toupper(tissue), " ---\n", sep = "")
  cat("Recomendación:", rec, "\n")
  cat("Confianza:", conf, "\n")
  cat("Diferencia DIC:", dic_diff, "\n\n")
  
  tissue_data <- final_comparison %>% filter(Tissue == tissue)
  
  cat("Modelo 7 fuentes:\n")
  cat("  DIC =", tissue_data$DIC[tissue_data$N_sources == 7], "\n")
  cat("  Sigma =", tissue_data$Mean_sigma[tissue_data$N_sources == 7], "\n")
  cat("  IC width =", tissue_data$Mean_credible_width[tissue_data$N_sources == 7], "\n\n")
  
  cat("Modelo 5 fuentes:\n")
  cat("  DIC =", tissue_data$DIC[tissue_data$N_sources == 5], "\n")
  cat("  Sigma =", tissue_data$Mean_sigma[tissue_data$N_sources == 5], "\n")
  cat("  IC width =", tissue_data$Mean_credible_width[tissue_data$N_sources == 5], "\n\n")
}

cat("=============================================================================\n")
cat("6. CONSIDERACIONES PARA LA DECISIÓN FINAL\n")
cat("=============================================================================\n\n")

cat("CRITERIOS ESTADÍSTICOS:\n")
cat("  1. DIC: criterio principal (balancea ajuste y complejidad)\n")
cat("  2. Delta_DIC: magnitud de diferencia\n")
cat("  3. Sigma: bondad de ajuste puro\n")
cat("  4. Intervalos de credibilidad: incertidumbre\n\n")

cat("CRITERIOS BIOLÓGICOS:\n")
cat("  1. ¿Todas las fuentes son ecológicamente relevantes?\n")
cat("  2. ¿Las fuentes son isotópicamente distinguibles?\n")
cat("  3. ¿Hay redundancia entre fuentes?\n")
cat("  4. Interpretabilidad del modelo\n\n")

cat("PRINCIPIO DE PARSIMONIA:\n")
cat("  Si DIC < 2: Usar modelo más simple (5 fuentes)\n")
cat("  Si DIC 2-7: Evaluar contexto biológico\n")
cat("  Si DIC > 10: Usar modelo con mejor DIC\n\n")

cat("=============================================================================\n")
cat("7. ARCHIVOS GENERADOS\n")
cat("=============================================================================\n\n")

cat("TABLAS:\n")
cat("  - model_comparison_basic_metrics.csv: métricas básicas\n")
cat("  - DIC_comparison.csv: análisis DIC\n")
cat("  - final_model_comparison_complete.csv: tabla completa\n")
cat("  - decision_summary.csv: resumen de decisiones\n\n")

cat("GRÁFICOS:\n")
cat("  - diet_by_weight_*.png: cambios dietarios por peso\n")
cat("  - comparison_sigma.png: bondad de ajuste\n")
cat("  - comparison_uncertainty.png: incertidumbre\n")
cat("  - comparison_DIC.png: criterio de información\n")
cat("  - comparison_Delta_DIC.png: diferencias DIC\n\n")

cat("=============================================================================\n")
cat("FIN DEL REPORTE\n")
cat("=============================================================================\n")

sink()

cat("✅ Reporte final guardado: REPORTE_FINAL_COMPARACION_MODELOS.txt\n\n")

