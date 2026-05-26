# =============================================================================
# ANÁLISIS EXPLORATORIO COMPLETO: ISÓTOPOS ESTABLES
# Sistema pelágico: Salmón, fuentes dietarias y base de la cadena trófica
# =============================================================================

library(tidyverse)
library(ggplot2)
library(corrplot)
library(vegan)
library(gridExtra)
library(patchwork)
library(RColorBrewer)

# =============================================================================
# 1. CARGAR Y PREPARAR DATOS
# =============================================================================

# Cargar datos
raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")


# Crear columna de grupo/especie
raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^h-") ~ "Salmon_liver", 
  str_detect(raw_data$Sample, "^Sar") ~ "Sardina",
  str_detect(raw_data$Sample, "^Zoo") ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito") ~ "Fitoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchoveta",
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

cat("\n=== CONTROL DE CALIDAD ===\n")

# Crear tabla de flags QC
qc_flags <- raw_data %>%
  mutate(
    flag_pct_sum = ifelse(pct_sum > 100, "suma_%>100", ""),
    flag_CN_extreme = case_when(
      group %in% c("Salmon_muscle", "Salmon_liver") & CtoN > 10 ~ "C/N_extremo",
      group %in% c("Sardina", "Anchoveta") & CtoN > 8 ~ "C/N_extremo", 
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

# Dataset limpio: MANTENER SOLO MÚSCULO DE SALMÓN + fuentes
clean_data <- qc_flags %>%
  filter(!ID %in% exclude_samples) %>%
  filter(group != "Salmon_liver")

cat("\nDataset después de QC (músculo + fuentes):", nrow(clean_data), "muestras\n")
print(table(clean_data$group))


# Función de corrección lipídica
lipid_correction_mcconnaughey <- function(d13C, CN_ratio, D = 6.0) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent <- 93 / (1 + (1/lipid_term))
  correction_factor <- -0.207 + (3.90 / (1 + (287 / L_percent)))
  correction <- D * correction_factor
  final_correction <- ifelse(CN_ratio > 3.5, correction, 0)
  return(d13C + final_correction)
}

# Función de corrección lipídica Hoffman & Sutton (2010)
lipid_correction_hoffman_sutton <- function(d13C_bulk, CN_bulk, delta_13C_lipid = -6.39, CN_protein = 3.76) {
  correction_factor <- (delta_13C_lipid * (CN_protein - CN_bulk)) / CN_bulk
  d13C_corrected <- d13C_bulk + correction_factor
  final_correction <- ifelse(CN_bulk > 3.5, d13C_corrected, d13C_bulk)
  return(final_correction)
}

# McConnaughey optimizado para Chinook salmon
lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent <- 93 / (1 + (1/lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  final_correction <- ifelse(CN_ratio > 3.5, correction_factor, 0)
  return(d13C + final_correction)
}

# Corrección lipídica para Clupeiformes (sardina/anchoveta)
lipid_correction_clupeiformes <- function(d13C_bulk, CN_bulk) {
  d13C_corrected <- d13C_bulk + 0.894 * CN_bulk - 2.377
  return(d13C_corrected)
}


# Aplicar ambas correcciones
# Aplicar correcciones específicas por grupo taxonómico
clean_data <- clean_data %>%
  mutate(
    needs_correction = CtoN > 3.5,
    d13C_corrected = case_when(
      # Salmón: parámetros optimizados para Chinook
      group == "Salmon_muscle" ~ lipid_correction_optimized_chinook(d13C, CtoN),
      
      # Sardina y Anchoveta: ecuación Clupeiformes (Sardenne et al. 2023)
      # group %in% c("Sardina", "Anchoveta") ~ lipid_correction_clupeiformes(d13C, CtoN),
      
      # Otros grupos: método general
      TRUE ~ lipid_correction_hoffman_sutton(d13C, CtoN)
    ),
    group_ordered = factor(group, levels = c("Fitoplankton", "Zooplankton", "Anchoveta", "Sardina", "Salmon_muscle"))
  )


# Datos específicos de salmones con edad
salmon_data <- clean_data %>%
  filter(group == "Salmon_muscle") %>%
  filter(!is.na(age)) %>%
  mutate(
    age_factor = factor(age, levels = c(1, 2, 3))
  )

cat("=== DATOS CARGADOS ===\n")
cat("Total de muestras:", nrow(clean_data), "\n")
print(table(clean_data$group))
cat("Salmones con edad:", nrow(salmon_data), "\n")
print(table(salmon_data$age_factor))

# =============================================================================
# 2. ESTADÍSTICAS DESCRIPTIVAS POR GRUPO
# =============================================================================

cat("\n=== ESTADÍSTICAS DESCRIPTIVAS ===\n")

# Estadísticas por grupo taxonómico
group_stats <- clean_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    # δ13C
    mean_d13C = round(mean(d13C_corrected, na.rm = TRUE), 2),
    sd_d13C = round(sd(d13C_corrected, na.rm = TRUE), 2),
    # δ15N  
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
    # δ34S
    mean_d34S = round(mean(d34S, na.rm = TRUE), 2),
    sd_d34S = round(sd(d34S, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  arrange(mean_d15N)  # Ordenar por nivel trófico

cat("ESTADÍSTICAS POR GRUPO TAXONÓMICO:\n")
print(group_stats)

# Estadísticas de salmones por edad
age_stats <- salmon_data %>%
  group_by(age_factor) %>%
  summarise(
    n = n(),
    # Morfometría
    mean_length = round(mean(fork_length, na.rm = TRUE), 1),
    sd_length = round(sd(fork_length, na.rm = TRUE), 1),
    # Isótopos
    mean_d13C = round(mean(d13C_corrected, na.rm = TRUE), 2),
    sd_d13C = round(sd(d13C_corrected, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
    mean_d34S = round(mean(d34S, na.rm = TRUE), 2),
    sd_d34S = round(sd(d34S, na.rm = TRUE), 2),
    .groups = 'drop'
  )

cat("\nESTADÍSTICAS DE SALMONES POR EDAD:\n")
print(age_stats)

# =============================================================================
# 3. VISUALIZACIÓN: ESPACIO ISOTÓPICO COMPLETO
# =============================================================================

cat("\n=== GENERANDO VISUALIZACIONES ===\n")

# Definir colores por grupo
group_colors <- c(
  "Fitoplankton" = "#2E8B57",     # Verde marino
  "Zooplankton" = "#4169E1",      # Azul real  
  "Anchoveta" = "#9370DB",        # Violeta medio
  "Sardina" = "#FF8C00",          # Naranja oscuro
  "Salmon_muscle" = "#DC143C"     # Rojo carmesí
)

age_colors <- c("1" = "#E31A1C", "2" = "#1F78B4", "3" = "#33A02C")

# Gráfico 1: Espacio isotópico δ13C vs δ15N - Todo el sistema
p1_system <- clean_data %>%
 mutate(
   # Cambiar nombre de grupo
   group_ordered = case_when(
     group_ordered == "Salmon_muscle" ~ "Chinook salmon",
     TRUE ~ group_ordered
   ),
   # Create size variable: age for salmon, fixed size for other species
   point_size = case_when(
     group_ordered == "Chinook salmon" & !is.na(age) ~ age,  # Use age directly
     group_ordered == "Chinook salmon" & is.na(age) ~ 2,     # Default for salmon without age
     TRUE ~ 2                                                # Fixed size for other species
   )
 ) %>%
 filter(!is.na(d15N) & !is.na(d13C_corrected)) %>%
 ggplot(aes(x = d13C_corrected, y = d15N, color = group_ordered)) +
 geom_point(aes(size = point_size), alpha = 0.8) +
 stat_ellipse(level = 0.95, linewidth = 1.2, alpha = 0.7) +
 scale_color_manual(values = group_colors, name = "Group") +
 scale_size_continuous(name = "Salmon age (years)", 
                      range = c(3, 7),  # Age 1=3, Age 2=5, Age 3=7
                      breaks = c(1, 2, 3),
                      labels = c("1", "2", "3")) +
 labs(
  #  title = "Isotopic Space of the Pelagic System",
  #  subtitle = "δ13C vs δ15N - Entire trophic chain (95% ellipses)",
   x = "\nδ13C corrected (‰)",
   y = "δ15N (‰)\n"
 ) +
 theme_classic() +
 theme(
   axis.text = element_text(size = 16),
   axis.title = element_text(size = 18, face = "bold"),
   legend.position = "right",
   legend.text = element_text(size = 18),
   legend.title = element_text(size = 19),
   panel.grid.minor = element_blank()
 ) +
 # Add separate legends
 guides(
   color = guide_legend(title = "Group", order = 1),
   size = guide_legend(title = "Salmon age (years)", order = 2,
                      override.aes = list(color = "red"))
 )

print(p1_system)


# OPCIÓN 3: Usar diferentes formas (shapes) además de tamaños - CORREGIDA
p1_system_4 <- clean_data %>%
  mutate(
    total_length_iso = as.numeric(total_length),
    group_ordered = case_when(
      group_ordered == "Salmon_muscle" ~ "Chinook salmon",
      TRUE ~ group_ordered
    ),
    # Crear categorías de tamaño con formas y ordenar los niveles
    size_shape = case_when(
      group_ordered == "Chinook salmon" & !is.na(total_length_iso) ~ case_when(
        total_length_iso < 50 ~ "Small (< 50 cm)",
        total_length_iso >= 50 ~ "Large (≥ 50 cm)",
        TRUE ~ "Unknown"
      ),
      TRUE ~ "Other species"
    ),
    # Convertir a factor con niveles ordenados
    size_shape = factor(size_shape, 
                       levels = c("Large (≥ 50 cm)", "Small (< 50 cm)", "Other species", "Unknown"))
  ) %>%
  filter(!is.na(d15N) & !is.na(d13C_corrected)) %>%
  ggplot(aes(x = d13C_corrected, y = d15N, color = group_ordered)) +
  geom_point(aes(size = size_shape, shape = size_shape), alpha = 0.7, stroke = 1) +
  stat_ellipse(level = 0.95, linewidth = 1.2, alpha = 0.7) +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_size_manual(name = "Salmon size",
                   values = c("Large (≥ 50 cm)" = 6, 
                             "Small (< 50 cm)" = 3,
                             "Other species" = 3, 
                             "Unknown" = 3)) +
  scale_shape_manual(name = "Salmon size",
                    values = c("Large (≥ 50 cm)" = 17,  # triángulo
                              "Small (< 50 cm)" = 16,   # círculo
                              "Other species" = 16,     # círculo
                              "Unknown" = 16)) +        # círculo
  labs(
    x = "\nδ13C corrected (‰)",
    y = "δ15N (‰)\n"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(title = "Group", order = 1),
    # Combinar size y shape en una sola leyenda
    size = "none",  # Ocultar leyenda de size separada
    shape = guide_legend(title = "Salmon size", order = 2,
                        override.aes = list(alpha = 1, 
                                          size = c(6, 3, 3)))  # Especificar tamaños manualmente
  )

print(p1_system_4)


# Gráfico 2: Espacio isotópico δ13C vs δ34S
p2_sulfur <- clean_data %>%
  filter(!is.na(d34S) & !is.na(d13C_corrected)) %>%
  ggplot(aes(x = d13C_corrected, y = d34S, color = group_ordered)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1.2, alpha = 0.7) +
  scale_color_manual(values = group_colors, name = "Grupo") +
  labs(
    title = "Espacio Isotópico: δ13C vs δ34S",
    subtitle = "Separación adicional con isótopos de azufre",
    x = "δ13C corregido (‰)",
    y = "δ34S (‰)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Gráfico 3: Salmones por edad en espacio isotópico
p3_salmon_age <- salmon_data %>%
  ggplot(aes(x = d13C_corrected, y = d15N)) +
  geom_point(aes(color = age_factor, size = fork_length), alpha = 0.8) +
  scale_color_manual(values = age_colors, name = "Edad") +
  scale_size_continuous(range = c(2, 6), name = "Longitud\n(cm)") +
  labs(
    title = "Ontogenia en el Espacio Isotópico",
    subtitle = "Salmones: tamaño = longitud, color = clase de edad",
    x = "δ13C corregido (‰)",
    y = "δ15N (‰)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.box = "horizontal"
  )

# Mostrar gráficos principales
print(p1_system)
print(p2_sulfur) 
print(p3_salmon_age)

# =============================================================================
# 4. ANÁLISIS DE POSICIÓN TRÓFICA
# =============================================================================

# Calcular posición trófica usando fitoplancton como baseline (si existe)
if("Fitoplankton" %in% clean_data$group && sum(!is.na(clean_data$d15N[clean_data$group == "Fitoplankton"])) > 0) {
  
  baseline_d15N <- mean(clean_data$d15N[clean_data$group == "Fitoplankton"], na.rm = TRUE)
  
  trophic_data <- clean_data %>%
    filter(!is.na(d15N)) %>%
    mutate(
      trophic_position = 1 + (d15N - baseline_d15N) / 3.4,  # Asumiendo discriminación de 3.4‰
      trophic_position = pmax(trophic_position, 1)  # Mínimo nivel 1
    )
  
  # Estadísticas de posición trófica
  tp_stats <- trophic_data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_TP = round(mean(trophic_position, na.rm = TRUE), 2),
      sd_TP = round(sd(trophic_position, na.rm = TRUE), 2),
      .groups = 'drop'
    ) %>%
    arrange(mean_TP)
  
  cat("\nPOSICIÓN TRÓFICA POR GRUPO:\n")
  print(tp_stats)
  
  # Gráfico de posición trófica
  p4_trophic <- trophic_data %>%
    ggplot(aes(x = reorder(group, trophic_position), y = trophic_position, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = "Posición Trófica por Grupo",
      subtitle = paste0("Baseline: Fitoplankton (δ15N = ", round(baseline_d15N, 1), "‰)"),
      x = "Grupo",
      y = "Posición Trófica"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(p4_trophic)
  
} else {
  cat("\nFitoplancton no disponible para cálculo de posición trófica\n")
  p4_trophic <- NULL
}

# =============================================================================
# 5. ANÁLISIS ONTOGENÉTICO DETALLADO
# =============================================================================

# Boxplots por edad para cada isótopo
p5_isotopes_age <- salmon_data %>%
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
  geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
  facet_wrap(~isotope_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = age_colors) +
  labs(
    title = "Cambios Ontogenéticos en Composición Isotópica",
    subtitle = "Salmones por clase de edad",
    x = "Clase de Edad",
    y = "Valor Isotópico",
    fill = "Edad"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

# Relaciones edad vs isótopos
p6_age_trends <- salmon_data %>%
  select(age, fork_length, d13C_corrected, d15N, d34S) %>%
  pivot_longer(cols = c(d13C_corrected, d15N, d34S), 
               names_to = "isotope", values_to = "value") %>%
  mutate(isotope_label = case_when(
    isotope == "d13C_corrected" ~ "δ13C (‰)",
    isotope == "d15N" ~ "δ15N (‰)",
    isotope == "d34S" ~ "δ34S (‰)"
  )) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = age, y = value)) +
  geom_point(aes(color = fork_length), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
  facet_wrap(~isotope_label, scales = "free_y", ncol = 3) +
  scale_color_viridis_c(name = "Longitud\n(cm)") +
  labs(
    title = "Tendencias Isotópicas con la Edad",
    subtitle = "Color = longitud corporal, línea = tendencia lineal",
    x = "Edad (años)",
    y = "Valor Isotópico"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  )

print(p5_isotopes_age)
print(p6_age_trends)

# =============================================================================
# 6. TESTS ESTADÍSTICOS
# =============================================================================

cat("\n=== ANÁLISIS ESTADÍSTICOS ===\n")

# Tests por edad en salmones
isotopes <- c("d13C_corrected", "d15N", "d34S")
cat("DIFERENCIAS POR EDAD (Kruskal-Wallis):\n")

for(iso in isotopes) {
  if(sum(!is.na(salmon_data[[iso]])) > 5) {
    kw_test <- kruskal.test(salmon_data[[iso]] ~ salmon_data$age_factor)
    cat(sprintf("%s: χ² = %.2f, p = %.4f", iso, kw_test$statistic, kw_test$p.value))
    if(kw_test$p.value < 0.05) cat(" ***")
    cat("\n")
  }
}

# Correlaciones edad-morfometría-isótopos
if(nrow(salmon_data) > 5) {
  vars_for_cor <- c("age", "fork_length", "d13C_corrected", "d15N")
  cor_matrix <- cor(salmon_data[vars_for_cor], use = "complete.obs", method = "spearman")
  
  cat("\nMATRIZ DE CORRELACIONES (Spearman):\n")
  print(round(cor_matrix, 3))
  
  # Gráfico de correlaciones
  corrplot(cor_matrix, method = "color", type = "upper", 
           order = "original", tl.col = "black", tl.srt = 45,
           addCoef.col = "black", number.cex = 0.8,
           title = "Correlaciones: Edad, Morfometría e Isótopos")
}

# PERMANOVA multivariado
if(nrow(salmon_data) > 10) {
  cat("\nPERMANOVA - Composición isotópica ~ Edad:\n")
  
  iso_matrix <- salmon_data %>%
    select(d13C_corrected, d15N, d34S) %>%
    na.omit()
  
  age_vector <- salmon_data$age_factor[complete.cases(salmon_data[c("d13C_corrected", "d15N")])]
  
  if(nrow(iso_matrix) > 5) {
    perm_age <- adonis2(iso_matrix ~ age_vector, permutations = 999, method = "euclidean")
    print(perm_age)
  }
}

# =============================================================================
# 7. PANEL FINAL COMBINADO
# =============================================================================

# Crear panel con los gráficos más importantes
if(!is.null(p4_trophic)) {
  panel_main <- (p1_system | p2_sulfur) / (p3_salmon_age | p4_trophic)
} else {
  panel_main <- (p1_system | p2_sulfur) / p3_salmon_age
}

panel_main <- panel_main + 
  plot_annotation(
    title = "Análisis Exploratorio: Isótopos Estables en Sistema Pelágico",
    subtitle = "Ecología trófica y cambios ontogenéticos en salmón",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(panel_main)

# Panel ontogenético
panel_ontogeny <- p5_isotopes_age / p6_age_trends
panel_ontogeny <- panel_ontogeny +
  plot_annotation(
    title = "Análisis Ontogenético Detallado",
    subtitle = "Cambios isotópicos durante el desarrollo",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold")
    )
  )

print(panel_ontogeny)

