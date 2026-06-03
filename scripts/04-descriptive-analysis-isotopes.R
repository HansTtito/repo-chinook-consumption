# Descriptive stable isotope analysis — Chinook salmon and dietary sources

library(tidyverse)
library(ggplot2)
library(vegan)

# ---- Cargar y preparar datos ----

raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")


raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^Sar") ~ "Sardina",
  str_detect(raw_data$Sample, "^Zoo") ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito") ~ "Fitoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchoveta",
  TRUE ~ "Unknown"
)

# Calcular suma de porcentajes
raw_data$pct_sum <- raw_data$pctN + raw_data$pctC + raw_data$pctS

print(table(raw_data$group))

# ---- Control de calidad ----

# Crear tabla de flags QC
qc_flags <- raw_data %>%
  mutate(
    flag_pct_sum = ifelse(pct_sum > 100, "suma_%>100", ""),
    flag_CN_extreme = case_when(
      group == "Salmon_muscle" & CtoN > 10 ~ "C/N_extremo",
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
  print(problems)
}

# Decisiones de exclusión
exclude_samples <- qc_flags %>%
  filter(str_detect(all_flags, "suma_%>100")) %>%
  pull(ID)

if(length(exclude_samples) > 0) {
  print(qc_flags[qc_flags$ID %in% exclude_samples, c("ID", "Sample", "group", "pct_sum", "all_flags")])
}

# Dataset limpio: músculo de salmón + fuentes
clean_data <- qc_flags %>%
  filter(!ID %in% exclude_samples)

print(table(clean_data$group))


# McConnaughey optimizado para Chinook salmon
lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent <- 93 / (1 + (1/lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  return(d13C + correction_factor)
}


clean_data <- clean_data %>%
  mutate(
    d13C_corrected = case_when(
      group == "Salmon_muscle" ~ lipid_correction_optimized_chinook(d13C, CtoN),
      TRUE ~ d13C
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

# ---- Estadísticas descriptivas ----

# Estadísticas por grupo taxonómico
group_stats <- clean_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_d13C = round(mean(d13C_corrected, na.rm = TRUE), 2),
    sd_d13C = round(sd(d13C_corrected, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
    mean_d34S = round(mean(d34S, na.rm = TRUE), 2),
    sd_d34S = round(sd(d34S, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  arrange(mean_d15N)

print(group_stats)

# Estadísticas de salmones por edad
age_stats <- salmon_data %>%
  group_by(age_factor) %>%
  summarise(
    n = n(),
    mean_length = round(mean(fork_length, na.rm = TRUE), 1),
    sd_length = round(sd(fork_length, na.rm = TRUE), 1),
    mean_d13C = round(mean(d13C_corrected, na.rm = TRUE), 2),
    sd_d13C = round(sd(d13C_corrected, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
    mean_d34S = round(mean(d34S, na.rm = TRUE), 2),
    sd_d34S = round(sd(d34S, na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(age_stats)

# ---- Visualizaciones ----

# Definir colores por grupo
group_colors <- c(
  "Fitoplankton" = "#2E8B57",     # Verde marino
  "Zooplankton" = "#4169E1",      # Azul real  
  "Anchoveta" = "#9370DB",        # Violeta medio
  "Sardina" = "#FF8C00",          # Naranja oscuro
  "Salmon_muscle" = "#DC143C"     # Rojo carmesí
)

# Gráfico 1: Espacio isotópico δ13C vs δ15N - Todo el sistema
p1_system <- clean_data %>%
 mutate(
   group_ordered = case_when(
     group_ordered == "Salmon_muscle" ~ "Chinook salmon",
     TRUE ~ group_ordered
   ),
   point_size = case_when(
     group_ordered == "Chinook salmon" & !is.na(age) ~ age,
     group_ordered == "Chinook salmon" & is.na(age) ~ 2,
     TRUE ~ 2
   )
 ) %>%
 filter(!is.na(d15N) & !is.na(d13C_corrected)) %>%
 ggplot(aes(x = d13C_corrected, y = d15N, color = group_ordered)) +
 geom_point(aes(size = point_size), alpha = 0.8) +
 stat_ellipse(level = 0.95, linewidth = 1.2, alpha = 0.7) +
 scale_color_manual(values = group_colors, name = "Group") +
 scale_size_continuous(name = "Salmon age (years)", 
                      range = c(3, 7),
                      breaks = c(1, 2, 3),
                      labels = c("1", "2", "3")) +
 labs(
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
 guides(
   color = guide_legend(title = "Group", order = 1),
   size = guide_legend(title = "Salmon age (years)", order = 2,
                      override.aes = list(color = "red"))
 )

print(p1_system)

# ---- Tests estadísticos ----

# Tests por edad en salmones
isotopes <- c("d13C_corrected", "d15N", "d34S")

for(iso in isotopes) {
  if(sum(!is.na(salmon_data[[iso]])) > 5) {
    kw_test <- kruskal.test(salmon_data[[iso]] ~ salmon_data$age_factor)
    cat(sprintf("%s: χ² = %.2f, p = %.4f", iso, kw_test$statistic, kw_test$p.value))
    if(kw_test$p.value < 0.05) cat(" ***")
    cat("\n")
  }
}

# ---- Correlaciones Spearman ----

salmon_complete <- salmon_data %>%
  select(age, fork_length, total_length, weight,
         d13C_corrected, d15N, d34S, CtoN, pctN, pctC) %>%
  mutate(total_length = as.numeric(total_length)) %>%
  filter(complete.cases(.))

cor_summary <- data.frame()
for(iso in c("d13C_corrected", "d15N", "d34S")) {
  if(sum(!is.na(salmon_complete[[iso]])) > 5) {
    test <- cor.test(salmon_complete$age, salmon_complete[[iso]], method = "spearman")
    cor_summary <- rbind(cor_summary, data.frame(
      Variable     = iso,
      Spearman_r   = round(as.numeric(test$estimate), 3),
      p_value      = round(test$p.value, 4),
      Significance = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      stringsAsFactors = FALSE
    ))
  }
}
print(cor_summary)
write.csv(cor_summary, "output/correlations_age_isotopes.csv", row.names = FALSE)

# PERMANOVA multivariado
if(nrow(salmon_data) > 10) {
  
  iso_matrix <- salmon_data %>%
    select(d13C_corrected, d15N, d34S) %>%
    na.omit()
  
  age_vector <- salmon_data$age_factor[complete.cases(salmon_data[c("d13C_corrected", "d15N")])]
  
  if(nrow(iso_matrix) > 5) {
    perm_age <- adonis2(iso_matrix ~ age_vector, permutations = 999, method = "euclidean")
    print(perm_age)
  }
}


