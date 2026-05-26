library(janitor)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

chinook_data <- read.csv(file = "data/data_raw/biological-data/data_chinook_cleaned.csv")

# ANÁLISIS CORRECTO - Incluir todos los individuos con datos de contenido estomacal (incluyendo 0s)
chinook_diet_corrected <- chinook_data |>
  clean_names() |>
  dplyr::select(region, location, species, year, month,
                day, fl_mm, hl_mm, tl_mm, tw_g, sex, age,
                stomach_weight_g, total_content_g, sardina_g, anchoveta_g,
                teleost_digested_g, bolus_g, p_sc, p_an, p_tni, p_rni, p_is,
                p_ca) |>
  dplyr::filter(species == "Oncorhynchus tshawytscha" | is.na(species), 
                !is.na(age),
                age >= 2 & age <= 5,  # Solo edades 2-5
                # Incluir individuos que tengan al menos UNA medida de contenido estomacal (no NA)
                (!is.na(sardina_g) | !is.na(anchoveta_g) | !is.na(teleost_digested_g) | 
                 !is.na(bolus_g) | !is.na(p_sc) | !is.na(p_an) | 
                 !is.na(p_tni) | !is.na(p_is) | !is.na(p_ca))) |>
  # Crear las categorías de alimento (0 cuando es NA)
  mutate(
    sardina_total = coalesce(sardina_g, 0) + coalesce(p_sc, 0),
    anchoveta_total = coalesce(anchoveta_g, 0) + coalesce(p_an, 0),
    otros_total = coalesce(teleost_digested_g, 0) + coalesce(bolus_g, 0) +
                  coalesce(p_tni, 0) + coalesce(p_is, 0) + coalesce(p_ca, 0),
    # Calcular peso total por individuo
    peso_total_individuo = sardina_total + anchoveta_total + otros_total
  ) |>
  # Solo filtrar individuos que tengan peso_total > 0 (para evitar divisiones por 0)
  filter(peso_total_individuo > 0) |>
  mutate(
    # Calcular proporciones por individuo
    prop_sardina = sardina_total / peso_total_individuo,
    prop_anchoveta = anchoveta_total / peso_total_individuo,
    prop_otros = otros_total / peso_total_individuo
  )

# Número de individuos únicos
cat("Número de individuos únicos con datos de dieta:", nrow(chinook_diet_corrected), "\n")

# ANÁLISIS POR AÑO Y EDAD - Proporción basada en pesos totales
chinook_diet_year_age <- chinook_diet_corrected |>
  group_by(year, age) |>
  summarise(
    n_individuos = n(),
    sardina_peso_total = sum(sardina_total, na.rm = TRUE),
    anchoveta_peso_total = sum(anchoveta_total, na.rm = TRUE),
    otros_peso_total = sum(otros_total, na.rm = TRUE),
    peso_total_edad = sum(peso_total_individuo, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    prop_sardina_mean = sardina_peso_total / peso_total_edad,
    prop_anchoveta_mean = anchoveta_peso_total / peso_total_edad,
    prop_otros_mean = otros_peso_total / peso_total_edad,
    sardina_porcentaje = round(prop_sardina_mean * 100, 2),
    anchoveta_porcentaje = round(prop_anchoveta_mean * 100, 2),
    otros_porcentaje = round(prop_otros_mean * 100, 2),
    # Promedios por individuo para referencia
    sardina_peso_mean = sardina_peso_total / n_individuos,
    anchoveta_peso_mean = anchoveta_peso_total / n_individuos,
    otros_peso_mean = otros_peso_total / n_individuos
  )

# ANÁLISIS SOLO POR EDAD - Proporción basada en pesos totales
chinook_diet_age <- chinook_diet_corrected |>
  group_by(age) |>
  summarise(
    n_individuos = n(),
    sardina_peso_total = sum(sardina_total, na.rm = TRUE),
    anchoveta_peso_total = sum(anchoveta_total, na.rm = TRUE),
    otros_peso_total = sum(otros_total, na.rm = TRUE),
    peso_total_edad = sum(peso_total_individuo, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    prop_sardina_mean = sardina_peso_total / peso_total_edad,
    prop_anchoveta_mean = anchoveta_peso_total / peso_total_edad,
    prop_otros_mean = otros_peso_total / peso_total_edad,
    sardina_porcentaje = round(prop_sardina_mean * 100, 2),
    anchoveta_porcentaje = round(prop_anchoveta_mean * 100, 2),
    otros_porcentaje = round(prop_otros_mean * 100, 2),
    # Mantener promedios individuales para referencia
    sardina_peso_mean = sardina_peso_total / n_individuos,
    anchoveta_peso_mean = anchoveta_peso_total / n_individuos,
    otros_peso_mean = otros_peso_total / n_individuos
  )

# FORMATO LARGO PARA VISUALIZACIÓN
chinook_diet_long <- chinook_diet_age |>
  pivot_longer(
    cols = c(sardina_porcentaje, anchoveta_porcentaje, otros_porcentaje),
    names_to = "tipo_alimento",
    values_to = "porcentaje"
  ) |>
  mutate(tipo_alimento = case_when(
    tipo_alimento == "sardina_porcentaje" ~ "Sardina",
    tipo_alimento == "anchoveta_porcentaje" ~ "Anchoveta",
    tipo_alimento == "otros_porcentaje" ~ "Otros"
  ))

# Añadir datos de edad 1 (asumiendo 100% otros)
age_1 <- data.frame(
  age = 1, 
  n_individuos = 1,
  sardina_peso_total = 0,
  anchoveta_peso_total = 0,
  otros_peso_total = 1,
  peso_total_edad = 1,
  prop_sardina_mean = 0,
  prop_anchoveta_mean = 0,
  prop_otros_mean = 1,
  sardina_porcentaje = 0,
  anchoveta_porcentaje = 0,
  otros_porcentaje = 100,
  sardina_peso_mean = 0,
  anchoveta_peso_mean = 0,
  otros_peso_mean = 1
)

chinook_diet_age_final <- rbind(age_1, chinook_diet_age)

# Mostrar resultados
print("=== ANÁLISIS POR AÑO Y EDAD ===")
print(chinook_diet_year_age |> as.data.frame())

print("\n=== ANÁLISIS POR EDAD (incluyendo edad 1) ===")
print(chinook_diet_age_final)

print("\n=== FORMATO LARGO PARA VISUALIZACIÓN ===")
print(chinook_diet_long)

# Guardar archivos
write.csv(chinook_diet_age_final, "data/data_raw/biological-data/diet_proportion_by_age.csv", row.names = FALSE)

# Crear tabla resumen por edad
diet_summary_table <- chinook_diet_age_final %>%
  select(age, n_individuos, sardina_porcentaje, anchoveta_porcentaje, otros_porcentaje) %>%
  rename(
    Age = age,
    n = n_individuos,
    `Sardine (%)` = sardina_porcentaje,
    `Anchovy (%)` = anchoveta_porcentaje,
    `Others (%)` = otros_porcentaje
  ) %>%
  arrange(Age)

# Mostrar la tabla
print("=== TABLA RESUMEN ===")
print(diet_summary_table)

# Datos de peso por edad
diet_with_fish_weight <- chinook_diet_corrected %>%
  group_by(age) %>%
  summarise(
    mean_fish_weight = mean(tw_g, na.rm = TRUE),
    se_fish_weight = sd(tw_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

weight_age_1 <- data.frame(age = 1, mean_fish_weight = 700, se_fish_weight = NA)
diet_with_fish_weight_final <- rbind(weight_age_1, diet_with_fish_weight)

# Combinar con tabla de dieta
diet_table_complete <- diet_summary_table %>%
  left_join(diet_with_fish_weight_final, by = c("Age" = "age")) %>%
  mutate(
    `Mean fish weight (g)` = ifelse(is.na(se_fish_weight), 
                                   sprintf("%.1f", mean_fish_weight),
                                   sprintf("%.1f ± %.1f", mean_fish_weight, se_fish_weight))
  ) %>%
  select(-mean_fish_weight, -se_fish_weight)

write.csv(diet_table_complete, "data/data_raw/biological-data/diet_weight_age.csv", row.names = FALSE)

# GRÁFICOS

# Preparar datos para gráficos - traducir al inglés
chinook_diet_age_eng <- chinook_diet_long %>%
  mutate(tipo_alimento = case_when(
    tipo_alimento == "Sardina" ~ "Sardine",
    tipo_alimento == "Anchoveta" ~ "Anchovy",
    tipo_alimento == "Otros" ~ "Others"
  ))

# Crear edad 1 en formato largo con TODAS las columnas necesarias
age_1_long <- data.frame(
  age = rep(1, 3),
  n_individuos = rep(1, 3),
  sardina_peso_total = rep(0, 3),
  anchoveta_peso_total = rep(0, 3),
  otros_peso_total = rep(1, 3),
  peso_total_edad = rep(1, 3),
  prop_sardina_mean = rep(0, 3),
  prop_anchoveta_mean = rep(0, 3),
  prop_otros_mean = rep(1, 3),
  # sardina_porcentaje = c(0, 0, 100),
  # anchoveta_porcentaje = c(0, 0, 100),
  # otros_porcentaje = c(0, 0, 100),
  sardina_peso_mean = rep(0, 3),
  anchoveta_peso_mean = rep(0, 3),
  otros_peso_mean = rep(1, 3),
  tipo_alimento = c("Sardine", "Anchovy", "Others"),
  porcentaje = c(0, 0, 100)
)

# Corregir los porcentajes para edad 1
# age_1_long$porcentaje <- c(0, 0, 100)

chinook_diet_age_eng <- rbind(age_1_long, chinook_diet_age_eng)

# Paleta de colores
diet_colors <- c("Sardine" = "#4472C4", "Anchovy" = "#E76F51", "Others" = "#8D6E63")

# 1. Gráfico principal de composición por edad
p_diet_age <- ggplot(chinook_diet_age_eng, aes(x = factor(age), y = porcentaje, fill = tipo_alimento)) +
  geom_col(alpha = 0.9, color = "white", size = 0.3) +
  scale_fill_manual(values = diet_colors, name = "Diet") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Age",
    y = "% by weight in the diet"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  geom_text(aes(label = ifelse(porcentaje > 5, paste0(round(porcentaje, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 5.5, fontface = "bold", color = "white")

print(p_diet_age)

# 2. Análisis temporal (solo si hay datos de múltiples años)
if(nrow(chinook_diet_year_age) > 0) {
  variacion_anual <- chinook_diet_year_age %>%
    group_by(year) %>%
    summarise(
      sardina_prom = mean(sardina_porcentaje, na.rm = TRUE),
      anchoveta_prom = mean(anchoveta_porcentaje, na.rm = TRUE),
      otros_prom = mean(otros_porcentaje, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(sardina_prom, anchoveta_prom, otros_prom),
      names_to = "tipo_alimento",
      values_to = "porcentaje_promedio"
    ) %>%
    mutate(tipo_alimento = case_when(
      tipo_alimento == "sardina_prom" ~ "Sardine",
      tipo_alimento == "anchoveta_prom" ~ "Anchovy",
      tipo_alimento == "otros_prom" ~ "Others"
    ))

  p_temporal <- ggplot(variacion_anual, aes(x = year, y = porcentaje_promedio, color = tipo_alimento)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = diet_colors, name = "Diet type") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = unique(variacion_anual$year)) +
    labs(
      title = "Temporal Evolution of Diet",
      subtitle = "Average percentage per year (all ages)",
      x = "Year",
      y = "Average percentage (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  print(p_temporal)
}

# Resumen final
cat("\nRESUMEN DE COMPOSICIÓN POR EDAD:\n")
print(diet_table_complete)
