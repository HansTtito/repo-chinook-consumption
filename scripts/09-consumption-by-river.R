# Population-level consumption scaling — isotope-derived diet, ages 1-2

library(tidyverse)
library(ggplot2)

# ---- Cargar datos ----

poblacion_global  <- read.csv("data/data_raw/river-population/poblacion_mar_global.csv")
poblacion_por_rio <- read.csv("data/data_raw/river-population/poblacion_total_mar_por_rio.csv")

consumo_individual <- read.csv("data/data_raw/bioenergetic-model/resultados_total_consumption_by_age_ISOTOPES.csv")

# ---- Análisis de compatibilidad temporal ----

años_poblacion <- sort(unique(poblacion_global$year))
años_consumo   <- sort(unique(consumo_individual$year)) - 1
años_overlap   <- intersect(años_poblacion, años_consumo)

if (length(años_overlap) == 0) {
  stop("No hay overlap temporal entre datos de población y consumo")
}

# ---- Filtrar datos al período de overlap ----

poblacion_global_filtrada <- poblacion_global %>%
  filter(year %in% años_overlap)

poblacion_rio_filtrada <- poblacion_por_rio %>%
  filter(year %in% años_overlap)

consumo_filtrado <- consumo_individual %>%
  mutate(year_consumo = year - 1) %>%
  filter(year_consumo %in% años_overlap)

# ---- Calcular consumo poblacional GLOBAL con propagación de incertidumbre ----

consumo_poblacional_global <- consumo_filtrado %>%
  mutate(edad_consumo = age_inicial) %>%
  left_join(
    poblacion_global_filtrada,
    by = c("year_consumo" = "year", "edad_consumo" = "age")
  ) %>%
  filter(!is.na(N_mar_total)) %>%
  mutate(
    consumo_total_poblacional_t     = (consumption_kg * N_mar_total) / 1000,
    consumo_anchoveta_poblacional_t = (consumo_anchoveta_kg * N_mar_total) / 1000,
    consumo_sardina_poblacional_t   = (consumo_sardina_kg * N_mar_total) / 1000,
    consumo_otros_poblacional_t     = (consumo_otros_kg * N_mar_total) / 1000,
    consumo_se_poblacional_t        = ifelse(!is.na(consumption_se_kg),
                                             (consumption_se_kg * N_mar_total) / 1000, NA),
    consumo_lower_poblacional_t     = ifelse(!is.na(consumo_min_kg),
                                             (consumo_min_kg * N_mar_total) / 1000, NA),
    consumo_upper_poblacional_t     = ifelse(!is.na(consumo_max_kg),
                                             (consumo_max_kg * N_mar_total) / 1000, NA),
    consumo_cv_poblacional          = ifelse(!is.na(consumption_se_kg) & consumption_kg > 0,
                                             consumption_se_kg / consumption_kg, NA)
  ) %>%
  select(
    year_consumo, edad_consumo, transicion_label,
    N_mar_total,
    consumption_kg, consumption_se_kg, consumo_min_kg, consumo_max_kg,
    consumo_anchoveta_kg, consumo_sardina_kg, consumo_otros_kg,
    consumo_total_poblacional_t, consumo_se_poblacional_t,
    consumo_lower_poblacional_t, consumo_upper_poblacional_t, consumo_cv_poblacional,
    consumo_anchoveta_poblacional_t, consumo_sardina_poblacional_t, consumo_otros_poblacional_t,
    prop_anchoveta, prop_sardina, prop_otros
  )

# ---- Calcular consumo poblacional POR RÍO con propagación de incertidumbre ----

consumo_poblacional_rio <- consumo_filtrado %>%
  mutate(edad_consumo = age_inicial) %>%
  left_join(
    poblacion_rio_filtrada,
    by = c("year_consumo" = "year", "edad_consumo" = "age")
  ) %>%
  filter(!is.na(N_mar)) %>%
  mutate(
    consumo_total_poblacional_t     = (consumption_kg * N_mar) / 1000,
    consumo_anchoveta_poblacional_t = (consumo_anchoveta_kg * N_mar) / 1000,
    consumo_sardina_poblacional_t   = (consumo_sardina_kg * N_mar) / 1000,
    consumo_otros_poblacional_t     = (consumo_otros_kg * N_mar) / 1000,
    consumo_se_poblacional_t        = ifelse(!is.na(consumption_se_kg),
                                             (consumption_se_kg * N_mar) / 1000, NA),
    consumo_lower_poblacional_t     = ifelse(!is.na(consumo_min_kg),
                                             (consumo_min_kg * N_mar) / 1000, NA),
    consumo_upper_poblacional_t     = ifelse(!is.na(consumo_max_kg),
                                             (consumo_max_kg * N_mar) / 1000, NA),
    consumo_cv_poblacional          = ifelse(!is.na(consumption_se_kg) & consumption_kg > 0,
                                             consumption_se_kg / consumption_kg, NA)
  ) %>%
  select(
    year_consumo, Rio, edad_consumo, transicion_label,
    N_mar,
    consumption_kg, consumption_se_kg, consumo_min_kg, consumo_max_kg,
    consumo_anchoveta_kg, consumo_sardina_kg, consumo_otros_kg,
    consumo_total_poblacional_t, consumo_se_poblacional_t,
    consumo_lower_poblacional_t, consumo_upper_poblacional_t, consumo_cv_poblacional,
    consumo_anchoveta_poblacional_t, consumo_sardina_poblacional_t, consumo_otros_poblacional_t,
    prop_anchoveta, prop_sardina, prop_otros
  )

# ---- Crear resúmenes anuales con propagación de incertidumbre ----

resumen_anual_global <- consumo_poblacional_global %>%
  group_by(year_consumo) %>%
  summarise(
    poblacion_total       = sum(N_mar_total, na.rm = TRUE),
    consumo_total_t       = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_t   = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_t     = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_t       = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    n_edades              = n_distinct(edad_consumo),
    consumo_se_total_t    = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    consumo_lower_total_t = sum(consumo_lower_poblacional_t, na.rm = TRUE),
    consumo_upper_total_t = sum(consumo_upper_poblacional_t, na.rm = TRUE),
    consumo_cv_total      = ifelse(consumo_total_t > 0, consumo_se_total_t / consumo_total_t, NA),
    consumo_ci95_lower_t  = consumo_total_t - (1.96 * consumo_se_total_t),
    consumo_ci95_upper_t  = consumo_total_t + (1.96 * consumo_se_total_t),
    n_casos_con_se        = sum(!is.na(consumo_se_poblacional_t)),
    .groups = 'drop'
  ) %>%
  arrange(year_consumo)

resumen_anual_rio <- consumo_poblacional_rio %>%
  group_by(year_consumo, Rio) %>%
  summarise(
    poblacion_total       = sum(N_mar, na.rm = TRUE),
    consumo_total_t       = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_t   = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_t     = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_t       = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    n_edades              = n_distinct(edad_consumo),
    consumo_se_total_t    = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    consumo_lower_total_t = sum(consumo_lower_poblacional_t, na.rm = TRUE),
    consumo_upper_total_t = sum(consumo_upper_poblacional_t, na.rm = TRUE),
    consumo_cv_total      = ifelse(consumo_total_t > 0, consumo_se_total_t / consumo_total_t, NA),
    consumo_ci95_lower_t  = consumo_total_t - (1.96 * consumo_se_total_t),
    consumo_ci95_upper_t  = consumo_total_t + (1.96 * consumo_se_total_t),
    n_casos_con_se        = sum(!is.na(consumo_se_poblacional_t)),
    .groups = 'drop'
  ) %>%
  arrange(year_consumo, Rio)

resumen_por_edad_global <- consumo_poblacional_global %>%
  group_by(edad_consumo, transicion_label) %>%
  summarise(
    años_disponibles                  = n_distinct(year_consumo),
    poblacion_promedio                = mean(N_mar_total, na.rm = TRUE),
    consumo_individual_promedio_kg    = mean(consumption_kg, na.rm = TRUE),
    consumo_individual_se_promedio_kg = mean(consumption_se_kg, na.rm = TRUE),
    consumo_total_promedio_t          = mean(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_total_se_promedio_t       = mean(consumo_se_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_promedio_t      = mean(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_promedio_t        = mean(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_promedio_t          = mean(consumo_otros_poblacional_t, na.rm = TRUE),
    consumo_cv_promedio               = mean(consumo_cv_poblacional, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(edad_consumo)

# ---- Mostrar resultados clave ----

print("\n=== RESUMEN ANUAL GLOBAL ===")
print(resumen_anual_global)

print("\n=== RESUMEN POR EDAD (GLOBAL) ===")
print(resumen_por_edad_global)

print("\n=== RESUMEN ANUAL POR RÍO (primeros registros) ===")
print(head(resumen_anual_rio, 10))

# ---- Crear gráficos ----

colores_edad_transicion <- c("1 \u2192 2" = "grey85", "2 \u2192 3" = "grey65",
                              "3 \u2192 4" = "grey45", "4 \u2192 5" = "grey25")
colores_presa <- c("Anchoveta" = "#E76F51", "Sardina" = "#4472C4", "Otros" = "#8D6E63")
colores_rio   <- c("Bueno" = "#1f77b4", "Tolten" = "#ff7f0e",
                   "Imperial" = "#2ca02c", "Valdivia" = "#d62728")

p1 <- ggplot(resumen_anual_global, aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "steelblue", alpha = 0.8, color = "white") +
  geom_errorbar(aes(ymin = pmax(0, consumo_ci95_lower_t), ymax = consumo_ci95_upper_t),
                width = 0.3, color = "darkred", size = 1, alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title   = "Annual population-level consumption",
    x       = "",
    y = "Annual consumption (t)\n",
    caption = "Barras de error: IC 95% | Ages 1-2 diet: muscle isotope model, 5 sources"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text   = element_text(color = "black", size = 18),
    axis.title  = element_text(size = 20)
  ) +
  geom_text(aes(label = paste0(scales::comma(round(consumo_total_t, 0)),
                               "\n(CV:", round(consumo_cv_total * 100, 1), "%)")),
            vjust = -2, size = 3.5)

print(p1)
ggsave("output/consumo_poblacional_total_ISOTOPES.png", plot = p1, width = 14, height = 10, dpi = 300)

p2 <- consumo_poblacional_global %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_poblacional_t,
             fill = factor(transicion_label))) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_edad_transicion, name = "Edad") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title   = "Annual consumption by age transition",
    x       = "",
    y = "Annual consumption (t)\n",
    caption = "Ages 1-2 diet: muscle isotope model, 5 sources"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text   = element_text(color = "black", size = 18),
    axis.title  = element_text(size = 20),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p2)
ggsave("output/consumo_poblacional_por_edad_ISOTOPES.png", plot = p2, width = 14, height = 10, dpi = 300)

consumo_presa_anual <- resumen_anual_global %>%
  select(year_consumo, consumo_anchoveta_t, consumo_sardina_t, consumo_otros_t) %>%
  pivot_longer(cols = starts_with("consumo_"), names_to = "presa", values_to = "consumo_t") %>%
  mutate(presa = case_when(
    presa == "consumo_anchoveta_t" ~ "Anchoveta",
    presa == "consumo_sardina_t"   ~ "Sardina",
    presa == "consumo_otros_t"     ~ "Otros"
  ))

p3 <- ggplot(consumo_presa_anual, aes(x = factor(year_consumo), y = consumo_t, fill = presa)) +
  geom_col(alpha = 0.9, color = "white", size = 0.3) +
  scale_fill_manual(values = colores_presa, name = "Tipo de presa") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title   = "Annual consumption by prey type",
    x       = "",
    y = "Annual consumption (t)\n",
    caption = "Ages 1-2 diet: muscle isotope model, 5 sources"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(angle = 90, hjust = 1),
    legend.position = "bottom",
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text    = element_text(color = "black", size = 18),
    axis.title   = element_text(size = 20),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p3)
ggsave("output/consumo_poblacional_por_presa_ISOTOPES.png", plot = p3, width = 14, height = 10, dpi = 300)

p4 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t, fill = Rio)) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_rio, name = "Río") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title   = "Annual consumption by river",
    x       = "",
    y = "Annual consumption (t)\n",
    caption = "Ages 1-2 diet: muscle isotope model, 5 sources"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text   = element_text(color = "black", size = 18),
    axis.title  = element_text(size = 20),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p4)
ggsave("output/consumo_poblacional_por_rio_ISOTOPES.png", plot = p4, width = 14, height = 10, dpi = 300)

p5 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "darkgreen", alpha = 0.8, color = "white") +
  facet_wrap(~ Rio, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title   = "Temporal consumption trends by river",
    x       = "",
    y = "Annual consumption (t)\n",
    caption = "Ages 1-2 diet: muscle isotope model, 5 sources"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text   = element_text(color = "black", size = 18),
    axis.title  = element_text(size = 20),
    strip.text  = element_text(size = 18)
  )

print(p5)
ggsave("output/consumo_poblacional_facetado_ISOTOPES.png", plot = p5, width = 14, height = 10, dpi = 300)

# ---- Guardar resultados ----

dir.create("data/data_raw/population-consumption", showWarnings = FALSE, recursive = TRUE)

write.csv(consumo_poblacional_global, "data/data_raw/population-consumption/consumo_poblacional_global_detallado_ISOTOPES.csv", row.names = FALSE)
write.csv(consumo_poblacional_rio,    "data/data_raw/population-consumption/consumo_poblacional_por_rio_detallado_ISOTOPES.csv",    row.names = FALSE)
write.csv(resumen_anual_global,       "data/data_raw/population-consumption/resumen_consumo_anual_global_ISOTOPES.csv",              row.names = FALSE)
write.csv(resumen_anual_rio,          "data/data_raw/population-consumption/resumen_consumo_anual_por_rio_ISOTOPES.csv",             row.names = FALSE)
write.csv(resumen_por_edad_global,    "data/data_raw/population-consumption/resumen_consumo_por_edad_global_ISOTOPES.csv",           row.names = FALSE)

tabla_consolidada <- resumen_anual_global %>%
  mutate(
    consumo_total_miles_t       = round(consumo_total_t / 1000, 2),
    consumo_se_miles_t          = round(consumo_se_total_t / 1000, 2),
    consumo_cv_pct              = round(consumo_cv_total * 100, 1),
    consumo_ci95_lower_miles_t  = round(pmax(0, consumo_ci95_lower_t / 1000), 2),
    consumo_ci95_upper_miles_t  = round(consumo_ci95_upper_t / 1000, 2),
    consumo_anchoveta_miles_t   = round(consumo_anchoveta_t / 1000, 2),
    consumo_sardina_miles_t     = round(consumo_sardina_t / 1000, 2),
    consumo_otros_miles_t       = round(consumo_otros_t / 1000, 2)
  ) %>%
  select(year_consumo, poblacion_total,
         consumo_total_miles_t, consumo_se_miles_t, consumo_cv_pct,
         consumo_ci95_lower_miles_t, consumo_ci95_upper_miles_t,
         consumo_anchoveta_miles_t, consumo_sardina_miles_t, consumo_otros_miles_t,
         n_casos_con_se)

write.csv(tabla_consolidada, "data/data_raw/population-consumption/tabla_consolidada_con_incertidumbre_ISOTOPES.csv", row.names = FALSE)

# ---- Análisis adicional de incertidumbre ----

cv_stats <- consumo_poblacional_global %>%
  filter(!is.na(consumo_cv_poblacional)) %>%
  group_by(edad_consumo) %>%
  summarise(
    cv_promedio_pct  = round(mean(consumo_cv_poblacional) * 100, 1),
    cv_mediano_pct   = round(median(consumo_cv_poblacional) * 100, 1),
    cv_min_pct       = round(min(consumo_cv_poblacional) * 100, 1),
    cv_max_pct       = round(max(consumo_cv_poblacional) * 100, 1),
    n_observaciones  = n(),
    .groups = 'drop'
  ) %>%
  arrange(edad_consumo)

print("Estadísticas de Coeficiente de Variación por Edad:")
print(cv_stats)
