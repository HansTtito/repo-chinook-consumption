# ================================================================
# EXTRAPOLACIÓN DE CONSUMO INDIVIDUAL A POBLACIONAL - ISOTOPES
# Calcula consumo total y por presas a nivel de río y global
# Usando dieta con isótopos para edades 1-2
# ================================================================

library(tidyverse)
library(ggplot2)

# ---- Cargar datos ----

print("=== CARGANDO DATOS ===")

# Datos de población en el mar
poblacion_global <- read.csv("data/data_raw/river-population/poblacion_mar_global.csv")
poblacion_por_rio <- read.csv("data/data_raw/river-population/poblacion_total_mar_por_rio.csv")

# Datos de consumo individual CON ISÓTOPOS
consumo_individual <- read.csv("data/data_raw/bioenergetic-model/resultados_total_consumption_by_age_ISOTOPES.csv")

print(paste("Población global:", nrow(poblacion_global), "registros"))
print(paste("Población por río:", nrow(poblacion_por_rio), "registros"))
print(paste("Consumo individual (ISOTOPES):", nrow(consumo_individual), "registros"))

# ---- Análisis de compatibilidad temporal ----

años_poblacion <- sort(unique(poblacion_global$year))
años_consumo <- sort(unique(consumo_individual$year)) - 1
años_overlap <- intersect(años_poblacion, años_consumo)

print(paste("Años población:", paste(range(años_poblacion), collapse = "-")))
print(paste("Años consumo:", paste(range(años_consumo), collapse = "-")))
print(paste("Overlap disponible:", paste(años_overlap, collapse = ", ")))

if(length(años_overlap) == 0) {
  stop("No hay overlap temporal entre datos de población y consumo")
}

# ---- Filtrar datos al período de overlap ----

poblacion_global_filtrada <- poblacion_global %>%
  filter(year %in% años_overlap)

poblacion_rio_filtrada <- poblacion_por_rio %>%
  filter(year %in% años_overlap)

consumo_filtrado <- consumo_individual %>%
  mutate(year_consumo = year - 1) %>%   # año real de consumo
  filter(year_consumo %in% años_overlap)

print(paste("Datos filtrados - Población global:", nrow(poblacion_global_filtrada)))
print(paste("Datos filtrados - Población por río:", nrow(poblacion_rio_filtrada)))
print(paste("Datos filtrados - Consumo:", nrow(consumo_filtrado)))

# ---- Función para mapear transición a edad ----

mapear_transicion_a_edad <- function(age_inicial) {
  # Transición X→Y corresponde al consumo de individuos de edad X
  return(age_inicial)
}

# ---- Calcular consumo poblacional GLOBAL con propagación de incertidumbre ----

print("\n=== CALCULANDO CONSUMO POBLACIONAL GLOBAL CON INCERTIDUMBRE ===")

consumo_poblacional_global <- consumo_filtrado %>%
  mutate(edad_consumo = mapear_transicion_a_edad(age_inicial)) %>%
  left_join(
    poblacion_global_filtrada,
    by = c("year_consumo" = "year", "edad_consumo" = "age")
  ) %>%
  filter(!is.na(N_mar_total)) %>%  # Eliminar casos sin población
  mutate(
    # Consumo total poblacional (kg → toneladas)
    consumo_total_poblacional_t = (consumption_kg * N_mar_total) / 1000,
    
    # Consumo poblacional por presa (kg → toneladas)
    consumo_anchoveta_poblacional_t = (consumo_anchoveta_kg * N_mar_total) / 1000,
    consumo_sardina_poblacional_t = (consumo_sardina_kg * N_mar_total) / 1000,
    consumo_otros_poblacional_t = (consumo_otros_kg * N_mar_total) / 1000,
    
    # PROPAGACIÓN DE INCERTIDUMBRE
    # Error estándar poblacional (kg → toneladas)
    consumo_se_poblacional_t = ifelse(!is.na(consumption_se_kg), 
                                      (consumption_se_kg * N_mar_total) / 1000, 
                                      NA),
    
    # Intervalos de confianza poblacionales (kg → toneladas)
    consumo_lower_poblacional_t = ifelse(!is.na(consumo_min_kg), 
                                         (consumo_min_kg * N_mar_total) / 1000, 
                                         NA),
    consumo_upper_poblacional_t = ifelse(!is.na(consumo_max_kg), 
                                         (consumo_max_kg * N_mar_total) / 1000, 
                                         NA),
    
    # Coeficiente de variación poblacional (igual al individual)
    consumo_cv_poblacional = ifelse(!is.na(consumption_se_kg) & consumption_kg > 0,
                                    consumption_se_kg / consumption_kg,
                                    NA)
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

print(paste("Registros de consumo poblacional global:", nrow(consumo_poblacional_global)))

# ---- Calcular consumo poblacional POR RÍO con propagación de incertidumbre ----

print("\n=== CALCULANDO CONSUMO POBLACIONAL POR RÍO CON INCERTIDUMBRE ===")

consumo_poblacional_rio <- consumo_filtrado %>%
  mutate(edad_consumo = mapear_transicion_a_edad(age_inicial)) %>%
  left_join(
    poblacion_rio_filtrada,
    by = c("year_consumo" = "year", "edad_consumo" = "age")
  ) %>%
  filter(!is.na(N_mar)) %>%  # Eliminar casos sin población
  mutate(
    # Consumo total poblacional por río (kg → toneladas)
    consumo_total_poblacional_t = (consumption_kg * N_mar) / 1000,
    
    # Consumo poblacional por presa y río (kg → toneladas)
    consumo_anchoveta_poblacional_t = (consumo_anchoveta_kg * N_mar) / 1000,
    consumo_sardina_poblacional_t = (consumo_sardina_kg * N_mar) / 1000,
    consumo_otros_poblacional_t = (consumo_otros_kg * N_mar) / 1000,
    
    # PROPAGACIÓN DE INCERTIDUMBRE POR RÍO
    # Error estándar poblacional por río (kg → toneladas)
    consumo_se_poblacional_t = ifelse(!is.na(consumption_se_kg), 
                                      (consumption_se_kg * N_mar) / 1000, 
                                      NA),
    
    # Intervalos de confianza poblacionales por río (kg → toneladas)
    consumo_lower_poblacional_t = ifelse(!is.na(consumo_min_kg), 
                                         (consumo_min_kg * N_mar) / 1000, 
                                         NA),
    consumo_upper_poblacional_t = ifelse(!is.na(consumo_max_kg), 
                                         (consumo_max_kg * N_mar) / 1000, 
                                         NA),
    
    # Coeficiente de variación poblacional por río
    consumo_cv_poblacional = ifelse(!is.na(consumption_se_kg) & consumption_kg > 0,
                                    consumption_se_kg / consumption_kg,
                                    NA)
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

print(paste("Registros de consumo poblacional por río:", nrow(consumo_poblacional_rio)))

# ---- Crear resúmenes anuales con propagación de incertidumbre ----

print("\n=== CREANDO RESÚMENES CON INCERTIDUMBRE ===")

# Resumen anual global (suma de todas las edades) CON PROPAGACIÓN
resumen_anual_global <- consumo_poblacional_global %>%
  group_by(year_consumo) %>%
  summarise(
    poblacion_total = sum(N_mar_total, na.rm = TRUE),
    consumo_total_t = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_t = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_t = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_t = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    n_edades = n_distinct(edad_consumo),
    
    # PROPAGACIÓN DE INCERTIDUMBRE (suma en cuadratura)
    consumo_se_total_t = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    
    # Intervalos de confianza agregados (método conservador)
    consumo_lower_total_t = sum(consumo_lower_poblacional_t, na.rm = TRUE),
    consumo_upper_total_t = sum(consumo_upper_poblacional_t, na.rm = TRUE),
    
    # Coeficiente de variación del total
    consumo_cv_total = ifelse(consumo_total_t > 0, consumo_se_total_t / consumo_total_t, NA),
    
    # Intervalos de confianza basados en SE propagado (95% CI asumiendo normalidad)
    consumo_ci95_lower_t = consumo_total_t - (1.96 * consumo_se_total_t),
    consumo_ci95_upper_t = consumo_total_t + (1.96 * consumo_se_total_t),
    
    # Número de casos con incertidumbre disponible
    n_casos_con_se = sum(!is.na(consumo_se_poblacional_t)),
    .groups = 'drop'
  ) %>%
  arrange(year_consumo)

# Resumen anual por río (suma de todas las edades por río) CON PROPAGACIÓN
resumen_anual_rio <- consumo_poblacional_rio %>%
  group_by(year_consumo, Rio) %>%
  summarise(
    poblacion_total = sum(N_mar, na.rm = TRUE),
    consumo_total_t = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_t = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_t = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_t = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    n_edades = n_distinct(edad_consumo),
    
    # PROPAGACIÓN DE INCERTIDUMBRE POR RÍO
    consumo_se_total_t = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    
    # Intervalos de confianza por río
    consumo_lower_total_t = sum(consumo_lower_poblacional_t, na.rm = TRUE),
    consumo_upper_total_t = sum(consumo_upper_poblacional_t, na.rm = TRUE),
    
    # Coeficiente de variación por río
    consumo_cv_total = ifelse(consumo_total_t > 0, consumo_se_total_t / consumo_total_t, NA),
    
    # IC 95% basado en SE propagado
    consumo_ci95_lower_t = consumo_total_t - (1.96 * consumo_se_total_t),
    consumo_ci95_upper_t = consumo_total_t + (1.96 * consumo_se_total_t),
    
    n_casos_con_se = sum(!is.na(consumo_se_poblacional_t)),
    .groups = 'drop'
  ) %>%
  arrange(year_consumo, Rio)

# Resumen por edad (global) CON INCERTIDUMBRE
resumen_por_edad_global <- consumo_poblacional_global %>%
  group_by(edad_consumo, transicion_label) %>%
  summarise(
    años_disponibles = n_distinct(year_consumo),
    poblacion_promedio = mean(N_mar_total, na.rm = TRUE),
    consumo_individual_promedio_kg = mean(consumption_kg, na.rm = TRUE),
    consumo_individual_se_promedio_kg = mean(consumption_se_kg, na.rm = TRUE),
    consumo_total_promedio_t = mean(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_total_se_promedio_t = mean(consumo_se_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_promedio_t = mean(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_promedio_t = mean(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_promedio_t = mean(consumo_otros_poblacional_t, na.rm = TRUE),
    
    # Coeficiente de variación promedio por edad
    consumo_cv_promedio = mean(consumo_cv_poblacional, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(edad_consumo)

print("Resúmenes creados:")
print(paste("- Anual global:", nrow(resumen_anual_global), "años"))
print(paste("- Anual por río:", nrow(resumen_anual_rio), "registros"))
print(paste("- Por edad global:", nrow(resumen_por_edad_global), "edades"))

# ---- Mostrar resultados clave ----

print("\n=== RESUMEN ANUAL GLOBAL ===")
print(resumen_anual_global)

print("\n=== RESUMEN POR EDAD (GLOBAL) ===")
print(resumen_por_edad_global)

print("\n=== RESUMEN ANUAL POR RÍO (primeros registros) ===")
print(head(resumen_anual_rio, 10))

# ---- Crear gráficos ----

print("\n=== GENERANDO GRÁFICOS ===")

# Colores
colores_edad <- c("1" = "grey85", "2" = "grey65", "3" = "grey45", "4" = "grey25")
colores_edad_transicion <- c("1 → 2" = "grey85", "2 → 3" = "grey65", "3 → 4" = "grey45", "4 → 5" = "grey25")

colores_presa <- c("Anchoveta" = "#E76F51", "Sardina" = "#4472C4", "Otros" = "#8D6E63")
colores_rio <- c("Bueno" = "#1f77b4", "Tolten" = "#ff7f0e", "Imperial" = "#2ca02c", "Valdivia" = "#d62728")

# Gráfico 1: Evolución temporal consumo total global CON INCERTIDUMBRE
p1 <- ggplot(resumen_anual_global, aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "steelblue", alpha = 0.8, color = "white") +
  geom_errorbar(aes(ymin = pmax(0, consumo_ci95_lower_t), ymax = consumo_ci95_upper_t),
                width = 0.3, color = "darkred", size = 1, alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Consumo Total Poblacional por Año - ISOTOPES",
    subtitle = "Suma de todas las edades con intervalos de confianza 95% - Nivel global",
    x = "",
    y = "Consumo total (toneladas/año)\n",
    caption = "Barras de error: IC 95% basado en propagación de incertidumbre | Dieta edades 1-2: isótopos músculo 7F"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 20)
  ) +
  geom_text(aes(label = paste0(scales::comma(round(consumo_total_t, 0)), 
                              "\n(CV:", round(consumo_cv_total*100, 1), "%)")),
            vjust = -2, size = 3.5)

print(p1)
ggsave("output/consumo_poblacional_total_ISOTOPES.png", plot = p1, width = 14, height = 10, dpi = 300)

# Gráfico 2: Consumo por edad y año (global) CON INCERTIDUMBRE
p2 <- consumo_poblacional_global %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_poblacional_t, fill = factor(transicion_label))) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_edad_transicion, name = "Edad") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Consumo Poblacional por Edad y Año - ISOTOPES",
    subtitle = "Desglose etario del consumo total - Nivel global",
    x = "",
    y = "Consumo (toneladas/año)\n",
    caption = "Dieta edades 1-2: isótopos músculo 7F"
  ) +
  theme_minimal()  +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p2)
ggsave("output/consumo_poblacional_por_edad_ISOTOPES.png", plot = p2, width = 14, height = 10, dpi = 300)

# Gráfico 3: Consumo por tipo de presa (global anual)
consumo_presa_anual <- resumen_anual_global %>%
  select(year_consumo, consumo_anchoveta_t, consumo_sardina_t, consumo_otros_t) %>%
  pivot_longer(cols = starts_with("consumo_"), names_to = "presa", values_to = "consumo_t") %>%
  mutate(presa = case_when(
    presa == "consumo_anchoveta_t" ~ "Anchoveta",
    presa == "consumo_sardina_t" ~ "Sardina",
    presa == "consumo_otros_t" ~ "Otros"
  ))

p3 <- ggplot(consumo_presa_anual, aes(x = factor(year_consumo), y = consumo_t, fill = presa)) +
  geom_col(alpha = 0.9, color = "white", size = 0.3) +
  scale_fill_manual(values = colores_presa, name = "Tipo de presa") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Consumo Poblacional por Tipo de Presa - ISOTOPES",
    subtitle = "Distribución del consumo según composición de dieta - Nivel global",
    x = "",
    y = "Consumo (toneladas/año)\n",
    caption = "Dieta edades 1-2: isótopos músculo 7F"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom",
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p3)
ggsave("output/consumo_poblacional_por_presa_ISOTOPES.png", plot = p3, width = 14, height = 10, dpi = 300)

# Gráfico 4: Consumo por río CON INCERTIDUMBRE
p4 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t, fill = Rio)) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_rio, name = "Río") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Consumo Poblacional por Río y Año - ISOTOPES",
    subtitle = "Contribución de cada río al consumo total",
    x = "",
    y = "Consumo (toneladas/año)\n",
    caption = "Dieta edades 1-2: isótopos músculo 7F"
  ) +
  theme_minimal()  +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(p4)
ggsave("output/consumo_poblacional_por_rio_ISOTOPES.png", plot = p4, width = 14, height = 10, dpi = 300)

# Gráfico 5: Facetado por río - evolución temporal
p5 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "darkgreen", alpha = 0.8, color = "white") +
  facet_wrap(~ Rio, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Evolución del Consumo por Río - ISOTOPES",
    subtitle = "Tendencias temporales independientes por río",
    x = "",
    y = "Consumo (toneladas/año)\n",
    caption = "Dieta edades 1-2: isótopos músculo 7F"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    plot.caption = element_text(size = 10, hjust = 1, color = "darkgray"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 18)
  )

print(p5)
ggsave("output/consumo_poblacional_facetado_ISOTOPES.png", plot = p5, width = 14, height = 10, dpi = 300)

# ---- Guardar resultados ----

print("\n=== GUARDANDO RESULTADOS ===")

# Crear directorio si no existe
dir.create("data/data_raw/population-consumption", showWarnings = FALSE, recursive = TRUE)

# Guardar datos detallados CON SUFIJO ISOTOPES
write.csv(consumo_poblacional_global, "data/data_raw/population-consumption/consumo_poblacional_global_detallado_ISOTOPES.csv", row.names = FALSE)
write.csv(consumo_poblacional_rio, "data/data_raw/population-consumption/consumo_poblacional_por_rio_detallado_ISOTOPES.csv", row.names = FALSE)

# Guardar resúmenes CON SUFIJO ISOTOPES
write.csv(resumen_anual_global, "data/data_raw/population-consumption/resumen_consumo_anual_global_ISOTOPES.csv", row.names = FALSE)
write.csv(resumen_anual_rio, "data/data_raw/population-consumption/resumen_consumo_anual_por_rio_ISOTOPES.csv", row.names = FALSE)
write.csv(resumen_por_edad_global, "data/data_raw/population-consumption/resumen_consumo_por_edad_global_ISOTOPES.csv", row.names = FALSE)

# Crear tabla final consolidada CON INCERTIDUMBRE
tabla_consolidada <- resumen_anual_global %>%
  mutate(
    consumo_total_miles_t = round(consumo_total_t / 1000, 2),
    consumo_se_miles_t = round(consumo_se_total_t / 1000, 2),
    consumo_cv_pct = round(consumo_cv_total * 100, 1),
    consumo_ci95_lower_miles_t = round(pmax(0, consumo_ci95_lower_t / 1000), 2),
    consumo_ci95_upper_miles_t = round(consumo_ci95_upper_t / 1000, 2),
    consumo_anchoveta_miles_t = round(consumo_anchoveta_t / 1000, 2),
    consumo_sardina_miles_t = round(consumo_sardina_t / 1000, 2),
    consumo_otros_miles_t = round(consumo_otros_t / 1000, 2)
  ) %>%
  select(year_consumo, poblacion_total, 
         consumo_total_miles_t, consumo_se_miles_t, consumo_cv_pct,
         consumo_ci95_lower_miles_t, consumo_ci95_upper_miles_t,
         consumo_anchoveta_miles_t, consumo_sardina_miles_t, consumo_otros_miles_t,
         n_casos_con_se)

write.csv(tabla_consolidada, "data/data_raw/population-consumption/tabla_consolidada_con_incertidumbre_ISOTOPES.csv", row.names = FALSE)

print("Archivos guardados (ISOTOPES):")
print("- consumo_poblacional_global_detallado_ISOTOPES.csv")
print("- consumo_poblacional_por_rio_detallado_ISOTOPES.csv")
print("- resumen_consumo_anual_global_ISOTOPES.csv")
print("- resumen_consumo_anual_por_rio_ISOTOPES.csv")
print("- resumen_consumo_por_edad_global_ISOTOPES.csv")
print("- tabla_consolidada_con_incertidumbre_ISOTOPES.csv")

print("\nGráficos guardados (ISOTOPES):")
print("- consumo_poblacional_total_ISOTOPES.png")
print("- consumo_poblacional_por_edad_ISOTOPES.png")
print("- consumo_poblacional_por_presa_ISOTOPES.png")
print("- consumo_poblacional_por_rio_ISOTOPES.png")
print("- consumo_poblacional_facetado_ISOTOPES.png")

# ---- Análisis adicional de incertidumbre ----

print("\n=== ANÁLISIS DE INCERTIDUMBRE ===")

# Estadísticas de coeficientes de variación
cv_stats <- consumo_poblacional_global %>%
  filter(!is.na(consumo_cv_poblacional)) %>%
  group_by(edad_consumo) %>%
  summarise(
    cv_promedio_pct = round(mean(consumo_cv_poblacional) * 100, 1),
    cv_mediano_pct = round(median(consumo_cv_poblacional) * 100, 1),
    cv_min_pct = round(min(consumo_cv_poblacional) * 100, 1),
    cv_max_pct = round(max(consumo_cv_poblacional) * 100, 1),
    n_observaciones = n(),
    .groups = 'drop'
  ) %>%
  arrange(edad_consumo)

print("Estadísticas de Coeficiente de Variación por Edad:")
print(cv_stats)

# Estadísticas globales de incertidumbre
print("\nEstadísticas Globales de Incertidumbre:")
incert_global <- resumen_anual_global %>%
  summarise(
    cv_promedio_anual_pct = round(mean(consumo_cv_total, na.rm = TRUE) * 100, 1),
    cv_mediano_anual_pct = round(median(consumo_cv_total, na.rm = TRUE) * 100, 1),
    prop_años_con_cv_bajo_20pct = round(sum(consumo_cv_total < 0.20, na.rm = TRUE) / 
                                        sum(!is.na(consumo_cv_total)) * 100, 1),
    años_con_incertidumbre = sum(n_casos_con_se > 0)
  )

print(paste("CV promedio anual:", incert_global$cv_promedio_anual_pct, "%"))
print(paste("CV mediano anual:", incert_global$cv_mediano_anual_pct, "%"))
print(paste("Proporción años con CV < 20%:", incert_global$prop_años_con_cv_bajo_20pct, "%"))
print(paste("Años con datos de incertidumbre:", incert_global$años_con_incertidumbre, "de", nrow(resumen_anual_global)))

# ---- Resumen ejecutivo CON INCERTIDUMBRE ----

print("\n=== RESUMEN EJECUTIVO CON INCERTIDUMBRE (ISOTOPES) ===")
print(paste("Período analizado:", paste(range(años_overlap), collapse = "-")))
print(paste("Años de datos:", length(años_overlap)))
print("Fuente dieta edades 1-2: Isótopos músculo 7 fuentes")
print("Fuente dieta edades 3-5: Contenido estomacal")

consumo_total_periodo <- sum(resumen_anual_global$consumo_total_t)
consumo_promedio_anual <- mean(resumen_anual_global$consumo_total_t)

# Incertidumbre agregada
se_total_periodo <- sqrt(sum(resumen_anual_global$consumo_se_total_t^2, na.rm = TRUE))
cv_promedio_periodo <- se_total_periodo / consumo_total_periodo * 100

print(paste("Consumo total período:", scales::comma(round(consumo_total_periodo)), "toneladas"))
print(paste("Consumo promedio anual:", scales::comma(round(consumo_promedio_anual)), "±", 
            scales::comma(round(mean(resumen_anual_global$consumo_se_total_t, na.rm = TRUE))), 
            "toneladas/año"))
print(paste("Coeficiente de variación del período:", round(cv_promedio_periodo, 1), "%"))

# Por presa (promedio anual) con incertidumbre
anchoveta_promedio <- mean(resumen_anual_global$consumo_anchoveta_t)
sardina_promedio <- mean(resumen_anual_global$consumo_sardina_t)
otros_promedio <- mean(resumen_anual_global$consumo_otros_t)

print(paste("Consumo promedio anual por presa:"))
print(paste("- Anchoveta:", scales::comma(round(anchoveta_promedio)), "t/año"))
print(paste("- Sardina:", scales::comma(round(sardina_promedio)), "t/año"))
print(paste("- Otros:", scales::comma(round(otros_promedio)), "t/año"))

# Rango de incertidumbre
cv_min <- min(resumen_anual_global$consumo_cv_total, na.rm = TRUE) * 100
cv_max <- max(resumen_anual_global$consumo_cv_total, na.rm = TRUE) * 100

print(paste("Rango de coeficientes de variación anual:", round(cv_min, 1), "% -", round(cv_max, 1), "%"))

# Intervalos de confianza del consumo total
ic_lower_promedio <- mean(resumen_anual_global$consumo_ci95_lower_t, na.rm = TRUE)
ic_upper_promedio <- mean(resumen_anual_global$consumo_ci95_upper_t, na.rm = TRUE)

print(paste("Intervalo de confianza 95% promedio anual:", 
            scales::comma(round(ic_lower_promedio)), "-", 
            scales::comma(round(ic_upper_promedio)), "toneladas/año"))

print("\n¡Extrapolación poblacional CON INCERTIDUMBRE completada (ISOTOPES)!")
print("NOTA: La incertidumbre se propagó asumiendo independencia entre edades")
print("IMPORTANTE: Edades 1-2 usan dieta isotópica (músculo 7F), edades 3-5 usan estómagos")




# ==============================================================
# ESTADÍSTICAS PARA ACTUALIZAR TEXTO DE CONSUMO POBLACIONAL
# ==============================================================

# 1. Rango de consumo total anual (con años)
resumen_anual_global %>%
  select(year_consumo, consumo_total_t) %>%
  arrange(consumo_total_t) %>%
  slice(1, n()) %>%
  mutate(consumo_total_t = round(consumo_total_t, 0))

# Fold-difference
max_consumo <- max(resumen_anual_global$consumo_total_t)
min_consumo <- min(resumen_anual_global$consumo_total_t)
cat("\nFold-difference:", round(max_consumo / min_consumo, 1), "x\n")

# 2. Top 3 años de peak consumption
resumen_anual_global %>%
  select(year_consumo, consumo_total_t) %>%
  arrange(desc(consumo_total_t)) %>%
  head(3) %>%
  mutate(consumo_total_t = round(consumo_total_t, 0))

# 3. Consumo acumulado total y por presa
cat("\n=== CONSUMO ACUMULADO PERÍODO ===\n")
resumen_anual_global %>%
  summarise(
    period = paste(min(year_consumo), "-", max(year_consumo)),
    total_t = round(sum(consumo_total_t), 0),
    anchovy_t = round(sum(consumo_anchoveta_t), 0),
    sardine_t = round(sum(consumo_sardina_t), 0),
    others_t = round(sum(consumo_otros_t), 0)
  ) %>%
  mutate(
    anchovy_pct = round(anchovy_t / total_t * 100, 1),
    sardine_pct = round(sardine_t / total_t * 100, 1),
    others_pct = round(others_t / total_t * 100, 1)
  ) %>%
  print()

# 4. Rango de CV anual
cat("\n=== RANGO DE CV ANUAL ===\n")
resumen_anual_global %>%
  summarise(
    cv_min = round(min(consumo_cv_total, na.rm = TRUE) * 100, 1),
    cv_max = round(max(consumo_cv_total, na.rm = TRUE) * 100, 1)
  ) %>%
  print()

# 5. Años con mayor CV
resumen_anual_global %>%
  select(year_consumo, consumo_total_t, consumo_cv_total) %>%
  mutate(
    consumo_total_t = round(consumo_total_t, 0),
    cv_pct = round(consumo_cv_total * 100, 1)
  ) %>%
  arrange(desc(cv_pct)) %>%
  head(3)

# 6. Últimos 2 años (2021-2022 o equivalente)
resumen_anual_global %>%
  arrange(desc(year_consumo)) %>%
  head(2) %>%
  select(year_consumo, consumo_total_t) %>%
  mutate(consumo_total_t = round(consumo_total_t, 0))
