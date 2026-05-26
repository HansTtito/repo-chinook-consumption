
# ================================================================
# ESTIMACIÓN POBLACIONAL DE SALMONES CHINOOK
# Comparación de métodos: Ratios Anuales vs Factor Promedio
# ================================================================

library(tidyverse)
library(ggplot2)

# ---- Funciones auxiliares ----

interpolar_ratios <- function(años_objetivo, ratios_conocidos) {
  ratios_completos <- data.frame(year = años_objetivo)
  
  ratios_completos <- ratios_completos %>%
    left_join(ratios_conocidos %>% dplyr::select(year, ratio_Nb_N), by = "year")
  
  # Para años faltantes, usar el promedio (método original)
  ratio_promedio <- mean(ratios_conocidos$ratio_Nb_N, na.rm = TRUE)
  ratios_completos$ratio_Nb_N[is.na(ratios_completos$ratio_Nb_N)] <- ratio_promedio
  
  return(ratios_completos)
}

# ---- Lectura y procesamiento de datos ----

print("=== CARGANDO DATOS ===")

# Cargar archivos
nb_rivers <- read.csv("data/data_raw/river-population/nb_rivers.csv", fileEncoding = "latin1")
poblacion_espinoza <- read.csv("data/data_raw/river-population/poblacion_espinoza.csv", fileEncoding = "latin1")

# Procesar retornantes Toltén (de Espinoza)
retornantes_tolten <- poblacion_espinoza %>%
  slice(5) %>%  # Fila 5 = número de retornantes
  pivot_longer(cols = -Temporadas, names_to = "temporada", values_to = "retornantes") %>%
  mutate(
    temporada = str_replace_all(temporada, "X", ""),
    temporada = str_replace(temporada, "\\.", "/")
  ) %>%
  separate(temporada, into = c("y_ini", "y_fin"), sep = "/", convert = TRUE) %>%
  mutate(
    y_fin = ifelse(y_fin < 50, y_fin + 2000, y_fin + 1900),
    year = y_fin
  ) %>%
  dplyr::select(year, retornantes) %>%
  filter(!is.na(retornantes), retornantes > 0) %>%
  arrange(year)

print(paste("Retornantes Toltén:", nrow(retornantes_tolten), "años"))

# Procesar datos Nb
nb_procesado <- nb_rivers %>%
  mutate(
    Rio = case_when(
      Site == "La Barra" ~ "Tolten",
      Site == "Puerto Saavedra" ~ "Imperial", 
      Site == "Valdivia" ~ "Valdivia",
      Site == "Bueno" ~ "Bueno",
      TRUE ~ Site
    )
  ) %>%
  dplyr::select(Rio, year, Nb, N_real)

# Calcular ratios Nb/N para Toltén
ratios_tolten <- nb_procesado %>%
  filter(Rio == "Tolten") %>%
  left_join(retornantes_tolten, by = "year") %>%
  filter(!is.na(Nb), !is.na(retornantes)) %>%
  mutate(ratio_Nb_N = Nb / retornantes) %>%
  arrange(year)

print("Ratios Nb/N calculados para Toltén:")
print(ratios_tolten)


# Cargar librerías necesarias
library(ggplot2)
library(ggrepel)  # Para etiquetas que no se superpongan

# Gráfico 1: Relación entre retornantes y ratio_Nb_N
plot1 <- ratios_tolten |>
  ggplot(aes(x = retornantes, y = ratio_Nb_N)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.7) +
  geom_text_repel(aes(label = year), 
                  size = 3.5,
                  box.padding = 0.5,
                  point.padding = 0.3) +
  labs(title = "Relación entre Retornantes y Ratio Nb/N",
       subtitle = "Río Toltén (2015-2021)",
       x = "Número de Retornantes",
       y = "Ratio Nb/N") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, face = "bold"))

# Gráfico 2: Relación entre retornantes y Nb
plot2 <- ratios_tolten |>
  ggplot(aes(x = retornantes, y = Nb)) +
  geom_point(size = 3, color = "coral", alpha = 0.7) +
  geom_text_repel(aes(label = year), 
                  size = 3.5,
                  box.padding = 0.5,
                  point.padding = 0.3) +
  labs(title = "Relación entre Retornantes y Nb",
       subtitle = "Río Toltén (2015-2021)",
       x = "Número de Retornantes",
       y = "Nb") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, face = "bold"))

# Mostrar los gráficos
print(plot1)
print(plot2)

# Opcional: Combinar ambos gráficos en una sola figura
library(patchwork)
combined_plot <- plot1 / plot2
print(combined_plot)

# Factor promedio
factor_promedio <- mean(ratios_tolten$ratio_Nb_N, na.rm = TRUE)
print(paste("Factor promedio Nb/N:", round(factor_promedio, 6)))

# ---- Extraer valores base automáticamente ----

print("\n=== VALORES BASE EXTRAÍDOS ===")

valores_base <- nb_procesado %>%
  group_by(Rio) %>%
  summarise(
    year_base = ifelse(any(!is.na(N_real)), 
                      year[which(!is.na(N_real))[1]], 
                      year[which(!is.na(Nb))[1]]),
    N_base = ifelse(any(!is.na(N_real)), 
                   N_real[which(!is.na(N_real))[1]], 
                   NA),
    Nb_base = ifelse(any(!is.na(Nb)), 
                    Nb[which(!is.na(Nb))[1]], 
                    NA),
    .groups = 'drop'
  )

# Agregar Toltén con datos reales
valores_base <- valores_base %>%
  filter(Rio != "Tolten") %>%
  bind_rows(data.frame(
    Rio = "Tolten",
    year_base = NA,
    N_base = NA,
    Nb_base = NA
  ))

print(valores_base)

# ---- Estimar poblaciones con ratios anuales ----

print("\n=== ESTIMACIONES CON RATIOS ANUALES ===")

# Interpolar ratios para todos los años
años_objetivo <- sort(unique(retornantes_tolten$year))
ratios_completos <- interpolar_ratios(años_objetivo, ratios_tolten)

# Índice poblacional
indice_poblacional <- retornantes_tolten %>%
  left_join(ratios_completos, by = "year") %>%
  mutate(indice_normalizado = retornantes / max(retornantes, na.rm = TRUE))

# Estimaciones ratios anuales
estimaciones_ratios <- list()

for(i in 1:nrow(valores_base)) {
  rio_info <- valores_base[i, ]
  rio_name <- rio_info$Rio
  
  if(rio_name == "Tolten") {
    # Toltén usa datos reales
    est <- indice_poblacional %>%
      mutate(
        Rio = "Tolten",
        N_estimado = retornantes,
        metodo = "Datos_Reales"
      )
    
  } else if(!is.na(rio_info$N_base)) {
    # Si tiene N_base real (Imperial)
    year_base <- rio_info$year_base
    N_base <- rio_info$N_base
    
    indice_base <- indice_poblacional$indice_normalizado[indice_poblacional$year == year_base]
    if(length(indice_base) == 0) indice_base <- 1
    
    est <- indice_poblacional %>%
      mutate(
        Rio = rio_name,
        N_estimado = N_base * (indice_normalizado / indice_base),
        metodo = "N_base_tendencia"
      )
    
  } else if(!is.na(rio_info$Nb_base)) {
    # Si solo tiene Nb_base (Valdivia, Bueno)
    Nb_base <- rio_info$Nb_base
    
    est <- indice_poblacional %>%
      mutate(
        Rio = rio_name,
        N_estimado = Nb_base / ratio_Nb_N,  # Ratios variables
        metodo = "Nb_ratios_anuales"
      )
  }
  
  estimaciones_ratios[[rio_name]] <- est %>%
    dplyr::select(year, Rio, N_estimado, metodo)
}

resultado_ratios <- bind_rows(estimaciones_ratios) %>%
  arrange(Rio, year)

# ---- Estimar poblaciones con factor promedio ----

print("\n=== ESTIMACIONES CON FACTOR PROMEDIO ===")

estimaciones_promedio <- list()

for(i in 1:nrow(valores_base)) {
  rio_info <- valores_base[i, ]
  rio_name <- rio_info$Rio
  
  if(rio_name == "Tolten") {
    # Toltén usa datos reales
    est <- retornantes_tolten %>%
      mutate(
        Rio = "Tolten",
        N_estimado = retornantes,
        metodo = "Datos_Reales"
      )
    
  } else if(!is.na(rio_info$N_base)) {
    # Si tiene N_base real (Imperial)
    year_base <- rio_info$year_base
    N_base <- rio_info$N_base
    
    indice_base <- retornantes_tolten$retornantes[retornantes_tolten$year == year_base]
    if(length(indice_base) == 0) indice_base <- max(retornantes_tolten$retornantes)
    
    est <- retornantes_tolten %>%
      mutate(
        Rio = rio_name,
        N_estimado = N_base * (retornantes / indice_base),
        metodo = "N_base_tendencia"
      )
    
  } else if(!is.na(rio_info$Nb_base)) {
    # Si solo tiene Nb_base (Valdivia, Bueno) - usar factor promedio
    Nb_base <- rio_info$Nb_base
    year_base <- rio_info$year_base
    N_base_estimado <- Nb_base / factor_promedio
    
    indice_base <- retornantes_tolten$retornantes[retornantes_tolten$year == year_base]
    if(length(indice_base) == 0) indice_base <- max(retornantes_tolten$retornantes)
    
    est <- retornantes_tolten %>%
      mutate(
        Rio = rio_name,
        N_estimado = N_base_estimado * (retornantes / indice_base),
        metodo = "Nb_factor_promedio"
      )
  }
  
  estimaciones_promedio[[rio_name]] <- est %>%
    dplyr::select(year, Rio, N_estimado, metodo)
}

resultado_promedio <- bind_rows(estimaciones_promedio) %>%
  arrange(Rio, year)

# ---- Comparación de métodos ----

print("\n=== COMPARACIÓN DE MÉTODOS ===")

comparacion <- resultado_ratios %>%
  rename(N_ratios_anuales = N_estimado) %>%
  dplyr::select(year, Rio, N_ratios_anuales) %>%
  left_join(
    resultado_promedio %>% 
      rename(N_factor_promedio = N_estimado) %>%
      dplyr::select(year, Rio, N_factor_promedio),
    by = c("year", "Rio")
  ) %>%
  mutate(
    diferencia_absoluta = N_ratios_anuales - N_factor_promedio,
    diferencia_pct = round((N_ratios_anuales / N_factor_promedio - 1) * 100, 1)
  )

# Mostrar estadísticas de comparación
stats_comparacion <- comparacion %>%
  filter(Rio != "Tolten") %>%  # Excluir Toltén (datos reales)
  group_by(Rio) %>%
  summarise(
    diferencia_promedio_pct = round(mean(abs(diferencia_pct), na.rm = TRUE), 1),
    diferencia_max_pct = round(max(abs(diferencia_pct), na.rm = TRUE), 1),
    n_años = n(),
    .groups = 'drop'
  )

print("Diferencias porcentuales entre métodos:")
print(stats_comparacion)

# Mostrar muestra de la comparación
print("\nMuestra de comparación (primeros años):")
print(head(comparacion, 12))

# ---- Gráficos de barras ----

print("\n=== GENERANDO GRÁFICOS DE BARRAS ===")

# Seleccionar años para mostrar (cada 2 años para no saturar)
# años_mostrar <- seq(min(comparacion$year), max(comparacion$year), by = 2)
años_mostrar <- sort(unique(comparacion$year))
# GRÁFICO 1: Ratios Anuales
datos_barras_ratios <- resultado_ratios %>%
  filter(year %in% años_mostrar) %>%
  mutate(year = factor(year))

p_barras_ratios <- ggplot(datos_barras_ratios, aes(x = year, y = N_estimado, fill = Rio)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Estimaciones con Ratios Anuales",
    subtitle = "Retornantes estimados por río usando ratios Nb/N variables",
    x = "Año",
    y = "Número de retornantes estimados",
    fill = "Río"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_barras_ratios)

# GRÁFICO 2: Factor Promedio  
datos_barras_promedio <- resultado_promedio %>%
  filter(year %in% años_mostrar) %>%
  mutate(year = factor(year))

p_barras_promedio <- ggplot(datos_barras_promedio, aes(x = year, y = N_estimado, fill = Rio)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    # title = "Estimaciones con Factor Promedio",
    # subtitle = "Retornantes estimados por río usando factor Nb/N fijo",
    x = "",
    y = "Estimated number of returning adults\n",
    fill = "River"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

print(p_barras_promedio)

# GRÁFICO 3: Comparación lado a lado
datos_barras_combinados <- bind_rows(
  resultado_ratios %>%
    filter(year %in% años_mostrar) %>%
    mutate(Método = "Ratios Anuales"),
  
  resultado_promedio %>%
    filter(year %in% años_mostrar) %>%
    mutate(Método = "Factor Promedio")
) %>%
  mutate(
    year = factor(year),
    Método = factor(Método, levels = c("Ratios Anuales", "Factor Promedio"))
  )

p_barras_combinado <- ggplot(datos_barras_combinados, aes(x = year, y = N_estimado, fill = Rio)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  facet_wrap(~ Método, ncol = 2) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Comparación: Ratios Anuales vs Factor Promedio",
    subtitle = "Los 4 ríos lado a lado para cada método",
    x = "Año",
    y = "Número de retornantes estimados",
    fill = "Río"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgray")
  )

print(p_barras_combinado)

# GRÁFICO 4: Un año específico detallado (el más reciente)
año_detalle <- max(comparacion$year)

datos_año_detalle <- bind_rows(
  resultado_ratios %>%
    filter(year == año_detalle) %>%
    mutate(Método = "Ratios Anuales"),
  
  resultado_promedio %>%
    filter(year == año_detalle) %>%
    mutate(Método = "Factor Promedio")
) %>%
  mutate(Método = factor(Método, levels = c("Ratios Anuales", "Factor Promedio")))

p_detalle_año <- ggplot(datos_año_detalle, aes(x = Rio, y = N_estimado, fill = Método)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("Ratios Anuales" = "#2E86AB", "Factor Promedio" = "#A23B72")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = paste("Comparación Detallada -", año_detalle),
    subtitle = "Estimaciones por río usando ambos métodos",
    x = "Río",
    y = "Número de retornantes estimados",
    fill = "Método"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  geom_text(aes(label = scales::comma(round(N_estimado, 0))), 
            position = position_dodge(width = 0.6), 
            vjust = -0.3, size = 3.5, fontface = "bold")

print(p_detalle_año)

# ---- Guardar resultado final (solo factor promedio) ----

print("\n=== GUARDANDO RESULTADO FINAL ===")

# Guardar solo estimaciones con factor promedio
write.csv(resultado_promedio, "data/data_raw/river-population/total_population_from_nb.csv", row.names = FALSE)

print("Archivo guardado: total_population_from_nb.csv")

# Mostrar resumen final
print("\n=== RESUMEN FINAL ===")
print(paste("Factor promedio usado:", round(factor_promedio, 6)))
print("Estimaciones guardadas para 4 ríos:")

resumen_final <- resultado_promedio %>%
  group_by(Rio) %>%
  summarise(
    años = n(),
    min_N = round(min(N_estimado)),
    max_N = round(max(N_estimado)),
    promedio_N = round(mean(N_estimado)),
    .groups = 'drop'
  )

print(resumen_final)

print("\n¡Análisis completado!")
print("- Comparación mostrada en gráfico")
print("- Datos guardados en outputs/poblacion_factor_promedio.csv")


datos_usados <- datos_barras_promedio 






# Ver los datos de Nb por río
print(nb_procesado)

# Y los valores base extraídos
print(valores_base)
