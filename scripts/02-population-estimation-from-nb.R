# Population estimation for Chinook salmon from Nb data

library(tidyverse)
library(ggplot2)

# ---- Lectura y procesamiento de datos ----

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

# Factor promedio
factor_promedio <- mean(ratios_tolten$ratio_Nb_N, na.rm = TRUE)

# ---- Extraer valores base automáticamente ----

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

# ---- Estimar poblaciones con factor promedio ----

estimaciones_promedio <- list()

for(i in 1:nrow(valores_base)) {
  rio_info <- valores_base[i, ]
  rio_name <- rio_info$Rio
  
  if(rio_name == "Tolten") {
    est <- retornantes_tolten %>%
      mutate(
        Rio = "Tolten",
        N_estimado = retornantes,
        metodo = "Datos_Reales"
      )
    
  } else if(!is.na(rio_info$N_base)) {
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

# ---- Gráfico ----

años_mostrar <- sort(unique(resultado_promedio$year))
datos_barras_promedio <- resultado_promedio %>%
  filter(year %in% años_mostrar) %>%
  mutate(year = factor(year))

p_barras_promedio <- ggplot(datos_barras_promedio, aes(x = year, y = N_estimado, fill = Rio)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::comma) +
  labs(
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

# ---- Exportación ----

write.csv(resultado_promedio, "data/data_raw/river-population/total_population_from_nb.csv", row.names = FALSE)

