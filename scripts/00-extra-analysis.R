# CÓDIGO LIMPIO - ANÁLISIS POR TEMPORADAS
library(readr)
library(janitor)
library(ggplot2)
library(tidyr)
library(dplyr)

# Cargar datos
biologico_salmones <- read.csv("data/data_raw/biological-data/data_chinook_limpia-Kevin.csv", fileEncoding = "latin1") |>
  clean_names() |> as.data.frame()

# CREAR TEMPORADAS (octubre-abril)
biologico_temporadas <- biologico_salmones |>
  dplyr::filter(species == "Oncorhynchus tshawytscha", 
                !is.na(year), !is.na(month)) |>
  mutate(
    temporada = case_when(
      month >= 10 ~ paste0(year, "-", year + 1),
      month <= 4 ~ paste0(year - 1, "-", year),
      TRUE ~ NA_character_
    ),
    temporada_inicio = case_when(
      month >= 10 ~ year,
      month <= 4 ~ year - 1,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::filter(!is.na(temporada), !is.na(month),
                location %in% c("La Barra", "Tolten","Puerto Saavedra", "Queule","ToltãN")) |> 
  distinct() |> 
  select(year, month, species, temporada, temporada_inicio, tw_g, 
         fl_mm, tl_mm, age, sex, region, location,
         day, hl_mm, stomach_weight_g, total_content_g, sardina_g, anchoveta_g,
         teleost_digested_g, bolus_g, p_sc, p_an, p_tni, p_rni, p_is,
         p_ca)

# PREPARAR DATOS DE DIETA
chinook_diet <- biologico_temporadas |>
  clean_names() |>
  dplyr::select(temporada, region, location, species, year, month,
                day, fl_mm, hl_mm, tl_mm, tw_g, sex, age,
                stomach_weight_g, total_content_g, sardina_g, anchoveta_g,
                teleost_digested_g, bolus_g, p_sc, p_an, p_tni, p_rni, p_is,
                p_ca) |>
  dplyr::filter(species == "Oncorhynchus tshawytscha" | is.na(species), 
                !is.na(age),
                age >= 2 & age <= 5,
                (!is.na(sardina_g) | !is.na(anchoveta_g) | 
                 !is.na(teleost_digested_g) | !is.na(bolus_g) | 
                 !is.na(p_sc) | !is.na(p_an) | !is.na(p_tni) | 
                 !is.na(p_is) | !is.na(p_ca))) |>
  mutate(
    sardina_total = coalesce(sardina_g, 0) + coalesce(p_sc, 0),
    anchoveta_total = coalesce(anchoveta_g, 0) + coalesce(p_an, 0),
    otros_total = coalesce(teleost_digested_g, 0) + coalesce(bolus_g, 0) +
                  coalesce(p_tni, 0) + coalesce(p_is, 0) + coalesce(p_ca, 0),
    peso_total_individuo = sardina_total + anchoveta_total + otros_total
  ) |>
  filter(peso_total_individuo > 0)

# CALCULAR TAMAÑOS DE MUESTRA POR TEMPORADA

# 1. Individuos con edad por temporada
edad_por_temporada <- biologico_temporadas |>
  filter(!is.na(age), age > 0) |>
  group_by(temporada) |>
  summarise(n_edad = n(), .groups = "drop") |>
  arrange(temporada)

# 2. Individuos con peso por temporada
peso_por_temporada <- biologico_temporadas |>
  filter(!is.na(tw_g)) |>
  group_by(temporada) |>
  summarise(n_peso = n(), .groups = "drop") |>
  arrange(temporada)

# 3. Individuos con edad y peso por temporada
peso_edad_por_temporada <- biologico_temporadas |>
  filter(!is.na(age), !is.na(tw_g), age > 0) |>
  group_by(temporada) |>
  summarise(n_peso_edad = n(), .groups = "drop") |>
  arrange(temporada)

# 4. Individuos con datos de dieta por temporada
dieta_por_temporada <- chinook_diet |>
  group_by(temporada) |>
  summarise(n_dieta = n(), .groups = "drop") |>
  arrange(temporada)

# COMBINAR TODOS LOS TAMAÑOS DE MUESTRA
tamaños_muestra <- edad_por_temporada |>
  full_join(peso_por_temporada, by = "temporada") |>
  full_join(peso_edad_por_temporada, by = "temporada") |>
  full_join(dieta_por_temporada, by = "temporada") |>
  replace_na(list(n_edad = 0, n_peso = 0, n_peso_edad = 0, n_dieta = 0)) |>
  arrange(temporada)

# ANÁLISIS POR EDAD Y TEMPORADA
temporada_edad <- biologico_temporadas |>
  filter(!is.na(age), age > 0) |>
  group_by(temporada, age) |>
  summarise(n_individuos = n(), .groups = "drop") |> 
  arrange(temporada)

# ESTADÍSTICAS POR TEMPORADA
estadisticas_temporada <- biologico_temporadas |>
  filter(!is.na(age), !is.na(tw_g), age > 0) |>
  group_by(temporada) |>
  summarise(
    n_total = n(),
    edad_promedio = mean(age),
    edad_mediana = median(age),
    peso_promedio = mean(tw_g),
    peso_mediana = median(tw_g),
    peso_sd = sd(tw_g),
    edad_min = min(age),
    edad_max = max(age),
    peso_min = min(tw_g),
    peso_max = max(tw_g),
    .groups = "drop"
  ) |>
  arrange(temporada)

# GUARDAR DATOS LIMPIOS
write.csv(biologico_temporadas, "data/data_raw/biological-data/data_chinook_cleaned.csv", row.names = FALSE)

# GRÁFICOS DE TAMAÑOS DE MUESTRA


biologico_temporadas |> 
  filter(!is.na(age), !is.na(tw_g), age > 0) |>
  ggplot(aes(x=factor(age), y = tw_g/1000, fill = factor(age))) +
  geom_boxplot() +
  facet_wrap(~temporada) +
  scale_fill_viridis_d(option = "plasma") +
  theme_classic() +
  labs(x = "\nAge",
       y = "Total Weight (kg)\n",
       fill = "Age") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))




# 1. Gráfico de barras - Tamaños de muestra por temporada
datos_temporada_largo <- tamaños_muestra |>
  pivot_longer(cols = c(n_edad, n_peso, n_peso_edad, n_dieta),
               names_to = "tipo_dato",
               values_to = "n_individuos") |>
  mutate(tipo_dato = case_when(
    tipo_dato == "n_edad" ~ "Con edad",
    tipo_dato == "n_peso" ~ "Con peso",
    tipo_dato == "n_peso_edad" ~ "Edad y peso",
    tipo_dato == "n_dieta" ~ "Con dieta"
  ))

ggplot(datos_temporada_largo, aes(x = temporada, y = n_individuos, fill = tipo_dato)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = n_individuos), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  labs(title = "Tamaño de muestra por temporada - Oncorhynchus tshawytscha",
       subtitle = "Temporadas: Octubre-Abril del año siguiente",
       x = "Temporada", y = "Número de individuos", fill = "Datos disponibles") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("Con edad" = "#2E86AB",
                               "Con peso" = "#A23B72",
                               "Edad y peso" = "#F18F01",
                               "Con dieta" = "#3BA55C"))

# 2. Gráfico de líneas - Evolución temporal de tamaños de muestra
ggplot(datos_temporada_largo, aes(x = temporada, y = n_individuos, color = tipo_dato, group = tipo_dato)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Evolución temporal de tamaños de muestra",
       subtitle = "Número de individuos con diferentes tipos de datos por temporada",
       x = "Temporada", y = "Número de individuos", color = "Tipo de datos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(values = c("Con edad" = "#2E86AB",
                                "Con peso" = "#A23B72",
                                "Edad y peso" = "#F18F01",
                                "Con dieta" = "#3BA55C"))

# 3. Heatmap - Individuos por temporada y edad
ggplot(temporada_edad, aes(x = temporada, y = factor(age), fill = n_individuos)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_individuos), color = "black", size = 4.5) +
  scale_fill_gradient(low = "white", high = "#2E86AB", na.value = "grey90") +
  labs(title = "Distribución de individuos por temporada y edad",
       x = "Temporada", y = "Edad (años)", fill = "N° Individuos") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))


# 3. Heatmap - Individuos por temporada y edad (todos con edad)
ggplot(temporada_edad, aes(x = temporada, y = factor(age), fill = n_individuos)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_individuos), color = "black", size = 4.5) +
  scale_fill_gradient(low = "white", high = "#2E86AB", na.value = "grey90") +
  labs(title = "Distribución de individuos por temporada y edad",
       subtitle = "Todos los individuos con edad registrada",
       x = "Temporada", y = "Edad (años)", fill = "N° Individuos") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))

# 3b. Heatmap - Individuos CON PESO por temporada y edad
temporada_edad_peso <- biologico_temporadas |>
  filter(!is.na(age), !is.na(tw_g), age > 0) |>
  group_by(temporada, age) |>
  summarise(n_individuos = n(), .groups = "drop") |> 
  arrange(temporada)

ggplot(temporada_edad_peso, aes(x = temporada, y = factor(age), fill = n_individuos)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_individuos), color = "black", size = 4.5) +
  scale_fill_gradient(low = "white", high = "#A23B72", na.value = "grey90") +
  labs(title = "Individuos con peso por temporada y edad",
       subtitle = "Solo individuos con datos de peso corporal registrado",
       x = "Temporada", y = "Edad (años)", fill = "N° Individuos") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))


# 4. Gráfico de barras apiladas - Composición por edad en cada temporada
ggplot(temporada_edad |> filter(!(temporada == "2023-2024")), aes(x = temporada, y = n_individuos, fill = factor(age))) +
  geom_col(position = "stack", alpha = 0.8) +
  labs(title = "Composición etaria por temporada",
       subtitle = "Distribución de individuos según edad en cada temporada",
       x = "Temporada", y = "Número de individuos", fill = "Edad (años)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_fill_viridis_d(option = "plasma")


# 4. Gráfico de barras apiladas - Composición por edad en cada temporada
ggplot(temporada_edad |> filter(!(temporada == "2023-2024")), aes(x = temporada, y = n_individuos, fill = factor(age))) +
  geom_col(position = "fill", alpha = 0.8) +
  labs(title = "Composición etaria por temporada",
       subtitle = "Distribución de individuos según edad en cada temporada",
       x = "", y = "Número de individuos\n", fill = "Edad (años)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold")) +
  scale_fill_viridis_d(option = "plasma")

# 5. Resumen de tamaños de muestra
print("RESUMEN DE TAMAÑOS DE MUESTRA POR TEMPORADA")
print("==========================================")
print(tamaños_muestra)

print("\nESTADÍSTICAS GENERALES:")
cat("Total temporadas:", nrow(tamaños_muestra), "\n")
cat("Total individuos con edad:", sum(tamaños_muestra$n_edad), "\n")
cat("Total individuos con peso:", sum(tamaños_muestra$n_peso), "\n")
cat("Total individuos con edad y peso:", sum(tamaños_muestra$n_peso_edad), "\n")
cat("Total individuos con datos de dieta:", sum(tamaños_muestra$n_dieta), "\n")
cat("Temporada con más individuos (edad):", tamaños_muestra$temporada[which.max(tamaños_muestra$n_edad)], "\n")
cat("Temporada con más individuos (dieta):", tamaños_muestra$temporada[which.max(tamaños_muestra$n_dieta)], "\n")