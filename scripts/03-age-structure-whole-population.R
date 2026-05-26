# # ================================================================
# # DISTRIBUCIÓN DE POBLACIONES POR ESTRUCTURA DE EDADES
# # Toma estimaciones poblacionales y las distribuye por edades
# # ================================================================

# library(tidyverse)
# library(ggplot2)
# library(janitor)
# box::use(
#   dplyr = dplyr[filter, select, case_when, mutate, reframe, n, arrange],
#   janitor = janitor[clean_names],
#   readr = readr[read_csv],
#   ggplot2 = ggplot2[ggplot, aes, geom_point, geom_line, labs, theme_minimal, 
#                      scale_x_continuous, theme, element_text, ggtitle]
# )
# # ---- Cargar datos ----

# print("=== CARGANDO DATOS ===")

# # 1. Poblaciones estimadas (del script anterior)
# poblacion_estimada <- read.csv("data/data_raw/river-population/total_population_from_nb.csv")

# # 2. Datos biológicos chinook
# biologico_salmones <- read.csv("data/data_raw/biological-data/data_chinook_cleaned.csv", fileEncoding = "latin1")

# print(paste("Poblaciones estimadas:", nrow(poblacion_estimada), "registros"))
# print(paste("Años disponibles:", paste(sort(unique(poblacion_estimada$year)), collapse = ", ")))
# print(paste("Ríos:", paste(unique(poblacion_estimada$Rio), collapse = ", ")))

# # ---- Procesar datos biológicos por temporada ----

# print("\n=== PROCESANDO DATOS BIOLÓGICOS POR TEMPORADA ===")

# chinook_cleaned <- biologico_salmones |> 
#   janitor$clean_names() |> 
#   dplyr$select(age, tw_g, location, species, year, month, fl_mm) |> 
#   dplyr$filter(
#     species == "Oncorhynchus tshawytscha",
#     !is.na(age), 
#     !is.na(year),
#     !is.na(month)
#   )

# print(paste("Datos biológicos con mes:", nrow(chinook_cleaned), "individuos"))

# # Asignar temporada de retorno según mes de captura
# chinook_temporada <- chinook_cleaned |>
#   dplyr$mutate(
#     temporada = case_when(
#       # Nov-Dic del año X corresponde a temporada X/(X+1)
#       month %in% c(10, 11, 12) ~ paste0(year, "/", substr(year + 1, 3, 4)),
#       # Ene-Mar del año Y corresponde a temporada (Y-1)/Y  
#       month %in% c(1, 2, 3, 4) ~ paste0(year - 1, "/", substr(year, 3, 4)),
#       # Otros meses no corresponden a temporada de retorno
#       TRUE ~ NA_character_
#     ),
#     # Año final de la temporada (para hacer join)
#     year_temporada = case_when(
#       month %in% c(10, 11, 12) ~ year + 1,
#       month %in% c(1, 2, 3, 4) ~ year,
#       TRUE ~ NA_real_
#     )
#   ) |>
#   dplyr$filter(!is.na(temporada)) |>  # Solo muestras de temporada de retorno
#   dplyr$select(-year, -month) |>  # Remover year original
#   dplyr$rename(year = year_temporada)  # Usar año de temporada

# print(paste("Muestras en temporada de retorno:", nrow(chinook_temporada), "individuos"))

# # Resumen de muestras por temporada
# muestras_por_temporada <- chinook_temporada |>
#   dplyr$group_by(temporada, year) |>
#   dplyr$summarise(
#     n_individuos = n(),
#     edades_min_max = paste(min(age, na.rm = TRUE), max(age, na.rm = TRUE), sep = "-"),
#     meses_muestreo = "Nov-Mar",
#     .groups = 'drop'
#   ) |>
#   dplyr$arrange(year)

# print("Muestras biológicas por temporada de retorno:")
# print(muestras_por_temporada)

# # ---- Calcular proporciones de edad por temporada ----

# print("\n=== CALCULANDO PROPORCIONES DE EDAD POR TEMPORADA ===")

# estructura_edades <- chinook_temporada |>
#   dplyr$group_by(year, age) |>
#   dplyr$summarise(n_individuos = n(), .groups = "drop") |>
#   dplyr$group_by(year) |>
#   dplyr$mutate(
#     total_year = sum(n_individuos),
#     proporcion = n_individuos / total_year
#   ) |>
#   dplyr$select(year, age, proporcion, n_individuos) |>
#   dplyr$arrange(year, age) |> 
#   as.data.frame()

# print("Estructura de edades calculada por temporada:")
# print(estructura_edades)

# tabla_edades_export <- estructura_edades %>%
#  # Convertir proporciones a porcentajes
#  mutate(proporcion_pct = round(proporcion * 100, 1)) %>%
#  # Seleccionar y renombrar columnas
#  select(year, age, proporcion_pct, n_individuos) %>%
#  # Pivotar para tener edades como columnas
#  pivot_wider(
#    names_from = age, 
#    values_from = c(proporcion_pct, n_individuos),
#    names_sep = "_"
#  ) %>%
#  # Reordenar columnas
#  select(year, 
#         proporcion_pct_2, n_individuos_2,
#         proporcion_pct_3, n_individuos_3, 
#         proporcion_pct_4, n_individuos_4,
#         proporcion_pct_5, n_individuos_5) %>%
#  # Reemplazar NA con 0
#  replace(is.na(.), 0) %>%
#  # Calcular total por año
#  rowwise() %>%
#  mutate(
#    n_total = sum(n_individuos_2, n_individuos_3, n_individuos_4, n_individuos_5, na.rm = TRUE)
#  ) %>%
#  ungroup() %>%
#  # Reordenar para poner total al principio
#  select(year, n_total, everything())

# # Exportar a CSV
# write.csv(tabla_edades_export, "data/data_raw/biological-data/age_composition_by_year.csv", row.names = FALSE)


# # Mostrar temporadas disponibles vs años poblacionales
# temporadas_disponibles <- sort(unique(estructura_edades$year))
# años_poblacion <- sort(unique(poblacion_estimada$year))

# print(paste("\nTemporadas con estructura de edades:", paste(temporadas_disponibles, collapse = ", ")))
# print(paste("Años con poblaciones estimadas:", paste(años_poblacion, collapse = ", ")))

# # Verificar que proporciones suman 1 por temporada
# verificacion <- estructura_edades |>
#   dplyr$group_by(year) |>
#   dplyr$summarise(suma_proporciones = round(sum(proporcion), 3), .groups = 'drop')

# print("\nVerificación por temporada (deben sumar 1.000):")
# print(verificacion)

# # ---- Manejo de años faltantes ----

# años_poblacion <- sort(unique(poblacion_estimada$year))
# años_estructura <- sort(unique(estructura_edades$year))

# años_faltantes <- setdiff(años_poblacion, temporadas_disponibles)

# if(length(años_faltantes) > 0) {
#   print(paste("\nTemporadas sin estructura de edades:", paste(años_faltantes, collapse = ", ")))
#   print("Usando promedio histórico de temporadas disponibles")
  
#   # Calcular estructura promedio
#   estructura_promedio <- estructura_edades |>
#     dplyr$group_by(age) |>
#     dplyr$summarise(proporcion = mean(proporcion, na.rm = TRUE), .groups = 'drop')
  
#   # Agregar temporadas faltantes con estructura promedio
#   for(temporada_faltante in años_faltantes) {
#     estructura_temporada_faltante <- estructura_promedio |>
#       dplyr$mutate(
#         year = temporada_faltante,
#         n_individuos = NA
#       ) |>
#       dplyr$select(year, age, proporcion, n_individuos)
    
#     estructura_edades <- rbind(estructura_edades, estructura_temporada_faltante)
#   }
  
#   print("Estructura promedio aplicada a temporadas faltantes:")
#   print(estructura_promedio)
# } else {
#   print("\n✅ Todas las temporadas poblacionales tienen estructura de edades")
# }

# # ---- Aplicar estructura de edades a poblaciones estimadas ----

# print("\n=== DISTRIBUYENDO POBLACIONES POR EDAD (TEMPORADAS CORREGIDAS) ===")

# poblacion_por_edad <- poblacion_estimada |>
#   dplyr$left_join(estructura_edades, by = "year") |>
#   dplyr$filter(!is.na(proporcion)) |>  # Eliminar combinaciones sin estructura
#   dplyr$mutate(
#     N_estimado_por_edad = N_estimado * proporcion
#   ) |>
#   dplyr$select(year, Rio, age, N_estimado_total = N_estimado, 
#                proporcion, N_estimado_por_edad, metodo) |>
#   dplyr$arrange(Rio, year, age)

# print(paste("Registros generados:", nrow(poblacion_por_edad)))

# # Mostrar ejemplo
# print("\nEjemplo de resultados:")
# ejemplo <- poblacion_por_edad |>
#   dplyr$filter(year == max(year), Rio == "Imperial") |>
#   dplyr$mutate(
#     N_estimado_por_edad = round(N_estimado_por_edad, 0),
#     proporcion_pct = round(proporcion * 100, 1)
#   )
# print(ejemplo)

# # ---- Validación: verificar que suma por año/río da el total ----

# print("\n=== VALIDACIÓN ===")

# validacion <- poblacion_por_edad |>
#   dplyr$group_by(year, Rio) |>
#   dplyr$summarise(
#     N_total_original = first(N_estimado_total),
#     N_suma_edades = sum(N_estimado_por_edad),
#     diferencia = abs(N_total_original - N_suma_edades),
#     .groups = 'drop'
#   ) |>
#   dplyr$mutate(diferencia_relativa = diferencia / N_total_original)

# max_diferencia <- max(validacion$diferencia_relativa, na.rm = TRUE)
# print(paste("Máxima diferencia relativa:", round(max_diferencia * 100, 4), "%"))

# if(max_diferencia < 0.001) {
#   print("✅ Validación exitosa - Las sumas coinciden")
# } else {
#   print("⚠️  Hay diferencias en las sumas - revisar cálculos")
#   print(validacion[validacion$diferencia_relativa > 0.001, ])
# }

# # ---- Crear tabla resumen ----

# print("\n=== CREANDO TABLA RESUMEN ===")

# resumen_por_rio_edad <- poblacion_por_edad |>
#   dplyr$group_by(Rio, age) |>
#   dplyr$summarise(
#     años_disponibles = n(),
#     N_promedio = round(mean(N_estimado_por_edad, na.rm = TRUE), 0),
#     N_min = round(min(N_estimado_por_edad, na.rm = TRUE), 0),
#     N_max = round(max(N_estimado_por_edad, na.rm = TRUE), 0),
#     .groups = 'drop'
#   ) |>
#   dplyr$arrange(Rio, age)

# print("Resumen por río y edad:")
# print(resumen_por_rio_edad)

# # ---- Gráficos ----

# print("\n=== GENERANDO GRÁFICOS ===")

# # Colores azules oceánicos por edad
# colores_edad <- c(
#  "1" = "#E3F2FD",  # azul muy claro
#  "2" = "#90CAF9",  # azul claro
#  "3" = "#42A5F5",  # azul medio
#  "4" = "#1976D2",  # azul oscuro
#  "5" = "#0D47A1"   # azul muy oscuro
# )

# # Gráfico 1: Estructura de edades por año
# p1 <- ggplot(estructura_edades, aes(x = factor(year), y = proporcion, fill = factor(age))) +
#   geom_col(alpha = 0.8) +
#   scale_fill_manual(values = colores_edad, name = "Edad") +
#   labs(
#     title = "Estructura de Edades por Año",
#     subtitle = "Proporciones observadas en muestras biológicas",
#     x = "",
#     y = "Proporción (Individuos por Edad)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 24, face = "bold"),
#     plot.subtitle = element_text(size = 20, margin = margin(b = 20)),
#     axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 20),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(color = "black", size = 22),
#     legend.text = element_text(size = 20),
#     legend.title = element_text(size = 22)
#   )

# print(p1)

# # Gráfico 1: Estructura de edades por año
# p1.5 <- ggplot(poblacion_por_edad |> reframe(N_age = sum(N_estimado_por_edad)/1000, .by = c("year","age")), aes(x = factor(year), y = N_age, fill = factor(age))) +
#   geom_col(alpha = 0.8) +
#   scale_fill_viridis_d(name = "Edad") +
#   labs(
#     title = "Estructura de Edades por Año",
#     subtitle = "Población total estimada por Edad ",
#     x = "Año",
#     y = "Número de individuos (miles)"
#   ) +
#   theme_minimal() +
#   ylim(0, 300) +
#   theme(
#     plot.title = element_text(size = 14, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )

# print(p1.5)

# # Gráfico 2: Poblaciones por edad - años recientes
# # años_recientes <- tail(sort(unique(poblacion_por_edad$year)), 5)
# años_recientes <- sort(unique(poblacion_por_edad$year))

# p2 <- poblacion_por_edad |>
#   dplyr$filter(year %in% años_recientes) |>
#   ggplot(aes(x = factor(year), y = N_estimado_por_edad, fill = factor(age))) +
#   geom_col(alpha = 0.8) +
#   facet_wrap(~ Rio, scales = "free_y") +
#   scale_fill_viridis_d(name = "Edad") +
#   scale_y_continuous(labels = scales::comma) +
#   labs(
#     title = "Población Estimada por Edad",
#     subtitle = paste("Años:", paste(años_recientes, collapse = ", ")),
#     x = "Año",
#     y = "Número de individuos"
#   ) +
#   theme_minimal() +
#   ylim(0,100000)+
#   theme(
#     plot.title = element_text(size = 14, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.background = element_rect(fill = "lightgray")
#   )

# print(p2)

# # Gráfico 3: Composición por edades para un año específico
# año_ejemplo <- max(poblacion_por_edad$year)

# p3 <- poblacion_por_edad |>
#   dplyr$filter(year == año_ejemplo) |>
#   ggplot(aes(x = Rio, y = N_estimado_por_edad, fill = factor(age))) +
#   geom_col(alpha = 0.8) +
#   scale_fill_viridis_d(name = "Edad") +
#   scale_y_continuous(labels = scales::comma) +
#   labs(
#     title = paste("Composición por Edades -", año_ejemplo),
#     subtitle = "Distribución de edades por río",
#     x = "Río",
#     y = "Número de individuos"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold")
#   )

# print(p3)

# # ---- Guardar resultados ----

# print("\n=== GUARDANDO RESULTADOS ===")

# # Guardar archivo principal
# write.csv(poblacion_por_edad, "data/data_raw/river-population/age-structure-population-by-river.csv", row.names = FALSE)


































# library(janitor)

# df <- read.csv("data/data_raw/biological-data/data_chinook_cleaned.csv", 
#                fileEncoding = "latin1") |>
#   clean_names() |>
#   dplyr::filter(species == "Oncorhynchus tshawytscha",
#          !is.na(age), !is.na(year), !is.na(month))

# # Ver ubicaciones únicas con edad registrada
# unique(df$location)

# # Y con conteo por ubicación
# table(df$location)



