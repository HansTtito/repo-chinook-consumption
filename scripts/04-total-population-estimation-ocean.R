# ================================================================
# RETROCÁLCULO POBLACIONAL POR RÍO
# Estima población en el mar antes del retorno
# ================================================================

library(tidyverse)
library(ggplot2)

# ---- Cargar datos ----

print("=== CARGANDO DATOS ===")

# Población por edad, río y año (del script anterior)
poblacion_por_edad <- read.csv("data/data_raw/river-population/age-structure-population-by-river.csv")

print(paste("Datos cargados:", nrow(poblacion_por_edad), "registros"))
print(paste("Ríos:", paste(unique(poblacion_por_edad$Rio), collapse = ", ")))
print(paste("Años:", paste(range(poblacion_por_edad$year), collapse = "-")))

# ---- Parámetros ----

# Mortalidades naturales por edad
# 1.221069147
# 0.677896542
# 0.488438297
# 0.396434143
# 0.343833221
# 0.310838767

# mortalidades <- c(`1` = 1.221069147, `2` = 0.677896542, `3` = 0.488438297, `4` = 0.396434143, `5` = 0.343833221)

mortalidades <- c(`1` = 0.2, `2` = 0.2, `3` = 0.2, `4` = 0.2, `5` = 0.2)

print("Mortalidades naturales:")
print(mortalidades)

# ---- Función retrocálculo por río ----

calcular_retrocalculo_por_rio <- function(datos_poblacion, mortalidades) {
  
  # Preparar datos
  datos_base <- datos_poblacion %>%
    filter(!is.na(N_estimado_por_edad), N_estimado_por_edad > 0) %>%
    mutate(age = as.numeric(age))
  
  # Lista para almacenar resultados
  poblacion_mar_completa <- data.frame()
  
  # Para cada registro de retornantes
  for(i in 1:nrow(datos_base)) {
    rio <- datos_base$Rio[i]
    edad_retorno <- datos_base$age[i]
    anio_retorno <- datos_base$year[i]
    retornantes <- datos_base$N_estimado_por_edad[i]
    
    # Reconstruir cohorte hacia atrás
    poblacion_actual <- retornantes
    
    for(edad_actual in edad_retorno:1) {
      anio_actual <- anio_retorno - (edad_retorno - edad_actual)
      
      # Aplicar mortalidad hacia atrás (excepto en edad de retorno)
      if(edad_actual < edad_retorno) {
        M_edad <- mortalidades[as.character(edad_actual)]
        if(!is.na(M_edad)) {
          poblacion_actual <- poblacion_actual * exp(M_edad)
        }
      }
      
      # Guardar solo si NO es retornante (población en el mar)
      if(edad_actual < edad_retorno) {
        poblacion_mar_completa <- rbind(
          poblacion_mar_completa,
          data.frame(
            year = anio_actual,
            Rio = rio,
            age = edad_actual,
            N_mar = poblacion_actual,
            cohorte_origen = paste(anio_retorno, edad_retorno, sep = "_")
          )
        )
      }
    }
  }
  
  # Sumar contribuciones por año-río-edad
  poblacion_mar_final <- poblacion_mar_completa %>%
    group_by(year, Rio, age) %>%
    summarise(
      N_mar = sum(N_mar, na.rm = TRUE),
      n_cohortes = n_distinct(cohorte_origen),
      .groups = "drop"
    ) %>%
    arrange(Rio, year, age)
  
  return(poblacion_mar_final)
}

# ---- Ejecutar retrocálculo ----

print("\n=== EJECUTANDO RETROCÁLCULO ===")

poblacion_mar_por_rio <- calcular_retrocalculo_por_rio(poblacion_por_edad, mortalidades)

print(paste("Población en mar calculada:", nrow(poblacion_mar_por_rio), "registros"))

# ---- Crear tabla global ----

poblacion_mar_global <- poblacion_mar_por_rio %>%
  group_by(year, age) %>%
  summarise(
    N_mar_total = sum(N_mar, na.rm = TRUE),
    n_rios = n_distinct(Rio),
    .groups = "drop"
  ) %>%
  arrange(year, age)

print("\n=== RESUMEN POR RÍO ===")
resumen_por_rio <- poblacion_mar_por_rio %>%
  group_by(Rio) %>%
  summarise(
    años_disponibles = n_distinct(year),
    total_poblacion = sum(N_mar, na.rm = TRUE),
    poblacion_promedio_anual = round(total_poblacion / años_disponibles, 0),
    .groups = "drop"
  )
print(resumen_por_rio)

print("\n=== RESUMEN GLOBAL ===")
resumen_global <- poblacion_mar_global %>%
  group_by(year) %>%
  summarise(
    total_año = sum(N_mar_total, na.rm = TRUE),
    .groups = "drop"
  )
print(head(resumen_global))

# ---- Gráficos ----

print("\n=== GENERANDO GRÁFICOS ===")

# # Colores por edad (escala de grises)
# colores_edad <- c(
#   "1" = "grey90",
#   "2" = "grey70", 
#   "3" = "grey50",
#   "4" = "grey30",
#   "5" = "grey10"
# )

# # Colores azules oceánicos por edad
# colores_edad <- c(
#  "1" = "#E3F2FD",  # azul muy claro
#  "2" = "#90CAF9",  # azul claro
#  "3" = "#42A5F5",  # azul medio
#  "4" = "#1976D2",  # azul oscuro
#  "5" = "#0D47A1"   # azul muy oscuro
# )

# Opción 2: Viridis
colores_edad <- c(
 "1" = "#FDE725",  # amarillo claro
 "2" = "#5DC863",  # verde
 "3" = "#21908C",  # verde-azul
 "4" = "#3B528B",  # azul-morado
 "5" = "#440154"   # morado oscuro
)

# Gráfico 1: Por río
p1 <- poblacion_mar_por_rio %>%
 filter(year >= 2014 & year <= 2022) %>%
 mutate(age = factor(age)) %>%
 ggplot(aes(x = factor(year), y = N_mar, fill = age)) +
 geom_col(alpha = 0.9, color = "white", size = 0.2) +
 facet_wrap(~ Rio, scales = "free_y") +
 scale_fill_manual(values = colores_edad, name = "Age") +
 scale_y_continuous(labels = scales::comma) +
 labs(
  #  title = "At-sea Population by River",
  #  subtitle = "Backward cohort reconstruction based on estimated spawners",
   x = "",
   y = "At-sea population\n"
 ) +
 theme_classic() +
 theme(
   axis.text = element_text(size = 18, color = "black"),
   strip.background = element_rect(fill = "lightgray"),
   axis.title = element_text(size = 20),
   strip.text = element_text(size = 20),
   legend.title = element_text(size = 20),
   legend.text = element_text(size = 18)
 )

print(p1)

# Gráfico 2: Global
p2 <- poblacion_mar_global  %>%
 filter(year >= 2014 & year <= 2022) %>%
 mutate(age = factor(age)) %>%
  mutate(age = factor(age)) %>%
  ggplot(aes(x = factor(year), y = N_mar_total, fill = age)) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_edad, name = "Age") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    # title = "Población Total en el Mar (Todos los Ríos)",
    # subtitle = "Suma de retrocálculos de todos los ríos",
    x = "", 
    y = "At-sea population\n"
  ) +
 theme_classic() +
 theme(
   axis.text = element_text(size = 18, color = "black"),
   strip.background = element_rect(fill = "lightgray"),
   axis.title = element_text(size = 20),
   strip.text = element_text(size = 20),
   legend.title = element_text(size = 20),
   legend.text = element_text(size = 18)
 )
print(p2)

# ---- Guardar resultados ----

print("\n=== GUARDANDO RESULTADOS ===")


# Guardar tablas
write.csv(poblacion_mar_por_rio, "data/data_raw/river-population/poblacion_total_mar_por_rio.csv", row.names = FALSE)
write.csv(poblacion_mar_global, "data/data_raw/river-population/poblacion_mar_global.csv", row.names = FALSE)

print("Archivos guardados:")
print("- poblacion_total_mar_por_rio.csv")
print("- poblacion_mar_global.csv")

# ---- Resumen final ----

print("\n=== RESUMEN FINAL ===")
print("Población en el mar calculada por río y globalmente")
print("Parámetros usados:")
print(paste("- Mortalidades:", paste(names(mortalidades), "=", mortalidades, collapse = ", ")))
print(paste("- Ríos procesados:", length(unique(poblacion_mar_por_rio$Rio))))
print(paste("- Años:", paste(range(poblacion_mar_por_rio$year), collapse = "-")))
print(paste("- Total registros por río:", nrow(poblacion_mar_por_rio)))
print(paste("- Total registros globales:", nrow(poblacion_mar_global)))
