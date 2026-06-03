# Backward cohort reconstruction — at-sea population by river

library(tidyverse)
library(ggplot2)

# ---- Cargar datos ----

# Población por edad, río y año (del script anterior)
poblacion_por_edad <- read.csv("data/data_raw/river-population/age-structure-population-by-river.csv")

# ---- Parámetros ----

mortalidades <- c(`1` = 0.2, `2` = 0.2, `3` = 0.2, `4` = 0.2, `5` = 0.2)

# ---- Función retrocálculo por río ----

calcular_retrocalculo_por_rio <- function(datos_poblacion, mortalidades) {
  
  datos_base <- datos_poblacion %>%
    filter(!is.na(N_estimado_por_edad), N_estimado_por_edad > 0) %>%
    mutate(age = as.numeric(age))
  
  poblacion_mar_completa <- data.frame()
  
  for(i in 1:nrow(datos_base)) {
    rio <- datos_base$Rio[i]
    edad_retorno <- datos_base$age[i]
    anio_retorno <- datos_base$year[i]
    retornantes <- datos_base$N_estimado_por_edad[i]
    
    poblacion_actual <- retornantes
    
    for(edad_actual in edad_retorno:1) {
      anio_actual <- anio_retorno - (edad_retorno - edad_actual)
      
      if(edad_actual < edad_retorno) {
        M_edad <- mortalidades[as.character(edad_actual)]
        if(!is.na(M_edad)) {
          poblacion_actual <- poblacion_actual * exp(M_edad)
        }
      }
      
      if(edad_actual < edad_retorno) {
        poblacion_mar_completa <- rbind(
          poblacion_mar_completa,
          data.frame(
            year = anio_actual,
            Rio = rio,
            age = edad_actual,
            N_mar = poblacion_actual
          )
        )
      }
    }
  }
  
  poblacion_mar_final <- poblacion_mar_completa %>%
    group_by(year, Rio, age) %>%
    summarise(
      N_mar = sum(N_mar, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Rio, year, age)
  
  return(poblacion_mar_final)
}

# ---- Ejecutar retrocálculo ----

poblacion_mar_por_rio <- calcular_retrocalculo_por_rio(poblacion_por_edad, mortalidades)

# ---- Crear tabla global ----

poblacion_mar_global <- poblacion_mar_por_rio %>%
  group_by(year, age) %>%
  summarise(
    N_mar_total = sum(N_mar, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year, age)

# ---- Gráficos ----

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
  ggplot(aes(x = factor(year), y = N_mar_total, fill = age)) +
  geom_col(alpha = 0.9, color = "white", size = 0.2) +
  scale_fill_manual(values = colores_edad, name = "Age") +
  scale_y_continuous(labels = scales::comma) +
  labs(
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

write.csv(poblacion_mar_por_rio, "data/data_raw/river-population/poblacion_total_mar_por_rio.csv", row.names = FALSE)
write.csv(poblacion_mar_global, "data/data_raw/river-population/poblacion_mar_global.csv", row.names = FALSE)
