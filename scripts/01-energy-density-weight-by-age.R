# Energy density analysis for Chinook salmon
# NLS models: SSasymp and exponential by age and weight

# ---- Carga de paquetes ----
box::use(
  dplyr = dplyr[filter, select, case_when, mutate, reframe, n, arrange, left_join],
  janitor = janitor[clean_names],
  readr = readr[read_csv],
  tidyr = tidyr[pivot_longer, pivot_wider]
)

library(nlstools)
library(ggplot2)

# ---- Lectura de datos ----
densidad_energetica_salmones <- readr$read_csv("data/data_raw/energy-density/DE_salmones.csv")
factor_wet_dry <- read.csv("data/data_raw/energy-density/wet_dry_weight.csv")
biologico_salmones <- readr$read_csv("data/data_raw/biological-data/chinook_model_inputs.csv")

# ---- Procesamiento de datos ----

# Factor de conversión húmedo-seco
fc <- factor_wet_dry |> 
  janitor$clean_names() |> 
  dplyr$mutate(factor = (peso_humedo - peso_seco)/peso_humedo) |> 
  dplyr$reframe(factor = mean(factor))


# Datos de densidad energética procesados
densidad_energetica_resumen <- densidad_energetica_salmones |> 
  janitor$clean_names() |> 
  dplyr$mutate(
    de_humedo = de * (1 - as.numeric(fc$factor))
  ) |> 
  dplyr$reframe(
    de_humedo = mean(de_humedo),
    peso = mean(peso, na.rm = TRUE),
    talla = mean(talla, na.rm = TRUE),
    .by = c(individuo, age)
  )


densidad_energetica_resumen |> reframe(n = n(), peso_medio = mean(peso), sd = sd(peso), se = sd/sqrt(n), .by = age) |> as.data.frame()

tabla_export_de <- densidad_energetica_resumen  |> 
  dplyr$reframe(
    n_individuos = n(),
    DE_media = mean(de_humedo),
    peso_mean = mean(peso, na.rm = TRUE),
    talla_mean = mean(talla, na.rm = TRUE),
    DE_sd = sd(de_humedo),
    DE_se = DE_sd / sqrt(n_individuos),
    peso_sd = sd(peso),
    peso_se = peso_sd/sqrt(n_individuos),
    .by = c(age)
  ) |> as.data.frame()

write.csv(tabla_export_de, "data/data_raw/energy-density/de_weight_by_age.csv", row.names = FALSE)

# ---- Ajuste de los 4 modelos NLS ----


# Modelo 1: SSasymp por edad
modelo1_ssasymp_edad <- nls(de_humedo ~ SSasymp(age, Asym, R0, lrc),
                           data = densidad_energetica_resumen)

# Modelo 2: SSasymp por peso
modelo2_ssasymp_peso <- nls(de_humedo ~ SSasymp(peso, Asym, R0, lrc),
                           data = densidad_energetica_resumen)

# Modelo 3: Exponencial por edad
start_vals_edad <- list(a = mean(densidad_energetica_resumen$de_humedo), b = 1)
modelo3_exp_edad <- nls(de_humedo ~ a * age^b,
                       data = densidad_energetica_resumen,
                       start = start_vals_edad)

# Modelo 4: Exponencial por peso
start_vals_peso <- list(a = mean(densidad_energetica_resumen$de_humedo), b = 0.1)
modelo4_exp_peso <- nls(de_humedo ~ a * peso^b,
                       data = densidad_energetica_resumen,
                       start = start_vals_peso)

# Modelos lineales adicionales para comparar
modelo_linear_peso <- lm(de_humedo ~ peso, data = densidad_energetica_resumen)
modelo_linear_edad <- lm(de_humedo ~ age, data = densidad_energetica_resumen)
modelo_poly_peso <- lm(de_humedo ~ poly(peso, 2), data = densidad_energetica_resumen)
modelo_poly_edad <- lm(de_humedo ~ poly(age, 2), data = densidad_energetica_resumen)


# ---- Función para contar parámetros significativos ----
contar_params_sig <- function(modelo) {
  coefs <- summary(modelo)$coefficients
  sig <- sum(coefs[, "Pr(>|t|)"] < 0.05)
  total <- nrow(coefs)
  return(paste0(sig, "/", total))
}

# ---- Función para calcular R² ----
calc_r2 <- function(modelo, datos_reales) {
  return(cor(fitted(modelo), datos_reales)^2)
}

# ---- Función para R² de modelos lineales ----
get_r2_linear <- function(modelo) {
  return(summary(modelo)$r.squared)
}

# ---- Crear tabla resumen EXPANDIDA de modelos ----

tabla_resumen_modelos_expandida <- data.frame(
  Modelo = c(
    "Linear_edad", "Linear_peso", 
    "Polynomial_edad", "Polynomial_peso",
    "NLS_SSasymp_edad", "NLS_SSasymp_peso", 
    "NLS_exponencial_edad", "NLS_exponencial_peso"
  ),
  R2 = round(c(
    get_r2_linear(modelo_linear_edad),
    get_r2_linear(modelo_linear_peso),
    get_r2_linear(modelo_poly_edad),
    get_r2_linear(modelo_poly_peso),
    calc_r2(modelo1_ssasymp_edad, densidad_energetica_resumen$de_humedo),
    calc_r2(modelo2_ssasymp_peso, densidad_energetica_resumen$de_humedo),
    calc_r2(modelo3_exp_edad, densidad_energetica_resumen$de_humedo),
    calc_r2(modelo4_exp_peso, densidad_energetica_resumen$de_humedo)
  ), 3),
  RMSE_J_g = round(c(
    sqrt(mean(residuals(modelo_linear_edad)^2)),
    sqrt(mean(residuals(modelo_linear_peso)^2)),
    sqrt(mean(residuals(modelo_poly_edad)^2)),
    sqrt(mean(residuals(modelo_poly_peso)^2)),
    sqrt(mean(residuals(modelo1_ssasymp_edad)^2)),
    sqrt(mean(residuals(modelo2_ssasymp_peso)^2)),
    sqrt(mean(residuals(modelo3_exp_edad)^2)),
    sqrt(mean(residuals(modelo4_exp_peso)^2))
  ), 1),
  AIC = round(c(
    AIC(modelo_linear_edad),
    AIC(modelo_linear_peso),
    AIC(modelo_poly_edad),
    AIC(modelo_poly_peso),
    AIC(modelo1_ssasymp_edad),
    AIC(modelo2_ssasymp_peso),
    AIC(modelo3_exp_edad),
    AIC(modelo4_exp_peso)
  ), 1),
  BIC = round(c(
    BIC(modelo_linear_edad),
    BIC(modelo_linear_peso),
    BIC(modelo_poly_edad),
    BIC(modelo_poly_peso),
    BIC(modelo1_ssasymp_edad),
    BIC(modelo2_ssasymp_peso),
    BIC(modelo3_exp_edad),
    BIC(modelo4_exp_peso)
  ), 1),
  Params_sig = c(
    contar_params_sig(modelo_linear_edad),
    contar_params_sig(modelo_linear_peso),
    contar_params_sig(modelo_poly_edad),
    contar_params_sig(modelo_poly_peso),
    contar_params_sig(modelo1_ssasymp_edad),
    contar_params_sig(modelo2_ssasymp_peso),
    contar_params_sig(modelo3_exp_edad),
    contar_params_sig(modelo4_exp_peso)
  ),
  Ecuacion = c(
    "DE ~ a + b×edad",
    "DE ~ a + b×peso",
    "DE ~ a + b×edad + c×edad²",
    "DE ~ a + b×peso + c×peso²",
    "DE ~ SSasymp(edad)",
    "DE ~ SSasymp(peso)",
    "DE ~ a × edad^b",
    "DE ~ a × peso^b"
  ),
  Tipo = c(
    "Linear", "Linear", "Polynomial", "Polynomial",
    "NLS", "NLS", "NLS", "NLS"
  )
)


# ---- Identificar mejor modelo por categoría ----

# Mejor modelo lineal
modelos_lineales <- tabla_resumen_modelos_expandida[tabla_resumen_modelos_expandida$Tipo %in% c("Linear", "Polynomial"), ]
mejor_lineal_idx <- which.min(modelos_lineales$AIC)
mejor_modelo_lineal <- modelos_lineales$Modelo[mejor_lineal_idx]

# Mejor modelo NLS
modelos_nls <- tabla_resumen_modelos_expandida[tabla_resumen_modelos_expandida$Tipo == "NLS", ]
mejor_nls_idx <- which.min(modelos_nls$AIC)
mejor_modelo_nls <- modelos_nls$Modelo[mejor_nls_idx]

# Mejor modelo global
mejor_global_idx <- which.min(tabla_resumen_modelos_expandida$AIC)
mejor_modelo_global <- tabla_resumen_modelos_expandida$Modelo[mejor_global_idx]

# Mejor modelo por edad (solo NLS)
modelos_edad_nls <- modelos_nls[grepl("edad", modelos_nls$Modelo), ]
mejor_edad_idx <- which.min(modelos_edad_nls$AIC)
mejor_modelo_edad <- modelos_edad_nls$Modelo[mejor_edad_idx]

# Mejor modelo por peso (solo NLS)
modelos_peso_nls <- modelos_nls[grepl("peso", modelos_nls$Modelo), ]
mejor_peso_idx <- which.min(modelos_peso_nls$AIC)
mejor_modelo_peso <- modelos_peso_nls$Modelo[mejor_peso_idx]


# Asignar modelos seleccionados basándose en la selección automática
if(mejor_modelo_edad == "NLS_exponencial_edad") {
  modelo_edad_final <- modelo3_exp_edad
} else {
  modelo_edad_final <- modelo1_ssasymp_edad
}

if(mejor_modelo_peso == "NLS_exponencial_peso") {
  modelo_peso_final <- modelo4_exp_peso
} else {
  modelo_peso_final <- modelo2_ssasymp_peso
}

# ---- Obtener pesos promedio del archivo biológico ----


chinook_cleaned <- biologico_salmones |> 
  janitor$clean_names() |> 
  dplyr$select("age", "tw_g")


pesos_promedio_edad <- chinook_cleaned |> 
  dplyr$filter(!is.na(tw_g), !is.na(age)) |>
  dplyr$reframe(peso_promedio = mean(tw_g, na.rm = TRUE),
                n = dplyr$n(), 
                peso_se = sd(tw_g, na.rm = TRUE)/sqrt(n), 
                .by = "age") |> 
  dplyr$arrange(age) |> 
  as.data.frame()


pesos_completos <- pesos_promedio_edad

# ---- Funciones de predicción ----


predecir_por_edad <- function(edad, modelo) {
  return(predict(modelo, newdata = data.frame(age = edad)))
}

predecir_por_peso <- function(peso, modelo) {
  if(is.na(peso)) return(NA)
  return(predict(modelo, newdata = data.frame(peso = peso)))
}

# ---- Predicciones para edades 4 y 5 (modelo por peso) ----

pesos_4_5_valores <- pesos_completos$peso_promedio[pesos_completos$age %in% c(4, 5)]
if(length(pesos_4_5_valores) == 0) {
  pesos_4_5_valores <- c(12000, 15000)
  warning("No weights found for ages 4-5, using estimated values: ", paste(pesos_4_5_valores, collapse = ", "))
}

pred_peso_4_5 <- sapply(pesos_4_5_valores, function(x) predecir_por_peso(x, modelo_peso_final))



# ---- Gráficos ----

# Datos observados
datos_obs <- densidad_energetica_resumen |> 
  dplyr$select(age, de_humedo, peso) |> 
  dplyr$rename(edad = age, densidad = de_humedo)

# Ecuación del modelo por peso
if(mejor_modelo_peso == "NLS_exponencial_peso") {
  params_peso <- coef(modelo4_exp_peso)
  ecuacion_peso <- paste0("ED = ", round(params_peso['a'], 0), " × Weight^", round(params_peso['b'], 3))
} else {
  params_peso <- coef(modelo2_ssasymp_peso)
  ecuacion_peso <- paste0("ED = ", round(params_peso['Asym'], 0), " - (", round(params_peso['Asym'] - params_peso['R0'], 0), ") × exp(-exp(", round(params_peso['lrc'], 2), ") × Weight)")
}

peso_seq <- seq(min(datos_obs$peso, na.rm = TRUE), 
                max(c(pesos_completos$peso_promedio, pesos_4_5_valores), na.rm = TRUE), length.out = 100)
pred_peso_seq <- sapply(peso_seq, function(x) predecir_por_peso(x, modelo_peso_final))

puntos_extrap_peso <- data.frame(
  peso = pesos_4_5_valores,
  densidad = pred_peso_4_5
)

# p2: modelo por peso
p2 <- ggplot() +
  geom_point(data = datos_obs, aes(x = peso, y = densidad),
             color = "black", size = 3, alpha = 0.7) +
  annotate("rect", xmin = 9500, xmax = 15000, ymin = 5500, ymax = 8500,
           alpha = 0.15, fill = "purple") +
  geom_line(data = data.frame(peso = peso_seq, densidad = pred_peso_seq),
            aes(x = peso, y = densidad), color = "red", size = 1.5) +
  geom_point(data = puntos_extrap_peso, aes(x = peso, y = densidad),
             color = "blue", size = 4, shape = 15) +
  geom_text(data = puntos_extrap_peso,
            aes(x = peso, y = densidad,
                label = paste0(round(peso, 0), "g\n", round(densidad, 0), " J/g")),
            vjust = -1, hjust = 0.5, size = 6, fontface = "bold") +
  annotate("text", x = min(peso_seq) + 200, y = 8200, 
           label = ecuacion_peso, 
           size = 6, fontface = "bold", hjust = 0, vjust = 1,
           color = "black") +
  labs(title = paste("Best Weight Model:", mejor_modelo_peso),
       subtitle = "Blue squares = Predictions for ages 4-5 using extrapolated weights",
       x = "Weight (g)", y = "Energy Density (J/g)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 15, color = "black", face = "italic"),
    axis.text = element_text(size = 17, color = "black"),
    axis.title = element_text(size = 18, color = "black", face = "bold")
  ) +
  ylim(5500, 8500)

p2

# ---- Exportación de archivos ----

write.csv(tabla_resumen_modelos_expandida, "data/data_raw/energy-density/tabla_resumen_modelos_nls.csv", row.names = FALSE)




# ---- Incertidumbre predicciones edades 4-5 (Delta Method) ----

library(propagate)

peso_age4 <- pesos_4_5_valores[1]
peso_age5 <- pesos_4_5_valores[2]

pred_4 <- predictNLS(modelo4_exp_peso,
                     newdata = data.frame(peso = peso_age4),
                     interval = "confidence",
                     level = 0.95)

pred_5 <- predictNLS(modelo4_exp_peso,
                     newdata = data.frame(peso = peso_age5),
                     interval = "confidence",
                     level = 0.95)

se_delta <- data.frame(
  edad     = c(4, 5),
  peso     = c(peso_age4, peso_age5),
  ED_mean  = round(c(pred_4$summary[1, "Prop.Mean.1"],
                     pred_5$summary[1, "Prop.Mean.1"]), 0),
  ED_sd    = round(c(pred_4$summary[1, "Prop.sd.1"],
                     pred_5$summary[1, "Prop.sd.1"]), 0),
  CI_lower = round(c(pred_4$summary[1, "Prop.2.5%"],
                     pred_5$summary[1, "Prop.2.5%"]), 0),
  CI_upper = round(c(pred_4$summary[1, "Prop.97.5%"],
                     pred_5$summary[1, "Prop.97.5%"]), 0),
  Source   = "Predicted"
)

print(se_delta)
write.csv(se_delta, 
          "data/data_raw/energy-density/se_delta_pred_45.csv", 
          row.names = FALSE)


# ---- Figure S2.1: Energy density model fit ----

library(tidyverse)

datos_obs_resumen <- densidad_energetica_resumen %>%
  group_by(age) %>%
  summarise(
    peso_mean = mean(peso, na.rm = TRUE),
    de_mean   = mean(de_humedo, na.rm = TRUE),
    de_se     = sd(de_humedo, na.rm = TRUE) / sqrt(sum(!is.na(de_humedo))),
    .groups   = "drop"
  )

peso_seq <- seq(
  min(densidad_energetica_resumen$peso, na.rm = TRUE),
  max(pesos_4_5_valores) * 1.05,
  length.out = 200
)

pred_seq <- sapply(peso_seq, function(x) predecir_por_peso(x, modelo_peso_final))

curva <- data.frame(
  peso = peso_seq,
  de   = pred_seq
)

extrap <- data.frame(
  peso = pesos_4_5_valores,
  de   = pred_peso_4_5,
  age  = c(4, 5)
)

x_obs_max <- max(datos_obs_resumen$peso_mean, na.rm = TRUE)
x_extra_min <- min(pesos_4_5_valores, na.rm = TRUE) - 500
x_extra_max <- max(pesos_4_5_valores, na.rm = TRUE) + 500

# ---- Figura ----
fig_S2_1 <- ggplot() +

  annotate(
    "rect",
    xmin = x_extra_min,
    xmax = x_extra_max,
    ymin = 5500,
    ymax = 8500,
    alpha = 0.08,
    fill = "grey40"
  ) +

  geom_line(
    data = curva,
    aes(x = peso, y = de),
    linewidth = 1,
    color = "black"
  ) +

  geom_vline(
    xintercept = x_obs_max,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +

  geom_point(
    data = densidad_energetica_resumen,
    aes(x = peso, y = de_humedo),
    shape = 21,
    fill = "white",
    color = "black",
    size = 2.8,
    alpha = 0.7,
    stroke = 0.5
  ) +

  geom_errorbar(
    data = datos_obs_resumen,
    aes(
      x = peso_mean,
      ymin = de_mean - de_se,
      ymax = de_mean + de_se
    ),
    width = 100,
    linewidth = 0.8
  ) +

  geom_point(
    data = datos_obs_resumen,
    aes(x = peso_mean, y = de_mean),
    shape = 19,
    size = 3.5,
    color = "black"
  ) +

  geom_point(
    data = extrap,
    aes(x = peso, y = de),
    shape = 17,
    size = 3.8,
    color = "black"
  ) +

  geom_text(
    data = extrap,
    aes(
      x = peso,
      y = de,
      label = paste0("Age ", age)
    ),
    vjust = -1,
    size = 4
  ) +

  annotate(
    "text",
    x = min(densidad_energetica_resumen$peso, na.rm = TRUE) + 100,
    y = 8300,
    label = ecuacion_peso,
    hjust = 0,
    size = 4.5,
    fontface = "italic"
  ) +

  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.04))
  ) +

  scale_y_continuous(
    limits = c(5500, 8500),
    expand = expansion(mult = c(0, 0.02))
  ) +

  labs(
    x = "Body weight (g)",
    y = "Energy density (J g\u207B\u00B9 wet mass)"
  ) +

  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold", color = "black")
  )

print(fig_S2_1)

# ---- Exportar figura S2.1 ----
ggsave(
  "output/graficos/Figure_S2_1_energy_density_model.png",
  plot = fig_S2_1,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300
)
ggsave(
  "output/graficos/Figure_S2_1_energy_density_model.pdf",
  plot = fig_S2_1,
  width = 180,
  height = 120,
  units = "mm"
)
ggsave(
  "output/graficos/Figure_S2_1_energy_density_model.tiff",
  plot = fig_S2_1,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


