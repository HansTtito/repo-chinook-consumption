# Diet proportions for bioenergetic model — isotope-derived ages 1 and 2

library(tidyverse)
library(cosimmr)

# =============================================================================
# 1. CARGAR DATOS
# =============================================================================

model_muscle_5f <- readRDS("output/cosimmr_model_8_muscle_5sources_cov.rds")
chinook_diet_original <- read.csv("data/data_raw/biological-data/diet_proportion_by_age.csv")

peso_edad1 <- 700.0
peso_edad2 <- 3271.3

# =============================================================================
# 2. EXTRACT DIET PROPORTIONS BY BODY MASS
# =============================================================================

extraer_proporciones <- function(model, peso_objetivo, n_pred = 500) {
  p <- plot(model, type = 'covariates_plot', cov_name = 'Mass_g', 
            one_plot = TRUE, n_pred = n_pred)
  
  datos <- p$data
  pesos_unicos <- unique(datos$cov)
  idx_cercano <- which.min(abs(pesos_unicos - peso_objetivo))
  peso_cercano <- pesos_unicos[idx_cercano]
  
  resultado <- datos %>%
    filter(cov == peso_cercano) %>%
    select(Source, mean) %>%
    mutate(Proportion = round(mean, 4))
  
  props <- setNames(resultado$Proportion, resultado$Source)
  return(props)
}

props_edad1_5f <- extraer_proporciones(model_muscle_5f, peso_edad1)
props_edad2_5f <- extraer_proporciones(model_muscle_5f, peso_edad2)

cat("Proporciones isotópicas edad 1 (5 fuentes):\n")
print(props_edad1_5f)
cat("\nProporciones isotópicas edad 2 (5 fuentes):\n")
print(props_edad2_5f)

# =============================================================================
# 3. MAPEAR 5 FUENTES → 3 CATEGORÍAS BIOENERGÉTICAS
# =============================================================================

# Fuentes incluidas: Anchovy, Sardine, Zooplankton, Fish, Invertebrates
# Phytoplankton y Primary_producers NO incluidos (línea base isotópica únicamente)

mapear_5_a_3 <- function(props_5f) {
  prop_anchoveta <- as.numeric(props_5f["Anchovy"])
  prop_sardina   <- as.numeric(props_5f["Sardine"])
  prop_otros     <- as.numeric(props_5f["Zooplankton"]) + 
                    as.numeric(props_5f["Fish"]) + 
                    as.numeric(props_5f["Invertebrates"])
  
  total <- prop_anchoveta + prop_sardina + prop_otros
  
  return(c(
    sardina   = prop_sardina,
    anchoveta = prop_anchoveta,
    otros     = prop_otros,
    total     = total
  ))
}

props_edad1_3cat <- mapear_5_a_3(props_edad1_5f)
props_edad2_3cat <- mapear_5_a_3(props_edad2_5f)

cat("Edad 1 (3 categorías):\n")
print(round(props_edad1_3cat, 4))
cat("\nEdad 2 (3 categorías):\n")
print(round(props_edad2_3cat, 4))

# Validar
if(abs(props_edad1_3cat["total"] - 1.0) > 0.01) warning("Edad 1 no suma 1.0")
if(abs(props_edad2_3cat["total"] - 1.0) > 0.01) warning("Edad 2 no suma 1.0")

# =============================================================================
# 4. CREAR TABLA DE DIETA MODIFICADA
# =============================================================================

dieta_modificada <- chinook_diet_original

dieta_modificada$prop_sardina_mean[dieta_modificada$age == 1]    <- props_edad1_3cat["sardina"]
dieta_modificada$prop_anchoveta_mean[dieta_modificada$age == 1]  <- props_edad1_3cat["anchoveta"]
dieta_modificada$prop_otros_mean[dieta_modificada$age == 1]      <- props_edad1_3cat["otros"]
dieta_modificada$sardina_porcentaje[dieta_modificada$age == 1]   <- round(props_edad1_3cat["sardina"] * 100, 2)
dieta_modificada$anchoveta_porcentaje[dieta_modificada$age == 1] <- round(props_edad1_3cat["anchoveta"] * 100, 2)
dieta_modificada$otros_porcentaje[dieta_modificada$age == 1]     <- round(props_edad1_3cat["otros"] * 100, 2)

dieta_modificada$prop_sardina_mean[dieta_modificada$age == 2]    <- props_edad2_3cat["sardina"]
dieta_modificada$prop_anchoveta_mean[dieta_modificada$age == 2]  <- props_edad2_3cat["anchoveta"]
dieta_modificada$prop_otros_mean[dieta_modificada$age == 2]      <- props_edad2_3cat["otros"]
dieta_modificada$sardina_porcentaje[dieta_modificada$age == 2]   <- round(props_edad2_3cat["sardina"] * 100, 2)
dieta_modificada$anchoveta_porcentaje[dieta_modificada$age == 2] <- round(props_edad2_3cat["anchoveta"] * 100, 2)
dieta_modificada$otros_porcentaje[dieta_modificada$age == 2]     <- round(props_edad2_3cat["otros"] * 100, 2)

# Mostrar comparación
cat("\n--- COMPARACIÓN EDADES 1-2 ---\n\n")
cat("ORIGINAL (estómagos):\n")
print(chinook_diet_original[chinook_diet_original$age %in% 1:2, 
                            c("age", "prop_sardina_mean", "prop_anchoveta_mean", "prop_otros_mean")])

cat("\nMODIFICADA (isótopos músculo 5f):\n")
print(dieta_modificada[dieta_modificada$age %in% 1:2, 
                       c("age", "prop_sardina_mean", "prop_anchoveta_mean", "prop_otros_mean")])

# =============================================================================
# 5. GUARDAR TABLA
# (nombre mantenido para compatibilidad con scripts posteriores 08, 09, 10)
# =============================================================================

write.csv(dieta_modificada, 
          "data/data_raw/biological-data/diet_proportion_by_age_ISOTOPES.csv", 
          row.names = FALSE)


# =============================================================================
# 6. GRÁFICO COMPARATIVO
# =============================================================================

datos_original <- chinook_diet_original %>%
  filter(age %in% 1:5) %>%
  select(age, prop_sardina_mean, prop_anchoveta_mean, prop_otros_mean) %>%
  mutate(Scenario = "Original (Stomachs)")

datos_modificado <- dieta_modificada %>%
  filter(age %in% 1:5) %>%
  select(age, prop_sardina_mean, prop_anchoveta_mean, prop_otros_mean) %>%
  mutate(Scenario = "Modified (Isotopes)")

datos_comparacion <- bind_rows(datos_original, datos_modificado) %>%
  pivot_longer(cols = c(prop_sardina_mean, prop_anchoveta_mean, prop_otros_mean),
               names_to = "Prey", values_to = "Proportion") %>%
  mutate(
    Prey = case_when(
      Prey == "prop_sardina_mean"   ~ "Sardine",
      Prey == "prop_anchoveta_mean" ~ "Anchovy",
      Prey == "prop_otros_mean"     ~ "Others"
    ),
    Prey = factor(Prey, levels = c("Others", "Sardine", "Anchovy"))
  )

colores <- c("Sardine" = "#5B9BD5", "Anchovy" = "#ED7D31", "Others" = "#A5A5A5")

p_comparison <- ggplot(datos_comparacion, aes(x = factor(age), y = Proportion, fill = Prey)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
  facet_wrap(~Scenario, ncol = 2) +
  scale_fill_manual(values = colores) +
  scale_y_continuous(
    name = "Diet Proportion",
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  scale_x_discrete(name = "Age") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  labs(
    title = "Diet Composition Comparison: Original vs Isotope-Modified",
    subtitle = "Ages 1-2 modified using muscle 5-source isotope model (Anchovy, Sardine, Zooplankton, Fish, Invertebrates)",
    fill = "Prey Type"
  )

print(p_comparison)

ggsave("output/diet_comparison_isotopes_vs_stomachs.png", 
       plot = p_comparison, width = 12, height = 6, dpi = 300)
