library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(cowplot)

# ---------------------------
# DATOS
# ---------------------------

poligono_shp <- st_read("data/spatial/poligono_costero_60mn.shp")

river_points <- read.csv(
  "data/spatial/river_points.csv",
  fileEncoding = "latin1"
) %>%
  filter(Puerto %in% c("Imperial","Toltén","Valdivia","Río bueno")) %>%
  mutate(marker = if_else(marker == "Study",
                          "Direct census",
                          "Estimated population"))

world <- ne_countries(scale = "large", returnclass = "sf")

# ---------------------------
# MAPA PRINCIPAL (B/N)
# ---------------------------

main_map <- ggplot() +

  # Fondo mundial tenue
  geom_sf(data = world,
          fill = "grey92",
          color = "grey60",
          linewidth = 0.2) +

  # Área de estudio 60 nm (sombreado)
  geom_sf(data = poligono_shp,
          fill = "grey70",
          alpha = 0.25,
          color = "black",
          linewidth = 0.7) +

  # Puntos ríos
  geom_point(
    data = river_points,
    aes(Longitud, Latitud,
        shape = marker,
        fill = marker),
    size = 3,
    stroke = 0.8,
    color = "black"
  ) +

  # Etiquetas
  geom_label(
    data = river_points,
    aes(Longitud, Latitud, label = Puerto),
    nudge_x = c(1.2, 1.2, 1.5, 1.2),
    nudge_y = c(0, 0, 0, -0.3),
    size = 2.8,
    label.size = 0,
    fill = "white",
    alpha = 0.9
  ) +

  # Escala
  annotation_scale(
    location = "br",
    width_hint = 0.25,
    height = unit(0.15, "cm")
  ) +

  # Norte
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(
      line_col = "black",
      fill = c("white", "black")
    )
  ) +

  coord_sf(
    xlim = c(-83, -60),
    ylim = c(-44, -31),
    expand = FALSE
  ) +

  scale_shape_manual(
    name = "River population data",
    values = c(
      "Direct census" = 21,
      "Estimated population" = 24
    )
  ) +

  scale_fill_manual(
    name = "River population data",
    values = c(
      "Direct census" = "black",
      "Estimated population" = "white"
    )
  ) +

  labs(x = "Longitude", y = "Latitude") +

  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = c(0.83, 0.83),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

inset_map <- ggplot(world) +
  geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_rect(
    aes(xmin = -83, xmax = -60,
        ymin = -44, ymax = -31),
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  coord_sf(
    xlim = c(-100, -20),
    ylim = c(-60, 15),
    expand = FALSE
  ) +
  theme_void()

final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map,
            x = 0.06,
            y = 0.66,
            width = 0.23,
            height = 0.23)

final_map

# ---------------------------
# EXPORTAR
# ---------------------------

# Opción 1: Para revistas - formato estándar (recomendado)
# Single column: ancho ~3.5", two-column: ~7"
ggsave("outputs/graficos/Study_area_final.pdf",
       final_map,
       width = 180, height = 140, units = "mm",  # ~7" x 5.5" en mm
       dpi = 300)  # 300 dpi es suficiente para PDF vectorial

# Opción 2: Si la revista pide tamaño específico (ejemplo Nature/Science)
ggsave("outputs/graficos/Study_area_final.pdf",
       final_map,
       width = 183, height = 130, units = "mm",  # Full page width común
       dpi = 300)

# Opción 3: Para presentaciones o posters (más grande)
ggsave("outputs/graficos/Study_area_final.pdf",
       final_map,
       width = 10, height = 7, units = "in",
       dpi = 300)

# Opción 4: Guardar TAMBIÉN en PNG de alta calidad (backup)
ggsave("outputs/graficos/Study_area_final.png",
       final_map,
       width = 180, height = 140, units = "mm",
       dpi = 600,  # Para PNG sí importa más el DPI
       bg = "white")




# Diet -------------------------------------------------------------------

library(tidyverse)
library(ggpattern)  

# Leer datos
diet_data <- read.csv("data/data_raw/biological-data/diet_proportion_by_age_ISOTOPES.csv")

# Preparar datos para gráfico
plot_data <- diet_data %>%
  filter(age %in% 1:5) %>%
  select(age, prop_sardina_mean, prop_anchoveta_mean, prop_otros_mean) %>%
  pivot_longer(
    cols = c(prop_sardina_mean, prop_anchoveta_mean, prop_otros_mean),
    names_to = "Prey",
    values_to = "Proportion"
  ) %>%
  mutate(
    Prey = case_when(
      Prey == "prop_sardina_mean" ~ "Sardine",
      Prey == "prop_anchoveta_mean" ~ "Anchovy",
      Prey == "prop_otros_mean" ~ "Others"
    ),
    Prey = factor(Prey, levels = c("Others", "Sardine", "Anchovy")),
    Proportion = Proportion * 100
  )

plot_data <- plot_data %>%
  mutate(Prey = factor(Prey,
                       levels = c("Others", "Anchovy", "Sardine")))


# Gráfico en escala de grises
p <- ggplot(plot_data,
            aes(x = factor(age),
                y = Proportion,
                fill = Prey)) +

  geom_bar(stat = "identity",
           position = "stack",
           color = "black",
           linewidth = 0.25) +

  scale_fill_manual(values = c(
    "Others" = "white",
    "Anchovy" = "grey60",
    "Sardine" = "grey20"
  )) +

  scale_y_continuous(
    name = "Diet proportion (%)",
    limits = c(0, 101),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +

  scale_x_discrete(name = "Age (years)") +

  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +

  labs(fill = "Prey type") 

# Mostrar
print(p)

# Guardar en alta resolución
ggsave("outputs/graficos/diet_composition_grayscale.pdf",
       plot = p,
       width = 180, height = 120, units = "mm",
       dpi = 300)

# También en PNG
ggsave("outputs/graficos/diet_composition_grayscale.png",
       plot = p,
       width = 180, height = 120, units = "mm",
       dpi = 600,
       bg = "white")

# Tabla resumen
cat("\n=== RESUMEN POR EDAD ===\n")
plot_data %>%
  pivot_wider(names_from = Prey, values_from = Proportion) %>%
  arrange(age) %>%
  print()



# Isotopes ---------------------------------------------------------------

library(tidyverse)

# =============================================================================
# BIPLOT COMPLETO Y CLARO: TODAS LAS FUENTES VISIBLES
# =============================================================================

# [Mismo código de preparación de datos hasta sources_all]
raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")

names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-") ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^h-") ~ "Salmon_liver", 
  str_detect(raw_data$Sample, "^Sar") ~ "Sardine",
  str_detect(raw_data$Sample, "^Zoo") ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito") ~ "Phytoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchovy",
  TRUE ~ "Unknown"
)

lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent <- 93 / (1 + (1/lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  final_correction <- ifelse(CN_ratio > 1, correction_factor, 0)
  return(d13C + final_correction)
}

raw_data <- raw_data %>%
  mutate(d13C_corrected = case_when(
    group %in% c("Salmon_muscle", "Salmon_liver") ~ lipid_correction_optimized_chinook(d13C, CtoN),
    TRUE ~ d13C
  ))

consumers <- raw_data %>%
  filter(group == "Salmon_muscle", !is.na(d13C_corrected), !is.na(d15N), !is.na(age), age %in% 1:5) %>%
  select(ID, Sample, age, d13C_corrected, d15N, d34S)

sources_marine <- raw_data %>%
  filter(group %in% c("Sardine", "Anchovy", "Zooplankton", "Phytoplankton"), 
         !is.na(d13C_corrected), !is.na(d15N)) %>%
  mutate(source_name = recode(group, 
    "Sardine" = "Sardine", "Anchovy" = "Anchovy",
    "Zooplankton" = "Zooplankton", "Phytoplankton" = "Phytoplankton"),
    habitat = "Marine") %>%
  select(source_name, d13C_corrected, d15N, habitat)

sources_freshwater <- read.csv("data/data_raw/diet-isotopes/datos_rios.csv") %>%
  rename(source_name = Source) %>%
  mutate(d13C_corrected = d13C_mean, d15N = d15N_mean, habitat = "Freshwater") %>%
  select(source_name, d13C_corrected, d15N, habitat)

sources_all <- bind_rows(sources_marine, sources_freshwater)

# Calcular elipses y centroides
calc_ellipse <- function(data, level = 0.95) {
  if(nrow(data) < 3) return(NULL)
  sigma <- cov(data[, c("d13C_corrected", "d15N")])
  mu <- colMeans(data[, c("d13C_corrected", "d15N")])
  eigenvalues <- eigen(sigma)$values
  eigenvectors <- eigen(sigma)$vectors
  chisq_val <- qchisq(level, df = 2)
  radii <- sqrt(chisq_val * eigenvalues)
  angle <- seq(0, 2 * pi, length.out = 100)
  data.frame(
    x = mu[1] + radii[1] * cos(angle) * eigenvectors[1, 1] + radii[2] * sin(angle) * eigenvectors[1, 2],
    y = mu[2] + radii[1] * cos(angle) * eigenvectors[2, 1] + radii[2] * sin(angle) * eigenvectors[2, 2]
  )
}

ellipses_list <- consumers %>%
  group_by(age) %>%
  group_split() %>%
  map(~{ ellipse <- calc_ellipse(.x); if(!is.null(ellipse)) ellipse$age <- .x$age[1]; ellipse }) %>%
  bind_rows()

centroides <- consumers %>%
  group_by(age) %>%
  summarise(d13C_mean = mean(d13C_corrected), d15N_mean = mean(d15N), .groups = 'drop')

# =============================================================================
# PREPARAR FUENTES CON SÍMBOLOS Y COLORES DIFERENCIADOS
# =============================================================================

age_greys <- c("1" = "grey85", "2" = "grey65", "3" = "grey45", "4" = "grey25", "5" = "grey5")

# Símbolos y configuración por fuente
source_config <- data.frame(
  source_name = c("Sardine", "Anchovy", "Zooplankton", "Phytoplankton", 
                  "Fish", "Invertebrates", "Primary_producers"),
  shape = c(15, 17, 19, 18, 4, 3, 8),  # Diferentes símbolos
  color = c("grey70", "grey70", "grey70", "grey70",  # Marinas en negro/gris oscuro
            "grey10", "grey10", "grey10"),         # Agua dulce en gris
  size = c(4.5, 4.5, 4.5, 4.5, 5, 5, 5),  # Agua dulce un poco más grandes
  habitat = c("Marine", "Marine", "Marine", "Marine",
              "Freshwater", "Freshwater", "Freshwater")
)

# Combinar con datos isotópicos
sources_plot <- sources_all %>%
  left_join(source_config, by = c("source_name","habitat"))

# Calcular posiciones promedio para etiquetas (para múltiples puntos por fuente)
source_labels <- sources_plot %>%
  group_by(source_name, shape, color, size, habitat) %>%
  summarise(
    d13C = mean(d13C_corrected),
    d15N = mean(d15N),
    .groups = 'drop'
  ) %>%
  mutate(
    # Ajustar posición de etiquetas para que no tapen símbolos
    nudge_x = case_when(
      source_name == "Sardine" ~ 0.5,
      source_name == "Anchovy" ~ 1.3,
      source_name == "Zooplankton" ~ 0.8,
      source_name == "Phytoplankton" ~ 0.8,
      TRUE ~ 0.5
    ),
    nudge_y = case_when(
      source_name == "Sardine" ~ -0.5,
      source_name == "Anchovy" ~ 0.5,
      source_name == "Zooplankton" ~ -0.5,
      source_name == "Phytoplankton" ~ 0.3,
      TRUE ~ 0.5
    )
  )

# =============================================================================
# GRÁFICO FINAL
# =============================================================================

p_biplot <- ggplot() +
  # Elipses por edad
  geom_path(data = ellipses_list,
            aes(x = x, y = y, group = age, color = factor(age)),
            linewidth = 1.2) +
  # Puntos individuales salmón (muy transparentes)
  geom_point(data = consumers,
             aes(x = d13C_corrected, y = d15N, fill = factor(age)),
             shape = 21, size = 2, color = "black", stroke = 0.3, alpha = 0.25) +
  # TODAS las fuentes - puntos individuales
  geom_point(data = sources_plot,
             aes(x = d13C_corrected, y = d15N, shape = shape),
             size = sources_plot$size,
             color = sources_plot$color,
             stroke = 1) +
  # Etiquetas de fuentes (ajustadas para no tapar)
  geom_text(data = source_labels,
            aes(x = d13C + nudge_x, y = d15N + nudge_y, 
                label = source_name,
                fontface = ifelse(habitat == "Marine", "bold", "italic")),
            size = 3, hjust = 0) +
  # Centroides GRANDES con números
  geom_point(data = centroides,
             aes(x = d13C_mean, y = d15N_mean, fill = factor(age)),
             size = 12, shape = 21, color = "black", stroke = 1.5) +
  geom_text(data = centroides,
            aes(x = d13C_mean, y = d15N_mean, label = age),
            size = 6, fontface = "bold", color = "white") +
  # Escalas
  scale_fill_manual(name = "Salmon Age", values = age_greys) +
  scale_color_manual(name = "Salmon Age", values = age_greys) +
  scale_shape_identity() +  # Usar símbolos directamente
  # Ejes
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030"))
  ) +
  # Tema
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = c(0.12, 0.88),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white")
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 6, alpha = 1)),
    color = "none"
  )

print(p_biplot)

ggsave("outputs/graficos/isotope_biplot_all_sources.pdf",
       plot = p_biplot, width = 210, height = 170, units = "mm", dpi = 300)

ggsave("outputs/graficos/isotope_biplot_all_sources.png",
       plot = p_biplot, width = 210, height = 170, units = "mm", dpi = 600, bg = "white")




# Consumption by age -----------------------------------------------------

# Paleta de grises para las transiciones
GRISES_TRANSICION <- c(
  "1 → 2" = "gray90",
  "2 → 3" = "gray65", 
  "3 → 4" = "gray40",
  "4 → 5" = "gray15"
)

# Boxplot en escala de grises
plot_consumo_boxplot_por_transicion <- function(consumo_individual, titulo = NULL, verbose = TRUE) {
  if(verbose) cat("=== CREANDO BOXPLOT DE CONSUMO POR TRANSICIÓN ===\n")
  
  p <- ggplot(consumo_individual, aes(x = transicion_label, y = consumption_kg, fill = transicion_label)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = GRISES_TRANSICION) +
    scale_y_continuous(
      name = "Annual prey consumption (kg)",
      labels = scales::comma_format()
    ) +
    scale_x_discrete(name = "Age Transition") +
    # theme_minimal() +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5)
    ) 
    # labs(
    #   title = titulo %||% "Distribution of Annual Consumption by Age Transition",
    #   subtitle = "Boxplots show median, quartiles and individual data points across years"
    # )
  
  return(p)
}

# Barras en escala de grises
plot_consumo_barras_promedio_transicion <- function(consumo_individual, titulo = NULL, verbose = TRUE) {
  if(verbose) cat("=== CREANDO BARRAS DE CONSUMO PROMEDIO POR TRANSICIÓN ===\n")
  
  datos_resumen <- consumo_individual %>%
    group_by(transicion_label, age_inicial) %>%
    summarise(
      consumo_promedio = mean(consumption_kg, na.rm = TRUE),
      consumo_se = sd(consumption_kg, na.rm = TRUE) / sqrt(n()),
      n_observaciones = n(),
      .groups = 'drop'
    ) %>%
    arrange(age_inicial)
  
  p <- ggplot(datos_resumen, aes(x = transicion_label, y = consumo_promedio, fill = transicion_label)) +
    geom_col(alpha = 0.8, color = "black", size = 0.5) +
    geom_errorbar(aes(ymin = consumo_promedio - consumo_se, 
                      ymax = consumo_promedio + consumo_se),
                  width = 0.3, alpha = 0.7, size = 0.8) +
    geom_text(aes(label = paste("n =", n_observaciones)),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = GRISES_TRANSICION) +
    scale_y_continuous(
      name = "Mean Annual Consumption (kg)",
      labels = scales::comma_format(),
      expand = expansion(mult = c(0, 0.15))
    ) +
    scale_x_discrete(name = "Age Transition") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5)
    ) +
    labs(
      title = titulo %||% "Mean Annual Consumption by Age Transition",
      subtitle = "Error bars show standard error; numbers indicate sample size"
    )
  
  if(verbose) {
    cat("Consumo promedio por transición (kg):\n")
    print(round(datos_resumen, 1))
  }
  
  return(p)
}

# Violin en escala de grises
plot_consumo_violin_por_transicion <- function(consumo_individual, titulo = NULL, verbose = TRUE) {
  if(verbose) cat("=== CREANDO VIOLIN PLOT DE CONSUMO POR TRANSICIÓN ===\n")
  
  p <- ggplot(consumo_individual, aes(x = transicion_label, y = consumption_kg, fill = transicion_label)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                 fill = "white", color = "black", stroke = 1) +
    scale_fill_manual(values = GRISES_TRANSICION) +
    scale_y_continuous(
      name = "Annual Consumption (kg)",
      labels = scales::comma_format()
    ) +
    scale_x_discrete(name = "Age Transition") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5)
    ) +
    labs(
      title = titulo %||% "Distribution Density of Annual Consumption by Age Transition",
      subtitle = "Violin plots show distribution density; diamonds show means; boxes show quartiles"
    )
  
  return(p)
}

# Líneas en escala de grises
plot_consumo_lineas_por_ano_transicion <- function(consumo_individual, titulo = NULL, verbose = TRUE) {
  if(verbose) cat("=== CREANDO GRÁFICO DE LÍNEAS POR AÑO Y TRANSICIÓN ===\n")
  
  p <- ggplot(consumo_individual, aes(x = year, y = consumption_kg, 
                                       color = transicion_label, 
                                       shape = transicion_label,
                                       linetype = transicion_label)) +
    geom_line(aes(group = transicion_label), size = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(name = "Age Transition", values = GRISES_TRANSICION) +
    scale_shape_manual(
      name = "Age Transition",
      values = c("1 → 2" = 19, "2 → 3" = 17, "3 → 4" = 15, "4 → 5" = 18)
    ) +
    scale_linetype_manual(
      name = "Age Transition",
      values = c("1 → 2" = "solid", "2 → 3" = "dashed", 
                 "3 → 4" = "dotted", "4 → 5" = "dotdash")
    ) +
    scale_x_continuous(
      name = "Year",
      breaks = scales::pretty_breaks(n = 8)
    ) +
    scale_y_continuous(
      name = "Annual Consumption (kg)",
      labels = scales::comma_format()
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5)
    ) +
    labs(
      title = titulo %||% "Temporal Trends in Annual Consumption by Age Transition",
      subtitle = "Lines show consumption patterns across years for each age transition"
    )
  
  return(p)
}

# Panel combinado (sin cambios, usa las funciones actualizadas)
plot_panel_consumo_edades <- function(consumo_individual, ncol = 2, titulo_general = NULL) {
  cat("=== CREANDO PANEL COMBINADO DE GRÁFICOS ===\n")
  
  p1 <- plot_consumo_boxplot_por_transicion(consumo_individual, "A) Distribution by Transition", verbose = FALSE)
  p2 <- plot_consumo_barras_promedio_transicion(consumo_individual, "B) Mean Consumption", verbose = FALSE)
  p3 <- plot_consumo_violin_por_transicion(consumo_individual, "C) Density Distribution", verbose = FALSE)
  p4 <- plot_consumo_lineas_por_ano_transicion(consumo_individual, "D) Temporal Trends", verbose = FALSE)
  
  panel_combined <- wrap_plots(p1, p2, p3, p4, ncol = ncol) +
    plot_annotation(
      title = titulo_general %||% "Comprehensive Analysis of Annual Consumption by Age Transition",
      subtitle = "Multiple perspectives on consumption patterns across age transitions and years",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40")
      )
    )
  
  return(panel_combined)
}

generar_tabla_consumo_individual <- function(resultados_modelos, verbose = TRUE) {
  if(verbose) cat("=== GENERANDO TABLA CONSUMO INDIVIDUAL ===\n")
  
  # Obtener métricas básicas
  modelos_metrics <- compare_scenarios(resultados_modelos$modelos, metrics = c("consumption", "p_value"))
  
  # Procesar datos básicos
  datos_base <- modelos_metrics$scenario_data %>%
    filter(converged == TRUE) %>%
    mutate(
      year = as.numeric(gsub("T(\\d{4})_.*", "\\1", scenario)),
      transicion = gsub("T\\d{4}_(.+)", "\\1", scenario),
      age_inicial = as.numeric(gsub("(\\d+)to\\d+", "\\1", transicion)),
      age_final = as.numeric(gsub("\\d+to(\\d+)", "\\1", transicion)),
      transicion_label = paste(age_inicial, "→", age_final),
      consumption_kg = consumption_est / 1000,
      consumption_se_kg = ifelse(!is.na(consumption_se), consumption_se / 1000, NA)
    )
  
  # Extraer metadatos de consumo por presa CON intervalos de confianza
  modelos_metadata <- data.frame()
  for(i in 1:nrow(datos_base)) {
    scenario_name <- datos_base$scenario[i]
    if(scenario_name %in% names(resultados_modelos$modelos)) {
      metadata <- resultados_modelos$modelos[[scenario_name]]$metadata
      
      fila_metadata <- data.frame(
        scenario = scenario_name,
        data_source = metadata$data_source %||% NA,
        consumo_min_kg = metadata$consumo_min_kg %||% NA,
        consumo_max_kg = metadata$consumo_max_kg %||% NA,
        consumo_se_kg = metadata$consumo_se_kg %||% NA,
        consumo_anchoveta_kg = metadata$consumo_anchoveta_kg %||% NA,
        consumo_anchoveta_min = metadata$consumo_anchoveta_min %||% NA,
        consumo_anchoveta_max = metadata$consumo_anchoveta_max %||% NA,
        consumo_sardina_kg = metadata$consumo_sardina_kg %||% NA,
        consumo_sardina_min = metadata$consumo_sardina_min %||% NA,
        consumo_sardina_max = metadata$consumo_sardina_max %||% NA,
        consumo_otros_kg = metadata$consumo_otros_kg %||% NA,
        consumo_otros_min = metadata$consumo_otros_min %||% NA,
        consumo_otros_max = metadata$consumo_otros_max %||% NA,
        prop_anchoveta = metadata$prop_anchoveta %||% NA,
        prop_sardina = metadata$prop_sardina %||% NA,
        prop_otros = metadata$prop_otros %||% NA,
        stringsAsFactors = FALSE
      )
      modelos_metadata <- rbind(modelos_metadata, fila_metadata)
    }
  }
  
  # Combinar datos CON intervalos de confianza completos
  consumo_individual <- datos_base %>%
    left_join(modelos_metadata, by = "scenario") %>%
    select(scenario, year, transicion_label, age_inicial, age_final,
          consumption_kg, consumption_se_kg,
          consumo_min_kg, consumo_max_kg,
          p_value_est, p_value_se, data_source,
          consumo_anchoveta_kg, consumo_anchoveta_min, consumo_anchoveta_max,
          consumo_sardina_kg, consumo_sardina_min, consumo_sardina_max,
          consumo_otros_kg, consumo_otros_min, consumo_otros_max,
          prop_anchoveta, prop_sardina, prop_otros) %>%
    arrange(year, age_inicial)
  
  if(verbose) {
    cat("Tabla final:", nrow(consumo_individual), "filas\n")
    cat("Años:", paste(sort(unique(consumo_individual$year)), collapse = ", "), "\n\n")
  }
  
  return(consumo_individual)
}


resultados_modelos <- readRDS("data/data_raw/bioenergetic-model/modelos_total_by_age_ISOTOPES.rds")

consumo_individual <- generar_tabla_consumo_individual(resultados_modelos)


individual_consumption_plot <- plot_consumo_boxplot_por_transicion(consumo_individual)


ggsave("outputs/graficos/individual_consumption_plot.pdf",
       plot = individual_consumption_plot, width = 210, height = 170, units = "mm", dpi = 300)

ggsave("outputs/graficos/individual_consumption_plot.png",
       plot = individual_consumption_plot, width = 210, height = 170, units = "mm", dpi = 600, bg = "white")






# ================================================================
# FIGURE 5 — Reconstructed at-sea abundance by age
# ================================================================

library(tidyverse)

# ---- Leer archivo generado previamente ----
poblacion_mar_global <- read.csv(
  "data/data_raw/river-population/poblacion_mar_global.csv"
)

# Filtrar período del paper
datos_fig5 <- poblacion_mar_global %>%
  filter(year >= 2014 & year <= 2022) %>%
  mutate(
    age = factor(age, levels = sort(unique(age)))
  )

# Paleta B/N con contraste claro (print-friendly)
colores_edad_bw <- c(
  "1" = "grey85",
  "2" = "grey65",
  "3" = "grey45",
  "4" = "grey25",
  "5" = "black"
)

# ---- Crear figura ----
fig5 <- ggplot(datos_fig5,
               aes(x = factor(year),
                   y = N_mar_total,
                   fill = age)) +

  geom_col(color = "black", linewidth = 0.3) +

  scale_fill_manual(
    values = colores_edad_bw,
    name = "Age class",
    breaks = c("1","2","3","4","5")
  ) +

  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +

  labs(
    x = NULL,
    y = "Reconstructed at-sea abundance"
  ) +

  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +

  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

print(fig5)

# ---- Guardar en formato revista ----
ggsave(
  "outputs/graficos/Figure_5_at_sea_abundance.tiff",
  plot = fig5,
  dpi = 300,
  width = 14,
  height = 10
)



# ================================================================
# ESTADÍSTICAS PARA RESULTADOS - POBLACIÓN AT-SEA
# ================================================================

# Cargar datos de población en el mar
poblacion_global <- read.csv("data/data_raw/river-population/poblacion_mar_global.csv") |> filter(year >= 2014)

cat("=== ABUNDANCIA AT-SEA TOTAL POR AÑO ===\n")
resumen_anual <- poblacion_global %>%
  group_by(year) %>%
  summarise(N_total = sum(N_mar_total, na.rm = TRUE)) %>%
  arrange(year)
print(resumen_anual)

cat("\n=== RANGO TOTAL ===\n")
cat(sprintf("Mínimo: %s individuos (año %d)\n",
    scales::comma(min(resumen_anual$N_total)),
    resumen_anual$year[which.min(resumen_anual$N_total)]))
cat(sprintf("Máximo: %s individuos (año %d)\n",
    scales::comma(max(resumen_anual$N_total)),
    resumen_anual$year[which.max(resumen_anual$N_total)]))
cat(sprintf("Fold-difference: %.0f-fold\n",
    max(resumen_anual$N_total) / min(resumen_anual$N_total)))

cat("\n=== ABUNDANCIA POR EDAD Y AÑO ===\n")
resumen_edad <- poblacion_global %>%
  group_by(year, age) %>%
  summarise(N = sum(N_mar_total, na.rm = TRUE), .groups = "drop") %>%
  arrange(year, age)
print(resumen_edad)

cat("\n=== PICO POR EDAD ===\n")
for(a in 1:4) {
  d <- resumen_edad %>% filter(age == a)
  cat(sprintf("Age-%d: máximo %s individuos en %d\n",
      a, scales::comma(max(d$N)), d$year[which.max(d$N)]))
}

cat("\n=== DOMINANCIA POR AÑO (edad dominante) ===\n")
resumen_edad %>%
  group_by(year) %>%
  slice_max(N, n = 1) %>%
  select(year, age, N) %>%
  print()




# Total Consumption ------------------------------------------------------


library(tidyverse)

# Cargar datos
poblacion_global <- read.csv("data/data_raw/river-population/poblacion_mar_global.csv")
poblacion_por_rio <- read.csv("data/data_raw/river-population/poblacion_total_mar_por_rio.csv")
consumo_individual <- read.csv("data/data_raw/bioenergetic-model/resultados_total_consumption_by_age_ISOTOPES.csv")

# Verificar overlap temporal
años_poblacion <- sort(unique(poblacion_global$year))
años_consumo <- sort(unique(consumo_individual$year)) - 1
años_overlap <- intersect(años_poblacion, años_consumo)

# Filtrar datos
poblacion_global_filtrada <- poblacion_global %>% filter(year %in% años_overlap)
poblacion_rio_filtrada <- poblacion_por_rio %>% filter(year %in% años_overlap)
consumo_filtrado <-consumo_individual %>%
  mutate(year_consumo = year - 1) %>%
  filter(year_consumo %in% años_overlap)
# Función auxiliar
mapear_transicion_a_edad <- function(age_inicial) {
  return(age_inicial)
}

# Crear consumo_poblacional_global
consumo_poblacional_global <- consumo_filtrado %>%
  mutate(edad_consumo = mapear_transicion_a_edad(age_inicial)) %>%
  left_join(poblacion_global_filtrada, by = c("year_consumo" = "year", "edad_consumo" = "age")) %>%
  filter(!is.na(N_mar_total)) %>%
  mutate(
    consumo_total_poblacional_t = (consumption_kg * N_mar_total) / 1000,
    consumo_anchoveta_poblacional_t = (consumo_anchoveta_kg * N_mar_total) / 1000,
    consumo_sardina_poblacional_t = (consumo_sardina_kg * N_mar_total) / 1000,
    consumo_otros_poblacional_t = (consumo_otros_kg * N_mar_total) / 1000,
    consumo_se_poblacional_t = ifelse(!is.na(consumption_se_kg), 
                                      (consumption_se_kg * N_mar_total) / 1000, NA),
    consumo_lower_poblacional_t = ifelse(!is.na(consumo_min_kg), 
                                         (consumo_min_kg * N_mar_total) / 1000, NA),
    consumo_upper_poblacional_t = ifelse(!is.na(consumo_max_kg), 
                                         (consumo_max_kg * N_mar_total) / 1000, NA),
    consumo_cv_poblacional = ifelse(!is.na(consumption_se_kg) & consumption_kg > 0,
                                    consumption_se_kg / consumption_kg, NA)
  )

# Crear resumen_anual_global
resumen_anual_global <- consumo_poblacional_global %>%
  group_by(year_consumo) %>%
  summarise(
    poblacion_total = sum(N_mar_total, na.rm = TRUE),
    consumo_total_t = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_anchoveta_t = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_t = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_t = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    consumo_se_total_t = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    consumo_lower_total_t = sum(consumo_lower_poblacional_t, na.rm = TRUE),
    consumo_upper_total_t = sum(consumo_upper_poblacional_t, na.rm = TRUE),
    consumo_cv_total = ifelse(consumo_total_t > 0, consumo_se_total_t / consumo_total_t, NA),
    consumo_ci95_lower_t = consumo_total_t - (1.96 * consumo_se_total_t),
    consumo_ci95_upper_t = consumo_total_t + (1.96 * consumo_se_total_t),
    .groups = 'drop'
  )

# Crear consumo_poblacional_rio y resumen_anual_rio
consumo_poblacional_rio <- consumo_filtrado %>%
  mutate(edad_consumo = mapear_transicion_a_edad(age_inicial)) %>%
  left_join(poblacion_rio_filtrada, by = c("year_consumo" = "year", "edad_consumo" = "age")) %>%
  filter(!is.na(N_mar)) %>%
  mutate(
    consumo_total_poblacional_t = (consumption_kg * N_mar) / 1000,
    consumo_se_poblacional_t = ifelse(!is.na(consumption_se_kg), 
                                      (consumption_se_kg * N_mar) / 1000, NA)
  )

resumen_anual_rio <- consumo_poblacional_rio %>%
  group_by(year_consumo, Rio) %>%
  summarise(
    consumo_total_t = sum(consumo_total_poblacional_t, na.rm = TRUE),
    consumo_se_total_t = sqrt(sum(consumo_se_poblacional_t^2, na.rm = TRUE)),
    .groups = 'drop'
  )

print("\n=== GENERATING PLOTS ===")

# Grayscale colors
colores_edad <- c("1" = "grey85", "2" = "grey65", "3" = "grey45", "4" = "grey25")
colores_edad_transicion <- c("1 → 2" = "grey85", "2 → 3" = "grey65", "3 → 4" = "grey45", "4 → 5" = "grey25")

colores_prey <- c("Anchovy" = "grey80", "Sardine" = "grey50", "Others" = "grey20")
colores_river <- c("Bueno" = "grey85", "Tolten" = "grey65", "Imperial" = "grey45", "Valdivia" = "grey25")

# Plot 1: Temporal evolution of total global consumption WITH UNCERTAINTY
p1 <- ggplot(resumen_anual_global, aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "grey60", alpha = 0.9, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, consumo_ci95_lower_t),
                    ymax = consumo_ci95_upper_t),
                width = 0.25, color = "black", linewidth = 0.9) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Total population consumption by year",
    x = NULL,
    y = "Total consumption (tonnes year⁻¹)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(face = "bold")
  ) +
  geom_text(aes(label = paste0(scales::comma(round(consumo_total_t, 0)),
                               "\n(CV:", round(consumo_cv_total*100, 1), "%)")),
            vjust = -1.5, size = 4)


print(p1)
ggsave("outputs/graficos/population_consumption_total_ISOTOPES.png", plot = p1, width = 14, height = 10, dpi = 300)

# Plot 2: Consumption by age and year (global) WITH UNCERTAINTY
p2 <- consumo_poblacional_global %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_poblacional_t, fill = factor(transicion_label))) +
  geom_col(alpha = 0.9, color = "black", size = 0.2) +
  scale_fill_manual(values = colores_edad_transicion, name = "Age") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Population Consumption by Age and Year - ISOTOPES",
    subtitle = "Age breakdown of total consumption - Global level",
    x = "",
    y = "Consumption (tonnes/year)\n",
    caption = "Diet ages 1-2: muscle isotopes 7F"
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
ggsave("outputs/graficos/population_consumption_by_age_ISOTOPES.png", plot = p2, width = 14, height = 10, dpi = 300)

# Plot 3: Consumption by prey type (annual global)

consumo_prey_anual <- resumen_anual_global %>%
  select(year_consumo, consumo_anchoveta_t, consumo_sardina_t, consumo_otros_t) %>%
  pivot_longer(
    cols = starts_with("consumo_"),
    names_to = "prey",
    values_to = "consumo_t"
  ) %>%
  mutate(
    prey = case_when(
      prey == "consumo_anchoveta_t" ~ "Anchovy",
      prey == "consumo_sardina_t"   ~ "Sardine",
      prey == "consumo_otros_t"     ~ "Others"
    ),
    # Force a consistent order for stacking + legend
    prey = factor(prey, levels = c("Anchovy", "Sardine", "Others"))
  )

# B/W palette with strong contrast (print-friendly)
colores_prey_bw <- c(
  "Anchovy" = "grey20",
  "Sardine" = "grey55",
  "Others"  = "grey85"
)

p3 <- ggplot(consumo_prey_anual, aes(x = factor(year_consumo), y = consumo_t, fill = prey)) +
  geom_col(alpha = 1, color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = colores_prey_bw,
    name = "Prey type",
    breaks = c("Anchovy", "Sardine", "Others") # legend order = visual order
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Population consumption by prey type",
    x = NULL,
    y = "Consumption (tonnes year\u207B\u00B9)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

print(p3)

ggsave(
  "outputs/graficos/population_consumption_by_prey.png",
  plot = p3, width = 14, height = 10, dpi = 300
)

# Plot 4: Consumption by river WITH UNCERTAINTY
p4 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t, fill = Rio)) +
  geom_col(alpha = 0.9, color = "black", size = 0.2) +
  scale_fill_manual(values = colores_river, name = "River") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Population Consumption by River and Year - ISOTOPES",
    subtitle = "Contribution of each river to total consumption",
    x = "",
    y = "Consumption (tonnes/year)\n",
    caption = "Diet ages 1-2: muscle isotopes 7F"
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
ggsave("outputs/graficos/population_consumption_by_river_ISOTOPES.png", plot = p4, width = 14, height = 10, dpi = 300)

# Plot 5: Faceted by river - temporal evolution
p5 <- resumen_anual_rio %>%
  ggplot(aes(x = factor(year_consumo), y = consumo_total_t)) +
  geom_col(fill = "grey50", alpha = 0.8, color = "black", size = 0.3) +
  facet_wrap(~ Rio, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Consumption Evolution by River - ISOTOPES",
    subtitle = "Independent temporal trends by river",
    x = "",
    y = "Consumption (tonnes/year)\n",
    caption = "Diet ages 1-2: muscle isotopes 7F"
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
ggsave("outputs/graficos/population_consumption_faceted_ISOTOPES.png", plot = p5, width = 14, height = 10, dpi = 300)





# ================================================================
# ESTADÍSTICAS PARA RESULTADOS - CONSUMO INDIVIDUAL
# ================================================================

# Cargar datos actualizados
resultados_modelos <- readRDS("data/data_raw/bioenergetic-model/modelos_total_by_age_ISOTOPES.rds")
consumo_individual <- generar_tabla_consumo_individual(resultados_modelos)

# Agregar year_consumo
consumo_individual <- consumo_individual %>%
  mutate(year_consumo = year - 1)

cat("=== P-VALUES ===\n")
cat(sprintf("Rango global: %.4f – %.4f\n", 
    min(consumo_individual$p_value_est), 
    max(consumo_individual$p_value_est)))

consumo_individual %>%
  group_by(age_inicial) %>%
  summarise(p_mean = round(mean(p_value_est), 3)) %>%
  arrange(age_inicial) %>%
  print()

cat("\n=== CONSUMO INDIVIDUAL POR EDAD ===\n")
consumo_individual %>%
  group_by(age_inicial) %>%
  summarise(
    mean_kg  = round(mean(consumption_kg), 1),
    min_kg   = round(min(consumption_kg, na.rm = TRUE), 1),
    max_kg   = round(max(consumption_kg, na.rm = TRUE), 1),
    cv_pct   = round(sd(consumption_kg) / mean(consumption_kg) * 100, 1),
    n        = n()
  ) %>%
  arrange(age_inicial) %>%
  print()

cat("\n=== RATIO CONSUMO/PESO INICIAL (age-1) ===\n")
peso_inicial_age1 <- 700  # g
consumo_age1 <- consumo_individual %>% 
  filter(age_inicial == 1) %>% 
  pull(consumption_kg) %>% mean()
cat(sprintf("%.1f kg / %.3f kg = %.1f-times initial body mass\n", 
    consumo_age1, peso_inicial_age1/1000, consumo_age1/(peso_inicial_age1/1000)))

cat("\n=== FOURFOLD INCREASE ===\n")
medias <- consumo_individual %>%
  group_by(age_inicial) %>%
  summarise(mean_kg = mean(consumption_kg)) %>%
  arrange(age_inicial)
cat(sprintf("Age-1: %.1f kg | Age-4: %.1f kg | Ratio: %.1f-fold\n",
    medias$mean_kg[1], medias$mean_kg[4], medias$mean_kg[4]/medias$mean_kg[1]))













cat("=== CONSUMO POBLACIONAL - ESTADÍSTICAS PARA RESULTS ===\n")

# Rango anual
resumen_anual_global %>%
  select(year_consumo, consumo_total_t, consumo_cv_total) %>%
  mutate(
    consumo_total_t = round(consumo_total_t, 0),
    cv_pct = round(consumo_cv_total * 100, 1)
  ) %>%
  arrange(year_consumo) %>%
  print()

cat("\n=== RANGO Y FOLD ===\n")
cat(sprintf("Mínimo: %s t (%d)\n", 
    scales::comma(round(min(resumen_anual_global$consumo_total_t))),
    resumen_anual_global$year_consumo[which.min(resumen_anual_global$consumo_total_t)]))
cat(sprintf("Máximo: %s t (%d)\n",
    scales::comma(round(max(resumen_anual_global$consumo_total_t))),
    resumen_anual_global$year_consumo[which.max(resumen_anual_global$consumo_total_t)]))
cat(sprintf("Fold: %.0f-fold\n",
    max(resumen_anual_global$consumo_total_t)/min(resumen_anual_global$consumo_total_t)))
cat(sprintf("CV range: %.1f%% – %.1f%%\n",
    min(resumen_anual_global$consumo_cv_total*100, na.rm=TRUE),
    max(resumen_anual_global$consumo_cv_total*100, na.rm=TRUE)))

cat("\n=== PORCENTAJES DE PRESA ACUMULADOS ===\n")
cat(sprintf("Sardine: %.1f%%\n", sum(resumen_anual_global$consumo_sardina_t)/sum(resumen_anual_global$consumo_total_t)*100))
cat(sprintf("Anchovy: %.1f%%\n", sum(resumen_anual_global$consumo_anchoveta_t)/sum(resumen_anual_global$consumo_total_t)*100))
cat(sprintf("Others:  %.1f%%\n", sum(resumen_anual_global$consumo_otros_t)/sum(resumen_anual_global$consumo_total_t)*100))









# ================================================================
# TABLE S4.2: Weight inputs for bioenergetic simulations
# ================================================================
library(tidyverse)
library(flextable)
library(officer)
library(janitor)

# Cargar datos
chinook_biologico <- read.csv("data/data_raw/biological-data/data_chinook_cleaned.csv") |>
  clean_names()

densidad_energetica <- read.csv("data/data_raw/energy-density/energy-density-weight-by-age.csv") |>
  clean_names()

# Procesar temporadas
chinook_temporadas <- chinook_biologico |>
  filter(!is.na(age), !is.na(year), !is.na(tw_g), !is.na(month),
         species == "Oncorhynchus tshawytscha",
         month %in% c(10, 11, 12, 1, 2, 3, 4)) |>
  mutate(
    year_temporada = case_when(
      month %in% c(10, 11, 12) ~ year + 1,
      month %in% c(1, 2, 3, 4) ~ year
    )
  )

umbral_datos <- 5
temporadas  <- 2015:2023
transiciones <- list(c(1,2), c(2,3), c(3,4), c(4,5))

tabla_s42 <- purrr::map_dfr(temporadas, function(temp) {
  purrr::map_dfr(transiciones, function(trans) {
    age_ini <- trans[1]
    age_fin <- trans[2]
    year_consumo <- temp - 1

    # Pesos finales
    pesos_fin <- chinook_temporadas |>
      filter(year_temporada == temp, age == age_fin) |> pull(tw_g)

    # Pesos iniciales
    pesos_ini <- chinook_temporadas |>
      filter(year_temporada == temp, age == age_ini) |> pull(tw_g)

    # Peso final y fuente
    if (length(pesos_fin) >= umbral_datos) {
      w_final      <- round(mean(pesos_fin), 0)
      n_final      <- length(pesos_fin)
      source_final <- "Year-specific"
    } else {
      pesos_hist   <- chinook_temporadas |> filter(age == age_fin) |> pull(tw_g)
      w_final      <- round(mean(pesos_hist), 0)
      n_final      <- length(pesos_fin)
      source_final <- "Historical mean"
    }

    # Peso inicial y fuente
    if (age_ini == 1) {
      params       <- densidad_energetica[densidad_energetica$edad == 1, ]
      w_initial    <- round(params$peso_inicial, 0)
      n_initial    <- NA_integer_
      source_ini   <- "Fixed (2025 sampling)"
    } else if (length(pesos_ini) >= umbral_datos) {
      w_initial    <- round(mean(pesos_ini), 0)
      n_initial    <- length(pesos_ini)
      source_ini   <- "Year-specific"
    } else if (length(pesos_ini) > 0) {
      w_initial    <- round(mean(pesos_ini), 0)
      n_initial    <- length(pesos_ini)
      source_ini   <- "Year-specific (n < 5)"
    } else {
      pesos_hist_i <- chinook_temporadas |> filter(age == age_ini) |> pull(tw_g)
      w_initial    <- round(mean(pesos_hist_i), 0)
      n_initial    <- 0L
      source_ini   <- "Historical mean"
    }

    data.frame(
      `Calendar year`    = year_consumo,
      `Age transition`   = paste0(age_ini, " → ", age_fin),
      `W_initial (g)`    = w_initial,
      `Source (initial)` = source_ini,
      `n (initial)`      = ifelse(is.na(n_initial), "—", as.character(n_initial)),
      `W_final (g)`      = w_final,
      `Source (final)`   = source_final,
      `n (final)`        = as.character(n_final),
      check.names = FALSE
    )
  })
})




# Crear flextable
ft <- flextable(tabla_s42) |>
  bold(part = "header") |>
  align(align = "center", part = "all") |>
  align(j = c("Age transition", "Source (initial)", "Source (final)"),
        align = "left", part = "all") |>
  bg(i = ~ `Source (initial)` == "Historical mean" | 
          `Source (final)`   == "Historical mean",
     bg = "#FFF3CD") |>  # resaltar filas con datos históricos
  autofit() |>
  set_caption("Table S4.2. Initial and final body mass inputs for each age–year bioenergetic simulation. Year-specific values are means from samples with n ≥ 5 individuals collected during the return season. Historical mean values were used when year-specific sample sizes were insufficient (n < 5). Age-1 initial weight was fixed at 700 g based on 2025 field sampling, applied uniformly across all simulation years. Rows highlighted in yellow indicate use of historical mean weights.")

# Exportar
doc <- read_docx() |> body_add_flextable(ft)
print(doc, target = "outputs/tablas/Table_S4_2_weight_inputs.docx")
cat("✅ Table S4.2 guardada\n")

# Preview en consola
cat("\n=== PREVIEW TABLE S4.2 ===\n")
print(tabla_s42)











# Tablas -----------------------------------------------------------------

library(tidyverse)
library(knitr)
library(kableExtra)

# ==============================================================
# TABLA 1: Consumo anual total con incertidumbre
# ==============================================================

tabla1 <- resumen_anual_global %>%
  mutate(
    Year = year_consumo,
    `Total Consumption (t)` = round(consumo_total_t, 0),
    `95% CI Lower (t)` = round(pmax(0, consumo_ci95_lower_t), 0),
    `95% CI Upper (t)` = round(consumo_ci95_upper_t, 0),
    `CV (%)` = round(consumo_cv_total * 100, 1)
  ) %>%
  select(Year, `Total Consumption (t)`, `95% CI Lower (t)`, `95% CI Upper (t)`, `CV (%)`)

print("=== TABLA 1: Annual Total Consumption ===")
print(tabla1)

# Guardar tabla 1
write.csv(tabla1, "outputs/tablas/Table1_Annual_Total_Consumption.csv", row.names = FALSE)

# Versión formateada para Word/LaTeX
kable(tabla1, format = "markdown", align = "c")

# ==============================================================
# TABLA 2: Composición por presa (% por año)
# ==============================================================

tabla2_year <- resumen_anual_global %>%
  mutate(
    Year = year_consumo,
    `Anchovy (%)` = round((consumo_anchoveta_t / consumo_total_t) * 100, 1),
    `Sardine (%)` = round((consumo_sardina_t / consumo_total_t) * 100, 1),
    `Others (%)` = round((consumo_otros_t / consumo_total_t) * 100, 1)
  ) %>%
  select(Year, `Anchovy (%)`, `Sardine (%)`, `Others (%)`)

print("=== TABLA 2A: Prey Composition by Year ===")
print(tabla2_year)

# Guardar tabla 2A
write.csv(tabla2_year, "outputs/tablas/Table2A_Prey_Composition_by_Year.csv", row.names = FALSE)

# Versión formateada
kable(tabla2_year, format = "markdown", align = "c")

# ==============================================================
# TABLA 2 ALTERNATIVA: Composición por presa (% por edad)
# ==============================================================

# Primero necesitamos calcular el resumen por edad
tabla2_age <- consumo_poblacional_global %>%
  group_by(transicion_label) %>%
  summarise(
    consumo_anchoveta_total = sum(consumo_anchoveta_poblacional_t, na.rm = TRUE),
    consumo_sardina_total = sum(consumo_sardina_poblacional_t, na.rm = TRUE),
    consumo_otros_total = sum(consumo_otros_poblacional_t, na.rm = TRUE),
    consumo_total = sum(consumo_total_poblacional_t, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    `Age Transition` = transicion_label,
    `Anchovy (%)` = round((consumo_anchoveta_total / consumo_total) * 100, 1),
    `Sardine (%)` = round((consumo_sardina_total / consumo_total) * 100, 1),
    `Others (%)` = round((consumo_otros_total / consumo_total) * 100, 1)
  ) %>%
  select(`Age Transition`, `Anchovy (%)`, `Sardine (%)`, `Others (%)`)

print("=== TABLA 2B: Prey Composition by Age Transition ===")
print(tabla2_age)

# Guardar tabla 2B
write.csv(tabla2_age, "outputs/tablas/Table2B_Prey_Composition_by_Age.csv", row.names = FALSE)

# Versión formateada
kable(tabla2_age, format = "markdown", align = "c")





library(tidyverse)

# Cargar datos
raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length","weight",
                     "Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

# Asignar grupos
raw_data$group <- case_when(
  str_detect(raw_data$Sample, "^m-")   ~ "Salmon_muscle",
  str_detect(raw_data$Sample, "^h-")   ~ "Salmon_liver",
  str_detect(raw_data$Sample, "^Sar")  ~ "Sardine",
  str_detect(raw_data$Sample, "^Zoo")  ~ "Zooplankton",
  str_detect(raw_data$Sample, "Fito")  ~ "Phytoplankton",
  str_detect(raw_data$Sample, "^Anch") ~ "Anchovy",
  TRUE ~ "Unknown"
)

# Corrección lipídica Chinook (sin umbral, consistente con Lerner et al., 2021)
lipid_correction_optimized_chinook <- function(d13C, CN_ratio, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN_ratio - 0.775, 0.001)
  L_percent  <- 93 / (1 + (1/lipid_term))
  correction_factor <- D * (I + (3.90 / (1 + (287 / L_percent))))
  d13C + correction_factor
}

raw_data <- raw_data %>%
  mutate(d13C_corrected = case_when(
    group == "Salmon_muscle" ~ lipid_correction_optimized_chinook(d13C, CtoN),
    TRUE ~ d13C
  ))

# Excluir Sar-2-2 (suma elemental > 100%, identificada en QC)
raw_data <- raw_data %>% filter(Sample != "Sar-2 - 2")

# --- TABLA ISOTÓPICA ---

# Parte 1: Salmón por edad (músculo)
salmon_table <- raw_data %>%
  filter(group == "Salmon_muscle", !is.na(age), age %in% 1:3) %>%
  group_by(age) %>%
  summarise(
    Group     = paste0("Chinook salmon (age ", first(age), ")"),
    n         = n(),
    `δ13C (‰)` = paste0(round(mean(d13C_corrected, na.rm=TRUE), 2), " ± ",
                         round(sd(d13C_corrected,   na.rm=TRUE), 2)),
    `δ15N (‰)` = paste0(round(mean(d15N, na.rm=TRUE), 2), " ± ",
                         round(sd(d15N,   na.rm=TRUE), 2)),
    `δ34S (‰)` = paste0(round(mean(d34S, na.rm=TRUE), 2), " ± ",
                         round(sd(d34S,   na.rm=TRUE), 2)),
    .groups = "drop"
  ) %>%
  select(Group, n, `δ13C (‰)`, `δ15N (‰)`, `δ34S (‰)`)

# Parte 2: Fuentes de presa (marinas)
sources_table <- raw_data %>%
  filter(group %in% c("Sardine","Anchovy","Zooplankton","Phytoplankton")) %>%
  mutate(Group = group) %>%
  group_by(Group) %>%
  summarise(
    n         = n(),
    `δ13C (‰)` = paste0(round(mean(d13C_corrected, na.rm=TRUE), 2), " ± ",
                         round(sd(d13C_corrected,   na.rm=TRUE), 2)),
    `δ15N (‰)` = paste0(round(mean(d15N, na.rm=TRUE), 2), " ± ",
                         round(sd(d15N,   na.rm=TRUE), 2)),
    `δ34S (‰)` = paste0(round(mean(d34S, na.rm=TRUE), 2), " ± ",
                         round(sd(d34S,   na.rm=TRUE), 2)),
    .groups = "drop"
  )

# Combinar y guardar
isotope_table <- bind_rows(salmon_table, sources_table)

write.csv(isotope_table, "outputs/tablas/Table3_Isotopic_Composition.csv", row.names = FALSE)

print(isotope_table)





# ---- TABLE 2: Energy Density by Age and Species ----

library(tidyverse)

# Cargar datos
de_salmones    <- read.csv("data/data_raw/energy-density/DE_salmones.csv")
de_sardinas    <- read.csv("data/data_raw/energy-density/sardina_densidad_energetica.csv")
factor_wet_dry <- read.csv("data/data_raw/energy-density/wet_dry_weight.csv")
se_delta       <- read.csv("data/data_raw/energy-density/se_delta_pred_45.csv")

# Factor húmedo-seco
fc <- factor_wet_dry %>%
  mutate(factor = (peso_humedo - peso_seco) / peso_humedo) %>%
  summarise(factor = mean(factor)) %>%
  pull(factor)

# Salmones edades 1-3: empírico
de_salmon_tabla <- de_salmones %>%
  group_by(Individuo, AGE, peso) %>%
  summarise(DE_ind = mean(DE) * (1 - fc), .groups = "drop") %>%
  group_by(AGE) %>%
  summarise(
    Species       = "Chinook salmon",
    n             = n(),
    `Mean weight (g)` = round(mean(peso), 0),
    `Mean ED (J/g)`   = round(mean(DE_ind), 0),
    `SD`              = round(sd(DE_ind), 0),
    Source        = "Empirical",
    .groups       = "drop"
  ) %>%
  rename(Age = AGE)

# Salmones edades 4-5: predicho con delta method
de_salmon_pred <- se_delta %>%
  transmute(
    Age               = edad,
    Species           = "Chinook salmon",
    n                 = NA_integer_,
    `Mean weight (g)` = round(peso, 0),
    `Mean ED (J/g)`   = ED_mean,
    `SD`              = ED_sd,
    Source            = "Predicted"
  )

# Sardinas: empírico
de_sardina_tabla <- de_sardinas %>%
  filter(!is.na(DE)) |> 
  group_by(Individuo) %>%
  summarise(DE_ind = mean(DE) * 0.2, .groups = "drop") %>%
  summarise(
    Age               = NA_integer_,
    Species           = "Sardine",
    n                 = n(),
    `Mean weight (g)` = NA_integer_,
    `Mean ED (J/g)`   = round(mean(DE_ind), 0),
    `SD`              = round(sd(DE_ind), 0),
    Source            = "Empirical"
  )

# Combinar
tabla_DE_final <- bind_rows(de_salmon_tabla, de_salmon_pred, de_sardina_tabla) %>%
  dplyr::select(Species, Age, n, `Mean weight (g)`, `Mean ED (J/g)`, SD, Source)

write.csv(tabla_DE_final, "outputs/tablas/Table2_Energy_Density.csv", row.names = FALSE)
print(tabla_DE_final)














# =============================================================================
# TABLA S1: Composición isotópica por grupo
# =============================================================================

library(tidyverse)
library(flextable)
library(officer)

# Cargar datos
raw_data <- read.csv("data/data_raw/diet-isotopes/stable-isotopes.csv", 
                     fileEncoding = "latin1")
names(raw_data) <- c("ID","Sample","date","age","fork_length","total_length",
                     "weight","Mass_mg","d15N","d13C","d34S","pctN","pctC","pctS","CtoN")

# Grupos
raw_data <- raw_data %>%
  mutate(group = case_when(
    str_detect(Sample, "^m-") ~ "Salmon_muscle",
    str_detect(Sample, "^h-") ~ "Salmon_liver",
    str_detect(Sample, "^Sar") ~ "Sardine",
    str_detect(Sample, "^Zoo") ~ "Zooplankton",
    str_detect(Sample, "Fito") ~ "Phytoplankton",
    str_detect(Sample, "^Anch") ~ "Anchovy",
    TRUE ~ "Unknown"
  ))

# Excluir Sar-2-2
raw_data <- raw_data %>% filter(Sample != "Sar-2 - 2")

# Corrección lipídica (McConnaughey optimizado Chinook)
lipid_correction <- function(d13C, CN, D = 6.31, I = 0.103) {
  lipid_term <- pmax(0.246 * CN - 0.775, 0.001)
  L <- 93 / (1 + (1 / lipid_term))
  d13C + D * (I + 3.90 / (1 + 287 / L))
}

raw_data <- raw_data %>%
  mutate(
    d13C_corrected = case_when(
      group %in% c("Salmon_muscle", "Salmon_liver") ~ lipid_correction(d13C, CtoN),
      TRUE ~ d13C
    )
  )

# Calcular estadísticas
table_s1 <- raw_data %>%
  filter(group %in% c("Salmon_muscle", "Salmon_liver",
                      "Sardine", "Anchovy", "Zooplankton", "Phytoplankton")) %>%
  mutate(
    group_label = case_when(
      group == "Salmon_muscle" ~ paste0("Salmon muscle age-", age),
      group == "Salmon_liver"  ~ paste0("Salmon liver age-",  age),
      TRUE ~ group
    )
  ) %>%
  group_by(group_label) %>%
  summarise(
    n        = n(),
    d13C     = sprintf("%.2f \u00b1 %.2f", mean(d13C_corrected, na.rm=T), sd(d13C_corrected, na.rm=T)),
    d15N     = sprintf("%.2f \u00b1 %.2f", mean(d15N,           na.rm=T), sd(d15N,           na.rm=T)),
    n_d34S   = sum(!is.na(d34S)),
    d34S     = ifelse(n_d34S > 0,
                      sprintf("%.2f \u00b1 %.2f", mean(d34S, na.rm=T), sd(d34S, na.rm=T)),
                      "\u2014"),
    .groups  = "drop"
  ) %>%
  mutate(
    notes = case_when(
      group_label == "Phytoplankton" ~ paste0("n=", n_d34S, " for \u03b434S"),
      TRUE ~ ""
    )
  ) %>%
  select(group_label, n, d13C, d15N, d34S, notes)

# Ordenar filas
order_groups <- c(
  "Salmon muscle age-1", "Salmon muscle age-2", "Salmon muscle age-3",
  "Salmon liver age-1",  "Salmon liver age-2",  "Salmon liver age-3",
  "Sardine", "Anchovy", "Zooplankton", "Phytoplankton"
)

table_s1 <- table_s1 %>%
  mutate(group_label = factor(group_label, levels = order_groups)) %>%
  arrange(group_label)

# Renombrar columnas
names(table_s1) <- c("Group", "n", "\u03b4\u00b9\u00b3C (\u2030)", 
                      "\u03b4\u00b9\u2075N (\u2030)", "\u03b4\u00b3\u2074S (\u2030)", "Notes")

print(table_s1)

# Exportar como Word con flextable
ft <- flextable(table_s1) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  autofit() %>%
  set_caption("Table S1. Isotopic composition (mean \u00b1 SD) of Chinook salmon tissues and dietary source groups. \u03b4\u00b9\u00b3C values for salmon are lipid-corrected. Sardine n=14 after exclusion of Sar-2-2. One phytoplankton sample (Fito 1-1) lacked \u03b4\u00b3\u2074S.")

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "output/Table_S1_isotopes.docx")
cat("\u2705 Guardado: output/Table_S1_isotopes.docx\n")


# Cargar fuentes de agua dulce
freshwater <- read.csv("data/data_raw/diet-isotopes/datos_rios.csv") %>%
  rename(group_label = Source) %>%
  mutate(
    n     = NA_integer_,
    d13C  = sprintf("%.2f \u00b1 %.2f", d13C_mean, d13C_sd),
    d15N  = sprintf("%.2f \u00b1 %.2f", d15N_mean, d15N_sd),
    d34S  = sprintf("%.2f \u00b1 %.2f", d34S_mean, d34S_sd),
    notes = "Literature"
  ) %>%
  select(group_label, n, d13C, d15N, d34S, notes)

names(freshwater) <- names(table_s1)

# Combinar con tabla principal
table_s1_complete <- bind_rows(table_s1, freshwater)

print(table_s1_complete)






# ================================================================
# TABLA: CV de p-values y consumo por clase de edad
# Para verificar el CV = 27% reportado en Discussion
# ================================================================

tabla_cv_pvalues <- consumo_individual |>
  group_by(age_inicial) |>
  summarise(
    n                = n(),
    # P-values
    p_mean           = round(mean(p_value_est, na.rm = TRUE), 3),
    p_sd             = round(sd(p_value_est, na.rm = TRUE), 3),
    p_min            = round(min(p_value_est, na.rm = TRUE), 3),
    p_max            = round(max(p_value_est, na.rm = TRUE), 3),
    p_cv_pct         = round(sd(p_value_est, na.rm = TRUE) / 
                             mean(p_value_est, na.rm = TRUE) * 100, 1),
    # Consumo (kg)
    consumo_mean_kg  = round(mean(consumption_kg, na.rm = TRUE), 1),
    consumo_sd_kg    = round(sd(consumption_kg, na.rm = TRUE), 1),
    consumo_cv_pct   = round(sd(consumption_kg, na.rm = TRUE) / 
                             mean(consumption_kg, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) |>
  rename(
    "Age class"      = age_inicial,
    "N"              = n,
    "P mean"         = p_mean,
    "P SD"           = p_sd,
    "P min"          = p_min,
    "P max"          = p_max,
    "P CV (%)"       = p_cv_pct,
    "Consumption mean (kg)" = consumo_mean_kg,
    "Consumption SD (kg)"   = consumo_sd_kg,
    "Consumption CV (%)"    = consumo_cv_pct
  )

print(tabla_cv_pvalues)

# Exportar como CSV
write.csv(tabla_cv_pvalues, 
          "outputs/tabla_cv_pvalues_consumo_por_edad.csv", 
          row.names = FALSE)

















# === TABLE S3: Individual annual consumption by age ===
library(tidyverse)
library(flextable)
library(officer)

# Cargar resultados del modelo bioenergetico
resultados_modelos <- readRDS("data/data_raw/bioenergetic-model/modelos_total_by_age_ISOTOPES.rds")
consumo_individual <- generar_tabla_consumo_individual(resultados_modelos, verbose = FALSE)

# Formatear para publicación
table_s42 <- consumo_individual %>%
  group_by(transicion_label, age_inicial) %>%
  summarise(
    `Mean (kg/year)`   = round(mean(consumption_kg, na.rm = TRUE), 1),
    `SD`               = round(sd(consumption_kg, na.rm = TRUE), 1),
    `Min (kg/year)`    = round(min(consumption_kg, na.rm = TRUE), 1),
    `Max (kg/year)`    = round(max(consumption_kg, na.rm = TRUE), 1),
    anchovy_kg         = round(mean(consumo_anchoveta_kg, na.rm = TRUE), 1),
    sardine_kg         = round(mean(consumo_sardina_kg, na.rm = TRUE), 1),
    others_kg          = round(mean(consumo_otros_kg, na.rm = TRUE), 1),
    `n years`          = n(),
    .groups = "drop"
  ) %>%
  mutate(
    total_kg           = anchovy_kg + sardine_kg + others_kg,
    `Anchovy (kg/year)` = anchovy_kg,
    `Anchovy (%)`       = round(anchovy_kg / total_kg * 100, 1),
    `Sardine (kg/year)` = sardine_kg,
    `Sardine (%)`       = round(sardine_kg / total_kg * 100, 1),
    `Others (kg/year)`  = others_kg,
    `Others (%)`        = round(others_kg / total_kg * 100, 1),
    `Age transition`    = transicion_label
  ) %>%
  arrange(age_inicial) %>%
  select(`Age transition`, `Mean (kg/year)`, SD, `Min (kg/year)`, `Max (kg/year)`,
         `Anchovy (kg/year)`, `Anchovy (%)`,
         `Sardine (kg/year)`, `Sardine (%)`,
         `Others (kg/year)`, `Others (%)`,
         `n years`)

table_s42 |> as.data.frame()

table_s42 |> write.csv("outputs/tablas/Table_S42_individual_consumption.csv", row.names = FALSE)

# Exportar Word
ft_s42 <- flextable(table_s42) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  autofit() %>%
  set_caption("Table S3. Individual annual consumption estimates (kg year⁻¹) per age transition, averaged across simulation years (2014–2022). Values represent consumption during the transition year (e.g., age 1→2 = fish in their first year at sea).")

doc_s42 <- read_docx() %>% body_add_flextable(ft_s42)
print(doc_s42, target = "outputs/tablas/Table_S42_individual_consumption.docx")
cat("✅ Table S3 guardada\n")


















# === TABLE S5.2: Annual age structure of returning spawners ===
library(tidyverse)
library(flextable)
library(officer)

# Cargar datos (objeto generado en script 03)
estructura_edades <- read.csv("data/data_raw/biological-data/age_composition_by_year.csv")

table_s52 <- estructura_edades %>%
  transmute(
    `Season`   = year,
    `n total`  = n_total,
    `Age-2 (%)` = proporcion_pct_2,
    `n age-2`  = n_individuos_2,
    `Age-3 (%)` = proporcion_pct_3,
    `n age-3`  = n_individuos_3,
    `Age-4 (%)` = proporcion_pct_4,
    `n age-4`  = n_individuos_4,
    `Age-5 (%)` = proporcion_pct_5,
    `n age-5`  = n_individuos_5
  )

ft_s52 <- flextable(table_s52) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit() %>%
  set_caption("Table S5.2. Annual age structure of returning Chinook salmon spawners. Proportions (%) and sample sizes (n) by age class and return season. Seasons without biological sampling received mean age structure across available seasons (see Section S5.3).")

table_s52 |> write.csv("outputs/tablas/Table_S5_2_age_structure.csv", row.names = FALSE)


doc_s52 <- read_docx() %>% body_add_flextable(ft_s52)
print(doc_s52, target = "outputs/tablas/Table_S5_2_age_structure.docx")
cat("✅ Table S5.2 guardada\n")






# === TABLE S5.3: Annual returning spawner abundance by river ===
library(tidyverse)
library(flextable)
library(officer)

# Cargar datos generados por script 02
poblacion <- read.csv("data/data_raw/river-population/total_population_from_nb.csv")

# Pivotar a formato ancho (una columna por río)
table_s53 <- poblacion %>%
  mutate(
    N_estimado = round(N_estimado, 0),
    # Identificar método por río para footnote
    metodo_label = case_when(
      metodo == "Datos_Reales"      ~ "a",
      metodo == "N_base_tendencia"  ~ "b",
      metodo == "Nb_factor_promedio" ~ "c"
    )
  ) %>%
  # Crear columna con N + superíndice método
  mutate(N_label = paste0(scales::comma(N_estimado))) %>%
  select(year, Rio, N_label) %>%
  pivot_wider(names_from = Rio, values_from = N_label) %>%
  arrange(year) %>%
  # Convertir year a temporada
  mutate(Season = paste0(year - 1, "-", year)) %>%
  select(Season, Tolten, Imperial, Valdivia, Bueno)

# Renombrar columnas
names(table_s53) <- c("Season", "Toltén (N)", "Imperial (N)", 
                       "Valdivia (N)", "Bueno (N)")

table_s53 |> write.csv("outputs/tablas/Table_S5_3_spawner_abundance.csv", row.names = FALSE)


# Crear flextable
ft_s53 <- flextable(table_s53) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  autofit() %>%
  add_footer_lines("ᵃ Direct census (hydroacoustic monitoring; Espinoza Henríquez, 2023).") %>%
  add_footer_lines("ᵇ Scaled from single reference census (N = 21,724 in 2020; Gómez-Uchida et al., 2021).") %>%
  add_footer_lines("ᶜ Estimated from Nb using mean conversion factor f̄ = 0.001007, scaled by Toltén River trend.") %>%
  set_caption("Table S5.3. Annual returning spawner abundance (N) for each river system across the study period (2015–2023). Values represent total census size used as input for backward cohort reconstruction.")

# Exportar
doc_s53 <- read_docx() %>% body_add_flextable(ft_s53)
print(doc_s53, target = "outputs/tablas/Table_S5_3_spawner_abundance.docx")
cat("✅ Table S5.3 guardada\n")




















# =============================================================================
# VERIFICACIÓN PROPORCIONES DIETARIAS - MIXING MODEL
# Agregar a Graficos.R
# =============================================================================

library(cosimmr)
library(tidyverse)

# Cargar modelo seleccionado (5 fuentes, músculo)
model_muscle_5f <- readRDS("output/cosimmr_moel_8_muscle_5sources.rds")

# Función para extraer proporciones a un peso dado
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

# Extraer proporciones brutas (5 fuentes)
props_edad1 <- extraer_proporciones(model_muscle_5f, 700)
props_edad2 <- extraer_proporciones(model_muscle_5f, 3271)

cat("=== PROPORCIONES 5 FUENTES ===\n")
cat("\nEdad 1 (700g):\n"); print(round(props_edad1 * 100, 1))
cat("\nEdad 2 (3271g):\n"); print(round(props_edad2 * 100, 1))

# Agregar a 3 categorías bioenergéticas
agregar_3cat <- function(props) {
  c(
    Sardine   = as.numeric(props["Sardine"]),
    Anchovy   = as.numeric(props["Anchovy"]),
    Others    = as.numeric(props["Zooplankton"]) +
                as.numeric(props["Fish"]) +
                as.numeric(props["Invertebrates"])
  )
}

cat1 <- round(agregar_3cat(props_edad1) * 100, 1)
cat2 <- round(agregar_3cat(props_edad2) * 100, 1)

cat("\n=== PROPORCIONES 3 CATEGORÍAS BIOENERGÉTICAS ===\n")
cat("\nEdad 1:\n")
cat(sprintf("  Sardine: %.1f%%\n  Anchovy: %.1f%%\n  Others:  %.1f%%\n  Total:   %.1f%%\n",
    cat1["Sardine"], cat1["Anchovy"], cat1["Others"], sum(cat1)))
cat("\nEdad 2:\n")
cat(sprintf("  Sardine: %.1f%%\n  Anchovy: %.1f%%\n  Others:  %.1f%%\n  Total:   %.1f%%\n",
    cat2["Sardine"], cat2["Anchovy"], cat2["Others"], sum(cat2)))
