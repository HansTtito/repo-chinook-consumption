# ================================================================
# MODELOS BIOENERGÉTICOS - ANÁLISIS CON ISÓTOPOS
# Dieta edades 1-2: Modelo músculo 7 fuentes
# Dieta edades 3-5: Contenido estomacal
#
# CHECKPOINTS disponibles:
#   1. datos_base + datos_procesados  → "checkpoint_01_datos.rds"
#   2. resultados_modelos             → "checkpoint_02_modelos.rds"
#   3. consumo_individual             → "checkpoint_03_consumo.rds"
#   4. sensitivity_results            → "checkpoint_04_sensitivity.rds"
# ================================================================

library(fb4package)
library(dplyr)
library(janitor)
library(ggplot2)
library(stringr)
library(tidyr)
library(patchwork)
library(purrr)

# ── Helpers ────────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper para guardar/cargar checkpoints
checkpoint_save <- function(obj, path, verbose = TRUE) {
  saveRDS(obj, path)
  if (verbose) cat("✓ Checkpoint guardado:", path, "\n")
}

checkpoint_load <- function(path, verbose = TRUE) {
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (verbose) cat("✓ Checkpoint cargado:", path, "\n")
  obj
}

# Evalúa si recalcular: si force = TRUE o el checkpoint no existe, corre expr
# De lo contrario, carga desde disco
with_checkpoint <- function(path, expr, force = FALSE, verbose = TRUE) {
  if (!force && file.exists(path)) {
    return(checkpoint_load(path, verbose))
  }
  result <- expr
  checkpoint_save(result, path, verbose)
  result
}

# ================================================================
# 1. FUNCIONES DE CARGA Y VALIDACIÓN
# ================================================================

cargar_datos_base <- function(verbose = TRUE) {
  if (verbose) cat("=== CARGANDO DATOS BASE ===\n")

  data("fish4_parameters", package = "fb4package")
  chinook_params <- fish4_parameters[["Oncorhynchus tshawytscha"]]

  temperature_base  <- read.csv("data/data_raw/temperature/temperature_time_serie_daily.csv") |>
    dplyr::select(1, 2, 4)

  chinook_biologico <- read.csv("data/data_raw/biological-data/data_chinook_cleaned.csv") |>
    clean_names()

  densidad_energetica <- read.csv("data/data_raw/energy-density/energy-density-weight-by-age.csv") |>
    clean_names()

  chinook_diet_age <- read.csv("data/data_raw/biological-data/diet_proportion_by_age_ISOTOPES.csv") |>
    clean_names()

  if (verbose) {
    cat("Archivos cargados:\n")
    cat("- Parámetros especies: ✓\n")
    cat("- Temperatura:", nrow(temperature_base), "registros\n")
    cat("- Datos biológicos:", nrow(chinook_biologico), "registros\n")
    cat("- Densidad energética:", nrow(densidad_energetica), "edades\n")
    cat("- Dieta por edad (ISOTOPES):", nrow(chinook_diet_age), "registros\n\n")
  }

  list(
    chinook_params      = chinook_params,
    temperature_base    = temperature_base,
    chinook_biologico   = chinook_biologico,
    densidad_energetica = densidad_energetica,
    chinook_diet_age    = chinook_diet_age
  )
}

validar_datos <- function(datos_base, verbose = TRUE) {
  if (verbose) cat("=== VALIDANDO DATOS ===\n")

  errores <- c()
  if (nrow(datos_base$chinook_biologico) == 0)   errores <- c(errores, "No hay datos biológicos")
  if (nrow(datos_base$densidad_energetica) == 0) errores <- c(errores, "No hay datos de densidad energética")
  if (nrow(datos_base$chinook_diet_age) == 0)    errores <- c(errores, "No hay datos de dieta")

  cols_requeridas <- c("year", "age", "tw_g", "month", "species")
  cols_faltantes  <- cols_requeridas[!cols_requeridas %in% colnames(datos_base$chinook_biologico)]
  if (length(cols_faltantes) > 0)
    errores <- c(errores, paste("Columnas faltantes:", paste(cols_faltantes, collapse = ", ")))

  if (length(errores) > 0) stop(paste("Errores de validación:", paste(errores, collapse = ", ")))

  if (verbose) cat("✓ Validación completada\n\n")
  invisible(TRUE)
}

# ================================================================
# 2. PROCESAMIENTO DE DATOS
# ================================================================

procesar_datos_temporadas <- function(chinook_biologico, verbose = TRUE) {
  if (verbose) cat("=== PROCESANDO DATOS POR TEMPORADAS ===\n")

  chinook_temporadas <- chinook_biologico |>
    filter(!is.na(age), !is.na(year), !is.na(tw_g), !is.na(month),
           species == "Oncorhynchus tshawytscha",
           month %in% c(10, 11, 12, 1, 2, 3, 4)) |>
    mutate(
      temporada = case_when(
        month %in% c(10, 11, 12) ~ paste0(year, "/", substr(year + 1, 3, 4)),
        month %in% c(1, 2, 3, 4) ~ paste0(year - 1, "/", substr(year, 3, 4))
      ),
      year_temporada = case_when(
        month %in% c(10, 11, 12) ~ year + 1,
        month %in% c(1, 2, 3, 4) ~ year
      )
    )

  temporadas_disponibles <- sort(unique(chinook_temporadas$year_temporada))

  if (verbose) {
    cat("Datos procesados:", nrow(chinook_temporadas), "de", nrow(chinook_biologico), "originales\n")
    cat("Temporadas disponibles:", paste(temporadas_disponibles, collapse = ", "), "\n\n")
  }

  list(
    chinook_temporadas      = chinook_temporadas,
    temporadas_disponibles  = temporadas_disponibles
  )
}

# ================================================================
# 3. FUNCIONES AUXILIARES PARA MODELOS
# ================================================================

obtener_temperatura_por_año <- function(año_temporada, temperature_data) {
  temperature_data |>
    filter(year == año_temporada - 1) |>
    select(Day, Temperature)
}

crear_datos_dieta_por_edad <- function(chinook_diet_age, age_inicial, duracion = 365, verbose = TRUE) {
  dieta_edad <- chinook_diet_age |> filter(age == age_inicial)

  if (nrow(dieta_edad) == 0) stop(paste("No hay datos de dieta para edad", age_inicial))

  anchoveta_prop <- dieta_edad$prop_anchoveta_mean %||% 0
  sardina_prop   <- dieta_edad$prop_sardina_mean   %||% 0
  otros_prop     <- dieta_edad$prop_otros_mean     %||% 0

  if (verbose)
    cat(sprintf("Dieta edad %d - Anchoveta: %.1f%%, Sardina: %.1f%%, Otros: %.1f%%\n",
                age_inicial, anchoveta_prop * 100, sardina_prop * 100, otros_prop * 100))

  list(
    diet_data = data.frame(Day = 1:duracion, anchoveta = anchoveta_prop,
                           sardina = sardina_prop, otros = otros_prop),
    prey_energy_data = data.frame(Day = 1:duracion, anchoveta = 5553,
                                  sardina = 5636, otros = 4670),
    indigestible_data = data.frame(Day = 1:duracion, anchoveta = 0.089,
                                   sardina = 0.089, otros = 0.089),
    prey_names   = c("anchoveta", "sardina", "otros"),
    proportions  = list(anchoveta = anchoveta_prop, sardina = sardina_prop, otros = otros_prop)
  )
}

configurar_parametros_modelo <- function(chinook_temporadas, densidad_energetica,
                                         temporada_year, age_inicial, age_final,
                                         umbral_datos = 5, verbose = TRUE) {
  if (verbose)
    cat(sprintf("=== Configurando modelo: %d - Transición %d → %d ===\n",
                temporada_year, age_inicial, age_final))

  pesos_observados <- chinook_temporadas |>
    filter(year_temporada == temporada_year, age == age_final) |> pull(tw_g)

  pesos_iniciales  <- chinook_temporadas |>
    filter(year_temporada == temporada_year, age == age_inicial) |> pull(tw_g)

  if (length(pesos_observados) >= umbral_datos) {
    pesos_para_modelo <- pesos_observados
    data_source       <- "específica"
  } else {
    pesos_para_modelo <- chinook_temporadas |> filter(age == age_final) |> pull(tw_g)
    data_source       <- "histórica"
  }

  if (length(pesos_para_modelo) == 0) return(NULL)

  params_energia <- densidad_energetica[densidad_energetica$edad == age_inicial, ]
  if (nrow(params_energia) == 0) return(NULL)

  if (age_inicial == 1) {
    peso_inicial        <- params_energia$peso_inicial
    peso_inicial_source <- "parámetro fijo"
  } else if (length(pesos_iniciales) > 0) {
    peso_inicial        <- mean(pesos_iniciales)
    peso_inicial_source <- "específico temporal"
  } else {
    pesos_historicos    <- chinook_temporadas |> filter(age == age_inicial) |> pull(tw_g)
    peso_inicial        <- ifelse(length(pesos_historicos) > 0,
                                  mean(pesos_historicos), params_energia$peso_inicial)
    peso_inicial_source <- ifelse(length(pesos_historicos) > 0, "histórico", "parámetro fijo")
  }

  if (verbose) {
    cat(paste("Datos:", data_source, "- N:", length(pesos_para_modelo), "\n"))
    cat(paste("Peso inicial:", round(peso_inicial), "g (", peso_inicial_source, ")\n"))
    cat(paste("Target promedio:", round(mean(pesos_para_modelo)), "g\n\n"))
  }

  list(
    pesos_para_modelo   = pesos_para_modelo,
    peso_inicial        = peso_inicial,
    peso_inicial_source = peso_inicial_source,
    params_energia      = params_energia,
    data_source         = data_source
  )
}

# ================================================================
# 4. EJECUCIÓN DE MODELOS
# ================================================================

ejecutar_modelo_individual <- function(chinook_params, temperature_base, parametros_modelo,
                                       chinook_diet_age, temporada_year, age_inicial, age_final,
                                       verbose = TRUE) {
  pesos_para_modelo <- parametros_modelo$pesos_para_modelo
  peso_inicial      <- parametros_modelo$peso_inicial
  params_energia    <- parametros_modelo$params_energia

  densidad_ini <- params_energia$densidad_mix_ini
  densidad_end <- ifelse(is.na(params_energia$densidad_mix_end),
                         params_energia$pred_edad_end,
                         params_energia$densidad_mix_end)

  if (any(is.na(c(peso_inicial, densidad_ini, densidad_end)))) {
    if (verbose) cat("Parámetros energéticos incompletos\n")
    return(NULL)
  }

  if (verbose) {
    cat(sprintf("Ejecutando: %d transición %d → %d\n", temporada_year, age_inicial, age_final))
    cat(sprintf("Peso inicial: %d g - Densidades: %d / %d J/g\n",
                round(peso_inicial), round(densidad_ini), round(densidad_end)))
  }

  temperature_year    <- obtener_temperatura_por_año(temporada_year, temperature_base)
  duracion_year       <- nrow(temperature_year)
  datos_dieta_year    <- crear_datos_dieta_por_edad(chinook_diet_age, age_inicial, duracion_year, verbose)

  bio_obj <- Bioenergetic(
    species_info       = chinook_params$species_info,
    species_params     = chinook_params$life_stages$adult,
    environmental_data = list(temperature = temperature_year),
    diet_data          = list(
      proportions  = datos_dieta_year$diet_data,
      energies     = datos_dieta_year$prey_energy_data,
      indigestible = datos_dieta_year$indigestible_data,
      prey_names   = datos_dieta_year$prey_names
    ),
    simulation_settings = list(
      initial_weight = peso_inicial,
      duration       = duracion_year,
      oxycal         = 13560
    )
  )

  bio_obj$species_params <- set_parameter_value(bio_obj$species_params, "ED_ini", densidad_ini)
  bio_obj$species_params <- set_parameter_value(bio_obj$species_params, "ED_end", densidad_end)
  bio_obj$species_params <- set_parameter_value(bio_obj$species_params, "RK5",    0.00001)

  tryCatch({
    modelo <- run_fb4(
      x                = bio_obj,
      fit_to           = "Weight",
      strategy         = "mle",
      observed_weights = pesos_para_modelo,
      first_day        = 1,
      last_day         = duracion_year,
      backend          = "tmb",
      verbose          = FALSE
    )

    consumption_total    <- get_consumption_uncertainty(modelo)
    consumo_total_kg     <- consumption_total$estimate / 1000
    consumo_min_kg       <- (consumption_total$ci_lower %||% NA) / 1000
    consumo_max_kg       <- (consumption_total$ci_upper %||% NA) / 1000
    consumo_se_kg        <- (consumption_total$se      %||% NA) / 1000

    props <- datos_dieta_year$proportions

    modelo$metadata <- list(
      temporada           = temporada_year,
      transicion          = paste(age_inicial, age_final, sep = "→"),
      data_source         = parametros_modelo$data_source,
      peso_inicial        = peso_inicial,
      peso_inicial_source = parametros_modelo$peso_inicial_source,
      densidad_ini        = densidad_ini,
      densidad_end        = densidad_end,
      consumo_total_kg    = consumo_total_kg,
      consumo_min_kg      = consumo_min_kg,
      consumo_max_kg      = consumo_max_kg,
      consumo_se_kg       = consumo_se_kg,
      consumo_anchoveta_kg  = consumo_total_kg * props$anchoveta,
      consumo_anchoveta_min = (consumo_min_kg  %||% NA) * props$anchoveta,
      consumo_anchoveta_max = (consumo_max_kg  %||% NA) * props$anchoveta,
      consumo_sardina_kg    = consumo_total_kg * props$sardina,
      consumo_sardina_min   = (consumo_min_kg  %||% NA) * props$sardina,
      consumo_sardina_max   = (consumo_max_kg  %||% NA) * props$sardina,
      consumo_otros_kg      = consumo_total_kg * props$otros,
      consumo_otros_min     = (consumo_min_kg  %||% NA) * props$otros,
      consumo_otros_max     = (consumo_max_kg  %||% NA) * props$otros,
      prop_anchoveta        = props$anchoveta,
      prop_sardina          = props$sardina,
      prop_otros            = props$otros
    )

    if (verbose) {
      cat(sprintf("✓ Consumo total: %.1f kg\n  - Anchoveta: %.1f kg\n  - Sardina: %.1f kg\n  - Otros: %.1f kg\n\n",
                  consumo_total_kg,
                  modelo$metadata$consumo_anchoveta_kg,
                  modelo$metadata$consumo_sardina_kg,
                  modelo$metadata$consumo_otros_kg))
    }

    modelo

  }, error = function(e) {
    if (verbose) cat(paste("✗ Error:", e$message, "\n"))
    NULL
  })
}

ejecutar_modelos_batch <- function(datos_procesados, datos_base,
                                   temporadas_seleccionadas  = NULL,
                                   transiciones_seleccionadas = 1:4,
                                   umbral_datos = 5, verbose = TRUE) {
  if (verbose) cat("=== EJECUTANDO MODELOS POR LOTES ===\n")

  if (is.null(temporadas_seleccionadas))
    temporadas_seleccionadas <- datos_procesados$temporadas_disponibles

  todos_los_modelos  <- list()
  contador_exitosos  <- 0
  contador_fallidos  <- 0

  for (temporada in temporadas_seleccionadas) {
    if (verbose) cat(paste("--- TEMPORADA", temporada, "---\n"))

    for (transicion in transiciones_seleccionadas) {
      age_inicial <- transicion
      age_final   <- transicion + 1

      parametros <- configurar_parametros_modelo(
        datos_procesados$chinook_temporadas, datos_base$densidad_energetica,
        temporada, age_inicial, age_final, umbral_datos, verbose = FALSE
      )

      if (!is.null(parametros)) {
        modelo <- ejecutar_modelo_individual(
          datos_base$chinook_params, datos_base$temperature_base,
          parametros, datos_base$chinook_diet_age,
          temporada, age_inicial, age_final, verbose = verbose
        )

        if (!is.null(modelo)) {
          nombre_modelo <- paste0("T", temporada, "_", transicion, "to", transicion + 1)
          todos_los_modelos[[nombre_modelo]] <- modelo
          contador_exitosos <- contador_exitosos + 1
        } else {
          contador_fallidos <- contador_fallidos + 1
        }
      } else {
        contador_fallidos <- contador_fallidos + 1
        if (verbose) cat(paste("✗ Sin datos para transición", transicion, "\n"))
      }
    }
  }

  if (verbose) {
    tasa <- round(contador_exitosos / (contador_exitosos + contador_fallidos) * 100, 1)
    cat(sprintf("\n=== RESUMEN EJECUCIÓN ===\nExitosos: %d | Fallidos: %d | Tasa: %.1f%%\n\n",
                contador_exitosos, contador_fallidos, tasa))
  }

  list(
    modelos       = todos_los_modelos,
    configuracion = list(
      n_exitosos            = contador_exitosos,
      n_fallidos            = contador_fallidos,
      temporadas_ejecutadas = temporadas_seleccionadas,
      transiciones_ejecutadas = transiciones_seleccionadas
    )
  )
}

# ================================================================
# 5. TABLAS DE RESULTADOS
# ================================================================

generar_tabla_consumo_individual <- function(resultados_modelos, verbose = TRUE) {
  if (verbose) cat("=== GENERANDO TABLA CONSUMO INDIVIDUAL ===\n")

  modelos_metrics <- compare_scenarios(resultados_modelos$modelos,
                                       metrics = c("consumption", "p_value"))

  datos_base_df <- modelos_metrics$scenario_data |>
    filter(converged == TRUE) |>
    mutate(
      year            = as.numeric(gsub("T(\\d{4})_.*", "\\1", scenario)),
      transicion      = gsub("T\\d{4}_(.+)", "\\1", scenario),
      age_inicial     = as.numeric(gsub("(\\d+)to\\d+", "\\1", transicion)),
      age_final       = as.numeric(gsub("\\d+to(\\d+)", "\\1", transicion)),
      transicion_label = paste(age_inicial, "→", age_final),
      consumption_kg  = consumption_est / 1000,
      consumption_se_kg = ifelse(!is.na(consumption_se), consumption_se / 1000, NA)
    )

  modelos_metadata <- purrr::map_dfr(datos_base_df$scenario, function(scenario_name) {
    if (!scenario_name %in% names(resultados_modelos$modelos)) return(NULL)
    m <- resultados_modelos$modelos[[scenario_name]]$metadata
    data.frame(
      scenario              = scenario_name,
      data_source           = m$data_source           %||% NA,
      consumo_min_kg        = m$consumo_min_kg        %||% NA,
      consumo_max_kg        = m$consumo_max_kg        %||% NA,
      consumo_se_kg         = m$consumo_se_kg         %||% NA,
      consumo_anchoveta_kg  = m$consumo_anchoveta_kg  %||% NA,
      consumo_anchoveta_min = m$consumo_anchoveta_min %||% NA,
      consumo_anchoveta_max = m$consumo_anchoveta_max %||% NA,
      consumo_sardina_kg    = m$consumo_sardina_kg    %||% NA,
      consumo_sardina_min   = m$consumo_sardina_min   %||% NA,
      consumo_sardina_max   = m$consumo_sardina_max   %||% NA,
      consumo_otros_kg      = m$consumo_otros_kg      %||% NA,
      consumo_otros_min     = m$consumo_otros_min     %||% NA,
      consumo_otros_max     = m$consumo_otros_max     %||% NA,
      prop_anchoveta        = m$prop_anchoveta        %||% NA,
      prop_sardina          = m$prop_sardina          %||% NA,
      prop_otros            = m$prop_otros            %||% NA,
      stringsAsFactors = FALSE
    )
  })

  consumo_individual <- datos_base_df |>
    left_join(modelos_metadata, by = "scenario") |>
    select(scenario, year, transicion_label, age_inicial, age_final,
           consumption_kg, consumption_se_kg,
           consumo_min_kg, consumo_max_kg,
           p_value_est, p_value_se, data_source,
           consumo_anchoveta_kg, consumo_anchoveta_min, consumo_anchoveta_max,
           consumo_sardina_kg, consumo_sardina_min, consumo_sardina_max,
           consumo_otros_kg, consumo_otros_min, consumo_otros_max,
           prop_anchoveta, prop_sardina, prop_otros) |>
    arrange(year, age_inicial)

  if (verbose) {
    cat("Tabla final:", nrow(consumo_individual), "filas\n")
    cat("Años:", paste(sort(unique(consumo_individual$year)), collapse = ", "), "\n\n")
  }

  consumo_individual
}

# ================================================================
# 6. VISUALIZACIONES
# ================================================================

# ── Paletas compartidas ──────────────────────────────────────────
.colores_transicion <- c(
  "1 → 2" = "#8c9aa5ff",
  "2 → 3" = "#87CEEB",
  "3 → 4" = "#4682B4",
  "4 → 5" = "#191970"
)

.colores_presas <- c(
  "Anchovy"  = "#E76F51",
  "Sardine"  = "#4472C4",
  "Others"   = "#8D6E63"
)

.theme_base <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size - 4, hjust = 0.5, color = "gray40"),
      axis.title       = element_text(size = base_size, face = "bold", color = "black"),
      axis.text        = element_text(size = base_size - 2, color = "black"),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "gray30", fill = NA, linewidth = 0.5)
    )
}

# ── Consumo por edad (barras apiladas) ───────────────────────────
plot_consumo_barras_apiladas <- function(resultados_modelos, titulo = NULL) {
  modelos_metrics <- compare_scenarios(resultados_modelos$modelos,
                                       metrics = c("consumption", "p_value"))

  datos <- modelos_metrics$scenario_data |>
    filter(converged == TRUE) |>
    mutate(
      año            = as.numeric(gsub("T(\\d{4})_.*", "\\1", scenario)) - 1,
      age_inicial    = as.numeric(gsub(".*_(\\d+)to\\d+", "\\1", scenario)),
      transicion_label = paste(age_inicial, "→", age_inicial + 1),
      consumption_kg = consumption_est / 1000
    ) |>
    arrange(año, age_inicial)

  ggplot(datos, aes(x = factor(año), y = consumption_kg, fill = transicion_label)) +
    geom_col(position = "stack", alpha = 0.85, color = "white", linewidth = 0.3) +
    scale_fill_manual(name = "Age Transition", values = .colores_transicion) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Annual Consumption (kg)\n",
                       labels = scales::comma_format(),
                       expand = expansion(mult = c(0, 0.05))) +
    .theme_base() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()) +
    labs(title    = titulo %||% "Annual Consumption by Age Transition",
         subtitle = "Stacked bars show contribution of each age transition to total consumption")
}

# ── Consumo por presa (barras apiladas) ─────────────────────────
plot_consumo_presas_barras_apiladas <- function(consumo_individual, titulo = NULL) {
  datos_largo <- consumo_individual |>
    group_by(year) |>
    summarise(Anchovy = sum(consumo_anchoveta_kg, na.rm = TRUE),
              Sardine = sum(consumo_sardina_kg,   na.rm = TRUE),
              Others  = sum(consumo_otros_kg,     na.rm = TRUE),
              .groups = "drop") |>
    pivot_longer(cols = c(Anchovy, Sardine, Others),
                 names_to = "Presa", values_to = "Consumo_kg") |>
    mutate(Presa = factor(Presa, levels = c("Others", "Sardine", "Anchovy")))

  ggplot(datos_largo, aes(x = factor(year), y = Consumo_kg, fill = Presa)) +
    geom_col(position = "stack", alpha = 0.85, color = "white", linewidth = 0.3) +
    scale_fill_manual(name = "Prey Type", values = .colores_presas) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Annual Consumption (kg)\n",
                       labels = scales::comma_format(),
                       expand = expansion(mult = c(0, 0.05))) +
    .theme_base() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()) +
    labs(title    = titulo %||% "Annual Consumption by Prey Type",
         subtitle = "Distribution of total consumption among anchovy, sardine and other prey")
}

# ── Consumo: boxplot por transición ─────────────────────────────
plot_consumo_boxplot_por_transicion <- function(consumo_individual, titulo = NULL) {
  ggplot(consumo_individual,
         aes(x = transicion_label, y = consumption_kg, fill = transicion_label)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = .colores_transicion) +
    scale_y_continuous(name = "Annual Consumption (kg)", labels = scales::comma_format()) +
    scale_x_discrete(name = "Age Transition") +
    .theme_base() +
    theme(legend.position = "none") +
    labs(title    = titulo %||% "Distribution of Annual Consumption by Age Transition",
         subtitle = "Boxplots show median, quartiles and individual data points across years")
}

# ── Consumo: barras de promedio por transición ───────────────────
plot_consumo_barras_promedio_transicion <- function(consumo_individual, titulo = NULL) {
  datos_resumen <- consumo_individual |>
    group_by(transicion_label, age_inicial) |>
    summarise(consumo_promedio  = mean(consumption_kg, na.rm = TRUE),
              consumo_se        = sd(consumption_kg, na.rm = TRUE) / sqrt(n()),
              n_observaciones   = n(),
              .groups = "drop") |>
    arrange(age_inicial)

  ggplot(datos_resumen, aes(x = transicion_label, y = consumo_promedio, fill = transicion_label)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_errorbar(aes(ymin = consumo_promedio - consumo_se,
                      ymax = consumo_promedio + consumo_se),
                  width = 0.3, alpha = 0.7, linewidth = 0.8) +
    geom_text(aes(label = paste("n =", n_observaciones)), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = .colores_transicion) +
    scale_y_continuous(name = "Mean Annual Consumption (kg)",
                       labels = scales::comma_format(),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_x_discrete(name = "Age Transition") +
    .theme_base() +
    theme(legend.position = "none", panel.grid.major.x = element_blank()) +
    labs(title    = titulo %||% "Mean Annual Consumption by Age Transition",
         subtitle = "Error bars show standard error; numbers indicate sample size")
}

# ── Consumo: violin por transición ───────────────────────────────
plot_consumo_violin_por_transicion <- function(consumo_individual, titulo = NULL) {
  ggplot(consumo_individual,
         aes(x = transicion_label, y = consumption_kg, fill = transicion_label)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", color = "black", stroke = 1) +
    scale_fill_manual(values = .colores_transicion) +
    scale_y_continuous(name = "Annual Consumption (kg)", labels = scales::comma_format()) +
    scale_x_discrete(name = "Age Transition") +
    .theme_base() +
    theme(legend.position = "none") +
    labs(title    = titulo %||% "Distribution Density of Annual Consumption by Age Transition",
         subtitle = "Violin plots show distribution density; diamonds show means; boxes show quartiles")
}

# ── Consumo: líneas temporales por transición ────────────────────
plot_consumo_lineas_por_ano_transicion <- function(consumo_individual, titulo = NULL) {
  ggplot(consumo_individual,
         aes(x = year, y = consumption_kg, color = transicion_label, shape = transicion_label)) +
    geom_line(aes(group = transicion_label), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(name = "Age Transition", values = .colores_transicion) +
    scale_shape_manual(name = "Age Transition", values = c(19, 17, 15, 18)) +
    scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(name = "Annual Consumption (kg)", labels = scales::comma_format()) +
    .theme_base() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title    = titulo %||% "Temporal Trends in Annual Consumption by Age Transition",
         subtitle = "Lines show consumption patterns across years for each age transition")
}

# ── Consumo: panel combinado ─────────────────────────────────────
plot_panel_consumo_edades <- function(consumo_individual, ncol = 2, titulo_general = NULL) {
  wrap_plots(
    plot_consumo_boxplot_por_transicion(consumo_individual, "A) Distribution by Transition"),
    plot_consumo_barras_promedio_transicion(consumo_individual, "B) Mean Consumption"),
    plot_consumo_violin_por_transicion(consumo_individual, "C) Density Distribution"),
    plot_consumo_lineas_por_ano_transicion(consumo_individual, "D) Temporal Trends"),
    ncol = ncol
  ) +
    plot_annotation(
      title    = titulo_general %||% "Comprehensive Analysis of Annual Consumption by Age Transition",
      subtitle = "Multiple perspectives on consumption patterns across age transitions and years",
      theme    = theme(plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
                       plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"))
    )
}

# ── Presas: boxplot ───────────────────────────────────────────────
plot_consumo_presas_boxplot <- function(consumo_individual, titulo = NULL) {
  datos_largo <- consumo_individual |>
    select(year, consumo_anchoveta_kg, consumo_sardina_kg, consumo_otros_kg) |>
    pivot_longer(cols = starts_with("consumo_"),
                 names_to = "tipo_presa", values_to = "consumo_kg") |>
    mutate(tipo_presa = recode(tipo_presa,
                               "consumo_anchoveta_kg" = "Anchovy",
                               "consumo_sardina_kg"   = "Sardine",
                               "consumo_otros_kg"     = "Others"),
           tipo_presa = factor(tipo_presa, levels = c("Anchovy", "Sardine", "Others")))

  ggplot(datos_largo, aes(x = tipo_presa, y = consumo_kg, fill = tipo_presa)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = .colores_presas) +
    scale_y_continuous(name = "Annual Consumption (kg)", labels = scales::comma_format()) +
    scale_x_discrete(name = "Prey Type") +
    .theme_base() +
    theme(legend.position = "none") +
    labs(title    = titulo %||% "Distribution of Annual Consumption by Prey Type",
         subtitle = "Boxplots show total consumption across all age transitions")
}

# ── P-values: serie temporal ─────────────────────────────────────
plot_pvalues_por_transicion_tiempo <- function(consumo_individual, titulo = NULL) {
  if (!"p_value_est" %in% colnames(consumo_individual))
    stop("No se encontraron datos de p-value en consumo_individual")

  datos_pvalue <- consumo_individual |>
    filter(!is.na(p_value_est)) |>
    mutate(p_value_lower = pmax(0, p_value_est - p_value_se, na.rm = TRUE),
           p_value_upper = pmin(1, p_value_est + p_value_se, na.rm = TRUE))

  ggplot(datos_pvalue,
         aes(x = year, y = p_value_est, color = transicion_label, shape = transicion_label)) +
    geom_line(aes(group = transicion_label), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbar(aes(ymin = p_value_lower, ymax = p_value_upper), width = 0.2, alpha = 0.6) +
    scale_color_manual(name = "Age Transition", values = .colores_transicion) +
    scale_shape_manual(name = "Age Transition", values = c(19, 17, 15, 18)) +
    scale_x_continuous(name = "", breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(name = "P-value (Feeding Level)", limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    .theme_base() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title    = titulo %||% "Feeding Level (P-value) Trends by Age Transition",
         subtitle = "Lines show temporal patterns in estimated feeding intensity; error bars show ± SE")
}

# ── P-values: boxplot por transición ────────────────────────────
plot_pvalues_boxplot_por_transicion <- function(consumo_individual, titulo = NULL) {
  datos_pvalue <- consumo_individual |> filter(!is.na(p_value_est))

  ggplot(datos_pvalue, aes(x = transicion_label, y = p_value_est, fill = transicion_label)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = .colores_transicion) +
    scale_y_continuous(name = "P-value (Feeding Level)", limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(name = "Age Transition") +
    .theme_base() +
    theme(legend.position = "none") +
    labs(title    = titulo %||% "Distribution of Feeding Levels (P-values) by Age Transition",
         subtitle = "Boxplots show median, quartiles and individual feeding intensity estimates")
}

# ================================================================
# 7. ANÁLISIS DE SENSITIVIDAD TEMPERATURA
# ================================================================

ejecutar_sensitivity_analysis_batch <- function(resultados_modelos,
                                                temperatures   = seq(8, 18, by = 1),
                                                p_values_base  = c(0.15, 0.25, 0.5, 0.75),
                                                verbose = TRUE) {
  cat("=== ANÁLISIS DE SENSITIVIDAD TEMPERATURA-CRECIMIENTO ===\n")

  modelos_metrics  <- compare_scenarios(resultados_modelos$modelos,
                                        metrics = c("consumption", "p_value"))
  modelos_exitosos <- names(resultados_modelos$modelos)

  if (verbose) {
    cat("Modelos a analizar:", length(modelos_exitosos), "\n")
    cat("Rango de temperatura:", min(temperatures), "°C a", max(temperatures), "°C\n\n")
  }

  sensitivity_results <- purrr::imap(
    setNames(modelos_exitosos, modelos_exitosos),
    function(modelo_name, i) {
      cat("[", which(modelos_exitosos == modelo_name), "/", length(modelos_exitosos),
          "] Procesando:", modelo_name, "\n")

      p_value_modelo <- modelos_metrics$scenario_data$p_value_est[
        modelos_metrics$scenario_data$scenario == modelo_name
      ]

      if (length(p_value_modelo) == 0 || is.na(p_value_modelo)) {
        cat("  ⚠️ P-value no encontrado, saltando\n\n")
        return(NULL)
      }

      p_values_completo <- sort(unique(c(p_values_base, p_value_modelo)))
      modelo    <- resultados_modelos$modelos[[modelo_name]]
      bio_obj   <- modelo$bioenergetic_object
      metadata  <- modelo$metadata

      if (is.null(bio_obj)) {
        cat("  ⚠️ Objeto bioenergético no encontrado\n\n")
        return(NULL)
      }

      tryCatch({
        sens_data <- analyze_growth_temperature_sensitivity(
          bio_obj      = bio_obj,
          p_values     = p_values_completo,
          temperatures = temperatures,
          verbose      = FALSE
        )
        sens_data$modelo        <- modelo_name
        sens_data$temporada     <- metadata$temporada
        sens_data$transicion    <- metadata$transicion
        sens_data$p_value_modelo <- p_value_modelo
        cat("  ✓ Análisis completado\n\n")
        sens_data
      }, error = function(e) {
        cat("  ✗ Error:", e$message, "\n\n")
        NULL
      })
    }
  )

  sensitivity_results <- Filter(Negate(is.null), sensitivity_results)

  cat("=== RESUMEN FINAL ===\n")
  cat("Modelos analizados:", length(sensitivity_results), "de", length(modelos_exitosos), "\n")

  list(
    sensitivity_data = sensitivity_results,
    modelos_metrics  = modelos_metrics,
    configuracion    = list(
      temperatures_range     = range(temperatures),
      p_values_base          = p_values_base,
      n_modelos_procesados   = length(sensitivity_results),
      fecha_analisis         = Sys.time()
    )
  )
}

plot_growth_temperature_sensitivity_ggplot <- function(sensitivity_data,
                                                       ylim   = c(-0.002, 0.008),
                                                       xlim   = c(1, 20),
                                                       titulo = NULL) {
  valid_results <- sensitivity_data[!is.na(sensitivity_data$daily_growth_rate), ]
  if (nrow(valid_results) == 0) stop("No valid sensitivity results to plot")

  feeding_levels     <- sort(unique(valid_results$feeding_pct), decreasing = TRUE)
  valid_results$p_label  <- paste("P =", round(valid_results$feeding_pct, 2))
  valid_results$p_factor <- factor(valid_results$p_label,
                                   levels = paste("P =", round(feeding_levels, 2)))

  optimal_temp <- valid_results[valid_results$feeding_pct == max(valid_results$feeding_pct), ] |>
    (\(d) d$temperature[which.max(d$daily_growth_rate)])()

  ggplot(valid_results, aes(x = temperature, y = daily_growth_rate)) +
    geom_line(aes(color = p_factor, group = p_factor), linewidth = 1.2, alpha = 0.8) +
    geom_point(aes(color = p_factor, shape = p_factor), size = 2, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = optimal_temp, linetype = "dashed", color = "red", alpha = 0.7) +
    annotate("text", x = optimal_temp + 1,
             y = max(valid_results$daily_growth_rate, na.rm = TRUE) * 0.8,
             label = paste("Optimal\n", optimal_temp, "°C"),
             color = "red", size = 3, hjust = 0) +
    scale_color_manual(name = "Feeding Level", values = rainbow(length(feeding_levels))) +
    scale_shape_manual(name = "Feeding Level",
                       values = rep(c(19, 1, 2, 0, 17), length.out = length(feeding_levels))) +
    scale_x_continuous(name = "Temperature (°C)", breaks = seq(0, 20, by = 5)) +
    scale_y_continuous(name = "Daily Growth Rate (g/g/day)") +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    theme_bw(base_size = 14) +
    theme(legend.title = element_text(size = 14, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "lightgray", linetype = "dotted")) +
    labs(title = titulo %||% "Temperature & Feeding Effects on Growth")
}

# Helper interno para extraer y graficar subconjuntos de sensitivity
.plot_sensitivity_subset <- function(sensitivity_results, modelo_names, ncol, titulo_general) {
  plot_list <- purrr::compact(purrr::imap(
    setNames(modelo_names, modelo_names),
    function(nombre, key) {
      sens_data <- sensitivity_results$sensitivity_data[[nombre]]
      if (is.null(sens_data)) return(NULL)
      tryCatch(
        plot_growth_temperature_sensitivity_ggplot(sens_data, titulo = key),
        error = function(e) NULL
      )
    }
  ))

  if (length(plot_list) == 0) stop("No se pudieron crear gráficos")

  wrap_plots(plot_list, ncol = ncol) +
    plot_annotation(
      title    = titulo_general,
      subtitle = "Sensitivity analysis",
      theme    = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
}

plot_por_anos_una_transicion_ggplot <- function(sensitivity_results,
                                                transicion = "1to2",
                                                ylim = c(-0.002, 0.008),
                                                ncol = 2) {
  modelos <- grep(paste0("_", transicion), names(sensitivity_results$sensitivity_data), value = TRUE)
  if (length(modelos) == 0) stop("No se encontraron modelos para transición: ", transicion)

  anos        <- sort(as.numeric(gsub("T(\\d{4})_.*", "\\1", modelos)))
  nombres_ord <- paste0("T", anos, "_", transicion)
  etiquetas   <- setNames(as.character(anos), nombres_ord)

  trans_nums  <- strsplit(transicion, "to")[[1]]
  titulo_gral <- paste("Temperature & Feeding Effects - Transition", trans_nums[1], "->", trans_nums[2])

  .plot_sensitivity_subset(sensitivity_results, nombres_ord[nombres_ord %in% names(sensitivity_results$sensitivity_data)],
                           ncol, titulo_gral) |>
    (\(p) p + plot_annotation(subtitle = "By year"))()
}

plot_por_transiciones_un_ano_ggplot <- function(sensitivity_results,
                                                ano  = 2015,
                                                ylim = c(-0.002, 0.008),
                                                ncol = 2) {
  trans_orden <- c("1to2", "2to3", "3to4", "4to5")
  modelos     <- grep(paste0("T", ano, "_"), names(sensitivity_results$sensitivity_data), value = TRUE)
  if (length(modelos) == 0) stop("No se encontraron modelos para año: ", ano)

  trans_disponibles <- gsub(paste0("T", ano, "_(.+)"), "\\1", modelos)
  trans_ordenadas   <- trans_orden[trans_orden %in% trans_disponibles]
  nombres_ord       <- paste0("T", ano, "_", trans_ordenadas)

  etiquetas <- purrr::set_names(
    purrr::map_chr(trans_ordenadas, ~ paste("Transition", gsub("(\\d)to(\\d)", "\\1 -> \\2", .x))),
    nombres_ord
  )

  .plot_sensitivity_subset(sensitivity_results,
                           nombres_ord[nombres_ord %in% names(sensitivity_results$sensitivity_data)],
                           ncol,
                           paste("Temperature & Feeding Effects - Year", ano))
}

# ================================================================
# 8. ESTADÍSTICAS DESCRIPTIVAS
# ================================================================

resumir_estadisticas <- function(consumo_individual) {
  cat("\n=== P-values por edad ===\n")
  consumo_individual |>
    reframe(p_mean = mean(p_value_est),
            p_min  = min(p_value_est),
            p_max  = max(p_value_est),
            .by    = age_inicial) |>
    arrange(age_inicial) |> print()

  cat(sprintf("\nRango global p-values: %.4f - %.4f\n",
              min(consumo_individual$p_value_est),
              max(consumo_individual$p_value_est)))

  cat("\n=== Consumo por edad (kg) ===\n")
  consumo_individual |>
    reframe(n        = n(),
            mean_kg  = round(mean(consumption_kg), 1),
            min_kg   = round(min(consumo_min_kg),  1),
            max_kg   = round(max(consumo_max_kg),  1),
            cv_pct   = round(sd(consumption_kg) / mean(consumption_kg) * 100, 1),
            .by      = age_inicial) |>
    arrange(age_inicial) |> print()

  cat("\n=== Proporción de presas por edad ===\n")
  consumo_individual |>
    reframe(anchovy_pct = round(mean(prop_anchoveta) * 100, 1),
            sardine_pct = round(mean(prop_sardina)   * 100, 1),
            others_pct  = round(mean(prop_otros)     * 100, 1),
            .by         = age_inicial) |>
    arrange(age_inicial) |> print()
}

# ================================================================
# 9. FLUJO PRINCIPAL
# ================================================================
# Cambiar force = TRUE para forzar recálculo desde ese paso en adelante.
# Si el checkpoint existe y force = FALSE, se carga directo desde disco.

FORCE_DATOS       <- TRUE
FORCE_MODELOS     <- TRUE
FORCE_CONSUMO     <- TRUE
FORCE_SENSITIVITY <- TRUE

PATHS <- list(
  cp1  = "data/data_raw/bioenergetic-model/checkpoint_01_datos.rds",
  cp2  = "data/data_raw/bioenergetic-model/checkpoint_02_modelos.rds",
  cp3  = "data/data_raw/bioenergetic-model/checkpoint_03_consumo.rds",
  cp4  = "data/data_raw/bioenergetic-model/checkpoint_04_sensitivity.rds",
  csv  = "data/data_raw/bioenergetic-model/resultados_total_consumption_by_age_ISOTOPES.csv"
)

# ── CHECKPOINT 1: Datos base + temporadas ────────────────────────
cp1 <- with_checkpoint(PATHS$cp1, force = FORCE_DATOS, {
  datos_base      <- cargar_datos_base()
  validar_datos(datos_base)
  datos_procesados <- procesar_datos_temporadas(datos_base$chinook_biologico)
  list(datos_base = datos_base, datos_procesados = datos_procesados)
})

datos_base       <- cp1$datos_base
datos_procesados <- cp1$datos_procesados

# ── CHECKPOINT 2: Modelos bioenergéticos ─────────────────────────
resultados_modelos <- with_checkpoint(PATHS$cp2, force = FORCE_MODELOS, {
  ejecutar_modelos_batch(
    datos_procesados, datos_base,
    temporadas_seleccionadas   = seq(2015, 2023),
    transiciones_seleccionadas = seq(1, 4),
    umbral_datos = 5
  )
})


saveRDS(resultados_modelos, "data/data_raw/bioenergetic-model/modelos_total_by_age_ISOTOPES.rds")


# ── CHECKPOINT 3: Tabla de consumo individual ────────────────────
consumo_individual <- with_checkpoint(PATHS$cp3, force = FORCE_CONSUMO, {
  ci <- generar_tabla_consumo_individual(resultados_modelos)
  write.csv(ci, PATHS$csv, row.names = FALSE)
  ci
})

# ── CHECKPOINT 4: Análisis de sensitividad ───────────────────────
sensitivity_results <- with_checkpoint(PATHS$cp4, force = FORCE_SENSITIVITY, {
  ejecutar_sensitivity_analysis_batch(
    resultados_modelos,
    temperatures  = seq(1, 20, by = 1),
    p_values_base = c(0.15, 0.25, 0.5, 0.75, 1)
  )
})

# ================================================================
# 10. VISUALIZACIONES FINALES
# ================================================================

dir.create("output", showWarnings = FALSE)

# Barras apiladas por edad
p1 <- plot_consumo_barras_apiladas(resultados_modelos)
ggsave("output/consumo_por_edad_ISOTOPES.png", p1, width = 12, height = 8, dpi = 300)

# Barras apiladas por presa
p2 <- plot_consumo_presas_barras_apiladas(consumo_individual)
ggsave("output/consumo_por_presa_ISOTOPES.png", p2, width = 12, height = 8, dpi = 300)

# P-values temporales
p3 <- plot_pvalues_por_transicion_tiempo(consumo_individual)
ggsave("output/pvalues_temporal_ISOTOPES.png", p3, width = 12, height = 8, dpi = 300)

# Panel completo de consumo
p_panel <- plot_panel_consumo_edades(consumo_individual)
ggsave("output/panel_consumo_ISOTOPES.png", p_panel, width = 16, height = 12, dpi = 300)

# Boxplots adicionales
p_presas    <- plot_consumo_presas_boxplot(consumo_individual)
p_pv_box    <- plot_pvalues_boxplot_por_transicion(consumo_individual)

# Sensitividad por transición (ejemplo)
p_sens_trans_1_2 <- plot_por_anos_una_transicion_ggplot(sensitivity_results, transicion = "1to2", ncol = 2)
p_sens_trans_2_3 <- plot_por_anos_una_transicion_ggplot(sensitivity_results, transicion = "2to3", ncol = 2)
p_sens_trans_3_4 <- plot_por_anos_una_transicion_ggplot(sensitivity_results, transicion = "3to4", ncol = 2)
p_sens_trans_4_5 <- plot_por_anos_una_transicion_ggplot(sensitivity_results, transicion = "4to5", ncol = 2)
p_sens_ano   <- plot_por_transiciones_un_ano_ggplot(sensitivity_results, ano = 2015, ncol = 2)

# Estadísticas descriptivas
resumir_estadisticas(consumo_individual)

cat("\n=== ANÁLISIS COMPLETADO ===\n")
cat("Checkpoints guardados en: data/data_raw/bioenergetic-model/\n")
cat("Gráficos guardados en: output/\n")





# ================================================================
# FIGURAS PARA MATERIAL SUPLEMENTARIO S4
# ================================================================

# ── Figura S4.1: Consumo vs temperatura (sensibilidad) ───────────
# Promedio de consumo anual por clase de edad a temperaturas fijas

# ── Figura S4.1: Consumo vs temperatura (sensibilidad) ───────────

sens_summary <- purrr::map_dfr(
  names(sensitivity_results$sensitivity_data),
  function(nombre) {
    d <- sensitivity_results$sensitivity_data[[nombre]]
    if (is.null(d)) return(NULL)
    age <- as.numeric(gsub("T\\d{4}_(\\d)to\\d+", "\\1", nombre))
    
    # Solo el p-value estimado por el modelo (no los otros)
    d |>
      filter(p_value == p_value_modelo) |>
      group_by(temperature) |>
      summarise(
        consumption_kg = mean(total_consumption / 1000, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(age_transition = paste0(age, " \u2192 ", age + 1))
  }
) |>
  group_by(temperature, age_transition) |>
  summarise(consumption_kg = mean(consumption_kg, na.rm = TRUE), .groups = "drop")

# ================================================================
# FIGURE S4.1 — Temperature sensitivity (B/W journal style)
# ================================================================

# ================================================================
# FIGURE S4.1 — Temperature sensitivity (B/W journal style)
# ================================================================

library(tidyverse)

# ---- Asegurar orden correcto de transiciones ----
sens_summary_plot <- sens_summary %>%
  mutate(
    age_transition = factor(
      age_transition,
      levels = c("1 → 2", "2 → 3", "3 → 4", "4 → 5")
    )
  )

# ---- Escala de grises ----
colores_transicion_bw <- c(
  "1 → 2" = "grey80",
  "2 → 3" = "grey60",
  "3 → 4" = "grey40",
  "4 → 5" = "grey20"
)

# ---- Tipos de línea ----
tipos_linea <- c(
  "1 → 2" = "solid",
  "2 → 3" = "dashed",
  "3 → 4" = "dotted",
  "4 → 5" = "dotdash"
)

# ---- Figura ----
fig_S4_1 <- ggplot(
  sens_summary_plot,
  aes(
    x = temperature,
    y = consumption_kg,
    color = age_transition,
    linetype = age_transition,
    group = age_transition
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = colores_transicion_bw,
    name = "Age transition"
  ) +
  scale_linetype_manual(
    values = tipos_linea,
    name = "Age transition"
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, 5),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Annual consumption (kg)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print(fig_S4_1)

# ---- Guardar PNG (rápido, buena calidad) ----
ggsave(
  "outputs/graficos/Figure_S4_1_temperature_sensitivity.png",
  plot = fig_S4_1,
  dpi = 300,
  width = 180,
  height = 120,
  units = "mm"
)

# ---- Guardar PDF (vectorial, editable) ----
ggsave(
  "outputs/graficos/Figure_S4_1_temperature_sensitivity.pdf",
  plot = fig_S4_1,
  width = 180,
  height = 120,
  units = "mm"
)

# ---- Guardar TIFF (formato journal) ----
ggsave(
  "outputs/graficos/Figure_S4_1_temperature_sensitivity.tiff",
  plot = fig_S4_1,
  dpi = 600,
  width = 180,
  height = 120,
  units = "mm",
  compression = "lzw"
)


# ── Figura S4.2: P-values estimados por transición y año ─────────
fig_S4_2 <- plot_pvalues_por_transicion_tiempo(consumo_individual,
               titulo = "Figure S4.2. Estimated p-values by age transition and year")

ggsave("output/FigS4_2_pvalues_temporal.png",
       fig_S4_2, width = 10, height = 6, dpi = 300)

# ── Figura S4.3: Pesos predichos vs observados ───────────────────
modelos_metrics <- compare_scenarios(resultados_modelos$modelos,
                                     metrics = c("consumption", "p_value"))

pred_obs <- modelos_metrics$scenario_data |>
  filter(converged == TRUE) |>
  mutate(
    age_transition = paste0(gsub("T\\d{4}_(\\d)to\\d", "\\1", scenario),
                            " \u2192 ",
                            gsub("T\\d{4}_\\dto(\\d)", "\\1", scenario))
  )

# fig_S4_3 <- ggplot(pred_obs,
#                    aes(x = initial_weight, y = final_weight_pred,
#                        color = age_transition)) +
#   geom_point(size = 3, alpha = 0.8) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
#   scale_color_manual(name = "Age Transition", values = .colores_transicion) +
#   scale_x_continuous(name = "Initial Weight (g)", labels = scales::comma_format()) +
#   scale_y_continuous(name = "Predicted Final Weight (g)", labels = scales::comma_format()) +
#   theme_classic(base_size = 14) +
#   labs(title    = "Figure S4.3",
#        subtitle = "Predicted vs. initial weight for all 36 age–year combinations")

# ggsave("output/FigS4_3_predicted_weights.png",
#        fig_S4_3, width = 8, height = 6, dpi = 300)

# cat("\n✅ Figuras S4 guardadas en output/\n")