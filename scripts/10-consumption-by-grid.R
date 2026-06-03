# Exploratory spatial analysis — individual consumption by temperature grid cell
#
# NOTE: This script is not part of the main analysis pipeline.
# It re-runs the bioenergetic model for each spatial grid cell using fixed
# p-values from script 08, producing consumption maps by age and year.
# Output is not used in downstream scripts.

library(fb4package)
library(dplyr)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# --- configuración ---

config <- list(
  años_procesar   = 2015:2023,
  edades_procesar = 1:4,
  max_grillas     = NULL,
  verbose         = TRUE
)

cat("Años a procesar:", paste(config$años_procesar, collapse = ", "), "\n")
cat("Edades a procesar:", paste(config$edades_procesar, collapse = ", "), "\n")
cat("Límite de grillas:", ifelse(is.null(config$max_grillas), "Sin límite", config$max_grillas), "\n")
cat("Verbose:", config$verbose, "\n\n")

# --- cargar datos base ---

# Cargar modelos bioenergéticos CON ISÓTOPOS
modelos <- readRDS("data/data_raw/bioenergetic-model/modelos_total_by_age_ISOTOPES.rds")
cat("Modelos ISOTOPES cargados:", length(modelos$modelos), "modelos\n")

# Cargar datos de temperatura por grilla
datos_temp <- readRDS("data/data_raw/temperature/total_data_temperature.rds")
cat("Datos temperatura cargados para años:", paste(names(datos_temp$series_grillas), collapse = ", "), "\n")

# Verificar disponibilidad de años
años_modelos <- unique(gsub("T(\\d{4})_.*", "\\1", names(modelos$modelos)))
años_temperatura <- names(datos_temp$series_grillas)
# años_disponibles <- intersect(config$años_procesar, intersect(años_modelos, años_temperatura))

años_retorno_disponibles <- intersect(años_modelos, 
                                       as.character(as.numeric(años_temperatura) + 1))
años_disponibles <- intersect(as.character(config$años_procesar), 
                               años_retorno_disponibles)

cat("Años disponibles para procesamiento:", paste(años_disponibles, collapse = ", "), "\n\n")

if(length(años_disponibles) == 0) {
  stop("No hay años disponibles que coincidan entre modelos y temperatura")
}

# --- extraer valores de P por año/edad ---

p_values_table <- data.frame()

for(año in años_disponibles) {
  for(edad in config$edades_procesar) {
    
    modelo_name <- paste0("T", año, "_", edad, "to", edad + 1)
    
    if(modelo_name %in% names(modelos$modelos)) {
      p_value <- modelos$modelos[[modelo_name]]$summary$p_estimate
      converged <- modelos$modelos[[modelo_name]]$summary$converged
      
      p_values_table <- rbind(p_values_table, data.frame(
        año = as.numeric(año),
        edad = edad,
        modelo_name = modelo_name,
        p_value = p_value,
        converged = converged,
        stringsAsFactors = FALSE
      ))
      
      if(config$verbose) {
        cat("  ", modelo_name, ": p =", round(p_value, 4), 
            "(convergencia:", converged, ")\n")
      }
    } else {
      if(config$verbose) cat("  ", modelo_name, ": NO ENCONTRADO\n")
    }
  }
}

cat("P-values extraídos:", nrow(p_values_table), "modelos\n")
cat("Modelos convergidos:", sum(p_values_table$converged), "\n\n")

# Filtrar solo modelos convergidos
p_values_table <- p_values_table[p_values_table$converged, ]

if(nrow(p_values_table) == 0) {
  stop("No hay modelos convergidos disponibles")
}

# --- preparar estructura de resultados ---

resultado_final <- list()

for(año in años_disponibles) {
  resultado_final[[paste0("año_", año)]] <- list()
  
  for(edad in config$edades_procesar) {
    resultado_final[[paste0("año_", año)]][[paste0("edad_", edad)]] <- data.frame()
  }
}

# --- procesamiento principal ---

total_combinaciones <- nrow(p_values_table)
contador_completados <- 0
contador_errores <- 0

for(i in 1:nrow(p_values_table)) {
  
  año_actual <- p_values_table$año[i]
  edad_actual <- p_values_table$edad[i]
  modelo_name <- p_values_table$modelo_name[i]
  p_value_actual <- p_values_table$p_value[i]
  
  if(config$verbose) {
    cat("\n--- PROCESANDO:", modelo_name, "(", i, "/", total_combinaciones, ") ---\n")
  }
  
  bio_obj_base <- modelos$modelos[[modelo_name]]$bioenergetic_object
  
  if(is.null(bio_obj_base)) {
    cat("ERROR: No se encontró objeto bioenergético para", modelo_name, "\n")
    contador_errores <- contador_errores + 1
    next
  }
  
  # año real de consumo
  año_str <- as.character(año_actual - 1)

  if(!año_str %in% names(datos_temp$series_grillas)) {
    cat("ERROR: No hay datos de temperatura para el año", año_actual, "\n")
    contador_errores <- contador_errores + 1
    next
  }
  
  grillas_año <- datos_temp$series_grillas[[año_str]]
  n_grillas_total <- length(grillas_año)
  
  # Aplicar límite de grillas si está configurado
  if(!is.null(config$max_grillas)) {
    grillas_procesar <- names(grillas_año)[1:min(config$max_grillas, n_grillas_total)]
    grillas_año <- grillas_año[grillas_procesar]
  }
  
  n_grillas_procesar <- length(grillas_año)
  
  if(config$verbose) {
    cat("Grillas disponibles:", n_grillas_total, "\n")
    cat("Grillas a procesar:", n_grillas_procesar, "\n")
  }
  
  resultados_grillas <- data.frame()
  
  for(j in 1:n_grillas_procesar) {
    
    grilla_name <- names(grillas_año)[j]
    grilla_data <- grillas_año[[j]]
    
    if(config$verbose && (j %% 50 == 0 || j == 1 || j == n_grillas_procesar)) {
      cat("  Procesando grilla", j, "/", n_grillas_procesar, "(", grilla_name, ")\n")
    }
    
    # Manejar datos faltantes (completar si es necesario)
    if(nrow(grilla_data) < 365) {
      if(config$verbose) cat("    Completando datos faltantes (", nrow(grilla_data), " días disponibles)\n")
      
      # Crear secuencia completa de fechas del año
      año_consumo <- año_actual - 1
      fechas_completas <- seq(as.Date(paste0(año_consumo, "-01-01")), 
                              as.Date(paste0(año_consumo, "-12-31")), 
                              by = "day")
      
      # Identificar fechas faltantes
      fechas_faltantes <- fechas_completas[!fechas_completas %in% grilla_data$fecha]
      
      if(length(fechas_faltantes) > 0) {
        # Calcular promedios por día del año de los datos existentes
        grilla_data$dia_del_año <- as.numeric(format(grilla_data$fecha, "%j"))
        promedios_por_dia <- aggregate(temperatura ~ dia_del_año, grilla_data, mean, na.rm = TRUE)
        
        # Crear filas para fechas faltantes
        for(fecha_faltante in fechas_faltantes) {
          dia_del_año_faltante <- as.numeric(format(as.Date(fecha_faltante), "%j"))
          temp_promedio <- promedios_por_dia$temperatura[promedios_por_dia$dia_del_año == dia_del_año_faltante]
          
          # Si no hay promedio para ese día, usar promedio general
          if(length(temp_promedio) == 0) {
            temp_promedio <- mean(grilla_data$temperatura, na.rm = TRUE)
          }
          
          nueva_fila <- data.frame(
            grilla_id = grilla_data$grilla_id[1],
            longitud = grilla_data$longitud[1], 
            latitud = grilla_data$latitud[1],
            fecha = as.Date(fecha_faltante),
            temperatura = temp_promedio,
            dia_del_año = dia_del_año_faltante,
            stringsAsFactors = FALSE
          )
          
          grilla_data <- rbind(grilla_data, nueva_fila)
        }
        
        # Reordenar por fecha
        grilla_data <- grilla_data[order(grilla_data$fecha), ]
        
        if(config$verbose) cat("    Agregadas", length(fechas_faltantes), "fechas faltantes\n")
      }
    }
    
    # Determinar duración real del año (365 o 366 días para años bisiestos)
    duracion_real <- nrow(grilla_data)
    
    if(config$verbose && duracion_real == 366) {
      cat("    Año bisiesto detectado (366 días)\n")
    }
    
    # Preparar temperatura para esta grilla
    temperatura_grilla <- data.frame(
      Day = 1:duracion_real,
      Temperature = grilla_data$temperatura
    )
    
    # Verificar que no haya NAs críticos
    if(any(is.na(temperatura_grilla$Temperature))) {
      if(config$verbose) cat("    SALTANDO: temperaturas NA\n")
      next
    }
    
    # Crear copia del objeto bioenergético y reemplazar temperatura
    bio_obj_grilla <- bio_obj_base
    bio_obj_grilla$environmental_data$temperature <- temperatura_grilla
    
    # Si es año bisiesto, expandir también los datos de dieta
    if(duracion_real == 366) {
      # Expandir proporciones de dieta (duplicar el día 59 = 28 feb)
      diet_props <- bio_obj_grilla$diet_data$proportions
      fila_extra <- diet_props[59, ]  # 28 de febrero
      fila_extra$Day <- 60            # Será el 29 de febrero
      bio_obj_grilla$diet_data$proportions <- rbind(
        diet_props[1:59, ],
        fila_extra,
        transform(diet_props[60:365, ], Day = Day + 1)
      )
      
      # Expandir energías de presas
      diet_energies <- bio_obj_grilla$diet_data$energies
      fila_extra_energy <- diet_energies[59, ]
      fila_extra_energy$Day <- 60
      bio_obj_grilla$diet_data$energies <- rbind(
        diet_energies[1:59, ],
        fila_extra_energy,
        transform(diet_energies[60:365, ], Day = Day + 1)
      )
      
      # Expandir datos indigestibles
      diet_indigest <- bio_obj_grilla$diet_data$indigestible
      fila_extra_indigest <- diet_indigest[59, ]
      fila_extra_indigest$Day <- 60
      bio_obj_grilla$diet_data$indigestible <- rbind(
        diet_indigest[1:59, ],
        fila_extra_indigest,
        transform(diet_indigest[60:365, ], Day = Day + 1)
      )
      
      if(config$verbose) cat("    Datos de dieta expandidos a 366 días\n")
    }
    
    tryCatch({
      
      # Ejecutar modelo con p-value fijo y duración adaptativa
      modelo_grilla <- run_fb4(
        x = bio_obj_grilla,
        fit_to = "p_value",
        fit_value = p_value_actual,
        strategy = "direct_p_value",
        first_day = 1,
        last_day = duracion_real,  # 365 para año normal, 366 para bisiesto
        verbose = FALSE
      )
      
      # Consumo total anual (sumar todos los días)
      consumo_total_g <- modelo_grilla$summary$total_consumption_g
      consumo_total_kg <- consumo_total_g / 1000
      
      # Obtener proporciones de dieta de este modelo específico
      diet_props <- bio_obj_grilla$diet_data$proportions
      prop_anchoveta <- mean(diet_props$anchoveta, na.rm = TRUE)
      prop_sardina <- mean(diet_props$sardina, na.rm = TRUE)
      prop_otros <- mean(diet_props$otros, na.rm = TRUE)
      
      # Calcular consumo por presas
      consumo_anchoveta_kg <- consumo_total_kg * prop_anchoveta
      consumo_sardina_kg <- consumo_total_kg * prop_sardina
      consumo_otros_kg <- consumo_total_kg * prop_otros
      
      # Extraer información de la grilla
      grilla_id <- grilla_data$grilla_id[1]
      longitud <- grilla_data$longitud[1]
      latitud <- grilla_data$latitud[1]
      
      resultado_grilla <- data.frame(
        grilla_id = grilla_id,
        longitud = longitud,
        latitud = latitud,
        consumo_total_kg = consumo_total_kg,
        consumo_anchoveta_kg = consumo_anchoveta_kg,
        consumo_sardina_kg = consumo_sardina_kg,
        consumo_otros_kg = consumo_otros_kg,
        prop_anchoveta = prop_anchoveta,
        prop_sardina = prop_sardina,
        prop_otros = prop_otros,
        p_value_usado = p_value_actual,
        stringsAsFactors = FALSE
      )
      
      resultados_grillas <- rbind(resultados_grillas, resultado_grilla)
      
    }, error = function(e) {
      if(config$verbose) cat("    ERROR en", grilla_name, ":", e$message, "\n")
    })
  }
  
  año_key <- paste0("año_", año_actual)
  edad_key <- paste0("edad_", edad_actual)
  
  resultado_final[[año_key]][[edad_key]] <- resultados_grillas
  
  contador_completados <- contador_completados + 1
  
  if(config$verbose) {
    cat("Completado:", nrow(resultados_grillas), "grillas exitosas\n")
    if(nrow(resultados_grillas) > 0) {
      cat("Consumo promedio:", round(mean(resultados_grillas$consumo_total_kg), 2), "kg\n")
    }
  }
}

# --- guardar resultado final ---

dir.create("data/data_raw/spatial-consumption", showWarnings = FALSE, recursive = TRUE)

saveRDS(resultado_final, "data/data_raw/spatial-consumption/resultado_consumo_grillas_ISOTOPES.rds")

cat("Combinaciones procesadas:", contador_completados, "/", total_combinaciones, "\n")
cat("Errores:", contador_errores, "\n")

# --- funciones de visualización ---

plot_consumo_mapa <- function(datos, año, edad, tipo = "total", 
                             titulo = NULL, escala_colores = "viridis",
                             mostrar_stats = TRUE, verbose = TRUE) {
  
  # Construir claves para acceder a los datos
  año_key <- paste0("año_", año)
  edad_key <- paste0("edad_", edad)
  
  if(!año_key %in% names(datos)) {
    stop("Año ", año, " no encontrado en los datos. Años disponibles: ", 
         paste(gsub("año_", "", names(datos)), collapse = ", "))
  }
  
  if(!edad_key %in% names(datos[[año_key]])) {
    stop("Edad ", edad, " no encontrada para año ", año, ". Edades disponibles: ",
         paste(gsub("edad_", "", names(datos[[año_key]])), collapse = ", "))
  }
  
  datos_grillas <- datos[[año_key]][[edad_key]]
  
  if(nrow(datos_grillas) == 0) {
    stop("No hay datos disponibles para año ", año, " y edad ", edad)
  }
  
  if(verbose) {
    cat("Datos extraídos:", nrow(datos_grillas), "grillas\n")
    cat("Año:", año, "| Edad:", edad, "| Tipo:", tipo, "\n")
  }
  
  if(tipo == "presas") {
    mapas_presas <- list()
    
    for(presa in c("anchoveta", "sardina", "otros")) {
      mapa_presa <- crear_mapa_individual(datos_grillas, año, edad, presa, 
                                         titulo, escala_colores, mostrar_stats, verbose)
      mapas_presas[[presa]] <- mapa_presa
    }
    
    if(verbose) cat("Mapas de presas creados: anchoveta, sardina, otros\n")
    return(mapas_presas)
    
  } else {
    mapa_individual <- crear_mapa_individual(datos_grillas, año, edad, tipo, 
                                           titulo, escala_colores, mostrar_stats, verbose)
    return(mapa_individual)
  }
}

crear_mapa_individual <- function(datos_grillas, año, edad, tipo, 
                                 titulo = NULL, escala_colores = "viridis", 
                                 mostrar_stats = TRUE, verbose = FALSE) {
  
  config_mapa <- switch(tipo,
    "total" = list(
      var = "consumo_total_kg",
      nombre_var = "Consumo Total\n(kg/año)",
      color_option = escala_colores,
      titulo_base = "Consumo Total Anual"
    ),
    "anchoveta" = list(
      var = "consumo_anchoveta_kg", 
      nombre_var = "Consumo Anchoveta\n(kg/año)",
      color_option = "plasma",
      titulo_base = "Consumo de Anchoveta"
    ),
    "sardina" = list(
      var = "consumo_sardina_kg",
      nombre_var = "Consumo Sardina\n(kg/año)", 
      color_option = "viridis",
      titulo_base = "Consumo de Sardina"
    ),
    "otros" = list(
      var = "consumo_otros_kg",
      nombre_var = "Consumo Otras Presas\n(kg/año)",
      color_option = "magma", 
      titulo_base = "Consumo de Otras Presas"
    ),
    stop("Tipo '", tipo, "' no válido. Opciones: 'total', 'anchoveta', 'sardina', 'otros', 'presas'")
  )
  
  if(!config_mapa$var %in% colnames(datos_grillas)) {
    stop("Variable '", config_mapa$var, "' no encontrada en los datos")
  }
  
  # Cargar mapa de Chile
  chile <- ne_countries(country = "Chile", returnclass = "sf", scale = "medium")
  
  datos_mapa <- data.frame(
    longitud = datos_grillas$longitud,
    latitud = datos_grillas$latitud,
    consumo = datos_grillas[[config_mapa$var]],
    grilla_id = datos_grillas$grilla_id
  )
  
  datos_mapa <- datos_mapa[!is.na(datos_mapa$consumo), ]
  
  if(nrow(datos_mapa) == 0) {
    stop("No hay datos válidos para graficar")
  }
  
  stats <- list(
    n_grillas = nrow(datos_mapa),
    consumo_min = min(datos_mapa$consumo, na.rm = TRUE),
    consumo_max = max(datos_mapa$consumo, na.rm = TRUE),
    consumo_medio = mean(datos_mapa$consumo, na.rm = TRUE),
    consumo_total = sum(datos_mapa$consumo, na.rm = TRUE)
  )
  
  xlim <- range(datos_mapa$longitud) + c(-0.3, 0.3)
  ylim <- range(datos_mapa$latitud) + c(-0.3, 0.3)
  
  año_consumo <- año - 1

  if(is.null(titulo)) {
      titulo <- paste(config_mapa$titulo_base, "- Chinook Edad", edad, "(ISOTOPES)")
  }

  if(mostrar_stats) {
      subtitle <- paste0(
          "Año ", año_consumo, " | ", stats$n_grillas, " grillas | ",
          "Promedio: ", round(stats$consumo_medio, 1), " kg | ",
          "Total: ", round(stats$consumo_total, 0), " kg"
      )
  } else {
      subtitle <- paste("Año", año_consumo, "|", stats$n_grillas, 
                        "grillas | Dieta edades 1-2: muscle isotope model, 5 sources")
  }

  escala_color <- switch(config_mapa$color_option,
    "viridis" = function(...) scale_fill_viridis_c(direction = -1, ...),
    "plasma" = function(...) scale_fill_viridis_c(option = "plasma", direction = -1, ...),
    "magma" = function(...) scale_fill_viridis_c(option = "magma", direction = -1, ...),
    "inferno" = function(...) scale_fill_viridis_c(option = "inferno", direction = -1, ...),
    scale_fill_viridis_c  # default
  )
  
  mapa <- ggplot() +
    
    # Mapa de Chile como fondo
    geom_sf(data = chile, fill = "grey95", color = "grey80", size = 0.3) +
    
    # Consumo como grillas (tiles)
    geom_tile(data = datos_mapa, 
              aes(x = longitud, y = latitud, fill = consumo),
              alpha = 0.9) +
    
    escala_color(
      name = config_mapa$nombre_var,
      labels = function(x) format(round(x, 1), nsmall = 1),
      na.value = "transparent",
      trans = if(stats$consumo_max > stats$consumo_min * 10) "sqrt" else "identity"
    ) +
    
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    
    labs(
      title = titulo,
      subtitle = subtitle,
      x = "Longitud",
      y = "Latitud",
      caption = "Dieta edades 1-2: muscle isotope model, 5 sources | Edades 3-5: estómagos"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      plot.caption = element_text(hjust = 1, size = 8, color = "grey50"),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, size = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.background = element_rect(fill = "white")
    )
  
  attr(mapa, "estadisticas") <- stats
  
  if(verbose) {
    cat("Mapa creado - Rango:", round(stats$consumo_min, 2), "a", 
        round(stats$consumo_max, 2), "kg\n")
  }
  
  return(mapa)
}

# --- comparar edades ---

plot_comparar_edades <- function(datos, año, edades = 1:4, tipo = "total",
                                ncol = 2, titulo_general = NULL) {
  
  if(tipo == "presas") {
    mapas_por_presa <- list()
    
    for(presa in c("anchoveta", "sardina", "otros")) {
      mapas_edad_presa <- list()
      
      for(edad in edades) {
        año_key <- paste0("año_", año)
        edad_key <- paste0("edad_", edad)
        
        if(año_key %in% names(datos) && edad_key %in% names(datos[[año_key]])) {
          
          mapa <- plot_consumo_mapa(datos, año, edad, presa, 
                                   titulo = paste("Edad", edad),
                                   mostrar_stats = FALSE, verbose = FALSE)
          
          mapas_edad_presa[[paste0("edad_", edad)]] <- mapa
          
        }
      }
      
      if(length(mapas_edad_presa) > 0) {
        if(is.null(titulo_general)) {
          titulo_presa <- paste("Consumo de", tools::toTitleCase(presa), "- Año", año - 1, "(ISOTOPES)")
        } else {
          titulo_presa <- paste(titulo_general, "-", tools::toTitleCase(presa))
        }
        
        combined_presa <- wrap_plots(mapas_edad_presa, ncol = ncol) +
          plot_annotation(
            title = titulo_presa,
            caption = "Dieta edades 1-2: muscle isotope model, 5 sources | Edades 3-5: estómagos",
            theme = theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.caption = element_text(size = 9, hjust = 1, color = "grey50")
            )
          )
        
        mapas_por_presa[[presa]] <- combined_presa
        cat("  ", tools::toTitleCase(presa), ":", length(mapas_edad_presa), "edades\n")
      }
    }
    
    cat("Mapas por presas completados:", length(mapas_por_presa), "presas\n")
    return(mapas_por_presa)
    
  } else {
    mapas_edades <- list()
    
    for(edad in edades) {
      
      año_key <- paste0("año_", año)
      edad_key <- paste0("edad_", edad)
      
      if(año_key %in% names(datos) && edad_key %in% names(datos[[año_key]])) {
        
        mapa <- plot_consumo_mapa(datos, año, edad, tipo, 
                                 titulo = paste("Edad", edad),
                                 mostrar_stats = FALSE, verbose = FALSE)
        
        mapas_edades[[paste0("edad_", edad)]] <- mapa
        cat("  Edad", edad, "\n")
        
      } else {
        cat("  Edad", edad, ": Sin datos\n")
      }
    }
    
    if(length(mapas_edades) == 0) {
      stop("No hay datos disponibles para ninguna edad")
    }
    
    if(is.null(titulo_general)) {
      titulo_general <- paste("Comparación de Consumo por Edades - Año", año - 1, "(ISOTOPES)")
    }
    
    combined_plot <- wrap_plots(mapas_edades, ncol = ncol) +
      plot_annotation(
        title = titulo_general,
        caption = "Dieta edades 1-2: muscle isotope model, 5 sources | Edades 3-5: estómagos",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.caption = element_text(size = 9, hjust = 1, color = "grey50")
        )
      )
    
    cat("Mapas combinados:", length(mapas_edades), "edades\n")
    
    return(combined_plot)
  }
}

# --- extraer estadísticas ---

extraer_estadisticas_consumo <- function(datos, año, edad, tipo = "total") {
  
  año_key <- paste0("año_", año)
  edad_key <- paste0("edad_", edad)
  
  if(!año_key %in% names(datos) || !edad_key %in% names(datos[[año_key]])) {
    return(NULL)
  }
  
  datos_grillas <- datos[[año_key]][[edad_key]]
  
  var_consumo <- switch(tipo,
    "total" = "consumo_total_kg",
    "anchoveta" = "consumo_anchoveta_kg", 
    "sardina" = "consumo_sardina_kg",
    "otros" = "consumo_otros_kg",
    stop("Tipo no válido")
  )
  
  if(!var_consumo %in% colnames(datos_grillas)) {
    return(NULL)
  }
  
  consumo_values <- datos_grillas[[var_consumo]][!is.na(datos_grillas[[var_consumo]])]
  
  if(length(consumo_values) == 0) {
    return(NULL)
  }
  
  return(list(
    año = año,
    edad = edad,
    tipo = tipo,
    n_grillas = length(consumo_values),
    min = min(consumo_values),
    max = max(consumo_values),
    media = mean(consumo_values),
    mediana = median(consumo_values),
    total = sum(consumo_values),
    sd = sd(consumo_values)
  ))
}

# ---- Ejemplos de uso ----
# resultado_final <- readRDS("data/data_raw/spatial-consumption/resultado_consumo_grillas_ISOTOPES.rds")
# plot_consumo_mapa(resultado_final, año = 2015, edad = 1, tipo = "total")
# plot_comparar_edades(resultado_final, 2015, edad = 1:4, tipo = "total")
# plot_comparar_edades(resultado_final, 2015, edad = 1:4, tipo = "presas")
