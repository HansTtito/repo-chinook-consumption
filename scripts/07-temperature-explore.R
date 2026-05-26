# FUNCIONES PARA PROCESAMIENTO DE DATOS DE TEMPERATURA NETCDF
# Autor: Procesamiento de datos satelitales de temperatura marina
# Objetivo: Extraer series temporales de temperatura por grilla para modelos

library(terra)
library(sf)
library(reshape2)
library(dplyr)
library(ncdf4)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# =============================================================================
# FUNCIÓN 1: CARGAR Y ANALIZAR ARCHIVO NETCDF
# =============================================================================

# Helper: convierte tiempo CF-compliant a POSIXct
convert_cf_time <- function(vals, units, calendar = "gregorian", tz = "UTC") {
  # Si no hay valores o units, retornar vacío
  if (length(vals) == 0 || is.null(units)) return(as.POSIXct(character(0), tz = tz))
  
  # Extraer la unidad base y el origen
  # Ej: "hours since 1950-01-01 00:00:00"
  unit_base <- tolower(sub("\\s+since.*$", "", units))
  origin_str <- sub("^.*since\\s*", "", units, ignore.case = TRUE)
  origin_str <- gsub("T", " ", origin_str)  # algunos NetCDF usan T
  
  # Asegurar que tenga hora
  if (!grepl("\\d{2}:\\d{2}", origin_str)) origin_str <- paste0(origin_str, " 00:00:00")
  
  origin <- as.POSIXct(origin_str, tz = tz)
  if (is.na(origin)) stop("No se pudo parsear el origen temporal desde: ", units)
  
  # Convertir vals a segundos según la unidad
  mult <- switch(unit_base,
                 "seconds" = 1,
                 "second"  = 1,
                 "minutes" = 60,
                 "minute"  = 60,
                 "hours"   = 3600,
                 "hour"    = 3600,
                 "days"    = 86400,
                 "day"     = 86400,
                 stop("Unidad de tiempo no soportada: ", unit_base))
  
  # Ignorar calendarios raros por ahora (solo gregorian, standard, proleptic)
  if (!tolower(calendar) %in% c("gregorian", "standard", "proleptic_gregorian", "", NA)) {
    warning("Calendario '", calendar, "' no completamente soportado. Para '360_day'/'noleap', considera PCICt::as.PCICt().")
  }
  
  # Convertir a POSIXct
  origin + vals * mult
}

analizar_estructura_netcdf <- function(archivo_path, poligono_path = NULL, variable = NULL) {
  cat("=== ANALIZANDO ESTRUCTURA DEL ARCHIVO NETCDF ===\n")
  
  # Abrir con ncdf4 para ver dimensiones
  nc <- nc_open(archivo_path)
  on.exit(nc_close(nc)) # cerrar cuando acabe
  
  # Información sobre variables disponibles
  variables <- names(nc$var)
  cat("Variables disponibles:", paste(variables, collapse = ", "), "\n")
  
  dims <- lapply(nc$dim, function(x) {
    list(
      name = x$name,
      size = x$len,
      units = x$units,
      vals = if (x$len <= 20) x$vals else paste0("(", x$len, " valores)")
    )
  })
  
  # Extraer profundidad si existe
  depth_vals <- if ("depth" %in% names(nc$dim)) nc$dim$depth$vals else NA
  
  # Tiempo bruto (valores + unidades + calendario)
  if ("time" %in% names(nc$dim)) {
    time_raw <- if (!is.null(nc$dim$time$vals)) nc$dim$time$vals else ncvar_get(nc, "time")
    time_units <- nc$dim$time$units %||% ncatt_get(nc, "time", "units")$value
    time_cal   <- nc$dim$time$calendar %||% ncatt_get(nc, "time", "calendar")$value %||% "gregorian"
    time_posix <- convert_cf_time(time_raw, units = time_units, calendar = time_cal, tz = "UTC")
  } else {
    time_raw <- numeric(0)
    time_units <- NA_character_
    time_cal <- NA_character_
    time_posix <- as.POSIXct(character(0), tz="UTC")
  }
  
  # Si se especifica una variable, extraer solo esa
  if (!is.null(variable)) {
    if (variable %in% variables) {
      r <- rast(archivo_path, subds = variable)
      cat("Variable seleccionada:", variable, "\n")
    } else {
      stop("Variable '", variable, "' no encontrada. Disponibles: ", paste(variables, collapse = ", "))
    }
  } else {
    r <- tryCatch(rast(archivo_path), error = function(e) NULL)
  }
  
  # Extraer información temporal y de profundidad de terra si existe
  if (!is.null(r)) {
    d <- tryCatch(depth(r), error = function(e) rep(NA, nlyr(r)))
    t <- tryCatch(time(r), error = function(e) rep(as.Date(NA), nlyr(r)))
  } else {
    d <- numeric(0)
    t <- as.Date(character(0))
  }
  
  # Análisis detallado como en la función original
  profundidades_unicas <- if (!all(is.na(depth_vals))) depth_vals else if (!all(is.na(d))) unique(d) else numeric(0)
  fechas_unicas <- if (length(time_posix) > 0) unique(as.Date(time_posix)) else if (length(t) > 0 && !all(is.na(t))) unique(t) else as.Date(character(0))
  
  n_profundidades <- length(profundidades_unicas)
  n_dias <- length(fechas_unicas)
  
  estructura <- list(
    # Información básica del archivo
    archivo = basename(archivo_path),
    dimensiones = if (!is.null(r)) c(filas = nrow(r), columnas = ncol(r)) else c(filas = NA, columnas = NA),
    total_capas = if (!is.null(r)) nlyr(r) else NA,
    n_layers = if (!is.null(r)) nlyr(r) else NA, # alias para compatibilidad
    
    # Información espacial
    resolucion = if (!is.null(r)) res(r) else c(NA, NA),
    extension = if (!is.null(r)) as.vector(ext(r)) else c(NA, NA, NA, NA),
    crs = if (!is.null(r)) crs(r) else NA,
    
    # Información de profundidad
    profundidades_unicas = profundidades_unicas,
    depth_vals = depth_vals,
    n_profundidades = n_profundidades,
    rango_profundidades = if (length(profundidades_unicas) > 0) range(profundidades_unicas) else c(NA, NA),
    
    # Información temporal detallada
    fechas_unicas = fechas_unicas,
    n_dias = n_dias,
    rango_fechas = if (length(fechas_unicas) > 0) range(fechas_unicas) else as.Date(c(NA, NA)),
    time_raw = time_raw,
    time_units = time_units,
    time_calendar = time_cal,
    time = time_posix,
    
    # Información de dimensiones NetCDF y variables
    dimensiones_netcdf = dims,
    variables_disponibles = variables,
    variable_seleccionada = variable,
    
    # Objeto raster original
    raster_original = r
  )
  
  # Mostrar información detallada
  cat("Archivo:", estructura$archivo, "\n")
  if (!is.null(r)) {
    cat("Dimensiones:", estructura$dimensiones[1], "x", estructura$dimensiones[2], "\n")
    cat("Total de capas:", estructura$total_capas, "\n")
    cat("Resolución:", estructura$resolucion, "\n")
    cat("Extensión:", estructura$extension, "\n")
  }
  
  if (n_profundidades > 0) {
    cat("Profundidades únicas:", estructura$n_profundidades, "\n")
    cat("Rango de profundidades:", round(estructura$rango_profundidades, 2), "m\n")
    cat("Profundidades disponibles:", paste(round(profundidades_unicas, 2), collapse = ", "), "m\n")
  }
  
  if (length(time_posix) > 0) {
    cat("Información temporal:\n")
    cat("  - Unidades CF:", time_units, "\n")
    cat("  - Calendario:", time_cal, "\n")
    cat("  - Total de tiempos:", length(time_posix), "\n")
    cat("  - Días únicos:", estructura$n_dias, "\n")
    cat("  - Rango de fechas:", as.character(estructura$rango_fechas), "\n")
    cat("  - Primeras fechas:", paste(head(time_posix, 3), collapse = ", "), "\n")
    if (length(time_posix) > 3) {
      cat("  - Últimas fechas:", paste(tail(time_posix, 3), collapse = ", "), "\n")
    }
  }
  
  # Cargar polígono si se proporciona
  if (!is.null(poligono_path)) {
    poligono <- st_read(poligono_path, quiet = TRUE)
    estructura$poligono <- vect(poligono)
    cat("Polígono cargado exitosamente\n")
  }
  
  return(estructura)
}

# =============================================================================
# FUNCIÓN 2: FILTRAR Y PROCESAR PROFUNDIDADES
# =============================================================================

procesar_profundidades <- function(estructura, prof_min = 20, prof_max = 40, 
                                   recortar_espacial = TRUE) {
  "
  Procesa las profundidades objetivo y calcula promedios diarios
  
  Args:
    estructura: Output de analizar_estructura_netcdf()
    prof_min: Profundidad mínima (metros)
    prof_max: Profundidad máxima (metros)
    recortar_espacial: Si recortar al polígono
    
  Returns:
    Raster stack con promedio diario de profundidades objetivo
  "
  
  cat("=== PROCESANDO PROFUNDIDADES", prof_min, "-", prof_max, "m ===\n")
  
  r <- estructura$raster_original
  profundidades_unicas <- estructura$profundidades_unicas
  n_profundidades <- estructura$n_profundidades
  n_dias <- estructura$n_dias
  
  # Filtrar profundidades objetivo
  indices_prof_objetivo <- which(profundidades_unicas >= prof_min & 
                                   profundidades_unicas <= prof_max)
  valores_prof_objetivo <- profundidades_unicas[indices_prof_objetivo]
  
  if (length(valores_prof_objetivo) == 0) {
    stop("No se encontraron profundidades en el rango especificado")
  }
  
  cat("Profundidades encontradas:", length(valores_prof_objetivo), "\n")
  cat("Valores:", round(valores_prof_objetivo, 2), "m\n")
  
  # Recorte espacial si se requiere
  if (recortar_espacial && !is.null(estructura$poligono)) {
    cat("Recortando espacialmente...\n")
    r <- crop(r, estructura$poligono)
    r <- mask(r, estructura$poligono)
    cat("Recorte completado\n")
  }
  
  # Calcular índices de capas por día
  indices_por_dia <- list()
  for (dia in 1:n_dias) {
    indices_dia <- sapply(indices_prof_objetivo, function(prof_idx) {
      (dia - 1) * n_profundidades + prof_idx
    })
    indices_por_dia[[dia]] <- indices_dia
  }
  
  # Crear rasters promedio por día
  cat("Calculando promedios diarios...\n")
  rasters_diarios <- list()
  
  for (dia in 1:n_dias) {
    if (dia %% 50 == 0 || dia == 1) {
      cat("Procesando día", dia, "de", n_dias, "\n")
    }
    
    indices <- indices_por_dia[[dia]]
    capas_dia <- r[[indices]]
    promedio_dia <- mean(capas_dia, na.rm = TRUE)
    rasters_diarios[[dia]] <- promedio_dia
  }
  
  # Combinar en stack
  r_avg <- rast(rasters_diarios)
  fechas_ordenadas <- sort(estructura$fechas_unicas)
  names(r_avg) <- format(fechas_ordenadas, "%Y-%m-%d")
  
  cat("Procesamiento completado:", nlyr(r_avg), "días\n")
  
  return(r_avg)
}

# =============================================================================
# FUNCIÓN 3: EXTRAER SERIES TEMPORALES POR GRILLA
# =============================================================================

extraer_series_temporales <- function(raster_stack, formato = "largo") {
  "
  Extrae series temporales de temperatura para cada grilla
  
  Args:
    raster_stack: Raster stack con temperatura diaria
    formato: 'largo' o 'ancho'
    
  Returns:
    DataFrame con series temporales por grilla
  "
  
  cat("=== EXTRAYENDO SERIES TEMPORALES ===\n")
  
  # Obtener coordenadas de grillas
  coords <- xyFromCell(raster_stack, 1:ncell(raster_stack))
  
  # Extraer valores de temperatura
  temperaturas_matrix <- values(raster_stack)
  
  # Identificar grillas válidas (sin NA en el primer día)
  grillas_validas <- which(!is.na(temperaturas_matrix[, 1]))
  
  cat("Total de grillas:", nrow(coords), "\n")
  cat("Grillas válidas:", length(grillas_validas), "\n")
  
  # Filtrar datos válidos con IDs continuos
  coords_validas <- data.frame(
    grilla_id = 1:length(grillas_validas),  # IDs continuos 1,2,3...
    longitud = coords[grillas_validas, 1],
    latitud = coords[grillas_validas, 2]
  )
  
  temperaturas_validas <- temperaturas_matrix[grillas_validas, ]
  
  # Crear DataFrame según formato
  if (formato == "ancho") {
    # Formato ancho: una fila por grilla, una columna por día
    df_resultado <- data.frame(
      grilla_id = coords_validas$grilla_id,  # Ya son continuos
      longitud = coords_validas$longitud,
      latitud = coords_validas$latitud,
      temperaturas_validas
    )
    
    # Nombres de columnas con fechas
    colnames(df_resultado)[4:ncol(df_resultado)] <- names(raster_stack)
    
    cat("Formato ancho:", nrow(df_resultado), "grillas x", ncol(df_resultado), "columnas\n")
    
  } else {
    # Formato largo: una fila por grilla-día
    df_ancho <- data.frame(
      grilla_id = coords_validas$grilla_id,  # Ya son continuos
      longitud = coords_validas$longitud,
      latitud = coords_validas$latitud,
      temperaturas_validas
    )
    
    colnames(df_ancho)[4:ncol(df_ancho)] <- names(raster_stack)
    
    df_resultado <- melt(df_ancho,
                         id.vars = c("grilla_id", "longitud", "latitud"),
                         variable.name = "fecha",
                         value.name = "temperatura")
    
    # Convertir fecha
    df_resultado$fecha <- as.Date(gsub("X", "", as.character(df_resultado$fecha)))
    
    cat("Formato largo:", nrow(df_resultado), "observaciones\n")
  }
  
  return(df_resultado)
}

# =============================================================================
# FUNCIÓN 4: GENERAR SERIES INDIVIDUALES POR GRILLA
# =============================================================================

generar_series_individuales <- function(datos_temperatura, directorio_salida = NULL) {
  "
  Genera series temporales individuales para cada grilla
  
  Args:
    datos_temperatura: DataFrame en formato largo de extraer_series_temporales()
    directorio_salida: Directorio donde guardar las series (opcional)
    
  Returns:
    Lista con series temporales por grilla
  "
  
  cat("=== GENERANDO SERIES INDIVIDUALES ===\n")
  
  # Obtener IDs únicos de grillas
  grillas_unicas <- unique(datos_temperatura$grilla_id)
  n_grillas <- length(grillas_unicas)
  
  cat("Generando series para", n_grillas, "grillas\n")
  
  # Lista para almacenar series
  series_por_grilla <- list()
  
  # Generar serie para cada grilla
  for (i in 1:n_grillas) {
    grilla_id <- grillas_unicas[i]
    
    if (i %% 100 == 0 || i == 1) {
      cat("Procesando grilla", i, "de", n_grillas, "(ID:", grilla_id, ")\n")
    }
    
    # Filtrar datos de esta grilla
    datos_grilla <- datos_temperatura[datos_temperatura$grilla_id == grilla_id, ]
    
    # Ordenar por fecha
    datos_grilla <- datos_grilla[order(datos_grilla$fecha), ]
    
    # Crear serie temporal
    serie <- data.frame(
      grilla_id = grilla_id,
      longitud = datos_grilla$longitud[1],
      latitud = datos_grilla$latitud[1],
      fecha = datos_grilla$fecha,
      temperatura = datos_grilla$temperatura,
      dia_del_año = as.numeric(format(datos_grilla$fecha, "%j"))
    )
    
    # Agregar a lista
    series_por_grilla[[paste0("grilla_", grilla_id)]] <- serie
    
    # Guardar individualmente si se especifica directorio
    if (!is.null(directorio_salida)) {
      if (!dir.exists(directorio_salida)) {
        dir.create(directorio_salida, recursive = TRUE)
      }
      
      archivo_salida <- file.path(directorio_salida, paste0("serie_grilla_", grilla_id, ".csv"))
      write.csv(serie, archivo_salida, row.names = FALSE)
    }
  }
  
  cat("Series generadas exitosamente\n")
  
  return(series_por_grilla)
}

# =============================================================================
# FUNCIÓN 5: PREPARAR DATOS PARA MODELADO
# =============================================================================

preparar_datos_modelo <- function(series_por_grilla, incluir_coordenadas = TRUE) {
  "
  Prepara datos en formato matriz para modelado
  
  Args:
    series_por_grilla: Output de generar_series_individuales()
    incluir_coordenadas: Si incluir coordenadas como features
    
  Returns:
    Lista con matrices de datos y metadatos
  "
  
  cat("=== PREPARANDO DATOS PARA MODELO ===\n")
  
  n_grillas <- length(series_por_grilla)
  
  # Obtener longitud de serie (asumiendo todas iguales)
  longitud_serie <- nrow(series_por_grilla[[1]])
  
  cat("Número de grillas:", n_grillas, "\n")
  cat("Longitud de serie temporal:", longitud_serie, "días\n")
  
  # Crear matriz de temperaturas (grillas x días)
  matriz_temperaturas <- matrix(NA, nrow = n_grillas, ncol = longitud_serie)
  
  # Metadatos de grillas
  metadatos_grillas <- data.frame(
    indice = 1:n_grillas,
    grilla_id = numeric(n_grillas),
    longitud = numeric(n_grillas),
    latitud = numeric(n_grillas)
  )
  
  # Llenar matriz y metadatos
  for (i in 1:n_grillas) {
    serie <- series_por_grilla[[i]]
    
    # Temperaturas
    matriz_temperaturas[i, ] <- serie$temperatura
    
    # Metadatos
    metadatos_grillas$grilla_id[i] <- serie$grilla_id[1]
    metadatos_grillas$longitud[i] <- serie$longitud[1]
    metadatos_grillas$latitud[i] <- serie$latitud[1]
  }
  
  # Crear nombres para filas y columnas
  rownames(matriz_temperaturas) <- paste0("grilla_", metadatos_grillas$grilla_id)
  colnames(matriz_temperaturas) <- paste0("dia_", 1:longitud_serie)
  
  # Fechas
  fechas <- series_por_grilla[[1]]$fecha
  
  # Preparar resultado
  resultado <- list(
    X = matriz_temperaturas,  # Matriz principal (grillas x días)
    metadatos = metadatos_grillas,
    fechas = fechas,
    n_grillas = n_grillas,
    n_dias = longitud_serie,
    estadisticas = list(
      temp_min = min(matriz_temperaturas, na.rm = TRUE),
      temp_max = max(matriz_temperaturas, na.rm = TRUE),
      temp_media = mean(matriz_temperaturas, na.rm = TRUE),
      temp_sd = sd(matriz_temperaturas, na.rm = TRUE)
    )
  )
  
  # Agregar coordenadas como features si se requiere
  if (incluir_coordenadas) {
    matriz_coords <- cbind(metadatos_grillas$longitud, metadatos_grillas$latitud)
    colnames(matriz_coords) <- c("longitud", "latitud")
    rownames(matriz_coords) <- rownames(matriz_temperaturas)
    resultado$coords = matriz_coords
  }
  
  cat("Datos preparados:\n")
  cat("- Matriz X:", nrow(resultado$X), "x", ncol(resultado$X), "\n")
  cat("- Rango temperaturas:", round(resultado$estadisticas$temp_min, 2), "a", 
      round(resultado$estadisticas$temp_max, 2), "°C\n")
  cat("- Temperatura media:", round(resultado$estadisticas$temp_media, 2), "°C\n")
  
  return(resultado)
}


# =============================================================================
# FUNCIÓN 7: UTILIDADES ADICIONALES
# =============================================================================

# Función para visualizar una grilla específica
visualizar_serie_grilla <- function(series_por_grilla, grilla_id, titulo = NULL) {
  "Visualiza la serie temporal de una grilla específica"
  
  nombre_grilla <- paste0("grilla_", grilla_id)
  
  if (nombre_grilla %in% names(series_por_grilla)) {
    serie <- series_por_grilla[[nombre_grilla]]
    
    if (is.null(titulo)) {
      titulo <- paste("Serie Temporal - Grilla", grilla_id)
    }
    
    plot(serie$fecha, serie$temperatura, type = "l", 
         main = titulo,
         xlab = "Fecha", ylab = "Temperatura (°C)",
         col = "blue", lwd = 2)
    
    grid(TRUE)
    
    # Agregar información
    subtitle <- paste("Coords:", round(serie$longitud[1], 3), ",", round(serie$latitud[1], 3))
    mtext(subtitle, side = 3, line = 0.5, cex = 0.8)
    
  } else {
    cat("Grilla", grilla_id, "no encontrada\n")
  }
}


# =============================================================================
# FUNCIÓN PARA OBTENER PROMEDIO DIARIO TOTAL
# =============================================================================

obtener_promedio_diario_total <- function(datos_modelo) {
  "
 Calcula el promedio diario de temperatura de todas las grillas
 
 Args:
   datos_modelo: Resultado de preparar_datos_modelo()
   
 Returns:
   DataFrame con dos columnas: Day y Temperature
 "
  
  # Calcular promedio por día (columna) de todas las grillas (filas)
  promedios_diarios <- colMeans(datos_modelo$X, na.rm = TRUE)
  
  # Crear DataFrame resultado
  serie_promedio <- data.frame(
    Day = 1:length(promedios_diarios),
    Temperature = promedios_diarios
  )
  
  # Agregar fechas si están disponibles
  if (!is.null(datos_modelo$fechas)) {
    serie_promedio$Date <- datos_modelo$fechas
  }
  
  cat("Serie promedio diaria creada:", nrow(serie_promedio), "días\n")
  cat("Rango de temperatura:", round(range(serie_promedio$Temperature), 2), "°C\n")
  
  return(serie_promedio)
}


# =============================================================================
# FUNCIÓN PARA CREAR MAPA DE TEMPERATURA
# =============================================================================

mapear_temperatura_chile <- function(datos_modelo, fecha_index = 1, titulo = NULL) {
  "
  Crea mapa de temperatura en grillas sobre mapa de Chile
  
  Args:
    datos_modelo: Resultado de preparar_datos_modelo()
    fecha_index: Índice del día a visualizar (1-365)
    titulo: Título personalizado (opcional)
    
  Returns:
    Gráfico ggplot
  "
  
  # Validar entrada
  if (fecha_index < 1 || fecha_index > ncol(datos_modelo$X)) {
    stop("fecha_index debe estar entre 1 y ", ncol(datos_modelo$X))
  }
  
  # Cargar mapa de Chile
  chile <- ne_countries(country = "Chile", returnclass = "sf", scale = "medium")
  
  # Preparar datos de temperatura para el día específico
  datos_temp <- data.frame(
    longitud = datos_modelo$metadatos$longitud,
    latitud = datos_modelo$metadatos$latitud,
    temperatura = datos_modelo$X[, fecha_index]
  )
  
  # Remover NAs
  datos_temp <- datos_temp[!is.na(datos_temp$temperatura), ]
  
  # Calcular límites espaciales
  xlim <- range(datos_temp$longitud) + c(-0.3, 0.3)
  ylim <- range(datos_temp$latitud) + c(-0.3, 0.3)
  
  # Calcular tamaño de grilla
  lon_unique <- sort(unique(datos_temp$longitud))
  lat_unique <- sort(unique(datos_temp$latitud))
  
  if (length(lon_unique) > 1 && length(lat_unique) > 1) {
    lon_spacing <- mean(diff(lon_unique), na.rm = TRUE)
    lat_spacing <- mean(diff(lat_unique), na.rm = TRUE)
  } else {
    lon_spacing <- 0.05
    lat_spacing <- 0.05
  }
  
  # Crear título si no se proporciona
  if (is.null(titulo)) {
    if (!is.null(datos_modelo$fechas)) {
      fecha_str <- format(datos_modelo$fechas[fecha_index], "%Y-%m-%d")
      titulo <- paste("Temperatura Superficial del Mar -", fecha_str)
    } else {
      titulo <- paste("Temperatura Superficial del Mar - Día", fecha_index)
    }
  }
  
  # Crear mapa con grillas
  mapa <- ggplot() +
    # Mapa de Chile como fondo
    geom_sf(data = chile, fill = "grey95", color = "black", size = 0.3) +
    # Temperatura como grillas (tiles)
    geom_tile(data = datos_temp, 
              aes(x = longitud, y = latitud, fill = temperatura),
              # width = lon_spacing * 0.95, 
              # height = lat_spacing * 0.95,
              alpha = 0.9) +
    # Escala de colores para temperatura
    scale_fill_viridis_c(
      name = "Temperatura\n(°C)",
      option = "plasma",
      direction = -1,  # Más oscuro = más caliente
      labels = function(x) format(round(x, 1), nsmall = 1),
      na.value = "transparent"
    ) +
    # Límites espaciales
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    # Etiquetas
    labs(
      title = "",
      subtitle = paste("Área de estudio: Costa de Chile |", nrow(datos_temp), "grillas"),
      x = "Longitud",
      y = "Latitud"
    ) +
    # Tema
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "grey40"),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm"),
      # panel.grid.major = element_line(color = "grey98", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, size = 0.5),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      panel.background = element_rect(fill = "white")
    )
  
  return(mapa)
}

# =============================================================================
# FUNCIÓN PARA CREAR MAPA DE TEMPERATURA
# =============================================================================

mapear_temperatura_promedio_chile <- function(datos_modelo, titulo = NULL) {
  "
  Crea mapa de temperatura promedio anual en grillas sobre mapa de Chile
  
  Args:
    datos_modelo: Resultado de preparar_datos_modelo()
    titulo: Título personalizado (opcional)
    
  Returns:
    Gráfico ggplot
  "
  
  # Cargar mapa de Chile
  chile <- ne_countries(country = "Chile", returnclass = "sf", scale = "medium")
  
  # Calcular temperatura promedio anual para cada grilla
  temp_promedio <- rowMeans(datos_modelo$X, na.rm = TRUE)
  
  # Preparar datos de temperatura promedio
  datos_temp <- data.frame(
    longitud = datos_modelo$metadatos$longitud,
    latitud = datos_modelo$metadatos$latitud,
    temperatura_promedio = temp_promedio
  )
  
  # Remover NAs
  datos_temp <- datos_temp[!is.na(datos_temp$temperatura_promedio), ]
  
  # Calcular límites espaciales
  xlim <- range(datos_temp$longitud) + c(-0.3, 0.3)
  ylim <- range(datos_temp$latitud) + c(-0.3, 0.3)
  
  # Calcular tamaño de grilla
  lon_unique <- sort(unique(datos_temp$longitud))
  lat_unique <- sort(unique(datos_temp$latitud))
  
  if (length(lon_unique) > 1 && length(lat_unique) > 1) {
    lon_spacing <- mean(diff(lon_unique), na.rm = TRUE)
    lat_spacing <- mean(diff(lat_unique), na.rm = TRUE)
  } else {
    lon_spacing <- 0.05
    lat_spacing <- 0.05
  }
  
  # Crear título si no se proporciona
  if (is.null(titulo)) {
    if (!is.null(datos_modelo$fechas)) {
      year <- format(datos_modelo$fechas[1], "%Y")
      titulo <- paste("Temperatura Promedio Anual del Mar -", year)
    } else {
      titulo <- "Temperatura Promedio Anual del Mar"
    }
  }
  
  # Crear mapa con grillas
  mapa <- ggplot() +
    # Mapa de Chile como fondo
    geom_sf(data = chile, fill = "grey95", color = "black", size = 0.3) +
    # Temperatura como grillas (tiles)
    geom_tile(data = datos_temp, 
              aes(x = longitud, y = latitud, fill = temperatura_promedio),
              # width = lon_spacing * 0.95, 
              # height = lat_spacing * 0.95,
              alpha = 0.9) +
    # Escala de colores para temperatura
    scale_fill_viridis_c(
      name = "Temperatura\nPromedio (°C)",
      option = "plasma",
      direction = -1,  # Más oscuro = más caliente
      labels = function(x) format(round(x, 1), nsmall = 1),
      na.value = "transparent"
    ) +
    # Límites espaciales
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    # Etiquetas
    labs(
      title = titulo,
      subtitle = paste("Área de estudio: Costa de Chile |", nrow(datos_temp), "grillas | Promedio de", ncol(datos_modelo$X), "días"),
      x = "Longitud",
      y = "Latitud"
    ) +
    # Tema
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm"),
      # panel.grid.major = element_line(color = "grey98", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, size = 0.5),
      axis.text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 22, color = "black"),
      panel.background = element_rect(fill = "white")
    )
  return(mapa)
}


archivos_temperatura <- list.files(path = "data/data_raw/temperature/nc_files/",  pattern = "mercator", full.names = TRUE)
poligono_shp <- "data/spatial/poligono_costero_60mn.shp"

library(lubridate)

procesar_todos_los_años <- function(archivos, poligono_path, prof_min = 15, prof_max = 45) {
  
  # Listas para almacenar resultados por año
  rasters_por_año <- list()
  series_por_año <- list()
  datos_modelo_por_año <- list()
  promedios_diarios_por_año <- list()  # NUEVO
  
  for (i in seq_along(archivos)) {
    
    año <- gsub(".*mercator_(\\d{4})\\.nc", "\\1", basename(archivos[i]))
    cat("Procesando año:", año, "(", i, "/", length(archivos), ")\n")
    
    # Procesar archivo
    estructura <- analizar_estructura_netcdf(archivos[i], 
                                            variable = "thetao",
                                            poligono_path = poligono_path)
    
    # Obtener raster procesado
    raster_año <- procesar_profundidades(estructura, prof_min, prof_max)
    
    # Extraer series temporales
    datos_largo <- extraer_series_temporales(raster_año, formato = "largo")
    series_individuales <- generar_series_individuales(datos_largo)
    datos_modelo <- preparar_datos_modelo(series_individuales)
    
    # Obtener promedio diario total
    promedio_diario <- obtener_promedio_diario_total(datos_modelo)
    promedio_diario$año <- año  # Agregar columna de año
    
    # Guardar en listas con nombre del año
    rasters_por_año[[año]] <- raster_año
    series_por_año[[año]] <- series_individuales
    datos_modelo_por_año[[año]] <- datos_modelo
    promedios_diarios_por_año[[año]] <- promedio_diario  # NUEVO
    
    cat("Año", año, "completado:", nlyr(raster_año), "días,", 
        length(series_individuales), "grillas\n\n")
  }
  
  # Resultado final
  resultado <- list(
    rasters = rasters_por_año,
    series_grillas = series_por_año,
    datos_modelo = datos_modelo_por_año,
    promedios_diarios = promedios_diarios_por_año,  # NUEVO
    años_procesados = names(rasters_por_año)
  )
  
  cat("=== PROCESAMIENTO COMPLETO ===\n")
  cat("Años procesados:", paste(resultado$años_procesados, collapse = ", "), "\n")
  
  return(resultado)
}

# Procesar todos los archivos
resultados_completos <- procesar_todos_los_años(archivos_temperatura, poligono_shp)

# Juntar todos los promedios diarios
promedio_diario_completo <- do.call(rbind, resultados_completos$promedios_diarios)

rownames(promedio_diario_completo) <- NULL

# Ver resultado
head(promedio_diario_completo)

names(promedio_diario_completo) <- c("Day","Temperature","Date","year")

write.csv(promedio_diario_completo, "data/data_raw/temperature/temperature_time_serie_daily.csv", row.names = FALSE)


saveRDS(resultados_completos, "data/data_raw/temperature/total_data_temperature.rds")



promedio_diario_completo |> 
  ggplot(aes(x = Date, y = Temperature)) +
  geom_line(color = "blue") +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Serie Temporal Diaria de Temperatura Superficial del Mar",
       x = "",
       y = "Temperatura (°C)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 22, color = "black")
  )



figura_2018 <- resultados_completos$rasters$`2018`

for(i in seq(1,365, 16)){
  plot(figura_2018[[i:(i + 15)]])
}

names(figura_2018)

fechas_2018 <- as.character(seq(as.Date("2018-01-01"), as.Date("2018-12-31"), by = "day"))

setdiff(fechas_2018, names(figura_2018))



mapear_temperatura_promedio_chile(datos_modelo = resultados_completos$datos_modelo$`2015`)





