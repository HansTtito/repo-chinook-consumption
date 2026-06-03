# Temperature data processing from Copernicus GLORYS12 NetCDF files
#
# NOTE: This script requires raw NetCDF files from Copernicus Marine Service
# (GLORYS12 reanalysis, variable thetao) which are not included in this
# repository due to size and redistribution restrictions.
# The two output files are already provided:
#   - data/data_raw/temperature/temperature_time_serie_daily.csv  (used downstream)
#   - data/data_raw/temperature/total_data_temperature.rds
# This script is kept for transparency and reproducibility of data processing.

library(terra)
library(sf)
library(reshape2)
library(dplyr)
library(ncdf4)
library(ggplot2)
library(viridis)

# =============================================================================
# FUNCIÓN 1: CARGAR Y ANALIZAR ARCHIVO NETCDF
# =============================================================================

# Helper: convierte tiempo CF-compliant a POSIXct
convert_cf_time <- function(vals, units, calendar = "gregorian", tz = "UTC") {
  if (length(vals) == 0 || is.null(units)) return(as.POSIXct(character(0), tz = tz))

  unit_base  <- tolower(sub("\\s+since.*$", "", units))
  origin_str <- sub("^.*since\\s*", "", units, ignore.case = TRUE)
  origin_str <- gsub("T", " ", origin_str)

  if (!grepl("\\d{2}:\\d{2}", origin_str)) origin_str <- paste0(origin_str, " 00:00:00")

  origin <- as.POSIXct(origin_str, tz = tz)
  if (is.na(origin)) stop("No se pudo parsear el origen temporal desde: ", units)

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

  if (!tolower(calendar) %in% c("gregorian", "standard", "proleptic_gregorian", "", NA)) {
    warning("Calendario '", calendar, "' no completamente soportado.")
  }

  origin + vals * mult
}

analizar_estructura_netcdf <- function(archivo_path, poligono_path = NULL, variable = NULL) {
  nc <- nc_open(archivo_path)
  on.exit(nc_close(nc))

  variables <- names(nc$var)

  dims <- lapply(nc$dim, function(x) {
    list(
      name  = x$name,
      size  = x$len,
      units = x$units,
      vals  = if (x$len <= 20) x$vals else paste0("(", x$len, " valores)")
    )
  })

  depth_vals <- if ("depth" %in% names(nc$dim)) nc$dim$depth$vals else NA

  if ("time" %in% names(nc$dim)) {
    time_raw   <- if (!is.null(nc$dim$time$vals)) nc$dim$time$vals else ncvar_get(nc, "time")
    time_units <- nc$dim$time$units %||% ncatt_get(nc, "time", "units")$value
    time_cal   <- nc$dim$time$calendar %||% ncatt_get(nc, "time", "calendar")$value %||% "gregorian"
    time_posix <- convert_cf_time(time_raw, units = time_units, calendar = time_cal, tz = "UTC")
  } else {
    time_raw   <- numeric(0)
    time_units <- NA_character_
    time_cal   <- NA_character_
    time_posix <- as.POSIXct(character(0), tz = "UTC")
  }

  if (!is.null(variable)) {
    if (variable %in% variables) {
      r <- rast(archivo_path, subds = variable)
    } else {
      stop("Variable '", variable, "' no encontrada. Disponibles: ", paste(variables, collapse = ", "))
    }
  } else {
    r <- tryCatch(rast(archivo_path), error = function(e) NULL)
  }

  if (!is.null(r)) {
    d <- tryCatch(depth(r), error = function(e) rep(NA, nlyr(r)))
    t <- tryCatch(time(r),  error = function(e) rep(as.Date(NA), nlyr(r)))
  } else {
    d <- numeric(0)
    t <- as.Date(character(0))
  }

  profundidades_unicas <- if (!all(is.na(depth_vals))) depth_vals else if (!all(is.na(d))) unique(d) else numeric(0)
  fechas_unicas        <- if (length(time_posix) > 0) unique(as.Date(time_posix)) else if (length(t) > 0 && !all(is.na(t))) unique(t) else as.Date(character(0))

  n_profundidades <- length(profundidades_unicas)
  n_dias          <- length(fechas_unicas)

  estructura <- list(
    archivo              = basename(archivo_path),
    dimensiones          = if (!is.null(r)) c(filas = nrow(r), columnas = ncol(r)) else c(filas = NA, columnas = NA),
    total_capas          = if (!is.null(r)) nlyr(r) else NA,
    n_layers             = if (!is.null(r)) nlyr(r) else NA,
    resolucion           = if (!is.null(r)) res(r) else c(NA, NA),
    extension            = if (!is.null(r)) as.vector(ext(r)) else c(NA, NA, NA, NA),
    crs                  = if (!is.null(r)) crs(r) else NA,
    profundidades_unicas = profundidades_unicas,
    depth_vals           = depth_vals,
    n_profundidades      = n_profundidades,
    rango_profundidades  = if (length(profundidades_unicas) > 0) range(profundidades_unicas) else c(NA, NA),
    fechas_unicas        = fechas_unicas,
    n_dias               = n_dias,
    rango_fechas         = if (length(fechas_unicas) > 0) range(fechas_unicas) else as.Date(c(NA, NA)),
    time_raw             = time_raw,
    time_units           = time_units,
    time_calendar        = time_cal,
    time                 = time_posix,
    dimensiones_netcdf   = dims,
    variables_disponibles = variables,
    variable_seleccionada = variable,
    raster_original      = r
  )

  if (!is.null(poligono_path)) {
    poligono          <- st_read(poligono_path, quiet = TRUE)
    estructura$poligono <- vect(poligono)
  }

  return(estructura)
}

# =============================================================================
# FUNCIÓN 2: FILTRAR Y PROCESAR PROFUNDIDADES
# =============================================================================

procesar_profundidades <- function(estructura, prof_min = 20, prof_max = 40,
                                   recortar_espacial = TRUE) {
  r                    <- estructura$raster_original
  profundidades_unicas <- estructura$profundidades_unicas
  n_profundidades      <- estructura$n_profundidades
  n_dias               <- estructura$n_dias

  indices_prof_objetivo <- which(profundidades_unicas >= prof_min &
                                   profundidades_unicas <= prof_max)
  valores_prof_objetivo <- profundidades_unicas[indices_prof_objetivo]

  if (length(valores_prof_objetivo) == 0) {
    stop("No se encontraron profundidades en el rango especificado")
  }

  if (recortar_espacial && !is.null(estructura$poligono)) {
    r <- crop(r, estructura$poligono)
    r <- mask(r, estructura$poligono)
  }

  indices_por_dia <- list()
  for (dia in 1:n_dias) {
    indices_por_dia[[dia]] <- sapply(indices_prof_objetivo, function(prof_idx) {
      (dia - 1) * n_profundidades + prof_idx
    })
  }

  rasters_diarios <- list()
  for (dia in 1:n_dias) {
    capas_dia          <- r[[indices_por_dia[[dia]]]]
    rasters_diarios[[dia]] <- mean(capas_dia, na.rm = TRUE)
  }

  r_avg           <- rast(rasters_diarios)
  fechas_ordenadas <- sort(estructura$fechas_unicas)
  names(r_avg)    <- format(fechas_ordenadas, "%Y-%m-%d")

  return(r_avg)
}

# =============================================================================
# FUNCIÓN 3: EXTRAER SERIES TEMPORALES POR GRILLA
# =============================================================================

extraer_series_temporales <- function(raster_stack, formato = "largo") {
  coords             <- xyFromCell(raster_stack, 1:ncell(raster_stack))
  temperaturas_matrix <- values(raster_stack)
  grillas_validas    <- which(!is.na(temperaturas_matrix[, 1]))

  coords_validas <- data.frame(
    grilla_id = 1:length(grillas_validas),
    longitud  = coords[grillas_validas, 1],
    latitud   = coords[grillas_validas, 2]
  )

  temperaturas_validas <- temperaturas_matrix[grillas_validas, ]

  if (formato == "ancho") {
    df_resultado <- data.frame(
      grilla_id = coords_validas$grilla_id,
      longitud  = coords_validas$longitud,
      latitud   = coords_validas$latitud,
      temperaturas_validas
    )
    colnames(df_resultado)[4:ncol(df_resultado)] <- names(raster_stack)

  } else {
    df_ancho <- data.frame(
      grilla_id = coords_validas$grilla_id,
      longitud  = coords_validas$longitud,
      latitud   = coords_validas$latitud,
      temperaturas_validas
    )
    colnames(df_ancho)[4:ncol(df_ancho)] <- names(raster_stack)

    df_resultado <- melt(df_ancho,
                         id.vars      = c("grilla_id", "longitud", "latitud"),
                         variable.name = "fecha",
                         value.name   = "temperatura")
    df_resultado$fecha <- as.Date(gsub("X", "", as.character(df_resultado$fecha)))
  }

  return(df_resultado)
}

# =============================================================================
# FUNCIÓN 4: GENERAR SERIES INDIVIDUALES POR GRILLA
# =============================================================================

generar_series_individuales <- function(datos_temperatura, directorio_salida = NULL) {
  grillas_unicas  <- unique(datos_temperatura$grilla_id)
  n_grillas       <- length(grillas_unicas)
  series_por_grilla <- list()

  for (i in 1:n_grillas) {
    grilla_id    <- grillas_unicas[i]
    datos_grilla <- datos_temperatura[datos_temperatura$grilla_id == grilla_id, ]
    datos_grilla <- datos_grilla[order(datos_grilla$fecha), ]

    serie <- data.frame(
      grilla_id  = grilla_id,
      longitud   = datos_grilla$longitud[1],
      latitud    = datos_grilla$latitud[1],
      fecha      = datos_grilla$fecha,
      temperatura = datos_grilla$temperatura,
      dia_del_año = as.numeric(format(datos_grilla$fecha, "%j"))
    )

    series_por_grilla[[paste0("grilla_", grilla_id)]] <- serie

    if (!is.null(directorio_salida)) {
      if (!dir.exists(directorio_salida)) dir.create(directorio_salida, recursive = TRUE)
      write.csv(serie, file.path(directorio_salida, paste0("serie_grilla_", grilla_id, ".csv")),
                row.names = FALSE)
    }
  }

  return(series_por_grilla)
}

# =============================================================================
# FUNCIÓN 5: PREPARAR DATOS PARA MODELADO
# =============================================================================

preparar_datos_modelo <- function(series_por_grilla, incluir_coordenadas = TRUE) {
  n_grillas      <- length(series_por_grilla)
  longitud_serie <- nrow(series_por_grilla[[1]])

  matriz_temperaturas <- matrix(NA, nrow = n_grillas, ncol = longitud_serie)

  metadatos_grillas <- data.frame(
    indice    = 1:n_grillas,
    grilla_id = numeric(n_grillas),
    longitud  = numeric(n_grillas),
    latitud   = numeric(n_grillas)
  )

  for (i in 1:n_grillas) {
    serie <- series_por_grilla[[i]]
    matriz_temperaturas[i, ]        <- serie$temperatura
    metadatos_grillas$grilla_id[i]  <- serie$grilla_id[1]
    metadatos_grillas$longitud[i]   <- serie$longitud[1]
    metadatos_grillas$latitud[i]    <- serie$latitud[1]
  }

  rownames(matriz_temperaturas) <- paste0("grilla_", metadatos_grillas$grilla_id)
  colnames(matriz_temperaturas) <- paste0("dia_", 1:longitud_serie)

  fechas <- series_por_grilla[[1]]$fecha

  resultado <- list(
    X          = matriz_temperaturas,
    metadatos  = metadatos_grillas,
    fechas     = fechas,
    n_grillas  = n_grillas,
    n_dias     = longitud_serie,
    estadisticas = list(
      temp_min   = min(matriz_temperaturas, na.rm = TRUE),
      temp_max   = max(matriz_temperaturas, na.rm = TRUE),
      temp_media = mean(matriz_temperaturas, na.rm = TRUE),
      temp_sd    = sd(matriz_temperaturas, na.rm = TRUE)
    )
  )

  if (incluir_coordenadas) {
    matriz_coords <- cbind(metadatos_grillas$longitud, metadatos_grillas$latitud)
    colnames(matriz_coords) <- c("longitud", "latitud")
    rownames(matriz_coords) <- rownames(matriz_temperaturas)
    resultado$coords <- matriz_coords
  }

  return(resultado)
}

# =============================================================================
# FUNCIÓN PARA OBTENER PROMEDIO DIARIO TOTAL
# =============================================================================

obtener_promedio_diario_total <- function(datos_modelo) {
  promedios_diarios <- colMeans(datos_modelo$X, na.rm = TRUE)

  serie_promedio <- data.frame(
    Day         = 1:length(promedios_diarios),
    Temperature = promedios_diarios
  )

  if (!is.null(datos_modelo$fechas)) {
    serie_promedio$Date <- datos_modelo$fechas
  }

  return(serie_promedio)
}


archivos_temperatura <- list.files(path = "data/data_raw/temperature/nc_files/",
                                   pattern = "mercator", full.names = TRUE)
poligono_shp <- "data/spatial/poligono_costero_60mn.shp"

procesar_todos_los_años <- function(archivos, poligono_path, prof_min = 20, prof_max = 40) {

  promedios_diarios_por_año <- list()

  for (i in seq_along(archivos)) {

    año <- gsub(".*mercator_(\\d{4})\\.nc", "\\1", basename(archivos[i]))
    cat("Procesando año:", año, "(", i, "/", length(archivos), ")\n")

    estructura  <- analizar_estructura_netcdf(archivos[i],
                                              variable      = "thetao",
                                              poligono_path = poligono_path)
    raster_año  <- procesar_profundidades(estructura, prof_min, prof_max)
    datos_largo <- extraer_series_temporales(raster_año, formato = "largo")

    series_individuales <- generar_series_individuales(datos_largo)
    datos_modelo        <- preparar_datos_modelo(series_individuales)
    promedio_diario     <- obtener_promedio_diario_total(datos_modelo)
    promedio_diario$año <- año

    promedios_diarios_por_año[[año]] <- promedio_diario

    cat("Año", año, "completado:", nlyr(raster_año), "días,",
        length(series_individuales), "grillas\n\n")
  }

  resultado <- list(
    promedios_diarios = promedios_diarios_por_año,
    años_procesados   = names(promedios_diarios_por_año)
  )

  cat("Años procesados:", paste(resultado$años_procesados, collapse = ", "), "\n")

  return(resultado)
}

resultados_completos <- procesar_todos_los_años(archivos_temperatura, poligono_shp)

promedio_diario_completo <- do.call(rbind, resultados_completos$promedios_diarios)
rownames(promedio_diario_completo) <- NULL

names(promedio_diario_completo) <- c("Day", "Temperature", "Date", "year")

# Imputar valor faltante 08 mayo 2018 (S3.4)
idx_may08 <- which(promedio_diario_completo$Date == as.Date("2018-05-08"))
if (length(idx_may08) == 1) {
  promedio_diario_completo$Temperature[idx_may08] <- mean(c(
    promedio_diario_completo$Temperature[idx_may08 - 1],
    promedio_diario_completo$Temperature[idx_may08 + 1]
  ))
}

write.csv(promedio_diario_completo,
          "data/data_raw/temperature/temperature_time_serie_daily.csv",
          row.names = FALSE)

saveRDS(resultados_completos, "data/data_raw/temperature/total_data_temperature.rds")
