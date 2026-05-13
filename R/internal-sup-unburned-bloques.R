# Block 5: legacy unburned helpers migrated from
# 2_SCRIPTS/00_FUNCTIONS/00_GENERAL_FUNCTIONS/METHODS_FUNCTIONS_BLOQUES.R
# (lines 1..897). Only the three helpers used by the legacy unburned
# pipeline are kept: polygonize_Otsu(), coverage_by_patch_raster(),
# run_scenarios(). The rest of the legacy file (tag_fire_id,
# coverage_by_fire, run_scenarios_fireonly_new, BLOQUE 3 helpers, ...)
# is intentionally NOT migrated - those are not part of the supervised
# one-year chain.
#
# These functions are package-internal: the supervised legacy unburned
# dispatcher (internal-sup-unburned-legacy.R) calls them via lexical
# scoping from the package namespace.
#
# NOT EXPORTED. NOT documented in roxygen.

#### PRIMER BLOQUE ####

##### POLIGONISATION #####

polygonize_Otsu <- function(
    burn_raster,
    python_exe = NULL,
    gdal_polygonize_script = NULL,
    tile = TRUE,
    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    dissolve_tiles = TRUE,    # <-- NUEVO: disolver costuras entre tiles por DN
    min_pixels = NULL,        # numero minimo de pixeles por poligono
    out_path = NULL,          # ruta opcional a .shp o .gpkg
    ogr2ogr_exe = NULL        # ruta a ogr2ogr.exe (Anaconda / OSGeo / QGIS)
) {
  # 1) Entrada: aceptar ruta o SpatRaster
  if (inherits(burn_raster, "character")) {
    if (length(burn_raster) != 1L) {
      stop("'burn_raster' debe ser una sola ruta a un raster (.tif).")
    }
    if (!file.exists(burn_raster)) {
      stop("El archivo especificado en 'burn_raster' no existe: ", burn_raster)
    }
    burn_raster <- terra::rast(burn_raster)
  } else if (!inherits(burn_raster, "SpatRaster")) {
    stop("'burn_raster' debe ser o bien una ruta a un .tif o un SpatRaster de terra.")
  }
  
  # Asegurar 1 banda
  if (terra::nlyr(burn_raster) > 1) {
    burn_raster <- burn_raster[[1]]
  }
  
  # CRS de salida
  cr_out <- sf::st_crs(terra::crs(burn_raster, proj = TRUE))
  if (is.na(cr_out)) {
    stop("El raster 'burn_raster' no tiene CRS definido.")
  }
  
  # Area de una celda (en unidades del CRS, en tu caso m2)
  res_xy       <- terra::res(burn_raster)
  cell_area_u2 <- abs(res_xy[1] * res_xy[2])   # p.ej. 90*90 = 8100 m2
  
  # 2) POLIGONIZAR: tiles+GDAL o fallback terra 
  polys <- NULL
  
  use_tiles <- isTRUE(tile) &&
    !is.null(python_exe) &&
    !is.null(gdal_polygonize_script) &&
    file.exists(python_exe) &&
    file.exists(gdal_polygonize_script)
  
  if (use_tiles) {
    message("Poligonizando con tiles + gdal_polygonize.py ...")
    
    # Directorio temporal para tiles
    if (is.null(out_path)) {
      tile_dir <- file.path(tempdir(), "polygonize_tiles")
    } else {
      tile_dir <- file.path(dirname(out_path), "polygonize_tiles")
    }
    if (!dir.exists(tile_dir)) dir.create(tile_dir, recursive = TRUE)
    
    ext_r <- terra::ext(burn_raster)
    tile_width  <- (ext_r[2] - ext_r[1]) / n_cols
    tile_height <- (ext_r[4] - ext_r[3]) / n_rows
    count <- 1L
    tile_shapefiles <- character(0)
    label <- "burn"
    
    for (i in 0:(n_rows - 1)) {
      for (j in 0:(n_cols - 1)) {
        xmin_tile <- max(ext_r[1] + j * tile_width  - tile_overlap, ext_r[1])
        xmax_tile <- min(ext_r[1] + (j + 1) * tile_width + tile_overlap, ext_r[2])
        ymin_tile <- max(ext_r[3] + i * tile_height - tile_overlap, ext_r[3])
        ymax_tile <- min(ext_r[3] + (i + 1) * tile_height + tile_overlap, ext_r[4])
        
        tile_crop <- terra::crop(
          burn_raster,
          terra::ext(xmin_tile, xmax_tile, ymin_tile, ymax_tile)
        )
        tile_path <- file.path(tile_dir, sprintf("tile_%s_%d.tif", label, count))
        terra::writeRaster(
          tile_crop,
          tile_path,
          overwrite = TRUE,
          datatype = "INT1U",
          NAflag   = 0,
          gdal     = c("COMPRESS=LZW")
        )
        
        shp_path <- file.path(tile_dir, sprintf("tile_%s_%d.shp", label, count))
        system(glue::glue(
          '"{python_exe}" "{gdal_polygonize_script}" "{tile_path}" -f "ESRI Shapefile" "{shp_path}" DN'
        ))
        if (file.exists(tile_path)) file.remove(tile_path)
        
        tile_shapefiles <- c(tile_shapefiles, shp_path)
        count <- count + 1L
      }
    }
    
    if (!length(tile_shapefiles)) {
      stop("No se generaron shapefiles de tiles en la poligonizacion.")
    }
    
    #  2A) Merge + dissolve en C++ con ogr2ogr (si existe) 
    if (!is.null(ogr2ogr_exe) && file.exists(ogr2ogr_exe) && isTRUE(dissolve_tiles)) {
      message("Uniendo y disolviendo tiles con ogr2ogr (GDAL, C++, usando GPKG intermedio) ...")
      
      merged_gpkg    <- file.path(tile_dir, "tiles_merged.gpkg")
      dissolved_gpkg <- file.path(tile_dir, "tiles_dissolved.gpkg")
      dissolved_shp  <- file.path(tile_dir, "tiles_dissolved.shp")
      
      # Borrar si existian
      if (file.exists(merged_gpkg))    file.remove(merged_gpkg)
      if (file.exists(dissolved_gpkg)) file.remove(dissolved_gpkg)
      if (file.exists(dissolved_shp)) {
        dissolved_base <- tools::file_path_sans_ext(dissolved_shp)
        side <- list.files(
          dirname(dissolved_shp),
          pattern = paste0("^", basename(dissolved_base),
                           "\\.(shp|dbf|shx|prj|cpg)$"),
          ignore.case = TRUE,
          full.names  = TRUE
        )
        if (length(side)) unlink(side)
      }
      
      # 1) MERGE de todos los tiles a un solo GPKG (capa 'merged')
      first <- TRUE
      for (shp in tile_shapefiles) {
        if (first) {
          cmd_merge <- sprintf(
            '"%s" -f "GPKG" "%s" "%s" -nln merged',
            ogr2ogr_exe, merged_gpkg, shp
          )
          first <- FALSE
        } else {
          cmd_merge <- sprintf(
            '"%s" -f "GPKG" -update -append "%s" "%s" -nln merged',
            ogr2ogr_exe, merged_gpkg, shp
          )
        }
        message("OGR MERGE CMD: ", cmd_merge)
        system(cmd_merge)
      }
      
      # 2) DISSOLVE por DN = 1 en el GPKG (dialecto sqlite)
      # OJO: si te da error "no such column: geom", cambia geom -> geometry
      sql <- "SELECT ST_Union(geom) AS geom, DN FROM merged WHERE DN = 1 GROUP BY DN"
      
      cmd_diss <- sprintf(
        '"%s" -f "GPKG" "%s" "%s" -dialect sqlite -sql "%s" -nln dissolved -explodecollections',
        ogr2ogr_exe, dissolved_gpkg, merged_gpkg, sql
      )
      message("OGR DISSOLVE CMD: ", cmd_diss)
      ogr_out <- try(system(cmd_diss, intern = TRUE), silent = TRUE)
      if (!inherits(ogr_out, "try-error")) {
        message(paste(ogr_out, collapse = "\n"))
      }
      
      # 3) Convertir el resultado disuelto (GPKG) a Shapefile
      cmd_to_shp <- sprintf(
        '"%s" -f "ESRI Shapefile" "%s" "%s" -nln dissolved',
        ogr2ogr_exe, dissolved_shp, dissolved_gpkg
      )
      message("OGR TO_SHP CMD: ", cmd_to_shp)
      system(cmd_to_shp)
      
      # Comprobar shapefile final
      if (!file.exists(dissolved_shp) || file.info(dissolved_shp)$size == 0) {
        stop("ogr2ogr no creo correctamente 'tiles_dissolved.shp'. Revisa los comandos OGR impresos arriba.")
      }
      
      polys <- sf::st_read(dissolved_shp, quiet = TRUE)
      polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON","MULTIPOLYGON"), ]
      polys <- sf::st_transform(polys, cr_out)
      
    } else {
      #  2B) Fallback / sin disolver: unir tiles en R 
      if (!isTRUE(dissolve_tiles)) {
        message("dissolve_tiles = FALSE -> no se disuelven costuras entre tiles. Veras lineas rectas en los bordes de teselas (puedes disolver luego en ArcGIS/QGIS).")
      } else {
        message("ogr2ogr_exe no proporcionado o no encontrado. Uniendo tiles en R sin disolver (se veran lineas de costura).")
      }
      
      polys <- do.call(
        rbind,
        lapply(tile_shapefiles, sf::st_read, quiet = TRUE)
      )
      polys <- polys[!duplicated(sf::st_as_binary(sf::st_geometry(polys))), ]
      polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON","MULTIPOLYGON"), ]
      polys <- sf::st_transform(polys, cr_out)
    }
  
    
    # Limpiar tiles intermedios (.shp, .dbf, etc.)
    for (shp in tile_shapefiles) {
      shp_base <- tools::file_path_sans_ext(shp)
      for (ext in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
        f <- paste0(shp_base, ext)
        if (file.exists(f)) file.remove(f)
      }
    }
    
  } else {
    # Fallback sin tiles (lento, pero simple) 
    message("Poligonizando con terra::as.polygons() (fallback, sin tiles/GDAL)...")
    patch_poly <- terra::as.polygons(burn_raster, dissolve = TRUE, values = TRUE)
    
    val_col <- names(patch_poly)[1]
    patch_poly <- patch_poly[!is.na(patch_poly[[val_col]]), ]
    patch_poly <- patch_poly[patch_poly[[val_col]] == 1, ]
    
    names(patch_poly)[names(patch_poly) == val_col] <- "DN"
    
    polys <- sf::st_as_sf(patch_poly)
    sf::st_crs(polys) <- cr_out
  }
  
  #  3) Limpiar DN y comprobar si hay poligonos 
  if ("DN" %in% names(polys)) {
    polys <- polys[!is.na(polys$DN) & polys$DN != 0, , drop = FALSE]
  }
  
  if (!nrow(polys)) {
    warning("No hay poligonos resultantes despues de la poligonizacion.")
    polys$patch_id <- integer(0)
    return(polys)
  }
  
  # 4) Filtro por min_pixels + metricas de area 
  area_u2 <- as.numeric(sf::st_area(polys))
  
  if (!is.null(min_pixels)) {
    if (!is.numeric(min_pixels) || length(min_pixels) != 1L || min_pixels <= 0) {
      stop("'min_pixels' debe ser un numero positivo de longitud 1.")
    }
    min_area_u2 <- min_pixels * cell_area_u2
    keep <- area_u2 >= min_area_u2
    polys <- polys[keep, , drop = FALSE]
    area_u2 <- area_u2[keep]
    
    if (!nrow(polys)) {
      warning("Tras aplicar 'min_pixels' no queda ningun poligono.")
      polys$patch_id <- integer(0)
      return(polys)
    }
  }
  
  polys$n_pix   <- round(area_u2 / cell_area_u2)
  polys$area_ha <- area_u2 / 1e4
  
  polys$patch_id <- seq_len(nrow(polys))
  
  # 5) Escritura opcional a disco 
  if (!is.null(out_path)) {
    message("Escribiendo shapefile/salida en: ", out_path)
    
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(out_path))
    
    if (ext == "shp") {
      stem <- tools::file_path_sans_ext(out_path)
      side <- list.files(
        dirname(out_path),
        pattern = paste0("^", basename(stem), "\\.(shp|dbf|shx|prj|cpg|qix|fix|sbn|sbx)$"),
        ignore.case = TRUE,
        full.names = TRUE
      )
      if (length(side)) {
        message("Borrando shapefile previo: ",
                paste(basename(side), collapse = ", "))
        unlink(side)
      }
      sf::st_write(polys, out_path, quiet = TRUE)
      
    } else if (ext == "gpkg") {
      if (file.exists(out_path)) {
        message("Borrando GPKG previo: ", basename(out_path))
        file.remove(out_path)
      }
      layer_name <- tools::file_path_sans_ext(basename(out_path))
      sf::st_write(polys, dsn = out_path, layer = layer_name, quiet = TRUE)
      
    } else {
      warning("Extension de 'out_path' no reconocida: '", ext,
              "'. Usa .shp o .gpkg. No se escribe a disco.")
    }
  }
  
  #  6) Devolver sf con patch_id 
  return(polys)
}


##### COBERTURA POR PARCHE RASTER#####
coverage_by_patch_raster <- function(
    patches,            # sf o ruta a shapefile/gpkg/geojson con patch_id (o se crea)
    ref_raster,         # SpatRaster (1 o varias capas) o lista de SpatRaster
    ref_names = NULL,   # nombres para las columnas de cada referencia
    out_path  = NULL    # ruta opcional para guardar el resultado
) {
  #  0) Preparar patches_sf 
  if (inherits(patches, "sf")) {
    patches_sf <- patches
    
  } else if (is.character(patches) && length(patches) == 1L) {
    ext <- tolower(tools::file_ext(patches))
    
    if (ext %in% c("shp", "gpkg", "geojson", "json")) {
      message("Leyendo parches desde fichero vectorial: ", patches)
      patches_sf <- sf::st_read(patches, quiet = TRUE)
      
    } else if (ext %in% c("tif", "tiff")) {
      stop(
        "Has pasado un raster (.tif) en 'patches'. ",
        "Esta funcion espera parches ya poligonizados (sf).\n",
        "-> Usa primero tu funcion de poligonizar (por ej. polygonize_Otsu())\n",
        "  y pasale aqui el shapefile resultante o el objeto sf."
      )
    } else {
      stop(
        "Extension de fichero no reconocida para 'patches': ", ext, "\n",
        "Soporto: .shp, .gpkg, .geojson, .json (o directamente un objeto sf)."
      )
    }
    
  } else {
    stop("El argumento 'patches' debe ser un objeto 'sf' o una ruta a un fichero vectorial.")
  }
  
  #  Preparar ref_raster (multi-capa) 
  # Puede ser:
  #  - SpatRaster (1 o varias capas)
  #  - lista de SpatRaster (se convierte a stack)
  if (inherits(ref_raster, "SpatRaster")) {
    ref_stack <- ref_raster
  } else if (is.list(ref_raster) &&
             all(vapply(ref_raster, inherits, logical(1), "SpatRaster"))) {
    ref_stack <- terra::rast(ref_raster)
  } else {
    stop(
      "'ref_raster' debe ser un SpatRaster (una o varias capas) ",
      "o una lista de SpatRaster con misma resolucion/extent."
    )
  }
  
  n_refs <- terra::nlyr(ref_stack)
  if (n_refs < 1) {
    stop("El SpatRaster de referencia no tiene capas.")
  }
  
  # Nombres de columnas para cada referencia
  if (is.null(ref_names)) {
    ref_names <- names(ref_stack)
    if (is.null(ref_names) || any(ref_names == "")) {
      ref_names <- paste0("ref", seq_len(n_refs))
    }
  }
  if (length(ref_names) != n_refs) {
    stop("'ref_names' debe tener la misma longitud que el numero de capas de ref_raster (",
         n_refs, ").")
  }
  
  # 2) Comprobaciones basicas
  # 2.1 Asegurar columna patch_id
  if (!"patch_id" %in% names(patches_sf)) {
    patches_sf$patch_id <- seq_len(nrow(patches_sf))
  }
  
  # 2.2 Asegurar mismo CRS (proyectamos los parches al CRS del raster)
  cr_ref <- sf::st_crs(terra::crs(ref_stack, proj = TRUE))
  if (is.na(cr_ref)) stop("El raster de referencia no tiene CRS definido.")
  if (is.na(sf::st_crs(patches_sf))) stop("Los parches no tienen CRS definido.")
  
  if (sf::st_crs(patches_sf)$wkt != cr_ref$wkt) {
    message("Transformando parches al CRS del raster de referencia...")
    patches_sf <- sf::st_transform(patches_sf, cr_ref)
  }
  
  # 3) Rasterizar patch_id sobre plantilla
  message("Rasterizando patch_id sobre el raster de referencia...")
  template  <- ref_stack[[1]]  # usamos la primera capa como plantilla
  patches_v <- terra::vect(patches_sf)
  
  patch_id_r <- terra::rasterize(
    patches_v,
    template,
    field      = "patch_id",
    background = NA_real_
  )
  
  # Extraer patch_id una sola vez
  pid_vals <- as.vector(patch_id_r[])
  valid    <- !is.na(pid_vals)
  pid_vals <- pid_vals[valid]
  
  if (!length(pid_vals)) {
    stop("No hay celdas con patch_id (?no solapan parches y raster?).")
  }
  
  pid_vals <- as.integer(round(pid_vals))
  
  #  4) Tabla de celdas totales por parche
  message("Contando celdas por parche (totales)...")
  
  tab_total <- as.data.frame(table(pid_vals), stringsAsFactors = FALSE)
  names(tab_total) <- c("patch_id", "n_total")
  tab_total$patch_id <- as.integer(tab_total$patch_id)
  
  # Esta sera la tabla acumulada
  cov_tab <- tab_total
  
  #5) Bucle sobre cada raster de referencia
  message("Calculando cobertura para cada raster de referencia...")
  
  for (k in seq_len(n_refs)) {
    rk_name <- ref_names[k]
    message("  - Procesando referencia: ", rk_name, " (capa ", k, "/", n_refs, ")")
    
    ref_k <- ref_stack[[k]]
    
    # Normalizar a 0/1 (NA -> 0, cualquier valor != 0 -> 1)
    ref_clean <- ref_k
    ref_clean[is.na(ref_clean)] <- 0
    ref_clean[ref_clean != 0]   <- 1
    
    ref_vals_all <- as.vector(ref_clean[])
    ref_vals     <- ref_vals_all[valid]  # solo donde hay patch_id
    
    # Sumar celdas "1" por parche
    tab_ref_k <- aggregate(ref_vals, by = list(patch_id = pid_vals), FUN = sum)

    n_col     <- paste0("n_", rk_name)
    names(tab_ref_k)[2] <- n_col
    tab_ref_k$patch_id  <- as.integer(tab_ref_k$patch_id)
    
    # Unir a cov_tab
    cov_tab <- dplyr::full_join(cov_tab, tab_ref_k, by = "patch_id")
    
    # Rellenar NA con 0 y calcular cobertura
    cov_col <- rk_name                 # <<--- AQUI el truco: usar el nombre tal cual
    cov_tab[[n_col]]   <- ifelse(is.na(cov_tab[[n_col]]), 0, cov_tab[[n_col]])
    cov_tab[[cov_col]] <- cov_tab[[n_col]] / cov_tab[["n_total"]]
    
    
  }
  
  # 6) Unir de vuelta a los poligonos 
  message("Uniendo resultados a los parches originales...")
  
  out <- dplyr::left_join(
    patches_sf,
    cov_tab,
    by = "patch_id"
  )
  
  #  7) Escritura opcional a disco 
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    
    ext_out <- tolower(tools::file_ext(out_path))
    
    if (ext_out == "shp") {
      stem <- tools::file_path_sans_ext(out_path)
      side <- list.files(dirname(out_path),
                         pattern = paste0("^", basename(stem),
                                          "\\.(shp|dbf|shx|prj|cpg)$"),
                         ignore.case = TRUE, full.names = TRUE)
      if (length(side)) unlink(side)
      sf::st_write(out, out_path, quiet = TRUE)
      
    } else if (ext_out == "gpkg") {
      if (file.exists(out_path)) unlink(out_path)
      layer_name <- tools::file_path_sans_ext(basename(out_path))
      sf::st_write(out, out_path, layer = layer_name, quiet = TRUE)
      
    } else if (ext_out %in% c("geojson", "json")) {
      if (file.exists(out_path)) file.remove(out_path)
      sf::st_write(out, out_path, quiet = TRUE)
      
    } else {
      warning(
        "No escribo a disco porque la extension de 'out_path' no es soportada: ",
        ext_out,
        " (usa .shp, .gpkg o .geojson)"
      )
    }
  }
  
  return(out)
}


##### NUCLEOS  #####
run_scenarios <- function(
    patches_path,
    out_dir,
    cv_fields      = c("cv_ge100"),
    buffers_m      = c(90,150),
    # --- Nucleo / boost ---
    core_thr       = 0.60,
    alpha_boost    = NULL,
    min_base_boost = NULL,
    dist_power     = 1,
    keep_hi        = 0.60,
    # --- Rescate debil opcional ---
    weak_field     = NULL,
    weak_thr       = 0.60,
    weak_wmin      = 0.50,
    # --- Descarte bajo ---
    drop_lo        = 0.20,
    target_epsg    = 3035,
    # --- IO / rendimiento ---
    keep_all_attrs = TRUE,
    keep_aux       = TRUE,
    write_audit    = FALSE,
    write_all      = TRUE,
    write_keep     = TRUE,
    write_review   = TRUE,
    driver         = c("GPKG","ESRI Shapefile"),
    # --- Distancias: fiabilidad vs velocidad ---
    dist_mode      = c("edge","centroid"),
    near_mode      = c("edge","centroid"),   # NUEVO: vecindad con poligonos o centroides
    # --- Geometria (opcional) ---
    simplify_tol   = NULL,                   # NUEVO: tolerancia en metros para simplificar (NULL = no)
    revalidate_on_write = FALSE
){
  # Block 7: require() calls removed. Imports declared in DESCRIPTION.
  driver    <- match.arg(driver)
  dist_mode <- match.arg(dist_mode)
  near_mode <- match.arg(near_mode)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  msg <- function(...) message(sprintf(...))
  
  # -------- IO helpers (igual que tu version, compacto)
  .dedup_case_insensitive <- function(nms){
    lower <- tolower(nms)
    uniq_lower <- make.unique(lower, sep = "_")
    add <- sub("^[^_]+", "", uniq_lower)
    out <- nms
    need <- add != ""
    out[need] <- paste0(nms[need], add[need])
    out
  }
  .sanitize_df_cols <- function(x, driver){
    geom_col <- attr(x, "sf_column")
    nm <- names(x); idxg <- which(nm == geom_col)
    nm_ng <- if (length(idxg)) nm[-idxg] else nm
    nm_ng <- .dedup_case_insensitive(nm_ng)
    if (driver == "ESRI Shapefile") {
      nm_ng <- toupper(substr(nm_ng, 1, 10))
      nm_ng <- make.unique(nm_ng, sep = "_")
    } else nm_ng <- make.unique(nm_ng, sep = "_")
    if (length(idxg)) nm[-idxg] <- nm_ng else nm <- nm_ng
    names(x) <- nm
    x
  }
  .harmonize_geom <- function(x){
    g <- sf::st_geometry(x)
    g <- tryCatch(suppressWarnings(sf::st_zm(g, drop = TRUE, what = "ZM")), error = function(e) g)
    inv <- !sf::st_is_valid(g)
    if (any(inv, na.rm = TRUE)) {
      g[inv] <- suppressWarnings(if (requireNamespace("lwgeom", quietly = TRUE)) sf::st_make_valid(g[inv]) else sf::st_make_valid(g[inv]))
    }
    types <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    if (any(types %in% c("GEOMETRYCOLLECTION","MULTISURFACE"))) {
      g <- suppressWarnings(sf::st_collection_extract(g, "POLYGON"))
    }
    g <- suppressWarnings(sf::st_cast(g, "MULTIPOLYGON", warn = FALSE))
    sf::st_geometry(x) <- g
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    x
  }
  .try_write <- function(expr, tries = 3, sleep0 = 0.4){
    last_err <- NULL
    for (i in seq_len(tries)){
      ok <- tryCatch({ force(expr); TRUE }, error = function(e){ last_err <<- e; FALSE })
      if (ok) return(TRUE)
      gc(); Sys.sleep(sleep0 * i)
    }
    if (!is.null(last_err)) stop(last_err$message) else FALSE
  }
  .write_vect <- function(x, path, driver){
    x <- .sanitize_df_cols(x, driver)
    if (revalidate_on_write) x <- .harmonize_geom(x)
    if (driver == "GPKG") {
      path <- sub("\\.shp$", ".gpkg", path, ignore.case = TRUE)
      lyr  <- tools::file_path_sans_ext(basename(path))
      if (file.exists(path)) unlink(path)
      ok <- .try_write(sf::st_write(x, dsn = path, layer = lyr, quiet = TRUE))
      if (!ok) {
        alt <- sub("\\.gpkg$", ".shp", path, ignore.case = TRUE)
        msg("[!]  GPKG bloqueado -> fallback a SHP: %s", basename(alt))
        path <- .write_vect(x, alt, "ESRI Shapefile")
      }
    } else {
      path <- sub("\\.gpkg$", ".shp", path, ignore.case = TRUE)
      stem <- tools::file_path_sans_ext(path)
      side <- list.files(dirname(path),
                         pattern = paste0("^", basename(stem), "\\.(shp|dbf|shx|prj|cpg|qix|fix|sbn|sbx)$"),
                         ignore.case = TRUE, full.names = TRUE)
      if (length(side)) unlink(side)
      ok <- .try_write(sf::st_write(x, dsn = path, delete_layer = TRUE, quiet = TRUE,
                                    layer_options = c("ENCODING=UTF-8")))
      if (!ok) stop("No pude escribir SHP: ", path)
    }
    path
  }
  
  # ---- s2 off (GEOS)
  old_s2 <- sf::sf_use_s2(); on.exit(sf::sf_use_s2(old_s2), add = TRUE); sf::sf_use_s2(FALSE)
  
  # ---- lectura base
  p0 <- sf::st_read(patches_path, quiet = TRUE)
  if (is.na(sf::st_crs(p0))) stop("La capa no tiene CRS. Asigna EPSG:3035 o ajusta target_epsg.")
  if (is.na(sf::st_crs(p0)$epsg) || sf::st_crs(p0)$epsg != target_epsg) p0 <- sf::st_transform(p0, target_epsg)
  
  # make_valid una vez
  inv <- !sf::st_is_valid(p0)
  if (any(inv, na.rm = TRUE)) {
    if (requireNamespace("lwgeom", quietly = TRUE)) p0[inv,] <- suppressWarnings(sf::st_make_valid(p0[inv,])) else p0[inv,] <- suppressWarnings(sf::st_make_valid(p0[inv,]))
  }
  
  # simplificacion opcional (ojo: cambia geometria, pero acelera mucho st_is_within_distance en poligonos complejos)
  if (!is.null(simplify_tol) && is.finite(simplify_tol) && simplify_tol > 0) {
    msg("Simplificando geometria (tol=%.2f m)...", simplify_tol)
    # st_simplify es rapido; preserve_topology si tienes lwgeom
    if (requireNamespace("lwgeom", quietly = TRUE)) {
      .lw_simplify <- getExportedValue("lwgeom",
                                       "st_simplify_preserve_topology")
      p0 <- suppressWarnings(.lw_simplify(p0, simplify_tol))
    } else {
      p0 <- suppressWarnings(sf::st_simplify(p0, dTolerance = simplify_tol, preserveTopology = TRUE))
    }
  }
  
  # ---- precomputes (CLAVE)
  geom0      <- sf::st_geometry(p0)
  area_m2    <- as.numeric(sf::st_area(geom0))
  area_ha    <- area_m2 / 1e4
  
  # centroides (si se usan en near o dist)
  need_cent <- (dist_mode == "centroid") || (near_mode == "centroid")
  p_cent_all <- if (need_cent) sf::st_centroid(p0) else NULL
  
  # punto representativo (siempre util para nucleos, evita st_point_on_surface por escenario)
  msg("Precalculando punto representativo por parche (st_point_on_surface) UNA vez...")
  p_rep_all <- sf::st_point_on_surface(geom0)  # sfc POINT
  
  # ---- normalizacion rapida 0-1 (cache)
  .norm01 <- function(x){
    ok <- is.finite(x)
    mx <- if (any(ok)) max(x[ok]) else 0
    prop_gt1 <- if (any(ok)) mean(x[ok] > 1) else 0
    is_pct <- (mx > 1.01) || (prop_gt1 > 0.05)
    out <- numeric(length(x))
    out[ok] <- x[ok]
    if (is_pct) out[ok] <- out[ok] / 100
    out[!ok] <- 0
    out <- pmin(pmax(out, 0), 1)
    list(v = out, is_pct = is_pct)
  }
  
  # cache weak_field (una sola vez)
  weak_cache <- NULL
  if (!is.null(weak_field) && weak_field %in% names(p0)) {
    tmp <- .norm01(p0[[weak_field]])
    weak_cache <- list(v = tmp$v, is_pct = tmp$is_pct)
  }
  
  # combos
  combos <- expand.grid(cv_field = cv_fields, buffer = buffers_m, stringsAsFactors = FALSE)
  
  results <- list()
  
  # ===== loop eficiente: 1) por cv_field (cachea cv y nucleos) 2) por buffer
  for (cvf in unique(combos$cv_field)) {
    
    if (!cvf %in% names(p0)) {
      msg("[!]  Campo %s no existe. Omito.", cvf)
      next
    }
    
    tmp <- .norm01(p0[[cvf]])
    cv_base <- tmp$v
    cv_is_pct <- tmp$is_pct
    
    nucleo_log <- cv_base >= core_thr
    nuclei_idx <- which(nucleo_log)
    nuclei_pts <- if (length(nuclei_idx)) p_rep_all[nuclei_idx] else sf::st_sfc(crs = sf::st_crs(p0))
    
    # si no hay nucleos, todo ira a review/drop sin vecindad (solo core no existe)
    for (buf in unique(combos$buffer[combos$cv_field == cvf])) {
      
      id   <- gsub("\\D","", cvf); if (identical(id, "")) id <- cvf
      slug <- sprintf("cv%s_b%dm_t%02d", id, buf, round(core_thr*100))
      
      # --- vecindad (punto clave de rendimiento)
      base_geom_for_near <- if (near_mode == "centroid") p_cent_all else p0
      vecindad <- rep(0L, nrow(p0))
      dist_m   <- rep(Inf, nrow(p0))
      cand_idx <- integer(0)
      
      if (length(nuclei_pts) > 0) {
        near <- sf::st_is_within_distance(base_geom_for_near, nuclei_pts, dist = buf)
        near_log <- lengths(near) > 0
        vecindad[near_log] <- 1L
        cand_idx <- which(near_log)
        
        if (length(cand_idx)) {
          if (dist_mode == "edge") {
            nn <- sf::st_nearest_feature(p0[cand_idx, ], nuclei_pts)
            dist_m[cand_idx] <- as.numeric(sf::st_distance(p0[cand_idx, ], nuclei_pts[nn], by_element = TRUE))
          } else {
            # centroid (mucho mas rapido)
            nn <- sf::st_nearest_feature(p_cent_all[cand_idx, ], nuclei_pts)
            dist_m[cand_idx] <- as.numeric(sf::st_distance(p_cent_all[cand_idx, ], nuclei_pts[nn], by_element = TRUE))
          }
        }
      }
      
      w_dist <- pmax(0, 1 - (dist_m / buf)^dist_power)
      w_dist[!is.finite(w_dist)] <- 0
      
      # --- min_base_boost (auto)
      min_base_boost_iter <- min_base_boost
      if (is.null(min_base_boost_iter)) {
        cand <- if (any(nucleo_log, na.rm = TRUE)) cv_base[vecindad == 1L] else cv_base
        cand <- cand[is.finite(cand)]
        if (length(cand) >= 30) {
          mb <- as.numeric(stats::quantile(cand, 0.50, na.rm = TRUE))
          min_base_boost_iter <- min(0.55, max(0.35, mb))
        } else min_base_boost_iter <- 0.45
      }
      
      # --- alpha auto
      keep_hi_prov <- if (!is.null(keep_hi)) keep_hi else 0.60
      alpha_boost_iter <- alpha_boost
      if (is.null(alpha_boost_iter)) {
        w0 <- 0.8
        alpha_boost_iter <- if (min_base_boost_iter < 1 && w0 > 0)
          (keep_hi_prov - min_base_boost_iter) / (w0 * (1 - min_base_boost_iter)) else 0.25
        alpha_boost_iter <- max(0.15, min(0.35, alpha_boost_iter))
      }
      
      # --- boost preliminar + score
      boost <- rep(0, nrow(p0))
      if (length(cand_idx)) {
        do_boost <- (cv_base[cand_idx] >= min_base_boost_iter) & (cv_base[cand_idx] < keep_hi_prov)
        idxb <- cand_idx[do_boost]
        if (length(idxb)) {
          boost[idxb] <- alpha_boost_iter * (1 - cv_base[idxb]) * w_dist[idxb]
        }
      }
      S_patch <- pmin(1, cv_base + boost)
      
      # --- keep_hi auto (si NULL)
      keep_hi_iter <- keep_hi
      if (is.null(keep_hi_iter)) {
        candS <- S_patch[(vecindad == 1L) & (!nucleo_log)]
        candS <- candS[is.finite(candS)]
        if (length(candS) < 30) candS <- S_patch[(!nucleo_log) & is.finite(S_patch)]
        if (length(candS) < 10) candS <- S_patch[is.finite(S_patch)]
        keep_hi_iter <- if (length(candS) >= 10) as.numeric(stats::quantile(candS, 0.80, na.rm = TRUE)) else 0.60
        keep_hi_iter <- max(0.50, min(0.85, keep_hi_iter))
      }
      
      # --- boost definitivo con keep_hi_iter
      if (length(cand_idx)) {
        do_boost <- (cv_base[cand_idx] >= min_base_boost_iter) & (cv_base[cand_idx] < keep_hi_iter)
        idxb <- cand_idx[do_boost]
        if (length(idxb)) boost[idxb] <- alpha_boost_iter * (1 - cv_base[idxb]) * w_dist[idxb]
      }
      S_patch <- pmin(1, cv_base + boost)
      
      # --- t_star
      den <- 1 - alpha_boost_iter * w_dist
      den[den <= 1e-9] <- 1e-9
      t_star <- rep(keep_hi_iter, nrow(p0))
      has_w <- (alpha_boost_iter * w_dist) > 0
      t_star[has_w] <- (keep_hi_iter - alpha_boost_iter * w_dist[has_w]) / den[has_w]
      
      # --- rescate debil (usa cache)
      weak_ok <- rep(FALSE, nrow(p0))
      weak_is_pct <- NA_integer_
      if (!is.null(weak_cache)) {
        weak_ok <- (vecindad == 1L) & (w_dist >= weak_wmin) & (weak_cache$v >= weak_thr)
        weak_is_pct <- as.integer(weak_cache$is_pct)
      }
      
      # --- drop_lo auto si NULL
      drop_lo_iter <- drop_lo
      if (is.null(drop_lo_iter)) {
        sp_nv <- S_patch[is.finite(S_patch) & (vecindad == 0L)]
        n_nv  <- length(sp_nv)
        if (n_nv >= 30) {
          qs <- stats::quantile(sp_nv, probs = c(0.10, 0.20), na.rm = TRUE)
          drop_lo_iter <- as.numeric(qs[1] + 0.5 * (qs[2] - qs[1]))
        } else if (n_nv >= 10) {
          drop_lo_iter <- as.numeric(stats::quantile(sp_nv, 0.15, na.rm = TRUE))
        } else drop_lo_iter <- 0.20
        drop_lo_iter <- max(0.05, min(0.30, drop_lo_iter))
      }
      
      # --- decision
      is_core   <- nucleo_log
      via_boost <- (S_patch >= keep_hi_iter) & (vecindad == 1L) & (cv_base >= min_base_boost_iter)
      via_weak  <- weak_ok
      
      keep_p <- as.integer(is_core | via_boost | via_weak)
      decision <- ifelse(keep_p == 1L, "keep",
                         ifelse(S_patch <= drop_lo_iter, "drop", "review"))
      
      # ========= construir outputs SOLO si se escriben (menos copias)
      build_out <- function(sel_idx){
        x <- p0[sel_idx, , drop = FALSE]
        x$area_ha   <- area_ha[sel_idx]
        x$cv_base   <- cv_base[sel_idx]
        x$nucleo    <- as.integer(nucleo_log[sel_idx])
        x$vecindad  <- as.integer(vecindad[sel_idx])
        x$dist_m    <- dist_m[sel_idx]
        x$w_dist    <- w_dist[sel_idx]
        x$boost     <- boost[sel_idx]
        x$S_patch   <- S_patch[sel_idx]
        x$t_star    <- t_star[sel_idx]
        x$BUF_M     <- buf
        x$THR_CORE  <- core_thr
        x$ALPHA     <- alpha_boost_iter
        x$MB_MIN    <- min_base_boost_iter
        x$DPWR      <- dist_power
        x$KEEPHI    <- keep_hi_iter
        x$DROP_LO   <- drop_lo_iter
        x$keep_p    <- keep_p[sel_idx]
        x$decision  <- decision[sel_idx]
        if (keep_aux) {
          x$CVPCT <- as.integer(cv_is_pct)
          x$WKPCT <- weak_is_pct
        } else {
          # si no quieres auxiliares, tambien puedes borrar dist/boost/etc aqui
        }
        
        if (!keep_all_attrs) {
          extra_keep <- intersect(c("CORINE_CLA","COR_ECO_LA","CORINE_CLASS","COR_ECO_LABEL",
                                    "ECO_CLASS","unit_id2","Label"), names(x))
          keep_cols <- c(extra_keep,
                         "area_ha","cv_base","nucleo","vecindad","dist_m","w_dist","boost","S_patch","t_star",
                         "BUF_M","THR_CORE","ALPHA","MB_MIN","DPWR","KEEPHI","DROP_LO","keep_p","decision")
          x <- x[, c(keep_cols, attr(x, "sf_column")), drop = FALSE]
        }
        x
      }
      
      idx_all    <- seq_len(nrow(p0))
      idx_keep   <- which(keep_p == 1L)
      idx_review <- which(keep_p == 0L & decision == "review")
      
      part_p     <- if (write_all)   build_out(idx_all)    else NULL
      part_keep  <- if (write_keep && length(idx_keep))   build_out(idx_keep)   else NULL
      part_rev   <- if (write_review && length(idx_review)) build_out(idx_review) else NULL
      
      # ========= escritura
      stem <- paste0("sc_", slug)
      ext  <- ifelse(driver=="GPKG","gpkg","shp")
      p_all_path  <- file.path(out_dir, paste0(stem, "_patches.", ext))
      p_keep_path <- file.path(out_dir, paste0(stem, "_keep.",    ext))
      p_rev_path  <- file.path(out_dir, paste0(stem, "_review.",  ext))
      
      if (!is.null(part_p))    p_all_path  <- .write_vect(part_p,    p_all_path,  driver) else p_all_path  <- NA
      if (!is.null(part_keep)) p_keep_path <- .write_vect(part_keep, p_keep_path, driver) else p_keep_path <- NA
      if (!is.null(part_rev))  p_rev_path  <- .write_vect(part_rev,  p_rev_path,  driver) else p_rev_path  <- NA
      
      if (write_audit && length(nuclei_idx)) {
        nuclei_poly <- if (!is.null(part_p)) part_p[part_p$nucleo==1L, , drop=FALSE] else build_out(nuclei_idx)
        if (nrow(nuclei_poly)) {
          # Block 7: sf::st_unary_union no longer exists in current sf;
          # fall back to lwgeom::st_unary_union when available, otherwise
          # sf::st_union.
          .unary_u <- if (requireNamespace("lwgeom", quietly = TRUE))
            getExportedValue("lwgeom", "st_unary_union") else sf::st_union
          nuc_sf <- sf::st_sf(geometry = .unary_u(nuclei_poly))
          buf_sf <- sf::st_buffer(nuc_sf, buf)
          .write_vect(nuc_sf, file.path(out_dir, paste0(stem, "_nuclei.", ext)), driver)
          .write_vect(buf_sf,  file.path(out_dir, paste0(stem, "_buf.",    ext)), driver)
        }
      }
      
      msg(" %s: %s%s%s",
          slug,
          if (!is.na(p_all_path)) basename(p_all_path) else "<no_all>",
          if (!is.na(p_keep_path)) paste0(", ", basename(p_keep_path)) else "",
          if (!is.na(p_rev_path))  paste0(", ", basename(p_rev_path)) else "")
      
      results[[slug]] <- list(all=p_all_path, keep=p_keep_path, review=p_rev_path)
    }
  }
  
  invisible(results)
}
