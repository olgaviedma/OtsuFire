#' Process Burned Area Rasters Using Otsu Thresholding or Percentile Clipping
#'
#' @description
#' This function computes binary burned area masks from a severity index raster (RBR or dNBR)
#' using Otsu's thresholding or percentile clipping. The segmentation can be applied
#' to the full raster or stratified by:
#' - CORINE land cover classes (`corine_raster_path`)
#' - WWF ecoregions (`ecoregion_shapefile_path`)
#' - or the intersection of both (CORINE ? Ecoregion)
#'
#' ## Key features:
#' - Supports Otsu thresholding after pre-filtering by minimum index values (`otsu_thresholds`) or using original values
#' - Supports percentile-based thresholding via `trim_percentiles`
#' - Enforces a minimum threshold (`min_otsu_threshold_value`) when Otsu yields low values
#' - Allows stratification by CORINE, ecoregions, or their intersection
#' - Enables CORINE reclassification using a custom `reclass_matrix`
#' - Automatically computes RBR or dNBR if `nbr_pre_path` and `nbr_post_path` are provided
#' - Outputs:
#'   - Burned area binary rasters
#'   - Polygon shapefiles (dissolved by threshold label)
#'   - Histogram and inter-class variance plots per threshold
#'   - Threshold log files
#'
#' @details
#' - Raster values are internally rescaled to [0, 255] before Otsu thresholding.
#' - Histogram smoothing and variance curves are used for enhanced threshold detection.
#' - Output shapefiles can be in ESRI Shapefile or GeoJSON format.
#' - Intermediate tile shapefiles and rasters are cleaned after processing.
#'
#' @name process_otsu_rasters
#' @rdname process_otsu_rasters
#'
#' @param raster_path Path to a single-band RBR or dNBR raster.
#' @param nbr_pre_path,nbr_post_path Optional. Paths to pre- and post-fire NBR rasters for RBR/dNBR calculation.
#' @param output_dir Directory to save output files.
#' @param year Optional year label.
#' @param otsu_thresholds Numeric vector of minimum values to filter raster before applying Otsu.
#' @param trim_percentiles Optional data.frame with `min` and `max` columns to define percentile clipping thresholds.
#' @param use_original Logical. Use raw values without filtering (ignored if `trim_percentiles` is set).
#' @param corine_raster_path Path to CORINE raster for stratified segmentation.
#' @param reclassify_corine Logical. If TRUE, reclassifies CORINE using `reclass_matrix`.
#' @param reclass_matrix Matrix with original and new values for CORINE reclassification.
#' @param peninsula_shapefile Shapefile used to crop CORINE before reclassification.
#' @param output_corine_raster_dir Directory to save reclassified CORINE raster.
#' @param output_corine_vector_dir Directory to save reclassified CORINE shapefile.
#' @param reproject Logical. Reproject reclassified CORINE raster to EPSG:3035.
#' @param resolution Numeric. Target resolution for reprojected CORINE raster.
#' @param corine_classes Optional vector of reclassified CORINE classes to keep.
#' @param ecoregion_shapefile_path Path to WWF ecoregions shapefile.
#' @param ecoregion_field Column in ecoregion shapefile used for class labels.
#' @param ecoregion_classes Optional. Vector of selected ecoregion class names.
#' @param segment_by_intersection Logical. If TRUE, intersects CORINE and ecoregions.
#' @param min_otsu_threshold_value Minimum acceptable threshold. Used if Otsu value is lower.
#' @param python_exe Path to Python executable.
#' @param gdal_polygonize_script Path to `gdal_polygonize.py` script.
#' @param gdalwarp_path Path to `gdalwarp` executable (default = "gdalwarp").
#' @param n_rows,n_cols Number of rows/columns for raster tiling.
#' @param tile_overlap Overlap buffer (in map units) for tiles.
#' @param tile Logical. Apply raster tiling before polygonization.
#' @param index_type "RBR" or "dNBR". Used if `nbr_pre_path`/`nbr_post_path` are provided.
#' @param output_format Output vector file format: "shp" or "geojson".
#'
#' @return Named list of `sf` polygons grouped by segmentation. Output files saved to `output_dir`.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic Otsu thresholds (global)
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   otsu_thresholds = c(0, 50, 100),
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 2: Global with percentile trimming
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   trim_percentiles = data.frame(min = 0.01, max = 0.99),
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py",
#'   output_format = "geojson"
#' )
#'
#' # Example 3: Stratified by reclassified CORINE
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   corine_raster_path = "data/corine_1990_etrs89.tif",
#'   reclassify_corine = TRUE,
#'   reclass_matrix = matrix(c(
#'     22, 1,  # grasslands ? 1
#'     26, 1,
#'     32, 1,
#'     27, 2,  # shrublands ? 2
#'     29, 2,
#'     33, 3,  # burnt areas ? 3
#'     23, 4,  # broadleaved forest ? 4
#'     25, 5,  # mixed forest ? 5
#'     24, 6   # coniferous ? 6
#'   ), ncol = 2, byrow = TRUE),
#'   corine_classes = c(1, 2, 4, 5, 6),
#'   peninsula_shapefile = "data/peninsula_clip.shp",
#'   output_corine_raster_dir = "output/corine",
#'   output_corine_vector_dir = "output/corine",
#'   otsu_thresholds = c(50),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 4: Stratified by WWF ecoregions (default field: "EnS_name")
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   ecoregion_shapefile_path = "data/ecoregions_olson.shp",
#'   otsu_thresholds = c(50),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 5: Stratified by WWF ecoregions with custom field
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   ecoregion_shapefile_path = "data/ecoregions_olson.shp",
#'   ecoregion_field = "EnZ_name",
#'   ecoregion_classes = c("Iberian Conifer Forests", "Mediterranean Shrublands"),
#'   otsu_thresholds = c(50),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 6: Stratified by CORINE ? Ecoregion intersection
#' process_otsu_rasters(
#'   raster_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   corine_raster_path = "data/corine_1990_etrs89.tif",
#'   reclassify_corine = TRUE,
#'   reclass_matrix = matrix(c(
#'     22, 1, 26, 1, 32, 1,  # grasslands
#'     27, 2, 29, 2,         # shrublands
#'     23, 4, 25, 5, 24, 6   # forests
#'   ), ncol = 2, byrow = TRUE),
#'   corine_classes = c(1, 2, 4, 5, 6),
#'   peninsula_shapefile = "data/peninsula_clip.shp",
#'   output_corine_raster_dir = "output/corine",
#'   output_corine_vector_dir = "output/corine",
#'   ecoregion_shapefile_path = "data/ecoregions_olson.shp",
#'   ecoregion_field = "EnZ_name",
#'   segment_by_intersection = TRUE,
#'   otsu_thresholds = c(50),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#' }
#' @importFrom terra rast crop mask ext resample ifel values freq writeRaster ncell
#' @importFrom sf st_read st_write st_transform st_geometry_type st_as_binary st_make_valid
#' @importFrom stats quantile
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom dplyr all_of any_of where first
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom stringr str_detect str_replace str_extract
#' @importFrom OtsuSeg otsu_threshold_smoothed
#' @export

utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":=","ecoregion_classes_internal", "unit_id2", "CORINE_CLASS",
  "CORINE_YEAR", "ECO_CLASS", ".data"
))

process_otsu_rasters <- function(
    raster_path = NULL,
    nbr_pre_path = NULL,
    nbr_post_path = NULL,
    output_dir,
    year = NULL,
    otsu_thresholds = c(0, 25, 50, 75, 100),
    trim_percentiles = NULL,
    use_original = FALSE,
    corine_raster_path = NULL,
    reclassify_corine = FALSE,
    reclass_matrix = NULL,
    peninsula_shapefile = NULL,
    output_corine_raster_dir = NULL,
    output_corine_vector_dir = NULL,
    reproject = TRUE,
    resolution = 100,
    ecoregion_shapefile_path = NULL,
    ecoregion_field = NULL,
    ecoregion_classes = NULL,
    min_otsu_threshold_value = NULL,
    python_exe,
    gdal_polygonize_script,
    gdalwarp_path = "gdalwarp",
    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    tile = TRUE,
    index_type = "RBR",
    segment_by_intersection = FALSE,
    corine_classes = NULL,
    vectorize = FALSE,
    output_format = c("shp", "geojson")){

  output_format <- match.arg(output_format, choices = c("shp", "geojson"))

  if (is.null(output_dir)) stop("'output_dir' must be provided.")
  if (is.null(python_exe)) stop("'python_exe' must be provided.")

  if (reclassify_corine && is.null(peninsula_shapefile)) {
    stop("You set 'reclassify_corine = TRUE' but did not provide 'peninsula_shapefile'.")
  }
  if (reclassify_corine && is.null(output_corine_raster_dir)) {
    stop("You must specify 'output_corine_raster_dir' to save the reclassified CORINE raster.")
  }

  tile <- as.logical(tile)[1]
  if (is.na(tile)) stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")

  # Infer year if not provided
  if (is.null(year)) {
    if (!is.null(raster_path)) {
      raster_base <- tools::file_path_sans_ext(basename(raster_path))
      year_match <- stringr::str_extract(raster_base, "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    } else if (!is.null(output_dir)) {
      year_match <- stringr::str_extract(basename(normalizePath(output_dir)), "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    } else {
      stop("You must provide either 'raster_path' or 'output_dir' to infer the year.")
    }
  }


  # Ensure only one thresholding method is active
  if (!is.null(trim_percentiles) && !is.null(otsu_thresholds)) {
    stop("You cannot use 'trim_percentiles' and 'otsu_thresholds' at the same time. Use only one.")
  }

  # Determine raster source
  if (!is.null(nbr_pre_path) && !is.null(nbr_post_path)) {
    message("Calculating burn severity index from NBR pre and post rasters...")

    if (!file.exists(nbr_pre_path)) stop("Pre-fire NBR raster not found: ", nbr_pre_path)
    if (!file.exists(nbr_post_path)) stop("Post-fire NBR raster not found: ", nbr_post_path)

    nbr_pre <- terra::rast(nbr_pre_path)
    nbr_post <- terra::rast(nbr_post_path)

    # Validate dimensions
    if (!terra::compareGeom(nbr_pre, nbr_post, stopOnError = FALSE)) {
      stop("NBR pre and post rasters must have the same extent, resolution, and CRS.")
    }

    if (tolower(index_type) == "dnbr") {
      r <- (nbr_pre - nbr_post) * 1000
    } else if (tolower(index_type) == "rbr") {
      r <- (nbr_pre - nbr_post) * 1000 / (nbr_pre + 1.001)
    } else {
      stop("`index_type` must be either 'RBR' or 'dNBR'")
    }

  } else if (!is.null(raster_path)) {
    r <- terra::rast(raster_path)[[1]]  # ensure single band
  } else {
    stop("You must provide either 'raster_path' or both 'nbr_pre_path' and 'nbr_post_path'.")
  }


  # Asegurar carpeta principal de salida
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


  # Initialize results
  results <- list()
  threshold_log <- data.frame(Label = character(), RealThreshold = numeric(), stringsAsFactors = FALSE)

  # === Cargar corine si se requiere ===
  if (!is.null(corine_raster_path)) {
    corine_raster_input_original <- corine_raster_path  # proteger la ruta original
    corine_raster <- terra::rast(corine_raster_path)
    message("Unique values in CORINE raster: ", paste(unique(terra::values(corine_raster)), collapse = ", "))

    corine_filename <- basename(corine_raster_path)
    corine_year <- stringr::str_extract(corine_filename, "\\d{4}")

    if (is.na(corine_year)) {
      corine_year_short <- stringr::str_extract(corine_filename, "corine(\\d{2})")
      if (!is.na(corine_year_short)) {
        year_digits <- stringr::str_extract(corine_year_short, "\\d{2}")
        corine_year <- paste0("19", year_digits)
      }
    }

    if (is.na(corine_year)) {
      warning("Could not extract CORINE year from filename.")
    } else {
      message("Extracted CORINE year: ", corine_year)
    }

    if (reclassify_corine) {
      message("Cropping CORINE by mask...")

      if (is.null(peninsula_shapefile)) stop("Missing 'peninsula_shapefile'")
      if (is.null(output_corine_raster_dir)) stop("Missing 'output_corine_raster_dir'")

      peninsula <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula_proj <- sf::st_transform(peninsula, crs = sf::st_crs(corine_raster)) |> sf::st_make_valid()

      corine_r <- terra::rast(corine_raster_path)
      cropped <- terra::crop(corine_r, peninsula_proj)
      masked <- terra::mask(cropped, peninsula_proj)

      if (is.null(reclass_matrix)) stop("Provide 'reclass_matrix' for reclassification")

      message("Reclassifying CORINE using the provided matrix...")
      reclass_matrix <- as.matrix(reclass_matrix)

      original_vals <- unique(masked[])
      missing_classes <- setdiff(original_vals, reclass_matrix[, 1])
      if (length(missing_classes) > 0) {
        warning("Unreclassified CORINE classes: ", paste(missing_classes, collapse = ", "))
      }

      reclassed_r <- terra::classify(masked, rcl = reclass_matrix, include.lowest = TRUE, right = NA)

      valid_classes <- unique(reclass_matrix[, 2])
      vals <- reclassed_r[]
      vals[!vals %in% valid_classes] <- NA
      reclassed_r[] <- vals

      class_r <- reclassed_r
      message("Unique values in reclassified CORINE raster: ", paste(unique(terra::values(class_r)), collapse = ", "))

      # GENERATE A NEW OUTPUT NAME AUTOMATICALLY
      input_name <- tools::file_path_sans_ext(basename(corine_raster_path))
      out_filename <- paste0("reclassified_", input_name, "_", format(Sys.Date(), "%Yj"), ".tif")
      out_reproj <- file.path(output_corine_raster_dir, out_filename)

      # Verificar que no se sobrescriba el raster original
      if (normalizePath(out_reproj, mustWork = FALSE) == normalizePath(corine_raster_path, mustWork = FALSE)) {
        stop("ERROR: Output raster path would overwrite the original CORINE raster.")
      }

      if (reproject) {
        tmp_unproj <- tempfile(fileext = ".tif")
        terra::writeRaster(class_r, tmp_unproj, overwrite = TRUE)
        if (file.exists(out_reproj)) file.remove(out_reproj)

        cmd <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:3035" -tr {resolution} {resolution} -r near -co "COMPRESS=LZW" "{tmp_unproj}" "{out_reproj}"')
        system(cmd)

        if (!file.exists(out_reproj)) stop("ERROR: gdalwarp did not create output: ", out_reproj)
      } else {
        terra::writeRaster(class_r, out_reproj, overwrite = TRUE,
                           datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))
        if (!file.exists(out_reproj)) stop("ERROR: writeRaster failed: ", out_reproj)
      }

      tryCatch({
        test_r <- terra::rast(out_reproj)
      }, error = function(e) {
        stop("ERROR: Cannot read reclassified raster: ", out_reproj)
      })

      message("saved reclassified CORINE raster: ", out_reproj)

      # VECTORIZE
      if (vectorize && file.exists(python_exe) && file.exists(gdal_polygonize_script)) {
        out_shp <- file.path(output_corine_vector_dir, paste0(tools::file_path_sans_ext(out_filename), ".shp"))

        tryCatch({
          existing_files <- list.files(output_corine_vector_dir,
                                       pattern = paste0("^", tools::file_path_sans_ext(out_filename), "\\.(shp|shx|dbf|prj|cpg)$"),
                                       full.names = TRUE)
          if (length(existing_files) > 0) file.remove(existing_files)
        }, error = function(e) warning("Could not clean old vector files."))

        cmd_vec <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{out_reproj}" -f "ESRI Shapefile" "{out_shp}" DN')
        system(cmd_vec)

        if (!file.exists(out_shp)) stop("ERROR: Vector shapefile not created: ", out_shp)

        corine_vect <- tryCatch({
          sf::st_read(out_shp, quiet = TRUE)
        }, error = function(e) stop("ERROR: Could not read vector file: ", conditionMessage(e)))

        if (!"DN" %in% names(corine_vect)) stop("ERROR: Field 'DN' missing in shapefile")
        names(corine_vect)[names(corine_vect) == "DN"] <- "CORINE_CLASS"

        sf::st_write(corine_vect, out_shp, delete_layer = TRUE, quiet = TRUE)
        message("vectorized reclassified CORINE: ", out_shp)
        message("Classes: ", paste(unique(corine_vect$CORINE_CLASS), collapse = ", "))
      }

      # assign updated raster
      corine_raster_reclassed <- terra::rast(out_reproj)
    }

    # Post-process: resample and mask
    corine_resampled <- terra::resample(corine_raster_reclassed, r, method = "near")
    corine_masked <- terra::mask(corine_resampled, r)
    corine_masked[!corine_masked[] %in% corine_classes] <- NA

    if (all(is.na(corine_masked[]))) stop("All CORINE values removed after filtering.")

    # Extract final classes
    freq_table <- terra::freq(corine_masked)
    vals <- terra::values(corine_masked)
    corine_classes <- sort(unique(vals[!is.na(vals)]))
  } else {
    corine_classes <- NULL
  }

  # === Cargar ecorregiones si se requiere ===
  if (!is.null(ecoregion_shapefile_path)) {
    message("Leyendo shapefile de ecorregiones...")
    ecoregions_sf <- sf::st_read(ecoregion_shapefile_path, quiet = TRUE)

    # Recorte por la mascara de peninsula si esta definida
    if (!is.null(peninsula_shapefile)) {
      message("Recortando ecorregiones por mascara de peninsula...")
      peninsula_sf <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula_sf <- sf::st_transform(peninsula_sf, sf::st_crs(ecoregions_sf))
      ecoregions_sf <- sf::st_intersection(ecoregions_sf, peninsula_sf)
    }

    # Validar que haya geometrias y transformar al CRS del CORINE
    if (nrow(ecoregions_sf) == 0) stop("Ecoregions shapefile no contiene geometrias.")
    ecoregions_sf <- sf::st_make_valid(ecoregions_sf)
    ecoregions_sf <- sf::st_transform(ecoregions_sf, sf::st_crs(corine_vect))
    message("Ecoregiones cargadas correctamente.")

    # Si no se especifica el campo, intentar usar 'EnS_name' por defecto
    if (is.null(ecoregion_field)) {
      if ("EnS_name" %in% names(ecoregions_sf)) {
        ecoregion_field <- "EnS_name"
        message("No 'ecoregion_field' specified. Using default field: 'EnS_name'")
      } else {
        stop("You must specify 'ecoregion_field'. The default field 'EnS_name' was not found in the shapefile. Available fields: ",
             paste(names(ecoregions_sf), collapse = ", "))
      }
    }

    # Validar que el campo proporcionado existe
    if (!(ecoregion_field %in% names(ecoregions_sf))) {
      stop(paste0("Field '", ecoregion_field, "' not found in the ecoregion shapefile. Available fields: ",
                  paste(names(ecoregions_sf), collapse = ", ")))
    }

    # Extraer clases unicas si no se han definido
    if (!exists("ecoregion_classes") || is.null(ecoregion_classes)) {
      ecoregion_classes <- unique(ecoregions_sf[[ecoregion_field]])
      message("ecoregion_classes no estaba definido; usando todas las clases unicas del campo ", ecoregion_field)
    }

    ecoregion_classes_internal <- ecoregion_classes

  if (is.null(ecoregion_classes)) {
    ecoregion_classes_internal <- unique(ecoregions_sf[[ecoregion_field]])
    ecoregion_classes <- ecoregion_classes_internal  # ? importante para activar el bloque
  } else {
    ecoregion_classes_internal <- ecoregion_classes
  }
  }


  process_single_threshold <- function(r_input,
                                       pmin = NULL,
                                       pmax = NULL,
                                       otsu_min = NULL,  # umbral de entrada para filtrar el raster antes de calcular el Otsu
                                       label_suffix = "",
                                       corine_class = NULL,
                                       corine_year = NULL,
                                       ecoregion_name = NULL,
                                       ecoregion_field = NULL,
                                       cor_eco_name = NULL,           # <- NUEVO
                                       cor_eco_field = NULL,          # <- NUEVO
                                       min_otsu_threshold_value = NULL,
                                       output_dir,
                                       year = NULL,
                                       tile = TRUE,
                                       n_rows = 2,
                                       n_cols = 3,
                                       tile_overlap = 1000,
                                       python_exe,
                                       gdal_polygonize_script,
                                       output_format = c("shp", "geojson"))
  {


    # Validar y crear output_dir
    if (missing(output_dir) || is.null(output_dir)) {
      stop("'output_dir' must be provided.")
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    tile_dir <- file.path(output_dir, "temp_tiles")
    if (!dir.exists(tile_dir)) dir.create(tile_dir, recursive = TRUE)

    figures_dir <- file.path(output_dir, "FIGURES")
    if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)


    # Derivar 'year' si no esta definido
    if (is.null(year)) {
      year_match <- stringr::str_extract(basename(output_dir), "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    }

    # Derivar label base
    label <- if (label_suffix != "") paste0(label_suffix) else "default"

    # Preparar nombre de archivo raster de salida
    message("Building raster name with: year = ", year, ", label = ", label, ", output_dir = ", output_dir)
    raster_name <- file.path(output_dir, paste0("raster_", year, "_otsu_", label, ".tif"))

    # === Preparar raster filtrado y etiqueta ===

    if (!is.null(otsu_min)) {
      # Filtrado por umbral minimo
      r_filtered <- terra::ifel(r_input >= otsu_min, r_input, NA)
      label <- if (nzchar(label_suffix)) label_suffix else paste0("ge", otsu_min)

    } else if (!is.null(pmin) && !is.null(pmax)) {
      # Filtrado por percentiles
      vals_all <- terra::values(r_input)
      vals_all <- vals_all[!is.na(vals_all)]

      min_val <- quantile(vals_all, probs = pmin, na.rm = TRUE)
      max_val <- quantile(vals_all, probs = pmax, na.rm = TRUE)

      r_filtered <- terra::ifel(r_input >= min_val & r_input <= max_val, r_input, NA)
      if (nzchar(label_suffix)) {
        label <- label_suffix
      } else {
        label <- paste0("P", formatC(pmin * 1000, width = 3, flag = "0"),
                        "toP", formatC(pmax * 1000, width = 3, flag = "0"))
      }


      message(sprintf("Percentile range: pmin = %.3f (%.2f), pmax = %.3f (%.2f)",
                      pmin, min_val, pmax, max_val))

    } else {
      # Sin filtro
      r_filtered <- r_input
      label <- if (nzchar(label_suffix)) label_suffix else "original_values"
    }

    # === Validar y normalizar ===

    vals <- terra::values(r_filtered)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) stop("No valid values remaining after filtering.")

    min_val <- min(vals)
    max_val <- max(vals)
    range_val <- max_val - min_val
    if (range_val == 0) stop("Cannot normalize: max_val equals min_val.")

    r_rescaled <- (r_filtered - min_val) / range_val * 255
    vals_rescaled <- terra::values(r_rescaled)
    vals_rescaled <- vals_rescaled[!is.na(vals_rescaled)]
    if (length(vals_rescaled) == 0) stop("No valid values after normalization.")

    # === Histograma y umbral Otsu ===

    hist_values <- graphics::hist(vals_rescaled, breaks = 256, plot = FALSE)
    smoothed_counts <- OtsuSeg::smooth_histogram(hist_values$counts)
    smoothed_counts[is.na(smoothed_counts)] <- 0

    threshold_value_smoothed <- OtsuSeg::otsu_threshold_smoothed(smoothed_counts, hist_values$mids)
    real_threshold <- (threshold_value_smoothed / 255) * range_val + min_val

    message(sprintf("Threshold (%s): real threshold value = %.4f", label, real_threshold))

    # === Aplicar umbral minimo (si se define) ===

    min_applied <- FALSE
    if (!is.null(min_otsu_threshold_value)) {
      if (real_threshold < min_otsu_threshold_value) {
        message(sprintf(" ? Threshold %.2f less than minimum %.2f. The minimum value is applied.",
                        real_threshold, min_otsu_threshold_value))
        real_threshold <- min_otsu_threshold_value
        min_applied <- TRUE
      }
    }

    label <- gsub("_minapplied$", "", label)
    if (min_applied) {
      label <- paste0(label, "_minapplied")
    }

    # === Registrar umbral ===

    threshold_row <- data.frame(
      Label = label,
      RealThreshold = real_threshold,
      CORINE_CLASS = if (!is.null(corine_class)) corine_class else NA,
      ECOREGION_NAME = if (!is.null(ecoregion_name)) ecoregion_name else NA,
      COR_ECO_LABEL = if (!is.null(cor_eco_name)) cor_eco_name else NA
    )


    # Binarizacion
    binary_raster <- terra::ifel(r_filtered > real_threshold, 1, NA)
    raster_name <- file.path(output_dir, paste0("raster_", year, "_otsu_", label, ".tif"))
    terra::writeRaster(binary_raster, raster_name, overwrite = TRUE,
                       datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=-9999"))

   # if (file.exists(raster_name)) file.remove(raster_name)

    # Vectorizacion
    if (tile) {
      ext_r <- terra::ext(binary_raster)
      tile_width <- (ext_r[2] - ext_r[1]) / n_cols
      tile_height <- (ext_r[4] - ext_r[3]) / n_rows
      count <- 1
      tile_shapefiles <- c()
      for (i in 0:(n_rows - 1)) {
        for (j in 0:(n_cols - 1)) {
          xmin_tile <- max(ext_r[1] + j * tile_width - tile_overlap, ext_r[1])
          xmax_tile <- min(ext_r[1] + (j + 1) * tile_width + tile_overlap, ext_r[2])
          ymin_tile <- max(ext_r[3] + i * tile_height - tile_overlap, ext_r[3])
          ymax_tile <- min(ext_r[3] + (i + 1) * tile_height + tile_overlap, ext_r[4])
          tile_crop <- terra::crop(binary_raster, terra::ext(xmin_tile, xmax_tile, ymin_tile, ymax_tile))
          tile_path <- file.path(tile_dir, sprintf("tile_%s_%d.tif", label, count))
          terra::writeRaster(tile_crop, tile_path, overwrite = TRUE,
                             datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
          shp_path <- file.path(output_dir, sprintf("tile_%s_%d.shp", label, count))
          system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tile_path}" -f "ESRI Shapefile" "{shp_path}" DN'))
          if (file.exists(tile_path)) file.remove(tile_path)
          tile_shapefiles <- c(tile_shapefiles, shp_path)
          count <- count + 1
        }
      }
      polys <- do.call(rbind, lapply(tile_shapefiles, sf::st_read, quiet = TRUE))
      polys <- polys[!duplicated(sf::st_as_binary(sf::st_geometry(polys))), ]
      polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON", "MULTIPOLYGON"), ]
      polys <- sf::st_transform(polys, 3035)

      for (shp in tile_shapefiles) {
        shp_base <- tools::file_path_sans_ext(shp)
        extensions <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        for (ext in extensions) {
          f <- paste0(shp_base, ext)
          if (file.exists(f)) file.remove(f)
        }
      }
    } else {
      temp_raster <- file.path(tile_dir, sprintf("binary_%s.tif", label))
      terra::writeRaster(binary_raster, temp_raster, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
      shp_path <- file.path(output_dir, sprintf("burned_areas_%s_otsu_%s.shp", year, label))
      system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{temp_raster}" -f "ESRI Shapefile" "{shp_path}" DN'))
      if (file.exists(temp_raster)) file.remove(temp_raster)
      polys <- sf::st_read(shp_path, quiet = TRUE)
      polys <- sf::st_transform(polys, 3035)
    }

    if (!is.null(corine_class)) {
      polys$CORINE_CLASS <- corine_class
    }
    if (!is.null(ecoregion_name) && !is.null(ecoregion_field)) {
      polys[[ecoregion_field]] <- ecoregion_name
    }

    if (!is.null(corine_year)) {
      polys$CORINE_YEAR <- corine_year
    }

    if (!is.null(cor_eco_name) && !is.null(cor_eco_field)) {
      if (cor_eco_field %in% names(polys)) {
        warning("Overwriting existing field: ", cor_eco_field)
      }
      polys[[cor_eco_field]] <- cor_eco_name
    }


    # Histograma y curva de varianza
    interclass_variance_curve <- sapply(1:(length(hist_values$counts) - 1), function(t) {
      w_b <- sum(hist_values$counts[1:t]) / sum(hist_values$counts)
      w_f <- 1 - w_b
      if (w_b == 0 || w_f == 0) return(NA)
      mu_b <- sum(hist_values$mids[1:t] * hist_values$counts[1:t]) / sum(hist_values$counts[1:t])
      mu_f <- sum(hist_values$mids[(t + 1):length(hist_values$mids)] * hist_values$counts[(t + 1):length(hist_values$mids)]) / sum(hist_values$counts[(t + 1):length(hist_values$mids)])
      w_b * w_f * (mu_b - mu_f)^2
    })

    figures_dir <- file.path(output_dir, "FIGURES")
    if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

    plot_path <- file.path(figures_dir, paste0("otsu_plot_", year, "_", label, ".png"))

    png(plot_path, width = 1600, height = 800, res = 150)
    par(mfrow = c(1, 2))
    par(mar = c(5, 4, 4, 5) + 0.1)
    plot(hist_values$mids, smoothed_counts, type = "h", col = "grey", lwd = 2,
         xlab = "Intensity", ylab = "Frequency", main = "")
    par(new = TRUE)
    plot(hist_values$mids[-1], interclass_variance_curve, type = "l", col = "blue", lwd = 2,
         xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
    axis(side = 4, col.axis = "blue", col = "blue", lwd = 1)
    mtext("Inter-Class Variance", side = 4, line = 3, col = "blue")
    abline(v = threshold_value_smoothed, col = "red", lty = 2)
    legend("topright",
           legend = c("Histogram", "Inter-Class", paste("Real Threshold:", round(real_threshold, 2))),
           col = c("grey", "blue", "red"), lty = c(1, 1, 2), bty = "n")
    par(mar = c(5, 4, 4, 2) + 0.1)
    plot(hist_values$mids[-1], interclass_variance_curve, type = "l", col = "blue", lwd = 2,
         xlab = "Threshold", ylab = "Variance")
    dev.off()

    # === Limpiar shapefile DN.* si existe ===
    dn_base <- file.path(output_dir, "DN")
    dn_extensions <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
    for (ext in dn_extensions) {
      f <- paste0(dn_base, ext)
      if (file.exists(f)) file.remove(f)
    }

    # === Limpiar raster intermedio principal ===
    if (file.exists(raster_name)) {
      message("Deleting intermediate raster: ", raster_name)
      file.remove(raster_name)
    }

    # Tambien eliminar el raster temporal si no es tileado
    if (!tile && exists("temp_raster") && file.exists(temp_raster)) {
      message("Deleting temporary raster: ", temp_raster)
      file.remove(temp_raster)
    }

    return(list(polys = polys, log = threshold_row, tile_dir = tile_dir))
  }

  if (segment_by_intersection && !is.null(corine_classes) && !is.null(ecoregion_classes)) {

    message("Extracting CORINE year from filename...")
    corine_filename <- basename(corine_raster_path)
    corine_year <- stringr::str_extract(corine_filename, "\\d{4}")
    if (is.na(corine_year)) {
      corine_year_short <- stringr::str_extract(corine_filename, "corine(\\d{2})")
      if (!is.na(corine_year_short)) {
        corine_year <- paste0("19", stringr::str_extract(corine_year_short, "\\d{2}"))
      }
    }
    if (is.na(corine_year)) {
      warning("Could not extract CORINE year from filename.")
    } else {
      message("Extracted CORINE year: ", corine_year)
    }

    corine_vect$CORINE_CLASS <- as.character(corine_vect$CORINE_CLASS)
    corine_vect$CORINE_YEAR <- corine_year
    corine_vect$ID <- seq_len(nrow(corine_vect))

    # Reproject ecoregions
    ecoregions_sf <- sf::st_transform(ecoregions_sf, sf::st_crs(corine_vect))

    message("Intersecting CORINE and ECOREGIONS...")

    # 1. Validar geometrias rapidamente (solo si hay errores potenciales)
    corine_vect <- sf::st_make_valid(corine_vect)
    ecoregions_sf <- sf::st_make_valid(ecoregions_sf)

    # Crear campo numerico para rasterizar
    ecoregions_sf$ECO_ID_INTERNAL <- as.numeric(as.factor(ecoregions_sf[[ecoregion_field]]))

    message("Rasterizing full ecoregions map...")
    # Rasterizar todas las ecoregiones en el dominio completo
    ecoregions_r_full <- terra::rasterize(
      terra::vect(ecoregions_sf),
      terra::rast(corine_masked),  # Asegura alineacion
      field = "ECO_ID_INTERNAL"
    )

    message("Cropping + masking ecoregions to CORINE extent...")
    # Cortar y enmascarar a la zona valida de CORINE
    ecoregions_r <- terra::mask(
      terra::crop(ecoregions_r_full, corine_masked),
      corine_masked
    )

    message("Extracting dominant ecoregion in each CORINE polygon using centroids...")

    # Usar centroides de CORINE para extraccion mas rapida
    corine_centroids <- sf::st_centroid(corine_vect)

    # Extraer valor de raster de ecorregiones en centroides
    eco_vals <- terra::extract(ecoregions_r, terra::vect(corine_centroids), ID = FALSE)

    # Anadir ID de los poligonos CORINE
    eco_vals$ID <- corine_vect$ID
    eco_vals$ECO_ID_INTERNAL <- eco_vals[[1]]

    # Asociar nombres reales de ecoregion
    eco_id_to_name <- levels(as.factor(ecoregions_sf[[ecoregion_field]]))
    eco_vals$ECO_CLASS <- eco_id_to_name[eco_vals$ECO_ID_INTERNAL]

    # Limpiar NAs
    eco_vals <- eco_vals[!is.na(eco_vals$ECO_CLASS), ]

    # Dominant ecoregion por CORINE
    dominant_eco <- eco_vals[, c("ID", "ECO_CLASS")]

    message("Joining dominant ecoregion with each CORINE polygon...")
    corine_vect <- dplyr::left_join(corine_vect, dominant_eco, by = "ID")

    if (any(is.na(corine_vect$ECO_CLASS))) {
      message("Rescuing missing ECO_CLASS values...")

      # Extraer los poligonos sin clase asignada
      missing_na <- corine_vect[is.na(corine_vect$ECO_CLASS), ]
      message("Found ", nrow(missing_na), " polygons with missing ECO_CLASS")

      # Calcular centroides
      message("Computing centroids of missing polygons...")
      centroids_na <- sf::st_centroid(missing_na)
      centroids_na_spat <- terra::vect(centroids_na)

      # Extraer valores del raster sin enmascarar (mas seguro)
      message("Extracting ECO_ID_INTERNAL from raster at centroid locations (terra::extract)...")
      vals <- terra::extract(ecoregions_r, centroids_na_spat)

      # Crear tabla de recuperacion
      recovered <- data.frame(ID = missing_na$ID, ECO_ID_INTERNAL = vals[[2]])
      recovered <- recovered[!is.na(recovered$ECO_ID_INTERNAL), ]

      # Mapear codigos internos a nombres reales
      eco_id_to_name <- levels(as.factor(ecoregions_sf[[ecoregion_field]]))
      recovered$ECO_CLASS <- eco_id_to_name[recovered$ECO_ID_INTERNAL]

      # Join con corine_vect
      if (nrow(recovered) > 0) {
        message("Rescuing ", nrow(recovered), " values via raster centroid extraction.")
        corine_vect <- dplyr::left_join(
          corine_vect,
          recovered[, c("ID", "ECO_CLASS")],
          by = "ID",
          suffix = c("", ".rescued")
        )

        # Reemplazar valores NA
        corine_vect$ECO_CLASS <- ifelse(
          is.na(corine_vect$ECO_CLASS),
          corine_vect$ECO_CLASS.rescued,
          corine_vect$ECO_CLASS
        )
        corine_vect$ECO_CLASS.rescued <- NULL
      } else {
        message("No ECO_CLASS could be rescued via raster centroids.")
      }

      # Si aun hay NAs, hacer ultimo intento con st_join vectorial (lento pero robusto)
      if (any(is.na(corine_vect$ECO_CLASS))) {
        message("Attempting final rescue with spatial join (st_join)...")

        missing_na2 <- corine_vect[is.na(corine_vect$ECO_CLASS), ]
        joined <- sf::st_join(missing_na2, ecoregions_sf[, ecoregion_field], left = FALSE)

        if (nrow(joined) > 0) {
          recovered2 <- joined[, c("ID")]
          recovered2$ECO_CLASS <- joined[[ecoregion_field]]

          corine_vect <- dplyr::left_join(
            corine_vect,
            sf::st_drop_geometry(recovered2),
            by = "ID",
            suffix = c("", ".rescued")
          )

          corine_vect$ECO_CLASS <- ifelse(
            is.na(corine_vect$ECO_CLASS),
            corine_vect$ECO_CLASS.rescued,
            corine_vect$ECO_CLASS
          )
          corine_vect$ECO_CLASS.rescued <- NULL

          message("Successfully rescued ", nrow(recovered2), " polygons via st_join().")
        } else {
          message("No ECO_CLASS could be rescued with st_join either.")
        }
      }

      # Finalmente eliminar los que aun estan vacios
      missing_final <- sum(is.na(corine_vect$ECO_CLASS))
      if (missing_final > 0) {
        warning(missing_final, " CORINE polygons had no ecoregion assigned and will be removed.")
        corine_vect <- corine_vect[!is.na(corine_vect$ECO_CLASS), ]
      } else {
        message("All CORINE polygons have been assigned an ECO_CLASS.")
      }
    }

    message("Creating combined labels ECO_CLASS..")
    # Crear etiquetas combinadas tipo "1_MDS"
    corine_vect$unit_id2 <- paste0(
      corine_vect$CORINE_CLASS, "_",
      gsub("[^A-Za-z0-9]", "", as.character(corine_vect$ECO_CLASS))
    )

    message("Ensuring individual polygons (only casting to POLYGON)...")
    intersected_units <- sf::st_cast(corine_vect, "POLYGON")

    # Guardar shapefile final
    out_name <- file.path(output_corine_vector_dir, paste0("corine_with_ecoregions_", corine_year, ".shp"))
    message("Saving CORINE + ecoregion unit shapefile: ", out_name)

    if (file.exists(out_name)) {
      unlink(list.files(dirname(out_name), pattern = tools::file_path_sans_ext(basename(out_name)), full.names = TRUE))
    }

    sf::st_write(intersected_units, out_name, delete_layer = TRUE, quiet = TRUE)


  message("Dissolving individual polygons...")

  # Agrupar por combinaciones CORINE x Ecoregion (unit_id2)

    # Paso 1: Juntar atributos primero (sobre el objeto con geometria)
    intersected_units_sf <- sf::st_as_sf(intersected_units)

    # Paso 2: Agrupar y disolver en un solo paso con sf (mas estable)
    sf_agg <- intersected_units_sf |>
      dplyr::group_by(unit_id2) |>
      dplyr::summarise(
        CORINE_CLASS = first(CORINE_CLASS),
        CORINE_YEAR  = first(CORINE_YEAR),
        ECO_CLASS    = first(ECO_CLASS),
        .groups = "drop"
      )

    # Geometry validation and cleaning for sf_agg

    # 1. Ensure object is an sf
    sf_agg <- sf::st_as_sf(sf_agg)

    # 2. Check and fix invalid geometries
    if (any(!sf::st_is_valid(sf_agg))) {
      message("Fixing invalid geometries using st_make_valid()...")
      sf_agg <- sf::st_make_valid(sf_agg)
    } else {
      message("All geometries are valid.")
    }

    # 3. Check geometry types
    geom_types <- unique(sf::st_geometry_type(sf_agg))
    message("Detected geometry types: ", paste(geom_types, collapse = ", "))

    # 4. Keep only POLYGON and MULTIPOLYGON types
    sf_agg <- sf_agg[sf::st_geometry_type(sf_agg) %in% c("POLYGON", "MULTIPOLYGON"), ]

    # 5. Convert all to MULTIPOLYGON (optional but recommended)
    sf_agg <- sf::st_cast(sf_agg, "MULTIPOLYGON")

    # 6. Final confirmation
    message("Geometry validation completed. Total features: ", nrow(sf_agg))

  units_grouped <- sf::st_as_sf(sf_agg)


  # Guardar shapefile final
  out_name1 <- file.path(output_corine_vector_dir, paste0("corine_with_ecoregions_", corine_year, "_dissolve.shp"))
  message("Saving CORINE + ecoregion unit shapefile: ", out_name1)

  if (file.exists(out_name1)) {
    unlink(list.files(dirname(out_name1), pattern = tools::file_path_sans_ext(basename(out_name1)), full.names = TRUE))
  }

  sf::st_write(units_grouped, out_name1, delete_layer = TRUE, quiet = TRUE)

  }

  if (!is.null(trim_percentiles)) {
      if (!all(c("min", "max") %in% names(trim_percentiles))) {
        stop("El data frame 'trim_percentiles' debe tener columnas llamadas 'min' y 'max'.")
      }
      if (any(trim_percentiles$min >= trim_percentiles$max)) {
        stop("Cada valor 'min' en 'trim_percentiles' debe ser menor que su correspondiente 'max'.")
      }

    # === CASO 1: Interseccion CORINE ? Ecoregiones ===
    if (segment_by_intersection && !is.null(corine_classes) && !is.null(ecoregion_classes)) {

        message("Processing percentiles by CORINE ? Ecoregions intersection...")

        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]

          threshold_log <- NULL
          all_polys <- list()

          for (j in seq_len(nrow(units_grouped))) {
            unit <- units_grouped[j, ]
            corine_class <-unit$CORINE_CLASS
            cor_eco_name <- unit$unit_id2
            ecoregion_name<- unit$ECO_CLASS
            label_suffix <- paste0(
              unit$unit_id2,
              "_P", pmin * 100, "_toP", pmax * 100
            )

            mask_r <- terra::mask(r, terra::vect(unit))
            vals   <- terra::values(mask_r)
            vals   <- vals[!is.na(vals)]


            res_output <- process_single_threshold(
              r_input = mask_r,
              pmin = pmin,
              pmax = pmax,
              label_suffix = label_suffix,
              corine_class = corine_class,
              cor_eco_name = cor_eco_name,
              ecoregion_name = ecoregion_name,
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_format = output_format,
              output_dir = output_dir,
              year = year,
              tile = tile,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script
            )

            if (!is.null(res_output$polys) && nrow(res_output$polys) > 0) {
              label_real <- res_output$log$Label  # pimero definir

              # asegurarte de que todas las columnas esten
              res_output$polys$CORINE_CLASS <- corine_class
              res_output$polys$ECO_CLASS <- ecoregion_name
              res_output$polys$Label <- label_real
              res_output$polys$unit_id2 <- unit$unit_id2
              res_output$log$unit_id2 <- unit$unit_id2

              # seleccionar columnas comunes
              res_output$polys <- dplyr::select(
                res_output$polys,
                any_of(c("DN", "geometry", "CORINE_CLASS", "ECO_CLASS", "Label", "unit_id2"))
              )

             # all_polys[[label_real]] <- res_output$polys
              all_polys[[length(all_polys) + 1]] <- res_output$polys
              threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)

            }

          }

          combined <- do.call(rbind, all_polys)

          filename_base <- sprintf("burned_areas_%s_percentiles_corine_ecoregion_P%02d_toP%02d",
                                   year, round(pmin * 100), round(pmax * 100))

          ext <- if (output_format == "geojson") ".geojson" else ".shp"
          out_file <- file.path(output_dir, paste0(filename_base, ext))

          if (output_format == "shp") {
            shp_base <- tools::file_path_sans_ext(out_file)
            for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
              f <- paste0(shp_base, ext_i)
              if (file.exists(f)) file.remove(f)
            }
          } else {
            if (file.exists(out_file)) file.remove(out_file)
          }

          sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
          message("Saved: ", out_file)

          log_file <- file.path(
            output_dir,
            sprintf("burned_areas_%s_corine_ecoregion_percentiles_P%02d_toP%02d.txt",
                    year, round(pmin * 100), round(pmax * 100))
          )
          if (file.exists(log_file)) file.remove(log_file)

          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)


          # === Limpiar carpeta temporal tile_dir ===

          tile_dir <- file.path(output_dir, "temp_tiles")
          if (dir.exists(tile_dir)) {
            message("Cleaning up tile_dir: ", tile_dir)
            tryCatch({
              unlink(tile_dir, recursive = TRUE, force = TRUE)
            }, error = function(e) {
              warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
            })
          }
          else {
            message("No polygons generated for percentile ", pmin, "-", pmax, ". Skipping output.")
          }
          message("Segmenting by CORINE ? Ecoregion intersection completed.")
        }
    }

    # === CASO 2: Percentiles por clase CORINE ===
    else if (!is.null(corine_classes)&& !segment_by_intersection) {

        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]

          all_polys <- list()
          threshold_log <- data.frame()

          for (cls in corine_classes) {
            mask_r <- terra::ifel(corine_masked == cls, r, NA)

            label_suffix <- paste0(
              "corine_", cls,
              "_P", pmin * 100, "_toP", pmax * 100
            )

            res_output <- process_single_threshold(
              r_input = mask_r,
              pmin = pmin,
              pmax = pmax,
              label_suffix = label_suffix,
              corine_class = cls,
              corine_year = corine_year,
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_dir = output_dir,
              year = year,
              tile = tile,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script,
              output_format = output_format
            )

            res_output$polys$CORINE_CLASS <- cls
            res_output$polys$Label <- res_output$log$Label
            all_polys[[length(all_polys) + 1]] <- res_output$polys # nuevo
            threshold_log <- rbind(threshold_log, res_output$log)
          }

          combined <- do.call(rbind, all_polys)
          combined <- sf::st_make_valid(combined)
          names(combined) <- make.unique(substr(names(combined), 1, 10))

    ###########  SAVING #######

          # Construir nombre base del archivo
          filename_base <- sprintf("burned_areas_%s_percentiles_corine_P%02d_toP%02d",
                                   year, round(pmin * 100), round(pmax * 100))

          # Definir la extension segun el formato
          ext <- if (output_format == "geojson") ".geojson" else ".shp"
          out_file <- file.path(output_dir, paste0(filename_base, ext))

          # Eliminar archivo(s) previos si existen
          if (output_format == "shp") {
            shp_base <- tools::file_path_sans_ext(out_file)
            for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
              f <- paste0(shp_base, ext_i)
              if (file.exists(f)) file.remove(f)
            }
          } else {
            if (file.exists(out_file)) file.remove(out_file)
          }

          # Guardar el shapefile o GeoJSON
          sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)

          # Mensaje opcional
          message("Saved: ", out_file)

          log_file <- file.path(
            output_dir,
            sprintf("burned_areas_%s_corine_percentiles_P%02d_toP%02d.txt",
                    year, round(pmin * 100), round(pmax * 100))
          )

          if (file.exists(log_file)) file.remove(log_file)

          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

          # === Limpiar carpeta temporal tile_dir ===
          tile_dir <- file.path(output_dir, "temp_tiles")
          if (dir.exists(tile_dir)) {
            message("Cleaning up tile_dir: ", tile_dir)
            tryCatch({
              unlink(tile_dir, recursive = TRUE, force = TRUE)
            }, error = function(e) {
              warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
            })
          }


        }
    }

    # === CASO 3: Percentiles por clase de ecoregion ===
    else if (!is.null(ecoregion_classes)&& !segment_by_intersection) {
        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]

          threshold_log <- data.frame()
          all_polys <- list()

          for (eco in ecoregion_classes) {
            eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
            eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
            clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)

            label_suffix <- paste0(
              "ecoregion_", clean_eco,
              "_P", pmin * 100, "_toP", pmax * 100
            )

            res_output <- process_single_threshold(
              r_input = eco_raster_mask,
              pmin = pmin,
              pmax = pmax,
              label_suffix = label_suffix,
              ecoregion_name = eco,
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_dir = output_dir,
              year = year,
              tile = tile,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script,
              output_format = output_format
            )

            label_real <- res_output$log$Label
            res_output$polys$ECO_CLASS <- eco
            res_output$polys$Label <- label_real
            #all_polys[[label_real]] <- res_output$polys
            all_polys[[length(all_polys) + 1]] <- res_output$polys

            threshold_log <- rbind(threshold_log, res_output$log)
          }

          combined <- do.call(rbind, all_polys)
          combined <- sf::st_make_valid(combined)
          names(combined) <- make.unique(substr(names(combined), 1, 10))

          # Construir nombre base del archivo
          filename_base <- sprintf("burned_areas_%s_percentiles_ecoregion_P%02d_toP%02d",
                                   year, round(pmin * 100), round(pmax * 100))


          # Definir la extension segun el formato
          ext <- if (output_format == "geojson") ".geojson" else ".shp"
          out_file <- file.path(output_dir, paste0(filename_base, ext))

          # Eliminar archivo(s) previos si existen
          if (output_format == "shp") {
            shp_base <- tools::file_path_sans_ext(out_file)
            for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
              f <- paste0(shp_base, ext_i)
              if (file.exists(f)) file.remove(f)
            }
          } else {
            if (file.exists(out_file)) file.remove(out_file)
          }

          # Guardar el shapefile o GeoJSON
          sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)

          log_file <- file.path(
            output_dir,
            sprintf("burned_areas_%s_ecoregion_percentiles_P%02d_toP%02d.txt",
                    year, round(pmin * 100), round(pmax * 100))
          )

          if (file.exists(log_file)) file.remove(log_file)

          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

          # === Limpiar carpeta temporal tile_dir ===
          tile_dir <- file.path(output_dir, "temp_tiles")
          if (dir.exists(tile_dir)) {
            message("Cleaning up tile_dir: ", tile_dir)
            tryCatch({
              unlink(tile_dir, recursive = TRUE, force = TRUE)
            }, error = function(e) {
              warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
            })
          }

          # Mensaje opcional
          message("Saved: ", out_file)
        }
    }

    # === CASO 4: Percentiles sobre raster global ===
    else {

        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]

          threshold_log <- data.frame()

          res_output <- process_single_threshold(
            r_input = r,
            pmin = pmin,
            pmax = pmax,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            output_format = output_format
          )


          label_real <- res_output$log$Label
          threshold_log <- rbind(threshold_log, res_output$log)

          combined <- sf::st_make_valid(res_output$polys)


          # Construir nombre base del archivo
          filename_base <- sprintf("burned_areas_%s_percentiles_original_P%02d_toP%02d",
                                   year, round(pmin * 100), round(pmax * 100))


          # Definir la extension segun el formato
          ext <- if (output_format == "geojson") ".geojson" else ".shp"
          out_file <- file.path(output_dir, paste0(filename_base, ext))

          # Eliminar archivo(s) previos si existen
          if (output_format == "shp") {
            shp_base <- tools::file_path_sans_ext(out_file)
            for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
              f <- paste0(shp_base, ext_i)
              if (file.exists(f)) file.remove(f)
            }
          } else {
            if (file.exists(out_file)) file.remove(out_file)
          }

          # Guardar el shapefile o GeoJSON
          sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)

          # Mensaje opcional
          message("Saved: ", out_file)


          log_file <- file.path(
            output_dir,
            sprintf("burned_areas_%s_percentiles_P%02d_toP%02d.txt",
                                year, round(pmin * 100), round(pmax * 100))
            )

          if (file.exists(log_file)) file.remove(log_file)

          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

          # === Limpiar carpeta temporal tile_dir ===
          tile_dir <- file.path(output_dir, "temp_tiles")
          if (dir.exists(tile_dir)) {
            message("Cleaning up tile_dir: ", tile_dir)
            tryCatch({
              unlink(tile_dir, recursive = TRUE, force = TRUE)
            }, error = function(e) {
              warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
            })
          }
        }
    }

      }

  else if (!is.null(corine_classes) && !segment_by_intersection) {
    if (use_original) {

      threshold_log <- data.frame()
      all_polys <- list()

      for (cls in corine_classes) {
        mask_r <- terra::ifel(corine_masked == cls, r, NA)
        label_suffix <- paste0("corine_class_", cls, "_original")

        res_output <- process_single_threshold(
          r_input = mask_r,
          otsu_min = NULL,
          label_suffix = label_suffix,
          corine_class = cls,
          corine_year = corine_year,
          min_otsu_threshold_value = min_otsu_threshold_value,
          output_dir = output_dir,
          year = year,
          tile = tile,
          python_exe = python_exe,
          gdal_polygonize_script = gdal_polygonize_script,
          output_format = output_format
        )


        label_real <- res_output$log$Label
        res_output$polys$CORINE_CLASS <- cls
        res_output$polys$Label <- label_real
        #all_polys[[label_real]] <- res_output$polys
        all_polys[[length(all_polys) + 1]] <- res_output$polys
        threshold_log <- rbind(threshold_log, res_output$log)
      }

      combined <- do.call(rbind, all_polys)
      combined <- sf::st_make_valid(combined)
      names(combined) <- make.unique(substr(names(combined), 1, 10))


      # Definir extension segun el formato
      ext <- if (output_format == "geojson") ".geojson" else ".shp"

      # Crear nombre de archivo con extension adecuada
      out_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_original%s", year, ext))

      # Eliminar archivo existente segun formato
      if (output_format == "shp") {
        shp_base <- tools::file_path_sans_ext(out_file)
        shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        for (ext_i in shp_exts) {
          f <- paste0(shp_base, ext_i)
          if (file.exists(f)) file.remove(f)
        }
      } else {
        if (file.exists(out_file)) file.remove(out_file)
      }

      # Guardar salida
      sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)

      # Mensaje opcional
      message("Saved: ", out_file)

      log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_original_log.txt", year, otsu_min))

      if (file.exists(log_file)) file.remove(log_file)

      write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      # === Limpiar carpeta temporal tile_dir ===
      tile_dir <- file.path(output_dir, "temp_tiles")
      if (dir.exists(tile_dir)) {
        message("Cleaning up tile_dir: ", tile_dir)
        tryCatch({
          unlink(tile_dir, recursive = TRUE, force = TRUE)
        }, error = function(e) {
          warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
        })
      }


    } else {
      for (otsu_min in otsu_thresholds) {

        all_polys <- list()
        threshold_log <- data.frame()


        for (cls in corine_classes) {
          mask_r <- terra::ifel(corine_masked == cls, r, NA)
          clean_cls <- gsub("[^a-zA-Z0-9]", "_", as.character(cls))
          label_suffix <- paste0("corine_", clean_cls, "_ge", otsu_min)

          res_output <- process_single_threshold(
            r_input = mask_r,
            otsu_min = otsu_min,
            label_suffix = label_suffix,
            corine_class = cls,
            corine_year = corine_year,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            output_format = output_format
          )


          label_real <- res_output$log$Label
          res_output$polys$CORINE_CLASS <- cls
          res_output$polys$Label <- label_real
          #all_polys[[label_real]] <- res_output$polys
          all_polys[[length(all_polys) + 1]] <- res_output$polys
          threshold_log <- rbind(threshold_log, res_output$log)
        }

        combined <- do.call(rbind, all_polys)
        combined <- sf::st_make_valid(combined)
        names(combined) <- make.unique(substr(names(combined), 1, 10))

        # Validar formato de salida
        output_format <- match.arg(output_format, choices = c("shp", "geojson"))

        # Determinar extension segun el formato
        ext <- if (output_format == "geojson") ".geojson" else ".shp"

        # Construir nombre de archivo con extension correcta
        out_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_ge%d%s", year, otsu_min, ext))

        # Eliminar archivos previos si existen
        if (output_format == "shp") {
          shp_base <- tools::file_path_sans_ext(out_file)
          shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
          for (ext_i in shp_exts) {
            f <- paste0(shp_base, ext_i)
            if (file.exists(f)) file.remove(f)
          }
        } else {
          if (file.exists(out_file)) file.remove(out_file)
        }

        # Guardar archivo
        sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)

        # Mensaje de confirmacion
        message("Saved: ", out_file)

        log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_ge%d_log.txt", year, otsu_min))

        if (file.exists(log_file)) file.remove(log_file)

        write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

        # === Limpiar carpeta temporal tile_dir ===
        tile_dir <- file.path(output_dir, "temp_tiles")
        if (dir.exists(tile_dir)) {
          message("Cleaning up tile_dir: ", tile_dir)
          tryCatch({
            unlink(tile_dir, recursive = TRUE, force = TRUE)
          }, error = function(e) {
            warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
          })
        }

      }
    }
  }

  else if (!is.null(ecoregion_classes) && !segment_by_intersection) {
    message("Processing Otsu by ecoregion classes...")

    if (use_original) {

      all_polys <- list()
      threshold_log <- data.frame()

      for (eco in ecoregion_classes) {
        eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
        eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
        clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)
        label_suffix <- paste0("ecoregion_", clean_eco, "_original")

        res_output <- process_single_threshold(
          r_input = eco_raster_mask,
          otsu_min = NULL,
          label_suffix = label_suffix,
          ecoregion_name = eco,
          output_format = output_format,
          output_dir = output_dir,
          year = year,
          tile = tile,
          python_exe = python_exe,
          gdal_polygonize_script = gdal_polygonize_script
        )


        label_real <- res_output$log$Label
        res_output$polys[[ecoregion_field]] <- eco  # <- Dinamico
        res_output$polys$Label <- label_real
        #all_polys[[label_real]] <- res_output$polys
        all_polys[[length(all_polys) + 1]] <- res_output$polys
        threshold_log <- rbind(threshold_log, res_output$log)
      }

      combined <- do.call(rbind, all_polys)
      combined <- sf::st_make_valid(combined)
      names(combined) <- make.unique(substr(names(combined), 1, 10))

      ext <- if (output_format == "geojson") ".geojson" else ".shp"
      out_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_ecoregion_original%s", year, ext))

      if (output_format == "shp") {
        shp_base <- tools::file_path_sans_ext(out_file)
        for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
          f <- paste0(shp_base, ext_i)
          if (file.exists(f)) file.remove(f)
        }
      } else {
        if (file.exists(out_file)) file.remove(out_file)
      }

      sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
      message("Saved: ", out_file)


      log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_ecoregion_original_log.txt", year))

      if (file.exists(log_file)) file.remove(log_file)

      write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      # === Limpiar carpeta temporal tile_dir ===
      tile_dir <- file.path(output_dir, "temp_tiles")
      if (dir.exists(tile_dir)) {
        message("Cleaning up tile_dir: ", tile_dir)
        tryCatch({
          unlink(tile_dir, recursive = TRUE, force = TRUE)
        }, error = function(e) {
          warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
        })
      }

    } else {
      for (otsu_min in otsu_thresholds) {

        threshold_log <- data.frame()
        all_polys <- list()

        for (eco in ecoregion_classes) {
          eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
          eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
          clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)
          label_suffix <- paste0("ecoregion_", clean_eco, "_ge", otsu_min)

          res_output <- process_single_threshold(
            r_input = eco_raster_mask,
            otsu_min = otsu_min,
            label_suffix = label_suffix,
            ecoregion_name = eco,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_format = output_format,
            output_dir = output_dir,
            year = year,
            tile = tile,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script
          )


          label_real <- res_output$log$Label
          res_output$polys[[ecoregion_field]] <- eco  # <- Dinamico
          res_output$polys$Label <- label_real
          #all_polys[[label_real]] <- res_output$polys
          all_polys[[length(all_polys) + 1]] <- res_output$polys
          threshold_log <- rbind(threshold_log, res_output$log)
        }

        combined <- do.call(rbind, all_polys)
        combined <- sf::st_make_valid(combined)
        names(combined) <- make.unique(substr(names(combined), 1, 10))

        filename <- sprintf("burned_areas_%s_otsu_ecoregion_ge%d", year, otsu_min)

        ext <- if (output_format == "geojson") ".geojson" else ".shp"
        out_file <- file.path(output_dir, paste0(filename, ext))

        if (output_format == "shp") {
          shp_base <- tools::file_path_sans_ext(out_file)
          for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
            f <- paste0(shp_base, ext_i)
            if (file.exists(f)) file.remove(f)
          }
        } else {
          if (file.exists(out_file)) file.remove(out_file)
        }

        sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
        message("Saved: ", out_file)

        log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_ecoregion_ge%d_log.txt", year, otsu_min))

        if (file.exists(log_file)) file.remove(log_file)
        write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

        # === Limpiar carpeta temporal tile_dir ===
        tile_dir <- file.path(output_dir, "temp_tiles")
        if (dir.exists(tile_dir)) {
          message("Cleaning up tile_dir: ", tile_dir)
          tryCatch({
            unlink(tile_dir, recursive = TRUE, force = TRUE)
          }, error = function(e) {
            warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
          })
        }

      }
    }
  }

  else if (segment_by_intersection && !is.null(corine_classes) && !is.null(ecoregion_classes)) {

    message("Processing ", if (use_original) "original values" else "Otsu thresholding", " by CORINE ? Ecoregions intersection...")

    if (use_original) {
      all_polys <- list()
      threshold_log <- NULL

      for (j in seq_len(nrow(units_grouped))) {
        unit <- units_grouped[j, ]
        corine_class   <- unit$CORINE_CLASS
        cor_eco_name   <- unit$unit_id2
        ecoregion_name <- unit$ECO_CLASS
        label_suffix   <- paste0(cor_eco_name, "_original")

        mask_r <- terra::mask(r, terra::vect(unit))
        vals   <- terra::values(mask_r)
        vals   <- vals[!is.na(vals)]

        res_output <- process_single_threshold(
          r_input = mask_r,
          otsu_min = NULL,
          label_suffix = label_suffix,
          corine_class = corine_class,
          ecoregion_name = ecoregion_name,
          cor_eco_name = cor_eco_name,
          cor_eco_field = "COR_ECO_LABEL",  # o el nombre que definas arriba
          min_otsu_threshold_value = min_otsu_threshold_value,
          output_dir = output_dir,
          output_format = output_format,
          year = year,
          tile = tile,
          python_exe = python_exe,
          gdal_polygonize_script = gdal_polygonize_script
        )

        label_real <- res_output$log$Label
        res_output$polys$CORINE_CLASS <- corine_class
        res_output$polys$ECO_CLASS    <- ecoregion_name
        res_output$polys$unit_id2     <- cor_eco_name
        res_output$polys$Label        <- label_real
        res_output$polys$COR_ECO_LABEL <- cor_eco_name  # campo fijo y valido

       # all_polys[[label_real]] <- res_output$polys
        all_polys[[length(all_polys) + 1]] <- res_output$polys
        threshold_log <- rbind(threshold_log, res_output$log)
      }


      combined <- do.call(rbind, all_polys)
      combined <- sf::st_make_valid(combined)
      names(combined) <- make.unique(substr(names(combined), 1, 10))

      filename <- sprintf("burned_areas_%s_original_corine_ecoregion.shp", year)
      out_file <- file.path(output_dir, filename)

      shp_base <- tools::file_path_sans_ext(out_file)
      for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }

      sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
      message("Saved: ", out_file)

      log_file <- file.path(output_dir, sprintf("burned_areas_%s_original_corine_ecoregion_log.txt", year))
      if (file.exists(log_file)) file.remove(log_file)
      write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      # === Limpiar carpeta temporal tile_dir ===
      tile_dir <- file.path(output_dir, "temp_tiles")
      if (dir.exists(tile_dir)) {
        message("Cleaning up tile_dir: ", tile_dir)
        tryCatch({
          unlink(tile_dir, recursive = TRUE, force = TRUE)
        }, error = function(e) {
          warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
        })
      }

    } else {
      for (otsu_min in otsu_thresholds) {
        all_polys <- list()
        threshold_log <- NULL

        for (j in seq_len(nrow(units_grouped))) {
          unit <- units_grouped[j, ]
          label_suffix    <- paste0(unit$unit_id2, "_ge", otsu_min)
          corine_class    <- unit$CORINE_CLASS
          ecoregion_name  <- unit$ECO_CLASS
          cor_eco_name    <- unit$unit_id2  # combinacion CORINE ? ECOREGION

          mask_r <- terra::mask(r, terra::vect(unit))
          vals   <- terra::values(mask_r)
          vals   <- vals[!is.na(vals)]

          res_output <- process_single_threshold(
            r_input = mask_r,
            otsu_min = otsu_min,
            label_suffix = label_suffix,
            corine_class = corine_class,
            ecoregion_name = ecoregion_name,
            cor_eco_name = cor_eco_name,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            output_format = output_format,
            year = year,
            tile = tile,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script
          )

          label_real <- res_output$log$Label
          res_output$polys$CORINE_CLASS <- corine_class
          res_output$polys[[ecoregion_field]] <- ecoregion_name
          res_output$polys$unit_id2 <- cor_eco_name
          res_output$polys$Label <- label_real

          #all_polys[[label_real]] <- res_output$polys
          all_polys[[length(all_polys) + 1]] <- res_output$polys
          threshold_log <- rbind(threshold_log, res_output$log)
        }

        combined <- do.call(rbind, all_polys)
        combined <- sf::st_make_valid(combined)
        names(combined) <- make.unique(substr(names(combined), 1, 10))

        # Definir extension segun el formato
        ext <- if (output_format == "geojson") ".geojson" else ".shp"

        # Crear nombre de archivo con extension adecuada
        out_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_ecoregion_ge%d%s", year, otsu_min, ext))

        # Eliminar archivo existente segun formato
        if (output_format == "shp") {
          shp_base <- tools::file_path_sans_ext(out_file)
          shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
          for (ext_i in shp_exts) {
            f <- paste0(shp_base, ext_i)
            if (file.exists(f)) file.remove(f)
          }
        } else {
          if (file.exists(out_file)) file.remove(out_file)
        }

        sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
        message("Saved: ", out_file)

        log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_ecoregion_ge%d_log.txt", year, otsu_min))
        if (file.exists(log_file)) file.remove(log_file)
        write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

        # === Limpiar carpeta temporal tile_dir ===
        tile_dir <- file.path(output_dir, "temp_tiles")
        if (dir.exists(tile_dir)) {
          message("Cleaning up tile_dir: ", tile_dir)
          tryCatch({
            unlink(tile_dir, recursive = TRUE, force = TRUE)
          }, error = function(e) {
            warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
          })
        }
      }
    }
  }

  else if (use_original) {
    threshold_log <- data.frame()
    label_suffix <- "original_values"

    res_output <- process_single_threshold(
      r_input = r,
      label_suffix = label_suffix,
      min_otsu_threshold_value = min_otsu_threshold_value,  # se puede aplicar si es relevante
      output_format = output_format,
      output_dir = output_dir,
      year = year,
      tile = tile,
      python_exe = python_exe,
      gdal_polygonize_script = gdal_polygonize_script
    )

    label_real <- res_output$log$Label
    threshold_log <- rbind(threshold_log, res_output$log)
    combined <- sf::st_make_valid(res_output$polys)

    ext <- if (output_format == "geojson") ".geojson" else ".shp"
    out_file <- file.path(output_dir, sprintf("burned_areas_%s_%s%s", year, label_real, ext))

    if (output_format == "shp") {
      shp_base <- tools::file_path_sans_ext(out_file)
      for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }
    } else {
      if (file.exists(out_file)) file.remove(out_file)
    }

    sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
    message("Saved: ", out_file)

    log_file <- file.path(output_dir, sprintf("burned_areas_%s_%s_log.txt", year, label_real))
    if (file.exists(log_file)) file.remove(log_file)
    write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

    tile_dir <- file.path(output_dir, "temp_tiles")
    if (dir.exists(tile_dir)) {
      message("Cleaning up tile_dir: ", tile_dir)
      tryCatch({
        unlink(tile_dir, recursive = TRUE, force = TRUE)
      }, error = function(e) {
        warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
      })
    }

  }

  # === CASO 3: Umbralizacion clasica con valores calculados (Otsu puro) ===
  else {
    for (otsu_min in otsu_thresholds) {
      threshold_log <- data.frame()
      label_suffix <- paste0("ge", otsu_min)

      res_output <- process_single_threshold(
        r_input = r,
        otsu_min = otsu_min,
        label_suffix = label_suffix,
        min_otsu_threshold_value = min_otsu_threshold_value,
        output_format = output_format,
        output_dir = output_dir,
        year = year,
        tile = tile,
        python_exe = python_exe,
        gdal_polygonize_script = gdal_polygonize_script
      )

      label_real <- res_output$log$Label
      threshold_log <- rbind(threshold_log, res_output$log)
      combined <- sf::st_make_valid(res_output$polys)

      ext <- if (output_format == "geojson") ".geojson" else ".shp"
      out_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_%s%s", year, label_real, ext))

      if (output_format == "shp") {
        shp_base <- tools::file_path_sans_ext(out_file)
        for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
          f <- paste0(shp_base, ext_i)
          if (file.exists(f)) file.remove(f)
        }
      } else {
        if (file.exists(out_file)) file.remove(out_file)
      }

      sf::st_write(combined, out_file, append = FALSE, quiet = TRUE)
      message("Saved: ", out_file)

      log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_%s_log.txt", year, label_real))
      if (file.exists(log_file)) file.remove(log_file)
      write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      tile_dir <- file.path(output_dir, "temp_tiles")
      if (dir.exists(tile_dir)) {
        message("Cleaning up tile_dir: ", tile_dir)
        tryCatch({
          unlink(tile_dir, recursive = TRUE, force = TRUE)
        }, error = function(e) {
          warning("Could not delete tile_dir: ", tile_dir, "\nReason: ", e$message)
        })
      }

    }
  }


  return(results)
}

