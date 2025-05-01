#' Validate burned area maps against reference polygons
#'
#' @description
#' The `validate_fire_maps()` function evaluates the spatial accuracy of burned area detection shapefiles by comparing them against independent reference fire polygons (e.g., from Focclim or national databases).
#'
#' It calculates both **pixel-based metrics** and **area-based metrics**. Reference polygons are masked to retain only burnable areas based on a CORINE-derived raster and filtered by a minimum area threshold if specified.
#'
#' **Important**: If `min_area_reference_ha` is changed, you must set `force_reprocess_ref = TRUE` to recalculate the masked reference polygons.
#'
#' @param input_shapefile Character vector. One or several paths to shapefiles containing the burned polygons to validate.
#' @param ref_shapefile Character. Path to the shapefile with reference burned area polygons.
#' @param mask_shapefile Character. Path to the shapefile defining the study area boundary.
#' @param burnable_raster Character. Path to the raster file indicating burnable areas (binary or categorical).
#' @param year_target Numeric. Target year for filtering reference polygons.
#' @param validation_dir Character. Output directory to save validation results.
#' @param binary_burnable Logical. TRUE if the burnable raster is binary (1 = burnable, 0 = non-burnable). Default is TRUE.
#' @param burnable_classes Optional vector of raster values considered burnable if `binary_burnable = FALSE`.
#' @param buffer Numeric. Optional buffer distance (meters) around reference polygons when rasterizing for pixel validation. Default is 0.
#' @param threshold_completely_detected Numeric. Minimum percent (e.g., 90) of a reference polygon area that must be intersected to be considered completely detected. Default is 90.
#' @param min_area_reference_ha Numeric. Minimum area (in hectares) to retain reference polygons after masking. Default is NULL (no filtering).
#' @param use_gdal Logical. Whether to use external GDAL (via Python) for polygonizing rasters (faster for large datasets). Default is TRUE.
#' @param python_exe Character. Path to the Python executable.
#' @param gdal_polygonize_script Character. Path to `gdal_polygonize.py` script.
#' @param force_reprocess_ref Logical. If TRUE, forces the recalculation of masked reference polygons even if they already exist. Default is FALSE.
#'
#' @details
#' ## Pixel-based metrics:
#' - **True Positives (TP)**: Pixels correctly detected as burned.
#' - **False Positives (FP)**: Pixels wrongly detected as burned.
#' - **False Negatives (FN)**: Pixels missed (burned in reference but not detected).
#' - **True Negatives (TN)**: Pixels correctly identified as unburned.
#'
#' Derived indicators:
#' - **Precision** = TP / (TP + FP)
#' - **Recall** = TP / (TP + FN)
#' - **F1 Score** = 2 ? (Precision ? Recall) / (Precision + Recall)
#' - **Intersection over Union (IoU)** = TP / (TP + FP + FN)
#' - **Error Rate** = (FP + FN) / (TP + FP + FN + TN)
#'
#' ## Area-based metrics:
#' - **N_Reference_Polygons**: Number of reference polygons considered after masking.
#' - **N_Completely_Detected**: Number of reference polygons detected over `threshold_completely_detected` % of their area.
#' - **N_Detected_Polygons**: Number of polygons partially or fully detected (>0% intersection).
#' - **N_Not_Detected**: Number of reference polygons without any intersection.
#' - **Perc_Detected_Polygons**: Percentage of reference polygons detected.
#' - **Area_Reference_ha**: Total area (ha) of reference polygons.
#' - **Area_Detected_ha**: Total area (ha) of burned polygons provided as input.
#' - **Area_Intersection_ha**: Total intersected area (ha) between detected and reference polygons.
#' - **Area_Reference_NotDetected_ha**: Area (ha) of reference polygons not intersected.
#' - **Perc_Reference_Area_NotDetected**: Percentage of total reference area not intersected.
#' - **Recall_Area_percent** = (Area_Intersection_ha / Area_Reference_ha) ? 100
#' - **Precision_Area_percent** = (Area_Intersection_ha / Area_Detected_ha) ? 100
#'
#' ## Additional outputs:
#' - Shapefiles of undetected reference polygons and unmatched detection polygons are saved in the `VALIDATION` subdirectory.
#'
#' ## Notes:
#' - Recall and Precision at area scale are **not penalized** if detected polygons cover more area than the reference.
#' - Always set `force_reprocess_ref = TRUE` when changing `min_area_reference_ha`.
#'
#' @return
#' A list containing:
#' - `metrics`: data.table with pixel-based metrics.
#' - `polygon_summary`: data.table with area-based metrics.
#'
#' Two CSV files are saved automatically inside the `VALIDATION` subdirectory of `validation_dir`.
#'
#' @examples
#' \dontrun{
#' burned_dir <- "C:/MY_PROJECT/FIRE_MAPPING/2012/OTSU"
#' polygons_sf <- list.files(
#'     burned_dir,
#'     pattern = glob2rx("burned_areas_2012_otsu_corine*.shp"),
#'     full.names = TRUE
#' )
#' validate_fire_maps(
#'     input_shapefile = polygons_sf,
#'     ref_shapefile = "C:/MY_PROJECT/REFERENCE/fires_2012.shp",
#'     mask_shapefile = "C:/MY_PROJECT/STUDY_AREA/portugal_mask.shp",
#'     burnable_raster = "C:/MY_PROJECT/BURNABLE_MASK/burnable_binary.tif",
#'     year_target = 2012,
#'     validation_dir = "C:/MY_PROJECT/FIRE_MAPPING/2012",
#'     binary_burnable = TRUE,
#'     burnable_classes = NULL,
#'     min_area_reference_ha = 2,
#'     buffer = 0,
#'     threshold_completely_detected = 90,
#'     use_gdal = TRUE,
#'     python_exe = "C:/ProgramData/anaconda3/python.exe",
#'     gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py",
#'     force_reprocess_ref = TRUE
#' )
#' }
#' @importFrom sf st_read st_write st_crs st_transform st_make_valid st_intersection st_area st_buffer st_is_empty
#' @importFrom terra rast crop mask project rasterize as.polygons vect writeRaster global
#' @importFrom data.table data.table fwrite rbindlist
#' @importFrom glue glue
#' @importFrom tools file_path_sans_ext
#' @importFrom dplyr filter
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom dplyr filter first
#' @importFrom data.table :=
#' @importFrom magrittr %>%
#' @export

validate_fire_maps <- function(input_shapefile,
                               ref_shapefile,
                               mask_shapefile,
                               burnable_raster,
                               year_target,
                               validation_dir,
                               binary_burnable = TRUE,
                               burnable_classes = NULL,
                               buffer = 0,
                               threshold_completely_detected = 90,
                               min_area_reference_ha = NULL,
                               use_gdal = TRUE,
                               python_exe = "C:/ProgramData/anaconda3/python.exe",
                               gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py",
                               force_reprocess_ref = FALSE) {

  validation_output_dir <- file.path(validation_dir, "VALIDATION")
  if (!dir.exists(validation_output_dir)) {
    dir.create(validation_output_dir, recursive = TRUE)
  }

  # Load mask and burnable raster
  mask_geom <- sf::st_read(mask_shapefile, quiet = TRUE) |> sf::st_make_valid()
  burnable <- terra::rast(burnable_raster)

  ref_masked_path <- file.path(validation_output_dir, paste0("ref_polygons_processed_", year_target, ".shp"))

  # Forced deletion if requested
  if (force_reprocess_ref) {
    message("Force reprocessing reference polygons: deleting cached shapefile...")
    unlink(file.path(validation_dir, "VALIDATION", paste0("ref_polygons_processed_", year_target, ".shp")))
    unlink(list.files(file.path(validation_dir, "VALIDATION"),
                      pattern = paste0("ref_polygons_processed_", year_target, "\\..*"),
                      full.names = TRUE))
  }

  if (file.exists(ref_masked_path)) {
    message("Loading cached masked reference polygons...")
    ref_polygons <- sf::st_read(ref_masked_path, quiet = TRUE)
  } else {
    message("Processing reference polygons...")
    ref_polygons <- sf::st_read(ref_shapefile, quiet = TRUE) |> sf::st_make_valid()

    if ("year" %in% names(ref_polygons)) {
      ref_polygons <- dplyr::filter(ref_polygons, year == year_target)
    }
    if (nrow(ref_polygons) == 0) stop("No reference polygons found for the specified year.")

    if (sf::st_crs(ref_polygons) != sf::st_crs(mask_geom)) {
      ref_polygons <- sf::st_transform(ref_polygons, sf::st_crs(mask_geom))
    }
    if (terra::crs(burnable) != sf::st_crs(mask_geom)$wkt) {
      burnable <- terra::project(burnable, sf::st_crs(mask_geom)$wkt)
    }

    ref_polygons <- suppressWarnings(sf::st_intersection(ref_polygons, mask_geom))
    ref_raster <- terra::rasterize(terra::vect(ref_polygons), burnable, field = 1, background = NA)
    masked_ref_raster <- ref_raster * burnable
    masked_ref_raster[masked_ref_raster != 1] <- NA

    raster_temp_path <- file.path(tempdir(), "masked_ref_raster.tif")
    terra::writeRaster(masked_ref_raster, raster_temp_path, overwrite = TRUE,
                       datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))

    if (use_gdal) {
      output_shapefile_path <- file.path(tempdir(), "masked_reference_polygons.shp")

      # Lanzar GDAL
      system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{raster_temp_path}" -f "ESRI Shapefile" "{output_shapefile_path}" DN'))

      # Verificar que se haya creado el shapefile
      if (!file.exists(output_shapefile_path)) {
        stop("GDAL polygonization failed: masked_reference_polygons.shp was not created.")
      }

      # Leer shapefile generado
      masked_ref_polygons <- sf::st_read(output_shapefile_path, quiet = TRUE)

      # Borrar shapefiles temporales
      unlink(list.files(tempdir(), pattern = "masked_reference_polygons\\..*", full.names = TRUE))

      # Verificar que tenga geometrias
      if (nrow(masked_ref_polygons) == 0) {
        stop("Polygonized reference is empty after GDAL. Check burnable mask and inputs.")
      }

      # Eliminar geometrias vacias si las hubiera
      masked_ref_polygons <- masked_ref_polygons[!sf::st_is_empty(masked_ref_polygons), ]


    } else {
      masked_ref_polygons <- terra::as.polygons(masked_ref_raster, dissolve = FALSE) |> sf::st_as_sf()
      masked_ref_polygons <- masked_ref_polygons[!sf::st_is_empty(masked_ref_polygons), ]
    }


    # Filtrado por area minima de referencia (nuevo)
    if (!is.null(min_area_reference_ha)) {
      area_ref_polygons <- as.numeric(sf::st_area(ref_polygons)) / 10000  # Area en hectareas
      n_before <- nrow(ref_polygons)
      ref_polygons <- ref_polygons[area_ref_polygons >= min_area_reference_ha, ]
      n_after <- nrow(ref_polygons)

      cat(sprintf("Filtered small reference polygons: %d ? %d polygons (Area ? %.2f ha)\n", n_before, n_after, min_area_reference_ha))

      if (nrow(ref_polygons) == 0) {
        stop("No reference polygons remain after filtering by minimum area.")
      }
    }

    sf::st_write(ref_polygons, ref_masked_path, delete_layer = TRUE, quiet = TRUE)

  }

  if (!is.list(input_shapefile)) input_shapefile <- as.list(input_shapefile)

  metrics_list <- list()
  polygon_summary_list <- list()

  for (shp in input_shapefile) {
    input_name <- tools::file_path_sans_ext(basename(shp))
    raster_path <- sub("\\.shp$", ".tif", shp)

    if (!file.exists(raster_path)) {
      detection_polygons <- sf::st_read(shp, quiet = TRUE) |> sf::st_make_valid()
      if (sf::st_crs(detection_polygons) != sf::st_crs(ref_polygons)) {
        detection_polygons <- sf::st_transform(detection_polygons, sf::st_crs(ref_polygons))
      }

      detection_polygons <- suppressWarnings(sf::st_intersection(detection_polygons, mask_geom))

      base_raster <- terra::crop(burnable, terra::vect(detection_polygons))
      burned_raster <- terra::rasterize(terra::vect(detection_polygons), base_raster, field = 1, background = NA)
      terra::writeRaster(burned_raster, raster_path, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
    }

    pred_raster <- terra::rast(raster_path)

    detection_polygons <- sf::st_read(shp, quiet = TRUE)
    if (sf::st_crs(detection_polygons) != sf::st_crs(ref_polygons)) detection_polygons <- sf::st_transform(detection_polygons, sf::st_crs(ref_polygons))
    detection_polygons <- suppressWarnings(sf::st_intersection(detection_polygons, mask_geom))

    ## Pixel metrics
    compute_metrics <- function(pred_raster, ref_polygons, buffer = 0) {
      if (buffer > 0) ref_polygons <- sf::st_buffer(ref_polygons, dist = buffer)
      ref_raster <- terra::rasterize(terra::vect(ref_polygons), pred_raster, field = 1, background = 0)
      mask_na <- is.na(pred_raster) | is.na(ref_raster)
      pred_raster[mask_na] <- NA
      ref_raster[mask_na] <- NA

      tp <- terra::global((pred_raster == 1) & (ref_raster == 1), "sum", na.rm = TRUE)[[1]]
      fp <- terra::global((pred_raster == 1) & (ref_raster == 0), "sum", na.rm = TRUE)[[1]]
      fn <- terra::global((pred_raster == 0) & (ref_raster == 1), "sum", na.rm = TRUE)[[1]]
      tn <- terra::global((pred_raster == 0) & (ref_raster == 0), "sum", na.rm = TRUE)[[1]]

      precision <- tp / (tp + fp)
      recall <- tp / (tp + fn)
      f1 <- 2 * precision * recall / (precision + recall)
      iou <- tp / (tp + fp + fn)
      error_rate <- (fp + fn) / (tp + fp + fn + tn)

      list(TP = tp, FP = fp, FN = fn, TN = tn, Precision = precision, Recall = recall, F1 = f1, IoU = iou, ErrorRate = error_rate)
    }

    message("Computing pixel-based metrics...")

    pixel_metrics_dt <- data.table::as.data.table(compute_metrics(pred_raster, ref_polygons, buffer))
    pixel_metrics_dt[, `:=`(InputName = input_name, Year = year_target)]


    # Area-based metrics
    ref_polygons <- sf::st_make_valid(ref_polygons)
    detection_polygons <- sf::st_make_valid(detection_polygons)

    intersec <- suppressWarnings(sf::st_intersection(ref_polygons, detection_polygons))

    area_ref <- as.numeric(sf::st_area(ref_polygons)) / 10000  # ha
    area_intersect <- numeric(length(area_ref))

    # calcular ids para deteccion
    ids <- sf::st_intersects(ref_polygons, detection_polygons)

    n_detected <- sum(lengths(ids) > 0)  # Aqui definimos explicitamente cuantos poligonos de referencia estan detectados (intersectados)

    if (nrow(intersec) > 0) {
      ids_intersec <- sf::st_intersects(ref_polygons, intersec)
      for (j in seq_along(ids_intersec)) {
        if (length(ids_intersec[[j]]) > 0) {
          parts <- intersec[ids_intersec[[j]], ]
          area_intersect[j] <- sum(as.numeric(sf::st_area(parts))) / 10000
        }
      }
    }

    percent_detected <- (area_intersect / area_ref) * 100
    percent_detected[is.nan(percent_detected)] <- 0

    n_total <- length(area_ref)
    n_not_detected <- n_total - n_detected
    perc_detected_polygons <- round((n_detected / n_total) * 100, 2)

    completely_detected <- sum(percent_detected >= threshold_completely_detected)

    area_reference_total <- round(sum(area_ref, na.rm = TRUE), 2)
    area_detected_total <- round(sum(as.numeric(sf::st_area(detection_polygons)), na.rm = TRUE) / 10000, 2)
    area_detected_in_reference <- round(sum(area_intersect, na.rm = TRUE), 2)

    recall_area <- ifelse(area_reference_total > 0, (area_detected_in_reference / area_reference_total) * 100, NA)
    precision_area <- ifelse(area_detected_total > 0, (area_detected_in_reference / area_detected_total) * 100, NA)

    # Ahora si el resumen:
    polygon_summary_dt <- data.table::data.table(
      InputName = input_name,
      Year = year_target,
      N_Reference_Polygons = n_total,
      N_Completely_Detected = completely_detected,
      N_Detected_Polygons = n_detected,
      N_Not_Detected = n_not_detected,
      Perc_Detected_Polygons = perc_detected_polygons,
      Area_Reference_ha = area_reference_total,
      Area_Detected_ha = area_detected_total,
      Area_Intersection_ha = area_detected_in_reference,
      Recall_Area_percent = round(recall_area, 2),
      Precision_Area_percent = round(precision_area, 2)
    )

    # Guardar shapefiles de no detectados y calcular areas

    # Identificar poligonos de referencia no detectados
    ref_not_detected_idx <- lengths(ids) == 0
    ref_not_detected_polygons <- ref_polygons[ref_not_detected_idx, ]

    # Identificar poligonos de deteccion no intersectados con referencia
    det_ids <- sf::st_intersects(detection_polygons, ref_polygons)
    det_not_matched_idx <- lengths(det_ids) == 0
    det_not_matched_polygons <- detection_polygons[det_not_matched_idx, ]

    # Guardar shapefiles de salida
    ref_not_detected_path <- file.path(validation_output_dir, paste0("ref_polygons_not_detected_", input_name, ".shp"))
    det_not_matched_path <- file.path(validation_output_dir, paste0("input_polygons_not_matched_", input_name, ".shp"))

    if (nrow(ref_not_detected_polygons) > 0) {
      sf::st_write(ref_not_detected_polygons, ref_not_detected_path, delete_layer = TRUE, quiet = TRUE)
    }

    if (nrow(det_not_matched_polygons) > 0) {
      sf::st_write(det_not_matched_polygons, det_not_matched_path, delete_layer = TRUE, quiet = TRUE)
    }

    # Calcular area no detectada
    area_not_detected_total <- round(sum(as.numeric(sf::st_area(ref_not_detected_polygons)), na.rm = TRUE) / 10000, 2)
    perc_area_not_detected <- ifelse(area_reference_total > 0, (area_not_detected_total / area_reference_total) * 100, NA)

    # Anadir al resumen de poligonos
    polygon_summary_dt[, `:=`(
      Area_Reference_NotDetected_ha = area_not_detected_total,
      Perc_Reference_Area_NotDetected = round(perc_area_not_detected, 2)
    )]


    metrics_list[[input_name]] <- pixel_metrics_dt
    polygon_summary_list[[input_name]] <- polygon_summary_dt
  }

  all_metrics <- data.table::rbindlist(metrics_list)
  all_polygon_summary <- data.table::rbindlist(polygon_summary_list)

  data.table::fwrite(all_metrics, file.path(validation_output_dir, paste0("metrics_summary_", year_target, ".csv")))
  data.table::fwrite(all_polygon_summary, file.path(validation_output_dir, paste0("polygon_summary_", year_target, ".csv")))

  return(list(metrics = all_metrics, polygon_summary = all_polygon_summary))
}

