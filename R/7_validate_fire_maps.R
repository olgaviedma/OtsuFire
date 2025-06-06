#' Validate burned area maps against reference polygons
#'
#' @description
#' The `validate_fire_maps()` function evaluates the spatial accuracy of burned area detection shapefiles by comparing them against independent reference fire polygons (e.g., from Focclim or national databases).
#'
#' It calculates both **pixel-based metrics** and **area-based metrics**, unless a specific mode is selected. Reference polygons are masked to retain only burnable areas based on a CORINE-derived raster and filtered by a minimum area threshold if specified.
#'
#' **Important**: If `min_area_reference_ha` is changed, you must set `force_reprocess_ref = TRUE` to recalculate the masked reference polygons.
#'
#' @name validate_fire_maps
#' @rdname validate_fire_maps
#'
#' @param input_shapefile Character vector. One or more paths to shapefiles containing the burned polygons to validate.
#' @param ref_shapefile Character. Path to the shapefile with reference burned area polygons.
#' @param mask_shapefile Character. Path to the shapefile defining the study area boundary.
#' @param burnable_raster Character. Path to the raster file indicating burnable areas (binary or categorical).
#' @param year_target Numeric. Target year for filtering reference polygons.
#' @param validation_dir Character. Output directory to save validation results.
#' @param binary_burnable Logical. TRUE if the burnable raster is binary (1 = burnable, 0 = non-burnable). Default is TRUE.
#' @param burnable_classes Optional numeric vector of raster values considered burnable if `binary_burnable = FALSE`.
#' @param class_shape Optional character. Path to a shapefile containing class information (e.g., CORINE or ecoregions) used for class-wise error breakdown.
#' @param class_field Optional character. Name of the field in `class_shape` to group omission and commission errors by class (e.g., "CORINE", "eco_name", "cor_eco").
#' @param buffer Numeric. Optional buffer distance (in meters) applied around reference polygons for pixel-based validation. Default is 0.
#' @param threshold_completely_detected Numeric. Minimum percentage (e.g., 90) of a reference polygon area that must be intersected to be considered completely detected. Default is 90.
#' @param min_area_reference_ha Numeric. Minimum area (in hectares) to retain reference polygons after masking. Default is NULL (no filtering).
#' @param use_gdal Logical. Whether to use external GDAL (via Python) for polygonizing rasters (faster for large datasets). Default is TRUE.
#' @param python_exe Character. Path to the Python executable used for GDAL polygonization.
#' @param gdal_polygonize_script Character. Path to the `gdal_polygonize.py` script.
#' @param force_reprocess_ref Logical. If TRUE, forces recalculation of masked reference polygons even if cached versions exist. Default is FALSE.
#' @param metrics_type Character. Type of metrics to compute. One of `"all"` (default), `"pixel"`, or `"area"`.
#'
#' @details
#' ## Pixel-based metrics (when `metrics_type = "pixel"` or `"all"`):
#' - **True Positives (TP)**: Pixels correctly detected as burned.
#' - **False Positives (FP)**: Pixels wrongly detected as burned.
#' - **False Negatives (FN)**: Burned pixels missed by the detection.
#' - **True Negatives (TN)**: Pixels correctly identified as unburned.
#'
#' Derived indicators:
#' - **Precision** = TP / (TP + FP)
#' - **Recall** = TP / (TP + FN)
#' - **F1 Score** = 2 * (Precision * Recall) / (Precision + Recall)
#' - **Intersection over Union (IoU)** = TP / (TP + FP + FN)
#' - **Error Rate** = (FP + FN) / (TP + FP + FN + TN)
#'
#' ## Area-based metrics (when `metrics_type = "area"` or `"all"`):
#' - **N_Reference_Polygons**: Number of reference polygons after masking and filtering.
#' - **N_Completely_Detected**: Count of reference polygons where detected area ? `threshold_completely_detected`.
#' - **N_Detected_Polygons**: Reference polygons partially or fully detected (>0% overlap).
#' - **N_Not_Detected**: Reference polygons without any detected overlap.
#' - **Perc_Detected_Polygons**: Share of reference polygons detected.
#' - **Area_Reference_ha**: Total area of reference polygons (ha).
#' - **Area_Detected_ha**: Total area of detected burned polygons (ha).
#' - **Area_Intersection_ha**: Area of intersection between detection and reference polygons (ha).
#' - **Area_Reference_NotDetected_ha**: Area of reference polygons not intersected (ha).
#' - **Perc_Reference_Area_NotDetected**: Share of reference area missed (%).
#' - **Recall_Area_percent** = (Area_Intersection_ha / Area_Reference_ha) * 100
#' - **Precision_Area_percent** = (Area_Intersection_ha / Area_Detected_ha) * 100
#'
#' ## Class-wise error breakdown:
#' If `class_shape` and `class_field` are provided and valid:
#' - `ref_not_detected` and `det_not_matched` polygons are joined to `class_shape`.
#' - Output CSV files:
#'   - `omission_by_<class_field>_<input>.csv`
#'   - `commission_by_<class_field>_<input>.csv`
#'
#' ## Output files:
#' - `metrics_summary_<year>.csv` and/or `polygon_summary_<year>.csv`
#' - Shapefiles of undetected and unmatched polygons
#' - Class-wise omission/commission CSVs (if class information is provided)
#'
#' @return A list with:
#' \describe{
#'   \item{metrics}{A data.table of pixel-based accuracy metrics, or NULL if `metrics_type` excludes them.}
#'   \item{polygon_summary}{A data.table of area-based metrics, or NULL if `metrics_type` excludes them.}
#' }
#'
#' @examples
#' \dontrun{
#' validate_fire_maps(
#'   input_shapefile = list.files("shapefiles", pattern = "\\.shp$", full.names = TRUE),
#'   ref_shapefile = "ref_polygons_2022.shp",
#'   mask_shapefile = "mask_region.shp",
#'   burnable_raster = "burnable.tif",
#'   year_target = 2022,
#'   validation_dir = "validation_results",
#'   binary_burnable = TRUE,
#'   min_area_reference_ha = 1,
#'   buffer = 30,
#'   threshold_completely_detected = 90,
#'   use_gdal = TRUE,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py",
#'   force_reprocess_ref = TRUE,
#'   metrics_type = "all",
#'   class_shape = "corine_eco_regions.shp",
#'   class_field = "cor_eco"
#' )
#' }
#'
#' @importFrom sf st_read st_write st_crs st_transform st_make_valid st_intersection st_area st_buffer st_is_empty
#' @importFrom terra rast crop mask project rasterize as.polygons vect writeRaster global ncell
#' @importFrom data.table data.table fwrite rbindlist
#' @importFrom glue glue
#' @importFrom tools file_path_sans_ext
#' @importFrom dplyr filter group_by summarise mutate ungroup n
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @export

utils::globalVariables(c(
  ".data", "CORINE_CLASS", "CORINE_YEAR", "ECO_CLASS",
  "unit_id2", "geometry"
))

validate_fire_maps <- function(input_shapefile,
                                ref_shapefile,
                                mask_shapefile,
                                burnable_raster,
                                year_target,
                                validation_dir,
                                binary_burnable = TRUE,
                                burnable_classes = NULL,
                                class_shape,
                                class_field = NULL,
                                buffer = 0,
                                threshold_completely_detected = 90,
                                min_area_reference_ha = NULL,
                                use_gdal = TRUE,
                                python_exe = "C:/ProgramData/anaconda3/python.exe",
                                gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py",
                                force_reprocess_ref = FALSE,
                                metrics_type = c("all", "pixel", "area")) {

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
    cat("Original reference polygons", nrow(ref_polygons), "\n")

    if ("year" %in% names(ref_polygons)) {
      ref_polygons <- dplyr::filter(ref_polygons, year == year_target)
      cat("reference polygons filtered by (", year_target, "):", nrow(ref_polygons), "\n")
    }
    if (nrow(ref_polygons) == 0) stop("No reference polygons found for the specified year.")

    # Reproyectar ref_polygons y burnable si es necesario
    if (sf::st_crs(ref_polygons) != sf::st_crs(mask_geom)) {
      ref_polygons <- sf::st_transform(ref_polygons, sf::st_crs(mask_geom))
    }
    if (terra::crs(burnable) != sf::st_crs(mask_geom)$wkt) {
      burnable <- terra::project(burnable, sf::st_crs(mask_geom)$wkt)
    }

    # Validar y armonizar geometrias
    ref_polygons <- sf::st_make_valid(ref_polygons)
    mask_geom <- sf::st_make_valid(mask_geom)

    # Prefiltrar por BBOX e intersectar con mascara
    ref_polygons <- sf::st_filter(ref_polygons, mask_geom, .predicate = sf::st_intersects)
    cat("reference polygons filtered by mask (st_filter):", nrow(ref_polygons), "\n")

    if (nrow(ref_polygons) > 0) {
      ref_polygons <- suppressWarnings(sf::st_intersection(ref_polygons, mask_geom))
    }

    # Rasterizar referencia y multiplicar por mascara
    ref_raster <- terra::rasterize(terra::vect(ref_polygons), burnable, field = 1, background = NA)
    masked_ref_raster <- ref_raster * burnable
    masked_ref_raster[masked_ref_raster < 0.99] <- NA

    raster_temp_path <- file.path(tempdir(), "masked_ref_raster.tif")
    terra::writeRaster(masked_ref_raster, raster_temp_path, overwrite = TRUE,
                       datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))

    if (use_gdal) {
      output_shapefile_path <- file.path(tempdir(), "masked_reference_polygons.shp")
      system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{raster_temp_path}" -f "ESRI Shapefile" "{output_shapefile_path}" DN'))

      if (!file.exists(output_shapefile_path)) {
        stop("GDAL polygonization failed: masked_reference_polygons.shp was not created.")
      }

      masked_ref_polygons <- sf::st_read(output_shapefile_path, quiet = TRUE)
      unlink(list.files(tempdir(), pattern = "masked_reference_polygons\\..*", full.names = TRUE))

      if (nrow(masked_ref_polygons) == 0) {
        stop("Polygonized reference is empty after GDAL. Check burnable mask and inputs.")
      }

      masked_ref_polygons <- masked_ref_polygons[!sf::st_is_empty(masked_ref_polygons), ]

    } else {
      masked_ref_polygons <- terra::as.polygons(masked_ref_raster, dissolve = FALSE) |> sf::st_as_sf()
      masked_ref_polygons <- masked_ref_polygons[!sf::st_is_empty(masked_ref_polygons), ]
    }

    # Filtrar por area minima si se especifica
    if (!is.null(min_area_reference_ha)) {
      area_ref_polygons <- as.numeric(sf::st_area(masked_ref_polygons)) / 10000
      n_before <- nrow(masked_ref_polygons)
      masked_ref_polygons <- masked_ref_polygons[area_ref_polygons >= min_area_reference_ha, ]
      n_after <- nrow(masked_ref_polygons)

      cat(sprintf("Filtered small reference polygons: %d ? %d polygons (Area ? %.2f ha)\n",
                  n_before, n_after, min_area_reference_ha))

      if (nrow(masked_ref_polygons) == 0) {
        stop("No reference polygons remain after filtering by minimum area.")
      }
    }

    # Guardar y asignar a ref_polygons
    sf::st_write(masked_ref_polygons, ref_masked_path, delete_layer = TRUE, quiet = TRUE)
    ref_polygons <- masked_ref_polygons
  }


  if (!is.list(input_shapefile)) input_shapefile <- as.list(input_shapefile)

  metrics_list <- list()
  polygon_summary_list <- list()

  for (shp in input_shapefile) {
      input_name <- tools::file_path_sans_ext(basename(shp))
      raster_path <- sub("\\.shp$", ".tif", shp)
      message("Processing input file: ", input_name)

    if (!file.exists(raster_path)) {
      detection_polygons <- sf::st_read(shp, quiet = TRUE) |> sf::st_make_valid()
      if (sf::st_crs(detection_polygons) != sf::st_crs(ref_polygons)) {
        detection_polygons <- sf::st_transform(detection_polygons, sf::st_crs(ref_polygons))
      }

      if (sf::st_crs(detection_polygons) != terra::crs(burnable)) {
        detection_polygons <- sf::st_transform(detection_polygons, terra::crs(burnable))
      }

      # Filtrar primero solo los que intersectan con la mascara
      detection_polygons <- sf::st_filter(detection_polygons, mask_geom, .predicate = sf::st_intersects)

      # Solo intersectar si hay elementos dentro del area de mascara
      if (nrow(detection_polygons) > 0) {
        detection_polygons <- suppressWarnings(sf::st_intersection(detection_polygons, mask_geom))
      } else {
        detection_polygons <- detection_polygons[0, ]
      }


      # Verificar si detection_polygons tiene geometria valida
      if (nrow(detection_polygons) == 0 || is.null(sf::st_geometry(detection_polygons))) {
        warning(paste("Skipping", input_name, "- no valid polygons for cropping."))
        next
      }

      # Verificar si los extents se solapan
      if (is.null(terra::intersect(terra::ext(burnable), terra::ext(terra::vect(detection_polygons))))){
        warning(paste("Skipping", input_name, "- burnable raster and detection polygons do not overlap."))
        next
      }

      base_raster <- terra::crop(burnable, terra::vect(detection_polygons))
      burned_raster <- terra::rasterize(terra::vect(detection_polygons), base_raster, field = 1, background = NA)
      terra::writeRaster(burned_raster, raster_path, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
    }

    pred_raster <- terra::rast(raster_path)

    detection_polygons <- sf::st_read(shp, quiet = TRUE)
    message("Processing detection polygons - ", nrow(detection_polygons), " polygons")

    if (sf::st_crs(detection_polygons) != sf::st_crs(ref_polygons)) detection_polygons <- sf::st_transform(detection_polygons, sf::st_crs(ref_polygons))
    detection_polygons <- suppressWarnings(sf::st_intersection(detection_polygons, mask_geom))

    if (metrics_type %in% c("all", "pixel")) {
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

      pixel_metrics_dt <- data.table::as.data.table(compute_metrics(pred_raster, ref_polygons, buffer))
      pixel_metrics_dt[, `:=`(InputName = input_name, Year = year_target)]
      metrics_list[[input_name]] <- pixel_metrics_dt
    }

    message("Computing metrics:")
    if (metrics_type %in% c("all", "area")) {
      ref_polygons <- sf::st_make_valid(ref_polygons)
      detection_polygons <- sf::st_make_valid(detection_polygons)

      # Prefiltrar detecciones que potencialmente intersectan para acelerar la interseccion
      detection_polygons_filtered <- sf::st_filter(detection_polygons, ref_polygons, .predicate = sf::st_intersects)

      if (nrow(detection_polygons_filtered) > 0) {
        intersec <- suppressWarnings(sf::st_intersection(ref_polygons, detection_polygons_filtered))
      } else {
        intersec <- detection_polygons[0, ]  # objeto vacio
      }


      area_ref <- as.numeric(sf::st_area(ref_polygons)) / 10000
      area_intersect <- numeric(length(area_ref))
      ids <- sf::st_intersects(ref_polygons, detection_polygons)
      n_detected <- sum(lengths(ids) > 0)

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

      ref_not_detected_idx <- lengths(ids) == 0
      ref_not_detected_polygons <- ref_polygons[ref_not_detected_idx, ]
      det_ids <- sf::st_intersects(detection_polygons, ref_polygons)
      det_not_matched_idx <- lengths(det_ids) == 0
      det_not_matched_polygons <- detection_polygons[det_not_matched_idx, ]

      ref_not_detected_path <- file.path(validation_output_dir, paste0("ref_polygons_not_detected_", input_name, ".shp"))
      det_not_matched_path <- file.path(validation_output_dir, paste0("input_polygons_not_matched_", input_name, ".shp"))

      message("Saving not detected polygons:")
      if (nrow(ref_not_detected_polygons) > 0) {
        sf::st_write(ref_not_detected_polygons, ref_not_detected_path, delete_layer = TRUE, quiet = TRUE)
      }
      if (nrow(det_not_matched_polygons) > 0) {
        sf::st_write(det_not_matched_polygons, det_not_matched_path, delete_layer = TRUE, quiet = TRUE)
      }

      area_not_detected_total <- round(sum(as.numeric(sf::st_area(ref_not_detected_polygons)), na.rm = TRUE) / 10000, 2)
      perc_area_not_detected <- ifelse(area_reference_total > 0, (area_not_detected_total / area_reference_total) * 100, NA)

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
        Area_Reference_NotDetected_ha = area_not_detected_total,
        Perc_Reference_Area_NotDetected = round(perc_area_not_detected, 2),
        Recall_Area_percent = round(recall_area, 2),
        Precision_Area_percent = round(precision_area, 2)
      )

      polygon_summary_list[[input_name]] <- polygon_summary_dt

      message("Analysing errors by classes:")


      # === ANALIZAR ERRORES POR CLASE ===
      # Load class map only if needed
      class_map <- NULL
      if (!is.null(class_shape) && file.exists(class_shape)) {
        class_map <- sf::st_read(class_shape, quiet = TRUE)
        class_map <- sf::st_make_valid(class_map)
        class_map <- sf::st_transform(class_map, sf::st_crs(ref_polygons))
      }

      # Only proceed with class-based omission/commission if class_field is valid
      if (!is.null(class_field) && !is.null(class_map) && class_field %in% names(class_map)) {
        # Spatial join with class map
        ref_not_detected_polygons <- sf::st_join(ref_not_detected_polygons, class_map[, class_field], left = TRUE)
        det_not_matched_polygons <- sf::st_join(det_not_matched_polygons, class_map[, class_field], left = TRUE)

        clean_class_field <- function(sf_obj, class_field) {
          if (paste0(class_field, ".y") %in% names(sf_obj)) {
            sf_obj[[class_field]] <- sf_obj[[paste0(class_field, ".y")]]
            sf_obj[[paste0(class_field, ".y")]] <- NULL
            if (paste0(class_field, ".x") %in% names(sf_obj)) {
              sf_obj[[paste0(class_field, ".x")]] <- NULL
            }
          }
          return(sf_obj)
        }

        ref_not_detected_polygons <- clean_class_field(ref_not_detected_polygons, class_field)
        det_not_matched_polygons <- clean_class_field(det_not_matched_polygons, class_field)

        if (class_field %in% names(ref_not_detected_polygons)) {
          omission_by_class <- ref_not_detected_polygons |>
            dplyr::mutate(area_ha = as.numeric(sf::st_area(geometry)) / 10000) |>
            dplyr::group_by(!!rlang::sym(class_field)) |>
            dplyr::summarise(
              N_Omitted = dplyr::n(),
              Area_Omitted_ha = sum(area_ha, na.rm = TRUE),
              .groups = "drop"
            ) |>
            sf::st_drop_geometry()

          omission_by_class[[class_field]] <- as.character(omission_by_class[[class_field]])

          data.table::fwrite(
            omission_by_class,
            file.path(validation_output_dir, paste0("omission_by_", class_field, "_", input_name, ".csv"))
          )
        }

        if (class_field %in% names(det_not_matched_polygons)) {
          commission_by_class <- det_not_matched_polygons |>
            dplyr::mutate(area_ha = as.numeric(sf::st_area(geometry)) / 10000) |>
            dplyr::group_by(!!rlang::sym(class_field)) |>
            dplyr::summarise(
              N_Commission = dplyr::n(),
              Area_Commission_ha = sum(area_ha, na.rm = TRUE),
              .groups = "drop"
            ) |>
            sf::st_drop_geometry()

          commission_by_class[[class_field]] <- as.character(commission_by_class[[class_field]])

          data.table::fwrite(
            commission_by_class,
            file.path(validation_output_dir, paste0("commission_by_", class_field, "_", input_name, ".csv"))
          )
        }
      } else if (!is.null(class_field)) {
        warning(paste0("Field '", class_field, "' not found or class map not provided. Skipping class-wise validation."))
      }

    }
    }

  if (metrics_type %in% c("all", "pixel")) {
    all_metrics <- data.table::rbindlist(metrics_list)
    data.table::fwrite(all_metrics, file.path(validation_output_dir, paste0("metrics_summary_", year_target, ".csv")))
  }

  if (metrics_type %in% c("all", "area")) {
    all_polygon_summary <- data.table::rbindlist(polygon_summary_list)
    data.table::fwrite(all_polygon_summary, file.path(validation_output_dir, paste0("polygon_summary_", year_target, ".csv")))
  }

  return(list(
    metrics = if (metrics_type %in% c("all", "pixel")) all_metrics else NULL,
    polygon_summary = if (metrics_type %in% c("all", "area")) all_polygon_summary else NULL
  ))
}
