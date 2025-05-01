#' Calculate DOY-Based Flags and Summary Statistics from Raster and Fire Polygons
#'
#' This function extracts Day-of-Year (DOY) values from a raster (e.g., from a Landsat composite),
#' masked by burned area polygons, and computes per-polygon statistics: mode, median,
#' and optionally percentiles (e.g., 5th, 25th, 95th). These DOY values are converted to calendar dates,
#' and their components (day, month, year) are also included.
#'
#' Optionally, the function can threshold the DOY band around the median and/or mode of each polygon
#' using user-defined DOY windows (e.g., ±10 days), create binary raster masks, and vectorize them using GDAL.
#'
#' The function supports either a single shapefile (original burned areas) or a named list
#' of shapefiles or `sf` objects generated with different Otsu thresholds.
#'
#' @param raster A multi-band `SpatRaster` object. One of the bands must contain DOY values.
#' @param doy_band Integer. The index of the DOY band in the raster (default is 2).
#' @param polygons_sf A single shapefile path, an `sf` object, or a named list of shapefiles or `sf` objects.
#'                    Names are used to label outputs (e.g., `"ge50"`, `"ge100"`).
#' @param output_dir Directory where output raster and shapefile files will be saved.
#' @param year Calendar year (numeric or character) used for DOY-to-date conversion and filenames.
#' @param doy_thresholds Numeric vector of DOY windows (e.g., `c(10, 15)`) to apply around the selected statistic(s).
#' @param stats Statistic(s) to use for thresholding. One of `"median"`, `"mode"`, or `"both"` (default `"both"`).
#' @param calc_percentiles Logical. If `TRUE`, compute additional DOY percentiles. Default is `TRUE`.
#' @param percentiles Numeric vector of percentiles to calculate (default: `c(0.05, 0.25, 0.95)`).
#' @param polygonize Logical. If `TRUE`, vectorize the DOY binary masks using GDAL (`gdal_polygonize.py`). Default is `TRUE`.
#' @param python_exe Path to the Python executable. Required if `polygonize = TRUE`.
#' @param gdal_polygonize_script Path to `gdal_polygonize.py`. Required if `polygonize = TRUE`.
#' @param get_mode Optional custom function to calculate the mode. Default function provided.
#' @param get_median Optional custom function to calculate the median. Default function provided.
#' @param get_quantile Optional custom function to calculate percentiles. Default function provided.
#'@param keep_all_polygons Logical. If `TRUE`, all polygons are kept and a flag column is added indicating which polygons meet the threshold condition.
#'       If `FALSE`, only polygons meeting the condition are saved. Default: `TRUE`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{sf_objects}{Named list of `sf` outputs containing DOY statistics and date components for each input polygon set.}
#'   \item{shp_paths}{Named list of output shapefile paths, including DOY summary shapefiles and any polygonized threshold layers.}
#' }
#'
#' @examples
#' \dontrun{
#' # Case 1: Multiple Otsu-based burned area shapefiles
#' rast <- terra::rast("rbr_doy_stack_1987.tif")
#' polys <- list(
#'   ge50 = "burned_areas_1987_otsu_ge50.shp",
#'   ge100 = "burned_areas_1987_otsu_ge100.shp"
#' )
#' calculate_doy_flags(
#'   raster = rast,
#'   doy_band = 2,
#'   polygons_sf = polys,
#'   output_dir = "output/",
#'   year = 1987,
#'   stats = "both",
#'   doy_thresholds = c(10, 15),
#'   calc_percentiles = TRUE,
#'   percentiles = c(0.05, 0.25, 0.95),
#'   polygonize = TRUE,
#'   python_exe = "/usr/bin/python",
#'   gdal_polygonize_script = "/usr/bin/gdal_polygonize.py"
#' )
#'
#' # Case 2: Single shapefile without polygonization
#' calculate_doy_flags(
#'   raster = rast,
#'   polygons_sf = "burned_areas_1987_origval.shp",
#'   output_dir = "output/",
#'   year = 1987,
#'   stats = "median",
#'   doy_thresholds = 10,
#'   calc_percentiles = FALSE,
#'   polygonize = FALSE
#' )
#' }
#'
#' @note The DOY raster band should be derived from a Landsat composite generated in Google Earth Engine (GEE),
#' following the approach described by Quintero et al. (2025).
#'
#' @references
#' Quintero, N., Viedma, O., Veraverbeke, S., & Moreno, J. M. (2025).
#' *Optimising Regional Fire Severity Mapping Using Pixel-Based Image Compositing*. SSRN 4929831.
#'
#' @importFrom sf st_read st_write st_transform st_crs
#' @importFrom terra crop mask extract classify writeRaster rasterize vect
#' @importFrom data.table as.data.table setnames
#' @importFrom stats median quantile
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom dplyr filter first
#' @importFrom data.table :=
#' @importFrom magrittr %>%
#' @export

utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":="
))

calculate_doy_flags <- function(
    raster,
    doy_band = 2,
    polygons_sf,
    output_dir,
    year,
    doy_thresholds = c(10, 15, 20, 30),
    stats = c("both", "median", "mode"),
    calc_percentiles = TRUE,
    percentiles = c(0.05, 0.25, 0.95),
    polygonize = TRUE,
    keep_all_polygons = TRUE,
    python_exe = NULL,
    gdal_polygonize_script = NULL,
    get_mode = function(x) {
      if (length(x) == 0 || all(is.na(x))) return(NA_real_)
      ux <- unique(na.omit(x))
      ux[which.max(tabulate(match(x, ux)))]
    },
    get_median = function(x) {
      if (length(x) == 0 || all(is.na(x))) return(NA_real_)
      as.numeric(median(x, na.rm = TRUE))
    },
    get_quantile = function(x, prob) {
      if (length(x) == 0 || all(is.na(x))) return(NA_real_)
      as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
    }
) {
  stats <- match.arg(stats)
  sf_outputs <- list(stats = list())
  shp_paths  <- list(stats = list())

  # Organizar entrada
  if (is.character(polygons_sf) && length(polygons_sf) > 1) {
    names(polygons_sf) <- tools::file_path_sans_ext(basename(polygons_sf))
    polygons_sf <- as.list(polygons_sf)
  } else if (!is.list(polygons_sf)) {
    label <- tools::file_path_sans_ext(basename(polygons_sf))
    polygons_sf <- setNames(list(polygons_sf), label)
  } else if (is.null(names(polygons_sf)) || any(names(polygons_sf) == "")) {
    names(polygons_sf) <- sapply(polygons_sf, function(x) {
      if (is.character(x)) tools::file_path_sans_ext(basename(x)) else "sf_input"
    })
  }

  # Procesar cada poligono
  for (label in names(polygons_sf)) {
    shp <- polygons_sf[[label]]
    if (is.character(shp) && file.exists(shp)) {
      shp <- sf::st_read(shp, quiet = TRUE)
    }
    if (!inherits(shp, "sf")) {
      warning("Input invalido para: ", label)
      next
    }

    shp$fire_id <- seq_len(nrow(shp))
    rbr_crop <- terra::mask(terra::crop(raster, shp), shp)
    doy_raster <- rbr_crop[[doy_band]]

    if (!sf::st_crs(shp) == terra::crs(doy_raster)) {
      shp <- sf::st_transform(shp, terra::crs(doy_raster))
    }

    vals <- terra::extract(doy_raster, terra::vect(shp), bind = FALSE)
    vals_dt <- data.table::as.data.table(vals)
    names(vals_dt)[2] <- "doy"
    vals_dt <- vals_dt[!is.na(doy)]
    if (!"doy" %in% names(vals_dt) || nrow(vals_dt) == 0) next

    # Calcular estadisticas
    mode_by_fire   <- vals_dt[, .(mode_doy   = get_mode(doy)), by = ID]
    median_by_fire <- vals_dt[, .(median_doy = get_median(doy)), by = ID]
    merge_list <- list(mode_by_fire, median_by_fire)

    if (calc_percentiles && length(percentiles) > 0) {
      for (p in percentiles) {
        pname <- paste0("q", round(p * 100))
        q_by_fire <- vals_dt[, .(val = get_quantile(doy, p)), by = ID]
        data.table::setnames(q_by_fire, "val", paste0(pname, "_doy"))
        merge_list <- c(merge_list, list(q_by_fire))
      }
    }

    doy_stats <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), merge_list)
    shp_stats <- merge(shp, doy_stats, by.x = "fire_id", by.y = "ID", all.x = TRUE)

    if (!"area_m2" %in% names(shp_stats)) {
      shp_stats$area_m2 <- as.numeric(sf::st_area(shp_stats))
    }
    shp_stats$area_m2 <- round(shp_stats$area_m2, 2)

        # Convertir DOY a fecha
    date_vars <- c("mode_doy", "median_doy")
    if (calc_percentiles && length(percentiles) > 0) {
      q_names <- paste0("q", round(percentiles * 100), "_doy")
      date_vars <- c(date_vars, q_names)
    }

    for (var in date_vars) {
      shp_stats[[paste0(var, "_date")]]  <- as.Date(shp_stats[[var]] - 1, origin = paste0(year, "-01-01"))
      shp_stats[[paste0(var, "_day")]]   <- as.integer(format(shp_stats[[paste0(var, "_date")]], "%d"))
      shp_stats[[paste0(var, "_month")]] <- as.integer(format(shp_stats[[paste0(var, "_date")]], "%m"))
      shp_stats[[paste0(var, "_year")]]  <- as.integer(format(shp_stats[[paste0(var, "_date")]], "%Y"))
    }

    label_clean <- sub(".*otsu_", "", label)
    short_label <- label_clean

    # LOOP de umbrales
    for (doy_thresh in doy_thresholds) {
      label_clean <- sub(".*otsu_", "", label)
      short_label <- label_clean

      # --- MEDIAN ---
      if (stats %in% c("median", "both")) {
        stat_name <- "median"
        flag_mask <- doy_raster >= (terra::rasterize(terra::vect(shp_stats), doy_raster, field = "median_doy") - doy_thresh)
        flag_mask1 <- terra::classify(flag_mask, cbind(0, NA))

        out_tif <- file.path(output_dir, sprintf("DOY_%d_%s_%s.tif", doy_thresh, short_label, stat_name))
        if (file.exists(out_tif)) file.remove(out_tif)
        terra::writeRaster(flag_mask1, out_tif, overwrite = TRUE, datatype = "INT1U")

        if (polygonize && !is.null(python_exe) && !is.null(gdal_polygonize_script)) {
          shp_mask <- sub("\\.tif$", ".shp", out_tif)

          # --- Asignar flag_doy ---
          vals_shp <- terra::extract(flag_mask1, terra::vect(shp_stats), fun = max, na.rm = TRUE, ID = FALSE)
          shp_stats$flag_doy <- ifelse(vals_shp == 1, 1L, 0L)

          # --- Filtrar si el usuario quiere ---
          shp_stats_out <- if (keep_all_polygons) shp_stats else shp_stats[shp_stats$flag_doy == 1, ]

          if (nrow(shp_stats_out) > 0) {
            suppressWarnings(sf::st_write(shp_stats_out, shp_mask, delete_dsn = TRUE, quiet = TRUE))
            sf_outputs$stats[[paste0("thresh_", doy_thresh, "_", stat_name, "_", label)]] <- shp_stats_out
            shp_paths$stats[[paste0("thresh_", doy_thresh, "_", stat_name, "_", label)]] <- shp_mask
          }
        }
      }

      # --- MODE ---
      if (stats %in% c("mode", "both")) {
        stat_name <- "mode"
        flag_mask <- doy_raster >= (terra::rasterize(terra::vect(shp_stats), doy_raster, field = "mode_doy") - doy_thresh)
        flag_mask1 <- terra::classify(flag_mask, cbind(0, NA))

        out_tif <- file.path(output_dir, sprintf("DOY_%d_%s_%s.tif", doy_thresh, short_label, stat_name))
        if (file.exists(out_tif)) file.remove(out_tif)
        terra::writeRaster(flag_mask1, out_tif, overwrite = TRUE, datatype = "INT1U")

        if (polygonize && !is.null(python_exe) && !is.null(gdal_polygonize_script)) {
          shp_mask <- sub("\\.tif$", ".shp", out_tif)

          # --- Asignar flag_doy ---
          vals_shp <- terra::extract(flag_mask1, terra::vect(shp_stats), fun = max, na.rm = TRUE, ID = FALSE)
          shp_stats$flag_doy <- ifelse(vals_shp == 1, 1L, 0L)

          # --- Filtrar si el usuario quiere ---
          shp_stats_out <- if (keep_all_polygons) shp_stats else shp_stats[shp_stats$flag_doy == 1, ]

          if (nrow(shp_stats_out) > 0) {
            suppressWarnings(sf::st_write(shp_stats_out, shp_mask, delete_dsn = TRUE, quiet = TRUE))
            sf_outputs$stats[[paste0("thresh_", doy_thresh, "_", stat_name, "_", label)]] <- shp_stats_out
            shp_paths$stats[[paste0("thresh_", doy_thresh, "_", stat_name, "_", label)]] <- shp_mask
          }
        }
      }
    }

  }

  return(list(
    sf_objects = sf_outputs,
    shp_paths  = shp_paths
  ))
}

