#' Detect Post-Fire Regeneration Using Otsu and Percentile-Based Clipping
#'
#' This function detects vegetation regeneration signals using negative RBR or dNBR values.
#' It applies percentile-based clipping to the raster, rescales values, and uses a smoothed
#' Otsu threshold to extract regeneration areas. Optionally, the raster can be tiled before
#' polygonization using GDAL. All outputs are saved as binary rasters and shapefiles.
#'
#' You can also mask out previously burned areas by providing a historical RBR raster
#' (`rbr_date`). If `bind_all = TRUE`, the function merges all resulting shapefiles
#' across thresholds into one dissolved shapefile with unique polygon IDs.
#'
#' If `use_fixed_threshold = TRUE`, a fixed threshold is applied instead of calculating it from Otsu.
#'
#' @param rbr_post Path to post-fire RBR or dNBR raster. Alternatively, provide `nbr_pre_path` and `nbr_post_path`.
#'                 If a named list is provided, names should match the format `"P1"`, `"P2"`, etc.
#' @param nbr_pre_path Path to pre-fire NBR raster (for index computation).
#' @param nbr_post_path Path to post-fire NBR raster (for index computation).
#' @param rbr_date Optional path to an RBR raster used to mask previously burned areas (RBR < 0).
#' @param output_dir Directory where outputs (rasters, shapefiles, plots) will be saved.
#' @param python_exe Path to the Python executable (used to call GDAL).
#' @param gdal_polygonize_script Path to `gdal_polygonize.py` script.
#' @param n_rows,n_cols Number of rows and columns for tiling. Ignored if `tile = FALSE`.
#' @param tile_overlap Overlap size (in meters) between tiles. Used for seamless polygonization.
#' @param tile Logical. If `TRUE`, raster is tiled before polygonization.
#' @param index_type Index type to compute: either `"RBR"` or `"dNBR"`. Ignored if `rbr_post` is provided.
#' @param trim_percentiles Data frame with columns `min` and `max`, defining percentile ranges
#'                         to clip negative RBR values before thresholding (only used if `use_fixed_threshold = FALSE`).
#' @param bind_all Logical. If `TRUE`, all shapefiles from each threshold are merged into one combined shapefile.
#' @param regen_year Integer vector specifying the number of years after the fire to consider for regeneration (e.g., `c(2)`).
#' @param fire_year Integer indicating the fire year. Used for output naming and labeling.
#' @param use_fixed_threshold Logical. If `TRUE`, applies a fixed threshold instead of Otsu.
#' @param fixed_threshold_value Numeric threshold to use when `use_fixed_threshold = TRUE`. Default is `-150`.
#'
#' @return A named list with entries for each regeneration year (`"P1"`, `"P2"`, etc.), each containing:
#' \describe{
#'   \item{raster}{Path to binary raster with detected regeneration.}
#'   \item{shapefile}{Path to polygonized shapefile of regeneration areas.}
#'   \item{combined}{(if `bind_all = TRUE`) Path to the merged shapefile.}
#' }
#'
#' @examples
#' \dontrun{
#' # Regeneration detection using precomputed RBR and masking out prior burns
#' process_otsu_regenera(
#'   rbr_post = "data/RBR_1986.tif",
#'   rbr_date = "data/RBR_1985.tif",
#'   output_dir = "output/regenera",
#'   fire_year = 1984,
#'   regen_year = 2,
#'   tile = TRUE,
#'   n_rows = 2,
#'   n_cols = 3,
#'   tile_overlap = 1000,
#'   bind_all = TRUE,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py",
#'   trim_percentiles = data.frame(
#'     min = c(0.01, 0.005),
#'     max = c(0.99, 0.995)
#'   )
#' )
#' }
#'
#' @importFrom terra rast mask crop writeRaster ext ifel values
#' @importFrom sf st_read st_write st_geometry_type st_transform st_as_binary st_make_valid st_union st_cast st_as_sf
#' @importFrom tools file_path_sans_ext
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom dplyr filter first
#' @importFrom data.table :=
#' @export

utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":="
))

process_otsu_regenera <- function(
    rbr_post = NULL,
    nbr_pre_path = NULL,
    nbr_post_path = NULL,
    rbr_date = NULL,
    output_dir,
    python_exe,
    gdal_polygonize_script,
    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    tile = TRUE,
    index_type = "RBR",
    trim_percentiles = data.frame(min = c(0.01, 0.005), max = c(0.99, 0.995)),
    bind_all = FALSE,
    regen_year = c(2),
    fire_year = NULL,
    use_fixed_threshold = FALSE,
    fixed_threshold_value = -100
) {
  tile <- as.logical(tile)[1]
  if (is.na(tile)) stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")

  if (!use_fixed_threshold && (!is.data.frame(trim_percentiles) || !all(c("min", "max") %in% names(trim_percentiles)))) {
    stop("`trim_percentiles` must be a data.frame with columns `min` and `max`.")
  }

  results_all <- list()

  for (year_offset in regen_year) {
    label_key <- paste0("P", year_offset)

    regen_year_value <- if (!is.null(rbr_post)) {
      rbr_path <- if (is.list(rbr_post)) {
        if (!label_key %in% names(rbr_post)) {
          stop(paste("No se encontro raster para", label_key, "en rbr_post"))
        }
        rbr_post[[label_key]]
      } else {
        rbr_post
      }
      as.numeric(gsub("\\D", "", basename(rbr_path)))
    } else if (!is.null(nbr_post_path)) {
      as.numeric(basename(dirname(nbr_post_path)))
    } else {
      stop("You must provide either 'rbr_post' or both 'nbr_pre_path' and 'nbr_post_path'.")
    }

    regeneration_year <- fire_year + year_offset

    r <- if (!is.null(rbr_post)) {
      terra::rast(rbr_path)[[1]]
    } else {
      nbr_pre <- terra::rast(nbr_pre_path)
      nbr_post <- terra::rast(nbr_post_path)
      if (tolower(index_type) == "dnbr") {
        (nbr_pre - nbr_post) * 1000
      } else if (tolower(index_type) == "rbr") {
        (nbr_pre - nbr_post) * 1000 / (nbr_pre + 1.001)
      } else {
        stop("`index_type` must be either 'RBR' or 'dNBR'")
      }
    }


    # 1. CREAR LA MÁSCARA si hay rbr_date
    if (!is.null(rbr_date)) {
      r_prev <- terra::rast(rbr_date)[[1]]

      # Alinear antes de aplicar mascara
      r_prev <- terra::resample(r_prev, r, method = "near")

      # Crear mascara: 1 si rbr_date >= 0 (area valida para regenerar)
      mask_r <- terra::ifel(r_prev >= 0, 1, NA)

      # Aplicar la mascara
      r <- terra::mask(r, mask_r)
    }

      # 2. Binarizamos directamente el r enmascarado
    binary_raster <- terra::ifel(r < fixed_threshold_value, 1, NA)

    folder_path <- output_dir
    threshold_log <- data.frame(Label = character(), RealThreshold = numeric(), stringsAsFactors = FALSE)
    all_polys_list <- list()
    results <- list()

    label_suffix <- if (use_fixed_threshold)
      paste0(regeneration_year, "_P", year_offset, "_thresh", abs(fixed_threshold_value))
    else
      paste0(regeneration_year, "_P", year_offset, "_otsu")

        real_threshold <- fixed_threshold_value

    raster_name <- file.path(output_dir, paste0("raster_", fire_year, "_regenera_", label_suffix, ".tif"))
    terra::writeRaster(binary_raster, raster_name, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=-9999"))

    tile_shapefiles <- c()
    if (tile) {
      ext_r <- terra::ext(binary_raster)
      tile_width <- (ext_r[2] - ext_r[1]) / n_cols
      tile_height <- (ext_r[4] - ext_r[3]) / n_rows
      count <- 1
      for (i in 0:(n_rows - 1)) {
        for (j in 0:(n_cols - 1)) {
          xmin_tile <- max(ext_r[1] + j * tile_width - tile_overlap, ext_r[1])
          xmax_tile <- min(ext_r[1] + (j + 1) * tile_width + tile_overlap, ext_r[2])
          ymin_tile <- max(ext_r[3] + i * tile_height - tile_overlap, ext_r[3])
          ymax_tile <- min(ext_r[3] + (i + 1) * tile_height + tile_overlap, ext_r[4])
          tile_crop <- terra::crop(binary_raster, terra::ext(xmin_tile, xmax_tile, ymin_tile, ymax_tile))
          if (all(is.na(terra::values(tile_crop)))) {
            message("Skipping tile ", count, ": all values NA")
            count <- count + 1
            next
          }
          tile_path <- file.path(folder_path, sprintf("tile_regenera_%s_%d.tif", label_suffix, count))
          shp_path <- file.path(folder_path, sprintf("tile_regenera_%s_%d.shp", label_suffix, count))
          terra::writeRaster(tile_crop, tile_path, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
          system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tile_path}" -f "ESRI Shapefile" "{shp_path}" DN'))
          if (file.exists(tile_path)) file.remove(tile_path)
          if (file.exists(shp_path)) {
            tile_shapefiles <- c(tile_shapefiles, shp_path)
          }
          count <- count + 1
        }
      }

      if (length(tile_shapefiles) > 0) {
        polys <- do.call(rbind, lapply(tile_shapefiles, function(x) {
          p <- sf::st_read(x, quiet = TRUE)
          p[p$DN == 1, ]
        }))
        polys <- polys[!duplicated(sf::st_as_binary(sf::st_geometry(polys))), ]
        polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON", "MULTIPOLYGON"), ]
        polys$otsu_threshold <- real_threshold
        polys$regen_year <- paste0("P", year_offset)
        polys <- sf::st_transform(polys, 3035)
        names(polys) <- abbreviate(names(polys), minlength = 10)
        final_shp <- file.path(folder_path, sprintf("burned_areas_%s_regenera_%s.shp", fire_year, label_suffix))

        sf::st_write(polys, final_shp, append = FALSE)
        if (bind_all) {
          polys$label <- label_suffix
          all_polys_list[[label_suffix]] <- polys
        }
        results[[label_suffix]] <- list(raster = raster_name, shapefile = final_shp)
        tile_prefix <- sprintf("tile_regenera_%s_", label_suffix)
        tile_files <- list.files(folder_path, pattern = paste0("^", tile_prefix, ".*"), full.names = TRUE)
        sapply(tile_files, function(f) tryCatch(file.remove(f), error = function(e) NULL))

      }
    }

    log_file <- file.path(output_dir, sprintf("otsu_thresholds_burned_%s_regen_%s.txt", fire_year, label_suffix))
    write_header <- !file.exists(log_file) || file.info(log_file)$size == 0
    suppressWarnings(write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE, append = TRUE, col.names = write_header))

    if (bind_all && length(all_polys_list) > 0) {
      message("\U0001F517 Binding all regenera polygons into a single shapefile...")
      combined_sf <- do.call(rbind, all_polys_list)
      combined_sf <- sf::st_make_valid(combined_sf)
      combined_union <- sf::st_union(combined_sf)
      combined_sf_final <- sf::st_cast(combined_union, "POLYGON")
      combined_sf_final <- sf::st_as_sf(combined_sf_final)
      combined_sf_final$ID <- seq_len(nrow(combined_sf_final))
      combined_sf_final$otsu_threshold <- real_threshold
      combined_sf_final$regen_year <- paste0("P", year_offset)
      combined_path <- file.path(folder_path, sprintf("regenera_combined_burned_%s_P%s_thresh%s.shp", fire_year, year_offset, abs(fixed_threshold_value)))
      sf::st_write(combined_sf_final, combined_path, append = FALSE)
      results$combined <- combined_path
    }

    results_all[[paste0("P", year_offset)]] <- results
  }

  return(results_all)
}



