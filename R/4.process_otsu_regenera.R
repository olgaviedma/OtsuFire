#' Detect Post-Fire Regeneration Using Otsu or Fixed Thresholds
#'
#' @description
#' This function detects vegetation regeneration signals using negative RBR or dNBR values.
#' By default (`use_fixed_threshold = FALSE`), the function identifies regeneration areas where
#' RBR or dNBR values are below 0. If `use_fixed_threshold = TRUE`, a fixed threshold value
#' (e.g., -100) is applied instead.
#'
#' Optionally, percentile-based clipping can be applied before thresholding, and tiling
#' can be used to handle large rasters. The function outputs binary rasters and shapefiles,
#' optionally merged into a single shapefile if `bind_all = TRUE`.
#'
#' If a historical burn raster (`rbr_date`) is provided, it will be used to mask out areas
#' previously burned (i.e., where RBR \\ge 0), so that regeneration is only detected in those
#' affected by the most recent fire.
#'
#' @name process_otsu_regenera
#' @rdname process_otsu_regenera
#'
#' @param rbr_post Path to post-fire RBR or dNBR raster. Alternatively, provide `nbr_pre_path` and `nbr_post_path`.
#'                 If a named list is provided, names should match the format `"P1"`, `"P2"`, etc.
#' @param nbr_pre_path Path to pre-fire NBR raster (for index computation).
#' @param nbr_post_path Path to post-fire NBR raster (for index computation).
#' @param rbr_date Optional path to an RBR raster used to mask previously burned areas (RBR \\ge 0).
#' @param output_dir Directory where outputs (rasters, shapefiles, plots) will be saved.
#' @param python_exe Path to the Python executable (used to call GDAL).
#' @param gdal_polygonize_script Path to `gdal_polygonize.py` script.
#' @param n_rows,n_cols Number of rows and columns for tiling. Ignored if `tile = FALSE`.
#' @param tile_overlap Overlap size (in meters) between tiles. Used for seamless polygonization.
#' @param tile Logical. If `TRUE`, raster is tiled before polygonization.
#' @param index_type Index type to compute: either `"RBR"` or `"dNBR"`. Ignored if `rbr_post` is provided.
#' @param trim_percentiles Data frame with columns `min` and `max`, defining percentile ranges
#'                         to clip negative RBR values before thresholding.
#'                         **Note:** This parameter is currently not used unless percentile-based thresholding is implemented.
#' @param bind_all Logical. If `TRUE`, all shapefiles from each threshold are merged into one combined shapefile.
#' @param regen_year Integer vector specifying the number of years after the fire to consider for regeneration (e.g., `c(2)`).
#' @param fire_year Integer indicating the fire year. Used for output naming and labeling.
#' @param use_fixed_threshold Logical. If `TRUE`, applies `fixed_threshold_value` instead of using the default threshold (< 0).
#' @param fixed_threshold_value Numeric threshold to use when `use_fixed_threshold = TRUE`. Default is `-100`.
#' @param output_format Character. Output format for saved files: `"shp"` for ESRI Shapefile (default) or `"geojson"` for GeoJSON.
#'
#' @return A named list of results per regeneration year (`P1`, `P2`, etc.) with:
#' \describe{
#'   \item{`raster`}{Path to the binary regeneration raster for that year.}
#'   \item{`shapefile`}{Path to the shapefile with polygons derived from the raster tiles for that year.}
#' }
#'
#' If `bind_all = TRUE`, an additional item named `combined` is returned, which contains the path to a single shapefile:
#' \item{`combined`}{Path to the shapefile with all valid regeneration polygons combined across all years. The filename includes all `P` labels used (e.g., `_P1P2_`).}
#'
#' @examples
#' \dontrun{
#' # Regeneration detection using RBR raster and default threshold (< 0)
#' process_otsu_regenera(
#'   rbr_post = list(P2 = "data/RBR_1986.tif"),
#'   output_dir = "output/regenera",
#'   python_exe = "/usr/bin/python3",
#'   gdal_polygonize_script = "/usr/bin/gdal_polygonize.py",
#'   fire_year = 1984,
#'   regen_year = c(2),
#'   use_fixed_threshold = FALSE,
#'   output_format = c("geojson")
#' )
#'
#' # Same detection but using a fixed threshold of -150
#' and using tiling if the input raster is very large
#' process_otsu_regenera(
#'   rbr_post = list(P2 = "data/RBR_1986.tif"),
#'   output_dir = "output/regenera",
#'   python_exe = "/usr/bin/python3",
#'   gdal_polygonize_script = "/usr/bin/gdal_polygonize.py",
#'   fire_year = 1984,
#'   regen_year = c(2),
#'   use_fixed_threshold = TRUE,
#'   fixed_threshold_value = -150
#'   bind_all = FALSE,
#'   n_rows = 2,
#'   n_cols = 3,
#'   tile_overlap = 1000,
#'   tile = TRUE,
#'   output_format = c("shp")
#' )
#' }
#'
#' @importFrom terra rast mask crop writeRaster ext ifel values ncell
#' @importFrom sf st_read st_write st_geometry_type st_transform st_as_binary st_make_valid st_union st_cast st_as_sf
#' @importFrom tools file_path_sans_ext
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom dplyr filter first
#' @importFrom data.table :=
#' @importFrom stringr str_detect str_replace
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
    fixed_threshold_value = -100,
    output_format = c("shp", "geojson")
) {

  # Validar formato
  output_format <- match.arg(output_format, choices = c("shp", "geojson"))


  tile <- as.logical(tile)[1]
  if (is.na(tile)) stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")

  if (!use_fixed_threshold && (!is.data.frame(trim_percentiles) || !all(c("min", "max") %in% names(trim_percentiles)))) {
    stop("`trim_percentiles` must be a data.frame with columns `min` and `max`.")
  }

  if (use_fixed_threshold && is.null(fixed_threshold_value)) {
    stop("`fixed_threshold_value` must be provided when `use_fixed_threshold = TRUE`.")
  }

  if (use_fixed_threshold && !missing(trim_percentiles)) {
    message("Note: 'trim_percentiles' is ignored when 'use_fixed_threshold = TRUE'.")
  }


  results_all <- list()

  all_polys_list_total <- list()
  all_labels <- c()

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

    if (!is.null(rbr_date)) {
      r_prev <- terra::rast(rbr_date)[[1]]
      r_prev <- terra::resample(r_prev, r, method = "near")
      mask_r <- terra::ifel(r_prev >= 0, 1, NA)
      r <- terra::mask(r, mask_r)
    }

    if (use_fixed_threshold) {
      binary_raster <- terra::ifel(r < fixed_threshold_value, 1, NA)
      real_threshold <- fixed_threshold_value
    } else {
      binary_raster <- terra::ifel(r < 0, 1, NA)
      real_threshold <- 0
    }

    folder_path <- output_dir
    threshold_log <- data.frame(Label = character(), RealThreshold = numeric(), stringsAsFactors = FALSE)
    results <- list()

    label_suffix <- if (use_fixed_threshold)
      paste0(regeneration_year, "_P", year_offset, "_thresh", abs(fixed_threshold_value))
    else
      paste0(regeneration_year, "_P", year_offset, "_otsu")

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
          if (file.exists(shp_path)) tile_shapefiles <- c(tile_shapefiles, shp_path)
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


        # Determinar extension segun formato
        ext <- if (output_format == "geojson") ".geojson" else ".shp"

        # Construir ruta de salida
        final_file <- file.path(folder_path, sprintf("burned_areas_%s_regenera_%s%s", fire_year, label_suffix, ext))

        # Eliminar archivos previos si es necesario
        if (output_format == "shp") {
          shp_base <- tools::file_path_sans_ext(final_file)
          shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
          for (ext_i in shp_exts) {
            f <- paste0(shp_base, ext_i)
            if (file.exists(f)) file.remove(f)
          }
        } else {
          if (file.exists(final_file)) file.remove(final_file)
        }

        # Guardar archivo
        sf::st_write(polys, final_file, append = FALSE, quiet = TRUE)

        # Mensaje de confirmacion
        message("Saved regeneration layer: ", final_file)


        results[[label_suffix]] <- list(raster = raster_name, shapefile = final_file)

        if (bind_all) {
          polys$label <- label_suffix
          all_polys_list_total[[label_suffix]] <- polys
          all_labels <- union(all_labels, paste0("P", year_offset))
        }

        tile_prefix <- sprintf("tile_regenera_%s_", label_suffix)
        tile_files <- list.files(folder_path, pattern = paste0("^", tile_prefix, ".*"), full.names = TRUE)
        sapply(tile_files, function(f) tryCatch(file.remove(f), error = function(e) NULL))
      }
    }

    log_file <- file.path(output_dir, sprintf("otsu_thresholds_burned_%s_regen_%s.txt", fire_year, label_suffix))
    write_header <- !file.exists(log_file) || file.info(log_file)$size == 0
    suppressWarnings(write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE, append = TRUE, col.names = write_header))

    results_all[[paste0("P", year_offset)]] <- results
  }

  # BLOQUE COMBINADO FINAL
  if (bind_all && length(all_polys_list_total) > 0) {
    message("\U0001F517 Binding all regenera polygons into a single shapefile...")
    combined_sf <- do.call(rbind, all_polys_list_total)
    combined_sf <- sf::st_make_valid(combined_sf)
    combined_union <- sf::st_union(combined_sf)
    combined_sf_final <- sf::st_cast(combined_union, "POLYGON")
    combined_sf_final <- sf::st_as_sf(combined_sf_final)
    combined_sf_final$ID <- seq_len(nrow(combined_sf_final))
    combined_sf_final$otsu_threshold <- real_threshold
    combined_sf_final$regen_year <- paste(all_labels, collapse = "")

    # Validar formato de salida
    output_format <- match.arg(output_format, choices = c("shp", "geojson"))

    # Determinar extension
    ext <- if (output_format == "geojson") ".geojson" else ".shp"

    # Construir nombre de archivo
    combined_path <- file.path(
      output_dir,
      sprintf("regenera_combined_burned_%s_%s_thresh%s%s",
              fire_year,
              paste(all_labels, collapse = ""),
              abs(real_threshold),
              ext)
    )

    # Eliminar archivo(s) existente(s) si corresponde
    if (output_format == "shp") {
      shp_base <- tools::file_path_sans_ext(combined_path)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext_i in shp_exts) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }
    } else {
      if (file.exists(combined_path)) file.remove(combined_path)
    }

    # Guardar archivo
    sf::st_write(combined_sf_final, combined_path, append = FALSE, quiet = TRUE)

    # Guardar ruta en resultado
    results_all$combined <- combined_path

    # Mensaje de confirmacion
    message("Saved combined regeneration layer: ", combined_path)

  }


  return(results_all)
}



