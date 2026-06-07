#' Process Burned Area Rasters Using Otsu Thresholding or Percentile Clipping
#'
#' @description
#' This function computes binary burned area masks from a severity index raster (RBR or dNBR)
#' using Otsu's thresholding or percentile clipping. The segmentation can be applied
#' to the full raster or stratified by:
#' - CORINE land cover classes (`corine_raster_path`)
#' - WWF ecoregions (`ecoregion_shapefile_path`)
#' - or the intersection of both (CORINE × Ecoregion)
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
#'   - Polygon shapefiles (no dissolve at write-time; attributes preserved)
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
#' @param dissolve Ignored at write-time. Writing never dissolves; if you need to dissolve, do it explicitly before writing.
#' @param vectorize Logical. If TRUE, polygonize the burned mask; if FALSE, keep only raster output.
#' @param write_raster Logical. If TRUE, writes the binary burned-area raster to disk.
#' @param burnable_mask Optional `terra::SpatRaster` (0/1). Cells with 0/NA are treated as non-burnable and masked out before binarization.

#'
#'
#' @return Invisibly returns a named list with small summaries per output written.
#'
#' @importFrom terra rast crop mask ext resample ifel values freq writeRaster ncell classify vect rasterize compareGeom
#' @importFrom sf st_read st_write st_transform st_geometry_type st_as_binary st_make_valid st_is_valid st_drop_geometry st_centroid st_intersection st_as_sf st_crs st_cast st_join
#' @importFrom stats quantile
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par plot
#' @importFrom dplyr all_of any_of where first group_by summarise select bind_rows left_join
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom stringr str_detect str_replace str_extract
#' @importFrom glue glue
process_otsu_rasters_ <- function(
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
    corine_classes = NULL,
    ecoregion_shapefile_path = NULL,
    ecoregion_field = NULL,
    ecoregion_classes = NULL,
    segment_by_intersection = FALSE,
    min_otsu_threshold_value = NULL,
    python_exe,
    gdal_polygonize_script,
    gdalwarp_path = "gdalwarp",
    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    tile = TRUE,
    index_type = "RBR",
    vectorize = FALSE,
    write_raster = TRUE,
    output_format = c("shp", "geojson"),
    dissolve = FALSE,          # ignorado al escribir (siempre NO dissolve)
    burnable_mask = NULL       # terra::SpatRaster 0/1
) {
  
  # AS07 (0.5.0): planar geometry. This function relies on s2 = FALSE for its
  # st_make_valid / st_intersection work. The historical top-level
  # `sf::sf_use_s2(FALSE)` in R/*.R is DEAD in installed-package mode (top-level
  # R/ expressions run at build time, not at library() load time), so the
  # setting only took effect under pkgload::load_all(). Scope it locally and
  # restore the caller's previous value on exit so the planar setting is
  # guaranteed in BOTH modes without permanently contaminating the session.
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  output_format <- match.arg(output_format, choices = c("shp", "geojson"))

  if (missing(output_dir) || is.null(output_dir)) stop("'output_dir' must be provided.")
  if (missing(python_exe) || is.null(python_exe)) stop("'python_exe' must be provided.")
  if (missing(gdal_polygonize_script) || is.null(gdal_polygonize_script)) stop("'gdal_polygonize_script' must be provided.")
  
  tile <- as.logical(tile)[1]
  if (is.na(tile)) stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  figures_dir <- file.path(output_dir, "FIGURES")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  
  # ---- Infer year if not provided ----
  if (is.null(year)) {
    if (!is.null(raster_path)) {
      raster_base <- tools::file_path_sans_ext(basename(raster_path))
      year_match <- stringr::str_extract(raster_base, "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    } else {
      year_match <- stringr::str_extract(basename(normalizePath(output_dir)), "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    }
  }
  
  # ---- Normalize thresholds ----
  if (!is.null(otsu_thresholds)) {
    otsu_thresholds <- suppressWarnings(as.numeric(otsu_thresholds))
    otsu_thresholds <- otsu_thresholds[is.finite(otsu_thresholds)]
    if (!length(otsu_thresholds)) otsu_thresholds <- NULL
  }
  
  do_orig <- isTRUE(use_original)
  do_otsu <- !is.null(otsu_thresholds) && length(otsu_thresholds) > 0
  
  # ---- Ensure only one thresholding method ----
  if (!is.null(trim_percentiles) && !is.null(otsu_thresholds)) {
    stop("You cannot use 'trim_percentiles' and 'otsu_thresholds' at the same time. Use only one.")
  }
  
  # ---- Determine raster source ----
  if (!is.null(nbr_pre_path) && !is.null(nbr_post_path)) {
    if (!file.exists(nbr_pre_path)) stop("Pre-fire NBR raster not found: ", nbr_pre_path)
    if (!file.exists(nbr_post_path)) stop("Post-fire NBR raster not found: ", nbr_post_path)
    
    nbr_pre  <- terra::rast(nbr_pre_path)
    nbr_post <- terra::rast(nbr_post_path)
    
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
    if (!file.exists(raster_path)) stop("Raster not found: ", raster_path)
    r <- terra::rast(raster_path)[[1]]
  } else {
    stop("You must provide either 'raster_path' or both 'nbr_pre_path' and 'nbr_post_path'.")
  }
  
  # ---------------------- helper: write vectors (NO DISSOLVE) ----------------------
  .write_vector <- function(sfobj, outfile, format = c("shp", "geojson")) {
    format <- match.arg(format)
    if (!inherits(sfobj, "sf")) stop("sfobj must be an sf object.")
    if (nrow(sfobj) == 0) stop("sfobj has 0 rows.")
    
    if (format == "shp") {
      names(sfobj) <- make.unique(substr(names(sfobj), 1, 10))
      shp_base <- tools::file_path_sans_ext(outfile)
      for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }
    } else {
      if (file.exists(outfile)) file.remove(outfile)
    }
    
    sf::st_write(sfobj, outfile, append = FALSE, quiet = TRUE)
    invisible(TRUE)
  }
  
  # ====================== HELPER: single threshold ====================
  process_single_threshold <- function(
    r_input,
    pmin = NULL,
    pmax = NULL,
    otsu_min = NULL,
    label_suffix = "",
    corine_class = NULL,
    corine_year = NULL,
    ecoregion_name = NULL,
    ecoregion_field = NULL,
    cor_eco_name = NULL,
    cor_eco_field = NULL,
    min_otsu_threshold_value = NULL,
    output_dir,
    year = NULL,
    tile = TRUE,
    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    python_exe,
    gdal_polygonize_script,
    burnable_mask = NULL,
    write_raster = TRUE,
    vectorize = TRUE,
    output_format = c("shp", "geojson")
  ) {
    
    output_format <- match.arg(output_format, choices = c("shp", "geojson"))
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    fig_dir <- file.path(output_dir, "FIGURES")
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
    
    tile_dir <- file.path(output_dir, "temp_tiles")
    if (!dir.exists(tile_dir)) dir.create(tile_dir, recursive = TRUE)
    
    if (is.null(year)) {
      year_match <- stringr::str_extract(basename(output_dir), "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    }
    
    label <- if (nzchar(label_suffix)) label_suffix else "default"
    
    # --------- Filter (otsu_min / percentiles / original) ----------
    if (!is.null(otsu_min)) {
      r_filtered <- terra::ifel(r_input >= otsu_min, r_input, NA)
      if (!nzchar(label_suffix)) label <- paste0("ge", otsu_min)
      
    } else if (!is.null(pmin) && !is.null(pmax)) {
      # (FIX) Percentiles via terra::quantile (block-wise; no terra::values)
      qq <- terra::quantile(r_input, probs = c(pmin, pmax), na.rm = TRUE)
      min_val <- as.numeric(qq[1, 1])
      max_val <- as.numeric(qq[2, 1])
      if (!is.finite(min_val) || !is.finite(max_val) || min_val >= max_val) {
        return(list(raster = NULL, polys = NULL, log = NULL, tile_dir = tile_dir))
      }
      
      r_filtered <- terra::ifel(r_input >= min_val & r_input <= max_val, r_input, NA)
      
      if (!nzchar(label_suffix)) {
        label <- paste0(
          "P", formatC(pmin * 100, width = 2, flag = "0"),
          "_toP", formatC(pmax * 100, width = 2, flag = "0")
        )
      }
      
    } else {
      r_filtered <- r_input
      if (!nzchar(label_suffix)) label <- "original_values"
    }
    
    # ---- (FIX) Apply burnable mask BEFORE minmax/hist/otsu and before binarization ----
    if (!is.null(burnable_mask)) {
      if (!terra::compareGeom(burnable_mask, r_filtered, stopOnError = FALSE)) {
        burnable_mask <- terra::resample(burnable_mask, r_filtered, method = "near")
      }
      burn01 <- terra::ifel(!is.na(burnable_mask) & burnable_mask == 1, 1, NA)
      r_filtered <- terra::mask(r_filtered, burn01)
      r_filtered <- terra::trim(r_filtered)
    }
    
    # --------- Validate & normalize (SIN values(): por bloques) ----------
    terra::minmax(r_filtered)
    mm <- terra::minmax(r_filtered)
    
    min_val <- mm[1, 1]
    max_val <- mm[2, 1]
    range_val <- max_val - min_val
    
    if (!is.finite(min_val) || !is.finite(max_val) || !is.finite(range_val) || range_val == 0) {
      return(list(raster = NULL, polys = NULL, log = NULL, tile_dir = tile_dir))
    }
    
    r_rescaled <- (r_filtered - min_val) / range_val * 255
    
    # --------- Otsu (histograma por bloques con terra::hist) ----------
    # Block 9b: OtsuSeg moved to Suggests; resolved at call time.
    if (!requireNamespace("OtsuSeg", quietly = TRUE)) {
      stop("Package 'OtsuSeg' is required for the Otsu-smoothed-histogram ",
           "step in process_otsu_rasters_(). Install it with ",
           "install.packages('OtsuSeg').", call. = FALSE)
    }
    .smooth_histogram_fun <- getExportedValue("OtsuSeg", "smooth_histogram")
    .otsu_threshold_smoothed_fun <- getExportedValue("OtsuSeg",
                                                     "otsu_threshold_smoothed")
    hist_values <- terra::hist(r_rescaled, breaks = 256, plot = FALSE)
    smoothed_counts <- .smooth_histogram_fun(hist_values$counts)
    smoothed_counts[is.na(smoothed_counts)] <- 0

    threshold_value_smoothed <- .otsu_threshold_smoothed_fun(
      smoothed_counts, hist_values$mids
    )
    real_threshold <- (threshold_value_smoothed / 255) * range_val + min_val
    
    # --------- Apply min threshold ----------
    min_applied <- FALSE
    if (!is.null(min_otsu_threshold_value) && is.finite(real_threshold) && real_threshold < min_otsu_threshold_value) {
      real_threshold <- min_otsu_threshold_value
      min_applied <- TRUE
    }
    label <- gsub("_minapplied$", "", label)
    if (min_applied) label <- paste0(label, "_minapplied")
    
    threshold_row <- data.frame(
      Label = label,
      RealThreshold = as.numeric(real_threshold),
      CORINE_CLASS = if (!is.null(corine_class)) corine_class else NA,
      ECOREGION_NAME = if (!is.null(ecoregion_name)) ecoregion_name else NA,
      COR_ECO_LABEL = if (!is.null(cor_eco_name)) cor_eco_name else NA,
      stringsAsFactors = FALSE
    )
    
    # --------- Binarize ----------
    binary_raster <- terra::ifel(r_filtered > real_threshold, 1, NA)
    
    # ---- Save raster burned mask if requested ----
    out_raster_path <- NULL
    if (isTRUE(write_raster)) {
      raster_name <- sprintf("BA_%s_%s_binary.tif", year, label)
      out_raster_path <- file.path(output_dir, raster_name)
      terra::writeRaster(
        binary_raster, out_raster_path,
        overwrite = TRUE,
        datatype = "INT1U",
        NAflag = 0,
        gdal = c("COMPRESS=LZW")
      )
    }
    
    # ---- Optional vectorization ----
    if (!isTRUE(vectorize)) {
      return(list(
        raster = binary_raster,
        polys = NULL,
        log = threshold_row,
        tile_dir = tile_dir,
        raster_path = out_raster_path
      ))
    }
    
    # --------- Vectorize ----------
    if (tile) {
      
      ext_r <- terra::ext(binary_raster)
      tile_width  <- (ext_r[2] - ext_r[1]) / n_cols
      tile_height <- (ext_r[4] - ext_r[3]) / n_rows
      
      count <- 1
      tile_vecs <- character(0)
      
      for (i in 0:(n_rows - 1)) {
        for (j in 0:(n_cols - 1)) {
          
          xmin_tile <- max(ext_r[1] + j * tile_width  - tile_overlap, ext_r[1])
          xmax_tile <- min(ext_r[1] + (j + 1) * tile_width + tile_overlap, ext_r[2])
          ymin_tile <- max(ext_r[3] + i * tile_height - tile_overlap, ext_r[3])
          ymax_tile <- min(ext_r[3] + (i + 1) * tile_height + tile_overlap, ext_r[4])
          
          tile_crop <- terra::crop(binary_raster, terra::ext(xmin_tile, xmax_tile, ymin_tile, ymax_tile))
          tile_path <- file.path(tile_dir, sprintf("tile_%s_%d.tif", label, count))
          
          terra::writeRaster(tile_crop, tile_path, overwrite = TRUE,
                             datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
          
          # polygonize -> SHP tile
          shp_path <- file.path(tile_dir, sprintf("tile_%s_%d.shp", label, count))
          cmd <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tile_path}" -f "ESRI Shapefile" "{shp_path}" DN')
          system(cmd)
          
          if (file.exists(tile_path)) file.remove(tile_path)
          tile_vecs <- c(tile_vecs, shp_path)
          count <- count + 1
        }
      }
      
      if (!length(tile_vecs)) {
        return(list(raster = binary_raster, polys = NULL, log = threshold_row, tile_dir = tile_dir, raster_path = out_raster_path))
      }
      
      polys <- do.call(rbind, lapply(tile_vecs, sf::st_read, quiet = TRUE))
      polys <- polys[!duplicated(sf::st_as_binary(sf::st_geometry(polys))), ]
      polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON", "MULTIPOLYGON"), ]
      polys <- sf::st_make_valid(polys)
      polys <- sf::st_transform(polys, 3035)
      
      # clean tile files
      for (shp in tile_vecs) {
        base <- tools::file_path_sans_ext(shp)
        for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
          f <- paste0(base, ext_i)
          if (file.exists(f)) file.remove(f)
        }
      }
      
    } else {
      
      tmp_raster <- file.path(tile_dir, sprintf("binary_%s.tif", label))
      terra::writeRaster(binary_raster, tmp_raster, overwrite = TRUE,
                         datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
      
      filename_base <- sprintf("BA_%s_%s", year, label)
      
      if (output_format == "geojson") {
        out_vec <- file.path(output_dir, paste0(filename_base, ".geojson"))
        if (file.exists(out_vec)) file.remove(out_vec)
        system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tmp_raster}" -f "GeoJSON" "{out_vec}" DN'))
      } else {
        out_vec <- file.path(output_dir, paste0(filename_base, ".shp"))
        base <- tools::file_path_sans_ext(out_vec)
        for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
          f <- paste0(base, ext_i)
          if (file.exists(f)) file.remove(f)
        }
        system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tmp_raster}" -f "ESRI Shapefile" "{out_vec}" DN'))
      }
      
      if (file.exists(tmp_raster)) file.remove(tmp_raster)
      polys <- sf::st_read(out_vec, quiet = TRUE)
      polys <- sf::st_make_valid(polys)
      polys <- sf::st_transform(polys, 3035)
    }
    
    # Keep only DN==1
    if (!is.null(polys) && "DN" %in% names(polys)) {
      polys <- polys[polys$DN == 1, , drop = FALSE]
    }
    
    if (is.null(polys) || !inherits(polys, "sf") || nrow(polys) == 0 || all(sf::st_is_empty(polys))) {
      return(list(raster = binary_raster, polys = NULL, log = threshold_row, tile_dir = tile_dir, raster_path = out_raster_path))
    }
    
    # Attributes
    if (!is.null(corine_class)) polys$CORINE_CLASS <- corine_class
    if (!is.null(ecoregion_name) && !is.null(ecoregion_field)) polys[[ecoregion_field]] <- ecoregion_name
    if (!is.null(corine_year)) polys$CORINE_YEAR <- corine_year
    if (!is.null(cor_eco_name) && !is.null(cor_eco_field)) polys[[cor_eco_field]] <- cor_eco_name
    polys$Label <- threshold_row$Label
    
    # Figures (hist + interclass)
    interclass_variance_curve <- sapply(1:(length(hist_values$counts) - 1), function(t) {
      w_b <- sum(hist_values$counts[1:t]) / sum(hist_values$counts)
      w_f <- 1 - w_b
      if (w_b == 0 || w_f == 0) return(NA_real_)
      mu_b <- sum(hist_values$mids[1:t] * hist_values$counts[1:t]) / sum(hist_values$counts[1:t])
      mu_f <- sum(hist_values$mids[(t + 1):length(hist_values$mids)] * hist_values$counts[(t + 1):length(hist_values$mids)]) /
        sum(hist_values$counts[(t + 1):length(hist_values$mids)])
      w_b * w_f * (mu_b - mu_f)^2
    })
    
    plot_path <- file.path(output_dir, "FIGURES", paste0("otsu_plot_", year, "_", label, ".png"))
    grDevices::png(plot_path, width = 1600, height = 800, res = 150)
    oldpar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    graphics::par(mar = c(5, 4, 4, 5) + 0.1)
    graphics::plot(hist_values$mids, smoothed_counts, type = "h", lwd = 2,
                   xlab = "Intensity", ylab = "Frequency", main = "")
    graphics::par(new = TRUE)
    graphics::plot(hist_values$mids[-1], interclass_variance_curve, type = "l", lwd = 2,
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
    graphics::axis(side = 4, lwd = 1)
    graphics::mtext("Inter-Class Variance", side = 4, line = 3)
    graphics::abline(v = threshold_value_smoothed, lty = 2)
    graphics::legend("topright",
                     legend = c("Histogram", "Inter-Class", paste("Real Threshold:", round(real_threshold, 2))),
                     lty = c(1, 1, 2), bty = "n")
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
    graphics::plot(hist_values$mids[-1], interclass_variance_curve, type = "l",
                   xlab = "Threshold", ylab = "Variance")
    grDevices::dev.off()
    
    return(list(
      raster = binary_raster,
      polys = polys,
      log = threshold_row,
      tile_dir = tile_dir,
      raster_path = out_raster_path
    ))
  }
  
  # ============================ CORINE PREP ============================
  corine_vect <- NULL
  corine_year <- NA_character_
  corine_masked <- NULL
  corine_raster_reclassed <- NULL
  
  if (!is.null(corine_raster_path)) {
    
    if (!file.exists(corine_raster_path)) stop("CORINE raster not found: ", corine_raster_path)
    corine_raster <- terra::rast(corine_raster_path)
    
    # Extract CORINE year from filename
    corine_filename <- basename(corine_raster_path)
    corine_year <- stringr::str_extract(corine_filename, "\\d{4}")
    if (is.na(corine_year)) {
      corine_year_short <- stringr::str_extract(corine_filename, "corine(\\d{2})")
      if (!is.na(corine_year_short)) {
        year_digits <- stringr::str_extract(corine_year_short, "\\d{2}")
        corine_year <- paste0("19", year_digits)
      }
    }
    
    if (reclassify_corine) {
      if (is.null(peninsula_shapefile)) stop("You set 'reclassify_corine = TRUE' but did not provide 'peninsula_shapefile'.")
      if (is.null(output_corine_raster_dir)) stop("You must specify 'output_corine_raster_dir' to save the reclassified CORINE raster.")
      if (is.null(reclass_matrix)) stop("Provide 'reclass_matrix' for reclassification.")
      
      peninsula <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula_proj <- sf::st_transform(peninsula, crs = sf::st_crs(corine_raster)) |> sf::st_make_valid()
      
      cropped <- terra::crop(corine_raster, peninsula_proj)
      masked  <- terra::mask(cropped, peninsula_proj)
      
      reclass_matrix <- as.matrix(reclass_matrix)
      reclassed_r <- terra::classify(masked, rcl = reclass_matrix, include.lowest = TRUE, right = NA)
      
      valid_classes <- unique(reclass_matrix[, 2])
      vals <- reclassed_r[]; vals[!vals %in% valid_classes] <- NA; reclassed_r[] <- vals
      
      if (!dir.exists(output_corine_raster_dir)) dir.create(output_corine_raster_dir, recursive = TRUE)
      input_name <- tools::file_path_sans_ext(basename(corine_raster_path))
      out_filename <- paste0("reclassified_", input_name, "_", format(Sys.Date(), "%Yj"), ".tif")
      out_reproj <- file.path(output_corine_raster_dir, out_filename)
      
      if (normalizePath(out_reproj, mustWork = FALSE) == normalizePath(corine_raster_path, mustWork = FALSE)) {
        stop("ERROR: Output raster path would overwrite the original CORINE raster.")
      }
      
      if (reproject) {
        tmp_unproj <- tempfile(fileext = ".tif")
        terra::writeRaster(reclassed_r, tmp_unproj, overwrite = TRUE)
        if (file.exists(out_reproj)) file.remove(out_reproj)
        cmd <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:3035" -tr {resolution} {resolution} -r near -co "COMPRESS=LZW" "{tmp_unproj}" "{out_reproj}"')
        system(cmd)
        if (!file.exists(out_reproj)) stop("ERROR: gdalwarp did not create output: ", out_reproj)
      } else {
        terra::writeRaster(reclassed_r, out_reproj, overwrite = TRUE,
                           datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))
        if (!file.exists(out_reproj)) stop("ERROR: writeRaster failed: ", out_reproj)
      }
      
      corine_raster_reclassed <- terra::rast(out_reproj)
      
      # Optional: vectorize CORINE (for intersection workflows)
      if ((vectorize || segment_by_intersection) && !is.null(output_corine_vector_dir)) {
        if (!dir.exists(output_corine_vector_dir)) dir.create(output_corine_vector_dir, recursive = TRUE)
        out_shp <- file.path(output_corine_vector_dir, paste0(tools::file_path_sans_ext(out_filename), ".shp"))
        
        # clean previous
        try({
          existing_files <- list.files(output_corine_vector_dir,
                                       pattern = paste0("^", tools::file_path_sans_ext(out_filename), "\\.(shp|shx|dbf|prj|cpg)$"),
                                       full.names = TRUE)
          if (length(existing_files) > 0) file.remove(existing_files)
        }, silent = TRUE)
        
        cmd_vec <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{out_reproj}" -f "ESRI Shapefile" "{out_shp}" DN')
        system(cmd_vec)
        if (!file.exists(out_shp)) stop("ERROR: CORINE vector shapefile not created: ", out_shp)
        
        corine_vect <- sf::st_read(out_shp, quiet = TRUE)
        if ("DN" %in% names(corine_vect)) names(corine_vect)[names(corine_vect) == "DN"] <- "CORINE_CLASS"
        corine_vect <- sf::st_make_valid(corine_vect)
        corine_vect <- sf::st_transform(corine_vect, 3035)
        sf::st_write(corine_vect, out_shp, delete_layer = TRUE, quiet = TRUE)
      }
      
    } else {
      corine_raster_reclassed <- corine_raster
    }
    
    # Align CORINE to r
    corine_resampled <- terra::resample(corine_raster_reclassed, r, method = "near")
    corine_masked <- terra::mask(corine_resampled, r)
    
    if (!is.null(corine_classes)) {
      corine_masked[!corine_masked[] %in% corine_classes] <- NA
      if (all(is.na(corine_masked[]))) stop("All CORINE values removed after filtering.")
      corine_classes <- sort(unique(terra::values(corine_masked)[!is.na(terra::values(corine_masked))]))
    } else {
      vals <- terra::values(corine_masked)
      corine_classes <- sort(unique(vals[!is.na(vals)]))
    }
  }
  
  # ============================ ECOREGIONS PREP ============================
  ecoregions_sf <- NULL
  
  if (!is.null(ecoregion_shapefile_path)) {
    if (!file.exists(ecoregion_shapefile_path)) stop("Ecoregion shapefile not found: ", ecoregion_shapefile_path)
    ecoregions_sf <- sf::st_read(ecoregion_shapefile_path, quiet = TRUE)
    
    if (!is.null(peninsula_shapefile)) {
      peninsula_sf <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula_sf <- sf::st_transform(peninsula_sf, sf::st_crs(ecoregions_sf))
      ecoregions_sf <- sf::st_intersection(ecoregions_sf, peninsula_sf)
    }
    
    if (nrow(ecoregions_sf) == 0) stop("Ecoregions shapefile has no geometries after crop.")
    ecoregions_sf <- sf::st_make_valid(ecoregions_sf)
    ecoregions_sf <- sf::st_transform(ecoregions_sf, 3035)
    
    if (is.null(ecoregion_field)) {
      if ("EnS_name" %in% names(ecoregions_sf)) {
        ecoregion_field <- "EnS_name"
      } else {
        stop("You must specify 'ecoregion_field'. Available fields: ", paste(names(ecoregions_sf), collapse = ", "))
      }
    }
    if (!(ecoregion_field %in% names(ecoregions_sf))) {
      stop(sprintf("Field '%s' not found in ecoregions. Available: %s",
                   ecoregion_field, paste(names(ecoregions_sf), collapse = ", ")))
    }
    
    if (is.null(ecoregion_classes)) {
      ecoregion_classes <- unique(ecoregions_sf[[ecoregion_field]])
    }
  }
  
  # ===================== INTERSECTION PREP (CORINE × ECO) =====================
  units_grouped <- NULL
  if (isTRUE(segment_by_intersection)) {
    
    if (is.null(corine_raster_path) || is.null(ecoregion_shapefile_path)) {
      stop("segment_by_intersection=TRUE requires both 'corine_raster_path' and 'ecoregion_shapefile_path'.")
    }
    if (is.null(corine_masked) || is.null(ecoregions_sf)) {
      stop("Intersection prep failed: missing corine_masked/ecoregions_sf.")
    }
    
    # Ensure CORINE vector exists
    if (is.null(corine_vect)) {
      tmp_corine <- tempfile(fileext = ".tif")
      terra::writeRaster(corine_masked, tmp_corine, overwrite = TRUE,
                         datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))
      tmp_shp <- tempfile(fileext = ".shp")
      system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{tmp_corine}" -f "ESRI Shapefile" "{tmp_shp}" DN'))
      corine_vect <- sf::st_read(tmp_shp, quiet = TRUE)
      if ("DN" %in% names(corine_vect)) names(corine_vect)[names(corine_vect) == "DN"] <- "CORINE_CLASS"
      corine_vect <- sf::st_make_valid(corine_vect)
      corine_vect <- sf::st_transform(corine_vect, 3035)
      corine_vect <- corine_vect[!is.na(corine_vect$CORINE_CLASS), ]
    }
    
    corine_vect$CORINE_CLASS <- as.character(corine_vect$CORINE_CLASS)
    corine_vect$CORINE_YEAR <- corine_year
    corine_vect$ID <- seq_len(nrow(corine_vect))
    
    # Rasterize ecoregions IDs on CORINE grid
    ecoregions_sf$ECO_ID_INTERNAL <- as.numeric(as.factor(ecoregions_sf[[ecoregion_field]]))
    ecoregions_r_full <- terra::rasterize(terra::vect(ecoregions_sf),
                                          terra::rast(corine_masked),
                                          field = "ECO_ID_INTERNAL")
    ecoregions_r <- terra::mask(terra::crop(ecoregions_r_full, corine_masked), corine_masked)
    
    # Assign dominant ecoregion to each CORINE polygon using centroids
    corine_centroids <- sf::st_centroid(corine_vect)
    eco_vals <- terra::extract(ecoregions_r, terra::vect(corine_centroids), ID = FALSE)
    eco_vals$ID <- corine_vect$ID
    eco_vals$ECO_ID_INTERNAL <- eco_vals[[1]]
    
    eco_id_to_name <- levels(as.factor(ecoregions_sf[[ecoregion_field]]))
    eco_vals$ECO_CLASS <- eco_id_to_name[eco_vals$ECO_ID_INTERNAL]
    eco_vals <- eco_vals[!is.na(eco_vals$ECO_CLASS), ]
    dominant_eco <- eco_vals[, c("ID", "ECO_CLASS")]
    
    corine_vect <- dplyr::left_join(corine_vect, dominant_eco, by = "ID")
    
    # Drop those still missing ECO_CLASS
    corine_vect <- corine_vect[!is.na(corine_vect$ECO_CLASS), ]
    
    # Combined label
    corine_vect$unit_id2 <- paste0(
      corine_vect$CORINE_CLASS, "_",
      gsub("[^A-Za-z0-9]", "", as.character(corine_vect$ECO_CLASS))
    )
    
    # Dissolve by unit_id2 (this is for defining units; NOT for final burned polygons writing)
    intersected_units <- sf::st_cast(corine_vect, "POLYGON")
    sf_agg <- intersected_units |>
      dplyr::group_by(unit_id2) |>
      dplyr::summarise(
        CORINE_CLASS = dplyr::first(CORINE_CLASS),
        CORINE_YEAR  = dplyr::first(CORINE_YEAR),
        ECO_CLASS    = dplyr::first(ECO_CLASS),
        .groups = "drop"
      )
    
    sf_agg <- sf::st_make_valid(sf_agg)
    sf_agg <- sf_agg[sf::st_geometry_type(sf_agg) %in% c("POLYGON", "MULTIPOLYGON"), ]
    sf_agg <- sf::st_cast(sf_agg, "MULTIPOLYGON")
    
    units_grouped <- sf::st_as_sf(sf_agg)
  }
  
  # ============================ MAIN MODES ============================
  # ---------------- Percentile trimming ----------------
  if (!is.null(trim_percentiles)) {
    
    if (!all(c("min", "max") %in% names(trim_percentiles))) {
      stop("The data frame 'trim_percentiles' must have columns named 'min' and 'max'.")
    }
    if (any(trim_percentiles$min >= trim_percentiles$max)) {
      stop("Each 'min' in 'trim_percentiles' must be smaller than its corresponding 'max'.")
    }
    
    # ---- CASE 1: CORINE × ECOREGIONS ----
    if (isTRUE(segment_by_intersection) && !is.null(units_grouped)) {
      
      for (i in seq_len(nrow(trim_percentiles))) {
        pmin <- trim_percentiles$min[i]
        pmax <- trim_percentiles$max[i]
        
        threshold_log <- NULL
        all_polys <- list()
        
        for (j in seq_len(nrow(units_grouped))) {
          unit <- units_grouped[j, ]
          u <- terra::vect(unit)
          
          mask_r <- terra::crop(r, u)
          mask_r <- terra::mask(mask_r, u)
          mask_r <- terra::trim(mask_r)
          
          label_suffix <- paste0(unit$unit_id2, "_P", pmin * 100, "_toP", pmax * 100)
          
          res_output <- process_single_threshold(
            r_input = mask_r,
            pmin = pmin, pmax = pmax,
            label_suffix = label_suffix,
            corine_class = unit$CORINE_CLASS,
            corine_year = unit$CORINE_YEAR,
            ecoregion_name = unit$ECO_CLASS,
            ecoregion_field = "ECO_CLASS",
            cor_eco_name = unit$unit_id2,
            cor_eco_field = "unit_id2",
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = FALSE,     # NO por unidad
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            res_output$polys$CORINE_CLASS <- unit$CORINE_CLASS
            res_output$polys$ECO_CLASS <- unit$ECO_CLASS
            res_output$polys$unit_id2 <- unit$unit_id2
            all_polys[[length(all_polys) + 1]] <- res_output$polys
          }
        }
        
        if (isTRUE(vectorize)) {
          if (!length(all_polys)) next
          combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
          filename_base <- sprintf("BA_%s_PERC_CORI_ECOREG_P%02d_toP%02d", year, round(pmin * 100), round(pmax * 100))
          out_file <- file.path(output_dir, paste0(filename_base, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_CORI_ECOREG_PERC_P%02d_toP%02d.txt", year, round(pmin * 100), round(pmax * 100)))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
    } else if (!is.null(corine_masked)) {
      # ---- CASE 2: Percentiles by CORINE class ----
      for (i in seq_len(nrow(trim_percentiles))) {
        pmin <- trim_percentiles$min[i]
        pmax <- trim_percentiles$max[i]
        
        threshold_log <- NULL
        all_polys <- list()
        
        for (cls in corine_classes) {
          mask_r <- terra::ifel(corine_masked == cls, r, NA)
          
          label_suffix <- paste0("corine_", gsub("[^a-zA-Z0-9]", "_", as.character(cls)), "_P", pmin * 100, "_toP", pmax * 100)
          
          res_output <- process_single_threshold(
            r_input = mask_r,
            pmin = pmin, pmax = pmax,
            label_suffix = label_suffix,
            corine_class = cls,
            corine_year = corine_year,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = write_raster,
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            res_output$polys$CORINE_CLASS <- cls
            all_polys[[length(all_polys) + 1]] <- res_output$polys
          }
        }
        
        if (isTRUE(vectorize) && length(all_polys) > 0) {
          combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
          filename <- sprintf("BA_%s_PERC_CORI_P%02d_toP%02d", year, round(pmin * 100), round(pmax * 100))
          out_file <- file.path(output_dir, paste0(filename, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_CORI_PERC_P%02d_toP%02d.txt", year, round(pmin * 100), round(pmax * 100)))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
    } else if (!is.null(ecoregions_sf)) {
      # ---- CASE 3: Percentiles by ECOREGION class ----
      for (i in seq_len(nrow(trim_percentiles))) {
        pmin <- trim_percentiles$min[i]
        pmax <- trim_percentiles$max[i]
        
        threshold_log <- NULL
        all_polys <- list()
        
        for (eco in ecoregion_classes) {
          eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
          if (nrow(eco_mask) == 0) next
          
          eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
          label_suffix <- paste0("ecoregion_", gsub("[^a-zA-Z0-9]", "_", eco), "_P", pmin * 100, "_toP", pmax * 100)
          
          res_output <- process_single_threshold(
            r_input = eco_raster_mask,
            pmin = pmin, pmax = pmax,
            label_suffix = label_suffix,
            ecoregion_name = eco,
            ecoregion_field = ecoregion_field,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = write_raster,
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            res_output$polys[[ecoregion_field]] <- eco
            all_polys[[length(all_polys) + 1]] <- res_output$polys
          }
        }
        
        if (isTRUE(vectorize) && length(all_polys) > 0) {
          combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
          filename <- sprintf("BA_%s_PERC_ECOREG_P%02d_toP%02d", year, round(pmin * 100), round(pmax * 100))
          out_file <- file.path(output_dir, paste0(filename, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_ECOREG_PERC_P%02d_toP%02d.txt", year, round(pmin * 100), round(pmax * 100)))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
    } else {
      # ---- CASE 4: Global percentiles ----
      for (i in seq_len(nrow(trim_percentiles))) {
        pmin <- trim_percentiles$min[i]
        pmax <- trim_percentiles$max[i]
        
        res_output <- process_single_threshold(
          r_input = r,
          pmin = pmin, pmax = pmax,
          min_otsu_threshold_value = min_otsu_threshold_value,
          output_dir = output_dir,
          year = year,
          tile = tile,
          n_rows = n_rows, n_cols = n_cols,
          tile_overlap = tile_overlap,
          python_exe = python_exe,
          gdal_polygonize_script = gdal_polygonize_script,
          burnable_mask = burnable_mask,
          write_raster = write_raster,
          vectorize = vectorize,
          output_format = output_format
        )
        
        if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
          combined <- sf::st_make_valid(res_output$polys)
          filename <- sprintf("BA_%s_PERC_ORIG_P%02d_toP%02d", year, round(pmin * 100), round(pmax * 100))
          out_file <- file.path(output_dir, paste0(filename, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(res_output$log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_PERC_P%02d_toP%02d.txt", year, round(pmin * 100), round(pmax * 100)))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(res_output$log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    }
    
  } else {
    # ---------------- OTSU / ORIGINAL ----------------
    
    # ============ CASE: CORINE × ECOREGION ============
    if (isTRUE(segment_by_intersection) && !is.null(units_grouped)) {
      
      # ---- ORIGINAL scenario ----
      if (do_orig) {
        threshold_log <- NULL
        all_polys <- list()
        mosaic_r <- NULL  # (FIX) mosaic accumulator for raster mode
        
        for (j in seq_len(nrow(units_grouped))) {
          unit <- units_grouped[j, ]
          u <- terra::vect(unit)
          
          mask_r <- terra::crop(r, u)
          mask_r <- terra::mask(mask_r, u)
          mask_r <- terra::trim(mask_r)
          
          res_output <- process_single_threshold(
            r_input = mask_r,
            otsu_min = NULL,
            label_suffix = paste0(unit$unit_id2, "_original"),
            corine_class = unit$CORINE_CLASS,
            corine_year = unit$CORINE_YEAR,
            ecoregion_name = unit$ECO_CLASS,
            ecoregion_field = "ECO_CLASS",
            cor_eco_name = unit$unit_id2,
            cor_eco_field = "unit_id2",
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = FALSE,    # NO por unidad (luego mosaico)
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          
          if (isTRUE(vectorize)) {
            if (!is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
              all_polys[[length(all_polys) + 1]] <- res_output$polys
            }
          } else {
            if (!is.null(res_output$raster)) {
              mosaic_r <- if (is.null(mosaic_r)) res_output$raster else terra::mosaic(mosaic_r, res_output$raster, fun = "max")
            }
          }
        } # (FIX) cierre correcto del for
        
        # --- OUTPUT (ORIGINAL CORI×ECO) ---
        if (isTRUE(vectorize)) {
          if (length(all_polys) > 0) {
            combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
            out_file <- file.path(output_dir, sprintf("BA_%s_ORIG_CORI_ECOREG%s", year, if (output_format == "geojson") ".geojson" else ".shp"))
            .write_vector(combined, out_file, format = output_format)
          }
        } else {
          if (isTRUE(write_raster) && !is.null(mosaic_r)) {
            out_raster <- file.path(output_dir, sprintf("BA_%s_ORIG_CORI_ECOREG_binary.tif", year))
            terra::writeRaster(mosaic_r, out_raster, overwrite = TRUE, datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
          }
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_ORIG_CORI_ECOREG_log.txt", year))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
      # ---- OTSU scenarios geX ----
      if (do_otsu) {
        for (otsu_min in otsu_thresholds) {
          
          threshold_log <- NULL
          all_polys <- list()
          mosaic_r <- NULL  # (FIX) reset mosaic per scenario
          
          for (j in seq_len(nrow(units_grouped))) {
            unit <- units_grouped[j, ]
            u <- terra::vect(unit)
            
            mask_r <- terra::crop(r, u)
            mask_r <- terra::mask(mask_r, u)
            mask_r <- terra::trim(mask_r)
            
            res_output <- process_single_threshold(
              r_input = mask_r,
              otsu_min = otsu_min,
              label_suffix = paste0(unit$unit_id2, "_ge", otsu_min),
              corine_class = unit$CORINE_CLASS,
              corine_year = unit$CORINE_YEAR,
              ecoregion_name = unit$ECO_CLASS,
              ecoregion_field = "ECO_CLASS",
              cor_eco_name = unit$unit_id2,
              cor_eco_field = "unit_id2",
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_dir = output_dir,
              year = year,
              tile = tile,
              n_rows = n_rows, n_cols = n_cols,
              tile_overlap = tile_overlap,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script,
              burnable_mask = burnable_mask,
              write_raster = FALSE,   # NO por unidad
              vectorize = vectorize,
              output_format = output_format
            )
            
            if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
            
            if (isTRUE(vectorize)) {
              if (!is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
                all_polys[[length(all_polys) + 1]] <- res_output$polys
              }
            } else {
              if (!is.null(res_output$raster)) {
                mosaic_r <- if (is.null(mosaic_r)) res_output$raster else terra::mosaic(mosaic_r, res_output$raster, fun = "max")
              }
            }
          } # cierre for units
          
          # --- OUTPUT per geX ---
          if (isTRUE(vectorize)) {
            if (length(all_polys) > 0) {
              combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
              out_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ECOREG_ge%d%s", year, otsu_min, if (output_format == "geojson") ".geojson" else ".shp"))
              .write_vector(combined, out_file, format = output_format)
            }
          } else {
            if (isTRUE(write_raster) && !is.null(mosaic_r)) {
              out_raster <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ECOREG_ge%d_binary.tif", year, otsu_min))
              terra::writeRaster(mosaic_r, out_raster, overwrite = TRUE, datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
            }
          }
          
          if (!is.null(threshold_log)) {
            log_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ECOREG_ge%d_log.txt", year, otsu_min))
            if (file.exists(log_file)) file.remove(log_file)
            write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
          }
        }
      }
      
    } else if (!is.null(corine_masked)) {
      # ============ CASE: CORINE ============
      if (do_orig) {
        threshold_log <- NULL
        all_polys <- list()
        
        for (cls in corine_classes) {
          mask_r <- terra::ifel(corine_masked == cls, r, NA)
          res_output <- process_single_threshold(
            r_input = mask_r,
            otsu_min = NULL,
            label_suffix = paste0("corine_", gsub("[^a-zA-Z0-9]", "_", as.character(cls)), "_original"),
            corine_class = cls,
            corine_year = corine_year,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = write_raster,
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            res_output$polys$CORINE_CLASS <- cls
            all_polys[[length(all_polys) + 1]] <- res_output$polys
          }
        }
        
        if (isTRUE(vectorize) && length(all_polys) > 0) {
          combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
          out_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ORIG%s", year, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ORIG_log.txt", year))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
      if (do_otsu) {
        for (otsu_min in otsu_thresholds) {
          threshold_log <- NULL
          all_polys <- list()
          
          for (cls in corine_classes) {
            mask_r <- terra::ifel(corine_masked == cls, r, NA)
            res_output <- process_single_threshold(
              r_input = mask_r,
              otsu_min = otsu_min,
              label_suffix = paste0("corine_", gsub("[^a-zA-Z0-9]", "_", as.character(cls)), "_ge", otsu_min),
              corine_class = cls,
              corine_year = corine_year,
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_dir = output_dir,
              year = year,
              tile = tile,
              n_rows = n_rows, n_cols = n_cols,
              tile_overlap = tile_overlap,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script,
              burnable_mask = burnable_mask,
              write_raster = write_raster,
              vectorize = vectorize,
              output_format = output_format
            )
            
            if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
            if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
              res_output$polys$CORINE_CLASS <- cls
              all_polys[[length(all_polys) + 1]] <- res_output$polys
            }
          }
          
          if (isTRUE(vectorize) && length(all_polys) > 0) {
            combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
            out_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ge%d%s", year, otsu_min, if (output_format == "geojson") ".geojson" else ".shp"))
            .write_vector(combined, out_file, format = output_format)
          }
          
          if (!is.null(threshold_log)) {
            log_file <- file.path(output_dir, sprintf("BA_%s_otsu_CORI_ge%d_log.txt", year, otsu_min))
            if (file.exists(log_file)) file.remove(log_file)
            write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
          }
        }
      }
      
    } else if (!is.null(ecoregions_sf)) {
      # ============ CASE: ECOREGION ============
      if (do_orig) {
        threshold_log <- NULL
        all_polys <- list()
        
        for (eco in ecoregion_classes) {
          eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
          if (nrow(eco_mask) == 0) next
          
          eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
          
          res_output <- process_single_threshold(
            r_input = eco_raster_mask,
            otsu_min = NULL,
            label_suffix = paste0("ecoregion_", gsub("[^a-zA-Z0-9]", "_", eco), "_original"),
            ecoregion_name = eco,
            ecoregion_field = ecoregion_field,
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = write_raster,
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            res_output$polys[[ecoregion_field]] <- eco
            all_polys[[length(all_polys) + 1]] <- res_output$polys
          }
        }
        
        if (isTRUE(vectorize) && length(all_polys) > 0) {
          combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
          out_file <- file.path(output_dir, sprintf("BA_%s_otsu_ECOREG_ORIG%s", year, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(threshold_log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_otsu_ECOREG_ORIG_log.txt", year))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
      if (do_otsu) {
        for (otsu_min in otsu_thresholds) {
          threshold_log <- NULL
          all_polys <- list()
          
          for (eco in ecoregion_classes) {
            eco_mask <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] == eco, ]
            if (nrow(eco_mask) == 0) next
            
            eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
            
            res_output <- process_single_threshold(
              r_input = eco_raster_mask,
              otsu_min = otsu_min,
              label_suffix = paste0("ecoregion_", gsub("[^a-zA-Z0-9]", "_", eco), "_ge", otsu_min),
              ecoregion_name = eco,
              ecoregion_field = ecoregion_field,
              min_otsu_threshold_value = min_otsu_threshold_value,
              output_dir = output_dir,
              year = year,
              tile = tile,
              n_rows = n_rows, n_cols = n_cols,
              tile_overlap = tile_overlap,
              python_exe = python_exe,
              gdal_polygonize_script = gdal_polygonize_script,
              burnable_mask = burnable_mask,
              write_raster = write_raster,
              vectorize = vectorize,
              output_format = output_format
            )
            
            if (!is.null(res_output$log)) threshold_log <- dplyr::bind_rows(threshold_log, res_output$log)
            if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
              res_output$polys[[ecoregion_field]] <- eco
              all_polys[[length(all_polys) + 1]] <- res_output$polys
            }
          }
          
          if (isTRUE(vectorize) && length(all_polys) > 0) {
            combined <- do.call(rbind, all_polys) |> sf::st_make_valid()
            out_file <- file.path(output_dir, sprintf("BA_%s_otsu_ECOREG_ge%d%s", year, otsu_min, if (output_format == "geojson") ".geojson" else ".shp"))
            .write_vector(combined, out_file, format = output_format)
          }
          
          if (!is.null(threshold_log)) {
            log_file <- file.path(output_dir, sprintf("BA_%s_otsu_ECOREG_ge%d_log.txt", year, otsu_min))
            if (file.exists(log_file)) file.remove(log_file)
            write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
          }
        }
      }
      
    } else {
      # ============ CASE: GLOBAL ============
      if (do_orig) {
        res_output <- process_single_threshold(
          r_input = r,
          otsu_min = NULL,
          label_suffix = "original_values",
          min_otsu_threshold_value = min_otsu_threshold_value,
          output_dir = output_dir,
          year = year,
          tile = tile,
          n_rows = n_rows, n_cols = n_cols,
          tile_overlap = tile_overlap,
          python_exe = python_exe,
          gdal_polygonize_script = gdal_polygonize_script,
          burnable_mask = burnable_mask,
          write_raster = write_raster,
          vectorize = vectorize,
          output_format = output_format
        )
        
        if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
          combined <- sf::st_make_valid(res_output$polys)
          out_file <- file.path(output_dir, sprintf("BA_%s_%s%s", year, res_output$log$Label, if (output_format == "geojson") ".geojson" else ".shp"))
          .write_vector(combined, out_file, format = output_format)
        }
        
        if (!is.null(res_output$log)) {
          log_file <- file.path(output_dir, sprintf("BA_%s_%s_log.txt", year, res_output$log$Label))
          if (file.exists(log_file)) file.remove(log_file)
          write.table(res_output$log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
      
      if (do_otsu) {
        for (otsu_min in otsu_thresholds) {
          res_output <- process_single_threshold(
            r_input = r,
            otsu_min = otsu_min,
            label_suffix = paste0("ge", otsu_min),
            min_otsu_threshold_value = min_otsu_threshold_value,
            output_dir = output_dir,
            year = year,
            tile = tile,
            n_rows = n_rows, n_cols = n_cols,
            tile_overlap = tile_overlap,
            python_exe = python_exe,
            gdal_polygonize_script = gdal_polygonize_script,
            burnable_mask = burnable_mask,
            write_raster = write_raster,
            vectorize = vectorize,
            output_format = output_format
          )
          
          if (isTRUE(vectorize) && !is.null(res_output$polys) && nrow(res_output$polys) > 0 && !all(sf::st_is_empty(res_output$polys))) {
            combined <- sf::st_make_valid(res_output$polys)
            out_file <- file.path(output_dir, sprintf("BA_%s_otsu_%s%s", year, res_output$log$Label, if (output_format == "geojson") ".geojson" else ".shp"))
            .write_vector(combined, out_file, format = output_format)
          }
          
          if (!is.null(res_output$log)) {
            log_file <- file.path(output_dir, sprintf("BA_%s_otsu_%s_log.txt", year, res_output$log$Label))
            if (file.exists(log_file)) file.remove(log_file)
            write.table(res_output$log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
          }
        }
      }
    }
  }
  
  # ---- best-effort cleanup of tiles dir ----
  tile_dir <- file.path(output_dir, "temp_tiles")
  if (dir.exists(tile_dir)) {
    try(unlink(tile_dir, recursive = TRUE, force = TRUE), silent = TRUE)
  }
  
  invisible(list(ok = TRUE, year = year))
}
