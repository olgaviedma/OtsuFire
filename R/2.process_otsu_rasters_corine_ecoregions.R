#' Process Burned Area Rasters Using Otsu Thresholding or Percentile Clipping
#' @description
#' This function computes binary burned area masks from a severity index raster (RBR or dNBR)
#' using either Otsu's thresholding or percentile-based clipping. The segmentation can be applied
#' to the full raster or stratified by CORINE land cover classes (raster) or WWF ecoregions (shapefile).
#'
#' ## Key features:
#' - Supports Otsu thresholding after filtering by minimum index values (`otsu_thresholds`)
#'   or after percentile-based clipping (`trim_percentiles`)
#' - Allows applying a minimum threshold (`min_otsu_threshold_value`) if Otsu yields low values
#' - Can stratify processing by CORINE land cover (`corine_raster_path`) or ecoregions (`ecoregion_shapefile_path`)
#' - Optionally tiles the raster and polygonizes each tile using GDAL (`gdal_polygonize.py`)
#'
#' For each segmentation, the function generates:
#' - A binary burned area raster (1 = burned, NA = unburned)
#' - A polygon shapefile of burned patches
#' - A histogram + inter-class variance plot
#' - A log file with actual thresholds used
#'
#' @name process_otsu_rasters
#' @rdname process_otsu_rasters
#'
#' @param folder_path Path to a directory or a single-band RBR or dNBR raster (required if `nbr_pre_path`/`nbr_post_path` are not used).
#' @param raster_path Path to a single-band raster file (optional alternative to `folder_path`).
#' @param nbr_pre_path Path to pre-fire NBR raster. Required if computing RBR/dNBR.
#' @param nbr_post_path Path to post-fire NBR raster. Required if computing RBR/dNBR.
#' @param output_dir Directory to save output rasters, shapefiles, plots, and logs.
#' @param year Integer or character. Year label used for naming output files. If NULL, defaults to folder name.
#' @param otsu_thresholds Numeric vector of lower bounds to pre-filter index values before applying Otsu (e.g., `c(0, 50, 100)`).
#' @param trim_percentiles Optional `data.frame` with columns `min` and `max` (e.g., 0.01 to 0.99) to clip index values before applying Otsu.
#' @param use_original Logical. If `TRUE`, uses all values from the raster without filtering (except NA). Ignored if `trim_percentiles` is set.
#' @param corine_raster_path Optional path to CORINE land cover raster for stratified segmentation.
#' @param ecoregion_shapefile_path Optional path to WWF ecoregions shapefile (must include column `EnS_name`).
#' @param min_otsu_threshold_value Numeric. Minimum acceptable threshold value. If Otsu yields a lower threshold, this value is used instead. Default: 100.
#' @param python_exe Path to the Python executable (used to call GDAL).
#' @param gdal_polygonize_script Path to `gdal_polygonize.py` script.
#' @param n_rows,n_cols Number of rows and columns for tiling the raster if `tile = TRUE`.
#' @param tile_overlap Size of buffer (in map units) around each tile to ensure seamless polygonization.
#' @param tile Logical. If `TRUE`, raster is split into tiles before polygonization.
#' @param index_type Character. One of `"RBR"` or `"dNBR"` if computing the index from pre/post NBR rasters.
#'
#' @return A named list of `sf` polygons for each segmentation. Additionally, a shapefile and log are saved to `output_dir`.
#'
#' @details
#' When stratifying by CORINE or ecoregion classes, a single combined shapefile is created for each threshold,
#' and intermediate shapefiles per class are removed. The function ensures geometry validity and avoids duplicate polygons.
#'
#' If both `trim_percentiles` and `otsu_thresholds` are provided, the function raises an error-only one method can be used.
#'
#' The function internally rescales raster values to [0, 255] before applying Otsu, and plots the histogram
#' along with the inter-class variance curve and the selected threshold.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic Otsu thresholds
#' process_otsu_rasters(
#'   folder_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   otsu_thresholds = c(0, 50, 100),
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 2: Percentile clipping (global)
#' process_otsu_rasters(
#'   folder_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   trim_percentiles = data.frame(min = 0.01, max = 0.99),
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 3: Stratified by CORINE with min threshold enforcement
#' process_otsu_rasters(
#'   folder_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   corine_raster_path = "data/corine.tif",
#'   otsu_thresholds = c(0, 50, 100),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#'
#' # Example 4: Stratified by WWF ecoregions
#' process_otsu_rasters(
#'   folder_path = "data/rbr_1985.tif",
#'   output_dir = "output/1985",
#'   ecoregion_shapefile_path = "data/ecoregions_olson.shp",
#'   otsu_thresholds = c(0, 50, 100),
#'   min_otsu_threshold_value = 150,
#'   python_exe = "C:/Python/python.exe",
#'   gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
#' )
#' }
#'
#' @importFrom terra rast crop mask ext resample ifel values freq writeRaster ncell
#' @importFrom sf st_read st_write st_transform st_geometry_type st_as_binary st_make_valid
#' @importFrom stats quantile
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom dplyr first
#' @importFrom data.table :=
#' @importFrom stringr str_detect str_replace
#' @examplesIf requireNamespace("otsuSeg", quietly = TRUE)
#' @export

utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":="
))

process_otsu_rasters <- function(folder_path = NULL,
                                 raster_path = NULL,
                                 nbr_pre_path = NULL,
                                 nbr_post_path = NULL,
                                 output_dir,
                                 year = NULL,
                                 otsu_thresholds = c(0, 25, 50, 75, 100),
                                 trim_percentiles = NULL,
                                 use_original = FALSE,
                                 corine_raster_path = NULL,
                                 ecoregion_shapefile_path = NULL,
                                 min_otsu_threshold_value = NULL,
                                 python_exe,
                                 gdal_polygonize_script,
                                 n_rows = 2,
                                 n_cols = 3,
                                 tile_overlap = 1000,
                                 tile = TRUE,
                                 index_type = "RBR") {

  tile <- as.logical(tile)[1]
  if (is.na(tile)) {
    stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")
  }

  if (is.null(year)) {
    if (!is.null(folder_path)) {
      year <- basename(folder_path)
    } else if (!is.null(raster_path)) {
      year <- tools::file_path_sans_ext(basename(raster_path))
    } else {
      stop("You must provide either 'folder_path' or 'raster_path'")
    }
  }

  if (!is.null(trim_percentiles) && !is.null(otsu_thresholds)) {
    stop("No puedes usar 'trim_percentiles' y 'otsu_thresholds' al mismo tiempo. Usa solo uno de los dos.")
  }

  if (!is.null(raster_path)) {
    raster <- terra::rast(raster_path)
    r <- raster[[1]]
    year <- year
  } else if (!is.null(nbr_pre_path) && !is.null(nbr_post_path)) {
    nbr_pre <- terra::rast(nbr_pre_path)
    nbr_post <- terra::rast(nbr_post_path)

    if (tolower(index_type) == "dnbr") {
      r <- (nbr_pre - nbr_post) * 1000
    } else if (tolower(index_type) == "rbr") {
      r <- (nbr_pre - nbr_post) * 1000 / (nbr_pre + 1.001)
    } else {
      stop("`index_type` must be either 'RBR' or 'dNBR'")
    }

    year <- year
  } else {
    stop("You must provide 'raster_path' or both 'nbr_pre_path' and 'nbr_post_path'.")
  }

  folder_path <- output_dir
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

  results <- list()
  threshold_log <- data.frame(Label = character(), RealThreshold = numeric(), stringsAsFactors = FALSE)


  if (!is.null(corine_raster_path)) {
    corine_raster <- terra::rast(corine_raster_path)
    corine_resampled <- terra::resample(corine_raster, r, method = "near")
    corine_masked <- terra::mask(corine_resampled, r)
    freq_table <- terra::freq(corine_masked)
    vals <- terra::values(corine_masked)
    corine_classes <- sort(unique(vals[!is.na(vals)]))


  } else {
    corine_classes <- NULL
  }


  if (!is.null(ecoregion_shapefile_path)) {
    ecoregions_sf <- sf::st_read(ecoregion_shapefile_path, quiet = TRUE)
    ecoregions_sf <- sf::st_transform(ecoregions_sf, sf::st_crs(r))
    ecoregion_classes <- unique(ecoregions_sf$EnS_name)
  } else {
    ecoregion_classes <- NULL
  }



  process_single_threshold <- function(r_input,
                                       pmin = NULL,
                                       pmax = NULL,
                                       otsu_min = NULL,
                                       label_suffix = "",
                                       corine_class = NULL,
                                       ecoregion_name = NULL,
                                       min_otsu_threshold_value = NULL) {
    # Subconjunto por percentiles
    if (!is.null(pmin) && !is.null(pmax)) {
      vals <- terra::values(r_input)
      vals <- vals[!is.na(vals)]
      min_val <- quantile(vals, probs = pmin, na.rm = TRUE)
      max_val <- quantile(vals, probs = pmax, na.rm = TRUE)
      vals <- vals[vals >= min_val & vals <= max_val]
      min_val <- min(vals)
      max_val <- max(vals)
      r_filtered <- terra::ifel(r_input >= min_val & r_input <= max_val, r_input, NA)
      label <- paste0("P", formatC(pmin * 1000, width = 3, flag = "0"),
                      "toP", formatC(pmax * 1000, width = 3, flag = "0"))
    } else if (!is.null(otsu_min)) {
      r_filtered <- terra::ifel(r_input >= otsu_min, r_input, NA)
      vals <- terra::values(r_filtered)
      vals <- vals[!is.na(vals)]
      min_val <- min(vals)
      max_val <- max(vals)
      label <- if (nzchar(label_suffix)) label_suffix else paste0("ge", otsu_min)
    } else {
      # Caso nuevo: usar todos los valores no NA
      vals <- terra::values(r_input)
      vals <- vals[!is.na(vals)]
      min_val <- min(vals)
      max_val <- max(vals)
      r_filtered <- terra::ifel(r_input >= min_val & r_input <= max_val, r_input, NA)
      label <- if (nzchar(label_suffix)) label_suffix else "original_values"
    }

    # Normalizacion para Otsu
    r_rescaled <- (r_filtered - min_val) / (max_val - min_val) * 255
    hist_values <- graphics::hist(terra::values(r_rescaled), breaks = 256, plot = FALSE)

    # Verificar si el paquete otsuSeg está instalado
    if (!requireNamespace("otsuSeg", quietly = TRUE)) {
      stop("Package 'otsuSeg' is required for Otsu thresholding. Please install it from GitHub:\n  remotes::install_github('olgaviedma/otsuSeg')", call. = FALSE)
    }

    smoothed_counts <- otsuSeg::smooth_histogram(hist_values$counts)
    smoothed_counts[is.na(smoothed_counts)] <- 0
    threshold_value_smoothed <- otsuSeg::otsu_threshold_smoothed(smoothed_counts, hist_values$mids)
    real_threshold <- (threshold_value_smoothed / 255) * (max_val - min_val) + min_val

    message(sprintf("Threshold (%s): real threshold value = %.4f", label, real_threshold))

    # Reglas de umbral minimo
    #if (!is.null(corine_class) && corine_class %in% c(22, 23, 26) && real_threshold < 300) { # 22: agroforestry; 23: broad-leaved; 26: grasslands
     # message(sprintf(" ? CORINE class %s: Otsu threshold %.2f sustituido por minimo 300.", corine_class, real_threshold))
     # real_threshold <- 300
   # } else

    # === Aplicar umbral minimo si corresponde ===

    min_applied <- FALSE  # Por defecto, no se aplica

    if (!is.null(min_otsu_threshold_value)) {
      if (real_threshold < min_otsu_threshold_value) {
        message(sprintf(" ? Threshold %.2f menor que el minimo %.2f. Se aplica el minimo.",
                        real_threshold, min_otsu_threshold_value))
        real_threshold <- min_otsu_threshold_value
        min_applied <- TRUE  # Se aplico el umbral minimo
      }
    }

    label <- gsub("_minapplied$", "", label)  # limpiar por si ya lo tiene
    if (min_applied) {
      label <- paste0(label, "_minapplied")
    }


    threshold_row <- data.frame(Label = label, RealThreshold = real_threshold)

    # Binarizacion
    binary_raster <- terra::ifel(r_filtered > real_threshold, 1, NA)
    raster_name <- file.path(output_dir, paste0("raster_", year, "_otsu_", label, ".tif"))
    terra::writeRaster(binary_raster, raster_name, overwrite = TRUE,
                       datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=-9999"))

    if (file.exists(raster_name)) file.remove(raster_name)

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
          tile_path <- file.path(folder_path, sprintf("tile_%s_%d.tif", label, count))
          terra::writeRaster(tile_crop, tile_path, overwrite = TRUE,
                             datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
          shp_path <- file.path(folder_path, sprintf("tile_%s_%d.shp", label, count))
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
      temp_raster <- file.path(folder_path, sprintf("binary_%s.tif", label))
      terra::writeRaster(binary_raster, temp_raster, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=0"))
      shp_path <- file.path(folder_path, sprintf("burned_areas_%s_otsu_%s.shp", year, label))
      system(glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{temp_raster}" -f "ESRI Shapefile" "{shp_path}" DN'))
      if (file.exists(temp_raster)) file.remove(temp_raster)
      polys <- sf::st_read(shp_path, quiet = TRUE)
      polys <- sf::st_transform(polys, 3035)
    }

    if (!is.null(corine_class)) {
      polys$CORINE_CLASS <- corine_class
    }
    if (!is.null(ecoregion_name)) {
      polys$EnS_name <- ecoregion_name
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

    figures_dir <- file.path(folder_path, "FIGURES")
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

    return(list(polys = polys, log = threshold_row))
  }

  if (!is.null(trim_percentiles)) {
    # Validaciones de estructura
    if (!all(c("min", "max") %in% names(trim_percentiles))) {
      stop("El data frame 'trim_percentiles' debe tener columnas llamadas 'min' y 'max'.")
    }
    if (any(trim_percentiles$min >= trim_percentiles$max)) {
      stop("Cada valor 'min' en 'trim_percentiles' debe ser menor que su correspondiente 'max'.")
    }

    if (!is.null(trim_percentiles)) {
      if (!all(c("min", "max") %in% names(trim_percentiles))) {
        stop("El data frame 'trim_percentiles' debe tener columnas llamadas 'min' y 'max'.")
      }
      if (any(trim_percentiles$min >= trim_percentiles$max)) {
        stop("Cada valor 'min' en 'trim_percentiles' debe ser menor que su correspondiente 'max'.")
      }

      # === CASO 1: Percentiles por clase CORINE ===
      if (!is.null(corine_classes)) {
        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]
          all_polys <- list()

          for (cls in corine_classes) {
            mask_r <- terra::ifel(corine_masked == cls, r, NA)
            label_suffix <- paste0("corine_", cls, "_P", formatC(pmin * 100, digits = 0), "_toP", formatC(pmax * 100, digits = 0))

            res <- process_single_threshold(
              r_input = mask_r,
              pmin = pmin,
              pmax = pmax,
              label_suffix = label_suffix,
              corine_class = cls,
              min_otsu_threshold_value = min_otsu_threshold_value
            )

            res$polys$CORINE_CLASS <- cls
            res$polys$Label <- res$log$Label
            all_polys[[length(all_polys) + 1]] <- res$polys
            threshold_log <- rbind(threshold_log, res$log)
          }

          combined <- do.call(rbind, all_polys)
          combined <- sf::st_make_valid(combined)
          names(combined) <- make.unique(substr(names(combined), 1, 10))

          # Verificar si alguna clase aplico el minimo
          any_minapplied <- any(grepl("_minapplied$", sapply(all_polys, function(p) unique(p$Label))))
          filename_base <- sprintf("burned_areas_%s_percentiles_corine_P%d_toP%d", year, pmin * 100, pmax * 100)
          filename <- if (any_minapplied) paste0(filename_base, "_minapplied") else filename_base
          out_shp <- file.path(folder_path, paste0(filename, ".shp"))

          shp_base <- tools::file_path_sans_ext(out_shp)
          for (ext in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
            f <- paste0(shp_base, ext)
            if (file.exists(f)) file.remove(f)
          }
          sf::st_write(combined, out_shp, append = FALSE)
        }

        # === CASO 2: Percentiles por clase de ecoregion ===
      } else if (!is.null(ecoregion_classes)) {
        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]
          all_polys <- list()

          for (eco in ecoregion_classes) {
            eco_mask <- ecoregions_sf[ecoregions_sf$EnS_name == eco, ]
            eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
            clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)
            label_suffix <- paste0("ecoregion_", clean_eco, "_P", formatC(pmin * 100, digits = 0), "_toP", formatC(pmax * 100, digits = 0))

            res <- process_single_threshold(
              r_input = eco_raster_mask,
              pmin = pmin,
              pmax = pmax,
              label_suffix = label_suffix,
              ecoregion_name = eco,
              min_otsu_threshold_value = min_otsu_threshold_value
            )

            res$polys$EnS_name <- eco
            res$polys$Label <- res$log$Label
            all_polys[[length(all_polys) + 1]] <- res$polys
            threshold_log <- rbind(threshold_log, res$log)
          }

          combined <- do.call(rbind, all_polys)
          combined <- sf::st_make_valid(combined)
          names(combined) <- make.unique(substr(names(combined), 1, 10))

          # Verificar si alguna clase aplico el minimo
          any_minapplied <- any(grepl("_minapplied$", sapply(all_polys, function(p) unique(p$Label))))
          filename_base <- sprintf("burned_areas_%s_percentiles_corine_P%d_toP%d", year, pmin * 100, pmax * 100)
          filename <- if (any_minapplied) paste0(filename_base, "_minapplied") else filename_base
          out_shp <- file.path(folder_path, paste0(filename, ".shp"))

          shp_base <- tools::file_path_sans_ext(out_shp)
          for (ext in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
            f <- paste0(shp_base, ext)
            if (file.exists(f)) file.remove(f)
          }
          sf::st_write(combined, out_shp, append = FALSE)
        }

        # === CASO 3: Percentiles sobre raster global ===
      } else {
        for (i in seq_len(nrow(trim_percentiles))) {
          pmin <- trim_percentiles$min[i]
          pmax <- trim_percentiles$max[i]

          res <- process_single_threshold(
            r_input = r,
            pmin = pmin,
            pmax = pmax,
            min_otsu_threshold_value = min_otsu_threshold_value
          )

          label_real <- res$log$Label
          results[[label_real]] <- res$polys
          threshold_log <- rbind(threshold_log, res$log)

          # Verificar si alguna clase aplico el minimo
          any_minapplied <- any(grepl("_minapplied$", sapply(all_polys, function(p) unique(p$Label))))
          filename_base <- sprintf("burned_areas_%s_percentiles_corine_P%d_toP%d", year, pmin * 100, pmax * 100)
          filename <- if (any_minapplied) paste0(filename_base, "_minapplied") else filename_base
          out_shp <- file.path(folder_path, paste0(filename, ".shp"))

          shp_base <- tools::file_path_sans_ext(out_shp)
          for (ext in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
            f <- paste0(shp_base, ext)
            if (file.exists(f)) file.remove(f)
          }
          sf::st_write(res$polys, out_shp, append = FALSE)
        }
      }
    }


  } else if (!is.null(corine_classes)) {
    if (use_original) {
      all_polys <- list()
      for (cls in corine_classes) {
        mask_r <- terra::ifel(corine_masked == cls, r, NA)
        label_suffix <- paste0("corine_class_", cls, "_original")

        res <- process_single_threshold(
          r_input = mask_r,
          otsu_min = NULL,
          label_suffix = label_suffix,
          corine_class = cls,
          min_otsu_threshold_value = min_otsu_threshold_value
        )

        label_real <- res$log$Label
        res$polys$CORINE_CLASS <- cls
        res$polys$Label <- label_real
        all_polys[[label_real]] <- res$polys
        threshold_log <- rbind(threshold_log, res$log)
      }

      combined <- do.call(rbind, all_polys)
      combined <- sf::st_make_valid(combined)
      names(combined) <- make.unique(substr(names(combined), 1, 10))

      out_shp <- file.path(folder_path, sprintf("burned_areas_%s_otsu_corine_original.shp", year))
      shp_base <- tools::file_path_sans_ext(out_shp)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext in shp_exts) {
        f <- paste0(shp_base, ext)
        if (file.exists(f)) file.remove(f)
      }
      sf::st_write(combined, out_shp, append=FALSE)

    } else {
      for (otsu_min in otsu_thresholds) {
        all_polys <- list()
        for (cls in corine_classes) {
          mask_r <- terra::ifel(corine_masked == cls, r, NA)
          clean_cls <- gsub("[^a-zA-Z0-9]", "_", as.character(cls))
          label_suffix <- paste0("corine_", clean_cls, "_ge", otsu_min)

          res <- process_single_threshold(
            r_input = mask_r,
            otsu_min = otsu_min,
            label_suffix = label_suffix,
            corine_class = cls,
            min_otsu_threshold_value = min_otsu_threshold_value
          )

          label_real <- res$log$Label
          res$polys$CORINE_CLASS <- cls
          res$polys$Label <- label_real
          all_polys[[label_real]] <- res$polys
          threshold_log <- rbind(threshold_log, res$log)
        }

        combined <- do.call(rbind, all_polys)
        combined <- sf::st_make_valid(combined)
        names(combined) <- make.unique(substr(names(combined), 1, 10))

        out_shp <- file.path(folder_path, sprintf("burned_areas_%s_otsu_corine_ge%d.shp", year, otsu_min))
        shp_base <- tools::file_path_sans_ext(out_shp)
        shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        for (ext in shp_exts) {
          f <- paste0(shp_base, ext)
          if (file.exists(f)) file.remove(f)
        }
        sf::st_write(combined, out_shp, append=FALSE)

        log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_corine_ge%d_log.txt", year, otsu_min))
        labels_used <- unlist(lapply(all_polys, function(p) {
          if (!is.null(p$Label)) unique(p$Label) else NA
        }))
        labels_used <- labels_used[!is.na(labels_used)]

        if (!is.null(threshold_log) && "Label" %in% colnames(threshold_log)) {
          log_to_write <- threshold_log[threshold_log$Label %in% labels_used, ]
        } else {
          log_to_write <- threshold_log
        }

        write.table(log_to_write,
                    file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      }
    }

  } else if (!is.null(ecoregion_classes)) {
    if (use_original) {
      all_polys <- list()
      for (eco in ecoregion_classes) {
        eco_mask <- ecoregions_sf[ecoregions_sf$EnS_name == eco, ]
        eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
        clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)
        label_suffix <- paste0("ecoregion_", clean_eco, "_original")

        res <- process_single_threshold(
          r_input = eco_raster_mask,
          otsu_min = NULL,
          label_suffix = label_suffix,
          ecoregion_name = eco,
          min_otsu_threshold_value = min_otsu_threshold_value
        )

        label_real <- res$log$Label
        res$polys$EnS_name <- eco
        res$polys$Label <- label_real
        all_polys[[label_real]] <- res$polys
        threshold_log <- rbind(threshold_log, res$log)
      }

      combined <- do.call(rbind, all_polys)
      combined <- sf::st_make_valid(combined)
      names(combined) <- make.unique(substr(names(combined), 1, 10))

      out_shp <- file.path(folder_path, sprintf("burned_areas_%s_otsu_ecoregion_original.shp", year))
      shp_base <- tools::file_path_sans_ext(out_shp)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext in shp_exts) {
        f <- paste0(shp_base, ext)
        if (file.exists(f)) file.remove(f)
      }
      sf::st_write(combined, out_shp, append=FALSE)

    } else {
      for (otsu_min in otsu_thresholds) {
        all_polys <- list()
        for (eco in ecoregion_classes) {
          eco_mask <- ecoregions_sf[ecoregions_sf$EnS_name == eco, ]
          eco_raster_mask <- terra::mask(r, terra::vect(eco_mask))
          clean_eco <- gsub("[^a-zA-Z0-9]", "_", eco)
          label_suffix <- paste0("ecoregion_", clean_eco, "_ge", otsu_min)

          res <- process_single_threshold(
            r_input = eco_raster_mask,
            otsu_min = otsu_min,
            label_suffix = label_suffix,
            ecoregion_name = eco,
            min_otsu_threshold_value = min_otsu_threshold_value
          )

          label_real <- res$log$Label
          res$polys$EnS_name <- eco
          res$polys$Label <- label_real
          all_polys[[label_real]] <- res$polys
          threshold_log <- rbind(threshold_log, res$log)
        }

        combined <- do.call(rbind, all_polys)
        combined <- sf::st_make_valid(combined)
        names(combined) <- make.unique(substr(names(combined), 1, 10))

        out_shp <- file.path(folder_path, sprintf("burned_areas_%s_otsu_ecoregion_ge%d.shp", year, otsu_min))
        shp_base <- tools::file_path_sans_ext(out_shp)
        shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        for (ext in shp_exts) {
          f <- paste0(shp_base, ext)
          if (file.exists(f)) file.remove(f)
        }
        sf::st_write(combined, out_shp, append=FALSE)

        log_file <- file.path(output_dir, sprintf("burned_areas_%s_otsu_ecoregion_ge%d_log.txt", year, otsu_min))
        labels_used <- unlist(lapply(all_polys, function(p) {
          if (!is.null(p$Label)) unique(p$Label) else NA
        }))
        labels_used <- labels_used[!is.na(labels_used)]

        if (!is.null(threshold_log) && "Label" %in% colnames(threshold_log)) {
          log_to_write <- threshold_log[threshold_log$Label %in% labels_used, ]
        } else {
          log_to_write <- threshold_log
        }

        write.table(log_to_write,
                    file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

      }
    }

  } else if (use_original) {
    all_polys <- list()
    for (otsu_min in otsu_thresholds) {
      label_suffix <- paste0("original_ge", otsu_min)
      r_filtered <- terra::ifel(r > otsu_min, r, NA)

      if (all(is.na(terra::values(r_filtered)))) {
        message(sprintf("No hay valores mayores que %.2f en el raster. Se omite '%s'.", otsu_min, label_suffix))
        next
      }

      res <- process_single_threshold(
        r_input = r_filtered,
        otsu_min = otsu_min,
        label_suffix = label_suffix,
        min_otsu_threshold_value = min_otsu_threshold_value
      )

      label_real <- res$log$Label
      results[[label_real]] <- res$polys
      threshold_log <- rbind(threshold_log, res$log)

      out_shp <- file.path(folder_path, sprintf("burned_areas_%s_%s.shp", year, label_real))
      shp_base <- tools::file_path_sans_ext(out_shp)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext in shp_exts) {
        f <- paste0(shp_base, ext)
        if (file.exists(f)) file.remove(f)
      }
      sf::st_write(res$polys, out_shp, append = FALSE)
    }

  } else {
    for (otsu_min in otsu_thresholds) {
      label_suffix <- paste0("ge", otsu_min)

      res <- process_single_threshold(
        r_input = r,
        otsu_min = otsu_min,
        label_suffix = label_suffix,
        min_otsu_threshold_value = min_otsu_threshold_value
      )

      label_real <- res$log$Label
      results[[label_real]] <- res$polys
      threshold_log <- rbind(threshold_log, res$log)

      out_shp <- file.path(folder_path, sprintf("burned_areas_%s_otsu_%s.shp", year, label_real))
      shp_base <- tools::file_path_sans_ext(out_shp)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext in shp_exts) {
        f <- paste0(shp_base, ext)
        if (file.exists(f)) file.remove(f)
      }
      sf::st_write(res$polys, out_shp, append = FALSE)
    }
  }

  # Ruta de salida del log
  # === Crear nombre de log coherente con shapefile de salida ===
  if (exists("out_shp")) {
    log_base <- tools::file_path_sans_ext(basename(out_shp))
    log_file <- file.path(output_dir, paste0(log_base, "_log.txt"))
  } else {
    # Caso para percentiles u otros casos sin shapefile combinado
    log_file <- file.path(output_dir, sprintf("otsu_thresholds_%s_percentiles_log.txt", year))
  }

  write.table(threshold_log, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)

  return(results)
}
