#' Mask and mosaic raster tiles, with optional one-step (robust) warp to target CRS/resolution
#'
#' Applies an optional polygon mask to raster tiles in a folder, mosaics them using
#' a maximum-value composite, and (optionally) warps the result to a target CRS and
#' output resolution using GDAL (\code{gdalwarp}).
#'
#' Input tiles are selected via flexible glob patterns (e.g., \code{"MinMin_*.tif"})
#' and can be excluded using additional globs (by default, already excludes mosaics
#' and derived outputs such as \code{*_proj.tif} and \code{*_res*m.tif}).
#'
#' \strong{Output stages (\code{return_stage})}
#' \itemize{
#'   \item \code{"mosaic"}: returns the raw mosaic GeoTIFF created from the VRT.
#'   If \code{mosaic_res_m} is provided, the mosaic filename includes
#'   \code{_mosaic_res{mosaic_res_m}m}.
#'   \item \code{"mosaic_clean"}: returns a cleaned mosaic with band names restored and
#'   a suffix indicating masking: \code{_mosaic_masked.tif} or \code{_mosaic_nomask.tif}.
#'   \item \code{"resampled"} (or alias \code{"resample"}): returns the final warped output
#'   with suffix \code{_res{out_res_m}m.tif}.
#' }
#'
#' \strong{Robust (recommended) warp mode}
#' By default (\code{proj_res_m = NULL}), the function uses a \emph{single} \code{gdalwarp}
#' call that reprojects and resamples directly to \code{out_res_m}. This avoids creating
#' very large intermediate rasters for big AOIs and reduces the risk of BigTIFF write errors.
#' If \code{proj_res_m} is provided, a two-step warp is performed (\code{proj_res_m} then
#' \code{out_res_m}).
#'
#' \strong{Resampling method and value handling}
#' For multiband outputs (e.g., RBR + DOY), the warp uses nearest-neighbor resampling to
#' preserve integer/categorical bands. For single-band rasters, bilinear resampling is used.
#' The final \code{"resampled"} stage is re-written as integer (values are \code{floor()}'d)
#' and saved as \code{INT2S} with \code{NAflag = -9999}. All GDAL outputs are written with
#' \code{BIGTIFF=YES}, \code{TILED=YES}, and \code{COMPRESS=LZW}.
#'
#' \strong{Intermediates and cleanup}
#' Even when requesting \code{"resampled"}, the function necessarily creates intermediate
#' files (mosaic and mosaic_clean) before warping. These are automatically deleted depending
#' on \code{keep_intermediates}. Masked tiles are removed unless \code{keep_masked_tiles=TRUE}.
#' Note that pre-existing outputs from older runs are not automatically removed unless they
#' are overwritten in the current run.
#'
#' @param folder_path Character. Folder containing input raster tiles.
#' @param mask_shapefile Character or NULL. Path to polygon shapefile used as mask.
#' If NULL/invalid and \code{apply_mask=TRUE}, the function will run without masking.
#' @param gdalwarp_path Character. Full path to the GDAL \code{gdalwarp} executable.
#' @param crs_target Integer. Target EPSG code for the final warp (default 3035; ETRS89 / LAEA).
#' @param tile_buffer Numeric. Buffer distance (map units) used if mask and raster do not overlap
#' initially (default 1000).
#' @param year Character/integer/NULL. Optional year used in output naming. If NULL, the function
#' attempts to infer a year from \code{folder_path} or the first raster filename.
#'
#' @param tif_glob Character vector. Glob(s) used to select input tiles (default \code{"*.tif"}).
#' Examples: \code{"MinMin_*.tif"}, \code{c("MinMin_*.tif","MaxMax_*.tif")}.
#' @param exclude_glob Character vector. Glob(s) used to exclude files (default excludes mosaics,
#' projections, and resampled outputs).
#'
#' @param apply_mask Logical. If TRUE (default) and \code{mask_shapefile} is valid, tiles are masked
#' before mosaicking.
#' @param return_stage Character. One of \code{"resampled"}, \code{"mosaic_clean"}, \code{"mosaic"}.
#' The alias \code{"resample"} is accepted and treated as \code{"resampled"}.
#'
#' @param keep_intermediates Character. Controls which intermediate products are kept:
#' \code{"none"} (default), \code{"mosaic"}, \code{"mosaic_clean"}, or \code{"all"}.
#' @param keep_masked_tiles Logical. If TRUE, masked tile rasters under \code{masked_rasters/} are
#' retained for debugging (default FALSE).
#'
#' @param mosaic_res_m Numeric or NULL. Optional resolution (meters) enforced at the mosaic step.
#' If provided, the mosaic filename includes \code{_mosaic_res{mosaic_res_m}m}.
#' @param proj_res_m Numeric or NULL. If NULL (default), performs a one-step robust warp directly to
#' \code{out_res_m}. If provided, performs a two-step warp: \code{proj_res_m} then \code{out_res_m}.
#' @param out_res_m Numeric. Output resolution (meters) for the final warped product (default 90).
#'
#' @param preserve_band_names Logical. If TRUE (default), attempts to preserve band names from the first
#' input tile and re-apply them to outputs.
#' @param band_names Character vector or NULL. If provided, forces these names to be used (must match
#' the number of layers in the input rasters).
#'
#' @return Character. Normalized path to the output raster corresponding to \code{return_stage}.
#'
#' @examples
#' \dontrun{
#' # 1) Mosaic only (no mask), forcing mosaic resolution and returning mosaic
#' out_mosaic <- mask_mosaic_raster(
#'   folder_path    = "C:/ZENODO/exdata",
#'   mask_shapefile = NULL,
#'   apply_mask     = FALSE,
#'   tif_glob       = "MinMin_*.tif",
#'   exclude_glob   = c("*_mosaic*.tif", "*_proj.tif", "*_res*m.tif"),
#'   mosaic_res_m   = 90,
#'   year           = 2012,
#'   gdalwarp_path  = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe",
#'   return_stage   = "mosaic"
#' )
#'
#' # 2) Full workflow to final output (robust one-step warp: proj_res_m=NULL)
#' out_res <- mask_mosaic_raster(
#'   folder_path    = "C:/ZENODO/exdata",
#'   mask_shapefile = "C:/ZENODO/exdata/iberian_peninsula_proj_final.shp",
#'   apply_mask     = TRUE,
#'   tif_glob       = "MinMin_*.tif",
#'   exclude_glob   = c("*_mosaic*.tif", "*_proj.tif", "*_res*m.tif"),
#'   year           = 2012,
#'   crs_target     = 3035,
#'   out_res_m      = 90,
#'   proj_res_m     = NULL,
#'   band_names     = c("rbr","doy"),
#'   gdalwarp_path  = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe",
#'   return_stage   = "resampled",
#'   keep_intermediates = "none"
#' )
#' }
#'
#' @importFrom terra vect rast nlyr crs project ext buffer intersect crop mask writeRaster NAflag
#' @importFrom gdalUtilities gdalbuildvrt
#' @importFrom stringr str_extract str_match
#' @importFrom tools file_path_sans_ext
#' @importFrom utils glob2rx
#'
#' @export

mask_mosaic_raster <- function(
    folder_path,
    mask_shapefile = NULL,
    gdalwarp_path,
    crs_target = 3035,
    tile_buffer = 1000,
    year = NULL,

    # file selection
    tif_glob = "*.tif",
    exclude_glob = c("*_mosaic*.tif", "*_proj.tif", "*_res*m.tif"),

    # mask / workflow control
    apply_mask = TRUE,
    return_stage = c("resampled", "mosaic_clean", "mosaic", "resample"),

    # cleanup control
    keep_intermediates = c("none", "mosaic", "mosaic_clean", "all"),
    keep_masked_tiles = FALSE,

    # mosaic resolution control
    mosaic_res_m = NULL,

    # reprojection / output resolution control
    # If proj_res_m is NULL -> one-step warp directly to out_res_m (recommended for large AOIs)
    proj_res_m = NULL,
    out_res_m  = 90,

    # band name handling
    preserve_band_names = TRUE,
    band_names = NULL
) {

  # ---------------------------
  # args
  # ---------------------------
  return_stage <- tolower(as.character(return_stage)[1])
  if (return_stage == "resample") return_stage <- "resampled"
  if (!return_stage %in% c("resampled", "mosaic_clean", "mosaic")) {
    stop("return_stage must be one of: 'resampled', 'mosaic_clean', 'mosaic' (or 'resample').")
  }

  keep_intermediates <- match.arg(keep_intermediates)

  if (missing(gdalwarp_path) || is.null(gdalwarp_path) || !file.exists(gdalwarp_path)) {
    stop("gdalwarp_path must be defined and point to a valid executable.")
  }
  if (!dir.exists(folder_path)) stop("folder_path does not exist: ", folder_path)

  message("Processing folder: ", folder_path)

  # ---------------------------
  # select rasters by glob
  # ---------------------------
  resolve_globs <- function(globs) {
    globs <- as.character(globs)
    out <- character(0)
    for (g in globs) {
      out <- c(out, list.files(folder_path, pattern = glob2rx(g), full.names = TRUE))
    }
    unique(out)
  }

  raster_files <- resolve_globs(tif_glob)

  if (!is.null(exclude_glob) && length(exclude_glob) > 0) {
    excl <- resolve_globs(exclude_glob)
    raster_files <- setdiff(raster_files, excl)
  }

  raster_files <- raster_files[file.exists(raster_files)]
  if (length(raster_files) == 0) {
    stop("No rasters found in ", folder_path, " after applying tif_glob/exclude_glob.")
  }

  # ---------------------------
  # band names
  # ---------------------------
  infer_band_names <- function(files, check_n = 5L) {
    r1 <- terra::rast(files[1])
    n1 <- terra::nlyr(r1)

    nm1 <- names(r1)
    if (is.null(nm1) || length(nm1) != n1 || any(!nzchar(nm1))) {
      nm1 <- paste0("band", seq_len(n1))
    }

    k <- min(length(files), as.integer(check_n))
    if (k > 1) {
      for (i in 2:k) {
        ri <- terra::rast(files[i])
        if (terra::nlyr(ri) != n1) {
          warning("[mask_mosaic_raster] Inputs have different number of layers. Using generic names.")
          return(paste0("band", seq_len(n1)))
        }
        nmi <- names(ri)
        if (is.null(nmi) || length(nmi) != n1 || any(!nzchar(nmi))) next
        if (!identical(nmi, nm1)) {
          warning("[mask_mosaic_raster] Inputs have different band names. Using names from the first raster: ",
                  paste(nm1, collapse = ", "))
          break
        }
      }
    }
    nm1
  }

  band_names_use <- NULL
  if (!is.null(band_names)) {
    band_names <- as.character(band_names)
    r0 <- terra::rast(raster_files[1])
    if (length(band_names) != terra::nlyr(r0)) {
      stop("band_names length (", length(band_names),
           ") must match number of layers in inputs (", terra::nlyr(r0), ").")
    }
    band_names_use <- band_names
  } else if (isTRUE(preserve_band_names)) {
    band_names_use <- infer_band_names(raster_files)
  }

  apply_names_if_possible <- function(r) {
    if (!is.null(band_names_use) && terra::nlyr(r) == length(band_names_use)) {
      names(r) <- band_names_use
    } else if (is.null(names(r)) || any(!nzchar(names(r)))) {
      names(r) <- paste0("band", seq_len(terra::nlyr(r)))
    }
    r
  }

  # ---------------------------
  # infer year for naming
  # ---------------------------
  infer_year <- function(x) {
    m <- stringr::str_extract(as.character(x), "(19|20)[0-9]{2}")
    if (is.na(m) || !nzchar(m)) return(NULL)
    m
  }

  if (is.null(year)) {
    year <- infer_year(basename(folder_path))
    if (is.null(year)) year <- infer_year(basename(raster_files[1]))
    if (is.null(year)) year <- "unknown"
  } else {
    year <- as.character(year)[1]
    if (!nzchar(year)) year <- "unknown"
  }

  # ---------------------------
  # masking decision
  # ---------------------------
  do_mask <- isTRUE(apply_mask) &&
    !is.null(mask_shapefile) &&
    is.character(mask_shapefile) &&
    length(mask_shapefile) == 1 &&
    nzchar(mask_shapefile) &&
    file.exists(mask_shapefile)

  if (isTRUE(apply_mask) && !do_mask) {
    message("[INFO] apply_mask=TRUE but mask_shapefile is NULL/invalid -> running WITHOUT mask.")
  }

  masked_files <- raster_files
  masked_output_folder <- NULL

  if (do_mask) {
    mask_vect <- terra::vect(mask_shapefile)

    masked_output_folder <- file.path(folder_path, "masked_rasters")
    dir.create(masked_output_folder, showWarnings = FALSE, recursive = TRUE)

    example_raster <- terra::rast(raster_files[1])
    if (terra::crs(example_raster) != terra::crs(mask_vect)) {
      message("Reprojecting mask to match raster CRS...")
      mask_vect <- terra::project(mask_vect, terra::crs(example_raster))
    }

    mask_raster <- function(raster_path, mask_vect, output_path, buffer_width) {
      r <- terra::rast(raster_path)
      if (terra::crs(r) != terra::crs(mask_vect)) {
        mask_vect <- terra::project(mask_vect, terra::crs(r))
      }

      if (is.null(terra::intersect(terra::ext(r), terra::ext(mask_vect)))) {
        message("No overlap. Applying buffer...")
        mask_vect <- terra::buffer(mask_vect, width = buffer_width)
      }
      if (is.null(terra::intersect(terra::ext(r), terra::ext(mask_vect)))) {
        message("No overlap even after buffering. Skipping: ", basename(raster_path))
        return(NA_character_)
      }

      r_masked <- terra::mask(terra::crop(r, mask_vect), mask_vect)
      terra::writeRaster(r_masked, output_path, overwrite = TRUE)
      output_path
    }

    masked_files <- vapply(raster_files, function(r) {
      outp <- file.path(masked_output_folder, paste0("masked_", tools::file_path_sans_ext(basename(r)), ".tif"))
      mask_raster(r, mask_vect, outp, tile_buffer)
    }, character(1))

    masked_files <- masked_files[!is.na(masked_files)]
    if (length(masked_files) == 0) stop("No masked rasters generated - all skipped due to lack of overlap.")
  }

  # ---------------------------
  # mosaic (max)
  # ---------------------------
  mosaic_tmp <- file.path(folder_path, "_mosaic_tmp")
  dir.create(mosaic_tmp, showWarnings = FALSE, recursive = TRUE)

  vrt_path <- file.path(mosaic_tmp, paste0("temp_mosaic_", year, ".vrt"))
  gdalUtilities::gdalbuildvrt(masked_files, vrt_path, overwrite = TRUE)
  message("VRT created: ", vrt_path)

  prefix_match <- stringr::str_match(basename(raster_files[1]), "(.*?)([12][0-9]{3})")[, 2]
  prefix <- if (!is.na(prefix_match)) prefix_match else "composite_"

  mosaic_tag <- if (!is.null(mosaic_res_m)) paste0("_mosaic_res", as.numeric(mosaic_res_m)[1], "m") else "_mosaic"
  tif_output <- file.path(folder_path, paste0(prefix, year, mosaic_tag, ".tif"))
  if (file.exists(tif_output)) file.remove(tif_output)

  dst_nodata <- "-9999"
  cmd_mosaic_args <- c(
    "-overwrite",
    "-r", "max",
    "-dstnodata", dst_nodata,
    "-co", "COMPRESS=LZW",
    "-co", "BIGTIFF=YES",
    "-co", "TILED=YES"
  )
  if (!is.null(mosaic_res_m)) {
    mosaic_res_m <- as.numeric(mosaic_res_m)[1]
    if (!is.finite(mosaic_res_m) || mosaic_res_m <= 0) stop("mosaic_res_m must be a positive number.")
    cmd_mosaic_args <- c(cmd_mosaic_args, "-tr", as.character(mosaic_res_m), as.character(mosaic_res_m))
  }
  cmd_mosaic_args <- c(cmd_mosaic_args, shQuote(vrt_path), shQuote(tif_output))

  gdal_out <- system2(gdalwarp_path, args = cmd_mosaic_args, stdout = TRUE, stderr = TRUE)
  status <- attr(gdal_out, "status"); if (is.null(status)) status <- 0L
  if (status != 0L || !file.exists(tif_output)) {
    stop("Mosaic creation failed or mosaic file not created: ", tif_output,
         "\n--- GDAL output (last lines) ---\n", paste(tail(gdal_out, 40), collapse = "\n"))
  }
  message("Final mosaic saved: ", tif_output)

  # early return: mosaic
  if (return_stage == "mosaic") {
    if (!isTRUE(keep_masked_tiles) && !is.null(masked_output_folder) && dir.exists(masked_output_folder)) {
      unlink(masked_output_folder, recursive = TRUE, force = TRUE)
    }
    unlink(mosaic_tmp, recursive = TRUE, force = TRUE)
    return(normalizePath(tif_output, winslash = "/", mustWork = FALSE))
  }

  # ---------------------------
  # mosaic_clean (restore band names)
  # ---------------------------
  r <- terra::rast(tif_output)
  r <- apply_names_if_possible(r)

  final_clean <- file.path(folder_path, paste0(prefix, year, if (do_mask) "_mosaic_masked.tif" else "_mosaic_nomask.tif"))
  terra::writeRaster(r, final_clean, overwrite = TRUE, gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"))

  if (!isTRUE(keep_masked_tiles) && !is.null(masked_output_folder) && dir.exists(masked_output_folder)) {
    unlink(masked_output_folder, recursive = TRUE, force = TRUE)
  }
  unlink(mosaic_tmp, recursive = TRUE, force = TRUE)

  # cleanup depending on keep_intermediates
  if (keep_intermediates %in% c("none", "mosaic_clean")) {
    # do not keep raw mosaic unless user asked for it
    if (file.exists(tif_output)) unlink(tif_output, force = TRUE)
  }

  if (return_stage == "mosaic_clean") {
    return(normalizePath(final_clean, winslash = "/", mustWork = FALSE))
  }

  # ---------------------------
  # resampled: warp to target CRS + out_res_m
  # If proj_res_m is NULL -> one-step (recommended). If not NULL -> two-step.
  # ---------------------------
  out_res_m <- as.numeric(out_res_m)[1]
  if (!is.finite(out_res_m) || out_res_m <= 0) stop("out_res_m must be a positive number.")

  resampled <- file.path(folder_path, paste0(tools::file_path_sans_ext(basename(final_clean)), "_res", out_res_m, "m.tif"))

  gdal_co <- c("-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES", "-co", "TILED=YES")

  # Choose a safe method for multiband products
  nlyr_in <- terra::nlyr(terra::rast(final_clean))
  warp_method <- if (nlyr_in > 1) "near" else "bilinear"

  if (is.null(proj_res_m)) {
    # one-step: reproject+resample directly to out_res_m
    cmd_out <- c(
      "-overwrite",
      "-t_srs", paste0("EPSG:", crs_target),
      "-tr", as.character(out_res_m), as.character(out_res_m),
      "-r", warp_method,
      "-dstnodata", dst_nodata,
      gdal_co,
      shQuote(final_clean),
      shQuote(resampled)
    )

    out_txt <- system2(gdalwarp_path, args = cmd_out, stdout = TRUE, stderr = TRUE)
    st <- attr(out_txt, "status"); if (is.null(st)) st <- 0L
    if (st != 0L || !file.exists(resampled)) {
      stop("Resample/warp failed or output not created: ", resampled,
           "\n--- GDAL output (last lines) ---\n", paste(tail(out_txt, 40), collapse = "\n"))
    }

  } else {
    # two-step (only if user really wants it)
    proj_res_m <- as.numeric(proj_res_m)[1]
    if (!is.finite(proj_res_m) || proj_res_m <= 0) stop("proj_res_m must be a positive number.")

    reprojected <- file.path(folder_path, paste0(tools::file_path_sans_ext(basename(final_clean)), "_proj.tif"))

    cmd_proj <- c(
      "-overwrite",
      "-t_srs", paste0("EPSG:", crs_target),
      "-tr", as.character(proj_res_m), as.character(proj_res_m),
      "-r", warp_method,
      "-dstnodata", dst_nodata,
      gdal_co,
      shQuote(final_clean),
      shQuote(reprojected)
    )

    proj_txt <- system2(gdalwarp_path, args = cmd_proj, stdout = TRUE, stderr = TRUE)
    stp <- attr(proj_txt, "status"); if (is.null(stp)) stp <- 0L
    if (stp != 0L || !file.exists(reprojected)) {
      stop("Reprojection failed or output not created: ", reprojected,
           "\n--- GDAL output (last lines) ---\n", paste(tail(proj_txt, 40), collapse = "\n"))
    }

    cmd_res <- c(
      "-overwrite",
      "-tr", as.character(out_res_m), as.character(out_res_m),
      "-r", warp_method,
      "-dstnodata", dst_nodata,
      gdal_co,
      shQuote(reprojected),
      shQuote(resampled)
    )

    res_txt <- system2(gdalwarp_path, args = cmd_res, stdout = TRUE, stderr = TRUE)
    strr <- attr(res_txt, "status"); if (is.null(strr)) strr <- 0L
    if (strr != 0L || !file.exists(resampled)) {
      stop("Resample failed or output not created: ", resampled,
           "\n--- GDAL output (last lines) ---\n", paste(tail(res_txt, 40), collapse = "\n"))
    }

    if (keep_intermediates == "none" && file.exists(reprojected)) unlink(reprojected, force = TRUE)
  }

  message("Resampled to ", out_res_m, "m: ", resampled)

  # read + rewrite to enforce integer + restore names
  r_final <- terra::rast(resampled)
  r_final <- apply_names_if_possible(r_final)

  r_final <- floor(r_final)
  terra::NAflag(r_final) <- -9999
  terra::writeRaster(
    r_final, resampled, overwrite = TRUE, datatype = "INT2S",
    gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES", "NAflag=-9999")
  )

  # final cleanup: keep only final file by default
  if (keep_intermediates == "none") {
    if (file.exists(final_clean)) unlink(final_clean, force = TRUE)
  }

  normalizePath(resampled, winslash = "/", mustWork = FALSE)
}
