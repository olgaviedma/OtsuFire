#' Build an annual change-index mosaic
#'
#' @description
#' Build a masked mosaic from raster tiles representing a burned-area change
#' index, such as dNBR or RBR. The function verifies that all input tiles have
#' matching coordinate reference systems, spatial resolutions, grid alignment,
#' and band counts. It then merges the tiles, clips and masks the mosaic to the
#' area of interest, optionally applies a user-defined lower bound to valid
#' index values, and writes the resulting GeoTIFF. It returns a compact,
#' structured summary of the output.
#'
#' @param folder_path Character. Directory containing the raster tiles to mosaic.
#' @param mask An `sf` object, a `terra::SpatVector`, or a path to a vector file
#'   used to crop and mask the mosaic to the area of interest.
#' @param year Optional integer or character scalar identifying the target year.
#'   If `NULL`, the function tries to infer it from `folder_path` or from the
#'   first tile name.
#' @param raster_pattern Character. Pattern used to list candidate tiles inside
#'   `folder_path`.
#' @param output_dir Character. Directory where the final mosaic will be written.
#'   Defaults to `folder_path`.
#' @param mosaic_name Optional character scalar. Stable output name to use for
#'   the mosaic. If `NULL`, a name is derived from the first tile.
#' @param nodata_value Numeric scalar used as NoData when reading and writing
#'   rasters.
#' @param tol Numeric scalar. Tolerance used when checking tile resolution and
#'   grid origin compatibility.
#' @param cap_below Logical. Whether values below `lower_cap` should be replaced
#'   by `lower_cap` before writing the mosaic.
#' @param lower_cap Optional numeric scalar. Lower bound applied when
#'   `cap_below = TRUE`.
#' @param method Character. Mosaic method. Supported values are `"merge"` and
#'   `"vrt"`.
#' @param overwrite Logical. Whether an existing output mosaic may be replaced.
#' @param keep_temporary_files Logical. Whether to keep the temporary VRT when
#'   `method = "vrt"`.
#'
#' @return A named list with class `change_index_mosaic_result`, containing
#'   `mosaic_path`, `tile_inventory`, `method_used`, `mask_summary`,
#'   `mosaic_diagnostics`, and `temporary_vrt`.
#'
#' @importFrom sf st_crs st_read st_transform
#' @importFrom terra NAflag crs crop datatype ifel mask merge nlyr origin project
#' @importFrom terra rast res vect writeRaster
#' @importFrom tools file_path_sans_ext
#' @importFrom utils glob2rx
#' @family workflow
#' @export
change_index_mosaic <- function(
    folder_path,
    mask,
    year = NULL,
    raster_pattern = "IBERIAN_MinMin_all_year_*.tif",
    output_dir = folder_path,
    mosaic_name = NULL,
    nodata_value = -9999,
    tol = 0.0000001,
    cap_below = TRUE,
    lower_cap = -1000,
    method = c("merge", "vrt"),
    overwrite = FALSE,
    keep_temporary_files = FALSE
) {
  method <- match.arg(method)

  if (!dir.exists(folder_path)) {
    stop("'folder_path' does not exist: ", folder_path, call. = FALSE)
  }
  if (!is.logical(cap_below) || length(cap_below) != 1L) {
    stop("'cap_below' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(keep_temporary_files) || length(keep_temporary_files) != 1L) {
    stop("'keep_temporary_files' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(lower_cap) && (!is.numeric(lower_cap) || length(lower_cap) != 1L)) {
    stop("'lower_cap' must be NULL or a single numeric value.", call. = FALSE)
  }

  output_dir <- output_dir %||% folder_path
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mask_source <- if (is.character(mask)) normalizePath(mask, winslash = "/", mustWork = FALSE) else class(mask)[1]
  mask_vect <- .of_as_mask_vector(mask)

  raster_files <- list.files(folder_path, pattern = utils::glob2rx(raster_pattern), full.names = TRUE)
  raster_basenames <- basename(raster_files)
  raster_files <- raster_files[!grepl("_mosaic|_res\\d+m|_proj", raster_basenames, ignore.case = TRUE)]
  raster_files <- sort(raster_files)

  if (length(raster_files) == 0L) {
    stop("No rasters found in ", folder_path, " using pattern ", raster_pattern, call. = FALSE)
  }

  ref_rast <- terra::rast(raster_files[1])
  terra::NAflag(ref_rast) <- nodata_value

  ref_crs <- terra::crs(ref_rast)
  ref_res <- terra::res(ref_rast)
  ref_org <- terra::origin(ref_rast)
  ref_nlyr <- terra::nlyr(ref_rast)
  ref_names <- names(ref_rast)
  ref_dtype <- terra::datatype(ref_rast)

  mask_reprojected <- FALSE
  if (terra::crs(mask_vect) != ref_crs) {
    mask_vect <- terra::project(mask_vect, ref_crs)
    mask_reprojected <- TRUE
  }

  if (length(raster_files) > 1L) {
    for (f in raster_files[-1]) {
      r <- terra::rast(f)
      terra::NAflag(r) <- nodata_value

      if (terra::nlyr(r) != ref_nlyr) {
        stop("Band count mismatch in ", basename(f), ".", call. = FALSE)
      }
      if (terra::crs(r) != ref_crs) {
        stop("CRS mismatch in ", basename(f), ".", call. = FALSE)
      }
      if (!isTRUE(all.equal(terra::res(r), ref_res, tolerance = tol))) {
        stop("Resolution mismatch in ", basename(f), ".", call. = FALSE)
      }
      if (!isTRUE(all.equal(terra::origin(r), ref_org, tolerance = tol))) {
        stop("Grid origin mismatch in ", basename(f), ".", call. = FALSE)
      }
    }
  }

  year_label <- .of_resolve_mosaic_year(year = year, folder_path = folder_path, raster_files = raster_files)
  base_name <- .of_resolve_mosaic_name(
    mosaic_name = mosaic_name,
    raster_file = raster_files[1],
    year_label = year_label
  )
  mosaic_output <- file.path(output_dir, paste0(base_name, "_mosaic.tif"))

  if (file.exists(mosaic_output) && !isTRUE(overwrite)) {
    stop("Output mosaic already exists and overwrite = FALSE: ", mosaic_output, call. = FALSE)
  }

  vrt_path <- NULL
  if (identical(method, "merge")) {
    rasters <- lapply(raster_files, function(f) {
      r <- terra::rast(f)
      terra::NAflag(r) <- nodata_value
      r
    })

    mos <- rasters[[1]]
    if (length(rasters) > 1L) {
      for (i in 2:length(rasters)) {
        mos <- terra::merge(mos, rasters[[i]])
      }
    }
  } else {
    # method = "vrt" needs gdalUtilities (Suggests). Guard so the merge
    # path keeps working on machines without the optional dependency.
    if (!requireNamespace("gdalUtilities", quietly = TRUE)) {
      stop(
        "method = \"vrt\" requires the 'gdalUtilities' package. ",
        "Install it or use method = \"merge\".",
        call. = FALSE
      )
    }
    vrt_path <- file.path(output_dir, paste0("temp_mosaic_", year_label, ".vrt"))
    gdalUtilities::gdalbuildvrt(
      gdalfile = raster_files,
      output.vrt = vrt_path,
      overwrite = TRUE,
      srcnodata = nodata_value,
      vrtnodata = nodata_value,
      hidenodata = TRUE
    )
    mos <- terra::rast(vrt_path)
    terra::NAflag(mos) <- nodata_value
  }

  out_crop <- terra::crop(mos, mask_vect)
  out_mask <- terra::mask(out_crop, mask_vect)

  if (isTRUE(cap_below) && !is.null(lower_cap)) {
    out_mask <- terra::ifel(!is.na(out_mask) & out_mask < lower_cap, lower_cap, out_mask)
  }

  if (!is.null(ref_names) && length(ref_names) == terra::nlyr(out_mask)) {
    names(out_mask) <- ref_names
  }

  terra::NAflag(out_mask) <- nodata_value

  force_dtype <- NULL
  if (any(grepl("U$", ref_dtype))) {
    force_dtype <- "INT4S"
  }

  wopt <- list(
    NAflag = nodata_value,
    gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES")
  )
  if (!is.null(force_dtype)) {
    wopt$datatype <- force_dtype
  }

  terra::writeRaster(out_mask, mosaic_output, overwrite = overwrite, wopt = wopt)

  temporary_vrt_out <- NULL
  if (!is.null(vrt_path)) {
    if (isTRUE(keep_temporary_files)) {
      temporary_vrt_out <- vrt_path
    } else if (file.exists(vrt_path)) {
      unlink(vrt_path, force = TRUE)
    }
  }

  tile_inventory <- data.frame(
    tile_index = seq_along(raster_files),
    file = normalizePath(raster_files, winslash = "/", mustWork = FALSE),
    basename = basename(raster_files),
    stringsAsFactors = FALSE
  )

  result <- list(
    mosaic_path = normalizePath(mosaic_output, winslash = "/", mustWork = FALSE),
    tile_inventory = tile_inventory,
    method_used = method,
    mask_summary = list(
      mask_source = mask_source,
      mask_reprojected = mask_reprojected,
      target_crs = ref_crs
    ),
    mosaic_diagnostics = list(
      n_tiles = length(raster_files),
      year = year_label,
      raster_pattern = raster_pattern,
      nodata_value = nodata_value,
      tolerance = tol,
      cap_below = cap_below,
      lower_cap = lower_cap,
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
    ),
    temporary_vrt = temporary_vrt_out
  )

  class(result) <- c("change_index_mosaic_result", class(result))
  result
}

# --- internal helpers (mosaic stage) -----------------------------------

.of_as_mask_vector <- function(mask) {
  if (inherits(mask, "SpatVector")) {
    return(mask)
  }
  if (inherits(mask, "sf")) {
    return(terra::vect(mask))
  }
  if (is.character(mask) && length(mask) == 1L) {
    if (!file.exists(mask)) {
      stop("'mask' path does not exist: ", mask, call. = FALSE)
    }
    return(terra::vect(mask))
  }
  stop("'mask' must be an sf object, a SpatVector, or a valid vector path.", call. = FALSE)
}

.of_resolve_mosaic_year <- function(year, folder_path, raster_files) {
  if (!is.null(year) && length(year) == 1L) {
    return(as.character(year))
  }

  from_folder <- regmatches(basename(folder_path), regexpr("(19|20)\\d{2}", basename(folder_path)))
  if (length(from_folder) == 1L && nzchar(from_folder)) {
    return(from_folder)
  }

  first_name <- basename(raster_files[1])
  from_file <- regmatches(first_name, regexpr("(19|20)\\d{2}", first_name))
  if (length(from_file) == 1L && nzchar(from_file)) {
    return(from_file)
  }

  "unknown_year"
}

.of_resolve_mosaic_name <- function(mosaic_name, raster_file, year_label) {
  if (!is.null(mosaic_name) && length(mosaic_name) == 1L && nzchar(mosaic_name)) {
    return(mosaic_name)
  }

  stem <- tools::file_path_sans_ext(basename(raster_file))
  pat <- paste0("^(.*?", year_label, ")_.*$")
  derived <- sub(pat, "\\1", stem)
  if (!identical(derived, stem)) {
    return(derived)
  }

  derived_generic <- sub("^(.*?\\d{4})_.*$", "\\1", stem)
  if (!identical(derived_generic, stem)) {
    return(derived_generic)
  }

  stem
}
