#' Build a mosaic from raster tiles
#'
#' @title Build a mosaic from raster tiles
#' @description
#' Build a masked mosaic from raster tiles stored in a folder. The function
#' checks that all tiles share CRS, resolution, origin alignment, and band
#' count, merges them with `terra::merge()`, crops and masks the result to the
#' input mask, optionally caps values below a lower bound, writes the final
#' GeoTIFF, and returns its path.
#' @param folder_path Character scalar. Folder containing the raster tiles to
#'   mosaic.
#' @param mask_path Character scalar. Path to the vector mask used to crop and
#'   mask the mosaic.
#' @param raster_pattern Character scalar. Pattern used to list candidate raster
#'   tiles inside `folder_path`.
#' @param nodata_value Numeric scalar used as NoData when reading and writing
#'   rasters.
#' @param tol Numeric scalar. Tolerance used when checking tile resolution and
#'   origin alignment.
#' @param cap_below Logical scalar. Whether values below `lower_cap` should be
#'   replaced by `lower_cap`.
#' @param lower_cap Numeric scalar. Lower bound applied when `cap_below = TRUE`.
#' @param method Character scalar. Mosaic method. Currently only `"merge"` is
#'   implemented.
#' @return Character scalar with the path to the written mosaic GeoTIFF.
#' @examples
#' \dontrun{
#' mosaic_from_tiles(
#'   folder_path = "C:/path/to/tiles",
#'   mask_path = "C:/path/to/mask.shp"
#' )
#' }
#' @noRd
mosaic_from_tiles <- function(
    folder_path,
    mask_path,
    raster_pattern = "IBERIAN_MinMin_all_year_*.tif",
    nodata_value = -9999,
    tol = 1e-7,
    cap_below = TRUE,
    lower_cap = -1000,
    method = "merge" # Reserved for future expansion; only merge implemented.
) {
  method <- match.arg(method)
  stopifnot(
    is.character(folder_path), length(folder_path) == 1L, dir.exists(folder_path),
    is.character(mask_path), length(mask_path) == 1L, file.exists(mask_path)
  )

  if (!is.logical(cap_below) || length(cap_below) != 1) {
    stop("'cap_below' must be TRUE or FALSE.")
  }

  if (!is.null(lower_cap) && (!is.numeric(lower_cap) || length(lower_cap) != 1)) {
    stop("'lower_cap' must be a single numeric value.")
  }

  message("[Step 1] Loading mask...")
  mask_vect <- terra::vect(mask_path)

  message("[Step 2] Listing raster tiles...")
  raster_files <- list.files(folder_path, pattern = utils::glob2rx(raster_pattern), full.names = TRUE)
  raster_basenames <- basename(raster_files)
  raster_files <- raster_files[!grepl("_mosaic|_res\\d+m|_proj", raster_basenames)]
  raster_files <- sort(raster_files)

  if (length(raster_files) == 0) stop("No rasters found in ", folder_path)

  # Reference = first tile
  ref_rast <- terra::rast(raster_files[1])
  terra::NAflag(ref_rast) <- nodata_value

  ref_crs <- terra::crs(ref_rast)
  ref_res <- terra::res(ref_rast)
  ref_org <- terra::origin(ref_rast)
  ref_nlyr <- terra::nlyr(ref_rast)
  ref_names <- names(ref_rast)
  ref_dtype <- terra::datatype(ref_rast)

  # Reproject mask to the raster CRS
  if (terra::crs(mask_vect) != ref_crs) {
    message("[Step 1b] Reprojecting mask to match rasters CRS...")
    mask_vect <- terra::project(mask_vect, ref_crs)
  }

  message("[Step 2b] Validating that all tiles share CRS/resolution/origin alignment/band count...")
  for (f in raster_files[-1]) {
    r <- terra::rast(f)
    terra::NAflag(r) <- nodata_value

    if (terra::nlyr(r) != ref_nlyr) {
      stop("Band count mismatch.\nRef has ", ref_nlyr, " bands, but ", basename(f),
           " has ", terra::nlyr(r), " bands.")
    }

    if (terra::crs(r) != ref_crs) {
      stop("CRS mismatch in: ", basename(f))
    }

    if (!isTRUE(all.equal(terra::res(r), ref_res, tolerance = tol))) {
      stop("Resolution mismatch in: ", basename(f),
           "\nRef res: ", paste(ref_res, collapse = ", "),
           "\nThis:    ", paste(terra::res(r), collapse = ", "))
    }

    if (!isTRUE(all.equal(terra::origin(r), ref_org, tolerance = tol))) {
      stop("Grid origin (alignment) mismatch in: ", basename(f),
           "\nRef origin: ", paste(ref_org, collapse = ", "),
           "\nThis:       ", paste(terra::origin(r), collapse = ", "),
           "\n(This indicates that the tiles do not fall on the same grid.)")
    }
  }

  # Output name
  base_full <- tools::file_path_sans_ext(basename(raster_files[1]))
  base_name <- sub("^(.*?\\d{4})_.*$", "\\1", base_full)
  mosaic_output <- file.path(folder_path, paste0(base_name, "_mosaic.tif"))

  message("[Step 3] Merging tiles with terra::merge() (recommended)...")

  rasters <- lapply(raster_files, function(f) {
    r <- terra::rast(f)
    terra::NAflag(r) <- nodata_value
    r
  })

  mos <- rasters[[1]]
  if (length(rasters) > 1) {
    for (i in 2:length(rasters)) {
      mos <- terra::merge(mos, rasters[[i]])
    }
  }

  message("[Step 4] Cropping + masking mosaic...")
  out_crop <- terra::crop(mos, mask_vect)
  out_mask <- terra::mask(out_crop, mask_vect)

  # Lower cap
  if (isTRUE(cap_below) && !is.null(lower_cap)) {
    message("[Step 4b] Capping values below ", lower_cap, " ...")
    out_mask <- terra::ifel(!is.na(out_mask) & out_mask < lower_cap, lower_cap, out_mask)
  }

  # Preserve band names
  if (!is.null(ref_names) && length(ref_names) == terra::nlyr(out_mask)) {
    names(out_mask) <- ref_names
  }

  terra::NAflag(out_mask) <- nodata_value

  # If any data type is unsigned integer (U), -9999 does not fit
  force_dtype <- NULL
  if (any(grepl("U$", ref_dtype))) force_dtype <- "INT4S"

  message("[Step 5] Writing final GeoTIFF (NoData always = ", nodata_value, ") ...")
  wopt <- list(
    NAflag = nodata_value,
    gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES")
  )
  if (!is.null(force_dtype)) wopt$datatype <- force_dtype

  terra::writeRaster(out_mask, mosaic_output, overwrite = TRUE, wopt = wopt)

  message("[Done] Final mosaic saved at: ", mosaic_output)
  return(mosaic_output)
}
