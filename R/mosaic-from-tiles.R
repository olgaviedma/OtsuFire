#' @title Build a mosaic from raster tiles
#' @description
#' Build a mosaic from raster tiles stored in a folder and mask it to a study
#' area. The function checks that the tiles have matching coordinate reference
#' systems, spatial resolutions, grid alignment, and band counts. It then
#' merges the tiles, crops and masks the mosaic, optionally applies a lower
#' bound to pixel values, and writes the resulting GeoTIFF.
#' @param folder_path Character scalar. Path to the folder containing the
#'   raster tiles.
#' @param mask_path Character scalar. Path to the vector file defining the
#'   study area used to crop and mask the mosaic.
#' @param raster_pattern Character scalar. Filename pattern used to select
#'   raster tiles within `folder_path`. Default:
#'   `"*.tif"` (all GeoTIFF files in the folder); set a narrower pattern if
#'   the folder contains other rasters.
#' @param nodata_value Numeric scalar. Value used to identify missing data when
#'   reading the tiles and to encode missing data in the output. Default:
#'   `-9999`.
#' @param tol Numeric scalar. Numerical tolerance used when checking tile
#'   resolution and grid alignment. Default: `1e-07`.
#' @param cap_below Logical scalar. Whether to replace values below
#'   `lower_cap` with `lower_cap`. Default: `TRUE`.
#' @param lower_cap Numeric scalar. Lower bound applied when
#'   `cap_below = TRUE`. Choose a value appropriate for the index and its
#'   scaling. Default: `-1000`.
#' @param method Character scalar. Method used to combine the tiles.
#'   Currently, only `"merge"` is supported.
#' @details
#' The function checks the compatibility of the selected tiles before
#' combining them with `terra::merge()`. Tiles must have matching coordinate
#' reference systems, spatial resolutions, grid alignment, and band counts.
#'
#' The merged raster is cropped to the extent of the study area and masked to
#' its geometry. Cells outside the mask are assigned missing values.
#'
#' When `cap_below = TRUE`, non-missing pixel values below `lower_cap` are
#' replaced with the lower bound. For example, with `lower_cap = -1000`, a
#' value of `-1500` becomes `-1000`. Set `cap_below = FALSE` to disable this
#' operation.
#'
#' The filename pattern determines which tiles are included. For an annual
#' mosaic, select tiles from a single year, for example with
#' `"RBR_2020_*.tif"`, or store each year's tiles in a separate folder.
#'
#' Although the default filename pattern refers to RBR, the function can also
#' be used with other raster indices by changing `raster_pattern` and, where
#' appropriate, `lower_cap`. All selected tiles should represent the same
#' index and use the same units and scaling.
#' @return A character scalar containing the path to the output mosaic
#'   GeoTIFF.
#' @examples
#' # Create a temporary folder for two example tiles
#' tmp_dir <- tempfile("mosaic_example_")
#' dir.create(tmp_dir)
#'
#' # Create adjacent raster tiles with matching grids
#' tile1 <- terra::rast(
#'   nrows = 2, ncols = 2,
#'   xmin = 0, xmax = 2,
#'   ymin = 0, ymax = 2,
#'   crs = "EPSG:3035",
#'   vals = 1:4
#' )
#' tile2 <- terra::rast(
#'   nrows = 2, ncols = 2,
#'   xmin = 2, xmax = 4,
#'   ymin = 0, ymax = 2,
#'   crs = "EPSG:3035",
#'   vals = 5:8
#' )
#'
#' terra::writeRaster(
#'   tile1,
#'   file.path(tmp_dir, "RBR_2020_tile_a.tif"),
#'   NAflag = -9999,
#'   overwrite = TRUE
#' )
#' terra::writeRaster(
#'   tile2,
#'   file.path(tmp_dir, "RBR_2020_tile_b.tif"),
#'   NAflag = -9999,
#'   overwrite = TRUE
#' )
#'
#' # Create a polygon mask covering both tiles
#' mask_vect <- terra::vect(
#'   matrix(
#'     c(
#'       0, 0,
#'       4, 0,
#'       4, 2,
#'       0, 2,
#'       0, 0
#'     ),
#'     ncol = 2,
#'     byrow = TRUE
#'   ),
#'   type = "polygons",
#'   crs = "EPSG:3035"
#' )
#'
#' mask_path <- file.path(tmp_dir, "study_area.gpkg")
#' terra::writeVector(
#'   mask_vect,
#'   mask_path,
#'   overwrite = TRUE
#' )
#'
#' # Build the mosaic
#' result_path <- mosaic_from_tiles(
#'   folder_path = tmp_dir,
#'   mask_path = mask_path,
#'   raster_pattern = "RBR_2020_*.tif",
#'   lower_cap = -1000
#' )
#'
#' # Inspect the output
#' result_path
#' terra::rast(result_path)
#'
#' \dontrun{
#' # Build an annual mosaic from your own raster tiles
#' result_path <- mosaic_from_tiles(
#'   folder_path = "path/to/raster_tiles",
#'   mask_path = "path/to/study_area.gpkg",
#'   raster_pattern = "RBR_2020_*.tif",
#'   lower_cap = -1000
#' )
#'
#' # Build a mosaic without applying a lower bound
#' result_path <- mosaic_from_tiles(
#'   folder_path = "path/to/raster_tiles",
#'   mask_path = "path/to/study_area.gpkg",
#'   raster_pattern = "RBR_2020_*.tif",
#'   cap_below = FALSE
#' )
#' }
#' @export
mosaic_from_tiles <- function(
    folder_path,
    mask_path,
    raster_pattern = "*.tif",
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
