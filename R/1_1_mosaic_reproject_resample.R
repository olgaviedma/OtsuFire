#' Mosaic, Reproject, and Resample Burned Area Tiles
#'
#' @description
#' This function creates a seamless mosaic from multiple raster tiles (e.g., RBR and DOY composites),
#' masks them with a shapefile (e.g., Iberian Peninsula), reprojects the mosaic to a specified CRS
#' (default: EPSG:3035), and resamples it to a target resolution (default: 90m).
#'
#' It assumes the input rasters are two-band GeoTIFFs:
#' - Band 1: RBR (Relative Burn Ratio)
#' - Band 2: DOY (Day of Year of fire)
#'
#' @param folder_path Character. Path to the folder containing the raster tiles.
#' @param mask_path Character. Path to a polygon shapefile used to mask the rasters (e.g., burnable area mask).
#' @param year Integer or character. Year label used for naming output files. If NULL, uses folder name.
#' @param raster_pattern Character. Regex pattern to identify input raster tiles. Default: 'IBERIAN_MinMin_all_year_*.tif'.
#' @param crs_target Character. EPSG code string for the target CRS to reproject (default: "EPSG:3035").
#' @param res_target Numeric. Target spatial resolution in meters. Default: 90.
#' @param nodata_value Numeric. NoData value to assign to missing values. Default: -9999.
#'
#' @return Character. Full path to the final resampled GeoTIFF mosaic.
#'
#' @details
#' The function uses `gdalbuildvrt` and `gdalwarp` for efficient VRT-based mosaicking and reprojection.
#' Steps:
#' 1. Apply spatial mask to each raster tile.
#' 2. Create VRT mosaic.
#' 3. Reproject and resample the mosaic with `gdalwarp`.
#' 4. Clean NoData values and assign band names ('rbr', 'doy').
#' 5. Save compressed GeoTIFF output.
#'
#' @importFrom terra rast vect crop mask crs NAflag writeRaster ext ncell
#' @importFrom gdalUtilities gdalbuildvrt gdalwarp
#' @importFrom stringr str_detect str_replace
#'
#' @examples
#' \dontrun{
#' mosaic_reproject_resample(
#'   folder_path = "ZENODO/exdata",
#'   mask_path = "ZENODO/exdata/iberian_peninsula_proj_final.shp",
#'   year = 2012,
#'   raster_pattern = "IBERIAN_MinMin_all_year_2012_*.tif",
#'   crs_target = "EPSG:3035",
#'   res_target = 90
#' )
#' }
#'
#' @export
mosaic_reproject_resample <- function(
    folder_path,
    mask_path,
    year = NULL,
    raster_pattern = "IBERIAN_MinMin_all_year_*.tif",
    crs_target = "EPSG:3035",
    res_target = 90,
    nodata_value = -9999
) {
  if (is.null(year)) year <- basename(folder_path)

  message("[Step 1] Loading and projecting mask...")
  mask_vect <- terra::vect(mask_path)

  message("[Step 2] Listing raster tiles...")
  raster_files <- list.files(folder_path, pattern = glob2rx(raster_pattern), full.names = TRUE)
  raster_files <- raster_files[!grepl("_mosaic|_res\\d+m|_proj", raster_files)]
  if (length(raster_files) == 0) stop("No rasters found in ", folder_path)

  # Project mask only once to match first raster
  ref_rast <- terra::rast(raster_files[1])
  if (!terra::crs(mask_vect) == terra::crs(ref_rast)) {
    mask_vect <- terra::project(mask_vect, terra::crs(ref_rast))
  }

  masked_output_folder <- file.path(folder_path, "masked_rasters")
  dir.create(masked_output_folder, showWarnings = FALSE)

  # Optimized masking
  message("[Step 3] Masking rasters...")
  masked_files <- vapply(raster_files, function(r) {
    r_rast <- terra::rast(r)

    # Verifica solapamiento espacial
    if (is.null(terra::intersect(terra::ext(r_rast), terra::ext(mask_vect)))) {
      return(NA_character_)
    }

    # Optimizado: crop + mask por separado
    r_crop <- terra::crop(r_rast, mask_vect)
    r_mask <- terra::mask(r_crop, mask_vect)

    output_path <- file.path(masked_output_folder, paste0("masked_", basename(r)))
    terra::writeRaster(
      r_mask, output_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "TILED=YES")
    )
    return(output_path)
  }, FUN.VALUE = character(1))


  masked_files <- masked_files[!is.na(masked_files)]

  message("[Step 4] Creating VRT...")
  vrt_path <- file.path(masked_output_folder, paste0("temp_mosaic_", year, ".vrt"))
  gdalUtilities::gdalbuildvrt(
    gdalfile = masked_files,
    output.vrt = vrt_path,
    overwrite = TRUE
  )

  message("[Step 5] Warping mosaic to ", crs_target, " with resolution ", res_target, "m...")
  mosaic_output <- file.path(folder_path, paste0("IBERIAN_MinMin_all_year_", year, "_mosaic_res", res_target, "m.tif"))

  gdalUtilities::gdalwarp(
    srcfile = vrt_path,
    dstfile = mosaic_output,
    t_srs = crs_target,
    tr = c(res_target, res_target),
    tap = TRUE,
    r = "max",
    co = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"),
    srcnodata = nodata_value,
    dstnodata = nodata_value,
    overwrite = TRUE
  )

  message("[Step 6] Naming bands (if 2 layers)...")
  mosaic_raster <- terra::rast(mosaic_output)

  # Reemplazar -Inf y Inf por NA
  mosaic_raster[is.infinite(mosaic_raster)] <- NA

  # Detectar o asignar NAflag
  na_flag <- terra::NAflag(mosaic_raster)
  if (is.na(na_flag)) na_flag <- -9999  # asignar si no existe

  if (terra::nlyr(mosaic_raster) == 2) {
    names(mosaic_raster) <- c("rbr", "doy")
    tmp_named <- file.path(folder_path, paste0("tmp_", basename(mosaic_output)))
    terra::writeRaster(
      mosaic_raster, tmp_named,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "TILED=YES", paste0("NAflag=", na_flag))
    )
    file.rename(tmp_named, mosaic_output)
  }



  message("[Step 7] Cleaning temporary files...")
  unlink(vrt_path, force = TRUE)
  unlink(masked_output_folder, recursive = TRUE)

  message("[Done] Final mosaic saved at: ", mosaic_output)
  return(mosaic_output)
}
