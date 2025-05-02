#' Mosaic Rasters with Burnable Mask and Resample to 90m
#' @description
#' This fucntion applies a shapefile mask to all raster files (tiles) in a single folder, mosaics them by maximum value,
#' and resamples the final mosaic to 90m resolution using GDAL. Output filenames preserve the core identifier
#' from the original raster name, right before the year digits. Intermediate masked rasters are deleted after use.
#'
#' @name mask_mosaic_raster
#' @rdname mask_mosaic_raster
#' @param folder_path Path to the folder containing input rasters to process.
#' @param mask_shapefile Path to the shapefile mask (burnable area).
#' @param gdalwarp_path Path to the GDAL `gdalwarp` executable (for reprojection/resampling).
#' @param crs_target EPSG code for reprojection (default is 3035; ETRS89 LAEA).
#' @param tile_buffer Buffer in meters applied if mask and raster do not initially overlap.
#'
#' @return The path to the final resampled raster.
#'
#' @importFrom terra vect rast crs project ext buffer crop mask writeRaster nlyr NAflag
#' @importFrom gdalUtilities gdalbuildvrt gdalwarp
#' @importFrom glue glue
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_match
#'
#' @examples
#' \dontrun{
#' mask_mosaic_raster(
#'   folder_path = "path/to/your/folder/2020",
#'   mask_shapefile = "data/burneable_mask.shp",
#'   gdalwarp_path = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe"
#' )
#' }
#'
#' @export
mask_mosaic_raster <- function(folder_path,
                                         mask_shapefile,
                                         gdalwarp_path,
                                         crs_target = 3035,
                                         tile_buffer = 1000) {

  if (missing(gdalwarp_path) || is.null(gdalwarp_path) || !file.exists(gdalwarp_path)) {
    stop("gdalwarp_path debe estar definido y apuntar a un ejecutable valido")
  }

  message("Processing folder:", folder_path)

  mask_vect <- terra::vect(mask_shapefile)
  masked_output_folder <- file.path(folder_path, "masked_rasters")
  dir.create(masked_output_folder, showWarnings = FALSE, recursive = TRUE)

  raster_files <- list.files(folder_path, pattern = glob2rx("*.tif"), full.names = TRUE)
  if (length(raster_files) == 0) {
    stop("No rasters found in ", folder_path)
  }

  example_raster <- terra::rast(raster_files[1])
  if (terra::crs(example_raster) != terra::crs(mask_vect)) {
    message("Reprojecting mask to match raster CRS...")
    mask_vect <- terra::project(mask_vect, terra::crs(example_raster))
  }

  mask_raster <- function(raster_path, mask_vect, output_path, buffer_width) {
    r <- terra::rast(raster_path)
    if (!terra::crs(r) == terra::crs(mask_vect)) {
      mask_vect <- terra::project(mask_vect, terra::crs(r))
    }
    if (is.null(terra::intersect(terra::ext(r), terra::ext(mask_vect)))) {
      message("No overlap. Applying buffer...")
      mask_vect <- terra::buffer(mask_vect, width = buffer_width)
    }

    # Verifica el solape tras aplicar buffer
    if (is.null(terra::intersect(terra::ext(r), terra::ext(mask_vect)))) {
      message("No overlap even after buffering. Skipping: ", basename(raster_path))
      return(NA_character_)  # Devuelve NA si no hay solape
    }

    r_masked <- terra::mask(terra::crop(r, mask_vect), mask_vect)
    terra::writeRaster(r_masked, output_path, overwrite = TRUE)
    return(output_path)
  }

  masked_files <- vapply(raster_files, function(r) {
    output_path <- file.path(masked_output_folder, paste0("masked_", tools::file_path_sans_ext(basename(r)), ".tif"))
    mask_raster(r, mask_vect, output_path, tile_buffer)
  }, character(1))

  masked_files <- masked_files[!is.na(masked_files)]
  if (length(masked_files) == 0) stop("No masked rasters generated - all skipped due to lack of overlap.")


  year <- basename(folder_path)
  vrt_path <- file.path(masked_output_folder, paste0("temp_mosaic_", year, ".vrt"))
  gdalUtilities::gdalbuildvrt(masked_files, vrt_path, overwrite = TRUE)
  message("VRT created: ", vrt_path)

  prefix_match <- stringr::str_match(basename(raster_files[1]), "(.*?)([12][0-9]{3})")[, 2]
  prefix <- if (!is.na(prefix_match)) prefix_match else paste0("composite_")

  tif_output <- file.path(folder_path, paste0(prefix, year, "_mosaic.tif"))
  if (file.exists(tif_output)) file.remove(tif_output)

  gdalUtilities::gdalwarp(vrt_path, tif_output, r = "max", dstnodata = "NA",
                          co = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"), overwrite = TRUE)
  message("Final mosaic saved: ", tif_output)

  r <- terra::rast(tif_output)
  if (terra::nlyr(r) == 2) names(r) <- c("rbr", "doy")

  final_clean <- file.path(folder_path, paste0(prefix, year, "_mosaic_masked.tif"))
  terra::writeRaster(r, final_clean, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  unlink(masked_output_folder, recursive = TRUE, force = TRUE)

  reprojected <- file.path(folder_path, paste0(tools::file_path_sans_ext(basename(final_clean)), "_proj.tif"))
  resampled <- file.path(folder_path, paste0(tools::file_path_sans_ext(basename(final_clean)), "_res90m.tif"))

  # Reproyeccion
  cmd_proj_args <- c(
    "-t_srs", paste0("EPSG:", crs_target),
    "-tr", "30", "30",
    "-r", "bilinear",
    "-dstnodata", "0",
    "-co", "COMPRESS=LZW",
    shQuote(final_clean),
    shQuote(reprojected)
  )

  proj_status <- tryCatch({
    system2(gdalwarp_path, args = cmd_proj_args, stdout = NULL, stderr = NULL)
  }, error = function(e) {
    warning("gdalwarp reprojection failed with error: ", conditionMessage(e))
    return(1)
  })

  if (proj_status != 0) {
    warning("gdalwarp reprojection failed with status ", proj_status)
  } else if (!file.exists(reprojected)) {
    warning("Reprojected file was not created.")
  } else {
    message("Reprojected: ", reprojected)
  }

  # Remuestreo
  cmd_resample_args <- c(
    "-tr", "90", "90",
    "-r", "near",
    "-dstnodata", "0",
    "-co", "COMPRESS=LZW",
    shQuote(reprojected),
    shQuote(resampled)
  )

  resample_status <- tryCatch({
    system2(gdalwarp_path, args = cmd_resample_args, stdout = NULL, stderr = NULL)
  }, error = function(e) {
    warning("gdalwarp resample failed with error: ", conditionMessage(e))
    return(1)
  })

  if (resample_status != 0) {
    warning("gdalwarp resample failed with status ", resample_status)
  } else if (!file.exists(resampled)) {
    warning("Resampled file was not created.")
  } else {
    message("Resampled to 90m: ", resampled)
  }


  r_final <- terra::rast(resampled)
  r_final <- floor(r_final)
  terra::NAflag(r_final) <- -9999
  terra::writeRaster(r_final, resampled, overwrite = TRUE, datatype = "INT2S", gdal = c("COMPRESS=LZW", "NAflag=-9999"))

  return(resampled)
}
