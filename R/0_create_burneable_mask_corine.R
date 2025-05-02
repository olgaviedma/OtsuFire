#' Generate Burnable Mask from CORINE Rasters
#' @description
#' Reclassifies CORINE land cover rasters into binary burnable masks based on predefined vegetation classes,
#' masks them using an Iberian Peninsula shapefile, reprojects them to EPSG:3035 and/or EPSG:4326, and optionally vectorizes the result.
#'
#' @name create_burneable_mask_corine
#' @rdname create_burneable_mask_corine
#' @param corine_rasters Named list of CORINE raster paths.
#' @param peninsula_shapefile Path to the Iberian Peninsula shapefile.
#' @param output_raster_dir Output directory for raster files.
#' @param output_vector_dir Output directory for vector shapefiles.
#' @param gdalwarp_path Path to GDAL gdalwarp executable.
#' @param python_exe Path to Python executable.
#' @param gdal_polygonize_script Path to gdal_polygonize.py script.
#' @param reproject Logical. Whether to reproject to EPSG:3035. Default is TRUE.
#' @param res Resolution (meters for EPSG:3035). Default is 90.
#' @param to_wgs84 Logical. Whether to create an EPSG:4326 version using gdalwarp. Default is TRUE.
#' @param vectorize Logical. Whether to polygonize final rasters using GDAL. Default is TRUE.
#'
#' @return Nothing is returned, but output files are written to disk.
#'
#' @importFrom terra rast crop mask project freq writeRaster
#' @importFrom raster raster extent
#' @importFrom sf st_read st_transform st_crs st_make_valid
#' @importFrom glue glue
#' @export

generate_burnable_mask <- function(
    corine_rasters,
    peninsula_shapefile,
    output_raster_dir,
    output_vector_dir,
    gdalwarp_path,
    python_exe,
    gdal_polygonize_script,
    reproject = TRUE,
    res = 90,
    to_wgs84 = TRUE,
    vectorize = TRUE
) {
  reclass_matrix <- matrix(c(
    1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0,
    11, 0, 12, 0, 13, 0, 14, 0, 15, 0, 16, 0, 17, 0, 18, 0, 19, 0, 20, 0,
    21, 0, 22, 1, 23, 1, 24, 1, 25, 1, 26, 1, 27, 1, 28, 1, 29, 1, 30, 1,
    31, 1, 32, 1, 33, 1, 34, 0, 35, 0, 36, 0, 37, 0, 38, 0, 39, 0, 40, 0,
    41, 0, 42, 0, 43, 0, 44, 0, 48, 0
  ), ncol = 2, byrow = TRUE)

  selected_original_classes <- reclass_matrix[reclass_matrix[, 2] == 1, 1]
  peninsula <- sf::st_read(peninsula_shapefile, quiet = TRUE)
  for (name in names(corine_rasters)) {
    raster_path <- corine_rasters[[name]]
    corine_r <- raster::raster(raster_path)
    peninsula_proj <- sf::st_transform(peninsula, crs = sf::st_crs(corine_r)) |> sf::st_make_valid()

    cropped <- raster::crop(corine_r, raster::extent(peninsula_proj))
    masked <- raster::mask(cropped, peninsula_proj)
    if (all(is.na(masked[]))) {
      message("No valid values in selected CORINE classes for: ", name)
      next
    }
    vals <- masked[]
    vals[!vals %in% selected_original_classes] <- NA
    binary_vals <- ifelse(!is.na(vals), 1, NA)
    masked[] <- binary_vals

    binary_r <- terra::rast(masked)
    if (reproject) {
      binary_r <- terra::project(binary_r, "EPSG:3035", method = "near", res = res)
    }

    output_r <- file.path(output_raster_dir, paste0("burneable_mask_binary_", name, "_ETRS89.tif"))
    terra::writeRaster(binary_r, output_r, overwrite = TRUE,
                       datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
    message("Saved binary burnable mask: ", output_r)

    if (to_wgs84) {
      output_wgs <- file.path(output_raster_dir, paste0("burneable_mask_binary_", name, "_wgs84.tif"))
      cmd <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:4326" -tr 0.0003 0.0003 -r near -dstnodata 0 -co "COMPRESS=LZW" "{output_r}" "{output_wgs}"')
      system(cmd)
      message("Reprojected to WGS84: ", output_wgs)
    }

    final_r <- terra::project(terra::rast(masked), "EPSG:3035", method = "near", res = res)
    out_reproj <- file.path(output_raster_dir, paste0("burneable_classes_def1_", name, "_ETRS89.tif"))
    terra::writeRaster(final_r, out_reproj, overwrite = TRUE, datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))

    message("Saved reprojected raster: ", out_reproj)
    if (vectorize) {
      out_shp <- file.path(output_vector_dir, paste0("burneable_classes_def1_", name, "_ETRS89.shp"))
      cmd_vec <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{out_reproj}" -f "ESRI Shapefile" "{out_shp}" DN')
      system(cmd_vec)
      message("Vectorized: ", out_shp)
    }
  }
}
