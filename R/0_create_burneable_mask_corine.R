#' Generate Burnable Mask from CORINE Rasters
#'
#' @description
#' This function reclassifies CORINE land cover rasters into binary burnable masks using a predefined set of vegetation-related land cover classes.
#' It masks the rasters using an Iberian Peninsula shapefile, optionally reprojects them to EPSG:3035 and/or EPSG:4326, and polygonizes the result using GDAL tools.
#' The raster masks can be used in GEE for masking the Landsat composites.
#'
#' The reclassification uses the following CORINE classes (corine_code / description):
#' - 244: Agro-forestry areas (22)
#' - 311: Broad-leaved forest (23)
#' - 312: Coniferous forest (24)
#' - 313: Mixed forest (25)
#' - 321: Natural grasslands (26)
#' - 322: Moors and heathland (27)
#' - 323: Sclerophyllous vegetation (28)
#' - 324: Transitional woodland-shrub (29)
#' - 333: Sparsely vegetated areas (32)
#' - 334: Burnt areas (33)
#'
#' These classes are assigned a gridcode from 22 to 33 in the final output shapefile.
#' The .csv file with gridcode and CORINE classes correspondences is provided in the README folder (clc_legend_gridcode.csv)
#' (if `vectorize = TRUE`).
#'
#' @name generate_burnable_mask
#' @rdname generate_burnable_mask
#'
#' @param corine_rasters Named list of CORINE raster file paths (e.g., one per year).
#' @param peninsula_shapefile Path to the Iberian Peninsula shapefile used to mask the extent.
#' @param output_raster_dir Directory where binary and class rasters will be written.
#' @param output_vector_dir Directory where vector shapefiles will be saved (if `vectorize = TRUE`).
#' @param gdalwarp_path Full path to the GDAL `gdalwarp` executable.
#' @param python_exe Full path to the Python executable.
#' @param gdal_polygonize_script Full path to `gdal_polygonize.py`.
#' @param reproject Logical. Whether to reproject rasters to EPSG:3035. Default is TRUE.
#' @param res Numeric. Output resolution in meters (used when reprojecting). Default is 90.
#' @param to_wgs84 Logical. If TRUE, generates an additional raster version reprojected to EPSG:4326.
#' @param vectorize Logical. Whether to convert the reclassified raster to polygon shapefile. Default is TRUE.
#'
#' @return No R object is returned, but the following files are written:
#' - Binary raster mask (`burneable_mask_binary_<year>_ETRS89.tif`)
#' - Optional EPSG:4326 version (`..._wgs84.tif`)
#' - Multiclass raster with selected CORINE codes (`burneable_classes_def1_<year>_ETRS89.tif`)
#' - Polygon shapefile with selected CORINE classes (`burneable_classes_def1_<year>_ETRS89.shp`)
#'
#' @examples
#' \dontrun{
#' corine_rasters <- list(
#'   "corine06" = "data/CORINE_2006.tif"
#' )
#'
#' generate_burnable_mask(
#'   corine_rasters = corine_rasters,
#'   peninsula_shapefile = "data/peninsula_ib.shp",
#'   output_raster_dir = "output/rasters",
#'   output_vector_dir = "output/vectors",
#'   gdalwarp_path = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe",
#'   python_exe = "C:/ProgramData/anaconda3/python.exe",
#'   gdal_polygonize_script =
#'   "C:/ProgramData/anaconda3/Lib/site-packages/osgeo/scripts/gdal_polygonize.py",
#'   reproject = TRUE,
#'   res = 90,
#'   to_wgs84 = TRUE,
#'   vectorize = TRUE
#' )
#' }
#'
#' @importFrom terra rast crop mask project freq writeRaster ncell
#' @importFrom raster raster extent
#' @importFrom sf st_read st_transform st_crs st_make_valid
#' @importFrom glue glue
#' @importFrom stringr str_detect str_replace
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
  # Reclasificación CORINE → 1 (vegetación seleccionada) o 0
  reclass_matrix <- matrix(c(
    1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0,
    11, 0, 12, 0, 13, 0, 14, 0, 15, 0, 16, 0, 17, 0, 18, 0, 19, 0, 20, 0,
    21, 0, 22, 1, 23, 1, 24, 1, 25, 1, 26, 1, 27, 1, 28, 1, 29, 1,
    30, 0, 31, 0, 32, 1, 33, 1, 34, 0, 35, 0, 36, 0, 37, 0,
    38, 0, 39, 0, 40, 0, 41, 0, 42, 0, 43, 0, 44, 0, 48, 0
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
    masked[] <- vals

    class_r <- terra::rast(masked)
    out_reproj <- file.path(output_raster_dir, paste0("burneable_classes_def1_", name, "_ETRS89.tif"))
    if (reproject) {
      tmp_unproj <- tempfile(fileext = ".tif")
      terra::writeRaster(class_r, tmp_unproj, overwrite = TRUE)
      cmd <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:3035" -tr {res} {res} -r near -co "COMPRESS=LZW" "{tmp_unproj}" "{out_reproj}"')
      system(cmd)
    } else {
      terra::writeRaster(class_r, out_reproj, overwrite = TRUE, datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))
    }
    message("Saved reprojected raster with selected CORINE classes: ", out_reproj)

    if (vectorize && file.exists(python_exe) && file.exists(gdal_polygonize_script)) {
      out_shp <- file.path(output_vector_dir, paste0("burneable_classes_def1_", name, "_ETRS89.shp"))
      cmd_vec <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{out_reproj}" -f "ESRI Shapefile" "{out_shp}" DN')
      system(cmd_vec)
      message("Vectorized: ", out_shp)
    }

    binary_vals <- ifelse(!is.na(vals), 1, NA)
    masked[] <- binary_vals
    binary_r <- terra::rast(masked)
    output_r <- file.path(output_raster_dir, paste0("burneable_mask_binary_", name, "_ETRS89.tif"))
    if (reproject) {
      tmp_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(binary_r, tmp_bin, overwrite = TRUE)
      cmd_bin <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:3035" -tr {res} {res} -r near -dstnodata 0 -co "COMPRESS=LZW" "{tmp_bin}" "{output_r}"')
      system(cmd_bin)
    } else {
      terra::writeRaster(binary_r, output_r, overwrite = TRUE, datatype = "INT1U", NAflag = 0, gdal = c("COMPRESS=LZW"))
    }
    message("Saved binary burnable mask: ", output_r)

    if (vectorize && file.exists(python_exe) && file.exists(gdal_polygonize_script)) {
      out_shp_bin <- file.path(output_vector_dir, paste0("burneable_mask_binary_", name, "_ETRS89.shp"))
      cmd_bin <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{output_r}" -f "ESRI Shapefile" "{out_shp_bin}" DN')
      system(cmd_bin)
      message("Vectorized binary mask: ", out_shp_bin)
    }

    if (to_wgs84) {
      output_wgs <- file.path(output_raster_dir, paste0("burneable_mask_binary_", name, "_wgs84.tif"))
      if (file.exists(output_wgs)) file.remove(output_wgs)
      cmd_wgs <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:4326" -tr 0.0003 0.0003 -r near -dstnodata 0 -co "COMPRESS=LZW" "{output_r}" "{output_wgs}"')
      system(cmd_wgs)
      message("Reprojected to WGS84: ", output_wgs)

      if (vectorize && file.exists(python_exe) && file.exists(gdal_polygonize_script)) {
        output_shp_wgs <- file.path(output_vector_dir, paste0("burneable_mask_binary_", name, "_wgs84.shp"))
        if (file.exists(output_shp_wgs)) {
          shp_base <- tools::file_path_sans_ext(output_shp_wgs)
          file.remove(paste0(shp_base, c(".shp", ".shx", ".dbf", ".prj", ".cpg")))
        }
        cmd_shp <- glue::glue('"{python_exe}" "{gdal_polygonize_script}" "{output_wgs}" -f "ESRI Shapefile" "{output_shp_wgs}" DN')
        system(cmd_shp)
        message("Vectorized binary mask in WGS84: ", output_shp_wgs)
      }
    }
  }
}
