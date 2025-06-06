#' Apply a Binary Raster Mask (and Optional Geographic Clip) to a Two-Band Mosaic
#'
#' @description
#' This function applies a binary raster mask to a 2-band burned area mosaic (RBR and DOY).
#' All pixels in the mosaic where the mask does not match `valid_mask_value` are assigned `NA`.
#' Optionally, a shapefile can be provided to crop and mask the mosaic to a specific geographic area
#' (e.g., the Iberian Peninsula).
#'
#' The function ensures that the mask raster is projected and aligned to the mosaic before applying.
#' The output is saved with the same resolution as the input mosaic.
#'
#' @name mask_to_mosaic
#' @rdname mask_to_mosaic
#'
#' @param mosaic_path Character. Path to the input 2-band raster mosaic (e.g., with RBR and DOY).
#' @param mask_raster_path Character. Path to the binary mask raster (e.g., burnable areas).
#' @param output_path Character (optional). Output destination. There are three behaviors:
#'   \itemize{
#'     \item If `NULL`, the output is saved in the same folder as the input mosaic, using the same
#'           base name with suffix `_mosaic_masked_res90m.tif`.
#'     \item If a directory path is provided, the output file will be saved inside that folder,
#'           using the same base name as the input mosaic, with suffix `_mosaic_masked_res90m.tif`.
#'     \item If a full file path ending in `.tif` is provided, that name is used directly as the output.
#'   }
#' @param valid_mask_value Numeric. Value in the mask that defines valid pixels (default is 1).
#' @param crs_target Integer or character. EPSG code for reprojection (only used if needed). Default: 3035.
#' @param shapefile_clip Character (optional). Path to a shapefile to further crop and mask the mosaic.
#'
#' @return Character. Full path to the final masked (and optionally clipped) raster.
#'
#' @details
#' This function is typically used after mosaicking burned area products derived from satellite imagery.
#' It filters areas outside a burnable mask (e.g., based on land cover) and optionally restricts output
#' to a defined geographic region (e.g., the Iberian Peninsula).
#'
#' Input mosaic must have:
#' - Band 1: Relative Burn Ratio (RBR)
#' - Band 2: Day of Year (DOY)
#'
#' Output raster is written in GeoTIFF format with LZW compression.
#'
#' @importFrom terra rast mask crop compareGeom project vect crs writeRaster ext nlyr ncell
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' \dontrun{
#' # Case 1: Let the function generate the output name in the same folder
#' mask_to_mosaic(
#'   mosaic_path = "data/IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif",
#'   mask_raster_path = "data/burneable_mask_corine_ETRS89.tif"
#' )
#'
#' # Case 2: Save the output to a specific folder, name generated automatically
#' mask_to_mosaic(
#'   mosaic_path = "data/IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif",
#'   mask_raster_path = "data/burneable_mask_corine_ETRS89.tif",
#'   output_path = "outputs/"
#' )
#'
#' # Case 3: Define the full output file name explicitly
#' mask_to_mosaic(
#'   mosaic_path = "data/IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif",
#'   mask_raster_path = "data/burneable_mask_corine_ETRS89.tif",
#'   output_path = "outputs/custom_masked_output_2012.tif"
#' )
#' }
#'
#' @export

mask_to_mosaic <- function(
    mosaic_path,
    mask_raster_path,
    output_path = NULL,
    valid_mask_value = 1,
    crs_target = 3035,
    shapefile_clip = NULL
) {
  message("[Step 1] Reading mosaic and mask...")
  rast <- terra::rast(mosaic_path)
  mask <- terra::rast(mask_raster_path)

  # === Lógica robusta para output_path ===
  if (is.null(output_path)) {
    base_name <- tools::file_path_sans_ext(basename(mosaic_path))
    out_dir <- dirname(mosaic_path)
    masked_name <- sub("(.*)_res", "\\1_masked_res", base_name)
    output_name <- file.path(out_dir, paste0(masked_name, ".tif"))

  } else if (dir.exists(output_path)) {
    base_name <- tools::file_path_sans_ext(basename(mosaic_path))
    masked_name <- sub("(.*)_res", "\\1_masked_res", base_name)
    output_name <- file.path(output_path, paste0(masked_name, ".tif"))
  } else {
    output_name <- normalizePath(output_path, mustWork = FALSE)
  }

  message("[Step 2] Aligning CRS and resolution of mask if needed...")
  if (!terra::compareGeom(rast, mask, stopOnError = FALSE)) {
    mask <- terra::project(mask, rast, method = "near")
  }
  mask[mask != valid_mask_value] <- NA

  message("[Step 3] Applying mask to each band...")
  masked_bands <- lapply(1:terra::nlyr(rast), function(i) {
    terra::mask(rast[[i]], mask)
  })
  masked_stack <- terra::rast(masked_bands)
  if (terra::nlyr(masked_stack) == 2) names(masked_stack) <- c("rbr", "doy")

  if (!is.null(shapefile_clip)) {
    message("[Step 4] Cropping to Iberian Peninsula shape...")
    vect_mask <- terra::vect(shapefile_clip)
    if (!terra::crs(vect_mask) == terra::crs(masked_stack)) {
      vect_mask <- terra::project(vect_mask, terra::crs(masked_stack))
    }
    masked_stack <- terra::crop(masked_stack, vect_mask)
    masked_stack <- terra::mask(masked_stack, vect_mask)
  }

  message("[Step 5] Writing output raster...")
  if (file.exists(output_name)) {
    message("[Warning] Output file already exists. Deleting...")
    unlink(output_name)
  }
  terra::writeRaster(masked_stack, output_name, overwrite = TRUE, gdal = c("COMPRESS=LZW"))


  message("[Done] Masked mosaic saved at: ", output_name)
  return(output_name)
}

