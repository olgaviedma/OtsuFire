#' Calculate Geometric Metrics and Filter Fire Polygons (Batch Mode)
#'
#' This function calculates geometric descriptors for fire-affected polygons from
#' a list of input shapefiles. Metrics include area, perimeter, bounding box dimensions,
#' and elongation. Optionally, spatial filters can be applied to retain only polygons
#' that meet certain geometric criteria.
#'
#' For each input shapefile, the function saves:
#' - A shapefile with all polygons and computed metrics.
#' - A filtered shapefile with only polygons passing all active filters.
#'
#' If the input polygons contain a `CORINE_CLA` column, geometries are dissolved before metrics are calculated.
#'
#' ## Metrics calculated per polygon:
#' - `area_ha`: Area in hectares.
#' - `bbox_wx`: Width of the axis-aligned bounding box (east-west).
#' - `bbox_hy`: Height of the axis-aligned bounding box (north-south).
#' - `perim_m`: Polygon perimeter in meters.
#' - `p_w_ratio`: Perimeter-to-width ratio (`perim_m / bbox_wx`), proxy for complexity.
#' - `h_w_ratio`: Height-to-width ratio (`bbox_hy / bbox_wx`), proxy for anisotropy.
#' - `mnbbx_wd`: Width of the minimum rotated bounding box.
#' - `mnbbx_ln`: Length of the minimum rotated bounding box.
#' - `mnbbx_el`: Elongation as length-to-width ratio (`mnbbx_ln / mnbbx_wd`).
#'
#' ## Filtering behavior:
#' All filter thresholds are optional. If a parameter is `NULL`, its corresponding filter is skipped.
#' A polygon must satisfy **all active filters** to appear in the filtered output.
#'
#' @param shapefile_paths Character vector of paths to input polygon shapefiles.
#' @param output_dir Directory to save the output shapefiles. If `NULL`, outputs are saved next to input files.
#' @param area_min_ha Minimum area in hectares (`area_ha`). Set to `NULL` to disable. Default: 10.
#' @param bbox_h_min Minimum height of axis-aligned bounding box (`bbox_hy`). Default: 630.
#' @param mnbbx_wd_min Minimum width of rotated bounding box (`mnbbx_wd`). Default: 800.
#' @param p_w_ratio_min Minimum perimeter-to-width ratio (`p_w_ratio`). Default: 4.0.
#' @param h_w_ratio_min Minimum height-to-width ratio (`h_w_ratio`). Default: 0.35.
#'
#' @return A named list with one entry per input shapefile. Each entry contains:
#' \describe{
#'   \item{metrics}{Path to shapefile with all polygons and computed metrics.}
#'   \item{filtered}{Path to shapefile with only filtered polygons.}
#'   \item{polygons_all}{`sf` object with all input polygons and metrics.}
#'   \item{polygons_filtered}{`sf` object with only polygons that passed filters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load example fire polygons and apply geometric filtering
#' burned_dir <- system.file("extdata", package = "OtsuFire")
#' burned_files <- list.files(
#'   burned_dir, pattern = glob2rx("burned_areas_*.shp$"), full.names = TRUE
#' )
#'
#' # Define filtering thresholds
#' result_list <- calculate_polygon_metrics(
#'   shapefile_paths = burned_files,
#'   output_dir = burned_dir,
#'   area_min_ha = 10,
#'   bbox_h_min = 630,
#'   mnbbx_wd_min = 800,
#'   p_w_ratio_min = 4.0,
#'   h_w_ratio_min = 0.35
#' )
#'
#' # Access paths and results
#' print(result_list[[1]]$filtered)
#' }
#'
#' @importFrom sf st_read st_write st_geometry st_bbox st_length st_transform st_as_binary st_minimum_rotated_rectangle st_coordinates st_cast st_area st_make_valid
#' @importFrom purrr map_dbl
#' @importFrom dplyr filter first
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom data.table :=
#' @importFrom magrittr %>%
#' @export
#'
utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":="
))

calculate_polygon_metrics <- function(
    shapefile_paths,
    output_dir = NULL,
    area_min_ha = 10,
    bbox_h_min = 630,
    mnbbx_wd_min = 800,
    p_w_ratio_min = 4.0,
    h_w_ratio_min = 0.35
) {
  results <- list()

  for (shapefile_path in shapefile_paths) {
    message("Processing: ", shapefile_path)

    # Set output directory if not provided
    out_dir <- if (is.null(output_dir)) dirname(shapefile_path) else output_dir
    year <- basename(dirname(shapefile_path))
    shp_name <- tools::file_path_sans_ext(basename(shapefile_path))

    # Load shapefile with error handling
    polys <- tryCatch(
      {
        sf::st_read(shapefile_path, quiet = TRUE)
      },
      error = function(e) {
        warning("Failed to read shapefile: ", shapefile_path, " (", e$message, ")")
        return(NULL)
      }
    )

    # If reading failed, skip to next
    if (is.null(polys)) next

    # Validate geometries
    if (!all(sf::st_is_valid(polys))) {
      polys <- sf::st_make_valid(polys)
    }

    # Detect if it is a CORINE or ECOREGION shapefile and dissolve if needed
    is_corine <- "CORINE_CLA" %in% names(polys)
    is_ecoregion <- "EnS_name" %in% names(polys)

    if (is_corine) {
      message("Dissolving CORINE polygons by geometry...")

      union_geom <- sf::st_union(polys)
      cast_polys <- sf::st_cast(union_geom, "POLYGON", warn = FALSE)

      if (length(cast_polys) == 0) {
        warning("The result of the CORINE dissolve is empty. Skipping shapefile.")
        next
      }

      polys <- sf::st_sf(geometry = cast_polys)

    } else if (is_ecoregion) {
      message("Dissolving ECOREGIONS polygons by geometry...")

      union_geom <- sf::st_union(polys)
      cast_polys <- sf::st_cast(union_geom, "POLYGON", warn = FALSE)

      if (length(cast_polys) == 0) {
        warning("The result of the ECOREGIONS dissolve is empty. Skipping shapefile.")
        next
      }

      polys <- sf::st_sf(geometry = cast_polys)

    } else {
      message("Dissolve not applied.")
    }

    # Area and ID
    polys$area_m2 <- sf::st_area(polys)
    polys$area_ha <- round(as.numeric(polys$area_m2) / 10000, 3)
    polys$uid <- seq_len(nrow(polys))

    # Bounding box metrics
    polys$bbox_wx <- purrr::map_dbl(polys$geometry, ~ sf::st_bbox(.x)["xmax"] - sf::st_bbox(.x)["xmin"])
    polys$bbox_hy <- purrr::map_dbl(polys$geometry, ~ sf::st_bbox(.x)["ymax"] - sf::st_bbox(.x)["ymin"])

    # Perimeter and elongation
    polys$perim_m <- purrr::map_dbl(polys$geometry, ~ as.numeric(sf::st_length(sf::st_cast(.x, "MULTILINESTRING"))))
    polys$p_w_ratio <- polys$perim_m / polys$bbox_wx
    polys$h_w_ratio <- polys$bbox_hy / polys$bbox_wx

    # Rotated bounding box
    polys$ornt_bbox <- sf::st_minimum_rotated_rectangle(polys$geometry)
    bbox_dims <- lapply(polys$ornt_bbox, function(poly) {
      coords <- sf::st_coordinates(sf::st_cast(poly, "LINESTRING"))
      dx <- max(coords[, 1]) - min(coords[, 1])
      dy <- max(coords[, 2]) - min(coords[, 2])
      c(wd = min(dx, dy), ln = max(dx, dy))
    })
    bbox_df <- do.call(rbind, bbox_dims)
    polys$mnbbx_wd <- bbox_df[, "wd"]
    polys$mnbbx_ln <- bbox_df[, "ln"]
    polys$mnbbx_el <- polys$mnbbx_ln / polys$mnbbx_wd

    # Save metrics shapefile
    metrics_path <- file.path(out_dir, paste0(shp_name, "_metrics.shp"))
    polys <- polys %>% dplyr::select(-area_m2, -ornt_bbox)
    sf::st_write(polys, metrics_path, append = FALSE)
    message("Saved metrics shapefile (unfiltered): ", metrics_path)

    # Apply filters
    filtered <- polys
    filter_labels <- c()

    if (!is.null(area_min_ha)) {
      filtered <- filtered[filtered$area_ha >= area_min_ha, ]
      filter_labels <- c(filter_labels, "area")
    }
    if (!is.null(bbox_h_min)) {
      filtered <- filtered[filtered$bbox_hy >= bbox_h_min, ]
      filter_labels <- c(filter_labels, "bbox_h")
    }
    if (!is.null(mnbbx_wd_min)) {
      filtered <- dplyr::filter(filtered, mnbbx_wd >= mnbbx_wd_min + 1e-6)
      filter_labels <- c(filter_labels, "mnbbx_wd")
    }
    if (!is.null(p_w_ratio_min)) {
      filtered <- dplyr::filter(filtered, p_w_ratio >= p_w_ratio_min + 1e-6)
      filter_labels <- c(filter_labels, "p_w_ratio")
    }
    if (!is.null(h_w_ratio_min)) {
      filtered <- dplyr::filter(filtered, h_w_ratio >= h_w_ratio_min + 1e-6)
      filter_labels <- c(filter_labels, "h_w_ratio")
    }

    # Suffix for filtered file
    suffix <- if (length(filter_labels) == 0) {
      "nofilter"
    } else if (length(filter_labels) == 5) {
      "all"
    } else {
      paste(filter_labels, collapse = "_")
    }

    # Save filtered shapefile
    filtered_path <- file.path(out_dir, paste0(shp_name, "_metrics_filt_", suffix, ".shp"))
    sf::st_write(filtered, filtered_path, append = FALSE)
    message("Saved filtered metrics shapefile: ", filtered_path)

    # Store result
    results[[shapefile_path]] <- list(
      metrics = metrics_path,
      filtered = filtered_path,
      polygons_all = polys,
      polygons_filtered = filtered
    )
  }

  return(results)
}
