#' Calculate Geometric Metrics and Filter Fire Polygons (Batch Mode)
#'
#' @description
#' This function calculates geometric descriptors for fire-affected polygons from
#' a list of input shapefiles. Metrics include area, perimeter, bounding box dimensions,
#' and elongation. Optionally, spatial filters can be applied to retain only polygons
#' that meet certain geometric criteria.
#'
#' For each input shapefile, the function saves:
#' - A shapefile with all polygons and computed metrics.
#' - A filtered shapefile with only polygons passing all active filters.
#' - Optionally, joined versions of the above with attributes from the original input shapefile.
#'
#' If the input fire events are subdivided by CORINE or ECOREGION classes, adjacent geometries can be dissolved before metric calculation to reconstruct contiguous burned areas.
#' In that case, a numeric `burned_id` is assigned to each disjoint polygon in the dissolved geometry. This ID is retained in the spatial joins to allow tracking the disaggregated components of each event.
#'
#' ## Metrics calculated per polygon:
#' - `area_ha`: Area in hectares.
#' - `bbox_wx`: Width of the axis-aligned bounding box (east-west).
#' - `bbox_hy`: Height of the axis-aligned bounding box (north-south).
#' - `perim_m`: Polygon perimeter in meters.
#' - `p_w_ratio`: Perimeter-to-width ratio (`perim_m / bbox_wx`), proxy for complexity.
#' - `h_w_ratio`: Height-to-width ratio (`bbox_hy / bbox_wx`), proxy for anisotropy.
#' - `burned_id`: Integer identifier (1 to N) assigned to each polygon in the dissolved output.
#'
#' ## Filtering behavior:
#' All filter thresholds are optional. If a parameter is `NULL`, its corresponding filter is skipped.
#' A polygon must satisfy **all active filters** to appear in the filtered output, unless `filter_logic = "OR"` is used.
#'
#' ## Spatial Join (optional):
#' If the input shapefiles contain attribute data (e.g., CORINE codes, ecoregion classes), these attributes
#' are preserved in the output by performing a spatial join **from the original shapefile to the metrics output**.
#' - `joined_metrics`: polygons from the original shapefile with appended metric fields and `event_id`.
#' - `joined_filtered`: filtered version of the above.
#'
#' @name calculate_polygon_metrics
#' @rdname calculate_polygon_metrics
#'
#' @param shapefile_paths Character vector of paths to input polygon shapefiles.
#' @param output_dir Directory to save the output shapefiles. If `NULL`, outputs are saved next to input files.
#' @param area_min_ha Minimum area in hectares (`area_ha`). Set to `NULL` to disable. Default: 10.
#' @param bbox_h_min Minimum height of axis-aligned bounding box (`bbox_hy`). Default: 630.
#' @param mnbbx_wd_min Minimum width of rotated bounding box (`mnbbx_wd`). Default: 800.
#' @param p_w_ratio_min Minimum perimeter-to-width ratio (`p_w_ratio`). Default: 4.0.
#' @param h_w_ratio_min Minimum height-to-width ratio (`h_w_ratio`). Default: 0.35.
#' @param output_format Character. Output format for saved files: `"shp"` (default) or `"geojson"`.
#' @param filter_logic Logical combination of filters: `"AND"` (default) requires all filters to be met, `"OR"` passes polygons meeting any filter.
#' @param dissolve Logical. If TRUE, dissolve adjacent polygons into contiguous shapes before computing metrics. Default: TRUE.
#' @param overlay_polygons Optional. `sf` object or path to shapefile to be joined by spatial intersection. If `NULL`, no spatial join is performed.
#'
#' @return A named list with one entry per input shapefile. Each entry contains:
#' \describe{
#'   \item{metrics}{Path to shapefile with all polygons and computed metrics.}
#'   \item{filtered}{Path to shapefile with only filtered polygons.}
#'   \item{joined_metrics}{Path to joined shapefile: original polygons + metrics + event ID.}
#'   \item{joined_filtered}{Path to joined filtered shapefile: original polygons + metrics + event ID.}
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
#'   h_w_ratio_min = 0.35,
#'   output_format = "geojson",
#'   filter_logic = "AND",
#'   dissolve = TRUE
#' )
#'
#' # Access filtered and joined outputs
#' print(result_list[[1]]$filtered)
#' print(result_list[[1]]$joined_filtered)
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
#' @importFrom stringr str_detect str_replace
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
    h_w_ratio_min = 0.35,
    output_format = c("shp", "geojson"),
    filter_logic = c("AND", "OR"),
    dissolve = TRUE,
    join_attributes = TRUE,
    overlay_polygons = NULL) {

  output_format <- match.arg(output_format, choices = c("shp", "geojson"))
  filter_logic <- match.arg(filter_logic, choices = c("AND", "OR"))

  results <- list()

  for (shapefile_path in shapefile_paths) {
    message("Processing: ", shapefile_path)

    out_dir <- if (is.null(output_dir)) dirname(shapefile_path) else output_dir
    year <- basename(dirname(shapefile_path))
    shp_name <- tools::file_path_sans_ext(basename(shapefile_path))

    polys <- tryCatch(
      {
        sf::st_read(shapefile_path, quiet = TRUE)
      },
      error = function(e) {
        warning("Failed to read shapefile: ", shapefile_path, " (", e$message, ")")
        return(NULL)
      }
    )

    if (is.null(polys)) next

    if (!all(sf::st_is_valid(polys))) {
      polys <- sf::st_make_valid(polys)
    }

    if (isTRUE(dissolve)) {
      message("Dissolving adjacent polygons...")
      union_geom <- sf::st_union(polys)
      cast_polys <- sf::st_cast(union_geom, "POLYGON", warn = FALSE)
      if (length(cast_polys) == 0) {
        warning("Dissolve result is empty. Skipping.")
        next
      }
      polys <- sf::st_sf(geometry = cast_polys)
      polys$burned_id <- seq_len(nrow(polys))

    } else {
      message("Dissolve not applied.")
    }


    polys$area_m2 <- sf::st_area(polys)
    polys$area_ha <- round(as.numeric(polys$area_m2) / 10000, 3)
    polys$uid <- seq_len(nrow(polys))

    polys$bbox_wx <- purrr::map_dbl(polys$geometry, ~ sf::st_bbox(.x)["xmax"] - sf::st_bbox(.x)["xmin"])
    polys$bbox_hy <- purrr::map_dbl(polys$geometry, ~ sf::st_bbox(.x)["ymax"] - sf::st_bbox(.x)["ymin"])

    polys$perim_m <- purrr::map_dbl(polys$geometry, ~ as.numeric(sf::st_length(sf::st_cast(.x, "MULTILINESTRING"))))
    polys$p_w_ratio <- polys$perim_m / polys$bbox_wx
    polys$h_w_ratio <- polys$bbox_hy / polys$bbox_wx

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

    ext <- if (output_format == "geojson") ".geojson" else ".shp"
    metrics_path <- file.path(out_dir, paste0(shp_name, "_metrics", ext))

    if (output_format == "shp") {
      shp_base <- tools::file_path_sans_ext(metrics_path)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext_i in shp_exts) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }
    } else {
      if (file.exists(metrics_path)) file.remove(metrics_path)
    }

    polys <- polys %>% dplyr::select(-area_m2, -ornt_bbox)
    sf::st_write(polys, metrics_path, append = FALSE, quiet = TRUE)
    message("Saved metrics file (unfiltered): ", metrics_path)

    # FILTERING LOGIC
    m_area     <- if (!is.null(area_min_ha)) polys$area_ha >= area_min_ha else NA
    m_bbox_h   <- if (!is.null(bbox_h_min)) polys$bbox_hy >= bbox_h_min else NA
    m_mnbbx_wd <- if (!is.null(mnbbx_wd_min)) polys$mnbbx_wd >= (mnbbx_wd_min + 1e-6) else NA
    m_p_w      <- if (!is.null(p_w_ratio_min)) polys$p_w_ratio >= (p_w_ratio_min + 1e-6) else NA
    m_h_w      <- if (!is.null(h_w_ratio_min)) polys$h_w_ratio >= (h_w_ratio_min + 1e-6) else NA

    masks <- list(m_area, m_bbox_h, m_mnbbx_wd, m_p_w, m_h_w)
    filter_labels <- c("area", "bbox_h", "mnbbx_wd", "p_w_ratio", "h_w_ratio")
    valid <- !sapply(masks, function(x) all(is.na(x)))
    masks <- masks[valid]
    filter_labels <- filter_labels[valid]

    if (length(masks) == 0) {
      filtered <- polys
      suffix <- "nofilter"
    } else {
      combined_mask <- Reduce(
        if (filter_logic == "AND") `&` else `|`,
        masks
      )
      filtered <- polys[combined_mask, ]
      suffix <- if (length(filter_labels) == 5) {
        "all"
      } else {
        paste(filter_labels, collapse = "_")
      }
    }

    filtered_path <- file.path(out_dir, paste0(shp_name, "_metrics_filt_", suffix, ext))

    if (output_format == "shp") {
      shp_base <- tools::file_path_sans_ext(filtered_path)
      shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
      for (ext_i in shp_exts) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) file.remove(f)
      }
    } else {
      if (file.exists(filtered_path)) file.remove(filtered_path)
    }

    sf::st_write(filtered, filtered_path, append = FALSE, quiet = TRUE)
    message("Saved filtered metrics file: ", filtered_path)


    # Unión espacial con shapefile original solo si join_attributes = TRUE
    if (isTRUE(join_attributes)) {

      message("Joining attributes from original shapefile: ", shapefile_path)

      overlay_polygons <- tryCatch({
        sf::st_read(shapefile_path, quiet = TRUE)
      }, error = function(e) {
        warning("Could not read original shapefile: ", e$message)
        return(NULL)
      })

      if (!is.null(overlay_polygons)) {
        overlay_polygons <- sf::st_make_valid(overlay_polygons)
        overlay_polygons <- sf::st_transform(overlay_polygons, sf::st_crs(polys))

        # Unión espacial sobre todos los polígonos con métricas
        joined <- tryCatch({
          sf::st_join(overlay_polygons, polys, join = sf::st_intersects, left = FALSE)
        }, error = function(e) {
          warning("st_join() failed for original polygons: ", e$message)
          return(NULL)
        })


        if (!is.null(joined)) {
          message("Join result over all polygons: ", nrow(joined), " rows.")
          polys <- joined
        }

        # Unión espacial sobre polígonos filtrados
        joined_filter <- tryCatch({
          sf::st_join(overlay_polygons, filtered, join = sf::st_intersects, left = FALSE)
        }, error = function(e) {
          warning("st_join() failed for filtered polygons: ", e$message)
          return(NULL)
        })


        if (!is.null(joined_filter)) {
          message("Join result over filtered polygons: ", nrow(joined_filter), " rows.")
          filtered <- joined_filter
        }

        # Guardar las versiones con unión espacial
        joined_metrics_path <- file.path(out_dir, paste0(shp_name, "_metrics_joined", ext))
        joined_filtered_path <- file.path(out_dir, paste0(shp_name, "_metrics_filt_", suffix, "_joined", ext))

        # Eliminar si ya existen
        for (f in c(joined_metrics_path, joined_filtered_path)) {
          if (output_format == "shp") {
            shp_base <- tools::file_path_sans_ext(f)
            shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
            for (ext_i in shp_exts) {
              file_i <- paste0(shp_base, ext_i)
              if (file.exists(file_i)) file.remove(file_i)
            }
          } else {
            if (file.exists(f)) file.remove(f)
          }
        }

        # Guardar los shapefiles con unión
        sf::st_write(polys, joined_metrics_path, append = FALSE, quiet = TRUE)
        sf::st_write(filtered, joined_filtered_path, append = FALSE, quiet = TRUE)
        message("Saved joined metrics file: ", joined_metrics_path)
        message("Saved joined filtered file: ", joined_filtered_path)

      }
    }

    # Resultado actualizado
    results[[shapefile_path]] <- list(
      metrics = metrics_path,
      filtered = filtered_path,
      joined_metrics = joined_metrics_path,
      joined_filtered = joined_filtered_path,
      polygons_all = polys,
      polygons_filtered = filtered
    )

  }

  return(results)
}

