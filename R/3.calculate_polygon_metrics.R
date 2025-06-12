#' Calculate Geometric Metrics and Filter Fire Polygons (Batch Mode)
#'
#' @description
#' Calculates geometric descriptors for fire-affected polygons from a list of input shapefiles.
#' Metrics include area, perimeter, bounding box dimensions, rotated bounding box, and shape ratios.
#' Optionally applies spatial filters and performs attribute joins to enrich outputs.
#'
#' For each input shapefile, the function:
#' - Computes metrics per polygon and saves a new shapefile.
#' - Applies spatial filters and saves a filtered version.
#' - Optionally joins metrics back to the original input polygons (if `join_attributes = TRUE`).
#'
#' If `dissolve = TRUE`, adjacent polygons are merged to reconstruct contiguous burned areas.
#' A numeric `burned_id` is then assigned to each resulting polygon and used in subsequent joins.
#'
#' ## Metrics computed per polygon:
#' - `area_ha`: Area in hectares.
#' - `bbox_wx`: Width of axis-aligned bounding box (x).
#' - `bbox_hy`: Height of axis-aligned bounding box (y).
#' - `perim_m`: Perimeter in meters.
#' - `p_w_ratio`: Perimeter-to-width ratio.
#' - `h_w_ratio`: Height-to-width ratio.
#' - `burned_id`: Polygon ID if `dissolve = TRUE`.
#'
#' ## Filter behavior:
#' All thresholds are optional. If a threshold is `NULL`, the corresponding filter is skipped.
#' If `filter_logic = "AND"`, all active filters must be satisfied. If `"OR"`, any filter match is sufficient.
#'
#' ## Spatial join (`join_attributes = TRUE`):
#' - Attributes from the original input shapefile are preserved via intersection-based joins.
#' - Only input polygons with area ? `min_input_area_ha` are used.
#' - The join selects, for each input polygon, the intersecting output polygon with the largest shared area.
#' - Only polygons that intersect filtered outputs are retained in the joined outputs.
#' - Optionally, you can restrict which variables from the original shapefile to retain using `columns_to_keep`.
#'
#' @name calculate_polygon_metrics
#' @rdname calculate_polygon_metrics
#'
#' @param shapefile_paths Character vector. Paths to input shapefiles.
#' @param output_dir Optional. Directory to save outputs. Defaults to the input file's folder.
#' @param area_min_ha Minimum area in hectares (`area_ha`). Default: 10.
#' @param bbox_h_min Minimum height of axis-aligned bounding box (`bbox_hy`). Default: 630.
#' @param mnbbx_wd_min Minimum width of rotated bounding box. Default: 800.
#' @param p_w_ratio_min Minimum perimeter-to-width ratio. Default: 4.0.
#' @param h_w_ratio_min Minimum height-to-width ratio. Default: 0.35.
#' @param output_format Output format for saved files: `"shp"` or `"geojson"`. Default: `"shp"`.
#' @param filter_logic Logical combination of filters: `"AND"` (default) or `"OR"`.
#' @param dissolve Logical. If `TRUE`, dissolve adjacent polygons before computing metrics. Default: TRUE.
#' @param join_attributes Logical. If `TRUE`, joins filtered/dissolved polygons back to the original shapefile. Default: TRUE.
#' @param min_input_area_ha Minimum area (in hectares) for input polygons used in spatial joins. Default: 10.
#' @param overlay_polygons_path Optional. `sf` object or path to shapefile to be used for join instead of the input shapefile. Default: `NULL`.
#' @param columns_to_keep Optional. Character vector of variable names to retain from the original shapefile during spatial join. Default: `NULL` (keep all).
#'
#' @return A named list (one element per input file), each containing:
#' \describe{
#'   \item{metrics}{Path to shapefile with all polygons and computed metrics.}
#'   \item{filtered}{Path to shapefile with only polygons passing the filters.}
#'   \item{joined_metrics}{(If `join_attributes = TRUE`) Path to joined shapefile with original attributes + metrics.}
#'   \item{joined_filtered}{(If `join_attributes = TRUE`) Path to joined filtered shapefile.}
#'   \item{polygons_all}{`sf` object of all polygons with computed metrics.}
#'   \item{polygons_filtered}{`sf` object of polygons that passed the filters.}
#' }
#'
#' @note Examples require large external raster files (hosted on Zenodo).
#' Therefore, they are wrapped#' in dontrun{} to avoid errors during R CMD check
#' and to ensure portability.
#'
#' @examples
#' \dontrun{
#' burned_files <- list.files("path/to/shapefiles", pattern = "\\.shp$", full.names = TRUE)
#'
#' results <- calculate_polygon_metrics(
#'   shapefile_paths = burned_files,
#'   output_dir = "outputs/",
#'   area_min_ha = 10,
#'   bbox_h_min = 600,
#'   mnbbx_wd_min = 800,
#'   p_w_ratio_min = 4.0,
#'   h_w_ratio_min = 0.25,
#'   output_format = "geojson",
#'   filter_logic = "AND",
#'   dissolve = TRUE,
#'   join_attributes = TRUE,
#'   min_input_area_ha = 5,
#'   columns_to_keep = c("CLC_CODE", "ECOR_REGION")
#' )
#' }
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
  "area_orig",
  "sub_uid",
  "bbox_wx",
  "bbox_hy",
  "perim_m",
  "mnbbx_ln",
  "mnbbx_el",
  "inter_area",
  "logic_tag"
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
    min_input_area_ha = 10,
    overlay_polygons_path = NULL,
    columns_to_keep = NULL
) {

  output_format <- match.arg(output_format, choices = c("shp", "geojson"))

  if (!is.null(filter_logic)) {
    filter_logic <- match.arg(filter_logic, choices = c("AND", "OR"))
  }

  join_attributes <- isTRUE(join_attributes)

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

    # Lista de patrones de metricas a recalcular
    metric_patterns <- c("area_", "bbox_", "perim_", "p_w_", "h_w_", "mnbbx_", "uid", "ornt_bbox")

    # Detectar columnas a eliminar
    cols_to_remove <- names(polys)[sapply(metric_patterns, function(pat) {
      any(grepl(pat, names(polys)))
    })]

    cols_to_remove <- unique(unlist(lapply(metric_patterns, function(pat) {
      grep(pat, names(polys), value = TRUE)
    })))

    if (length(cols_to_remove) > 0) {
      message("Removing existing metric columns (lax match): ", paste(cols_to_remove, collapse = ", "))
      polys <- polys[, !(names(polys) %in% cols_to_remove)]
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


    if (!all(sf::st_is_valid(polys))) {
      polys <- sf::st_make_valid(polys)
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

    # FILTERING LOGIC (revisado)
    if (!is.null(filter_logic)) {
      safe_ge <- function(x, threshold) {
        !is.na(x) & x >= threshold
      }

      masks <- list()
      filter_labels <- c()

      # Filtros activos
      if (!is.null(area_min_ha)) {
        m_area <- safe_ge(polys$area_ha, area_min_ha)
        masks <- c(masks, list(m_area))
        filter_labels <- c(filter_labels, "area")
        message("- Area ? ", area_min_ha, ": ", sum(m_area), " / ", length(m_area))
      }

      if (!is.null(bbox_h_min)) {
        m_bbox_h <- safe_ge(polys$bbox_hy, bbox_h_min)
        masks <- c(masks, list(m_bbox_h))
        filter_labels <- c(filter_labels, "bbox_h")
        message("- BBox Height ? ", bbox_h_min, ": ", sum(m_bbox_h))
      }

      if (!is.null(mnbbx_wd_min)) {
        m_mnbbx_wd <- safe_ge(polys$mnbbx_wd, mnbbx_wd_min + 1e-6)
        masks <- c(masks, list(m_mnbbx_wd))
        filter_labels <- c(filter_labels, "mnbbx_wd")
        message("- MinBBox Width ? ", mnbbx_wd_min, ": ", sum(m_mnbbx_wd))
      }

      if (!is.null(p_w_ratio_min)) {
        m_p_w <- safe_ge(polys$p_w_ratio, p_w_ratio_min + 1e-6)
        masks <- c(masks, list(m_p_w))
        filter_labels <- c(filter_labels, "p_w_ratio")
        message("- Perim/Width Ratio ? ", p_w_ratio_min, ": ", sum(m_p_w))
      }

      if (!is.null(h_w_ratio_min)) {
        m_h_w <- safe_ge(polys$h_w_ratio, h_w_ratio_min + 1e-6)
        masks <- c(masks, list(m_h_w))
        filter_labels <- c(filter_labels, "h_w_ratio")
        message("- Height/Width Ratio ? ", h_w_ratio_min, ": ", sum(m_h_w))
      }

      if (length(masks) == 0) {
        warning("No filters specified, skipping filtering.")
        combined_mask <- rep(TRUE, nrow(polys))
      } else {
        combined_mask <- if (filter_logic == "AND") {
          Reduce(`&`, masks)
        } else {
          Reduce(`|`, masks)
        }
      }

      filtered <- polys[combined_mask, ]

      suffix <- if (length(filter_labels) == 5) "all" else paste(filter_labels, collapse = "_")
      suffix_tag <- paste0(suffix, "_", tolower(filter_logic))

      base_name <- paste0(shp_name, "_metrics_filt_", suffix_tag)
      filtered_path <- file.path(out_dir, paste0(base_name, ext))

      # Cambiar a GPKG si es muy largo
      if (nchar(filtered_path) > 240 && output_format == "shp") {
        message("Path too long for SHP format. Switching to GPKG.")
        ext <- ".gpkg"
        filtered_path <- file.path(out_dir, paste0(base_name, ext))
      }

      # Eliminar archivos previos si existen
      if (ext == ".shp") {
        shp_base <- tools::file_path_sans_ext(filtered_path)
        shp_exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        for (ext_i in shp_exts) {
          f <- paste0(shp_base, ext_i)
          tryCatch({ if (file.exists(f)) file.remove(f) }, error = function(e) {})
        }
      } else {
        if (file.exists(filtered_path)) file.remove(filtered_path)
      }

      sf::st_write(filtered, filtered_path, append = FALSE, quiet = TRUE)
      message("saved filtered metrics file (", filter_logic, "): ", filtered_path)

    } else {
      message("No filtering applied because 'filter_logic' is NULL.")
      filtered <- NULL
      filtered_path <- NULL
    }


 # Perform spatial join to propagate metrics from dissolved polygons back to original sub-polygons

    if (!is.null(join_attributes) && isTRUE(join_attributes)) {
      message("Joining attributes from original shapefile: ", shapefile_path)

      overlay_polygons <- if (is.null(overlay_polygons_path)) {
        sf::st_read(shapefile_path, quiet = TRUE)
      } else {
        sf::st_read(overlay_polygons_path, quiet = TRUE)
      }

      overlay_polygons <- sf::st_make_valid(overlay_polygons)
      overlay_polygons <- sf::st_transform(overlay_polygons, sf::st_crs(polys))
      overlay_polygons$area_orig <- round(as.numeric(sf::st_area(overlay_polygons)) / 10000, 3)

      if (!is.null(min_input_area_ha)) {
        initial_n <- nrow(overlay_polygons)
        overlay_polygons <- overlay_polygons %>% dplyr::filter(area_orig >= min_input_area_ha)
        filtered_n <- nrow(overlay_polygons)
        message("Filtered input polygons: ", initial_n - filtered_n,
                " polygons removed (area < ", min_input_area_ha, " ha).")
      }

      overlay_polygons$sub_uid <- seq_len(nrow(overlay_polygons))

      if (!is.null(columns_to_keep)) {
        valid_cols <- intersect(columns_to_keep, names(overlay_polygons))
        overlay_polygons <- overlay_polygons %>%
          dplyr::select(all_of(c("sub_uid", valid_cols)))
      } else {
        overlay_polygons <- overlay_polygons %>% dplyr::select(sub_uid)
      }

      # ==== JOINED_ALL ====
      message("Computing intersection-based join for full metrics (preserving original geometry)...")
      joined_all <- tryCatch({
        inter_geom <- sf::st_intersection(
          overlay_polygons,
          polys %>% dplyr::select(burned_id, area_ha, bbox_wx, bbox_hy, perim_m,
                                  p_w_ratio, h_w_ratio, mnbbx_wd, mnbbx_ln, mnbbx_el)
        )
        inter_geom$inter_area <- sf::st_area(inter_geom)

        best_match <- inter_geom %>%
          dplyr::group_by(sub_uid) %>%
          dplyr::slice_max(order_by = inter_area, n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(inter_area = as.numeric(inter_area)) %>%
          dplyr::select(-inter_area)

        best_match
      }, error = function(e) {
        warning("Intersection-based join failed for full polygons: ", e$message)
        return(NULL)
      })

      # ==== JOINED_FILTERED ====
      message("Computing intersection-based join for filtered metrics...")
      joined_filtered <- tryCatch({
        inter_geom_filt <- sf::st_intersection(
          overlay_polygons,
          filtered %>% dplyr::select(burned_id, area_ha, bbox_wx, bbox_hy, perim_m,
                                     p_w_ratio, h_w_ratio, mnbbx_wd, mnbbx_ln, mnbbx_el)
        )
        inter_geom_filt$inter_area <- sf::st_area(inter_geom_filt)

        best_match_filt <- inter_geom_filt %>%
          dplyr::group_by(sub_uid) %>%
          dplyr::slice_max(order_by = inter_area, n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(inter_area = as.numeric(inter_area)) %>%
          dplyr::select(-inter_area)

        best_match_filt
      }, error = function(e) {
        warning("Intersection-based join failed for filtered polygons: ", e$message)
        return(NULL)
      })


      # ==== GUARDAR ====
      if (!exists("result_list")) result_list <- list()
      if (is.null(result_list[[shapefile_path]])) result_list[[shapefile_path]] <- list()

      ext <- if (output_format[1] == "geojson") ".geojson" else ".shp"

      # Guardar joined_all
      if (!is.null(joined_all)) {
        joined_all <- sf::st_make_valid(joined_all)
        joined_all <- joined_all[!sf::st_is_empty(joined_all), ]
        joined_all <- joined_all[sf::st_geometry_type(joined_all) %in% c("POLYGON", "MULTIPOLYGON"), ]

        # Abreviar nombres si es SHP
        if (ext == ".shp") names(joined_all) <- abbreviate(names(joined_all), minlength = 8, strict = TRUE)

        joined_all_path <- sub("\\.shp$", paste0("_metricsJ", ext), shapefile_path)

        # Convertir a .gpkg si el nombre es muy largo
        if (nchar(joined_all_path) > 240 && ext == ".shp") {
          message("Path too long for SHP (joined_all). Switching to GPKG.")
          ext <- ".gpkg"
          joined_all_path <- sub("\\.shp$", paste0("_metricsJ", ext), shapefile_path)
        }

        message("Saving joined full metrics: ", joined_all_path)
        tryCatch({
          suppressWarnings(
            suppressMessages(
              sf::st_write(joined_all, joined_all_path, delete_dsn = TRUE, quiet = TRUE)
            )
          )
          result_list[[shapefile_path]]$joined_metrics <- joined_all_path
        }, error = function(e) {
          message("error while writing joined_all: ", e$message)
        })

      }


      # Guardar joined_filtered
      if (!is.null(joined_filtered)) {
        joined_filtered <- sf::st_make_valid(joined_filtered)
        joined_filtered <- joined_filtered[!sf::st_is_empty(joined_filtered), ]
        joined_filtered <- joined_filtered[sf::st_geometry_type(joined_filtered) %in% c("POLYGON", "MULTIPOLYGON"), ]

        if (nrow(joined_filtered) == 0) {
          warning("No valid geometries in joined_filtered; skipping save.")
        } else {
          if (ext == ".shp") names(joined_filtered) <- abbreviate(names(joined_filtered), minlength = 8, strict = TRUE)

          # Usar mismo sufijo de filtrado que en filtered_path
          suffix_tag <- if (length(filter_labels) == 5) {
            "all"
          } else {
            paste(filter_labels, collapse = "_")
          }
          suffix_tag <- paste0(suffix_tag, "_", tolower(filter_logic))

          joined_filtered_name <- paste0("_metrics_filtJ_", suffix_tag, ext)
          joined_filtered_path <- sub("\\.shp$", joined_filtered_name, shapefile_path)

          if (nchar(joined_filtered_path) > 240 && ext == ".shp") {
            message("Path too long for SHP (joined_filtered). Switching to GPKG.")
            ext <- ".gpkg"
            joined_filtered_path <- sub("\\.shp$", paste0("_metrics_filtJ", logic_tag, ext), shapefile_path)
          }

          message("Saving joined filtered metrics: ", joined_filtered_path)
          tryCatch({
            suppressWarnings(
              suppressMessages(
                sf::st_write(joined_filtered, joined_filtered_path, delete_dsn = TRUE, quiet = TRUE)
              )
            )
            result_list[[shapefile_path]]$joined_filtered <- joined_filtered_path
          }, error = function(e) {
            message("error while writing joined_filtered: ", e$message)
          })
        }
      }

    }

    # Resultado actualizado
    result_list <- list(
      metrics = metrics_path,
      filtered = filtered_path,
      joined_metrics = if (join_attributes && exists("joined_all_path")) joined_all_path else NULL,
      joined_filtered = if (join_attributes && exists("joined_filtered_path")) joined_filtered_path else NULL,
      polygons_all = polys,
      polygons_filtered = filtered
    )


    results[[shapefile_path]] <- result_list

  }

  return(results)
}

