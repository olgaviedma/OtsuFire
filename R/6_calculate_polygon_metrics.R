#' Calculate geometric metrics for polygons, optionally filter, and propagate back to original geometry
#'
#' @description
#' This function computes geometric metrics on a (optionally) **dissolved** representation
#' of polygons (contiguous patches), and then propagates those metrics back to the
#' **original** polygons (or to an optional overlay geometry). It can also create a
#' filtered subset based on user thresholds.
#'
#' The implementation is optimized to avoid expensive polygon intersections by using a
#' fast point-in-polygon join (point-on-surface + within). If the join is ambiguous
#' (a point falls inside multiple dissolved polygons), the largest dissolved polygon
#' by area is selected.
#'
#' By default, outputs store only **clean metric names** (no legacy aliases). Set
#' \code{keep_legacy = TRUE} to also store legacy fields for backward compatibility.
#'
#' @section Clean metrics (saved by default):
#' \describe{
#'   \item{\code{burned_id}}{ID of the dissolved polygon (patch) used to attach metrics.}
#'   \item{\code{area_ha}}{Area (ha) of the dissolved polygon.}
#'   \item{\code{bbox_w_m}, \code{bbox_h_m}}{Axis-aligned bounding box width/height (m).}
#'   \item{\code{perimeter_m}}{Perimeter (m).}
#'   \item{\code{mrr_w_m}, \code{mrr_l_m}}{Minimum rotated rectangle short/long side (m).}
#'   \item{\code{elong_mrr}}{Elongation from MRR: \code{mrr_l_m / mrr_w_m} (rotation-invariant).}
#'   \item{\code{compact_pp}}{Polsby-Popper compactness: \code{4*pi*A/P^2} (0-1).}
#'   \item{\code{shape_index}}{Shape index: \code{P/(2*sqrt(pi*A))} (>= 1).}
#'   \item{\code{rect_fill}}{Rectangular fill: \code{A / (mrr_w_m*mrr_l_m)} (0-1).}
#'   \item{\code{join_ambiguous}}{TRUE if the join was ambiguous (multiple dissolved hits).}
#' }
#'
#' @section Filtering:
#' Filtering is applied on dissolved polygons using the provided thresholds. If
#' \code{filter_logic = "AND"}, all active filters must pass; if \code{"OR"}, any active
#' filter is sufficient. Set \code{filter_logic = NULL} to disable filtering entirely.
#'
#' Note: Some filter thresholds are defined on *legacy* metrics (\code{p_w_ratio_min},
#' \code{h_w_ratio_min}). These are computed internally if needed for filtering, but are
#' not written unless \code{keep_legacy = TRUE}.
#'
#' @section Attribute retention:
#' If \code{columns_to_keep = NULL} (default), all attributes in the original/overlay layer
#' are retained, and metrics are appended. If \code{columns_to_keep} is provided, only
#' those fields are retained (plus \code{sub_uid}).
#'
#' @param shapefile_paths Character vector. Paths to vector files (SHP/GeoJSON/GPKG).
#' @param output_dir Character or NULL. Output directory. Defaults to each input folder.
#'
#' @param input_layer Character or NULL. If \code{shapefile_paths} contains GeoPackages,
#'   this selects the layer to read.
#' @param overlay_layer Character or NULL. If \code{overlay_polygons_path} is a GeoPackage path,
#'   this selects the overlay layer to read.
#'
#' @param area_min_ha Numeric or NULL. Minimum \code{area_ha} threshold for filtering.
#' @param bbox_h_min Numeric or NULL. Minimum bbox height (m) threshold for filtering.
#' @param mnbbx_wd_min Numeric or NULL. Minimum MRR short side (m) threshold for filtering.
#' @param p_w_ratio_min Numeric or NULL. Minimum \code{perimeter / bbox width} (legacy) for filtering.
#' @param h_w_ratio_min Numeric or NULL. Minimum \code{bbox height / bbox width} (legacy) for filtering.
#'
#' @param compute_all Logical. If TRUE, compute all metrics and prerequisites.
#' @param compute_area Logical. Compute \code{area_ha}. Default TRUE.
#' @param compute_bbox Logical. Compute bbox metrics. Default FALSE.
#' @param compute_perim Logical. Compute perimeter. Default FALSE.
#' @param compute_mrr Logical. Compute MRR metrics. Default FALSE.
#' @param compute_compact Logical. Compute compactness metrics. Default FALSE.
#' @param compute_rectfill Logical. Compute rectangular fill. Default FALSE.
#'
#' @param output_format Character. One of \code{"shp"}, \code{"geojson"}, \code{"gpkg"}.
#'   For ArcGIS, \code{"gpkg"} is recommended.
#' @param gpkg_layer Character. Layer name used when writing GeoPackage outputs. Default "polys".
#'
#' @param filter_logic Character or NULL. "AND" or "OR". Use NULL to disable filtering.
#' @param dissolve Logical. If TRUE, dissolve adjacent polygons before computing metrics.
#' @param save_dissolved Logical. If TRUE, also write dissolved outputs.
#'
#' @param keep_legacy Logical. If TRUE, also write legacy alias columns (\code{bbox_wx},
#'   \code{bbox_hy}, \code{perim_m}, \code{mnbbx_wd}, \code{mnbbx_ln}, \code{mnbbx_el},
#'   \code{p_w_ratio}, \code{h_w_ratio}). Default FALSE.
#'
#' @param overlay_polygons_path NULL, character path, or sf. Geometry to receive metrics.
#'   Defaults to the input polygons.
#' @param columns_to_keep Character or NULL. Original/overlay columns to retain before adding metrics.
#'   NULL keeps all columns.
#'
#' @return Named list (one element per input) with output paths and in-memory sf objects.
#'
#' @export
calculate_polygon_metrics <- function(
  shapefile_paths,
  output_dir = NULL,

  input_layer = NULL,
  overlay_layer = NULL,

  # ---- OPTIONAL FILTER THRESHOLDS ----
  area_min_ha   = NULL,
  bbox_h_min    = NULL,
  mnbbx_wd_min  = NULL,
  p_w_ratio_min = NULL,
  h_w_ratio_min = NULL,

  # ---- METRIC COMPUTATION CONTROL ----
  compute_all     = FALSE,
  compute_area    = TRUE,
  compute_bbox    = FALSE,
  compute_perim   = FALSE,
  compute_mrr     = FALSE,
  compute_compact = FALSE,
  compute_rectfill = FALSE,

  # ---- I/O options ----
  output_format  = c("gpkg", "shp", "geojson"),
  gpkg_layer     = "polys",
  filter_logic   = c("AND", "OR"),  # can be NULL to disable filtering
  dissolve       = TRUE,
  save_dissolved = FALSE,

  keep_legacy = FALSE,

  # --- Back-compat: accepted but not used for writing ---
  overlay_polygons_path = NULL,
  columns_to_keep       = NULL
) {

  # ---------------------------
  # Helpers
  # ---------------------------
  msg <- function(...) message(...)

  output_format <- match.arg(output_format, choices = c("gpkg", "shp", "geojson"))
  if (!is.null(filter_logic)) filter_logic <- match.arg(filter_logic, choices = c("AND", "OR"))

  is_longlat_safe <- function(x) {
    crs <- sf::st_crs(x)
    if (is.na(crs)) return(NA)
    isTRUE(sf::st_is_longlat(x))
  }

  normalize_crs <- function(x) {
    crs <- sf::st_crs(x)
    if (!is.na(crs$epsg)) {
      # force canonical CRS object from EPSG
      sf::st_crs(x) <- sf::st_crs(crs$epsg)
    }
    x
  }

  safe_make_valid_polys <- function(x) {
    if (is.null(x)) return(NULL)
    x <- sf::st_as_sf(x)

    if (is.na(sf::st_crs(x))) stop("Input has NA CRS. Please assign a CRS before running.")

    # drop empty
    emp <- sf::st_is_empty(sf::st_geometry(x))
    if (any(emp)) x <- x[!emp, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)

    # make valid
    if (any(!sf::st_is_valid(x))) x <- sf::st_make_valid(x)

    # keep polygonal
    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
    if (nrow(x) == 0) return(NULL)

    # homogenize
    x$geometry <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)
    sf::st_as_sf(x)
  }

  # Drop legacy + clean metric fields to avoid duplication across reruns
  drop_existing_metrics <- function(sfobj) {
    metric_patterns <- c(
      # clean
      "area_ha", "bbox_w_m", "bbox_h_m", "perimeter_m",
      "mrr_w_m", "mrr_l_m", "elong_mrr",
      "compact_pp", "shape_index", "rect_fill",
      "burned_id", "join_ambiguous", "sub_uid",
      # internal temp
      "area_m2", "ornt_bbox",
      # legacy
      "bbox_wx", "bbox_hy", "perim_m",
      "mnbbx_wd", "mnbbx_ln", "mnbbx_el",
      "p_w_ratio", "h_w_ratio"
    )
    cols_to_remove <- unique(unlist(lapply(metric_patterns, function(pat) grep(pat, names(sfobj), value = TRUE, fixed = TRUE))))
    if (length(cols_to_remove) > 0) {
      sfobj <- sfobj[, !(names(sfobj) %in% cols_to_remove), drop = FALSE]
    }
    sfobj
  }

  # Make SHP-safe, unique field names (<=10 chars)
  shp_safe_names <- function(nms) {
    nms2 <- abbreviate(nms, minlength = 8, strict = TRUE)
    nms2 <- substr(nms2, 1, 10)
    nmsu <- make.unique(nms2, sep = "_")
    substr(nmsu, 1, 10)
  }

  shorten_path_if_needed <- function(path, max_chars = 240) {
    if (nchar(path) <= max_chars) return(path)
    dirn <- dirname(path)
    ext  <- tools::file_ext(path)
    base <- tools::file_path_sans_ext(basename(path))
    h <- sum(utf8ToInt(base)) %% 4294967295
    h <- sprintf("%08x", h)
    base_short <- substr(base, 1, 60)
    file.path(dirn, paste0(base_short, "_", h, ".", ext))
  }

  delete_existing_vector <- function(path) {
    if (!file.exists(path)) return(invisible(TRUE))
    ext <- tolower(tools::file_ext(path))
    if (ext == "shp") {
      base <- tools::file_path_sans_ext(path)
      for (e in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
        f <- paste0(base, e)
        if (file.exists(f)) tryCatch(file.remove(f), error = function(e) {})
      }
    } else {
      tryCatch(file.remove(path), error = function(e) {})
    }
    invisible(TRUE)
  }

  write_vector <- function(x, path, output_format, gpkg_layer) {
    x <- normalize_crs(x)
    path <- shorten_path_if_needed(path)
    delete_existing_vector(path)

    if (output_format == "shp") {
      names(x) <- shp_safe_names(names(x))
      suppressWarnings(sf::st_write(x, path, append = FALSE, quiet = TRUE))
    } else if (output_format == "gpkg") {
      suppressWarnings(sf::st_write(x, path, layer = gpkg_layer, driver = "GPKG",
                                    delete_dsn = TRUE, quiet = TRUE))
    } else { # geojson
      suppressWarnings(sf::st_write(x, path, driver = "GeoJSON", append = FALSE, quiet = TRUE))
    }
    path
  }

  # Compute min/max side length of a rotated rectangle polygon
  get_wd_ln_from_mrr <- function(poly) {
    tryCatch({
      ls <- sf::st_cast(poly, "LINESTRING", warn = FALSE)
      coords <- sf::st_coordinates(ls)[, 1:2, drop = FALSE]
      if (nrow(coords) < 2) return(c(wd = NA_real_, ln = NA_real_))
      if (all(coords[1, ] == coords[nrow(coords), ])) coords <- coords[-nrow(coords), , drop = FALSE]
      if (nrow(coords) < 2) return(c(wd = NA_real_, ln = NA_real_))
      prev <- rbind(coords[nrow(coords), , drop = FALSE], coords[-nrow(coords), , drop = FALSE])
      segs <- sqrt(rowSums((coords - prev) ^ 2))
      segs <- segs[is.finite(segs) & segs > 0]
      if (length(segs) == 0) return(c(wd = NA_real_, ln = NA_real_))
      c(wd = min(segs, na.rm = TRUE), ln = max(segs, na.rm = TRUE))
    }, error = function(e) c(wd = NA_real_, ln = NA_real_))
  }

  # ---------------------------
  # Main
  # ---------------------------
  if (length(shapefile_paths) == 0) return(list())
  shapefile_paths <- as.character(shapefile_paths)

  results <- list()

  for (shapefile_path in shapefile_paths) {
    msg("Processing: ", shapefile_path)

    out_dir <- if (is.null(output_dir)) dirname(shapefile_path) else output_dir
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    shp_name <- tools::file_path_sans_ext(basename(shapefile_path))

    # --- Read original polygons once (support GPKG layer) ---
    orig <- tryCatch(
      {
        if (!is.null(input_layer)) {
          sf::st_read(shapefile_path, layer = input_layer, quiet = TRUE)
        } else {
          sf::st_read(shapefile_path, quiet = TRUE)
        }
      },
      error = function(e) {
        warning("Failed to read: ", shapefile_path, " (", e$message, ")")
        NULL
      }
    )
    if (is.null(orig)) next

    orig <- safe_make_valid_polys(orig)
    if (is.null(orig) || nrow(orig) == 0) {
      warning("No valid polygon geometries after validation. Skipping: ", shapefile_path)
      next
    }

    # CRS sanity: metrics assume projected units
    ll <- is_longlat_safe(orig)
    if (isTRUE(ll)) {
      stop(
        "Detected lon/lat CRS (degrees). Metrics require a projected CRS in meters.\n",
        "Reproject your data first (e.g., sf::st_transform) and rerun."
      )
    }

    # Remove any old metric fields from original
    orig_clean <- drop_existing_metrics(orig)

    # ---------------------------
    # DISSOLVE (optional)
    # ---------------------------
    if (isTRUE(dissolve)) {
      msg("Dissolving adjacent polygons...")
      union_geom <- sf::st_union(sf::st_geometry(orig_clean))
      cast_polys <- sf::st_cast(union_geom, "POLYGON", warn = FALSE)

      if (length(cast_polys) == 0) {
        warning("Dissolve result is empty. Skipping: ", shapefile_path)
        next
      }

      polys_d <- sf::st_sf(geometry = sf::st_sfc(cast_polys, crs = sf::st_crs(orig_clean)))
      polys_d$geometry <- sf::st_cast(polys_d$geometry, "MULTIPOLYGON", warn = FALSE)
      polys_d <- sf::st_as_sf(polys_d)
    } else {
      msg("Dissolve not applied.")
      polys_d <- orig_clean
    }

    polys_d <- safe_make_valid_polys(polys_d)
    if (is.null(polys_d) || nrow(polys_d) == 0) {
      warning("No valid dissolved polygons. Skipping: ", shapefile_path)
      next
    }

    # IDs for linking
    polys_d$burned_id <- seq_len(nrow(polys_d))

    # ---------------------------
    # Decide which metrics are needed
    # ---------------------------
    any_filter <- !is.null(area_min_ha) || !is.null(bbox_h_min) || !is.null(mnbbx_wd_min) ||
      !is.null(p_w_ratio_min) || !is.null(h_w_ratio_min)

    # prerequisites driven by requested metrics or filtering
    need_area <- isTRUE(compute_all) || isTRUE(compute_area) || !is.null(area_min_ha) ||
      isTRUE(compute_compact) || isTRUE(compute_rectfill)

    # bbox needed for bbox outputs OR filtering OR legacy filters
    need_bbox <- isTRUE(compute_all) || isTRUE(compute_bbox) || !is.null(bbox_h_min) ||
      !is.null(p_w_ratio_min) || !is.null(h_w_ratio_min)

    # perimeter needed for perimeter output OR compactness OR legacy p_w filter
    need_perim <- isTRUE(compute_all) || isTRUE(compute_perim) || isTRUE(compute_compact) ||
      !is.null(p_w_ratio_min)

    # MRR needed for MRR outputs OR rect_fill OR filter on mnbbx_wd
    need_mrr <- isTRUE(compute_all) || isTRUE(compute_mrr) || isTRUE(compute_rectfill) ||
      !is.null(mnbbx_wd_min)

    # ---------------------------
    # Compute internal primitives on dissolved polygons
    # ---------------------------
    if (need_area) {
      polys_d$area_m2 <- sf::st_area(polys_d)
      polys_d$area_ha <- round(as.numeric(polys_d$area_m2) / 10000, 3)
    }

    if (need_bbox) {
      bb_list <- lapply(sf::st_geometry(polys_d), sf::st_bbox)
      polys_d$bbox_wx <- vapply(bb_list, function(b) unname(b["xmax"] - b["xmin"]), numeric(1))
      polys_d$bbox_hy <- vapply(bb_list, function(b) unname(b["ymax"] - b["ymin"]), numeric(1))
    }

    if (need_perim) {
      polys_d$perim_m <- vapply(
        sf::st_geometry(polys_d),
        function(g) as.numeric(sf::st_length(sf::st_cast(g, "MULTILINESTRING"))),
        numeric(1)
      )
    }

    if (need_mrr) {
      mrr_list <- lapply(sf::st_geometry(polys_d), function(g) {
        tryCatch(
          sf::st_minimum_rotated_rectangle(g),
          error = function(e1) {
            g2 <- tryCatch(sf::st_make_valid(g), error = function(e2) NULL)
            if (is.null(g2)) return(NULL)
            tryCatch(sf::st_minimum_rotated_rectangle(g2), error = function(e3) NULL)
          }
        )
      })

      ok <- vapply(mrr_list, function(x) !is.null(x), logical(1))
      ornt <- vector("list", length(mrr_list))
      for (i in seq_along(mrr_list)) {
        ornt[[i]] <- if (ok[i]) mrr_list[[i]] else sf::st_geometrycollection()
      }
      polys_d$ornt_bbox <- sf::st_sfc(ornt, crs = sf::st_crs(polys_d))

      bbox_df <- do.call(rbind, lapply(polys_d$ornt_bbox, function(x) {
        if (inherits(x, "sfc") || inherits(x, "sf")) x <- sf::st_geometry(x)
        if (length(x) == 0) return(c(wd = NA_real_, ln = NA_real_))
        get_wd_ln_from_mrr(x)
      }))

      polys_d$mnbbx_wd <- bbox_df[, "wd"]
      polys_d$mnbbx_ln <- bbox_df[, "ln"]
      polys_d$mnbbx_el <- polys_d$mnbbx_ln / pmax(polys_d$mnbbx_wd, 1e-6)
    }

    # legacy ratios (only needed for filtering or keep_legacy)
    if (!is.null(p_w_ratio_min) || isTRUE(keep_legacy)) {
      if (!("bbox_wx" %in% names(polys_d))) {
        bb_list <- lapply(sf::st_geometry(polys_d), sf::st_bbox)
        polys_d$bbox_wx <- vapply(bb_list, function(b) unname(b["xmax"] - b["xmin"]), numeric(1))
      }
      if (!("perim_m" %in% names(polys_d))) {
        polys_d$perim_m <- vapply(
          sf::st_geometry(polys_d),
          function(g) as.numeric(sf::st_length(sf::st_cast(g, "MULTILINESTRING"))),
          numeric(1)
        )
      }
      polys_d$p_w_ratio <- polys_d$perim_m / pmax(as.numeric(polys_d$bbox_wx), 1e-6)
    }

    if (!is.null(h_w_ratio_min) || isTRUE(keep_legacy)) {
      if (!all(c("bbox_wx", "bbox_hy") %in% names(polys_d))) {
        bb_list <- lapply(sf::st_geometry(polys_d), sf::st_bbox)
        polys_d$bbox_wx <- vapply(bb_list, function(b) unname(b["xmax"] - b["xmin"]), numeric(1))
        polys_d$bbox_hy <- vapply(bb_list, function(b) unname(b["ymax"] - b["ymin"]), numeric(1))
      }
      polys_d$h_w_ratio <- polys_d$bbox_hy / pmax(as.numeric(polys_d$bbox_wx), 1e-6)
    }

    # compactness metrics
    if (isTRUE(compute_all) || isTRUE(compute_compact)) {
      if (!("area_m2" %in% names(polys_d))) polys_d$area_m2 <- sf::st_area(polys_d)
      if (!("perim_m" %in% names(polys_d))) {
        polys_d$perim_m <- vapply(
          sf::st_geometry(polys_d),
          function(g) as.numeric(sf::st_length(sf::st_cast(g, "MULTILINESTRING"))),
          numeric(1)
        )
      }
      A <- as.numeric(polys_d$area_m2)
      P <- as.numeric(polys_d$perim_m)
      polys_d$compact_pp  <- (4 * pi * A) / pmax(P^2, 1e-12)
      polys_d$shape_index <- P / pmax(2 * sqrt(pi * A), 1e-12)
    }

    # rectangular fill
    if (isTRUE(compute_all) || isTRUE(compute_rectfill)) {
      if (!("area_m2" %in% names(polys_d))) polys_d$area_m2 <- sf::st_area(polys_d)
      if (!all(c("mnbbx_wd", "mnbbx_ln") %in% names(polys_d))) {
        # ensure MRR dims exist
        mrr_list <- lapply(sf::st_geometry(polys_d), function(g) {
          tryCatch(sf::st_minimum_rotated_rectangle(g), error = function(e) NULL)
        })
        ok <- vapply(mrr_list, function(x) !is.null(x), logical(1))
        ornt <- vector("list", length(mrr_list))
        for (i in seq_along(mrr_list)) {
          ornt[[i]] <- if (ok[i]) mrr_list[[i]] else sf::st_geometrycollection()
        }
        polys_d$ornt_bbox <- sf::st_sfc(ornt, crs = sf::st_crs(polys_d))

        bbox_df <- do.call(rbind, lapply(polys_d$ornt_bbox, function(x) {
          if (inherits(x, "sfc") || inherits(x, "sf")) x <- sf::st_geometry(x)
          if (length(x) == 0) return(c(wd = NA_real_, ln = NA_real_))
          get_wd_ln_from_mrr(x)
        }))
        polys_d$mnbbx_wd <- bbox_df[, "wd"]
        polys_d$mnbbx_ln <- bbox_df[, "ln"]
      }
      A <- as.numeric(polys_d$area_m2)
      polys_d$rect_fill <- A / pmax(polys_d$mnbbx_wd * polys_d$mnbbx_ln, 1e-12)
    }

    # ---------------------------
    # Filtering on dissolved polygons (optional)
    # ---------------------------
    filtered_d <- NULL
    filter_labels <- character(0)

    if (!is.null(filter_logic) && isTRUE(any_filter)) {
      safe_ge <- function(x, thr) !is.na(x) & x >= thr
      masks <- list()

      if (!is.null(area_min_ha)) {
        masks <- c(masks, list(safe_ge(polys_d$area_ha, area_min_ha)))
        filter_labels <- c(filter_labels, "area")
      }
      if (!is.null(bbox_h_min)) {
        masks <- c(masks, list(safe_ge(polys_d$bbox_hy, bbox_h_min)))
        filter_labels <- c(filter_labels, "bbox_h")
      }
      if (!is.null(mnbbx_wd_min)) {
        masks <- c(masks, list(safe_ge(polys_d$mnbbx_wd, mnbbx_wd_min + 1e-6)))
        filter_labels <- c(filter_labels, "mrr_w")
      }
      if (!is.null(p_w_ratio_min)) {
        masks <- c(masks, list(safe_ge(polys_d$p_w_ratio, p_w_ratio_min + 1e-6)))
        filter_labels <- c(filter_labels, "perim_over_bboxw")
      }
      if (!is.null(h_w_ratio_min)) {
        masks <- c(masks, list(safe_ge(polys_d$h_w_ratio, h_w_ratio_min + 1e-6)))
        filter_labels <- c(filter_labels, "bbox_aspect")
      }

      combined_mask <- if (filter_logic == "AND") Reduce(`&`, masks) else Reduce(`|`, masks)
      filtered_d <- polys_d[combined_mask, , drop = FALSE]

      msg("Filtering (", filter_logic, "): kept ", nrow(filtered_d), " / ", nrow(polys_d), " dissolved polygons.")
    } else {
      if (!is.null(filter_logic) && !isTRUE(any_filter)) {
        msg("No filter thresholds provided -> filtered output disabled.")
      } else if (is.null(filter_logic)) {
        msg("filter_logic is NULL -> filtering disabled.")
      }
    }

    # ---------------------------
    # Build CLEAN metric table for dissolved polygons (to propagate + to write)
    # ---------------------------
    clean_cols <- c(
      "burned_id",
      "area_ha",
      "bbox_w_m", "bbox_h_m",
      "perimeter_m",
      "mrr_w_m", "mrr_l_m",
      "elong_mrr",
      "compact_pp", "shape_index",
      "rect_fill"
    )

    # map internal names -> clean names
    dissolved_metrics <- polys_d[, c("burned_id"), drop = FALSE]

    if ("area_ha" %in% names(polys_d)) dissolved_metrics$area_ha <- polys_d$area_ha

    if ("bbox_wx" %in% names(polys_d)) dissolved_metrics$bbox_w_m <- polys_d$bbox_wx
    if ("bbox_hy" %in% names(polys_d)) dissolved_metrics$bbox_h_m <- polys_d$bbox_hy

    if ("perim_m" %in% names(polys_d)) dissolved_metrics$perimeter_m <- polys_d$perim_m

    if ("mnbbx_wd" %in% names(polys_d)) dissolved_metrics$mrr_w_m <- polys_d$mnbbx_wd
    if ("mnbbx_ln" %in% names(polys_d)) dissolved_metrics$mrr_l_m <- polys_d$mnbbx_ln
    if (all(c("mnbbx_ln", "mnbbx_wd") %in% names(polys_d))) {
      dissolved_metrics$elong_mrr <- polys_d$mnbbx_ln / pmax(polys_d$mnbbx_wd, 1e-6)
    }

    if ("compact_pp" %in% names(polys_d)) dissolved_metrics$compact_pp <- polys_d$compact_pp
    if ("shape_index" %in% names(polys_d)) dissolved_metrics$shape_index <- polys_d$shape_index
    if ("rect_fill" %in% names(polys_d)) dissolved_metrics$rect_fill <- polys_d$rect_fill

    # attach geometry (for dissolved write if needed)
    dissolved_metrics <- sf::st_as_sf(dissolved_metrics)
    sf::st_geometry(dissolved_metrics) <- sf::st_geometry(polys_d)

    # optional: include legacy fields in dissolved output
    if (isTRUE(keep_legacy)) {
      # add legacy names as stored fields
      if ("bbox_wx" %in% names(polys_d)) dissolved_metrics$bbox_wx <- polys_d$bbox_wx
      if ("bbox_hy" %in% names(polys_d)) dissolved_metrics$bbox_hy <- polys_d$bbox_hy
      if ("perim_m" %in% names(polys_d)) dissolved_metrics$perim_m <- polys_d$perim_m
      if ("p_w_ratio" %in% names(polys_d)) dissolved_metrics$p_w_ratio <- polys_d$p_w_ratio
      if ("h_w_ratio" %in% names(polys_d)) dissolved_metrics$h_w_ratio <- polys_d$h_w_ratio
      if ("mnbbx_wd" %in% names(polys_d)) dissolved_metrics$mnbbx_wd <- polys_d$mnbbx_wd
      if ("mnbbx_ln" %in% names(polys_d)) dissolved_metrics$mnbbx_ln <- polys_d$mnbbx_ln
      if ("mnbbx_el" %in% names(polys_d)) dissolved_metrics$mnbbx_el <- polys_d$mnbbx_el
    }

    # ---------------------------
    # Propagate dissolved metrics back to ORIGINAL / overlay geometry (fast)
    # ---------------------------
    overlay <- NULL
    if (is.null(overlay_polygons_path)) {
      overlay <- orig_clean
    } else if (inherits(overlay_polygons_path, "sf")) {
      overlay <- overlay_polygons_path
    } else {
      overlay <- tryCatch({
        if (!is.null(overlay_layer)) {
          sf::st_read(overlay_polygons_path, layer = overlay_layer, quiet = TRUE)
        } else {
          sf::st_read(overlay_polygons_path, quiet = TRUE)
        }
      }, error = function(e) {
        warning("Failed to read overlay_polygons_path: ", overlay_polygons_path, " (", e$message, ")")
        NULL
      })
    }

    overlay <- safe_make_valid_polys(overlay)
    if (is.null(overlay) || nrow(overlay) == 0) {
      warning("Overlay polygons are empty/invalid. Skipping: ", shapefile_path)
      next
    }

    overlay <- drop_existing_metrics(overlay)

    # Ensure same CRS as dissolved
    if (sf::st_crs(overlay) != sf::st_crs(polys_d)) {
      overlay <- sf::st_transform(overlay, sf::st_crs(polys_d))
    }

    overlay$sub_uid <- seq_len(nrow(overlay))

    # retain all attributes by default; restrict only if user requested
    if (!is.null(columns_to_keep)) {
      keep_cols <- intersect(columns_to_keep, names(overlay))
      keep_cols <- unique(c("sub_uid", keep_cols))
      overlay <- overlay[, keep_cols, drop = FALSE]
    }

    # point-in-polygon join
    current_s2 <- sf::sf_use_s2()
    sf::sf_use_s2(FALSE)

    pts  <- sf::st_point_on_surface(overlay)
    hits <- sf::st_within(pts, polys_d, sparse = TRUE)

    idx <- rep(NA_integer_, length(hits))
    amb <- rep(FALSE, length(hits))

    lens <- lengths(hits)
    if (any(lens == 1)) {
      idx[lens == 1] <- vapply(hits[lens == 1], function(i) i[[1]], integer(1))
    }

    multi <- which(lens > 1)
    if (length(multi) > 0) {
      amb[multi] <- TRUE
      if (!("area_m2" %in% names(polys_d))) polys_d$area_m2 <- sf::st_area(polys_d)
      area_num <- as.numeric(polys_d$area_m2)
      for (k in multi) {
        cand <- hits[[k]]
        if (length(cand) == 0) next
        idx[k] <- cand[which.max(area_num[cand])]
      }
    }

    sf::sf_use_s2(current_s2)

    overlay_all <- overlay
    overlay_all$join_ambiguous <- amb

    # attach clean metrics by burned_id index
    overlay_all$burned_id <- ifelse(is.na(idx), NA_integer_, polys_d$burned_id[idx])

    # join remaining metrics from dissolved_metrics using idx
    # (avoid matching by burned_id string; this is faster)
    # Note: dissolved_metrics is aligned with polys_d row order.
    attach_cols <- setdiff(names(dissolved_metrics), c("geometry", "burned_id"))
    for (nm in attach_cols) {
      overlay_all[[nm]] <- ifelse(is.na(idx), NA, dissolved_metrics[[nm]][idx])
    }

    # optional: store legacy fields too
    if (isTRUE(keep_legacy)) {
      legacy_cols <- intersect(
        c("bbox_wx","bbox_hy","perim_m","p_w_ratio","h_w_ratio","mnbbx_wd","mnbbx_ln","mnbbx_el"),
        names(polys_d)
      )
      for (nm in legacy_cols) {
        overlay_all[[nm]] <- ifelse(is.na(idx), NA, polys_d[[nm]][idx])
      }
    }

    # Filtered overlay (cheap: by dissolved burned_id membership)
    overlay_filt <- NULL
    if (!is.null(filtered_d) && nrow(filtered_d) > 0) {
      keep <- !is.na(overlay_all$burned_id) & overlay_all$burned_id %in% filtered_d$burned_id
      if (any(keep)) overlay_filt <- overlay_all[keep, , drop = FALSE]
    }

    # ---------------------------
    # Write outputs (ORIGINAL geometry)
    # ---------------------------
    ext_main <- if (output_format == "geojson") ".geojson" else if (output_format == "gpkg") ".gpkg" else ".shp"

    out1 <- file.path(out_dir, paste0(shp_name, "_orig_metrics", ext_main))
    out1 <- write_vector(overlay_all, out1, output_format, gpkg_layer)
    msg("Saved ORIGINAL metrics: ", out1)

    out2 <- NULL
    if (!is.null(overlay_filt) && nrow(overlay_filt) > 0) {
      filt_labels <- if (length(filter_labels) > 0) paste(filter_labels, collapse = "_") else "none"
      suffix_tag  <- paste0(filt_labels, "_", tolower(filter_logic))
      out2 <- file.path(out_dir, paste0(shp_name, "_orig_metrics_filt_", suffix_tag, ext_main))
      out2 <- write_vector(overlay_filt, out2, output_format, gpkg_layer)
      msg("Saved ORIGINAL filtered metrics: ", out2)
    } else {
      msg("No ORIGINAL polygons matched filtered dissolved set (or filtering disabled).")
    }

    # ---------------------------
    # Optionally write dissolved outputs
    # ---------------------------
    dissolved_metrics_path  <- NULL
    dissolved_filtered_path <- NULL

    if (isTRUE(save_dissolved)) {
      out_d1 <- file.path(out_dir, paste0(shp_name, "_diss_metrics", ext_main))
      dissolved_metrics_path <- write_vector(dissolved_metrics, out_d1, output_format, gpkg_layer)
      msg("Saved DISSOLVED metrics: ", dissolved_metrics_path)

      if (!is.null(filtered_d) && nrow(filtered_d) > 0) {
        # build dissolved filtered sf using clean table
        keep_d <- polys_d$burned_id %in% filtered_d$burned_id
        dissolved_metrics_filt <- dissolved_metrics[keep_d, , drop = FALSE]

        filt_labels <- if (length(filter_labels) > 0) paste(filter_labels, collapse = "_") else "none"
        suffix_tag  <- paste0(filt_labels, "_", tolower(filter_logic))

        out_d2 <- file.path(out_dir, paste0(shp_name, "_diss_metrics_filt_", suffix_tag, ext_main))
        dissolved_filtered_path <- write_vector(dissolved_metrics_filt, out_d2, output_format, gpkg_layer)
        msg("Saved DISSOLVED filtered metrics: ", dissolved_filtered_path)
      }
    }

    # ---------------------------
    # Return per input
    # ---------------------------
    results[[shapefile_path]] <- list(
      metrics  = out1,
      filtered = out2,

      dissolved_metrics  = dissolved_metrics_path,
      dissolved_filtered = dissolved_filtered_path,

      polygons_all                = overlay_all,
      polygons_filtered           = overlay_filt,
      polygons_dissolved          = dissolved_metrics,
      polygons_dissolved_filtered = if (!is.null(filtered_d) && nrow(filtered_d) > 0) dissolved_metrics[polys_d$burned_id %in% filtered_d$burned_id, , drop = FALSE] else NULL,

      filter_labels = filter_labels
    )
  }

  results
}
