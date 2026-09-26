#' Canonical empty-geometry-safe polygon sanitiser (probabilistic refinement phase).
#'
#' @description
#' The ONE canonical polygon-sanitisation routine used across the supervised
#' phase. It is empty-geometry-safe BY CONSTRUCTION and version-independent --
#' it never relies on `sf::st_collection_extract()` swallowing or zeroing a
#' layer that merely contains empty geometries (the sf 1.0-20 defect that the
#' real-2017 smoke surfaced: an indiscriminate `st_collection_extract()` on a
#' whole layer that contained 4 empty geometries among 1278 features returned 0
#' rows, zeroing the layer and aborting `audit_deterministic_pools()` with
#' "internal_sf is empty").
#'
#' Steps, IN ORDER:
#' \enumerate{
#'   \item drop PRE-EXISTING empty geometries;
#'   \item repair via [sf::st_make_valid()];
#'   \item drop geometries that BECAME empty after make_valid;
#'   \item extract polygonal components ONLY from rows that are a REAL
#'     `GEOMETRYCOLLECTION` -- `POLYGON` / `MULTIPOLYGON` rows are kept INTACT;
#'     a `GEOMETRYCOLLECTION` with NO polygonal component is DISCARDED and
#'     COUNTED (never vanishes silently);
#'   \item cast the kept polygonal rows to `MULTIPOLYGON` for a homogeneous
#'     output type, preserving source attributes / identifiers and the original
#'     row order of the kept rows.
#' }
#' If AFTER sanitisation NO usable geometry remains, an INFORMATIVE error is
#' raised (not a bare "empty").
#'
#' @param x An `sf` object.
#' @param geom_name Target active geometry column name (default `"geometry"`).
#' @param id_col Optional name of a stable identifier column; when present the
#'   IDs of discarded rows are recorded in the audit.
#' @param normalize_geom_col Optional function `(x, geom_name) -> x` used to
#'   normalise the active geometry column (the safety kit passes its own); when
#'   `NULL` a built-in fallback is used.
#' @param error_on_empty Logical. When `TRUE` (default) raise an informative
#'   error if no usable geometry remains; when `FALSE` return the (empty) layer
#'   plus the audit so a DRY-RUN caller can inspect counts without aborting.
#'
#' @return A `list` with:
#'   \describe{
#'     \item{`sf`}{the sanitised `sf` (MULTIPOLYGON, deterministic order);}
#'     \item{`audit`}{a one-row-style `list` with `n_input`,
#'       `n_empty_before`, `n_after_drop_preexisting_empty`,
#'       `n_invalid_before_make_valid`, `n_empty_after_make_valid`,
#'       `n_geometrycollection`, `n_without_polygon_component`, `n_output`, and
#'       `discarded_ids` (character; only when `id_col` resolves to a column).}
#'   }
#'
#' @keywords internal
#' @noRd
.of_sanitize_supervised_polygons <- function(x,
                                             geom_name = "geometry",
                                             id_col = NULL,
                                             normalize_geom_col = NULL,
                                             error_on_empty = TRUE) {
  if (!inherits(x, "sf")) {
    stop(".of_sanitize_supervised_polygons(): 'x' must be an sf object.",
         call. = FALSE)
  }

  norm <- normalize_geom_col
  if (is.null(norm)) {
    norm <- function(x, geom_name = "geometry") {
      g <- attr(x, "sf_column")
      if (!is.null(g) && g != geom_name) {
        # Rename the active geometry column to `geom_name`. Detach the geometry
        # first (so the data frame no longer carries an sfc under the old name),
        # then re-attach it under the new name -- this keeps the sf_column
        # attribute consistent (a bare names<- would leave it dangling).
        geom <- sf::st_geometry(x)
        df <- sf::st_drop_geometry(x)
        df[[geom_name]] <- geom
        x <- sf::st_as_sf(df, sf_column_name = geom_name)
      }
      x
    }
  }

  # Auto-detect a stable id column (default convention: "fire_uid") so the
  # audit's discarded_ids are populated even when the caller does not name one.
  if (is.null(id_col) && "fire_uid" %in% names(x)) id_col <- "fire_uid"

  resolve_ids <- function(x, idx) {
    if (is.null(id_col) || !(id_col %in% names(x))) return(character(0))
    as.character(sf::st_drop_geometry(x)[[id_col]][idx])
  }

  x <- norm(x, geom_name)

  n_input <- nrow(x)
  audit <- list(
    n_input                          = n_input,
    n_empty_before                   = 0L,
    n_after_drop_preexisting_empty   = n_input,
    n_invalid_before_make_valid      = 0L,
    n_empty_after_make_valid         = 0L,
    n_geometrycollection             = 0L,
    n_without_polygon_component      = 0L,
    n_output                         = 0L,
    discarded_ids                    = character(0)
  )

  # --- (1) drop PRE-EXISTING empties --------------------------------------
  empty_before <- sf::st_is_empty(sf::st_geometry(x))
  empty_before[is.na(empty_before)] <- TRUE
  audit$n_empty_before <- sum(empty_before)
  if (audit$n_empty_before > 0L) {
    audit$discarded_ids <- c(audit$discarded_ids,
                             resolve_ids(x, which(empty_before)))
    x <- x[!empty_before, , drop = FALSE]
  }
  audit$n_after_drop_preexisting_empty <- nrow(x)

  if (nrow(x) > 0L) {
    # --- (2) repair via make_valid ----------------------------------------
    invalid_before <- !sf::st_is_valid(x)
    invalid_before[is.na(invalid_before)] <- TRUE
    audit$n_invalid_before_make_valid <- sum(invalid_before)
    x <- sf::st_make_valid(x)

    # --- (3) drop geometries that BECAME empty after make_valid -----------
    empty_after <- sf::st_is_empty(sf::st_geometry(x))
    empty_after[is.na(empty_after)] <- TRUE
    audit$n_empty_after_make_valid <- sum(empty_after)
    if (audit$n_empty_after_make_valid > 0L) {
      audit$discarded_ids <- c(audit$discarded_ids,
                               resolve_ids(x, which(empty_after)))
      x <- x[!empty_after, , drop = FALSE]
    }
  }

  if (nrow(x) > 0L) {
    # --- (4) extract polygonal components ONLY from real GEOMETRYCOLLECTION
    gtype <- as.character(sf::st_geometry_type(x, by_geometry = TRUE))
    is_gc <- gtype == "GEOMETRYCOLLECTION"
    audit$n_geometrycollection <- sum(is_gc)

    if (any(is_gc)) {
      gc_idx <- which(is_gc)
      crs0 <- sf::st_crs(x)
      g_all <- sf::st_geometry(x)
      no_poly <- logical(length(gc_idx))
      # Extract polygonal parts ROW-WISE: st_collection_extract on a whole sfc
      # may flatten a GC into MULTIPLE geometries (breaking the 1:1 row mapping),
      # so we extract per row and re-union the polygonal parts into ONE
      # (multi)polygon for that row. A GC with NO polygonal part yields an EMPTY
      # geometry, which we detect + discard + count.
      for (k in seq_along(gc_idx)) {
        one <- sf::st_sfc(g_all[[gc_idx[k]]], crs = crs0)
        poly_parts <- suppressWarnings(
          sf::st_collection_extract(one, "POLYGON", warn = FALSE))
        if (length(poly_parts) == 0L ||
            all(sf::st_is_empty(poly_parts))) {
          no_poly[k] <- TRUE
        } else {
          merged <- suppressWarnings(sf::st_union(poly_parts))
          g_all[[gc_idx[k]]] <- merged[[1L]]
        }
      }
      audit$n_without_polygon_component <- sum(no_poly)
      if (audit$n_without_polygon_component > 0L) {
        audit$discarded_ids <- c(
          audit$discarded_ids, resolve_ids(x, gc_idx[no_poly]))
      }
      # Write the extracted polygonal geometries back into their ORIGINAL rows
      # (preserving deterministic order); POLYGON / MULTIPOLYGON rows untouched.
      sf::st_geometry(x) <- g_all
      # Discard the GC rows that had no polygonal component.
      if (audit$n_without_polygon_component > 0L) {
        drop_rows <- gc_idx[no_poly]
        keep <- setdiff(seq_len(nrow(x)), drop_rows)
        x <- x[keep, , drop = FALSE]
      }
    }
  }

  if (nrow(x) > 0L) {
    # Guard: st_make_valid can collapse a degenerate polygon to a NON-empty,
    # NON-polygonal geometry (a POINT / LINESTRING) that survives the empty
    # check yet has no polygonal component. Keep ONLY POLYGON / MULTIPOLYGON
    # rows so the homogenising cast below is well-defined; drop + count any
    # stray non-polygonal survivor (same semantic class as a GC without a
    # polygonal part: a row with no usable polygon).
    ftype <- as.character(sf::st_geometry_type(x, by_geometry = TRUE))
    keep_poly <- ftype %in% c("POLYGON", "MULTIPOLYGON")
    if (!all(keep_poly)) {
      n_strays <- sum(!keep_poly)
      audit$n_without_polygon_component <-
        audit$n_without_polygon_component + n_strays
      audit$discarded_ids <- c(audit$discarded_ids,
                               resolve_ids(x, which(!keep_poly)))
      x <- x[keep_poly, , drop = FALSE]
    }
  }

  if (nrow(x) > 0L) {
    # Homogenise to MULTIPOLYGON; keep attributes + deterministic row order.
    x <- suppressWarnings(sf::st_cast(x, "MULTIPOLYGON", warn = FALSE))
    x <- norm(x, geom_name)
  }

  audit$n_output <- nrow(x)

  if (audit$n_output == 0L && isTRUE(error_on_empty)) {
    stop(sprintf(
      paste0(".of_sanitize_supervised_polygons(): no usable polygon geometry ",
             "remains after sanitisation. n_input=%d, n_empty_before=%d, ",
             "n_empty_after_make_valid=%d, n_geometrycollection=%d, ",
             "n_without_polygon_component=%d. The layer was NOT silently ",
             "zeroed by st_collection_extract -- it genuinely contains no ",
             "polygonal geometry after dropping empties and repairing."),
      audit$n_input, audit$n_empty_before, audit$n_empty_after_make_valid,
      audit$n_geometrycollection, audit$n_without_polygon_component),
      call. = FALSE)
  }

  list(sf = x, audit = audit)
}

make_sf_safety_kit <- function(dirs = NULL,
                               result_dir = NULL,
                               target_year = NULL) {
  
  # ============================== GENERIC HELPERS ==========================
  msg <- function(...) {
    cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf(...), "\n")
  }
  
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }
  
  announce_start <- function(title = "MASTER pipeline") {
    cat(title, "started at: ", format(Sys.time()), "\n")
    invisible(NULL)
  }
  
  announce_end <- function(title = "MASTER pipeline") {
    cat(title, "finished at: ", format(Sys.time()), "\n")
    invisible(NULL)
  }
  
  get_corine_year <- function(y) {
    if (y >= 1984 && y <= 1999) {
      "1990"
    } else if (y <= 2005) {
      "2000"
    } else if (y <= 2011) {
      "2006"
    } else if (y <= 2017) {
      "2012"
    } else {
      "2018"
    }
  }
  
  check_exists <- function(x, label = deparse(substitute(x))) {
    if (!file.exists(x)) {
      stop(sprintf("Does not exist: %s\n%s", label, x), call. = FALSE)
    }
  }
  
  stop_if_missing <- function(paths) {
    miss <- paths[!file.exists(paths)]
    if (length(miss)) {
      stop("Missing file(s):\n- ", paste(miss, collapse = "\n- "), call. = FALSE)
    }
  }
  
  safe_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }
  
  # ========================= SF / GEOMETRY HELPERS =========================
  normalize_geom_col <- function(x, geom_name = "geometry") {
    stopifnot(inherits(x, "sf"))
    g <- attr(x, "sf_column")
    
    if (!is.null(g) && g != geom_name && g %in% names(x)) {
      x <- dplyr::rename(x, !!geom_name := dplyr::all_of(g))
      x <- sf::st_as_sf(x, sf_column_name = geom_name)
    } else {
      if (geom_name %in% names(x)) sf::st_geometry(x) <- geom_name
    }
    x
  }
  
  sanitize_polygons <- function(x, id_col = "fire_uid") {
    stopifnot(inherits(x, "sf"))
    res <- .of_sanitize_supervised_polygons(
      x,
      geom_name          = "geometry",
      id_col             = if (!is.null(id_col) && id_col %in% names(x)) id_col else NULL,
      normalize_geom_col = normalize_geom_col,
      error_on_empty     = TRUE)
    a <- res$audit
    # Only log when sanitisation actually changed something (avoid noise on the
    # common all-clean case).
    if (a$n_empty_before > 0L || a$n_empty_after_make_valid > 0L ||
        a$n_without_polygon_component > 0L ||
        a$n_invalid_before_make_valid > 0L) {
      msg(paste0("[sanitize_polygons] n_input=%d -> n_output=%d | ",
                 "empty_before=%d | invalid_before=%d | empty_after_valid=%d | ",
                 "geomcollection=%d | gc_without_polygon=%d"),
          a$n_input, a$n_output, a$n_empty_before, a$n_invalid_before_make_valid,
          a$n_empty_after_make_valid, a$n_geometrycollection,
          a$n_without_polygon_component)
    }
    res$sf
  }
  
  ensure_area_ha <- function(x, area_col = "area_ha") {
    stopifnot(inherits(x, "sf"))
    x[[area_col]] <- as.numeric(sf::st_area(x) / 1e4)
    x
  }
  
  drop_empty <- function(x, tag, dump_dir = NULL) {
    stopifnot(inherits(x, "sf"))
    x <- normalize_geom_col(x, "geometry")
    idx <- which(sf::st_is_empty(sf::st_geometry(x)))
    
    if (length(idx) > 0) {
      msg("[drop_empty] %s: removing %d EMPTY geometries out of %d", tag, length(idx), nrow(x))
      
      if (!is.null(dump_dir)) {
        dir.create(dump_dir, recursive = TRUE, showWarnings = FALSE)
        dump_path <- file.path(dump_dir, paste0(tag, "_EMPTY.gpkg"))
        sf::st_write(x[idx, ], dump_path, layer = "empties", delete_layer = TRUE, quiet = TRUE)
      }
      
      x <- x[-idx, ]
    }
    x
  }
  
  check_sf <- function(x, tag) {
    stopifnot(inherits(x, "sf"))
    x <- normalize_geom_col(x, "geometry")
    ne <- sum(sf::st_is_empty(sf::st_geometry(x)))
    nv <- sum(!sf::st_is_valid(x), na.rm = TRUE)
    
    msg("[check_sf] %s | n=%d | empty=%d | invalid=%d | geom_col=%s",
        tag, nrow(x), ne, nv, attr(x, "sf_column"))
    invisible(x)
  }
  
  check_sf_if_nonempty <- function(x, tag) {
    stopifnot(inherits(x, "sf"))
    if (nrow(x) > 0) check_sf(x, tag)
    invisible(x)
  }
  
  empty_sf <- function(crs_obj = sf::NA_crs_) {
    sf::st_sf(geometry = sf::st_sfc(crs = crs_obj))[0, ]
  }
  
  to_crs_safe <- function(x, crs_target) {
    stopifnot(inherits(x, "sf"))
    if (nrow(x) == 0) return(x)
    if (is.na(sf::st_crs(x))) stop("Object has NA CRS; cannot transform safely.")
    
    if (!is.null(crs_target) && !is.na(crs_target) && sf::st_crs(x) != crs_target) {
      x <- sf::st_transform(x, crs_target)
    }
    x
  }
  
  safe_write_gpkg <- function(x, gpkg_path, layer, tag) {
    dump_dir <- if (!is.null(dirs) && "99_LOGS_EMPTY" %in% names(dirs)) dirs$`99_LOGS_EMPTY` else NULL
    
    x <- sanitize_polygons(x)
    x <- drop_empty(x, tag = tag, dump_dir = dump_dir)
    
    g_active <- attr(x, "sf_column")
    sfc_cols <- names(x)[vapply(x, inherits, logical(1), "sfc")]
    drop_cols <- setdiff(sfc_cols, g_active)
    if (length(drop_cols) > 0) x <- dplyr::select(x, -dplyr::all_of(drop_cols))
    
    sf::st_write(x, gpkg_path, layer = layer, delete_layer = TRUE, quiet = TRUE)
    invisible(gpkg_path)
  }
  
  safe_read_gpkg <- function(gpkg_path, layer, tag) {
    dump_dir <- if (!is.null(dirs) && "99_LOGS_EMPTY" %in% names(dirs)) dirs$`99_LOGS_EMPTY` else NULL
    
    x <- sf::st_read(gpkg_path, layer = layer, quiet = TRUE)
    x <- sanitize_polygons(x)
    x <- drop_empty(x, tag = paste0(tag, "_read"), dump_dir = dump_dir)
    check_sf(x, paste0(tag, " (after read)"))
    x
  }
  
  read_layer_if_exists <- function(gpkg_path, layer, crs_target = NULL, dump_tag = layer) {
    dump_dir <- if (!is.null(dirs) && "99_LOGS_EMPTY" %in% names(dirs)) dirs$`99_LOGS_EMPTY` else NULL
    
    if (!file.exists(gpkg_path)) {
      return(empty_sf(crs_target %||% sf::NA_crs_))
    }
    
    lyr <- tryCatch(sf::st_layers(gpkg_path)$name, error = function(e) character())
    if (!(layer %in% lyr)) {
      return(empty_sf(crs_target %||% sf::NA_crs_))
    }
    
    x <- sf::st_read(gpkg_path, layer = layer, quiet = TRUE)
    x <- sanitize_polygons(x)
    x <- ensure_area_ha(x)
    x <- drop_empty(x, tag = dump_tag, dump_dir = dump_dir)
    
    if (!is.null(crs_target) && !is.na(crs_target)) {
      x <- to_crs_safe(x, crs_target)
    }
    
    x
  }
  
  .pick_one <- function(x, what = "file") {
    if (length(x) == 0) return(NULL)
    if (length(x) > 1) stop("Found more than one ", what, " candidate.\n", paste(" -", x, collapse = "\n"))
    x[[1]]
  }
  
  .find_in_dir <- function(dir, pattern, what = "file") {
    if (!dir.exists(dir)) return(NULL)
    .pick_one(list.files(dir, pattern = pattern, full.names = TRUE), what = what)
  }
  
  # ============================== RASTER HELPERS ===========================
  read_single_raster <- function(path) {
    r <- terra::rast(path)
    if (terra::nlyr(r) > 1) r <- r[[1]]
    r
  }
  
  align_to_template <- function(x, template, method = "near") {
    if (terra::same.crs(x, template)) {
      terra::resample(x, template, method = method)
    } else {
      terra::project(x, template, method = method)
    }
  }
  
  # ============================== WHITEBOX HELPERS =========================
  init_whitebox_safe <- function(exe_path) {
    check_exists(exe_path, "whitebox_exe")
    whitebox::wbt_init(exe_path = exe_path)
    msg("WhiteboxTools OK: %s", exe_path)
    invisible(TRUE)
  }
  
  # ============================== TIMING HELPERS ===========================
  timelog <- data.frame(
    step = character(),
    start = as.POSIXct(character()),
    end = as.POSIXct(character()),
    elapsed_sec = numeric(),
    stringsAsFactors = FALSE
  )
  
  timing_dir <- NULL
  timing_csv <- NULL
  
  if (!is.null(result_dir) && !is.null(target_year)) {
    timing_dir <- file.path(result_dir, "99_LOGS_TIMING")
    timing_csv <- file.path(timing_dir, sprintf("%d_timing_steps.csv", target_year))
    dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  time_step <- function(step, expr) {
    t0 <- Sys.time()
    out <- eval.parent(substitute(expr))
    t1 <- Sys.time()
    
    timelog <<- rbind(
      timelog,
      data.frame(
        step = step,
        start = t0,
        end = t1,
        elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")),
        stringsAsFactors = FALSE
      )
    )
    out
  }
  
  save_timelog <- function(top_n = 5) {
    if (nrow(timelog) == 0) return(invisible(NULL))
    
    if (is.null(timing_csv)) {
      msg("Timing log exists, but result_dir/target_year were not provided; skipping CSV write.")
      return(invisible(timelog))
    }
    
    total_sec <- sum(timelog$elapsed_sec)
    
    out <- timelog |>
      dplyr::mutate(
        elapsed_min = elapsed_sec / 60,
        pct_total   = if (total_sec > 0) 100 * elapsed_sec / total_sec else NA_real_
      ) |>
      dplyr::arrange(dplyr::desc(elapsed_sec))
    
    total_row <- data.frame(
      step = "TOTAL",
      start = min(timelog$start),
      end = max(timelog$end),
      elapsed_sec = total_sec,
      elapsed_min = total_sec / 60,
      pct_total   = 100,
      stringsAsFactors = FALSE
    )
    
    write.csv(rbind(out, total_row), timing_csv, row.names = FALSE)
    
    msg("Timing summary (top %d):", top_n)
    print(utils::head(out[, c("step", "elapsed_sec", "elapsed_min", "pct_total")], top_n))
    
    invisible(out)
  }
  
  get_timelog <- function() timelog
  
  reset_timelog <- function() {
    timelog <<- timelog[0, ]
    invisible(NULL)
  }
  
  # ============================== RETURN KIT ===============================
  list(
    msg                  = msg,
    `%||%`               = `%||%`,
    announce_start       = announce_start,
    announce_end         = announce_end,
    
    get_corine_year      = get_corine_year,
    check_exists         = check_exists,
    stop_if_missing      = stop_if_missing,
    safe_dir             = safe_dir,
    
    normalize_geom_col   = normalize_geom_col,
    sanitize_polygons    = sanitize_polygons,
    ensure_area_ha       = ensure_area_ha,
    drop_empty           = drop_empty,
    check_sf             = check_sf,
    check_sf_if_nonempty = check_sf_if_nonempty,
    empty_sf             = empty_sf,
    to_crs_safe          = to_crs_safe,
    
    safe_write_gpkg      = safe_write_gpkg,
    safe_read_gpkg       = safe_read_gpkg,
    read_layer_if_exists = read_layer_if_exists,
    .pick_one            = .pick_one,
    .find_in_dir         = .find_in_dir,
    
    read_single_raster   = read_single_raster,
    align_to_template    = align_to_template,
    
    init_whitebox_safe   = init_whitebox_safe,
    
    time_step            = time_step,
    save_timelog         = save_timelog,
    get_timelog          = get_timelog,
    reset_timelog        = reset_timelog,
    timing_dir           = timing_dir,
    timing_csv           = timing_csv
  )
}
