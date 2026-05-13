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
      stop(sprintf("No existe: %s\n%s", label, x), call. = FALSE)
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
  
  sanitize_polygons <- function(x) {
    stopifnot(inherits(x, "sf"))
    x <- normalize_geom_col(x, "geometry")
    x <- sf::st_make_valid(x)
    x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON", warn = FALSE))
    x <- suppressWarnings(sf::st_cast(x, "MULTIPOLYGON", warn = FALSE))
    normalize_geom_col(x, "geometry")
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
    if (length(x) > 1) stop("Encontre mas de un ", what, " candidato.\n", paste(" -", x, collapse = "\n"))
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
