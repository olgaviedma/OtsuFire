# Internal engine dispatcher. NOT exported.
#
# Implements the full validated deterministic chain in v0.2 contract shape:
#
#   (detect) process_otsu_rasters_grow
#         -> segmentation_refinement
#         -> merge_aoi_shapefiles
#   (score ) scoring_burned_area_stage2
#         -> [build_rbr_keep_pool_from_registry]     (optional, registry-based)
#         -> score_rbr_keep_classes
#         -> build_internal_decisions
#         -> write_canonical_decision_output
#
# Engine source is loaded from config$engine_root into an isolated env.
# Nothing here shadows the public API or reintroduces the 0.1.x bridge.

# Block 4A: deterministic engine code is migrated into the package as
# internal helpers (see R/internal-det-*.R). The dispatcher no longer
# sources external files; engine functions live in the package namespace.
# .engine_env_cache and .of_load_engine() are retained as no-ops for any
# legacy call sites, but require_engine_root is now optional.

# Cached engine environments — kept as a noop placeholder for backward
# compatibility with any code that may still call .of_load_engine(); the
# real engine functions are now package-internal.
.engine_env_cache <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.of_engine_files <- function() {
  # Self-contained: no external file dependencies remain. Returned for
  # backward compatibility / introspection only.
  character(0)
}

#' @keywords internal
#' @noRd
.of_require_engine_root <- function(config) {
  if (!inherits(config, "otsufire_burned_mapping_config")) {
    stop("Internal error: dispatcher expects otsufire_burned_mapping_config.",
         call. = FALSE)
  }
  # No external dependency: engine_root is optional.
  config$engine_root %||% ""
}

#' @keywords internal
#' @noRd
.of_load_engine <- function(engine_root) {
  # Returns an env that resolves engine functions through the package
  # namespace. This indirection lets the rest of the dispatcher keep using
  # `engine$process_otsu_rasters_grow` without changes.
  ns <- topenv(parent.frame())  # falls back below if topenv is empty
  if (!exists("process_otsu_rasters_grow", envir = ns, inherits = TRUE)) {
    ns <- asNamespace("OtsuFire")
  }
  env <- new.env(parent = ns)
  # Bind canonical names that the dispatcher uses to the package internals.
  assign("process_otsu_rasters_grow",
         get("process_otsu_rasters_grow", envir = ns), envir = env)
  assign("segmentation_refinement",
         get("segmentation_refinement", envir = ns), envir = env)
  assign("merge_aoi_shapefiles",
         get("merge_aoi_shapefiles", envir = ns), envir = env)
  assign("scoring_burned_area_stage2",
         get("scoring_burned_area_stage2", envir = ns), envir = env)
  assign("score_rbr_keep_classes",
         get("score_rbr_keep_classes", envir = ns), envir = env)
  if (exists("build_rbr_keep_pool_from_registry", envir = ns, inherits = TRUE)) {
    assign("build_rbr_keep_pool_from_registry",
           get("build_rbr_keep_pool_from_registry", envir = ns), envir = env)
  }
  if (exists("build_internal_decisions", envir = ns, inherits = TRUE)) {
    assign("build_internal_decisions",
           get("build_internal_decisions", envir = ns), envir = env)
  }
  if (exists("write_canonical_decision_output", envir = ns, inherits = TRUE)) {
    assign("write_canonical_decision_output",
           get("write_canonical_decision_output", envir = ns), envir = env)
  }
  if (exists("build_reference_decisions", envir = ns, inherits = TRUE)) {
    assign("build_reference_decisions",
           get("build_reference_decisions", envir = ns), envir = env)
  }
  env
}

#' @keywords internal
#' @noRd
.of_input_to_path <- function(spec, label, tmpdir = tempdir()) {
  if (is.null(spec)) return(NULL)
  if (!is.null(spec$path) && !is.na(spec$path) && nzchar(spec$path)) {
    if (!file.exists(spec$path)) {
      stop(sprintf("'%s' path does not exist: %s", label, spec$path),
           call. = FALSE)
    }
    return(spec$path)
  }
  val <- spec$value
  if (is.null(val)) {
    stop(sprintf("'%s' has neither a path nor an in-memory value.", label),
         call. = FALSE)
  }
  if (inherits(val, "SpatRaster")) {
    out <- tempfile(pattern = paste0(label, "_"), tmpdir = tmpdir,
                    fileext = ".tif")
    terra::writeRaster(val, out, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    return(out)
  }
  if (inherits(val, c("sf", "SpatVector"))) {
    out <- tempfile(pattern = paste0(label, "_"), tmpdir = tmpdir,
                    fileext = ".gpkg")
    sf::st_write(if (inherits(val, "SpatVector")) sf::st_as_sf(val) else val,
                 out, delete_dsn = TRUE, quiet = TRUE)
    return(out)
  }
  stop(sprintf("'%s' has an unsupported in-memory type.", label),
       call. = FALSE)
}

# ---------- canonical CORINE reclass matrix (from validated pipeline) -----
#' @keywords internal
#' @noRd
.of_corine_reclass_matrix <- function() {
  matrix(c(
    22,  1,   26, 2,   32, 3,   33, 3,
    27,  4,   28, 4,   29, 4,   23, 5,
    25,  6,   24, 7,
    1, 8,  2, 8,  3, 8,  4, 8,  5, 8,  6, 8,  7, 8,  8, 8,  9, 8, 10, 8, 11, 8,
    12, 9, 13, 9, 14, 9, 15, 9, 16, 9, 17, 9, 18, 9, 19, 9, 20, 9, 21, 9,
    35, 10, 36, 10, 37, 10, 38, 10, 39, 10, 40, 10, 41, 10, 42, 10, 43, 10, 44, 10,
    30, 11, 31, 11, 34, 11
  ), ncol = 2, byrow = TRUE)
}

# ---------- scenario presets (mirrors DETERMINISTIC_PIPELINE.R) -----------
#' @keywords internal
#' @noRd
.of_to_corine_list <- function(x) {
  # Mirror legacy to_corine_list(): wrap a class->value map in $corine so
  # the engine's resolve_*_local() helpers hit `map$corine[[key]]`.
  stopifnot(!is.null(names(x)))
  inner <- as.list(as.numeric(x))
  names(inner) <- names(x)
  list(corine = inner)
}
#' @keywords internal
#' @noRd
.of_to_corine_vec <- function(x) {
  # Mirror legacy to_corine_vec(): wrap a named numeric vector in $corine.
  stopifnot(!is.null(names(x)))
  list(corine = setNames(as.numeric(x), names(x)))
}

# Build the engine-shaped params (params_otsu / params_refine) from the
# resolved public methodological blocks stored on the config. The public
# surface (detect_params / refine_params / scoring_params) is resolved in
# build_burned_mapping_config(); here we only attach the fixed internal
# parameters and apply the engine-side $corine wrapping. No closed
# scenario presets remain.
#' @keywords internal
#' @noRd
.of_engine_params_from_config <- function(config) {
  det <- config$detect_params
  ref <- config$refine_params
  rsc <- config$rescue_params
  if (is.null(det) || is.null(ref) || is.null(rsc)) {
    stop("Internal error: config is missing resolved detect/refine/rescue ",
         "params. Rebuild it with build_burned_mapping_config().",
         call. = FALSE)
  }

  params_otsu <- list(
    otsu_value_range = c(0, 1500),
    trim_percentiles = list(min = 0.05, max = 0.99),
    otsu_thresholds = det$otsu_thresholds,
    grow_delta = det$grow_delta,
    min_grow_threshold_value = det$min_grow_threshold_value,
    otsu_min_by_class = .of_to_corine_list(det$otsu_min_by_class),
    grow_delta_by_class = .of_to_corine_list(det$grow_delta_by_class),
    min_grow_threshold_by_class = .of_to_corine_vec(det$min_grow_threshold_by_class),
    mask_edge_n_pixels = 1L,
    min_seed_pixels_per_component = det$min_seed_pixels_per_component,
    grow_engine = "whitebox", wbt_diag = FALSE,
    crop_to_units = TRUE, save_debug_rasters = FALSE,
    segment_by_intersection = TRUE, ecoregion_touches = TRUE
  )

  params_refine <- list(
    n_workers = 4, target_crs = "EPSG:3035",
    aoi_buffer_m = ref$aoi_buffer_m, omission_buffer_m = ref$omission_buffer_m,
    top_n_if_no_ref = Inf, min_detected_area_m2 = ref$min_detected_area_m2,
    clip_aois_to_raster_extent = TRUE, merge_overlaps = ref$merge_overlaps,
    merge_overlaps_buffer_m = ref$merge_overlaps_buffer_m, verbose = TRUE,
    rescue_args = list(
      trim_percentiles = list(min = 0.05, max = 0.99),
      otsu_value_range = c(0, 1500),
      otsu_thresholds = det$otsu_thresholds,
      otsu_min_by_class = .of_to_corine_list(rsc$otsu_min_by_class),
      grow_delta = det$grow_delta,
      grow_delta_by_class = .of_to_corine_list(rsc$grow_delta_by_class),
      min_grow_threshold_value = det$min_grow_threshold_value,
      min_grow_threshold_by_class = .of_to_corine_vec(rsc$min_grow_threshold_by_class),
      mask_edge_n_pixels = 1L,
      min_seed_pixels_per_component = det$min_seed_pixels_per_component,
      grow_engine = "whitebox", wbt_diag = FALSE,
      ecoregion_touches = TRUE, segment_by_intersection = TRUE,
      crop_to_units = TRUE, save_debug_rasters = FALSE,
      tile = FALSE, resolution = 90
    )
  )

  list(params_otsu = params_otsu, params_refine = params_refine)
}

#' @keywords internal
#' @noRd
.of_make_grow_tag <- function(otsu_threshold, grow_delta, seed_pixels,
                              grow_engine) {
  # Mirrors legacy make_grow_tag: e.g. "otsu310_d90_seed30_whitebox"
  fmt <- function(x) {
    s <- formatC(x, format = "fg", digits = 12, flag = "#")
    s <- gsub("\\s+", "", s); s <- sub("\\.?0+$", "", s)
    s <- gsub("-", "m", s);   s <- gsub("\\.", "p", s)
    s
  }
  paste0("otsu", fmt(otsu_threshold),
         "_d",  fmt(grow_delta),
         "_seed", as.integer(seed_pixels),
         "_",   grow_engine)
}

# ---------- DETECTION adapter --------------------------------------------
#' @keywords internal
#' @noRd
.of_run_detection <- function(config, aoi = NULL, write_outputs = TRUE,
                              overwrite = FALSE) {
  if (!inherits(config, "otsufire_burned_mapping_config")) {
    stop("Internal error: dispatcher expects otsufire_burned_mapping_config.",
         call. = FALSE)
  }
  engine_root <- config$engine_root %||% ""
  engine <- .of_load_engine(engine_root)

  change_index_path <- .of_input_to_path(config$inputs$change_index,
                                          "change_index")
  if (is.null(change_index_path)) {
    stop("detect_burned_patches(): config$inputs$change_index is required.",
         call. = FALSE)
  }
  corine_path <- .of_input_to_path(config$inputs$vegetation_map,
                                    "vegetation_map")
  if (is.null(corine_path)) {
    stop("detect_burned_patches(): config$inputs$vegetation_map (CORINE) is required.",
         call. = FALSE)
  }

  ecoregion_path <- config$options$ecoregion_shapefile_path
  if (is.null(ecoregion_path) || !nzchar(ecoregion_path)) {
    stop("detect_burned_patches(): supply config$options$ecoregion_shapefile_path.",
         call. = FALSE)
  }
  if (!file.exists(ecoregion_path)) {
    stop("ecoregion_shapefile_path does not exist: ", ecoregion_path,
         call. = FALSE)
  }

  peninsula_path <- config$options$peninsula_shapefile_path
  gdalwarp_path  <- config$tool_paths$gdalwarp_path %||%
                    config$options$gdalwarp_path

  reclassify_corine <- isTRUE(config$options$reclassify_corine %||% TRUE)
  if (reclassify_corine) {
    if (is.null(peninsula_path) || !file.exists(peninsula_path)) {
      stop("reclassify_corine=TRUE requires options$peninsula_shapefile_path to an existing file.",
           call. = FALSE)
    }
    if (is.null(gdalwarp_path) || !file.exists(gdalwarp_path)) {
      stop("reclassify_corine=TRUE requires tool_paths$gdalwarp_path (or options$gdalwarp_path) to an existing gdalwarp executable.",
           call. = FALSE)
    }
  }

  preset <- .of_engine_params_from_config(config)
  params_otsu   <- preset$params_otsu
  params_refine <- preset$params_refine

  grow_tag <- .of_make_grow_tag(
    otsu_threshold = params_otsu$otsu_thresholds,
    grow_delta     = params_otsu$grow_delta,
    seed_pixels    = params_otsu$min_seed_pixels_per_component,
    grow_engine    = params_otsu$grow_engine
  )

  base_dir        <- config$output_routes$base
  grow_dir        <- config$output_routes$grow_dir
  refine_dir      <- config$output_routes$refine_dir
  corine_reclass_dir <- file.path(base_dir, "00_CORINE_RECLASS")
  if (isTRUE(write_outputs)) {
    for (d in c(grow_dir, refine_dir, corine_reclass_dir)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  grow_shp_path <- file.path(
    grow_dir, paste0("BA_", config$target_year,
                     "_OTSUGROW_CORI_ECOREG_", grow_tag, ".shp")
  )
  refine_merged_path <- file.path(
    refine_dir, paste0("BA_", config$target_year,
                       "_REFINE_MERGED_", grow_tag, ".gpkg")
  )

  # ---- Stage 1: Otsu + grow -----------------------------------------
  grow_args <- list(
    raster_path = change_index_path,
    output_dir  = grow_dir,
    year        = config$target_year,
    peninsula_shapefile = peninsula_path,
    corine_raster_path = corine_path,
    reclassify_corine  = reclassify_corine,
    reclass_matrix     = .of_corine_reclass_matrix(),
    output_corine_raster_dir = corine_reclass_dir,
    corine_classes     = NULL,
    gdalwarp_path      = gdalwarp_path,
    ecoregion_shapefile_path = ecoregion_path,
    ecoregion_field    = "EnZ_name",
    segment_by_intersection = params_otsu$segment_by_intersection,
    ecoregion_touches  = params_otsu$ecoregion_touches,
    otsu_value_range   = params_otsu$otsu_value_range,
    trim_percentiles   = params_otsu$trim_percentiles,
    otsu_thresholds    = params_otsu$otsu_thresholds,
    otsu_min_by_class  = params_otsu$otsu_min_by_class,
    grow_delta         = params_otsu$grow_delta,
    grow_delta_by_class = params_otsu$grow_delta_by_class,
    min_grow_threshold_value    = params_otsu$min_grow_threshold_value,
    min_grow_threshold_by_class = params_otsu$min_grow_threshold_by_class,
    mask_edge_n_pixels = params_otsu$mask_edge_n_pixels,
    min_seed_pixels_per_component = params_otsu$min_seed_pixels_per_component,
    grow_engine        = params_otsu$grow_engine,
    wbt_diag           = params_otsu$wbt_diag,
    crop_to_units      = params_otsu$crop_to_units,
    save_debug_rasters = params_otsu$save_debug_rasters,
    random_seed        = config$deterministic_seed
  )

  det_over <- config$options$detect_overrides
  if (is.list(det_over)) {
    for (nm in names(det_over)) grow_args[[nm]] <- det_over[[nm]]
  }

  do.call(engine$process_otsu_rasters_grow, grow_args)
  if (!file.exists(grow_shp_path)) {
    stop("Stage 1 did not produce the expected grown shapefile: ",
         grow_shp_path, call. = FALSE)
  }

  # ---- Stage 2: segmentation_refinement (per-AOI) -------------------
  detected_polys <- sf::read_sf(grow_shp_path, quiet = TRUE)

  refine_args <- list(
    n_workers = params_refine$n_workers,
    raster_path = change_index_path,
    base_output_dir = refine_dir,
    year = config$target_year,
    detected_polys = detected_polys,
    reference_polys = NULL,
    binary_final_raster = NULL,
    target_crs = params_refine$target_crs,
    aoi_buffer_m = params_refine$aoi_buffer_m,
    omission_buffer_m = params_refine$omission_buffer_m,
    top_n_if_no_ref = params_refine$top_n_if_no_ref,
    min_detected_area_m2 = params_refine$min_detected_area_m2,
    clip_aois_to_raster_extent = params_refine$clip_aois_to_raster_extent,
    merge_overlaps = params_refine$merge_overlaps,
    merge_overlaps_buffer_m = params_refine$merge_overlaps_buffer_m,
    verbose = params_refine$verbose,
    corine_raster_path = corine_path,
    corine_classes = NULL,
    ecoregion_shapefile_path = ecoregion_path,
    ecoregion_field = "EnZ_name",
    process_fun = engine$process_otsu_rasters_grow,
    write_aoi_per_folder_shp = TRUE,
    standardize_ba_output = TRUE,
    keep_original_ba = FALSE,
    stop_on_unused_args = TRUE,
    random_seed = as.integer(config$deterministic_seed) + 1000L,
    rescue_args = params_refine$rescue_args
  )

  refine_over <- config$options$refine_overrides
  if (is.list(refine_over)) {
    for (nm in names(refine_over)) refine_args[[nm]] <- refine_over[[nm]]
  }

  # The refinement stage uses future.apply in parallel workers. F4b pre-loads
  # the ecoregion sf into the closure, which inflates serialized globals size
  # beyond the default 500 MiB ceiling. The validated pipeline sets 8 GB; we
  # match that for the duration of this call and restore afterwards.
  globals_budget <- config$options$future_globals_maxsize %||% (8 * 1024^3)
  old_globals <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = globals_budget)
  on.exit(options(future.globals.maxSize = old_globals), add = TRUE)

  do.call(engine$segmentation_refinement, refine_args)

  # ---- Stage 3: merge_aoi_shapefiles --------------------------------
  engine$merge_aoi_shapefiles(
    base_output_dir = refine_dir,
    pattern         = "^BA_AOI_\\d+\\.shp$",
    out_path        = refine_merged_path,
    dissolve        = FALSE,
    verbose         = TRUE
  )

  if (!file.exists(refine_merged_path)) {
    stop("Stage 3 did not produce the expected refined merged gpkg: ",
         refine_merged_path, call. = FALSE)
  }

  list(
    grow_result = NULL,
    grown_patches_path   = grow_shp_path,
    refined_patches_path = refine_merged_path,
    grow_tag             = grow_tag,
    detection_diagnostics = list(
      engine_root      = engine_root,
      scenario         = config$scenario,
      grow_tag         = grow_tag,
      params_otsu      = params_otsu,
      ran_refinement   = TRUE
    )
  )
}

#' @keywords internal
#' @noRd
.of_read_candidates <- function(burned_candidates) {
  if (inherits(burned_candidates, "sf")) return(burned_candidates)
  if (inherits(burned_candidates, "SpatVector"))
    return(sf::st_as_sf(burned_candidates))
  if (is.character(burned_candidates) && length(burned_candidates) == 1L &&
      file.exists(burned_candidates))
    return(sf::st_read(burned_candidates, quiet = TRUE))
  stop("burned_candidates could not be read.", call. = FALSE)
}

# ---------- SCORING adapter ----------------------------------------------
#' @keywords internal
#' @noRd
.of_run_scoring <- function(burned_candidates, config, keep_pool = NULL,
                            write_outputs = TRUE, overwrite = FALSE) {
  if (!inherits(config, "otsufire_burned_mapping_config")) {
    stop("Internal error: dispatcher expects otsufire_burned_mapping_config.",
         call. = FALSE)
  }
  engine_root <- config$engine_root %||% ""
  engine <- .of_load_engine(engine_root)

  change_index_path <- .of_input_to_path(config$inputs$change_index,
                                          "change_index")
  burnable_path <- .of_input_to_path(config$inputs$burnable_mask,
                                      "burnable_mask")
  preyear_path  <- .of_input_to_path(config$inputs$previous_year_burned,
                                      "previous_year_burned")
  refmap_path   <- .of_input_to_path(config$inputs$reference_burned_map,
                                      "reference_burned_map")
  hotspots_path <- .of_input_to_path(config$inputs$hotspots, "hotspots")

  base_dir     <- config$output_routes$base
  phase1_dir   <- file.path(base_dir, "03_SCORE_PHASE1")
  phase2_dir   <- file.path(base_dir, "04_SCORE_PHASE2")
  scoring_dir  <- config$output_routes$scoring_dir
  if (isTRUE(write_outputs)) {
    for (d in c(phase1_dir, phase2_dir, scoring_dir)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Read polys_stage2 (refined merged) and polys_stage1 (grown). The
  # public API passes `burned_candidates` as the stage2 layer; we
  # derive stage1 from the detection output routes.
  polys_stage2 <- .of_read_candidates(burned_candidates)

  grown_shp <- NULL
  det_paths <- config$options$.detect_paths
  if (is.list(det_paths) && !is.null(det_paths$grown_patches_path)) {
    grown_shp <- det_paths$grown_patches_path
  } else {
    # Infer from the grow_tag recorded in polys_stage2 attributes if present.
    grown_shp <- NULL
  }
  if (is.null(grown_shp) || !file.exists(grown_shp)) {
    # Fallback: try <grow_dir>/BA_<year>_OTSUGROW_CORI_ECOREG_*.shp
    cand <- list.files(
      config$output_routes$grow_dir,
      pattern = paste0("^BA_", config$target_year,
                       "_OTSUGROW_CORI_ECOREG_.*\\.shp$"),
      full.names = TRUE
    )
    if (length(cand) > 0L) grown_shp <- cand[1L]
  }
  if (is.null(grown_shp) || !file.exists(grown_shp)) {
    stop("score_burned_patches(): could not locate stage1 (grown) polygons.\n",
         "Supply options$.detect_paths$grown_patches_path explicitly, or run\n",
         "the full pipeline via run_deterministic_pipeline().", call. = FALSE)
  }
  polys_stage1 <- sf::read_sf(grown_shp, quiet = TRUE)

  preyear_polys <- NULL
  if (!is.null(preyear_path) && file.exists(preyear_path)) {
    preyear_polys <- sf::st_read(preyear_path, quiet = TRUE)
  }
  ref_polys <- NULL
  if (!is.null(refmap_path) && file.exists(refmap_path)) {
    ref_polys <- sf::st_read(refmap_path, quiet = TRUE)
  }

  burnable_corine <- terra::rast(burnable_path)
  rbr_rast <- terra::rast(change_index_path)
  if (terra::nlyr(rbr_rast) > 1L) rbr_rast <- rbr_rast[[1L]]

  # ---- Stage 4: scoring_burned_area_stage2 --------------------------
  # Public scoring surface comes from config$scoring_params (resolved in
  # build_burned_mapping_config()); the remaining constants stay fixed.
  scoring <- config$scoring_params
  if (is.null(scoring)) {
    stop("Internal error: config is missing resolved scoring_params. ",
         "Rebuild it with build_burned_mapping_config().", call. = FALSE)
  }
  stage2_args <- list(
    polys_stage2 = polys_stage2,
    polys_stage1 = polys_stage1,
    burnable_corine = burnable_corine,
    burnable_classes = NULL,
    rbr_rast = rbr_rast,
    corine_na_scope = "all",
    max_corine_na_frac = 0.65,
    support_buffer_m = scoring$support_buffer_m,
    preyear_polys = preyear_polys,
    erase_mask_buffer_m = scoring$erase_mask_buffer_m,
    erase_post_shave_m = scoring$erase_post_shave_m,
    erase_min_area_m2 = scoring$erase_min_area_m2,
    ref_polys = ref_polys,
    score_ref_polys = !is.null(ref_polys),
    ref_buffer_m = scoring$ref_buffer_m,
    internal_keep_value = "keep",
    validate_use_clean = TRUE,
    validate_ref_use_clean = TRUE,
    validate_ref_use_kept = TRUE,
    save_outputs = isTRUE(write_outputs),
    out_dir = phase1_dir,
    prefix = paste0("BA_", config$target_year, "_", config$scenario),
    driver = "GPKG", overwrite = TRUE, quiet = TRUE,
    save_splits = FALSE, suppress_sf_warnings = TRUE
  )
  stage2_over <- config$options$stage2_overrides
  if (is.list(stage2_over)) {
    for (nm in names(stage2_over)) stage2_args[[nm]] <- stage2_over[[nm]]
  }
  res_phase1 <- do.call(engine$scoring_burned_area_stage2, stage2_args)

  stage2_keep_review <- res_phase1$internal$stage2_keep_review
  stage2_flagged     <- res_phase1$internal$stage2_flagged
  stage2_clean       <- res_phase1$internal$stage2_keep_review_clean

  if (is.null(stage2_keep_review) || nrow(stage2_keep_review) == 0L) {
    stop("Stage 4 (scoring_burned_area_stage2) produced an empty internal$stage2_keep_review.",
         call. = FALSE)
  }

  # ---- Stage 5a: optional registry-based keep pool ------------------
  registry_keep_pool <- NULL
  if (!is.null(keep_pool)) {
    registry_keep_pool <- keep_pool
  } else if (!is.null(config$registry_path) &&
             file.exists(config$registry_path)) {
    registry_keep_pool <- tryCatch(
      engine$build_rbr_keep_pool_from_registry(
        registry_path = config$registry_path,
        rbr_rast = rbr_rast,
        registry_layer = "burned_high_conf_registry",
        target_year = config$target_year,
        verbose = TRUE
      ),
      error = function(e) {
        warning("build_rbr_keep_pool_from_registry failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  # ---- Stage 5b: score_rbr_keep_classes (internal) ------------------
  rbr_args_internal <- list(
    polys_sf = stage2_keep_review,
    rbr_rast = rbr_rast,
    out_dir  = phase2_dir,
    prefix   = paste0("internal_", config$target_year, "_",
                      config$scenario, "_p10_p5"),
    keep_pool = if (!is.null(registry_keep_pool))
                  registry_keep_pool$keep_pool else NULL,
    use_keep_common_pool = FALSE,
    keep_fallback_col = "flag_internal",
    keep_fallback_value = "keep",
    max_keep_samples = 5e6,
    above_keep_prob = 0.05,
    promote_p_above_ref = 0.50,
    promote_percentile = 0.10,
    min_area_ha = 10, min_pix = NULL,
    fill_na_corine = TRUE, fill_na_corine_value = 0,
    chunk_size = 500, verbose = TRUE,
    save_outputs = isTRUE(write_outputs),
    driver = "GPKG", overwrite = TRUE, quiet = TRUE,
    save_splits = FALSE, arcgis_fix = TRUE, output_epsg = 3035,
    random_seed = as.integer(config$deterministic_seed) + 2000L
  )
  rbr_over <- config$options$rbr_overrides
  if (is.list(rbr_over)) {
    for (nm in names(rbr_over)) rbr_args_internal[[nm]] <- rbr_over[[nm]]
  }
  res_internal_rbr <- do.call(engine$score_rbr_keep_classes,
                              rbr_args_internal)

  # ---- Stage 5c: score_rbr_keep_classes (external/reference) --------
  res_external_rbr <- NULL
  ref_flagged <- res_phase1$reference$ref_flagged
  if (!is.null(ref_flagged) && inherits(ref_flagged, "sf") &&
      nrow(ref_flagged) > 0L) {
    res_external_rbr <- tryCatch(
      engine$score_rbr_keep_classes(
        polys_sf  = ref_flagged,
        rbr_rast  = rbr_rast,
        out_dir   = phase2_dir,
        prefix    = paste0("reference_", config$target_year, "_",
                           config$scenario, "_p10_p5"),
        keep_pool = res_internal_rbr$keep_pool,
        save_outputs = isTRUE(write_outputs),
        driver = "GPKG", overwrite = TRUE, quiet = TRUE,
        random_seed = as.integer(config$deterministic_seed) + 3000L
      ),
      error = function(e) {
        warning("score_rbr_keep_classes external failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  # ---- Stage 6: canonical decisions ---------------------------------
  hotspots_sf <- NULL
  if (!is.null(hotspots_path) && file.exists(hotspots_path)) {
    hotspots_sf <- suppressWarnings(sf::st_make_valid(
      sf::st_read(hotspots_path, quiet = TRUE)
    ))
  }
  preyear_available <- !is.null(preyear_polys) &&
                        inherits(preyear_polys, "sf") &&
                        nrow(preyear_polys) > 0L

  internal_decisions_sf <- engine$build_internal_decisions(
    internal_flagged = stage2_flagged,
    internal_clean   = stage2_clean,
    internal_rbr     = res_internal_rbr$scored_sf,
    hotspots_sf      = hotspots_sf,
    target_year      = config$target_year,
    scenario_name    = config$scenario,
    preyear_available = preyear_available
  )

  internal_decisions_path <- NA_character_
  if (isTRUE(write_outputs) && !is.null(internal_decisions_sf) &&
      nrow(internal_decisions_sf) > 0L) {
    target <- config$output_routes$internal_decisions
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(target)) {
      if (!isTRUE(overwrite)) {
        stop("internal_decisions already exists and overwrite = FALSE: ",
             target, call. = FALSE)
      }
      try(unlink(target, force = TRUE), silent = TRUE)
    }
    internal_decisions_path <- engine$write_canonical_decision_output(
      sf_obj = internal_decisions_sf,
      out_dir = dirname(target),
      gpkg_name = basename(target),
      layer_name = "internal_decisions",
      overwrite = TRUE, quiet = TRUE,
      arcgis_fix = TRUE, output_epsg = 3035
    )
  }

  reference_decisions_sf <- NULL
  reference_decisions_path <- NA_character_
  if (!is.null(res_external_rbr) && !is.null(ref_flagged) &&
      !is.null(res_phase1$validation$ref_flagged)) {
    reference_decisions_sf <- tryCatch(
      engine$build_reference_decisions(
        ref_flagged    = ref_flagged,
        ref_validation = res_phase1$validation$ref_flagged,
        ref_rbr        = res_external_rbr$scored_sf,
        target_year    = config$target_year,
        scenario_name  = config$scenario
      ),
      error = function(e) {
        warning("build_reference_decisions failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (isTRUE(write_outputs) && !is.null(reference_decisions_sf) &&
        nrow(reference_decisions_sf) > 0L) {
      target <- config$output_routes$reference_decisions
      reference_decisions_path <- engine$write_canonical_decision_output(
        sf_obj = reference_decisions_sf,
        out_dir = dirname(target),
        gpkg_name = basename(target),
        layer_name = "reference_decisions",
        overwrite = TRUE, quiet = TRUE,
        arcgis_fix = TRUE, output_epsg = 3035
      )
    }
  }

  list(
    internal_decisions = internal_decisions_sf,
    internal_decisions_path = internal_decisions_path,
    reference_decisions = reference_decisions_sf,
    reference_decisions_path = reference_decisions_path,
    keep_pool_summary = res_internal_rbr$keep_pool,
    scoring_diagnostics = list(
      engine_root = engine_root,
      scenario = config$scenario,
      n_stage1 = nrow(polys_stage1),
      n_stage2 = nrow(polys_stage2),
      n_stage2_keep_review = nrow(stage2_keep_review),
      n_decisions = if (is.null(internal_decisions_sf)) 0L
                    else nrow(internal_decisions_sf),
      registry_keep_pool_used = !is.null(registry_keep_pool)
    )
  )
}
