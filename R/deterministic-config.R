#' Build a deterministic burned-mapping configuration object
#'
#' @description
#' Public entry point for the deterministic workflow. Collects the spatial
#' inputs, the methodological parameters, the output location, and the
#' technical/runtime options needed by the rest of the deterministic
#' stages into a single normalized object of class
#' `otsufire_burned_mapping_config`. The function does not run any heavy
#' processing: it validates arguments, resolves the methodological
#' parameter blocks, and computes registry and output routes.
#'
#' Pass the returned object to [detect_burned_patches()],
#' [score_burned_patches()], or [run_deterministic_pipeline()].
#'
#' The public methodological surface is exposed through three named
#' lists: `detect_params`, `refine_params`, and `scoring_params`. Each
#' accepts a small, fixed set of environmental-user-facing keys (see
#' the parameter sections below). Any other internal/engine parameter is
#' fixed to its validated default and is not configurable through this
#' public API.
#'
#' @details
#' The public parameter surface is organized as ecological controls rather
#' than engine internals.
#'
#' **Vegetation classes used by `*_by_vegetation` vectors**
#' - `1` Agroforestry
#' - `2` Grassland
#' - `3` Sparse / burned
#' - `4` Shrubland
#' - `5` Broadleaved forest
#' - `6` Mixed forest
#' - `7` Conifer forest
#' - `8` Artificial areas
#' - `9` Agricultural areas
#' - `10` Water / wetlands
#' - `11` Bare
#'
#' **How to interpret `detect_params`**
#' - `seed_threshold` (change-index units): global minimum change-index
#'   value required for a pixel to become a burned seed. Increasing it
#'   creates fewer, more conservative seeds with lower commission but
#'   higher omission; decreasing it creates more seeds and higher
#'   sensitivity but raises false-positive risk. Adjust when the whole
#'   map is systematically over- or under-detecting burned cores. Acts as
#'   a global floor that class-specific bounds can override.
#' - `seed_threshold_by_vegetation` (change-index units): class-specific
#'   lower bound for seed generation. Increasing a class value makes seed
#'   detection more conservative in that land-cover class; decreasing it
#'   makes seed generation more permissive there. Adjust when errors are
#'   concentrated in specific vegetation or land-cover types. The
#'   validated defaults are highest in water / wetlands and artificial /
#'   bare classes, and lower in vegetated burnable classes.
#' - `growth_delta` (change-index units): global reduction applied from
#'   seed threshold to growth threshold. Increasing it allows more
#'   permissive expansion from seeds, larger patches, and higher
#'   commission risk; decreasing it makes expansion more restrictive and
#'   patches smaller, at the cost of higher omission. Adjust when mapped
#'   patch boundaries are consistently too small or too large. Acts as a
#'   general fallback when no class-specific delta is supplied.
#' - `growth_delta_by_vegetation` (change-index units): class-specific
#'   amount by which the growth threshold is relaxed relative to the seed
#'   threshold. Increasing a class value allows stronger expansion within
#'   that class and increases connectivity; decreasing it limits expansion
#'   there and produces tighter patch boundaries. Adjust when over- or
#'   under-expansion is vegetation-specific. The validated defaults are
#'   larger in shrubland and forest classes, and very small in water /
#'   wetlands, artificial, agricultural, and bare classes.
#' - `minimum_growth_threshold` (change-index units): global minimum
#'   change-index floor allowed during region growing. Increasing it
#'   restricts growth to stronger change-index values and reduces
#'   over-expansion; decreasing it allows growth into weaker change-index
#'   values and can absorb unburned pixels. Adjust when growing
#'   systematically leaks into weak-change backgrounds or misses
#'   low-severity burn edges. Acts as a global floor that class-specific
#'   floors usually override.
#' - `minimum_growth_threshold_by_vegetation` (change-index units):
#'   class-specific minimum floor constraining region growing. Increasing
#'   a class value prevents growth into low-change pixels within that
#'   class and reduces commission; decreasing it allows growth into lower
#'   change pixels and can improve boundary completeness at the cost of
#'   more false positives. Adjust when boundary errors are class-specific.
#'   The validated defaults are highest in water / wetlands, artificial,
#'   agricultural, and bare classes, and lower in shrubland and forested
#'   classes.
#' - `minimum_seed_pixels` (pixels): minimum number of core seed pixels
#'   required to retain a grown patch. Increasing it removes small
#'   isolated detections and reduces salt-and-pepper false positives;
#'   decreasing it keeps smaller candidate patches and improves small-fire
#'   sensitivity at the cost of more noise. Adjust according to minimum
#'   fire size and raster resolution. Interpret this parameter together
#'   with pixel size: 30 pixels correspond to different areas at 20 m and
#'   30 m resolution.
#'
#' **How to interpret `refine_params`**
#' - `aoi_buffer_m` (meters): buffer around the area of interest used
#'   during refinement. Increasing it reduces edge artefacts near the AOI
#'   boundary but increases processing area; decreasing it reduces
#'   processing area but increases truncation risk at the edges. Adjust
#'   when fires near AOI boundaries are clipped or when processing cost is
#'   excessive. Not vegetation-specific.
#' - `minimum_detected_area_m2` (m^2): minimum area of detected polygons
#'   retained after segmentation. Increasing it removes tiny patches
#'   earlier in the workflow; decreasing it retains very small detections
#'   for later filtering. Adjust when segmentation produces excessive tiny
#'   fragments. This threshold is usually kept permissive because later
#'   filters handle patch reliability.
#' - `merge_overlaps` (logical): whether overlapping candidate geometries
#'   are merged during refinement. Setting it to `TRUE` produces more
#'   consolidated patches and avoids duplicate candidates; setting it to
#'   `FALSE` preserves overlapping detections separately. Adjust only if
#'   overlap structure must be preserved for diagnostics. Not
#'   vegetation-specific.
#'
#' **How to interpret `scoring_params`**
#' - `support_buffer_m` (meters): buffer used to evaluate contextual
#'   support around candidate patches. Increasing it allows nearby
#'   support evidence to influence the patch; decreasing it restricts
#'   support evidence to the patch geometry itself. Adjust when evidence
#'   is spatially misaligned with patch boundaries. Usually conservative;
#'   hotspot and reference buffers are handled separately.
#' - `previous_fire_exclusion_buffer_m` (meters): buffer around
#'   previous-year burned areas used to identify temporal conflicts.
#'   Increasing it more aggressively removes or downgrades patches
#'   overlapping recent burns; decreasing it is more permissive with
#'   repeated or adjacent burns. Adjust when previous-year scars are being
#'   repeatedly detected as new burns. The validated default of 90 m is
#'   roughly one 90 m analysis pixel.
#' - `previous_fire_cleanup_buffer_m` (meters): buffer used to clean
#'   residual overlap with previous-year burned polygons. Increasing it
#'   removes more temporal-cleanup residue; decreasing it retains more
#'   boundary-adjacent area near old burns. Adjust when temporal-cleanup
#'   artefacts are visible along old fire boundaries. Not
#'   vegetation-specific.
#' - `minimum_remaining_area_m2` (m^2): minimum area required after
#'   previous-fire overlap removal. Increasing it drops more small
#'   remnants after temporal cleanup; decreasing it retains more residual
#'   fragments. Adjust when previous-fire erasure leaves many meaningless
#'   slivers or removes valid reburns. It should reflect the minimum
#'   meaningful patch size after cleanup.
#' - `reference_buffer_m` (meters): buffer used when constructing or
#'   comparing spectral support references. Increasing it allows more
#'   spatial tolerance around reference or support evidence; decreasing it
#'   makes the comparison stricter and more spatially exact. Adjust when
#'   spatial mismatch exists between candidate patches and support /
#'   reference layers. Not vegetation-specific.
#'
#' @param change_index Raster path, `terra::SpatRaster`, or `NULL`.
#'   Annual change-index raster (RBR, dNBR, RdNBR, ...). Required before
#'   detection runs; may be `NULL` at config-build time for configs that
#'   will be used only with pre-built candidates (e.g. direct calls to
#'   [score_burned_patches()]).
#' @param vegetation_map Raster, sf object, or path. Optional vegetation or
#'   land-cover stratification layer (CORINE-compatible in 0.2.x).
#' @param burnable_mask Raster, sf object, or path. Required: binary mask
#'   of burnable area (1 burnable, 0 non-burnable).
#' @param hotspots sf POINT layer or path. Optional hotspot layer for the
#'   target year.
#' @param previous_year_burned sf POLYGON layer or path. Optional
#'   previous-year burned map for temporal conflict assessment.
#' @param reference_burned_map sf POLYGON, raster, or path. Optional burned
#'   reference used by [validate_fire_maps()] when validation is requested.
#' @param target_year Integer scalar. Target year. Required: used for
#'   naming outputs and filtering hotspots/previous-year layers.
#' @param output_dir Character scalar. Root output directory.
#'   Defaults to `tempdir()`.
#' @param run_name Character scalar. The single visible experiment
#'   identifier. Used to name output folders/files, the
#'   `.../DETERMINISTIC/<run_name>/...` route level, and the
#'   `Results/<toupper(run_name)>/RBR_TRAINING_REGISTRY.gpkg` registry
#'   path. Must be a single non-empty string and must not contain `/` or
#'   `\\`. Defaults to `"deterministic_burned_map"`.
#' @param detect_params Named list. Detection parameters. Recognized keys:
#'   `seed_threshold` (default `310`), `growth_delta` (default `90`),
#'   `minimum_growth_threshold` (default `240`), `minimum_seed_pixels`
#'   (default `30`), and the per-vegetation-class vectors
#'   `seed_threshold_by_vegetation`, `growth_delta_by_vegetation`,
#'   `minimum_growth_threshold_by_vegetation` (named numeric vectors over
#'   classes `"1"`..`"11"`). If a global scalar is given without its
#'   matching `*_by_vegetation` vector, the global value is broadcast to
#'   all 11 classes. A partial `*_by_vegetation` vector is completed with
#'   the user global (if given) or otherwise the legacy default for the
#'   missing class. Unknown keys raise an error.
#' @param refine_params Named list. Refinement parameters. Recognized
#'   keys: `aoi_buffer_m` (default `5000`), `omission_buffer_m`
#'   (default `0`), `minimum_detected_area_m2` (default `10`),
#'   `merge_overlaps` (default `TRUE`), `merge_overlaps_buffer_m`
#'   (default `0`). Unknown keys raise an error.
#' @param scoring_params Named list. Scoring parameters. Recognized keys:
#'   `support_buffer_m` (default `0`),
#'   `previous_fire_exclusion_buffer_m` (default `90`),
#'   `previous_fire_cleanup_buffer_m` (default `90`),
#'   `minimum_remaining_area_m2` (default `20000`),
#'   `reference_buffer_m` (default `90`). Unknown keys raise an error.
#' @param options Named list. Technical/runtime container only (NOT a
#'   methodological surface). Recognized keys: `deterministic_seed`,
#'   `engine_root`, `registry_path`, `whitebox_exe`, `gdalwarp_path`,
#'   `tool_paths`, `ecoregion_shapefile_path`, `peninsula_shapefile_path`,
#'   `future_globals_maxsize`, `change_index_validation`.
#'
#'   `options$change_index_validation` is an optional advanced sub-list
#'   that tunes the in-memory sanity check applied to `change_index` by
#'   [detect_burned_patches()] (it is never applied to `vegetation_map`
#'   or `burnable_mask`, and no raster is read at config-build time).
#'   Recognized sub-keys (all optional; defaults shown):
#'   `expected_nodata` (`-9999`), `lower_cap` (`-1000`), `cap_below`
#'   (`TRUE`), `check_finite` (`TRUE`). `lower_cap` is configurable on
#'   purpose: the validated value `-1000` suits RBR, but other change
#'   indices have different valid ranges, so it must not be hard-coded.
#'   If supplied, the block must be a named list using only these keys.
#'
#' @return An S3 object of class `otsufire_burned_mapping_config` (a named
#'   list). Stable fields: `target_year`, `inputs`, `run_name`,
#'   `output_routes`, `detect_params`, `refine_params`, `scoring_params`,
#'   `rescue_params`, `options`, `engine_root`, `deterministic_seed`,
#'   `tool_paths`, `registry_path`. The `detect_params`, `refine_params`,
#'   and `scoring_params` fields hold the fully resolved parameter blocks
#'   in their internal engine-name shape.
#'
#' @family workflow
#' @export
build_burned_mapping_config <- function(
    change_index = NULL,
    vegetation_map = NULL,
    burnable_mask,
    hotspots = NULL,
    previous_year_burned = NULL,
    reference_burned_map = NULL,
    target_year,
    output_dir = tempdir(),
    run_name = "deterministic_burned_map",
    detect_params = list(),
    refine_params = list(),
    scoring_params = list(),
    options = list()
) {
  if (missing(target_year) || is.null(target_year)) {
    stop("'target_year' is required.", call. = FALSE)
  }
  target_year <- suppressWarnings(as.integer(target_year)[1L])
  if (is.na(target_year) || target_year < 1900L || target_year > 2100L) {
    stop("'target_year' must be an integer in [1900, 2100].", call. = FALSE)
  }

  if (missing(burnable_mask) || is.null(burnable_mask)) {
    stop("'burnable_mask' is required.", call. = FALSE)
  }

  if (!is.character(run_name) || length(run_name) != 1L || is.na(run_name) ||
      !nzchar(run_name)) {
    stop("'run_name' must be a single non-empty character string.", call. = FALSE)
  }
  if (grepl("[/\\\\]", run_name)) {
    stop("'run_name' must not contain '/' or '\\\\' (it is used in output paths and the registry path).",
         call. = FALSE)
  }

  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("'output_dir' must be a single non-empty character string.", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  if (!is.list(options)) {
    stop("'options' must be a list().", call. = FALSE)
  }
  if (length(options) > 0L && (is.null(names(options)) || any(!nzchar(names(options))))) {
    stop("'options' must be a named list (all elements must have names).", call. = FALSE)
  }

  # Advanced change-index validation block: only structural validation
  # here (named list, recognized keys). It is consumed at detection time
  # by detect_burned_patches(); no raster is read here.
  ci_val <- options[["change_index_validation"]]
  if (!is.null(ci_val)) {
    .of_check_named_list(ci_val, "options$change_index_validation")
    .of_check_unknown_keys(
      ci_val, names(.of_change_index_validation_defaults()),
      "options$change_index_validation")
  }

  # ---- methodological parameter resolution --------------------------
  detect_resolved  <- .of_resolve_detect_params(detect_params)
  refine_resolved  <- .of_resolve_refine_params(refine_params)
  scoring_resolved <- .of_resolve_scoring_params(scoring_params)
  rescue_resolved  <- .of_derive_rescue_params(detect_resolved)

  inputs <- list(
    change_index = .of_normalize_input_spec(change_index, "change_index",
                                            allow_null = TRUE),
    vegetation_map = .of_normalize_input_spec(vegetation_map, "vegetation_map",
                                              allow_null = TRUE),
    burnable_mask = .of_normalize_input_spec(burnable_mask, "burnable_mask",
                                              allow_null = FALSE),
    hotspots = .of_normalize_input_spec(hotspots, "hotspots",
                                         allow_null = TRUE),
    previous_year_burned = .of_normalize_input_spec(previous_year_burned,
                                                     "previous_year_burned",
                                                     allow_null = TRUE),
    reference_burned_map = .of_normalize_input_spec(reference_burned_map,
                                                     "reference_burned_map",
                                                     allow_null = TRUE)
  )

  deterministic_seed <- .of_coerce_seed(options$deterministic_seed,
                                        default = 12345L)

  tool_paths <- list(
    whitebox_exe  = options$whitebox_exe  %||% options$tool_paths$whitebox_exe,
    gdalwarp_path = options$gdalwarp_path %||% options$tool_paths$gdalwarp_path
  )

  # Block 4A: deterministic engine is now in the package itself; engine_root
  # is optional and only used for legacy diagnostics. Resolution is best-effort.
  engine_root <- options$engine_root %||% .of_find_engine_root()
  if (!is.null(engine_root)) {
    engine_root <- normalizePath(engine_root, winslash = "/", mustWork = FALSE)
    if (!dir.exists(engine_root)) engine_root <- NULL
  }

  output_routes <- .of_build_output_routes(
    output_dir = output_dir,
    target_year = target_year,
    run_name = run_name
  )

  registry_path <- options$registry_path %||% .of_resolve_registry_path(
    output_dir = output_dir,
    run_name = run_name,
    explicit = options$registry_path
  )

  cfg <- list(
    target_year = target_year,
    inputs = inputs,
    run_name = run_name,
    output_dir = output_dir,
    output_routes = output_routes,
    registry_path = registry_path,
    engine_root = engine_root,
    deterministic_seed = deterministic_seed,
    tool_paths = tool_paths,
    detect_params = detect_resolved,
    refine_params = refine_resolved,
    scoring_params = scoring_resolved,
    rescue_params = rescue_resolved,
    options = options
  )

  class(cfg) <- c("otsufire_burned_mapping_config", "list")
  cfg
}

#' @export
print.otsufire_burned_mapping_config <- function(x, ...) {
  cat("<otsufire_burned_mapping_config>\n")
  cat("  target_year  :", x$target_year, "\n")
  cat("  run_name     :", x$run_name, "\n")
  cat("  output_dir   :", x$output_dir, "\n")
  cat("  engine_root  :", x$engine_root %||% "<unresolved>", "\n")
  cat("  deterministic_seed:", x$deterministic_seed, "\n")
  cat("  registry_path:", x$registry_path %||% "<auto>", "\n")
  cat("  detect_params: seed_threshold=", x$detect_params$otsu_thresholds,
      " growth_delta=", x$detect_params$grow_delta,
      " minimum_growth_threshold=", x$detect_params$min_grow_threshold_value,
      " minimum_seed_pixels=", x$detect_params$min_seed_pixels_per_component,
      "\n", sep = "")
  cat("  refine_params: aoi_buffer_m=", x$refine_params$aoi_buffer_m,
      " minimum_detected_area_m2=", x$refine_params$min_detected_area_m2,
      "\n", sep = "")
  cat("  scoring_params: support_buffer_m=", x$scoring_params$support_buffer_m,
      " reference_buffer_m=", x$scoring_params$ref_buffer_m,
      "\n", sep = "")
  cat("  inputs       :\n")
  for (nm in names(x$inputs)) {
    v <- x$inputs[[nm]]
    cat(sprintf("    - %-22s: %s\n", nm,
                if (is.null(v)) "<NULL>" else
                  sprintf("%s [%s]",
                          if (!is.null(v$path)) v$path else "<in-memory>",
                          v$type)))
  }
  invisible(x)
}

# --- methodological parameter resolution -------------------------------

# Vegetation classes the per-class vectors are defined over.
.OF_VEG_CLASSES <- as.character(1:11)

#' @keywords internal
#' @noRd
.of_detect_param_defaults <- function() {
  list(
    seed_threshold = 310,
    growth_delta = 90,
    minimum_growth_threshold = 240,
    minimum_seed_pixels = 30L,
    seed_threshold_by_vegetation = c(
      "1" = 310, "2" = 420, "3" = 380, "4" = 320, "5" = 320, "6" = 320,
      "7" = 320, "8" = 500, "9" = 470, "10" = 670, "11" = 480),
    growth_delta_by_vegetation = c(
      "1" = 80, "2" = 20, "3" = 40, "4" = 100, "5" = 85, "6" = 85,
      "7" = 85, "8" = 10, "9" = 10, "10" = 5, "11" = 10),
    minimum_growth_threshold_by_vegetation = c(
      "1" = 230, "2" = 360, "3" = 300, "4" = 230, "5" = 240, "6" = 240,
      "7" = 240, "8" = 430, "9" = 400, "10" = 560, "11" = 400)
  )
}

#' @keywords internal
#' @noRd
.of_check_named_list <- function(x, label) {
  if (!is.list(x)) {
    stop(sprintf("'%s' must be a list().", label), call. = FALSE)
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop(sprintf("'%s' must be a named list (all elements must have names).",
                 label), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.of_check_unknown_keys <- function(x, allowed, label) {
  unknown <- setdiff(names(x), allowed)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown %s key(s): %s. Allowed: %s.",
                 label, paste(unknown, collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}

# Default tuning for the in-memory change-index sanity check. lower_cap
# is intentionally a default, not a hard-coded constant: RBR uses -1000
# but other change indices have different valid ranges. Shared by the
# builder (key validation) and detect_burned_patches() (resolution).
#' @keywords internal
#' @noRd
.of_change_index_validation_defaults <- function() {
  list(
    expected_nodata = -9999,
    lower_cap       = -1000,
    cap_below       = TRUE,
    check_finite    = TRUE
  )
}

#' @keywords internal
#' @noRd
.of_check_scalar_num <- function(x, label) {
  if (is.null(x)) return(invisible(NULL))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric value.", label),
         call. = FALSE)
  }
  invisible(NULL)
}

# Resolve one global/by-vegetation triple into a complete named numeric
# vector over all 11 classes, following the closed resolution rules.
#' @keywords internal
#' @noRd
.of_resolve_by_veg <- function(user_global, user_byveg,
                               default_global, default_byveg, label) {
  classes <- .OF_VEG_CLASSES
  if (!is.null(user_byveg)) {
    if (!is.numeric(user_byveg) || is.null(names(user_byveg)) ||
        any(!nzchar(names(user_byveg)))) {
      stop(sprintf("detect_params$%s must be a named numeric vector.",
                   label), call. = FALSE)
    }
    unknown <- setdiff(names(user_byveg), classes)
    if (length(unknown) > 0L) {
      stop(sprintf(
        "detect_params$%s has unknown vegetation class(es): %s (expected '1'..'11').",
        label, paste(unknown, collapse = ", ")), call. = FALSE)
    }
  }
  out <- stats::setNames(numeric(length(classes)), classes)
  for (cl in classes) {
    if (!is.null(user_byveg) && cl %in% names(user_byveg)) {
      out[cl] <- as.numeric(user_byveg[[cl]])
    } else if (!is.null(user_global)) {
      out[cl] <- as.numeric(user_global)
    } else {
      out[cl] <- as.numeric(default_byveg[[cl]])
    }
  }
  out
}

#' @keywords internal
#' @noRd
.of_resolve_detect_params <- function(detect_params) {
  .of_check_named_list(detect_params, "detect_params")
  allowed <- c(
    "seed_threshold", "growth_delta", "minimum_growth_threshold",
    "minimum_seed_pixels", "seed_threshold_by_vegetation",
    "growth_delta_by_vegetation", "minimum_growth_threshold_by_vegetation")
  .of_check_unknown_keys(detect_params, allowed, "detect_params")

  d <- .of_detect_param_defaults()

  # Exact extraction only: `$` would partial-match e.g. `growth_delta`
  # against `growth_delta_by_vegetation`.
  seed_g  <- detect_params[["seed_threshold"]]
  delta_g <- detect_params[["growth_delta"]]
  floor_g <- detect_params[["minimum_growth_threshold"]]
  .of_check_scalar_num(seed_g,  "detect_params$seed_threshold")
  .of_check_scalar_num(delta_g, "detect_params$growth_delta")
  .of_check_scalar_num(floor_g, "detect_params$minimum_growth_threshold")

  msp <- detect_params[["minimum_seed_pixels"]]
  if (is.null(msp)) {
    msp <- d$minimum_seed_pixels
  } else {
    .of_check_scalar_num(msp, "detect_params$minimum_seed_pixels")
    msp <- as.integer(round(msp))
  }

  list(
    otsu_thresholds = if (!is.null(seed_g)) as.numeric(seed_g)
                      else d$seed_threshold,
    grow_delta = if (!is.null(delta_g)) as.numeric(delta_g)
                 else d$growth_delta,
    min_grow_threshold_value = if (!is.null(floor_g)) as.numeric(floor_g)
                               else d$minimum_growth_threshold,
    min_seed_pixels_per_component = msp,
    otsu_min_by_class = .of_resolve_by_veg(
      seed_g, detect_params[["seed_threshold_by_vegetation"]],
      d$seed_threshold, d$seed_threshold_by_vegetation,
      "seed_threshold_by_vegetation"),
    grow_delta_by_class = .of_resolve_by_veg(
      delta_g, detect_params[["growth_delta_by_vegetation"]],
      d$growth_delta, d$growth_delta_by_vegetation,
      "growth_delta_by_vegetation"),
    min_grow_threshold_by_class = .of_resolve_by_veg(
      floor_g, detect_params[["minimum_growth_threshold_by_vegetation"]],
      d$minimum_growth_threshold, d$minimum_growth_threshold_by_vegetation,
      "minimum_growth_threshold_by_vegetation")
  )
}

#' @keywords internal
#' @noRd
.of_resolve_refine_params <- function(refine_params) {
  .of_check_named_list(refine_params, "refine_params")
  allowed <- c("aoi_buffer_m", "omission_buffer_m",
               "minimum_detected_area_m2", "merge_overlaps",
               "merge_overlaps_buffer_m")
  .of_check_unknown_keys(refine_params, allowed, "refine_params")

  g <- function(key, default) {
    if (is.null(refine_params[[key]])) default else refine_params[[key]]
  }

  merge_overlaps <- g("merge_overlaps", TRUE)
  if (!is.logical(merge_overlaps) || length(merge_overlaps) != 1L ||
      is.na(merge_overlaps)) {
    stop("refine_params$merge_overlaps must be a single TRUE/FALSE value.",
         call. = FALSE)
  }
  out <- list(
    aoi_buffer_m            = g("aoi_buffer_m", 5000),
    omission_buffer_m       = g("omission_buffer_m", 0),
    min_detected_area_m2    = g("minimum_detected_area_m2", 10),
    merge_overlaps          = merge_overlaps,
    merge_overlaps_buffer_m = g("merge_overlaps_buffer_m", 0)
  )
  .of_check_scalar_num(out$aoi_buffer_m, "refine_params$aoi_buffer_m")
  .of_check_scalar_num(out$omission_buffer_m,
                       "refine_params$omission_buffer_m")
  .of_check_scalar_num(out$min_detected_area_m2,
                       "refine_params$minimum_detected_area_m2")
  .of_check_scalar_num(out$merge_overlaps_buffer_m,
                       "refine_params$merge_overlaps_buffer_m")
  out
}

#' @keywords internal
#' @noRd
.of_resolve_scoring_params <- function(scoring_params) {
  .of_check_named_list(scoring_params, "scoring_params")
  allowed <- c("support_buffer_m", "previous_fire_exclusion_buffer_m",
               "previous_fire_cleanup_buffer_m",
               "minimum_remaining_area_m2", "reference_buffer_m")
  .of_check_unknown_keys(scoring_params, allowed, "scoring_params")

  g <- function(key, default) {
    if (is.null(scoring_params[[key]])) default else scoring_params[[key]]
  }
  out <- list(
    support_buffer_m    = g("support_buffer_m", 0),
    erase_mask_buffer_m = g("previous_fire_exclusion_buffer_m", 90),
    erase_post_shave_m  = g("previous_fire_cleanup_buffer_m", 90),
    erase_min_area_m2   = g("minimum_remaining_area_m2", 20000),
    ref_buffer_m        = g("reference_buffer_m", 90)
  )
  .of_check_scalar_num(out$support_buffer_m,
                       "scoring_params$support_buffer_m")
  .of_check_scalar_num(out$erase_mask_buffer_m,
                       "scoring_params$previous_fire_exclusion_buffer_m")
  .of_check_scalar_num(out$erase_post_shave_m,
                       "scoring_params$previous_fire_cleanup_buffer_m")
  .of_check_scalar_num(out$erase_min_area_m2,
                       "scoring_params$minimum_remaining_area_m2")
  .of_check_scalar_num(out$ref_buffer_m,
                       "scoring_params$reference_buffer_m")
  out
}

# Rescue thresholds are NOT public: they are derived automatically from
# the resolved detection by-class vectors using the legacy relaxation
# logic (classes 1,4,5,6,7: seed -20 / delta +10 / floor -20; class 3:
# seed -10 / delta +5 / floor -10; all other classes unchanged).
#' @keywords internal
#' @noRd
.of_derive_rescue_params <- function(detect_resolved) {
  rescue_seed  <- detect_resolved$otsu_min_by_class
  rescue_delta <- detect_resolved$grow_delta_by_class
  rescue_floor <- detect_resolved$min_grow_threshold_by_class

  relax_classes <- c("1", "4", "5", "6", "7")
  rescue_seed[relax_classes]  <- pmax(0, rescue_seed[relax_classes] - 20)
  rescue_delta[relax_classes] <- rescue_delta[relax_classes] + 10
  rescue_floor[relax_classes] <- pmax(0, rescue_floor[relax_classes] - 20)
  rescue_seed["3"]  <- pmax(0, rescue_seed["3"] - 10)
  rescue_delta["3"] <- rescue_delta["3"] + 5
  rescue_floor["3"] <- pmax(0, rescue_floor["3"] - 10)

  list(
    otsu_min_by_class           = rescue_seed,
    grow_delta_by_class         = rescue_delta,
    min_grow_threshold_by_class = rescue_floor
  )
}

# --- internal helpers --------------------------------------------------

.of_normalize_input_spec <- function(x, name, allow_null) {
  if (is.null(x)) {
    if (!allow_null) {
      stop(sprintf("'%s' is required and cannot be NULL.", name), call. = FALSE)
    }
    return(NULL)
  }

  if (inherits(x, "SpatRaster")) {
    return(list(type = "SpatRaster", path = NA_character_, value = x))
  }
  if (inherits(x, "SpatVector")) {
    return(list(type = "SpatVector", path = NA_character_, value = x))
  }
  if (inherits(x, "sf")) {
    return(list(type = "sf", path = NA_character_, value = x))
  }
  if (is.character(x) && length(x) == 1L && nzchar(x)) {
    # Existence is not enforced at config-build time: users may build a config
    # ahead of a batch that materializes inputs later. Downstream stages
    # re-validate file existence before use.
    return(list(
      type = "path",
      path = normalizePath(x, winslash = "/", mustWork = FALSE),
      value = NULL
    ))
  }

  stop(sprintf(
    "'%s' must be a SpatRaster, SpatVector, sf object, or a single character path.",
    name
  ), call. = FALSE)
}

.of_coerce_seed <- function(x, default = 12345L) {
  if (is.null(x)) return(as.integer(default))
  y <- suppressWarnings(as.integer(x)[1L])
  if (is.na(y)) {
    stop("options$deterministic_seed must be coercible to an integer.", call. = FALSE)
  }
  y
}

.of_find_engine_root <- function() {
  candidates <- unique(c(
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd()))
  ))
  for (cand in candidates) {
    cur <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    repeat {
      probe <- file.path(cur, "2_SCRIPTS", "00_FUNCTIONS",
                         "02_DETERMINISTIC_FUNCTIONS")
      if (dir.exists(probe)) return(probe)
      # Also accept being launched from inside 2_SCRIPTS/
      probe2 <- file.path(cur, "00_FUNCTIONS", "02_DETERMINISTIC_FUNCTIONS")
      if (dir.exists(probe2)) return(probe2)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

.of_build_output_routes <- function(output_dir, target_year, run_name) {
  # Base layout: <output_dir>/<year>/DETERMINISTIC/<run_name>/...
  # run_name appears exactly once (no duplicate segment).
  base <- file.path(output_dir, as.character(target_year),
                    "DETERMINISTIC", run_name)
  list(
    base          = base,
    grow_dir      = file.path(base, "01_GROW"),
    refine_dir    = file.path(base, "02_REFINE"),
    scoring_dir   = file.path(base, "05_DECISIONS"),
    validation_dir = file.path(base, "06_VALIDATION"),
    timing_dir    = file.path(base, "99_LOGS_TIMING"),
    # Canonical file names per DETERMINISTIC_OUTPUTS_FINAL.csv.
    grow_vector           = file.path(base, "01_GROW",
                                       paste0("BA_", target_year,
                                              "_OTSUGROW_CORI_ECOREG.shp")),
    refine_merged_gpkg    = file.path(base, "02_REFINE",
                                       paste0("BA_", target_year,
                                              "_REFINE_MERGED.gpkg")),
    internal_decisions    = file.path(base, "05_DECISIONS",
                                       "internal_decisions.gpkg"),
    reference_decisions   = file.path(base, "05_DECISIONS",
                                       "reference_decisions.gpkg"),
    validation_workbook   = file.path(base, "06_VALIDATION",
                                       paste0("validation_ALL_", target_year,
                                              "_", run_name, "_res30.xlsx")),
    timing_csv            = file.path(base, "99_LOGS_TIMING",
                                       paste0(target_year, "_timing_steps.csv"))
  )
}

.of_resolve_registry_path <- function(output_dir, run_name, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = FALSE))
  }
  registry_dir <- toupper(run_name)

  # Run-specific registry convention:
  #   <results_root>/<toupper(run_name)>/RBR_TRAINING_REGISTRY.gpkg
  # If the caller's output_dir is already under a "Results" folder (either
  # ending in .../Results or already containing a .../Results/... segment),
  # we reuse that segment instead of appending a second "Results/".
  out_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  parts <- strsplit(out_norm, "/", fixed = TRUE)[[1L]]

  results_root <- NULL
  res_idx <- which(parts == "Results")
  if (length(res_idx) > 0L) {
    # Use the deepest "Results" segment as the anchor.
    anchor <- tail(res_idx, 1L)
    results_root <- paste(parts[seq_len(anchor)], collapse = "/")
  } else {
    # No "Results" segment — append one under output_dir.
    results_root <- file.path(out_norm, "Results")
  }

  file.path(results_root, registry_dir, "RBR_TRAINING_REGISTRY.gpkg")
}
