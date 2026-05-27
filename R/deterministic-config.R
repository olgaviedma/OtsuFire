#' @title Build a deterministic burned-area mapping configuration object
#'
#' @description
#' Build a complete configuration object for the deterministic burned-area
#' workflow.
#'
#' This function is the main public entry point used to configure
#' deterministic burned-area mapping in OtsuFire. It gathers:
#' \itemize{
#'   \item the spatial inputs,
#'   \item the methodological parameters,
#'   \item the output structure,
#'   \item and the technical runtime options
#' }
#' into a single validated configuration object that can later be passed
#' to:
#' \itemize{
#'   \item \code{\link[=detect_burned_patches]{detect_burned_patches()}},
#'   \item \code{\link[=score_burned_patches]{score_burned_patches()}},
#'   \item \code{\link[=run_deterministic_pipeline]{run_deterministic_pipeline()}}.
#' }
#'
#' The function itself does not execute segmentation, polygonization, or
#' burned-area detection. It performs validation, normalization, and
#' workflow preparation only.
#'
#' The public methodological surface is intentionally restricted to three
#' parameter blocks:
#' \itemize{
#'   \item \code{detect_params},
#'   \item \code{refine_params},
#'   \item \code{scoring_params}.
#' }
#' Each block exposes a small set of user-facing methodological controls.
#' Internal engine parameters remain fixed to validated defaults and are
#' not configurable through the public API.
#'
#' \strong{Deterministic scoring references}
#'
#' The deterministic workflow is fully self-contained and does not depend
#' on external burned-area registries or external burned-like archives
#' during scoring.
#'
#' Spectral-support references used during scoring are obtained either:
#' \itemize{
#'   \item from an explicit \code{keep_pool} supplied to
#'     \code{\link[=score_burned_patches]{score_burned_patches()}} and
#'     built from trusted \emph{external} years,
#'   \item or from a local same-year reference automatically constructed
#'     from high-confidence current-year candidate patches.
#' }
#'
#' Explicit same-year \code{keep_pool} references are not a supported
#' public mode; same-year scoring should use the native local fallback
#' (\code{keep_pool = NULL}). These references are used internally only
#' to compute spectral-support metrics during scoring.
#'
#' @param change_index Raster path, `terra::SpatRaster`, or `NULL`.
#'   Annual change-index raster (e.g. RBR, dNBR, RdNBR). Required before
#'   detection runs but may be `NULL` when the configuration will only be
#'   used with pre-built candidates, for example direct calls to
#'   \code{\link[=score_burned_patches]{score_burned_patches()}}.
#' @param vegetation_map Raster, `sf` object, or path. Optional vegetation
#'   or land-cover stratification layer (CORINE-compatible in OtsuFire
#'   0.2.x).
#' @param burnable_mask Raster, `sf` object, or path. Required binary
#'   burnable-area mask (`1 = burnable`, `0 = non-burnable`).
#' @param hotspots `sf` POINT layer or path. Optional hotspot layer for
#'   the target year.
#' @param previous_year_burned `sf` POLYGON layer or path. Optional
#'   previous-year burned map used for temporal conflict assessment.
#' @param reference_burned_map `sf` POLYGON, raster, or path. Optional
#'   burned reference used by
#'   \code{\link[=validate_fire_maps]{validate_fire_maps()}} when
#'   validation is requested.
#' @param target_year Integer scalar. Target year used for naming outputs
#'   and filtering temporal layers.
#' @param output_dir Character scalar. Root output directory. Defaults to
#'   `tempdir()`.
#' @param run_name Character scalar. Visible experiment identifier used to
#'   name output folders and files. Must be a single non-empty string and
#'   must not contain `/` or `\\`. Defaults to
#'   `"deterministic_burned_map"`.
#' @param detect_params Named list containing detection-stage parameters.
#'   Recognized keys:
#'   \itemize{
#'     \item `seed_threshold`
#'     \item `growth_delta`
#'     \item `minimum_growth_threshold`
#'     \item `minimum_seed_pixels`
#'     \item `seed_threshold_by_vegetation`
#'     \item `growth_delta_by_vegetation`
#'     \item `minimum_growth_threshold_by_vegetation`
#'   }
#'   Unknown keys raise an error.
#' @param refine_params Named list containing refinement-stage parameters.
#'   Recognized keys:
#'   \itemize{
#'     \item `aoi_buffer_m`
#'     \item `omission_buffer_m`
#'     \item `minimum_detected_area_m2`
#'     \item `merge_overlaps`
#'     \item `merge_overlaps_buffer_m`
#'   }
#'   Unknown keys raise an error.
#' @param scoring_params Named list containing scoring-stage parameters.
#'   Recognized keys:
#'   \itemize{
#'     \item `support_buffer_m`
#'     \item `previous_fire_exclusion_buffer_m`
#'     \item `previous_fire_cleanup_buffer_m`
#'     \item `minimum_remaining_area_m2`
#'     \item `reference_buffer_m`
#'   }
#'   Unknown keys raise an error.
#' @param options Named list containing technical/runtime options only
#'   (not methodological controls). Recognized keys include:
#'   \itemize{
#'     \item `deterministic_seed`
#'     \item `engine_root`
#'     \item `whitebox_exe`
#'     \item `gdalwarp_path`
#'     \item `tool_paths`
#'     \item `ecoregion_shapefile_path`
#'     \item `peninsula_shapefile_path`
#'     \item `future_globals_maxsize`
#'     \item `change_index_validation`
#'   }
#'
#'   The deterministic workflow does not consume external burned-area
#'   registries during scoring. Consequently, supplying
#'   `options$registry_path` raises an error.
#'
#'   \strong{\code{options$change_index_validation}}
#'
#'   Optional advanced sub-list controlling the sanity check applied to
#'   `change_index` before detection. Recognized sub-keys:
#'   \itemize{
#'     \item `expected_nodata`
#'     \item `lower_cap`
#'     \item `cap_below`
#'     \item `check_finite`
#'   }
#'
#'   Defaults:
#'   \preformatted{list(
#'   expected_nodata = -9999,
#'   lower_cap = -1000,
#'   cap_below = TRUE,
#'   check_finite = TRUE
#' )}
#'
#'   The validated default `lower_cap = -1000` is appropriate for RBR
#'   products. Other change indices may require different valid ranges.
#'
#' @details
#' \strong{Workflow structure}
#'
#' The deterministic workflow is organised into three methodological
#' stages:
#' \enumerate{
#'   \item detection,
#'   \item refinement,
#'   \item scoring.
#' }
#'
#' Each stage exposes a restricted public parameter surface designed to
#' preserve methodological reproducibility while still allowing
#' ecological adaptation.
#'
#' \strong{Vegetation classes used by \code{*_by_vegetation}}
#' \itemize{
#'   \item \code{1} Agroforestry
#'   \item \code{2} Grassland
#'   \item \code{3} Sparse / burned
#'   \item \code{4} Shrubland
#'   \item \code{5} Broadleaved forest
#'   \item \code{6} Mixed forest
#'   \item \code{7} Conifer forest
#'   \item \code{8} Artificial areas
#'   \item \code{9} Agricultural areas
#'   \item \code{10} Water / wetlands
#'   \item \code{11} Bare
#' }
#'
#' \strong{Detection parameters}
#'
#' \itemize{
#'   \item \code{seed_threshold}: global minimum change-index value
#'   required for a pixel to become a burned seed. Higher values produce
#'   fewer and more conservative seeds; lower values increase sensitivity
#'   but also false-positive risk.
#'
#'   \item \code{seed_threshold_by_vegetation}: class-specific lower
#'   bounds for seed generation. Use when omission or commission errors
#'   are concentrated in specific vegetation types.
#'
#'   \item \code{growth_delta}: global reduction applied from the seed
#'   threshold to the growth threshold. Higher values allow more
#'   permissive expansion and larger patches; lower values produce
#'   tighter and more conservative boundaries.
#'
#'   \item \code{growth_delta_by_vegetation}: class-specific relaxation
#'   applied during region growing. Controls vegetation-specific patch
#'   expansion behaviour.
#'
#'   \item \code{minimum_growth_threshold}: global minimum floor allowed
#'   during region growing. Higher values restrict expansion into weak
#'   spectral change; lower values allow growth into lower-severity
#'   pixels.
#'
#'   \item \code{minimum_growth_threshold_by_vegetation}: class-specific
#'   minimum floor constraining region growing. Useful when over- or
#'   under-expansion is concentrated in particular land-cover classes.
#'
#'   \item \code{minimum_seed_pixels}: minimum number of seed pixels
#'   required to retain a candidate patch. Higher values remove small
#'   isolated detections; lower values increase small-fire sensitivity.
#'   Interpret together with raster resolution.
#' }
#'
#' \strong{Refinement parameters}
#'
#' \itemize{
#'   \item \code{aoi_buffer_m}: buffer around the area of interest used
#'   during refinement. Controls edge artefacts near AOI boundaries.
#'
#'   \item \code{minimum_detected_area_m2}: minimum polygon area retained
#'   after segmentation. Usually kept permissive because later filters
#'   evaluate patch reliability.
#'
#'   \item \code{merge_overlaps}: logical flag controlling whether
#'   overlapping polygons are merged during refinement.
#' }
#'
#' \strong{Scoring parameters}
#'
#' \itemize{
#'   \item \code{support_buffer_m}: buffer used to evaluate contextual
#'   support around candidate patches.
#'
#'   \item \code{previous_fire_exclusion_buffer_m}: buffer around
#'   previous-year burned areas used to identify temporal conflicts.
#'   Increasing this value makes repeated-burn detection more
#'   conservative.
#'
#'   \item \code{previous_fire_cleanup_buffer_m}: buffer used to clean
#'   residual overlap with previous-year burned polygons.
#'
#'   \item \code{minimum_remaining_area_m2}: minimum polygon area
#'   retained after temporal cleanup. Useful for removing residual
#'   slivers after overlap removal.
#'
#'   \item \code{reference_buffer_m}: buffer used when constructing or
#'   comparing spectral-support references.
#' }
#'
#' @return Returns an S3 object of class
#'   `otsufire_burned_mapping_config`.
#'
#'   Stable public fields include:
#'   \itemize{
#'     \item `target_year`
#'     \item `inputs`
#'     \item `run_name`
#'     \item `output_routes`
#'     \item `detect_params`
#'     \item `refine_params`
#'     \item `scoring_params`
#'     \item `rescue_params`
#'     \item `options`
#'     \item `engine_root`
#'     \item `deterministic_seed`
#'     \item `tool_paths`
#'   }
#'
#'   The configuration object stores fully resolved parameter blocks in
#'   their internal engine-ready structure.
#'
#'   Internal implementation details beyond these stable public fields
#'   should not be relied upon by downstream user code.
#'
#' @examples
#' \dontrun{
#' config <- build_burned_mapping_config(
#'   change_index = "MinMin_2022_mosaic_res90m.tif",
#'   vegetation_map = "CLC_2018_peninsula.tif",
#'   burnable_mask = "burnable_mask_binary_corine_2018_ETRS89.tif",
#'   hotspots = "hotspots_2022.gpkg",
#'   target_year = 2022,
#'   output_dir = "results/",
#'   run_name = "balanced_2022"
#' )
#'
#' final_map <- run_deterministic_pipeline(config)
#' }
#'
#' @seealso
#' Related workflow functions:
#' \itemize{
#'   \item \code{\link[=build_supervised_burned_config]{build_supervised_burned_config()}}
#'   \item \code{\link[=change_index_mosaic]{change_index_mosaic()}}
#'   \item \code{\link[=detect_burned_patches]{detect_burned_patches()}}
#'   \item \code{\link[=score_burned_patches]{score_burned_patches()}}
#'   \item \code{\link[=run_deterministic_pipeline]{run_deterministic_pipeline()}}
#'   \item \code{\link[=run_oneyear_supervised_pipeline]{run_oneyear_supervised_pipeline()}}
#' }
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
    stop("'run_name' must not contain '/' or '\\\\' (it is used in output paths).",
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

  if ("registry_path" %in% names(options)) {
    stop(
      "'options$registry_path' is not accepted: the deterministic ",
      "workflow is decoupled from any external burned-like registry and ",
      "no longer reads or resolves a registry path. When an external ",
      "reference is needed for scoring, supply it later as an explicit ",
      "'keep_pool' to score_burned_patches(); there is no implicit ",
      "registry. Remove 'registry_path' from 'options'.",
      call. = FALSE
    )
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

  # The deterministic workflow is fully decoupled from any external
  # burned-like registry: no registry_path is resolved, stored, or read.
  # Keep-pool resolution at scoring time is strictly explicit keep_pool
  # or local current-year fallback.

  cfg <- list(
    target_year = target_year,
    inputs = inputs,
    run_name = run_name,
    output_dir = output_dir,
    output_routes = output_routes,
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
