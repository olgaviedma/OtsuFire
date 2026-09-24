#' @title Detect candidate burned patches from an annual change-index raster
#'
#' @description
#' Run the detection stage of the deterministic OtsuFire workflow on an
#' annual change-index raster.
#'
#' This function takes an
#' [otsufire_burned_mapping_config][build_burned_mapping_config] object
#' and produces the candidate burned-patch outputs later used by
#' [score_burned_patches()].
#' It is the "find candidate fires" step, not the final burned-map step.
#'
#' Conceptually, this stage includes:
#' \enumerate{
#'   \item Otsu-based threshold estimation,
#'   \item seed generation,
#'   \item seed-supported region growing,
#'   \item segmentation refinement,
#'   \item candidate merging.
#' }
#'
#' In practice, the wrapper validates inputs, prepares the runtime
#' context, rejects unsupported `aoi` input explicitly, and then calls
#' the internal detection engine.
#'
#' Given the same inputs and parameter settings, this stage is
#' deterministic and reproducible.
#'
#' @param config An `otsufire_burned_mapping_config` object generated with
#'   [build_burned_mapping_config()]. This object contains validated
#'   workflow inputs, methodological parameters, output routes, and
#'   runtime options.
#' @param aoi NOT YET IMPLEMENTED — must be `NULL`. The parameter is
#'   reserved for a future release in which detection can be spatially
#'   restricted to an AOI. The current detection engine does not propagate
#'   the AOI to the underlying grow/refine stages, so a non-`NULL` value
#'   is rejected with an explicit error rather than silently ignored.
#'   Pass `aoi = NULL` (the default) to run detection over the full valid
#'   extent of the change-index raster.
#' @param write_outputs Logical scalar. Whether the detection products
#'   should be written to disk. Defaults to `TRUE` because this stage is
#'   commonly inspected visually during workflow development and QA/QC.
#' @param overwrite Logical scalar. Whether existing outputs may be
#'   replaced when `write_outputs = TRUE`. Defaults to `FALSE`.
#'
#' @details
#' \strong{Workflow overview}
#'
#' The detection stage transforms a continuous annual change-index raster
#' into spatially coherent candidate burned patches.
#'
#' The workflow proceeds conceptually as follows:
#' \enumerate{
#'   \item estimate vegetation-aware seed thresholds from the
#'     change-index distribution,
#'   \item identify high-confidence burned seed pixels,
#'   \item expand seeds into neighbouring burned-like pixels using
#'     constrained region growing,
#'   \item refine and filter candidate polygons,
#'   \item merge overlapping or fragmented detections where appropriate.
#' }
#'
#' The resulting outputs are candidate burned patches, not final burned
#' decisions. Reliability assessment and final `keep` / `review` /
#' `drop` assignment happen later in [score_burned_patches()].
#'
#' The detection stage operates exclusively on the deterministic
#' candidate-generation workflow and does not perform probabilistic
#' classification or supervised scoring.
#'
#' Detection behaviour is controlled through the parameter blocks stored
#' inside `config`, particularly:
#' \itemize{
#'   \item `detect_params`
#'   \item `refine_params`
#' }
#'
#' Key methodological controls include:
#' \itemize{
#'   \item seed-generation thresholds,
#'   \item vegetation-specific constraints,
#'   \item region-growing permissiveness,
#'   \item minimum patch size,
#'   \item overlap-merging behaviour.
#' }
#'
#' Before delegating to the internal engine, the function also validates
#' that the required inputs exist on disk and runs the configured
#' `change_index` integrity check (the same rules as [clean_raster_file()]),
#' stopping with an error if the raster is not clean.
#'
#' \strong{Current AOI limitation and what to do instead}
#'
#' At the moment, `aoi` is intentionally rejected rather than partially
#' applied. The current detection engine still runs through the full grow
#' and refine machinery without a true AOI-aware execution path, so
#' allowing `aoi` now would be misleading.
#'
#' If you need an AOI-limited run today, you must guarantee that the
#' workflow inputs have already been prepared externally for that AOI
#' before calling `detect_burned_patches()`. In practice, this means that
#' the files on disk should already be cropped or masked to the intended
#' AOI for at least:
#' \itemize{
#'   \item `change_index`,
#'   \item `vegetation_map`,
#'   \item `burnable_mask`,
#'   \item and any supporting spatial layers used later in the workflow,
#'     especially the study-area boundary and ecoregion layers referenced
#'     through `config$options`.
#' }
#'
#' If later deterministic stages should also remain AOI-limited, the same
#' external preparation principle should be applied consistently to other
#' relevant spatial inputs, such as hotspot layers or previous-year
#' burned maps.
#'
#' \strong{Why this stage matters}
#'
#' Many burned-area mapping errors start during candidate generation,
#' before scoring ever begins.
#'
#' The deterministic detection stage is designed to:
#' \itemize{
#'   \item maximise spatial coherence,
#'   \item reduce obvious false positives early,
#'   \item preserve low-severity burned edges where possible,
#'   \item produce a traceable candidate universe for later scoring.
#' }
#'
#' Keeping candidate generation separate from scoring makes the workflow
#' easier to inspect, explain, and reproduce.
#'
#' @return Returns a named list containing the main detection products and
#'   their associated file paths.
#'
#'   Stable public fields currently include:
#'   \itemize{
#'     \item `otsu_raster`
#'     \item `seed_raster`
#'     \item `grown_patches`
#'     \item `grown_patches_path`
#'     \item `refined_patches`
#'     \item `refined_patches_path`
#'     \item `detection_diagnostics`
#'   }
#'
#'   In the current public wrapper, `grown_patches_path`,
#'   `refined_patches_path`, and `detection_diagnostics` are the main
#'   populated outputs. The object placeholders `otsu_raster`,
#'   `seed_raster`, `grown_patches`, and `refined_patches` are kept for
#'   compatibility and currently return `NULL`.
#'
#'   The exact internal implementation should not be relied upon beyond
#'   these stable public outputs.
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
#' detection <- detect_burned_patches(config)
#'
#' detection$refined_patches_path
#' refined <- sf::st_read(detection$refined_patches_path, quiet = TRUE)
#' plot(sf::st_geometry(refined))
#' }
#'
#' @seealso
#' Related deterministic workflow functions:
#' \itemize{
#'   \item [build_burned_mapping_config()]
#'   \item [score_burned_patches()]
#'   \item [run_deterministic_pipeline()]
#'   \item [validate_fire_maps()]
#'   \item [build_supervised_burned_config()]
#' }
#'
#' @family modular
#' @export
detect_burned_patches <- function(config, aoi = NULL, write_outputs = TRUE,
                                  overwrite = FALSE) {

  if (!inherits(config, "otsufire_burned_mapping_config")) {
    stop("'config' must be created by build_burned_mapping_config().",
         call. = FALSE)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  # Contract requirement: change_index and burnable_mask must be present
  # and resolvable before detection runs.
  if (is.null(config$inputs$change_index)) {
    stop(
      "detect_burned_patches(): config$inputs$change_index is NULL.\n",
      "Rebuild the config with an explicit 'change_index' argument.",
      call. = FALSE
    )
  }
  if (is.null(config$inputs$burnable_mask)) {
    stop(
      "detect_burned_patches(): config$inputs$burnable_mask is NULL.\n",
      "Rebuild the config with an explicit 'burnable_mask' argument.",
      call. = FALSE
    )
  }

  # Validate that on-disk inputs exist (paths were stored at config-build
  # time without existence enforcement).
  .of_check_input_file(config$inputs$change_index, "change_index")
  .of_check_input_file(config$inputs$burnable_mask, "burnable_mask")
  if (!is.null(config$inputs$vegetation_map)) {
    .of_check_input_file(config$inputs$vegetation_map, "vegetation_map")
  }

  # P1-DET-01: `aoi` is documented but the engine adapters
  # (.of_run_detection -> grow_args / refine_args) do not propagate it to
  # the underlying stages. Until that wiring exists we reject any
  # non-NULL value explicitly so users do not silently believe they ran
  # a regional restriction.
  if (!is.null(aoi)) {
    stop(
      "aoi parameter is documented but not yet implemented in the ",
      "current detection engine; pass aoi=NULL.",
      call. = FALSE
    )
  }

  # Sanity-check ONLY the change_index raster in memory before detection.
  # Never rewrites disk, never cleans silently: a dirty raster aborts with
  # a clear error (clean_raster_inmem(action = "fail")). Tunable via
  # config$options$change_index_validation (see build_burned_mapping_config()).
  .of_validate_change_index(config)

  res <- .of_run_detection(config, aoi = aoi,
                           write_outputs = isTRUE(write_outputs),
                           overwrite = isTRUE(overwrite))

  # Contract-shaped return. Diagnostic rasters (otsu/seed) are only
  # written when the engine's save_debug_rasters flag is on; here we
  # surface the paths when the engine declared them.
  list(
    otsu_raster           = NULL,
    seed_raster           = NULL,
    grown_patches         = NULL,
    grown_patches_path    = res$grown_patches_path %||% NA_character_,
    refined_patches       = NULL,
    refined_patches_path  = res$refined_patches_path %||% NA_character_,
    detection_diagnostics = res$detection_diagnostics
  )
}

# --- internal validation helpers ---------------------------------------

#' @keywords internal
#' @noRd
.of_check_input_file <- function(spec, label) {
  if (is.null(spec)) return(invisible(NULL))
  if (!is.null(spec$path) && !is.na(spec$path) && nzchar(spec$path) &&
      spec$type == "path") {
    if (!file.exists(spec$path)) {
      stop(sprintf("%s: file does not exist: %s", label, spec$path),
           call. = FALSE)
    }
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.of_resolve_change_index_validation <- function(config) {
  defaults <- .of_change_index_validation_defaults()
  blk <- config$options[["change_index_validation"]]
  if (is.null(blk)) return(defaults)
  # Structural validity (named list, recognized keys) was already enforced
  # by build_burned_mapping_config(); merge user values over the defaults.
  utils::modifyList(defaults, blk)
}

#' @keywords internal
#' @noRd
.of_validate_change_index <- function(config) {
  spec <- config$inputs$change_index
  ci_rast <-
    if (!is.null(spec$value) && inherits(spec$value, "SpatRaster")) {
      spec$value
    } else if (!is.null(spec$path) && !is.na(spec$path) &&
               nzchar(spec$path)) {
      terra::rast(spec$path)
    } else {
      stop("detect_burned_patches(): change_index could not be read for ",
           "validation.", call. = FALSE)
    }

  v <- .of_resolve_change_index_validation(config)

  # action = "fail": clean raster -> returns silently; dirty raster ->
  # stop() with a message that includes the lower_cap value actually used.
  clean_raster_inmem(
    ci_rast,
    name            = "change_index",
    expected_nodata = v$expected_nodata,
    lower_cap       = v$lower_cap,
    cap_below       = v$cap_below,
    check_finite    = v$check_finite,
    action          = "fail",
    verbose         = TRUE
  )
  invisible(NULL)
}

#' @keywords internal
#' @noRd
# P1-DET-01: kept only for backward namespace compatibility. After the
# explicit aoi-rejection in detect_burned_patches() this helper is
# unreachable from the public surface; do not call it from new code.
.of_validate_aoi <- function(aoi) {
  if (inherits(aoi, c("sf", "SpatVector"))) return(aoi)
  if (is.character(aoi) && length(aoi) == 1L && file.exists(aoi)) return(aoi)
  stop("'aoi' must be an sf, SpatVector, or a readable file path.",
       call. = FALSE)
}
