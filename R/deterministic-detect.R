#' Detect candidate burned patches from an annual change-index raster
#'
#' @description
#' Detect candidate burned patches from an annual change-index raster using a
#' configuration created with [build_burned_mapping_config()].
#'
#' The function estimates Otsu-based thresholds, identifies seed pixels,
#' grows candidate patches, and refines and merges the resulting polygons.
#' These candidates can then be evaluated with [score_burned_patches()] to
#' assign final `"keep"`, `"review"`, or `"drop"` decisions.
#'
#' @param config An object of class `otsufire_burned_mapping_config`, created
#'   with [build_burned_mapping_config()]. Contains the spatial inputs,
#'   detection and refinement settings, output locations, and runtime
#'   options. A change-index raster is required for detection.
#' @param aoi Reserved argument. Must be `NULL`, the default. Supplying
#'   another value raises an error. To restrict detection to a study area,
#'   prepare the spatial inputs for that area before running the function.
#'   See \strong{Spatial extent} below.
#' @param write_outputs Logical scalar. Whether to write detection products to
#'   disk. Default: `TRUE`.
#' @param overwrite Logical scalar. Whether existing output files may be
#'   replaced when `write_outputs = TRUE`. Default: `FALSE`.
#'
#' @section Detection workflow:
#' The function converts a continuous change-index raster into candidate
#' burned patches through five steps:
#' 1. Estimate seed thresholds using Otsu-based thresholding and
#'    vegetation-specific settings where applicable.
#' 2. Identify seed pixels showing strong spectral change.
#' 3. Expand seed-supported patches into neighbouring pixels that meet the
#'    growth criteria.
#' 4. Refine candidate polygons and apply the configured filtering and size
#'    constraints.
#' 5. Merge detections according to the refinement settings.
#'
#' The outputs represent potential burned areas. Their reliability and final
#' decisions are assessed separately by [score_burned_patches()].
#'
#' @section Detection and refinement settings:
#' Detection behaviour is controlled through two parameter lists supplied to
#' [build_burned_mapping_config()]:
#'
#' | Parameter list | Controls |
#' |---|---|
#' | `detect_params` | Seed thresholds, vegetation-specific constraints, region-growing thresholds, and minimum seed size. |
#' | `refine_params` | Refinement buffers, minimum detected polygon area, and overlap merging. |
#'
#' See [build_burned_mapping_config()] for the accepted parameters and their
#' defaults.
#'
#' Change-index thresholds must match the index and its numerical scaling.
#' The ground area represented by a minimum pixel count depends on raster
#' resolution.
#'
#' @section Input validation:
#' Before detection, the function checks the required inputs and evaluates
#' the change-index raster using [clean_raster_inmem()] with
#' `action = "fail"`.
#'
#' The raster checks are controlled by `config$options$change_index_validation`.
#'
#' Detection stops if the raster fails an enabled check. Resolve the reported
#' issue before running the function again.
#'
#' @section Spatial extent:
#' The `aoi` argument is not currently supported. Detection uses the valid
#' extent of the change-index raster, subject to the configured masks and
#' spatial constraints.
#'
#' To restrict the workflow to a study area, crop or mask the relevant inputs
#' before running detection:
#' * the change-index raster;
#' * the vegetation layer, when supplied;
#' * the burnable mask;
#' * supporting study-area and ecoregion layers referenced through
#'   `config$options`.
#'
#' Prepare these layers consistently for the intended study area. If
#' subsequent scoring must use the same spatial restriction, prepare its
#' supporting inputs consistently, including hotspots and previous-year
#' burned-area maps where applicable.
#'
#' Use `aoi = NULL` with the prepared inputs.
#'
#' @return A named list containing detection output paths, diagnostics, and
#'   reserved object fields.
#'
#' | Field | Current contents |
#' |---|---|
#' | `otsu_raster` | Reserved field. Currently `NULL`. |
#' | `seed_raster` | Reserved field. Currently `NULL`. |
#' | `grown_patches` | Reserved field. Currently `NULL`. |
#' | `grown_patches_path` | Path to the grown-patch output. |
#' | `refined_patches` | Reserved field. Currently `NULL`. |
#' | `refined_patches_path` | Path to the refined candidate-patch output. |
#' | `detection_diagnostics` | Diagnostics returned by the detection stage. |
#'
#' The main populated outputs are `grown_patches_path`,
#' `refined_patches_path`, and `detection_diagnostics`. When products are
#' written, use the returned paths to load the spatial layers.
#'
#' The reserved object fields do not currently provide in-memory alternatives
#' to the written products.
#'
#' @seealso [build_burned_mapping_config()], [score_burned_patches()],
#'   [run_deterministic_pipeline()], [clean_raster_inmem()].
#'
#' @examples
#' \dontrun{
#' # Configure detection from an annual RBR raster
#' config <- build_burned_mapping_config(
#'   change_index = "data/RBR_2022.tif",
#'   vegetation_map = "data/vegetation_classes.tif",
#'   burnable_mask = "data/burnable_mask.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022"
#' )
#'
#' # Detect candidate patches and write the outputs
#' detection <- detect_burned_patches(
#'   config = config,
#'   write_outputs = TRUE,
#'   overwrite = FALSE
#' )
#'
#' # Inspect detection diagnostics
#' detection$detection_diagnostics
#'
#' # Load and plot the refined candidate patches
#' refined <- sf::st_read(
#'   detection$refined_patches_path,
#'   quiet = TRUE
#' )
#' plot(sf::st_geometry(refined))
#'
#' # Pass the refined candidates to the scoring stage
#' scored <- score_burned_patches(
#'   burned_candidates = detection$refined_patches_path,
#'   config = config
#' )
#'
#' # Summarise final decisions
#' table(
#'   scored$internal_decisions$class_final,
#'   useNA = "ifany"
#' )
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
