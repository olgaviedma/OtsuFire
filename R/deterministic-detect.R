#' Detect candidate burned patches from the annual change index
#'
#' @description
#' Public modular detection function. Consumes an
#' [otsufire_burned_mapping_config][build_burned_mapping_config]
#' and produces the intermediate and final detection products needed by
#' [score_burned_patches()]. Conceptually this stage covers
#' Otsu thresholding, seed-supported growth, segmentation refinement, and
#' candidate merging.
#'
#' This function validates inputs and delegates the heavy numerical work
#' to the validated internal engine (not exported). Its public signature
#' and return shape are frozen by
#' `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv` and
#' `DETERMINISTIC_OUTPUTS_FINAL.csv`.
#'
#' @param config An `otsufire_burned_mapping_config` object from
#'   [build_burned_mapping_config()].
#' @param aoi Optional sf POLYGON, SpatVector, or path. Restricts detection
#'   to the supplied geometry. When `NULL`, detection runs over the full
#'   valid extent of the change-index raster.
#' @param write_outputs Logical scalar. Whether to write detection products
#'   to disk. Default `TRUE` because running this stage alone is usually
#'   for inspection.
#' @param overwrite Logical scalar. Whether existing outputs may be
#'   replaced when `write_outputs = TRUE`. Default `FALSE`.
#'
#' @return A named list with fields:
#'   `otsu_raster`, `seed_raster`, `grown_patches`, `refined_patches`,
#'   `detection_diagnostics`, and the corresponding `*_path` entries.
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

  if (!is.null(aoi)) {
    aoi <- .of_validate_aoi(aoi)
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
.of_validate_aoi <- function(aoi) {
  if (inherits(aoi, c("sf", "SpatVector"))) return(aoi)
  if (is.character(aoi) && length(aoi) == 1L && file.exists(aoi)) return(aoi)
  stop("'aoi' must be an sf, SpatVector, or a readable file path.",
       call. = FALSE)
}
