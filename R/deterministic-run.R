#' Run the full deterministic burned-area pipeline for one year
#'
#' @description
#' High-level public wrapper that orchestrates the complete deterministic
#' burned-area workflow for a single target year by sequentially chaining
#' \code{\link[=detect_burned_patches]{detect_burned_patches()}},
#' \code{\link[=score_burned_patches]{score_burned_patches()}}, and — when
#' requested — the shared
#' \code{\link[=validate_fire_maps]{validate_fire_maps()}} utility, using
#' a single
#' \code{\link[=build_burned_mapping_config]{otsufire_burned_mapping_config}}.
#'
#' The pipeline implements the full unsupervised OtsuFire decision
#' workflow, including candidate patch delineation through adaptive
#' seed-and-grow segmentation, sequential rule-based filtering,
#' temporal-conflict assessment, spectral-support scoring, and final
#' deterministic confidence assignment into `keep`, `review`, or `drop`.
#'
#' Deterministic scoring within the pipeline can reuse an explicit
#' spectral-support \code{keep_pool} supplied by the caller from trusted
#' external years, or fall back to the same-year local reference
#' automatically when \code{keep_pool} is omitted. The deterministic
#' pipeline does not read or depend on implicit external burned-like
#' registries.
#'
#' This orchestrator is new in 0.2.x and is \strong{not} a thin copy of
#' the legacy 0.1.x `8_deterministic_bridge.R`. It consumes the new
#' public configuration object, calls the modular public functions, and
#' assembles a contract-shaped return object together with timing
#' diagnostics and output routes.
#'
#' Contract source: `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `DETERMINISTIC_OUTPUTS_FINAL.csv`.
#'
#' @details
#' \strong{Deterministic workflow structure}
#'
#' The pipeline executes the deterministic burned-area workflow in three
#' main stages:
#'
#' \enumerate{
#'   \item \strong{Detection stage.}
#'     \code{\link[=detect_burned_patches]{detect_burned_patches()}}
#'     applies adaptive Otsu thresholding and seed-and-grow segmentation
#'     to the annual change-index raster in order to generate candidate
#'     burned patches.
#'   \item \strong{Scoring stage.}
#'     \code{\link[=score_burned_patches]{score_burned_patches()}}
#'     evaluates candidate patches through a sequential rule-based system
#'     including:
#'     \itemize{
#'       \item Filter 1: seed support and burnable-context plausibility;
#'       \item Filter 2: optional active-fire hotspot corroboration;
#'       \item Temporal consistency assessment against previous-year
#'         burned areas (a distinct step, not a numbered filter);
#'       \item Filter 3: spectral-support evaluation relative to a
#'         high-confidence keep-like reference distribution.
#'     }
#'     The outputs of these components are integrated into a final
#'     deterministic confidence decision assigning each patch to `keep`,
#'     `review`, or `drop`, while preserving full decision traceability.
#'
#'     Within the current wrapper, this scoring stage can use either an
#'     explicit \code{keep_pool} built from trusted external years or the
#'     local same-year spectral-support fallback when \code{keep_pool} is
#'     omitted. Explicit same-year \code{keep_pool} references are not a
#'     supported public mode.
#'   \item \strong{Validation stage.} When validation is enabled and a
#'     reference burned-area layer is available,
#'     \code{\link[=validate_fire_maps]{validate_fire_maps()}} evaluates
#'     the resulting deterministic burned-area outputs against external
#'     burned-area references.
#' }
#'
#' @param config An `otsufire_burned_mapping_config` object returned by
#'   \code{\link[=build_burned_mapping_config]{build_burned_mapping_config()}}.
#' @param write_outputs Logical scalar. Whether the pipeline writes its
#'   outputs to disk. Default `TRUE`.
#' @param overwrite Logical scalar. Whether existing outputs from the same
#'   run may be replaced. Default `FALSE`.
#' @param keep_pool Optional explicit keep-like reference pool passed
#'   through to
#'   \code{\link[=score_burned_patches]{score_burned_patches()}}. This
#'   should be a pre-built keep-pool summary object such as the output of
#'   \code{\link[=build_keep_pool_from_samples]{build_keep_pool_from_samples()}}.
#'   When `NULL` (default), the pipeline uses the deterministic scoring
#'   local fallback for the target year. Explicit pools that include the
#'   same target year are intentionally rejected when helper metadata
#'   reveal that overlap; for same-year scoring, use `keep_pool = NULL`.
#' @param run_validation Logical scalar or `"auto"`. Controls validation
#'   behaviour.
#'   \itemize{
#'     \item `TRUE`: always attempts validation;
#'     \item `FALSE`: skips validation;
#'     \item `"auto"` (default): runs validation only when
#'       `config$inputs$reference_burned_map` is non-NULL.
#'   }
#'
#' @return A named list of class `otsufire_deterministic_run` with fields:
#' \describe{
#'   \item{`detection`}{Output object returned by
#'     \code{\link[=detect_burned_patches]{detect_burned_patches()}}.}
#'   \item{`scoring`}{Output object returned by
#'     \code{\link[=score_burned_patches]{score_burned_patches()}}.}
#'   \item{`validation`}{Validation outputs returned by
#'     \code{\link[=validate_fire_maps]{validate_fire_maps()}} when
#'     validation is executed; otherwise `NULL`.}
#'   \item{`result_paths`}{Named list of important written output paths
#'     produced during the workflow.}
#'   \item{`timing_log`}{Timing diagnostics summarising execution time
#'     for each deterministic stage.}
#'   \item{`config`}{The input configuration object used to run the
#'     workflow.}
#' }
#'
#' @family workflow
#' @export
run_deterministic_pipeline <- function(config, write_outputs = TRUE,
                                       overwrite = FALSE,
                                       keep_pool = NULL,
                                       run_validation = "auto") {

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
  keep_pool <- .of_normalize_keep_pool_arg(keep_pool)
  keep_pool <- .of_validate_keep_pool_target_year(
    keep_pool,
    config$target_year
  )
  if (!((is.logical(run_validation) && length(run_validation) == 1L) ||
        (is.character(run_validation) && length(run_validation) == 1L &&
         identical(run_validation, "auto")))) {
    stop("'run_validation' must be TRUE, FALSE, or 'auto'.", call. = FALSE)
  }

  do_validate <- .of_resolve_run_validation(run_validation, config)

  timing <- list()
  record_step <- function(name, started) {
    timing[[name]] <<- list(
      step = name,
      elapsed_sec = as.numeric(Sys.time() - started, units = "secs")
    )
  }

  # --- detection ------------------------------------------------------
  t0 <- Sys.time()
  detection <- detect_burned_patches(
    config = config,
    aoi = NULL,
    write_outputs = write_outputs,
    overwrite = overwrite
  )
  record_step("detection", t0)

  # Candidates to feed the scoring stage. Prefer refined_patches, fall back
  # to grown_patches.
  candidates_path <- detection$refined_patches_path
  if (is.na(candidates_path) || !nzchar(candidates_path) ||
      !file.exists(candidates_path)) {
    candidates_path <- detection$grown_patches_path
  }
  if (is.na(candidates_path) || !nzchar(candidates_path) ||
      !file.exists(candidates_path)) {
    stop(
      "run_deterministic_pipeline(): detection stage did not produce a usable ",
      "candidate layer. Inspect the detection_diagnostics for details.",
      call. = FALSE
    )
  }

  # Thread detection output paths into the scoring call so stage 4 can
  # locate the stage-1 (grown) layer without rediscovery heuristics.
  config_scoring <- config
  config_scoring$options$.detect_paths <- list(
    grown_patches_path   = detection$grown_patches_path,
    refined_patches_path = detection$refined_patches_path
  )

  # --- scoring --------------------------------------------------------
  t0 <- Sys.time()
  scoring <- score_burned_patches(
    burned_candidates = candidates_path,
    config = config_scoring,
    keep_pool = keep_pool,
    write_outputs = write_outputs,
    overwrite = overwrite
  )
  record_step("scoring", t0)

  # --- validation (optional, shared) ----------------------------------
  validation <- NULL
  if (isTRUE(do_validate)) {
    t0 <- Sys.time()
    validation <- tryCatch(
      .of_run_shared_validation(scoring, config, write_outputs = write_outputs),
      error = function(e) {
        warning("Shared validation failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    record_step("validation", t0)
  }

  # --- timing log -----------------------------------------------------
  timing_df <- .of_timing_to_df(timing)
  timing_csv <- NA_character_
  if (isTRUE(write_outputs) && !is.null(timing_df) && nrow(timing_df) > 0L) {
    timing_csv <- config$output_routes$timing_csv
    dir.create(dirname(timing_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(timing_df, timing_csv, row.names = FALSE)
  }

  result_paths <- list(
    grown_patches       = detection$grown_patches_path,
    refined_patches     = detection$refined_patches_path,
    internal_decisions  = scoring$internal_decisions_path,
    reference_decisions = scoring$reference_decisions_path,
    validation_workbook = if (!is.null(validation))
      validation$workbook_path %||% NA_character_ else NA_character_,
    timing_csv          = timing_csv
  )

  out <- list(
    config       = config,
    detection    = detection,
    scoring      = scoring,
    validation   = validation,
    result_paths = result_paths,
    timing_log   = timing_df
  )
  class(out) <- c("otsufire_deterministic_run", "list")
  out
}

# --- internal helpers --------------------------------------------------

#' @keywords internal
#' @noRd
.of_resolve_run_validation <- function(run_validation, config) {
  if (is.logical(run_validation)) return(isTRUE(run_validation))
  # "auto": run only if a reference burned map was provided.
  !is.null(config$inputs$reference_burned_map)
}

#' @keywords internal
#' @noRd
.of_run_shared_validation <- function(scoring, config, write_outputs) {
  ref_spec <- config$inputs$reference_burned_map
  if (is.null(ref_spec)) return(NULL)

  burnable_spec <- config$inputs$burnable_mask

  # The shared validator expects paths. If the scoring produced an on-disk
  # internal_decisions layer, use that directly. Otherwise we can't wire
  # the shared validator without materializing a temporary file.
  input_shp <- scoring$internal_decisions_path
  if (is.na(input_shp) || !nzchar(input_shp) || !file.exists(input_shp)) {
    return(NULL)
  }

  validation_dir <- config$output_routes$validation_dir
  if (isTRUE(write_outputs)) {
    dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ref_path <- .of_input_to_path(ref_spec, "reference_burned_map")
  burnable_path <- .of_input_to_path(burnable_spec, "burnable_mask")
  mask_path <- .of_input_to_path(burnable_spec, "burnable_mask")

  res <- validate_fire_maps(
    input_shapefile = input_shp,
    ref_shapefile   = ref_path,
    mask_shapefile  = mask_path,
    burnable_raster = burnable_path,
    year_target     = config$target_year,
    validation_dir  = dirname(validation_dir)
  )

  list(
    metrics         = res$metrics,
    polygon_summary = res$polygon_summary,
    workbook_path   = NA_character_
  )
}

#' @keywords internal
#' @noRd
.of_timing_to_df <- function(timing) {
  if (length(timing) == 0L) return(NULL)
  do.call(rbind, lapply(timing, function(r)
    data.frame(step = r$step, elapsed_sec = r$elapsed_sec,
               stringsAsFactors = FALSE)
  ))
}

#' @export
print.otsufire_deterministic_run <- function(x, ...) {
  cat("<otsufire_deterministic_run>\n")
  cat("  run_name   :", x$config$run_name, "\n")
  cat("  target_year:", x$config$target_year, "\n")
  cat("  result_paths:\n")
  for (nm in names(x$result_paths)) {
    v <- x$result_paths[[nm]]
    cat(sprintf("    - %-20s: %s\n", nm,
                if (is.null(v) || is.na(v) || !nzchar(v)) "<none>" else v))
  }
  invisible(x)
}
