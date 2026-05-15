#' Run the full deterministic burned-area pipeline for one year
#'
#' @description
#' High-level public wrapper that chains
#' [detect_burned_patches()], [score_burned_patches()], and — when
#' requested — the shared [validate_fire_maps()] utility, using a single
#' [otsufire_burned_mapping_config][build_burned_mapping_config].
#'
#' This orchestrator is new in 0.2.x and is **not** a thin copy of the
#' 0.1.x `8_deterministic_bridge.R`. It consumes the new public config
#' object, calls the new modular public functions, and assembles a
#' contract-shaped return value plus a timing log.
#'
#' Contract source: `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `DETERMINISTIC_OUTPUTS_FINAL.csv`.
#'
#' @param config An `otsufire_burned_mapping_config` object from
#'   [build_burned_mapping_config()].
#' @param write_outputs Logical scalar. Whether the pipeline writes its
#'   outputs to disk. Default `TRUE`.
#' @param overwrite Logical scalar. Whether existing outputs from the same
#'   run may be replaced. Default `FALSE`.
#' @param run_validation Logical scalar or `"auto"`. Controls validation.
#'   `TRUE` always attempts validation, `FALSE` skips, `"auto"` (default)
#'   runs only when `config$inputs$reference_burned_map` is non-NULL.
#'
#' @return A named list of class `otsufire_deterministic_run` with fields:
#'   `detection`, `scoring`, `validation`, `result_paths`, `timing_log`,
#'   and `config`.
#'
#' @family workflow
#' @export
run_deterministic_pipeline <- function(config, write_outputs = TRUE,
                                       overwrite = FALSE,
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
    keep_pool = NULL,
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
