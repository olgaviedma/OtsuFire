#' Assign deterministic filters and final decisions to candidate burned patches
#'
#' @description
#' Public modular scoring function. Receives candidate burned patches
#' produced by [detect_burned_patches()] and applies the canonical
#' deterministic decision logic: change-index scoring against the keep-like
#' reference distribution, support-based filters (filter_1, filter_2,
#' filter_3), temporal overlap assessment against `previous_year_burned`,
#' and final class integration (`class_final`).
#'
#' Contract source: `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `DETERMINISTIC_OUTPUTS_FINAL.csv`, `DETERMINISTIC_DATA_MODEL_FINAL.csv`.
#'
#' Validation of the decision layer against an external reference is NOT
#' part of this function; users call [validate_fire_maps()] on the
#' resulting `internal_decisions` (or supervised thresholded outputs).
#' There is no `validate_burned_maps()` wrapper in 0.2.x.
#'
#' @param burned_candidates sf POLYGON, SpatVector, or path to a readable
#'   vector file. Candidate burned patches from
#'   [detect_burned_patches()].
#' @param config An `otsufire_burned_mapping_config` object from
#'   [build_burned_mapping_config()].
#' @param keep_pool Optional keep-like reference pool (sf POLYGON or a
#'   pre-built keep-pool summary list). When `NULL`, the scoring engine
#'   builds the pool from the burned-like registry or the local fallback.
#' @param write_outputs Logical scalar. Whether to write scoring products.
#'   Default `TRUE`.
#' @param overwrite Logical scalar. Whether existing scoring outputs may be
#'   replaced. Default `FALSE`.
#'
#' @return A named list with fields: `internal_decisions`,
#'   `internal_decisions_path`, `reference_decisions`,
#'   `reference_decisions_path`, `keep_pool_summary`, `scoring_diagnostics`.
#'
#' @family modular
#' @export
score_burned_patches <- function(burned_candidates, config, keep_pool = NULL,
                                 write_outputs = TRUE, overwrite = FALSE) {

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

  if (missing(burned_candidates) || is.null(burned_candidates)) {
    stop("'burned_candidates' is required.", call. = FALSE)
  }

  # Validate config content (change_index) before touching burned_candidates
  # on disk — if the config itself is not usable, we fail fast with a
  # config-level error instead of a misleading path-existence error.
  if (is.null(config$inputs$change_index)) {
    stop(
      "score_burned_patches(): config$inputs$change_index is NULL.\n",
      "Rebuild the config with an explicit 'change_index' argument.",
      call. = FALSE
    )
  }
  .of_check_input_file(config$inputs$change_index, "change_index")

  if (is.character(burned_candidates) && length(burned_candidates) == 1L) {
    if (!file.exists(burned_candidates)) {
      stop("'burned_candidates' path does not exist: ", burned_candidates,
           call. = FALSE)
    }
  } else if (!inherits(burned_candidates, c("sf", "SpatVector"))) {
    stop(
      "'burned_candidates' must be an sf POLYGON, SpatVector, or a readable path.",
      call. = FALSE
    )
  }

  if (!is.null(keep_pool) &&
      !(inherits(keep_pool, c("sf", "SpatVector")) || is.list(keep_pool))) {
    stop("'keep_pool' must be NULL, an sf / SpatVector, or a keep-pool list.",
         call. = FALSE)
  }

  res <- .of_run_scoring(burned_candidates, config, keep_pool = keep_pool,
                         write_outputs = isTRUE(write_outputs),
                         overwrite = isTRUE(overwrite))

  list(
    internal_decisions       = res$internal_decisions,
    internal_decisions_path  = res$internal_decisions_path %||% NA_character_,
    reference_decisions      = res$reference_decisions,
    reference_decisions_path = res$reference_decisions_path %||% NA_character_,
    keep_pool_summary        = res$keep_pool_summary,
    scoring_diagnostics      = res$scoring_diagnostics
  )
}
