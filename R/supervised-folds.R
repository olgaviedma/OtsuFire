#' Create spatial folds for out-of-fold diagnostics
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' Fold-construction entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.3.x this is a
#' contract-shaped wrapper; the engine execution is performed by
#' [run_oneyear_supervised_pipeline()]. Calling this function standalone
#' validates inputs and surfaces the canonical output paths resolved from
#' `config$output_routes`.
#'
#' @param train_labelled sf or GPKG path. Labelled training layer.
#' @param split_unit Character scalar. `"fire"` (default) or `"poly"`.
#' @param block_sizes_m Numeric vector. Candidate block sizes in meters.
#' @param k_candidates Integer vector. Candidate fold counts.
#' @param out_dir Character. Output folder (defaults to config route).
#' @param config Optional `otsufire_supervised_burned_config`; required
#'   when `out_dir` is not supplied.
#' @param write_outputs Logical.
#' @param overwrite Logical.
#'
#' @keywords internal
#' @noRd
make_spatial_folds <- function(train_labelled,
                               split_unit = c("fire", "poly"),
                               block_sizes_m = c(5000, 3000, 2000),
                               k_candidates = c(5L, 4L, 3L),
                               out_dir = NULL,
                               config = NULL,
                               write_outputs = TRUE,
                               overwrite = FALSE) {
  split_unit <- match.arg(split_unit)
  if (missing(train_labelled) || is.null(train_labelled)) {
    stop("'train_labelled' is required.", call. = FALSE)
  }
  if (!(inherits(train_labelled, c("sf", "SpatVector")) ||
        (is.character(train_labelled) && length(train_labelled) == 1L &&
         file.exists(train_labelled)))) {
    stop("'train_labelled' must be sf, SpatVector, or an existing path.",
         call. = FALSE)
  }
  if (!is.numeric(block_sizes_m) || length(block_sizes_m) < 1L ||
      any(block_sizes_m <= 0)) {
    stop("'block_sizes_m' must be a positive numeric vector.", call. = FALSE)
  }
  if (!is.numeric(k_candidates) || length(k_candidates) < 1L ||
      any(k_candidates < 2L)) {
    stop("'k_candidates' must be an integer-like vector with values >= 2.",
         call. = FALSE)
  }
  if (is.null(out_dir)) {
    if (is.null(config) ||
        !inherits(config, "otsufire_supervised_burned_config")) {
      stop("Supply either 'out_dir' or a valid 'config'.", call. = FALSE)
    }
    out_dir <- config$output_routes$folds_dir
  }

  base_dir <- out_dir
  year <- if (!is.null(config)) config$target_year else NA_integer_
  list(
    blocks            = NULL,
    folds_table       = NULL,
    train_with_folds  = NULL,
    blocks_gpkg       = if (!is.na(year))
      file.path(base_dir, sprintf("%d_blocks_5000m.gpkg", year)) else NA_character_,
    folds_csv         = if (!is.na(year))
      file.path(base_dir, sprintf("%d_folds_5000m_k5_r1_unit-%s.csv",
                                   year, split_unit)) else NA_character_,
    train_with_folds_gpkg = if (!is.null(config))
      config$output_routes$train_with_folds_gpkg else NA_character_,
    note = paste0(
      "In 0.2.x spatial folds are materialized by ",
      "run_oneyear_supervised_pipeline(); calling make_spatial_folds() ",
      "standalone validates inputs and returns canonical paths."
    )
  )
}
