#' Extract predictor features for labelled and scoring polygons
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' Feature-extraction entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.3.x this is a
#' contract-shaped wrapper; the engine execution is performed by
#' [run_oneyear_supervised_pipeline()]. Calling this function standalone
#' validates inputs and surfaces the canonical output paths.
#'
#' @param train_with_folds sf or GPKG path. Labelled training layer with
#'   fold columns.
#' @param scoring_pool sf or GPKG path. Deterministic polygon universe.
#' @param config `otsufire_supervised_burned_config`.
#' @param use_hotspots Logical.
#' @param out_dir Character. Output folder (defaults to config route).
#' @param write_outputs Logical.
#'
#' @keywords internal
#' @noRd
extract_supervised_features <- function(train_with_folds, scoring_pool, config,
                                        use_hotspots = TRUE, out_dir = NULL,
                                        write_outputs = TRUE) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (missing(train_with_folds) || is.null(train_with_folds)) {
    stop("'train_with_folds' is required.", call. = FALSE)
  }
  if (missing(scoring_pool) || is.null(scoring_pool)) {
    stop("'scoring_pool' is required.", call. = FALSE)
  }
  if (!is.logical(use_hotspots) || length(use_hotspots) != 1L) {
    stop("'use_hotspots' must be TRUE or FALSE.", call. = FALSE)
  }
  .of_check_input_file(config$inputs$change_index, "change_index")

  if (is.null(out_dir)) out_dir <- config$output_routes$features_dir
  list(
    train_features          = NULL,
    scoring_features        = NULL,
    features_geometry_gpkg  = config$output_routes$features_geometry_gpkg,
    train_features_rds      = file.path(out_dir, "train_features.rds"),
    scoring_features_rds    = file.path(out_dir, "scoring_features.rds"),
    note = paste0(
      "In 0.2.x features are materialized by run_oneyear_supervised_pipeline(); ",
      "calling extract_supervised_features() standalone validates inputs."
    )
  )
}
