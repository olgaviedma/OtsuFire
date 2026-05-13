#' Score the deterministic polygon universe and export supervised final maps
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' Scoring entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.3.x this is a
#' contract-shaped wrapper; real execution happens inside
#' [run_oneyear_supervised_pipeline()]. Thresholded `p_burned` or
#' `p_burned_current_year` outputs may then be passed to
#' [validate_fire_maps()] for external validation.
#'
#' @param scoring_features sf or GPKG path.
#' @param model Fitted XGBoost model or RDS path.
#' @param recipe Named list or RDS path.
#' @param oof_summary Optional sf or GPKG path.
#' @param export_burned_like Logical.
#' @param out_map_dir Output folder; defaults to config route.
#' @param config Optional `otsufire_supervised_burned_config`.
#'
#' @keywords internal
#' @noRd
score_supervised_burned_map <- function(scoring_features, model, recipe,
                                        oof_summary = NULL,
                                        export_burned_like = TRUE,
                                        out_map_dir = NULL, config = NULL) {
  if (missing(scoring_features) || is.null(scoring_features)) {
    stop("'scoring_features' is required.", call. = FALSE)
  }
  if (missing(model) || is.null(model)) {
    stop("'model' is required (xgboost model object or RDS path).",
         call. = FALSE)
  }
  if (missing(recipe) || is.null(recipe)) {
    stop("'recipe' is required (named list or RDS path).", call. = FALSE)
  }
  if (!is.logical(export_burned_like) || length(export_burned_like) != 1L) {
    stop("'export_burned_like' must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(out_map_dir)) {
    if (is.null(config) ||
        !inherits(config, "otsufire_supervised_burned_config")) {
      stop("Supply either 'out_map_dir' or a valid 'config'.", call. = FALSE)
    }
    out_map_dir <- config$output_routes$final_map_dir
  }
  prefix <- if (!is.null(config))
    sprintf("%d_%s_patch_certified", config$target_year, config$scenario) else
    "scored"
  list(
    deterministic_scored      = NULL,
    final_map_full            = NULL,
    final_map                 = NULL,
    burned_like_scored        = NULL,
    final_map_gpkg            = file.path(out_map_dir,
                                   paste0(prefix, "_final_map.gpkg")),
    final_map_counts_csv      = file.path(out_map_dir,
                                   paste0(prefix, "_final_map_counts.csv")),
    burned_like_gpkg          = file.path(out_map_dir,
                                   paste0(prefix, "_burned_like_scored.gpkg")),
    burned_like_counts_csv    = file.path(out_map_dir,
                                   paste0(prefix, "_burned_like_scored_counts.csv")),
    note = paste0(
      "In 0.2.x scoring is performed by run_oneyear_supervised_pipeline(); ",
      "calling score_supervised_burned_map() standalone validates inputs."
    )
  )
}
