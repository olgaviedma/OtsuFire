#' Run out-of-fold diagnostics on supervised features
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' OOF-diagnostic entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.3.x this is a
#' contract-shaped wrapper; real execution happens inside
#' [run_oneyear_supervised_pipeline()].
#'
#' @param train_features data.frame / sf / RDS path.
#' @param scoring_features data.frame / sf / RDS path.
#' @param fold_cols Character vector of fold column names.
#' @param params Optional named list of xgboost params.
#' @param out_dir Output folder; defaults to config route.
#' @param config Optional `otsufire_supervised_burned_config`.
#'
#' @keywords internal
#' @noRd
run_oof_diagnostics <- function(train_features, scoring_features,
                                fold_cols = c("fold_rep1", "fold_rep2"),
                                params = NULL, out_dir = NULL,
                                config = NULL) {
  if (missing(train_features) || is.null(train_features)) {
    stop("'train_features' is required.", call. = FALSE)
  }
  if (missing(scoring_features) || is.null(scoring_features)) {
    stop("'scoring_features' is required.", call. = FALSE)
  }
  if (!is.character(fold_cols) || length(fold_cols) < 1L) {
    stop("'fold_cols' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.null(params) && !is.list(params)) {
    stop("'params' must be NULL or a named list.", call. = FALSE)
  }
  if (is.null(out_dir)) {
    if (is.null(config) ||
        !inherits(config, "otsufire_supervised_burned_config")) {
      stop("Supply either 'out_dir' or a valid 'config'.", call. = FALSE)
    }
    out_dir <- config$output_routes$oof_dir
  }
  prefix <- if (!is.null(config))
    sprintf("%d_%s_patch", config$target_year, config$scenario) else "oof"
  list(
    design_bundle            = NULL,
    oof_agg                  = NULL,
    oof_long                 = NULL,
    labeled_oof_summary      = NULL,
    oof_agg_csv              = file.path(out_dir, paste0(prefix, "_oof_agg.csv")),
    oof_long_csv             = file.path(out_dir, paste0(prefix, "_oof_long.csv")),
    labeled_oof_summary_gpkg = file.path(out_dir,
                                 paste0(prefix, "_labeled_oof_summary.gpkg")),
    oof_metrics_by_threshold = file.path(out_dir,
                                 paste0(prefix, "_oof_metrics_by_threshold.csv")),
    oof_metrics_summary      = file.path(out_dir,
                                 paste0(prefix, "_oof_metrics_summary.txt")),
    oof_best_thresholds      = file.path(out_dir,
                                 paste0(prefix, "_oof_best_thresholds.csv")),
    note = paste0(
      "In 0.2.x OOF is materialized by run_oneyear_supervised_pipeline(); ",
      "calling run_oof_diagnostics() standalone validates inputs."
    )
  )
}
