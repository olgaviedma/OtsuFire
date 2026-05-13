#' Train the final supervised burned-area model
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' Final-model training entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.3.x this is a
#' contract-shaped wrapper; real execution happens inside
#' [run_oneyear_supervised_pipeline()].
#'
#' The legacy `hard_negative_source` parameter is deprecated (Phase 2A
#' four-bucket taxonomy; HANDOFF Sections 25-30). It remains in the
#' signature for backward API compatibility and is ignored.
#'
#' @param train_features sf or GPKG path.
#' @param oof_agg data.frame or CSV path (optional).
#' @param hard_negative_source DEPRECATED. Ignored. Kept for signature
#'   compatibility.
#' @param random_negative_source Character vector.
#' @param out_dir Output folder; defaults to config route.
#' @param config Optional `otsufire_supervised_burned_config`.
#'
#' @keywords internal
#' @noRd
train_final_burned_model <- function(
    train_features, oof_agg = NULL,
    hard_negative_source = NULL,
    random_negative_source = c("random_burnable_background",
                                "legacy_review_residual",
                                "legacy_keep_residual"),
    out_dir = NULL, config = NULL) {
  if (missing(train_features) || is.null(train_features)) {
    stop("'train_features' is required.", call. = FALSE)
  }
  if (!is.null(hard_negative_source)) {
    warning(
      "'hard_negative_source' is deprecated and ignored. Phase 2A four-bucket ",
      "negative-pool taxonomy replaces it (see HANDOFF Sections 25-30).",
      call. = FALSE
    )
  }
  if (!is.character(random_negative_source) ||
      length(random_negative_source) < 1L) {
    stop("'random_negative_source' must be a non-empty character vector.",
         call. = FALSE)
  }
  if (is.null(out_dir)) {
    if (is.null(config) ||
        !inherits(config, "otsufire_supervised_burned_config")) {
      stop("Supply either 'out_dir' or a valid 'config'.", call. = FALSE)
    }
    out_dir <- config$output_routes$final_model_dir
  }
  prefix <- if (!is.null(config))
    sprintf("%d_%s_patch_certified", config$target_year, config$scenario) else
    "final_model"
  list(
    training_ok            = NULL,
    model                  = NULL,
    recipe                 = NULL,
    split_idx              = NULL,
    model_rds              = file.path(out_dir,
                                        paste0(prefix, "_final_model.rds")),
    recipe_rds             = file.path(out_dir,
                                        paste0(prefix, "_recipe.rds")),
    training_ok_gpkg       = file.path(out_dir,
                                        paste0(prefix,
                                               "_certified_training_ok.gpkg")),
    feature_importance_csv = file.path(out_dir,
                                        paste0(prefix,
                                               "_feature_importance.csv")),
    meta_txt               = file.path(out_dir, paste0(prefix, "_meta.txt")),
    model_summary_txt      = file.path(out_dir,
                                        paste0(prefix, "_model_summary.txt")),
    note = paste0(
      "In 0.2.x the final model is trained by run_oneyear_supervised_pipeline(); ",
      "calling train_final_burned_model() standalone validates inputs."
    )
  )
}
