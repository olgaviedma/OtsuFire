#' Run the complete one-year supervised burned-area workflow
#'
#' @description
#' Runs the whole supervised workflow for one target year from a single
#' configuration object: training pools, spatial folds, feature extraction,
#' out-of-fold (OOF) diagnostics, final model training, scoring, and the
#' optional deterministic-vs-supervised consistency check. Each stage can
#' also be run on its own with the corresponding function
#' ([build_supervised_training_pools()], [make_spatial_folds()],
#' [extract_supervised_features()], [run_oof_diagnostics()],
#' [train_final_burned_model()], [score_supervised_burned_map()],
#' [check_supervised_consistency()]); the results are the same.
#'
#' The final map returned in `final_map_gpkg` carries the per-patch
#' `p_burned` score. Threshold it (for example with
#' [write_thresholded_burned()]) and pass it to [validate_fire_maps()] for
#' external validation.
#'
#' @param config `otsufire_supervised_burned_config` object from
#'   [build_supervised_burned_config()].
#' @param run_consistency Logical. Whether to run the
#'   deterministic-vs-supervised consistency check. Default `TRUE`.
#' @param overwrite Logical. Default `FALSE`: stage outputs that already
#'   exist on disk are reused, so an interrupted run can resume without
#'   recomputing completed stages. Set `TRUE` to regenerate every stage
#'   output.
#' @param random_to_burned_ratio,otsu_unburned_to_burned_ratio Numeric or
#'   `NULL`. Deprecated; set `negative_pool_params$caps` in
#'   [build_supervised_burned_config()] instead. Caps on the
#'   burnable-background and moderately burned-like negatives, as multiples
#'   of the number of burned labels. `NULL` (default) uses the configured
#'   caps (package default `1.0` each).
#' @param feature_whitelist_override Character vector or `NULL`. Deprecated;
#'   set it in [build_supervised_burned_config()] instead. Restricts the
#'   features the OOF and final models may use to a subset of the package's
#'   supervised feature list; names outside that list are an error, so the
#'   list can only be narrowed, never extended. `NULL` (default) uses the
#'   configured feature set.
#' @param feature_weights Named numeric vector or `NULL`. Deprecated; set it
#'   in [build_supervised_burned_config()] instead. Per-feature weights passed
#'   to XGBoost for the OOF and final models. Names must be model features (or
#'   their `_isNA` companions); unknown names are dropped with a warning and
#'   unnamed features keep weight 1. Values must be finite and `>= 0`: values
#'   in (0, 1) attenuate a feature, values above 1 boost it. `NULL` (default)
#'   uses the configured weights (uniform 1 by default).
#' @param reuse_upstream Logical. When `TRUE`, the pools, folds and features
#'   stages are skipped and their outputs are read from an existing run
#'   (`01_POOLS/<year>_<run_label>_pools.gpkg`,
#'   `02_FOLDS/<year>_train_with_folds_<size>m.gpkg` and
#'   `03_FEATURES/features_geometry.gpkg` must already exist). Default
#'   `FALSE`.
#' @param oof_nrounds_max,oof_early_stop,oof_seed_base Integer or `NULL`.
#'   Deprecated; set them in [build_supervised_burned_config()] instead.
#'   Maximum boosting rounds, early-stopping patience and base RNG seed of the
#'   OOF models. `NULL` (default) uses the configured values (package defaults
#'   `4000`, `80` and `42`).
#' @param final_sampling_seed,final_seed Integer or `NULL`. Deprecated; set
#'   them in [build_supervised_burned_config()] instead. RNG seeds of the
#'   final model's negative sampling and of its train/validation split and
#'   training. `NULL` (default) uses the configured values (package default
#'   `42`).
#' @param final_val_frac Numeric in (0, 1) or `NULL`. Deprecated. Validation
#'   fraction of the final model's train/validation split. `NULL` (default)
#'   uses the configured value (package default `0.15`).
#' @param final_group_col Character or `NULL`. Deprecated. Grouping column of
#'   the final model's grouped train/validation split. `NULL` (default) uses
#'   the configured value (package default `"block_id"`).
#' @param final_nrounds_max,final_early_stopping_rounds Integer or `NULL`.
#'   Deprecated. Maximum boosting rounds and early-stopping patience of the
#'   final model. `NULL` (default) uses the configured values (package
#'   defaults `4000` and `80`, the same as the OOF models).
#' @param final_impute_numeric Character or `NULL`. Deprecated.
#'   Numeric-imputation rule of the final model, `"median"` or `"zero"`.
#'   `NULL` (default) uses the configured value (package default
#'   `"median"`).
#' @param final_impute_factor_missing Character or `NULL`. Deprecated. Level
#'   used for missing factor/character values in the final model. `NULL`
#'   (default) uses the configured value (package default `"MISSING"`).
#' @param ... Not used. Any extra named argument is an error, so a misspelt
#'   or removed argument is never silently ignored. The removed arguments
#'   `additional_drop_cols` / `extra_drop_cols` give an error that points to
#'   `feature_whitelist_override`.
#'
#' @details
#' The OOF and final models share one training procedure: imputation
#' medians and `scale_pos_weight` are fit on the training rows only, an inner
#' validation split selects the number of boosting rounds by early stopping,
#' and a fresh model is refit on all training rows at that round count. OOF
#' applies the same capped negative sampling as the final model,
#' independently within each training fold.
#'
#' Training eligibility is defined by explicit class, never by negation: only
#' burned rows (positives) and unburned rows that belong to a valid negative
#' pool (`random`, `otsu`, and `artifact_hard` when enabled) enter training;
#' review, keep, `NA` or unknown rows never become negatives.
#'
#' @section Setting the methodological parameters:
#' Every methodological and training-control parameter (negative-pool caps,
#' feature set and weights, XGBoost settings, seeds, imputation) is set in
#' [build_supervised_burned_config()] and read from the configuration by every
#' stage, so the OOF and final models always use the same values.
#'
#' The corresponding arguments of this function are deprecated and will be
#' removed in a future version. They all default to `NULL` (use the
#' configured value). A non-`NULL` value:
#' \itemize{
#'   \item replaces a value the configuration left at its default, with a
#'     deprecation warning of class `"otsufire_deprecated_param"`;
#'   \item is an error if it contradicts a value set explicitly in the
#'     configuration;
#'   \item is applied silently if it equals the configured value.
#' }
#' Every resolution is recorded in the `resolved_params_provenance` field of
#' the returned object, and the per-field XGBoost settings in
#' `model_params_provenance`.
#'
#' @section Memory:
#' While it runs, the function raises the terra memory limit
#' (`terra::terraOptions(memmax)`) to 16 GB, so that large annual mosaics are
#' not processed feature by feature, and restores the previous value on exit.
#' The machine should have at least 16 GB of free RAM.
#'
#' @return A named list of class `otsufire_supervised_run` with fields:
#'   `config`, `result_dir`, `pools_gpkg`, `train_with_folds_gpkg`,
#'   `features_geometry_gpkg`, `oof_agg_csv`, `final_model_rds`,
#'   `final_map_gpkg`, `burned_like_gpkg`, `timing_csv`,
#'   `consistency_*`, `resolved_params_provenance`,
#'   `model_params_provenance` (per-field XGBoost settings and their origin),
#'   `legacy_run_summary`, `legacy_consistency_summary`.
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L,
#'   options = list(data_base = "D:/FIRE", composite_base = "D:/FIRE/Composites")
#' )
#' run <- run_oneyear_supervised_pipeline(cfg, run_consistency = TRUE)
#' run$final_map_gpkg
#' }
#'
#' @family workflow
#' @export
run_oneyear_supervised_pipeline <- function(config, run_consistency = TRUE,
                                            overwrite = FALSE,
                                            # Gate 1B (2026-06-07): all
                                            # methodological knobs DEFAULT TO
                                            # NULL = "use cfg$train_control /
                                            # cfg$model_params" (the single
                                            # source of truth). A non-NULL value
                                            # OVERRIDES the cfg value and is
                                            # written through to BOTH OOF and
                                            # FINAL identically. The canonical
                                            # numeric defaults live ONLY in
                                            # build_supervised_burned_config()
                                            # (.of_canonical_train_control /
                                            # .of_canonical_model_params); they
                                            # are no longer duplicated here.
                                            random_to_burned_ratio                 = NULL,
                                            otsu_unburned_to_burned_ratio          = NULL,
                                            feature_whitelist_override             = NULL,
                                            feature_weights                        = NULL,
                                            reuse_upstream                         = FALSE,
                                            oof_nrounds_max                        = NULL,
                                            oof_early_stop                         = NULL,
                                            oof_seed_base                          = NULL,
                                            final_sampling_seed                    = NULL,
                                            final_seed                             = NULL,
                                            final_val_frac                         = NULL,
                                            final_group_col                        = NULL,
                                            final_nrounds_max                      = NULL,
                                            final_early_stopping_rounds            = NULL,
                                            final_impute_numeric                   = NULL,
                                            final_impute_factor_missing            = NULL,
                                            ...) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }

  # 0.5.0 hard removal: `additional_drop_cols` was replaced by
  # `feature_whitelist_override`. Catch via `...` so old callers get a clear
  # migration message. Checked BEFORE the Gate 1B cfg resolution so the
  # migration message wins over the train_control validation for stale calls.
  .dots <- list(...)
  if ("training_protocol" %in% names(.dots)) {
    stop("training_protocol is no longer an argument; OtsuFire always uses ",
         "inner-early-stopping selection + full-data refit.", call. = FALSE)
  }
  if ("oof_sampling" %in% names(.dots)) {
    stop("oof_sampling is no longer an argument; OtsuFire OOF always uses the ",
         "same capped negative-sampling policy as the final model, applied ",
         "independently within each training fold.", call. = FALSE)
  }
  if ("additional_drop_cols" %in% names(.dots) ||
      "extra_drop_cols" %in% names(.dots)) {
    stop(
      "`additional_drop_cols` / `extra_drop_cols` were removed in ",
      "OtsuFire 0.5.0. Use `feature_whitelist_override` instead. ",
      "Pass a subset of `.supervised_feature_cols` (e.g., the ",
      "current whitelist minus the columns you want to drop). ",
      "See NEWS.md and the migration note in HANDOFF Section N+20.",
      call. = FALSE
    )
  }
  if ("cap_contextual" %in% names(.dots) ||
      "contextual_exclusion_to_burned_ratio" %in% names(.dots)) {
    stop(
      "`cap_contextual` / `contextual_exclusion_to_burned_ratio` were removed in ",
      "OtsuFire (GATE 6.5): the contextual (deterministic-drop) negative bucket ",
      "no longer exists: drop patches enter training only through the ",
      "moderately burned-like (`otsu`) and strongly burned-like ",
      "(`artifact_hard`) pools. Set the caps with ",
      "negative_pool_params$caps = c(random = , otsu = ).",
      call. = FALSE
    )
  }
  if (length(.dots) > 0L) {
    stop("Unused arguments passed to run_oneyear_supervised_pipeline(): ",
         paste(names(.dots), collapse = ", "), call. = FALSE)
  }

  # Gate 1B + Precision 1: resolve every methodological param from cfg at THIS
  # public boundary ONLY. The CANONICAL way to set these is
  # build_supervised_burned_config(); the function-level arguments are
  # DEPRECATED COMPATIBILITY SHIMS. .of_resolve_methodological_shim() folds each
  # override into the resolved value here and (a) errors on a conflict with an
  # explicit builder-user value, (b) warns (class "otsufire_deprecated_param")
  # when an override supersedes a canonical default, (c) records a provenance
  # row. Downstream (dispatcher / orchestrator / OOF / FINAL / engines) consume
  # ONLY the resolved scalars below, never the raw arg + cfg in parallel.
  .tc <- config$train_control
  if (is.null(.tc) || !is.list(.tc)) {
    stop("'config' has no resolved train_control; rebuild it with ",
         "build_supervised_burned_config().", call. = FALSE)
  }
  .prov  <- config$resolved_params_provenance$train_control %||% list()
  .rec   <- .of_new_shim_record()
  .shim  <- function(override, cfg_value, param, arg_name = param) {
    .of_resolve_methodological_shim(
      override = override, cfg_value = cfg_value, param = param,
      arg_name = arg_name, cfg_provenance = .prov[[param]] %||% "default",
      record = .rec)
  }
  random_to_burned_ratio                 <- .shim(random_to_burned_ratio,                 .tc$caps$random,     "cap_random",     "random_to_burned_ratio")
  otsu_unburned_to_burned_ratio          <- .shim(otsu_unburned_to_burned_ratio,          .tc$caps$otsu,       "cap_otsu",       "otsu_unburned_to_burned_ratio")
  feature_whitelist_override             <- .shim(feature_whitelist_override,             .tc$feature_whitelist_override, "feature_whitelist_override", "feature_whitelist_override")
  feature_weights                        <- .shim(feature_weights,                        .tc$feature_weights, "feature_weights", "feature_weights")
  oof_nrounds_max                        <- .shim(oof_nrounds_max,                        .tc$nrounds_max,     "nrounds_max",    "oof_nrounds_max")
  oof_early_stop                         <- .shim(oof_early_stop,                          .tc$early_stop,      "early_stop",     "oof_early_stop")
  oof_seed_base                          <- .shim(oof_seed_base,                          .tc$seeds$oof_seed_base,       "oof_seed_base",       "oof_seed_base")
  final_sampling_seed                    <- .shim(final_sampling_seed,                    .tc$seeds$final_sampling_seed, "final_sampling_seed", "final_sampling_seed")
  final_seed                             <- .shim(final_seed,                             .tc$seeds$final_seed,          "final_seed",          "final_seed")
  final_val_frac                         <- .shim(final_val_frac,                         .tc$val_frac,        "val_frac",       "final_val_frac")
  final_group_col                        <- .shim(final_group_col,                        .tc$group_col,       "group_col",      "final_group_col")
  final_nrounds_max                      <- .shim(final_nrounds_max,                      .tc$nrounds_max,     "nrounds_max",    "final_nrounds_max")
  final_early_stopping_rounds            <- .shim(final_early_stopping_rounds,            .tc$early_stop,      "early_stop",     "final_early_stopping_rounds")
  final_impute_numeric                   <- .shim(final_impute_numeric,                   .tc$impute_numeric,  "impute_numeric", "final_impute_numeric")
  final_impute_factor_missing            <- .shim(final_impute_factor_missing,            .tc$impute_factor_missing, "impute_factor_missing", "final_impute_factor_missing")
  # Precision 1 condition 6: the resolved-params provenance record (cfg value /
  # requested override / resolved value / provenance label per methodological
  # field). Attached to the run summary the pipeline emits below.
  .resolved_params_provenance <- .of_shim_record_to_df(.rec)
  # Gate 1B (2026-06-07): the PER-FIELD cfg$model_params provenance (canonical /
  # requested / resolved / provenance per xgb field). The builder folds a PARTIAL
  # model_params override onto the canonical block, so this is tracked per field
  # exactly like train_control. Surfaced on the run object so a downstream
  # manifest can render the model_params provenance table.
  .model_params_provenance <- .of_model_params_provenance_to_df(
    config$resolved_params_provenance$model_params)
  if (!is.logical(run_consistency) || length(run_consistency) != 1L) {
    stop("'run_consistency' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(reuse_upstream) || length(reuse_upstream) != 1L) {
    stop("'reuse_upstream' must be TRUE or FALSE.", call. = FALSE)
  }

  # terra memory ceiling: lift from default (~1 GB) to 16 GB for the
  # whole one-year supervised run, so large mosaics (e.g. ~16k polygons
  # in 2025) do not force terra into per-feature mode inside
  # extract_features() and downstream raster stages. Restored on exit
  # so the user's R session is not contaminated.
  .prev_memmax <- tryCatch(
    terra::terraOptions(print = FALSE)$memmax,
    error = function(e) NA_real_
  )
  terra::terraOptions(memmax = 16)
  on.exit({
    if (is.finite(.prev_memmax)) {
      terra::terraOptions(memmax = .prev_memmax)
    } else {
      terra::terraOptions(default = TRUE)
    }
  }, add = TRUE)

  .of_check_input_file(config$inputs$internal_decisions, "internal_decisions")
  .of_check_input_file(config$inputs$change_index,       "change_index")

  t0 <- Sys.time()
  res <- .of_run_supervised_oneyear(config,
                                     run_consistency = run_consistency,
                                     overwrite = overwrite,
                                     random_to_burned_ratio                 = random_to_burned_ratio,
                                     otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
                                     feature_whitelist_override             = feature_whitelist_override,
                                     feature_weights                        = feature_weights,
                                     reuse_upstream                         = reuse_upstream,
                                     oof_nrounds_max                        = oof_nrounds_max,
                                     oof_early_stop                         = oof_early_stop,
                                     oof_seed_base                          = oof_seed_base,
                                     final_sampling_seed                    = final_sampling_seed,
                                     final_seed                             = final_seed,
                                     final_val_frac                         = final_val_frac,
                                     final_group_col                        = final_group_col,
                                     final_nrounds_max                      = final_nrounds_max,
                                     final_early_stopping_rounds            = final_early_stopping_rounds,
                                     final_impute_numeric                   = final_impute_numeric,
                                     final_impute_factor_missing            = final_impute_factor_missing)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  legacy <- res$legacy_run_summary %||% list()
  result_dir <- legacy$result_dir %||%
    normalizePath(config$output_routes$base, winslash = "/", mustWork = FALSE)

  # B3 (2026-06-06): resolve prefixes from the SAME config knobs the
  # orchestrator and modular wrappers use (config$options$prefix_base /
  # prefix_oof_base, defaults "patch_certified" / "patch"), so the advertised
  # run-result paths match the files actually written when the option is
  # overridden. Defaults unchanged -> byte-identical for the default config.
  prefix_base     <- config$options$prefix_base %||% "patch_certified"
  prefix_oof_base <- config$options$prefix_oof_base %||% "patch"
  prefix     <- sprintf("%d_%s_%s",
                         config$target_year, config$scenario, prefix_base)
  prefix_oof <- sprintf("%d_%s_%s",
                         config$target_year, config$scenario, prefix_oof_base)

  timing_csv <- legacy$timing_csv %||% config$output_routes$timing_csv
  # §N+27 (2026-06-05): supervised burned-like registry removed (abandoned
  # research line). The run result no longer carries registry path fields.

  consistency_out <- res$consistency
  structure(
    list(
      config                  = config,
      result_dir              = result_dir,
      pools_gpkg              = file.path(result_dir, "01_POOLS",
                                 sprintf("%d_%s_pools.gpkg",
                                         config$target_year, config$scenario)),
      train_with_folds_gpkg   = file.path(result_dir, "02_FOLDS",
                                 sprintf("%d_train_with_folds_5000m.gpkg",
                                         config$target_year)),
      features_geometry_gpkg  = file.path(result_dir, "03_FEATURES",
                                           "features_geometry.gpkg"),
      oof_agg_csv             = file.path(result_dir, "05_OOF",
                                           paste0(prefix_oof, "_oof_agg.csv")),
      final_model_rds         = file.path(result_dir, "07_FINAL_MODEL_V2",
                                           paste0(prefix, "_final_model.rds")),
      final_map_gpkg          = file.path(result_dir, "09_FINAL_MAP",
                                           paste0(prefix, "_final_map.gpkg")),
      burned_like_gpkg        = file.path(result_dir, "09_FINAL_MAP",
                                           paste0(prefix, "_burned_like_scored.gpkg")),
      timing_csv              = timing_csv,
      consistency_issues_gpkg = consistency_out$issues_gpkg %||% NA_character_,
      consistency_summary_csv = consistency_out$summary_csv %||% NA_character_,
      consistency_summary_txt = consistency_out$summary_txt %||% NA_character_,
      elapsed_sec             = elapsed,
      # Precision 1 condition 6: resolved-params provenance written into the
      # run-level artifact the pipeline emits (cfg value / requested override /
      # resolved value / provenance label per methodological field). Also
      # carried on cfg$resolved_params_provenance (builder-time field flags).
      resolved_params_provenance = .resolved_params_provenance,
      # Gate 1B (2026-06-07): per-field cfg$model_params provenance table
      # (param / canonical / requested / resolved / provenance), e.g. a partial
      # build_supervised_burned_config(model_params = list(eta = 0.03)) yields
      # eta=user (requested 0.03) and every other xgb field default.
      model_params_provenance = .model_params_provenance,
      legacy_run_summary      = legacy,
      legacy_consistency_summary = consistency_out
    ),
    class = c("otsufire_supervised_run", "list")
  )
}

#' @export
print.otsufire_supervised_run <- function(x, ...) {
  cat("<otsufire_supervised_run>\n")
  cat("  scenario         :", x$config$scenario, "\n")
  cat("  target_year      :", x$config$target_year, "\n")
  cat("  result_dir       :", x$result_dir, "\n")
  cat("  final_map_gpkg   :", x$final_map_gpkg, "\n")
  cat("  burned_like_gpkg :", x$burned_like_gpkg, "\n")
  cat("  elapsed_sec      :", sprintf("%.1f", x$elapsed_sec), "\n")
  invisible(x)
}
