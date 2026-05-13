#' Run the complete one-year supervised burned-area workflow
#'
#' @description
#' Top-level one-year supervised orchestrator defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. Chains pool
#' generation, spatial folds, feature extraction, OOF diagnostics, final
#' model training, scoring, scenario-specific registry update, and
#' optional consistency checks for a single target year and scenario.
#'
#' Delegates the heavy numerical chain to the validated supervised engine
#' (`run_supervised_pipeline()` inside `03_PATCH_LEVEL_PIPELINE_2.R`) via
#' the internal dispatcher, with configuration injected. Mirrors the
#' Block 2b deterministic dispatcher pattern.
#'
#' The generic multi-year `run_supervised_pipeline()` is intentionally
#' NOT exported in 0.2.x and is deferred until LOYO is rebuilt.
#'
#' Thresholded `p_burned` outputs from the returned `final_map_gpkg` may
#' be passed to [validate_fire_maps()] for external validation.
#'
#' @param config `otsufire_supervised_burned_config` object from
#'   [build_supervised_burned_config()].
#' @param run_consistency Logical. Whether to execute the
#'   deterministic-vs-supervised consistency block.
#' @param overwrite Logical. Forwarded to the engine.
#' @param contextual_exclusion_to_burned_ratio Numeric. Cap on the
#'   contextual-exclusion negative pool (deterministic-drop polygons whose
#'   `neg_type` is NOT `spectral_reject_medium`), expressed as a multiple
#'   of `n_burned`, applied when fitting the FINAL supervised model.
#'   Default `0.25` reproduces historical behaviour. Use `Inf` to disable
#'   the cap. See "Phase B sampling caps" below.
#' @param spectral_hard_negative_to_burned_ratio Numeric. Cap on the
#'   spectral-hard-negative pool (deterministic-drop polygons with
#'   `neg_type = "spectral_reject_medium"`) for the FINAL model, as a
#'   multiple of `n_burned`. Default `1.0` (historical behaviour).
#' @param random_to_burned_ratio Numeric. Cap on
#'   random-burnable-background negatives for the FINAL model, as a
#'   multiple of `n_burned`. Default `1.0` (historical behaviour).
#' @param otsu_unburned_to_burned_ratio Numeric. Cap on Otsu current-year
#'   unburned-patch negatives for the FINAL model, as a multiple of
#'   `n_burned`. Default `1.0` (historical behaviour).
#' @param feature_whitelist_override Character vector. Optional subset
#'   of the canonical `.supervised_feature_cols` whitelist that the
#'   supervised stages (OOF + final model) are allowed to see. Names
#'   not in `.supervised_feature_cols` are rejected with an error: the
#'   canonical list is fixed and this argument can only RESTRICT it,
#'   never extend it. Default `NULL` reproduces the canonical 51-feature
#'   behaviour. The OOF wrapper and the final model both consume the
#'   same active whitelist, so OOF metrics and final-model behaviour
#'   stay symmetric. Phase B Exp4b example: pass
#'   `setdiff(OtsuFire:::.supervised_feature_cols,
#'           c("hs_in_poly","hs_in_buffer","hs_min_dist_m","hs_support_present"))`
#'   to drop the four hotspot-presence features.
#' @param feature_weights Named numeric vector forwarded to xgboost via
#'   the `xgb.DMatrix` `feature_weights` info on every per-fold OOF
#'   dtrain/dtest and on the final-model dtrain/dval. Names must
#'   correspond to columns of the supervised design matrix
#'   (i.e. members of the active whitelist or their `_isNA`
#'   companions). Names that do not appear in the active feature space
#'   trigger a warning and are silently dropped; features the caller
#'   does not name implicitly receive weight 1.0. Values must be
#'   finite and >= 0; values in (0, 1) attenuate, > 1 boost. Default
#'   `NULL` reproduces the canonical uniform-1.0 behaviour. Phase B
#'   Exp5a example: pass
#'   `setNames(rep(0.05, length(hs_cols)), hs_cols)` for the 13
#'   hotspot columns to attenuate them by 20×.
#'
#' @param reuse_upstream Logical. When `TRUE`, skip STEP A (pools),
#'   STEP B1-B2 (folds), and STEP B3 (features), and consume pre-existing
#'   upstream artefacts from disk instead. Requires an existing baseline
#'   run with `01_POOLS/<year>_<scenario>_pools.gpkg`,
#'   `02_FOLDS/<year>_train_with_folds_<size>m.gpkg`, and
#'   `03_FEATURES/features_geometry.gpkg` already present. Default
#'   `FALSE` reproduces the historical behaviour.
#'
#' @section Phase B sampling caps and feature space:
#' The pass-through hooks above support the Phase B experiment matrix
#' (cap_contextual ∈ {0.25, 1.0, Inf} × cap_spectral ∈ {1.0, 2.0} ×
#' hotspot-feature variants {All, L1, L2, None}) and the eventual mass
#' re-training (10 years × 3 scenarios). All defaults reproduce the
#' historical behaviour byte-for-byte; a caller that does not pass
#' them obtains exactly the same final model as before.
#'
#' OtsuFire 0.5.0 (2026-05-09): the deprecated `additional_drop_cols`
#' deny-list was removed. Callers that pass it now get a hard error
#' with a migration message pointing at `feature_whitelist_override`.
#' Phase B Exp4b is reproduced by passing
#' `feature_whitelist_override = setdiff(.supervised_feature_cols,
#' c("hs_in_poly","hs_in_buffer","hs_min_dist_m","hs_support_present"))`;
#' Phase B Exp5a is reproduced by passing a `feature_weights` vector
#' that sets the 13 hotspot columns to 0.05.
#'
#' @return A named list of class `otsufire_supervised_run` with fields:
#'   `config`, `result_dir`, `pools_gpkg`, `train_with_folds_gpkg`,
#'   `features_geometry_gpkg`, `oof_agg_csv`, `final_model_rds`,
#'   `final_map_gpkg`, `burned_like_gpkg`, `timing_csv`,
#'   `burned_like_registry_path`, `consistency_*`, `legacy_run_summary`,
#'   `legacy_consistency_summary`.
#'
#' @family workflow
#' @export
run_oneyear_supervised_pipeline <- function(config, run_consistency = TRUE,
                                            overwrite = FALSE,
                                            contextual_exclusion_to_burned_ratio   = 0.25,
                                            spectral_hard_negative_to_burned_ratio = 1.0,
                                            random_to_burned_ratio                 = 1.0,
                                            otsu_unburned_to_burned_ratio          = 1.0,
                                            feature_whitelist_override             = NULL,
                                            feature_weights                        = NULL,
                                            reuse_upstream                         = FALSE,
                                            ...) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (!is.logical(run_consistency) || length(run_consistency) != 1L) {
    stop("'run_consistency' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(reuse_upstream) || length(reuse_upstream) != 1L) {
    stop("'reuse_upstream' must be TRUE or FALSE.", call. = FALSE)
  }

  # 0.5.0 hard removal: `additional_drop_cols` was replaced by
  # `feature_whitelist_override`. Catch via `...` so old callers get
  # a clear migration message.
  .dots <- list(...)
  if ("additional_drop_cols" %in% names(.dots) ||
      "extra_drop_cols" %in% names(.dots)) {
    stop(
      "`additional_drop_cols` / `extra_drop_cols` were removed in ",
      "OtsuFire 0.5.0. Use `feature_whitelist_override` instead. ",
      "Pass a subset of `.supervised_feature_cols` (e.g., the ",
      "current whitelist minus the columns you want to drop). ",
      "See NEWS.md and the migration note in HANDOFF §N+20.",
      call. = FALSE
    )
  }
  if (length(.dots) > 0L) {
    stop("Unused arguments passed to run_oneyear_supervised_pipeline(): ",
         paste(names(.dots), collapse = ", "), call. = FALSE)
  }

  .of_check_input_file(config$inputs$internal_decisions, "internal_decisions")
  .of_check_input_file(config$inputs$change_index,       "change_index")

  t0 <- Sys.time()
  res <- .of_run_supervised_oneyear(config,
                                     run_consistency = run_consistency,
                                     overwrite = overwrite,
                                     contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
                                     spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
                                     random_to_burned_ratio                 = random_to_burned_ratio,
                                     otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
                                     feature_whitelist_override             = feature_whitelist_override,
                                     feature_weights                        = feature_weights,
                                     reuse_upstream                         = reuse_upstream)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  legacy <- res$legacy_run_summary %||% list()
  result_dir <- legacy$result_dir %||%
    normalizePath(config$output_routes$base, winslash = "/", mustWork = FALSE)

  prefix     <- sprintf("%d_%s_patch_certified",
                         config$target_year, config$scenario)
  prefix_oof <- sprintf("%d_%s_patch",
                         config$target_year, config$scenario)

  timing_csv <- legacy$timing_csv %||% config$output_routes$timing_csv
  registry_path <- legacy$registry_path %||% config$burned_like_registry_path
  registry_ranked_csv <- if (is.character(registry_path) &&
                              nzchar(registry_path))
    file.path(dirname(registry_path), "burned_high_conf_registry_ranked.csv") else NA_character_
  registry_summary_txt <- if (is.character(registry_path) &&
                               nzchar(registry_path))
    file.path(dirname(registry_path), "burned_high_conf_registry_summary.txt") else NA_character_

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
      burned_like_registry_path = registry_path,
      burned_like_registry_ranked_csv = registry_ranked_csv,
      burned_like_registry_summary_txt = registry_summary_txt,
      consistency_issues_gpkg = consistency_out$issues_gpkg %||% NA_character_,
      consistency_summary_csv = consistency_out$summary_csv %||% NA_character_,
      consistency_summary_txt = consistency_out$summary_txt %||% NA_character_,
      elapsed_sec             = elapsed,
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
  cat("  registry_path    :", x$burned_like_registry_path %||% "<none>", "\n")
  cat("  elapsed_sec      :", sprintf("%.1f", x$elapsed_sec), "\n")
  invisible(x)
}
