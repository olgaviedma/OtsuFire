#' Run the complete one-year supervised burned-area workflow
#'
#' @description
#' Top-level one-year supervised orchestrator defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. Chains pool
#' generation, spatial folds, feature extraction, OOF diagnostics, final
#' model training, scoring, and optional consistency checks for a single
#' target year and scenario.
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
#' @param overwrite Logical. Default `TRUE` regenerates/clobbers stage
#'   outputs (historical behaviour). When `FALSE`, the engine forwards
#'   the flag to the unburned builder and the OOF / final-model wrappers,
#'   which skip or reuse existing on-disk stage outputs instead of
#'   overwriting them. (Before 2026-06-05 this flag was cosmetic: it was
#'   validated but never threaded, and the engine always behaved as if
#'   `overwrite = TRUE`.)
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
#'   hotspot columns to attenuate them by 20x.
#'
#' @param reuse_upstream Logical. When `TRUE`, skip STEP A (pools),
#'   STEP B1-B2 (folds), and STEP B3 (features), and consume pre-existing
#'   upstream artefacts from disk instead. Requires an existing baseline
#'   run with `01_POOLS/<year>_<scenario>_pools.gpkg`,
#'   `02_FOLDS/<year>_train_with_folds_<size>m.gpkg`, and
#'   `03_FEATURES/features_geometry.gpkg` already present. Default
#'   `FALSE` reproduces the historical behaviour.
#'
#' @param oof_nrounds_max Integer. Maximum xgboost boosting rounds for the OOF
#'   per-fold models. Default `4000`. FRENTE 1 (2026-06-05): UNIFIED with the
#'   FINAL-model default (was `3000`). Threaded to [run_oof_diagnostics()].
#' @param oof_early_stop Integer. xgboost early-stopping patience for the OOF
#'   per-fold models. Default `80`. FRENTE 1: UNIFIED with the FINAL stage
#'   (was `75`).
#' @param oof_seed_base Integer. Base RNG seed for the repeated block-CV OOF
#'   runs. Default `42`. The canonical seed shared with the FINAL stage.
#' @param final_sampling_seed Integer. RNG seed for the FINAL-model
#'   negative-pool sampling. Default `42`. FRENTE 1: UNIFIED with the OOF
#'   canonical seed (was `999`). Threaded to [train_final_burned_model()].
#' @param final_seed Integer. RNG seed for the FINAL-model train/val split and
#'   xgboost training. Default `42`. FRENTE 1: UNIFIED with the OOF canonical
#'   seed (was `999`).
#' @param final_val_frac Numeric in (0, 1). Validation fraction for the
#'   FINAL-model train/val split. Default `0.15`.
#' @param final_group_col Character or `NULL`. Grouping column for the
#'   FINAL-model grouped train/val split. Default `"block_id"`.
#' @param final_nrounds_max Integer. Maximum xgboost boosting rounds for the
#'   FINAL model. Default `4000`. FRENTE 1 (2026-06-05): now UNIFIED with
#'   `oof_nrounds_max`.
#' @param final_early_stopping_rounds Integer. xgboost early-stopping patience
#'   for the FINAL model. Default `80`. FRENTE 1: now UNIFIED with
#'   `oof_early_stop`.
#' @param final_impute_numeric Character. Numeric-imputation rule for the
#'   FINAL model: `"median"` (default) or `"zero"`.
#' @param final_impute_factor_missing Character. Sentinel level for missing
#'   factor/character values in the FINAL model. Default `"MISSING"`.
#' @param training_protocol Character. Supervised training protocol, one of
#'   `"legacy"` (default) or `"nested_refit"`. `"legacy"` reproduces the
#'   historical OOF + FINAL behaviour byte-for-byte. `"nested_refit"` opts into
#'   the leakage-free protocol (B1, 2026-06-07) used IDENTICALLY by the OOF and
#'   FINAL stages via one shared internal core: medians and `scale_pos_weight`
#'   are fit on the training rows only, an inner validation split is the sole
#'   early-stopping set, and a fresh model is refit on all training rows at the
#'   selected `best_iteration` before deployment. Opt-in; the production default
#'   is unchanged.
#' @param oof_sampling Character. OOF per-fold negative sampling mode, one of
#'   `"capped"` (default) or `"full"`. Only active under
#'   `training_protocol = "nested_refit"`: `"capped"` applies the SAME four
#'   bucket caps the FINAL model uses to the outer-train negatives of each fold;
#'   `"full"` uses all outer-train rows. Ignored on the legacy path.
#'
#' @section D1 training knobs (2026-06-05) + FRENTE 1 unification:
#' The `oof_*` and `final_*` arguments expose the OOF and final-model training
#' hyperparameters that were previously hardcoded at the engine defaults. They
#' remain fully user-overridable. FRENTE 1 (2026-06-05) UNIFIED the OOF and
#' FINAL hyperparameters to a single CANONICAL set: the XGBoost params block
#' (objective/eval_metric/eta/max_depth/min_child_weight/subsample/
#' colsample_bytree/gamma/lambda/alpha) now comes from one source of truth
#' [.of_canonical_xgb_params()] used by BOTH stages, and the training controls
#' were aligned (`oof_nrounds_max` 3000 -> 4000, `oof_early_stop` 75 -> 80,
#' `final_sampling_seed`/`final_seed` 999 -> 42 to match `oof_seed_base`). This
#' is RESULT-AFFECTING (the FINAL model changes vs the historical baseline) and
#' was signed off by Natalia; it is first evaluated in the 2017 run. Only
#' `scale_pos_weight` differs between the two stages, because each computes it
#' from its own training labels.
#'
#' @param ... Migration trap for removed arguments. Reserved for
#'   detecting calls that still pass the OtsuFire 0.4.x deny-list
#'   arguments `additional_drop_cols` / `extra_drop_cols`, which were
#'   removed in 0.5.0; such calls error with a migration message
#'   pointing at `feature_whitelist_override`. Any other named argument
#'   passed here is also rejected with an "Unused arguments" error
#'   (never silently ignored).
#'
#' @section Phase B sampling caps and feature space:
#' The pass-through hooks above support the Phase B experiment matrix
#' (cap_contextual  in  {0.25, 1.0, Inf} x cap_spectral  in  {1.0, 2.0} x
#' hotspot-feature variants {All, L1, L2, None}) and the eventual mass
#' re-training (10 years x 3 scenarios). All defaults reproduce the
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
#' @section terra memory ceiling:
#' On entry the function lifts `terra::terraOptions(memmax)` to 16 GB
#' and restores the previous value on exit (`on.exit()`), so the
#' caller's R session is not contaminated. The default terra ceiling
#' (~1 GB) forces per-feature processing on large mosaics and was
#' exhausting RAM inside `extract_features()` and downstream raster
#' stages when the annual mosaic carried ~16k polygons (e.g. 2025).
#' The 16 GB ceiling assumes the host has at least that much free
#' RAM available; the package was developed on a 64 GB machine.
#'
#' @return A named list of class `otsufire_supervised_run` with fields:
#'   `config`, `result_dir`, `pools_gpkg`, `train_with_folds_gpkg`,
#'   `features_geometry_gpkg`, `oof_agg_csv`, `final_model_rds`,
#'   `final_map_gpkg`, `burned_like_gpkg`, `timing_csv`,
#'   `consistency_*`, `legacy_run_summary`,
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
                                            # FRENTE 1 (2026-06-05): OOF + FINAL
                                            # training knobs UNIFIED to the
                                            # canonical set. oof_nrounds_max
                                            # 3000 -> 4000, oof_early_stop 75 ->
                                            # 80, final_sampling_seed / final_seed
                                            # 999 -> 42 (= oof_seed_base). All
                                            # user-overridable.
                                            oof_nrounds_max                        = 4000,
                                            oof_early_stop                         = 80,
                                            oof_seed_base                          = 42,
                                            final_sampling_seed                    = 42,
                                            final_seed                             = 42,
                                            final_val_frac                         = 0.15,
                                            final_group_col                        = "block_id",
                                            final_nrounds_max                      = 4000,
                                            final_early_stopping_rounds            = 80,
                                            final_impute_numeric                   = "median",
                                            final_impute_factor_missing            = "MISSING",
                                            # B1 (2026-06-07): nested-refit
                                            # training protocol. Default "legacy"
                                            # reproduces the production model
                                            # byte-for-byte; "nested_refit" opts
                                            # into the leakage-free protocol
                                            # (B1a/B1b/B1c/B1d) used identically
                                            # by OOF and FINAL. `oof_sampling`
                                            # ("capped"/"full") controls the OOF
                                            # per-fold negative bucketing and is
                                            # only active under nested_refit.
                                            training_protocol                      = c("legacy", "nested_refit"),
                                            oof_sampling                           = c("capped", "full"),
                                            ...) {
  training_protocol <- match.arg(training_protocol)
  oof_sampling <- match.arg(oof_sampling)
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
      "See NEWS.md and the migration note in HANDOFF Section N+20.",
      call. = FALSE
    )
  }
  if (length(.dots) > 0L) {
    stop("Unused arguments passed to run_oneyear_supervised_pipeline(): ",
         paste(names(.dots), collapse = ", "), call. = FALSE)
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
                                     contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
                                     spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
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
                                     final_impute_factor_missing            = final_impute_factor_missing,
                                     # B1 (2026-06-07): forward the protocol
                                     # toggle + OOF sampling mode.
                                     training_protocol                      = training_protocol,
                                     oof_sampling                           = oof_sampling)
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
