#' Run the probabilistic refinement workflow for one year
#'
#' @description
#' Run the probabilistic refinement workflow for one target year using a
#' configuration created with [build_supervised_burned_config()].
#'
#' The function builds training pools from classified Otsu-guided patches and
#' background samples, assigns spatial folds, extracts features, generates
#' out-of-fold diagnostics, trains the final model, and scores candidate
#' polygons.
#'
#' Optional consistency checks compare the classifications of Otsu-guided
#' patches with the probabilistic refinement results.
#'
#' The function returns output paths, configuration information, and run
#' summaries. External map validation is performed separately with
#' [validate_fire_maps()] after selecting a threshold for `p_burned`.
#'
#' @param config An object of class `otsufire_supervised_burned_config`,
#'   created with [build_supervised_burned_config()]. Contains the inputs,
#'   sampling settings, feature controls, model parameters, and output
#'   locations.
#' @param run_consistency Logical scalar. Whether to compare the
#'   classifications of Otsu-guided patches with the probabilistic refinement results.
#'   Default: `TRUE`. These checks do not replace external map validation.
#' @param overwrite Logical scalar. Whether existing stage outputs may be
#'   regenerated and replaced. With `FALSE`, existing outputs are reused or
#'   skipped where supported. Default: `FALSE`.
#' @param reuse_upstream Logical scalar. Whether to skip pool generation,
#'   spatial-fold assignment, and feature extraction and use existing outputs
#'   from those stages. Default: `FALSE`. See \strong{Reusing previous
#'   outputs}.
#' @param ... Reserved for detecting unsupported or removed arguments.
#'   Additional arguments are rejected. The removed arguments
#'   `additional_drop_cols` and `extra_drop_cols` produce a migration message
#'   directing users to `feature_whitelist_override`.
#' @param random_to_burned_ratio Compatibility argument: cap on background
#'   negatives relative to the burned-label count. Preferred configuration
#'   setting: `negative_pool_params = list(caps = c(random = ..., otsu = ...))`.
#'   See \strong{Compatibility arguments}.
#' @param otsu_unburned_to_burned_ratio Compatibility argument: cap on
#'   Otsu-based negatives relative to the burned-label count. Preferred
#'   configuration setting:
#'   `negative_pool_params = list(caps = c(random = ..., otsu = ...))`.
#' @param feature_whitelist_override Compatibility argument: restrict the
#'   permitted predictor set. Preferred configuration setting:
#'   `feature_whitelist_override`.
#' @param feature_weights Compatibility argument: supply named per-feature
#'   weights. Preferred configuration setting: `feature_weights`.
#' @param oof_nrounds_max Compatibility argument: maximum boosting rounds.
#'   Preferred configuration setting: `nrounds_max`.
#' @param oof_early_stop Compatibility argument: early-stopping patience.
#'   Preferred configuration setting: `early_stop`.
#' @param oof_seed_base Compatibility argument: random seed for out-of-fold
#'   training. Preferred configuration setting: `oof_seed_base`.
#' @param final_sampling_seed Compatibility argument: random seed for final
#'   negative-pool sampling. Preferred configuration setting:
#'   `final_sampling_seed`.
#' @param final_seed Compatibility argument: random seed for the final split
#'   and model training. Preferred configuration setting: `final_seed`.
#' @param final_val_frac Compatibility argument: inner validation fraction.
#'   Preferred configuration setting: `val_frac`.
#' @param final_group_col Compatibility argument: grouping column for the
#'   inner train/validation split. Preferred configuration setting:
#'   `group_col`.
#' @param final_nrounds_max Compatibility argument: maximum boosting rounds.
#'   Preferred configuration setting: `nrounds_max`.
#' @param final_early_stopping_rounds Compatibility argument: early-stopping
#'   patience. Preferred configuration setting: `early_stop`.
#' @param final_impute_numeric Compatibility argument: numeric imputation
#'   rule, `"median"` or `"zero"`. Preferred configuration setting:
#'   `impute_numeric`.
#' @param final_impute_factor_missing Compatibility argument: category used
#'   for missing factor values. Preferred configuration setting:
#'   `impute_factor_missing`.
#'
#' @section Compatibility arguments:
#' The arguments `random_to_burned_ratio`, `otsu_unburned_to_burned_ratio`,
#' `feature_whitelist_override`, `feature_weights` and the `oof_*` /
#' `final_*` arguments are retained for compatibility. Their default value is
#' `NULL`, meaning that the resolved configuration value is used.
#'
#' For new code, supply the corresponding setting to
#' [build_supervised_burned_config()].
#'
#' The shared configuration controls keep out-of-fold and final training
#' settings aligned. Historical argument prefixes such as `oof_` and `final_`
#' should not be used to define separate training procedures.
#'
#' See \strong{Compatibility and parameter conflicts} for how non-`NULL`
#' compatibility arguments are resolved.
#'
#' @section Workflow stages:
#' The function runs the following stages:
#'
#' | Stage | Purpose |
#' |---|---|
#' | Training pools | Build burned and unburned training pools and the candidate pool to be scored. |
#' | Spatial folds | Assign spatial cross-validation folds. |
#' | Feature extraction | Compute predictor variables for training and scoring polygons. |
#' | Out-of-fold diagnostics | Evaluate predictions for spatially held-out labelled observations. |
#' | Final model training | Select the number of boosting rounds and fit the final model. |
#' | Scoring | Apply the final model to candidate polygons. |
#' | Consistency checks | Optionally compare probabilistic refinement results with the original classifications of Otsu-guided patches. |
#'
#' This function runs a single-year workflow. It does not perform
#' leave-one-year-out evaluation or coordinate multi-year training.
#'
#' @section Training eligibility and negative sampling:
#' Training uses explicit burned and unburned labels. Otsu-guided patches
#' classified as `"keep"` provide the initial burned pool after label
#' preparation.
#'
#' Unburned training examples must belong to a supported negative pool. The
#' standard pools are background negatives (`random`) and moderately
#' burned-like negatives (`otsu`). Additional strongly burned-like negatives
#' (`artifact_hard`) are included when enabled and eligible under the
#' configuration.
#'
#' Review-class candidates and rows with missing or unknown labels do not
#' automatically become negatives.
#'
#' The same negative-sampling policy is used for out-of-fold and final
#' training. Background and Otsu-based caps are applied within the relevant
#' training subset. Optional hard-negative weighting follows the settings in
#' `config$negative_pool_params`.
#'
#' @section Out-of-fold and final training:
#' Both stages use an inner validation split to select the number of boosting
#' rounds through early stopping. A fresh model is then fitted on the
#' available training observations using the selected number of rounds.
#'
#' During out-of-fold evaluation, this procedure is repeated within each
#' spatial training fold. The outer held-out fold is reserved for prediction
#' and evaluation.
#'
#' Preprocessing statistics and class-balance settings are derived from the
#' relevant training observations. Consequently, quantities such as
#' imputation medians and `scale_pos_weight` may differ between folds and the
#' final model.
#'
#' Out-of-fold diagnostics measure agreement with the internally derived
#' labels. They do not provide an independent assessment of map accuracy.
#'
#' @section Feature controls:
#' Configure the permitted feature set through `feature_whitelist_override`
#' in [build_supervised_burned_config()].
#'
#' The whitelist restricts the package-supported feature set; it does not
#' introduce arbitrary new predictors. Optional shape-feature availability is
#' controlled separately through `include_shape_features` in the
#' configuration.
#'
#' `feature_weights` supplies named, finite, non-negative weights for
#' design-matrix columns. Unspecified features receive a weight of 1. Names
#' outside the active feature space generate a warning and are excluded.
#'
#' These are predictor weights passed to XGBoost, not training-example
#' weights or direct multipliers of feature values.
#'
#' @section Scores and external validation:
#' The scored map includes `p_burned`. This is a model score and is not
#' necessarily a calibrated probability.
#'
#' To create a binary burned-area map, select a score threshold. The
#' thresholded output can then be evaluated against an external reference
#' with [validate_fire_maps()].
#'
#' Setting `run_consistency = TRUE` enables internal comparisons with
#' Otsu-guided patch classifications. It does not run external map
#' validation.
#'
#' @section Reusing previous outputs:
#' `overwrite` and `reuse_upstream` control different aspects of execution.
#'
#' | Setting | Behaviour |
#' |---|---|
#' | `overwrite = FALSE` | Reuse or skip existing stage outputs where supported. Useful for resuming an interrupted run with compatible settings. |
#' | `overwrite = TRUE` | Allow stage outputs to be regenerated and replaced. |
#' | `reuse_upstream = TRUE` | Explicitly skip pool generation, fold assignment, and feature extraction, and load their existing outputs. |
#'
#' Explicit upstream reuse requires the expected files under the configured
#' output routes, including `01_POOLS/<year>_<run_label>_pools.gpkg`,
#' `02_FOLDS/<year>_train_with_folds_<size>m.gpkg` and
#' `03_FEATURES/features_geometry.gpkg`.
#'
#' Reused products must be compatible with the current inputs and settings.
#' Changes to labels, sampling rules, spatial folds, or extracted feature
#' requirements may require rebuilding the affected stages.
#'
#' With `reuse_upstream = TRUE`, upstream generation remains skipped even
#' when later outputs are regenerated.
#'
#' @section Compatibility and parameter conflicts:
#' Methodological and training settings should be defined in
#' [build_supervised_burned_config()]. The corresponding function-level
#' arguments are deprecated compatibility options.
#'
#' A non-`NULL` compatibility value is handled as follows:
#'
#' | Situation | Behaviour |
#' |---|---|
#' | The configuration field uses a package default and the requested value differs | Apply the compatibility value and emit an `otsufire_deprecated_param` warning. |
#' | The configuration contains an explicit user value that differs from the requested value | Stop with an error because two explicit settings conflict. |
#' | The requested value equals the configuration value | Accept it without a deprecation warning. |
#'
#' Parameter resolution is recorded in `resolved_params_provenance`.
#' Provenance for individual XGBoost parameters is returned separately in
#' `model_params_provenance`.
#'
#' Migrate calls to the configuration builder to avoid relying on these
#' compatibility arguments.
#'
#' @section Memory settings:
#' During execution, the function temporarily sets the terra memory ceiling,
#' `memmax`, to 16 GB and restores the previous setting on exit.
#'
#' This setting controls terra processing; it is not a limit on the total
#' memory used by the R process. Available memory must also accommodate
#' feature tables, model training, and other workflow objects.
#'
#' @return A named list of class `otsufire_supervised_run`.
#'
#' | Field | Contents |
#' |---|---|
#' | `config` | Configuration used for the run, including resolved settings. |
#' | `result_dir` | Main result directory. |
#' | `pools_gpkg` | Path to the pool GeoPackage. |
#' | `train_with_folds_gpkg` | Path to the training layer with spatial-fold assignments. |
#' | `features_geometry_gpkg` | Path to the feature layer with geometries. |
#' | `oof_agg_csv` | Path to the aggregated out-of-fold output. |
#' | `final_model_rds` | Path to the saved final-model RDS file. |
#' | `final_map_gpkg` | Path to the scored candidate map, including `p_burned`. |
#' | `burned_like_gpkg` | Path to the burned-like output layer. |
#' | `timing_csv` | Path to the execution-time summary. |
#' | `consistency_*` | Fields associated with the optional consistency checks. |
#' | `resolved_params_provenance` | Record of configuration and compatibility-argument resolution. |
#' | `model_params_provenance` | Per-parameter provenance for the resolved XGBoost settings. |
#' | `legacy_run_summary` | Run summary retained for compatibility. |
#' | `legacy_consistency_summary` | Consistency summary retained for compatibility. |
#'
#' The `final_map_gpkg` output is a scored layer. A score threshold is needed
#' to derive a binary burned-area map.
#'
#' @seealso [build_supervised_burned_config()],
#'   [build_supervised_training_pools()], [make_spatial_folds()],
#'   [extract_supervised_features()], [run_oof_diagnostics()],
#'   [train_final_burned_model()], [score_supervised_burned_map()],
#'   [check_supervised_consistency()], [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' # Configure a probabilistic refinement run from classified Otsu-guided patches
#' config <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   reference_burned_map = "data/reference_burned_2022.gpkg",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022",
#'   negative_pool_params = list(
#'     caps = c(random = 1.0, otsu = 1.0)
#'   ),
#'   nrounds_max = 4000L,
#'   early_stop = 80L,
#'   oof_seed_base = 42L,
#'   final_sampling_seed = 42L,
#'   final_seed = 42L,
#'   random_seed = 42L
#' )
#'
#' # Run the complete single-year workflow
#' result <- run_oneyear_supervised_pipeline(
#'   config = config,
#'   run_consistency = TRUE,
#'   overwrite = FALSE
#' )
#'
#' # Inspect output locations
#' result$result_dir
#' result$final_model_rds
#' result$final_map_gpkg
#' result$timing_csv
#'
#' # Load the scored candidate map
#' scored_map <- sf::st_read(
#'   result$final_map_gpkg,
#'   quiet = TRUE
#' )
#' summary(scored_map$p_burned)
#'
#' # Inspect how parameters were resolved
#' result$resolved_params_provenance
#' result$model_params_provenance
#'
#' # Rerun downstream stages using compatible existing pools,
#' # spatial folds, and extracted features
#' result_reused <- run_oneyear_supervised_pipeline(
#'   config = config,
#'   reuse_upstream = TRUE,
#'   overwrite = TRUE,
#'   run_consistency = TRUE
#' )
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
      "Pass the subset of package features you want to keep ",
      "(the current whitelist minus the columns to drop). ",
      "See ?build_supervised_burned_config.",
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
