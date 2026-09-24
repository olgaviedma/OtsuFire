#' Run spatial out-of-fold diagnostics for the supervised burned-area model
#'
#' @description
#' Evaluate the supervised burned-area model using repeated spatial
#' cross-validation and summarise its performance across score thresholds.
#'
#' Each held-out fold is scored by a model fitted without that fold. The
#' function returns out-of-fold predictions, diagnostic metrics, threshold
#' summaries, and the design-matrix bundle used by the supervised workflow.
#'
#' Run this stage after [extract_supervised_features()] and before
#' [train_final_burned_model()].
#'
#' The diagnostics measure agreement with the internal burned and unburned
#' labels, including labels derived from Otsu-guided patches. External map
#' accuracy must be assessed separately with suitable reference data.
#'
#' @param train_features `sf` object, `data.frame`, or GeoPackage path
#'   containing labelled training features. Must contain `class` and the
#'   columns listed in `fold_cols`. When a path is supplied, reads the
#'   `train_features` layer.
#' @param scoring_features `sf` object, `data.frame`, or GeoPackage path
#'   containing features for the candidate scoring pool. When a path is
#'   supplied, reads the `scoring_features` layer.
#' @param config Required object of class `otsufire_supervised_burned_config`,
#'   created with [build_supervised_burned_config()]. Supplies resolved
#'   methodological settings, run identifiers, and default output locations.
#' @param fold_cols Character vector naming the fold-assignment columns to
#'   evaluate. Default: `c("fold_rep1", "fold_rep2")`. These columns must
#'   exist in `train_features`.
#' @param params Optional named list of XGBoost parameters. Leave as `NULL` to
#'   use the shared configuration-based parameter resolution. For consistent
#'   OOF and final training, set model parameters through `model_params` in
#'   [build_supervised_burned_config()].
#' @param feature_whitelist_override Character vector restricting the
#'   permitted feature set, or `NULL` to use the configured setting.
#'   Deprecated as a function-level control; configure it in the builder.
#' @param feature_weights Named numeric vector of per-feature weights, or
#'   `NULL` to use the configured setting. Deprecated as a function-level
#'   control; configure it in the builder.
#' @param include_shape_features Logical scalar or `NULL`. Whether shape and
#'   size features are eligible for modelling. When `NULL`, uses
#'   `config$train_control$include_shape_features`. Required features must
#'   have been prepared during extraction.
#' @param nrounds_max Training control: maximum number of boosting rounds per
#'   fold. Builder argument: `nrounds_max`. See \strong{Training controls}.
#' @param early_stop Training control: early-stopping patience. Builder
#'   argument: `early_stop`.
#' @param seed_base Training control: base random seed for repeated OOF
#'   training. Builder argument: `oof_seed_base`.
#' @param random_to_burned_ratio Training control: background-negative cap
#'   relative to the burned-label count. Builder setting:
#'   `negative_pool_params$caps["random"]`.
#' @param otsu_unburned_to_burned_ratio Training control: Otsu-negative cap
#'   relative to the burned-label count. Builder setting:
#'   `negative_pool_params$caps["otsu"]`.
#' @param val_frac Training control: inner validation fraction used for early
#'   stopping. Builder argument: `val_frac`.
#' @param impute_numeric Training control: numeric imputation rule,
#'   `"median"` or `"zero"`. Builder argument: `impute_numeric`.
#' @param impute_factor_missing Training control: category used for missing
#'   factor values. Builder argument: `impute_factor_missing`.
#' @param group_col Training control: grouping column used for the inner
#'   train/validation split. Builder argument: `group_col`.
#' @param source_weights Optional named numeric vector assigning a row weight
#'   to each specified value of the `source` column. When `NULL`, uses
#'   `config$negative_pool_params$source_weights` if available. These are
#'   training-example weights, distinct from `feature_weights`.
#' @param artifact_hard_source Character vector identifying `source` values
#'   belonging to promoted hard negatives, or `NULL`. When `NULL`, resolves to
#'   `"artifact_hard"` if the mechanism is enabled, otherwise `character(0)`.
#' @param total_weight_ratio Numeric scalar or `NULL`. Target ratio of total
#'   hard-negative weight to total burned-pool weight. When `NULL`, uses
#'   `config$negative_pool_params$artifact_hard$total_weight_ratio`. With no
#'   source-weight override and no hard-negative rows, no explicit
#'   sample-weight vector is constructed.
#' @param out_dir Character scalar or `NULL`. Directory for OOF results.
#'   Defaults to `config$output_routes$oof_dir`, normally `05_OOF`.
#' @param matrix_dir Character scalar or `NULL`. Directory for the
#'   design-matrix bundle. Defaults to `config$output_routes$matrix_dir`,
#'   normally `04_MATRIX`.
#' @param labelled_gpkg Character scalar or `NULL`. Path to the training
#'   GeoPackage with fold assignments, used to attach geometry to aggregated
#'   OOF results. When `NULL`, the spatial OOF summary is not written.
#' @param labelled_layer Character scalar. Layer to read from
#'   `labelled_gpkg`. Default: `"train_with_folds"`.
#' @param overwrite Logical scalar. Passed to the underlying OOF workflow to
#'   control replacement of the design-matrix bundle. Default: `TRUE`.
#' @param .internal_resolved Internal argument used by the full pipeline to
#'   indicate that compatibility arguments have already been resolved. Leave
#'   as `FALSE` for direct calls.
#'
#' @section Training controls:
#' The arguments `nrounds_max`, `early_stop`, `seed_base`,
#' `random_to_burned_ratio`, `otsu_unburned_to_burned_ratio`, `val_frac`,
#' `impute_numeric`, `impute_factor_missing` and `group_col` default to
#' `NULL`, which uses the configured value. They are retained as
#' compatibility controls; set them through [build_supervised_burned_config()]
#' for new workflows.
#'
#' @section Workflow:
#' The function:
#' 1. Prepares the design-matrix bundle from the extracted features.
#' 2. Runs spatially blocked OOF training for each supplied repetition.
#' 3. Predicts the held-out observations.
#' 4. Aggregates OOF predictions and computes metrics across the evaluated
#'    thresholds.
#' 5. Reports selected thresholds under the implemented diagnostic criteria.
#' 6. Writes the OOF results and design-matrix bundle.
#' 7. Optionally joins aggregated predictions to polygon geometries.
#'
#' The supplied fold columns determine the outer evaluation partitions. This
#' function uses the assignments prepared by [make_spatial_folds()].
#'
#' @section Training within each fold:
#' For each outer fold, the remaining observations form the available
#' training partition.
#'
#' Within that partition, the workflow applies training eligibility and
#' negative-sampling rules, uses an inner validation split to select the
#' number of boosting rounds, and refits the preprocessing recipe and model on
#' the selected outer-training observations.
#'
#' The fitted model then predicts the outer held-out fold.
#'
#' The outer test observations must remain excluded from preprocessing
#' estimation, feature selection, sample-weight resolution, early stopping,
#' and model fitting. Preparing a shared matrix bundle does not make the
#' held-out observations eligible for training.
#'
#' @section Relationship to final-model training:
#' OOF and final training use a shared training procedure and configuration.
#' This keeps the feature policy, sampling rules, and training controls
#' aligned.
#'
#' The fitted models still differ because they use different training
#' subsets. Imputation statistics, class-balance weights, selected boosting
#' rounds, and predictions can therefore differ across folds and from the
#' final model.
#'
#' Configure shared settings through [build_supervised_burned_config()] to
#' keep the diagnostic and final training procedures comparable.
#'
#' @section Training labels and negative pools:
#' Training eligibility requires explicit burned or unburned labels. Missing,
#' unknown, or review labels are not automatically interpreted as unburned.
#'
#' Eligible negatives belong to supported pools:
#' * background negatives, identified as `random`;
#' * moderately burned-like negatives, identified as `otsu`;
#' * promoted strongly burned-like negatives, identified through
#'   `artifact_hard_source`, when enabled.
#'
#' Background and Otsu caps are applied within each training fold. Promoted
#' hard negatives use their configured eligibility and weighting rules.
#'
#' The labels derived from Otsu-guided patch classifications may contain
#' errors. OOF evaluation measures how well the model predicts those labels
#' under the selected spatial partition.
#'
#' @section Features and weights:
#' The active feature whitelist and per-feature weights follow the resolved
#' settings.
#'
#' When shape features are enabled, the permitted feature set includes
#' `area_ha`, `n_pix`, `log_area`, `perim_m`, `compactness` and
#' `elongation`. Associated missingness indicators may also be included in
#' the design matrix.
#'
#' Per-feature weights and training-example weights have different purposes:
#'
#' | Weight control | Applies to |
#' |---|---|
#' | `feature_weights` | Predictor columns. |
#' | `source_weights` | Training rows grouped by source. |
#' | `total_weight_ratio` | Total contribution of the hard-negative pool relative to the burned pool. |
#'
#' Use consistent weighting settings for OOF and final training.
#'
#' @section Threshold diagnostics:
#' The function evaluates OOF scores across a set of thresholds and writes
#' metric summaries and selected thresholds.
#'
#' These thresholds are derived from the internal labelled dataset. They
#' provide diagnostic choices for converting scores into classes, rather than
#' a universally optimal threshold for burned-area mapping.
#'
#' `p_burned` is a model score and is not necessarily a calibrated
#' probability. Threshold performance may differ when applied to the full
#' candidate population or evaluated against an independent reference.
#'
#' Metrics used to select a threshold should not also be presented as an
#' independent evaluation of that selection.
#'
#' @section Spatial OOF summaries:
#' Supply `labelled_gpkg` to join aggregated OOF results back to the training
#' geometries.
#'
#' Use the actual training file returned by [make_spatial_folds()], since the
#' selected block size can vary between runs. The full pipeline supplies this
#' file automatically.
#'
#' Without `labelled_gpkg`, tabular OOF outputs remain available, but
#' `labeled_oof_summary` is `NULL`.
#'
#' @section Compatibility arguments:
#' Function-level training controls and the feature-whitelist and
#' feature-weight arguments are deprecated compatibility options.
#'
#' For these arguments:
#' * `NULL` uses the configuration value;
#' * an override of a package-default setting emits an
#'   `otsufire_deprecated_param` warning;
#' * an override that conflicts with an explicit builder setting raises an
#'   error.
#'
#' Use [build_supervised_burned_config()] to define the methodological
#' settings for new workflows.
#'
#' @section Output files:
#' OOF filenames use the prefix `<year>_<run_label>_patch`.
#'
#' | File | Contents |
#' |---|---|
#' | `<prefix>_oof_agg.csv` | Aggregated OOF predictions by unit. |
#' | `<prefix>_oof_long.csv` | Detailed OOF predictions by fold and repetition. |
#' | `<prefix>_oof_metrics_by_threshold.csv` | Metrics across evaluated thresholds. |
#' | `<prefix>_oof_metrics_summary.txt` | Diagnostic summary. |
#' | `<prefix>_oof_best_thresholds.csv` | Thresholds selected by the diagnostic criteria. |
#' | `<prefix>_labeled_oof_summary.gpkg` | Aggregated OOF results with geometries, when requested. |
#'
#' The design-matrix bundle is written separately under `matrix_dir`. Use the
#' returned paths to locate the outputs.
#'
#' @return A named list containing diagnostic objects and output paths.
#'
#' | Field | Contents |
#' |---|---|
#' | `design_bundle` | Design-matrix bundle retained for downstream final-model training. |
#' | `oof_agg` | `data.frame` of aggregated OOF predictions by unit. |
#' | `oof_long` | `data.frame` of detailed OOF predictions by fold and repetition. |
#' | `labeled_oof_summary` | Spatial OOF summary as an `sf` object, or `NULL` when the labelled geometry source is unavailable. |
#' | `oof_agg_csv` | Path to aggregated OOF predictions. |
#' | `oof_long_csv` | Path to detailed OOF predictions. |
#' | `oof_metrics_by_threshold_csv` | Path to threshold-dependent metrics. |
#' | `oof_metrics_summary_txt` | Path to the diagnostic summary. |
#' | `oof_best_thresholds_csv` | Path to selected thresholds. |
#' | `labeled_oof_summary_gpkg` | Path to the spatial OOF summary when written. |
#' | `design_bundle_rds` | Path to the saved design-matrix bundle. |
#'
#' @seealso [build_supervised_burned_config()], [make_spatial_folds()],
#'   [extract_supervised_features()], [train_final_burned_model()],
#'   [score_supervised_burned_map()], [validate_supervised_execution()],
#'   [run_oneyear_supervised_pipeline()], [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' # Configure the supervised workflow
#' config <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022",
#'   nrounds_max = 4000L,
#'   early_stop = 80L,
#'   oof_seed_base = 42L
#' )
#'
#' # Build the pools and assign spatial folds
#' pools <- build_supervised_training_pools(config)
#'
#' folds <- make_spatial_folds(
#'   train_labelled = pools$pools_gpkg,
#'   config = config
#' )
#'
#' if (!isTRUE(folds$selected$ok)) {
#'   stop("Review the fallback partition before continuing.")
#' }
#'
#' # Extract features without hotspot predictors
#' features <- extract_supervised_features(
#'   train_with_folds = folds$train_with_folds_gpkg,
#'   scoring_pool = pools$pools_gpkg,
#'   config = config,
#'   use_hotspots = FALSE
#' )
#'
#' # Use the fold-assignment columns present in the feature table
#' fold_columns <- grep(
#'   "^fold_rep",
#'   names(features$train_features),
#'   value = TRUE
#' )
#'
#' # Run OOF diagnostics and attach the training geometries
#' oof <- run_oof_diagnostics(
#'   train_features = features$train_features,
#'   scoring_features = features$scoring_features,
#'   config = config,
#'   fold_cols = fold_columns,
#'   labelled_gpkg = folds$train_with_folds_gpkg
#' )
#'
#' # Inspect aggregated predictions
#' head(oof$oof_agg)
#'
#' # Inspect thresholds selected by the diagnostic criteria
#' thresholds <- read.csv(oof$oof_best_thresholds_csv)
#' thresholds
#'
#' # Locate the diagnostic and model-input outputs
#' oof$oof_metrics_by_threshold_csv
#' oof$labeled_oof_summary_gpkg
#' oof$design_bundle_rds
#' }
#'
#' @family workflow
#' @export
run_oof_diagnostics <- function(train_features, scoring_features,
                                config,
                                fold_cols = c("fold_rep1", "fold_rep2"),
                                params = NULL,
                                # Gate 1B (2026-06-07): all methodological knobs
                                # DEFAULT TO NULL = "read from cfg$train_control /
                                # cfg$model_params" (single source of truth). A
                                # non-NULL value overrides cfg for standalone use.
                                # No literal methodological numbers live here.
                                feature_whitelist_override = NULL,
                                feature_weights = NULL,
                                # OPTIONAL shape/size block (OtsuFire 0.12.0).
                                # NULL = read cfg$train_control$include_shape_features
                                # (single source of truth shared with FINAL). A
                                # non-NULL logical overrides cfg for standalone use.
                                # OFF -> byte-identical OOF.
                                include_shape_features = NULL,
                                nrounds_max = NULL,
                                early_stop = NULL,
                                seed_base = NULL,
                                random_to_burned_ratio                 = NULL,
                                otsu_unburned_to_burned_ratio          = NULL,
                                val_frac = NULL,
                                impute_numeric = NULL,
                                impute_factor_missing = NULL,
                                group_col = NULL,
                                # PHASE 2 (artifact_hard): optional per-row source
                                # weighting forwarded to the OOF engine.
                                # `source_weights` is a NAMED numeric
                                # (source -> weight); `artifact_hard_source` marks
                                # promoted artifact_hard rows. Both DEFAULT TO NULL
                                # = "read from cfg$negative_pool_params"
                                # (source_weights) / use the canonical
                                # "artifact_hard" tag when the feature is enabled.
                                # OFF by default -> byte-identical to today.
                                source_weights = NULL,
                                artifact_hard_source = NULL,
                                # PHASE 2 (artifact_hard): pool-level weight
                                # balance. NULL = read from
                                # cfg$negative_pool_params$artifact_hard$total_weight_ratio.
                                total_weight_ratio = NULL,
                                out_dir = NULL,
                                matrix_dir = NULL,
                                labelled_gpkg = NULL,
                                labelled_layer = "train_with_folds",
                                overwrite = TRUE,
                                # Precision 1 (2026-06-07): internal sentinel set
                                # TRUE by the orchestrator, which already resolved
                                # the methodological shims at the
                                # run_oneyear_supervised_pipeline() boundary
                                # (warning + conflict-error + provenance happen
                                # ONCE there). When TRUE this standalone boundary
                                # skips re-warning so the orchestrated path does
                                # not double-warn. Direct callers leave it FALSE
                                # and get the full deprecated-shim treatment.
                                .internal_resolved = FALSE) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(train_features) || is.null(train_features)) {
    stop("'train_features' is required.", call. = FALSE)
  }
  if (missing(scoring_features) || is.null(scoring_features)) {
    stop("'scoring_features' is required.", call. = FALSE)
  }
  if (missing(config) || is.null(config) ||
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  # Gate 1B + Precision 1: resolve methodological params from cfg at this public
  # boundary. The CANONICAL path is build_supervised_burned_config(); the
  # function-level args here are DEPRECATED COMPATIBILITY SHIMS. Each override is
  # folded via .of_resolve_methodological_shim() (conflict-error vs an explicit
  # builder-user value; deprecation warning of class "otsufire_deprecated_param"
  # over a canonical default; provenance recording). When invoked by the
  # orchestrator (.internal_resolved=TRUE) the shims were ALREADY resolved at the
  # run_oneyear_supervised_pipeline() boundary, so we use the plain cfg-precedence
  # %||% to avoid a second warning for the same override.
  .tc <- config$train_control
  if (is.null(.tc) || !is.list(.tc)) {
    stop("'config' has no resolved train_control; rebuild it with ",
         "build_supervised_burned_config().", call. = FALSE)
  }
  if (isTRUE(.internal_resolved)) {
    .shim <- function(override, cfg_value, param, arg_name = param) override %||% cfg_value
  } else {
    .prov <- config$resolved_params_provenance$train_control %||% list()
    .shim <- function(override, cfg_value, param, arg_name = param) {
      .of_resolve_methodological_shim(
        override = override, cfg_value = cfg_value, param = param,
        arg_name = arg_name, cfg_provenance = .prov[[param]] %||% "default",
        record = NULL)
    }
  }
  nrounds_max <- .shim(nrounds_max, .tc$nrounds_max, "nrounds_max")
  early_stop  <- .shim(early_stop,  .tc$early_stop,  "early_stop")
  seed_base   <- .shim(seed_base,   .tc$seeds$oof_seed_base, "oof_seed_base", "seed_base")
  feature_whitelist_override <- .shim(feature_whitelist_override, .tc$feature_whitelist_override, "feature_whitelist_override")
  feature_weights            <- .shim(feature_weights,            .tc$feature_weights, "feature_weights")
  random_to_burned_ratio                 <- .shim(random_to_burned_ratio,                 .tc$caps$random,     "cap_random",     "random_to_burned_ratio")
  otsu_unburned_to_burned_ratio          <- .shim(otsu_unburned_to_burned_ratio,          .tc$caps$otsu,       "cap_otsu",       "otsu_unburned_to_burned_ratio")
  val_frac              <- .shim(val_frac,              .tc$val_frac,        "val_frac")
  impute_numeric        <- .shim(impute_numeric,        .tc$impute_numeric,  "impute_numeric")
  impute_factor_missing <- .shim(impute_factor_missing, .tc$impute_factor_missing, "impute_factor_missing")
  # Gate 1B (2026-06-07): group_col is sourced from cfg$train_control the SAME
  # way as val_frac / impute_* etc. so OOF's grouped block-CV column matches
  # the FINAL stage's group_col (single source of truth). Previously OOF used a
  # hardcoded "block_id" literal in run_dm_oof_pipeline(); that literal is
  # removed and group_col is now threaded from here.
  group_col             <- .shim(group_col,            .tc$group_col,       "group_col")
  # OPTIONAL shape/size block (OtsuFire 0.12.0): resolve from cfg when the arg
  # is left NULL (single source of truth shared with FINAL via the SAME cfg
  # field, so OOF and FINAL cannot diverge). Plain cfg-precedence %||% (no
  # deprecation warning: the FALSE default is the canonical baseline).
  include_shape_features <- include_shape_features %||% .tc$include_shape_features
  include_shape_features <- isTRUE(include_shape_features)
  if (!is.character(fold_cols) || length(fold_cols) < 1L) {
    stop("'fold_cols' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.null(params) && !is.list(params)) {
    stop("'params' must be NULL or a named list.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve out_dir / matrix_dir / result_dir / prefix / year-scenario from
  #    config, mirroring the orchestrator (dirs$05_OOF == oof_dir,
  #    dirs$04_MATRIX == matrix_dir, result_dir == output_routes$base,
  #    prefix_oof == "<year>_<run_label>_patch").
  # ---------------------------------------------------------------------------
  if (is.null(out_dir)) {
    out_dir <- config$output_routes$oof_dir
  }
  if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L ||
      !nzchar(out_dir)) {
    stop("Could not resolve 'out_dir' (05_OOF) from config.", call. = FALSE)
  }
  if (is.null(matrix_dir)) {
    matrix_dir <- config$output_routes$matrix_dir
  }
  if (is.null(matrix_dir) || !is.character(matrix_dir) ||
      length(matrix_dir) != 1L || !nzchar(matrix_dir)) {
    stop("Could not resolve 'matrix_dir' (04_MATRIX) from config.",
         call. = FALSE)
  }
  result_dir  <- config$output_routes$base
  target_year <- config$target_year
  scenario    <- config$scenario
  # The orchestrator builds prefix_oof as
  # sprintf("%d_%s_%s", year, scenario, prefix_oof_base) where prefix_oof_base
  # comes from config$options$prefix_oof_base (default "patch"). B3
  # (2026-06-06): read the SAME knob here so overriding the option keeps the
  # advertised paths consistent across orchestrator / run-result / wrapper.
  # Default unchanged -> byte-identical for the default config.
  prefix_oof_base <- config$options$prefix_oof_base %||% "patch"
  prefix_oof  <- sprintf("%d_%s_%s", target_year, scenario, prefix_oof_base)
  save_prefix <- paste0(target_year, "_", scenario, "_", prefix_oof_base)

  # ---------------------------------------------------------------------------
  # 2) Resolve inputs to data.frames. Accept sf / data.frame / GPKG path. When
  #    a path is supplied, read the canonical feature layers (the same layers
  #    the orchestrator reads from features_geometry.gpkg). The OOF stage body
  #    historically starts at sf::st_drop_geometry(train_features) /
  #    (scoring_features), so we reproduce that exactly.
  # ---------------------------------------------------------------------------
  read_feature_layer <- function(x, layer) {
    if (inherits(x, "SpatVector")) {
      return(sf::st_as_sf(x))
    }
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      return(sf::read_sf(x, layer = layer))
    }
    x
  }
  # PHASE 2 (artifact_hard): resolve the per-row weighting inputs. source_weights
  # defaults to the cfg value; artifact_hard_source defaults to the canonical tag
  # when the feature is enabled, else character(0) (OFF -> NO weight vector ->
  # byte-identical to today).
  if (is.null(source_weights)) {
    source_weights <- config$negative_pool_params$source_weights
  }
  if (is.null(artifact_hard_source)) {
    .ah_enabled <- isTRUE(config$negative_pool_params$artifact_hard$enabled)
    artifact_hard_source <- if (.ah_enabled) "artifact_hard" else character(0)
  }
  if (is.null(total_weight_ratio)) {
    total_weight_ratio <-
      config$negative_pool_params$artifact_hard$total_weight_ratio %||% 0.10
  }

  train_features_sf   <- read_feature_layer(train_features,   "train_features")
  scoring_features_sf <- read_feature_layer(scoring_features, "scoring_features")

  labelled    <- if (inherits(train_features_sf, "sf")) {
    sf::st_drop_geometry(train_features_sf)
  } else {
    as.data.frame(train_features_sf)
  }
  burned_like <- if (inherits(scoring_features_sf, "sf")) {
    sf::st_drop_geometry(scoring_features_sf)
  } else {
    as.data.frame(scoring_features_sf)
  }
  labelled_df <- labelled

  # ---------------------------------------------------------------------------
  # 3) Build XGB params (faithful MOVE of orchestrator STEP C1). When params is
  #    supplied use it verbatim; otherwise rebuild inline EXACTLY as today.
  # ---------------------------------------------------------------------------
  if (is.null(params)) {
    n_pos <- sum(labelled$class == "burned",   na.rm = TRUE)
    n_neg <- sum(labelled$class == "unburned", na.rm = TRUE)
    spw   <- if (n_pos > 0) n_neg / n_pos else 1

    # Gate 1B (2026-06-07): build the xgb params from cfg$model_params (the
    # SINGLE SOURCE OF TRUTH), then merge the site-specific scale_pos_weight
    # (computed here from the OOF labels). cfg$model_params is itself sourced
    # from .of_canonical_model_params(), so OOF and FINAL can never diverge.
    params <- config$model_params
    params[["scale_pos_weight"]] <- spw
  }

  # ---------------------------------------------------------------------------
  # 4) DM -> OOF (faithful MOVE of orchestrator STEP C2). Every argument and
  #    value matches the historical inline run_dm_oof_pipeline(...) call.
  # ---------------------------------------------------------------------------
  pipe1 <- run_dm_oof_pipeline(
    labelled    = labelled,
    burned_like = burned_like,
    labelled_df = labelled_df,
    params      = params,

    result_dir  = result_dir,
    target_year = target_year,
    fold_cols   = fold_cols,
    # 2026-06-05 (D1 expose): forward the OOF training knobs. Defaults equal
    # the engine's historical hardcoded values -> byte-identical OOF outputs.
    nrounds_max = nrounds_max,
    early_stop  = early_stop,
    seed_base   = seed_base,
    prefix      = prefix_oof,

    id_cols     = c("fire_uid", "class", "source", "poly_id", "block_id",
                    fold_cols),   # fold-rep columns dynamic (was hard-coded fold_rep1/2)
    # 0.4.0 (Agent H): the orchestrator does not pass drop_regex. The OOF
    # wrapper applies the canonical whitelist filter directly; drop_regex is a
    # no-op since 0.4.0. 2026-06-05: ecoregions removed; no categorical
    # predictors remain.
    cat_cols    = character(0),
    hs_n_col    = "hs_used_n",
    hs_conf_col = "hs_conf_mean",
    hs_frp_col  = "hs_frp_max",
    median_from = "labelled",
    save_dir_dm = matrix_dir,
    save_prefix = save_prefix,
    overwrite   = overwrite,
    # 0.5.0: pass the active whitelist override and per-feature weights down to
    # the OOF stage so it trains under the same conditions as the final model
    # (KB1/KB2 symmetry).
    feature_whitelist_override = feature_whitelist_override,
    feature_weights            = feature_weights,
    # OPTIONAL shape/size block (OtsuFire 0.12.0): forward the SAME flag the
    # FINAL stage uses so the OOF whitelist universe matches FINAL exactly.
    include_shape_features     = include_shape_features,

    # The cap ratios (the SAME variables forwarded to FINAL) so the OOF chain
    # carries identical values. OOF always uses the capped negative-sampling
    # policy, applied independently within each training fold (no toggle). GATE
    # 6.5: contextual cap removed (random + otsu only).
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    val_frac          = val_frac,
    impute_numeric    = impute_numeric,
    impute_factor_missing = impute_factor_missing,
    # Gate 1B (2026-06-07): thread the cfg-sourced grouped-CV column down to the
    # OOF engine (was a hardcoded "block_id" literal in the wrapper).
    group_col         = group_col,
    # PHASE 2 (artifact_hard): forward the per-row weighting inputs + the
    # artifact_hard source tag (run_dm_oof_pipeline threads them to run_oof_xgb).
    # Default-off -> byte-identical.
    source_weights       = source_weights,
    artifact_hard_source = artifact_hard_source,
    total_weight_ratio   = total_weight_ratio,

    labelled_gpkg  = labelled_gpkg,
    labelled_layer = labelled_layer
  )

  # ---------------------------------------------------------------------------
  # 5) Assemble the documented return (objects + paths). The path fields mirror
  #    what the inner wrapper wrote; the kept objects (design_bundle / oof_agg /
  #    oof_summary) are exactly what the downstream final-model stage consumes.
  # ---------------------------------------------------------------------------
  files <- pipe1$files
  labeled_oof_summary <- NULL
  if (!is.null(labelled_gpkg) && is.character(labelled_gpkg) &&
      length(labelled_gpkg) == 1L &&
      !is.null(files$labeled_oof_summary_gpkg) &&
      file.exists(files$labeled_oof_summary_gpkg)) {
    labeled_oof_summary <- tryCatch(
      sf::read_sf(files$labeled_oof_summary_gpkg,
                  layer = "labeled_oof_summary"),
      error = function(e) NULL
    )
  }

  list(
    design_bundle                 = pipe1$dm,
    oof_agg                       = pipe1$oof$oof_agg,
    oof_long                      = pipe1$oof$oof_long,
    labeled_oof_summary           = labeled_oof_summary,
    oof_agg_csv                   = files$oof_agg_path,
    oof_long_csv                  = files$oof_long_path,
    oof_metrics_by_threshold_csv  = files$oof_metrics_path,
    oof_metrics_summary_txt       = files$oof_metrics_summary_path,
    oof_best_thresholds_csv       = files$oof_best_thresholds_path,
    labeled_oof_summary_gpkg      = files$labeled_oof_summary_gpkg,
    design_bundle_rds             = file.path(
      matrix_dir, paste0(save_prefix, "_design_bundle.rds")
    ),
    pipe = pipe1
  )
}
