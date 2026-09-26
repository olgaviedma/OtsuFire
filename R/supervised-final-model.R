#' Train the final supervised burned-area model
#'
#' @description
#' Train the final XGBoost burned-area model using the eligible labelled
#' examples selected under the configured sampling and weighting rules.
#'
#' The function selects the number of boosting rounds using an inner
#' validation split, then refits the model on the complete selected training
#' set. It saves the fitted model, preprocessing recipe, training records, and
#' diagnostic summaries.
#'
#' Run this stage after feature extraction and out-of-fold diagnostics. The
#' returned `model` and `recipe` are used by [score_supervised_burned_map()]
#' to score candidate polygons.
#'
#' @param train_features `sf` object, `data.frame`, or GeoPackage path
#'   containing labelled training features. Must contain `class`, `source`,
#'   `neg_type`, and the required feature columns. When a path is supplied,
#'   the layer specified by `labelled_layer` is read.
#' @param config Required object of class `otsufire_supervised_burned_config`,
#'   created with [build_supervised_burned_config()]. Supplies the resolved
#'   model parameters, sampling settings, feature controls, random seeds, and
#'   output locations.
#' @param oof_agg Character scalar or `NULL`. Path to the aggregated OOF CSV
#'   returned by [run_oof_diagnostics()]. Used to enrich the model summary
#'   with OOF threshold information. Does not affect model fitting. Default:
#'   `NULL`.
#' @param labelled_layer Character scalar. Layer read when `train_features` is
#'   a GeoPackage path. Default: `"train_features"`.
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
#' @param random_to_burned_ratio Training control: cap on background
#'   negatives as a multiple of the burned-label count. Builder setting:
#'   `negative_pool_params$caps["random"]`. See \strong{Training and sampling
#'   controls}.
#' @param otsu_unburned_to_burned_ratio Training control: cap on Otsu-based
#'   negatives as a multiple of the burned-label count. Builder setting:
#'   `negative_pool_params$caps["otsu"]`.
#' @param sampling_seed Training control: random seed for negative-pool
#'   sampling. Builder argument: `final_sampling_seed`.
#' @param seed Training control: random seed for the inner split and XGBoost
#'   training. Builder argument: `final_seed`.
#' @param val_frac Training control: inner validation fraction, between `0`
#'   and `1`, excluding the endpoints. Builder argument: `val_frac`.
#' @param group_col Training control: column identifying observations kept
#'   together in the inner train/validation split. Builder argument:
#'   `group_col`.
#' @param nrounds_max Training control: maximum number of boosting rounds.
#'   Builder argument: `nrounds_max`.
#' @param early_stopping_rounds Training control: early-stopping patience.
#'   Builder argument: `early_stop`.
#' @param impute_numeric Training control: numeric imputation rule,
#'   `"median"` or `"zero"`. Builder argument: `impute_numeric`.
#' @param impute_factor_missing Training control: category used for missing
#'   factor values. Builder argument: `impute_factor_missing`.
#' @param source_weights Optional named numeric vector assigning row weights
#'   by values of the `source` column, for example `c(artifact_hard = 0.3)`.
#'   When `NULL`, reads `config$negative_pool_params$source_weights` if
#'   available.
#' @param artifact_hard_source Character vector identifying source values for
#'   promoted hard negatives, or `NULL`. When `NULL`, resolves to
#'   `"artifact_hard"` if the mechanism is enabled, otherwise `character(0)`.
#' @param total_weight_ratio Numeric scalar or `NULL`. Target ratio of total
#'   hard-negative weight to total burned-pool weight. When `NULL`, uses
#'   `config$negative_pool_params$artifact_hard$total_weight_ratio`. Ignored
#'   when an explicit hard-negative source weight is supplied through
#'   `source_weights`. With no source-weight override and no hard-negative
#'   rows, no explicit sample-weight vector is constructed.
#' @param out_dir Character scalar or `NULL`. Output directory. Defaults to
#'   `config$output_routes$final_model_dir`, normally the run's
#'   `07_FINAL_MODEL_V2` folder.
#' @param canonical_oof_fingerprint Internal argument used by the full
#'   pipeline to compare the final feature schema with the OOF schema. A
#'   mismatch raises an error. Leave as `NULL` for standalone use, which
#'   computes and saves its own fingerprint.
#' @param overwrite Logical scalar. Controls replacement of the
#'   training-selection GeoPackage. Default: `TRUE`. This argument is not
#'   documented as a general reuse switch for all model files.
#' @param verbose Logical scalar. Whether to display progress messages.
#'   Default: `TRUE`.
#' @param .internal_resolved Internal argument indicating that compatibility
#'   settings have already been resolved by the pipeline. Leave as `FALSE` for
#'   direct calls.
#'
#' @section Training and sampling controls:
#' The arguments `random_to_burned_ratio`, `otsu_unburned_to_burned_ratio`,
#' `sampling_seed`, `seed`, `val_frac`, `group_col`, `nrounds_max`,
#' `early_stopping_rounds`, `impute_numeric` and `impute_factor_missing`
#' default to `NULL`, meaning that the resolved configuration value is used.
#' They are retained for compatibility; use the corresponding builder
#' settings for new workflows.
#'
#' @section Model-fitting procedure:
#' The function uses the following procedure:
#' 1. Identify eligible burned and unburned training examples.
#' 2. Apply the configured negative-pool caps and sampling rules.
#' 3. Resolve the active features and training weights.
#' 4. Create an inner train/validation split, respecting the configured
#'    grouping column.
#' 5. Use early stopping to select the number of boosting rounds.
#' 6. Refit the preprocessing recipe and a fresh model on all selected
#'    eligible training examples.
#' 7. Save the refitted model, recipe, training records, and summaries.
#'
#' The saved model is the full selected-set refit, not the intermediate model
#' used for early stopping.
#'
#' "All selected training examples" refers to the observations remaining
#' after eligibility checks and sampling. It does not necessarily include
#' every row supplied in `train_features`.
#'
#' @section Training labels and negative pools:
#' Training requires explicit burned or unburned labels.
#'
#' Positive labels originate from the labelled burned pool, including
#' Otsu-guided patches classified as `"keep"` during label preparation.
#' Review candidates and rows with missing or unknown labels do not
#' automatically become negatives.
#'
#' Eligible negative sources include:
#' * background negatives, identified as `random`;
#' * moderately burned-like negatives, identified as `otsu`;
#' * promoted strongly burned-like negatives, identified through
#'   `artifact_hard_source`, when enabled.
#'
#' Background and Otsu negatives are subject to their configured caps.
#' Hard-negative contribution is controlled through the applicable weighting
#' rules.
#'
#' @section Weighting:
#' Predictor weights and training-example weights serve different purposes.
#'
#' | Setting | Applies to |
#' |---|---|
#' | `feature_weights` | Predictor columns. |
#' | `source_weights` | Training rows belonging to named sources. |
#' | `total_weight_ratio` | Total hard-negative contribution relative to the burned pool. |
#'
#' Under pool-ratio weighting:
#' \preformatted{
#' total hard-negative weight =
#'   total_weight_ratio x total burned-pool weight
#' }
#'
#' An explicit hard-negative weight in `source_weights` takes precedence over
#' this ratio. For example, a source weight of `0.3` is a per-row setting and
#' does not imply that the hard-negative pool has 30% of the burned pool's
#' total weight.
#'
#' Use consistent sampling and weighting settings for OOF and final training.
#'
#' @section Relationship to OOF diagnostics:
#' The final model and OOF models use the same training procedure and shared
#' configuration.
#'
#' They are fitted on different training subsets, so their preprocessing
#' statistics, class-balance weights, selected boosting rounds, and
#' predictions may differ.
#'
#' The optional `oof_agg` file adds OOF threshold information to the model
#' summary. It does not change the selected training rows, fit a model,
#' calibrate scores, or alter the final fitted model.
#'
#' OOF diagnostics assess agreement with internal labels. Independent map
#' validation is still needed to assess burned-area mapping accuracy.
#'
#' @section Feature recipe:
#' The returned recipe records the transformations needed for prediction,
#' including imputation information, feature names and order, and applied
#' feature controls.
#'
#' Keep the model and its recipe together. Use [score_supervised_burned_map()]
#' to apply the saved transformations consistently to candidate features.
#'
#' When shape features are enabled, the permitted feature set includes
#' `area_ha`, `n_pix`, `log_area`, `perim_m`, `compactness` and
#' `elongation`. Associated missingness indicators may also be included.
#'
#' @section Feature-schema checks:
#' When the full pipeline supplies `canonical_oof_fingerprint`, the function
#' checks that the final structural feature schema matches the OOF schema
#' before training and saving.
#'
#' Standalone calls with `canonical_oof_fingerprint = NULL` create their own
#' fingerprint. They do not, through that argument alone, verify agreement
#' with a previous OOF run.
#'
#' @section Inner validation split:
#' `split_idx` records the inner train/validation partition used to select
#' the boosting-round count.
#'
#' The final refit subsequently uses both parts of the selected dataset. The
#' inner validation observations are therefore not an independent test set
#' for the saved final model.
#'
#' @section Compatibility arguments:
#' Methodological settings should be supplied through
#' [build_supervised_burned_config()].
#'
#' The function-level caps, feature whitelist and weights, seeds, validation
#' fraction, grouping, boosting controls, and imputation arguments are
#' deprecated compatibility options.
#'
#' For these arguments:
#' * `NULL` uses the resolved configuration value;
#' * an override of a package-default setting emits an
#'   `otsufire_deprecated_param` warning;
#' * an override that conflicts with an explicit builder setting raises an
#'   error.
#'
#' @return A named list containing fitted objects, training records, and
#'   output paths.
#'
#' \strong{Returned objects}
#'
#' | Field | Contents |
#' |---|---|
#' | `model` | Fitted XGBoost model from the final refit. |
#' | `recipe` | Preprocessing and feature specification associated with the fitted model. |
#' | `training_ok` | Selected eligible training polygons as an `sf` object. |
#' | `split_idx` | Indices defining the inner train/validation split used for round selection. |
#'
#' \strong{Output paths}
#'
#' | Field | Contents |
#' |---|---|
#' | `final_model_rds` | Path to the saved final model. |
#' | `recipe_rds` | Path to the saved recipe. |
#' | `feature_importance_csv` | Path to the feature-importance table. |
#' | `model_summary_txt` | Path to the model summary, including OOF information when supplied. |
#' | `meta_txt` | Path to the model metadata file. |
#' | `split_idx_rds` | Path to the saved split indices. |
#' | `training_ok_csv` | Path to the selected training table. |
#' | `training_ok_gpkg` | Path to the selected training polygons. |
#'
#' Outputs are written under `out_dir`, normally `07_FINAL_MODEL_V2`. Use the
#' returned paths to locate the files.
#'
#' @seealso [build_supervised_burned_config()],
#'   [extract_supervised_features()], [run_oof_diagnostics()],
#'   [score_supervised_burned_map()], [validate_supervised_execution()],
#'   [run_oneyear_supervised_pipeline()].
#'
#' @examples
#' \dontrun{
#' # Configure final training with the same settings used for OOF
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
#'   negative_pool_params = list(
#'     caps = c(random = 1.0, otsu = 1.0)
#'   ),
#'   nrounds_max = 4000L,
#'   early_stop = 80L,
#'   final_sampling_seed = 42L,
#'   final_seed = 42L,
#'   val_frac = 0.15,
#'   group_col = "block_id",
#'   impute_numeric = "median",
#'   impute_factor_missing = "MISSING"
#' )
#'
#' # Train from existing feature and OOF outputs
#' # Replace these paths with the outputs from the preceding stages.
#' final <- train_final_burned_model(
#'   train_features = "path/to/features_geometry.gpkg",
#'   config = config,
#'   oof_agg = "path/to/2022_balanced_patch_oof_agg.csv"
#' )
#'
#' # Inspect the fitted model and preprocessing recipe
#' final$model
#' final$recipe
#'
#' # Inspect the examples selected for training
#' table(final$training_ok$class, useNA = "ifany")
#' table(final$training_ok$source, useNA = "ifany")
#'
#' # Locate the saved model and recipe
#' final$final_model_rds
#' final$recipe_rds
#'
#' # Inspect feature importance
#' importance <- read.csv(final$feature_importance_csv)
#' head(importance)
#'
#' # The model and recipe are now available for use with
#' # score_supervised_burned_map().
#' }
#'
#' @family workflow
#' @export
train_final_burned_model <- function(
    train_features,
    config,
    oof_agg = NULL,
    # Gate 1B (2026-06-07): all methodological knobs DEFAULT TO NULL = "read
    # from cfg$train_control / cfg$model_params" (single source of truth). A
    # non-NULL value overrides cfg for standalone use. No literal methodological
    # numbers live here.
    random_to_burned_ratio                 = NULL,
    otsu_unburned_to_burned_ratio          = NULL,
    feature_whitelist_override = NULL,
    feature_weights = NULL,
    # OPTIONAL shape/size block (OtsuFire 0.12.0). NULL = read
    # cfg$train_control$include_shape_features (the SAME cfg field the OOF stage
    # reads, so OOF and FINAL cannot diverge). A non-NULL logical overrides cfg
    # for standalone use. OFF -> byte-identical FINAL model.
    include_shape_features = NULL,
    sampling_seed = NULL,
    seed = NULL,
    val_frac = NULL,
    group_col = NULL,
    nrounds_max = NULL,
    early_stopping_rounds = NULL,
    impute_numeric = NULL,
    impute_factor_missing = NULL,
    # PHASE 2 (artifact_hard): optional per-row source weighting. `source_weights`
    # is a NAMED numeric (source -> weight); `artifact_hard_source` is the source
    # value(s) marking promoted artifact_hard rows. Both DEFAULT TO NULL = "read
    # from cfg$negative_pool_params" (source_weights) / use the canonical
    # "artifact_hard" tag. With no source_weights override AND no artifact_hard
    # row present this is byte-identical to today (NULL weight vector).
    source_weights = NULL,
    artifact_hard_source = NULL,
    # PHASE 2 (artifact_hard): pool-level weight balance (artifact_hard total
    # weight / burned total weight). NULL = read from
    # cfg$negative_pool_params$artifact_hard$total_weight_ratio. Ignored when an
    # explicit artifact_hard pin is supplied via source_weights.
    total_weight_ratio = NULL,
    labelled_layer = "train_features",
    out_dir = NULL,
    # Gate 1E (2026-06-09): the CANONICAL OOF structural feature-schema
    # fingerprint, threaded by the orchestrator from the OOF stage. Forwarded to
    # the engine, which asserts the FINAL refit's structural fingerprint equals
    # it BEFORE training/saving (ERROR on mismatch). NULL = FINAL runs standalone
    # (still computes + persists its own fingerprint).
    canonical_oof_fingerprint = NULL,
    overwrite = TRUE,
    verbose = TRUE,
    # Precision 1 (2026-06-07): internal sentinel set TRUE by the orchestrator,
    # which already resolved the deprecated methodological shims at the
    # run_oneyear_supervised_pipeline() boundary. When TRUE this standalone
    # boundary skips re-warning (avoids double-warning on the orchestrated path).
    # Direct callers leave it FALSE and get the full deprecated-shim treatment.
    .internal_resolved = FALSE) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(train_features) || is.null(train_features)) {
    stop("'train_features' is required.", call. = FALSE)
  }
  if (!(inherits(train_features, c("sf", "data.frame", "SpatVector")) ||
        (is.character(train_features) && length(train_features) == 1L &&
         file.exists(train_features)))) {
    stop("'train_features' must be sf / data.frame, SpatVector, or an ",
         "existing GPKG path.", call. = FALSE)
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
  # over a canonical default). When invoked by the orchestrator
  # (.internal_resolved=TRUE) the shims were ALREADY resolved at the
  # run_oneyear_supervised_pipeline() boundary, so we use the plain
  # cfg-precedence %||% to avoid a second warning for the same override.
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
  random_to_burned_ratio                 <- .shim(random_to_burned_ratio,                 .tc$caps$random,     "cap_random",     "random_to_burned_ratio")
  otsu_unburned_to_burned_ratio          <- .shim(otsu_unburned_to_burned_ratio,          .tc$caps$otsu,       "cap_otsu",       "otsu_unburned_to_burned_ratio")
  feature_whitelist_override <- .shim(feature_whitelist_override, .tc$feature_whitelist_override, "feature_whitelist_override")
  feature_weights            <- .shim(feature_weights,            .tc$feature_weights, "feature_weights")
  sampling_seed         <- .shim(sampling_seed,         .tc$seeds$final_sampling_seed, "final_sampling_seed", "sampling_seed")
  seed                  <- .shim(seed,                  .tc$seeds$final_seed,          "final_seed",          "seed")
  val_frac              <- .shim(val_frac,              .tc$val_frac,        "val_frac")
  group_col             <- .shim(group_col,            .tc$group_col,       "group_col")
  nrounds_max           <- .shim(nrounds_max,           .tc$nrounds_max,     "nrounds_max")
  early_stopping_rounds <- .shim(early_stopping_rounds, .tc$early_stop,      "early_stop", "early_stopping_rounds")
  impute_numeric        <- .shim(impute_numeric,        .tc$impute_numeric,  "impute_numeric")
  impute_factor_missing <- .shim(impute_factor_missing, .tc$impute_factor_missing, "impute_factor_missing")
  # OPTIONAL shape/size block (OtsuFire 0.12.0): resolve from cfg when NULL (the
  # SAME cfg field run_oof_diagnostics() reads -> OOF/FINAL cannot diverge).
  # Plain cfg-precedence %||% (the FALSE default is the canonical baseline).
  include_shape_features <- include_shape_features %||% .tc$include_shape_features
  include_shape_features <- isTRUE(include_shape_features)
  if (!is.null(oof_agg) &&
      !(is.character(oof_agg) && length(oof_agg) == 1L)) {
    stop("'oof_agg' must be NULL or a single CSV path.", call. = FALSE)
  }
  for (nm in c("random_to_burned_ratio",
               "otsu_unburned_to_burned_ratio")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop(sprintf("'%s' must be a single non-negative number.", nm),
           call. = FALSE)
    }
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve out_dir + prefix from config
  # ---------------------------------------------------------------------------
  if (is.null(out_dir)) {
    out_dir <- config$output_routes$final_model_dir
  }
  if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L ||
      !nzchar(out_dir)) {
    stop("Could not resolve 'out_dir' (07_FINAL_MODEL_V2) from config.",
         call. = FALSE)
  }
  target_year <- config$target_year
  scenario    <- config$scenario
  # The orchestrator's final-stage prefix is
  # sprintf("%d_%s_%s", year, scenario, prefix_base) where prefix_base comes
  # from config$options$prefix_base (default "patch_certified"). NOT the OOF
  # "_patch". B3 (2026-06-06): read the SAME knob here so overriding the option
  # keeps the advertised paths consistent. Default unchanged -> byte-identical.
  prefix_base <- config$options$prefix_base %||% "patch_certified"
  prefix <- sprintf("%d_%s_%s", target_year, scenario, prefix_base)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # 2) Resolve train_features to a GPKG path the engine reads. The engine reads
  #    a labelled GPKG by `labelled_gpkg` + `labelled_layer`. When the caller
  #    supplies an in-memory sf / data.frame we materialise it to a temp GPKG
  #    under `labelled_layer` so the engine path is byte-identical.
  # ---------------------------------------------------------------------------
  # PHASE 2 (artifact_hard): resolve the per-row weighting inputs. source_weights
  # defaults to the cfg value (config$negative_pool_params$source_weights, NULL
  # off by default); a non-NULL argument overrides it. artifact_hard_source
  # defaults to the canonical "artifact_hard" source tag when the feature is
  # enabled, else character(0) (so the OFF path resolves NO weight vector and is
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

  if (is.character(train_features)) {
    labelled_gpkg <- train_features
  } else {
    tf_sf <- if (inherits(train_features, "SpatVector")) {
      sf::st_as_sf(train_features)
    } else if (inherits(train_features, "sf")) {
      train_features
    } else {
      sf::st_as_sf(as.data.frame(train_features))
    }
    labelled_gpkg <- tempfile(fileext = ".gpkg")
    sf::st_write(tf_sf, labelled_gpkg, layer = labelled_layer,
                 quiet = TRUE, delete_dsn = TRUE)
  }

  # ---------------------------------------------------------------------------
  # 3) Train (faithful MOVE of the fused wrapper's train_final_model_direct
  #    call). Every argument/value matches the historical wrapper -> orchestrator
  #    call: the two caps (random + Otsu residual), the whitelist override +
  #    weights, the `qa` OOF
  #    aggregate, overwrite/verbose. The engine's own sampling/training seeds
  #    (sampling_seed = 999, seed = 999) are its untouched defaults -> identical
  #    RNG stream -> byte-identical model.
  # ---------------------------------------------------------------------------
  m2 <- train_final_model_direct(
    qa             = oof_agg,
    labelled_gpkg  = labelled_gpkg,
    labelled_layer = labelled_layer,
    out_dir        = out_dir,
    prefix         = prefix,
    overwrite      = overwrite,
    verbose        = verbose,
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    feature_whitelist_override             = feature_whitelist_override,
    feature_weights                        = feature_weights,
    # OPTIONAL shape/size block (OtsuFire 0.12.0): forward the SAME flag the OOF
    # stage uses so the FINAL whitelist universe matches OOF exactly.
    include_shape_features                 = include_shape_features,
    # PHASE 2 (artifact_hard): forward the per-row weighting inputs. Default-off
    # (NULL source_weights + character(0) artifact_hard_source) -> byte-identical.
    source_weights                         = source_weights,
    artifact_hard_source                   = artifact_hard_source,
    total_weight_ratio                     = total_weight_ratio,
    # 2026-06-05 (D1 expose): forward the FINAL training knobs. Defaults equal
    # the engine's historical hardcoded values -> byte-identical RNG/model.
    sampling_seed                          = sampling_seed,
    seed                                   = seed,
    val_frac                               = val_frac,
    group_col                              = group_col,
    nrounds_max                            = nrounds_max,
    early_stopping_rounds                  = early_stopping_rounds,
    impute_numeric                         = impute_numeric,
    impute_factor_missing                  = impute_factor_missing,
    # Gate 1E (2026-06-09): the canonical OOF structural fingerprint to assert
    # the FINAL refit against (pre-train/save). NULL = standalone FINAL.
    canonical_oof_fingerprint              = canonical_oof_fingerprint,
    # Gate 1B (2026-06-07): hand cfg$model_params (the single source of truth,
    # WITHOUT scale_pos_weight) to the engine, which merges the site-specific
    # spw computed from its own training split. The engine no longer carries
    # its own methodological xgb defaults.
    model_params_base                      = config$model_params
  )

  # ---------------------------------------------------------------------------
  # 4) Assemble the documented return (objects + paths). The path fields mirror
  #    what the engine wrote (its `files` list), recovered from the canonical
  #    prefix when absent. The kept objects (model / recipe) are exactly what
  #    the scoring half consumes.
  # ---------------------------------------------------------------------------
  files <- m2$files %||% list()

  final_model_rds <- files$model_rds %||%
    file.path(out_dir, paste0(prefix, "_final_model.rds"))
  recipe_rds <- files$recipe_rds %||%
    file.path(out_dir, paste0(prefix, "_recipe.rds"))
  split_idx_rds <- files$split_rds %||%
    file.path(out_dir, paste0(prefix, "_split_idx.rds"))

  recipe <- if (!is.null(recipe_rds) && file.exists(recipe_rds)) {
    readRDS(recipe_rds)
  } else {
    NULL
  }

  list(
    model                  = m2$model,
    recipe                 = recipe,
    training_ok            = m2$training_ok_sf,
    split_idx              = m2$split,
    final_model_rds        = final_model_rds,
    recipe_rds             = recipe_rds,
    feature_importance_csv = files$feature_importance_csv %||%
      file.path(out_dir, paste0(prefix, "_feature_importance.csv")),
    model_summary_txt      = files$model_summary_txt %||%
      file.path(out_dir, paste0(prefix, "_model_summary.txt")),
    meta_txt               = files$meta_txt %||%
      file.path(out_dir, paste0(prefix, "_meta.txt")),
    split_idx_rds          = split_idx_rds,
    training_ok_csv        = files$training_ok_csv %||%
      file.path(out_dir, paste0(prefix, "_training_ok.csv")),
    training_ok_gpkg       = files$training_ok_gpkg %||%
      file.path(out_dir, paste0(prefix, "_training_ok.gpkg")),
    # Gate 1E (2026-06-09): the FINAL structural feature-schema fingerprint
    # (also persisted in the recipe under recipe$schema_fingerprint).
    schema_fingerprint     = m2$schema_fingerprint,
    # PHASE 2 (artifact_hard): the per-source sample-weight log (n / per-row /
    # total weight per source). Empty 0-row frame on the OFF path.
    sample_weight_log      = m2$sample_weight_log,
    direct                 = m2
  )
}
