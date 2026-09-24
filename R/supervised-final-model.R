#' Train the final supervised burned-area model
#'
#' @description
#' Trains the final XGBoost burned-area model on all labelled training rows and
#' writes the model artifacts. This is the training stage of the supervised
#' pipeline; the companion [score_supervised_burned_map()] then uses the
#' returned `model` and `recipe` to score the candidate polygons.
#'
#' Use it after the features and out-of-fold (OOF) stages: pass the labelled
#' training features, the configuration from [build_supervised_burned_config()],
#' and (optionally) the OOF aggregate CSV. The methodological settings (caps,
#' seeds, feature whitelist/weights, XGBoost rounds, imputation rules) are read
#' from the configuration, so calling this function directly produces the same
#' result as the equivalent stage of the full pipeline
#' [run_oneyear_supervised_pipeline()].
#'
#' It returns the in-memory objects (`model`, `recipe`, `training_ok`,
#' `split_idx`) together with the paths it wrote under the `07_FINAL_MODEL_V2`
#' folder.
#'
#' @section How the model is fit:
#' OtsuFire uses one training procedure with no protocol choice. The number of
#' boosting rounds is selected on an inner train/validation split with XGBoost
#' early stopping, then a fresh model is refit on \emph{all} labelled rows at
#' that round count. The saved model, recipe, feature order and parameters are
#' always the refit model. The OOF stage applies the same procedure per fold, so
#' OOF diagnostics and the final model share an identical training core.
#'
#' Training eligibility is defined by explicit class, never by negation: only
#' explicit burned rows (positives) and explicit unburned rows that resolve to a
#' valid negative bucket (`random`, `otsu`, and `artifact_hard` when enabled)
#' enter training; review / keep / `NA` /
#' unknown rows never become negatives. The OOF and final stages share one
#' eligibility resolver and one capping helper, so they cannot diverge on which
#' rows are used or how negatives are capped.
#'
#' The `oof_agg` CSV (when supplied) only enriches `model_summary.txt` with the
#' OOF recommended-threshold block; it does not change the fitted model.
#'
#' @section Configuration is the source of truth:
#' The methodological / training-control arguments here (the `*_to_burned_ratio`
#' caps, `feature_whitelist_override`, `feature_weights`, `sampling_seed`,
#' `seed`, `val_frac`, `group_col`, `nrounds_max`, `early_stopping_rounds`,
#' `impute_*`) are deprecated compatibility shims. Set them in
#' [build_supervised_burned_config()] instead (`cfg$train_control` /
#' `cfg$model_params`), which every stage reads. Each argument defaults to
#' `NULL`, meaning "use the configured value". A non-`NULL` override of a
#' canonical-default field emits a deprecation warning of class
#' `"otsufire_deprecated_param"`; an override that conflicts with an explicit
#' builder value errors. These shims are scheduled for removal in a future
#' minor version.
#'
#' @param train_features sf / data.frame OR a single GPKG path. The labelled
#'   training features produced by the features stage (the `train_features`
#'   layer of `03_FEATURES/features_geometry.gpkg`). Must carry the `class`,
#'   `source`, `neg_type` and feature columns. A path is read using
#'   `labelled_layer`; an in-memory object is written to a temporary GPKG and
#'   read back identically.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Provides the output folder, the file
#'   prefix, `target_year` / `scenario`, and the resolved methodological
#'   parameters.
#' @param oof_agg Character path to the OOF aggregate CSV (`_oof_agg.csv`) OR
#'   `NULL`. Used for summary enrichment only.
#' @param random_to_burned_ratio Numeric or `NULL`. Cap on
#'   random-burnable-background negatives, as a multiple of the burned-label
#'   count. `NULL` (default) reads the configured cap.
#' @param otsu_unburned_to_burned_ratio Numeric or `NULL`. Cap on Otsu
#'   current-year unburned-patch negatives. `NULL` (default) reads the
#'   configured cap.
#' @param feature_whitelist_override Character vector restricting the active
#'   feature space to a subset of the canonical whitelist, or `NULL` (default)
#'   to use the configured whitelist.
#' @param feature_weights Named numeric vector of per-feature weights, or `NULL`
#'   (default) to use the configured weights.
#' @param include_shape_features Logical or `NULL`. Optional shape/size feature
#'   block. `NULL` (default) reads `cfg$train_control$include_shape_features`
#'   (the same field the OOF stage reads, so the two cannot diverge). When
#'   effectively `TRUE` the feature universe gains the six geometric columns and
#'   their `_isNA` companions; when `FALSE` the model is unchanged.
#' @param source_weights Named numeric vector of per-row weights keyed by the
#'   `source` column (e.g. `c(artifact_hard = 0.3)`), or `NULL` (default) to
#'   read `config$negative_pool_params$source_weights` (off by default). With no
#'   override and no artifact_hard rows present, no weight vector is built.
#' @param artifact_hard_source Character vector of `source` value(s) marking
#'   promoted artifact_hard rows, or `NULL` (default). `NULL` resolves to
#'   `"artifact_hard"` when the artifact_hard feature is enabled, else
#'   `character(0)` (off). Used by the eligibility and per-row weight resolvers.
#' @param total_weight_ratio Numeric or `NULL`. Pool-level weight balance
#'   (artifact_hard total weight / burned total weight). `NULL` (default) reads
#'   `config$negative_pool_params$artifact_hard$total_weight_ratio`. Ignored
#'   when an explicit artifact_hard pin is supplied via `source_weights`.
#' @param sampling_seed Integer or `NULL`. RNG seed for the negative-pool
#'   sampling step. `NULL` (default) reads the configured seed.
#' @param seed Integer or `NULL`. RNG seed for the train/val split and the
#'   XGBoost training call. `NULL` (default) reads the configured seed.
#' @param val_frac Numeric in (0, 1) or `NULL`. Validation fraction for the
#'   final-model train/val split. `NULL` (default) reads the configured value.
#' @param group_col Character or `NULL`. Grouping column for the grouped
#'   train/val split (kept-together unit). `NULL` (default) reads the configured
#'   value.
#' @param nrounds_max Integer or `NULL`. Maximum XGBoost boosting rounds. `NULL`
#'   (default) reads the configured value.
#' @param early_stopping_rounds Integer or `NULL`. XGBoost early-stopping
#'   patience. `NULL` (default) reads the configured value.
#' @param impute_numeric Character or `NULL`. Numeric-imputation rule
#'   (`"median"` or `"zero"`). `NULL` (default) reads the configured value.
#' @param impute_factor_missing Character or `NULL`. Sentinel level for missing
#'   factor/character values. `NULL` (default) reads the configured value.
#' @param labelled_layer Character. Layer name read from `train_features` when
#'   it is a GPKG path. Default `"train_features"`.
#' @param out_dir Character or `NULL`. Output folder for the
#'   `07_FINAL_MODEL_V2` artifacts. Defaults to
#'   `config$output_routes$final_model_dir`.
#' @param canonical_oof_fingerprint For internal use by
#'   [run_oneyear_supervised_pipeline()], which passes the OOF feature-schema
#'   fingerprint so the final model is checked against it. Leave it `NULL`
#'   (the default).
#' @param overwrite Logical. Controls clobbering of the training-ok GPKG.
#'   Default `TRUE`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @param .internal_resolved For internal use by
#'   [run_oneyear_supervised_pipeline()]. Leave it `FALSE` (the default).
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `model` — the fitted xgboost model.
#'     \item `recipe` — the training recipe list (impute medians, feature/x
#'       columns, params, applied whitelist/weights).
#'     \item `training_ok` — the selected training-ok sf (the
#'       `training_ok.gpkg` contents).
#'     \item `split_idx` — the train/val split index list.
#'     \item `final_model_rds`, `recipe_rds`, `feature_importance_csv`,
#'       `model_summary_txt`, `meta_txt`, `split_idx_rds`, `training_ok_csv`,
#'       `training_ok_gpkg` — written paths.
#'   }
#'
#' @seealso
#' [build_supervised_burned_config()], [run_oof_diagnostics()],
#' [score_supervised_burned_map()], [validate_supervised_execution()],
#' [run_oneyear_supervised_pipeline()]
#'
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L
#' )
#' tm <- train_final_burned_model(
#'   train_features = "03_FEATURES/features_geometry.gpkg",
#'   config = cfg,
#'   oof_agg = "05_OOF/2017_balanced_patch_oof_agg.csv"
#' )
#' tm$final_model_rds
#' tm$model
#' }
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
