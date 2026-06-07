#' Train the final supervised burned-area model
#'
#' @description
#' Supervised-pipeline stage D (final-model training). This is the TRAINING
#' half of the formerly fused
#' `run_train_final_model_and_export_final_map()` wrapper. It wraps the
#' internal engine [train_final_model_direct()] with the exact arguments and
#' values the fused wrapper passed (the four sampling caps, the active feature
#' whitelist / per-feature weights, the OOF aggregate `qa` input, and the
#' historical sampling/training seeds), writes the canonical `07_FINAL_MODEL_V2`
#' artifacts under the `"<year>_<scenario>_patch_certified"` prefix, and returns
#' both the in-memory objects (`model`, `recipe`, `training_ok`, `split_idx`)
#' and the written paths. The `score_supervised_burned_map()` scoring half then
#' consumes the returned `model` + `recipe`.
#'
#' This is the standalone, exported implementation of the final-model stage that
#' the one-year orchestrator [run_oneyear_supervised_pipeline()] delegates to,
#' so calling it directly produces byte-identical final-model outputs to a full
#' run.
#'
#' @details
#' The four `*_to_burned_ratio` caps, `feature_whitelist_override`, and
#' `feature_weights` are forwarded verbatim to [train_final_model_direct()].
#' Their defaults reproduce the historical behaviour byte-for-byte. The
#' whitelist override / weights are part of the signature for the same reason
#' the OOF stage carries them: the final model must train under the SAME feature
#' space + per-feature weights as the OOF diagnostics (KB1/KB2 symmetry).
#'
#' The `oof_agg` argument is the OOF per-unit aggregate CSV (`_oof_agg.csv`,
#' the `qa` argument of the engine). It is consumed only to enrich the
#' `model_summary.txt` with the OOF recommended-threshold block; it does not
#' change the fitted model. The engine derives the
#' `_oof_metrics_summary.txt` / `_oof_best_thresholds.csv` paths from the
#' `oof_agg` path by name substitution, exactly as before.
#'
#' @param train_features sf / data.frame OR a single GPKG path. The labelled
#'   training features produced by the features stage (the `train_features`
#'   layer of `03_FEATURES/features_geometry.gpkg`). Must carry the `class`,
#'   `source`, `neg_type` and feature columns. When a path is supplied the
#'   `labelled_layer` layer is read by the engine; when an in-memory sf /
#'   data.frame is supplied it is written to a temporary GPKG so the engine can
#'   read it back identically.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Used to derive `out_dir` (the
#'   `07_FINAL_MODEL_V2` folder), the output prefix
#'   (`"<year>_<scenario>_patch_certified"`), and `target_year` / `scenario`.
#' @param oof_agg Character path to the OOF aggregate CSV (`_oof_agg.csv`) OR
#'   `NULL`. Forwarded as the engine's `qa` argument (summary enrichment only).
#' @param contextual_exclusion_to_burned_ratio Numeric. Cap on the
#'   contextual-exclusion negative pool. Default `0.25`.
#' @param spectral_hard_negative_to_burned_ratio Numeric. Cap on the
#'   spectral-hard-negative pool. Default `1.0`.
#' @param random_to_burned_ratio Numeric. Cap on random-burnable-background
#'   negatives. Default `1.0`.
#' @param otsu_unburned_to_burned_ratio Numeric. Cap on Otsu current-year
#'   unburned-patch negatives. Default `1.0`.
#' @param feature_whitelist_override Character vector subset of the canonical
#'   `.supervised_feature_cols` to use as the active feature space, or `NULL`
#'   (default) for the full whitelist. Forwarded to
#'   [train_final_model_direct()].
#' @param feature_weights Named numeric vector of per-feature weights, or
#'   `NULL` (default). Forwarded to [train_final_model_direct()].
#' @param sampling_seed Integer. RNG seed for the negative-pool sampling step
#'   in [train_final_model_direct()]. Default `42`. FRENTE 1 (2026-06-05):
#'   UNIFIED with the OOF canonical seed (was `999`). User-overridable.
#' @param seed Integer. RNG seed for the train/val split and the xgboost
#'   training call in [train_final_model_direct()]. Default `42`. FRENTE 1:
#'   UNIFIED with the OOF canonical seed (was `999`). User-overridable.
#' @param val_frac Numeric in (0, 1). Validation fraction for the final-model
#'   train/val split. Default `0.15`.
#' @param group_col Character or `NULL`. Grouping column for the grouped
#'   train/val split (kept-together unit). Default `"block_id"`.
#' @param nrounds_max Integer. Maximum xgboost boosting rounds for the FINAL
#'   model. Default `4000`. FRENTE 1 (2026-06-05): the OOF stage default was
#'   raised `3000` -> `4000` to MATCH this; the two are now unified.
#' @param early_stopping_rounds Integer. xgboost early-stopping patience for
#'   the FINAL model. Default `80`. FRENTE 1: the OOF stage default was raised
#'   `75` -> `80` to MATCH this; the two are now unified.
#' @param impute_numeric Character. Numeric-imputation rule passed to
#'   [train_final_model_direct()]: `"median"` (default) or `"zero"`.
#' @param impute_factor_missing Character. Sentinel level used for missing
#'   factor/character values. Default `"MISSING"`.
#' @param training_protocol Character. `"legacy"` (default) or `"nested_refit"`.
#'   Forwarded to [train_final_model_direct()]. `"legacy"` reproduces the
#'   historical FINAL model byte-for-byte; `"nested_refit"` routes through the
#'   shared leakage-free core (B1, 2026-06-07) and DEPLOYS the model refit on
#'   all training rows at the selected `best_iteration`.
#' @param labelled_layer Character. Layer name read from `train_features` when
#'   it is a GPKG path. Default `"train_features"`.
#' @param out_dir Character or `NULL`. Output folder for the
#'   `07_FINAL_MODEL_V2` artifacts. Defaults to
#'   `config$output_routes$final_model_dir`.
#' @param overwrite Logical. Forwarded to [train_final_model_direct()]
#'   (controls clobbering of the training-ok GPKG). Default `TRUE` (historical
#'   behaviour).
#' @param verbose Logical. Forwarded to [train_final_model_direct()]. Default
#'   `TRUE`.
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
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   scenario = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L
#' )
#' tm <- train_final_burned_model(
#'   train_features = "03_FEATURES/features_geometry.gpkg",
#'   config = cfg,
#'   oof_agg = "05_OOF/2017_balanced_patch_oof_agg.csv"
#' )
#' tm$final_model_rds
#' }
train_final_burned_model <- function(
    train_features,
    config,
    oof_agg = NULL,
    # Gate 1B (2026-06-07): all methodological knobs DEFAULT TO NULL = "read
    # from cfg$train_control / cfg$model_params" (single source of truth). A
    # non-NULL value overrides cfg for standalone use. No literal methodological
    # numbers live here.
    contextual_exclusion_to_burned_ratio   = NULL,
    spectral_hard_negative_to_burned_ratio = NULL,
    random_to_burned_ratio                 = NULL,
    otsu_unburned_to_burned_ratio          = NULL,
    feature_whitelist_override = NULL,
    feature_weights = NULL,
    sampling_seed = NULL,
    seed = NULL,
    val_frac = NULL,
    group_col = NULL,
    nrounds_max = NULL,
    early_stopping_rounds = NULL,
    impute_numeric = NULL,
    impute_factor_missing = NULL,
    training_protocol = NULL,
    labelled_layer = "train_features",
    out_dir = NULL,
    overwrite = TRUE,
    verbose = TRUE) {
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
  # Gate 1B: resolve methodological params from cfg (explicit arg > cfg).
  .tc <- config$train_control
  if (is.null(.tc) || !is.list(.tc)) {
    stop("'config' has no resolved train_control; rebuild it with ",
         "build_supervised_burned_config().", call. = FALSE)
  }
  contextual_exclusion_to_burned_ratio   <- contextual_exclusion_to_burned_ratio   %||% .tc$caps$contextual
  spectral_hard_negative_to_burned_ratio <- spectral_hard_negative_to_burned_ratio %||% .tc$caps$spectral
  random_to_burned_ratio                 <- random_to_burned_ratio                 %||% .tc$caps$random
  otsu_unburned_to_burned_ratio          <- otsu_unburned_to_burned_ratio          %||% .tc$caps$otsu
  feature_whitelist_override <- feature_whitelist_override %||% .tc$feature_whitelist_override
  feature_weights            <- feature_weights            %||% .tc$feature_weights
  sampling_seed         <- sampling_seed         %||% .tc$seeds$final_sampling_seed
  seed                  <- seed                  %||% .tc$seeds$final_seed
  val_frac              <- val_frac              %||% .tc$val_frac
  group_col             <- group_col             %||% .tc$group_col
  nrounds_max           <- nrounds_max           %||% .tc$nrounds_max
  early_stopping_rounds <- early_stopping_rounds %||% .tc$early_stop
  impute_numeric        <- impute_numeric        %||% .tc$impute_numeric
  impute_factor_missing <- impute_factor_missing %||% .tc$impute_factor_missing
  training_protocol     <- training_protocol     %||% .tc$training_protocol
  training_protocol <- match.arg(training_protocol, c("legacy", "nested_refit"))
  if (!is.null(oof_agg) &&
      !(is.character(oof_agg) && length(oof_agg) == 1L)) {
    stop("'oof_agg' must be NULL or a single CSV path.", call. = FALSE)
  }
  for (nm in c("contextual_exclusion_to_burned_ratio",
               "spectral_hard_negative_to_burned_ratio",
               "random_to_burned_ratio",
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
  #    call: the four caps, the whitelist override + weights, the `qa` OOF
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
    contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
    spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    feature_whitelist_override             = feature_whitelist_override,
    feature_weights                        = feature_weights,
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
    training_protocol                      = training_protocol,
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
    direct                 = m2
  )
}
