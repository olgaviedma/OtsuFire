#' Run out-of-fold (OOF) diagnostics for the supervised burned-area model
#'
#' @description
#' Cross-validates the burned-area model on the spatial folds and reports how
#' well it separates burned from unburned at every probability threshold. This
#' is the diagnostics stage of the supervised pipeline: it tells you which score
#' threshold to use and how the model is expected to perform before you train
#' and apply the final model.
#'
#' Run it after the features stage and before [train_final_burned_model()]. Pass
#' the labelled and scoring features and the configuration from
#' [build_supervised_burned_config()]; the methodological settings (XGBoost
#' parameters, negative-pool caps, feature whitelist/weights, seeds) are read
#' from the configuration, so calling this function directly produces the same
#' result as the equivalent stage of the full pipeline
#' [run_oneyear_supervised_pipeline()].
#'
#' @section What it does:
#' \enumerate{
#'   \item Builds the XGBoost design matrix from the extracted features.
#'   \item Runs repeated spatially blocked out-of-fold scoring (each fold is
#'     predicted by a model that never saw it).
#'   \item Computes per-threshold OOF metrics and picks the recommended and best
#'     thresholds.
#'   \item Writes the `05_OOF` metric outputs and the `04_MATRIX` design-matrix
#'     bundle.
#' }
#'
#' @section Outputs:
#' Files are prefixed `"<year>_<scenario>_patch"`:
#' \itemize{
#'   \item `<prefix>_oof_agg.csv`, `<prefix>_oof_long.csv` — per-unit and
#'     per-fold OOF scores.
#'   \item `<prefix>_oof_metrics_by_threshold.csv`,
#'     `<prefix>_oof_metrics_summary.txt`, `<prefix>_oof_best_thresholds.csv` —
#'     the threshold sweep and chosen thresholds.
#'   \item `<prefix>_labeled_oof_summary.gpkg` — OOF scores joined back to
#'     geometry (when a labelled GPKG is available).
#'   \item the design-matrix bundle under `04_MATRIX`.
#' }
#'
#' @details
#' The OOF stage and the final-model stage share one training core, so their
#' diagnostics and the fitted model cannot diverge. With `params = NULL` the
#' XGBoost parameter block is the canonical one used by both stages (booster
#' `gbtree`, objective `binary:logistic`, `eval_metric` of `logloss` first —
#' which drives early stopping — then `aucpr`, `eta = 0.05`, `max_depth = 5`,
#' `min_child_weight = 5`, `subsample = 0.8`, `colsample_bytree = 0.75`,
#' `gamma = 0`, `lambda = 1`, `alpha = 0`), with `scale_pos_weight` computed
#' here as `n_unburned / n_burned` from the OOF training labels.
#'
#' `feature_whitelist_override` and `feature_weights` are forwarded so the OOF
#' stage trains under the same feature space and per-feature weights as the
#' final model; this is why both are part of the signature.
#'
#' @param train_features sf / data.frame OR a single GPKG path. The labelled
#'   training features produced by the features stage (the `train_features`
#'   object / the `train_features` layer of
#'   `03_FEATURES/features_geometry.gpkg`). Must carry the `class` column and
#'   the fold columns (`fold_cols`). When a path is supplied the
#'   `train_features` layer is read.
#' @param scoring_features sf / data.frame OR a single GPKG path. The scoring
#'   (burned-like) features (`scoring_features` object / layer). When a path is
#'   supplied the `scoring_features` layer is read.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Used to derive `out_dir` (the
#'   `05_OOF` folder), `matrix_dir` (the `04_MATRIX` folder), `result_dir`
#'   (the SUPERVISED scenario base), the output prefix
#'   (`"<year>_<scenario>_patch"`), and `target_year` / `scenario`.
#' @param fold_cols Character vector of fold column names. Default
#'   `c("fold_rep1", "fold_rep2")`.
#' @param params Optional named list of XGBoost params. When `NULL` (default)
#'   the canonical parameter block is built (see Details).
#' @param feature_whitelist_override Character vector restricting the active OOF
#'   feature space to a subset of the canonical whitelist, or `NULL` (default)
#'   for the full whitelist. Forwarded to [run_dm_oof_pipeline()].
#' @param feature_weights Named numeric vector of per-feature weights, or
#'   `NULL` (default). Forwarded to [run_dm_oof_pipeline()].
#' @param include_shape_features Logical or `NULL`. Optional shape/size feature
#'   block. `NULL` (default) reads `cfg$train_control$include_shape_features` —
#'   the same field the final stage reads, so OOF and final cannot diverge on
#'   the active feature set. When effectively `TRUE` the OOF whitelist universe
#'   gains the six shape names and their `_isNA` companions; when `FALSE` the
#'   OOF output is unchanged.
#' @param source_weights Named numeric vector of per-row weights keyed by the
#'   `source` column, or `NULL` (default). `NULL` reads the configured value
#'   (`config$negative_pool_params$source_weights`, off by default). With no
#'   override and no artifact_hard rows present, no weight vector is built.
#' @param artifact_hard_source Character vector of `source` value(s) marking
#'   promoted artifact_hard rows, or `NULL` (default). `NULL` resolves to
#'   `"artifact_hard"` when that feature is enabled, else `character(0)` (off).
#'   Used by the per-fold eligibility resolver (promoted rows enter the uncapped
#'   artifact_hard bucket) and the per-row sample-weight resolver.
#' @param total_weight_ratio Numeric or `NULL`. Pool-level weight balance
#'   (artifact_hard total weight / burned total weight). `NULL` (default) reads
#'   `config$negative_pool_params$artifact_hard$total_weight_ratio`.
#' @param nrounds_max Integer or `NULL`. Maximum XGBoost boosting rounds for the
#'   OOF per-fold models. `NULL` (default) reads the configured value (matches
#'   the final-model default).
#' @param early_stop Integer or `NULL`. XGBoost early-stopping patience for the
#'   OOF per-fold models. `NULL` (default) reads the configured value (matches
#'   the final stage).
#' @param seed_base Integer or `NULL`. Base RNG seed for the repeated block-CV
#'   OOF runs. `NULL` (default) reads the configured value (matches the final
#'   stage).
#' @param out_dir Character or `NULL`. Output folder for the `05_OOF` outputs.
#'   Defaults to `config$output_routes$oof_dir`.
#' @param matrix_dir Character or `NULL`. Output folder for the `04_MATRIX`
#'   design-matrix bundle. Defaults to `config$output_routes$matrix_dir`.
#' @param labelled_gpkg Character or `NULL`. Path to the train-with-folds GPKG
#'   used to attach geometry to the OOF aggregate for the
#'   `labeled_oof_summary` layer. When `NULL` no labelled-summary GPKG is
#'   written by the inner wrapper (the orchestrator always passes the runtime
#'   `02_FOLDS` train-with-folds GPKG; its block size is not fixed, so it is an
#'   explicit argument rather than a config route).
#' @param labelled_layer Character. Layer name inside `labelled_gpkg`. Default
#'   `"train_with_folds"`.
#' @details
#' OOF uses the same capped negative-sampling policy as the final model, applied
#' independently within each training fold. Each outer fold runs through the
#' shared leakage-free core: per-fold medians and `scale_pos_weight` are fit on
#' the outer-train rows only, an inner-validation split is the sole
#' early-stopping set, the model is refit on all outer-train rows at the best
#' iteration, and the untouched outer-test fold is predicted.
#'
#' Training eligibility is defined by explicit class, never by negation: only
#' explicit burned rows (positives) and explicit unburned rows that resolve to a
#' valid negative bucket (random, otsu) enter training; review / keep / `NA` /
#' unknown rows never become negatives. The OOF and final stages share one
#' eligibility resolver and one capping helper, so they cannot diverge on which
#' rows are used or how negatives are capped.
#' @param random_to_burned_ratio,otsu_unburned_to_burned_ratio
#'   Numeric or `NULL`. The two negative-bucket caps forwarded to the OOF chain
#'   so OOF sees the same caps as the final model. `NULL` (default) reads the
#'   configured caps.
#' @param val_frac Numeric or `NULL`. Inner validation fraction for the
#'   per-fold split. `NULL` (default) reads the configured value.
#' @param impute_numeric Character or `NULL`. Numeric-imputation rule
#'   (`"median"` or `"zero"`). `NULL` (default) reads the configured value.
#' @param impute_factor_missing Character or `NULL`. Sentinel level for missing
#'   factor/character values. `NULL` (default) reads the configured value.
#' @param group_col Character or `NULL`. Grouping column for the grouped
#'   block-CV. `NULL` (default) reads `cfg$train_control$group_col` (the same
#'   field the final stage reads, so OOF and final never diverge).
#' @param overwrite Logical. Forwarded to [run_dm_oof_pipeline()] (controls
#'   whether the design-matrix bundle is recomputed/clobbered). Default `TRUE`.
#' @param .internal_resolved Internal use only; set by the orchestrator. When
#'   `TRUE` the deprecated methodological shims were already resolved upstream,
#'   so this boundary skips re-warning. Direct callers leave it `FALSE`.
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `design_bundle` — the design-matrix bundle (`dm`) returned by the
#'       inner wrapper (the object kept for the downstream final-model stage).
#'     \item `oof_agg` — the per-unit OOF aggregate data.frame.
#'     \item `oof_long` — the per-fold long OOF data.frame.
#'     \item `labeled_oof_summary` — the OOF aggregate joined back to geometry
#'       (sf), or `NULL` when `labelled_gpkg` is not available.
#'     \item `oof_agg_csv`, `oof_long_csv`, `oof_metrics_by_threshold_csv`,
#'       `oof_metrics_summary_txt`, `oof_best_thresholds_csv`,
#'       `labeled_oof_summary_gpkg`, `design_bundle_rds` — written paths.
#'   }
#'
#' @section Deprecated function-level parameter shims:
#' The methodological / training-control arguments here (`nrounds_max`,
#' `early_stop`, `seed_base`, the `*_to_burned_ratio` caps,
#' `feature_whitelist_override`, `feature_weights`, `val_frac`, `impute_*`,
#' `group_col`) are deprecated compatibility shims. Set them in
#' [build_supervised_burned_config()] instead (`cfg$train_control`), which every
#' stage reads. A non-`NULL` override of a canonical-default field emits a
#' deprecation warning of class `"otsufire_deprecated_param"`; an override that
#' conflicts with an explicit builder value errors. These shims are scheduled
#' for removal in a future minor version.
#'
#' @seealso
#' [build_supervised_burned_config()], [extract_supervised_features()],
#' [train_final_burned_model()], [score_supervised_burned_map()],
#' [make_spatial_folds()], [validate_supervised_execution()],
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
#' feats <- extract_supervised_features(
#'   train_with_folds = "02_FOLDS/2017_train_with_folds_2000m.gpkg",
#'   scoring_pool     = "01_POOLS/2017_balanced_pools.gpkg",
#'   config = cfg
#' )
#' oof <- run_oof_diagnostics(
#'   train_features   = feats$features_geometry_gpkg,
#'   scoring_features = feats$features_geometry_gpkg,
#'   config = cfg
#' )
#' oof$oof_agg_csv
#' }
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
  #    prefix_oof == "<year>_<scenario>_patch").
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
