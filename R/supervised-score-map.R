#' Score Otsu-guided patches and export probabilistic refinement burned-area outputs
#'
#' @description
#' Apply a trained supervised model to the Otsu-guided patches in the scoring
#' pool and export their scores and spatial outputs.
#'
#' The function uses the model and feature-processing recipe returned by
#' [train_final_burned_model()]. It predicts `p_burned`, joins available
#' out-of-fold (OOF) information, and applies the current-year temporal
#' adjustment.
#'
#' Outputs include the complete scored candidate layer and a public map layer
#' that excludes candidates flagged by the current-year temporal filter. An
#' additional burned-like subset can be exported.
#'
#' Run this function after final-model training. It performs the scoring
#' stage used by [run_oneyear_supervised_pipeline()].
#'
#' @param scoring_features `sf` object, `data.frame`, or GeoPackage path
#'   containing features for the Otsu-guided patches to score. When a path is
#'   supplied, `scoring_layer` identifies the layer to read. Spatial exports
#'   require polygon geometry.
#' @param model Fitted XGBoost model object or path to its RDS file. Use the
#'   final model returned by [train_final_burned_model()].
#' @param recipe Training recipe list or path to its RDS file. Required. Must
#'   contain `recipe$cols$feature_cols` and the feature-processing information
#'   associated with `model`. Use a model and recipe saved from the same
#'   training run.
#' @param config An `otsufire_supervised_burned_config` object created by
#'   [build_supervised_burned_config()]. Supplies run identifiers and output
#'   locations.
#' @param oof_summary Path to the OOF labelled-summary GeoPackage, or `NULL`.
#'   Supplies OOF information for the `p_oof_mean` join. When `NULL`, the
#'   function uses the conventional path under `05_OOF`.
#' @param labelled_features `sf` object, `data.frame`, GeoPackage path, or
#'   `NULL`. Labelled training features used to associate OOF predictions with
#'   scoring polygons through `source_poly_id`. When `NULL`, uses
#'   `config$output_routes$features_geometry_gpkg`.
#' @param preyear_overlap_threshold Numeric scalar. Previous-year overlap
#'   threshold used to identify temporal conflicts. Default: `0.70`.
#' @param hotspot_density_threshold Numeric scalar. Hotspot-density threshold
#'   used to assess current-year fire support during temporal adjustment.
#'   Default: `0.001`.
#' @param temporal_penalty_floor Numeric scalar. Floor parameter used by the
#'   temporal score-adjustment rule. Default: `0.10`. These three arguments
#'   control temporal adjustment. They do not specify the score threshold for
#'   producing a binary burned-area map.
#' @param export_burned_like Logical scalar. Whether to export the burned-like
#'   subset and its counts table. Default: `TRUE`.
#' @param qa_labelled_layer Layer name in `oof_summary`. Default:
#'   `"labeled_oof_summary"`.
#' @param labelled_features_layer Layer name read from `labelled_features` when
#'   supplied as a GeoPackage path. Default: `"train_features"`.
#' @param scoring_layer Layer name read from `scoring_features` when supplied
#'   as a GeoPackage path. Default: `"scoring_features"`.
#' @param out_map_dir Output directory for final-map products, or `NULL` to
#'   use `config$output_routes$final_map_dir`, normally the `09_FINAL_MAP`
#'   folder.
#' @param out_score_dir Output directory for scored products, or `NULL` to use
#'   `config$output_routes$scored_dir`, normally the `08_SCORED` folder.
#' @param overwrite Logical scalar. Whether existing outputs may be replaced.
#'   Default: `TRUE`.
#' @param verbose Logical scalar. Whether to print progress messages. Default:
#'   `TRUE`.
#'
#' @section Scoring workflow:
#' The function:
#' 1. validates the supplied inputs and resolves output locations;
#' 2. prepares scoring features using the saved training recipe;
#' 3. builds the scoring matrix and predicts `p_burned`;
#' 4. joins available OOF information using the labelled-feature identifiers;
#' 5. applies the current-year temporal adjustment;
#' 6. writes the complete scored layer, public map layer, and optional
#'    burned-like subset.
#'
#' `p_burned` is a model score and is not necessarily a calibrated
#' probability.
#'
#' @section Feature preparation:
#' The training recipe defines the feature names, predictor order, expected
#' types, missingness indicators, imputation values, and categorical levels
#' used for scoring.
#'
#' The function requires this saved schema rather than deriving a new one
#' from the scoring data. It runs [validate_supervised_execution()] in strict
#' mode with the model, recipe, and available feature names before
#' constructing the scoring matrix.
#'
#' Recoverable differences are reconciled using the recipe. These may
#' include:
#'
#' | Difference | Treatment |
#' |---|---|
#' | Different predictor order | Reorder predictors to match the training schema. |
#' | Missing predictor | Create the missing column and apply the recipe's missing-value handling. |
#' | Additional predictor | Exclude it from the model matrix when it is outside the saved schema. |
#' | Changed column type | Convert to the expected type where possible; otherwise handle the affected values as missing. |
#' | Unseen categorical level | Map to the recipe's missing or unknown-category representation. |
#' | Entirely missing numeric predictor | Apply the saved numeric imputation rule. |
#'
#' Corrections are recorded by the scoring workflow. An incompatible schema
#' raises an error.
#'
#' The initial check based on column names assesses feature availability.
#' Type conversion and missing-value handling occur when the scoring data are
#' processed.
#'
#' @section OOF information:
#' OOF predictions describe labelled examples scored while held out from
#' their corresponding training folds.
#'
#' The optional OOF join adds this information to matching polygons through
#' the labelled-feature identifiers. It does not make the final-model
#' predictions out-of-fold predictions.
#'
#' Keep OOF scores and final-model scores distinct when evaluating
#' performance, particularly for polygons that also contributed to
#' final-model training.
#'
#' @section Current-year temporal filtering:
#' The temporal rule identifies candidates with a previous-year conflict and
#' weak current-year hotspot support.
#'
#' A temporal conflict is identified when either:
#' * previous-year overlap meets or exceeds `preyear_overlap_threshold`; or
#' * `preyear_action == "drop"`.
#'
#' Candidates meeting the public-map exclusion rule are flagged through
#' `current_year_public_drop`. They remain available in the complete scored
#' layers but are excluded from the public `final_map` layer.
#'
#' This flag represents a rule-based temporal exclusion. It does not
#' independently establish that a polygon is an old fire scar.
#'
#' @section Output layers:
#' The final-map GeoPackage contains three layers:
#'
#' | Layer | Contents |
#' |---|---|
#' | `deterministic_scored` | Complete scored Otsu-guided patch layer, retaining all scored candidates and their identifiers and geometries. |
#' | `final_map_full` | Same contents as `deterministic_scored`. |
#' | `final_map` | Public-column subset after removing candidates flagged by `current_year_public_drop`. |
#'
#' The name `deterministic_scored` is retained for API compatibility.
#'
#' The complete layers contain one row per input scoring polygon. The public
#' layer may contain fewer rows because of temporal exclusions.
#'
#' When the exclusion flag is fully populated with logical values, the
#' expected relationship is:
#' \preformatted{
#' nrow(final_map) ==
#'   nrow(final_map_full) -
#'   sum(final_map_full$current_year_public_drop)
#' }
#'
#' The `final_map` layer is temporally filtered. Its name alone does not
#' imply that a binary threshold has been applied to `p_burned`. Apply the
#' selected score threshold when preparing a binary burned-area product for
#' external validation.
#'
#' @return A named list containing scored spatial objects and output paths.
#'
#' | Field | Contents |
#' |---|---|
#' | `deterministic_scored` | Complete scored Otsu-guided patch layer as an `sf` object. |
#' | `final_map_full` | Same spatial object contents as `deterministic_scored`. |
#' | `final_map` | Public map layer after temporal exclusions and public-column selection. |
#' | `burned_like_scored` | Burned-like subset as an `sf` object, or `NULL` when `export_burned_like = FALSE`. |
#' | `deterministic_scored_gpkg` | Path to the complete scored output. |
#' | `final_map_gpkg` | Path to the final-map GeoPackage. |
#' | `final_map_counts_csv` | Path to the final-map counts table. |
#' | `burned_like_gpkg` | Path to the optional burned-like output. |
#' | `burned_like_counts_csv` | Path to the optional burned-like counts table. |
#'
#' Burned-like export paths do not represent newly written products when
#' `export_burned_like = FALSE`.
#'
#' @seealso [build_supervised_burned_config()],
#'   [extract_supervised_features()], [run_oof_diagnostics()],
#'   [train_final_burned_model()], [validate_supervised_execution()],
#'   [validate_fire_maps()], [run_oneyear_supervised_pipeline()].
#'
#' @examples
#' \dontrun{
#' # Configure the run using the inputs used in the preceding stages
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced",
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022"
#' )
#'
#' # Use the feature and OOF files produced by the preceding stages
#' features_path <- "path/to/features_geometry.gpkg"
#' oof_agg_path <- "path/to/2022_balanced_patch_oof_agg.csv"
#' oof_summary_path <-
#'   "path/to/2022_balanced_patch_labeled_oof_summary.gpkg"
#'
#' # Train the final model
#' trained <- train_final_burned_model(
#'   train_features = features_path,
#'   config = cfg,
#'   oof_agg = oof_agg_path
#' )
#'
#' # Score the Otsu-guided patches
#' scored <- score_supervised_burned_map(
#'   scoring_features = features_path,
#'   model = trained$model,
#'   recipe = trained$recipe,
#'   config = cfg,
#'   oof_summary = oof_summary_path,
#'   labelled_features = features_path,
#'   export_burned_like = TRUE
#' )
#'
#' # Locate the exported map
#' scored$final_map_gpkg
#'
#' # Compare complete and public-layer counts
#' c(
#'   all_candidates = nrow(scored$final_map_full),
#'   public_candidates = nrow(scored$final_map)
#' )
#'
#' # Inspect temporal exclusions
#' table(
#'   scored$final_map_full$current_year_public_drop,
#'   useNA = "ifany"
#' )
#'
#' # Inspect the public-layer scores
#' summary(scored$final_map$p_burned)
#'
#' # Read the exported public layer
#' public_map <- sf::st_read(
#'   scored$final_map_gpkg,
#'   layer = "final_map",
#'   quiet = TRUE
#' )
#'
#' plot(public_map["p_burned"])
#' }
#'
#' @family workflow
#' @export
score_supervised_burned_map <- function(
    scoring_features, model, recipe,
    config,
    oof_summary = NULL,
    labelled_features = NULL,
    export_burned_like = TRUE,
    preyear_overlap_threshold = 0.70,
    hotspot_density_threshold = 0.001,
    temporal_penalty_floor = 0.10,
    qa_labelled_layer = "labeled_oof_summary",
    labelled_features_layer = "train_features",
    scoring_layer = "scoring_features",
    out_map_dir = NULL,
    out_score_dir = NULL,
    overwrite = TRUE,
    verbose = TRUE) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(scoring_features) || is.null(scoring_features)) {
    stop("'scoring_features' is required.", call. = FALSE)
  }
  if (missing(model) || is.null(model)) {
    stop("'model' is required (xgboost model object or RDS path).",
         call. = FALSE)
  }
  if (missing(recipe) || is.null(recipe)) {
    stop("'recipe' is required (named list or RDS path).", call. = FALSE)
  }
  if (missing(config) || is.null(config) ||
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (!is.logical(export_burned_like) || length(export_burned_like) != 1L ||
      is.na(export_burned_like)) {
    stop("'export_burned_like' must be TRUE or FALSE.", call. = FALSE)
  }
  for (nm in c("preyear_overlap_threshold", "hotspot_density_threshold",
               "temporal_penalty_floor")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be a single number.", nm), call. = FALSE)
    }
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve dirs + prefix from config (result_dir == output_routes$base,
  #    out_score_dir == 08_SCORED, out_map_dir == 09_FINAL_MAP, prefix ==
  #    "<year>_<run_label>_patch_certified").
  # ---------------------------------------------------------------------------
  result_dir <- config$output_routes$base
  if (is.null(result_dir) || !is.character(result_dir) ||
      length(result_dir) != 1L || !nzchar(result_dir)) {
    stop("Could not resolve 'result_dir' (output base) from config.",
         call. = FALSE)
  }
  if (is.null(out_map_dir)) {
    out_map_dir <- config$output_routes$final_map_dir
  }
  if (is.null(out_map_dir) || !is.character(out_map_dir) ||
      length(out_map_dir) != 1L || !nzchar(out_map_dir)) {
    stop("Could not resolve 'out_map_dir' (09_FINAL_MAP) from config.",
         call. = FALSE)
  }
  if (is.null(out_score_dir)) {
    out_score_dir <- config$output_routes$scored_dir
  }
  if (is.null(out_score_dir) || !is.character(out_score_dir) ||
      length(out_score_dir) != 1L || !nzchar(out_score_dir)) {
    stop("Could not resolve 'out_score_dir' (08_SCORED) from config.",
         call. = FALSE)
  }
  target_year <- config$target_year
  scenario    <- config$scenario
  # B3 (2026-06-06): the scoring/final-map stage shares the orchestrator's
  # `prefix_base` (config$options$prefix_base, default "patch_certified").
  # Read the knob so overriding the option keeps advertised paths consistent.
  # Default unchanged -> byte-identical.
  prefix_base <- config$options$prefix_base %||% "patch_certified"
  prefix <- sprintf("%d_%s_%s", target_year, scenario, prefix_base)
  year_tag <- as.character(target_year)

  dir.create(out_score_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_map_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # 2) Resolve model / recipe to RDS paths the engine reads. Accept objects
  #    (write to temp RDS so the engine's readRDS path is identical) or paths
  #    (forward directly).
  # ---------------------------------------------------------------------------
  resolve_rds <- function(x, what) {
    if (is.character(x) && length(x) == 1L) {
      if (!file.exists(x)) {
        stop(sprintf("'%s' path does not exist: %s", what, x), call. = FALSE)
      }
      return(x)
    }
    p <- tempfile(fileext = ".rds")
    saveRDS(x, p)
    p
  }
  model_rds  <- resolve_rds(model,  "model")
  recipe_rds <- resolve_rds(recipe, "recipe")

  # ---------------------------------------------------------------------------
  # 3) Resolve scoring universe + labelled features + OOF summary to GPKG paths.
  #    The orchestrator passed gpkg_features (features_geometry.gpkg) for BOTH
  #    `unlabeled_gpkg` (scoring) and `labelled_features_gpkg` (train_features
  #    layer). When objects are supplied they are materialised to temp GPKGs.
  # ---------------------------------------------------------------------------
  resolve_gpkg <- function(x, layer, what) {
    if (is.character(x) && length(x) == 1L) {
      if (!file.exists(x)) {
        stop(sprintf("'%s' path does not exist: %s", what, x), call. = FALSE)
      }
      return(x)
    }
    xs <- if (inherits(x, "SpatVector")) sf::st_as_sf(x) else x
    p <- tempfile(fileext = ".gpkg")
    sf::st_write(xs, p, layer = layer, quiet = TRUE, delete_dsn = TRUE)
    p
  }
  unlabeled_gpkg <- resolve_gpkg(scoring_features, scoring_layer,
                                 "scoring_features")

  if (is.null(labelled_features)) {
    labelled_features_gpkg <- config$output_routes$features_geometry_gpkg
  } else {
    labelled_features_gpkg <- resolve_gpkg(labelled_features,
                                           labelled_features_layer,
                                           "labelled_features")
  }
  if (is.null(labelled_features_gpkg) || !nzchar(labelled_features_gpkg)) {
    stop("Could not resolve 'labelled_features' GPKG from config.",
         call. = FALSE)
  }

  if (is.null(oof_summary)) {
    oof_summary <- file.path(
      config$output_routes$oof_dir,
      sprintf("%d_%s_patch_labeled_oof_summary.gpkg", target_year, scenario)
    )
  }
  if (!is.character(oof_summary) || length(oof_summary) != 1L) {
    stop("'oof_summary' must be NULL or a single GPKG path.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 3b) Gate 1D.2 — MANDATORY recipe-driven feature-schema enforcement (check 8)
  #     on the MAIN path, BEFORE the scoring matrix is built.
  #
  # The recipe saved during the FINAL refit is the CANONICAL schema source
  # (recipe$cols$feature_cols, incl. any _isNA companions; medians; factor
  # handling; the 1C.3 5-case reconciliation policy). It is REQUIRED above (a
  # missing recipe already errors), so there is NO silent schema reconstruction
  # from the scoring year. Here we make Gate 1C.4's feature-schema check (check
  # 8) VERIFIABLE and ENFORCED on the orchestrator main path by running the
  # SINGLE Gate 1D.1 validation engine with the model + recipe + the scoring
  # frame's column names — so a schema incompatibility fails fast with a clear
  # aggregated error, instead of being SKIPPED/NOT_VERIFIABLE (the Gate 1C.4
  # flag #1). The recipe stays the schema source: the candidate is the scoring
  # inputs, the expected set is read from the recipe, never re-derived.
  #
  # The reconciliation that follows in the engine (.of_reconcile_scoring_schema)
  # then HEALS the 5 recoverable cases (order / missing / extra / type / level /
  # all-NA) per the explicit 1C.3 policy and RECORDS every decision; check 8
  # only blocks the genuinely-incompatible case (a recipe feature that is
  # neither present nor creatable), which the reconciler cannot recover.
  # Resolve the loaded recipe object (recipe may be an object or an RDS path).
  recipe_obj <- if (is.character(recipe) && length(recipe) == 1L) {
    tryCatch(readRDS(recipe), error = function(e)
      stop(sprintf("Could not read 'recipe' RDS for schema enforcement: %s (%s)",
                   recipe, conditionMessage(e)), call. = FALSE))
  } else {
    recipe
  }
  model_obj <- if (is.character(model) && length(model) == 1L) {
    tryCatch(readRDS(model), error = function(e) NULL)
  } else {
    model
  }
  expected_schema <- .of_model_expected_features(model = model_obj,
                                                 recipe = recipe_obj)
  if (!length(expected_schema)) {
    stop("score_supervised_burned_map(): the FINAL-refit recipe does not carry ",
         "a feature schema (recipe$cols$feature_cols). The recipe is the ",
         "MANDATORY canonical schema source at scoring (Gate 1D.2); refusing to ",
         "reconstruct the schema from the scoring year. Re-run the FINAL refit ",
         "so a complete recipe is produced.", call. = FALSE)
  }
  # Read the scoring frame's column names (the candidate schema) WITHOUT loading
  # geometry/cells: read 0 rows of the resolved scoring layer.
  scoring_feature_names <- tryCatch({
    hd <- sf::st_read(unlabeled_gpkg, layer = scoring_layer, quiet = TRUE)
    nm <- names(sf::st_drop_geometry(hd))
    nm
  }, error = function(e)
    stop(sprintf("score_supervised_burned_map(): could not read scoring layer ",
                 "'%s' column names for schema enforcement (%s).",
                 scoring_layer, conditionMessage(e)), call. = FALSE))
  # Run the SINGLE validation engine in strict mode (check 8 now verifiable).
  # data_base/composite_base/result_name come from cfg$options so the engine
  # resolves the SAME input paths the run consumes (other checks stay PASS /
  # NOT_VERIFIABLE; only the feature-schema check is newly enforced here).
  .opts <- config$options %||% list()
  validate_supervised_execution(
    config                = config,
    strict                = TRUE,
    target_year           = target_year,
    model                 = model_obj,
    recipe                = recipe_obj,
    scoring_feature_names = scoring_feature_names,
    data_base             = .opts$data_base,
    composite_base        = .opts$composite_base,
    result_name           = .opts$result_name %||% "Min_Min")
  if (isTRUE(verbose)) {
    message(sprintf(paste0("[Gate 1D.2] feature-schema check PASS: recipe ",
                           "schema (%d features) is covered by the scoring ",
                           "inputs; reconciliation will heal recoverable cases."),
                    length(expected_schema)))
  }

  # ---------------------------------------------------------------------------
  # 4) Score + export final map (faithful MOVE of the fused wrapper's
  #    score_burnedlike_and_export_final_map call). Every argument/value matches
  #    the historical orchestrator call: result_dir, prefix, the qa/labelled/
  #    unlabeled inputs, model/recipe, the 08_SCORED/09_FINAL_MAP out dirs, the
  #    export flag, and the three CURRENTYEAR_* temporal thresholds.
  # ---------------------------------------------------------------------------
  out <- score_burnedlike_and_export_final_map(
    result_dir = result_dir,
    prefix     = prefix,
    year_tag   = year_tag,

    qa_labelled_gpkg  = oof_summary,
    qa_labelled_layer = qa_labelled_layer,

    labelled_features_gpkg  = labelled_features_gpkg,
    labelled_features_layer = labelled_features_layer,

    model_rds  = model_rds,
    recipe_rds = recipe_rds,

    unlabeled_gpkg  = unlabeled_gpkg,
    unlabeled_layer = scoring_layer,

    out_score_dir = out_score_dir,
    out_map_dir   = out_map_dir,

    export_burned_like = export_burned_like,
    preyear_overlap_threshold = preyear_overlap_threshold,
    hotspot_density_threshold = hotspot_density_threshold,
    temporal_penalty_floor = temporal_penalty_floor,
    overwrite = overwrite,
    verbose   = verbose
  )

  # ---------------------------------------------------------------------------
  # 5) Assemble the documented return (objects + paths). The path fields mirror
  #    what the engine wrote (its `files` list).
  # ---------------------------------------------------------------------------
  files <- out$files %||% list()

  list(
    deterministic_scored      = out$deterministic_scored,
    final_map_full            = out$final_map_full,
    final_map                 = out$final_map,
    burned_like_scored        = out$burned_like_scored,
    deterministic_scored_gpkg = files$scored_gpkg %||%
      file.path(out_score_dir, paste0(prefix, "_deterministic_scored.gpkg")),
    final_map_gpkg            = files$final_map_gpkg %||%
      file.path(out_map_dir, paste0(prefix, "_final_map.gpkg")),
    final_map_counts_csv      = files$counts_csv %||%
      file.path(out_map_dir, paste0(prefix, "_final_map_counts.csv")),
    burned_like_gpkg          = files$burned_like_gpkg,
    burned_like_counts_csv    = files$burned_like_counts_csv,
    # Gate 1E (2026-06-09): the scoring-matrix structural fingerprint (asserted
    # == the SAVED FINAL fingerprint before predicting). The orchestrator threads
    # this into the run manifest's scoring leg.
    scoring_schema_fingerprint = out$scoring_schema_fingerprint,
    scored                    = out
  )
}
