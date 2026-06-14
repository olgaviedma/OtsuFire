#' Score the deterministic polygon universe and export supervised final maps
#'
#' @description
#' Supervised-pipeline stage E (scoring + final-map export). This is the
#' SCORING half of the formerly fused
#' `run_train_final_model_and_export_final_map()` wrapper. It wraps the internal
#' engine `score_burnedlike_and_export_final_map()` with the exact arguments and
#' values the fused wrapper passed (the trained model + recipe, the labelled
#' features + OOF labelled-summary join inputs, the current-year temporal
#' adjustment thresholds from config, and the burned-like export flag), and
#' writes the canonical `08_SCORED` + `09_FINAL_MAP` outputs under the
#' `"<year>_<scenario>_patch_certified"` prefix.
#'
#' This is the standalone, exported implementation of the scoring stage that the
#' one-year orchestrator [run_oneyear_supervised_pipeline()] delegates to, so
#' calling it directly produces byte-identical scoring / final-map outputs to a
#' full run.
#'
#' @section Mandatory recipe + enforced feature-schema check (Gate 1D.2):
#' The `recipe` produced by the FINAL refit is a \strong{mandatory} input and is
#' the \strong{canonical} source of the scoring feature schema on the MAIN path:
#' feature names and order, types, `_isNA` companions, imputation medians,
#' categorical levels, and the missing/extra/unknown-level policies all come
#' from `recipe$cols$feature_cols` / `recipe$cols$x_cols` /
#' `recipe$impute$numeric_medians` / `recipe$impute$impute_factor_missing` —
#' the schema is \strong{never} reconstructed from the scoring-year data. Before
#' the scoring matrix is built, this function runs the SINGLE Gate 1D.1
#' validation engine ([validate_supervised_execution()], `strict = TRUE`) with
#' the model + recipe + the scoring frame's column names, so the
#' feature-schema check (check 8) is \strong{ENFORCED and verifiable} on the
#' main path (a genuinely incompatible schema fails fast; it is no longer
#' SKIPPED / NOT_VERIFIABLE). The check is recipe-driven and aligned with the
#' Gate 1C.3 5-case reconciliation policy: column-order altered → realigned to
#' recipe order; a feature missing → created as `NA_real_` (xgboost-missing) +
#' `_isNA`; an extra column → dropped (cannot shift the matrix); a type change →
#' coerced per the recipe or set NA; a new categorical level → recipe sentinel
#' (never invented); a feature all-NA → recipe-median impute / xgboost-missing.
#' Every reconciliation decision is RECORDED (no silent corrections). Proving:
#' `recipe of refit model -> final scoring -> exact same schema`.
#'
#' @details
#' The current-year temporal adjustment is applied inside the engine using the
#' three thresholds `preyear_overlap_threshold`, `hotspot_density_threshold` and
#' `temporal_penalty_floor`. Their defaults are the orchestrator's
#' `CURRENTYEAR_*` constants (`0.70`, `0.001`, `0.10`) so the adjustment is
#' byte-identical; they are exposed as parameters so standalone callers can
#' override them.
#'
#' `model` and `recipe` are accepted as in-memory objects (the
#' `train_final_burned_model()` return fields) OR as RDS paths. When objects are
#' supplied they are written to temporary RDS files so the engine's
#' `readRDS()`-based loading path is byte-identical; when paths are supplied they
#' are forwarded directly.
#'
#' @param scoring_features sf / data.frame OR a single GPKG path. The
#'   deterministic scoring universe (the `scoring_features` layer of
#'   `03_FEATURES/features_geometry.gpkg`). Forwarded as the engine's
#'   `unlabeled_gpkg` (read back by `scoring_layer`).
#' @param model Fitted xgboost model object OR an RDS path.
#' @param recipe Training recipe list OR an RDS path. \strong{Mandatory} — this
#'   is the canonical feature-schema source at scoring (Gate 1D.2). It must carry
#'   `recipe$cols$feature_cols` (the refit feature schema); the function errors
#'   rather than reconstruct the schema from the scoring year.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Used to derive `result_dir`
#'   (`output_routes$base`), `out_score_dir` (`08_SCORED`), `out_map_dir`
#'   (`09_FINAL_MAP`), the output prefix
#'   (`"<year>_<scenario>_patch_certified"`), and `target_year` / `scenario`.
#' @param oof_summary Character path to the OOF labelled-summary GPKG
#'   (`_labeled_oof_summary.gpkg`) OR `NULL`. Forwarded as the engine's
#'   `qa_labelled_gpkg`. Required for the OOF `p_oof_mean` join; when `NULL` the
#'   function falls back to the canonical `05_OOF` path
#'   (`"<year>_<scenario>_patch_labeled_oof_summary.gpkg"`).
#' @param labelled_features sf / data.frame OR a single GPKG path. The labelled
#'   training features (the `train_features` layer of
#'   `03_FEATURES/features_geometry.gpkg`) used to join OOF predictions by
#'   `source_poly_id`. When `NULL` the function falls back to
#'   `config$output_routes$features_geometry_gpkg`.
#' @param export_burned_like Logical. Whether to export the burned-like subset
#'   (`_burned_like_scored.gpkg` + `_counts.csv`). Default `TRUE`.
#' @param preyear_overlap_threshold Numeric. Current-year temporal-adjustment
#'   pre-year overlap threshold. Default `0.70` (orchestrator
#'   `CURRENTYEAR_PREYEAR_OVERLAP_THR`).
#' @param hotspot_density_threshold Numeric. Current-year hotspot-density
#'   threshold. Default `0.001` (orchestrator `CURRENTYEAR_HOTSPOT_DENSITY_THR`).
#' @param temporal_penalty_floor Numeric. Current-year temporal penalty floor.
#'   Default `0.10` (orchestrator `CURRENTYEAR_TEMPORAL_PENALTY_FLOOR`).
#' @param qa_labelled_layer Character. Layer in `oof_summary`. Default
#'   `"labeled_oof_summary"`.
#' @param labelled_features_layer Character. Layer in `labelled_features`.
#'   Default `"train_features"`.
#' @param scoring_layer Character. Layer in `scoring_features`. Default
#'   `"scoring_features"`.
#' @param out_map_dir Character or `NULL`. Output folder for the `09_FINAL_MAP`
#'   outputs. Defaults to `config$output_routes$final_map_dir`.
#' @param out_score_dir Character or `NULL`. Output folder for the `08_SCORED`
#'   outputs. Defaults to `config$output_routes$scored_dir`.
#' @param overwrite Logical. Forwarded to the engine. Default `TRUE`.
#' @param verbose Logical. Forwarded to the engine. Default `TRUE`.
#'
#' @section Output-layer contract (Gate 2):
#' The written `<prefix>_final_map.gpkg` carries three layers with DIFFERENT,
#' deliberate contracts, so they are NOT expected to have equal row counts:
#' \itemize{
#'   \item `deterministic_scored` and `final_map_full` (identical content) — the
#'     COMPLETE scored deterministic universe: every candidate the FINAL model
#'     scored, one row per input polygon, all IDs/geometries preserved (the
#'     authoritative "no candidate is lost" layers; `nrow ==` the scoring
#'     universe size, asserted).
#'   \item `final_map` — the CURRENT-YEAR PUBLIC burned map: `final_map_full`
#'     with the public column subset AND the current-year temporal filter
#'     applied. It DROPS exactly the rows flagged `current_year_public_drop`
#'     (a current-year temporal conflict — pre-year overlap
#'     `>= preyear_overlap_threshold` or `preyear_action == "drop"` — with WEAK
#'     current-year hotspot support: a polygon almost certainly re-detecting the
#'     PRE-YEAR fire with no current-year evidence). Such polygons are
#'     legitimately excluded from the current-year public map while remaining in
#'     the full scored layers, so
#'     `nrow(final_map) == nrow(final_map_full) -`
#'     `sum(final_map_full$current_year_public_drop %in% TRUE)`. This is a
#'     documented methodological exclusion, NOT a silent row loss.
#' }
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `deterministic_scored` — the full scored deterministic universe
#'       (sf, == `final_map_full`). ALL scored candidates (no row dropped).
#'     \item `final_map_full` — same as `deterministic_scored`.
#'     \item `final_map` — the current-year public final map (sf); the
#'       temporally-filtered, public-column subset (see the output-layer
#'       contract above): `nrow <= nrow(final_map_full)`.
#'     \item `burned_like_scored` — the burned-like subset (sf), or `NULL` when
#'       `export_burned_like = FALSE`.
#'     \item `deterministic_scored_gpkg`, `final_map_gpkg`,
#'       `final_map_counts_csv`, `burned_like_gpkg`,
#'       `burned_like_counts_csv` — written paths.
#'   }
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
#'   train_features = "03_FEATURES/features_geometry.gpkg", config = cfg,
#'   oof_agg = "05_OOF/2017_balanced_patch_oof_agg.csv"
#' )
#' sm <- score_supervised_burned_map(
#'   scoring_features  = "03_FEATURES/features_geometry.gpkg",
#'   model = tm$model, recipe = tm$recipe, config = cfg,
#'   oof_summary = "05_OOF/2017_balanced_patch_labeled_oof_summary.gpkg"
#' )
#' sm$final_map_gpkg
#' }
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
  #    "<year>_<scenario>_patch_certified").
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
