# =============================================================================
# CANONICAL versioned supervised usage script  (01 / 4)
# OtsuFire supervised burned-area mapping — ONE-YEAR pipeline
# =============================================================================
# *** THIS IS THE CANONICAL, VERSIONED COPY (lives in inst/scripts/ of the
#     OtsuFire package). Do NOT treat the external 00_USAGE/02_SUPERVISED_USAGE
#     copies as the source of truth — those are WORKING COPIES. See inst/scripts/
#     README.md and MANIFEST.csv for the sync contract. ***
#
# Targets OtsuFire >= 0.5.0.
#
# WHAT THIS SCRIPT DEMONSTRATES
#   The corrected common-recipe supervised pipeline for a single year, using the
#   PUBLIC API only and the cfg-as-single-source-of-truth pattern. The pipeline
#   is the chain of 6 exported modular functions:
#       build_supervised_burned_config()  -> cfg   (single source of truth)
#       build_supervised_training_pools() -> 01_POOLS
#       make_spatial_folds()              -> 02_FOLDS
#       extract_supervised_features()     -> 03_FEATURES
#       run_oof_diagnostics()             -> 04_MATRIX + 05_OOF  (SANITY signal)
#       train_final_burned_model()        -> 07_FINAL_MODEL_V2
#       score_supervised_burned_map()     -> 08_SCORED + 09_FINAL_MAP
#   OtsuFire uses a SINGLE training procedure: the number of boosting rounds is
#   selected with an inner validation split and early stopping, then the recipe
#   and model are refitted on all available training data before prediction.
#   There is no training-protocol choice.
#
# METHODOLOGICAL PARAMS LIVE IN THE CFG (NOT in function-level args)
#   Every methodological / training-control knob (the 4 negative-bucket caps,
#   oof_sampling, seeds, xgb params, feature whitelist) is set
#   ONCE in build_supervised_burned_config(). Passing those same knobs as
#   function-level arguments to the stage functions is a DEPRECATED shim that
#   emits an `otsufire_deprecated_param` warning — this script does NOT do that.
#
# p_burned IS A MODEL SCORE, NOT A CALIBRATED PROBABILITY
#   The exported `p_burned` / `p_burned_model` column is the raw XGBoost score in
#   [0, 1]. It is NOT a calibrated probability and must not be read as one (per
#   the package docs). Use it for ranking / thresholding, and pick the operating
#   threshold from external (EFFIS) validation, not from the score scale.
#
# NEGATIVE POOL (4 buckets, all_sources policy — the only supported policy):
#   1. contextual exclusions  (deterministic DROP, non-spectral)   cap 0.25 x n_burned
#   2. spectral hard negatives (deterministic DROP, spectral)      cap 1.0  x n_burned (pkg default; see 02_phaseB for 2.0)
#   3. random burnable bg      (burnable cells outside buffers)    cap 1.0  x n_burned
#   4. otsu unburned           (legacy Otsu residual patches)      cap 1.0  x n_burned
#
# HOW TO RUN
#   This script does NOT execute heavy work when sourced (for a parse/syntax
#   check). The real pipeline body is guarded behind `RUN <- FALSE`. Edit the
#   CONFIG block below (replace every <PATH_TO_...> placeholder), set RUN <- TRUE,
#   then source the file.
# =============================================================================


# -------------------------------------------------------------------
# CONFIG BLOCK  --  EDIT THESE  (replace every <PATH_TO_...> placeholder)
# -------------------------------------------------------------------
RUN <- FALSE   # <<< set TRUE to actually run the heavy pipeline

paths <- list(
  # Development package tree (loaded via pkgload::load_all, NOT library()).
  pkg_root           = "<PATH_TO_OtsuFire_v02_rebuild>",          # TODO

  # MAIN supervised input: the deterministic decision layer for the SAME
  # year + scenario (a GPKG with layer internal_decisions).
  internal_decisions = "<PATH_TO_internal_decisions.gpkg>",       # TODO

  # Main annual change-index composite (band1 = summer RBR, band2 = DOY_post).
  change_index       = "<PATH_TO_MinMin_YEAR_mosaic_res90m.tif>", # TODO

  # Target-year hotspot points (used post-2000). Set to NULL for pre-MODIS years
  # (see 06_no_hotspots_reduced_historical.R).
  hotspots           = "<PATH_TO_hotspots_iberia_YEAR.geojson>",  # TODO or NULL

  # Optional external burned reference (only consumed by validate_fire_maps()).
  reference_burned   = "<PATH_TO_Effis_CA_YEAR_maskKeep_summer.shp>", # TODO or NULL

  # Roots forwarded into cfg$options for the legacy all_sources Otsu builder.
  data_base          = "<PATH_TO_1_DATA>",                        # TODO
  composite_base     = "<PATH_TO_1_DATA/Imagery/Composites_90m>", # TODO
  output_dir         = "<PATH_TO_1_DATA/Results>"                 # root; package appends <year>/<run>/SUPERVISED/<scenario>
)

# External GDAL/Python tool paths (the all_sources legacy Otsu builder needs them).
tool_paths <- list(
  python_exe             = "<PATH_TO_python.exe>",                # TODO
  gdal_polygonize_script = "<PATH_TO_gdal_polygonize.py>",        # TODO
  gdalwarp_path          = "<PATH_TO_gdalwarp.exe>",              # TODO
  ogr2ogr_exe            = "<PATH_TO_ogr2ogr.exe>"                # TODO
)

# Run-level settings.
target_year  <- 2017L
scenario     <- "balanced"     # one of: balanced | original | lax | restrictive
run_name     <- "Min_Min"      # composite family / run identifier
use_hotspots <- target_year > 2000   # hotspots are used post-2000

# 4-bucket negative-pool caps (package defaults).
caps <- list(contextual = 0.25, spectral = 1.0, random = 1.0, otsu = 1.0)


# -------------------------------------------------------------------
# main() — the heavy pipeline. NOT called unless RUN is TRUE.
# -------------------------------------------------------------------
main <- function() {
  suppressPackageStartupMessages({
    library(pkgload); library(sf); library(terra)
  })

  # Load the development tree (NOT an installed copy that may lack the modular API).
  if ("OtsuFire" %in% loadedNamespaces()) {
    try(unloadNamespace("OtsuFire"), silent = TRUE)
  }
  suppressMessages(pkgload::load_all(paths$pkg_root, quiet = TRUE))
  ver <- utils::packageVersion("OtsuFire")
  cat("OtsuFire version:", as.character(ver), "\n")
  stopifnot(ver >= "0.5.0")

  # ---- STEP 1. Build the cfg (single source of truth). Nothing computed yet. --
  cfg <- build_supervised_burned_config(
    scenario             = scenario,
    internal_decisions   = paths$internal_decisions,
    change_index         = paths$change_index,
    hotspots             = paths$hotspots,
    reference_burned_map = paths$reference_burned,
    target_year          = target_year,
    output_dir           = paths$output_dir,
    run_name             = run_name,
    # ---- CANONICAL methodological params (set HERE, read by every stage) ----
    cap_contextual       = caps$contextual,
    cap_spectral         = caps$spectral,
    cap_random           = caps$random,
    cap_otsu             = caps$otsu,
    options = list(
      data_base                       = paths$data_base,
      composite_base                  = paths$composite_base,
      result_name                     = run_name,
      legacy_otsu_mode                = "burnable_only",
      legacy_otsu_threshold           = 0,
      legacy_reference_otsu_threshold = 100,
      legacy_sample_n                 = 2000L,
      legacy_reuse_existing           = TRUE,
      legacy_write_output             = TRUE,
      unb_verbose                     = TRUE
    )
  )
  # Attach external tool paths (consumed by the pools stage).
  cfg$tool_paths$python_exe             <- tool_paths$python_exe
  cfg$tool_paths$gdal_polygonize_script <- tool_paths$gdal_polygonize_script
  cfg$tool_paths$gdalwarp_path          <- tool_paths$gdalwarp_path
  cfg$tool_paths$ogr2ogr_exe            <- tool_paths$ogr2ogr_exe
  stopifnot(all(vapply(cfg$tool_paths, file.exists, logical(1))))

  # ---- PRE-RUN CHECK: validate_supervised_execution(cfg, strict = TRUE) ------
  # Fails fast on blocking inputs/schema problems BEFORE any heavy compute.
  validate_supervised_execution(cfg, strict = TRUE)

  # ---- STEP 2. Training pools (4-bucket all_sources negative pool). ----------
  pools <- build_supervised_training_pools(
    config = cfg, write_outputs = TRUE, overwrite = TRUE
  )

  # ---- STEP 3. Spatial CV folds (kept together by `fire`). -------------------
  folds <- make_spatial_folds(
    train_labelled = pools$train_labeled, config = cfg,
    split_unit = "fire", out_dir = NULL, write_outputs = TRUE, overwrite = TRUE
  )

  # ---- STEP 4. Features (raster support stack + optional hotspots). ----------
  feats <- extract_supervised_features(
    train_with_folds = folds$train_with_folds,
    scoring_pool     = pools$scoring_pool,
    config           = cfg, use_hotspots = use_hotspots,
    out_dir = NULL, write_outputs = TRUE
  )

  # ---- STEP 5. OOF diagnostics (SANITY signal, not the final metrics). -------
  oof <- run_oof_diagnostics(
    train_features   = feats$train_features,
    scoring_features = feats$scoring_features,
    config           = cfg,
    fold_cols        = c("fold_rep1", "fold_rep2"),
    out_dir          = NULL,
    labelled_gpkg    = folds$train_with_folds_gpkg,  # PATH (writes labeled_oof_summary)
    labelled_layer   = "train_with_folds"
  )

  # ---- STEP 6. Train the deployed FINAL model. -------------------------------
  model <- train_final_burned_model(
    train_features = feats$train_features,
    config         = cfg,
    oof_agg        = oof$oof_agg_csv,   # PATH (qa enrichment only)
    out_dir        = NULL, overwrite = TRUE, verbose = TRUE
  )

  # ---- STEP 7. Score + export the final map. p_burned is a MODEL SCORE. ------
  scored <- score_supervised_burned_map(
    scoring_features   = feats$scoring_features,
    model              = model$model,
    recipe             = model$recipe,
    config             = cfg,
    oof_summary        = oof$labeled_oof_summary_gpkg,  # PATH
    export_burned_like = TRUE, out_map_dir = NULL, verbose = TRUE
  )

  # ---- STEP 8. ACCEPTANCE SIGNAL (KB3): drops should NOT be rescued. ---------
  fm <- scored$final_map_full
  drops_med <- stats::median(fm$p_burned_model[fm$class_final == "drop"], na.rm = TRUE)
  keeps_med <- stats::median(fm$p_burned_model[fm$class_final == "keep"], na.rm = TRUE)
  cat(sprintf("[KB3] drops median p_burned_model=%.4f ; keeps median=%.4f\n",
              drops_med, keeps_med))

  invisible(list(cfg = cfg, pools = pools, folds = folds, feats = feats,
                 oof = oof, model = model, scored = scored))
}

if (isTRUE(RUN)) {
  main()
} else {
  message("01_supervised_oneyear_legacy.R sourced with RUN=FALSE (no compute). ",
          "Edit the CONFIG block, set RUN <- TRUE, then source again.")
}
# =============================================================================
# End of canonical one-year supervised script.
# (Filename retains the historical "legacy" suffix for the in-repo sync
#  contract; it no longer refers to any training-protocol choice — OtsuFire
#  always uses inner-early-stopping selection + full-data refit.)
# =============================================================================
