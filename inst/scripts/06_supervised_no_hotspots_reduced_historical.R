# =============================================================================
# CANONICAL versioned supervised usage script  (04 / 4)
# OtsuFire — no-hotspots / reduced-historical (pre-MODIS) example
# =============================================================================
# *** CANONICAL, VERSIONED COPY (inst/scripts/). External 00_USAGE copies are
#     WORKING COPIES. *** Targets OtsuFire >= 0.5.0.
#
# WHAT THIS SCRIPT DEMONSTRATES — the package runs WITHOUT hotspots:
#   * hotspots = NULL  and  use_hotspots = FALSE  (a pre-MODIS-style year, e.g.
#     before ~2000). The hs_* hotspot feature block is simply absent / NA, and
#     the engine treats those features as missing (not as a sentinel value).
#   * This is the configuration to use for reduced-historical years where no
#     target-year fire-detection points exist. It also exercises the all-NA
#     feature / no-hotspot scoring path the package was hardened against.
#
# p_burned IS A MODEL SCORE, NOT A CALIBRATED PROBABILITY (per the package docs).
# Heavy work is guarded behind RUN <- FALSE (this file parses without computing).
# =============================================================================


# -------------------------------------------------------------------
# CONFIG BLOCK  --  EDIT THESE  (replace every <PATH_TO_...> placeholder)
# -------------------------------------------------------------------
RUN <- FALSE

paths <- list(
  pkg_root           = "<PATH_TO_OtsuFire_v02_rebuild>",          # TODO
  internal_decisions = "<PATH_TO_internal_decisions.gpkg>",       # TODO (pre-MODIS year)
  change_index       = "<PATH_TO_MinMin_YEAR_mosaic_res90m.tif>", # TODO
  data_base          = "<PATH_TO_1_DATA>",                        # TODO
  composite_base     = "<PATH_TO_1_DATA/Imagery/Composites_90m>", # TODO
  output_dir         = "<PATH_TO_1_DATA/Results>"                 # root
)
tool_paths <- list(
  python_exe             = "<PATH_TO_python.exe>",                # TODO
  gdal_polygonize_script = "<PATH_TO_gdal_polygonize.py>",        # TODO
  gdalwarp_path          = "<PATH_TO_gdalwarp.exe>",              # TODO
  ogr2ogr_exe            = "<PATH_TO_ogr2ogr.exe>"                # TODO
)

# A pre-MODIS / reduced-historical year. NO hotspots exist for it.
target_year  <- 1994L
scenario     <- "balanced"
run_name     <- "Min_Min"

# THE KEY OF THIS EXAMPLE: no hotspots at all.
hotspots     <- NULL    # <<< no target-year fire-detection points
use_hotspots <- FALSE   # <<< the hs_* feature block is disabled

caps <- list(contextual = 0.25, spectral = 1.0, random = 1.0, otsu = 1.0)


main <- function() {
  suppressPackageStartupMessages({ library(pkgload); library(sf); library(terra) })
  if ("OtsuFire" %in% loadedNamespaces()) try(unloadNamespace("OtsuFire"), silent = TRUE)
  suppressMessages(pkgload::load_all(paths$pkg_root, quiet = TRUE))
  stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

  cfg <- build_supervised_burned_config(
    scenario           = scenario,
    internal_decisions = paths$internal_decisions,
    change_index       = paths$change_index,
    hotspots           = hotspots,            # NULL
    target_year        = target_year,
    output_dir         = paths$output_dir,
    run_name           = run_name,
    cap_contextual     = caps$contextual,
    cap_spectral       = caps$spectral,
    cap_random         = caps$random,
    cap_otsu           = caps$otsu,
    options = list(
      data_base = paths$data_base, composite_base = paths$composite_base,
      result_name = run_name, otsu_negative_mode = "burnable_only",
      otsu_negative_threshold = 0, otsu_negative_reference_threshold = 100,
      otsu_negative_reuse_existing = TRUE,
      otsu_negative_write_output = TRUE, unb_verbose = TRUE
    )
  )
  cfg$tool_paths$python_exe             <- tool_paths$python_exe
  cfg$tool_paths$gdal_polygonize_script <- tool_paths$gdal_polygonize_script
  cfg$tool_paths$gdalwarp_path          <- tool_paths$gdalwarp_path
  cfg$tool_paths$ogr2ogr_exe            <- tool_paths$ogr2ogr_exe

  validate_supervised_execution(cfg, strict = TRUE)

  pools <- build_supervised_training_pools(config = cfg, write_outputs = TRUE, overwrite = TRUE)
  folds <- make_spatial_folds(train_labelled = pools$train_labeled, config = cfg,
                              split_unit = "fire", write_outputs = TRUE, overwrite = TRUE)
  # use_hotspots = FALSE -> the hs_* feature block is not built.
  feats <- extract_supervised_features(train_with_folds = folds$train_with_folds,
                                       scoring_pool = pools$scoring_pool, config = cfg,
                                       use_hotspots = use_hotspots, write_outputs = TRUE)
  oof   <- run_oof_diagnostics(train_features = feats$train_features,
                               scoring_features = feats$scoring_features, config = cfg,
                               fold_cols = c("fold_rep1", "fold_rep2"),
                               labelled_gpkg = folds$train_with_folds_gpkg,
                               labelled_layer = "train_with_folds")
  model <- train_final_burned_model(train_features = feats$train_features, config = cfg,
                                    oof_agg = oof$oof_agg_csv, overwrite = TRUE, verbose = TRUE)
  scored <- score_supervised_burned_map(scoring_features = feats$scoring_features,
                                        model = model$model, recipe = model$recipe, config = cfg,
                                        oof_summary = oof$labeled_oof_summary_gpkg,
                                        export_burned_like = TRUE, verbose = TRUE)

  hs_cols <- grep("^hs_", names(feats$train_features), value = TRUE)
  cat(sprintf("[no-hotspots] hs_* feature columns present: %d (expected 0 or all-NA)\n",
              length(hs_cols)))
  invisible(list(cfg = cfg, scored = scored))
}

if (isTRUE(RUN)) main() else
  message("06_supervised_no_hotspots_reduced_historical.R sourced with RUN=FALSE (no compute).")
# =============================================================================
# End of no-hotspots / reduced-historical example.
# =============================================================================
