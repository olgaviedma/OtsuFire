# =============================================================================
# CANONICAL versioned supervised usage script  (03 / 4)
# OtsuFire — EFFIS external validation of a supervised burned map
# =============================================================================
# *** CANONICAL, VERSIONED COPY (inst/scripts/). External 00_USAGE copies are
#     WORKING COPIES. *** Targets OtsuFire >= 0.5.0.
#
# WHAT THIS SCRIPT DEMONSTRATES — TWO DISTINCT validators, do not confuse them:
#   * validate_supervised_execution(cfg, strict = TRUE)
#       PRE-RUN structural check of the cfg/inputs/schema (no external truth).
#       Used as the fail-fast gate before heavy compute (shown in scripts 01-02).
#   * validate_fire_maps(input_shapefile, ref_shapefile, ...)
#       POST-RUN ACCURACY validation of the thresholded supervised map against
#       an EFFIS reference burned layer (pixel/area metrics). This is what THIS
#       script focuses on.
#
# RECIPE:
#   1. Build cfg + run the pipeline (or load an already-scored final map).
#   2. Threshold the final map on p_burned_model (a MODEL SCORE, not a calibrated
#      probability — choose the operating threshold from validation, not the
#      score scale).
#   3. Write the thresholded polygons + a geometry-only validation mask.
#   4. Call validate_fire_maps() against the EFFIS reference.
#
# Heavy work is guarded behind RUN <- FALSE (this file parses without computing).
# =============================================================================


# -------------------------------------------------------------------
# CONFIG BLOCK  --  EDIT THESE  (replace every <PATH_TO_...> placeholder)
# -------------------------------------------------------------------
RUN <- FALSE

paths <- list(
  pkg_root           = "<PATH_TO_OtsuFire_v02_rebuild>",          # TODO
  internal_decisions = "<PATH_TO_internal_decisions.gpkg>",       # TODO
  change_index       = "<PATH_TO_MinMin_YEAR_mosaic_res90m.tif>", # TODO
  hotspots           = "<PATH_TO_hotspots_iberia_YEAR.geojson>",  # TODO or NULL
  data_base          = "<PATH_TO_1_DATA>",                        # TODO
  composite_base     = "<PATH_TO_1_DATA/Imagery/Composites_90m>", # TODO
  output_dir         = "<PATH_TO_1_DATA/Results>",                # root

  # ---- EFFIS validation inputs ----
  reference_burned   = "<PATH_TO_Effis_CA_YEAR_maskKeep_summer.shp>",   # EFFIS reference
  burnable_raster    = "<PATH_TO_burneable_mask_binary_corineXX_wgs84.tif>", # burnable mask (WGS84 ok; reprojected internally)
  validation_mask    = "<PATH_TO_mask_Peninsula_3035.shp>",            # study-area mask
  strata_raster      = "<PATH_TO_STRATA/strata_CLC_YYYY_res30.tif>",   # optional; NULL to skip strata
  strata_lut         = "<PATH_TO_LUT/lut_full_strata8_v1.csv>"         # optional; NULL to skip strata
)
tool_paths <- list(
  python_exe             = "<PATH_TO_python.exe>",                # TODO
  gdal_polygonize_script = "<PATH_TO_gdal_polygonize.py>",        # TODO
  gdalwarp_path          = "<PATH_TO_gdalwarp.exe>",              # TODO
  ogr2ogr_exe            = "<PATH_TO_ogr2ogr.exe>"                # TODO
)

target_year         <- 2017L
scenario            <- "balanced"
run_name            <- "Min_Min"
use_hotspots        <- target_year > 2000
OPERATING_THRESHOLD <- 0.50   # threshold on the MODEL SCORE p_burned_model
# GATE 6.7 (2026-06-12): caps via the typed PUBLIC negative_pool_params block.
caps <- c(random = 1.0, otsu = 1.0)


main <- function() {
  suppressPackageStartupMessages({ library(pkgload); library(sf); library(terra) })
  if ("OtsuFire" %in% loadedNamespaces()) try(unloadNamespace("OtsuFire"), silent = TRUE)
  suppressMessages(pkgload::load_all(paths$pkg_root, quiet = TRUE))
  stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

  # ---- Build cfg + run the pipeline to obtain a scored final map. ------------
  cfg <- build_supervised_burned_config(
    run_label            = scenario,
    internal_decisions   = paths$internal_decisions,
    change_index         = paths$change_index,
    hotspots             = paths$hotspots,
    reference_burned_map = paths$reference_burned,
    target_year          = target_year,
    output_dir           = paths$output_dir,
    run_name             = run_name,
    # GATE 6.7 (2026-06-12): typed PUBLIC negative_pool_params + runtime_options.
    negative_pool_params = list(
      random = list(n_cells = 1500L, rbr_quantile = 0.50),
      otsu   = list(candidate_threshold = 0, reference_threshold = 100),
      caps   = caps
    ),
    runtime_options = list(reuse_existing = TRUE, write_outputs = TRUE,
                           verbose = TRUE),
    options = list(
      data_base = paths$data_base, composite_base = paths$composite_base,
      result_name = run_name
    )
  )
  cfg$tool_paths$python_exe             <- tool_paths$python_exe
  cfg$tool_paths$gdal_polygonize_script <- tool_paths$gdal_polygonize_script
  cfg$tool_paths$gdalwarp_path          <- tool_paths$gdalwarp_path
  cfg$tool_paths$ogr2ogr_exe            <- tool_paths$ogr2ogr_exe

  validate_supervised_execution(cfg, strict = TRUE)   # PRE-RUN check (distinct validator)

  pools <- build_supervised_training_pools(config = cfg, write_outputs = TRUE, overwrite = TRUE)
  folds <- make_spatial_folds(train_labelled = pools$train_labeled, config = cfg,
                              split_unit = "fire", write_outputs = TRUE, overwrite = TRUE)
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

  # ---- POST-RUN: EFFIS validation via validate_fire_maps(). ------------------
  vdir <- file.path(paths$output_dir, as.character(target_year), run_name,
                    "SUPERVISED", scenario, "VALIDATION_EFFIS")
  dir.create(vdir, recursive = TRUE, showWarnings = FALSE)

  # Threshold the final map on the MODEL SCORE p_burned_model.
  fm <- scored$final_map_full
  p  <- suppressWarnings(as.numeric(sf::st_drop_geometry(fm)$p_burned_model))
  keep <- fm[is.finite(p) & p >= OPERATING_THRESHOLD, , drop = FALSE]
  keep <- sf::st_make_valid(sf::st_zm(keep, drop = TRUE, what = "ZM"))
  keep <- keep[!sf::st_is_empty(keep), , drop = FALSE]
  stopifnot(nrow(keep) > 0L)
  keep_path <- file.path(vdir, "thresholded_burned.gpkg")
  if (file.exists(keep_path)) unlink(keep_path, force = TRUE)
  sf::st_write(keep, keep_path, layer = "thresholded_burned", quiet = TRUE)

  # Geometry-only validation mask.
  mask_sf <- sf::st_make_valid(
    sf::st_zm(sf::st_read(paths$validation_mask, quiet = TRUE), drop = TRUE, what = "ZM"))
  mask_sf <- mask_sf[!sf::st_is_empty(mask_sf), 0, drop = FALSE]
  mask_sf$mask_id <- seq_len(nrow(mask_sf))
  mask_path <- file.path(vdir, "validation_mask_geometry_only.gpkg")
  if (file.exists(mask_path)) unlink(mask_path, force = TRUE)
  sf::st_write(mask_sf, mask_path, layer = "mask", quiet = TRUE)

  has_strata <- !is.null(paths$strata_raster) && !is.null(paths$strata_lut) &&
    file.exists(paths$strata_raster) && file.exists(paths$strata_lut)

  validation <- validate_fire_maps(
    input_shapefile      = keep_path,
    ref_shapefile        = paths$reference_burned,
    mask_shapefile       = mask_path,
    burnable_raster      = paths$burnable_raster,
    year_target          = target_year,
    validation_dir       = vdir,
    metrics_type         = "all",
    dissolve_ref_by      = "id",
    dissolve_input_by    = NULL,
    strata_raster        = if (has_strata) paths$strata_raster else NULL,
    strata_lut           = if (has_strata) paths$strata_lut else NULL,
    observability_raster = paths$change_index,
    ref_end_doy_col      = "end_doy",
    ref_start_doy_col    = "start_doy"
  )
  cat("[EFFIS] validation complete. Metrics:\n")
  print(validation$metrics)
  invisible(list(cfg = cfg, scored = scored, validation = validation))
}

if (isTRUE(RUN)) main() else
  message("05_supervised_effis_validation.R sourced with RUN=FALSE (no compute).")
# =============================================================================
# End of EFFIS validation example.
# =============================================================================
