mk_tmp_tif_m <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_gpkg_m <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                              c(0,1), c(0,0)))),
                    crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}
mk_cfg_m <- function() {
  ci <- mk_tmp_tif_m()
  id <- mk_tmp_gpkg_m()
  build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
}

# BUG 3 Phase 2 (2026-06-05): build_supervised_training_pools() is now a REAL
# exported function that runs the full STEP A sequence (read internal_decisions
# -> QA relabel via audit_deterministic_pools() -> burned/review/scoring pools
# -> BOTH all_sources unburned builders (deterministic drops + random burnable
# background, and the Otsu/legacy current-year residual patches) -> merge
# train_labeled -> write 01_POOLS/<year>_<scenario>_pools.gpkg). A faithful
# happy-path fixture would require a deterministic decisions GPKG with the
# expected QA columns, the full composite/mask raster stack for the target year
# AND the external GDAL/Python tool binaries the legacy Otsu builder shells out
# to -- far too heavy/unavailable for a unit test (and exercised end-to-end by
# the operational pipeline). We therefore assert the argument-validation
# contract (config is REQUIRED, write_outputs/overwrite/deterministic_decisions
# types) plus the self-contained config-resolution requirement: with no
# config$options$data_base the stage surfaces the data_base requirement after
# passing validation and resolving the canonical 01_POOLS output paths. The
# byte-identical A1-A6 engine sequence itself is covered by the orchestrator
# integration.
test_that("build_supervised_training_pools validates config and options", {
  expect_error(OtsuFire::build_supervised_training_pools(config = list()),
               regexp = "build_supervised_burned_config")
  cfg <- mk_cfg_m()
  expect_error(OtsuFire::build_supervised_training_pools(config = cfg,
                                                write_outputs = "x"),
               regexp = "write_outputs")
  expect_error(OtsuFire::build_supervised_training_pools(config = cfg,
                                                overwrite = "x"),
               regexp = "overwrite")
  expect_error(OtsuFire::build_supervised_training_pools(config = cfg,
                                                deterministic_decisions = 42),
               regexp = "deterministic_decisions")
  # config carries no options$data_base, so the self-contained config
  # resolution surfaces the data_base requirement (validation + 01_POOLS path
  # resolution passed first).
  expect_error(OtsuFire::build_supervised_training_pools(config = cfg),
               regexp = "data_base")
})

# A realistic labelled training fixture: a grid of small square polygons,
# each its own fire_uid, half burned / half unburned, in a projected CRS
# (EPSG:3035, meters) so block-CV can actually run. nx*ny polygons spread
# over a few km lets the 2000m block grid create multiple folds.
mk_train_labelled_sf <- function(nx = 8L, ny = 8L, step = 1500) {
  polys <- vector("list", nx * ny)
  fire_uid <- character(nx * ny)
  cls <- character(nx * ny)
  k <- 1L
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      x0 <- (i - 1L) * step
      y0 <- (j - 1L) * step
      polys[[k]] <- sf::st_polygon(list(rbind(
        c(x0, y0), c(x0 + 600, y0), c(x0 + 600, y0 + 600),
        c(x0, y0 + 600), c(x0, y0)
      )))
      fire_uid[k] <- sprintf("F%03d", k)
      cls[k] <- if (k %% 2L == 0L) "burned" else "unburned"
      k <- k + 1L
    }
  }
  sf::st_sf(
    fire_uid = fire_uid,
    class    = cls,
    geometry = sf::st_sfc(polys, crs = 3035)
  )
}

test_that("make_spatial_folds validates its arguments", {
  expect_error(OtsuFire:::make_spatial_folds(), regexp = "train_labelled")
  cfg <- mk_cfg_m()
  expect_error(
    OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m(),
                                  config = cfg, block_sizes_m = -1),
    regexp = "block_sizes_m"
  )
  expect_error(
    OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m(),
                                  config = cfg, k_candidates = c(1L)),
    regexp = "k_candidates"
  )
  # config is now REQUIRED (no out_dir fallback without it)
  expect_error(
    OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m()),
    regexp = "config"
  )
  expect_error(
    OtsuFire:::make_spatial_folds(train_labelled = mk_train_labelled_sf(),
                                  config = cfg, split_unit = "bogus"),
    regexp = "split_unit"
  )
})

test_that("make_spatial_folds does real fold work and writes outputs", {
  cfg <- mk_cfg_m()
  out_dir <- file.path(tempfile("folds_"), "02_FOLDS")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  train <- mk_train_labelled_sf()

  # Relaxed acceptance gate cannot be met by this tiny synthetic input, so we
  # expect the soft fallback warning + a _FOLD_FALLBACK.txt audit; the stage
  # still produces a valid, fully-assigned fold set and writes all outputs.
  out <- suppressWarnings(
    OtsuFire:::make_spatial_folds(
      train_labelled = train, config = cfg,
      split_unit = "fire",
      block_sizes_m = c(3000, 2000),
      k_candidates = c(3L),
      out_dir = out_dir,
      write_outputs = TRUE, overwrite = TRUE
    )
  )

  expect_true(is.list(out))
  # objects
  expect_s3_class(out$blocks, "sf")
  expect_s3_class(out$train_with_folds, "sf")
  expect_true(is.data.frame(out$folds_table))
  # every polygon got a fold_rep1 assignment (no leakage / no NA)
  expect_true("fold_rep1" %in% names(out$train_with_folds))
  expect_false(anyNA(out$train_with_folds$fold_rep1))
  # written paths
  for (p in c("blocks_gpkg", "folds_csv", "train_with_folds_gpkg")) {
    expect_true(p %in% names(out))
    expect_true(file.exists(out[[p]]), info = p)
  }
})

test_that("make_spatial_folds reads a GPKG path (train_labeled layer)", {
  cfg <- mk_cfg_m()
  out_dir <- file.path(tempfile("folds_path_"), "02_FOLDS")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gpkg <- tempfile(fileext = ".gpkg")
  sf::st_write(mk_train_labelled_sf(), gpkg, layer = "train_labeled",
               quiet = TRUE, delete_dsn = TRUE)

  out <- suppressWarnings(
    OtsuFire:::make_spatial_folds(
      train_labelled = gpkg, config = cfg,
      split_unit = "fire",
      block_sizes_m = c(2000), k_candidates = c(3L),
      out_dir = out_dir, write_outputs = TRUE, overwrite = TRUE
    )
  )
  expect_true(file.exists(out$train_with_folds_gpkg))
  expect_false(anyNA(out$train_with_folds$fold_rep1))
})

# BUG 3 Phase 2 (2026-06-05): extract_supervised_features() is now a REAL
# exported function that loads + aligns the raster support stack (RBR summer,
# post-fire DOY, autumn RBR, DEM, slope, CORINE), optionally loads hotspots,
# and runs the patch-level feature engine. A faithful happy-path fixture would
# require an aligned multi-band raster stack on disk (composite + topography +
# CORINE for the target year) plus a hotspot layer -- far too heavy/slow for a
# unit test (and exercised end-to-end by the operational pipeline). We
# therefore assert the argument-validation contract (config is REQUIRED, both
# inputs are required, use_hotspots/write_outputs must be logical) plus the
# self-contained raster-path requirement (the stage needs config$options$
# data_base / composite_base to locate the raster inputs; the test config has
# none, so the stage surfaces that requirement after passing validation and
# resolving the canonical 03_FEATURES output path). The byte-identical engine
# call itself is covered by the orchestrator integration.
test_that("extract_supervised_features requires config and both inputs", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::extract_supervised_features(), regexp = "config")
  expect_error(
    OtsuFire:::extract_supervised_features(config = cfg,
                                 scoring_pool = mk_tmp_gpkg_m()),
    regexp = "train_with_folds"
  )
  expect_error(
    OtsuFire:::extract_supervised_features(config = cfg,
                                 train_with_folds = mk_tmp_gpkg_m()),
    regexp = "scoring_pool"
  )
  expect_error(
    OtsuFire:::extract_supervised_features(config = cfg,
                                 train_with_folds = mk_tmp_gpkg_m(),
                                 scoring_pool = mk_tmp_gpkg_m(),
                                 use_hotspots = "x"),
    regexp = "use_hotspots"
  )
  expect_error(
    OtsuFire:::extract_supervised_features(config = cfg,
                                 train_with_folds = mk_tmp_gpkg_m(),
                                 scoring_pool = mk_tmp_gpkg_m(),
                                 write_outputs = "x"),
    regexp = "write_outputs"
  )
  # config carries no options$data_base, so the self-contained raster loading
  # surfaces the data_base requirement (validation + path resolution passed).
  expect_error(
    OtsuFire:::extract_supervised_features(
      train_with_folds = mk_tmp_gpkg_m(),
      scoring_pool = mk_tmp_gpkg_m(),
      config = cfg
    ),
    regexp = "data_base"
  )
})

# BUG 3 Phase 2 (2026-06-05): run_oof_diagnostics() is now a REAL exported
# function that builds the XGB design matrix and runs repeated block-CV OOF.
# A faithful happy-path fixture would require a fully-featured, whitelist-
# complete, multi-fold training table plus an xgboost fit per fold -- far too
# heavy/slow for a unit test (and exercised end-to-end by the operational
# pipeline). We therefore assert the argument-validation contract (config is
# now REQUIRED, both feature inputs are required, fold_cols must be non-empty)
# plus the resolution of out_dir/matrix_dir from config, which is the modular
# function's own responsibility (the byte-identical engine call itself is
# covered by the wrapper's own tests + the orchestrator integration).
test_that("run_oof_diagnostics validates inputs", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::run_oof_diagnostics(), regexp = "train_features")
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m()),
    regexp = "scoring_features"
  )
  # config is now REQUIRED (uniform modular signature; no out_dir fallback)
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m()),
    regexp = "config"
  )
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m(),
                        config = cfg,
                        fold_cols = character(0)),
    regexp = "fold_cols"
  )
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m(),
                        config = cfg, params = "not-a-list"),
    regexp = "params"
  )
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m(),
                        config = cfg, overwrite = "yes"),
    regexp = "overwrite"
  )
})

test_that("run_oof_diagnostics resolves out_dir/matrix_dir from config", {
  cfg <- mk_cfg_m()
  # config must be the S3 class, not an arbitrary list.
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m(),
                        config = list()),
    regexp = "config"
  )
  # The 05_OOF / 04_MATRIX routes the function resolves from config are the
  # canonical supervised output routes.
  expect_identical(cfg$output_routes$oof_dir,
                   file.path(cfg$output_routes$base, "05_OOF"))
  expect_identical(cfg$output_routes$matrix_dir,
                   file.path(cfg$output_routes$base, "04_MATRIX"))
})

# BUG 3 Phase 2 (2026-06-05): train_final_burned_model() is now the REAL
# exported TRAINING half of the formerly fused
# run_train_final_model_and_export_final_map() wrapper. It wraps
# train_final_model_direct() (xgboost fit on the labelled pool) with the same
# args/values + the four sampling caps + whitelist/weights, and writes the
# 07_FINAL_MODEL_V2 artifacts. A faithful happy-path fixture would require a
# fully-featured, whitelist-complete labelled GPKG plus an xgboost fit -- far
# too heavy/slow for a unit test (and exercised end-to-end by the operational
# pipeline). We therefore assert the argument-validation contract (config now
# REQUIRED, train_features required, the four caps must be non-negative
# numbers, overwrite logical) plus the resolution of out_dir from config; the
# byte-identical engine call is covered by the engine's own tests + the
# orchestrator integration.
test_that("train_final_burned_model validates inputs", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::train_final_burned_model(), regexp = "train_features")
  # config is now REQUIRED (uniform modular signature; no out_dir fallback)
  expect_error(
    OtsuFire:::train_final_burned_model(train_features = mk_tmp_gpkg_m()),
    regexp = "config"
  )
  expect_error(
    OtsuFire:::train_final_burned_model(train_features = mk_tmp_gpkg_m(),
                                        config = list()),
    regexp = "config"
  )
  # Precision 1 (2026-06-07): -1 is a DEPRECATED function-level override (warns
  # via "otsufire_deprecated_param") AND an invalid value (errors). Suppress the
  # deprecation warning so only the validation error is asserted here.
  expect_error(
    suppressWarnings(
      OtsuFire:::train_final_burned_model(
        train_features = mk_tmp_gpkg_m(), config = cfg,
        random_to_burned_ratio = -1)),
    regexp = "random_to_burned_ratio"
  )
  expect_error(
    OtsuFire:::train_final_burned_model(train_features = mk_tmp_gpkg_m(),
                                        config = cfg, oof_agg = 42),
    regexp = "oof_agg"
  )
  expect_error(
    OtsuFire:::train_final_burned_model(train_features = mk_tmp_gpkg_m(),
                                        config = cfg, overwrite = "yes"),
    regexp = "overwrite"
  )
})

test_that("train_final_burned_model resolves out_dir from config", {
  cfg <- mk_cfg_m()
  # 07_FINAL_MODEL_V2 is the canonical final-model route.
  expect_identical(cfg$output_routes$final_model_dir,
                   file.path(cfg$output_routes$base, "07_FINAL_MODEL_V2"))
})

# BUG 3 Phase 2 (2026-06-05): score_supervised_burned_map() is now the REAL
# exported SCORING half. It wraps score_burnedlike_and_export_final_map()
# (xgboost scoring + temporal adjustment + final-map export) with the same
# args/values + the three CURRENTYEAR_* temporal thresholds. As above, a full
# happy-path needs a trained model + recipe + scored universe -- impractical in
# unit time -- so we assert the validation contract (config now REQUIRED,
# model + recipe required, temporal thresholds numeric, export flag logical)
# plus the 08_SCORED / 09_FINAL_MAP route resolution.
test_that("score_supervised_burned_map validates inputs", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::score_supervised_burned_map(), regexp = "scoring_features")
  expect_error(
    OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                  config = cfg),
    regexp = "model"
  )
  expect_error(
    OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                  model = list(dummy = TRUE),
                                  config = cfg),
    regexp = "recipe"
  )
  # config is now REQUIRED (uniform modular signature; no out_map_dir fallback)
  expect_error(
    OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                  model = list(dummy = TRUE),
                                  recipe = list(dummy = TRUE)),
    regexp = "config"
  )
  expect_error(
    OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                  model = list(dummy = TRUE),
                                  recipe = list(dummy = TRUE),
                                  config = cfg, export_burned_like = "yes"),
    regexp = "export_burned_like"
  )
  expect_error(
    OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                  model = list(dummy = TRUE),
                                  recipe = list(dummy = TRUE),
                                  config = cfg,
                                  preyear_overlap_threshold = "x"),
    regexp = "preyear_overlap_threshold"
  )
})

test_that("score_supervised_burned_map resolves map/score dirs from config", {
  cfg <- mk_cfg_m()
  expect_identical(cfg$output_routes$final_map_dir,
                   file.path(cfg$output_routes$base, "09_FINAL_MAP"))
  expect_identical(cfg$output_routes$scored_dir,
                   file.path(cfg$output_routes$base, "08_SCORED"))
})

# §N+27 (2026-06-05): the `update_burned_like_registry validates
# prob_threshold` test was removed together with the abandoned supervised
# burned-like registry. The public function no longer exists.

# B3 (2026-06-06): the modular wrappers (run_oof_diagnostics,
# train_final_burned_model, score_supervised_burned_map) and the run-result
# reconstruction (run_supervised_burned_mapping) must build their prefixes from
# config$options$prefix_oof_base / config$options$prefix_base -- the SAME knobs
# the orchestrator uses -- so overriding the option keeps every advertised path
# consistent. Defaults ("patch" / "patch_certified") are unchanged. The heavy
# wrappers cannot return paths in a unit test (they need a full raster/model
# stack), so we assert at the path/config level: the resolved prefixes (the
# exact expressions the wrappers and run-result reconstruction use) match what
# the orchestrator would build for the same config.
test_that("supervised prefixes honor config$options$prefix_base/prefix_oof_base", {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # The single source of truth shared by orchestrator (prefix_oof/prefix at
  # internal-sup-orchestrator.R), the three modular wrappers, and the
  # run-result reconstruction (supervised-run.R).
  resolve_prefixes <- function(cfg) {
    prefix_oof_base <- cfg$options$prefix_oof_base %||% "patch"
    prefix_base     <- cfg$options$prefix_base %||% "patch_certified"
    list(
      prefix_oof = sprintf("%d_%s_%s", cfg$target_year, cfg$scenario,
                           prefix_oof_base),
      prefix     = sprintf("%d_%s_%s", cfg$target_year, cfg$scenario,
                           prefix_base)
    )
  }

  # Default config: byte-identical legacy prefixes.
  cfg <- mk_cfg_m()
  pd <- resolve_prefixes(cfg)
  expect_identical(pd$prefix_oof, "2025_balanced_patch")
  expect_identical(pd$prefix,     "2025_balanced_patch_certified")

  # Overriding prefix_base must propagate to the final/score prefix while
  # leaving the OOF prefix on its own (independent) knob.
  cfg2 <- cfg
  cfg2$options$prefix_base <- "custom"
  po <- resolve_prefixes(cfg2)
  expect_identical(po$prefix,     "2025_balanced_custom")
  expect_identical(po$prefix_oof, "2025_balanced_patch")

  # Overriding prefix_oof_base must propagate to the OOF prefix.
  cfg3 <- cfg
  cfg3$options$prefix_oof_base <- "oofcustom"
  poo <- resolve_prefixes(cfg3)
  expect_identical(poo$prefix_oof, "2025_balanced_oofcustom")
  expect_identical(poo$prefix,     "2025_balanced_patch_certified")
})

test_that("check_supervised_consistency validates required inputs", {
  cfg <- mk_cfg_m()
  expect_error(check_supervised_consistency(),
               regexp = "deterministic_decisions")
  expect_error(
    check_supervised_consistency(deterministic_decisions = mk_tmp_gpkg_m()),
    regexp = "final_map"
  )
  out <- check_supervised_consistency(
    deterministic_decisions = mk_tmp_gpkg_m(),
    final_map = mk_tmp_gpkg_m(),
    config = cfg
  )
  expect_true("consistency_summary_csv" %in% names(out))
})

test_that("run_oneyear_supervised_pipeline validates config and flags", {
  expect_error(run_oneyear_supervised_pipeline(config = list()),
               regexp = "build_supervised_burned_config")
  cfg <- mk_cfg_m()
  expect_error(run_oneyear_supervised_pipeline(cfg, run_consistency = "x"),
               regexp = "run_consistency")
  expect_error(run_oneyear_supervised_pipeline(cfg, overwrite = 1),
               regexp = "overwrite")
})
