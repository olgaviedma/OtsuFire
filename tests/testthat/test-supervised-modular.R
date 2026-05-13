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

test_that("OtsuFire:::build_supervised_training_pools validates config and options", {
  expect_error(OtsuFire:::build_supervised_training_pools(config = list()),
               regexp = "build_supervised_burned_config")
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::build_supervised_training_pools(config = cfg,
                                                write_outputs = "x"),
               regexp = "write_outputs")
  expect_error(OtsuFire:::build_supervised_training_pools(config = cfg,
                                                deterministic_decisions = 42),
               regexp = "deterministic_decisions")
  out <- OtsuFire:::build_supervised_training_pools(config = cfg)
  expect_true(is.list(out))
  expect_true("pools_gpkg" %in% names(out))
})

test_that("make_spatial_folds validates its arguments", {
  expect_error(OtsuFire:::make_spatial_folds(), regexp = "train_labelled")
  expect_error(OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m(),
                                   block_sizes_m = -1),
               regexp = "block_sizes_m")
  expect_error(OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m(),
                                   k_candidates = c(1L)),
               regexp = "k_candidates")
  cfg <- mk_cfg_m()
  out <- OtsuFire:::make_spatial_folds(train_labelled = mk_tmp_gpkg_m(), config = cfg)
  expect_true(is.list(out))
  expect_true("train_with_folds_gpkg" %in% names(out))
})

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
  out <- OtsuFire:::extract_supervised_features(
    train_with_folds = mk_tmp_gpkg_m(),
    scoring_pool = mk_tmp_gpkg_m(),
    config = cfg
  )
  expect_true(is.list(out))
  expect_true("features_geometry_gpkg" %in% names(out))
})

test_that("run_oof_diagnostics validates inputs", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::run_oof_diagnostics(), regexp = "train_features")
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m()),
    regexp = "scoring_features"
  )
  expect_error(
    OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                        scoring_features = mk_tmp_gpkg_m(),
                        fold_cols = character(0)),
    regexp = "fold_cols"
  )
  out <- OtsuFire:::run_oof_diagnostics(train_features = mk_tmp_gpkg_m(),
                              scoring_features = mk_tmp_gpkg_m(),
                              config = cfg)
  expect_true(is.list(out))
  expect_true("oof_agg_csv" %in% names(out))
})

test_that("train_final_burned_model honors deprecation and deps", {
  cfg <- mk_cfg_m()
  expect_error(OtsuFire:::train_final_burned_model(), regexp = "train_features")
  expect_warning(
    out <- OtsuFire:::train_final_burned_model(
      train_features = mk_tmp_gpkg_m(),
      hard_negative_source = "deterministic_drop_source",
      config = cfg
    ),
    regexp = "deprecated"
  )
  expect_true(is.list(out))
  expect_true("model_rds" %in% names(out))
})

test_that("score_supervised_burned_map requires model + recipe", {
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
  out <- OtsuFire:::score_supervised_burned_map(scoring_features = mk_tmp_gpkg_m(),
                                       model = list(dummy = TRUE),
                                       recipe = list(dummy = TRUE),
                                       config = cfg)
  expect_true("final_map_gpkg" %in% names(out))
})

test_that("update_burned_like_registry validates prob_threshold", {
  cfg <- mk_cfg_m()
  expect_error(update_burned_like_registry(),
               regexp = "burned_like_scored")
  expect_error(
    update_burned_like_registry(burned_like_scored = mk_tmp_gpkg_m(),
                                 burned_like_registry_path = tempfile(),
                                 prob_threshold = 1.5),
    regexp = "prob_threshold"
  )
  out <- update_burned_like_registry(
    burned_like_scored = mk_tmp_gpkg_m(),
    config = cfg
  )
  expect_true("registry_path" %in% names(out))
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
