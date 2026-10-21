# Tests for the optional DOY split in build_burned_mapping_config() and run_deterministic_pipeline().
# Detection and scoring are replaced by stubs; fixture in helper-doy-split.R.

# ---- configuration -------------------------------------------------------

mk_cfg_inputs <- function() {
  ci <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 8, nrow = 8, vals = 1:64), ci)
  bm <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32)), bm)
  list(ci = ci, bm = bm)
}

test_that("the DOY split is off by default and leaves refine_params unchanged", {
  i <- mk_cfg_inputs()
  mk <- function(...) build_burned_mapping_config(change_index = i$ci, burnable_mask = i$bm,
                                                  target_year = 2025L, output_dir = tempdir(),
                                                  run_name = "t", ...)
  base <- mk()
  expect_false(base$doy_split_params$split_by_doy)
  expect_identical(base$doy_split_params$gap_days, 30)
  expect_identical(base$doy_split_params$min_subpatch_px, 9L)
  on <- mk(refine_params = list(aoi_buffer_m = 1000, split_by_doy = TRUE, doy_gap_days = 15))
  expect_true(on$doy_split_params$split_by_doy)
  expect_identical(on$doy_split_params$gap_days, 15)
  # engine settings do not carry the DOY keys
  expect_identical(names(on$refine_params), names(base$refine_params))
  expect_identical(mk(refine_params = list(split_by_doy = FALSE))$refine_params, base$refine_params)
  # unknown keys are still rejected
  expect_error(mk(refine_params = list(bogus = 1)), "Unknown refine_params")
})

test_that("invalid DOY split settings are rejected by the configuration", {
  i <- mk_cfg_inputs()
  mk <- function(rp) build_burned_mapping_config(change_index = i$ci, burnable_mask = i$bm,
                                                 target_year = 2025L, refine_params = rp)
  expect_error(mk(list(split_by_doy = "yes")), "split_by_doy")
  expect_error(mk(list(split_by_doy = NA)), "split_by_doy")
  expect_error(mk(list(split_by_doy = TRUE, doy_gap_days = 0)), "doy_gap_days")
  expect_error(mk(list(doy_min_group_share = 2)), "min_group_share")
})

# ---- pipeline wiring (detection and scoring replaced by stubs) -------------

mk_pipeline_stubs <- function(fx) {
  captured <- new.env()
  stub_detect <- function(config, aoi = NULL, write_outputs = TRUE, overwrite = FALSE) {
    list(grown_patches_path = fx$gpkg, refined_patches_path = fx$gpkg,
         detection_diagnostics = NULL)
  }
  stub_score <- function(burned_candidates, config, keep_pool = NULL,
                         write_outputs = TRUE, overwrite = FALSE) {
    captured$candidates <- burned_candidates
    list(internal_decisions_path = NA_character_, reference_decisions_path = NA_character_)
  }
  list(detect = stub_detect, score = stub_score, captured = captured)
}

test_that("with split_by_doy = FALSE the pipeline scores exactly the detected candidates", {
  fx <- mk_doy_fixture(); i <- mk_cfg_inputs()
  st <- mk_pipeline_stubs(fx)
  mk <- function(...) build_burned_mapping_config(change_index = fx$ci, burnable_mask = i$bm,
                                                  target_year = 2006L, output_dir = tempfile("out"),
                                                  run_name = "t", ...)
  testthat::local_mocked_bindings(detect_burned_patches = st$detect,
                                  score_burned_patches = st$score)
  res_default <- run_deterministic_pipeline(mk(), run_validation = FALSE)
  cand_default <- st$captured$candidates
  res_off <- run_deterministic_pipeline(mk(refine_params = list(split_by_doy = FALSE)),
                                        run_validation = FALSE)
  expect_identical(cand_default, fx$gpkg)
  expect_identical(st$captured$candidates, fx$gpkg)
  expect_identical(names(res_off), names(res_default))
  expect_identical(names(res_off$result_paths), names(res_default$result_paths))
  expect_false("doy_split" %in% names(res_off))
  expect_false("doy_split" %in% res_off$timing_log$step)
})

test_that("with split_by_doy = TRUE the pipeline scores the split layer and audits it", {
  fx <- mk_doy_fixture(); i <- mk_cfg_inputs()
  st <- mk_pipeline_stubs(fx)
  cfg <- build_burned_mapping_config(change_index = fx$ci, burnable_mask = i$bm,
                                     target_year = 2006L, output_dir = tempfile("out"),
                                     run_name = "t",
                                     refine_params = list(split_by_doy = TRUE))
  testthat::local_mocked_bindings(detect_burned_patches = st$detect,
                                  score_burned_patches = st$score)
  res <- run_deterministic_pipeline(cfg, run_validation = FALSE)
  expect_false(identical(st$captured$candidates, fx$gpkg))
  expect_match(st$captured$candidates, "DOYSPLIT/dsplit_30d", fixed = TRUE)
  expect_equal(nrow(sf::st_read(st$captured$candidates, quiet = TRUE)), 3L)
  expect_equal(res$doy_split$n_mixed, 1L)
  expect_true(file.exists(res$result_paths$doy_split_audit))
  expect_true("doy_split" %in% res$timing_log$step)
})
