# Tests for build_burned_mapping_config(). Lightweight — no engine runs.

mk_tmp_tif <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mk_tmp_mask <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("config builder returns the expected S3 class and required fields", {
  ci <- mk_tmp_tif()
  bm <- mk_tmp_mask()

  cfg <- build_burned_mapping_config(
    scenario = "balanced",
    change_index = ci,
    burnable_mask = bm,
    target_year = 2025,
    output_dir = tempdir()
  )

  expect_s3_class(cfg, "otsufire_burned_mapping_config")
  expect_identical(cfg$scenario, "balanced")
  expect_identical(cfg$target_year, 2025L)
  expect_identical(cfg$run_name, "deterministic_burned_map")
  expect_type(cfg$inputs, "list")
  expect_named(cfg$inputs,
               c("change_index", "vegetation_map", "burnable_mask",
                 "hotspots", "previous_year_burned", "reference_burned_map"))
  expect_identical(cfg$deterministic_seed, 12345L)
  expect_type(cfg$output_routes, "list")
  expect_true(all(c("grow_vector", "refine_merged_gpkg",
                    "internal_decisions", "timing_csv") %in%
                  names(cfg$output_routes)))
})

test_that("config builder respects scenario choices and explicit overrides", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()

  cfg <- build_burned_mapping_config(
    scenario = "lax",
    change_index = ci,
    burnable_mask = bm,
    target_year = 2012L,
    run_name = "custom_run",
    output_dir = tempdir(),
    options = list(deterministic_seed = 7L)
  )

  expect_identical(cfg$scenario, "lax")
  expect_identical(cfg$run_name, "custom_run")
  expect_identical(cfg$deterministic_seed, 7L)
})

test_that("config builder rejects invalid scenario", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      scenario = "ultra",
      change_index = ci,
      burnable_mask = bm,
      target_year = 2025L
    ),
    regexp = "should be one of"
  )
})

test_that("config builder rejects missing target_year and invalid year values", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(change_index = ci, burnable_mask = bm),
    regexp = "target_year"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 1500L
    ),
    regexp = "target_year"
  )
})

test_that("config builder rejects missing burnable_mask", {
  ci <- mk_tmp_tif()
  expect_error(
    build_burned_mapping_config(change_index = ci, target_year = 2025L),
    regexp = "burnable_mask"
  )
})

test_that("config builder rejects unnamed options list entries", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      options = list(1, 2, 3)
    ),
    regexp = "named list"
  )
})

test_that("config builder accepts NULL change_index at build time", {
  bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = NULL, burnable_mask = bm, target_year = 2025L
  )
  expect_s3_class(cfg, "otsufire_burned_mapping_config")
  expect_null(cfg$inputs$change_index)
})

test_that("output routes respect run_name, year, and scenario", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    scenario = "original", change_index = ci, burnable_mask = bm,
    target_year = 1985L, run_name = "smoke_x",
    output_dir = tempdir()
  )
  expect_match(cfg$output_routes$base,
               "1985/smoke_x/DETERMINISTIC/original", fixed = FALSE)
  expect_match(cfg$output_routes$internal_decisions,
               "internal_decisions\\.gpkg$")
})
