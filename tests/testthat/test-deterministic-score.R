# Tests for score_burned_patches(). Input validation only — no engine runs.

mk_tmp_tif3 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_mask3 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("score_burned_patches exists and is no longer NotYetImplemented", {
  expect_true(is.function(score_burned_patches))
  err <- tryCatch(score_burned_patches(), error = function(e) conditionMessage(e))
  expect_false(grepl("not implemented yet", err, ignore.case = TRUE))
})

test_that("score_burned_patches rejects non-config", {
  expect_error(score_burned_patches(burned_candidates = "x.gpkg",
                                    config = list()),
               regexp = "build_burned_mapping_config")
})

test_that("score_burned_patches requires burned_candidates", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(score_burned_patches(config = cfg),
               regexp = "burned_candidates")
})

test_that("score_burned_patches fails when burned_candidates path is missing", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(
    score_burned_patches("nowhere/at/all.gpkg", cfg),
    regexp = "burned_candidates.*does not exist"
  )
})

test_that("score_burned_patches rejects bad burned_candidates types", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(
    score_burned_patches(burned_candidates = 42, config = cfg),
    regexp = "burned_candidates"
  )
})

test_that("score_burned_patches rejects config with NULL change_index", {
  bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(
    change_index = NULL, burnable_mask = bm, target_year = 2025L
  )
  # Any plausible burned_candidates placeholder will fail at the
  # change_index check first.
  expect_error(
    score_burned_patches(burned_candidates = "x.gpkg", config = cfg),
    regexp = "change_index"
  )
})

test_that("score_burned_patches rejects bad keep_pool types", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  # Pass an existing path so we get past the burned_candidates check and hit
  # the keep_pool type check.
  fake_candidates <- tempfile(fileext = ".gpkg"); file.create(fake_candidates)
  expect_error(
    score_burned_patches(fake_candidates, cfg, keep_pool = 1:10),
    regexp = "keep_pool"
  )
})
