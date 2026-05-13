# Tests for detect_burned_patches(). Input validation only — no engine runs.

mk_tmp_tif2 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mk_tmp_mask2 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("detect_burned_patches exists and is no longer NotYetImplemented", {
  expect_true(is.function(detect_burned_patches))
  err <- tryCatch(detect_burned_patches(), error = function(e) conditionMessage(e))
  expect_false(grepl("not implemented yet", err, ignore.case = TRUE))
})

test_that("detect_burned_patches rejects a non-config object", {
  expect_error(detect_burned_patches(list(scenario = "balanced")),
               regexp = "build_burned_mapping_config")
})

test_that("detect_burned_patches rejects a config with NULL change_index", {
  bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(
    change_index = NULL, burnable_mask = bm, target_year = 2025L
  )
  expect_error(detect_burned_patches(cfg), regexp = "change_index")
})

test_that("detect_burned_patches fails cleanly when change_index file is missing", {
  bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(
    change_index = "does/not/exist.tif",
    burnable_mask = bm,
    target_year = 2025L
  )
  expect_error(detect_burned_patches(cfg),
               regexp = "change_index.*does not exist")
})

test_that("detect_burned_patches rejects non-logical write_outputs/overwrite", {
  ci <- mk_tmp_tif2(); bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(detect_burned_patches(cfg, write_outputs = "yes"),
               regexp = "write_outputs")
  expect_error(detect_burned_patches(cfg, overwrite = 1),
               regexp = "overwrite")
})

test_that("detect_burned_patches rejects invalid aoi inputs", {
  ci <- mk_tmp_tif2(); bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L,
                                      options = list(engine_root = tempdir()))
  expect_error(detect_burned_patches(cfg, aoi = 42),
               regexp = "aoi")
})
