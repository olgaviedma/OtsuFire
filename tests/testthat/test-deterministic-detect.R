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

# ---- change_index in-memory validation -------------------------------

# A "clean" change index: NoData declared as -9999, all values finite and
# well above any sane lower_cap.
mk_clean_ci <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 100:163)
  terra::writeRaster(r, f, overwrite = TRUE,
                     datatype = "INT4S", NAflag = -9999)
  f
}

# A change index whose minimum is -500: clean under the default
# lower_cap (-1000) but dirty if lower_cap is raised above -500.
mk_lowmin_ci <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = c(-500L, rep(50L, 63)))
  terra::writeRaster(r, f, overwrite = TRUE,
                     datatype = "INT4S", NAflag = -9999)
  f
}

test_that("change_index validation uses defaults and passes a clean raster", {
  bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(
    change_index = mk_clean_ci(), burnable_mask = bm, target_year = 2025L)
  # No change_index_validation block -> defaults (lower_cap = -1000).
  # verbose = TRUE emits informational messages by design; only assert
  # that a clean raster does not error.
  expect_no_error(suppressMessages(OtsuFire:::.of_validate_change_index(cfg)))
  res <- OtsuFire:::.of_resolve_change_index_validation(cfg)
  expect_identical(res$lower_cap, -1000)
  expect_identical(res$expected_nodata, -9999)
})

test_that("change_index validation honours an options$ lower_cap override", {
  bm <- mk_tmp_mask2()
  ci <- mk_lowmin_ci()

  # Default lower_cap (-1000): min -500 is not below it -> clean.
  cfg_def <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L)
  expect_no_error(
    suppressMessages(OtsuFire:::.of_validate_change_index(cfg_def)))

  # Override lower_cap to -100: min -500 is now below it -> dirty, and the
  # error must show the actual lower_cap used.
  cfg_ovr <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    options = list(change_index_validation = list(lower_cap = -100)))
  expect_identical(
    OtsuFire:::.of_resolve_change_index_validation(cfg_ovr)$lower_cap, -100)
  expect_error(OtsuFire:::.of_validate_change_index(cfg_ovr),
               regexp = "lower_cap = -100")
})

test_that("detect_burned_patches() aborts when change_index fails validation", {
  ci <- mk_tmp_tif2()   # synthetic raster: no declared NoData -> dirty
  bm <- mk_tmp_mask2()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L)
  expect_error(detect_burned_patches(cfg),
               regexp = "change_index.*dirty")
})
