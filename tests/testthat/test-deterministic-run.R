# Tests for run_deterministic_pipeline(). Plumbing only — no engine run.

mk_tmp_tif4 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_mask4 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("run_deterministic_pipeline exists and is no longer NotYetImplemented", {
  expect_true(is.function(run_deterministic_pipeline))
  expect_true("keep_pool" %in% names(formals(run_deterministic_pipeline)))
  err <- tryCatch(run_deterministic_pipeline(),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("not implemented yet", err, ignore.case = TRUE))
})

test_that("run_deterministic_pipeline rejects non-config", {
  expect_error(run_deterministic_pipeline(config = list()),
               regexp = "build_burned_mapping_config")
})

test_that("run_deterministic_pipeline rejects invalid run_validation", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(run_deterministic_pipeline(cfg, run_validation = "maybe"),
               regexp = "run_validation")
  expect_error(run_deterministic_pipeline(cfg, run_validation = 1),
               regexp = "run_validation")
})

test_that("run_deterministic_pipeline rejects non-logical flags", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(run_deterministic_pipeline(cfg, write_outputs = "x"),
               regexp = "write_outputs")
  expect_error(run_deterministic_pipeline(cfg, overwrite = 7),
               regexp = "overwrite")
})

test_that("run_deterministic_pipeline validates keep_pool before engine work", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(run_deterministic_pipeline(cfg, keep_pool = 1:10),
               regexp = "keep-pool summary list")
  expect_error(run_deterministic_pipeline(cfg, keep_pool = list(foo = 1)),
               regexp = "keep-pool summary list")
})

test_that("run_deterministic_pipeline threads keep_pool into scoring", {
  body_src <- paste(deparse(body(OtsuFire::run_deterministic_pipeline)),
                    collapse = "\n")
  expect_true(grepl("keep_pool\\s*=\\s*keep_pool", body_src))
  expect_false(grepl("keep_pool\\s*=\\s*NULL", body_src))
})

test_that("run_deterministic_pipeline 'auto' resolves correctly from config", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg_no_ref <- build_burned_mapping_config(change_index = ci,
                                             burnable_mask = bm,
                                             target_year = 2025L)
  # We can call the internal resolver through the test env to check the
  # decision without running the engine.
  env <- asNamespace("OtsuFire")
  expect_false(
    env$.of_resolve_run_validation("auto", cfg_no_ref)
  )

  ref_f <- tempfile(fileext = ".gpkg")
  file.create(ref_f)
  cfg_ref <- build_burned_mapping_config(change_index = ci,
                                          burnable_mask = bm,
                                          reference_burned_map = ref_f,
                                          target_year = 2025L)
  expect_true(
    env$.of_resolve_run_validation("auto", cfg_ref)
  )
})
