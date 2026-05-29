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
mk_keep_pool_summary4 <- function(years_used = NULL) {
  md <- if (is.null(years_used)) NULL else
    list(years_used = as.integer(years_used))
  structure(
    list(
      min_area_ha = 10,
      min_pix = 13L,
      promote_percentile = 0.10,
      promote_p_above_ref = 0.50,
      qref_prob = 0.05,
      qref_keep = 250,
      q25_keep = 300,
      pixel_area_ha = 0.81,
      keep_medians = c(260, 310, 400),
      metadata = md
    ),
    class = c("otsufire_keep_pool", "list")
  )
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

test_that("run_deterministic_pipeline rejects same-year explicit keep_pool before detection", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(
    run_deterministic_pipeline(cfg, keep_pool = mk_keep_pool_summary4(2025L)),
    regexp = "cannot include the target year"
  )
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

# ---- P1-DET-03: (write_outputs=FALSE, run_validation=TRUE) rejected -------
# The shared validator (validate_fire_maps) consumes the on-disk decision
# layer produced when write_outputs = TRUE. The previous behavior silently
# returned NULL from .of_run_shared_validation() when that layer was not
# materialized. The pipeline now refuses the combination explicitly.

test_that("run_deterministic_pipeline rejects (write_outputs=FALSE, run_validation=TRUE) early", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L
  )
  err <- tryCatch(
    run_deterministic_pipeline(
      cfg, write_outputs = FALSE, run_validation = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "write_outputs = FALSE", fixed = TRUE)
  expect_match(err, "run_validation = TRUE", fixed = TRUE)
  # Must fire before any detection/scoring work.
  expect_false(grepl("change_index", err, fixed = TRUE))
})

test_that("run_deterministic_pipeline: 'auto' + reference + write_outputs=FALSE also rejected", {
  ci <- mk_tmp_tif4(); bm <- mk_tmp_mask4()
  ref_f <- tempfile(fileext = ".gpkg"); file.create(ref_f)
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm,
    reference_burned_map = ref_f, target_year = 2025L
  )
  # 'auto' resolves to TRUE because reference_burned_map is set; combined
  # with write_outputs=FALSE this should fail at the early-validation gate.
  err <- tryCatch(
    run_deterministic_pipeline(
      cfg, write_outputs = FALSE, run_validation = "auto"
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "write_outputs = FALSE", fixed = TRUE)
})

# ---- P2-DET-01: workbook_path is wired to validate_fire_maps' excel_path ---
# We verify two complementary things without spinning up the whole engine:
#   (a) the shared-validation helper does request write_excel = isTRUE(
#       write_outputs) from validate_fire_maps() and threads the
#       excel_filename derived from output_routes$validation_workbook,
#       so workbook_path is no longer a hardcoded NA;
#   (b) the assembled result_paths$validation_workbook field corresponds
#       to validation$workbook_path (no longer NA when validation runs).

test_that(".of_run_shared_validation body wires write_excel and excel_path", {
  body_src <- paste(
    deparse(body(OtsuFire:::.of_run_shared_validation)),
    collapse = "\n"
  )
  # write_excel must be propagated from isTRUE(write_outputs).
  expect_true(grepl("write_excel\\s*=\\s*isTRUE\\(write_outputs\\)",
                    body_src))
  # excel_filename must be derived from output_routes$validation_workbook.
  expect_true(grepl("validation_workbook", body_src))
  # workbook_path must come from the validator's excel_path, not a
  # hardcoded NA_character_.
  expect_true(grepl("excel_path", body_src))
  # Sanity: the old "workbook_path   = NA_character_" hardcode is gone.
  expect_false(grepl("workbook_path\\s*=\\s*NA_character_(?![^)]*excel)",
                     body_src, perl = TRUE))
})

test_that("run_deterministic_pipeline body threads workbook_path through validation", {
  body_src <- paste(deparse(body(OtsuFire::run_deterministic_pipeline)),
                    collapse = "\n")
  # result_paths$validation_workbook must come from validation$workbook_path.
  expect_true(grepl("validation_workbook\\s*=", body_src))
  expect_true(grepl("validation\\$workbook_path", body_src))
})

# §N+7.5 — surfaced by the 2005/balanced deterministic smoke (2026-05-29).
# The shared validator's `mask_shapefile` argument expects a study-area
# boundary polygon. The original wrapper passed the burnable RASTER for
# both `burnable_path` and `mask_path`, which broke validate_fire_maps()
# with "Cannot open ...tif" — the error got caught by the outer tryCatch
# and the warning was buried in R's end-of-run "50 or more warnings"
# batch, so the user saw `validation_workbook = NA` without any
# actionable signal. These tests pin the fix.

test_that(".of_run_shared_validation body sources mask_path from peninsula_shapefile_path", {
  body_src <- paste(
    deparse(body(OtsuFire:::.of_run_shared_validation)),
    collapse = "\n"
  )
  # mask_path must read from config$options$peninsula_shapefile_path.
  expect_true(grepl("peninsula_shapefile_path", body_src))
  expect_true(grepl("mask_path\\s*<-\\s*peninsula_path", body_src))
  # The old aliasing of mask_path to the burnable spec must be gone.
  expect_false(grepl(
    "mask_path\\s*<-\\s*\\.of_input_to_path\\(burnable_spec",
    body_src
  ))
  # When peninsula_shapefile_path is missing, validation must stop()
  # with a message mentioning the option key (the outer tryCatch then
  # turns it into a visible warning).
  expect_true(grepl("peninsula_shapefile_path", body_src))
  expect_true(grepl("stop\\(", body_src))
})

test_that("run_deterministic_pipeline body marks shared-validation warning as immediate", {
  body_src <- paste(deparse(body(OtsuFire::run_deterministic_pipeline)),
                    collapse = "\n")
  # The tryCatch around .of_run_shared_validation must emit its warning
  # with immediate. = TRUE so it survives R's batched end-of-run output
  # ("There were 50 or more warnings").
  expect_true(grepl("Shared validation failed", body_src))
  expect_true(grepl("immediate\\.\\s*=\\s*TRUE", body_src))
})
