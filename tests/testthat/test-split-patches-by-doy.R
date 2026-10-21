# Tests for split_patches_by_doy() (fixture in helper-doy-split.R).


test_that("a mixed patch is split into one sub-patch per date group", {
  fx <- mk_doy_fixture()
  res <- split_patches_by_doy(fx$cand, fx$ci, year = 2006L, verbose = FALSE)
  expect_equal(res$n_examined, 2L)
  expect_equal(res$n_mixed, 1L)
  expect_equal(res$n_subpatches, 2L)
  expect_equal(nrow(res$patches), 3L)            # 1 untouched + 2 sub-patches
  expect_identical(res$patches$CLUMP_ID, 1:3)
  pa <- res$patches_audit
  expect_setequal(pa$doy_med, c(168, 240))
  expect_true(all(pa$parent_id == "1"))
  # the two halves cover the parent (rasterised footprint, 162 ha each)
  expect_equal(sort(as.numeric(sf::st_area(res$patches[pa$CLUMP_ID, ])) / 1e4),
               c(162, 162), tolerance = 0.01)
  # no output is written for an in-memory input without out_dir
  expect_true(is.na(res$output_path))
})

test_that("a patch with homogeneous dates is not split and keeps its geometry", {
  fx <- mk_doy_fixture()
  res <- split_patches_by_doy(fx$cand[2, ], fx$ci, year = 2006L, verbose = FALSE)
  expect_equal(res$n_mixed, 0L)
  expect_equal(res$n_subpatches, 0L)
  expect_equal(nrow(res$patches), 1L)
  expect_true(sf::st_equals(sf::st_geometry(res$patches),
                            sf::st_geometry(fx$cand[2, ]), sparse = FALSE)[1, 1])
})

test_that("dates closer than doy_gap_days are not split", {
  fx <- mk_doy_fixture(doy_left = 200, doy_right = 225)   # 25 days apart
  expect_equal(split_patches_by_doy(fx$cand, fx$ci, doy_gap_days = 30,
                                    verbose = FALSE)$n_mixed, 0L)
  expect_equal(split_patches_by_doy(fx$cand, fx$ci, doy_gap_days = 15,
                                    verbose = FALSE)$n_mixed, 1L)
})

test_that("candidates below min_area_ha are not examined", {
  fx <- mk_doy_fixture()
  res <- split_patches_by_doy(fx$cand, fx$ci, min_area_ha = 400, verbose = FALSE)
  expect_equal(res$n_examined, 0L)
  expect_equal(nrow(res$patches), 2L)
})

test_that("outputs are written under DOYSPLIT, audited, and never overwritten", {
  fx <- mk_doy_fixture()
  res <- split_patches_by_doy(fx$gpkg, fx$ci, doy_gap_days = 30, verbose = FALSE)
  expect_match(res$output_path, "DOYSPLIT/dsplit_30d/BA_2006_REFINE_MERGED_dsplit30d.gpkg", fixed = TRUE)
  expect_true(all(file.exists(c(res$output_path, res$audit_path, res$patches_audit_path,
                                res$parents_audit_path, res$params_path))))
  audit <- utils::read.csv(res$audit_path)
  expect_equal(audit$n_mixed, 1L)
  expect_equal(audit$doy_gap_days, 30)
  expect_equal(nrow(sf::st_read(res$output_path, quiet = TRUE)), 3L)
  # the original refined layer is preserved
  expect_equal(nrow(sf::st_read(fx$gpkg, quiet = TRUE)), 2L)
  expect_error(split_patches_by_doy(fx$gpkg, fx$ci, verbose = FALSE),
               "overwrite = FALSE")
  expect_no_error(split_patches_by_doy(fx$gpkg, fx$ci, overwrite = TRUE, verbose = FALSE))
})

test_that("parameters are validated", {
  fx <- mk_doy_fixture()
  f <- function(...) split_patches_by_doy(fx$cand, fx$ci, verbose = FALSE, ...)
  expect_error(f(doy_gap_days = 0), "doy_gap_days")
  expect_error(f(doy_gap_days = -5), "doy_gap_days")
  expect_error(f(doy_gap_days = NA_real_), "doy_gap_days")
  expect_error(f(doy_gap_days = "30"), "doy_gap_days")
  expect_error(f(min_group_share = 1), "min_group_share")
  expect_error(f(min_group_share = -0.1), "min_group_share")
  expect_error(f(min_subpatch_px = 0), "min_subpatch_px")
  expect_error(f(min_subpatch_px = 2.5), "min_subpatch_px")
  expect_error(f(min_area_ha = -1), "min_area_ha")
  expect_error(f(overwrite = NA), "overwrite")
  expect_error(split_patches_by_doy(fx$cand, "no_such_file.tif"), "change_index")
  one_band <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(fx$ci)[[1]], one_band)
  expect_error(split_patches_by_doy(fx$cand, one_band), "at least 2 bands")
  expect_error(split_patches_by_doy("no_such.gpkg", fx$ci), "refine_path")
})

