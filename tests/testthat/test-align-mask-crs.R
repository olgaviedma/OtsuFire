# Gate 1C.1 (2026-06-08): the reusable burnable-mask/CRS alignment helper
# `.of_align_mask_to_template()`. These tests cover the core regression — a
# mask in a DIFFERENT CRS than the template used to be `resample()`d without
# `project()`, silently producing an all-NaN / zero-burnable raster. The helper
# now PROJECTS on a CRS mismatch and validates the result.

align_fn <- function() get(".of_align_mask_to_template", envir = asNamespace("OtsuFire"))

# A small binary burnable mask in EPSG:3035 (ETRS89-LAEA, metres) with a
# rectangular burnable block. Template shares geometry by default.
mk_mask_3035 <- function(ncol = 20, nrow = 20, crs = "EPSG:3035",
                         xmin = 3000000, ymin = 2000000, res = 100) {
  skip_if_not_installed("terra")
  r <- terra::rast(
    ncol = ncol, nrow = nrow,
    xmin = xmin, xmax = xmin + ncol * res,
    ymin = ymin, ymax = ymin + nrow * res,
    crs = crs
  )
  terra::values(r) <- 0
  # burnable block: a 10x10 interior square == 1
  m <- terra::as.matrix(r, wide = TRUE)
  m[6:15, 6:15] <- 1
  terra::values(r) <- as.vector(t(m))
  r
}

test_that("1. identical CRS + identical geometry returns mask unchanged", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035()
  template <- mk_mask_3035()  # same geometry

  out <- f(mask, template, allowed_values = c(0, 1), binary = TRUE)
  expect_true(out$report$non_empty)
  expect_true(out$report$crs_ok)
  expect_true(out$report$within_tol)
  # burnable preserved exactly (no reproject, no resample)
  expect_equal(out$report$n_burnable, 100)
  expect_equal(out$report$area_abs_diff, 0)
  expect_true(terra::compareGeom(out$aligned, template, stopOnError = FALSE))
})

test_that("2. identical CRS + different grid uses resample path, preserves burnable", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035(res = 100)
  # template: same CRS/extent but finer grid (50 m -> 40x40)
  template <- mk_mask_3035(ncol = 40, nrow = 40, res = 50)

  expect_false(terra::compareGeom(mask, template, stopOnError = FALSE))
  out <- f(mask, template, allowed_values = c(0, 1), binary = TRUE)
  expect_true(out$report$non_empty)
  expect_true(out$report$within_tol)
  expect_true(out$report$n_burnable > 0)
  # geometry now matches the (finer) template
  expect_true(terra::compareGeom(out$aligned, template, stopOnError = FALSE))
})

test_that("3. different CRS uses project path; burnable preserved; NOT all-NA (core regression)", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035()  # EPSG:3035

  # template covering the SAME ground footprint but in EPSG:25830 (UTM 30N, m).
  # Build it by projecting the mask's geometry, then use it purely as a grid.
  template_geom <- terra::project(mask, "EPSG:25830", method = "near")
  template <- terra::rast(template_geom)  # empty grid, same CRS/extent/res

  expect_false(terra::same.crs(mask, template))

  out <- f(mask, template, allowed_values = c(0, 1), binary = TRUE)
  expect_true(out$report$non_empty)
  expect_true(out$report$n_burnable > 0)        # NOT zeroed out
  expect_true(out$report$within_tol)
  expect_true(out$report$crs_ok)
  expect_true(terra::compareGeom(out$aligned, template, stopOnError = FALSE))

  # The whole point: a naive resample-without-project would have produced this:
  naive <- terra::resample(mask, template, method = "near")
  naive_burn <- terra::freq(naive, value = 1)
  naive_n <- if (is.null(naive_burn) || nrow(naive_burn) == 0L) 0 else sum(naive_burn[, "count"])
  expect_equal(naive_n, 0)                       # confirms the bug exists
  expect_true(out$report$n_burnable > 0)         # confirms the helper fixes it
})

test_that("4. no overlap -> stop with clear error", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035(xmin = 3000000, ymin = 2000000)
  # template far away in the SAME CRS (no extent intersection)
  template <- mk_mask_3035(xmin = 4000000, ymin = 5000000)

  expect_error(
    f(mask, template, allowed_values = c(0, 1), binary = TRUE),
    regexp = "overlap"
  )
})

test_that("5. empty mask / zero burnable cells -> stop with clear error", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035()
  terra::values(mask) <- 0  # all non-burnable, no value == 1
  template <- mk_mask_3035()

  expect_error(
    f(mask, template, allowed_values = c(0, 1), binary = TRUE),
    regexp = "ZERO burnable"
  )
})

test_that("6. invalid values in a binary mask -> stop (no silent coercion)", {
  skip_if_not_installed("terra")
  f <- align_fn()

  # value 2 present
  mask2 <- mk_mask_3035()
  m <- terra::as.matrix(mask2, wide = TRUE)
  m[1, 1] <- 2
  terra::values(mask2) <- as.vector(t(m))
  template <- mk_mask_3035()
  expect_error(
    f(mask2, template, allowed_values = c(0, 1), binary = TRUE),
    regexp = "outside the allowed set"
  )

  # fractional value 0.5 present
  maskf <- mk_mask_3035()
  mf <- terra::as.matrix(maskf, wide = TRUE)
  mf[1, 1] <- 0.5
  terra::values(maskf) <- as.vector(t(mf))
  expect_error(
    f(maskf, template, allowed_values = c(0, 1), binary = TRUE),
    regexp = "outside the allowed set"
  )
})

test_that("7. correct alignment result: dims/crs/extent match template", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035()
  # template in a DIFFERENT CRS covering the same footprint (project path);
  # near-resampling preserves burnable area within the default tolerance.
  template <- terra::rast(terra::project(mask, "EPSG:25830", method = "near"))

  out <- f(mask, template, allowed_values = c(0, 1), binary = TRUE)
  expect_true(terra::same.crs(out$aligned, template))
  expect_equal(unname(out$report$dims["nrow"]), terra::nrow(template))
  expect_equal(unname(out$report$dims["ncol"]), terra::ncol(template))
  expect_equal(out$report$res, terra::res(template))
  expect_equal(out$report$extent, as.vector(terra::ext(template)))
})

test_that("8. cross-CRS: no silent loss of burnable cells; count > 0 within tolerance", {
  skip_if_not_installed("terra")
  f <- align_fn()
  mask <- mk_mask_3035()
  area_before <- {
    cs <- terra::cellSize(mask, unit = "m")
    b <- terra::ifel(!is.na(mask) & mask == 1, 1, NA)
    terra::global(terra::ifel(!is.na(b), cs, NA), "sum", na.rm = TRUE)[1, 1]
  }

  template <- terra::rast(terra::project(mask, "EPSG:25830", method = "near"))
  out <- f(mask, template, allowed_values = c(0, 1), binary = TRUE,
           area_tol_rel = 0.05)

  expect_true(out$report$n_burnable > 0)
  expect_true(out$report$area_after > 0)
  # absolute and relative diff are recorded
  expect_true(is.finite(out$report$area_abs_diff))
  expect_true(is.finite(out$report$area_rel_diff))
  expect_true(out$report$within_tol)
  # area roughly preserved (within tolerance) vs the original footprint
  expect_lt(abs(out$report$area_after - area_before) / area_before, 0.05)
})
