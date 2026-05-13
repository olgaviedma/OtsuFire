# Tests for the metrics extension of validate_fire_maps() (0.2.1).
# Synthetic 10x10 grid (30 m cells) with two reference polygons and one
# detection polygon designed so that TP / FP / FN / TN and per-polygon
# coverage are known by hand.
#
# Layout (rows 1-5 = top half, rows 6-10 = bottom half):
#   ref1 = top-left  quadrant (rows 1-5 x cols 1-5)  -> 25 cells
#   ref2 = top-right quadrant (rows 1-5 x cols 6-10) -> 25 cells
#   pred = top-left  + first col of top-right        -> 30 cells
# So:
#   TP = 30 (pred fully inside ref1 + 5 cells inside ref2)
#   FP = 0  (pred is fully within ref1 + ref2 union)
#   FN = 20 (rest of ref2 not covered by pred)
#   TN = 50 (entire bottom half)
# And per-reference coverage:
#   ref1 -> 100%
#   ref2 ->  20%

build_synth_validation_fixture <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  burnable <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 0, xmax = 300, ymin = 0, ymax = 300,
    crs = "EPSG:3035"
  )
  terra::values(burnable) <- 1
  burnable_path <- file.path(dir, "burnable.tif")
  terra::writeRaster(burnable, burnable_path, overwrite = TRUE,
                     datatype = "FLT4S")

  mask_poly <- sf::st_polygon(list(matrix(c(
    0,   0,
    300, 0,
    300, 300,
    0,   300,
    0,   0
  ), ncol = 2, byrow = TRUE)))
  mask_sf <- sf::st_sf(
    id = 1L,
    geometry = sf::st_sfc(mask_poly, crs = 3035)
  )
  mask_path <- file.path(dir, "mask.shp")
  suppressWarnings(sf::st_write(mask_sf, mask_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  ref1_poly <- sf::st_polygon(list(matrix(c(
    1,   151,
    149, 151,
    149, 299,
    1,   299,
    1,   151
  ), ncol = 2, byrow = TRUE)))
  ref2_poly <- sf::st_polygon(list(matrix(c(
    151, 151,
    299, 151,
    299, 299,
    151, 299,
    151, 151
  ), ncol = 2, byrow = TRUE)))
  ref_sf <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(ref1_poly, ref2_poly, crs = 3035)
  )
  ref_path <- file.path(dir, "ref.shp")
  suppressWarnings(sf::st_write(ref_sf, ref_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  pred_poly <- sf::st_polygon(list(matrix(c(
    1,   151,
    179, 151,
    179, 299,
    1,   299,
    1,   151
  ), ncol = 2, byrow = TRUE)))
  pred_sf <- sf::st_sf(
    id = 1L,
    geometry = sf::st_sfc(pred_poly, crs = 3035)
  )
  pred_path <- file.path(dir, "pred.shp")
  suppressWarnings(sf::st_write(pred_sf, pred_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  list(
    burnable = burnable_path,
    mask     = mask_path,
    ref      = ref_path,
    pred     = pred_path
  )
}

run_validate_synth <- function(dir, threshold_min_detected = 10,
                                year_target = 2000L) {
  fx <- build_synth_validation_fixture(dir)
  suppressWarnings(suppressMessages(
    validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      year_target            = year_target,
      validation_dir         = dir,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = threshold_min_detected,
      metrics_type           = "all"
    )
  ))
}

test_that("pixel_global has Specificity, BalancedAccuracy, ErrorRate with finite values", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_pixglob_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  res <- run_validate_synth(tmp)
  m <- res$metrics
  expect_true(!is.null(m) && nrow(m) == 1L)

  expect_true(all(c("Specificity", "BalancedAccuracy", "ErrorRate") %in%
                    names(m)))
  expect_true(is.finite(m$Specificity))
  expect_true(is.finite(m$BalancedAccuracy))
  expect_true(is.finite(m$ErrorRate))

  # Hand-checked values for the synthetic fixture:
  # TP=30, FP=0, FN=20, TN=50  ->  Recall=0.6, Specificity=1.0,
  # BalancedAccuracy=0.8, ErrorRate=0.20.
  expect_equal(unname(m$TP), 30)
  expect_equal(unname(m$FP), 0)
  expect_equal(unname(m$FN), 20)
  expect_equal(unname(m$TN), 50)
  expect_equal(unname(m$Specificity),     1.0, tolerance = 1e-9)
  expect_equal(unname(m$BalancedAccuracy), 0.8, tolerance = 1e-9)
  expect_equal(unname(m$ErrorRate),        0.2, tolerance = 1e-9)
})

test_that("polygon_global has Coverage_mean/median/p10/p90 and Detected_Definition", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_polyglob_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  res <- run_validate_synth(tmp)
  p <- res$polygon_summary
  expect_true(!is.null(p) && nrow(p) == 1L)

  expect_true(all(c("Coverage_mean", "Coverage_median",
                    "Coverage_p10",  "Coverage_p90",
                    "Detected_Definition") %in% names(p)))
  expect_true(is.finite(p$Coverage_mean))
  expect_true(is.finite(p$Coverage_median))
  expect_true(is.finite(p$Coverage_p10))
  expect_true(is.finite(p$Coverage_p90))

  # Coverages: ref1 = 100, ref2 = 20  ->  mean = 60, median = 60
  expect_equal(unname(p$Coverage_mean),   60, tolerance = 1e-9)
  expect_equal(unname(p$Coverage_median), 60, tolerance = 1e-9)
  expect_equal(unname(p$Detected_Definition), "coverage_ref >= 10%")
})

test_that("threshold_min_detected = 50 yields fewer N_Detected_Polygons than default 10", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp_default <- tempfile("of_vfm_thr10_")
  tmp_strict  <- tempfile("of_vfm_thr50_")
  on.exit({
    unlink(tmp_default, recursive = TRUE)
    unlink(tmp_strict,  recursive = TRUE)
  }, add = TRUE)

  res_default <- run_validate_synth(tmp_default, threshold_min_detected = 10)
  res_strict  <- run_validate_synth(tmp_strict,  threshold_min_detected = 50)

  n_default <- res_default$polygon_summary$N_Detected_Polygons
  n_strict  <- res_strict$polygon_summary$N_Detected_Polygons

  # ref1 = 100% (kept under both), ref2 = 20% (kept at 10, dropped at 50).
  expect_equal(unname(n_default), 2L)
  expect_equal(unname(n_strict),  1L)
  expect_lt(n_strict, n_default)

  expect_equal(unname(res_default$polygon_summary$Detected_Definition),
               "coverage_ref >= 10%")
  expect_equal(unname(res_strict$polygon_summary$Detected_Definition),
               "coverage_ref >= 50%")
})
