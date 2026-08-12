# Frozen temporal-observability semantics + partial-observability audit (0.10.1).
#
# Contract (frozen): observability in validate_fire_maps() is TEMPORAL and
# WHOLE-FIRE. A reference fire with no valid observable DOY, or whose fire date
# is after the last observation, is excluded from the reference BEFORE
# rasterization, so it contributes 0 TP, 0 FN and 0 evaluated area and never
# enters omission / recall / F1 / IoU / per-fire detection. Partial spatial
# (pixel-level) observability inside KEPT fires is REPORTED for audit but never
# modifies the validation contract.

oa_ns <- asNamespace("OtsuFire")

# --- unit: the audit summarizer ---------------------------------------------

test_that(".vfm_observability_audit summarizes excluded + partial fires", {
  # 4 original fires: 2 excluded (one missing obs, one obs-before-fire), 1 fully
  # observable, 1 partially observable (kept whole).
  ro <- data.frame(
    id                    = c("F1", "F2", "F3", "F4"),
    reference_row         = 1:4,
    ref_area_domain_ha    = c(10, 20, 30, 40),
    observable_flag       = c(FALSE, FALSE, TRUE, TRUE),
    observable_reason     = c("no_observability_data", "obs_before_reference_doy",
                              "observable", "observable"),
    total_pixels          = c(100L, 200L, 300L, 400L),
    observable_pixels     = c(0L, 0L, 300L, 360L),
    non_observable_pixels = c(100L, 200L, 0L, 40L),
    observable_fraction   = c(0, 0, 1, 0.9),
    stringsAsFactors      = FALSE
  )
  au <- oa_ns$.vfm_observability_audit(ro, cell_area_ha = 0.81)

  s <- au$summary
  expect_equal(s$n_reference_fires_original, 4L)
  expect_equal(s$n_reference_fires_observable, 2L)
  expect_equal(s$n_reference_fires_excluded, 2L)
  expect_equal(s$n_reference_fires_partially_observable, 1L)   # only F4
  expect_equal(s$excluded_reference_area_ha, 30)               # 10 + 20
  expect_equal(s$partial_non_observable_pixels, 40L)           # F4 only (F3 = 0)
  expect_equal(s$partial_non_observable_area_ha, 40 * 0.81)

  expect_equal(nrow(au$excluded), 2L)
  expect_setequal(au$excluded$fire_id, c("F1", "F2"))
  expect_true(all(grepl("0 TP, 0 FN", au$excluded$confirmation_not_in_metrics)))

  expect_equal(nrow(au$partial), 1L)
  expect_equal(au$partial$fire_id, "F4")
  expect_equal(au$partial$temporally_valid_pixels, 360L)
  expect_equal(au$partial$non_observable_pixels, 40L)
  expect_equal(au$partial$observable_fraction, 0.9)
  expect_true(au$partial$contract_keeps_whole_fire)
})

test_that(".vfm_observability_audit handles the all-observable / no-partial case", {
  ro <- data.frame(
    id = c("A", "B"), reference_row = 1:2,
    ref_area_domain_ha = c(5, 6),
    observable_flag = c(TRUE, TRUE),
    observable_reason = c("observable", "observable"),
    total_pixels = c(50L, 60L), observable_pixels = c(50L, 60L),
    non_observable_pixels = c(0L, 0L), observable_fraction = c(1, 1),
    stringsAsFactors = FALSE
  )
  au <- oa_ns$.vfm_observability_audit(ro, cell_area_ha = 0.81)
  expect_equal(au$summary$n_reference_fires_excluded, 0L)
  expect_equal(au$summary$n_reference_fires_partially_observable, 0L)
  expect_equal(nrow(au$excluded), 0L)
  expect_equal(nrow(au$partial), 0L)
})

# --- end-to-end on the real validator ---------------------------------------

oa_synth <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  burnable <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 300,
                          ymin = 0, ymax = 300, crs = "EPSG:3035")
  terra::values(burnable) <- 1
  bp <- file.path(dir, "burnable.tif")
  terra::writeRaster(burnable, bp, overwrite = TRUE, datatype = "FLT4S")

  mask_poly <- sf::st_polygon(list(matrix(
    c(0, 0, 300, 0, 300, 300, 0, 300, 0, 0), ncol = 2, byrow = TRUE)))
  mp <- file.path(dir, "mask.shp")
  suppressWarnings(sf::st_write(
    sf::st_sf(id = 1L, geometry = sf::st_sfc(mask_poly, crs = 3035)),
    mp, quiet = TRUE, delete_dsn = TRUE))

  ref1 <- sf::st_polygon(list(matrix(
    c(1, 151, 149, 151, 149, 299, 1, 299, 1, 151), ncol = 2, byrow = TRUE)))
  ref2 <- sf::st_polygon(list(matrix(
    c(151, 151, 299, 151, 299, 299, 151, 299, 151, 151), ncol = 2, byrow = TRUE)))
  rp <- file.path(dir, "ref.shp")
  suppressWarnings(sf::st_write(sf::st_sf(
    id = 1:2, year = c(2000L, 2000L),
    end_doy = c(100, NA), start_doy = c(90, 150),
    geometry = sf::st_sfc(ref1, ref2, crs = 3035)),
    rp, quiet = TRUE, delete_dsn = TRUE))

  # ref1 observable (obs 120 >= ref 100); ref2 non-observable (obs 140 < ref 150)
  obs <- burnable; terra::values(obs) <- NA_real_
  obs <- terra::rasterize(
    terra::vect(sf::st_sf(obs_doy = c(120, 140),
                          geometry = sf::st_sfc(ref1, ref2, crs = 3035))),
    obs, field = "obs_doy", background = NA)
  op <- file.path(dir, "obs.tif")
  terra::writeRaster(obs, op, overwrite = TRUE, datatype = "FLT4S")

  pred <- sf::st_polygon(list(matrix(
    c(1, 151, 179, 151, 179, 299, 1, 299, 1, 151), ncol = 2, byrow = TRUE)))
  pp <- file.path(dir, "pred.shp")
  suppressWarnings(sf::st_write(
    sf::st_sf(id = 1L, geometry = sf::st_sfc(pred, crs = 3035)),
    pp, quiet = TRUE, delete_dsn = TRUE))
  list(burnable = bp, mask = mp, ref = rp, obs = op, pred = pp)
}

oa_run <- function(fx, dir, obs = TRUE) {
  suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile = fx$pred, ref_shapefile = fx$ref, mask_shapefile = fx$mask,
    burnable_raster = fx$burnable,
    observability_raster = if (obs) fx$obs else NULL,
    year_target = 2000L, validation_dir = dir,
    binary_burnable = TRUE, burnable_threshold = 0.5,
    threshold_completely_detected = 90, threshold_min_detected = 10,
    metrics_type = "all")))
}

test_that("a temporally non-observable fire contributes 0 TP / 0 FN / 0 area and reconciles", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d_obs <- tempfile("oa_obs_"); d_no <- tempfile("oa_no_")
  on.exit(unlink(c(d_obs, d_no), recursive = TRUE, force = TRUE), add = TRUE)

  res_obs <- oa_run(oa_synth(d_obs), d_obs, obs = TRUE)
  res_no  <- oa_run(oa_synth(d_no),  d_no,  obs = FALSE)

  # Reconciliation: original - excluded = observable (the 677-14=663 analogue).
  au <- res_obs$observability_audit
  expect_false(is.null(au))
  expect_equal(au$summary$n_reference_fires_original, 2L)
  expect_equal(au$summary$n_reference_fires_excluded, 1L)
  expect_equal(au$summary$n_reference_fires_observable, 1L)
  expect_equal(au$summary$n_reference_fires_original -
                 au$summary$n_reference_fires_excluded,
               au$summary$n_reference_fires_observable)

  # The excluded fire (ref2) is the sole source of FN without observability.
  # With the temporal filter it is removed BEFORE rasterization: its 20 FN
  # pixels and its reference area leave the confusion matrix entirely.
  expect_equal(unname(res_no$metrics$FN), 20)    # ref2 omitted -> 20 FN
  expect_equal(unname(res_obs$metrics$FN), 0)    # ref2 excluded -> 0 FN
  expect_equal(unname(res_obs$metrics$TP), 25)   # only ref1 remains
  expect_equal(unname(res_obs$metrics$FP), 5)

  # The excluded fire is absent from per-fire detection / omission accounting.
  expect_equal(unname(res_obs$polygon_summary$N_Reference_Polygons), 1L)
  expect_equal(unname(res_obs$polygon_summary$N_Not_Detected), 0L)

  # Excluded area is reported and equals ref2's burnable area (not in metrics).
  expect_gt(au$summary$excluded_reference_area_ha, 0)
  excl_area_from_table <-
    res_obs$reference_observability$ref_area_domain_ha[
      !res_obs$reference_observability$observable_flag]
  expect_equal(au$summary$excluded_reference_area_ha, round(excl_area_from_table, 4))
})

test_that("the three audit CSVs are written and never change the metrics", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- tempfile("oa_csv_"); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  res <- oa_run(oa_synth(d), d, obs = TRUE)

  obs_dir <- file.path(d, "VALIDATION", "02_OBSERVABILITY")
  expect_true(file.exists(file.path(obs_dir, "OBSERVABILITY_AUDIT_2000.csv")))
  expect_true(file.exists(file.path(obs_dir, "OBSERVABILITY_EXCLUDED_FIRES_2000.csv")))
  expect_true(file.exists(file.path(obs_dir, "OBSERVABILITY_PARTIAL_FIRES_2000.csv")))

  audit_csv <- read.csv(file.path(obs_dir, "OBSERVABILITY_AUDIT_2000.csv"))
  expect_equal(audit_csv$n_reference_fires_excluded, 1L)
  excl_csv <- read.csv(file.path(obs_dir, "OBSERVABILITY_EXCLUDED_FIRES_2000.csv"))
  expect_equal(nrow(excl_csv), 1L)
  expect_match(excl_csv$exclusion_reason, "obs_before_reference_doy")

  # Metrics are exactly the frozen hand-checked values (audit is side-effect free
  # w.r.t. the confusion matrix).
  expect_equal(unname(res$metrics$TP), 25)
  expect_equal(unname(res$metrics$FP), 5)
  expect_equal(unname(res$metrics$FN), 0)
  expect_equal(unname(res$metrics$TN), 70)
})
