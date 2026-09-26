# Tests for the UNDETERMINED_OBS_DATE policy (2.3.0).
#
# Missing authoritative observability metadata is NOT evidence of temporal
# non-observability, so it must not be treated like NOT_OBSERVABLE. These tests
# pin that separation: `observable_flag` is evidence, `in_evaluation_domain` is
# domain membership, and only the latter decides who is evaluated.
#
# The fixture builder `build_obs_fixture()` and `run_obs()` live in
# test-observability-wholefire.R (same testthat directory, sourced together).

# ref2 is fully observed and carries a perfectly usable end_doy, but its
# authoritative observability date is deliberately NA. The prediction covers
# BOTH quadrants, so ref2 can produce TP.
build_undet_fixture <- function(dir, pred_over_ref2 = TRUE) {
  build_obs_fixture(dir, ref2_obs_rows = 5L, ref2_end_doy = 100,
                    extra_ref_cols = list(obs_required_doy = c(100, NA_real_)),
                    pred_over_ref2 = pred_over_ref2)
}
run_undet <- function(fx, dir, policy) {
  suppressWarnings(run_obs(
    fx, dir, observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75,
    ref_obs_doy_col = "obs_required_doy",
    observability_undetermined_policy = policy))
}

test_that("T11 evaluate keeps UNDETERMINED as a positive reference feature", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_undet_fixture(file.path(tmp, "fx"))

  # the default IS evaluate
  expect_identical(
    eval(formals(validate_fire_maps)$observability_undetermined_policy)[1],
    "evaluate")

  res <- run_undet(fx, file.path(tmp, "ev"), "evaluate")
  ro <- res$reference_observability

  # status and flags: evidence FALSE, domain TRUE, and no date invented
  expect_identical(ro$observability_status, c("OBSERVABLE", "UNDETERMINED_OBS_DATE"))
  expect_identical(ro$observable_flag, c(TRUE, FALSE))
  expect_identical(ro$in_evaluation_domain, c(TRUE, TRUE))
  expect_true(is.na(ro$fire_required_doy[2]))
  expect_true(is.na(ro$required_doy_source[2]))
  expect_identical(ro$undetermined_cause,
                   c(NA_character_, "na_in_authoritative_obs_required_doy"))

  # stays a POSITIVE reference feature: both quadrants are reference-burned and
  # the prediction covers both, so all 50 cells are TP and nothing is removed.
  expect_equal(res$metrics$TP, 50)
  expect_equal(res$metrics$FP, 0)
  expect_equal(res$metrics$FN, 0)
  expect_equal(res$metrics$TN, 50)
  expect_equal(res$observability_settings$domain_removed_ha, 0)
  expect_equal(res$metrics$Evaluable_Area_ha, 100 * 0.09)
  expect_equal(res$polygon_summary$N_Reference_Polygons, 2L)

  # and it can produce FN too, when the prediction does NOT cover it
  fx2 <- build_undet_fixture(file.path(tmp, "fx2"), pred_over_ref2 = FALSE)
  res2 <- run_undet(fx2, file.path(tmp, "ev2"), "evaluate")
  expect_equal(res2$metrics$TP, 25)
  expect_equal(res2$metrics$FN, 25)   # the undetermined fire is omitted
  expect_equal(res2$metrics$FP, 0)    # and never becomes background
  expect_equal(res2$observability_settings$domain_removed_ha, 0)
})

test_that("T11b exclude removes UNDETERMINED from reference and domain symmetrically", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_undet_fixture(file.path(tmp, "fx"))
  res <- run_undet(fx, file.path(tmp, "ex"), "exclude")
  ro <- res$reference_observability

  expect_identical(ro$observability_status, c("OBSERVABLE", "UNDETERMINED_OBS_DATE"))
  expect_identical(ro$observable_flag, c(TRUE, FALSE))
  expect_identical(ro$in_evaluation_domain, c(TRUE, FALSE))

  # ref2 leaves the domain entirely: no TP, no FP, no FN, no TN there
  expect_equal(res$metrics$TP, 25)
  expect_equal(res$metrics$FP, 0)
  expect_equal(res$metrics$FN, 0)
  expect_equal(res$metrics$TN, 50)
  expect_equal(res$observability_settings$domain_removed_ha, 25 * 0.09)
  expect_equal(res$metrics$Evaluable_Area_ha, 75 * 0.09)
  expect_equal(res$polygon_summary$N_Reference_Polygons, 1L)
})

test_that("T11c the two policies differ ONLY in domain membership", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_undet_fixture(file.path(tmp, "fx"))
  ev <- run_undet(fx, file.path(tmp, "ev"), "evaluate")
  ex <- run_undet(fx, file.path(tmp, "ex"), "exclude")

  # accounting closure: everything that leaves a cell of the matrix leaves the
  # domain, exactly
  d <- (ex$metrics$TP - ev$metrics$TP) + (ex$metrics$FP - ev$metrics$FP) +
       (ex$metrics$FN - ev$metrics$FN) + (ex$metrics$TN - ev$metrics$TN)
  expect_equal(-d * 0.09, ex$observability_settings$domain_removed_ha)

  # the verdicts themselves are identical; only membership moved
  expect_identical(ev$reference_observability$observability_status,
                   ex$reference_observability$observability_status)
  expect_identical(ev$reference_observability$observable_flag,
                   ex$reference_observability$observable_flag)
  expect_false(identical(ev$reference_observability$in_evaluation_domain,
                         ex$reference_observability$in_evaluation_domain))
})

test_that("T12 no required DOY means NA fractions, never 0", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_undet_fixture(file.path(tmp, "fx"))
  ro <- run_undet(fx, file.path(tmp, "na"), "evaluate")$reference_observability

  # a 0 here would read as "observed nowhere", i.e. as EVIDENCE of
  # non-observability, which is precisely what an undetermined date is not
  expect_true(is.na(ro$temporal_observable_fraction[2]))
  expect_true(is.na(ro$post_fire_observable_fraction[2]))
  expect_true(is.na(ro$pre_fire_or_too_early_fraction[2]))
  expect_true(is.na(ro$n_temporally_observable[2]))
  # data coverage does not depend on the fire date and must stay valid
  expect_equal(ro$data_coverage_fraction[2], 1)
  expect_false(is.na(ro$obs_doy_max[2]))
  # the fire WITH a date keeps real numbers
  expect_equal(ro$temporal_observable_fraction[1], 1)
})

test_that("T13 OUTSIDE_BURNABLE_DOMAIN is separated from NO_OBSERVABILITY_DATA", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()

  # (a) burnable is 1 only on the top-left quadrant, so ref2 has no domain cell
  dir <- file.path(tmp, "fxout")
  fx <- build_obs_fixture(dir, ref2_obs_rows = 5L)
  bur <- terra::rast(fx$burnable)
  m <- matrix(0, nrow = 10, ncol = 10); m[1:5, 1:5] <- 1
  terra::values(bur) <- as.vector(t(m))
  terra::writeRaster(bur, fx$burnable, overwrite = TRUE, datatype = "FLT4S")
  ro <- suppressWarnings(run_obs(
    fx, file.path(tmp, "out"), observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75))$reference_observability
  expect_identical(ro$observability_status,
                   c("OBSERVABLE", "OUTSIDE_BURNABLE_DOMAIN"))
  expect_equal(ro$total_pixels[2], 0L)
  expect_false(ro$in_evaluation_domain[2])

  # (b) a TRUE NO_OBSERVABILITY_DATA: burnable surface exists and the fire has a
  # date, but the composite carries no observation anywhere inside it
  dir2 <- file.path(tmp, "fxnod")
  fx2 <- build_obs_fixture(dir2, ref2_obs_rows = 5L)
  obs <- terra::rast(fx2$obs)
  mo <- matrix(NA_real_, nrow = 10, ncol = 10); mo[1:5, 1:5] <- 120
  terra::values(obs) <- as.vector(t(mo))
  terra::writeRaster(obs, fx2$obs, overwrite = TRUE, datatype = "FLT4S")
  ro2 <- suppressWarnings(run_obs(
    fx2, file.path(tmp, "nod"), observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75))$reference_observability
  expect_identical(ro2$observability_status,
                   c("OBSERVABLE", "NO_OBSERVABILITY_DATA"))
  expect_gt(ro2$total_pixels[2], 0)          # it DOES have evaluable surface
  expect_equal(ro2$data_coverage_fraction[2], 0)
  expect_false(ro2$in_evaluation_domain[2])  # evidence of non-observation -> out
})

test_that("T14 NOT_OBSERVABLE still leaves the domain at TOF below the threshold", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # TOF = 3/5 = 0.60 < 0.75, with a perfectly valid authoritative date
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 3L,
                          extra_ref_cols = list(obs_required_doy = c(100, 100)),
                          pred_over_ref2 = TRUE)
  res <- suppressWarnings(run_obs(
    fx, file.path(tmp, "no"), observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75, ref_obs_doy_col = "obs_required_doy",
    observability_undetermined_policy = "evaluate"))
  ro <- res$reference_observability
  expect_identical(ro$observability_status, c("OBSERVABLE", "NOT_OBSERVABLE"))
  expect_equal(ro$temporal_observable_fraction[2], 0.6)
  expect_false(ro$in_evaluation_domain[2])
  # the evaluate policy must NOT rescue it: this is evidence, not missing metadata
  expect_equal(res$observability_settings$domain_removed_ha, 25 * 0.09)
  expect_equal(res$metrics$FP, 0)
})

test_that("T15 the undetermined policy moves the cache tag", {
  obs_raster <- withr::local_tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0,
                   ymax = 120, crs = "EPSG:3035")
  terra::values(r) <- 1:16; names(r) <- "doy"
  terra::writeRaster(r, obs_raster, overwrite = TRUE)
  fp <- function(...) .vfm_observability_fingerprint(
    observability_raster = obs_raster, ref_end_doy_col = "end_doy",
    ref_start_doy_col = "start_doy", ref_obs_doy_col = "obs_required_doy", ...)

  ev <- fp(observability_mode = "wholefire_fraction",
           observability_undetermined_policy = "evaluate")
  ex <- fp(observability_mode = "wholefire_fraction",
           observability_undetermined_policy = "exclude")
  expect_false(identical(ev, ex))
  # inert under legacy_any, which has no such policy
  expect_identical(
    fp(observability_mode = "legacy_any",
       observability_undetermined_policy = "evaluate"),
    fp(observability_mode = "legacy_any",
       observability_undetermined_policy = "exclude"))
  expect_identical(.VFM_OBS_RULE_VERSION,
                   "obsrule-3:wholefire{legacy_any|fraction}+undetpolicy")
  expect_identical(.VFM_REF_CACHE_SCHEMA, "refschema-5")
})

test_that("T16 the per-fire detection audit is written and filterable", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_undet_fixture(file.path(tmp, "fx"))
  d <- file.path(tmp, "aud")
  run_undet(fx, d, "evaluate")
  obs_dir <- file.path(d, "VALIDATION", "02_OBSERVABILITY")
  det <- list.files(obs_dir, pattern = "^REFERENCE_FIRES_DETECTION_2000_.*[.]csv$",
                    full.names = TRUE)
  und <- list.files(obs_dir, pattern = "^UNDETERMINED_FIRES_DETECTION_2000_.*[.]csv$",
                    full.names = TRUE)
  expect_length(det, 1L); expect_length(und, 1L)
  x <- utils::read.csv(det[1]); u <- utils::read.csv(und[1])
  expect_equal(nrow(x), 2L)
  expect_equal(nrow(u), 1L)
  expect_identical(unique(u$observability_status), "UNDETERMINED_OBS_DATE")

  # everything needed to revisit a fire without recomputing the validation
  for (nm in c("observability_status", "undetermined_cause", "obs_required_doy",
               "fire_required_doy", "observable_flag", "in_evaluation_domain",
               "data_coverage_fraction", "temporal_observable_fraction",
               "reference_area_ha", "detected_area_ha", "omitted_area_ha",
               "TP_cells", "TP_area_ha", "FN_cells", "FN_area_ha",
               "fire_recall_percent")) {
    expect_true(nm %in% names(x), info = nm)
  }
  # the undetermined fire is fully covered by the prediction here
  expect_equal(u$TP_cells, 25)
  expect_equal(u$FN_cells, 0)
  expect_equal(u$fire_recall_percent, 100)
  expect_true(is.na(u$temporal_observable_fraction))
  expect_true(is.na(u$fire_required_doy))
})
