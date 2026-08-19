# Tests for the whole-fire observability v2 rule.
#
# Fixture and runner: helper-observability-fixture.R (same directory).
#
# Rules under test
#   legacy_any         : observable <=> max(DOY) >= required DOY
#   wholefire_fraction : observable <=> TOF >= observability_min_fraction,
#                        plus removal of the excluded geometries from the
#                        evaluation domain.

# ---- T1: legacy regression ---------------------------------------------------

test_that("T1 legacy_any reproduces the historical rule exactly", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()

  # (a) default mode is the legacy one
  expect_identical(
    eval(formals(validate_fire_maps)$observability_mode)[1],
    "legacy_any"
  )

  # (b) ref2 observed on 1 of its 5 rows -> TOF = 0.2, but max(DOY) = 200 >= 100,
  #     so the legacy rule must KEEP it (this is the "one pixel rescues the
  #     fire" behaviour the legacy mode is defined by).
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 1L)
  res <- run_obs(fx, file.path(tmp, "legacy"), observability_mode = "legacy_any")

  ro <- res$reference_observability
  expect_equal(nrow(ro), 2L)
  expect_true(all(ro$observable_flag))
  expect_identical(ro$observable_reason, c("observable", "observable"))
  expect_identical(ro$observability_status, c("OBSERVABLE", "OBSERVABLE"))
  expect_equal(ro$obs_doy_max, c(120, 200))
  # legacy required DOY still uses the end -> start cascade
  expect_identical(ro$required_doy_source, c("end_doy", "end_doy"))
  # the diagnostic fraction is computed even in legacy mode and is NOT used
  expect_equal(ro$temporal_observable_fraction, c(1, 0.2))
  # legacy does NOT touch the evaluation domain
  expect_false(res$observability_settings$domain_masked)

  # (c) hand-computed confusion matrix, identical to pre-0.11.1 behaviour:
  #     pred covers ref1 entirely (25 TP), ref2 untouched (25 FN),
  #     bottom half unburned and unpredicted (50 TN), no FP.
  m <- res$metrics
  expect_equal(m$TP, 25); expect_equal(m$FP, 0)
  expect_equal(m$FN, 25); expect_equal(m$TN, 50)
  expect_equal(m$Precision, 1)
  expect_equal(m$Recall, 0.5)
  expect_equal(m$F1, 2 * 1 * 0.5 / 1.5)
  expect_equal(m$IoU, 0.5)

  # (d) the legacy verdict is exposed and agrees with observable_flag
  expect_identical(ro$observable_flag_legacy_any, ro$observable_flag)
})

test_that("T1b legacy_any and the ANY rule are mathematically identical", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # sweep every possible number of observable rows in ref2, including 0
  for (k in 0:5) {
    fx <- build_obs_fixture(file.path(tmp, paste0("fx", k)), ref2_obs_rows = k)
    res <- run_obs(fx, file.path(tmp, paste0("run", k)),
                   observability_mode = "legacy_any")
    ro <- res$reference_observability
    any_rule <- ro$temporal_observable_fraction > 0
    expect_identical(ro$observable_flag, any_rule,
                     info = sprintf("ref2_obs_rows = %d", k))
  }
})

# ---- T2/T3/T4: the 0.75 threshold -------------------------------------------

test_that("T2/T3/T4 the fraction rule cuts exactly at observability_min_fraction", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()

  # ref2 rows -> TOF: 3/5 = 0.60 (below), 4/5 = 0.80 (above), 5/5 = 1.
  # 0.74 / 0.75 are exercised by moving the THRESHOLD onto the observed TOF,
  # which is the same comparison with the same operator.
  fx4 <- build_obs_fixture(file.path(tmp, "fx4"), ref2_obs_rows = 4L)

  # T2: TOF = 0.80 but threshold 0.81 -> NOT OBSERVABLE (strictly below)
  r_below <- run_obs(fx4, file.path(tmp, "below"),
                     observability_mode = "wholefire_fraction",
                     observability_min_fraction = 0.81)
  ro <- r_below$reference_observability
  expect_equal(ro$temporal_observable_fraction, c(1, 0.8))
  expect_identical(ro$observability_status, c("OBSERVABLE", "NOT_OBSERVABLE"))
  expect_identical(ro$observable_flag, c(TRUE, FALSE))

  # T3: TOF = 0.80 with threshold exactly 0.80 -> OBSERVABLE (>= is inclusive)
  r_eq <- run_obs(fx4, file.path(tmp, "eq"),
                  observability_mode = "wholefire_fraction",
                  observability_min_fraction = 0.80)
  expect_identical(r_eq$reference_observability$observability_status,
                   c("OBSERVABLE", "OBSERVABLE"))

  # T2 again with the production threshold: TOF = 0.60 < 0.75 -> excluded
  fx3 <- build_obs_fixture(file.path(tmp, "fx3"), ref2_obs_rows = 3L)
  r075 <- run_obs(fx3, file.path(tmp, "t075"),
                  observability_mode = "wholefire_fraction",
                  observability_min_fraction = 0.75)
  ro3 <- r075$reference_observability
  expect_equal(ro3$temporal_observable_fraction, c(1, 0.6))
  expect_identical(ro3$observability_status, c("OBSERVABLE", "NOT_OBSERVABLE"))
  # and 0.80 >= 0.75 -> kept
  r075b <- run_obs(fx4, file.path(tmp, "t075b"),
                   observability_mode = "wholefire_fraction",
                   observability_min_fraction = 0.75)
  expect_identical(r075b$reference_observability$observability_status,
                   c("OBSERVABLE", "OBSERVABLE"))

  # T4: TOF = 1 -> OBSERVABLE
  fx5 <- build_obs_fixture(file.path(tmp, "fx5"), ref2_obs_rows = 5L)
  r_full <- run_obs(fx5, file.path(tmp, "full"),
                    observability_mode = "wholefire_fraction",
                    observability_min_fraction = 0.75)
  rof <- r_full$reference_observability
  expect_equal(rof$temporal_observable_fraction, c(1, 1))
  expect_true(all(rof$observable_flag))
  expect_equal(r_full$observability_settings$domain_removed_ha, 0)
})

# ---- T5/T6: effect on FP -----------------------------------------------------

test_that("T5 a non-observable fire covered by the prediction produces no FP", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # ref2 is only 20% observable and the prediction covers BOTH quadrants.
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 1L,
                          pred_over_ref2 = TRUE)

  legacy <- run_obs(fx, file.path(tmp, "legacy"),
                    observability_mode = "legacy_any")
  # under legacy ref2 is kept, so the prediction over it is TP, not FP
  expect_equal(legacy$metrics$TP, 50)
  expect_equal(legacy$metrics$FP, 0)

  # Now make ref2 non-observable under BOTH rules by removing every post-fire
  # pixel from it, and keep the prediction on top of it.
  fx0 <- build_obs_fixture(file.path(tmp, "fx0"), ref2_obs_rows = 0L,
                           pred_over_ref2 = TRUE)

  leg0 <- run_obs(fx0, file.path(tmp, "leg0"), observability_mode = "legacy_any")
  # LEGACY: ref2 is dropped from the reference but its territory stays in the
  # domain as reference = 0, so the 25 predicted cells become FALSE POSITIVES.
  expect_equal(leg0$metrics$TP, 25)
  expect_equal(leg0$metrics$FP, 25)
  expect_equal(leg0$metrics$FN, 0)
  expect_equal(leg0$metrics$TN, 50)

  new0 <- run_obs(fx0, file.path(tmp, "new0"),
                  observability_mode = "wholefire_fraction",
                  observability_min_fraction = 0.75)
  # NEW: ref2's geometry leaves the evaluation domain entirely -> no FP, no FN,
  # no TN there. Only ref1 (25 TP) and the unburned bottom half (50 TN) remain.
  expect_equal(new0$metrics$TP, 25)
  expect_equal(new0$metrics$FP, 0)
  expect_equal(new0$metrics$FN, 0)
  expect_equal(new0$metrics$TN, 50)
  expect_equal(new0$metrics$Precision, 1)
  expect_true(new0$observability_settings$domain_masked)
  # 25 cells of 30 m = 25 * 0.09 ha
  expect_equal(new0$observability_settings$domain_removed_ha, 25 * 0.09)
  expect_equal(new0$metrics$Evaluable_Area_ha, 75 * 0.09)
})

test_that("T5b the removed domain follows the cell-centre rule, not `touches`", {
  skip_if_not_installed("terra")
  # terra::mask(SpatRaster, SpatVector) defaults to touches = TRUE, which would
  # remove every cell the excluded perimeter merely grazes. The validator must
  # force touches = FALSE so the removed set matches the reference
  # rasterization exactly.
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100,
                   ymin = 0, ymax = 100, crs = "EPSG:3035")
  terra::values(r) <- 1
  p <- terra::vect(cbind(id = 1, x = c(48, 58, 58, 48, 48),
                         y = c(48, 48, 58, 58, 48)),
                   type = "polygons", crs = "EPSG:3035")
  removed <- function(...) terra::global(
    is.na(terra::mask(r, p, inverse = TRUE, updatevalue = NA, ...)),
    "sum", na.rm = TRUE)[[1]]
  # the default over-excludes; the cell-centre rule agrees with rasterize()
  expect_gt(removed(), removed(touches = FALSE))
  expect_equal(
    removed(touches = FALSE),
    terra::global(!is.na(terra::rasterize(p, r)), "sum", na.rm = TRUE)[[1]]
  )

  # end to end: the domain removed for a non-observable fire equals that fire's
  # burnable-domain area exactly (ref_area_domain_ha), with no dilation.
  tmp <- withr::local_tempdir()
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 0L)
  res <- run_obs(fx, file.path(tmp, "run"),
                 observability_mode = "wholefire_fraction",
                 observability_min_fraction = 0.75)
  ro <- res$reference_observability
  excl_area <- sum(ro$ref_area_domain_ha[!ro$observable_flag])
  expect_equal(res$observability_settings$domain_removed_ha, excl_area)
})

test_that("T6 a false positive outside every reference fire still counts", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 0L)

  res <- suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile      = fx$pred_outside,   # bottom-left quadrant, no ref there
    ref_shapefile        = fx$ref,
    mask_shapefile       = fx$mask,
    burnable_raster      = fx$burnable,
    year_target          = 2000L,
    validation_dir       = file.path(tmp, "out"),
    observability_raster = fx$obs,
    metrics_type         = "all",
    observability_mode   = "wholefire_fraction",
    observability_min_fraction = 0.75,
    force_reprocess_ref  = TRUE,
    force_reprocess_pred = TRUE
  )))
  # ref2 excluded (25 cells out of the domain). Remaining domain = 75 cells:
  # ref1 (25, unpredicted -> FN) + bottom half (50, of which 25 predicted).
  expect_equal(res$metrics$FP, 25)   # commission outside reference is preserved
  expect_equal(res$metrics$TP, 0)
  expect_equal(res$metrics$FN, 25)
  expect_equal(res$metrics$TN, 25)
})

# ---- T7: unknown end date ----------------------------------------------------

test_that("T7 end_doy = NA yields UNDETERMINED_OBS_DATE with no silent fallback", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # ref2 has no end date but a perfectly usable start date, and its pixels are
  # all observed at DOY 200 (i.e. after start_doy = 90). The legacy rule would
  # silently fall back to start_doy and call it observable.
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 5L,
                          ref2_end_doy = NA, ref2_start_doy = 90)

  legacy <- run_obs(fx, file.path(tmp, "legacy"), observability_mode = "legacy_any")
  rol <- legacy$reference_observability
  expect_identical(rol$required_doy_source, c("end_doy", "start_doy"))
  expect_equal(rol$fire_required_doy, c(100, 90))
  expect_true(all(rol$observable_flag))

  new <- run_obs(fx, file.path(tmp, "new"),
                 observability_mode = "wholefire_fraction",
                 observability_min_fraction = 0.75)
  ron <- new$reference_observability
  expect_identical(ron$observability_status,
                   c("OBSERVABLE", "UNDETERMINED_OBS_DATE"))
  expect_true(is.na(ron$fire_required_doy[2]))
  expect_true(is.na(ron$required_doy_source[2]))
  expect_false(ron$observable_flag[2])
  # 2.3.0: under the default "evaluate" policy it is NOT excluded. Missing
  # metadata is not evidence of non-observability, so it stays a positive
  # reference feature and the prediction (which does not cover ref2) omits it.
  expect_true(ron$in_evaluation_domain[2])
  expect_equal(new$observability_settings$domain_removed_ha, 0)
  expect_equal(new$metrics$FN, 25)
  expect_equal(new$metrics$Evaluable_Area_ha, 100 * 0.09)
  # the "exclude" policy is what reproduces the old exclusion, on request only
  excl <- suppressWarnings(run_obs(
    fx, file.path(tmp, "excl"), observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75,
    observability_undetermined_policy = "exclude"))
  expect_false(excl$reference_observability$in_evaluation_domain[2])
  expect_true(excl$observability_settings$domain_masked)
  expect_equal(excl$metrics$FN, 0)
  expect_equal(excl$metrics$Evaluable_Area_ha, 75 * 0.09)
  # reported separately, by source
  s <- new$observability_summary
  expect_equal(s$N_Undetermined_Obs_Date, 1L)
  expect_equal(s$N_Observable, 1L)
  und <- file.path(tmp, "new", "VALIDATION", "02_OBSERVABILITY",
                   "OBSERVABILITY_UNDETERMINED_2000.csv")
  expect_true(file.exists(und))
  expect_equal(nrow(utils::read.csv(und)), 1L)
})

# ---- T8: Neves uses DOY, not meanDOY ----------------------------------------

test_that("T8 ref_obs_doy_col overrides end_doy for observability only", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # ref2 pixels are observed at DOY 200 on 4 of 5 rows and at 50 on the rest.
  # end_doy = 100 -> TOF = 0.8 -> observable at 0.75.
  # DOY       = 250 -> nothing is post-fire -> TOF = 0 -> excluded.
  # Both rows carry a finite DOY so the override is exercised on both.
  fx <- build_obs_fixture(
    file.path(tmp, "fx"), ref2_obs_rows = 4L, ref2_end_doy = 100,
    extra_ref_cols = list(DOY = c(100, 250))
  )

  without <- run_obs(fx, file.path(tmp, "without"),
                     observability_mode = "wholefire_fraction",
                     observability_min_fraction = 0.75)
  expect_identical(without$reference_observability$required_doy_source,
                   c("end_doy", "end_doy"))
  expect_identical(without$reference_observability$observability_status,
                   c("OBSERVABLE", "OBSERVABLE"))

  with_doy <- run_obs(fx, file.path(tmp, "with"),
                      observability_mode = "wholefire_fraction",
                      observability_min_fraction = 0.75,
                      ref_obs_doy_col = "DOY")
  ro <- with_doy$reference_observability
  expect_identical(ro$required_doy_source, c("DOY", "DOY"))
  expect_equal(ro$fire_required_doy, c(100, 250))
  expect_equal(ro$temporal_observable_fraction, c(1, 0))
  expect_identical(ro$observability_status, c("OBSERVABLE", "NOT_OBSERVABLE"))
  expect_true(with_doy$observability_settings$ref_obs_doy_authoritative)

  # under legacy_any the column is ignored entirely (membership semantics
  # untouched: end_doy keeps driving the legacy rule)
  legacy <- run_obs(fx, file.path(tmp, "legacy"),
                    observability_mode = "legacy_any", ref_obs_doy_col = "DOY")
  expect_identical(legacy$reference_observability$required_doy_source,
                   c("end_doy", "end_doy"))
})

# ---- T10: the observability column is AUTHORITATIVE --------------------------

test_that("T10 an NA in ref_obs_doy_col means UNDETERMINED, with no silent fallback", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  # ref2 is fully observed (DOY 200 everywhere) and carries a perfectly usable
  # end_doy = 100. Its authoritative observability date is deliberately NA.
  # The authoritative rule must NOT reach for end_doy.
  fx <- build_obs_fixture(
    file.path(tmp, "fx"), ref2_obs_rows = 5L, ref2_end_doy = 100,
    extra_ref_cols = list(obs_required_doy = c(100, NA_real_))
  )

  auth <- suppressWarnings(run_obs(
    fx, file.path(tmp, "auth"),
    observability_mode = "wholefire_fraction",
    observability_min_fraction = 0.75,
    ref_obs_doy_col = "obs_required_doy"))
  ro <- auth$reference_observability
  expect_identical(ro$required_doy_source, c("obs_required_doy", NA_character_))
  expect_equal(ro$fire_required_doy, c(100, NA_real_))
  expect_identical(ro$observability_status,
                   c("OBSERVABLE", "UNDETERMINED_OBS_DATE"))
  expect_false(ro$observable_flag[2])
  # the cause is attributed to the authoritative column, not to a missing end
  expect_identical(ro$undetermined_cause,
                   c(NA_character_, "na_in_authoritative_obs_required_doy"))
  # 2.3.0: under the default "evaluate" policy it stays IN the domain as a
  # positive reference feature, so the uncovered ref2 becomes omission, not a
  # removed geometry.
  expect_true(ro$in_evaluation_domain[2])
  expect_equal(auth$observability_settings$domain_removed_ha, 0)
  expect_equal(auth$metrics$FN, 25)
  expect_equal(auth$metrics$Evaluable_Area_ha, 100 * 0.09)
  expect_true(auth$observability_settings$ref_obs_doy_authoritative)

  # the exclusion is announced, not silent (called directly: run_obs() wraps the
  # call in suppressWarnings(), which would swallow the very warning under test)
  # two warnings fire here under the default "evaluate" policy: the
  # authoritative-column one and the "kept for audit" one. Both are consumed so
  # the suite stays clean.
  expect_warning(expect_warning(
    suppressMessages(validate_fire_maps(
      input_shapefile      = fx$pred,
      ref_shapefile        = fx$ref,
      mask_shapefile       = fx$mask,
      burnable_raster      = fx$burnable,
      year_target          = 2000L,
      validation_dir       = file.path(tmp, "warn"),
      observability_raster = fx$obs,
      metrics_type         = "all",
      observability_mode   = "wholefire_fraction",
      observability_min_fraction = 0.75,
      ref_obs_doy_col      = "obs_required_doy",
      force_reprocess_ref  = TRUE,
      force_reprocess_pred = TRUE
    )),
    "Authoritative observability column 'obs_required_doy'"
  ), "are being EVALUATED as positive reference features")

  # the opt-in fallback restores the previous behaviour, explicitly
  fb <- run_obs(fx, file.path(tmp, "fb"),
                observability_mode = "wholefire_fraction",
                observability_min_fraction = 0.75,
                ref_obs_doy_col = "obs_required_doy",
                ref_obs_doy_fallback = "end_doy")
  rf <- fb$reference_observability
  expect_identical(rf$required_doy_source, c("obs_required_doy", "end_doy"))
  expect_equal(rf$fire_required_doy, c(100, 100))
  expect_identical(rf$observability_status, c("OBSERVABLE", "OBSERVABLE"))
  expect_false(fb$observability_settings$ref_obs_doy_authoritative)
})

test_that("T10c undetermined_cause names what was actually missing", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()

  # (a) no observability column at all, end_doy missing -> the end DOY was the
  #     only candidate. The cause must NOT mention an observability column.
  fx_noend <- build_obs_fixture(file.path(tmp, "fx1"), ref2_obs_rows = 5L,
                                ref2_end_doy = NA)
  r1 <- run_obs(fx_noend, file.path(tmp, "r1"),
                observability_mode = "wholefire_fraction",
                observability_min_fraction = 0.75)
  expect_identical(r1$reference_observability$undetermined_cause,
                   c(NA_character_, "no_end_doy"))

  # (b) observability column supplied WITH an explicit fallback, and neither
  #     column carries a date
  fx_both <- build_obs_fixture(file.path(tmp, "fx2"), ref2_obs_rows = 5L,
                               ref2_end_doy = NA,
                               extra_ref_cols = list(obs_required_doy = c(100, NA_real_)))
  r2 <- run_obs(fx_both, file.path(tmp, "r2"),
                observability_mode = "wholefire_fraction",
                observability_min_fraction = 0.75,
                ref_obs_doy_col = "obs_required_doy",
                ref_obs_doy_fallback = "end_doy")
  expect_identical(r2$reference_observability$undetermined_cause,
                   c(NA_character_, "no_obs_doy_and_no_end_doy"))

  # (c) legacy_any is the only mode that consults the start DOY
  fx_leg <- build_obs_fixture(file.path(tmp, "fx3"), ref2_obs_rows = 5L,
                              ref2_end_doy = NA, ref2_start_doy = NA)
  r3 <- run_obs(fx_leg, file.path(tmp, "r3"), observability_mode = "legacy_any")
  expect_identical(r3$reference_observability$undetermined_cause,
                   c(NA_character_, "no_end_or_start_doy"))

  # the summary reports the causes and the per-status areas
  s <- r1$observability_summary
  expect_identical(s$Undetermined_Causes, "no_end_doy=1")
  expect_equal(s$N_Undetermined_Obs_Date, 1L)
  expect_equal(s$Area_Undetermined_Obs_Date_ha, 25 * 0.09)
  expect_equal(s$Area_Not_Observable_ha, 0)
  expect_equal(s$Area_No_Observability_Data_ha, 0)
  # and the per-year CSV carries the new vocabulary
  und <- utils::read.csv(file.path(tmp, "r1", "VALIDATION", "02_OBSERVABILITY",
                                   "OBSERVABILITY_UNDETERMINED_2000.csv"))
  expect_true(all(c("observability_status", "area_domain_ha") %in% names(und)))
  expect_identical(unique(und$observability_status), "UNDETERMINED_OBS_DATE")
})

test_that("T10b the fallback policy moves the cache tag", {
  obs_raster <- withr::local_tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0,
                   ymax = 120, crs = "EPSG:3035")
  terra::values(r) <- 1:16
  names(r) <- "doy"
  terra::writeRaster(r, obs_raster, overwrite = TRUE)
  fp <- function(...) .vfm_observability_fingerprint(
    observability_raster = obs_raster, ref_end_doy_col = "end_doy",
    ref_start_doy_col = "start_doy",
    observability_mode = "wholefire_fraction", ...)

  auth <- fp(ref_obs_doy_col = "obs_required_doy", ref_obs_doy_fallback = "none")
  fb   <- fp(ref_obs_doy_col = "obs_required_doy", ref_obs_doy_fallback = "end_doy")
  expect_false(identical(auth, fb))
  # inert when no observability column is supplied
  expect_identical(
    fp(ref_obs_doy_col = NULL, ref_obs_doy_fallback = "none"),
    fp(ref_obs_doy_col = NULL, ref_obs_doy_fallback = "end_doy")
  )
})

# ---- T9: cache keys ----------------------------------------------------------

test_that("T9 the observability configuration changes the cache tags", {
  obs_raster <- withr::local_tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0,
                   ymax = 120, crs = "EPSG:3035")
  terra::values(r) <- 1:16
  names(r) <- "doy"
  terra::writeRaster(r, obs_raster, overwrite = TRUE)

  fp <- function(...) .vfm_observability_fingerprint(
    observability_raster = obs_raster,
    ref_end_doy_col = "end_doy", ref_start_doy_col = "start_doy", ...)

  legacy <- fp(observability_mode = "legacy_any")
  frac75 <- fp(observability_mode = "wholefire_fraction",
               observability_min_fraction = 0.75)
  frac90 <- fp(observability_mode = "wholefire_fraction",
               observability_min_fraction = 0.90)
  frac75_doy <- fp(observability_mode = "wholefire_fraction",
                   observability_min_fraction = 0.75, ref_obs_doy_col = "DOY")

  # mode, threshold and the observability DOY column all move the tag
  expect_false(identical(legacy, frac75))
  expect_false(identical(frac75, frac90))
  expect_false(identical(frac75, frac75_doy))
  expect_true(all(grepl("^obs-[A-F0-9]{8}$", c(legacy, frac75, frac90, frac75_doy))))

  # the threshold is inert under legacy_any and must NOT bust its cache
  expect_identical(
    fp(observability_mode = "legacy_any", observability_min_fraction = 0.75),
    fp(observability_mode = "legacy_any", observability_min_fraction = 0.10)
  )

  # the rule-family version is part of the tag
  expect_identical(.VFM_OBS_RULE_VERSION, "obsrule-3:wholefire{legacy_any|fraction}+undetpolicy")
  # bumped again in 2.3.0: the cached table gained in_evaluation_domain and the
  # OUTSIDE_BURNABLE_DOMAIN state, so no earlier table can be read back
  expect_identical(.VFM_REF_CACHE_SCHEMA, "refschema-5")

  # and the reference key inherits all of it
  k <- function(tag) .vfm_reference_cache_key(
    obs_cache_tag = tag, dom_cache_tag = "dom-DEADBEEF",
    ref_shapefile = NULL, min_area_reference_ha = NULL, dissolve_ref_by = NULL)
  expect_false(identical(k(legacy), k(frac75)))
})

test_that("T9c a cache HIT reproduces the fraction-mode metrics exactly", {
  skip_if_not_installed("terra")
  # The excluded-fire geometries are needed to rebuild the evaluation domain,
  # and on the cache-hit path they are read back from the GeoPackage instead of
  # being built in memory. A second call must therefore be identical to the
  # first, not silently lose the exclusion.
  tmp <- withr::local_tempdir()
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 0L,
                          pred_over_ref2 = TRUE)
  d <- file.path(tmp, "run")

  first <- run_obs(fx, d, observability_mode = "wholefire_fraction",
                   observability_min_fraction = 0.75)
  second <- suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile      = fx$pred,
    ref_shapefile        = fx$ref,
    mask_shapefile       = fx$mask,
    burnable_raster      = fx$burnable,
    year_target          = 2000L,
    validation_dir       = d,
    observability_raster = fx$obs,
    metrics_type         = "all",
    observability_mode   = "wholefire_fraction",
    observability_min_fraction = 0.75,
    force_reprocess_ref  = FALSE,   # <- cache HIT
    force_reprocess_pred = FALSE
  )))

  keep <- c("TP", "FP", "FN", "TN", "Precision", "Recall", "F1", "IoU",
            "Evaluable_Area_ha")
  expect_equal(as.list(first$metrics[, keep, with = FALSE]),
               as.list(second$metrics[, keep, with = FALSE]))
  expect_equal(first$observability_settings$domain_removed_ha,
               second$observability_settings$domain_removed_ha)
  expect_true(second$observability_settings$domain_masked)
  expect_equal(second$metrics$FP, 0)
})

test_that("T9b the prediction cache key depends on observability only when the domain is masked", {
  skip_if_not_installed("terra")
  tmp <- withr::local_tempdir()
  fx <- build_obs_fixture(file.path(tmp, "fx"), ref2_obs_rows = 0L)

  cache_files <- function(d) {
    list.files(file.path(d, "VALIDATION", "_CACHE"),
               pattern = "^predicted_fire_mask_")
  }

  d_leg <- file.path(tmp, "leg")
  run_obs(fx, d_leg, observability_mode = "legacy_any")
  f_leg <- cache_files(d_leg)
  expect_length(f_leg, 1L)
  # legacy keeps the historical shape: <input>_<dom>, no observability tag
  expect_false(grepl("obs-[A-F0-9]{8}", f_leg))

  d_new <- file.path(tmp, "new")
  run_obs(fx, d_new, observability_mode = "wholefire_fraction",
          observability_min_fraction = 0.75)
  f_new <- cache_files(d_new)
  expect_length(f_new, 1L)
  expect_true(grepl("obs-[A-F0-9]{8}", f_new))
  expect_false(identical(f_leg, f_new))

  # two different thresholds must not share a prediction cache file
  d_new2 <- file.path(tmp, "new2")
  run_obs(fx, d_new2, observability_mode = "wholefire_fraction",
          observability_min_fraction = 0.50)
  expect_false(identical(f_new, cache_files(d_new2)))
})
