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
    year = c(2000L, 2000L),
    end_doy = c(100, NA),
    start_doy = c(90, 150),
    geometry = sf::st_sfc(ref1_poly, ref2_poly, crs = 3035)
  )
  ref_path <- file.path(dir, "ref.shp")
  suppressWarnings(sf::st_write(ref_sf, ref_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  observability_sf <- sf::st_sf(
    obs_doy = c(120, 140),
    geometry = sf::st_sfc(ref1_poly, ref2_poly, crs = 3035)
  )
  observability <- burnable
  terra::values(observability) <- NA_real_
  observability <- terra::rasterize(terra::vect(observability_sf), observability,
                                    field = "obs_doy", background = NA)
  observability_path <- file.path(dir, "observability.tif")
  terra::writeRaster(observability, observability_path, overwrite = TRUE,
                     datatype = "FLT4S")

  observability_multiband <- c(burnable, observability)
  names(observability_multiband) <- c("rbr", "doy")
  observability_multiband_path <- file.path(dir, "observability_multiband.tif")
  terra::writeRaster(observability_multiband, observability_multiband_path,
                     overwrite = TRUE, datatype = "FLT4S")

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
    observability = observability_path,
    observability_multiband = observability_multiband_path,
    mask     = mask_path,
    ref      = ref_path,
    pred     = pred_path
  )
}

run_validate_synth <- function(dir, threshold_min_detected = 10,
                                year_target = 2000L,
                                observability_raster = NULL) {
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
      observability_raster   = observability_raster,
      metrics_type           = "all"
    )
  ))
}

validation_subdirs <- function(dir) {
  file.path(
    dir,
    "VALIDATION",
    c("01_SUMMARY", "02_OBSERVABILITY", "03_STRATA", "04_ERROR_LAYERS", "_CACHE")
  )
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

  expect_true(all(dir.exists(validation_subdirs(tmp))))
  expect_true(file.exists(file.path(
    tmp, "VALIDATION", "01_SUMMARY", "pixel_metrics_global_2000.csv"
  )))
  expect_true(file.exists(file.path(
    tmp, "VALIDATION", "01_SUMMARY", "fire_metrics_global_2000.csv"
  )))
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

test_that("observability filter excludes non-observable reference polygons and returns summary", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_obs_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_synth_validation_fixture(tmp)

  expect_warning(
    res <- suppressMessages(validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      observability_raster   = fx$observability,
      year_target            = 2000L,
      validation_dir         = tmp,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = 10,
      metrics_type           = "all"
    )),
    "Temporal observability filter excluded 1 of 2 reference polygons"
  )

  obs <- res$reference_observability
  expect_true(!is.null(obs) && nrow(obs) == 2L)
  expect_true(all(c("obs_ref_doy", "obs_doy_max", "obs_doy_margin",
                    "observable_flag", "observable_reason") %in% names(obs)))

  ref2 <- obs[obs$id == 2, ]
  expect_equal(ref2$obs_ref_doy, 150)
  expect_equal(ref2$obs_doy_max, 140)
  expect_equal(ref2$obs_doy_margin, -10)
  expect_false(ref2$observable_flag)
  expect_equal(ref2$observable_reason, "obs_before_reference_doy")

  expect_equal(unname(res$polygon_summary$N_Reference_Polygons), 1L)
  expect_equal(unname(res$polygon_summary$N_Detected_Polygons), 1L)
  expect_equal(unname(res$metrics$TP), 25)
  expect_equal(unname(res$metrics$FP), 5)
  expect_equal(unname(res$metrics$FN), 0)
  expect_equal(unname(res$metrics$TN), 70)
  expect_true(file.exists(file.path(
    tmp, "VALIDATION", "01_SUMMARY", "pixel_metrics_global_2000.csv"
  )))
  expect_true(file.exists(file.path(
    tmp, "VALIDATION", "01_SUMMARY", "fire_metrics_global_2000.csv"
  )))

  # Content-aware cache contract (0.10.1): observability-tagged caches carry the
  # full content-aware tag `obs-<8HEX>_dom-<8HEX>_ref-<8HEX>` (observability +
  # burnable-domain + reference identity) immediately before the extension.
  obs_dir <- file.path(tmp, "VALIDATION", "02_OBSERVABILITY")
  obs_csv <- list.files(
    obs_dir,
    pattern = "^reference_fires_observability_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.csv$",
    full.names = TRUE
  )
  expect_length(obs_csv, 1L)
  # The mandatory dom+ref suffix: a name WITHOUT it must NOT be produced.
  expect_length(
    list.files(obs_dir, pattern = "^reference_fires_observability_2000_obs-[A-F0-9]{8}\\.csv$"),
    0L
  )

  obs_gpkg <- list.files(
    obs_dir,
    pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(obs_gpkg, 1L)
  # OLD nameless-hash filename (obs tag directly before extension) must be absent.
  expect_length(
    list.files(obs_dir, pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}\\.gpkg$"),
    0L
  )

  dropped_sf <- suppressWarnings(sf::st_read(obs_gpkg, quiet = TRUE))
  expect_s3_class(dropped_sf, "sf")
  expect_equal(nrow(dropped_sf), 1L)
  expect_equal(dropped_sf$id, 2)
  expect_true(all(c("obs_ref_doy", "obs_ref_doy_source", "obs_doy_max",
                    "obs_doy_margin", "observable_reason") %in% names(dropped_sf)))
  expect_false(dropped_sf$observable_flag)
  expect_equal(dropped_sf$observable_reason, "obs_before_reference_doy")
  expect_false(any(sf::st_is_empty(dropped_sf)))
})

test_that("observability filter uses the doy layer from a multi-band raster", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_obs_multiband_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_synth_validation_fixture(tmp)

  expect_warning(
    res <- suppressMessages(validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      observability_raster   = fx$observability_multiband,
      year_target            = 2000L,
      validation_dir         = tmp,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = 10,
      metrics_type           = "all"
    )),
    "Temporal observability filter excluded 1 of 2 reference polygons"
  )

  ref2 <- res$reference_observability[res$reference_observability$id == 2, ]
  expect_equal(ref2$obs_doy_max, 140)
  expect_equal(ref2$obs_ref_doy_source, "start_doy")
  expect_false(ref2$observable_flag)
  expect_equal(unname(res$polygon_summary$N_Reference_Polygons), 1L)
})

test_that("observability filter skips the excluded-polygons gpkg when all references are observable", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_obs_all_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_synth_validation_fixture(tmp)

  observability_all <- terra::rast(fx$observability)
  observability_all[!is.na(observability_all)] <- 200
  observability_all_path <- file.path(tmp, "observability_all.tif")
  terra::writeRaster(observability_all, observability_all_path, overwrite = TRUE,
                     datatype = "FLT4S")

  res <- suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile        = fx$pred,
    ref_shapefile          = fx$ref,
    mask_shapefile         = fx$mask,
    burnable_raster        = fx$burnable,
    observability_raster   = observability_all_path,
    year_target            = 2000L,
    validation_dir         = tmp,
    binary_burnable        = TRUE,
    burnable_threshold     = 0.5,
    threshold_completely_detected = 90,
    threshold_min_detected = 10,
    metrics_type           = "all"
  )))

  expect_true(all(res$reference_observability$observable_flag))
  # No not-observable gpkg is produced when every reference is observable —
  # checked against the strict content-aware contract (obs hash + mandatory
  # `_dom-<8HEX>_ref-<8HEX>`) AND against the OLD nameless-hash pattern, so neither variant
  # is silently produced.
  obs_dir <- file.path(tmp, "VALIDATION", "02_OBSERVABILITY")
  obs_gpkg <- list.files(
    obs_dir,
    pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(obs_gpkg, 0L)
  expect_length(
    list.files(obs_dir, pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}\\.gpkg$"),
    0L
  )
})

test_that("reference cache is separated by observability tag", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_obs_cache_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_synth_validation_fixture(tmp)

  res_no_obs <- suppressMessages(suppressWarnings(
    validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      year_target            = 2000L,
      validation_dir         = tmp,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = 10,
      metrics_type           = "all"
    )
  ))
  expect_equal(unname(res_no_obs$polygon_summary$N_Reference_Polygons), 2L)

  expect_warning(
    res_obs <- suppressMessages(validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      observability_raster   = fx$observability,
      year_target            = 2000L,
      validation_dir         = tmp,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = 10,
      metrics_type           = "all"
    )),
    "Temporal observability filter excluded 1 of 2 reference polygons"
  )

  expect_equal(unname(res_obs$polygon_summary$N_Reference_Polygons), 1L)

  cache_dir <- file.path(tmp, "VALIDATION", "_CACHE")
  obs_subdir <- file.path(tmp, "VALIDATION", "02_OBSERVABILITY")

  # ---- Content-aware cache contract (f9130a6) ----------------------------
  # Point 3 + point 1: the no-observability run writes EXACTLY ONE processed
  # cache tagged `obs-none` followed by the MANDATORY `_dom-<8HEX>_ref-<8HEX>` ref-input
  # fingerprint immediately before the extension.
  no_obs_cache <- list.files(
    cache_dir,
    pattern = "^reference_fires_processed_2000_obs-none_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(no_obs_cache, 1L)
  expect_s3_class(suppressWarnings(sf::st_read(no_obs_cache, quiet = TRUE)), "sf")

  # Point 8: the OLD nameless-hash filename (obs tag directly before extension)
  # must NOT be produced and must NOT be silently accepted/reused. Assert absence
  # of ANY processed cache whose obs tag is followed directly by `.gpkg`
  # (obs-none OR obs-<hex>) with no dom/ref segment.
  expect_false(file.exists(file.path(
    cache_dir, "reference_fires_processed_2000_obs-none.gpkg"
  )))
  expect_length(
    list.files(
      cache_dir,
      pattern = "^reference_fires_processed_2000_obs-(none|[A-F0-9]{8})\\.gpkg$"
    ),
    0L
  )

  # Point 2 + point 1: the observability run writes EXACTLY ONE processed cache
  # carrying the observability hash `obs-<8HEX>` AND the mandatory `_dom-<8HEX>_ref-<8HEX>`
  # suffix. obs-none and obs-<hex> are distinct files (cache separated by tag).
  obs_cache <- list.files(
    cache_dir,
    pattern = "^reference_fires_processed_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(obs_cache, 1L)
  expect_s3_class(suppressWarnings(sf::st_read(obs_cache, quiet = TRUE)), "sf")
  expect_false(identical(basename(no_obs_cache), basename(obs_cache)))

  # Point 6: the SAME inputs produce the SAME cache name. Capture the obs-run
  # `_dom-<8HEX>_ref-<8HEX>` fingerprint; a second identical obs run must reuse it (no new
  # `ref-` tag appears, the cache count under obs-<hex> stays at exactly 1).
  in_hash_run1 <- sub(
    "^.*_(ref-[A-F0-9]{8})\\.gpkg$", "\\1", basename(obs_cache)
  )
  expect_match(in_hash_run1, "^ref-[A-F0-9]{8}$")

  obs_not_observable <- list.files(
    obs_subdir,
    pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(obs_not_observable, 1L)
  # Point 8 (observability dir): no OLD nameless-hash not-observable gpkg.
  expect_length(
    list.files(
      obs_subdir,
      pattern = "^reference_fires_not_observable_2000_obs-[A-F0-9]{8}\\.gpkg$"
    ),
    0L
  )

  unlink(obs_not_observable)
  expect_false(file.exists(obs_not_observable))

  expect_warning(
    suppressMessages(validate_fire_maps(
      input_shapefile        = fx$pred,
      ref_shapefile          = fx$ref,
      mask_shapefile         = fx$mask,
      burnable_raster        = fx$burnable,
      observability_raster   = fx$observability,
      year_target            = 2000L,
      validation_dir         = tmp,
      binary_burnable        = TRUE,
      burnable_threshold     = 0.5,
      threshold_completely_detected = 90,
      threshold_min_detected = 10,
      metrics_type           = "all"
    )),
    "Temporal observability filter excluded 1 of 2 reference polygons"
  )
  expect_true(file.exists(obs_not_observable))

  # Point 6 (deterministic name): the rerun used the SAME inputs, so the
  # processed cache name — including its `_dom-<8HEX>_ref-<8HEX>` fingerprint — is stable.
  # There is still EXACTLY ONE obs-tagged processed cache, under the SAME
  # `in-` hash captured from run 1 (no second `ref-` tag was minted).
  obs_cache_run2 <- list.files(
    cache_dir,
    pattern = "^reference_fires_processed_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  expect_length(obs_cache_run2, 1L)
  in_hash_run2 <- sub(
    "^.*_(ref-[A-F0-9]{8})\\.gpkg$", "\\1", basename(obs_cache_run2)
  )
  expect_identical(in_hash_run2, in_hash_run1)

  # Point 7 (content-aware busting): perturb a RELEVANT ref input — drop one
  # reference feature and rewrite ref.shp (changes nrow + .dbf size/mtime).
  # The `_dom-<8HEX>_ref-<8HEX>` fingerprint must change, so a SECOND processed cache file
  # under a DIFFERENT `ref-` tag appears and the old one is NOT reused.
  ref_perturbed <- sf::st_read(fx$ref, quiet = TRUE)
  ref_perturbed <- ref_perturbed[ref_perturbed$id == 1L, , drop = FALSE]
  suppressWarnings(sf::st_write(ref_perturbed, fx$ref, quiet = TRUE,
                                delete_dsn = TRUE))

  suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile        = fx$pred,
    ref_shapefile          = fx$ref,
    mask_shapefile         = fx$mask,
    burnable_raster        = fx$burnable,
    observability_raster   = fx$observability,
    year_target            = 2000L,
    validation_dir         = tmp,
    binary_burnable        = TRUE,
    burnable_threshold     = 0.5,
    threshold_completely_detected = 90,
    threshold_min_detected = 10,
    metrics_type           = "all"
  )))

  obs_caches_after <- list.files(
    cache_dir,
    pattern = "^reference_fires_processed_2000_obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}\\.gpkg$",
    full.names = TRUE
  )
  # A second, distinctly-tagged cache now exists alongside the first.
  expect_length(obs_caches_after, 2L)
  in_hashes_after <- sub(
    "^.*_(ref-[A-F0-9]{8})\\.gpkg$", "\\1", basename(obs_caches_after)
  )
  expect_length(unique(in_hashes_after), 2L)
  expect_true(in_hash_run1 %in% in_hashes_after)         # old not deleted...
  expect_true(any(in_hashes_after != in_hash_run1))      # ...and not reused
})

test_that("legacy flat outputs are moved out of the validation root", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  tmp <- tempfile("of_vfm_legacy_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_synth_validation_fixture(tmp)

  root_dir <- file.path(tmp, "VALIDATION")
  legacy_dir <- file.path(root_dir, "_LEGACY_FLAT_OUTPUTS")
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(legacy_dir, recursive = TRUE, showWarnings = FALSE)
  legacy_csv <- file.path(root_dir, "metrics_summary_2000.csv")
  archived_csv <- file.path(legacy_dir, "metrics_summary_2000.csv")
  legacy_csv_humanized <- file.path(root_dir, "pixel_metrics_global_2000.csv")
  legacy_gpkg <- file.path(root_dir, "ref_polygons_not_detected_keep.gpkg")
  writeLines("old-root", legacy_csv)
  writeLines("already-archived", archived_csv)
  writeLines("old-humanized-root", legacy_csv_humanized)
  writeLines("old-gpkg", legacy_gpkg)

  suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile        = fx$pred,
    ref_shapefile          = fx$ref,
    mask_shapefile         = fx$mask,
    burnable_raster        = fx$burnable,
    year_target            = 2000L,
    validation_dir         = tmp,
    binary_burnable        = TRUE,
    burnable_threshold     = 0.5,
    threshold_completely_detected = 90,
    threshold_min_detected = 10,
    metrics_type           = "all"
  )))

  expect_false(file.exists(legacy_csv))
  expect_false(file.exists(legacy_csv_humanized))
  expect_false(file.exists(legacy_gpkg))
  migrated_metrics <- list.files(
    legacy_dir,
    pattern = "^metrics_summary_2000(?:__legacy[0-9]{2})?\\.csv$",
    full.names = TRUE
  )
  expect_length(migrated_metrics, 2L)
  expect_equal(
    unname(sort(vapply(migrated_metrics, readLines, character(1), warn = FALSE))),
    sort(c("already-archived", "old-root"))
  )
  expect_true(file.exists(file.path(
    legacy_dir, "pixel_metrics_global_2000.csv"
  )))
  expect_true(file.exists(file.path(
    legacy_dir, "ref_polygons_not_detected_keep.gpkg"
  )))

  legacy_count_before <- length(list.files(legacy_dir, full.names = TRUE))
  suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile        = fx$pred,
    ref_shapefile          = fx$ref,
    mask_shapefile         = fx$mask,
    burnable_raster        = fx$burnable,
    year_target            = 2000L,
    validation_dir         = tmp,
    binary_burnable        = TRUE,
    burnable_threshold     = 0.5,
    threshold_completely_detected = 90,
    threshold_min_detected = 10,
    metrics_type           = "all"
  )))
  expect_equal(length(list.files(legacy_dir, full.names = TRUE)), legacy_count_before)
})

# §N+7.6 (2026-05-29) — surfaced by the 2005/balanced deterministic
# smoke. The Iberian-peninsula mask carries a column called `Id` and the
# EFFIS reference shapefile carries a column called `id`; after
# `sf::st_intersection(ref_polygons, mask_geom)` both names landed in
# ref_polygons unchanged, then `sf::st_write()` to the GPKG cache
# aborted with "Field count reached: duplicate names present?" because
# SQLite identifiers are case-insensitive. The previous F8b drop only
# removed the `FID` column. The fix at the top of `validate_fire_maps()`
# reduces mask_geom to its geometry column (we never consume mask
# attributes downstream); this test pins both the body and the
# functional behaviour.

test_that("validate_fire_maps() body keeps only mask geometry (§N+7.6)", {
  body_src <- paste(
    deparse(body(OtsuFire::validate_fire_maps)),
    collapse = "\n"
  )
  # Geometry-only narrowing must reach mask_geom.
  expect_true(grepl(
    "mask_geom\\s*<-\\s*mask_geom\\[\\s*,\\s*geom_col\\s*,\\s*drop\\s*=\\s*FALSE\\s*\\]",
    body_src
  ))
  expect_true(grepl("attr\\(mask_geom,\\s*\"sf_column\"\\)", body_src))
  # The narrow F8b drop ("FID" only) must be gone.
  expect_false(grepl(
    "!\\(names\\(mask_geom\\)\\s*%in%\\s*\"FID\"\\)",
    body_src
  ))
})

test_that("validate_fire_maps() handles case-insensitive id collision between mask and ref", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  dir <- file.path(tempdir(), "validate_n76")
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)

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
    0, 0, 300, 0, 300, 300, 0, 300, 0, 0
  ), ncol = 2, byrow = TRUE)))
  # Mask carries uppercase "Id" (mirrors Iberian_Peninsula_mask_3035.shp).
  mask_sf <- sf::st_sf(
    Id   = 1L,
    PAIS = "ESP",
    geometry = sf::st_sfc(mask_poly, crs = 3035)
  )
  mask_path <- file.path(dir, "mask_case.shp")
  suppressWarnings(sf::st_write(mask_sf, mask_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  ref_poly <- sf::st_polygon(list(matrix(c(
    50, 50, 200, 50, 200, 200, 50, 200, 50, 50
  ), ncol = 2, byrow = TRUE)))
  # Reference carries lowercase "id" (mirrors EFFIS shapefile schema).
  ref_sf <- sf::st_sf(
    id        = 1L,
    year      = 2005L,
    end_doy   = 200,
    start_doy = 100,
    geometry  = sf::st_sfc(ref_poly, crs = 3035)
  )
  ref_path <- file.path(dir, "ref.shp")
  suppressWarnings(sf::st_write(ref_sf, ref_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  pred_poly <- ref_poly
  pred_sf <- sf::st_sf(
    id       = 1L,
    geometry = sf::st_sfc(pred_poly, crs = 3035)
  )
  pred_path <- file.path(dir, "pred.shp")
  suppressWarnings(sf::st_write(pred_sf, pred_path, quiet = TRUE,
                                 delete_dsn = TRUE))

  # Without the §N+7.6 fix, the cache GPKG write at
  # `sf::st_write(ref_polygons, ref_vec_cache, ...)` aborts with
  # "Field count reached: duplicate names present?". With the fix the
  # call completes cleanly and the cache GPKG is written.
  expect_no_error(
    suppressWarnings(suppressMessages(
      validate_fire_maps(
        input_shapefile        = pred_path,
        ref_shapefile          = ref_path,
        mask_shapefile         = mask_path,
        burnable_raster        = burnable_path,
        year_target            = 2005L,
        validation_dir         = dir,
        binary_burnable        = TRUE,
        burnable_threshold     = 0.5,
        threshold_completely_detected = 90,
        threshold_min_detected = 10,
        metrics_type           = "all"
      )
    ))
  )
  # The cache GPKG must exist on disk and carry no mask-side attributes.
  cache_dir <- file.path(dir, "VALIDATION", "_CACHE")
  cache_gpkg <- list.files(cache_dir,
                            pattern = "^reference_fires_processed_2005.*\\.gpkg$",
                            full.names = TRUE)
  expect_length(cache_gpkg, 1L)
  cached <- sf::st_read(cache_gpkg, quiet = TRUE)
  # The reference attributes survive; mask "Id" and "PAIS" do not.
  expect_true("id" %in% names(cached))
  expect_false(any(c("Id", "PAIS") %in% names(cached)))
})
