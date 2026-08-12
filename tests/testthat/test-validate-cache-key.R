# Content-aware cache contract for validate_fire_maps() (0.10.1).
#
# The reference-side cache must be invalidated whenever the CONTENT or
# methodological identity of any input that shaped it changes: the
# observability raster content / layer / DOY columns / rule, the burnable mask,
# resolution, CRS, grid/extent, the EFFIS reference, the burnable thresholding
# and the reference-only options. It must NOT be invalidated by cosmetic run
# settings (logging), and `force_reprocess_ref = TRUE` must always rebuild.
#
# Tests 1-7 + 9 exercise the pure key helpers directly (fast, deterministic).
# Tests 8 + 10 exercise the real validator end-to-end on a tiny synthetic fixture.

vc_ns <- asNamespace("OtsuFire")

# --- tiny fixture builders ---------------------------------------------------

# Binary burnable-style raster. Extent fixed at 0..(ncell*base) so resolution
# can be varied (via `res`) while keeping the same footprint; CRS configurable.
vc_mk_rast <- function(path, res = 30, crs = "EPSG:3035", val = 1,
                       extent = 300, doy_name = FALSE, doy_vals = NULL) {
  n <- as.integer(extent / res)
  r <- terra::rast(nrows = n, ncols = n,
                   xmin = 0, xmax = extent, ymin = 0, ymax = extent, crs = crs)
  terra::values(r) <- if (is.null(doy_vals)) val else doy_vals
  if (doy_name) names(r) <- "doy"
  terra::writeRaster(r, path, overwrite = TRUE, datatype = "FLT4S")
  path
}

vc_mk_vec <- function(path, n_feat = 2L, shift = 0, crs = 3035,
                      with_doy = FALSE) {
  polys <- lapply(seq_len(n_feat), function(i) {
    x0 <- shift + (i - 1L) * 50
    sf::st_polygon(list(rbind(
      c(x0, 0), c(x0 + 40, 0), c(x0 + 40, 40), c(x0, 40), c(x0, 0))))
  })
  df <- data.frame(id = seq_len(n_feat))
  if (with_doy) {
    df$year <- 2000L
    df$end_doy <- 100
    df$start_doy <- 90
  }
  obj <- sf::st_sf(df, geometry = sf::st_sfc(polys, crs = crs))
  suppressWarnings(sf::st_write(obj, path, quiet = TRUE, delete_dsn = TRUE))
  path
}

# Full reference cache key from raw inputs.
vc_key <- function(burnable, mask, ref, obs = NULL,
                   end = "end_doy", start = "start_doy",
                   binary = TRUE, classes = NULL, threshold = 0.5,
                   minha = NULL, dissolve = NULL) {
  ot <- vc_ns$.vfm_observability_fingerprint(obs, end, start)
  dt <- vc_ns$.vfm_domain_fingerprint(burnable, mask, binary, classes, threshold)
  vc_ns$.vfm_reference_cache_key(ot, dt, ref, minha, dissolve)
}

vc_dir <- function() {
  d <- tempfile("vc_cache_")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# ---------------------------------------------------------------------------
# 1. Same file, same content -> cache reused (identical key).
# ---------------------------------------------------------------------------
test_that("1. identical inputs produce an identical (reusable) cache key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  bn <- vc_mk_rast(file.path(d, "burn.tif"))
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  ob <- vc_mk_rast(file.path(d, "obs.tif"), doy_name = TRUE, val = 200)

  k1 <- vc_key(bn, mk, rf, obs = ob)
  k2 <- vc_key(bn, mk, rf, obs = ob)
  expect_identical(k1, k2)
  expect_match(k1, "^obs-[A-F0-9]{8}_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}$")
})

# ---------------------------------------------------------------------------
# 2. Same name, different content -> cache invalidated (key changes).
# ---------------------------------------------------------------------------
test_that("2. same path, different observability content invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  bn <- vc_mk_rast(file.path(d, "burn.tif"))
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  ob <- file.path(d, "obs.tif")

  vc_mk_rast(ob, doy_name = TRUE, doy_vals = rep(120, 100))
  k1 <- vc_key(bn, mk, rf, obs = ob)
  # Overwrite the SAME path with different DOY content.
  vc_mk_rast(ob, doy_name = TRUE, doy_vals = rep(250, 100))
  k2 <- vc_key(bn, mk, rf, obs = ob)
  expect_false(identical(k1, k2))
})

# ---------------------------------------------------------------------------
# 3. Different burnable mask -> cache invalidated.
# ---------------------------------------------------------------------------
test_that("3. a different burnable mask invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  b1 <- vc_mk_rast(file.path(d, "burn1.tif"), val = 1)
  b2 <- vc_mk_rast(file.path(d, "burn2.tif"), val = 0)  # different content
  expect_false(identical(vc_key(b1, mk, rf), vc_key(b2, mk, rf)))

  # Also: a different study-area MASK invalidates.
  m2 <- vc_mk_vec(file.path(d, "mask2.shp"), n_feat = 2L, shift = 5)
  expect_false(identical(vc_key(b1, mk, rf), vc_key(b1, m2, rf)))
})

# ---------------------------------------------------------------------------
# 4. Different resolution -> cache invalidated.
# ---------------------------------------------------------------------------
test_that("4. a different burnable resolution invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  b30 <- vc_mk_rast(file.path(d, "b30.tif"), res = 30)  # 10x10
  b60 <- vc_mk_rast(file.path(d, "b60.tif"), res = 60)  # 5x5, same extent
  expect_false(identical(vc_key(b30, mk, rf), vc_key(b60, mk, rf)))
})

# ---------------------------------------------------------------------------
# 5. Different CRS / grid -> cache invalidated.
# ---------------------------------------------------------------------------
test_that("5. a different CRS / grid invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  b3035  <- vc_mk_rast(file.path(d, "b3035.tif"),  crs = "EPSG:3035")
  b25830 <- vc_mk_rast(file.path(d, "b25830.tif"), crs = "EPSG:25830")
  expect_false(identical(vc_key(b3035, mk, rf), vc_key(b25830, mk, rf)))
})

# ---------------------------------------------------------------------------
# 6. Different DOY band -> cache invalidated.
# ---------------------------------------------------------------------------
test_that("6. a different observability DOY band invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  bn <- vc_mk_rast(file.path(d, "burn.tif"))
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  o1 <- vc_mk_rast(file.path(d, "o1.tif"), doy_name = TRUE, doy_vals = rep(120, 100))
  o2 <- vc_mk_rast(file.path(d, "o2.tif"), doy_name = TRUE, doy_vals = rep(220, 100))
  expect_false(identical(vc_key(bn, mk, rf, obs = o1),
                         vc_key(bn, mk, rf, obs = o2)))

  # The DOY COLUMNS on the reference are part of the obs identity too.
  expect_false(identical(
    vc_key(bn, mk, rf, obs = o1, end = "end_doy"),
    vc_key(bn, mk, rf, obs = o1, end = "fire_doy")
  ))
  # Presence/absence of observability is distinguished (obs-none vs obs-hash).
  expect_match(vc_key(bn, mk, rf, obs = NULL),
               "^obs-none_dom-[A-F0-9]{8}_ref-[A-F0-9]{8}$")
})

# ---------------------------------------------------------------------------
# 7. Different EFFIS reference -> cache invalidated.
# ---------------------------------------------------------------------------
test_that("7. a different EFFIS reference invalidates the key", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  bn <- vc_mk_rast(file.path(d, "burn.tif"))
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  r1 <- vc_mk_vec(file.path(d, "ref1.shp"), n_feat = 2L, with_doy = TRUE)
  r2 <- vc_mk_vec(file.path(d, "ref2.shp"), n_feat = 3L, with_doy = TRUE)
  expect_false(identical(vc_key(bn, mk, r1), vc_key(bn, mk, r2)))

  # Reference-only options (min area, dissolve field) are part of the key too.
  expect_false(identical(vc_key(bn, mk, r1),
                         vc_key(bn, mk, r1, minha = 5)))
  expect_false(identical(vc_key(bn, mk, r1),
                         vc_key(bn, mk, r1, dissolve = "id")))
  # Burnable threshold is part of the domain identity.
  expect_false(identical(vc_key(bn, mk, r1, threshold = 0.5),
                         vc_key(bn, mk, r1, threshold = 0.9)))
})

# ---------------------------------------------------------------------------
# 9. Logging / cosmetic settings do NOT invalidate the cache key.
# ---------------------------------------------------------------------------
test_that("9. the cache key depends only on content/methodology, not on logging", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  # The key helpers take NO logging / verbose / quiet / message argument: their
  # only inputs are content + methodology, so cosmetic run settings cannot bust
  # the cache.
  for (fn in c(".vfm_reference_cache_key", ".vfm_domain_fingerprint",
               ".vfm_observability_fingerprint", ".vfm_raster_token",
               ".vfm_vector_token")) {
    fm <- names(formals(get(fn, envir = vc_ns)))
    expect_false(any(grepl("verbose|quiet|log|message|debug", fm,
                           ignore.case = TRUE)),
                 info = paste(fn, "must not take a logging-style argument"))
  }
  # And re-deriving the key from identical inputs is stable across calls.
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  bn <- vc_mk_rast(file.path(d, "burn.tif"))
  mk <- vc_mk_vec(file.path(d, "mask.shp"), n_feat = 1L)
  rf <- vc_mk_vec(file.path(d, "ref.shp"), with_doy = TRUE)
  expect_identical(vc_key(bn, mk, rf), vc_key(bn, mk, rf))
})

# --- behavioral tests on the real validator ---------------------------------

vc_synth <- function(dir) {
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

vc_run <- function(fx, dir) {
  suppressWarnings(suppressMessages(validate_fire_maps(
    input_shapefile = fx$pred, ref_shapefile = fx$ref, mask_shapefile = fx$mask,
    burnable_raster = fx$burnable, observability_raster = fx$obs,
    year_target = 2000L, validation_dir = dir,
    binary_burnable = TRUE, burnable_threshold = 0.5,
    threshold_completely_detected = 90, threshold_min_detected = 10,
    metrics_type = "all")))
}

# ---------------------------------------------------------------------------
# 8. force_reprocess_ref = TRUE ignores (rebuilds) the cache.
# ---------------------------------------------------------------------------
test_that("8. force_reprocess_ref = TRUE rebuilds the reference cache", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  fx <- vc_synth(d)

  suppressWarnings(vc_run(fx, d))  # build cache
  cache_dir <- file.path(d, "VALIDATION", "_CACHE")
  cache_files <- list.files(cache_dir,
                            pattern = "^reference_fires_processed_2000_.*\\.gpkg$",
                            full.names = TRUE)
  expect_length(cache_files, 1L)

  # A force_reprocess_ref run must announce the rebuild and re-derive the cache.
  expect_message(
    suppressWarnings(validate_fire_maps(
      input_shapefile = fx$pred, ref_shapefile = fx$ref, mask_shapefile = fx$mask,
      burnable_raster = fx$burnable, observability_raster = fx$obs,
      year_target = 2000L, validation_dir = d,
      binary_burnable = TRUE, burnable_threshold = 0.5,
      threshold_completely_detected = 90, threshold_min_detected = 10,
      metrics_type = "all", force_reprocess_ref = TRUE)),
    "Force reprocessing reference"
  )
  # Still exactly one processed cache (same deterministic key), rebuilt in place.
  expect_length(
    list.files(cache_dir, pattern = "^reference_fires_processed_2000_.*\\.gpkg$"),
    1L
  )
})

# ---------------------------------------------------------------------------
# 10. Identical inputs -> identical metrics (cache reuse is value-preserving).
# ---------------------------------------------------------------------------
test_that("10. identical inputs yield identical metrics (build vs cache reuse)", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  d <- vc_dir(); on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  fx <- vc_synth(d)

  res1 <- vc_run(fx, d)   # cold: builds cache
  res2 <- vc_run(fx, d)   # warm: reuses cache

  m1 <- res1$metrics; m2 <- res2$metrics
  for (col in c("TP", "FP", "FN", "TN", "Recall", "Precision", "F1", "IoU")) {
    if (col %in% names(m1)) {
      expect_identical(m1[[col]], m2[[col]],
                       info = paste("metric", col, "must be identical on cache reuse"))
    }
  }
  expect_identical(res1$polygon_summary$N_Reference_Polygons,
                   res2$polygon_summary$N_Reference_Polygons)
  expect_identical(res1$polygon_summary$N_Detected_Polygons,
                   res2$polygon_summary$N_Detected_Polygons)
})
