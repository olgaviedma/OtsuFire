# AS01 (OtsuFire 0.3.0): external-tool paths plumbed through
# config$tool_paths -> dispatcher bindings -> orchestrator ->
# build_unburned_from_legacy_pipeline.

mk_tmp_tif_tp <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_gpkg_tp <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                                c(0,1), c(0,0)))),
                    crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}

# --- T18 ------------------------------------------------------------
test_that("T18: cfg$tool_paths is populated and exposed via dispatcher", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    options = list(
      python_exe             = "/fake/python",
      gdal_polygonize_script = "/fake/gdal_polygonize.py",
      gdalwarp_path          = "/fake/gdalwarp",
      ogr2ogr_exe            = "/fake/ogr2ogr"
    )
  )

  expect_identical(cfg$tool_paths$python_exe, "/fake/python")
  expect_identical(cfg$tool_paths$gdal_polygonize_script, "/fake/gdal_polygonize.py")
  expect_identical(cfg$tool_paths$gdalwarp_path, "/fake/gdalwarp")
  expect_identical(cfg$tool_paths$ogr2ogr_exe, "/fake/ogr2ogr")

  # Dispatcher binding mirrors them.
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  bindings <- fn(cfg)
  expect_identical(bindings$python_exe, "/fake/python")
  expect_identical(bindings$gdal_polygonize_script, "/fake/gdal_polygonize.py")
  expect_identical(bindings$gdalwarp_path, "/fake/gdalwarp")
  expect_identical(bindings$ogr2ogr_exe, "/fake/ogr2ogr")
})

test_that("AS01: legacy pipeline defaults are NULL (no Olga.Viedma path leak)", {
  fn <- get("build_unburned_from_legacy_pipeline",
            envir = asNamespace("OtsuFire"))
  fmls <- formals(fn)
  expect_null(eval(fmls$python_exe))
  expect_null(eval(fmls$gdal_polygonize_script))
  expect_null(eval(fmls$gdalwarp_path))
  expect_null(eval(fmls$ogr2ogr_exe))
})

# --- AS03 ------------------------------------------------------------
# GATE 6.2 (2026-06-11): `UNB_LEGACY_RANDOM_SEED` (legacy_random_seed) binding
# was removed — it seeded ONLY the now-deleted Otsu generation-side pre-thinning.
# The drop max_s_patch cap binding is still exposed; the review/keep caps remain
# until GATE 6.4 removes them.
test_that("AS03: dispatcher exposes the drop max_s_patch cap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    options = list(
      legacy_drop_max_s_patch    = 0.10
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)
  expect_null(b$UNB_LEGACY_RANDOM_SEED)
  expect_equal(b$UNB_LEGACY_DROP_MAX_S_PATCH, 0.10)
})

# --- B2 (2026-06-05) -------------------------------------------------
# Plumb the previously-pinned UNB_* / UNB_LEGACY_* shared negative-pool
# params and the supervised run prefixes from config$options into the
# dispatcher bindings, defaulting to the orchestrator's hardcoded values.
test_that("B2: dispatcher binding defaults match the orchestrator hardcoded values", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir()
    # no options -> every key must fall back to its historical default
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)

  # shared deterministic-decision negative params
  expect_equal(b$UNB_EXCL_BUFFER_M, 500)
  expect_equal(b$UNB_N_RANDOM_CELLS, 1500)
  expect_equal(b$UNB_RANDOM_RBR_Q, 0.50)
  expect_equal(b$UNB_RANDOM_SEED, 42)
  expect_equal(b$UNB_RANDOM_PATCH_SIZE_CELLS, 3)

  # supervised run prefixes
  expect_identical(b$prefix_oof_base, "patch")
  expect_identical(b$prefix_base, "patch_certified")

  # remaining legacy Otsu params
  expect_equal(b$UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE, 0)
  expect_equal(b$UNB_LEGACY_MIN_PIXELS, 8)
  expect_equal(b$UNB_LEGACY_BUFFERS_M, 90)
  expect_equal(b$UNB_LEGACY_CORE_THR, 0.60)
  expect_equal(b$UNB_LEGACY_ALPHA_BOOST, 0.25)
  expect_equal(b$UNB_LEGACY_MIN_BASE_BOOST, 0.35)
  expect_equal(b$UNB_LEGACY_DIST_POWER, 1)
  expect_equal(b$UNB_LEGACY_KEEP_HI, 0.45)
  expect_equal(b$UNB_LEGACY_DROP_LO, 0.15)
  expect_equal(b$UNB_LEGACY_EXCL_BUFFER_M, 0)
  expect_equal(b$UNB_LEGACY_MIN_AREA_HA, 0)
})

test_that("B2: config$options overrides reach the dispatcher bindings", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    options = list(
      unb_excl_buffer_m           = 750,
      unb_n_random_cells          = 2000L,
      unb_random_rbr_q            = 0.60,
      unb_random_seed             = 7L,
      unb_random_patch_size_cells = 5L,
      prefix_oof_base             = "patchX",
      prefix_base                 = "patchX_certified",
      legacy_min_otsu_threshold_value = 1,
      legacy_min_pixels           = 12L,
      legacy_buffers_m            = 120,
      legacy_core_thr             = 0.70,
      legacy_alpha_boost          = 0.30,
      legacy_min_base_boost       = 0.40,
      legacy_dist_power           = 2,
      legacy_keep_hi              = 0.50,
      legacy_drop_lo              = 0.20,
      legacy_excl_buffer_m        = 30,
      legacy_min_area_ha          = 5
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)

  expect_equal(b$UNB_EXCL_BUFFER_M, 750)
  expect_equal(b$UNB_N_RANDOM_CELLS, 2000L)
  expect_equal(b$UNB_RANDOM_RBR_Q, 0.60)
  expect_equal(b$UNB_RANDOM_SEED, 7L)
  expect_equal(b$UNB_RANDOM_PATCH_SIZE_CELLS, 5L)
  expect_identical(b$prefix_oof_base, "patchX")
  expect_identical(b$prefix_base, "patchX_certified")
  expect_equal(b$UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE, 1)
  expect_equal(b$UNB_LEGACY_MIN_PIXELS, 12L)
  expect_equal(b$UNB_LEGACY_BUFFERS_M, 120)
  expect_equal(b$UNB_LEGACY_CORE_THR, 0.70)
  expect_equal(b$UNB_LEGACY_ALPHA_BOOST, 0.30)
  expect_equal(b$UNB_LEGACY_MIN_BASE_BOOST, 0.40)
  expect_equal(b$UNB_LEGACY_DIST_POWER, 2)
  expect_equal(b$UNB_LEGACY_KEEP_HI, 0.50)
  expect_equal(b$UNB_LEGACY_DROP_LO, 0.20)
  expect_equal(b$UNB_LEGACY_EXCL_BUFFER_M, 30)
  expect_equal(b$UNB_LEGACY_MIN_AREA_HA, 5)
})
