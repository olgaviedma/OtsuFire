# AS01 (OtsuFire 0.3.0): external-tool paths plumbed through
# config$tool_paths -> dispatcher bindings -> orchestrator ->
# build_otsu_negative_pipeline.

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
    run_label = "balanced",
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
  fn <- get("build_otsu_negative_pipeline",
            envir = asNamespace("OtsuFire"))
  fmls <- formals(fn)
  expect_null(eval(fmls$python_exe))
  expect_null(eval(fmls$gdal_polygonize_script))
  expect_null(eval(fmls$gdalwarp_path))
  expect_null(eval(fmls$ogr2ogr_exe))
})

# --- AS03 ------------------------------------------------------------
# GATE 6.2 (2026-06-11): `UNB_OTSU_NEG_RANDOM_SEED` (otsu_negative_random_seed) binding
# was removed — it seeded ONLY the now-deleted Otsu generation-side pre-thinning.
# The drop max_s_patch cap binding is still exposed; the review/keep caps remain
# until GATE 6.4 removes them.
test_that("AS03 / GATE 6.7: dispatcher exposes the FIXED drop max_s_patch cap (internal default, not user-settable)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  # GATE 6.7 (2026-06-12): otsu_negative_drop_max_s_patch is a FIXED internal
  # default (0.15), no longer a free-form option. A stale option is ignored
  # (not read); the binding always carries the internal-defaults value.
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    options = list(
      otsu_negative_drop_max_s_patch    = 0.10  # ignored now (internal default)
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)
  expect_null(b$UNB_OTSU_NEG_RANDOM_SEED)
  expect_equal(b$UNB_OTSU_NEG_DROP_MAX_S_PATCH, 0.15)
})

# --- B2 (2026-06-05) -------------------------------------------------
# Plumb the previously-pinned UNB_* / UNB_OTSU_NEG_* shared negative-pool
# params and the supervised run prefixes from config$options into the
# dispatcher bindings, defaulting to the orchestrator's hardcoded values.
test_that("B2: dispatcher binding defaults match the orchestrator hardcoded values", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
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
  expect_equal(b$UNB_OTSU_NEG_MIN_THRESHOLD_VALUE, 0)
  expect_equal(b$UNB_OTSU_NEG_MIN_PIXELS, 8)
  expect_equal(b$UNB_OTSU_NEG_BUFFERS_M, 90)
  expect_equal(b$UNB_OTSU_NEG_CORE_THR, 0.60)
  expect_equal(b$UNB_OTSU_NEG_ALPHA_BOOST, 0.25)
  expect_equal(b$UNB_OTSU_NEG_MIN_BASE_BOOST, 0.35)
  expect_equal(b$UNB_OTSU_NEG_DIST_POWER, 1)
  expect_equal(b$UNB_OTSU_NEG_KEEP_HI, 0.45)
  expect_equal(b$UNB_OTSU_NEG_DROP_LO, 0.15)
  expect_equal(b$UNB_OTSU_NEG_EXCL_BUFFER_M, 0)
  expect_equal(b$UNB_OTSU_NEG_MIN_AREA_HA, 0)
})

test_that("GATE 6.7: the PUBLIC negative_pool_params knobs + random_seed reach the dispatcher bindings (single source)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  # GATE 6.7 (2026-06-12): the user-settable negative-pool knobs are the typed
  # PUBLIC block (random n_cells / rbr_quantile, otsu candidate / reference
  # thresholds + caps) and the canonical random_seed. They flow into the
  # dispatcher bindings from cfg$negative_pool_params / train_control$seeds.
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    random_seed = 7L,
    negative_pool_params = list(
      random = list(n_cells = 2000L, rbr_quantile = 0.60),
      otsu   = list(candidate_threshold = 1, reference_threshold = 120)
    ),
    options = list(
      prefix_oof_base = "patchX",
      prefix_base     = "patchX_certified"
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)

  # PUBLIC random knobs + random_seed.
  expect_equal(b$UNB_N_RANDOM_CELLS, 2000L)
  expect_equal(b$UNB_RANDOM_RBR_Q, 0.60)
  expect_equal(b$UNB_RANDOM_SEED, 7L)
  # PUBLIC otsu thresholds.
  expect_equal(b$UNB_OTSU_NEG_THRESHOLD, 1)
  expect_equal(b$UNB_OTSU_NEG_REFERENCE_THRESHOLD, 120)
  # Output-naming prefixes (still plain options).
  expect_identical(b$prefix_oof_base, "patchX")
  expect_identical(b$prefix_base, "patchX_certified")
})

test_that("GATE 6.7: the FORMERLY free-form unb_* / otsu_negative_* options are NO LONGER read (fixed internal defaults)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_tp()
  id <- mk_tmp_gpkg_tp()
  # These knobs are now FIXED internal defaults
  # (.of_supervised_negative_pool_internal_defaults()), not user-settable. A
  # config still carrying them does NOT override the bindings: the single source
  # of truth is the internal-defaults block.
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir(),
    options = list(
      unb_excl_buffer_m                 = 750,   # ignored
      unb_random_patch_size_cells       = 5L,    # ignored
      otsu_negative_min_threshold_value = 1,     # ignored
      otsu_negative_min_pixels          = 12L,   # ignored
      otsu_negative_buffers_m           = 120,   # ignored
      otsu_negative_core_thr            = 0.70,  # ignored
      otsu_negative_alpha_boost         = 0.30,  # ignored
      otsu_negative_min_base_boost      = 0.40,  # ignored
      otsu_negative_dist_power          = 2,     # ignored
      otsu_negative_keep_hi             = 0.50,  # ignored
      otsu_negative_drop_lo             = 0.20,  # ignored
      otsu_negative_excl_buffer_m       = 30,    # ignored
      otsu_negative_min_area_ha         = 5      # ignored
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)
  ni <- get(".of_supervised_negative_pool_internal_defaults",
            envir = asNamespace("OtsuFire"))()

  # Every binding equals the FIXED internal default, NOT the stale option value.
  expect_equal(b$UNB_EXCL_BUFFER_M, ni$random_exclusion_buffer_m)
  expect_equal(b$UNB_RANDOM_PATCH_SIZE_CELLS, ni$random_patch_size_cells)
  expect_equal(b$UNB_OTSU_NEG_MIN_THRESHOLD_VALUE, ni$otsu_min_threshold_value)
  expect_equal(b$UNB_OTSU_NEG_MIN_PIXELS, ni$otsu_min_pixels)
  expect_equal(b$UNB_OTSU_NEG_BUFFERS_M, ni$otsu_buffers_m)
  expect_equal(b$UNB_OTSU_NEG_CORE_THR, ni$otsu_core_thr)
  expect_equal(b$UNB_OTSU_NEG_ALPHA_BOOST, ni$otsu_alpha_boost)
  expect_equal(b$UNB_OTSU_NEG_MIN_BASE_BOOST, ni$otsu_min_base_boost)
  expect_equal(b$UNB_OTSU_NEG_DIST_POWER, ni$otsu_dist_power)
  expect_equal(b$UNB_OTSU_NEG_KEEP_HI, ni$otsu_keep_hi)
  expect_equal(b$UNB_OTSU_NEG_DROP_LO, ni$otsu_drop_lo)
  expect_equal(b$UNB_OTSU_NEG_EXCL_BUFFER_M, ni$otsu_excl_buffer_m)
  expect_equal(b$UNB_OTSU_NEG_MIN_AREA_HA, ni$otsu_min_area_ha)
})
