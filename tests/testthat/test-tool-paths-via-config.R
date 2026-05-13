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
test_that("AS03: dispatcher exposes UNB_LEGACY_RANDOM_SEED + 3 max_s_patch caps", {
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
      legacy_random_seed         = 99L,
      legacy_drop_max_s_patch    = 0.10,
      legacy_review_max_s_patch  = 0.40,
      legacy_keep_max_s_patch    = 0.65
    )
  )
  fn <- get(".of_supervised_engine_bindings", envir = asNamespace("OtsuFire"))
  b <- fn(cfg)
  expect_equal(b$UNB_LEGACY_RANDOM_SEED, 99L)
  expect_equal(b$UNB_LEGACY_DROP_MAX_S_PATCH, 0.10)
  expect_equal(b$UNB_LEGACY_REVIEW_MAX_S_PATCH, 0.40)
  expect_equal(b$UNB_LEGACY_KEEP_MAX_S_PATCH, 0.65)
})
