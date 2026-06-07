# AS06 (OtsuFire 0.3.0): when sanitisation removes every legacy patch,
# the function must warn() and write `_LEGACY_POOL_EMPTY.txt`.

mk_rect_pool <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(matrix(
    c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
    ncol = 2, byrow = TRUE
  )))
}

# Shared fixture builder for the empty-Otsu-pool scenario (100% overlap).
.mk_empty_pool_inputs <- function() {
  legacy_geom <- sf::st_sfc(
    mk_rect_pool(0, 0, 10, 10),
    mk_rect_pool(20, 0, 30, 10),
    crs = 3035
  )
  legacy <- sf::st_sf(
    DECISION = c("drop", "drop"),
    S_PATCH_PA = c(0.05, 0.05),
    geometry = legacy_geom
  )
  internal_geom <- sf::st_sfc(
    mk_rect_pool(-5, -5, 35, 15),
    crs = 3035
  )
  internal <- sf::st_sf(geometry = internal_geom)

  # The reader defaults to layer "internal_decisions" for GPKG inputs;
  # use shapefile for the legacy patches (no layer concept) and GPKG
  # with the expected layer name for the internal-decisions layer.
  legacy_path <- tempfile(fileext = ".shp")
  internal_path <- tempfile(fileext = ".gpkg")
  sf::st_write(legacy, legacy_path, quiet = TRUE, delete_dsn = TRUE)
  sf::st_write(internal, internal_path, layer = "internal_decisions",
               quiet = TRUE, delete_dsn = TRUE)
  list(legacy_path = legacy_path, internal_path = internal_path)
}

# --- D4a ------------------------------------------------------------
test_that("D4a: empty Otsu pool ERRORS by default (no silent degradation)", {
  skip_if_not_installed("sf")

  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  inp <- .mk_empty_pool_inputs()
  out_gpkg <- tempfile(fileext = ".gpkg")

  expect_error(
    fn(
      legacy_patches_path = inp$legacy_path,
      internal_decisions_path = inp$internal_path,
      out_gpkg = out_gpkg,
      use_drop = TRUE,
      use_review = FALSE,
      use_keep = FALSE,
      drop_max_s_patch = 0.15,
      sample_n = 100L,
      sample_props = c(drop = 1.0),
      exclude_buffer_m = 50,
      # allow_empty_otsu_pool defaults to FALSE
      verbose = FALSE
    ),
    regexp = "Otsu legacy pool is EMPTY|allow_empty_otsu_pool"
  )

  # Default-FALSE path must NOT write the audit file (it aborts first).
  audit_path <- file.path(dirname(out_gpkg), "_LEGACY_POOL_EMPTY.txt")
  expect_false(file.exists(audit_path))
})

test_that("D4a: build_unburned_from_legacy_decisions defaults allow_empty_otsu_pool to FALSE", {
  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  expect_false(eval(formals(fn)$allow_empty_otsu_pool))
  fn2 <- get("build_unburned_from_legacy_pipeline",
             envir = asNamespace("OtsuFire"))
  expect_false(eval(formals(fn2)$allow_empty_otsu_pool))
})

# --- T19 ------------------------------------------------------------
test_that("T19: empty pool with opt-in -> warning + _LEGACY_POOL_EMPTY.txt", {
  skip_if_not_installed("sf")

  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  inp <- .mk_empty_pool_inputs()
  out_gpkg <- tempfile(fileext = ".gpkg")

  expect_warning(
    res <- fn(
      legacy_patches_path = inp$legacy_path,
      internal_decisions_path = inp$internal_path,
      out_gpkg = out_gpkg,
      use_drop = TRUE,
      use_review = FALSE,
      use_keep = FALSE,
      drop_max_s_patch = 0.15,
      sample_n = 100L,
      sample_props = c(drop = 1.0),
      exclude_buffer_m = 50,
      allow_empty_otsu_pool = TRUE,
      verbose = FALSE
    ),
    regexp = "all_sources degrading|legacy pool empty"
  )

  audit_path <- file.path(dirname(out_gpkg), "_LEGACY_POOL_EMPTY.txt")
  expect_true(file.exists(audit_path))
})

# --- AS12 -----------------------------------------------------------
test_that("AS12: build_unburned_from_legacy_decisions defaults use_review/use_keep to FALSE", {
  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  fmls <- formals(fn)
  expect_false(eval(fmls$use_review))
  expect_false(eval(fmls$use_keep))
  expect_true(eval(fmls$use_drop))
})
