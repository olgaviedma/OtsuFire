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
  # `tempfile()` returns a fresh, non-existent path, so `delete_dsn = TRUE`
  # has nothing to delete; for a `.shp` it makes GDAL probe the missing file
  # and emit a spurious "GDAL Error 1: ... does not appear to be a file"
  # message. Omit it (the GPKG keeps it harmlessly, but drop it for symmetry).
  sf::st_write(legacy, legacy_path, quiet = TRUE)
  sf::st_write(internal, internal_path, layer = "internal_decisions",
               quiet = TRUE)
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
      drop_max_s_patch = 0.15,
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
      drop_max_s_patch = 0.15,
      exclude_buffer_m = 50,
      allow_empty_otsu_pool = TRUE,
      verbose = FALSE
    ),
    regexp = "all_sources degrading|legacy pool empty"
  )

  audit_path <- file.path(dirname(out_gpkg), "_LEGACY_POOL_EMPTY.txt")
  expect_true(file.exists(audit_path))
})

# --- GATE 6.4 (2026-06-11): supersedes AS12 -------------------------
# Otsu review/keep are NEVER negatives, so the use_review / use_keep /
# review_max_s_patch / keep_max_s_patch parameters were removed entirely. Only
# the live `use_drop` path remains (effectively always TRUE).
test_that("GATE 6.4: build_unburned_from_legacy_decisions dropped the dead review/keep params; use_drop kept", {
  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  fmls <- names(formals(fn))
  expect_false("use_review" %in% fmls)
  expect_false("use_keep" %in% fmls)
  expect_false("review_max_s_patch" %in% fmls)
  expect_false("keep_max_s_patch" %in% fmls)
  expect_true("use_drop" %in% fmls)
  expect_true(eval(formals(fn)$use_drop))
})
