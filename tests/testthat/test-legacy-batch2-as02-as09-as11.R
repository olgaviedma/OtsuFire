# OPERATIONAL BATCH 2 (2026-06-06): low-level tests for the legacy
# Otsu-unburned stack robustness fixes.
#   AS02 - parameter fingerprint for honest `reuse_existing`.
#   AS09 - single canonical area column (no duplicate AREA_HA / AREA_HA_1).
#   AS11 - coverage_by_patch_raster() does not mutate its ref_raster in place.

.mk_rect_b2 <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(matrix(
    c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
    ncol = 2, byrow = TRUE
  )))
}

# --- AS02 -----------------------------------------------------------
test_that("AS02: fingerprint is deterministic and param-sensitive", {
  fp <- get("legacy_param_fingerprint_unb_legacy",
            envir = asNamespace("OtsuFire"))

  p1 <- list(target_year = 2017L, otsu_threshold = 0, buffers_m = 90,
             sample_props = c(drop = 0.7, review = 0.25, keep = 0.05),
             internal_decisions_path = "a/b.gpkg")
  # Same params, but sample_props supplied in a different element order:
  # named-vector sort makes the fingerprint order-stable.
  p2 <- list(target_year = 2017L, otsu_threshold = 0, buffers_m = 90,
             sample_props = c(keep = 0.05, drop = 0.7, review = 0.25),
             internal_decisions_path = "a/b.gpkg")
  # A genuinely different param (core_thr changed via buffers_m here).
  p3 <- modifyList(p1, list(buffers_m = 150))

  a <- fp(p1); b <- fp(p2); c <- fp(p3)

  expect_identical(a$checksum, b$checksum)
  expect_identical(a$text, b$text)
  expect_false(identical(a$checksum, c$checksum))
  expect_match(a$checksum, "^[0-9]{9}$")
})

# --- AS09 -----------------------------------------------------------
test_that("AS09: legacy pool keeps a single canonical area_ha column", {
  skip_if_not_installed("sf")

  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))

  # Two non-overlapping legacy patches carrying ESRI-truncated area dupes.
  geom <- sf::st_sfc(
    .mk_rect_b2(0, 0, 100, 100),
    .mk_rect_b2(500, 500, 700, 700),
    crs = 3035
  )
  legacy <- sf::st_sf(
    DECISION   = c("drop", "drop"),
    S_PATCH_PA = c(0.05, 0.05),
    AREA_HA    = c(-999, -999),   # stale truncated copy (wrong on purpose)
    AREA_HA_1  = c(-999, -999),   # make.unique collision copy
    geometry   = geom
  )
  internal <- sf::st_sf(
    geometry = sf::st_sfc(.mk_rect_b2(-100, -100, -50, -50), crs = 3035)
  )

  legacy_path   <- tempfile(fileext = ".shp")
  internal_path <- tempfile(fileext = ".gpkg")
  # Fresh `tempfile()` paths do not exist yet, so `delete_dsn = TRUE` has
  # nothing to delete; for a `.shp` it makes GDAL probe the missing file and
  # emit a spurious "GDAL Error 1: ... does not appear to be a file" message.
  sf::st_write(legacy, legacy_path, quiet = TRUE)
  sf::st_write(internal, internal_path, layer = "internal_decisions",
               quiet = TRUE)

  res <- fn(
    legacy_patches_path     = legacy_path,
    internal_decisions_path = internal_path,
    out_gpkg                = NULL,
    use_drop                = TRUE,
    exclude_buffer_m        = 0,
    verbose                 = FALSE
  )

  pool <- res$legacy_unburned_pool
  nm   <- names(pool)
  # Exactly one area column, lowercase, no truncated leftovers.
  expect_true("area_ha" %in% nm)
  expect_false("AREA_HA" %in% nm)
  expect_false("AREA_HA_1" %in% nm)
  # Value is geometry-derived (100m x 100m = 1 ha), NOT the stale -999.
  geom_ha <- as.numeric(sf::st_area(pool)) / 1e4
  expect_equal(pool$area_ha, geom_ha)
  expect_true(all(pool$area_ha > 0))
})

# --- AS11 -----------------------------------------------------------
test_that("AS11: coverage_by_patch_raster does not mutate ref_raster", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  fn <- get("coverage_by_patch_raster", envir = asNamespace("OtsuFire"))

  # Reference raster with NA and non-1 values so the 0/1 normalisation is
  # observable; snapshot its values before the call.
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                   ymin = 0, ymax = 4, crs = "EPSG:3035")
  terra::values(r) <- c(NA, 5, 0, 3,
                        NA, 5, 0, 3,
                        2, 2, 0, 0,
                        2, 2, 0, 0)
  before <- terra::values(r)[, 1]

  patches <- sf::st_sf(
    patch_id = 1L,
    geometry = sf::st_sfc(.mk_rect_b2(0, 0, 2, 2), crs = 3035)
  )

  out <- fn(patches = patches, ref_raster = r, ref_names = "ref",
            out_path = NULL)

  after <- terra::values(r)[, 1]
  # Caller's raster must be untouched (NAs still NA, 5 still 5, etc.).
  expect_equal(after, before)
  # And the coverage column was actually produced.
  expect_true("ref" %in% names(out))
})
