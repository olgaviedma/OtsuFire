# 2026-06-12 deterministic GPKG overwrite hotfix.
#
# The three stage-2 scoring writers (scoring_internal_burned_area,
# scoring_reference_burned_area, scoring_validation_burned_area) used to replace
# their first GPKG layer with a bare `unlink()` + `st_write(... )`. On Windows a
# bare unlink() can silently fail when a stale GPKG from a previous failed run is
# still present, leaving st_write to hit an existing dataset ("Creation
# failed."). The fix propagates `delete_dsn = isTRUE(overwrite)` to the first
# st_write so GDAL drops the datasource robustly.
#
# This is a pure I/O robustness change: it must NOT alter geometries, attributes
# or feature counts, and overwrite=FALSE must still neither delete nor overwrite
# an existing datasource.

scoring_internal_internal <- function() {
  get("scoring_internal_burned_area", envir = asNamespace("OtsuFire"))
}

# Two unit squares (EPSG:3035). stage2 == stage1 => every stage-2 polygon
# intersects a stage-1 seed and is flagged "keep", so the writer is reached with
# burnable_corine = NULL (no raster needed) and no "nohit_stage1" polygons.
mk_two_squares <- function(origin = 0) {
  geom <- lapply(0:1, function(i) {
    x0 <- origin + i * 10
    sf::st_polygon(list(rbind(
      c(x0, 0), c(x0 + 5, 0), c(x0 + 5, 5), c(x0, 5), c(x0, 0)
    )))
  })
  sf::st_sf(
    pid = 1:2,
    geometry = sf::st_sfc(geom, crs = 3035)
  )
}

run_scoring_internal <- function(out_dir, prefix, overwrite) {
  fn <- scoring_internal_internal()
  s1 <- mk_two_squares()
  s2 <- mk_two_squares()
  suppressWarnings(suppressMessages(fn(
    polys_stage2    = s2,
    polys_stage1    = s1,
    burnable_corine = NULL,
    save_outputs    = TRUE,
    out_dir         = out_dir,
    prefix          = prefix,
    overwrite       = overwrite,
    quiet           = TRUE
  )))
}

test_that("all three stage-2 writers propagate delete_dsn = isTRUE(overwrite)", {
  ns <- asNamespace("OtsuFire")
  for (fnm in c("scoring_internal_burned_area",
                "scoring_reference_burned_area",
                "scoring_validation_burned_area")) {
    src <- paste(deparse(body(get(fnm, envir = ns))), collapse = "\n")
    expect_true(
      grepl("delete_dsn = isTRUE(overwrite)", src, fixed = TRUE),
      info = paste0(fnm, " must propagate delete_dsn = isTRUE(overwrite)")
    )
    # The bare first-layer write without a datasource-deletion guard must be gone.
    expect_false(
      grepl('st_write(out, gpkg_path, layer = "stage2_flagged", quiet = quiet)',
            src, fixed = TRUE),
      info = paste0(fnm, " must not use a bare first-layer st_write")
    )
  }
})

test_that("overwrite=TRUE replaces a stale GPKG robustly (delete_dsn)", {
  skip_if_not_installed("sf")
  out_dir <- tempfile("deldsn_ow_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  prefix <- "T"
  gpkg <- file.path(out_dir, paste0(prefix, "_phase1.gpkg"))

  # Simulate a STALE GPKG left over from a previous (failed) run: a valid GPKG
  # at the exact target path carrying an unrelated layer.
  stale <- sf::st_sf(z = 99L, geometry = sf::st_sfc(
    sf::st_point(c(0, 0)), crs = 3035))
  suppressWarnings(sf::st_write(stale, gpkg, layer = "stale_layer", quiet = TRUE))
  expect_true(file.exists(gpkg))

  # overwrite=TRUE must succeed and produce the canonical layer set, with the
  # stale layer gone (datasource was dropped).
  expect_no_error(run_scoring_internal(out_dir, prefix, overwrite = TRUE))
  layers <- sf::st_layers(gpkg)$name
  expect_true("stage2_flagged" %in% layers)
  expect_false("stale_layer" %in% layers)

  flagged <- sf::st_read(gpkg, layer = "stage2_flagged", quiet = TRUE)
  expect_equal(nrow(flagged), 2L)
  expect_true(all(flagged$flag_internal == "keep"))
})

test_that("two overwrite=TRUE runs are byte-stable in counts/attrs/geometry", {
  skip_if_not_installed("sf")
  dir_a <- tempfile("deldsn_a_"); dir_b <- tempfile("deldsn_b_")
  dir.create(dir_a, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_b, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(c(dir_a, dir_b), recursive = TRUE, force = TRUE), add = TRUE)

  run_scoring_internal(dir_a, "T", overwrite = TRUE)
  # Second run into a directory that already holds the GPKG (the stale-path case)
  run_scoring_internal(dir_b, "T", overwrite = TRUE)
  run_scoring_internal(dir_b, "T", overwrite = TRUE)

  fa <- sf::st_read(file.path(dir_a, "T_phase1.gpkg"),
                    layer = "stage2_flagged", quiet = TRUE)
  fb <- sf::st_read(file.path(dir_b, "T_phase1.gpkg"),
                    layer = "stage2_flagged", quiet = TRUE)

  expect_equal(nrow(fa), nrow(fb))
  expect_identical(
    sf::st_drop_geometry(fa)[order(fa$source_poly_id), , drop = FALSE],
    sf::st_drop_geometry(fb)[order(fb$source_poly_id), , drop = FALSE]
  )
  # Geometry identical to the millimetre.
  expect_true(all(sf::st_equals(
    fa[order(fa$source_poly_id), ], fb[order(fb$source_poly_id), ],
    sparse = FALSE)[cbind(seq_len(nrow(fa)), seq_len(nrow(fb)))]))
})

test_that("overwrite=FALSE neither deletes nor overwrites an existing datasource", {
  skip_if_not_installed("sf")
  out_dir <- tempfile("deldsn_noow_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  prefix <- "T"
  gpkg <- file.path(out_dir, paste0(prefix, "_phase1.gpkg"))

  # Pre-existing datasource with a sentinel layer whose name does NOT collide
  # with the writer's layers.
  sentinel <- sf::st_sf(keep = 1L, geometry = sf::st_sfc(
    sf::st_point(c(1, 1)), crs = 3035))
  suppressWarnings(sf::st_write(sentinel, gpkg, layer = "sentinel_keep",
                                quiet = TRUE))
  before_mtime <- file.info(gpkg)$mtime

  # overwrite=FALSE => delete_dsn=FALSE: the datasource must survive intact
  # (the sentinel layer is preserved). New layers may be appended, but nothing
  # pre-existing is deleted or overwritten.
  suppressWarnings(run_scoring_internal(out_dir, prefix, overwrite = FALSE))

  expect_true(file.exists(gpkg))
  expect_true("sentinel_keep" %in% sf::st_layers(gpkg)$name)
  kept <- sf::st_read(gpkg, layer = "sentinel_keep", quiet = TRUE)
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$keep, 1L)
})
