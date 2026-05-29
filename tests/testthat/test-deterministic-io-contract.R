# Tests for honest I/O contracts in the deterministic phase.
# Currently covers P1-DET-02: merge_aoi_shapefiles() must refuse to
# replace an existing output when overwrite = FALSE and must replace it
# when overwrite = TRUE.

mk_tiny_aoi_polys <- function(dir, n_files = 2L) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(n_files)
  for (i in seq_len(n_files)) {
    poly <- sf::st_sf(
      id = i,
      geometry = sf::st_sfc(
        sf::st_polygon(list(rbind(
          c(0 + i, 0 + i), c(1 + i, 0 + i),
          c(1 + i, 1 + i), c(0 + i, 1 + i),
          c(0 + i, 0 + i)
        ))),
        crs = 3035
      )
    )
    p <- file.path(dir, sprintf("BA_AOI_%d.shp", i))
    sf::st_write(poly, p, quiet = TRUE, delete_layer = TRUE)
    paths[i] <- p
  }
  paths
}

test_that("merge_aoi_shapefiles refuses to overwrite when overwrite=FALSE", {
  base <- tempfile("merge_io_")
  mk_tiny_aoi_polys(base, n_files = 2L)

  out_path <- file.path(base, "merged.gpkg")
  # First call writes the destination.
  merge_aoi_shapefiles(
    base_output_dir = base,
    pattern = "^BA_AOI_\\d+\\.shp$",
    out_path = out_path,
    dissolve = FALSE,
    verbose = FALSE,
    overwrite = FALSE
  )
  expect_true(file.exists(out_path))
  first_mtime <- file.info(out_path)$mtime

  # Second call must refuse to clobber.
  err <- tryCatch(
    merge_aoi_shapefiles(
      base_output_dir = base,
      pattern = "^BA_AOI_\\d+\\.shp$",
      out_path = out_path,
      dissolve = FALSE,
      verbose = FALSE,
      overwrite = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "already exists", fixed = TRUE)
  expect_match(err, "overwrite = FALSE", fixed = TRUE)
  # The destination must NOT have been replaced.
  expect_identical(file.info(out_path)$mtime, first_mtime)
})

test_that("merge_aoi_shapefiles replaces destination when overwrite=TRUE", {
  base <- tempfile("merge_io_ok_")
  mk_tiny_aoi_polys(base, n_files = 2L)
  out_path <- file.path(base, "merged.gpkg")

  merge_aoi_shapefiles(
    base_output_dir = base,
    pattern = "^BA_AOI_\\d+\\.shp$",
    out_path = out_path,
    dissolve = FALSE,
    verbose = FALSE,
    overwrite = FALSE
  )
  expect_true(file.exists(out_path))

  # Touch the file in the past so we can observe replacement.
  Sys.setFileTime(out_path, Sys.time() - 60)
  before_mtime <- file.info(out_path)$mtime

  # overwrite=TRUE replaces the destination without error.
  expect_no_error(
    merge_aoi_shapefiles(
      base_output_dir = base,
      pattern = "^BA_AOI_\\d+\\.shp$",
      out_path = out_path,
      dissolve = FALSE,
      verbose = FALSE,
      overwrite = TRUE
    )
  )
  expect_true(file.exists(out_path))
  expect_true(file.info(out_path)$mtime >= before_mtime)
})

test_that("merge_aoi_shapefiles signature exposes overwrite parameter", {
  expect_true("overwrite" %in% names(formals(merge_aoi_shapefiles)))
  expect_identical(formals(merge_aoi_shapefiles)$overwrite, FALSE)
})
