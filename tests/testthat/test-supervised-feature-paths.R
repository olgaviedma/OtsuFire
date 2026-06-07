# Gate 1B PIECE 2 (2026-06-07): every FEATURE/CANDIDATE input path is a
# cfg$inputs field (no hardcoded absolute paths, no ecoregion_shapefile-style
# hardcode). PIECE 2 ONLY makes the paths configurable + consumed from cfg; it
# does NOT change fallback BEHAVIOUR (fail-fast on silent fallback is PIECE 5)
# and does NOT remove the ecoregion branch (PIECE 4).

ns <- asNamespace("OtsuFire")

mk_tmp_tif_fp <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}
mk_tmp_gpkg_fp <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}
mk_tmp_shp_fp <- function() {
  f <- tempfile(fileext = ".shp")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  suppressWarnings(sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f,
                                quiet = TRUE, delete_dsn = TRUE))
  f
}

# Build a complete fake data_base + composite_base so build_supervised_burned_config's
# config-time fail-fast loop (peninsula/topo/corine/burnable existence) passes,
# including the Ecoregion convention file.
mk_full_data_base_fp <- function(target_year, corine_year) {
  db <- file.path(tempdir(), paste0("ecodb_", as.integer(stats::runif(1, 1, 1e8))))
  for (d in c("Borders", "Topography", "Corine_Masks", "Ecoregion"))
    dir.create(file.path(db, d), recursive = TRUE, showWarnings = FALSE)
  cb <- file.path(db, "Composites")
  dir.create(file.path(cb, "Autumn"), recursive = TRUE, showWarnings = FALSE)
  r1 <- terra::rast(ncol = 4, nrow = 4, vals = 1:16)
  terra::writeRaster(c(r1, r1), file.path(db, "Topography", "elevation_slope.tif"),
                     overwrite = TRUE)
  terra::writeRaster(r1, file.path(db, "Corine_Masks",
                                   paste0("CLC_", corine_year, "_peninsula.tif")),
                     overwrite = TRUE)
  terra::writeRaster(r1, file.path(db, "Corine_Masks",
                                   paste0("burneable_mask_binary_corine_",
                                          corine_year, "_ETRS89.tif")),
                     overwrite = TRUE)
  terra::writeRaster(r1, file.path(cb, "Autumn",
                                   paste0("mean_mean_", target_year, "_mosaic.tif")),
                     overwrite = TRUE)
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  suppressWarnings(sf::st_write(sf::st_sf(id = 1L, geometry = sfc),
                                file.path(db, "Borders", "Iberian_peninsula.shp"),
                                quiet = TRUE, delete_dsn = TRUE))
  list(data_base = db, composite_base = cb)
}

# --------------------------------------------------------------------------
# 1) ecoregion_shapefile is a declared, validated cfg$inputs field with a
#    convention default and a user-override that is honoured.
# --------------------------------------------------------------------------
test_that("PIECE 2: ecoregion_shapefile defaults to the convention path", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_fp(); ci <- mk_tmp_tif_fp()
  paths <- mk_full_data_base_fp(2012L, "2012")
  # Create the Ecoregion convention file so the (path) spec is honoured as-is.
  eco_conv <- file.path(paths$data_base, "Ecoregion", "ecoregiones_olson.shp")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  suppressWarnings(sf::st_write(sf::st_sf(id = 1L, geometry = sfc), eco_conv,
                                quiet = TRUE, delete_dsn = TRUE))

  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2012L,
    options = list(data_base = paths$data_base,
                   composite_base = paths$composite_base)
  )
  expect_true("ecoregion_shapefile" %in% names(cfg$inputs))
  expect_identical(
    cfg$inputs$ecoregion_shapefile$path,
    normalizePath(eco_conv, winslash = "/", mustWork = FALSE)
  )
})

test_that("PIECE 2: a user-supplied ecoregion_shapefile overrides the convention", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_fp(); ci <- mk_tmp_tif_fp(); eco <- mk_tmp_shp_fp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2012L, ecoregion_shapefile = eco
  )
  expect_identical(cfg$inputs$ecoregion_shapefile$path,
                   normalizePath(eco, winslash = "/", mustWork = FALSE))
  # The shared accessor returns it (single source of truth).
  resolve <- get(".of_sup_input_path", envir = ns)
  expect_identical(resolve(cfg, "ecoregion_shapefile"),
                   normalizePath(eco, winslash = "/", mustWork = FALSE))
})

test_that("PIECE 2: data_base-less config leaves ecoregion_shapefile deferred (NULL)", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_fp(); ci <- mk_tmp_tif_fp()
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2012L  # no data_base -> convention cannot be built
  )
  expect_true("ecoregion_shapefile" %in% names(cfg$inputs))
  expect_null(cfg$inputs$ecoregion_shapefile)
})

# --------------------------------------------------------------------------
# 2) The legacy Otsu builder exposes ecoregion_shapefile_path and CONSUMES it
#    (no unconditional data_base/Ecoregion hardcode). The branch is retained.
# --------------------------------------------------------------------------
test_that("PIECE 2: legacy builder exposes + consumes ecoregion_shapefile_path", {
  fn <- get("build_unburned_from_legacy_pipeline", envir = ns)
  fmls <- formals(fn)
  expect_true("ecoregion_shapefile_path" %in% names(fmls))
  expect_null(eval(fmls$ecoregion_shapefile_path))  # NULL default (convention)

  src <- paste(deparse(fn), collapse = "\n")
  # The hardcode is now gated behind the configurable path (no unconditional
  # ecoregion_shapefile <- file.path(data_base, "Ecoregion", ...) assignment).
  expect_match(src, "ecoregion_shapefile_path")
  # The convention path is still REACHABLE as a fallback (branch retained for PIECE 4).
  expect_match(src, "ecoregiones_olson", fixed = TRUE)
  # The assignment must consult the configurable path first (it appears before
  # the convention file.path in the resolution expression).
  i_cfg  <- regexpr("ecoregion_shapefile_path", src)
  i_conv <- regexpr("ecoregiones_olson", src)
  expect_gt(i_cfg, 0); expect_gt(i_conv, 0)
  expect_lt(i_cfg, i_conv)
})

# --------------------------------------------------------------------------
# 3) The pool builder threads the cfg ecoregion path to the legacy builder.
# --------------------------------------------------------------------------
test_that("PIECE 2: pool builder threads cfg ecoregion_shapefile to the legacy builder", {
  src <- paste(deparse(get("build_supervised_training_pools", envir = ns)),
               collapse = "\n")
  expect_match(src, 'sup_input_path\\([^,]*,\\s*"ecoregion_shapefile"\\)|cfg_input_path\\("ecoregion_shapefile"\\)')
  expect_match(src, "ecoregion_shapefile_path")
})

# --------------------------------------------------------------------------
# 4) The feature extractor's standalone branch reads topo + corine from cfg.
# --------------------------------------------------------------------------
test_that("PIECE 2: feature extractor consumes topo + corine_raster from cfg", {
  src <- paste(deparse(get("extract_supervised_features", envir = ns)),
               collapse = "\n")
  expect_match(src, 'cfg_input_path\\("topo"\\)|sup_input_path\\([^,]*,\\s*"topo"\\)')
  expect_match(src, 'cfg_input_path\\("corine_raster"\\)|sup_input_path\\([^,]*,\\s*"corine_raster"\\)')
})
