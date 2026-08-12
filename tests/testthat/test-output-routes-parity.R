# Bug 1 (OtsuFire 0.3.0): config output_routes match the runtime
# convention.

mk_tmp_tif_routes <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_gpkg_routes <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                                c(0,1), c(0,0)))),
                    crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}

# --- T11 ------------------------------------------------------------
test_that("T11: cfg$output_routes folder names match runtime convention", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_routes()
  id <- mk_tmp_gpkg_routes()
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir()
  )

  # Final-model directory must end with the runtime convention.
  expect_match(cfg$output_routes$final_model_dir, "07_FINAL_MODEL_V2$")
  expect_match(cfg$output_routes$scored_dir,      "08_SCORED$")
  expect_match(cfg$output_routes$final_map_dir,   "09_FINAL_MAP$")
  expect_match(cfg$output_routes$consistency_dir, "11_CONSISTENCY_CHECKS$")

  # Final-map gpkg path lives under 09_FINAL_MAP, not the legacy 08.
  expect_match(cfg$output_routes$final_map_gpkg, "09_FINAL_MAP")
  expect_match(cfg$output_routes$burned_like_gpkg, "09_FINAL_MAP")
  expect_match(cfg$output_routes$final_model_rds, "07_FINAL_MODEL_V2")
})

# --- Change 3 (2026-06-05): output location is config-driven --------
# The config's output_routes$base (driven by output_dir + run_name) is the
# single source of truth for WHERE supervised outputs are written. For
# Natalia's default config (output_dir = <data_base>/Results and
# run_name = result_name), this must resolve to the SAME path the engine
# previously reconstructed from data_base + result_name, so outputs do NOT
# move for the existing setup.
test_that("output_routes$base coincides with the legacy engine reconstruction for the default config", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_routes()
  id <- mk_tmp_gpkg_routes()

  data_base   <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  result_name <- "Min_Min"
  target_year <- 2025L
  scenario    <- "balanced"

  cfg <- build_supervised_burned_config(
    run_label = scenario,
    internal_decisions = id,
    change_index = ci,
    target_year = target_year,
    # Natalia's default: output_dir = <data_base>/Results, run_name = result_name
    output_dir = file.path(data_base, "Results"),
    run_name = result_name
  )

  legacy_reconstruction <- file.path(
    data_base, "Results", as.character(target_year), result_name,
    "SUPERVISED", scenario
  )

  expect_identical(
    normalizePath(cfg$output_routes$base, winslash = "/", mustWork = FALSE),
    normalizePath(legacy_reconstruction, winslash = "/", mustWork = FALSE)
  )
})

test_that("config output_routes do NOT advertise the legacy folder names", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  ci <- mk_tmp_tif_routes()
  id <- mk_tmp_gpkg_routes()
  cfg <- build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir()
  )

  # No path should contain the pre-0.3.0 legacy folder names.
  paths <- unlist(cfg$output_routes, use.names = FALSE)
  expect_false(any(grepl("06_FINAL_MODEL(?!_V2)", paths, perl = TRUE)),
               info = "06_FINAL_MODEL must be replaced by 07_FINAL_MODEL_V2")
  expect_false(any(grepl("/07_SCORED/", paths)),
               info = "07_SCORED must become 08_SCORED")
  expect_false(any(grepl("/08_FINAL_MAP/", paths)),
               info = "08_FINAL_MAP must become 09_FINAL_MAP")
  expect_false(any(grepl("/09_CONSISTENCY/", paths)),
               info = "09_CONSISTENCY must become 11_CONSISTENCY_CHECKS")
})
