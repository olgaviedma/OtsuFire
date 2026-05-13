mk_tmp_tif_s <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_gpkg_s <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                              c(0,1), c(0,0)))),
                    crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}

test_that("supervised config returns expected S3 and required fields", {
  ci <- mk_tmp_tif_s()
  id <- mk_tmp_gpkg_s()
  cfg <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = id,
    change_index = ci,
    target_year = 2025,
    output_dir = tempdir()
  )
  expect_s3_class(cfg, "otsufire_supervised_burned_config")
  expect_identical(cfg$scenario, "balanced")
  expect_identical(cfg$target_year, 2025L)
  expect_identical(cfg$min_burned_pool_n, 5L)
  expect_type(cfg$output_routes, "list")
  expect_true(all(c("pools_gpkg", "final_map_gpkg", "burned_like_gpkg",
                    "timing_csv") %in% names(cfg$output_routes)))
  expect_identical(cfg$negative_pool_policy, "all_sources")
})

test_that("supervised config respects scenario-based negative-pool defaults and overrides", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()

  cfg_b <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_b$negative_pool_policy, "all_sources")

  cfg_o <- build_supervised_burned_config(
    scenario = "original", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_o$negative_pool_policy, "deterministic_direct")

  cfg_over <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L,
    options = list(negative_pool_policy = "deterministic_direct")
  )
  expect_identical(cfg_over$negative_pool_policy, "deterministic_direct")
})

test_that("supervised config rejects invalid inputs", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                    internal_decisions = id,
                                    change_index = ci),
    regexp = "target_year"
  )
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                    change_index = ci, target_year = 2025L),
    regexp = "internal_decisions"
  )
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                    internal_decisions = id,
                                    target_year = 2025L),
    regexp = "change_index"
  )
  expect_error(
    build_supervised_burned_config(scenario = "ultra",
                                    internal_decisions = id,
                                    change_index = ci, target_year = 2025L),
    regexp = "should be one of"
  )
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                    internal_decisions = id,
                                    change_index = ci, target_year = 2025L,
                                    min_burned_pool_n = -1),
    regexp = "min_burned_pool_n"
  )
})

test_that("supervised config deprecates otsu_unburned_generation", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_warning(
    cfg <- build_supervised_burned_config(
      scenario = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(negative_pool_policy = "otsu_unburned_generation")
    ),
    regexp = "deprecated"
  )
  expect_identical(cfg$negative_pool_policy, "all_sources")
})
