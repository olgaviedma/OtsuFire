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

test_that("supervised config is always all_sources (knob removed, honest guard)", {
  # §N+26 (2026-06-05): `negative_pool_policy` was FULLY REMOVED as a user knob.
  # The supervised pipeline is always `all_sources`:
  #   * a config built with no policy option -> cfg$negative_pool_policy ==
  #     "all_sources" (informational constant), regardless of scenario;
  #   * passing the historical canonical `negative_pool_policy = "all_sources"`
  #     is still accepted as a no-op (existing scripts don't break);
  #   * passing ANY other value (incl. the retired "deterministic_direct" and
  #     the old "otsu_unburned_generation" alias) now raises the honest guard
  #     error instead of being silently ignored / redirected.
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()

  # No policy option: always all_sources, all scenarios.
  cfg_b <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_b$negative_pool_policy, "all_sources")

  cfg_o <- build_supervised_burned_config(
    scenario = "original", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_o$negative_pool_policy, "all_sources")

  # Historical canonical value accepted as a no-op (no error, no warning).
  cfg_over <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L,
    options = list(negative_pool_policy = "all_sources")
  )
  expect_identical(cfg_over$negative_pool_policy, "all_sources")

  # Any stale value errors via the honest guard.
  expect_error(
    build_supervised_burned_config(
      scenario = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(negative_pool_policy = "deterministic_direct")
    ),
    regexp = "was removed in 2026-06"
  )
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

# ---------------------------------------------------------------------------
# §N+25 (2026-06-05): supervised-RUN inputs wired to the config.
# ---------------------------------------------------------------------------

# Build a fake data_base tree on disk so the convention paths resolve to real
# files, mirroring the package's path conventions for the RUN inputs.
mk_fake_data_base <- function(target_year, corine_year) {
  db <- file.path(tempdir(), paste0("db_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(file.path(db, "Borders"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(db, "Topography"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(db, "Corine_Masks"), recursive = TRUE, showWarnings = FALSE)
  cb <- file.path(db, "Composites"); dir.create(file.path(cb, "Autumn"),
                                                recursive = TRUE, showWarnings = FALSE)

  r1 <- terra::rast(ncol = 4, nrow = 4, vals = 1:16)
  topo <- c(r1, r1)
  terra::writeRaster(topo, file.path(db, "Topography", "elevation_slope.tif"),
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
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                              c(0,1), c(0,0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc),
               file.path(db, "Borders", "Iberian_peninsula.shp"),
               quiet = TRUE, delete_dsn = TRUE)
  list(data_base = db, composite_base = cb)
}

test_that("§N+25: RUN inputs default to the convention paths", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ty <- 2012L; cy <- "2012"
  paths <- mk_fake_data_base(ty, cy)

  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = ty,
    options = list(data_base = paths$data_base,
                   composite_base = paths$composite_base)
  )

  # The new RUN inputs are present in cfg$inputs.
  expect_true(all(c("peninsula_shapefile", "topo", "corine_raster",
                    "burnable_mask", "delayed_change_index") %in%
                  names(cfg$inputs)))

  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  expect_identical(cfg$inputs$topo$path,
                   norm(file.path(paths$data_base, "Topography",
                                  "elevation_slope.tif")))
  expect_identical(cfg$inputs$corine_raster$path,
                   norm(file.path(paths$data_base, "Corine_Masks",
                                  paste0("CLC_", cy, "_peninsula.tif"))))
  expect_identical(cfg$inputs$burnable_mask$path,
                   norm(file.path(paths$data_base, "Corine_Masks",
                                  paste0("burneable_mask_binary_corine_", cy,
                                         "_ETRS89.tif"))))
  expect_identical(cfg$inputs$peninsula_shapefile$path,
                   norm(file.path(paths$data_base, "Borders",
                                  "Iberian_peninsula.shp")))
  expect_identical(cfg$inputs$delayed_change_index$path,
                   norm(file.path(paths$composite_base, "Autumn",
                                  paste0("mean_mean_", ty, "_mosaic.tif"))))
})

test_that("§N+25: a user-supplied RUN-input path overrides the convention", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ty <- 2012L; cy <- "2012"
  paths <- mk_fake_data_base(ty, cy)
  my_corine <- mk_tmp_tif_s()

  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = ty, corine_raster = my_corine,
    options = list(data_base = paths$data_base,
                   composite_base = paths$composite_base)
  )
  expect_identical(cfg$inputs$corine_raster$path,
                   normalizePath(my_corine, winslash = "/", mustWork = FALSE))
})

test_that("§N+25: config construction FAILS FAST on a missing RUN input", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  missing_topo <- file.path(tempdir(), "does_not_exist_topo.tif")
  expect_error(
    build_supervised_burned_config(
      scenario = "balanced", internal_decisions = id, change_index = ci,
      target_year = 2012L, topo = missing_topo
    ),
    regexp = "topo.*does not exist|does not exist.*topo"
  )
})

test_that("§N+25: VALIDATION-only inputs are NOT part of cfg$inputs", {
  # strata / study-area mask / EFFIS reference belong to the separate
  # validate_fire_maps() call, not the supervised run config.
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2012L
  )
  expect_false(any(c("strata_tif", "strata_lut", "mask_tif", "mask_shp",
                     "ref_tif", "ref_shp") %in% names(cfg$inputs)))
})

test_that("§N+25: a data_base-less config leaves RUN inputs deferred (NULL)", {
  # Without options$data_base the convention cannot be built; the run inputs
  # stay NULL (deferred to runtime) and construction does not fail.
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2012L
  )
  expect_null(cfg$inputs$topo)
  expect_null(cfg$inputs$corine_raster)
  expect_null(cfg$inputs$burnable_mask)
  expect_null(cfg$inputs$peninsula_shapefile)
})

test_that("§N+25: the supervised run no longer builds/validates EFFIS/strata/study-area-mask files", {
  # The validation-only inputs were previously CONSTRUCTED and stopifnot-
  # validated in the orchestrator's COMMON INPUTS block but never consumed by
  # the run. They are removed from the run start (they belong to the separate
  # validate_fire_maps() call). Assert the run engine no longer references
  # them, so a run cannot fail/require those files. (Static introspection of
  # the in-package orchestrator function body keeps this cheap and avoids the
  # ~60-min pipeline.)
  ns <- asNamespace("OtsuFire")
  body_txt <- paste(
    deparse(body(get("run_supervised_pipeline", envir = ns))),
    collapse = "\n"
  )
  for (tok in c("strata_tif_path", "strata_lut_csv", "mask_tif_path",
                "mask_shp_path", "ref_tif_path", "ref_shp_path",
                "Effis_CA_", "Mask_StudyArea", "Validation_fires_burneable")) {
    expect_false(grepl(tok, body_txt, fixed = TRUE),
                 info = paste("orchestrator still references", tok))
  }
})

test_that("supervised config rejects the removed otsu_unburned_generation knob", {
  # §N+26 (2026-06-05): the `otsu_unburned_generation` alias is gone. It is no
  # longer redirected-with-a-warning; it now hits the honest guard, whether it
  # is passed as `negative_pool_policy` or as a bare `otsu_unburned_generation`
  # option. The pipeline is always all_sources.
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(
    build_supervised_burned_config(
      scenario = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(negative_pool_policy = "otsu_unburned_generation")
    ),
    regexp = "was removed in 2026-06"
  )
  expect_error(
    build_supervised_burned_config(
      scenario = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(otsu_unburned_generation = "otsu_unburned_generation")
    ),
    regexp = "was removed in 2026-06"
  )
})

# ===========================================================================
# Gate 1B PIECE 4 (2026-06-08): supervised ecoregion Otsu branch removed.
# The canonical supervised negative-pool Otsu mode is `burnable_only`; the
# only other supported mode is the ecoregion-free `corine`. The supervised
# pipeline neither exposes nor derives any ecoregion path. The DETERMINISTIC
# CORINE x ecoregion stratified Otsu is unaffected.
# ===========================================================================

test_that("PIECE 4: legacy unburned otsu_mode enum is burnable_only + corine only", {
  ns <- asNamespace("OtsuFire")
  fn <- get("build_unburned_from_legacy_pipeline", envir = ns)
  modes <- eval(formals(fn)$otsu_mode)
  expect_setequal(modes, c("burnable_only", "corine"))
  # `burnable_only` is the canonical default (match.arg picks the first).
  expect_identical(modes[[1L]], "burnable_only")
})

test_that("PIECE 4: dispatch + pools default the legacy otsu_mode to burnable_only", {
  ns <- asNamespace("OtsuFire")
  disp <- paste(deparse(get(".of_supervised_engine_bindings", envir = ns)),
                collapse = "\n")
  expect_match(disp, 'legacy_otsu_mode\\s*%\\|\\|%\\s*"burnable_only"')
  pools <- paste(deparse(get("build_supervised_training_pools", envir = ns)),
                 collapse = "\n")
  expect_match(pools, 'legacy_otsu_mode\\s*%\\|\\|%\\s*"burnable_only"')
})

test_that("PIECE 4: no ecoregion path is exposed or derived in the supervised path", {
  ns <- asNamespace("OtsuFire")
  # The cfg builder no longer has an ecoregion_shapefile argument.
  expect_false("ecoregion_shapefile" %in% names(formals(build_supervised_burned_config)))
  # The legacy builder no longer has an ecoregion_shapefile_path argument.
  fn <- get("build_unburned_from_legacy_pipeline", envir = ns)
  expect_false("ecoregion_shapefile_path" %in% names(formals(fn)))
  # No ecoregion path is hardcoded/derived in the legacy builder body.
  src <- paste(deparse(fn), collapse = "\n")
  expect_no_match(src, "ecoregiones_olson", fixed = TRUE)
  expect_no_match(src, "ecoregion_shapefile")
  # The supervised pool builder threads no ecoregion path.
  pools <- paste(deparse(get("build_supervised_training_pools", envir = ns)),
                 collapse = "\n")
  expect_no_match(pools, "ecoregion")
})

test_that("PIECE 4: a cfg without any ecoregion input still builds + keeps corine_raster", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2017L
  )
  # No ecoregion input field exists; corine_raster remains a supervised input.
  expect_false("ecoregion_shapefile" %in% names(cfg$inputs))
  expect_true("corine_raster" %in% names(cfg$inputs))
})

test_that("PIECE 4: supervised Otsu engine carries no ecoregion machinery", {
  ns <- asNamespace("OtsuFire")
  fn <- get("process_otsu_rasters_", envir = ns)
  # No ecoregion / intersection arguments survive in the engine signature.
  fmls <- names(formals(fn))
  for (a in c("ecoregion_shapefile_path", "ecoregion_field",
              "ecoregion_classes", "segment_by_intersection")) {
    expect_false(a %in% fmls)
  }
  body_src <- paste(deparse(body(fn)), collapse = "\n")
  expect_no_match(body_src, "ecoregions_sf")
  expect_no_match(body_src, "units_grouped")
  expect_no_match(body_src, "ECO_ID_INTERNAL")
})

test_that("PIECE 4: DETERMINISTIC CORINE x ecoregion Otsu wiring is preserved", {
  ns <- asNamespace("OtsuFire")
  # The deterministic stratified Otsu engine still accepts its ecoregion args.
  grow <- get("process_otsu_rasters_grow", envir = ns)
  gf <- names(formals(grow))
  for (a in c("ecoregion_shapefile_path", "ecoregion_field",
              "ecoregion_classes", "segment_by_intersection")) {
    expect_true(a %in% gf)
  }
  # The deterministic detection dispatcher still wires the ecoregion path and
  # the CORINE x ecoregion intersection through to the grow engine.
  det <- paste(deparse(get(".of_run_detection", envir = ns)),
               collapse = "\n")
  expect_match(det, "ecoregion_shapefile_path")
  expect_match(det, "segment_by_intersection")
})

