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
    run_label = "balanced",
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
    run_label = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_b$negative_pool_policy, "all_sources")

  cfg_o <- build_supervised_burned_config(
    run_label = "original", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_identical(cfg_o$negative_pool_policy, "all_sources")

  # Historical canonical value accepted as a no-op (no error, no warning).
  cfg_over <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L,
    options = list(negative_pool_policy = "all_sources")
  )
  expect_identical(cfg_over$negative_pool_policy, "all_sources")

  # Any stale value errors via the honest guard.
  expect_error(
    build_supervised_burned_config(
      run_label = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(negative_pool_policy = "deterministic_direct")
    ),
    regexp = "was removed in 2026-06"
  )
})

test_that("supervised config rejects invalid inputs", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(
    build_supervised_burned_config(run_label = "balanced",
                                    internal_decisions = id,
                                    change_index = ci),
    regexp = "target_year"
  )
  expect_error(
    build_supervised_burned_config(run_label = "balanced",
                                    change_index = ci, target_year = 2025L),
    regexp = "internal_decisions"
  )
  expect_error(
    build_supervised_burned_config(run_label = "balanced",
                                    internal_decisions = id,
                                    target_year = 2025L),
    regexp = "change_index"
  )
  expect_error(
    build_supervised_burned_config(run_label = "balanced",
                                    internal_decisions = id,
                                    change_index = ci, target_year = 2025L,
                                    min_burned_pool_n = -1),
    regexp = "min_burned_pool_n"
  )
})

test_that("run_label is a free-text label (any non-empty string accepted)", {
  # 2026-06-14 breaking change: `run_label` replaced the former `scenario`
  # enum. There is NO match.arg / enum check anymore; any non-empty character
  # string is a valid run label and is stored verbatim in cfg$scenario (the
  # internal field name is unchanged and drives all downstream naming/routing).
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()

  cfg_free <- build_supervised_burned_config(
    run_label = "anything_free", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )
  expect_s3_class(cfg_free, "otsufire_supervised_burned_config")
  expect_identical(cfg_free$scenario, "anything_free")

  # An empty string is rejected with a run_label-specific message.
  expect_error(
    build_supervised_burned_config(run_label = "",
                                    internal_decisions = id,
                                    change_index = ci, target_year = 2025L),
    regexp = "run_label"
  )
  # A non-character value is rejected too.
  expect_error(
    build_supervised_burned_config(run_label = 123,
                                    internal_decisions = id,
                                    change_index = ci, target_year = 2025L),
    regexp = "run_label"
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
  # `db` is a fresh temp dir, so the target `.shp` does not exist yet;
  # `delete_dsn = TRUE` would make GDAL probe the missing file and emit a
  # spurious "GDAL Error 1: ... does not appear to be a file" message.
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc),
               file.path(db, "Borders", "Iberian_peninsula.shp"),
               quiet = TRUE)
  list(data_base = db, composite_base = cb)
}

test_that("§N+25: RUN inputs default to the convention paths", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ty <- 2012L; cy <- "2012"
  paths <- mk_fake_data_base(ty, cy)

  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci,
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
    run_label = "balanced", internal_decisions = id, change_index = ci,
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
      run_label = "balanced", internal_decisions = id, change_index = ci,
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
    run_label = "balanced", internal_decisions = id, change_index = ci,
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
    run_label = "balanced", internal_decisions = id, change_index = ci,
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
      run_label = "balanced", internal_decisions = id,
      change_index = ci, target_year = 2025L,
      options = list(negative_pool_policy = "otsu_unburned_generation")
    ),
    regexp = "was removed in 2026-06"
  )
  expect_error(
    build_supervised_burned_config(
      run_label = "balanced", internal_decisions = id,
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

test_that("PIECE 4: Otsu residual negative otsu_mode enum is burnable_only + corine only", {
  ns <- asNamespace("OtsuFire")
  fn <- get("build_otsu_negative_pipeline", envir = ns)
  modes <- eval(formals(fn)$otsu_mode)
  expect_setequal(modes, c("burnable_only", "corine"))
  # `burnable_only` is the canonical default (match.arg picks the first).
  expect_identical(modes[[1L]], "burnable_only")
})

test_that("PIECE 4 / GATE 6.7: dispatch + pools source the Otsu negative otsu_mode from the fixed internal defaults (burnable_only)", {
  ns <- asNamespace("OtsuFire")
  # GATE 6.7 (2026-06-12): otsu_mode is a FIXED internal default
  # (.of_supervised_negative_pool_internal_defaults()$otsu_mode == "burnable_only"),
  # not a free-form option. Both the dispatch bindings and the pool builder read
  # UNB_OTSU_NEG_MODE from the internal-defaults block (ni$otsu_mode).
  expect_identical(
    get(".of_supervised_negative_pool_internal_defaults", envir = ns)()$otsu_mode,
    "burnable_only")
  disp <- paste(deparse(get(".of_supervised_engine_bindings", envir = ns)),
                collapse = "\n")
  expect_match(disp, "UNB_OTSU_NEG_MODE\\s*=\\s*ni\\$otsu_mode")
  pools <- paste(deparse(get("build_supervised_training_pools", envir = ns)),
                 collapse = "\n")
  expect_match(pools, "UNB_OTSU_NEG_MODE\\s*<-\\s*\\.ni\\$otsu_mode")
})

test_that("PIECE 4: no ecoregion path is exposed or derived in the supervised path", {
  ns <- asNamespace("OtsuFire")
  # The cfg builder no longer has an ecoregion_shapefile argument.
  expect_false("ecoregion_shapefile" %in% names(formals(build_supervised_burned_config)))
  # The Otsu negative builder no longer has an ecoregion_shapefile_path argument.
  fn <- get("build_otsu_negative_pipeline", envir = ns)
  expect_false("ecoregion_shapefile_path" %in% names(formals(fn)))
  # No ecoregion path is hardcoded/derived in the Otsu negative builder body.
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
    run_label = "balanced", internal_decisions = id,
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

# ===========================================================================
# GATE 6.7 (2026-06-12): typed PUBLIC negative_pool_params + runtime_options;
# caps single-source; random_seed in the canonical seeds block; runtime
# excluded from the methodological fingerprint.
# ===========================================================================

mk_g67_cfg <- function(ci, id, ...) {
  build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2025L, ...)
}

test_that("GATE 6.7: typed defaults are visible and stored on the cfg", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  cfg <- mk_g67_cfg(ci, id)
  np <- cfg$negative_pool_params
  expect_equal(np$random$n_cells, 1500L)
  expect_equal(np$random$rbr_quantile, 0.50)
  expect_equal(np$otsu$candidate_threshold, 0)
  expect_equal(np$otsu$reference_threshold, 100)
  expect_equal(unname(np$caps[c("random", "otsu")]), c(1.0, 1.0))
  expect_equal(cfg$negative_pool_runtime,
               list(reuse_existing = TRUE, write_outputs = TRUE, verbose = TRUE))
  # random_seed lives in the canonical seeds block, default 42.
  expect_equal(cfg$train_control$seeds$random_seed, 42L)
})

test_that("GATE 6.7: unknown keys ERROR in negative_pool_params / each sub-list / runtime_options", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(nope = 1)),
               regexp = "Unknown negative_pool_params key")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(random = list(nope = 1))),
               regexp = "Unknown negative_pool_params\\$random key")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(otsu = list(nope = 1))),
               regexp = "Unknown negative_pool_params\\$otsu key")
  expect_error(mk_g67_cfg(ci, id, runtime_options = list(nope = TRUE)),
               regexp = "Unknown runtime_options key")
})

test_that("GATE 6.7: type validation (n_cells / rbr_quantile / thresholds / caps / runtime flags)", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(random = list(n_cells = 0))),
               regexp = "n_cells")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(random = list(rbr_quantile = 1.5))),
               regexp = "rbr_quantile")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(otsu = list(candidate_threshold = "x"))),
               regexp = "candidate_threshold")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(caps = c(random = 1))),
               regexp = "caps")
  expect_error(mk_g67_cfg(ci, id, negative_pool_params = list(caps = c(random = -1, otsu = 1))),
               regexp = "caps")
  expect_error(mk_g67_cfg(ci, id, runtime_options = list(verbose = 1)),
               regexp = "verbose")
})

test_that("GATE 6.7: caps are single-source -- top-level cap_* ERRORS; resolved caps == negative_pool_params$caps", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  expect_error(mk_g67_cfg(ci, id, cap_random = 1.0), regexp = "unused argument")
  expect_error(mk_g67_cfg(ci, id, cap_otsu = 1.0), regexp = "unused argument")
  cfg <- mk_g67_cfg(ci, id,
                    negative_pool_params = list(caps = c(random = 0.75, otsu = 1.5)))
  # The resolved train_control caps equal negative_pool_params$caps (one source).
  expect_equal(cfg$train_control$caps$random, cfg$negative_pool_params$caps[["random"]])
  expect_equal(cfg$train_control$caps$otsu,   cfg$negative_pool_params$caps[["otsu"]])
  expect_equal(cfg$train_control$caps$random, 0.75)
  expect_equal(cfg$train_control$caps$otsu, 1.5)
})

test_that("GATE 6.7: random_seed lives in train_control$seeds with provenance + appears in the methodological fingerprint", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ns <- asNamespace("OtsuFire")
  fp <- get(".of_cfg_run_fingerprint", envir = ns)
  cfg_def  <- mk_g67_cfg(ci, id)
  cfg_seed <- mk_g67_cfg(ci, id, random_seed = 7L)
  expect_equal(cfg_seed$train_control$seeds$random_seed, 7L)
  # provenance recorded.
  expect_equal(cfg_def$resolved_params_provenance$train_control$random_seed, "default")
  expect_equal(cfg_seed$resolved_params_provenance$train_control$random_seed, "user")
  # changing random_seed CHANGES the methodological fingerprint (SAME inputs).
  expect_false(identical(fp(cfg_def)$checksum, fp(cfg_seed)$checksum))
})

test_that("GATE 6.7: a RUNTIME change does NOT alter the methodological fingerprint; a methodological change DOES", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ns <- asNamespace("OtsuFire")
  fp   <- get(".of_cfg_run_fingerprint", envir = ns)
  prev <- get(".of_cfg_neg_pool_fingerprint_preview", envir = ns)
  base <- mk_g67_cfg(ci, id)
  # Runtime-only change: fingerprint unchanged (both the run fp and the neg-pool
  # preview, which exclude runtime).
  rt   <- mk_g67_cfg(ci, id,
                     runtime_options = list(reuse_existing = FALSE,
                                            write_outputs = FALSE, verbose = FALSE))
  expect_identical(fp(base)$checksum, fp(rt)$checksum)
  expect_identical(prev(base, 2025L)$checksum, prev(rt, 2025L)$checksum)
  # Methodological change: each of n_cells / rbr_quantile / candidate_threshold /
  # a cap moves the fingerprint.
  for (np in list(
    list(random = list(n_cells = 999L)),
    list(random = list(rbr_quantile = 0.33)),
    list(otsu   = list(candidate_threshold = 5)),
    list(caps   = c(random = 2.0, otsu = 1.0))
  )) {
    cfg <- mk_g67_cfg(ci, id, negative_pool_params = np)
    expect_false(identical(fp(base)$checksum, fp(cfg)$checksum),
                 info = paste("np change:", names(np)[1]))
  }
})

test_that("GATE 6.7: reproducibility -- same random_seed yields the same resolved cfg negative-pool identity", {
  ci <- mk_tmp_tif_s(); id <- mk_tmp_gpkg_s()
  ns <- asNamespace("OtsuFire")
  prev <- get(".of_cfg_neg_pool_fingerprint_preview", envir = ns)
  c1 <- mk_g67_cfg(ci, id, random_seed = 123L)
  c2 <- mk_g67_cfg(ci, id, random_seed = 123L)
  c3 <- mk_g67_cfg(ci, id, random_seed = 456L)
  expect_identical(prev(c1, 2025L)$checksum, prev(c2, 2025L)$checksum)
  expect_false(identical(prev(c1, 2025L)$checksum, prev(c3, 2025L)$checksum))
})

test_that("GATE 6.7: the FIXED internal defaults are kept (not public, not settable)", {
  ns <- asNamespace("OtsuFire")
  ni <- get(".of_supervised_negative_pool_internal_defaults", envir = ns)()
  # Otsu mode is fixed burnable_only; use_drop fixed TRUE; the rest pinned.
  expect_identical(ni$otsu_mode, "burnable_only")
  expect_true(ni$otsu_use_drop)
  expect_equal(ni$random_exclusion_buffer_m, 500)
  expect_equal(ni$random_patch_size_cells, 3)
  expect_equal(ni$otsu_min_pixels, 8)
  expect_equal(ni$otsu_keep_hi, 0.45)
  expect_equal(ni$otsu_drop_lo, 0.15)
  expect_false(ni$allow_empty_otsu_pool)
  # None of these is a builder argument.
  fmls <- names(formals(build_supervised_burned_config))
  for (k in c("otsu_mode", "otsu_min_pixels", "otsu_keep_hi", "otsu_drop_lo",
              "random_patch_size_cells", "allow_empty_otsu_pool")) {
    expect_false(k %in% fmls, info = k)
  }
})

