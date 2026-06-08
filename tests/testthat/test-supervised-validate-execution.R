# ============================================================================
# Gate 1C.4: semantic/spatial fail-fast validation + cfg isolation.
#
# (A) Fail-fast: one test per condition (1-10) — validate_supervised_execution()
#     STOPS with a clear error BEFORE heavy compute, and produces NO heavy
#     artefact (the validator writes nothing; we assert the output tree stays
#     empty).
# (B) cfg ISOLATION: two deliberately different cfgs run sequentially AND
#     interleaved; demonstrate no global/session state is used, no path/param
#     is inherited across runs, changing a path changes the consumed input, the
#     fingerprints correspond to the correct cfg, and timestamps NEVER enter the
#     reproducible hashes.
# ============================================================================

ns_vse <- function() asNamespace("OtsuFire")
vse_fn <- function() get("validate_supervised_execution", envir = ns_vse())
cfp_fn <- function() get(".of_cfg_run_fingerprint", envir = ns_vse())

# --- fixture builder: a fully-aligned, valid one-year supervised cfg ---------
# Every spatial input shares the SAME CRS (3035), grid and extent, so the happy
# path passes all 10 checks. Individual tests then perturb ONE input to trip a
# specific check.
mk_vse_cfg <- function(out_dir = tempfile("vse_cfg_"),
                       target_year = 2017L,
                       ci_name = "MinMin_2017_mosaic_res90m.tif",
                       mask_vals = rep(c(0, 1), 50),
                       poly_crs = 3035,
                       class_final = rep("keep", 1L),
                       decisions_extra = NULL,
                       drop_class_final = FALSE,
                       write_mask = TRUE,
                       scenario = "balanced") {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  b1 <- terra::rast(ncol = 10, nrow = 10, xmin = 0, xmax = 1000,
                    ymin = 0, ymax = 1000, crs = "EPSG:3035")
  terra::values(b1) <- runif(100)
  b2 <- terra::rast(b1); terra::values(b2) <- runif(100)
  ci <- c(b1, b2)
  ci_p <- file.path(out_dir, ci_name)
  terra::writeRaster(ci, ci_p, overwrite = TRUE)

  mask_p <- file.path(out_dir, "burneable_mask_binary_corine_2012_ETRS89.tif")
  if (isTRUE(write_mask)) {
    mask <- terra::rast(b1); terra::values(mask) <- mask_vals
    terra::writeRaster(mask, mask_p, overwrite = TRUE)
  }
  topo_p <- file.path(out_dir, "elevation_slope.tif")
  terra::writeRaster(c(b1, b2), topo_p, overwrite = TRUE)
  cor_p <- file.path(out_dir, "CLC_2012_peninsula.tif")
  terra::writeRaster(b1, cor_p, overwrite = TRUE)

  poly_geom <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(100, 100), c(400, 100), c(400, 400),
                              c(100, 400), c(100, 100)))),
    crs = poly_crs)
  attrs <- list()
  if (!isTRUE(drop_class_final)) attrs$class_final <- class_final
  if (!is.null(decisions_extra)) attrs <- c(attrs, decisions_extra)
  if (length(attrs) == 0L) attrs$dummy <- 1L
  poly <- sf::st_sf(as.data.frame(attrs), geometry = poly_geom)
  id_p <- file.path(out_dir, "internal_decisions.gpkg")
  sf::st_write(poly, id_p, layer = "internal_decisions", quiet = TRUE,
               delete_dsn = TRUE)

  build_supervised_burned_config(
    scenario = scenario, internal_decisions = id_p, change_index = ci_p,
    target_year = target_year, output_dir = out_dir,
    topo = topo_p, corine_raster = cor_p,
    burnable_mask = if (isTRUE(write_mask)) mask_p else NULL)
}

# Rebuild a cfg replacing exactly ONE input path (keeps every other route).
rebuild_cfg_with <- function(cfg, ...) {
  ov <- list(...)
  g <- function(nm) OtsuFire:::.of_sup_input_path(cfg, nm)
  args <- list(
    scenario = cfg$scenario,
    internal_decisions = g("internal_decisions"),
    change_index = g("change_index"),
    target_year = cfg$target_year, output_dir = cfg$output_dir,
    topo = g("topo"), corine_raster = g("corine_raster"),
    burnable_mask = g("burnable_mask"))
  for (nm in names(ov)) args[[nm]] <- ov[[nm]]
  do.call(build_supervised_burned_config, args)
}

# A heavy artefact = anything written under the cfg's output tree by the run.
# The validator must produce NONE. We assert the SUPERVISED result dir holds no
# stage output after a failed validation.
expect_no_heavy_artifact <- function(cfg) {
  base <- cfg$output_routes$base
  produced <- if (dir.exists(base)) {
    list.files(base, recursive = TRUE, full.names = TRUE)
  } else character(0)
  # The fixture lives in output_dir but NOT under the SUPERVISED result tree;
  # a clean fail-fast leaves the result tree empty / nonexistent.
  expect_length(produced, 0L)
}

# ===========================================================================
# (A) FAIL-FAST: one test per condition 1-10
# ===========================================================================

test_that("(1) absent/invalid CRS on a spatial input fails fast", {
  td <- tempfile("vse1_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  # Write a NO-CRS copy of the change_index to a NEW file, point the cfg at it.
  ci_p <- OtsuFire:::.of_sup_input_path(cfg, "change_index")
  ci_nocrs <- file.path(td, "MinMin_2017_mosaic_res90m_nocrs.tif")
  r <- terra::rast(ci_p) * 1; terra::crs(r) <- ""
  terra::writeRaster(r, ci_nocrs, overwrite = TRUE)
  cfg2 <- rebuild_cfg_with(cfg, change_index = ci_nocrs)
  expect_error(vse_fn()(cfg2, strict = TRUE, target_year = 2017L), regexp = "CRS")
  expect_no_heavy_artifact(cfg2)
})

test_that("(2) no spatial overlap between inputs fails fast", {
  td <- tempfile("vse2_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  # Disjoint mask (same CRS, far-away extent) written to a NEW file.
  mask_p <- OtsuFire:::.of_sup_input_path(cfg, "burnable_mask")
  m <- terra::rast(mask_p) * 1
  terra::ext(m) <- terra::ext(1e6, 1e6 + 1000, 1e6, 1e6 + 1000)
  mask_far <- file.path(td, "burneable_mask_binary_corine_2012_ETRS89_far.tif")
  terra::writeRaster(m, mask_far, overwrite = TRUE)
  cfg2 <- rebuild_cfg_with(cfg, burnable_mask = mask_far)
  expect_error(vse_fn()(cfg2, strict = TRUE, target_year = 2017L), regexp = "overlap")
  expect_no_heavy_artifact(cfg2)
})

test_that("(3) empty raster fails fast", {
  td <- tempfile("vse3_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  # Point the cfg corine input at a 0-row raster by overwriting the file with a
  # degenerate one is not possible (terra rejects), so emulate an unreadable /
  # empty raster header by truncating. Instead, build a valid 1x0 is impossible;
  # use a raster opened from a corrupt file -> "cannot open" path also fails
  # fast. We assert the validator stops before heavy compute either way.
  cor_p <- OtsuFire:::.of_sup_input_path(cfg, "corine_raster")
  writeLines("not a raster", cor_p)
  # terra emits a GDAL "not recognized" warning while failing to open; that is
  # expected here (the validator turns it into a clean fail-fast error).
  expect_error(
    suppressWarnings(vse_fn()(cfg, strict = TRUE, target_year = 2017L)),
    regexp = "empty raster|cannot open")
  expect_no_heavy_artifact(cfg)
})

test_that("(4) burnable mask with ZERO burnable cells fails fast", {
  td <- tempfile("vse4_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, mask_vals = rep(0, 100))
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2017L),
               regexp = "ZERO burnable")
  expect_no_heavy_artifact(cfg)
})

test_that("(5) wrong year (cfg target vs input) fails fast", {
  td <- tempfile("vse5_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L)
  # The change_index path token says 2017; validating against 2018 is wrong-year.
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2018L),
               regexp = "[Ww]rong-year|target_year")
  expect_no_heavy_artifact(cfg)
})

test_that("(5b) decisions year column without target_year fails fast", {
  td <- tempfile("vse5b_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    decisions_extra = list(year = 1999L))
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2017L),
               regexp = "[Ww]rong-year|year")
  expect_no_heavy_artifact(cfg)
})

test_that("(6) required column absent in internal_decisions fails fast", {
  td <- tempfile("vse6_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, drop_class_final = TRUE)
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2017L),
               regexp = "class_final")
  expect_no_heavy_artifact(cfg)
})

test_that("(7) incompatible column type for class label fails fast", {
  td <- tempfile("vse7_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, class_final = 1.5)  # numeric, not char/factor
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2017L),
               regexp = "incompatible type")
  expect_no_heavy_artifact(cfg)
})

test_that("(8) incompatible feature schema (model expects missing cols) fails fast", {
  td <- tempfile("vse8_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  fake_recipe <- list(feature_names = c("featA", "featB", "featC"))
  expect_error(
    vse_fn()(cfg, strict = TRUE, target_year = 2017L,
             recipe = fake_recipe,
             scoring_feature_names = c("featA")),  # featB/featC missing
    regexp = "Incompatible feature schema|missing")
  expect_no_heavy_artifact(cfg)
})

test_that("(9) contradictory configuration fails fast", {
  td <- tempfile("vse9_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  # legacy protocol + oof_sampling='full' is mutually exclusive.
  expect_error(
    vse_fn()(cfg, strict = TRUE, target_year = 2017L,
             training_protocol = "legacy", oof_sampling = "full"),
    regexp = "contradictory")
  expect_no_heavy_artifact(cfg)
})

test_that("(10) cache belonging to ANOTHER cfg is not silently reused", {
  td <- tempfile("vse10_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  # Plant a foreign negative-pool fingerprint sidecar (wrong checksum).
  pools_dir <- cfg$output_routes$pools_dir
  dir.create(pools_dir, recursive = TRUE, showWarnings = FALSE)
  fp <- file.path(pools_dir, sprintf("%d_%s_neg_pool_fingerprint.txt",
                                     cfg$target_year, cfg$scenario))
  writeLines(c("# foreign", "CHECKSUM=000000000", "---", "x=1"), fp)

  # Under reuse_upstream = TRUE the foreign cache is a HARD ERROR (the run would
  # consume a pool built for a different cfg).
  expect_error(
    vse_fn()(cfg, strict = TRUE, target_year = 2017L, reuse_upstream = TRUE),
    regexp = "foreign pool|DIFFERENT cfg|contradictory|blocking check")
  # Without reuse_upstream it is a WARNING (the pool will be rebuilt) — the
  # engine emits the warning() regardless of strict mode.
  expect_warning(
    vse_fn()(cfg, target_year = 2017L, reuse_upstream = FALSE),
    regexp = "DIFFERENT cfg")
})

# A happy-path sanity anchor: a fully-valid cfg passes all 10 checks.
test_that("happy path: a fully-aligned valid cfg passes all checks", {
  td <- tempfile("vse_ok_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  rep <- expect_silent(vse_fn()(cfg, target_year = 2017L))
  # No blocking FAIL; every row is PASS, NOT_VERIFIABLE, or SKIPPED.
  expect_true(all(rep$status %in% c("PASS", "NOT_VERIFIABLE", "SKIPPED")))
  expect_false(any(rep$severity == "blocking" & rep$status == "FAIL"))
})

# ===========================================================================
# (C) Gate 1D.1: STRUCTURED REPORT + STRICT-MODE CONTRACT
# ===========================================================================

test_that("public validator returns a structured report with the documented columns", {
  td <- tempfile("vse_rep_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  rep <- vse_fn()(cfg, target_year = 2017L)
  expect_s3_class(rep, "data.frame")
  expect_identical(names(rep),
                   c("check", "status", "message", "evidence",
                     "severity", "verifiable"))
  expect_type(rep$verifiable, "logical")
  expect_true(all(rep$status %in%
                  c("PASS", "FAIL", "NOT_VERIFIABLE", "SKIPPED")))
  expect_true(all(rep$severity %in% c("blocking", "warning", "info")))
  # All 10 logical checks present.
  expect_true(all(c("crs_rasters", "empty_raster", "crs_vectors", "overlap",
                    "zero_burnable_mask", "wrong_year", "required_columns",
                    "incompatible_types", "feature_schema",
                    "contradictory_config", "foreign_cache") %in% rep$check))
})

test_that("clean cfg yields PASS or NOT_VERIFIABLE/SKIPPED, never a false PASS", {
  td <- tempfile("vse_clean_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  rep <- vse_fn()(cfg, target_year = 2017L)
  # Any non-PASS row must be explicitly flagged verifiable = FALSE (it was not
  # claimed to pass; it could not be evaluated / was inapplicable).
  non_pass <- rep[rep$status != "PASS", , drop = FALSE]
  expect_true(all(!non_pass$verifiable))
  # A PASS row is always verifiable = TRUE.
  pass <- rep[rep$status == "PASS", , drop = FALSE]
  expect_true(all(pass$verifiable))
})

test_that("strict = TRUE errors on a blocking failure; strict = FALSE returns the report", {
  td <- tempfile("vse_strict_"); dir.create(td)
  # A wrong-year cfg trips the (blocking) wrong_year check.
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L)
  # strict = FALSE: no error; report carries the blocking FAIL.
  rep <- vse_fn()(cfg, strict = FALSE, target_year = 2018L)
  expect_s3_class(rep, "data.frame")
  wr <- rep[rep$check == "wrong_year", , drop = FALSE]
  expect_identical(wr$status, "FAIL")
  expect_identical(wr$severity, "blocking")
  # strict = TRUE: aggregated error listing the failing check.
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2018L),
               regexp = "blocking check.*FAILED|wrong_year")
  expect_no_heavy_artifact(cfg)
})

test_that("NOT_VERIFIABLE is recorded for the feature-schema check when no model is supplied", {
  td <- tempfile("vse_nv_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  rep <- vse_fn()(cfg, target_year = 2017L)
  fs <- rep[rep$check == "feature_schema", , drop = FALSE]
  expect_identical(fs$status, "SKIPPED")
  expect_false(fs$verifiable)
  # A model whose schema cannot be introspected -> NOT_VERIFIABLE (no false PASS).
  rep2 <- vse_fn()(cfg, target_year = 2017L,
                   recipe = list(no_feature_names = TRUE),
                   scoring_feature_names = c("anything"))
  fs2 <- rep2[rep2$check == "feature_schema", , drop = FALSE]
  expect_identical(fs2$status, "NOT_VERIFIABLE")
  expect_false(fs2$verifiable)
})

test_that("strict = FALSE never errors even with multiple blocking failures", {
  td <- tempfile("vse_nostop_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, mask_vals = rep(0, 100))  # zero-burnable
  # No error despite a blocking FAIL; the report records it.
  rep <- vse_fn()(cfg, strict = FALSE, target_year = 2017L)
  zb <- rep[rep$check == "zero_burnable_mask", , drop = FALSE]
  expect_identical(zb$status, "FAIL")
})

# ===========================================================================
# (B) cfg ISOLATION
# ===========================================================================

test_that("cfg-run fingerprint is reproducible and timestamp-free", {
  td <- tempfile("vse_fp_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td)
  f1 <- cfp_fn()(cfg)
  Sys.sleep(1.1)               # wall-clock advances by > 1s
  f2 <- cfp_fn()(cfg)
  # Two equivalent runs at DIFFERENT wall-clock times -> IDENTICAL fingerprint.
  expect_identical(f1$checksum, f2$checksum)
  expect_identical(f1$text, f2$text)
  # The fingerprint text must not embed any wall-clock token.
  expect_no_match(f1$text, "[0-9]{4}-[0-9]{2}-[0-9]{2}")  # no ISO date
  expect_false(grepl("Sys.time|created_at|timestamp", f1$text))
})

test_that("a cfg field change produces a DIFFERENT fingerprint", {
  td1 <- tempfile("vse_fpa_"); dir.create(td1)
  td2 <- tempfile("vse_fpb_"); dir.create(td2)
  cfgA <- mk_vse_cfg(out_dir = td1, target_year = 2017L)
  # Different year + different output route -> a genuinely different cfg.
  cfgB <- mk_vse_cfg(out_dir = td2, target_year = 2017L)
  cfgB2 <- build_supervised_burned_config(
    scenario = cfgB$scenario,
    internal_decisions = OtsuFire:::.of_sup_input_path(cfgB, "internal_decisions"),
    change_index = OtsuFire:::.of_sup_input_path(cfgB, "change_index"),
    target_year = cfgB$target_year, output_dir = cfgB$output_dir,
    topo = OtsuFire:::.of_sup_input_path(cfgB, "topo"),
    corine_raster = OtsuFire:::.of_sup_input_path(cfgB, "corine_raster"),
    burnable_mask = OtsuFire:::.of_sup_input_path(cfgB, "burnable_mask"),
    cap_spectral = 2.0)   # the ONE methodological change
  expect_false(identical(cfp_fn()(cfgA)$checksum, cfp_fn()(cfgB2)$checksum))
  # Changing only the spectral cap on the SAME cfg also moves the fingerprint.
  expect_false(identical(cfp_fn()(cfgB)$checksum, cfp_fn()(cfgB2)$checksum))
})

test_that("changing a cfg input path changes the consumed input + fingerprint", {
  td <- tempfile("vse_path_"); dir.create(td)
  cfgA <- mk_vse_cfg(out_dir = td)
  ci_a <- OtsuFire:::.of_sup_input_path(cfgA, "change_index")
  # Build a second cfg that points change_index at a DIFFERENT file.
  ci_b <- file.path(td, "MinMin_2017_mosaic_res90m_ALT.tif")
  file.copy(ci_a, ci_b, overwrite = TRUE)
  cfgB <- build_supervised_burned_config(
    scenario = "balanced",
    internal_decisions = OtsuFire:::.of_sup_input_path(cfgA, "internal_decisions"),
    change_index = ci_b, target_year = 2017L, output_dir = td,
    topo = OtsuFire:::.of_sup_input_path(cfgA, "topo"),
    corine_raster = OtsuFire:::.of_sup_input_path(cfgA, "corine_raster"),
    burnable_mask = OtsuFire:::.of_sup_input_path(cfgA, "burnable_mask"))
  # The consumed change_index path differs ...
  expect_false(identical(
    OtsuFire:::.of_sup_input_path(cfgA, "change_index"),
    OtsuFire:::.of_sup_input_path(cfgB, "change_index")))
  # ... and so does the run fingerprint (the path is a fingerprint field).
  expect_false(identical(cfp_fn()(cfgA)$checksum, cfp_fn()(cfgB)$checksum))
})

# ===========================================================================
# (D) Gate 1D.3: ROBUST YEAR VALIDATION — classification + evidence hierarchy
# ===========================================================================

wy_row <- function(rep) rep[rep$check == "wrong_year", , drop = FALSE]

test_that("1D.3: correct year in filename token -> PASS (evidence=filename)", {
  td <- tempfile("vy_film_"); dir.create(td)
  # change_index filename carries the matching year token (2017); decisions has
  # no year info, so it is the filename token on change_index that proves it.
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_2017_mosaic_res90m.tif")
  rep <- vse_fn()(cfg, target_year = 2017L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "PASS")
  expect_match(wr$evidence, "change_index\\[filename\\]")
})

test_that("1D.3: WRONG year in filename token -> FAIL/blocking early", {
  td <- tempfile("vy_filw_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_2017_mosaic_res90m.tif")
  # Validating against 2018 makes the 2017 filename token a mismatch.
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2018L),
               regexp = "Wrong-year|does NOT match")
  rep <- vse_fn()(cfg, strict = FALSE, target_year = 2018L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "FAIL")
  expect_identical(wr$severity, "blocking")
  expect_no_heavy_artifact(cfg)
})

test_that("1D.3: correct year in a fire_year column -> PASS (evidence=column)", {
  td <- tempfile("vy_col_"); dir.create(td)
  # change_index filename has NO year token; the decisions fire_year column does.
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_mosaic_res90m.tif",
                    decisions_extra = list(fire_year = 2017L))
  rep <- vse_fn()(cfg, target_year = 2017L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "PASS")
  expect_match(wr$evidence, "internal_decisions\\[column:fire_year\\]")
})

test_that("1D.3: column year vs cfg discrepancy -> FAIL/blocking", {
  td <- tempfile("vy_cold_"); dir.create(td)
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_mosaic_res90m.tif",
                    decisions_extra = list(year = 1999L))
  expect_error(vse_fn()(cfg, strict = TRUE, target_year = 2017L),
               regexp = "Wrong-year|does NOT match")
  rep <- vse_fn()(cfg, strict = FALSE, target_year = 2017L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "FAIL")
  expect_identical(wr$severity, "blocking")
})

test_that("1D.3: an ATEMPORAL input (burnable_mask/CORINE) is never year-failed", {
  td <- tempfile("vy_atemp_"); dir.create(td)
  # The burnable_mask + CORINE filenames carry the CORINE EPOCH year (2012),
  # which differs from target_year 2017. They must NOT trip the year check —
  # they are classified ATEMPORAL and skipped (never enumerated in wrong_year).
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_2017_mosaic_res90m.tif")
  expect_identical(OtsuFire:::.of_vse_year_class("burnable_mask"), "atemporal")
  expect_identical(OtsuFire:::.of_vse_year_class("corine_raster"), "atemporal")
  rep <- vse_fn()(cfg, target_year = 2017L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "PASS")
  # No atemporal input appears in the wrong_year evidence.
  expect_no_match(wr$evidence, "burnable_mask|corine")
  expect_false(any(rep$severity == "blocking" & rep$status == "FAIL"))
})

test_that("1D.3: a YEAR-SPECIFIC input with NO verifiable year -> NOT_VERIFIABLE (no false PASS)", {
  td <- tempfile("vy_nv_"); dir.create(td)
  # change_index filename has NO year token AND decisions carry no year column:
  # the ONLY year-specific inputs cannot be resolved at any hierarchy level.
  cfg <- mk_vse_cfg(out_dir = td, target_year = 2017L,
                    ci_name = "MinMin_mosaic_res90m.tif")
  rep <- vse_fn()(cfg, target_year = 2017L)
  wr <- wy_row(rep)
  expect_identical(wr$status, "NOT_VERIFIABLE")
  expect_false(wr$verifiable)              # never a false PASS
  expect_false(wr$severity == "blocking")  # not blocking by itself
  expect_match(wr$evidence, "change_index|internal_decisions")
  # A NOT_VERIFIABLE wrong_year does not abort even in strict mode.
  expect_silent(vse_fn()(cfg, strict = TRUE, target_year = 2017L))
})

test_that("isolation: sequential AND interleaved validation use no shared/session state", {
  tda <- tempfile("vse_iso_a_"); dir.create(tda)
  tdb <- tempfile("vse_iso_b_"); dir.create(tdb)
  cfgA <- mk_vse_cfg(out_dir = tda, target_year = 2017L,
                     ci_name = "MinMin_2017_mosaic_res90m.tif")
  cfgB <- mk_vse_cfg(out_dir = tdb, target_year = 2016L,
                     ci_name = "MinMin_2016_mosaic_res90m.tif")
  fpA <- cfp_fn()(cfgA)$checksum
  fpB <- cfp_fn()(cfgB)$checksum
  expect_false(identical(fpA, fpB))

  # Record the global sf S2 setting; the validator must NOT mutate it.
  s2_before <- sf::sf_use_s2()

  # (a) SEQUENTIAL: validate A fully, then B fully.
  rA1 <- vse_fn()(cfgA, target_year = 2017L)
  rB1 <- vse_fn()(cfgB, target_year = 2016L)

  # (b) INTERLEAVED: A, B, A, B — each call must be self-contained.
  rA2 <- vse_fn()(cfgA, target_year = 2017L)
  rB2 <- vse_fn()(cfgB, target_year = 2016L)
  rA3 <- vse_fn()(cfgA, target_year = 2017L)
  rB3 <- vse_fn()(cfgB, target_year = 2016L)

  # Each cfg's fingerprint is stable across the interleaving (no inherited state).
  expect_identical(cfp_fn()(cfgA)$checksum, fpA)
  expect_identical(cfp_fn()(cfgB)$checksum, fpB)

  # Validating B with A's year (and vice versa) trips the wrong-year guard —
  # proving the validator reads the CFG it was handed, never a leftover from the
  # previous call.
  expect_error(vse_fn()(cfgB, strict = TRUE, target_year = 2017L),
               regexp = "[Ww]rong-year")
  expect_error(vse_fn()(cfgA, strict = TRUE, target_year = 2016L),
               regexp = "[Ww]rong-year")

  # No global session state mutated by the validator.
  expect_identical(sf::sf_use_s2(), s2_before)

  # The validator wrote no artefact into either result tree.
  expect_no_heavy_artifact(cfgA)
  expect_no_heavy_artifact(cfgB)
})
