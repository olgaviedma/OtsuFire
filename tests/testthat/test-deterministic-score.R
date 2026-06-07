# Tests for score_burned_patches(). Input validation only — no engine runs.

mk_tmp_tif3 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_tmp_mask3 <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}
mk_keep_pool_summary3 <- function(years_used = NULL) {
  md <- if (is.null(years_used)) NULL else
    list(years_used = as.integer(years_used))
  structure(
    list(
      min_area_ha = 10,
      min_pix = 13L,
      promote_percentile = 0.10,
      promote_p_above_ref = 0.50,
      qref_prob = 0.05,
      qref_keep = 250,
      q25_keep = 300,
      pixel_area_ha = 0.81,
      keep_medians = c(260, 310, 400),
      metadata = md
    ),
    class = c("otsufire_keep_pool", "list")
  )
}

mk_det_square_sf <- function(df, sizes, origin_x = 0) {
  stopifnot(nrow(df) == length(sizes))
  geom <- lapply(seq_along(sizes), function(i) {
    x0 <- origin_x + ((i - 1L) * 30)
    y0 <- 0
    s <- sizes[[i]]
    sf::st_polygon(list(rbind(
      c(x0, y0),
      c(x0 + s, y0),
      c(x0 + s, y0 + s),
      c(x0, y0 + s),
      c(x0, y0)
    )))
  })
  sf::st_sf(df, geometry = sf::st_sfc(geom, crs = 3035))
}

test_that("score_burned_patches exists and is no longer NotYetImplemented", {
  expect_true(is.function(score_burned_patches))
  err <- tryCatch(score_burned_patches(), error = function(e) conditionMessage(e))
  expect_false(grepl("not implemented yet", err, ignore.case = TRUE))
})

test_that("score_burned_patches rejects non-config", {
  expect_error(score_burned_patches(burned_candidates = "x.gpkg",
                                    config = list()),
               regexp = "build_burned_mapping_config")
})

test_that("score_burned_patches requires burned_candidates", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(score_burned_patches(config = cfg),
               regexp = "burned_candidates")
})

test_that("score_burned_patches fails when burned_candidates path is missing", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(
    score_burned_patches("nowhere/at/all.gpkg", cfg),
    regexp = "burned_candidates.*does not exist"
  )
})

test_that("score_burned_patches rejects bad burned_candidates types", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  expect_error(
    score_burned_patches(burned_candidates = 42, config = cfg),
    regexp = "burned_candidates"
  )
})

test_that("score_burned_patches rejects config with NULL change_index", {
  bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(
    change_index = NULL, burnable_mask = bm, target_year = 2025L
  )
  # Any plausible burned_candidates placeholder will fail at the
  # change_index check first.
  expect_error(
    score_burned_patches(burned_candidates = "x.gpkg", config = cfg),
    regexp = "change_index"
  )
})

test_that("deterministic scoring is decoupled from any external registry", {
  # The scoring adapter must not read config$registry_path nor call the
  # registry pool builder. Resolution is strictly explicit keep_pool or
  # local fallback (engine cannot run here, so assert on the adapter body).
  body_src <- paste(
    deparse(body(OtsuFire:::.of_run_scoring)), collapse = "\n")
  expect_false(grepl("registry_path", body_src, fixed = TRUE))
  expect_false(grepl("build_rbr_keep_pool_from_registry", body_src,
                      fixed = TRUE))
  expect_false(grepl("registry_keep_pool\\$keep_pool", body_src))
})

test_that("score_burned_patches rejects bad keep_pool types", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  # Pass an existing path so we get past the burned_candidates check and hit
  # the keep_pool type check.
  fake_candidates <- tempfile(fileext = ".gpkg"); file.create(fake_candidates)
  expect_error(
    score_burned_patches(fake_candidates, cfg, keep_pool = 1:10),
    regexp = "keep_pool"
  )
})

test_that("keep_pool normalization accepts flat summaries and legacy wrappers", {
  env <- asNamespace("OtsuFire")
  pool <- mk_keep_pool_summary3()

  expect_identical(env$.of_normalize_keep_pool_arg(pool), pool)
  expect_identical(env$.of_normalize_keep_pool_arg(list(keep_pool = pool)),
                   pool)
  expect_error(
    env$.of_normalize_keep_pool_arg(list(foo = 1)),
    regexp = "keep-pool summary list"
  )
})

test_that("same-year explicit keep_pool is rejected when metadata reveal target-year overlap", {
  env <- asNamespace("OtsuFire")
  pool_same  <- mk_keep_pool_summary3(2025L)
  pool_other <- mk_keep_pool_summary3(c(2017L, 2022L))
  pool_nomd  <- mk_keep_pool_summary3()

  expect_error(
    env$.of_validate_keep_pool_target_year(pool_same, 2025L),
    regexp = "cannot include the target year"
  )
  expect_identical(
    env$.of_validate_keep_pool_target_year(pool_other, 2025L),
    pool_other
  )
  expect_identical(
    env$.of_validate_keep_pool_target_year(pool_nomd, 2025L),
    pool_nomd
  )
})

test_that("score_burned_patches rejects malformed keep_pool lists early", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  fake_candidates <- tempfile(fileext = ".gpkg"); file.create(fake_candidates)

  expect_error(
    score_burned_patches(fake_candidates, cfg, keep_pool = list(foo = 1)),
    regexp = "keep-pool summary list"
  )
})

test_that("score_burned_patches rejects same-year explicit keep_pool before engine work", {
  ci <- mk_tmp_tif3(); bm <- mk_tmp_mask3()
  cfg <- build_burned_mapping_config(change_index = ci, burnable_mask = bm,
                                      target_year = 2025L)
  fake_candidates <- tempfile(fileext = ".gpkg"); file.create(fake_candidates)

  expect_error(
    score_burned_patches(fake_candidates, cfg,
                         keep_pool = mk_keep_pool_summary3(2025L)),
    regexp = "cannot include the target year"
  )
})

test_that("build_internal_decisions publishes cleaned geometry after previous-year erase", {
  env <- asNamespace("OtsuFire")

  internal_flagged <- mk_det_square_sf(
    data.frame(
      source_poly_id = c(1L, 2L),
      flag_internal = c("keep", "review"),
      why_flag = c("intersects_stage1_seed", "burnable_corine_ok"),
      stringsAsFactors = FALSE
    ),
    sizes = c(10, 10)
  )

  internal_clean <- mk_det_square_sf(
    data.frame(
      source_poly_id = 1L,
      flag_internal = "keep",
      stringsAsFactors = FALSE
    ),
    sizes = 6
  )

  internal_rbr <- mk_det_square_sf(
    data.frame(
      source_poly_id = 1L,
      area_ha = 0.0036,
      n_pix = 4L,
      median_rbr = 250,
      p_above_keep_q25 = 0.8,
      p_above_keep_ref = 0.8,
      percentile_in_keep = 0.9,
      conf_area = "ok",
      conf_pix = "ok",
      flag_rbr = "keep",
      stringsAsFactors = FALSE
    ),
    sizes = 6
  )

  out <- env$build_internal_decisions(
    internal_flagged = internal_flagged,
    internal_clean = internal_clean,
    internal_rbr = internal_rbr,
    target_year = 2025L,
    scenario_name = "test",
    preyear_available = TRUE
  )

  expect_s3_class(out, "sf")
  expect_equal(nrow(out), 2L)

  area_out <- as.numeric(sf::st_area(out))
  expect_equal(area_out[match(1L, out$source_poly_id)], 36, tolerance = 1e-6)
  expect_equal(area_out[match(2L, out$source_poly_id)], 0, tolerance = 1e-6)

  expect_equal(out$preyear_action[match(1L, out$source_poly_id)], "review")
  expect_equal(out$preyear_reason[match(1L, out$source_poly_id)], "overlap_previous_year_removed")
  expect_equal(out$preyear_action[match(2L, out$source_poly_id)], "drop")
  expect_equal(out$preyear_reason[match(2L, out$source_poly_id)], "previous_year_conflict")

  expect_equal(out$area_ha[match(1L, out$source_poly_id)], 0.0036)
  expect_true(sf::st_is_empty(out$geometry[match(2L, out$source_poly_id)]))
})
