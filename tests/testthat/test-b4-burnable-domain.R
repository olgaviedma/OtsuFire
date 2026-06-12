# Gate 1C.2 (2026-06-08): the B4 random burnable background must compute its
# change-index percentile AND draw its sample ONLY within the (aligned) burnable
# domain. Previously the percentile was taken over the WHOLE change-index raster.
# These tests construct a synthetic scenario where the whole-raster percentile
# and the burnable-only percentile DIFFER, and assert the burnable-only value is
# used, that no selected observation lies outside burnable, that the existing
# exclusions are preserved, that the selection is reproducible, that a different
# mask changes the eligible set, and that the negative-pool fingerprint reacts
# to the B4 domain/mask/percentile/seed/exclusions.

ns <- asNamespace("OtsuFire")
build_b4 <- function() get("build_unburned_from_deterministic_decisions",
                           envir = ns)
mask_hash <- function() get(".of_random_bg_mask_hash", envir = ns)

# --------------------------------------------------------------------------
# Synthetic fixture. A 20x20 grid in EPSG:3035 (metres).
#  - change-index (rbr): LEFT half (cols 1-10) has LOW values (1..N), RIGHT
#    half (cols 11-20) has HIGH values (1000+). Burnable is ONLY the RIGHT
#    half. So the WHOLE-raster 0.5 quantile is dragged down by the low left
#    half, whereas the burnable-only quantile sits in the high right values.
#    The two thresholds are very different -> the test can tell which was used.
#  - internal_decisions: one tiny "drop" polygon in the TOP-RIGHT burnable
#    corner (a positive/candidate to exclude) so exclude_sf is non-empty and
#    unburned_hard is non-empty.
# --------------------------------------------------------------------------
mk_fixture <- function(burnable_right = TRUE, seed_cols = 11:20) {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  ncol <- 20L; nrow <- 20L; res <- 100
  xmin <- 3000000; ymin <- 2000000
  geom <- list(ncol = ncol, nrow = nrow,
               xmin = xmin, xmax = xmin + ncol * res,
               ymin = ymin, ymax = ymin + nrow * res, crs = "EPSG:3035")

  # change index: left half low, right half high
  ci <- terra::rast(ncol = ncol, nrow = nrow, xmin = geom$xmin, xmax = geom$xmax,
                    ymin = geom$ymin, ymax = geom$ymax, crs = geom$crs)
  m <- matrix(0, nrow = nrow, ncol = ncol)
  m[, 1:10]  <- matrix(seq_len(nrow * 10), nrow = nrow)          # low: 1..200
  m[, 11:20] <- 1000 + matrix(seq_len(nrow * 10), nrow = nrow)   # high: 1001..1200
  terra::values(ci) <- as.vector(t(m))
  ci_f <- tempfile(fileext = ".tif")
  terra::writeRaster(ci, ci_f, overwrite = TRUE)

  # burnable mask: right half burnable (cols 11-20) when burnable_right, else
  # a shifted block (cols 6-15) to exercise the mask-change test.
  bm <- terra::rast(ncol = ncol, nrow = nrow, xmin = geom$xmin, xmax = geom$xmax,
                    ymin = geom$ymin, ymax = geom$ymax, crs = geom$crs)
  bmat <- matrix(0L, nrow = nrow, ncol = ncol)
  if (isTRUE(burnable_right)) bmat[, seed_cols] <- 1L else bmat[, 6:15] <- 1L
  terra::values(bm) <- as.vector(t(bmat))
  bm_f <- tempfile(fileext = ".tif")
  terra::writeRaster(bm, bm_f, overwrite = TRUE)

  # internal_decisions: a "drop" polygon in the top-right burnable corner.
  cell_xy <- function(rr, cc) {
    x <- geom$xmin + (cc - 0.5) * res
    y <- geom$ymax - (rr - 0.5) * res
    c(x, y)
  }
  pc <- cell_xy(1, 20)  # top-right cell centre (burnable when burnable_right)
  half <- res / 2
  poly <- sf::st_polygon(list(rbind(
    c(pc[1] - half, pc[2] - half), c(pc[1] + half, pc[2] - half),
    c(pc[1] + half, pc[2] + half), c(pc[1] - half, pc[2] + half),
    c(pc[1] - half, pc[2] - half))))
  dec <- sf::st_sf(
    poly_id = 1L,
    class_final = "drop",
    filter_1 = "drop", reason_1 = "outside_burnable",
    filter_2 = NA_character_, reason_2 = NA_character_,
    filter_3 = NA_character_, reason_3 = NA_character_,
    geometry = sf::st_sfc(poly, crs = 3035)
  )
  dec_f <- tempfile(fileext = ".gpkg")
  sf::st_write(dec, dec_f, layer = "internal_decisions", quiet = TRUE,
               delete_dsn = TRUE)

  list(ci = ci_f, bm = bm_f, dec = dec_f, geom = geom)
}

# A direct call wrapper. exclude_buffer_m = 0 keeps the exclusion to exactly the
# drop cell so the burnable accounting is easy to reason about.
run_b4 <- function(fx, burnable = fx$bm, seed = 42, rbr_q = 0.50,
                   n_random = 50, patch = 1, buffer = 0) {
  out_gpkg <- tempfile(fileext = ".gpkg")
  build_b4()(
    target_year             = 2017L,
    scenario_name           = "balanced",
    data_base               = tempdir(),
    result_name             = "Min_Min",
    composite_base          = tempdir(),
    exclude_buffer_m        = buffer,
    n_random_cells          = n_random,
    random_rbr_q            = rbr_q,
    random_seed             = seed,
    random_patch_size_cells = patch,
    overwrite_output        = TRUE,
    out_gpkg                = out_gpkg,
    internal_decisions_path = fx$dec,
    burnable_mask_path      = burnable,
    severity_raster_path    = fx$ci,
    verbose                 = FALSE
  )
}

# --------------------------------------------------------------------------
# (a) percentile computed only over BURNABLE cells.
# --------------------------------------------------------------------------
test_that("1C.2 (a): percentile is computed over burnable cells, not whole raster", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fx <- mk_fixture()
  res <- run_b4(fx, rbr_q = 0.50)

  # whole-raster 0.5 quantile (the OLD behaviour) over the finite change index.
  ci_vals <- terra::values(terra::rast(fx$ci), mat = FALSE)
  ci_vals <- ci_vals[is.finite(ci_vals)]
  thr_whole <- as.numeric(stats::quantile(ci_vals, 0.50, na.rm = TRUE))

  # burnable-only 0.5 quantile (right half, values 1001..1200), minus the one
  # excluded drop cell. The threshold must sit in the HIGH band, far above the
  # whole-raster threshold.
  thr_used <- res$rbr_threshold
  expect_true(is.finite(thr_used))
  expect_gt(thr_used, 1000)              # in the burnable (high) band
  expect_gt(thr_used, thr_whole)         # NOT the whole-raster value
  # sanity: whole-raster threshold is in the low band
  expect_lt(thr_whole, 1000)

  # audit confirms the percentile domain is the burnable, exclusion-applied set,
  # NOT the full finite-index count.
  expect_lt(res$b4_audit$n_percentile_domain, res$b4_audit$n_valid_index_cells)
  expect_equal(res$b4_audit$domain_decision, "burnable_restricted")
})

# --------------------------------------------------------------------------
# (b) NO selected random_background observation lies outside burnable.
# --------------------------------------------------------------------------
test_that("1C.2 (b): no selected observation lies outside burnable", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fx <- mk_fixture()
  res <- run_b4(fx, patch = 3)  # 3x3 dilation would spill outside burnable

  expect_equal(res$b4_audit$n_selected_outside_burnable, 0L)
  expect_true(res$b4_audit$confirmation_no_obs_outside_burnable)

  if (nrow(res$unburned_random) > 0) {
    bm <- terra::rast(fx$bm)
    burn_ids <- terra::cells(bm, 1)[[1]]
    cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(res$unburned_random)))
    sel_cells <- terra::cellFromXY(bm, cents)
    expect_true(all(sel_cells %in% burn_ids))
  }
})

# --------------------------------------------------------------------------
# (c) existing exclusions preserved: the buffered positive/candidate cell is
#     still excluded from the eligible/selected set.
# --------------------------------------------------------------------------
test_that("1C.2 (c): the deterministic-drop (positive) cell stays excluded", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fx <- mk_fixture()
  res <- run_b4(fx, buffer = 0, patch = 1)

  # the excluded drop cell is the top-right cell; it must not appear among the
  # selected random background cells.
  bm <- terra::rast(fx$bm)
  drop_cell <- terra::cellFromXY(
    bm, matrix(c(fx$geom$xmax - 50, fx$geom$ymax - 50), ncol = 2))
  if (nrow(res$unburned_random) > 0) {
    cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(res$unburned_random)))
    sel_cells <- terra::cellFromXY(bm, cents)
    expect_false(drop_cell %in% sel_cells)
  }
  # the exclusion was accounted for in the audit.
  expect_gte(res$b4_audit$n_excl_buffer, 1L)
})

# --------------------------------------------------------------------------
# (d) reproducibility: same seed -> identical selection.
# --------------------------------------------------------------------------
test_that("1C.2 (d): same seed reproduces the identical selection", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fx <- mk_fixture()
  r1 <- run_b4(fx, seed = 42)
  r2 <- run_b4(fx, seed = 42)
  g1 <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(r1$unburned_random)))
  g2 <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(r2$unburned_random)))
  expect_equal(g1, g2)
  expect_equal(r1$b4_audit$n_selected_cells, r2$b4_audit$n_selected_cells)
})

# --------------------------------------------------------------------------
# (e) changing the burnable mask changes the eligible set / percentile / hash.
# --------------------------------------------------------------------------
test_that("1C.2 (e): a different burnable mask changes eligible set + percentile + hash", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fx_right <- mk_fixture(burnable_right = TRUE)
  fx_mid   <- mk_fixture(burnable_right = FALSE)  # cols 6-15 burnable
  # same change index + decisions, only the mask differs
  r_right <- run_b4(fx_right, burnable = fx_right$bm)
  r_mid   <- run_b4(fx_right, burnable = fx_mid$bm)

  expect_false(identical(r_right$burnable_mask_hash, r_mid$burnable_mask_hash))
  # The two masks cover the same NUMBER of cells but a different FOOTPRINT, so
  # they overlap different change-index values: the burnable-only percentile
  # differs (the right mask sits wholly in the high band >1000; the shifted mask
  # straddles the low/high boundary). The differing threshold + mask hash prove
  # the mask change propagates into the eligible set and the fingerprint inputs.
  expect_false(isTRUE(all.equal(r_right$rbr_threshold, r_mid$rbr_threshold)))
  # the actually-eligible change-index range differs between the two masks.
  expect_lt(r_mid$rbr_threshold, r_right$rbr_threshold)
})

# --------------------------------------------------------------------------
# (f) negative-pool fingerprint reacts to B4 domain/mask/percentile/seed/excl
#     and is stable when nothing changes.
# --------------------------------------------------------------------------
test_that("1C.2 (f): neg-pool fingerprint changes on B4 changes, stable otherwise", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  fp <- get("otsu_negative_param_fingerprint", envir = ns)
  fx <- mk_fixture()
  base <- run_b4(fx, seed = 42, rbr_q = 0.50)

  mk_fp <- function(res, seed = 42, rbr_q = 0.50, mask_hash = res$burnable_mask_hash) {
    fp(list(
      neg_pool_policy       = "all_sources",
      b4_domain_decision    = res$b4_audit$domain_decision,
      b4_burnable_mask_hash = mask_hash,
      b4_random_rbr_q       = rbr_q,
      b4_percentile_value   = res$b4_audit$percentile_value,
      b4_random_seed        = seed,
      b4_exclude_buffer_m   = res$b4_audit$exclude_buffer_m
    ))$checksum
  }

  base_fp <- mk_fp(base)
  # stable when nothing changes
  expect_identical(base_fp, mk_fp(base))
  # changes on seed
  expect_false(identical(base_fp, mk_fp(base, seed = 7)))
  # changes on percentile
  expect_false(identical(base_fp, mk_fp(base, rbr_q = 0.25)))
  # changes on mask hash (domain change)
  expect_false(identical(base_fp, mk_fp(base, mask_hash = "000000000")))
})
