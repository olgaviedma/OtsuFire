# GATE 6.2 (2026-06-11):
#   PART A - cap_otsu (otsu_unburned_to_burned_ratio) is the SOLE Otsu selector;
#            the generation-side `legacy_sample_n` pre-thinning was removed, so
#            the FULL valid Otsu drop pool flows into the negative pool and the
#            per-bucket capping draws from it. Removing the pre-thinning changes
#            the Otsu AVAILABILITY (e.g. 2017: 2000 -> full pool) and therefore
#            the SELECTED Otsu ids, BY DESIGN; the selected COUNT stays
#            ceiling(n_burned * cap_otsu) capped by availability.
#   PART B - Otsu > random spatial dedup (Otsu has priority): random rows that
#            share POSITIVE area with any Otsu patch are removed before capping;
#            Otsu rows are never removed; no positive-area random-otsu duplicate
#            remains afterwards.

`%||%` <- function(x, y) if (is.null(x)) y else x

mk_sq_g6 <- function(xmin, ymin, side = 90) {
  sf::st_polygon(list(matrix(
    c(xmin, ymin,
      xmin + side, ymin,
      xmin + side, ymin + side,
      xmin, ymin + side,
      xmin, ymin),
    ncol = 2, byrow = TRUE
  )))
}

# =====================================================================
# PART A - legacy_sample_n removed; cap_otsu is the sole Otsu selector.
# =====================================================================

test_that("PART A: build_unburned_from_legacy_decisions no longer accepts sample_n / sample_props / random_seed", {
  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))
  fmls <- names(formals(fn))
  expect_false("sample_n" %in% fmls)
  expect_false("sample_props" %in% fmls)
  expect_false("random_seed" %in% fmls)
})

test_that("PART A: build_unburned_from_legacy_pipeline no longer accepts sample_n / sample_props / random_seed", {
  fn <- get("build_unburned_from_legacy_pipeline",
            envir = asNamespace("OtsuFire"))
  fmls <- names(formals(fn))
  expect_false("sample_n" %in% fmls)
  expect_false("sample_props" %in% fmls)
  expect_false("random_seed" %in% fmls)
})

test_that("PART A: the stratified pre-thinning helper was removed (no zombie)", {
  ns <- asNamespace("OtsuFire")
  expect_false(exists("sample_stratified_legacy_unburned",
                      envir = ns, inherits = FALSE))
})

test_that("PART A: full valid Otsu drop pool flows on (legacy_unburned_sampled == pool)", {
  skip_if_not_installed("sf")
  fn <- get("build_unburned_from_legacy_decisions",
            envir = asNamespace("OtsuFire"))

  # 5 non-overlapping drop patches, all passing the S_PATCH filter.
  geom <- sf::st_sfc(
    mk_sq_g6(0, 0), mk_sq_g6(500, 0), mk_sq_g6(1000, 0),
    mk_sq_g6(1500, 0), mk_sq_g6(2000, 0),
    crs = 3035
  )
  legacy <- sf::st_sf(DECISION = rep("drop", 5),
                      S_PATCH_PA = rep(0.05, 5),
                      geometry = geom)
  internal <- sf::st_sf(
    geometry = sf::st_sfc(mk_sq_g6(-1000, -1000, side = 50), crs = 3035)
  )
  lp <- tempfile(fileext = ".shp")
  ip <- tempfile(fileext = ".gpkg")
  sf::st_write(legacy, lp, quiet = TRUE)
  sf::st_write(internal, ip, layer = "internal_decisions", quiet = TRUE)

  res <- fn(
    legacy_patches_path = lp,
    internal_decisions_path = ip,
    out_gpkg = NULL,
    use_drop = TRUE, use_review = FALSE, use_keep = FALSE,
    drop_max_s_patch = 0.15, exclude_buffer_m = 0,
    verbose = FALSE
  )
  # No generation-side subsampling: the whole valid pool flows on, and the
  # `sampled` alias equals the full pool exactly.
  expect_equal(nrow(res$legacy_unburned_pool), 5L)
  expect_equal(nrow(res$legacy_unburned_sampled),
               nrow(res$legacy_unburned_pool))
})

test_that("PART A: cap_otsu is the SOLE Otsu selector; count == ceiling(n_burned*cap) capped by availability (2017 invariant)", {
  cap_fn <- get(".of_cap_negative_buckets", envir = asNamespace("OtsuFire"))

  # 2017-like: n_burned = 532, cap_otsu = 1.0. Old behaviour drew 532 from a
  # PRE-THINNED 2000; the new behaviour draws 532 from the FULL pool (e.g.
  # ~2481). The COUNT is identical (ceiling(532*1.0) = 532) but the selected ids
  # differ by design (different availability). We assert the count invariant and
  # that more availability does NOT change the selected count.
  n_burned <- 532L
  caps <- c(contextual = 1.0, random = 1.0, otsu = 1.0)

  draw_count <- function(n_avail) {
    negs <- list(
      contextual = integer(0),
      random     = integer(0),
      otsu       = seq_len(n_avail)
    )
    out <- cap_fn(
      positive_idx = integer(0),
      negatives_by_bucket = negs,
      n_burned = n_burned,
      caps = caps,
      seed = 42L,
      id = as.character(seq_len(n_avail)),
      context = "g6"
    )
    arow <- out$audit[out$audit$bucket == "otsu", ]
    arow$n_selected
  }

  # Pre-thinned availability (old) vs full pool (new): SAME selected count.
  expect_equal(draw_count(2000L), 532L)
  expect_equal(draw_count(2481L), 532L)
  # Availability below the cap -> take all (count == availability).
  expect_equal(draw_count(400L), 400L)
})

# =====================================================================
# PART B - Otsu > random spatial dedup.
# =====================================================================

test_that("PART B (WITH overlap): a random cell inside an Otsu patch is removed; Otsu kept", {
  skip_if_not_installed("sf")
  fn <- get("dedup_random_vs_otsu_unb_legacy", envir = asNamespace("OtsuFire"))

  # Otsu patch is a big square; random cells: one fully INSIDE it (duplicate),
  # one far away (kept), one sharing only an EDGE (kept - zero shared area).
  otsu <- sf::st_sf(geometry = sf::st_sfc(mk_sq_g6(0, 0, side = 300), crs = 3035))
  rnd <- sf::st_sf(
    rid = 1:3,
    geometry = sf::st_sfc(
      mk_sq_g6(100, 100, side = 90),   # inside  -> removed
      mk_sq_g6(5000, 5000, side = 90), # far     -> kept
      mk_sq_g6(300, 0, side = 90),     # edge-adjacent (shares x=300 edge) -> kept
      crs = 3035
    )
  )

  res <- fn(random_sf = rnd, otsu_sf = otsu)
  expect_equal(res$n_random_before, 3L)
  expect_equal(res$n_random_removed, 1L)        # only the truly-inside cell
  expect_equal(res$n_random_final, 2L)
  expect_true(res$area_removed > 0)
  expect_setequal(res$random_kept$rid, c(2L, 3L))

  # No positive-area duplicate remains; Otsu untouched (the function never edits
  # the otsu set, only returns the kept random set).
  ou <- sf::st_union(sf::st_geometry(otsu))
  gi <- suppressWarnings(sf::st_intersection(sf::st_geometry(res$random_kept), ou))
  areas <- if (length(gi)) as.numeric(sf::st_area(gi)) else numeric(0)
  expect_true(all(areas == 0))                  # any residual is boundary-only
})

test_that("PART B (WITHOUT overlap): nothing removed, random pool unchanged", {
  skip_if_not_installed("sf")
  fn <- get("dedup_random_vs_otsu_unb_legacy", envir = asNamespace("OtsuFire"))

  otsu <- sf::st_sf(geometry = sf::st_sfc(mk_sq_g6(0, 0, side = 90), crs = 3035))
  rnd <- sf::st_sf(
    rid = 1:3,
    geometry = sf::st_sfc(
      mk_sq_g6(10000, 0, side = 90),
      mk_sq_g6(20000, 0, side = 90),
      mk_sq_g6(30000, 0, side = 90),
      crs = 3035
    )
  )
  res <- fn(random_sf = rnd, otsu_sf = otsu)
  expect_equal(res$n_random_removed, 0L)
  expect_equal(res$n_random_final, 3L)
  expect_equal(res$area_removed, 0)
  expect_equal(nrow(res$random_kept), nrow(rnd))
})

test_that("PART B: empty inputs are safe no-ops", {
  skip_if_not_installed("sf")
  fn <- get("dedup_random_vs_otsu_unb_legacy", envir = asNamespace("OtsuFire"))
  rnd <- sf::st_sf(rid = 1L,
                   geometry = sf::st_sfc(mk_sq_g6(0, 0), crs = 3035))
  otsu0 <- rnd[0, , drop = FALSE]
  res <- fn(random_sf = rnd, otsu_sf = otsu0)
  expect_equal(res$n_random_removed, 0L)
  expect_equal(nrow(res$random_kept), 1L)

  res2 <- fn(random_sf = rnd[0, , drop = FALSE], otsu_sf = rnd)
  expect_equal(res2$n_random_before, 0L)
  expect_equal(res2$n_random_removed, 0L)
})

# Persisted 2017 proof (no pool rebuild): run the dedup DIRECTLY on the persisted
# random + Otsu geometries. The Gate-6 brief expected 0 overlap; the persisted
# geometry actually contains exactly 3 random cells fully inside Otsu patches
# (plus boundary-only touches that are correctly KEPT). The dedup removes exactly
# those 3 and leaves no positive-area duplicate. (Discrepancy flagged for Natalia.)
test_that("PART B (persisted 2017): dedup removes the genuine location duplicates, leaves no positive-area overlap, never edits Otsu", {
  skip_if_not_installed("sf")
  pools <- file.path(
    "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/00_LONG_RUNS",
    "LONG_2017_BALANCED_20260609_215851/SHARED/ENGINE_ROUTES/2017/Min_Min",
    "SUPERVISED/balanced/01_POOLS/2017_balanced_pools.gpkg"
  )
  skip_if_not(file.exists(pools), "persisted 2017 pools.gpkg not present")

  old_s2 <- sf::sf_use_s2(); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  rnd <- sf::st_read(pools, "unburned_random", quiet = TRUE)
  fr  <- sf::st_read(pools, "unburned_final_raw", quiet = TRUE)
  otsu <- fr[as.character(fr$source) == "otsu_patch_residual", , drop = FALSE]
  n_otsu_before <- nrow(otsu)

  fn <- get("dedup_random_vs_otsu_unb_legacy", envir = asNamespace("OtsuFire"))
  res <- fn(random_sf = rnd, otsu_sf = otsu)

  expect_equal(res$n_random_before, nrow(rnd))
  expect_equal(res$n_random_removed, 3L)          # genuine inside-otsu duplicates
  expect_equal(res$n_random_final, nrow(rnd) - 3L)
  expect_true(res$area_removed > 0)

  # Otsu set is NEVER modified by the dedup (priority direction).
  expect_equal(nrow(otsu), n_otsu_before)

  # No positive-area random-otsu duplicate remains.
  ou <- sf::st_union(sf::st_geometry(otsu))
  gi <- suppressWarnings(sf::st_intersection(sf::st_geometry(res$random_kept), ou))
  areas <- if (length(gi)) as.numeric(sf::st_area(gi)) else numeric(0)
  expect_true(all(areas == 0))
})
