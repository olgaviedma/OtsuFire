# Tests for Phase 2B ERA cross-year assembly (R/supervised-era-pool.R).
# Uses synthetic in-memory per-year inputs; reuses apply_visual_validation()
# as the VISUAL entry. Two years share the SAME raw ids to exercise the
# year-prefix collision avoidance.

mk_year <- function(y) {
  pts <- function(n) sf::st_sfc(lapply(seq_len(n), function(i) sf::st_point(c(i, i))), crs = 3035)
  tf <- sf::st_sf(
    fire_uid  = c("f1", "f2", "r1"),
    class     = c("burned", "burned", "unburned"),
    source    = c("internal_keep_qc", "internal_keep_qc", "random_burnable_background"),
    block_id  = c("b1", "b2", "b3"),
    fold_rep1 = c(1L, 2L, 3L),
    fold_rep2 = c(2L, 3L, 1L),
    geometry  = pts(3))
  sc <- sf::st_sf(
    poly_id  = c("c1", "c2"),
    fire_uid = c("cf1", "cf2"),
    block_id = c("sb1", "sb2"),
    geometry = pts(2))
  cl <- data.frame(
    fire_uid = c("f1", "f2", "r1", "cf1", "cf2"),
    poly_id  = c("pf1", "pf2", "pr1", "c1", "c2"),
    pool_source = c("high_confidence_keep", "high_confidence_keep", "random",
                    "artifact_hard", "artifact_hard"),
    artifact_hard_eligible = c(FALSE, FALSE, FALSE, TRUE, TRUE),
    VISUAL   = c(NA, 0L, NA, 1L, 0L),    # f2 dropped; c1 promoted; c2 omission
    stringsAsFactors = FALSE)
  list(year = y, train_features = tf, scoring_features = sc, consolidated_pool = cl)
}

two_years <- function() list(mk_year(1985), mk_year(1988))

test_that("1. multi-year assembly has no fire_uid collisions", {
  era <- assemble_era_training_pool(two_years())
  d <- sf::st_drop_geometry(era$pool)
  expect_false(any(duplicated(d$fire_uid)))
  expect_true(all(grepl("^(1985|1988)_", d$fire_uid)))
})

test_that("2. multi-year assembly has no block_id collisions across years", {
  era <- assemble_era_training_pool(two_years())
  d <- sf::st_drop_geometry(era$pool)
  # each block_id maps to exactly ONE source_year (no cross-year sharing)
  yr_per_block <- tapply(d$source_year, d$block_id, function(x) length(unique(x)))
  expect_true(all(yr_per_block == 1L))
  expect_true(all(grepl("^(1985|1988)_", d$block_id)))
})

test_that("5. VISUAL=0 real fires stay OUT of training and land in omissions", {
  era <- assemble_era_training_pool(two_years())
  d <- sf::st_drop_geometry(era$pool)
  om <- sf::st_drop_geometry(era$omissions)
  # dropped positive f2 never enters the pool (prefixed)
  expect_false(any(d$fire_uid %in% c("1985_f2", "1988_f2")))
  # the VISUAL=0 candidate c2 is omitted, not trained
  expect_true(all(om$omit_reason == "VISUAL0_real_fire"))
  expect_setequal(om$poly_id, c("c2", "c2"))
  expect_equal(nrow(om), 2L)
})

test_that("6. VISUAL=1 artifact_hard candidates are promoted to hard negatives", {
  era <- assemble_era_training_pool(two_years())
  d <- sf::st_drop_geometry(era$pool)
  hard <- d[d$source == "artifact_hard", , drop = FALSE]
  expect_equal(nrow(hard), 2L)                          # c1 from each year
  expect_true(all(hard$class == "unburned"))
  expect_true(all(hard$neg_type == "artifact_hard_negative"))
})

test_that("7. hard-negative folds are valid and reproducible", {
  era1 <- assemble_era_training_pool(two_years())
  era2 <- assemble_era_training_pool(two_years())
  d1 <- sf::st_drop_geometry(era1$pool); d2 <- sf::st_drop_geometry(era2$pool)
  h1 <- d1[d1$source == "artifact_hard", ]
  # folds drawn from the valid CV-fold set {1,3} (positives/background folds)
  expect_true(all(h1$fold_rep1 %in% c(1L, 3L)))
  expect_true(all(h1$fold_rep2 %in% c(1L, 3L)))
  # reproducible across runs (fixed hard_neg_fold_seed)
  h2 <- d2[d2$source == "artifact_hard", ]
  expect_identical(h1$fold_rep1, h2$fold_rep1)
  expect_identical(h1$fold_rep2, h2$fold_rep2)
})

test_that("assembled pool is an sf and the summary tallies per year", {
  era <- assemble_era_training_pool(two_years())
  expect_s3_class(era$pool, "sf")
  expect_equal(nrow(era$summary), 2L)
  expect_equal(sum(era$summary$hard_neg), 2L)
  expect_equal(sum(era$summary$positives), 2L)   # f1 from each year (f2 dropped)
})
