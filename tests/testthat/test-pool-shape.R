# Tests for Phase 2B shape features + SUBR05 random subsample
# (R/supervised-pool-shape.R).

mk_subsample_pool <- function() {
  data.frame(
    poly_id = 1:21,
    source  = c(rep("internal_keep_qc", 4),
                rep("random_burnable_background", 10),
                rep("otsu_patch_residual", 5),
                rep("artifact_hard", 2)),
    class   = c(rep("burned", 4), rep("unburned", 17)),
    stringsAsFactors = FALSE)
}

test_that("3. SUBR05 reduces ONLY the random background (ratio x positives), not others", {
  p <- mk_subsample_pool()                 # n_burned = 4 -> target = ceil(4*0.5) = 2
  r <- subsample_random_background(p, random_to_burned_ratio = 0.5, seed = 1985L)
  expect_equal(sum(r$source == "random_burnable_background"), 2L)  # 10 -> 2
  expect_equal(sum(r$class  == "burned"), 4L)                      # positives intact
  expect_equal(sum(r$source == "otsu_patch_residual"), 5L)         # otsu intact
  expect_equal(sum(r$source == "artifact_hard"), 2L)               # hard-negs intact
})

test_that("4. SUBR05 is reproducible under a fixed seed (and seed-sensitive)", {
  p <- mk_subsample_pool()
  r1 <- subsample_random_background(p, random_to_burned_ratio = 0.5, seed = 1985L)
  r2 <- subsample_random_background(p, random_to_burned_ratio = 0.5, seed = 1985L)
  expect_identical(r1$poly_id, r2$poly_id)                         # same rows kept
  kept1 <- sort(r1$poly_id[r1$source == "random_burnable_background"])
  kept9 <- sort(subsample_random_background(p, 0.5, seed = 9L)$poly_id[
                  subsample_random_background(p, 0.5, seed = 9L)$source == "random_burnable_background"])
  expect_false(identical(kept1, kept9))                            # different seed -> different subset
})

test_that("SUBR05 default seed is the canonical champion seed 42L", {
  p <- mk_subsample_pool()
  # default (no seed arg) must equal an explicit seed = 42L, and the default
  # must NOT silently behave like the old 1985L default.
  d  <- subsample_random_background(p, random_to_burned_ratio = 0.5)         # default seed
  e42 <- subsample_random_background(p, random_to_burned_ratio = 0.5, seed = 42L)
  expect_identical(d$poly_id, e42$poly_id)
  expect_identical(as.integer(formals(subsample_random_background)$seed), 42L)
})

test_that("SUBR05 is a no-op when random is already below target", {
  p <- mk_subsample_pool()
  expect_equal(nrow(subsample_random_background(p, random_to_burned_ratio = 5, seed = 1L)), nrow(p))
})

test_that("8. the 6 champion shape columns exist after add_pool_shape_features", {
  sq <- function(x0, s) sf::st_polygon(list(rbind(c(x0, 0), c(x0 + s, 0),
                                                  c(x0 + s, s), c(x0, s), c(x0, 0))))
  geom <- sf::st_sfc(sq(0, 100), sq(300, 150), sq(700, 80), sq(1000, 300), crs = 3035)
  p <- sf::st_sf(fire_uid = c("a", "b", "c", "d"), n_pix = c(10, 20, 5, 40),
                 source = c("internal_keep_qc", "internal_keep_qc",
                            "random_burnable_background", "otsu_patch_residual"),
                 class  = c("burned", "burned", "unburned", "unburned"),
                 geometry = geom)
  res <- add_pool_shape_features(p, seed = 1985L)
  expect_true(all(c("area_ha", "n_pix", "log_area", "perim_m", "compactness", "elongation")
                  %in% names(res)))
  expect_s3_class(res, "sf")
  # decorrelation: the random row 'c' inherits a burned donor's (a/b) area_ha
  rd <- sf::st_drop_geometry(res)
  burned_area <- rd$area_ha[rd$class == "burned"]
  expect_true(rd$area_ha[rd$fire_uid == "c"] %in% burned_area)
})
