# Tests for Phase 3A public threshold + write (R/supervised-threshold.R).
# Unit tests use synthetic in-memory sf (no local paths). A final integration
# test compares against the real champion 1985 thresholded_burned.gpkg and is
# SKIPPED when the file is absent, so the unit suite never depends on local data.

mk_scored <- function() {
  sq <- function(x0, s) sf::st_polygon(list(rbind(c(x0, 0), c(x0 + s, 0),
                                                  c(x0 + s, s), c(x0, s), c(x0, 0))))
  geom <- sf::st_sfc(sq(0, 100), sq(300, 100), sq(700, 100), sq(1000, 100), sq(1300, 100),
                     crs = 3035)
  sf::st_sf(
    poly_id        = paste0("p", 1:5),
    fire_uid       = paste0("f", 1:5),
    p_burned_model = c(0.95, 0.50, 0.4999, NA, 0.10),  # keep: p1,p2 ; drop: p3,p4(NA),p5
    geometry       = geom)
}

test_that("keeps exactly rows with finite score >= threshold (>= is inclusive)", {
  sc  <- mk_scored()
  out <- tempfile("thr_"); dir.create(out)
  res <- write_thresholded_burned(sc, out_dir = out, threshold = 0.50, overwrite = TRUE, verbose = FALSE)
  expect_equal(res$n_in, 5L)
  expect_equal(res$n_kept, 2L)              # 0.95 and 0.50 (inclusive); 0.4999 out
  expect_equal(res$n_nonfinite, 1L)         # the NA row
  kept <- sf::st_read(res$path, layer = "thresholded_burned", quiet = TRUE)
  expect_setequal(sf::st_drop_geometry(kept)$fire_uid, c("f1", "f2"))
})

test_that("CRS and attribute columns (incl. score) are preserved", {
  sc  <- mk_scored()
  out <- tempfile("thr_"); dir.create(out)
  res <- write_thresholded_burned(sc, out_dir = out, threshold = 0.50, overwrite = TRUE, verbose = FALSE)
  expect_equal(res$crs_epsg, 3035L)
  kept <- sf::st_read(res$path, layer = "thresholded_burned", quiet = TRUE)
  expect_true(all(c("poly_id", "fire_uid", "p_burned_model") %in% names(kept)))
  expect_equal(sf::st_crs(kept)$epsg, 3035L)
})

test_that("threshold and score_col are honoured", {
  sc  <- mk_scored()
  names(sc)[names(sc) == "p_burned_model"] <- "prob"
  out <- tempfile("thr_"); dir.create(out)
  res <- write_thresholded_burned(sc, out_dir = out, threshold = 0.90, score_col = "prob",
                                  overwrite = TRUE, verbose = FALSE)
  expect_equal(res$n_kept, 1L)              # only 0.95
})

test_that("overwrite=FALSE errors when the file exists; TRUE replaces it", {
  sc  <- mk_scored()
  out <- tempfile("thr_"); dir.create(out)
  write_thresholded_burned(sc, out_dir = out, overwrite = TRUE, verbose = FALSE)
  expect_error(write_thresholded_burned(sc, out_dir = out, overwrite = FALSE, verbose = FALSE),
               "overwrite=FALSE")
  expect_silent(write_thresholded_burned(sc, out_dir = out, overwrite = TRUE, verbose = FALSE))
})

test_that("missing score column and bad inputs are rejected", {
  sc  <- mk_scored(); out <- tempfile("thr_"); dir.create(out)
  expect_error(write_thresholded_burned(sc, out_dir = out, score_col = "nope", verbose = FALSE),
               "score column")
  expect_error(write_thresholded_burned(42, out_dir = out, verbose = FALSE), "path or an sf")
  expect_error(write_thresholded_burned(sc, out_dir = out, threshold = NA, verbose = FALSE),
               "threshold")
})

# ---- integration: parity vs the real champion 1985 thresholded_burned ----
test_that("reproduces the champion 1985 thresholded_burned from its scored_all", {
  base <- Sys.getenv(
    "OF_CHAMPION_1985_DIR",
    paste0("C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results/1985/SUPERVISED/",
           "1985_balanced_keepPool_2017_2022_2025/_CHAMPION_DROPSLHN20SUBR05_TEST/champion/balanced"))
  scored_fp <- file.path(base, "scored_all.gpkg")
  champ_fp  <- file.path(base, "thresholded_burned.gpkg")
  skip_if_not(file.exists(scored_fp), "champion scored_all.gpkg not on disk")
  skip_if_not(file.exists(champ_fp),  "champion thresholded_burned.gpkg not on disk")
  skip_if_not_installed("sf")

  out <- tempfile("thr_champ_"); dir.create(out)
  res <- write_thresholded_burned(scored_fp, out_dir = out, threshold = 0.50, overwrite = TRUE, verbose = FALSE)

  champ <- sf::st_read(champ_fp, layer = "thresholded_burned", quiet = TRUE)
  mine  <- sf::st_read(res$path,  layer = "thresholded_burned", quiet = TRUE)

  # (1) same polygon count
  expect_equal(res$n_kept, nrow(champ))
  # (2) same CRS
  expect_equal(sf::st_crs(mine)$epsg, sf::st_crs(champ)$epsg)
  # (3) same row identity by fire_uid (and poly_id)
  expect_setequal(as.character(sf::st_drop_geometry(mine)$fire_uid),
                  as.character(sf::st_drop_geometry(champ)$fire_uid))
  expect_setequal(as.character(sf::st_drop_geometry(mine)$poly_id),
                  as.character(sf::st_drop_geometry(champ)$poly_id))
  # (4) total area parity (metric parity; tolerant of writer rounding)
  a_mine  <- sum(as.numeric(sf::st_area(mine)))
  a_champ <- sum(as.numeric(sf::st_area(champ)))
  expect_equal(a_mine, a_champ, tolerance = 1e-6)
})
