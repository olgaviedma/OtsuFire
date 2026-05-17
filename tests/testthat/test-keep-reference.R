# Focused tests for the historical keep-reference helpers:
#   collect_keep_reference_samples() and build_keep_pool_from_samples().
# Minimal, controlled fixtures only. These functions are standalone
# (not integrated into score_burned_patches() / run_deterministic_pipeline()).

# --- fixtures --------------------------------------------------------------

# Build a tiny sf of square polygons in EPSG:3035 (meters).
mk_polys <- function(df) {
  n <- nrow(df)
  geom <- lapply(seq_len(n), function(i) {
    x0 <- (i - 1L) * 30
    y0 <- 0
    sf::st_polygon(list(rbind(
      c(x0,       y0),
      c(x0 + 20,  y0),
      c(x0 + 20,  y0 + 20),
      c(x0,       y0 + 20),
      c(x0,       y0)
    )))
  })
  sf::st_sf(df, geometry = sf::st_sfc(geom, crs = 3035))
}

# A "good" deterministic decision row that passes strict_hotspot.
good_row <- function(n, year = 2020L) {
  data.frame(
    year           = rep(year, n),
    class_final    = rep("keep", n),
    filter_1       = rep("keep", n),
    filter_2       = rep("keep", n),
    filter_3       = rep("keep", n),
    area_ha        = rep(250, n),
    conf_area      = rep("ok", n),
    conf_pix       = rep("ok", n),
    source_poly_id = seq_len(n),
    stringsAsFactors = FALSE
  )
}

write_decision_gpkg <- function(sf_obj, layer = "internal_decisions") {
  f <- tempfile(fileext = ".gpkg")
  sf::st_write(sf_obj, f, layer = layer, quiet = TRUE)
  f
}

# --- collect_keep_reference_samples() --------------------------------------

test_that("collect_keep_reference_samples errors when a path is missing", {
  expect_error(
    collect_keep_reference_samples(c("definitely/not/here.gpkg")),
    regexp = "do not exist"
  )
})

test_that("collect_keep_reference_samples errors on missing required column", {
  d <- good_row(3)
  d$class_final <- NULL          # drop a mandatory column
  p <- write_decision_gpkg(mk_polys(d))
  expect_error(
    collect_keep_reference_samples(p),
    regexp = "Missing required column.*class_final"
  )
})

test_that("collect_keep_reference_samples excludes target_year", {
  d <- rbind(good_row(3, 2020L), good_row(2, 2021L))
  d$source_poly_id <- seq_len(nrow(d))
  p <- write_decision_gpkg(mk_polys(d))

  res <- collect_keep_reference_samples(p, target_year = 2021L)
  expect_s3_class(res, "otsufire_keep_reference_samples")
  expect_equal(res$summary$target_year_excluded, 2021L)
  expect_true(all(res$samples$year == 2020L))
  expect_equal(res$summary$n_polygons_selected, 3L)

  ec <- res$exclusion_counts
  expect_equal(ec$n[ec$rule == "excluded_target_year"], 2L)
  expect_equal(ec$n[ec$rule == "selected"], 3L)
})

test_that("collect_keep_reference_samples applies strict_hotspot correctly", {
  d <- good_row(5)
  d$source_poly_id <- 1:5
  d$filter_3[2]   <- "review"   # fails filter_3
  d$area_ha[3]    <- 10         # fails min_area_ha (default 100)
  d$conf_pix[4]   <- "too_few_pix" # fails conf_pix
  d$preyear_action <- "keep"
  d$preyear_action[5] <- "drop" # fails preyear rule
  p <- write_decision_gpkg(mk_polys(d))

  res <- collect_keep_reference_samples(p)
  # Only row 1 passes every predicate.
  expect_equal(res$summary$n_polygons_selected, 1L)
  expect_equal(res$samples$source_poly_id, 1L)

  ec <- res$exclusion_counts
  expect_equal(ec$n[ec$rule == "failed_filter_3_keep"], 1L)
  expect_equal(ec$n[ec$rule == "failed_min_area_ha"], 1L)
  expect_equal(ec$n[ec$rule == "failed_conf_pix"], 1L)
  expect_equal(ec$n[ec$rule == "failed_preyear_drop"], 1L)
  expect_equal(ec$n[ec$rule == "input_total"], 5L)
})

test_that("collect_keep_reference_samples errors when nothing is selected", {
  d <- good_row(3)
  d$class_final <- "drop"        # nothing can pass
  p <- write_decision_gpkg(mk_polys(d))
  expect_error(
    collect_keep_reference_samples(p),
    regexp = "No polygons satisfy the strict_hotspot rule"
  )
})

test_that("collect_keep_reference_samples resolves run_label from run_name or scenario", {
  d1 <- good_row(2, 2020L); d1$run_name <- "RUN_A"
  d2 <- good_row(2, 2019L); d2$scenario <- "SCEN_B"
  d2$source_poly_id <- 3:4
  p1 <- write_decision_gpkg(mk_polys(d1))
  p2 <- write_decision_gpkg(mk_polys(d2))

  res <- collect_keep_reference_samples(c(p1, p2))
  expect_setequal(unique(res$samples$run_label), c("RUN_A", "SCEN_B"))
  expect_true(all(c("RUN_A", "SCEN_B") %in% res$summary$run_labels_used))
})

# --- build_keep_pool_from_samples() ----------------------------------------

mk_change_raster <- function(seed = 1L) {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 100,
                   resolution = 10, crs = "EPSG:3035")
  set.seed(seed)
  terra::values(r) <- runif(terra::ncell(r), 0, 1)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mk_samples <- function(n, year = 2020L) {
  d <- data.frame(
    year           = rep(year, n),
    source_poly_id = seq_len(n),
    stringsAsFactors = FALSE
  )
  mk_polys(d)
}

test_that("build_keep_pool_from_samples errors on missing year_col / id_col", {
  s <- mk_samples(5)
  r <- mk_change_raster()
  expect_error(
    build_keep_pool_from_samples(s, list(`2020` = r),
                                 year_col = "missing_year",
                                 min_samples = 3L),
    regexp = "Missing required column.*missing_year"
  )
  expect_error(
    build_keep_pool_from_samples(s, list(`2020` = r),
                                 sample_id_col = "missing_id",
                                 min_samples = 3L),
    regexp = "Missing required column.*missing_id"
  )
})

test_that("build_keep_pool_from_samples errors when a year raster is absent", {
  s <- rbind(mk_samples(3, 2020L), mk_samples(2, 2021L))
  s$source_poly_id <- seq_len(nrow(s))
  r <- mk_change_raster()
  expect_error(
    build_keep_pool_from_samples(s, list(`2020` = r), min_samples = 3L),
    regexp = "missing entries for year.*2021"
  )
})

test_that("build_keep_pool_from_samples errors when nrow(samples) < min_samples", {
  s <- mk_samples(4)
  r <- mk_change_raster()
  expect_error(
    build_keep_pool_from_samples(s, list(`2020` = r), min_samples = 25L),
    regexp = "Not enough samples"
  )
})

test_that("build_keep_pool_from_samples builds a score-compatible keep_pool", {
  s <- mk_samples(6, 2020L)
  r <- mk_change_raster()
  pool <- build_keep_pool_from_samples(s, list(`2020` = r),
                                       min_samples = 3L, quiet = TRUE)

  expect_s3_class(pool, "otsufire_keep_pool")

  # Exact contract consumed by score_burned_patches() -> score_rbr_keep_classes().
  req <- c("min_area_ha", "min_pix", "promote_percentile",
           "promote_p_above_ref", "qref_prob", "qref_keep",
           "q25_keep", "pixel_area_ha", "keep_medians")
  expect_length(setdiff(req, names(pool)), 0L)

  # Fixed V1 internal parameters.
  expect_equal(pool$qref_prob, 0.05)
  expect_equal(pool$promote_p_above_ref, 0.50)
  expect_equal(pool$promote_percentile, 0.10)
  expect_equal(pool$min_area_ha, 10)

  # min_pix must be DERIVED (not NULL) with the exact engine logic
  # max(ceiling(min_area_ha / pixel_area_ha), 13). The fixture raster is
  # 10 m x 10 m => pixel_area_ha = 0.01 => ceiling(10 / 0.01) = 1000.
  expect_false(is.null(pool$min_pix))
  expect_true(is.finite(pool$min_pix))
  expect_equal(pool$min_pix,
               as.integer(max(ceiling(pool$min_area_ha /
                                       pool$pixel_area_ha), 13)))
  expect_equal(pool$min_pix, 1000L)

  # Functional compatibility: the engine, when keep_pool is supplied,
  # does `min_pix <- keep_pool$min_pix` and then `n_pix >= min_pix`.
  # That comparison must yield a usable logical (length 1), which fails
  # if min_pix were NULL.
  expect_length(c(5L, 2000L) >= pool$min_pix, 2L)
  expect_identical(c(5L, 2000L) >= pool$min_pix, c(FALSE, TRUE))

  # Computed scalars / vectors.
  expect_true(is.numeric(pool$qref_keep) && length(pool$qref_keep) == 1L)
  expect_true(is.numeric(pool$q25_keep) && length(pool$q25_keep) == 1L)
  expect_true(is.numeric(pool$keep_medians) && length(pool$keep_medians) >= 1L)
  expect_true(all(is.finite(pool$keep_medians)))
  expect_true(is.numeric(pool$pixel_area_ha) && pool$pixel_area_ha > 0)

  # Metadata.
  expect_equal(pool$metadata$n_samples, 6L)
  expect_equal(pool$metadata$n_years, 1L)
  expect_equal(pool$metadata$years_used, 2020L)
  expect_true(pool$metadata$n_pixels_used >= 1L)
  expect_equal(pool$metadata$year_col, "year")
  expect_equal(pool$metadata$sample_id_col, "source_poly_id")
})
