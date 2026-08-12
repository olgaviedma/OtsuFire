# Tests for apply_chain_cleaning() / apply_chain()  -- optional de-chain stage.
# Mostly synthetic (no real data needed) so they run fast and anywhere.

skip_if_not_installed("sf")
skip_if_not_installed("terra")

# --- helpers: build a tiny synthetic REFINE_MERGED in a temp dir ----------------
make_refine <- function(geoms, dir = tempfile("run_")) {
  refdir <- file.path(dir, "02_REFINE")
  dir.create(refdir, recursive = TRUE, showWarnings = FALSE)
  sfo <- sf::st_sf(
    CLUMP_ID = seq_along(geoms), YEAR = 1985L, RUN = "test",
    LB = 285, DELTA = 115, MINSEED = 12, INDEX = "RBR",
    WORKFLOW = "wb", ENGINE = "test", SOURCE_AOI_FILE = "synthetic",
    geometry = sf::st_sfc(geoms, crs = 3035))
  p <- file.path(refdir, "BA_1985_REFINE_MERGED_test.gpkg")
  sf::st_write(sfo, p, quiet = TRUE)
  p
}
# two ~600 m squares joined by a thin 90 m neck -> highly inflated (a "chain")
dumbbell <- function() {
  sq <- function(x0, y0, s) sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
  a <- sq(0, 0, 600); b <- sq(1500, 0, 600)
  neck <- sf::st_polygon(list(rbind(c(600, 255), c(1500, 255), c(1500, 345), c(600, 345), c(600, 255))))
  sf::st_union(sf::st_union(a, b), neck)
}
compact_square <- function() sf::st_polygon(list(rbind(
  c(5000, 5000), c(5900, 5000), c(5900, 5900), c(5000, 5900), c(5000, 5000))))

test_that("function is exported under both names", {
  expect_true(exists("apply_chain_cleaning", mode = "function"))
  expect_true(exists("apply_chain", mode = "function"))
  expect_identical(apply_chain, apply_chain_cleaning)
})

test_that("dechain_distance_m / pixel_size_m -> opening_px (90/180/270 m)", {
  p <- make_refine(list(compact_square()))
  expect_equal(apply_chain(p, dechain_distance_m = 90,  pixel_size_m = 90)$opening_px, 1L)
  expect_equal(apply_chain(p, dechain_distance_m = 180, pixel_size_m = 90, overwrite = TRUE)$opening_px, 2L)
  expect_equal(apply_chain(p, dechain_distance_m = 270, pixel_size_m = 90, overwrite = TRUE)$opening_px, 3L)
  expect_error(apply_chain(p, dechain_distance_m = 45, pixel_size_m = 90), "opening_px")  # < 1 px
})

test_that("input REFINE_MERGED is never overwritten; output goes to CHAIN/dechain_<N>m/", {
  p <- make_refine(list(dumbbell()))
  before <- file.info(p)$mtime
  res <- apply_chain(p, dechain_distance_m = 180, verbose = FALSE)
  expect_true(file.exists(p))
  expect_identical(file.info(p)$mtime, before)            # untouched
  expect_true(grepl("CHAIN[\\\\/]dechain_180m", res$output_path))
  expect_true(file.exists(res$output_path))
  expect_false(identical(normalizePath(res$output_path), normalizePath(p)))
})

test_that("output schema is REFINE-compatible (CLUMP_ID + run metadata, valid MULTIPOLYGON)", {
  p <- make_refine(list(dumbbell()))
  res <- apply_chain(p, dechain_distance_m = 180, verbose = FALSE)
  out <- sf::st_read(res$output_path, quiet = TRUE)
  expect_true(all(c("CLUMP_ID","YEAR","RUN","SOURCE_AOI_FILE") %in% names(out)))
  expect_equal(out$CLUMP_ID, seq_len(nrow(out)))          # renumbered, unique
  expect_true(all(sf::st_is_valid(out)))
  expect_equal(sf::st_crs(out)$epsg, 3035L)
  expect_true(all(sf::st_geometry_type(out) == "MULTIPOLYGON"))
})

test_that("audit + params files are written", {
  p <- make_refine(list(dumbbell()))
  res <- apply_chain(p, dechain_distance_m = 180, verbose = FALSE)
  expect_true(file.exists(res$audit_path))
  expect_true(file.exists(res$params_path))
  a <- utils::read.csv(res$audit_path)
  expect_true(all(c("n_candidates_before","n_candidates_after","opening_px","area_removed_ha") %in% names(a)))
})

test_that("a thin neck between two blobs is split (do-the-work)", {
  # NB: a simple dumbbell has inflation ~2.0; real chains reach >22 (default thr).
  # Here we lower inflation_thr to deterministically flag the synthetic chain and
  # exercise the rasterise -> opening -> connected-components -> split mechanism.
  p <- make_refine(list(dumbbell()))
  res <- apply_chain(p, dechain_distance_m = 180, inflation_thr = 1.5, verbose = FALSE)  # cut necks up to 2 px (180 m)
  expect_equal(res$n_candidates_before, 1L)
  expect_gte(res$n_candidates_after, 2L)                  # one chain -> >=2 pieces
  expect_equal(res$n_chained_flagged, 1L)
  expect_gt(res$area_removed_ha, 0)
})

test_that("a compact blob (inflation < thr) passes through untouched (do-no-harm)", {
  p <- make_refine(list(compact_square()))
  res <- apply_chain(p, dechain_distance_m = 180, verbose = FALSE)
  expect_equal(res$n_chained_flagged, 0L)
  expect_equal(res$n_candidates_after, 1L)
  expect_equal(res$area_removed_ha, 0)
  out <- sf::st_read(res$output_path, quiet = TRUE)
  expect_equal(out$SOURCE_AOI_FILE[1], "untouched")
})

test_that("mixed set: compact survives, chain splits", {
  p <- make_refine(list(dumbbell(), compact_square()))
  res <- apply_chain(p, dechain_distance_m = 180, inflation_thr = 1.5, verbose = FALSE)
  expect_equal(res$n_chained_flagged, 1L)                 # only the dumbbell (infl~2.0); compact (~1.13) survives
  expect_gte(res$n_candidates_after, 3L)                  # >=2 split + 1 untouched compact
})

test_that("overwrite = FALSE protects an existing output", {
  p <- make_refine(list(dumbbell()))
  res <- apply_chain(p, dechain_distance_m = 180, verbose = FALSE)
  expect_error(apply_chain(p, dechain_distance_m = 180, verbose = FALSE), "overwrite")
  expect_silent(r2 <- apply_chain(p, dechain_distance_m = 180, overwrite = TRUE, verbose = FALSE))
})
