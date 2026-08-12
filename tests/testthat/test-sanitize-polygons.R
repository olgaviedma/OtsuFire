# ============================================================================
# Regression suite for the empty-geometry-safe supervised polygon sanitiser
# .of_sanitize_supervised_polygons() and the validate_supervised_execution()
# `sanitize_decisions` dry-run check.
#
# Root cause (real-2017 smoke, snapshot GATE1_FINAL_v0.6.0_20260609_141533,
# HEAD bf612a7): the former sanitize_polygons() called
# sf::st_collection_extract(<whole layer>, "POLYGON") indiscriminately. Under
# sf 1.0-20, if the layer contained ANY empty geometry the call returned 0 rows,
# zeroing the layer; audit_deterministic_pools() then aborted with
# "internal_sf is empty". 2017 internal_decisions has 4 empty geoms among 1278.
#
# These tests assert the fix is CORRECT BY CONSTRUCTION (drop empties FIRST;
# st_collection_extract ONLY on real GEOMETRYCOLLECTION rows), independent of
# whether the installed sf version reproduces the original defect.
# ============================================================================

san_fn <- function() get(".of_sanitize_supervised_polygons",
                         envir = asNamespace("OtsuFire"))

# --- geometry constructors --------------------------------------------------
.poly <- function(x0, y0, s = 100) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                            c(x0, y0 + s), c(x0, y0))))
}
.mpoly <- function(x0, y0, s = 100) {
  sf::st_multipolygon(list(
    list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
               c(x0, y0 + s), c(x0, y0))),
    list(rbind(c(x0 + 200, y0), c(x0 + 200 + s, y0), c(x0 + 200 + s, y0 + s),
               c(x0 + 200, y0 + s), c(x0 + 200, y0)))))
}
.empty_poly <- function() sf::st_polygon()
# A GEOMETRYCOLLECTION carrying a polygonal component (+ a point).
.gc_with_poly <- function(x0 = 500, y0 = 500) {
  sf::st_geometrycollection(list(.poly(x0, y0), sf::st_point(c(x0, y0))))
}
# A GEOMETRYCOLLECTION with NO polygonal component (only a point + a line).
.gc_without_poly <- function(x0 = 800, y0 = 800) {
  sf::st_geometrycollection(list(
    sf::st_point(c(x0, y0)),
    sf::st_linestring(rbind(c(x0, y0), c(x0 + 10, y0 + 10)))))
}
# A self-intersecting bow-tie that is INVALID but becomes valid under make_valid.
.bowtie <- function(x0 = 0, y0 = 0) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 100, y0 + 100),
                            c(x0 + 100, y0), c(x0, y0 + 100), c(x0, y0))))
}

mk_sf <- function(geoms, ids = NULL, crs = 3035, extra = NULL) {
  g <- sf::st_sfc(geoms, crs = crs)
  df <- data.frame(stringsAsFactors = FALSE)
  if (is.null(ids)) ids <- seq_along(geoms)
  df <- data.frame(fire_uid = as.character(ids), stringsAsFactors = FALSE)
  if (!is.null(extra)) for (nm in names(extra)) df[[nm]] <- extra[[nm]]
  sf::st_sf(df, geometry = g)
}

# ===========================================================================
# (1) mix of valid polygons + empty geometries: empties dropped, layer kept.
# ===========================================================================
test_that("(1) valid polygons + empty geometries: empties dropped, not zeroed", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.poly(0, 0), .empty_poly(), .poly(300, 300), .empty_poly()),
             ids = c("a", "b", "c", "d"))
  res <- san_fn()(x)
  expect_equal(res$audit$n_input, 4L)
  expect_equal(res$audit$n_empty_before, 2L)
  expect_equal(res$audit$n_output, 2L)
  expect_equal(nrow(res$sf), 2L)
  expect_setequal(as.character(res$sf$fire_uid), c("a", "c"))
})

# ===========================================================================
# (2) MULTIPOLYGON preserved (not corrupted / not dropped).
# ===========================================================================
test_that("(2) MULTIPOLYGON rows are preserved", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.poly(0, 0), .mpoly(0, 500)), ids = c("p", "mp"))
  res <- san_fn()(x)
  expect_equal(res$audit$n_output, 2L)
  expect_true(all(as.character(sf::st_geometry_type(res$sf)) == "MULTIPOLYGON"))
  expect_setequal(as.character(res$sf$fire_uid), c("p", "mp"))
})

# ===========================================================================
# (3) GEOMETRYCOLLECTION WITH polygonal components -> polygon extracted.
# ===========================================================================
test_that("(3) GEOMETRYCOLLECTION with a polygon: polygon component extracted", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.poly(0, 0), .gc_with_poly()), ids = c("plain", "gc"))
  res <- san_fn()(x)
  expect_equal(res$audit$n_geometrycollection, 1L)
  expect_equal(res$audit$n_without_polygon_component, 0L)
  expect_equal(res$audit$n_output, 2L)
  expect_true(all(as.character(sf::st_geometry_type(res$sf)) == "MULTIPOLYGON"))
  expect_setequal(as.character(res$sf$fire_uid), c("plain", "gc"))
})

# ===========================================================================
# (4) GEOMETRYCOLLECTION WITHOUT polygonal components -> discarded + counted.
# ===========================================================================
test_that("(4) GEOMETRYCOLLECTION without a polygon: discarded + counted", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.poly(0, 0), .gc_without_poly()), ids = c("plain", "gcnp"))
  res <- san_fn()(x)
  expect_equal(res$audit$n_geometrycollection, 1L)
  expect_equal(res$audit$n_without_polygon_component, 1L)
  expect_equal(res$audit$n_output, 1L)
  expect_identical(as.character(res$sf$fire_uid), "plain")
  # The discarded GC id is recorded, never vanishes silently.
  expect_true("gcnp" %in% res$audit$discarded_ids)
})

# ===========================================================================
# (5) invalid geometry that becomes valid via make_valid -> kept.
# ===========================================================================
test_that("(5) invalid geometry repaired by make_valid is kept", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.bowtie(), .poly(500, 500)), ids = c("bow", "ok"))
  # The bow-tie is invalid on input.
  expect_false(all(sf::st_is_valid(x)))
  res <- san_fn()(x)
  expect_gte(res$audit$n_invalid_before_make_valid, 1L)
  expect_equal(res$audit$n_output, 2L)
  expect_true(all(sf::st_is_valid(res$sf)))
  expect_setequal(as.character(res$sf$fire_uid), c("bow", "ok"))
})

# ===========================================================================
# (6) geometry that becomes empty after make_valid -> dropped + counted.
# ===========================================================================
test_that("(6) geometry empty after make_valid is dropped + counted", {
  skip_if_not_installed("sf")
  # A zero-area degenerate polygon (all identical vertices / collapsed ring):
  # st_make_valid collapses it to an EMPTY geometry.
  degenerate <- sf::st_polygon(list(rbind(c(0, 0), c(0, 0), c(0, 0), c(0, 0))))
  x <- mk_sf(list(degenerate, .poly(500, 500)), ids = c("deg", "ok"))
  res <- san_fn()(x)
  # The degenerate ring is removed by ONE of the safety steps: it is either
  # empty-before (step 1), empty-after-make_valid (step 3), or collapsed by
  # make_valid to a non-empty NON-polygonal geometry (POINT/LINESTRING) that the
  # step-4 stray-type guard drops + counts (n_without_polygon_component). All
  # three paths drop + count it; assert the survivor is the real polygon.
  expect_equal(res$audit$n_output, 1L)
  expect_identical(as.character(res$sf$fire_uid), "ok")
  n_removed_safely <- res$audit$n_empty_before +
    res$audit$n_empty_after_make_valid +
    res$audit$n_without_polygon_component
  expect_true(n_removed_safely >= 1L)
  # The dropped degenerate row's id is recorded, never vanishes silently.
  expect_true("deg" %in% res$audit$discarded_ids)
})

# ===========================================================================
# (7) completely empty / all-empty layer -> informative error.
# ===========================================================================
test_that("(7) all-empty layer raises an INFORMATIVE error (not a bare empty)", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.empty_poly(), .empty_poly()), ids = c("e1", "e2"))
  expect_error(san_fn()(x),
               regexp = "no usable polygon geometry remains|genuinely contains no")
  # With error_on_empty = FALSE the same input returns the audit (for dry-runs).
  res <- san_fn()(x, error_on_empty = FALSE)
  expect_equal(res$audit$n_input, 2L)
  expect_equal(res$audit$n_output, 0L)
  expect_equal(res$audit$n_empty_before, 2L)
})

# ===========================================================================
# (8) attribute + ID preservation across sanitisation.
# ===========================================================================
test_that("(8) attributes + identifiers are preserved", {
  skip_if_not_installed("sf")
  x <- mk_sf(list(.poly(0, 0), .empty_poly(), .poly(300, 300)),
             ids = c("k1", "k2", "k3"),
             extra = list(class_final = c("keep", "drop", "review"),
                          area_src = c(10.0, 20.0, 30.0)))
  res <- san_fn()(x)
  expect_true(all(c("fire_uid", "class_final", "area_src") %in% names(res$sf)))
  # The two surviving rows keep their EXACT attribute tuples.
  surv <- sf::st_drop_geometry(res$sf)
  expect_equal(surv$class_final[surv$fire_uid == "k1"], "keep")
  expect_equal(surv$class_final[surv$fire_uid == "k3"], "review")
  expect_equal(surv$area_src[surv$fire_uid == "k3"], 30.0)
})

# ===========================================================================
# (9) deterministic order preservation (kept rows keep original order).
# ===========================================================================
test_that("(9) deterministic order of kept rows is preserved", {
  skip_if_not_installed("sf")
  geoms <- list(.poly(0, 0), .empty_poly(), .poly(300, 0),
                .poly(600, 0), .empty_poly(), .poly(900, 0))
  x <- mk_sf(geoms, ids = c("r1", "r2", "r3", "r4", "r5", "r6"))
  res <- san_fn()(x)
  # r2 and r5 were empty; the survivors keep their input relative order.
  expect_identical(as.character(res$sf$fire_uid),
                   c("r1", "r3", "r4", "r6"))
})

# ===========================================================================
# (10) THE REAL SIMPLIFIED CASE: 4 empties among valid polygons must NOT zero
#      the whole layer. n_output == n_valid.
# ===========================================================================
test_that("(10) 4 empties among valid polygons do NOT zero the layer (n_output==n_valid)", {
  skip_if_not_installed("sf")
  set.seed(2017)
  n_valid <- 1274L
  n_empty <- 4L
  valid_geoms <- lapply(seq_len(n_valid), function(i) {
    .poly((i %% 100) * 200, (i %/% 100) * 200, s = 50)
  })
  # Interleave the 4 empties at scattered positions (as in the real layer).
  geoms <- valid_geoms
  empty_positions <- c(13L, 512L, 900L, 1278L)
  out <- vector("list", n_valid + n_empty)
  vi <- 1L
  for (pos in seq_len(n_valid + n_empty)) {
    if (pos %in% empty_positions) {
      out[[pos]] <- .empty_poly()
    } else {
      out[[pos]] <- valid_geoms[[vi]]; vi <- vi + 1L
    }
  }
  ids <- sprintf("F%04d", seq_len(n_valid + n_empty))
  x <- mk_sf(out, ids = ids)
  expect_equal(nrow(x), 1278L)

  res <- san_fn()(x)
  expect_equal(res$audit$n_input, 1278L)
  expect_equal(res$audit$n_empty_before, 4L)
  expect_equal(res$audit$n_output, 1274L)              # NOT zeroed
  expect_equal(nrow(res$sf), n_valid)
  # The 4 empty-position ids are exactly the dropped ones.
  expect_setequal(res$audit$discarded_ids, ids[empty_positions])
  # Every surviving id is a non-empty-position id, in deterministic order.
  expect_identical(as.character(res$sf$fire_uid), ids[-empty_positions])
})

# ===========================================================================
# validate_supervised_execution(): sanitize_decisions dry-run behaviour.
# ===========================================================================

# Reuse the same cfg fixture pattern as the validate-execution suite, but write
# the internal_decisions layer with controllable empty / all-empty content.
mk_cfg_with_decisions <- function(out_dir, poly_list, ids = NULL,
                                  class_final = NULL, target_year = 2017L) {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  b1 <- terra::rast(ncol = 10, nrow = 10, xmin = 0, xmax = 2000,
                    ymin = 0, ymax = 2000, crs = "EPSG:3035")
  terra::values(b1) <- runif(100)
  ci_p <- file.path(out_dir, "MinMin_2017_mosaic_res90m.tif")
  terra::writeRaster(c(b1, b1), ci_p, overwrite = TRUE)
  mask_p <- file.path(out_dir, "burneable_mask_binary_corine_2012_ETRS89.tif")
  m <- terra::rast(b1); terra::values(m) <- rep(c(0, 1), 50)
  terra::writeRaster(m, mask_p, overwrite = TRUE)
  topo_p <- file.path(out_dir, "elevation_slope.tif")
  terra::writeRaster(c(b1, b1), topo_p, overwrite = TRUE)
  cor_p <- file.path(out_dir, "CLC_2012_peninsula.tif")
  terra::writeRaster(b1, cor_p, overwrite = TRUE)

  if (is.null(ids)) ids <- seq_along(poly_list)
  if (is.null(class_final)) class_final <- rep("keep", length(poly_list))
  poly <- sf::st_sf(
    data.frame(fire_uid = as.character(ids),
               class_final = class_final, stringsAsFactors = FALSE),
    geometry = sf::st_sfc(poly_list, crs = 3035))
  id_p <- file.path(out_dir, "internal_decisions.gpkg")
  sf::st_write(poly, id_p, layer = "internal_decisions", quiet = TRUE,
               delete_dsn = TRUE)

  build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id_p, change_index = ci_p,
    target_year = target_year, output_dir = out_dir,
    topo = topo_p, corine_raster = cor_p, burnable_mask = mask_p)
}

sdec_row <- function(rep) rep[rep$check == "sanitize_decisions", , drop = FALSE]

test_that("validate: decisions with empties -> PASS_WITH_SANITIZATION + right counts", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("vse_san_emp_"); dir.create(td)
  # 4 valid polygons + 2 empties = mimic the real 1278 -> 1274 pattern at scale.
  polys <- list(.poly(0, 0), .empty_poly(), .poly(300, 300),
                .poly(600, 600), .empty_poly(), .poly(900, 900))
  cfg <- mk_cfg_with_decisions(td, polys,
                               ids = c("a", "b", "c", "d", "e", "f"))
  rep <- OtsuFire::validate_supervised_execution(cfg, target_year = 2017L)
  sd <- sdec_row(rep)
  expect_identical(sd$status, "PASS_WITH_SANITIZATION")
  expect_identical(sd$severity, "warning")          # non-blocking
  expect_false(sd$severity == "blocking" && sd$status == "FAIL")
  # Evidence carries the explicit counts (6 -> 4, 2 dropped).
  expect_match(sd$evidence, "n_input=6 -> n_output=4")
  expect_match(sd$evidence, "empty_before=2")
  # A PASS_WITH_SANITIZATION is NOT a blocking failure: strict mode does NOT abort.
  expect_silent(OtsuFire::validate_supervised_execution(
    cfg, strict = TRUE, target_year = 2017L))
})

test_that("validate: all-empty decisions layer is BLOCKING (would zero the layer)", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("vse_san_all_"); dir.create(td)
  polys <- list(.empty_poly(), .empty_poly(), .empty_poly())
  cfg <- mk_cfg_with_decisions(td, polys, ids = c("e1", "e2", "e3"))
  # strict = FALSE: report records a blocking FAIL on sanitize_decisions.
  rep <- OtsuFire::validate_supervised_execution(cfg, strict = FALSE,
                                                 target_year = 2017L)
  sd <- sdec_row(rep)
  expect_identical(sd$status, "FAIL")
  expect_identical(sd$severity, "blocking")
  expect_match(sd$message, "COMPLETELY\\s+EMPTY|no usable polygon")
  # strict = TRUE: aggregated error aborts before heavy compute.
  expect_error(
    OtsuFire::validate_supervised_execution(cfg, strict = TRUE,
                                            target_year = 2017L),
    regexp = "COMPLETELY EMPTY|sanitize_decisions|blocking check")
})

test_that("validate: clean decisions (no empties) -> plain PASS (blocking severity)", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("vse_san_clean_"); dir.create(td)
  polys <- list(.poly(0, 0), .poly(300, 300), .poly(600, 600))
  cfg <- mk_cfg_with_decisions(td, polys, ids = c("a", "b", "c"))
  rep <- OtsuFire::validate_supervised_execution(cfg, target_year = 2017L)
  sd <- sdec_row(rep)
  expect_identical(sd$status, "PASS")
  expect_identical(sd$severity, "blocking")
  expect_true(sd$verifiable)
})
