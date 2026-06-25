# Tests for export_visual_validation_pools() (R/supervised-visual-export.R).
# Unit tests use a synthetic in-memory pool (no local paths). The final
# integration test runs on the real 1985 consolidated pool and is SKIPPED when
# the file is absent.

mk_pool <- function() {
  sq <- function(x0, s) sf::st_polygon(list(rbind(c(x0, 0), c(x0 + s, 0),
                                                  c(x0 + s, s), c(x0, s), c(x0, 0))))
  geom <- sf::st_sfc(lapply(seq(0, by = 300, length.out = 8), sq, s = 100), crs = 3035)
  sf::st_sf(
    fire_uid               = paste0("f", 1:8),
    poly_id                = paste0("p", 1:8),
    pool_source            = c("high_confidence_keep", "high_confidence_keep",
                               "otsu", "otsu",
                               "random", "random",
                               "deterministic_drop", "deterministic_drop"),
    artifact_hard_eligible = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    VISUAL                 = c(NA, NA, NA, NA, NA, NA, 1L, 0L),  # one promoted, one omitted
    geometry               = geom)
}

test_that("returns the 5 expected layers with correct bucket assignment", {
  ex <- export_visual_validation_pools(mk_pool(), out_path = NULL, verbose = FALSE)
  expect_named(ex$layers, c("all_tagged", "keep", "otsu_residual",
                            "random_background", "artifact_hard"))
  expect_equal(nrow(ex$layers$keep), 2L)
  expect_equal(nrow(ex$layers$otsu_residual), 2L)
  expect_equal(nrow(ex$layers$random_background), 2L)
  expect_equal(nrow(ex$layers$artifact_hard), 2L)
  expect_equal(nrow(ex$layers$all_tagged), 8L)
  # artifact_hard kept as ONE family, distinguished by visual_action
  va <- sf::st_drop_geometry(ex$layers$artifact_hard)$visual_action
  expect_setequal(as.character(va), c("promoted_artifact_hard", "omitted_real_fire"))
})

test_that("all_tagged front columns and review_bucket are present", {
  ex <- export_visual_validation_pools(mk_pool(), verbose = FALSE)
  nm <- names(sf::st_drop_geometry(ex$layers$all_tagged))
  expect_true(all(c("VISUAL", "review_bucket", "visual_action", "used_for_training",
                    "pool_source", "artifact_hard_eligible", "fire_uid", "poly_id") %in% nm))
  # requested front fields come first (those that exist)
  expect_equal(nm[1], "VISUAL")
  expect_equal(nm[2], "review_bucket")
})

test_that("summary reports counts by bucket / VISUAL / action / training", {
  ex <- export_visual_validation_pools(mk_pool(), verbose = FALSE)
  expect_true(all(c("by_bucket", "by_visual", "by_action", "by_training",
                    "bucket_x_action") %in% names(ex$summary)))
  bb <- ex$summary$by_bucket
  expect_equal(bb$n[bb$review_bucket == "artifact_hard"], 2L)
})

test_that("writes a GPKG with exactly the expected layers", {
  out <- tempfile(fileext = ".gpkg")
  ex <- export_visual_validation_pools(mk_pool(), out_path = out, overwrite = TRUE, verbose = FALSE)
  on.exit(unlink(out, force = TRUE), add = TRUE)
  expect_true(file.exists(out))
  lyrs <- sf::st_layers(out)$name
  expect_setequal(lyrs, c("all_tagged", "keep", "otsu_residual",
                          "random_background", "artifact_hard"))
})

test_that("the all_tagged layer can be re-fed into apply_visual_validation", {
  out <- tempfile(fileext = ".gpkg")
  export_visual_validation_pools(mk_pool(), out_path = out, overwrite = TRUE, verbose = FALSE)
  on.exit(unlink(out, force = TRUE), add = TRUE)
  edited <- sf::st_read(out, layer = "all_tagged", quiet = TRUE)
  vv_re   <- apply_visual_validation(edited)
  vv_orig <- apply_visual_validation(mk_pool())
  # same training / omissions / dropped counts after the round-trip
  expect_equal(nrow(vv_re$training),  nrow(vv_orig$training))
  expect_equal(nrow(vv_re$omissions), nrow(vv_orig$omissions))
  expect_equal(nrow(vv_re$dropped),   nrow(vv_orig$dropped))
})

test_that("overwrite = FALSE refuses to clobber an existing file", {
  out <- tempfile(fileext = ".gpkg")
  export_visual_validation_pools(mk_pool(), out_path = out, overwrite = TRUE, verbose = FALSE)
  on.exit(unlink(out, force = TRUE), add = TRUE)
  expect_error(
    export_visual_validation_pools(mk_pool(), out_path = out, overwrite = FALSE, verbose = FALSE),
    "overwrite=FALSE")
})

test_that("the input object is not modified", {
  p <- mk_pool()
  export_visual_validation_pools(p, verbose = FALSE)
  expect_false("review_bucket" %in% names(p))   # caller's object untouched
})

# ---- integration: the real 1985 consolidated pool ----
test_that("buckets match the real 1985 consolidated pool", {
  cg <- Sys.getenv(
    "OF_CONSOLIDATED_1985",
    paste0("C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results/1985/SUPERVISED/",
           "1985_balanced_keepPool_2017_2022_2025/01_POOLS/supervised_training_pool.gpkg"))
  skip_if_not(file.exists(cg), "1985 consolidated pool not on disk")
  skip_if_not_installed("sf")

  out <- tempfile(fileext = ".gpkg")
  ex <- export_visual_validation_pools(cg, out_path = out, overwrite = TRUE, verbose = FALSE)
  on.exit(unlink(out, force = TRUE), add = TRUE)

  bb <- ex$summary$by_bucket
  getn <- function(b) { v <- bb$n[bb$review_bucket == b]; if (length(v)) v else 0L }
  expect_equal(getn("keep"),              342L)
  expect_equal(getn("otsu_residual"),     27646L)
  expect_equal(getn("random_background"), 12963L)
  expect_equal(getn("artifact_hard"),     255L)
  expect_setequal(sf::st_layers(out)$name,
                  c("all_tagged", "keep", "otsu_residual", "random_background", "artifact_hard"))

  # re-entry parity on real data
  edited  <- sf::st_read(out, layer = "all_tagged", quiet = TRUE)
  vv_re   <- apply_visual_validation(edited)
  vv_orig <- apply_visual_validation(cg)
  expect_equal(nrow(vv_re$training),  nrow(vv_orig$training))
  expect_equal(nrow(vv_re$omissions), nrow(vv_orig$omissions))
})
