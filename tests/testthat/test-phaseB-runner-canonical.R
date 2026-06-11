# =============================================================================
# Gate 1B (2026-06-07): the canonical experiment runner configures
# methodological params via the cfg BUILDER (no deprecated function-level
# shims), so the canonical run is WARNING-FREE.
#
# GATE 6.1 (2026-06-11): the dead spectral negative bucket was removed
# end-to-end. The former Phase-B spectral=2.0 axis no longer exists. This file
# now asserts the builder-path canonical contract on an OPERATIVE cap
# (cap_contextual): a builder-set cap is provenance "user", reaches BOTH the OOF
# and FINAL engines, and the canonical call sites stay warning-free.
# =============================================================================

ns <- asNamespace("OtsuFire")

mk_pb_tif <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}
mk_pb_gpkg <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}
mk_pb_oof_feats <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_sf(fire_uid = "a", class = "burned",
            fold_rep1 = 1L, fold_rep2 = 1L, geometry = sfc)
}
mk_pb_final_gpkg <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc),
               g, layer = "train_features", quiet = TRUE, delete_dsn = TRUE)
  g
}

# A canonical (builder-path) cfg with a USER contextual cap.
build_canonical_cfg <- function() {
  build_supervised_burned_config(
    scenario           = "balanced",
    internal_decisions = mk_pb_gpkg(),
    change_index       = mk_pb_tif(),
    target_year        = 2017L,
    cap_contextual     = 1.0
  )
}

# Capture the contextual cap each engine receives, engines mocked so no real fit.
capture_oof_contextual <- function(cfg, train_feats) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      a <- list(...)
      captured$contextual <- a$contextual_exclusion_to_burned_ratio
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL), files = list())
    },
    .package = "OtsuFire"
  )
  run_oof_diagnostics(train_features = train_feats,
                      scoring_features = train_feats, config = cfg,
                      out_dir = tempfile(), matrix_dir = tempfile())
  captured$contextual
}
capture_final_contextual <- function(cfg, train_gpkg) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      a <- list(...)
      captured$contextual <- a$contextual_exclusion_to_burned_ratio
      list(model = NULL, recipe = NULL, training_ok_sf = NULL, split = NULL,
           files = list())
    },
    .package = "OtsuFire"
  )
  train_final_burned_model(train_features = train_gpkg, config = cfg,
                           out_dir = tempfile())
  captured$contextual
}

test_that("canonical cfg: builder-set contextual cap 1.0 with provenance 'user'", {
  cfg <- build_canonical_cfg()
  expect_equal(cfg$train_control$caps$contextual, 1.0)
  expect_equal(cfg$train_control$caps$random, 1.0)
  expect_equal(cfg$train_control$caps$otsu, 1.0)
  # Only three buckets exist (GATE 6.1).
  expect_setequal(names(cfg$train_control$caps),
                  c("contextual", "random", "otsu"))
  expect_equal(cfg$train_control$training_protocol, "nested_refit")
  expect_equal(cfg$train_control$oof_sampling, "capped")
  prov <- cfg$resolved_params_provenance$train_control
  expect_equal(prov$cap_contextual, "user")
})

test_that("canonical OOF call is warning-free and the builder cap reaches the OOF engine", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_canonical_cfg()
  oof_ctx <- NULL
  expect_no_warning(
    oof_ctx <- capture_oof_contextual(cfg, mk_pb_oof_feats()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(oof_ctx, 1.0)
})

test_that("canonical FINAL call is warning-free and the builder cap reaches the FINAL engine", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_canonical_cfg()
  final_ctx <- NULL
  expect_no_warning(
    final_ctx <- capture_final_contextual(cfg, mk_pb_final_gpkg()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(final_ctx, 1.0)
})

test_that("the cap reaching OOF and FINAL is identical (single source of truth)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_canonical_cfg()
  oof_ctx   <- capture_oof_contextual(cfg, mk_pb_oof_feats())
  final_ctx <- capture_final_contextual(cfg, mk_pb_final_gpkg())
  expect_equal(oof_ctx, final_ctx)
  expect_equal(oof_ctx, 1.0)
})
