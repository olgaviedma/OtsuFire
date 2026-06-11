# =============================================================================
# Gate 1B (2026-06-07): the OFFICIAL B1 Phase 2 runner configures methodological
# params via the cfg BUILDER (no deprecated function-level shims), so the
# canonical Phase-B run is WARNING-FREE.
#
# This test constructs the Phase-B cfg EXACTLY as the runner does -- by sourcing
# the SAME cfg-construction helper the runner sources (B1_PHASE2_CFG_HELPER.R) --
# and asserts:
#   (1) build_b1_phase2_cfg("nested_refit", "capped") -> cfg$train_control$caps$
#       spectral == 2.0 with builder provenance "user";
#   (2) calling run_oof_diagnostics() / train_final_burned_model() with that cfg
#       and NO function-level methodological override produces NO warning of
#       class "otsufire_deprecated_param" (the canonical path is warning-free);
#   (3) the spectral cap reaching BOTH the OOF engine and the FINAL engine is
#       2.0 (the builder value flows through both stages identically).
#
# The deprecated-shim path (function-level overrides that DO warn) is tested
# separately in test-deprecated-shim-policy.R; this file is the canonical
# builder-path counterpart for the actual runner's cfg construction.
# =============================================================================

ns <- asNamespace("OtsuFire")

# The runner + helper live under the USAGE scripts tree (a sibling of the
# package). Resolve the helper from the canonical absolute path; skip cleanly if
# it is not present (e.g. the package is being tested outside the project tree).
.b1_helper <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/00_USAGE/02_SUPERVISED_USAGE",
  "B1_PHASE2_CFG_HELPER.R"
)

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

# Build the Phase-B (Config B: nested_refit / capped) cfg the SAME way the
# runner does. We source the runner's helper into a private environment and call
# build_b1_phase2_cfg() with throwaway inputs (the methodological knobs are what
# matters here, not the run inputs).
build_runner_phaseB_cfg <- function() {
  helper_env <- new.env(parent = environment())
  # build_supervised_burned_config must be visible to the sourced helper.
  assign("build_supervised_burned_config", build_supervised_burned_config,
         envir = helper_env)
  sys.source(.b1_helper, envir = helper_env)
  helper_env$build_b1_phase2_cfg(
    oof_sampling      = "capped",
    scenario          = "balanced",
    internal_decisions = mk_pb_gpkg(),
    change_index       = mk_pb_tif(),
    target_year        = 2017L
  )
}

# Capture the spectral cap each engine receives, with the engines mocked so no
# real fit runs (mirrors the lightweight pattern in test-cfg-single-source.R).
capture_oof_spectral <- function(cfg, train_feats) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      a <- list(...)
      captured$spectral <- a$spectral_hard_negative_to_burned_ratio
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL), files = list())
    },
    .package = "OtsuFire"
  )
  run_oof_diagnostics(train_features = train_feats,
                      scoring_features = train_feats, config = cfg,
                      out_dir = tempfile(), matrix_dir = tempfile())
  captured$spectral
}
capture_final_spectral <- function(cfg, train_gpkg) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      a <- list(...)
      captured$spectral <- a$spectral_hard_negative_to_burned_ratio
      list(model = NULL, recipe = NULL, training_ok_sf = NULL, split = NULL,
           files = list())
    },
    .package = "OtsuFire"
  )
  train_final_burned_model(train_features = train_gpkg, config = cfg,
                           out_dir = tempfile())
  captured$spectral
}

# ---------------------------------------------------------------------------
# (1) The runner's Phase-B cfg carries spectral cap 2.0, provenance "user".
# ---------------------------------------------------------------------------
test_that("Phase B runner cfg: spectral cap 2.0 with provenance 'user'", {
  skip_if_not(file.exists(.b1_helper),
              "B1_PHASE2_CFG_HELPER.R not found in project tree")
  cfg <- build_runner_phaseB_cfg()

  expect_equal(cfg$train_control$caps$spectral, 2.0)
  # the other Phase-B caps are configured via the builder too.
  expect_equal(cfg$train_control$caps$contextual, 0.25)
  expect_equal(cfg$train_control$caps$random, 1.0)
  expect_equal(cfg$train_control$caps$otsu, 1.0)
  # sampling is a Config B axis, set via the builder. training_protocol is now a
  # fixed internal constant (single training procedure, not a builder argument).
  expect_equal(cfg$train_control$training_protocol, "nested_refit")
  expect_equal(cfg$train_control$oof_sampling, "capped")
  # provenance: the cap came from the builder (user), not a default.
  prov <- cfg$resolved_params_provenance$train_control
  expect_equal(prov$cap_spectral, "user")
})

# ---------------------------------------------------------------------------
# (2)+(3) The canonical call sites (no function-level override) are warning-free
#         AND the spectral cap 2.0 reaches BOTH the OOF and FINAL engines.
# ---------------------------------------------------------------------------
test_that("Phase B runner: canonical OOF call is warning-free and 2.0 reaches the OOF engine", {
  skip_if_not(file.exists(.b1_helper),
              "B1_PHASE2_CFG_HELPER.R not found in project tree")
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_runner_phaseB_cfg()

  oof_spectral <- NULL
  expect_no_warning(
    oof_spectral <- capture_oof_spectral(cfg, mk_pb_oof_feats()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(oof_spectral, 2.0)
})

test_that("Phase B runner: canonical FINAL call is warning-free and 2.0 reaches the FINAL engine", {
  skip_if_not(file.exists(.b1_helper),
              "B1_PHASE2_CFG_HELPER.R not found in project tree")
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_runner_phaseB_cfg()

  final_spectral <- NULL
  expect_no_warning(
    final_spectral <- capture_final_spectral(cfg, mk_pb_final_gpkg()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(final_spectral, 2.0)
})

test_that("Phase B runner: the spectral cap reaching OOF and FINAL is identical (2.0) and passes the parity guard", {
  skip_if_not(file.exists(.b1_helper),
              "B1_PHASE2_CFG_HELPER.R not found in project tree")
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- build_runner_phaseB_cfg()

  oof_spectral   <- capture_oof_spectral(cfg, mk_pb_oof_feats())
  final_spectral <- capture_final_spectral(cfg, mk_pb_final_gpkg())
  expect_equal(oof_spectral, final_spectral)
  expect_equal(oof_spectral, 2.0)

  # The package-level parity guard agrees the capped path is consistent.
  guard <- get(".of_assert_spectral_cap_parity", envir = ns)
  expect_silent(guard(oof_spectral, final_spectral,
                      cfg$train_control$caps$spectral))
})
