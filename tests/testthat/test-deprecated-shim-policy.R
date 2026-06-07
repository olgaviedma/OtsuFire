# Precision 1 (2026-06-07): deprecated function-level methodological-parameter
# SHIMS. The canonical path is build_supervised_burned_config(); function-level
# overrides are deprecated compatibility shims that (a) propagate when used,
# (b) emit a deprecation warning of class "otsufire_deprecated_param" when they
# override a canonical default, (c) ERROR when they conflict with an explicit
# builder-user value, and (d) are recorded in a provenance record. The builder
# path stays warning-free.
#
# Precision 2 (2026-06-07): the package-level spectral-cap parity guard aborts
# when the spectral cap reaching OOF and FINAL is not the resolved cfg value.

ns <- asNamespace("OtsuFire")

mk_tif <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}
mk_gpkg <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}
mk_cfg <- function(...) {
  build_supervised_burned_config(
    scenario = "balanced", internal_decisions = mk_gpkg(),
    change_index = mk_tif(), target_year = 2025L, ...
  )
}
mk_oof_train_feats <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_sf(fire_uid = "a", class = "burned",
            fold_rep1 = 1L, fold_rep2 = 1L, geometry = sfc)
}
mk_final_train_gpkg <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc),
               g, layer = "train_features", quiet = TRUE, delete_dsn = TRUE)
  g
}

# Mock helpers that capture the spectral cap the engines receive (so we can prove
# a function-level override STILL propagates while it warns).
capture_oof_spectral <- function(cfg, train_feats, ...) {
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
                      out_dir = tempfile(), matrix_dir = tempfile(), ...)
  captured
}
capture_final_spectral <- function(cfg, train_gpkg, ...) {
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
                           out_dir = tempfile(), ...)
  captured
}

# ---------------------------------------------------------------------------
# (1) Builder-set overrides (canonical path) produce NO deprecation warning.
# ---------------------------------------------------------------------------
test_that("P1: builder-set override (canonical path) is warning-free and propagates", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg(cap_spectral = 2.0)
  # No deprecation warning when the override is set via the builder and consumed
  # without a function-level override.
  expect_no_warning(
    cap <- capture_oof_spectral(cfg, mk_oof_train_feats()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$spectral, 2.0)
  expect_no_warning(
    capf <- capture_final_spectral(cfg, mk_final_train_gpkg()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(capf$spectral, 2.0)
})

# ---------------------------------------------------------------------------
# (2) A function-level override emits the deprecation warning AND propagates.
# ---------------------------------------------------------------------------
test_that("P1: a function-level override warns (otsufire_deprecated_param) AND propagates (OOF)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg()  # canonical default spectral = 1.0
  cap <- NULL
  expect_warning(
    cap <- capture_oof_spectral(cfg, mk_oof_train_feats(),
                                spectral_hard_negative_to_burned_ratio = 2.0),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$spectral, 2.0)  # STILL propagates
})

test_that("P1: a function-level override warns AND propagates (FINAL)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg()
  cap <- NULL
  expect_warning(
    cap <- capture_final_spectral(cfg, mk_final_train_gpkg(),
                                  spectral_hard_negative_to_burned_ratio = 2.0),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$spectral, 2.0)
})

# ---------------------------------------------------------------------------
# (3) Explicit builder-user value + a CONFLICTING function-level override -> error.
# ---------------------------------------------------------------------------
test_that("P1: builder-user value + conflicting function-level override -> error", {
  skip_if_not_installed("sf")
  cfg <- mk_cfg(cap_spectral = 2.0)  # provenance "user"
  expect_error(
    run_oof_diagnostics(train_features = mk_oof_train_feats(),
                        scoring_features = mk_oof_train_feats(), config = cfg,
                        out_dir = tempfile(), matrix_dir = tempfile(),
                        spectral_hard_negative_to_burned_ratio = 3.0),
    regexp = "Conflicting supervised parameter"
  )
  expect_error(
    train_final_burned_model(train_features = mk_final_train_gpkg(), config = cfg,
                             out_dir = tempfile(),
                             spectral_hard_negative_to_burned_ratio = 3.0),
    regexp = "Conflicting supervised parameter"
  )
})

test_that("P1: builder-user value + EQUAL function-level override -> no error, no warning", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg(cap_spectral = 2.0)
  expect_no_warning(
    cap <- capture_oof_spectral(cfg, mk_oof_train_feats(),
                                spectral_hard_negative_to_burned_ratio = 2.0),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$spectral, 2.0)
})

# ---------------------------------------------------------------------------
# (4) Provenance record on the cfg + the resolver labels.
# ---------------------------------------------------------------------------
test_that("P1: cfg carries per-field provenance ('default' vs 'user')", {
  cfg_def  <- mk_cfg()
  cfg_user <- mk_cfg(cap_spectral = 2.0, nrounds_max = 10L)
  prov_def  <- cfg_def$resolved_params_provenance$train_control
  prov_user <- cfg_user$resolved_params_provenance$train_control
  expect_equal(prov_def$cap_spectral, "default")
  expect_equal(prov_def$nrounds_max, "default")
  expect_equal(prov_user$cap_spectral, "user")
  expect_equal(prov_user$nrounds_max, "user")
  expect_equal(prov_user$cap_random, "default")  # untouched stays default
})

test_that("P1: shim resolver labels the resolution provenance correctly", {
  f <- get(".of_resolve_methodological_shim", envir = ns)
  rec <- get(".of_new_shim_record", envir = ns)()
  # cfg-default, no override
  expect_equal(f(NULL, 1.0, "cap_spectral", cfg_provenance = "default",
                 record = rec), 1.0)
  # function-override over a default -> warns + applies
  expect_warning(
    v <- f(2.0, 1.0, "cap_spectral", cfg_provenance = "default", record = rec),
    class = "otsufire_deprecated_param")
  expect_equal(v, 2.0)
  # cfg-user, no override
  expect_equal(f(NULL, 2.0, "cap_spectral", cfg_provenance = "user",
                 record = rec), 2.0)
  df <- get(".of_shim_record_to_df", envir = ns)(rec)
  expect_true(all(c("cfg-default", "function-override", "cfg-user") %in%
                    df$provenance))
})

# ---------------------------------------------------------------------------
# (5) Precision 2: the package-level spectral-cap parity guard.
# ---------------------------------------------------------------------------
test_that("P2: spectral-cap parity guard passes when OOF=FINAL=cfg", {
  f <- get(".of_assert_spectral_cap_parity", envir = ns)
  expect_equal(f(2.0, 2.0, 2.0), 2.0)
  expect_equal(f(1.0, 1.0, 1.0), 1.0)
})

test_that("P2: a simulated silent reversion trips the parity guard", {
  f <- get(".of_assert_spectral_cap_parity", envir = ns)
  # OOF capped at 2.0 but FINAL silently reverted to 1.0.
  expect_error(f(2.0, 1.0, 2.0), regexp = "parity guard FAILED")
  # Both 1.0 but the resolved cfg expected 2.0 (silent reversion on both paths).
  expect_error(f(1.0, 1.0, 2.0), regexp = "parity guard FAILED")
})

test_that("P2: cfg with spectral=2.0 -> the value reaching BOTH OOF and FINAL is 2.0", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg(cap_spectral = 2.0)
  cap_o <- capture_oof_spectral(cfg, mk_oof_train_feats())
  cap_f <- capture_final_spectral(cfg, mk_final_train_gpkg())
  expect_equal(cap_o$spectral, 2.0)
  expect_equal(cap_f$spectral, 2.0)
  expect_equal(cap_o$spectral, cap_f$spectral)
  # And the guard agrees the capped path is consistent.
  f <- get(".of_assert_spectral_cap_parity", envir = ns)
  expect_silent(f(cap_o$spectral, cap_f$spectral,
                  cfg$train_control$caps$spectral))
})
