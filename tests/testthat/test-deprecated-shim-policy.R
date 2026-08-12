# Precision 1 (2026-06-07): deprecated function-level methodological-parameter
# SHIMS. The canonical path is build_supervised_burned_config(); function-level
# overrides are deprecated compatibility shims that (a) propagate when used,
# (b) emit a deprecation warning of class "otsufire_deprecated_param" when they
# override a canonical default, (c) ERROR when they conflict with an explicit
# builder-user value, and (d) are recorded in a provenance record. The builder
# path stays warning-free.
#
# GATE 6.1 (2026-06-11): the dead spectral negative bucket was removed.
# GATE 6.5 (2026-06-12): the contextual (deterministic-drop) bucket was removed.
# The shim policy is exercised here on the RANDOM cap, which remains an operative
# bucket with a deprecated function-level shim (random_to_burned_ratio).

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
    run_label = "balanced", internal_decisions = mk_gpkg(),
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

# Mock helpers that capture the RANDOM cap the engines receive (so we can prove a
# function-level override STILL propagates while it warns).
capture_oof_random <- function(cfg, train_feats, ...) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      a <- list(...)
      captured$random <- a$random_to_burned_ratio
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL), files = list())
    },
    .package = "OtsuFire"
  )
  run_oof_diagnostics(train_features = train_feats,
                      scoring_features = train_feats, config = cfg,
                      out_dir = tempfile(), matrix_dir = tempfile(), ...)
  captured
}
capture_final_random <- function(cfg, train_gpkg, ...) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      a <- list(...)
      captured$random <- a$random_to_burned_ratio
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
  # GATE 6.7 (2026-06-12): caps via the typed negative_pool_params block.
  cfg <- mk_cfg(negative_pool_params = list(caps = c(random = 0.5, otsu = 1.0)))
  expect_no_warning(
    cap <- capture_oof_random(cfg, mk_oof_train_feats()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$random, 0.5)
  expect_no_warning(
    capf <- capture_final_random(cfg, mk_final_train_gpkg()),
    class = "otsufire_deprecated_param"
  )
  expect_equal(capf$random, 0.5)
})

# ---------------------------------------------------------------------------
# (2) A function-level override emits the deprecation warning AND propagates.
# ---------------------------------------------------------------------------
test_that("P1: a function-level override warns (otsufire_deprecated_param) AND propagates (OOF)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg()  # canonical default random = 1.0
  cap <- NULL
  expect_warning(
    cap <- capture_oof_random(cfg, mk_oof_train_feats(),
                              random_to_burned_ratio = 0.5),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$random, 0.5)  # STILL propagates
})

test_that("P1: a function-level override warns AND propagates (FINAL)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg()
  cap <- NULL
  expect_warning(
    cap <- capture_final_random(cfg, mk_final_train_gpkg(),
                                random_to_burned_ratio = 0.5),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$random, 0.5)
})

# ---------------------------------------------------------------------------
# (3) Explicit builder-user value + a CONFLICTING function-level override -> error.
# ---------------------------------------------------------------------------
test_that("P1: builder-user value + conflicting function-level override -> error", {
  skip_if_not_installed("sf")
  cfg <- mk_cfg(negative_pool_params =
                  list(caps = c(random = 1.0, otsu = 1.0)))  # provenance "user"
  expect_error(
    run_oof_diagnostics(train_features = mk_oof_train_feats(),
                        scoring_features = mk_oof_train_feats(), config = cfg,
                        out_dir = tempfile(), matrix_dir = tempfile(),
                        random_to_burned_ratio = 0.5),
    regexp = "Conflicting supervised parameter"
  )
  expect_error(
    train_final_burned_model(train_features = mk_final_train_gpkg(), config = cfg,
                             out_dir = tempfile(),
                             random_to_burned_ratio = 0.5),
    regexp = "Conflicting supervised parameter"
  )
})

test_that("P1: builder-user value + EQUAL function-level override -> no error, no warning", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_cfg(negative_pool_params = list(caps = c(random = 1.0, otsu = 1.0)))
  expect_no_warning(
    cap <- capture_oof_random(cfg, mk_oof_train_feats(),
                              random_to_burned_ratio = 1.0),
    class = "otsufire_deprecated_param"
  )
  expect_equal(cap$random, 1.0)
})

# ---------------------------------------------------------------------------
# (4) Provenance record on the cfg + the resolver labels.
# ---------------------------------------------------------------------------
test_that("P1: cfg carries per-field provenance ('default' vs 'user')", {
  cfg_def  <- mk_cfg()
  cfg_user <- mk_cfg(
    negative_pool_params = list(caps = c(random = 0.5, otsu = 1.0)),
    nrounds_max = 10L)
  prov_def  <- cfg_def$resolved_params_provenance$train_control
  prov_user <- cfg_user$resolved_params_provenance$train_control
  expect_equal(prov_def$cap_random, "default")
  expect_equal(prov_def$cap_otsu, "default")
  expect_equal(prov_def$nrounds_max, "default")
  expect_equal(prov_user$cap_random, "user")
  expect_equal(prov_user$nrounds_max, "user")
  # GATE 6.7 (2026-06-12): caps are supplied as a single named vector via
  # negative_pool_params$caps, so providing it marks BOTH bucket caps "user".
  expect_equal(prov_user$cap_otsu, "user")
})

test_that("P1: shim resolver labels the resolution provenance correctly", {
  f <- get(".of_resolve_methodological_shim", envir = ns)
  rec <- get(".of_new_shim_record", envir = ns)()
  # cfg-default, no override
  expect_equal(f(NULL, 1.0, "cap_random", cfg_provenance = "default",
                 record = rec), 1.0)
  # function-override over a default -> warns + applies
  expect_warning(
    v <- f(0.5, 1.0, "cap_random", cfg_provenance = "default", record = rec),
    class = "otsufire_deprecated_param")
  expect_equal(v, 0.5)
  # cfg-user, no override
  expect_equal(f(NULL, 0.5, "cap_random", cfg_provenance = "user",
                 record = rec), 0.5)
  df <- get(".of_shim_record_to_df", envir = ns)(rec)
  expect_true(all(c("cfg-default", "function-override", "cfg-user") %in%
                    df$provenance))
})

# ---------------------------------------------------------------------------
# (5) GATE 6.1 / 6.5: the spectral + contextual caps / their guards no longer
#     exist.
# ---------------------------------------------------------------------------
test_that("GATE 6.1/6.5: spectral + contextual cap args (and the parity guard) are gone", {
  expect_false(exists(".of_assert_spectral_cap_parity", envir = ns,
                      inherits = FALSE))
  expect_error(mk_cfg(cap_spectral = 2.0), regexp = "cap_spectral|unused argument")
  expect_error(mk_cfg(cap_contextual = 1.0),
               regexp = "cap_contextual|unused argument")
})
