# ============================================================================
# Shared fixtures for the supervised CONTRACT-TEST layer (Gate 1C / 1D).
#
# Gate 1D.7 (2026-06-08): several Gate-1C/1D contract files independently
# re-built the SAME tiny fixtures (a 4x4 change-index .tif, a one-polygon
# decisions .gpkg, a minimal resolved cfg, a one-row burned train sf / gpkg,
# and the cfg->engine capture-mock pattern). This helper centralises those so
# the contract files share ONE definition without weakening any assertion.
#
# Naming convention: every exported fixture here is prefixed `sc_` (supervised
# contracts) so it never collides with a file-local fixture. Helpers are loaded
# automatically by testthat (helper-*.R) before any test file runs.
#
# IMPORTANT: these are *fixtures only*. They construct inputs; they assert
# nothing. The contracts themselves live in the individual test files.
# ============================================================================
#
# ----------------------------------------------------------------------------
# SUPERVISED CONTRACT INDEX (Gate 1C/1D). The 9 contracts from the Gate-1C/1D
# closure list and the test file that OWNS each (the authoritative stable
# assertion). Cross-cutting contracts list the primary owner first.
#
#   #  contract                          owning test file
#   -- --------------------------------- ------------------------------------
#   1  mask / CRS alignment              test-align-mask-crs.R
#   2  B4 sampling in burnable domain    test-b4-burnable-domain.R
#   3  all-NA features preserve rows     test-scoring-na-preserving.R
#   4  no-hotspots (hs_* absent) year    test-scoring-na-preserving.R
#   5  row / id / order preservation     test-scoring-na-preserving.R
#                                          (+ reproducibility Level C:
#                                           test-reproducibility-contract.R)
#   6  save / load / predict round-trip  test-scoring-na-preserving.R
#   7  cfg isolation (no shared state)   test-supervised-validate-execution.R
#   8  fingerprints WITHOUT timestamps   test-reproducibility-contract.R
#                                          (+ run/pool fp:
#                                           test-supervised-validate-execution.R;
#                                           B4 neg-pool fp: test-b4-burnable-domain.R)
#   9  content-aware cache               test-overwrite-contract.R
#                                          (+ vector-input cache:
#                                           test-validate-fire-maps.R)
#
# Adjacent supervised contract suites (not in the 9 but part of the layer):
#   * cfg single-source of truth         test-cfg-single-source.R
#   * negative-pool caps (cfg->engines)  test-caps-contract.R
#   * overwrite end-to-end               test-overwrite-contract.R
#   * OOF<->FINAL recipe/param symmetry  test-oof-final-symmetry.R
#   * OOF<->FINAL `_isNA` parity (1D.7)  test-oof-final-isna-parity.R
#   * fail-fast validation (1-10)        test-supervised-validate-execution.R
# ----------------------------------------------------------------------------

# A 4x4 single-band change-index raster on disk (EPSG-agnostic; values 1..16).
# Used by the cfg-builder contracts that only need a readable raster path.
sc_change_index_tif <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}

# A one-polygon decisions .gpkg on disk (EPSG:3035). Used by the cfg-builder
# contracts that only need a readable vector path.
sc_decisions_gpkg <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}

# A minimal resolved supervised cfg built from fresh on-disk inputs. Forwards
# `...` to build_supervised_burned_config() so a single override (e.g.
# cap_spectral = 2.0, nrounds_max = 10L) can be threaded through. target_year
# defaults to 2025L (the value the caps / single-source suites used).
sc_min_cfg <- function(...) {
  build_supervised_burned_config(
    scenario = "balanced", internal_decisions = sc_decisions_gpkg(),
    change_index = sc_change_index_tif(), target_year = 2025L, ...
  )
}

# A one-row burned training sf (single polygon, EPSG:3035) carrying the
# fold-assignment metadata the OOF wrapper preserves. This is the in-memory
# `train_features` fixture for the cfg->OOF capture contracts.
sc_burned_train_sf <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_sf(fire_uid = "a", class = "burned",
            fold_rep1 = 1L, fold_rep2 = 1L, geometry = sfc)
}

# The same one-row burned training frame written to a `train_features` layer in
# a fresh .gpkg. This is the on-disk fixture for the cfg->FINAL capture
# contracts (train_final_model_direct reads a gpkg).
sc_burned_train_gpkg <- function() {
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc),
               g, layer = "train_features", quiet = TRUE, delete_dsn = TRUE)
  g
}
