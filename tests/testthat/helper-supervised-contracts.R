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
# cap_contextual = 1.0, nrounds_max = 10L) can be threaded through. target_year
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

# ============================================================================
# Gate 1D.8 (2026-06-09): SHARED builders for the OOF/FINAL/scoring feature-recipe
# parity CONTRACT suite (test-oof-final-isna-parity.R + the 18-point contract
# test-oof-final-recipe-parity-18.R). Promoted here from the parity test file so
# BOTH contract files consume ONE definition (testthat 3e isolates each test
# file's top-level environment; helper-*.R is the only shared scope).
#
# Every builder takes a RAW upstream frame with NO `_isNA` columns pre-seeded
# (the analogue of one extract_supervised_features() GPKG feeding both stages),
# and resolves the DEPLOYED feature/column set each stage actually trains /
# scores on, through the SAME internal helpers the production pipeline uses.
# ============================================================================

# Internal-namespace accessors (the contract files reach into non-exported
# helpers / constants to drive each stage exactly as production does).
sc_ns      <- function() asNamespace("OtsuFire")
sc_ip      <- function(nm) get(nm, envir = sc_ns())

# `_isNA` / base column partition utilities, shared by every parity assertion.
sc_isna_of <- function(v) sort(grep("_isNA$", v, value = TRUE))
sc_base_of <- function(v) v[!grepl("_isNA$", v)]

# One upstream feature frame == the analogue of a single
# extract_supervised_features() GPKG feeding BOTH stages. Carries the canonical
# whitelist features (some with injected NA so `_isNA` companions are
# meaningful) + the admin / fold-assignment columns both stages preserve. NO
# `_isNA` columns are pre-written, so each stage's synthesis behaviour is what's
# under test.
#
# @param na_targets character; whitelist features into which NA is injected (so
#   the synthesised companions are non-degenerate). Defaults to a spread of
#   RBR / topo / hotspot features.
# @param all_na character; whitelist features forced ENTIRELY NA (point 7: an
#   all-NA feature -> `_isNA` all-1, base per missing policy, no row dropped).
sc_isna_upstream_df <- function(n = 60L, seed = 7L,
                                na_targets = c("rbr_med", "elev_med", "slope_med",
                                               "hs_conf_mean", "hs_frp_sum"),
                                all_na = character(0)) {
  set.seed(seed)
  wl <- sc_ip(".supervised_feature_cols")
  feat <- as.data.frame(
    matrix(stats::runif(n * length(wl)), nrow = n,
           dimnames = list(NULL, wl)))
  for (t in intersect(na_targets, names(feat))) {
    feat[[t]][sample.int(n, max(1L, floor(n / 5L)))] <- NA
  }
  for (t in intersect(all_na, names(feat))) feat[[t]] <- NA_real_
  admin <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = rep(c("burned", "unburned"), length.out = n),
    source    = rep(c("burned_truth", "random_burnable_background",
                      "deterministic_drop_hard", "otsu_patch_residual"),
                    length.out = n),
    neg_type  = rep(c(NA_character_, "background_cell",
                      "geo_excluded_hot", "otsu_patch_drop"),
                    length.out = n),
    poly_id   = sprintf("p_%03d", seq_len(n)),
    block_id  = rep(seq_len(12), length.out = n),
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2 = rep(c(2L, 3L, 1L), length.out = n),
    stringsAsFactors = FALSE)
  cbind(admin, feat)
}

# Resolve the OOF deployed design-matrix model_cols (what the OOF per-fold core
# trains on) from one upstream frame, replicating the OOF wrapper's whitelist
# filter + deferred design-matrix build.
sc_oof_model_cols_from <- function(df) {
  wl <- sc_ip(".supervised_feature_cols")
  allowed <- c(wl, paste0(wl, "_isNA"))
  id_cols <- c("fire_uid", "poly_id", "block_id", "fold_rep1", "fold_rep2",
               "class", "source", "neg_type")
  apply_wl <- function(d) {
    keep <- unique(c(intersect(names(d), id_cols),
                     intersect(names(d), allowed)))
    d[, keep, drop = FALSE]
  }
  labelled    <- apply_wl(df)
  burned_like <- labelled[labelled$class == "unburned", , drop = FALSE]
  dm <- sc_ip("build_design_matrix_patches")(
    labelled = labelled, burned_like = burned_like, id_cols = id_cols,
    defer_impute = TRUE, save_dir = NULL, verbose = FALSE)
  dm$model_cols
}

# Resolve the FINAL deployed design-matrix recipe (x_cols + recipe partition)
# from the SAME upstream frame, replicating train_final_model_direct()'s
# feat_cols derivation + the shared nested-refit core. Returns the full fit list.
sc_final_fit_from <- function(df) {
  wl <- sc_ip(".supervised_feature_cols")
  feat_cols <- sc_ip(".filter_to_supervised_whitelist")(names(df), whitelist = wl)
  params_fn <- function(spw) {
    p <- sc_ip(".of_canonical_xgb_params")(spw); p$nthread <- 1L; p
  }
  suppressMessages(suppressWarnings(sc_ip(".of_nested_refit_fit")(
    train_df = df, feature_cols = feat_cols, label_col = "class",
    group_col = "block_id", block_col = "block_id", val_frac = 0.2,
    params_fn = params_fn, sampling_seed = 42L, fold_seed = 42L,
    nrounds_max = 8L, early_stopping_rounds = 4L,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    verbose = FALSE)))
}
sc_final_x_cols_from <- function(df) sc_final_fit_from(df)$x_cols

# Build the SAVED-recipe shape (the persisted FINAL recipe schema) from a fit.
sc_recipe_from_fit <- function(fit) {
  medians <- fit$medians[!vapply(fit$medians,
                                 function(z) is.null(z) || is.na(z), logical(1))]
  list(
    cols = list(id_col = "fire_uid", class_col = "class", group_col = "block_id",
                feature_cols = fit$feature_cols,
                base_features = fit$base_features,
                missing_indicator_features = fit$missing_indicator_features,
                final_feature_order = fit$final_feature_order,
                x_cols = fit$x_cols),
    impute = list(impute_numeric = "median", numeric_medians = medians,
                  impute_factor_missing = "MISSING"))
}

# Resolve the SCORING deployed x_cols from the SAME upstream frame: build the
# saved-recipe shape from the FINAL fit, then drive the production scoring build
# (.of_reconcile_scoring_schema -> .of_nested_coerce_features ->
# .of_nested_apply_medians -> .of_nested_build_matrix(ref = x_cols)) on the raw
# (no `_isNA`) feature frame.
sc_scoring_x_cols_from <- function(df, fit) {
  recipe  <- sc_recipe_from_fit(fit)
  medians <- recipe$impute$numeric_medians
  wl <- sc_ip(".supervised_feature_cols")
  x_df <- df[, intersect(names(df), c("fire_uid", wl)), drop = FALSE]
  rec <- sc_ip(".of_reconcile_scoring_schema")(x_df, recipe)
  X_raw <- sc_ip(".of_nested_coerce_features")(rec$df, "MISSING")
  X_imp <- sc_ip(".of_nested_apply_medians")(X_raw, medians)
  M <- sc_ip(".of_nested_build_matrix")(X_imp, ref_cols = recipe$cols$x_cols)
  colnames(M)
}

# Faithful reproduction of the production scoring matrix build + predict
# (score_with_final_model does exactly this sequence on the recipe + sf data).
# Returns the predictions, the matrix and the reconciliation record.
sc_score_via_recipe <- function(model, recipe, x_df) {
  ip <- sc_ip
  `%or%` <- function(a, b) if (is.null(a)) b else a
  rec   <- ip(".of_reconcile_scoring_schema")(x_df, recipe)
  X_raw <- ip(".of_nested_coerce_features")(
    rec$df, recipe$impute$impute_factor_missing %or% "MISSING")
  X_imp <- ip(".of_nested_apply_medians")(
    X_raw, recipe$impute$numeric_medians %or% list())
  M <- ip(".of_nested_build_matrix")(X_imp, ref_cols = recipe$cols$x_cols)
  dmat <- xgboost::xgb.DMatrix(M, missing = NA)
  list(p = as.numeric(stats::predict(model, dmat)), M = M, record = rec$record)
}

# A REAL, tiny nested_refit FINAL model + recipe persisted to disk, built from a
# fresh on-disk GPKG whose train_features layer carries the whitelist features +
# admin columns (NO `_isNA` pre-seeded). Reuses the whitelist-fixture GPKG
# builder. Returns model/recipe paths + the scoring sf (admin + features).
# Keep it tiny (few rows, nthread via canonical params) for points 13/14/15/16.
sc_train_persist_fit <- function(prefix = "sc18",
                                 n_burned = 18L, n_neg_random = 9L,
                                 n_neg_drop = 9L, seed = 31L) {
  fixture_path <- testthat::test_path("test-final-model-uses-whitelist.R")
  source(fixture_path, local = TRUE)
  gpkg <- make_whitelist_fixture_gpkg(n_burned = n_burned,
                                      n_neg_random = n_neg_random,
                                      n_neg_drop = n_neg_drop, seed = seed)
  out_dir <- tempfile("sc18_train_")
  res <- suppressMessages(suppressWarnings(sc_ip("train_final_model_direct")(
    labelled_gpkg = gpkg, labelled_layer = "train_features",
    out_dir = out_dir, prefix = prefix, overwrite = TRUE, verbose = FALSE,
    nrounds_max = 12L, early_stopping_rounds = 6L,
    contextual_exclusion_to_burned_ratio = 1,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42, seed = 42, val_frac = 0.2, group_col = "block_id",
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = sc_ip(".of_canonical_model_params")())))
  scoring_sf <- sf::read_sf(gpkg, layer = "train_features", quiet = TRUE)
  list(gpkg = gpkg, out_dir = out_dir, res = res,
       model_rds = res$files$model_rds, recipe_rds = res$files$recipe_rds,
       scoring_sf = scoring_sf)
}
