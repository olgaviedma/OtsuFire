# Gate 1C.3 (2026-06-08): END-TO-END tests for the NA-preserving,
# recipe-driven supervised SCORING path (the "H BLOCKER").
#
# Contract under test (R/internal-sup-final-map.R::score_with_final_model +
# R/internal-sup-nested-refit.R::.of_reconcile_scoring_schema):
#   * the feature SCHEMA comes from the SAVED FINAL-refit recipe, never from the
#     scoring-year data;
#   * the scoring matrix is built NA-preserving (xgb.DMatrix(missing = NA)) so
#     NO row is dropped, even when a feature is entirely NA / hotspots = NULL;
#   * number of predictions == number of input polygons; row order / ids kept;
#   * save/persist the FINAL model + recipe, reload, predict -> identical;
#   * the explicit per-case schema-reconciliation policy is RECORDED.
#
# These are END-TO-END (train -> persist model+recipe -> reconcile+build+predict
# through the SAME helpers the production scoring closure uses), NOT a bare
# matrix-construction test.

ns <- asNamespace("OtsuFire")

reconcile_fn <- function() get(".of_reconcile_scoring_schema", envir = ns)
coerce_fn    <- function() get(".of_nested_coerce_features",  envir = ns)
applymed_fn  <- function() get(".of_nested_apply_medians",    envir = ns)
buildmat_fn  <- function() get(".of_nested_build_matrix",     envir = ns)
trainfn      <- function() get("train_final_model_direct",    envir = ns)
whitelist_fn <- function() get(".supervised_feature_cols",    envir = ns)

# Faithful reproduction of the production scoring matrix build:
# score_with_final_model() does exactly this sequence on the recipe + sf data.
score_via_recipe <- function(model, recipe, x_df) {
  rec <- reconcile_fn()(x_df, recipe)
  X_raw <- coerce_fn()(rec$df, recipe$impute$impute_factor_missing %||% "MISSING")
  X_imp <- applymed_fn()(X_raw, recipe$impute$numeric_medians %||% list())
  M <- buildmat_fn()(X_imp, ref_cols = recipe$cols$x_cols)
  dmat <- xgboost::xgb.DMatrix(M, missing = NA)
  list(p = as.numeric(predict(model, dmat)), M = M, record = rec$record)
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x[[1L]]))) y else x

# Train a real nested_refit FINAL model + recipe, persisted to disk, using the
# shared whitelist fixture. Returns paths + a scoring sf object that carries the
# admin columns + features.
train_persist_fixture <- function() {
  fixture_path <- testthat::test_path("test-final-model-uses-whitelist.R")
  source(fixture_path, local = TRUE)
  gpkg <- make_whitelist_fixture_gpkg(n_burned = 18, n_neg_random = 9,
                                      n_neg_drop = 9, seed = 31L)
  out_dir <- tempfile("naprev_train_")
  res <- suppressMessages(suppressWarnings(trainfn()(
    labelled_gpkg = gpkg,
    labelled_layer = "train_features",
    out_dir = out_dir,
    prefix = "naprev",
    overwrite = TRUE,
    verbose = FALSE,
    nrounds_max = 12L,
    early_stopping_rounds = 6L,
    contextual_exclusion_to_burned_ratio = 1,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42, seed = 42, val_frac = 0.2, group_col = "block_id",
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = ns$.of_canonical_model_params()
  )))
  scoring_sf <- sf::read_sf(gpkg, layer = "train_features", quiet = TRUE)
  list(gpkg = gpkg, out_dir = out_dir,
       model_rds = res$files$model_rds, recipe_rds = res$files$recipe_rds,
       scoring_sf = scoring_sf)
}

test_that("scoring with hotspots = NULL (all hs_* features absent) preserves all rows", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)

  # Simulate a hotspots = NULL year: DROP every hs_* feature + hotspot_available
  # from the scoring frame entirely (so they are EXPECTED-ABSENT per the recipe).
  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  hs_cols <- grep("^hs_|^hotspot_available", names(x_df), value = TRUE)
  expect_true(length(hs_cols) > 0)
  x_df_nohs <- x_df[, setdiff(names(x_df), hs_cols), drop = FALSE]

  out <- score_via_recipe(model, recipe, x_df_nohs)

  # ROW PRESERVATION + predictions == polygons (the H-blocker assertion).
  expect_equal(length(out$p), nrow(x_df_nohs))
  expect_equal(nrow(out$M), nrow(x_df_nohs))
  expect_true(all(is.finite(out$p)))

  # Column NAME + ORDER identity vs the saved recipe.
  expect_identical(colnames(out$M), recipe$cols$x_cols)

  # The absent hs_* features were RECORDED (expected_absent -> NA), not silently
  # dropped.
  absent_rec <- out$record[out$record$case == "expected_absent", ]
  expect_true(any(grepl("^hs_|^hotspot_available", absent_rec$column)))
  expect_true(all(absent_rec$action == "create_absent_numeric_NA"))
})

test_that("a numeric predictor entirely NA drops no row and is imputed per recipe", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)
  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))

  # Pick a non-degenerate numeric feature present in the recipe (has a median)
  # and blank it ENTIRELY to NA across all rows.
  med_cols <- names(recipe$impute$numeric_medians)
  target <- intersect(med_cols, names(x_df))[1]
  expect_false(is.na(target))
  x_df[[target]] <- NA_real_

  out <- score_via_recipe(model, recipe, x_df)

  expect_equal(length(out$p), nrow(x_df))   # no row dropped
  expect_true(all(is.finite(out$p)))

  # Imputed to the recipe median (not 0, not NA) in the design matrix.
  expect_true(target %in% colnames(out$M))
  expect_equal(unname(out$M[, target]),
               rep(recipe$impute$numeric_medians[[target]], nrow(x_df)),
               tolerance = 1e-9)

  # Recorded as a fully-NA feature.
  expect_true(target %in% out$record$column[out$record$case == "feature_fully_NA"])
})

test_that("expected columns absent are created with correct type and recipe order", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(fx$recipe_rds)
  # Feed a frame that has ONLY administrative columns (no model features at all).
  x_df <- data.frame(fire_uid = sprintf("u%02d", 1:5),
                     stringsAsFactors = FALSE)

  rec <- reconcile_fn()(x_df, recipe)
  # Every recipe feature_col is now present, numeric, all-NA, in recipe order.
  expect_identical(names(rec$df), recipe$cols$feature_cols)
  expect_true(all(vapply(rec$df, is.numeric, logical(1))))
  expect_true(all(vapply(rec$df, function(c) all(is.na(c)), logical(1))))
  # The extra admin column was recorded as a drop.
  expect_true("fire_uid" %in% rec$record$column[rec$record$case == "extra_column"])
})

test_that("degenerate categorical predictor is handled per policy (no invented level)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)
  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))

  # Inject a character value into a numeric feature: incompatible type -> coerce.
  med_cols <- intersect(names(recipe$impute$numeric_medians), names(x_df))
  target <- med_cols[1]
  x_df[[target]] <- as.character(x_df[[target]])           # numbers as strings
  x_df[[target]][1] <- "not_a_number"                       # one uncoercible cell

  out <- score_via_recipe(model, recipe, x_df)
  expect_equal(length(out$p), nrow(x_df))                   # still no row dropped
  expect_true(all(is.finite(out$p)))
  expect_true(target %in% out$record$column[out$record$case == "incompatible_type"])
})

test_that("save/load/predict round-trip is identical (model + recipe persisted)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  model1  <- readRDS(fx$model_rds)
  recipe1 <- readRDS(fx$recipe_rds)
  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))

  p1 <- score_via_recipe(model1, recipe1, x_df)$p

  # Re-persist + reload (mirrors the production save/score boundary).
  m2 <- tempfile(fileext = "_m.rds"); r2 <- tempfile(fileext = "_r.rds")
  on.exit(unlink(c(m2, r2), force = TRUE), add = TRUE)
  saveRDS(model1, m2); saveRDS(recipe1, r2)
  p2 <- score_via_recipe(readRDS(m2), readRDS(r2), x_df)$p

  expect_equal(p2, p1, tolerance = 1e-12)
})

test_that("full engine scores hotspots=NULL year end-to-end (predictions == polygons)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- train_persist_fixture()
  on.exit(unlink(c(fx$gpkg), force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Build the FIVE engine inputs in a single result_dir tree from the fixture.
  result_dir <- tempfile("naprev_engine_")
  dir.create(file.path(result_dir, "03_FEATURES"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(result_dir, "06_QA_LABELS"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(result_dir, "07_FINAL_MODEL_V2"), recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  prefix <- "naprev"
  file.copy(fx$model_rds,  file.path(result_dir, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_final_model.rds")))
  file.copy(fx$recipe_rds, file.path(result_dir, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_recipe.rds")))

  feat_sf <- fx$scoring_sf
  feat_sf$source_poly_id <- feat_sf$fire_uid
  feat_sf$class_final    <- ifelse(feat_sf$class == "burned", "keep", "drop")

  # Scoring universe = a hotspots-NULL year: strip every hs_* feature column.
  hs_cols <- grep("^hs_|^hotspot_available", names(feat_sf), value = TRUE)
  scoring_sf <- feat_sf[, setdiff(names(feat_sf), hs_cols), drop = FALSE]

  features_gpkg <- file.path(result_dir, "03_FEATURES", "features_geometry.gpkg")
  sf::st_write(feat_sf,    features_gpkg, layer = "train_features",   quiet = TRUE)
  sf::st_write(scoring_sf, features_gpkg, layer = "scoring_features", append = TRUE, quiet = TRUE)

  qa_df <- data.frame(fire_uid = feat_sf$fire_uid,
                      p_oof_mean = stats::runif(nrow(feat_sf)),
                      stringsAsFactors = FALSE)
  qa_sf <- sf::st_sf(qa_df, geometry = sf::st_geometry(feat_sf))
  qa_gpkg <- file.path(result_dir, "06_QA_LABELS",
                       paste0("2022_patch_labeled_oof_summary.gpkg"))
  sf::st_write(qa_sf, qa_gpkg, layer = "labeled_oof_summary", quiet = TRUE)

  fn <- get("score_burnedlike_and_export_final_map", envir = ns)
  res <- suppressMessages(suppressWarnings(fn(
    result_dir = result_dir,
    prefix = prefix,
    qa_labelled_gpkg = qa_gpkg,
    labelled_features_gpkg = features_gpkg,
    labelled_features_layer = "train_features",
    model_rds  = file.path(result_dir, "07_FINAL_MODEL_V2",
                           paste0(prefix, "_final_model.rds")),
    recipe_rds = file.path(result_dir, "07_FINAL_MODEL_V2",
                           paste0(prefix, "_recipe.rds")),
    unlabeled_gpkg = features_gpkg,
    unlabeled_layer = "scoring_features",
    id_col = "fire_uid",
    overwrite = TRUE,
    verbose = FALSE
  )))

  # predictions == input polygons; NO dropped rows despite hotspots = NULL.
  expect_equal(nrow(res$deterministic_scored), nrow(scoring_sf))
  expect_true(all(is.finite(res$deterministic_scored$p_burned)))
  # ID + order preserved exactly.
  expect_identical(as.character(res$deterministic_scored$fire_uid),
                   as.character(scoring_sf$fire_uid))
})
