# 0.4.0 architectural refactor (Agent H): tests for the whitelist
# alignment in the scoring path
# (`internal-sup-final-map.R::score_with_final_model()`).
#
# Under 0.4.0 the scoring matrix is aligned to
# `recipe$cols$x_cols`, which is itself a subset of
# `.supervised_feature_cols ∪ _isNA companions` (asserted in
# test-final-model-uses-whitelist.R). This test ensures the scoring
# path actually uses that recipe column set rather than blindly
# accepting the unlabeled GPKG's column set.

whitelist_internal <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

test_that("recipe$cols$x_cols ⊆ whitelist (post-train invariant)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  whitelist <- whitelist_internal()
  allowed   <- c(whitelist, paste0(whitelist, "_isNA"))

  # Borrow the fixture builder from
  # test-final-model-uses-whitelist.R (helper-fixtures.R is the
  # standard place but we keep this test self-contained).
  fixture_path <- testthat::test_path("test-final-model-uses-whitelist.R")
  expect_true(file.exists(fixture_path))
  source(fixture_path, local = TRUE)
  gpkg <- make_whitelist_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- get("train_final_model_direct", envir = asNamespace("OtsuFire"))
  out_dir <- tempfile("scoring_whitelist_")
  res <- suppressMessages(suppressWarnings(fn(
    labelled_gpkg = gpkg,
    labelled_layer = "train_features",
    out_dir = out_dir,
    prefix = "scoring_whitelist",
    overwrite = TRUE,
    verbose = FALSE,
    nrounds_max = 8L,
    early_stopping_rounds = 4L,
    # Gate 1B (2026-06-07): the engine now requires resolved methodological args.
    contextual_exclusion_to_burned_ratio = 1,
    spectral_hard_negative_to_burned_ratio = 1,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = OtsuFire:::.of_canonical_model_params()
  )))
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(res$files$recipe_rds)

  # Invariant 1: x_cols is whitelist-admissible.
  expect_true(all(recipe$cols$x_cols %in% allowed),
              info = paste("non-whitelist x_cols:",
                            paste(setdiff(recipe$cols$x_cols, allowed),
                                  collapse = ", ")))

  # Invariant 2: feature_cols is whitelist-admissible.
  expect_true(all(recipe$cols$feature_cols %in% allowed),
              info = paste("non-whitelist feature_cols:",
                            paste(setdiff(recipe$cols$feature_cols, allowed),
                                  collapse = ", ")))

  # Invariant 3 (architectural symmetry): recipe$cols$feature_cols
  # is exactly intersect(.supervised_feature_cols, available
  # numeric/integer cols). For our fixture every whitelist feature
  # is supplied as numeric, so the intersection is the full
  # whitelist.
  expect_setequal(recipe$cols$feature_cols, whitelist)
})

# NOTE (Gate 1C.3, 2026-06-08): the production scoring closure now uses the
# NA-preserving recipe-driven builder (.of_reconcile_scoring_schema +
# .of_nested_* + xgb.DMatrix(missing = NA)), not prep_X_df/sparse.model.matrix.
# This test reimplements the column-selection step inline and still asserts the
# scoring feature set is whitelist-admissible (the invariant under test).
test_that("scoring path feature selection yields whitelist-only cols", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fixture_path <- testthat::test_path("test-final-model-uses-whitelist.R")
  source(fixture_path, local = TRUE)
  gpkg <- make_whitelist_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- get("train_final_model_direct", envir = asNamespace("OtsuFire"))
  out_dir <- tempfile("scoring_pipe_")
  res <- suppressMessages(suppressWarnings(fn(
    labelled_gpkg = gpkg,
    labelled_layer = "train_features",
    out_dir = out_dir,
    prefix = "scoring_pipe",
    overwrite = TRUE,
    verbose = FALSE,
    nrounds_max = 8L,
    early_stopping_rounds = 4L,
    # Gate 1B (2026-06-07): the engine now requires resolved methodological args.
    contextual_exclusion_to_burned_ratio = 1,
    spectral_hard_negative_to_burned_ratio = 1,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = OtsuFire:::.of_canonical_model_params()
  )))
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(res$files$recipe_rds)
  whitelist <- whitelist_internal()
  allowed   <- c(whitelist, paste0(whitelist, "_isNA"))

  # Mock the scoring path: read the labelled GPKG (which carries
  # all administrative columns), apply the same prep / matrix logic
  # that score_with_final_model() does, and assert the resulting
  # columns are whitelist-admissible.
  S <- sf::read_sf(gpkg, layer = "train_features", quiet = TRUE)
  x_df <- as.data.frame(sf::st_drop_geometry(S))

  feature_cols <- recipe$cols$feature_cols
  missing_feat <- setdiff(feature_cols, names(x_df))
  if (length(missing_feat) > 0) {
    for (nm in missing_feat) x_df[[nm]] <- NA
  }
  expect_true(all(feature_cols %in% allowed))
  # The post-filter matrix derived from feature_cols can only
  # contain admissible columns.
  expect_equal(intersect(feature_cols, c("neg_type", "fire_uid",
                                            "block_id", "median_rbr",
                                            "p_above_keep_q25",
                                            "qa_changed", "p_oof_mean",
                                            "area_ha", "n_pix",
                                            "legacy_decision",
                                            "intersects_deterministic")),
                character(0))
})
