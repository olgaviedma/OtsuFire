# B1 Phase 2 (2026-06-07): save/load round-trip for the nested_refit core.
#
# Goal: prove the DEPLOYED/SCORED model == the refit core model. We train a
# small nested_refit model on a synthetic fixture, get the DIRECT core
# prediction on a held-out matrix, then SAVE the model + recipe (refit medians,
# x_cols order, _isNA handling via the transform, best_iteration) to tempfiles,
# RELOAD, re-apply the recipe to the SAME held-out rows, and assert the reloaded
# prediction EQUALS the direct prediction bit-for-bit.
#
# This is the unit-level analogue of the production save/score path: the FINAL
# stage persists model + recipe via saveRDS and the scoring stage reloads them;
# this test confirms that round-trip preserves predictions for the nested core.

ns <- asNamespace("OtsuFire")

nr_fit       <- function() get(".of_nested_refit_fit", envir = ns)
nr_transform <- function() get(".of_nested_transform", envir = ns)
canon_params <- function() get(".of_canonical_xgb_params", envir = ns)

# Build a synthetic train fixture + a held-out (outer-test-like) fixture whose
# first row has a MISSING feature (so the recipe's median imputation + the
# _isNA-style NA-preserving handling are exercised on reload).
make_saveload_fixture <- function(seed = 202L) {
  set.seed(seed)
  n_tr <- 80
  train_df <- data.frame(
    class    = rep(c("burned", "unburned"), each = n_tr / 2),
    block_id = rep(seq_len(16), length.out = n_tr),
    rbr_med  = c(rnorm(n_tr / 2, mean = 5, sd = 0.5),
                 rnorm(n_tr / 2, mean = 1, sd = 0.5)),
    elev_med = c(rnorm(n_tr / 2, 3), rnorm(n_tr / 2, -3)),
    slope_med = rnorm(n_tr, 10, 2),
    stringsAsFactors = FALSE
  )
  # Held-out rows: 12 rows; the FIRST has an NA rbr_med (-> imputed from the
  # train recipe median on both direct and reloaded paths).
  held_out <- data.frame(
    rbr_med  = c(NA, rnorm(11, mean = 3, sd = 2)),
    elev_med = rnorm(12),
    slope_med = rnorm(12, 10, 2),
    stringsAsFactors = FALSE
  )
  list(train_df = train_df, held_out = held_out,
       feature_cols = c("rbr_med", "elev_med", "slope_med"))
}

test_that("nested_refit model + recipe round-trip (save/reload) reproduces predictions", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fit_fn       <- nr_fit()
  transform_fn <- nr_transform()
  fx <- make_saveload_fixture()

  fit <- fit_fn(
    train_df = fx$train_df, feature_cols = fx$feature_cols,
    label_col = "class", group_col = "block_id", block_col = "block_id",
    val_frac = 0.2, params_fn = canon_params(),
    sampling_seed = 42L, fold_seed = 42L, nrounds_max = 60L,
    early_stopping_rounds = 15L
  )

  # ---- DIRECT core prediction on the held-out matrix ----
  X_direct <- transform_fn(
    df = fx$held_out, feature_cols = fit$feature_cols,
    med = fit$medians, ref_x_cols = fit$x_cols
  )
  d_direct <- xgboost::xgb.DMatrix(X_direct, missing = NA)
  p_direct <- predict(fit$model, d_direct)

  # ---- SAVE the deployable artifacts to tempfiles ----
  # The recipe captures exactly what the scoring path needs to reproduce the
  # design matrix: refit medians (NA = degenerate), x_cols order, feature_cols,
  # and best_iteration (provenance of the refit nrounds).
  recipe <- list(
    medians        = fit$medians,
    x_cols         = fit$x_cols,
    feature_cols   = fit$feature_cols,
    best_iteration = fit$best_iteration,
    impute_factor_missing = "MISSING"
  )
  model_path  <- tempfile(fileext = "_model.rds")
  recipe_path <- tempfile(fileext = "_recipe.rds")
  on.exit(unlink(c(model_path, recipe_path), force = TRUE), add = TRUE)
  # xgboost models must be saved via xgb.save.raw inside the list for full
  # fidelity; saveRDS of the booster object is the path the package uses
  # (train_final_model_direct saveRDS(model, rds_mod)), so mirror that here.
  saveRDS(fit$model, model_path)
  saveRDS(recipe, recipe_path)

  # ---- RELOAD + re-apply the recipe to the SAME held-out rows ----
  model_rl  <- readRDS(model_path)
  recipe_rl <- readRDS(recipe_path)

  # Recipe contents survived the round-trip.
  expect_equal(recipe_rl$x_cols, fit$x_cols)
  expect_equal(recipe_rl$feature_cols, fit$feature_cols)
  expect_equal(recipe_rl$best_iteration, fit$best_iteration)
  expect_equal(recipe_rl$medians, fit$medians)

  X_reload <- transform_fn(
    df = fx$held_out, feature_cols = recipe_rl$feature_cols,
    med = recipe_rl$medians, ref_x_cols = recipe_rl$x_cols,
    impute_factor_missing = recipe_rl$impute_factor_missing
  )
  d_reload <- xgboost::xgb.DMatrix(X_reload, missing = NA)
  p_reload <- predict(model_rl, d_reload)

  # ---- ASSERT: reloaded prediction == direct prediction on the same rows ----
  # The design matrices must be byte-identical (same column order, same imputed
  # values, same NA-preserving handling) ...
  expect_equal(colnames(X_reload), colnames(X_direct))
  expect_equal(as.matrix(X_reload), as.matrix(X_direct), tolerance = 0)
  # ... and therefore the predictions must be identical.
  expect_equal(p_reload, p_direct, tolerance = 1e-12)

  # The first held-out row's missing rbr_med was imputed with the TRAIN refit
  # median on BOTH paths (proves the recipe -- not the held-out data -- drove
  # imputation, and survived save/reload).
  expect_equal(unname(X_reload[1, "rbr_med"]), fit$medians$rbr_med,
               tolerance = 1e-12)
})
