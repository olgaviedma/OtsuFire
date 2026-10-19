# ============================================================================
# OOF <-> FINAL parameter symmetry: a changed `params` / `model_params` must
# reach EVERY OOF fold (selection fit AND refit).
#
# Regression test for the defect in which run_oof_xgb() accepted `params` but
# never used it: each fold trained with .of_canonical_xgb_params, so a
# `model_params` override only changed the FINAL model. test-cfg-single-source.R
# mocks run_dm_oof_pipeline(), so it only proves the argument is *passed*, never
# that it is *used*.
#
# NO MOCK OF THE FIT HERE. The only binding replaced is a pass-through spy on
# xgboost::xgb.train() that records the `params` handed to xgboost and then
# delegates to the REAL function: every fold still trains a real booster.
# ============================================================================

ns_ovr <- asNamespace("OtsuFire")

make_oof_override_fixture <- function(n = 240L, seed = 11L) {
  set.seed(seed)
  feat <- c("rbr_aw_med", "elev_sd", "slope_sd")
  X <- matrix(stats::rnorm(n * 3L), n, 3L, dimnames = list(NULL, feat))
  src <- sample(rep(c("burned", "random", "otsu"), length.out = n))
  X[, 1] <- X[, 1] + ifelse(src == "burned", 1.5, 0)
  y <- as.integer(src == "burned")
  blk <- rep(seq_len(24L), length.out = n)
  L <- data.frame(
    fire_uid  = sprintf("u%04d", seq_len(n)),
    class     = ifelse(y == 1L, "burned", "unburned"),
    source    = c(burned = "internal_keep_qc",
                  random = "random_burnable_background",
                  otsu   = "otsu_patch_residual")[src],
    neg_type  = c(burned = NA_character_, random = "background_cell",
                  otsu   = "otsu_patch_drop")[src],
    block_id  = blk,
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = 24L)[blk],
    fold_rep2 = rep(c(3L, 1L, 2L), length.out = 24L)[blk],
    stringsAsFactors = FALSE
  )
  list(X = Matrix::Matrix(X, sparse = TRUE), y = y, labelled_df = L,
       prepared = as.data.frame(X), model_cols = feat)
}

# Run run_oof_xgb() for real, recording every params list that reaches xgboost.
run_oof_recording_params <- function(params) {
  fx <- make_oof_override_fixture()
  real_train <- xgboost::xgb.train
  seen <- list()
  testthat::local_mocked_bindings(
    xgb.train = function(params, ...) {
      seen[[length(seen) + 1L]] <<- params
      real_train(params = params, ...)
    },
    .package = "xgboost"
  )
  res <- suppressWarnings(suppressMessages(get("run_oof_xgb", envir = ns_ovr)(
    XL_mat = fx$X, y = fx$y, labelled_df = fx$labelled_df,
    fold_cols = c("fold_rep1", "fold_rep2"), params = params,
    nrounds_max = 6L, early_stop = 3L, seed_base = 42L,
    group_col = "block_id", val_frac = 0.15,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    prepared_labelled = fx$prepared, model_cols = fx$model_cols,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    verbose = 0
  )))
  list(res = res, seen = seen)
}

test_that("a non-canonical params block reaches every OOF fold (selection + refit)", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat")))

  blk <- get(".of_canonical_model_params", envir = ns_ovr)()
  blk$max_depth <- 1
  blk$eta       <- 0.3
  out <- run_oof_recording_params(c(blk, list(scale_pos_weight = 2)))

  # 2 reps x 3 folds x (selection fit + refit) = 12 real xgboost trainings.
  expect_length(out$seen, 12L)
  for (i in seq_along(out$seen)) {
    expect_equal(out$seen[[i]]$max_depth, 1, info = paste("fit", i))
    expect_equal(out$seen[[i]]$eta, 0.3, info = paste("fit", i))
    # every other canonical field is carried unchanged
    for (nm in setdiff(names(blk), c("max_depth", "eta"))) {
      expect_identical(out$seen[[i]][[nm]], blk[[nm]],
                       info = paste("fit", i, "field", nm))
    }
  }
})

test_that("the canonical block still yields the canonical params in every fold", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat")))

  blk <- get(".of_canonical_model_params", envir = ns_ovr)()
  out <- run_oof_recording_params(c(blk, list(scale_pos_weight = 2)))
  expect_length(out$seen, 12L)
  for (p in out$seen) {
    expect_equal(p$max_depth, 5)
    expect_equal(p$eta, 0.05)
    expect_equal(p$min_child_weight, 5)
  }
})

test_that("scale_pos_weight is recomputed per fold, never taken from `params`", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat")))

  blk <- get(".of_canonical_model_params", envir = ns_ovr)()
  out <- run_oof_recording_params(c(blk, list(scale_pos_weight = 999)))
  spw <- vapply(out$seen, function(p) p$scale_pos_weight, numeric(1))
  expect_false(any(spw == 999))
  expect_true(all(is.finite(spw)) && all(spw > 0))
})

test_that("run_oof_xgb() validates `params`", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- make_oof_override_fixture()
  call_oof <- function(params) {
    get("run_oof_xgb", envir = ns_ovr)(
      XL_mat = fx$X, y = fx$y, labelled_df = fx$labelled_df,
      params = params, nrounds_max = 3L, early_stop = 2L, seed_base = 1L,
      group_col = "block_id", val_frac = 0.15, impute_numeric = "median",
      impute_factor_missing = "MISSING", prepared_labelled = fx$prepared,
      model_cols = fx$model_cols, random_to_burned_ratio = 1,
      otsu_unburned_to_burned_ratio = 1, verbose = 0)
  }
  expect_error(call_oof(NULL), "params")
  expect_error(call_oof("eta=0.3"), "params")
})
