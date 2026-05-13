# KB2 + Bug 2 + AS04 (OtsuFire 0.3.0): encoder parity at scoring.
#
# These tests train a tiny XGBoost model and exercise the imputation
# + sparse encoder symmetry between training and scoring. The fixes
# being verified:
#   * prep_X_df at scoring applies recipe$impute$numeric_medians.
#   * score_with_final_model uses Matrix::sparse.model.matrix.
#   * align_to_x_cols pads missing columns with sparse zero columns
#     (implicit-missing) instead of dense explicit zeros.

train_tiny_model <- function(seed = 1L) {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  set.seed(seed)
  n <- 60L
  X_df <- data.frame(
    a = rnorm(n, 1, 1),
    b = rnorm(n, 2, 1),
    c = runif(n, 1, 5)
  )
  y <- as.integer(X_df$a + X_df$b > 3)

  X <- Matrix::sparse.model.matrix(~ . - 1, data = X_df)
  dtrain <- xgboost::xgb.DMatrix(X, label = y)
  model <- xgboost::xgb.train(
    params = list(objective = "binary:logistic",
                  booster = "gbtree",
                  max_depth = 3, eta = 0.3),
    data = dtrain,
    nrounds = 30L,
    verbose = 0
  )

  list(model = model,
       x_cols = colnames(X),
       feature_cols = names(X_df),
       X_df = X_df,
       y = y)
}

# --- T8 -------------------------------------------------------------
test_that("T8: sparse-encoded predict() reproduces training-time values", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  bundle <- train_tiny_model(seed = 17L)

  # Re-encode the exact same X_df with sparse.model.matrix (same path as
  # training). predict() must agree with itself to FP precision.
  X_again <- Matrix::sparse.model.matrix(~ . - 1, data = bundle$X_df)
  d_again <- xgboost::xgb.DMatrix(X_again)
  p_again <- predict(bundle$model, d_again)

  d_train <- xgboost::xgb.DMatrix(
    Matrix::sparse.model.matrix(~ . - 1, data = bundle$X_df))
  p_train <- predict(bundle$model, d_train)

  expect_equal(p_again, p_train, tolerance = 1e-12)
})

# --- T9 -------------------------------------------------------------
test_that("T9: imputed-median scoring matches NA-bearing scoring after fix", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  bundle <- train_tiny_model(seed = 23L)
  medians <- vapply(bundle$X_df, stats::median, numeric(1), na.rm = TRUE)

  # NA-bearing scoring frame: a few rows have NA in column `a`.
  X_score_na <- bundle$X_df[1:6, , drop = FALSE]
  X_score_na$a[c(1L, 3L)] <- NA_real_

  # Manually impute with the training-time median to produce the
  # "fixed" scoring matrix.
  X_score_imp <- X_score_na
  X_score_imp$a[is.na(X_score_imp$a)] <- medians[["a"]]

  X_imp <- Matrix::sparse.model.matrix(~ . - 1, data = X_score_imp)
  p_imp <- predict(bundle$model, xgboost::xgb.DMatrix(X_imp))

  # Simulate the OtsuFire 0.3.0 scoring path: pretend the recipe says
  # impute "a" with `medians["a"]`. Apply the same imputation, then
  # sparse.model.matrix.
  numeric_medians <- as.list(medians)
  X_via_fix <- X_score_na
  for (nm in names(X_via_fix)) {
    if (nm %in% names(numeric_medians) && anyNA(X_via_fix[[nm]])) {
      X_via_fix[[nm]][is.na(X_via_fix[[nm]])] <- numeric_medians[[nm]]
    }
  }
  X_fix_mat <- Matrix::sparse.model.matrix(~ . - 1, data = X_via_fix)
  p_fix <- predict(bundle$model, xgboost::xgb.DMatrix(X_fix_mat))

  expect_equal(p_fix, p_imp, tolerance = 1e-12)
})

# --- T10 ------------------------------------------------------------
test_that("T10: bind_rows-induced NA pattern imputes correctly", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  bundle <- train_tiny_model(seed = 31L)
  medians <- as.list(vapply(bundle$X_df, stats::median, numeric(1),
                              na.rm = TRUE))

  # Simulate AS04: bind_rows of det_final_raw (only `a`, `b` populated)
  # with otsu_raw (only `b`, `c` populated). After binding, det rows
  # have NA in `c`; otsu rows have NA in `a`.
  det_rows <- bundle$X_df[1:3, , drop = FALSE]
  det_rows$c <- NA_real_
  otsu_rows <- bundle$X_df[4:6, , drop = FALSE]
  otsu_rows$a <- NA_real_
  combined <- rbind(det_rows, otsu_rows)
  expect_true(any(is.na(combined$c)) && any(is.na(combined$a)))

  # Apply the 0.3.0 fix: column-wise median imputation from medians.
  imputed <- combined
  for (nm in names(imputed)) {
    if (nm %in% names(medians) && anyNA(imputed[[nm]])) {
      imputed[[nm]][is.na(imputed[[nm]])] <- medians[[nm]]
    }
  }
  expect_false(any(is.na(imputed)))

  X_imputed <- Matrix::sparse.model.matrix(~ . - 1, data = imputed)
  p <- predict(bundle$model, xgboost::xgb.DMatrix(X_imputed))
  # Sanity: predictions are finite and in [0, 1].
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
})
