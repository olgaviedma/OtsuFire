# =============================================================================
# PHASE 2 (artifact_hard): per-row sample weights (spec E).
#
# Covers spec G item 3: total artifact_hard weight = 0.5*positives total;
# per-source logging; scale_pos_weight recomputed from WEIGHTED sums (no double
# weighting); and the OFF path (NULL weights) is byte-identical to today.
# =============================================================================

ns <- asNamespace("OtsuFire")
resolve_w <- ns$.of_resolve_sample_weights
nested_fit <- ns$.of_nested_refit_fit

test_that("OFF path: no source_weights + no artifact_hard rows -> NULL weights", {
  r <- resolve_w(class = c("burned", "unburned", "burned"),
                 source = c(NA, "random_burnable_background", NA))
  expect_null(r$weights)
  expect_equal(nrow(r$log), 0L)
})

test_that("artifact_hard total weight = total_weight_ratio * positives total", {
  cls <- c(rep("burned", 10), rep("unburned", 4), rep("unburned", 5))
  src <- c(rep(NA, 10), rep("random_burnable_background", 4),
           rep("artifact_hard", 5))
  is_pos <- cls == "burned"
  is_ah  <- src %in% "artifact_hard"
  # explicit ratio drives the pool-level balance exactly.
  for (twr in c(0.10, 0.25, 0.50)) {
    r <- resolve_w(class = cls, source = src, total_weight_ratio = twr)
    expect_equal(sum(r$weights[is_ah]), twr * sum(r$weights[is_pos]))
    # positives keep weight 1 (not double-weighted)
    expect_true(all(r$weights[is_pos] == 1))
    # random negatives keep weight 1
    expect_true(all(r$weights[src %in% "random_burnable_background"] == 1))
  }
})

test_that("total_weight_ratio defaults to 0.10 (the conservative config default)", {
  cls <- c(rep("burned", 10), rep("unburned", 5))
  src <- c(rep(NA, 10), rep("artifact_hard", 5))
  r <- resolve_w(class = cls, source = src)            # no ratio supplied
  is_pos <- cls == "burned"; is_ah <- src %in% "artifact_hard"
  expect_equal(sum(r$weights[is_ah]), 0.10 * sum(r$weights[is_pos]))
})

test_that("per-source weight log records n / per-row / total", {
  cls <- c(rep("burned", 6), rep("unburned", 3))
  src <- c(rep(NA, 6), rep("artifact_hard", 3))
  r <- resolve_w(class = cls, source = src, total_weight_ratio = 0.5)
  lg <- r$log
  ah_row <- lg[lg$source == "artifact_hard", ]
  expect_equal(ah_row$n, 3L)
  expect_equal(ah_row$total_weight, 0.5 * 6)
  expect_equal(ah_row$per_row_weight, (0.5 * 6) / 3)
})

test_that("source_weights override pins a source weight and disables auto-balance", {
  cls <- c(rep("burned", 4), rep("unburned", 2))
  src <- c(rep(NA, 4), rep("artifact_hard", 2))
  r <- resolve_w(class = cls, source = src,
                 source_weights = c(artifact_hard = 0.3))
  is_ah <- src %in% "artifact_hard"
  expect_true(all(r$weights[is_ah] == 0.3))   # pinned, not the 0.5-balance value
})

# ---- nested_refit weighting interaction -----------------------------------
# Build a tiny labelled frame with a couple of model features so the shared core
# can actually fit. We assert: (a) NULL sample_weights => identical model bytes
# to the no-arg call (OFF path byte-identical); (b) weighted spw == weighted sum
# ratio (no double weighting).

mk_train <- function(n = 80, seed = 7) {
  set.seed(seed)
  data.frame(
    class = c(rep("burned", n / 2), rep("unburned", n / 2)),
    block_id = sample(letters[1:8], n, TRUE),
    rbr_med = c(rnorm(n / 2, 0.6, 0.1), rnorm(n / 2, 0.2, 0.1)),
    elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE
  )
}
feat <- c("rbr_med", "elev_med")

test_that("sample_weights = NULL reproduces the unweighted model bytes", {
  td <- mk_train()
  f0 <- nested_fit(td, feature_cols = feat, label_col = "class",
                   block_col = "block_id", val_frac = 0.2,
                   sampling_seed = 1L, fold_seed = 1L,
                   nrounds_max = 30L, early_stopping_rounds = 10L)
  f1 <- nested_fit(td, feature_cols = feat, label_col = "class",
                   block_col = "block_id", val_frac = 0.2,
                   sampling_seed = 1L, fold_seed = 1L,
                   sample_weights = NULL,
                   nrounds_max = 30L, early_stopping_rounds = 10L)
  expect_identical(xgboost::xgb.save.raw(f0$model),
                   xgboost::xgb.save.raw(f1$model))
  # unweighted refit spw == count ratio
  expect_equal(f0$spw_refit, sum(td$class == "unburned") / sum(td$class == "burned"))
})

test_that("weighted refit spw uses WEIGHTED sums (no double weighting)", {
  td <- mk_train()
  y <- as.integer(td$class == "burned")
  # give the negatives weight 0.5, positives weight 1
  w <- ifelse(y == 1L, 1.0, 0.5)
  f <- nested_fit(td, feature_cols = feat, label_col = "class",
                  block_col = "block_id", val_frac = 0.2,
                  sampling_seed = 1L, fold_seed = 1L,
                  sample_weights = w,
                  nrounds_max = 30L, early_stopping_rounds = 10L)
  expect_equal(f$spw_refit, sum(w[y == 0L]) / sum(w[y == 1L]))
})

test_that("sample_weights length mismatch errors", {
  td <- mk_train()
  expect_error(
    nested_fit(td, feature_cols = feat, label_col = "class",
               block_col = "block_id", val_frac = 0.2,
               sampling_seed = 1L, fold_seed = 1L,
               sample_weights = rep(1, nrow(td) - 1L),
               nrounds_max = 10L, early_stopping_rounds = 5L),
    "align to train_df"
  )
})
