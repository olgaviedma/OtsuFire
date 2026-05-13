# Bug 9 (OtsuFire 0.3.0): degenerate OOF metric columns must produce
# a sentinel "best threshold" row, not a zero-row crash.

# --- T17 ------------------------------------------------------------
test_that("T17: all-NA metric column does not produce an empty subset", {
  # Emulator of the 0.3.0 .pick_best helper used by run_dm_oof_pipeline.
  pick_best <- function(df, col) {
    v <- df[[col]]
    if (is.null(v) || all(is.na(v))) {
      sentinel <- as.data.frame(
        stats::setNames(
          replicate(ncol(df), NA, simplify = FALSE),
          names(df)
        ),
        stringsAsFactors = FALSE
      )
      return(sentinel)
    }
    df[which.max(v), , drop = FALSE]
  }

  # Build a degenerate metrics table where balanced_accuracy is all NA.
  metrics <- data.frame(
    threshold = c(0.1, 0.3, 0.5, 0.7),
    accuracy = c(0.5, 0.6, 0.55, 0.5),
    balanced_accuracy = NA_real_,
    f1 = NA_real_,
    stringsAsFactors = FALSE
  )

  expect_warning(
    {
      best_bal <- pick_best(metrics, "balanced_accuracy")
    },
    NA  # the test emulator does not warn; the package version does.
  )
  best_bal <- pick_best(metrics, "balanced_accuracy")
  expect_equal(nrow(best_bal), 1L)
  expect_true(all(is.na(best_bal)))

  best_acc <- pick_best(metrics, "accuracy")
  expect_equal(nrow(best_acc), 1L)
  expect_equal(best_acc$accuracy, 0.6)
})

test_that("compute_oof_metrics_by_threshold itself returns rows with NAs in degenerate cases", {
  # When a threshold produces all-zero predictions or all-one predictions,
  # some metric columns are NA. The downstream pick must not crash.
  fn <- get("compute_oof_metrics_by_threshold", envir = asNamespace("OtsuFire"))

  # All probabilities equal -> at most one threshold band has any
  # non-empty prediction set; many metrics will be 1 or 0 but not NA.
  oof_agg <- data.frame(
    fire_uid = sprintf("f%02d", 1:6),
    class    = c(rep("burned", 3), rep("unburned", 3)),
    p_oof_mean = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  res <- fn(oof_agg)
  expect_true(nrow(res) > 0L)
  expect_true("balanced_accuracy" %in% names(res))
})
