# Tests for the dynamic fold-repeat column detection used by the orchestrator
# when handing fold_cols to run_oof_diagnostics().

test_that(".of_fold_rep_cols equals the old hard-coded value for 2 reps", {
  df <- data.frame(fold_rep1 = 1:3, fold_rep2 = 1:3, x = 1:3)
  expect_identical(OtsuFire:::.of_fold_rep_cols(df), c("fold_rep1", "fold_rep2"))
  # also works on a plain names() vector
  expect_identical(OtsuFire:::.of_fold_rep_cols(names(df)), c("fold_rep1", "fold_rep2"))
})

test_that(".of_fold_rep_cols detects a 3rd repeat and orders by repeat NUMBER", {
  df <- data.frame(fold_rep2 = 1:3, fold_rep1 = 1:3, fold_rep3 = 1:3)
  expect_identical(OtsuFire:::.of_fold_rep_cols(df),
                   c("fold_rep1", "fold_rep2", "fold_rep3"))
  # numeric, NOT lexicographic: fold_rep2 must come before fold_rep10
  nm <- c("fold_rep10", "fold_rep2", "fold_rep1")
  expect_identical(OtsuFire:::.of_fold_rep_cols(nm),
                   c("fold_rep1", "fold_rep2", "fold_rep10"))
})

test_that(".of_fold_rep_cols ignores non-matching names and returns empty when none", {
  df <- data.frame(fold_repA = 1, fold_rep = 1, foldrep1 = 1, x = 1)
  expect_identical(OtsuFire:::.of_fold_rep_cols(df), character(0))
})

test_that("the value passed to run_oof_diagnostics is the dynamic detection", {
  # The orchestrator computes fold_cols <- .of_fold_rep_cols(train_features) and
  # forwards it verbatim. Assert that, for the canonical 2-rep training frame, the
  # detection equals exactly what the OOF stage used to receive hard-coded.
  train_features <- data.frame(
    fire_uid = c("f1", "f2"), class = c("burned", "unburned"),
    fold_rep1 = c(1L, 2L), fold_rep2 = c(2L, 1L), rbr_med = c(100, 50))
  expect_identical(OtsuFire:::.of_fold_rep_cols(train_features),
                   c("fold_rep1", "fold_rep2"))
})
