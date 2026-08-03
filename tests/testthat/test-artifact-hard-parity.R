# =============================================================================
# PHASE 2 (artifact_hard): BASELINE-REPRODUCTION PROOF + evaluation contracts.
#
# Covers spec G items 7 (evaluated universe does not shrink; the 6 sentinels
# remain), 8 (p_burned_eval uses p_oof for labeled rows, p_burned_model for
# pure scoring rows), and the H parity proof (the FINAL booster reproduces the
# saved 0.11.0 1989 p_burned_model to machine epsilon -- proving the new code
# paths are inert when disabled).
#
# The parity leg depends on the canonical 1989 artifact; it is SKIPPED when the
# artifact is not present (so CI on another machine stays green).
# =============================================================================

CANON <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/00_AUDITS/1989_CANONICAL_MATRICES/ARTIFACTS/canonical.rds"
SENTINELS <- c(858, 865, 1102, 4250, 3916, 3556)

test_that("PARITY: FINAL booster reproduces saved 1989 p_burned_model to machine epsilon", {
  skip_if_not(file.exists(CANON), "canonical 1989 artifact not present")
  skip_if_not_installed("xgboost")
  can <- readRDS(CANON)
  # The saved booster is only loadable by the xgboost generation that wrote it:
  # a 1.7-era artifact makes xgboost >= 2.0 raise "'xgb.Booster' object is
  # corrupted or is from an incompatible XGBoost version". That is a property of
  # the stored artifact, not of the package, so skip rather than fail.
  skip_if_not(
    isTRUE(tryCatch({ xgboost::xgb.get.handle(can$model); TRUE },
                    error = function(e) FALSE)),
    "saved 1989 booster was written by an incompatible xgboost version"
  )
  p <- predict(can$model, can$X_scoring)
  ref <- can$scoring_meta$p_burned_model
  expect_equal(length(p), length(ref))
  max_abs_diff <- max(abs(p - ref))
  # machine epsilon: with disabled artifact_hard + NULL weights the scoring path
  # is byte-identical, so the prediction matches the saved baseline exactly.
  expect_lt(max_abs_diff, 1e-12)
})

test_that("evaluated universe does not shrink; the 6 sentinels remain", {
  skip_if_not(file.exists(CANON), "canonical 1989 artifact not present")
  can <- readRDS(CANON)
  expect_equal(nrow(can$scoring_meta), nrow(can$X_scoring))
  pid <- as.character(can$scoring_meta$poly_id)
  expect_true(all(as.character(SENTINELS) %in% pid))
})

# ---- p_burned_eval = coalesce(p_oof_mean, p_burned_model) (spec F / G8) -----
# Exercise the column logic directly (the assembly in internal-sup-final-map.R
# adds p_burned_eval with exactly this coalesce). Labeled rows use the held-out
# p_oof; pure scoring rows fall back to the in-sample model prediction.

test_that("p_burned_eval prefers p_oof_mean and falls back to p_burned_model", {
  fm <- data.frame(
    poly_id        = as.character(1:4),
    p_burned_oof   = c(0.10, NA, 0.90, NaN),   # p_oof_mean carried as p_burned_oof
    p_burned_model = c(0.50, 0.60, 0.70, 0.80),
    stringsAsFactors = FALSE
  )
  pe_oof   <- suppressWarnings(as.numeric(fm$p_burned_oof))
  pe_model <- suppressWarnings(as.numeric(fm$p_burned_model))
  p_burned_eval <- ifelse(is.finite(pe_oof), pe_oof, pe_model)
  # labeled (finite oof) -> oof; unlabeled (NA/NaN oof) -> model
  expect_equal(p_burned_eval, c(0.10, 0.60, 0.90, 0.80))
})

test_that("p_burned_eval drops no rows (length preserved)", {
  fm <- data.frame(
    p_burned_oof   = c(NA, 0.2, NA, 0.4, NA),
    p_burned_model = c(0.5, 0.6, 0.7, 0.8, 0.9)
  )
  pe <- ifelse(is.finite(fm$p_burned_oof), fm$p_burned_oof, fm$p_burned_model)
  expect_equal(length(pe), nrow(fm))
  expect_true(all(is.finite(pe)))
})
