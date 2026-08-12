# ============================================================================
# Gate 1E (2026-06-09): RUNTIME feature-schema PARITY GUARD.
#
# Proves the STRUCTURAL feature-space contract guard: identical structure PASSES;
# a reordered / missing / extra column or a dropped `_isNA` indicator ABORTS; and
# the fingerprint is INVARIANT to the legitimately-varying fitted statistics
# (imputation medians, best_iteration) and to wall-clock time. Also proves the
# saved->loaded scoring round-trip and the canonical-OOF / Phase-B cross-checks.
#
# These assert the CONTRACT (names, order, counts, encoding, weights policy,
# contract version), NOT the equality of fitted recipes (medians / spw / seeds).
# ============================================================================

ns <- asNamespace("OtsuFire")
g  <- function(nm) get(nm, envir = ns)

fsf      <- g("feature_schema_fingerprint")
assert_eq <- g(".of_assert_schema_fingerprints_equal")
oof_canon <- g(".of_oof_canonical_fingerprint")
build_mani <- g(".of_build_schema_parity_manifest")

# A canonical structural recipe shape (base + `_isNA` block) used as the
# reference schema across the structural tests.
ref_cols <- function() {
  base <- c("rbr_med", "elev_med", "slope_med", "hs_conf_mean")
  isna <- paste0(base, "_isNA")
  list(
    base_features              = base,
    missing_indicator_features = isna,
    final_feature_order        = c(base, isna),
    feature_cols               = c(base, isna),
    x_cols                     = c(base, isna)
  )
}

# ---------------------------------------------------------------------------
# (1) same schema -> PASS (identical fingerprint).
# ---------------------------------------------------------------------------
test_that("same structural schema yields the SAME fingerprint (guard PASSES)", {
  a <- fsf(ref_cols())
  b <- fsf(ref_cols())
  expect_identical(a$hash, b$hash)
  expect_silent(assert_eq("same", a, b))
  expect_true(assert_eq("same", a, b))
})

# ---------------------------------------------------------------------------
# (2) changed column ORDER -> FAIL (guard aborts).
# ---------------------------------------------------------------------------
test_that("changed column ORDER changes the fingerprint and ABORTS", {
  ref <- ref_cols()
  reordered <- ref
  # swap the first two base features (same SET, different ORDER).
  reordered$base_features <- ref$base_features[c(2, 1, 3, 4)]
  reordered$final_feature_order <-
    c(reordered$base_features, ref$missing_indicator_features)
  reordered$feature_cols <- reordered$final_feature_order
  reordered$x_cols       <- reordered$final_feature_order

  a <- fsf(ref)
  b <- fsf(reordered)
  expect_false(identical(a$hash, b$hash))
  expect_error(assert_eq("order", a, b),
               "Feature-schema parity guard ABORT")
})

# ---------------------------------------------------------------------------
# (3) a feature ABSENT -> FAIL.
# ---------------------------------------------------------------------------
test_that("a base feature ABSENT changes the fingerprint and ABORTS", {
  ref <- ref_cols()
  miss <- ref
  drop <- "elev_med"
  miss$base_features       <- setdiff(ref$base_features, drop)
  miss$final_feature_order <- setdiff(ref$final_feature_order, drop)
  miss$feature_cols        <- setdiff(ref$feature_cols, drop)
  miss$x_cols              <- setdiff(ref$x_cols, drop)

  expect_error(assert_eq("absent feature", fsf(ref), fsf(miss)),
               "STRUCTURAL feature-space")
})

# ---------------------------------------------------------------------------
# (4) an `_isNA` indicator ABSENT -> FAIL.
# ---------------------------------------------------------------------------
test_that("an `_isNA` indicator ABSENT changes the fingerprint and ABORTS", {
  ref <- ref_cols()
  miss <- ref
  drop <- "elev_med_isNA"
  miss$missing_indicator_features <- setdiff(ref$missing_indicator_features, drop)
  miss$final_feature_order        <- setdiff(ref$final_feature_order, drop)
  miss$feature_cols               <- setdiff(ref$feature_cols, drop)
  miss$x_cols                     <- setdiff(ref$x_cols, drop)

  a <- fsf(ref); b <- fsf(miss)
  expect_false(identical(a$hash, b$hash))
  # n_indicators must reflect the drop.
  expect_identical(a$payload$n_indicators, length(ref$missing_indicator_features))
  expect_identical(b$payload$n_indicators,
                   length(ref$missing_indicator_features) - 1L)
  expect_error(assert_eq("absent isNA", a, b),
               "Feature-schema parity guard ABORT")
})

# ---------------------------------------------------------------------------
# (5) an EXTRA column -> FAIL.
# ---------------------------------------------------------------------------
test_that("an EXTRA column changes the fingerprint and ABORTS", {
  ref <- ref_cols()
  extra <- ref
  extra$base_features       <- c(ref$base_features, "doy_med")
  extra$final_feature_order <- c(ref$base_features, "doy_med",
                                 ref$missing_indicator_features)
  extra$feature_cols        <- extra$final_feature_order
  extra$x_cols              <- extra$final_feature_order

  a <- fsf(ref); b <- fsf(extra)
  expect_false(identical(a$hash, b$hash))
  expect_identical(b$payload$n_total, a$payload$n_total + 1L)
  expect_error(assert_eq("extra", a, b),
               "Feature-schema parity guard ABORT")
})

# ---------------------------------------------------------------------------
# (6) same schema but DIFFERENT imputation medians -> PASS (medians not in the
#     structural contract).
# ---------------------------------------------------------------------------
test_that("DIFFERENT imputation medians do NOT change the fingerprint (PASS)", {
  r1 <- c(ref_cols(), list(
    numeric_medians = list(rbr_med = 0.1, elev_med = 100, slope_med = 5),
    impute_numeric = "median", impute_factor_missing = "MISSING"))
  r2 <- c(ref_cols(), list(
    # wildly different medians, same STRUCTURE.
    numeric_medians = list(rbr_med = 999, elev_med = -42, slope_med = 17),
    impute_numeric = "median", impute_factor_missing = "MISSING"))
  expect_identical(fsf(r1)$hash, fsf(r2)$hash)
})

# ---------------------------------------------------------------------------
# (7) same schema but DIFFERENT best_iteration / spw / seed -> PASS (fitted
#     statistics are not part of the structural contract).
# ---------------------------------------------------------------------------
test_that("DIFFERENT best_iteration / spw / seed do NOT change the fingerprint (PASS)", {
  r1 <- c(ref_cols(), list(best_iteration = 12L, spw_refit = 1.3, seed = 42L))
  r2 <- c(ref_cols(), list(best_iteration = 9999L, spw_refit = 7.7, seed = 7L))
  expect_identical(fsf(r1)$hash, fsf(r2)$hash)
})

# ---------------------------------------------------------------------------
# (8) scoring after save/load -> PASS (fingerprint round-trips through RDS).
#     Persist a recipe carrying the saved fingerprint, reload it, recompute the
#     scoring-leg fingerprint from the SAVED structural cols, and assert the
#     guard's round-trip assertion passes.
# ---------------------------------------------------------------------------
test_that("saved->loaded recipe fingerprint round-trips (scoring guard PASSES)", {
  ref <- ref_cols()
  saved_fp <- fsf(ref)
  recipe <- list(
    cols = list(base_features = ref$base_features,
                missing_indicator_features = ref$missing_indicator_features,
                final_feature_order = ref$final_feature_order,
                feature_cols = ref$feature_cols, x_cols = ref$x_cols),
    schema_fingerprint = list(hash = saved_fp$hash, payload = saved_fp$payload))
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds, force = TRUE), add = TRUE)
  saveRDS(recipe, rds)
  loaded <- readRDS(rds)

  # Recompute the scoring-leg fingerprint exactly as score_with_final_model does:
  # from the recipe's structural partition + the produced (here: identical)
  # x_cols.
  produced_fp <- fsf(list(
    base_features              = loaded$cols$base_features,
    missing_indicator_features = loaded$cols$missing_indicator_features,
    final_feature_order        = loaded$cols$final_feature_order,
    feature_cols               = loaded$cols$feature_cols,
    x_cols                     = loaded$cols$x_cols))
  expect_identical(loaded$schema_fingerprint$hash, produced_fp$hash)
  expect_true(assert_eq("round-trip", loaded$schema_fingerprint, produced_fp))
})

# ---------------------------------------------------------------------------
# (9) fingerprint INDEPENDENT of timestamps (two computations at different
#     wall-clock times are identical).
# ---------------------------------------------------------------------------
test_that("fingerprint is INDEPENDENT of wall-clock time", {
  a <- fsf(ref_cols())
  Sys.sleep(1.1)  # cross a wall-clock second boundary.
  b <- fsf(ref_cols())
  expect_identical(a$hash, b$hash)
  expect_identical(a$text, b$text)
  # And no timestamp leaks into the hashed body.
  expect_false(grepl(format(Sys.Date()), a$text, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# Cross-leg helpers: canonical OOF (all folds must agree) + Phase B abort.
# ---------------------------------------------------------------------------
test_that("canonical OOF requires ALL folds to share the contract; a divergent fold ABORTS", {
  ref <- ref_cols()
  same <- list(rep1_fold1 = fsf(ref), rep1_fold2 = fsf(ref),
               rep2_fold1 = fsf(ref))
  canon <- oof_canon(same)
  expect_identical(canon$hash, fsf(ref)$hash)

  # One fold drops an `_isNA` -> canonical establishment ABORTS.
  bad <- ref
  bad$missing_indicator_features <- ref$missing_indicator_features[-1]
  bad$final_feature_order <- setdiff(ref$final_feature_order,
                                     ref$missing_indicator_features[1])
  bad$feature_cols <- bad$final_feature_order
  bad$x_cols <- bad$final_feature_order
  diverged <- list(rep1_fold1 = fsf(ref), rep1_fold2 = fsf(bad))
  expect_error(oof_canon(diverged), "Feature-schema parity guard ABORT")
})

test_that("Phase B manifest ABORTS when OOF/FINAL/scoring fingerprints are incompatible", {
  ref <- ref_cols()
  fp_ok  <- fsf(ref)
  bad <- ref
  bad$base_features <- c(ref$base_features, "doy_med")
  bad$final_feature_order <- c(bad$base_features, ref$missing_indicator_features)
  bad$feature_cols <- bad$final_feature_order
  bad$x_cols <- bad$final_feature_order
  fp_bad <- fsf(bad)

  oof_guard <- list(canonical = fp_ok,
                    per_fold = list(rep1_fold1 = fp_ok, rep1_fold2 = fp_ok))

  # Compatible legs + phase_b -> pass manifest, no abort.
  mani <- build_mani(oof_guard = oof_guard, final_fp = fp_ok,
                     scoring_fp = fp_ok, phase_b = TRUE)
  expect_identical(mani$guard_result, "pass")
  expect_identical(mani$contract_version,
                   g(".OF_SUPERVISED_FEATURE_CONTRACT_VERSION"))
  expect_identical(mani$oof_canonical_fingerprint, fp_ok$hash)
  expect_length(mani$oof_per_fold_fingerprints, 2L)

  # FINAL diverges + phase_b -> ABORT the run.
  expect_error(
    build_mani(oof_guard = oof_guard, final_fp = fp_bad,
               scoring_fp = fp_ok, phase_b = TRUE),
    "Feature-schema parity guard ABORT \\(Phase B run\\)")

  # Same divergence but NOT phase_b -> manifest records "abort" without stopping.
  mani2 <- build_mani(oof_guard = oof_guard, final_fp = fp_bad,
                      scoring_fp = fp_ok, phase_b = FALSE)
  expect_identical(mani2$guard_result, "abort")
})

# ---------------------------------------------------------------------------
# END-TO-END parity on REAL data: drive the OOF / FINAL / scoring column
# derivations through the production helpers from ONE upstream frame and assert
# the structural fingerprints are IDENTICAL across the three legs (the guard's
# real-data PASS path). Reuses the shared sc_* contract builders.
# ---------------------------------------------------------------------------
test_that("OOF, FINAL and scoring produce the SAME structural fingerprint on one upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("dplyr")

  df <- sc_isna_upstream_df()
  oof_cols   <- sc_oof_model_cols_from(df)
  fit        <- sc_final_fit_from(df)
  final_cols <- fit$x_cols
  score_cols <- sc_scoring_x_cols_from(df, fit)

  fp_oof   <- fsf(list(final_feature_order = oof_cols, feature_cols = oof_cols,
                       x_cols = oof_cols))
  fp_final <- fsf(list(base_features = fit$base_features,
                       missing_indicator_features = fit$missing_indicator_features,
                       final_feature_order = fit$final_feature_order,
                       feature_cols = fit$feature_cols, x_cols = final_cols))
  fp_score <- fsf(list(base_features = fit$base_features,
                       missing_indicator_features = fit$missing_indicator_features,
                       final_feature_order = fit$final_feature_order,
                       feature_cols = fit$feature_cols, x_cols = score_cols))

  expect_identical(fp_oof$hash, fp_final$hash)
  expect_identical(fp_final$hash, fp_score$hash)
  # The structural counts are coherent and non-trivial.
  expect_gt(fp_final$payload$n_indicators, 0L)
  expect_identical(fp_final$payload$n_total,
                   fp_final$payload$n_base + fp_final$payload$n_indicators)
})
