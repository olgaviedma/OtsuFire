# ============================================================================
# Gate 1D.8 (2026-06-09): OOF <-> FINAL <-> SCORING `_isNA` COMPANION PARITY,
# END-TO-END FROM ONE RAW FEATURE FRAME (NO pre-seeded `_isNA` columns).
#
# This suite was the executable record of the KEYSTONE defect (Gate 1D.7): fed
# the SAME extract_supervised_features() output, OOF synthesised 51 `_isNA`
# companions (build_design_matrix_patches) while FINAL synthesised NONE (the
# nested core consumed feature_cols verbatim and extract_supervised_features()
# writes no `_isNA` to the GPKG). The deployed model + the OOF diagnostics ran on
# DIFFERENT feature spaces (51 vs 102 x_cols).
#
# Gate 1D.8 unified OOF, FINAL and scoring on ONE shared recipe / builder: the
# `_isNA` synthesis now lives in the shared core (.of_synthesize_isna_companions,
# invoked by .of_nested_coerce_features) and is the SAME for every path. This
# suite is now a NORMAL PASSING parity test: from ONE raw frame WITHOUT any
# pre-seeded `_isNA` columns it asserts identical column SET, ORDER, COUNT and
# hash of the deployed feature order across OOF, FINAL and SCORING, and that the
# legacy and nested_refit protocols yield the SAME feature space.
# ============================================================================

# Gate 1D.8 (2026-06-09): the upstream-frame builder + the per-stage column
# resolvers (OOF / FINAL / scoring) and the `_isNA` / base partition utilities
# now live in helper-supervised-contracts.R as `sc_*` builders, shared with the
# 18-point contract suite (test-oof-final-recipe-parity-18.R) so there is ONE
# definition. This file keeps short local aliases for readability of the
# historical assertions below.
ns_ip <- asNamespace("OtsuFire")
g_ip  <- function(nm) get(nm, envir = ns_ip)

make_isna_upstream_df <- sc_isna_upstream_df
oof_model_cols_from   <- sc_oof_model_cols_from
final_fit_from        <- sc_final_fit_from
final_x_cols_from     <- sc_final_x_cols_from
scoring_x_cols_from   <- sc_scoring_x_cols_from
isna_of               <- sc_isna_of
base_of               <- sc_base_of

# ---------------------------------------------------------------------------
# (1) BASE-FEATURE parity (load-bearing 1D.4 guarantee): stripped of `_isNA`,
#     OOF and FINAL train on the SAME canonical whitelist features in the SAME
#     canonical order, end-to-end from ONE upstream frame.
# ---------------------------------------------------------------------------
test_that("base feature set (whitelist, _isNA stripped) is identical OOF vs FINAL from one upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  wl <- g_ip(".supervised_feature_cols")

  df <- make_isna_upstream_df()
  oof_cols   <- oof_model_cols_from(df)
  final_cols <- final_x_cols_from(df)

  oof_base   <- base_of(oof_cols)
  final_base <- base_of(final_cols)

  expect_setequal(oof_base, final_base)
  expect_identical(oof_base,   intersect(wl, oof_base))
  expect_identical(final_base, intersect(wl, final_base))
  expect_identical(oof_base, final_base)
})

# ---------------------------------------------------------------------------
# (2) KEYSTONE PARITY (the formerly-xfail boundary, NOW PASSING): the `_isNA`
#     companion SET + the FULL deployed x_cols are IDENTICAL across OOF, FINAL
#     and SCORING when all three are fed the SAME raw frame with NO pre-seeded
#     `_isNA` columns. hash(OOF) == hash(FINAL) == hash(SCORING).
# ---------------------------------------------------------------------------
test_that("_isNA companion set + x_cols are IDENTICAL across OOF, FINAL and SCORING (one upstream frame)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  df <- make_isna_upstream_df()
  fit        <- final_fit_from(df)
  oof_cols   <- oof_model_cols_from(df)
  final_cols <- fit$x_cols
  score_cols <- scoring_x_cols_from(df, fit)

  # Every path actually synthesises `_isNA` companions now (the fix).
  expect_true(length(isna_of(oof_cols))   > 0L)
  expect_true(length(isna_of(final_cols)) > 0L)
  expect_true(length(isna_of(score_cols)) > 0L)

  # COUNT parity.
  expect_identical(length(oof_cols), length(final_cols))
  expect_identical(length(final_cols), length(score_cols))

  # SET parity.
  expect_setequal(oof_cols, final_cols)
  expect_setequal(final_cols, score_cols)
  expect_setequal(isna_of(oof_cols), isna_of(final_cols))
  expect_setequal(isna_of(final_cols), isna_of(score_cols))

  # ORDER parity (byte-identical column order).
  expect_identical(oof_cols, final_cols)
  expect_identical(final_cols, score_cols)

  # HASH parity: the single load-bearing acceptance assertion.
  h_oof   <- digest::digest(oof_cols)
  h_final <- digest::digest(final_cols)
  h_score <- digest::digest(score_cols)
  expect_identical(h_oof, h_final)
  expect_identical(h_final, h_score)
})

# ---------------------------------------------------------------------------
# (3) The single FINAL training procedure deploys the SAME shared recipe /
#     feature space (base + `_isNA`) as the standalone core. Driven through the
#     real FINAL engine on a GPKG built from the raw frame.
# ---------------------------------------------------------------------------
test_that("the single FINAL protocol deploys the shared feature space (one upstream frame)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("dplyr")

  df <- make_isna_upstream_df()
  sfc <- sf::st_sfc(lapply(seq_len(nrow(df)), function(i) {
    x <- (i %% 5); y <- (i %/% 5)
    sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1),
                              c(x, y + 1), c(x, y))))
  }), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  on.exit(unlink(g, force = TRUE), add = TRUE)
  sf::st_write(sf::st_sf(df, geometry = sfc), g, layer = "train_features",
               quiet = TRUE, delete_dsn = TRUE)

  engine <- g_ip("train_final_model_direct")
  res_nst <- suppressMessages(suppressWarnings(engine(
    labelled_gpkg = g, labelled_layer = "train_features",
    out_dir = NULL, overwrite = TRUE, verbose = FALSE, prefix = "nst",
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42L, group_col = "block_id", val_frac = 0.2,
    seed = 42L, nrounds_max = 8L, early_stopping_rounds = 4L,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = g_ip(".of_canonical_model_params")())))

  # The FINAL model deploys the SAME shared feature space (base + `_isNA`).
  expect_true(length(isna_of(res_nst$x_cols)) > 0L)
  # And it matches the standalone-core feature space derived above.
  expect_setequal(res_nst$x_cols, final_x_cols_from(df))
})

# ---------------------------------------------------------------------------
# (4) DETERMINISM of the boundary: re-running each stage's column derivation on
#     the SAME upstream frame is byte-identical.
# ---------------------------------------------------------------------------
test_that("each stage's deployed columns are deterministic for a fixed upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  df <- make_isna_upstream_df()
  expect_identical(oof_model_cols_from(df), oof_model_cols_from(df))
  expect_identical(final_x_cols_from(df),   final_x_cols_from(df))
})
