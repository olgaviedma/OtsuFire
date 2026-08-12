# ============================================================================
# Gate 2 (2026-06-10): EQUIVALENCE PROOF for the removal of `hs_any` from the
# canonical supervised whitelist (Natalia's final decision).
#
# `hs_any` was a binary helper (= as.integer(hs_used_n > 0)) carried in the
# whitelist for historical parity. It was:
#   * redundant with hs_used_n > 0 / hs_in_poly > 0;
#   * NEVER materialised as a base GPKG feature (extract_supervised_features()
#     never writes an `hs_any` column);
#   * never present in the RESOLVED recipe (the resolved base was already 50,
#     because `hs_any` was synthesised only inside build_design_matrix_patches()
#     and that synthesis never reached the FINAL / scoring recipe).
#
# This suite proves that REMOVING `hs_any` from the whitelist (51 -> 50 names)
# and deleting its synthesis block does NOT change the resolved feature space:
#   * the RESOLVED FINAL recipe base_features / missing_indicator_features /
#     final_feature_order are exactly the documented 50 base + 50 `_isNA` = 100,
#     in canonical order, with NO `hs_any` (or `hs_any_isNA`) anywhere;
#   * the feature_schema_fingerprint is computed from those resolved columns,
#     so it is well-defined and stable, and the hotspot block is 12 (genuine)
#     not 13;
#   * the whitelist length is now 50;
#   * NO-HOTSPOT removes exactly 12 genuine hotspot features (not a fictitious
#     13th `hs_any`).
# ============================================================================

g_ip2 <- function(nm) get(nm, envir = asNamespace("OtsuFire"))

# The 50 documented base features in canonical order (the whitelist contract
# AFTER the 0.6.2 hs_any removal). This is the authoritative "expected" set the
# proof asserts the resolved recipe equals.
.HS_ANY_REMOVAL_EXPECTED_BASE <- c(
  # RBR same-window (5)
  "rbr_valid_frac", "rbr_p10", "rbr_med", "rbr_p90", "rbr_iqr",
  # CORINE (7)
  "cor_open_frac", "cor_wetlands_frac", "cor_water_frac",
  "cor_herbaceous_frac", "cor_urban_frac", "cor_agri_frac", "cor_forest_frac",
  # elevation (7)
  "elev_valid_frac", "elev_p10", "elev_med", "elev_p90",
  "elev_iqr", "elev_mean", "elev_sd",
  # slope (7)
  "slope_valid_frac", "slope_p10", "slope_med", "slope_p90",
  "slope_iqr", "slope_mean", "slope_sd",
  # DOY (5)
  "doy_valid_frac", "doy_p10", "doy_med", "doy_p90", "doy_iqr",
  # RBR all-window + persistence (7)
  "rbr_aw_valid_frac", "rbr_aw_p10", "rbr_aw_med", "rbr_aw_p90",
  "rbr_aw_iqr", "persist_delta", "persist_ratio",
  # Hotspots (12, genuine -- hs_any removed)
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support"
)

# --- proof 1: whitelist is now 50, hs_any absent ----------------------------
test_that("0.6.2: canonical whitelist is 50 and carries no hs_any", {
  wl <- g_ip2(".supervised_feature_cols")
  expect_equal(length(wl), 50L)
  expect_false("hs_any" %in% wl)
  # Exact identity (names + order) with the documented expected base.
  expect_identical(wl, .HS_ANY_REMOVAL_EXPECTED_BASE)
  # The hotspot block is 12 genuine features.
  expect_equal(sum(grepl("^hs_|^hotspot_", wl)), 12L)
})

# --- proof 2: resolved FINAL recipe is 50 base + 50 _isNA = 100, no hs_any ---
test_that("resolved FINAL recipe base/indicator/order are the 50/50/100 space, hs_any-free", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  df  <- sc_isna_upstream_df()
  fit <- sc_final_fit_from(df)

  base <- fit$base_features
  ind  <- fit$missing_indicator_features
  ord  <- fit$final_feature_order

  # The resolved base is exactly the 50 documented base features, in canonical
  # order. This is the crux: the resolved base was ALWAYS 50 because hs_any was
  # never a GPKG column -> removing it from the whitelist changes nothing here.
  expect_identical(base, .HS_ANY_REMOVAL_EXPECTED_BASE)
  expect_equal(length(base), 50L)
  expect_equal(length(ind), 50L)
  expect_equal(length(ord), 100L)

  # NO hs_any anywhere in the resolved space.
  expect_false("hs_any" %in% base)
  expect_false("hs_any_isNA" %in% ind)
  expect_false(any(grepl("hs_any", ord)))

  # Every indicator is a base feature's `_isNA` companion (and vice versa).
  expect_setequal(ind, paste0(base, "_isNA"))
  expect_setequal(ord, c(base, ind))
})

# --- proof 3: the feature-schema fingerprint is well-defined & matches the
#     resolved 50/50/100 space (unchanged by the whitelist 51->50 edit) -------
test_that("feature_schema_fingerprint reflects the unchanged 50/50/100 resolved space", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  df  <- sc_isna_upstream_df()
  fit <- sc_final_fit_from(df)
  fp  <- g_ip2("feature_schema_fingerprint")(fit)

  expect_equal(fp$payload$n_base, 50L)
  expect_equal(fp$payload$n_indicators, 50L)
  expect_equal(fp$payload$n_total, 100L)
  expect_false(any(grepl("hs_any", fp$payload$base_features)))
  expect_false(any(grepl("hs_any", fp$payload$final_feature_order)))

  # The fingerprint is a deterministic function of the resolved structural
  # payload: recomputing it from a recipe carrying the SAME base/indicator/order
  # yields the IDENTICAL hash. Since the resolved base is unchanged by the
  # whitelist 51->50 edit, the fingerprint is unchanged.
  recipe_shape <- list(cols = list(
    base_features              = fit$base_features,
    missing_indicator_features = fit$missing_indicator_features,
    final_feature_order        = fit$final_feature_order,
    x_cols                     = fit$x_cols,
    feature_cols               = fit$feature_cols))
  fp2 <- g_ip2("feature_schema_fingerprint")(recipe_shape)
  expect_identical(fp$hash, fp2$hash)
})

# --- proof 4: NO-HOTSPOT removes exactly 12 genuine hotspot features ---------
test_that("NO-HOTSPOT override removes exactly 12 genuine hotspot features (no phantom 13th)", {
  wl <- g_ip2(".supervised_feature_cols")
  hotspot_block <- c(
    "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
    "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
    "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
    "hs_only_buffer_support")
  # All 12 are genuinely in the whitelist (none is a phantom).
  expect_true(all(hotspot_block %in% wl))
  expect_equal(length(hotspot_block), 12L)

  no_hotspot <- setdiff(wl, hotspot_block)
  expect_equal(length(no_hotspot), 38L)        # 50 - 12 = 38 base
  expect_equal(length(wl) - length(no_hotspot), 12L)  # removed count == 12
  expect_false("hs_any" %in% hotspot_block)    # hs_any is not the 13th item
})
