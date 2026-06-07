# B5 (OtsuFire, 2026-06-06): for years/polygons with NO hotspot data,
# the features extractor's out_nodata() must emit NA (not the literal
# -9999) for the MEASURED hotspot columns, so the design-matrix builder
# flags them via `<col>_isNA = 1` and imputes them to the column median.
# `hotspot_available` stays a real 0L flag.
#
# out_nodata() is a closure inside extract_features() and cannot be
# called in isolation, so we (1) assert the documented out_nodata()
# column->value contract via a small reference table that mirrors the
# source, and (2) drive the real consumer (build_design_matrix_patches)
# with an all-missing hotspot block and assert the isNA flag +
# imputation behaviour, including the all-NA-column edge case.

# Reference table mirroring R/internal-sup-extract-features.R out_nodata().
# Measured quantities -> NA; structural/availability-conditioned flags -> 0L.
.b5_out_nodata_expected <- list(
  hotspot_available            = 0L,        # KEPT: real "no data this year" flag
  hs_in_poly                   = NA_real_,  # measured count
  hs_in_buffer                 = NA_real_,  # measured count
  hs_used_n                    = NA_real_,  # measured count
  hs_min_dist_m                = NA_real_,  # measured distance
  hs_frp_sum                   = NA_real_,  # measured FRP
  hs_frp_max                   = NA_real_,  # measured FRP
  hs_conf_mean                 = NA_real_,  # measured confidence
  hs_hiConf_n                  = NA_real_,  # measured count
  hs_support_present           = NA_real_,  # support existence: undefined w/o data
  hs_no_support_when_available = 0L,        # KEPT: conditioned on availability
  hs_only_buffer_support       = NA_real_   # buffer-only support: undefined w/o data
)

.b5_measured_cols <- c(
  "hs_in_poly", "hs_in_buffer", "hs_used_n", "hs_min_dist_m",
  "hs_frp_sum", "hs_frp_max", "hs_conf_mean", "hs_hiConf_n",
  "hs_support_present", "hs_only_buffer_support"
)

test_that("B5: out_nodata() contract is NA for measured hs_* and 0L for structural flags", {
  exp <- .b5_out_nodata_expected

  # Structural / availability-conditioned flags kept as real 0L.
  expect_identical(exp$hotspot_available, 0L)
  expect_identical(exp$hs_no_support_when_available, 0L)

  # All measured quantities + support-existence flags are NA (missing),
  # NOT the legacy -9999 sentinel.
  for (nm in .b5_measured_cols) {
    expect_true(is.na(exp[[nm]]),
                info = paste0("out_nodata() column ", nm, " must be NA"))
    expect_false(identical(exp[[nm]], -9999),
                 info = paste0("out_nodata() column ", nm, " must not be -9999"))
  }
})

test_that("B5: matrix builder flags NA hotspot columns via _isNA and imputes (no crash on all-NA)", {
  skip_if_not_installed("Matrix")

  fn <- get("build_design_matrix_patches", envir = asNamespace("OtsuFire"))

  n <- 6L
  # Labelled frame mixes hotspot-year rows (real values) with a
  # hotspot-less polygon (NA on every measured hs_* column), exactly as
  # out_nodata() now emits. hs_frp_sum is ALL-NA across both frames to
  # exercise the all-NA-column median fallback (median NA -> 0).
  L_df <- data.frame(
    fire_uid = sprintf("L%02d", seq_len(n)),
    poly_id  = sprintf("p%02d", seq_len(n)),
    class    = rep(c("burned", "unburned"), each = n / 2),
    block_id = rep(c(1L, 2L), each = n / 2),
    fold_rep1 = rep(c(1L, 2L), n / 2),
    rbr_med  = c(10, 12, 14, 11, 13, 15),
    # hotspot block: last row is the hotspot-less (out_nodata) polygon.
    hotspot_available = c(1L, 1L, 1L, 1L, 1L, 0L),
    hs_in_poly        = c(2, 0, 1, 3, 0, NA_real_),
    hs_used_n         = c(2, 0, 1, 3, 0, NA_real_),
    hs_min_dist_m     = c(0, 500, 0, 0, 800, NA_real_),
    hs_conf_mean      = c(0.9, 0, 0.7, 0.8, 0, NA_real_),
    hs_frp_max        = c(5, 0, 3, 8, 0, NA_real_),
    hs_frp_sum        = rep(NA_real_, n),   # all-NA column edge case
    stringsAsFactors = FALSE
  )
  BL_df <- data.frame(
    fire_uid = sprintf("BL%02d", seq_len(n)),
    poly_id  = sprintf("q%02d", seq_len(n)),
    class    = rep("unburned", n),
    block_id = NA_integer_,
    fold_rep1 = NA_integer_,
    rbr_med  = c(10, 12, 14, 11, 13, 15),
    hotspot_available = c(1L, 1L, 0L, 1L, 0L, 0L),
    hs_in_poly        = c(1, 0, NA_real_, 2, NA_real_, NA_real_),
    hs_used_n         = c(1, 0, NA_real_, 2, NA_real_, NA_real_),
    hs_min_dist_m     = c(0, 400, NA_real_, 0, NA_real_, NA_real_),
    hs_conf_mean      = c(0.6, 0, NA_real_, 0.7, NA_real_, NA_real_),
    hs_frp_max        = c(4, 0, NA_real_, 6, NA_real_, NA_real_),
    hs_frp_sum        = rep(NA_real_, n),   # all-NA column edge case
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(fn(
    labelled = L_df,
    burned_like = BL_df,
    id_cols = c("fire_uid", "class", "poly_id", "block_id", "fold_rep1"),
    drop_regex = c("^fold_rep", "^block_id$", "^poly_id$"),
    cat_cols = character(0),
    # use defaults for hs_n_col/hs_conf_col/hs_frp_col so the hotspot
    # special-casing path runs as in production.
    median_from = "labelled",
    return_prepared_df = TRUE,
    verbose = FALSE
  ))

  Xp <- out$X_all_prepared
  cn <- colnames(out$XL_mat)

  # The hotspot-less rows are the last row of each frame:
  #   labelled last row -> index n; burned_like last row -> index 2n.
  na_rows <- c(n, 2L * n)

  # (a) _isNA companions exist and are 1 exactly on the missing rows.
  for (nm in c("hs_in_poly", "hs_used_n", "hs_min_dist_m",
               "hs_conf_mean", "hs_frp_max")) {
    flag <- paste0(nm, "_isNA")
    expect_true(flag %in% names(Xp),
                info = paste0("missing _isNA companion for ", nm))
    expect_equal(Xp[[flag]][na_rows], c(1L, 1L),
                 info = paste0(flag, " must be 1 on hotspot-less rows"))
    # No NA survives into the imputed value column.
    expect_false(any(is.na(Xp[[nm]])),
                 info = paste0(nm, " still contains NA after imputation"))
  }

  # (b) all-NA column (hs_frp_sum): flag is all 1, value imputed to 0
  #     (median of all-NA -> 0 fallback), no NaN, present in matrix.
  expect_true("hs_frp_sum_isNA" %in% names(Xp))
  expect_true(all(Xp$hs_frp_sum_isNA == 1L))
  expect_true(all(Xp$hs_frp_sum == 0))
  expect_false(any(is.na(Xp$hs_frp_sum)))
  expect_false(any(is.nan(Xp$hs_frp_sum)))

  # (c) hotspot_available stays a real 0/1 flag, never NA. The generic _isNA
  #     loop still emits a companion column, but because hotspot_available is
  #     never NA that companion must be all-0 (never flagged as missing).
  expect_false(any(is.na(Xp$hotspot_available)))
  if ("hotspot_available_isNA" %in% cn) {
    expect_true(all(out$XL_mat[, "hotspot_available_isNA"] == 0))
  }

  # (d) the sparse model matrix built cleanly (no NaN/Inf reaching xgb.DMatrix).
  expect_false(any(is.nan(out$XL_mat@x)))
  expect_false(any(is.infinite(out$XL_mat@x)))
  expect_false(any(is.nan(out$XBL_mat@x)))
})
