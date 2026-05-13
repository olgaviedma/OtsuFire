# 0.4.0 architectural refactor (Agent H): tests for the
# extract_features() ADD-only contract.
#
# Under 0.4.0, extract_features() preserves ALL input columns
# unchanged and only ADDS the 50 generated supervised features. The
# whitelist filter that controls what reaches the model matrix lives
# in train_final_model_direct() (see
# test-final-model-uses-whitelist.R). The legacy `strict_input` and
# `allow_passthrough_cols` arguments are no-ops since 0.4.0.

extract_features_internal <- function() {
  get("extract_features", envir = asNamespace("OtsuFire"))
}

residual_cols_internal <- function() {
  get(".deterministic_residual_cols", envir = asNamespace("OtsuFire"))
}

whitelist_internal <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

make_kb1_rect <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(matrix(
    c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
    ncol = 2, byrow = TRUE
  )))
}

make_kb1_fixture <- function(extra_train = NULL, extra_unl = NULL) {
  crs_epsg <- 3035

  train_folds <- sf::st_sf(
    fire_uid = c("t1", "t2", "t3"),
    class = c("burned", "unburned", "burned"),
    block_id = c(1L, 1L, 2L),
    fold_rep1 = c(1L, 1L, 2L),
    geometry = sf::st_sfc(
      make_kb1_rect(0, 0, 1, 1),
      make_kb1_rect(3, 0, 4, 1),
      make_kb1_rect(1, 3, 2, 4),
      crs = crs_epsg
    )
  )

  unlabeled <- sf::st_sf(
    fire_uid = c("u1", "u2", "u3"),
    geometry = sf::st_sfc(
      make_kb1_rect(0, 1, 1, 2),
      make_kb1_rect(3, 1, 4, 2),
      make_kb1_rect(4, 3, 5, 4),
      crs = crs_epsg
    )
  )

  if (!is.null(extra_train)) {
    for (nm in names(extra_train)) {
      train_folds[[nm]] <- extra_train[[nm]]
    }
  }
  if (!is.null(extra_unl)) {
    for (nm in names(extra_unl)) {
      unlabeled[[nm]] <- extra_unl[[nm]]
    }
  }

  template <- terra::rast(
    nrows = 5, ncols = 5,
    xmin = 0, xmax = 5, ymin = 0, ymax = 5,
    crs = "EPSG:3035"
  )
  rbr_summer <- template; terra::values(rbr_summer) <- seq_len(terra::ncell(template))
  dem <- template; terra::values(dem) <- seq(100, 124)
  slope <- template; terra::values(slope) <- seq(5, 29)
  corine_r <- template; terra::values(corine_r) <- rep(c(1, 1, 1, 2, 2), times = 5)

  list(
    train_folds = train_folds,
    unlabeled   = unlabeled,
    rbr_summer  = rbr_summer,
    dem = dem, slope = slope, corine_r = corine_r,
    cor_groups = list(natural = c(1), developed = c(2))
  )
}

run_kb1_extract <- function(strict_input = TRUE,
                            allow_passthrough_cols = character(0),
                            extra_train = NULL,
                            extra_unl = NULL) {
  fx <- make_kb1_fixture(extra_train = extra_train, extra_unl = extra_unl)
  ef <- extract_features_internal()
  ef(
    train_folds = fx$train_folds,
    unlabeled = fx$unlabeled,
    build_features = TRUE,
    rbr_summer = fx$rbr_summer,
    dem = fx$dem,
    slope = fx$slope,
    corine_r = fx$corine_r,
    cor_groups = fx$cor_groups,
    use_doy = FALSE,
    use_aw = FALSE,
    use_nbr = FALSE,
    use_hotspots = FALSE,
    use_ecoregions = FALSE,
    return_features = TRUE,
    return_features_geometry = FALSE,
    save_features_gpkg = FALSE,
    strict_input = strict_input,
    allow_passthrough_cols = allow_passthrough_cols,
    verbose = FALSE
  )
}

# --- T1 -------------------------------------------------------------
# 0.4.0 contract reversal: under strict_input = TRUE (default), every
# input column — including legacy deterministic-stage residuals — is
# preserved. The supervised model never sees them because they're
# filtered at train_final_model_direct() against the canonical
# whitelist, but extract_features() must NOT drop them.
test_that("T1: 0.4.0 strict_input preserves all input columns including residuals", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  residuals <- list(
    median_rbr = c(123, 234, 345),
    p_above_keep_q25 = c(0.1, 0.2, 0.3),
    p_above_keep_ref = c(0.5, 0.6, 0.7),
    percentile_in_keep = c(45, 55, 65),
    qa_changed = c(0L, 1L, 0L),
    n_pix = c(10L, 20L, 30L),
    area_ha = c(0.81, 1.62, 2.43)
  )
  out <- suppressMessages(run_kb1_extract(strict_input = TRUE,
                                            extra_train = residuals,
                                            extra_unl   = residuals))

  for (nm in names(residuals)) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("residual lost from train_feat:", nm))
    expect_true(nm %in% names(out$unl_feat),
                info = paste("residual lost from unl_feat:", nm))
  }
})

# --- T2 -------------------------------------------------------------
# 0.4.0: allow_passthrough_cols is a no-op since the function is
# add-only. Whatever is in the input survives.
test_that("T2: allow_passthrough_cols is a no-op since 0.4.0", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  residuals <- list(
    median_rbr = c(123, 234, 345),
    p_above_keep_q25 = c(0.1, 0.2, 0.3),
    n_pix = c(10L, 20L, 30L)
  )
  out <- suppressMessages(run_kb1_extract(
    strict_input = TRUE,
    allow_passthrough_cols = c("median_rbr"),
    extra_train = residuals,
    extra_unl   = residuals
  ))

  expect_true("median_rbr" %in% names(out$train_feat))
  expect_true("p_above_keep_q25" %in% names(out$train_feat))
  expect_true("n_pix" %in% names(out$train_feat))
})

# --- T3 -------------------------------------------------------------
# strict_input = FALSE has always been a passthrough mode; under
# 0.4.0 it remains a no-op (same as TRUE).
test_that("T3: strict_input = FALSE preserves all input columns", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  residuals <- list(
    median_rbr = c(123, 234, 345),
    p_above_keep_q25 = c(0.1, 0.2, 0.3),
    n_pix = c(10L, 20L, 30L),
    qa_changed = c(0L, 1L, 0L)
  )
  out <- suppressMessages(run_kb1_extract(
    strict_input = FALSE,
    extra_train = residuals,
    extra_unl   = residuals
  ))

  for (nm in names(residuals)) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("legacy mode dropped:", nm))
    expect_true(nm %in% names(out$unl_feat),
                info = paste("legacy mode dropped:", nm))
  }
})

# --- T5 -------------------------------------------------------------
# 0.4.0: input contract is symmetric — residuals on the unlabeled
# side also survive.
test_that("T5: residuals survive on the unlabeled side too", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  residuals <- list(median_rbr = c(123, 234, 345))
  out <- suppressMessages(run_kb1_extract(strict_input = TRUE,
                                            extra_unl = residuals))

  expect_false("median_rbr" %in% names(out$train_feat))  # not present in input
  expect_true("median_rbr" %in% names(out$unl_feat))     # preserved on unlabeled
})

# Bonus: the legacy deny-list constant retains the leak quartet as an
# audit log (no longer a runtime filter under 0.4.0, but useful for
# downstream auditors).
test_that(".deterministic_residual_cols still records the leak quartet", {
  rc <- residual_cols_internal()
  expect_true(all(c("median_rbr", "p_above_keep_q25",
                    "p_above_keep_ref", "percentile_in_keep") %in% rc))
})

# --- T_extra_a -----------------------------------------------------
# Administrative columns (added by the supervised pipeline itself)
# survive extract_features(). Was reversed in 0.3.1; remains correct
# under 0.4.0.
test_that("T_extra_a: administrative columns survive extract_features()", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  admin_cols <- list(
    neg_type                 = c("a", "b", "c"),
    cell_id                  = c(101L, 102L, 103L),
    training_selected        = c(TRUE, FALSE, TRUE),
    training_group           = c("g1", "g2", "g1"),
    training_reason          = c("r1", "r2", "r3"),
    intersects_deterministic = c(TRUE, TRUE, FALSE)
  )

  rc <- residual_cols_internal()
  for (nm in names(admin_cols)) {
    expect_false(nm %in% rc,
                 info = paste("admin col still in deny list:", nm))
  }

  out <- suppressWarnings(suppressMessages(run_kb1_extract(
    strict_input = TRUE,
    extra_train  = admin_cols,
    extra_unl    = admin_cols
  )))

  for (nm in names(admin_cols)) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("admin col dropped from train_feat:", nm))
    expect_true(nm %in% names(out$unl_feat),
                info = paste("admin col dropped from unl_feat:", nm))
  }
})

# --- T_extra_b -----------------------------------------------------
# 0.4.0 reversal: geometric columns are NOT in the supervised
# whitelist (sampling-induced size bias) but extract_features() is
# now ADD-only, so they SURVIVE through the feature-extraction stage.
# The whitelist filter at train_final_model_direct() drops them
# before the model matrix is built — that assertion lives in
# test-final-model-uses-whitelist.R.
test_that("T_extra_b: 0.4.0 — geometric columns survive extract_features()", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  geom_cols <- list(
    area_ha     = c(0.81, 1.62, 2.43),
    n_pix       = c(10L, 20L, 30L),
    log_area    = c(-0.21, 0.48, 0.89),
    perim_m     = c(40.0, 80.0, 120.0),
    compactness = c(0.7, 0.6, 0.5),
    elongation  = c(1.1, 1.3, 1.0),
    n_holes    = c(0L, 1L, 0L)
  )

  out <- suppressWarnings(suppressMessages(run_kb1_extract(
    strict_input = TRUE,
    extra_train  = geom_cols,
    extra_unl    = geom_cols
  )))

  for (nm in names(geom_cols)) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("geometric col lost from train_feat:", nm))
    expect_true(nm %in% names(out$unl_feat),
                info = paste("geometric col lost from unl_feat:", nm))
  }
})

# --- T_preserve_all ------------------------------------------------
# 0.4.0 architectural test (per Agent H brief §Phase 5a item a):
# build a fixture with a wide set of administrative + residual
# columns and assert ALL of them appear in the output, alongside the
# 50 generated features.
test_that("T_preserve_all: 0.4.0 — extract_features() preserves all input columns", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  bag <- list(
    # Administrative metadata (added by pool/fold/orchestrator).
    neg_type                 = c("a", "b", "c"),
    cell_id                  = c(101L, 102L, 103L),
    class_final              = c("keep", "drop", "keep"),
    legacy_decision          = c("keep", "drop", "keep"),
    intersects_deterministic = c(TRUE, TRUE, FALSE),
    raw_class                = c("burned", "unburned", "burned"),
    class_audited            = c("burned", "unburned", "burned"),
    training_selected        = c(TRUE, FALSE, TRUE),
    training_group           = c("g1", "g2", "g1"),
    training_reason          = c("r1", "r2", "r3"),
    poly_id                  = c("p1", "p2", "p3"),
    source                   = c("s1", "s2", "s3"),
    source_poly_id           = c("sp1", "sp2", "sp3"),
    fold_rep2                = c(1L, 2L, 1L),
    year                     = c(2005L, 2005L, 2005L),
    scenario                 = c("balanced", "balanced", "balanced"),
    # Geometric / sampling-biased columns (excluded from whitelist
    # but must survive extract_features()).
    area_ha                  = c(0.81, 1.62, 2.43),
    n_pix                    = c(10L, 20L, 30L),
    log_area                 = c(-0.21, 0.48, 0.89),
    perim_m                  = c(40.0, 80.0, 120.0),
    compactness              = c(0.7, 0.6, 0.5),
    elongation               = c(1.1, 1.3, 1.0),
    n_holes                 = c(0L, 1L, 0L),
    # Deterministic-stage residuals.
    median_rbr               = c(123, 234, 345),
    p_above_keep_q25         = c(0.1, 0.2, 0.3),
    p_above_keep_ref         = c(0.5, 0.6, 0.7),
    percentile_in_keep       = c(45, 55, 65),
    qa_changed               = c(0L, 1L, 0L),
    qa_final                 = c("yes", "no", "yes"),
    p_oof_mean               = c(0.5, 0.6, 0.7)
  )

  out <- suppressWarnings(suppressMessages(run_kb1_extract(
    strict_input = TRUE,
    extra_train  = bag,
    extra_unl    = bag
  )))

  # Every input column survives.
  for (nm in names(bag)) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("input col lost from train_feat:", nm))
    expect_true(nm %in% names(out$unl_feat),
                info = paste("input col lost from unl_feat:", nm))
  }

  # A sample of the generated features are present. CORINE column
  # names depend on the `cor_groups` argument passed to
  # extract_features(); the fixture uses {natural, developed}, so we
  # check those rather than the production whitelist names. The full
  # whitelist test lives in test-final-model-uses-whitelist.R.
  expected_added <- c("rbr_med", "rbr_p10", "rbr_p90", "rbr_iqr",
                      "elev_med", "slope_med",
                      "cor_natural_frac", "cor_developed_frac")
  for (nm in expected_added) {
    expect_true(nm %in% names(out$train_feat),
                info = paste("expected generated feature missing:", nm))
  }
})
