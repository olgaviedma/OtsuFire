# 0.4.0 architectural refactor (Agent H): tests for the whitelist
# filter inside train_final_model_direct().
#
# Under 0.4.0 the supervised model sees ONLY the columns in
# `.supervised_feature_cols` (plus their `_isNA` companions when
# present in the input). Administrative metadata, fold assignments,
# deterministic-stage residuals, geometric sampling-biased features,
# OOF outputs, and labels are all dropped at the model-matrix step
# regardless of what `extra_drop_cols` says.

whitelist_internal <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

train_final_internal <- function() {
  get("train_final_model_direct", envir = asNamespace("OtsuFire"))
}

residual_cols_internal <- function() {
  get(".deterministic_residual_cols", envir = asNamespace("OtsuFire"))
}

# Build a synthetic training-pool gpkg with the 50 whitelist features,
# administrative metadata, and known residuals. Returns the gpkg
# path; the caller is responsible for cleanup if desired.
make_whitelist_fixture_gpkg <- function(n_burned = 14, n_neg_random = 7,
                                          n_neg_drop = 7, seed = 17L) {
  set.seed(seed)
  whitelist <- whitelist_internal()
  n <- n_burned + n_neg_random + n_neg_drop

  feat_df <- as.data.frame(
    matrix(stats::runif(n * length(whitelist)),
           nrow = n, dimnames = list(NULL, whitelist))
  )

  # `hotspot_available` and the four hotspot-related flags are
  # integer 0/1 in the real pipeline; coerce for realism.
  for (nm in c("hotspot_available", "hs_in_poly", "hs_in_buffer",
                "hs_used_n", "hs_hiConf_n", "hs_support_present",
                "hs_no_support_when_available", "hs_only_buffer_support")) {
    if (nm %in% names(feat_df)) feat_df[[nm]] <- as.integer(round(feat_df[[nm]]))
  }

  # Administrative metadata that MUST survive in the GPKG and MUST
  # NOT enter the model matrix.
  admin_df <- data.frame(
    fire_uid                 = sprintf("uid_%03d", seq_len(n)),
    class                    = c(rep("burned", n_burned),
                                  rep("unburned", n_neg_random + n_neg_drop)),
    source                   = c(rep("burned_truth", n_burned),
                                  rep("random_burnable_background", n_neg_random),
                                  rep("deterministic_drop_hard", n_neg_drop)),
    neg_type                 = c(rep(NA_character_, n_burned),
                                  rep("background_cell", n_neg_random),
                                  rep("geo_excluded_hot", n_neg_drop)),
    block_id                 = rep(c(1L, 2L, 3L, 4L), length.out = n),
    fold_rep1                = rep(c(1L, 2L), length.out = n),
    fold_rep2                = rep(c(2L, 1L), length.out = n),
    poly_id                  = sprintf("p_%03d", seq_len(n)),
    cell_id                  = seq_len(n) + 1000L,
    intersects_deterministic = sample(c(TRUE, FALSE), n, replace = TRUE),
    legacy_decision          = sample(c("keep", "drop"), n, replace = TRUE),
    class_final              = c(rep("burned", n_burned),
                                  rep("drop", n_neg_random + n_neg_drop)),
    raw_class                = c(rep("burned", n_burned),
                                  rep("unburned", n_neg_random + n_neg_drop)),
    year                     = rep(2005L, n),
    scenario                 = rep("balanced", n),
    # Residuals from deterministic stage that must NOT enter model.
    median_rbr               = stats::runif(n) * 1000,
    p_above_keep_q25         = stats::runif(n),
    p_above_keep_ref         = stats::runif(n),
    percentile_in_keep       = stats::runif(n) * 100,
    qa_changed               = sample(0:1, n, replace = TRUE),
    qa_final                 = sample(c("yes", "no"), n, replace = TRUE),
    p_oof_mean               = stats::runif(n),
    # Geometric sampling-biased columns.
    area_ha                  = stats::runif(n, 0.5, 5),
    n_pix                    = sample(10:200, n, replace = TRUE),
    log_area                 = stats::runif(n, -1, 2),
    perim_m                  = stats::runif(n, 30, 500),
    compactness              = stats::runif(n),
    elongation               = stats::runif(n, 1, 3),
    n_holes                 = sample(0:2, n, replace = TRUE),
    stringsAsFactors         = FALSE
  )

  full_df <- cbind(admin_df, feat_df)

  # Build trivial polygons (all distinct, not overlapping).
  rect <- function(i) sf::st_polygon(list(matrix(c(
    i, 0, i + 0.5, 0,
    i + 0.5, 0.5, i, 0.5, i, 0
  ), ncol = 2, byrow = TRUE)))
  geom <- sf::st_sfc(lapply(seq_len(n), rect), crs = 3035)
  sf_obj <- sf::st_sf(full_df, geometry = geom)

  gpkg_path <- tempfile(fileext = ".gpkg")
  sf::st_write(sf_obj, gpkg_path, layer = "train_features", quiet = TRUE)
  gpkg_path
}

# --- T_whitelist_complete ----------------------------------------
# `.supervised_feature_cols` has the 50 entries Natalia approved on
# 2026-05-09 (50 in 0.4.0; hs_any restored in 0.4.1, then REMOVED again
# in 0.6.2 as a redundant, never-materialised helper). The constant is
# the contract; this test guards inadvertent edits.
test_that(".supervised_feature_cols has the 50 approved entries", {
  whitelist <- whitelist_internal()
  expect_equal(length(whitelist), 50L)

  # Block tally: 5 (RBR sw) + 7 (CORINE) + 7 (elev) + 7 (slope) +
  # 5 (DOY) + 7 (RBR aw + persistence) + 12 (hotspots) = 50.
  expect_true(all(c("rbr_valid_frac", "rbr_med", "rbr_iqr") %in% whitelist))
  expect_true(all(c("cor_open_frac", "cor_forest_frac",
                    "cor_water_frac") %in% whitelist))
  expect_true(all(c("elev_med", "elev_sd") %in% whitelist))
  expect_true(all(c("slope_med", "slope_sd") %in% whitelist))
  expect_true(all(c("doy_med", "doy_iqr") %in% whitelist))
  expect_true(all(c("rbr_aw_med", "persist_delta",
                    "persist_ratio") %in% whitelist))
  expect_true(all(c("hotspot_available", "hs_in_poly",
                    "hs_used_n") %in% whitelist))

  # 0.6.2: hs_any was removed from the whitelist (redundant with
  # hs_used_n > 0; never materialised as a base GPKG feature).
  expect_false("hs_any" %in% whitelist)

  # No ecoregion features in this build.
  expect_false(any(grepl("^eco_", whitelist)))

  # No fold-assignment / administrative columns.
  expect_false("block_id" %in% whitelist)
  expect_false("fold_rep1" %in% whitelist)
  expect_false("fold_rep2" %in% whitelist)
  expect_false("neg_type" %in% whitelist)
  expect_false("class" %in% whitelist)
  expect_false("fire_uid" %in% whitelist)
})

# --- T_extra_drop_removed_0_5_0 ----------------------------------
# 0.5.0: `extra_drop_cols` was removed entirely. The argument is no
# longer in the signature; passing it via `...` triggers the hard
# migration error.
test_that("extra_drop_cols is no longer in the train_final_model_direct signature", {
  fn <- train_final_internal()
  expect_false("extra_drop_cols" %in% names(formals(fn)))
})

# --- T_recipe_uses_whitelist -------------------------------------
# End-to-end assertion: train_final_model_direct() filters its
# feature set to the whitelist before constructing the model
# matrix. The recipe persists `feature_cols` and `x_cols`; both
# must be subsets of `.supervised_feature_cols ∪ _isNA companions`.
test_that("recipe$cols$feature_cols and x_cols are within the whitelist", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  whitelist <- whitelist_internal()
  rc <- residual_cols_internal()
  gpkg <- make_whitelist_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- train_final_internal()
  out_dir <- tempfile("whitelist_train_")
  res <- suppressMessages(suppressWarnings(fn(
    labelled_gpkg = gpkg,
    labelled_layer = "train_features",
    out_dir = out_dir,
    prefix = "whitelist_smoke",
    overwrite = TRUE,
    verbose = FALSE,
    nrounds_max = 8L,
    early_stopping_rounds = 4L,
    contextual_exclusion_to_burned_ratio = 1,
    spectral_hard_negative_to_burned_ratio = 1,
    random_to_burned_ratio = 1,
    otsu_unburned_to_burned_ratio = 1,
    # Gate 1B (2026-06-07): the engine now requires resolved methodological args.
    sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = OtsuFire:::.of_canonical_model_params()
  )))

  recipe <- readRDS(res$files$recipe_rds)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  feat_cols <- recipe$cols$feature_cols
  x_cols    <- recipe$cols$x_cols
  allowed   <- c(whitelist, paste0(whitelist, "_isNA"))

  # Every column the model trained on is whitelist-admissible.
  expect_true(all(feat_cols %in% allowed),
              info = paste("non-whitelist feature_cols:",
                            paste(setdiff(feat_cols, allowed),
                                  collapse = ", ")))
  expect_true(all(x_cols %in% allowed),
              info = paste("non-whitelist x_cols:",
                            paste(setdiff(x_cols, allowed),
                                  collapse = ", ")))

  # 0.6.2: hs_any was REMOVED from the whitelist (redundant with
  # hs_used_n > 0; never a base GPKG feature) and its synthesis block
  # was deleted, so it must NOT appear in the recipe column set.
  expect_false("hs_any" %in% feat_cols,
               info = "hs_any leaked into feature_cols after 0.6.2 removal")
  expect_false("hs_any" %in% x_cols,
               info = "hs_any leaked into x_cols after 0.6.2 removal")

  # Specific assertions: known administrative / residual columns
  # are NOT in the recipe column set.
  for (nm in c("neg_type", "fire_uid", "block_id", "fold_rep1",
                "fold_rep2", "class_final", "median_rbr",
                "p_above_keep_q25", "qa_changed", "p_oof_mean",
                "area_ha", "n_pix", "legacy_decision",
                "intersects_deterministic")) {
    expect_false(nm %in% feat_cols,
                  info = paste("admin/residual leaked into feature_cols:", nm))
    expect_false(nm %in% x_cols,
                  info = paste("admin/residual leaked into x_cols:", nm))
  }

  # The deterministic_residual_cols deny list and the recipe column
  # set are disjoint (the legacy audit-log invariant continues to
  # hold under 0.4.0).
  expect_equal(intersect(feat_cols, rc), character(0))

  # extra_drop_cols deprecation is silent for default callers (no
  # exception raised here) — the absence of a stop() is the test.
  expect_true(file.exists(res$files$model_rds))
})

# --- T_extra_drop_migration_error --------------------------------
# 0.5.0: passing `extra_drop_cols` is a hard error with the
# migration message pointing at `feature_whitelist_override`.
test_that("extra_drop_cols triggers the 0.5.0 migration error", {
  skip_if_not_installed("sf")

  gpkg <- make_whitelist_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- train_final_internal()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      extra_drop_cols = c("rbr_med", "elev_med")
    )),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      additional_drop_cols = c("rbr_med")
    )),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )
})
