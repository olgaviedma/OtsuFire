# 0.5.0: tests for the `feature_weights` argument added to
# `train_final_model_direct()` and `run_dm_oof_pipeline()`. The
# argument takes a NAMED numeric vector. Names not in the active
# feature space are silently dropped with a warning. Negative,
# non-finite, non-numeric, or unnamed values raise an error.

whitelist_internal_fw <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

train_final_internal_fw <- function() {
  get("train_final_model_direct", envir = asNamespace("OtsuFire"))
}

# Synthetic GPKG (mirrors the override fixture's shape).
make_weights_fixture_gpkg <- function(n_burned = 14, n_neg_random = 7,
                                       n_neg_drop = 7, seed = 41L) {
  set.seed(seed)
  whitelist <- whitelist_internal_fw()
  n <- n_burned + n_neg_random + n_neg_drop

  feat_df <- as.data.frame(
    matrix(stats::runif(n * length(whitelist)),
           nrow = n, dimnames = list(NULL, whitelist))
  )
  for (nm in c("hotspot_available", "hs_in_poly", "hs_in_buffer",
                "hs_used_n", "hs_hiConf_n", "hs_support_present",
                "hs_no_support_when_available", "hs_only_buffer_support")) {
    if (nm %in% names(feat_df)) feat_df[[nm]] <- as.integer(round(feat_df[[nm]]))
  }
  admin_df <- data.frame(
    fire_uid    = sprintf("uid_%03d", seq_len(n)),
    class       = c(rep("burned", n_burned),
                     rep("unburned", n_neg_random + n_neg_drop)),
    # GATE 6.5 (2026-06-12): deterministic drops are not training negatives; the
    # former "drop" negatives are now Otsu-residual negatives (valid bucket).
    source      = c(rep("burned_truth", n_burned),
                     rep("random_burnable_background", n_neg_random),
                     rep("otsu_patch_residual", n_neg_drop)),
    neg_type    = c(rep(NA_character_, n_burned),
                     rep("background_cell", n_neg_random),
                     rep("otsu_patch_drop", n_neg_drop)),
    block_id    = rep(c(1L, 2L, 3L, 4L), length.out = n),
    fold_rep1   = rep(c(1L, 2L), length.out = n),
    fold_rep2   = rep(c(2L, 1L), length.out = n),
    poly_id     = sprintf("p_%03d", seq_len(n)),
    stringsAsFactors = FALSE
  )
  full_df <- cbind(admin_df, feat_df)
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

run_train_with_weights <- function(gpkg, weights = NULL, override = NULL,
                                    prefix = "fw_smoke", out_dir = NULL) {
  if (is.null(out_dir)) out_dir <- tempfile("fw_train_")
  fn <- train_final_internal_fw()
  args <- list(
    labelled_gpkg              = gpkg,
    labelled_layer             = "train_features",
    out_dir                    = out_dir,
    prefix                     = prefix,
    overwrite                  = TRUE,
    verbose                    = FALSE,
    nrounds_max                = 8L,
    early_stopping_rounds      = 4L,
    random_to_burned_ratio                 = 1,
    otsu_unburned_to_burned_ratio          = 1,
    # Gate 1B (2026-06-07): the engine now requires resolved methodological args.
    sampling_seed              = 42,
    seed                       = 42,
    val_frac                   = 0.15,
    group_col                  = "block_id",
    impute_numeric             = "median",
    impute_factor_missing      = "MISSING",
    model_params_base          = OtsuFire:::.of_canonical_model_params(),
    feature_weights            = weights,
    feature_whitelist_override = override
  )
  res <- suppressMessages(do.call(fn, args))
  list(res = res, out_dir = out_dir)
}

# --- T6: default (no weights) --------------------------------------
test_that("default behaviour (no weights) leaves feature_weights NULL", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  gpkg <- make_weights_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  out <- run_train_with_weights(gpkg)
  on.exit(unlink(out$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(out$res$files$recipe_rds)
  expect_true(is.null(recipe$feature_weights) ||
                length(recipe$feature_weights) == 0L,
              info = "recipe$feature_weights should be empty when not used")

  meta_lines <- readLines(out$res$files$meta_txt)
  expect_true(any(grepl("feature_weights_applied: FALSE", meta_lines,
                          fixed = TRUE)))
})

# --- T7: valid named weight applied ---------------------------------
test_that("valid feature_weights applied and persisted in recipe + meta", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  gpkg <- make_weights_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  weights <- c(hs_in_poly = 0.05)
  out <- run_train_with_weights(gpkg, weights = weights)
  on.exit(unlink(out$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(out$res$files$recipe_rds)
  expect_true(is.list(recipe$feature_weights))
  expect_equal(recipe$feature_weights$applied$hs_in_poly, 0.05)
  expect_equal(length(recipe$feature_weights$applied), 1L)

  meta_lines <- readLines(out$res$files$meta_txt)
  expect_true(any(grepl("feature_weights_applied: TRUE", meta_lines,
                          fixed = TRUE)))
  expect_true(any(grepl("hs_in_poly=0.05", meta_lines, fixed = TRUE)))
})

# --- T8: negative weight ---------------------------------------------
test_that("negative weight errors", {
  gpkg <- "<unused>"
  fn <- train_final_internal_fw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      feature_weights = c(hs_in_poly = -0.1)
    )),
    regexp = "feature_weights.+finite and >= 0"
  )
})

# --- T9: non-numeric weight ------------------------------------------
test_that("non-numeric feature_weights errors", {
  fn <- train_final_internal_fw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = "<unused>",
      labelled_layer = "train_features",
      feature_weights = c(hs_in_poly = "low")
    )),
    regexp = "feature_weights.+named numeric"
  )
})

# --- T10: unnamed numeric vector errors -------------------------------
test_that("unnamed feature_weights vector errors", {
  fn <- train_final_internal_fw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = "<unused>",
      labelled_layer = "train_features",
      feature_weights = c(0.05, 0.10)
    )),
    regexp = "feature_weights.+named numeric"
  )
})

# --- T11: names not in the active feature space ----------------------
test_that("unknown weight names trigger warning + are dropped silently", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  gpkg <- make_weights_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  weights <- c(hs_in_poly = 0.05, this_does_not_exist = 0.5)
  out_dir <- tempfile("fw_unknown_")
  fn <- train_final_internal_fw()
  expect_warning(
    suppressMessages(fn(
      labelled_gpkg              = gpkg,
      labelled_layer             = "train_features",
      out_dir                    = out_dir,
      prefix                     = "fw_unknown",
      overwrite                  = TRUE,
      verbose                    = FALSE,
      nrounds_max                = 8L,
      early_stopping_rounds      = 4L,
      random_to_burned_ratio                 = 1,
      otsu_unburned_to_burned_ratio          = 1,
      sampling_seed              = 42,
      seed                       = 42,
      val_frac                   = 0.15,
      group_col                  = "block_id",
      impute_numeric             = "median",
      impute_factor_missing      = "MISSING",
      model_params_base          = OtsuFire:::.of_canonical_model_params(),
      feature_weights            = weights
    )),
    regexp = "this_does_not_exist"
  )
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(file.path(out_dir, "fw_unknown_recipe.rds"))
  # The unknown name is recorded as dropped; the valid one is applied.
  expect_true("this_does_not_exist" %in%
                recipe$feature_weights$unknown_dropped)
  expect_equal(recipe$feature_weights$applied$hs_in_poly, 0.05)
})

# --- T12: weights for features removed by the override --------------
test_that("weights for override-removed features warn + drop", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  gpkg <- make_weights_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  whitelist <- whitelist_internal_fw()
  override <- setdiff(whitelist,
                       c("hs_in_poly", "hs_in_buffer", "hs_min_dist_m",
                         "hs_support_present"))

  # Weight a feature that the override removes — should warn + drop.
  weights <- c(hs_in_poly = 0.05, rbr_med = 0.5)

  out_dir <- tempfile("fw_combined_")
  fn <- train_final_internal_fw()
  expect_warning(
    suppressMessages(fn(
      labelled_gpkg              = gpkg,
      labelled_layer             = "train_features",
      out_dir                    = out_dir,
      prefix                     = "fw_combined",
      overwrite                  = TRUE,
      verbose                    = FALSE,
      nrounds_max                = 8L,
      early_stopping_rounds      = 4L,
      random_to_burned_ratio                 = 1,
      otsu_unburned_to_burned_ratio          = 1,
      sampling_seed              = 42,
      seed                       = 42,
      val_frac                   = 0.15,
      group_col                  = "block_id",
      impute_numeric             = "median",
      impute_factor_missing      = "MISSING",
      model_params_base          = OtsuFire:::.of_canonical_model_params(),
      feature_whitelist_override = override,
      feature_weights            = weights
    )),
    regexp = "hs_in_poly"
  )
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(file.path(out_dir, "fw_combined_recipe.rds"))
  expect_true("hs_in_poly" %in% recipe$feature_weights$unknown_dropped)
  expect_equal(recipe$feature_weights$applied$rbr_med, 0.5)
  expect_false("hs_in_poly" %in% recipe$cols$feature_cols)
  expect_true("rbr_med" %in% recipe$cols$feature_cols)
})
