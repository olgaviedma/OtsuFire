# 0.5.0: tests for the `feature_whitelist_override` argument added to
# `train_final_model_direct()` and `run_dm_oof_pipeline()`. The override
# is strict: any name not in the canonical `.supervised_feature_cols`
# raises an error. The deprecated `additional_drop_cols` /
# `extra_drop_cols` arguments raise the migration error.

whitelist_internal_tw <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

train_final_internal_tw <- function() {
  get("train_final_model_direct", envir = asNamespace("OtsuFire"))
}

oof_pipeline_internal_tw <- function() {
  get("run_dm_oof_pipeline", envir = asNamespace("OtsuFire"))
}

# Synthetic GPKG with all canonical features + admin columns.
make_override_fixture_gpkg <- function(n_burned = 14, n_neg_random = 7,
                                        n_neg_drop = 7, seed = 31L) {
  set.seed(seed)
  whitelist <- whitelist_internal_tw()
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
    stringsAsFactors         = FALSE
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

run_train_smoke <- function(gpkg, override = NULL,
                              prefix = "wlist_smoke", out_dir = NULL) {
  if (is.null(out_dir)) out_dir <- tempfile("wlist_train_")
  fn <- train_final_internal_tw()
  res <- suppressMessages(suppressWarnings(fn(
    labelled_gpkg              = gpkg,
    labelled_layer             = "train_features",
    out_dir                    = out_dir,
    prefix                     = prefix,
    overwrite                  = TRUE,
    verbose                    = FALSE,
    nrounds_max                = 8L,
    early_stopping_rounds      = 4L,
    contextual_exclusion_to_burned_ratio = 1,
    random_to_burned_ratio = 1,
    otsu_unburned_to_burned_ratio = 1,
    # Gate 1B (2026-06-07): the engine now requires the resolved methodological
    # args (no defaults). Supply the canonical values it used to default to.
    sampling_seed              = 42,
    seed                       = 42,
    val_frac                   = 0.15,
    group_col                  = "block_id",
    impute_numeric             = "median",
    impute_factor_missing      = "MISSING",
    model_params_base          = OtsuFire:::.of_canonical_model_params(),
    feature_whitelist_override = override
  )))
  list(res = res, out_dir = out_dir)
}

# --- T1: default behaviour (no override) -------------------------------
test_that("default behaviour (no override) uses all canonical features", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  gpkg <- make_override_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  out <- run_train_smoke(gpkg)
  on.exit(unlink(out$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(out$res$files$recipe_rds)
  whitelist <- whitelist_internal_tw()

  # All canonical names that the matrix builder did not strip end up in
  # `feature_cols`. Allow the matrix builder to drop any rare-NA flags
  # that turn out to be all-zero — but the canonical set must dominate.
  expect_true(length(setdiff(whitelist, recipe$cols$feature_cols)) <= 1L,
              info = paste("missing canonical cols:",
                            paste(setdiff(whitelist, recipe$cols$feature_cols),
                                  collapse = ", ")))
  expect_null(recipe$feature_whitelist_override)
})

# --- T2: valid 47-element override -------------------------------------
test_that("valid override restricts feature_cols and x_cols", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  whitelist <- whitelist_internal_tw()
  drop_four <- c("hs_in_poly", "hs_in_buffer", "hs_min_dist_m",
                 "hs_support_present")
  override <- setdiff(whitelist, drop_four)
  expect_equal(length(override), length(whitelist) - 4L)

  gpkg <- make_override_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  out <- run_train_smoke(gpkg, override = override)
  on.exit(unlink(out$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  recipe <- readRDS(out$res$files$recipe_rds)

  # Allowed = override + their _isNA companions.
  allowed <- c(override, paste0(override, "_isNA"))

  expect_true(all(recipe$cols$feature_cols %in% allowed),
              info = paste("non-allowed feature_cols:",
                            paste(setdiff(recipe$cols$feature_cols, allowed),
                                  collapse = ", ")))
  expect_true(all(recipe$cols$x_cols %in% allowed))

  # The four dropped names must NOT appear.
  for (nm in drop_four) {
    expect_false(nm %in% recipe$cols$feature_cols)
    expect_false(nm %in% recipe$cols$x_cols)
  }

  # Persisted override.
  expect_equal(recipe$feature_whitelist_override, override)

  # meta.txt records the override.
  meta_lines <- readLines(out$res$files$meta_txt)
  expect_true(any(grepl("feature_whitelist_override_applied: TRUE",
                          meta_lines, fixed = TRUE)))
  expect_true(any(grepl(paste0("feature_whitelist_dropped_n: ",
                                 length(drop_four)),
                          meta_lines, fixed = TRUE)))
})

# --- T3: invalid override (unknown name) ------------------------------
test_that("override with unknown name errors clearly", {
  skip_if_not_installed("sf")

  gpkg <- make_override_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- train_final_internal_tw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      feature_whitelist_override = c("rbr_med", "this_is_not_canonical")
    )),
    regexp = "feature_whitelist_override.+not in the canonical"
  )
})

# --- T4: empty override --------------------------------------------------
test_that("empty override errors", {
  skip_if_not_installed("sf")

  gpkg <- make_override_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- train_final_internal_tw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      feature_whitelist_override = character(0)
    )),
    regexp = "non-empty character vector"
  )
})

# --- T5: migration error for additional_drop_cols / extra_drop_cols ---
test_that("additional_drop_cols / extra_drop_cols trigger migration error", {
  skip_if_not_installed("sf")

  gpkg <- make_override_fixture_gpkg()
  on.exit(unlink(gpkg, force = TRUE), add = TRUE)

  fn <- train_final_internal_tw()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      additional_drop_cols = c("rbr_med")
    )),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = gpkg,
      labelled_layer = "train_features",
      extra_drop_cols = c("rbr_med")
    )),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )

  # Same migration error in the OOF wrapper.
  oof_fn <- oof_pipeline_internal_tw()
  expect_error(
    suppressMessages(oof_fn(
      labelled = data.frame(),
      burned_like = data.frame(),
      labelled_df = data.frame(),
      params = list(),
      result_dir = tempdir(),
      additional_drop_cols = c("rbr_med")
    )),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )
})

# --- T_orchestrator_migration: top-level orchestrator also errors -----
test_that("run_oneyear_supervised_pipeline rejects additional_drop_cols", {
  expect_error(
    OtsuFire::run_oneyear_supervised_pipeline(
      config = structure(list(), class = "otsufire_supervised_burned_config"),
      additional_drop_cols = c("rbr_med")
    ),
    regexp = "0\\.5\\.0.+feature_whitelist_override"
  )
})
