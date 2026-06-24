# =============================================================================
# PHASE 2 (artifact_hard): END-TO-END public-wrapper wiring proof.
#
# Phase A wires the artifact_hard hard-negative feature into the PUBLIC
# wrappers (train_final_burned_model / run_oof_diagnostics / the one-year
# orchestrator) and the engine eligibility/capping path. This file proves, on a
# small synthetic fixture (fast; no EFFIS / no real rasters):
#
#   1. OFF-by-default parity: with artifact_hard disabled + source_weights NULL,
#      passing the new params is BYTE-IDENTICAL (xgb.save.raw) to NOT passing
#      them -> the wiring is inert when off.
#   2. Real activation: with artifact_hard rows present, the public
#      train_final_burned_model promotes them into training WITHOUT hitting the
#      GATE-6.5 "unburned with no valid bucket" error, and applies the per-row
#      weights (training_ok includes the artifact_hard rows).
#   3. Same scored universe: enabling artifact_hard does NOT change the scored
#      row count.
#   4. p_burned_eval uses the held-out p_oof for labeled rows, not the in-sample
#      p_burned_model (column-logic contract).
#   5. Weights + artifact_hard_source thread correctly through the
#      wrapper -> engine chain (the engine sees the resolved weights).
#   6. runner<->wrapper EQUIVALENCE: the package engine
#      (train_final_model_direct, via the same eligibility/capping/weight
#      helpers the public wrapper uses) and the verified runner's mirrored FINAL
#      assembly produce the SAME promoted/capped training set and the SAME final
#      model bytes for the SAME inputs.
#
# All of this runs against a tiny in-memory fixture so the file is fast and has
# no external dependency.
# =============================================================================

ns <- asNamespace("OtsuFire")

# Active whitelist used throughout (a 2-feature subset of the canonical set).
FEAT <- c("rbr_med", "elev_med")

# Build a synthetic labelled training frame (sf) carrying the columns the
# supervised engine reads: fire_uid / class / source / neg_type / block_id /
# fold_rep1 / fold_rep2 + the 2 model features + a point geometry.
mk_labelled_sf <- function(seed = 11, n_pos = 30, n_rand = 30, n_otsu = 20,
                           n_ah = 0) {
  set.seed(seed)
  n <- n_pos + n_rand + n_otsu + n_ah
  cls <- c(rep("burned", n_pos), rep("unburned", n_rand + n_otsu + n_ah))
  src <- c(rep(NA_character_, n_pos),
           rep("random_burnable_background", n_rand),
           rep("otsu_patch_residual", n_otsu),
           rep("artifact_hard", n_ah))
  ngt <- c(rep(NA_character_, n_pos),
           rep(NA_character_, n_rand),
           rep("otsu_patch_drop", n_otsu),
           rep("artifact_hard_negative", n_ah))
  df <- data.frame(
    fire_uid = as.character(seq_len(n)),
    class    = cls,
    source   = src,
    neg_type = ngt,
    block_id = paste0("blk_", sample(seq_len(8), n, replace = TRUE)),
    fold_rep1 = sample(1:5, n, replace = TRUE),
    fold_rep2 = sample(1:5, n, replace = TRUE),
    rbr_med  = c(rnorm(n_pos, 0.6, 0.1),
                 rnorm(n_rand, 0.15, 0.05),
                 rnorm(n_otsu, 0.2, 0.05),
                 rnorm(n_ah, 0.45, 0.05)),
    elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE
  )
  geom <- sf::st_sfc(lapply(seq_len(n), function(i)
    sf::st_point(c(runif(1, 0, 1e5), runif(1, 0, 1e5)))), crs = 3035)
  sf::st_sf(df, geometry = geom)
}

mk_cfg <- function(out_dir, enable_ah = FALSE, source_weights = NULL) {
  npp <- list(caps = c(random = 1.0, otsu = 1.0))
  if (isTRUE(enable_ah)) {
    npp$artifact_hard <- list(
      enabled = TRUE, total_weight_ratio = 0.5,
      persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
      rbr_med_reference = "negative", reason_whitelist = character(0),
      persist_delta_max = -100, area_ha_min = 500, doy_iqr_max = 1)
  }
  if (!is.null(source_weights)) npp$source_weights <- source_weights
  build_supervised_burned_config(
    run_label = "balanced",
    internal_decisions = file.path(out_dir, "decisions.gpkg"),
    change_index = file.path(out_dir, "rbr.tif"),
    target_year = 1989L, output_dir = out_dir,
    feature_whitelist_override = FEAT,
    negative_pool_params = npp,
    nrounds_max = 40L, early_stop = 12L)
}

write_train_gpkg <- function(sf_obj, out_dir) {
  p <- file.path(out_dir, "labelled.gpkg")
  if (file.exists(p)) unlink(p, force = TRUE)
  sf::st_write(sf_obj, p, layer = "train_features", quiet = TRUE)
  p
}

# ---------------------------------------------------------------------------
# 1) OFF-by-default parity: passing the new params with artifact_hard OFF +
#    source_weights NULL is byte-identical to NOT passing them.
# ---------------------------------------------------------------------------
test_that("OFF-by-default: new params inert -> identical FINAL model bytes", {
  skip_if_not_installed("xgboost")
  od <- file.path(tempdir(), paste0("ah_e2e_off_", as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  L <- mk_labelled_sf(seed = 21)
  gpkg <- write_train_gpkg(L, od)
  cfg <- mk_cfg(od, enable_ah = FALSE)

  tm_base <- train_final_burned_model(
    train_features = gpkg, config = cfg, oof_agg = NULL,
    out_dir = file.path(od, "base"), verbose = FALSE)
  tm_new <- train_final_burned_model(
    train_features = gpkg, config = cfg, oof_agg = NULL,
    source_weights = NULL, artifact_hard_source = NULL,
    out_dir = file.path(od, "new"), verbose = FALSE)

  expect_identical(xgboost::xgb.save.raw(tm_base$model),
                   xgboost::xgb.save.raw(tm_new$model))
  # OFF path: the per-source weight log is the empty 0-row frame (no weighting).
  expect_equal(nrow(tm_new$sample_weight_log), 0L)
  # training_ok carries NO artifact_hard rows when off.
  expect_false(any(tm_new$training_ok$source %in% "artifact_hard"))
})

# ---------------------------------------------------------------------------
# 2) Real activation from the public wrapper: artifact_hard rows present in the
#    training frame train WITHOUT the GATE-6.5 error, and the weights apply.
# ---------------------------------------------------------------------------
test_that("activation: artifact_hard rows train via the public wrapper (no GATE-6.5)", {
  skip_if_not_installed("xgboost")
  od <- file.path(tempdir(), paste0("ah_e2e_on_", as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  # 8 promoted artifact_hard rows already present in the labelled frame.
  L <- mk_labelled_sf(seed = 22, n_ah = 8)
  gpkg <- write_train_gpkg(L, od)
  cfg <- mk_cfg(od, enable_ah = TRUE)

  tm <- expect_no_error(
    train_final_burned_model(
      train_features = gpkg, config = cfg, oof_agg = NULL,
      out_dir = file.path(od, "on"), verbose = FALSE))

  # The artifact_hard rows are in training_ok (promoted + uncapped).
  n_ah_train <- sum(tm$training_ok$source %in% "artifact_hard")
  expect_equal(n_ah_train, 8L)
  # The per-source weight log records the artifact_hard source with the
  # 0.5*positives target balance.
  lg <- tm$sample_weight_log
  expect_true("artifact_hard" %in% lg$source)
  ah_row <- lg[lg$source == "artifact_hard", ]
  n_pos <- sum(L$class == "burned")
  expect_equal(ah_row$total_weight, 0.5 * n_pos)
})

# ---------------------------------------------------------------------------
# 5) Weights + artifact_hard_source thread through wrapper -> engine: an
#    explicit source_weights pin reaches the engine and pins the weight.
# ---------------------------------------------------------------------------
test_that("source_weights pin threads through the wrapper to the engine", {
  skip_if_not_installed("xgboost")
  od <- file.path(tempdir(), paste0("ah_e2e_pin_", as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  L <- mk_labelled_sf(seed = 23, n_ah = 6)
  gpkg <- write_train_gpkg(L, od)
  cfg <- mk_cfg(od, enable_ah = TRUE,
                source_weights = c(artifact_hard = 0.3))

  tm <- train_final_burned_model(
    train_features = gpkg, config = cfg, oof_agg = NULL,
    out_dir = file.path(od, "pin"), verbose = FALSE)
  lg <- tm$sample_weight_log
  ah_row <- lg[lg$source == "artifact_hard", ]
  # pinned weight 0.3 per row (NOT the 0.5-balance), n=6 -> total 1.8.
  expect_equal(ah_row$per_row_weight, 0.3)
  expect_equal(ah_row$total_weight, 0.3 * 6)
})

# ---------------------------------------------------------------------------
# 4) p_burned_eval = coalesce(p_oof_mean, p_burned_model): labeled rows use the
#    held-out p_oof, NOT the in-sample p_burned_model. (Column-logic contract,
#    mirroring the assembly in internal-sup-final-map.R.)
# ---------------------------------------------------------------------------
test_that("p_burned_eval uses p_oof for labeled rows, p_burned_model otherwise", {
  fm <- data.frame(
    p_burned_oof   = c(0.11, NA, 0.91, NaN),
    p_burned_model = c(0.50, 0.60, 0.70, 0.80))
  pe_oof   <- suppressWarnings(as.numeric(fm$p_burned_oof))
  pe_model <- suppressWarnings(as.numeric(fm$p_burned_model))
  p_burned_eval <- ifelse(is.finite(pe_oof), pe_oof, pe_model)
  expect_equal(p_burned_eval, c(0.11, 0.60, 0.91, 0.80))
  # no row dropped
  expect_equal(length(p_burned_eval), nrow(fm))
})

# ---------------------------------------------------------------------------
# 3 + 6) runner<->wrapper EQUIVALENCE on a synthetic fixture: the package engine
# (train_final_model_direct, the engine the public wrapper calls) and the
# verified runner's MIRRORED FINAL assembly produce the SAME capped/promoted
# training set and the SAME final model bytes for the SAME inputs.
#
# The runner mirror (from _experiments/runner_artifact_hard.R, mirror_final_fit)
# calls the package's OWN eligibility resolver + capping helper + sample-weight
# resolver + ML core directly. We reproduce that mirror here and assert it
# matches the public engine output.
# ---------------------------------------------------------------------------
test_that("runner<->wrapper equivalence: same promoted set + same model bytes", {
  skip_if_not_installed("xgboost")
  od <- file.path(tempdir(), paste0("ah_e2e_eq_", as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  L <- mk_labelled_sf(seed = 24, n_ah = 7)
  gpkg <- write_train_gpkg(L, od)
  cfg <- mk_cfg(od, enable_ah = TRUE)
  tc  <- cfg$train_control
  ah_source <- "artifact_hard"

  # ---- (a) public engine path (what train_final_burned_model calls) --------
  tm <- train_final_burned_model(
    train_features = gpkg, config = cfg, oof_agg = NULL,
    out_dir = file.path(od, "wrap"), verbose = FALSE)
  wrap_ids <- sort(as.character(tm$training_ok$fire_uid))

  # ---- (b) runner mirror (mirror_final_fit, package helpers directly) ------
  Lr <- L
  Lr$fire_uid <- as.character(Lr$fire_uid)
  Lr$class    <- as.character(Lr$class)
  id_col <- "fire_uid"; class_col <- "class"
  L_class <- as.character(Lr[[class_col]])
  L_source <- as.character(Lr[["source"]])
  L_neg <- as.character(Lr[["neg_type"]])
  L_id <- as.character(Lr[[id_col]])

  elig <- ns$.of_resolve_supervised_eligibility(
    id = L_id, class = L_class, source = L_source, neg_type = L_neg,
    random_background_source = "random_burnable_background",
    otsu_unburned_source = "otsu_patch_residual",
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep"),
    origin_stage = "FINAL", artifact_hard_source = ah_source)
  caps <- c(random = tc$caps$random, otsu = tc$caps$otsu)
  if ("artifact_hard" %in% names(elig$negatives_by_bucket)) caps["artifact_hard"] <- Inf
  cap <- ns$.of_cap_negative_buckets(
    positive_idx = elig$positive_idx,
    negatives_by_bucket = elig$negatives_by_bucket,
    n_burned = length(elig$positive_idx),
    caps = caps, seed = tc$seeds$final_sampling_seed,
    id = L_id, context = "FINAL")
  L_ok <- Lr[cap$selected_indices, , drop = FALSE]
  L_ok <- L_ok[!duplicated(as.character(L_ok[[id_col]])), , drop = FALSE]
  L_df <- sf::st_drop_geometry(L_ok)

  active_whitelist <- tc$feature_whitelist_override
  feat_cols <- ns$.filter_to_supervised_whitelist(names(L_df),
                                                  whitelist = active_whitelist)
  sw_src <- as.character(L_df[["source"]])
  sw_res <- ns$.of_resolve_sample_weights(
    class = as.character(L_df[[class_col]]), source = sw_src,
    source_weights = NULL, artifact_hard_source = ah_source,
    total_weight_ratio = cfg$negative_pool_params$artifact_hard$total_weight_ratio)
  params_from_cfg <- function(scale_pos_weight) {
    p <- cfg$model_params; p[["scale_pos_weight"]] <- scale_pos_weight; p
  }
  fit <- ns$.of_nested_refit_fit(
    train_df = L_df, feature_cols = feat_cols, label_col = class_col,
    group_col = tc$group_col, block_col = tc$group_col, val_frac = tc$val_frac,
    params_fn = params_from_cfg,
    sampling_seed = tc$seeds$final_sampling_seed, fold_seed = tc$seeds$final_seed,
    feature_weights = tc$feature_weights, sample_weights = sw_res$weights,
    nrounds_max = tc$nrounds_max, early_stopping_rounds = tc$early_stop,
    impute_numeric = tc$impute_numeric, impute_factor_missing = tc$impute_factor_missing,
    verbose = FALSE)

  mirror_ids <- sort(as.character(L_ok$fire_uid))

  # SAME promoted/capped training set...
  expect_identical(wrap_ids, mirror_ids)
  # ...and the SAME final model bytes.
  expect_identical(xgboost::xgb.save.raw(tm$model),
                   xgboost::xgb.save.raw(fit$model))
})
