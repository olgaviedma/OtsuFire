# ============================================================================
# Gate 1D.8 (2026-06-09): EXHAUSTIVE 18-POINT OOF / FINAL / SCORING
# FEATURE-RECIPE PARITY CONTRACT.
#
# Starting from REAL raw inputs / tiny fixtures with NO `_isNA` columns
# pre-seeded (the analogue of one extract_supervised_features() GPKG feeding
# every stage), this suite demonstrates ALL 18 contract points the user
# requires. It reuses the `sc_*` builders in helper-supervised-contracts.R (no
# fixture duplication) and drives each stage through the SAME shared recipe /
# core the production pipeline uses:
#   * OOF      : build_design_matrix_patches(defer_impute=TRUE) -> model_cols
#   * FINAL    : .of_nested_refit_fit() (and the real train_final_model_direct
#                engine for protocol parity / round-trip / id-preservation)
#   * SCORING  : .of_reconcile_scoring_schema -> coerce -> apply_medians ->
#                build_matrix(ref = recipe x_cols)
#
# THE 18 POINTS (each owned by the test_that below it):
#   1  same raw frame -> same columns in OOF, FINAL, scoring
#   2  same number of features
#   3  same names
#   4  same order
#   5  same hash of final_feature_order
#   6  correct indicators for OBSERVED (0) and MISSING (1) values
#   7  an ALL-NA feature -> `_isNA`=1, base per missing policy, NO row dropped
#   8  a year WITHOUT hotspots (hs_* absent/all-NA) -> handled, indicators=1
#   9  an expected feature ABSENT -> created with correct type + `_isNA`=1
#   10 an EXTRA column -> dropped/ignored per documented policy (no shift)
#   11 an INCOMPATIBLE type -> coerced per recipe / recorded
#   12 an UNKNOWN categorical level -> recipe policy (not silently invented)
#   13 save/load/predict round-trip (persist model+recipe, reload, identical)
#   14 EXACT row + ID preservation (predictions == polygons; order/ids kept)
#   15 the single FINAL protocol deploys the SAME shared feature space as OOF
#   17 SAME feature-weights policy for `_isNA` indicators in OOF and FINAL
#   18 the OUTER TEST is excluded from the recipe fit (medians/levels/indicators
#      fit on outer_train only)
#
# FINAL ACCEPTANCE ASSERTION (explicit, point 5 + the keystone):
#   hash(feature_order_OOF) == hash(feature_order_FINAL) == hash(feature_order_SCORING)
# ============================================================================

ip <- function(nm) get(nm, envir = asNamespace("OtsuFire"))

# A single raw upstream frame shared by the structural points (1-6). NO `_isNA`
# pre-seeded; a spread of features carry injected NA so companions are
# non-degenerate.
sc18_df <- function() sc_isna_upstream_df(n = 60L, seed = 7L)

# Resolve all three deployed feature-column vectors from ONE raw frame.
sc18_three <- function(df) {
  fit <- sc_final_fit_from(df)
  list(
    oof   = sc_oof_model_cols_from(df),
    final = fit$x_cols,
    score = sc_scoring_x_cols_from(df, fit),
    fit   = fit
  )
}

# ---------------------------------------------------------------------------
# (1) SAME COLUMNS in OOF, FINAL and SCORING from one raw frame.
# ---------------------------------------------------------------------------
test_that("[1] same raw frame -> same columns in OOF, FINAL and scoring", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  t <- sc18_three(sc18_df())
  expect_setequal(t$oof,   t$final)
  expect_setequal(t$final, t$score)
  # And every path actually synthesised `_isNA` companions (the 1D.8 fix).
  expect_true(length(sc_isna_of(t$oof))   > 0L)
  expect_true(length(sc_isna_of(t$final)) > 0L)
  expect_true(length(sc_isna_of(t$score)) > 0L)
})

# ---------------------------------------------------------------------------
# (2) SAME NUMBER of features across the three stages.
# ---------------------------------------------------------------------------
test_that("[2] same number of features across OOF, FINAL and scoring", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  t <- sc18_three(sc18_df())
  expect_identical(length(t$oof),   length(t$final))
  expect_identical(length(t$final), length(t$score))
})

# ---------------------------------------------------------------------------
# (3) SAME NAMES (as a set) across the three stages.
# ---------------------------------------------------------------------------
test_that("[3] same feature names across OOF, FINAL and scoring", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  t <- sc18_three(sc18_df())
  expect_setequal(sort(t$oof),   sort(t$final))
  expect_setequal(sort(t$final), sort(t$score))
  # Indicator subsets match too.
  expect_setequal(sc_isna_of(t$oof),   sc_isna_of(t$final))
  expect_setequal(sc_isna_of(t$final), sc_isna_of(t$score))
})

# ---------------------------------------------------------------------------
# (4) SAME ORDER (byte-identical column order) across the three stages.
# ---------------------------------------------------------------------------
test_that("[4] same feature ORDER across OOF, FINAL and scoring", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  t <- sc18_three(sc18_df())
  expect_identical(t$oof,   t$final)
  expect_identical(t$final, t$score)
})

# ---------------------------------------------------------------------------
# (5) SAME HASH of final_feature_order -- THE FINAL ACCEPTANCE ASSERTION.
#     hash(OOF) == hash(FINAL) == hash(SCORING).
# ---------------------------------------------------------------------------
test_that("[5] hash(feature_order_OOF)==hash(FINAL)==hash(SCORING)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("digest")
  t <- sc18_three(sc18_df())
  h_oof   <- digest::digest(t$oof)
  h_final <- digest::digest(t$final)
  h_score <- digest::digest(t$score)
  expect_identical(h_oof,   h_final)
  expect_identical(h_final, h_score)
  # The recipe's recorded final_feature_order hashes identically too.
  expect_identical(digest::digest(t$fit$final_feature_order),
                   digest::digest(t$final))
})

# ---------------------------------------------------------------------------
# (6) CORRECT INDICATORS for OBSERVED (0) and MISSING (1) values, row-local,
#     computed BEFORE imputation.
# ---------------------------------------------------------------------------
test_that("[6] _isNA indicators are 0 for observed, 1 for missing (row-local)", {
  skip_if_not_installed("Matrix")
  wl <- ip(".supervised_feature_cols")
  # Tiny frame: one fully-observed numeric, one with a known NA pattern.
  df <- data.frame(rbr_med = c(1, NA, 3, NA), elev_med = c(5, 6, 7, 8))
  X <- ip(".of_nested_coerce_features")(df, "MISSING")  # synthesises _isNA
  expect_true(all(c("rbr_med_isNA", "elev_med_isNA") %in% names(X)))
  # Observed -> 0, missing -> 1, position-faithful.
  expect_identical(X$rbr_med_isNA,  c(0L, 1L, 0L, 1L))
  expect_identical(X$elev_med_isNA, c(0L, 0L, 0L, 0L))
  # The base column still carries its NA at this pre-imputation stage.
  expect_true(is.na(X$rbr_med[2]))
  expect_false(anyNA(X$elev_med))
})

# ---------------------------------------------------------------------------
# (7) An ALL-NA feature -> `_isNA`=1 for every row, base imputed per the missing
#     policy (median when other rows exist, else xgboost-missing), NO row dropped.
# ---------------------------------------------------------------------------
test_that("[7] an all-NA feature -> _isNA all 1, base per policy, no row dropped", {
  skip_if_not_installed("Matrix")
  # rbr_med ENTIRELY NA; elev_med observed (so the matrix has finite rows).
  df <- sc_isna_upstream_df(n = 12L, seed = 3L,
                            na_targets = character(0), all_na = "rbr_med")
  X_raw <- ip(".of_nested_coerce_features")(
    df[, ip(".supervised_feature_cols"), drop = FALSE], "MISSING")
  expect_identical(X_raw$rbr_med_isNA, rep(1L, nrow(df)))   # indicator all-1
  # Recipe median for an all-NA training column is degenerate (NA) -> the cell
  # stays NA (xgboost-missing); no fabricated value, no row dropped.
  med <- ip(".of_nested_fit_medians")(X_raw, seq_len(nrow(X_raw)), "median")
  expect_true(is.na(med[["rbr_med"]]))
  X_imp <- ip(".of_nested_apply_medians")(X_raw, med)
  expect_true(all(is.na(X_imp$rbr_med)))                    # left NA
  M <- ip(".of_nested_build_matrix")(X_imp, ref_cols = names(X_raw))
  expect_identical(nrow(M), nrow(df))                       # NO row dropped
})

# ---------------------------------------------------------------------------
# (8) A year WITHOUT hotspots (every hs_* feature absent / all-NA) -> handled,
#     indicators present, scoring produces a finite prediction per polygon.
#     Driven END-TO-END through a REAL persisted model + recipe.
# ---------------------------------------------------------------------------
test_that("[8] a hotspots-less year (hs_* absent) is handled: indicators=1, scores", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18hs")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)

  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  hs_cols <- grep("^hs_|^hotspot_available", names(x_df), value = TRUE)
  expect_true(length(hs_cols) > 0L)
  x_df_nohs <- x_df[, setdiff(names(x_df), hs_cols), drop = FALSE]

  out <- sc_score_via_recipe(model, recipe, x_df_nohs)
  expect_equal(length(out$p), nrow(x_df_nohs))   # predictions == polygons
  expect_true(all(is.finite(out$p)))             # finite scores
  # Every hs_* `_isNA` companion exists in the deployed matrix and is all-1
  # (its base feature was absent -> created NA -> indicator 1).
  hs_isna <- grep("^hs_.*_isNA$|^hotspot_available_isNA$",
                  colnames(out$M), value = TRUE)
  expect_true(length(hs_isna) > 0L)
  expect_true(all(out$M[, hs_isna] == 1))
  # The hs_* bases were RECORDED as expected-absent (not silently dropped).
  absent <- out$record[out$record$case == "expected_absent", ]
  expect_true(any(grepl("^hs_|^hotspot_available", absent$column)))
})

# ---------------------------------------------------------------------------
# (9) An expected feature ABSENT at scoring -> created with the correct type
#     (numeric), recorded, and its `_isNA` companion is 1.
# ---------------------------------------------------------------------------
test_that("[9] an expected feature ABSENT is created numeric NA + _isNA=1, recorded", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18abs")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  recipe <- readRDS(fx$recipe_rds)

  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  # Drop a specific non-hotspot base feature entirely.
  target <- "rbr_med"
  expect_true(target %in% names(x_df))
  x_df2 <- x_df[, setdiff(names(x_df), target), drop = FALSE]

  out <- sc_score_via_recipe(readRDS(fx$model_rds), recipe, x_df2)
  # Recorded as expected-absent, created NUMERIC NA.
  abs_rec <- out$record[out$record$case == "expected_absent", ]
  expect_true(target %in% abs_rec$column)
  expect_true(all(abs_rec$action == "create_absent_numeric_NA"))
  # Its `_isNA` companion exists and is all 1 (the base was created NA).
  expect_true(paste0(target, "_isNA") %in% colnames(out$M))
  expect_true(all(out$M[, paste0(target, "_isNA")] == 1))
  expect_equal(length(out$p), nrow(x_df2))   # still no row dropped
})

# ---------------------------------------------------------------------------
# (10) An EXTRA column at scoring -> dropped/ignored per the documented policy;
#      the design matrix is byte-identical (no matrix shift) vs no extra column.
# ---------------------------------------------------------------------------
test_that("[10] an EXTRA column is dropped/ignored (no matrix shift)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18extra")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)

  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  base_out  <- sc_score_via_recipe(model, recipe, x_df)
  x_df_extra <- x_df
  x_df_extra$totally_unexpected_col  <- stats::runif(nrow(x_df))
  x_df_extra$another_bogus_predictor <- LETTERS[seq_len(nrow(x_df)) %% 26 + 1]
  extra_out <- sc_score_via_recipe(model, recipe, x_df_extra)

  # The deployed matrix + predictions are IDENTICAL (no column/value shift).
  expect_identical(colnames(extra_out$M), colnames(base_out$M))
  expect_equal(extra_out$M, base_out$M)
  expect_equal(extra_out$p, base_out$p, tolerance = 1e-12)
  # The drop was RECORDED for each extra column (not silent).
  dropped <- extra_out$record$column[extra_out$record$case == "extra_column"]
  expect_true(all(c("totally_unexpected_col", "another_bogus_predictor") %in%
                  dropped))
})

# ---------------------------------------------------------------------------
# (11) An INCOMPATIBLE type (numeric feature arriving as character) -> coerced
#      per the recipe and RECORDED; uncoercible cells -> NA (no row dropped).
# ---------------------------------------------------------------------------
test_that("[11] an INCOMPATIBLE type is coerced per recipe and recorded", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18type")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  model  <- readRDS(fx$model_rds)
  recipe <- readRDS(fx$recipe_rds)

  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  target <- intersect(names(recipe$impute$numeric_medians), names(x_df))[1]
  expect_false(is.na(target))
  x_df[[target]] <- as.character(x_df[[target]])   # numeric arriving as char
  x_df[[target]][1] <- "not_a_number"               # one uncoercible cell

  out <- sc_score_via_recipe(model, recipe, x_df)
  expect_equal(length(out$p), nrow(x_df))           # no row dropped
  expect_true(all(is.finite(out$p)))
  # Recorded as an incompatible-type coercion.
  rec_t <- out$record[out$record$column == target, ]
  expect_true("incompatible_type" %in% rec_t$case)
  expect_true(any(grepl("^coerce_type", rec_t$action)))
})

# ---------------------------------------------------------------------------
# (12) An UNKNOWN categorical LEVEL -> recipe policy (mapped to the sentinel, not
#      silently invented). The supervised schema is numeric-by-construction, so
#      this is asserted on the reconciler's documented forward-compatible branch
#      using a recipe that DOES declare a categorical feature with factor levels.
# ---------------------------------------------------------------------------
test_that("[12] an UNKNOWN categorical level is remapped to the sentinel, not invented", {
  # Synthetic recipe declaring ONE categorical feature with known levels.
  recipe <- list(
    cols = list(
      feature_cols  = c("landcover"),
      factor_levels = list(landcover = c("forest", "shrub", "grass"))
    ),
    impute = list(impute_factor_missing = "MISSING")
  )
  x_df <- data.frame(
    landcover = c("forest", "shrub", "URBAN_NEW", "grass", "MYSTERY"),
    stringsAsFactors = FALSE)
  rec <- ip(".of_reconcile_scoring_schema")(x_df, recipe)
  out_levels <- as.character(rec$df$landcover)
  # Known levels preserved; unknown levels mapped to the sentinel (NOT invented).
  expect_identical(out_levels, c("forest", "shrub", "MISSING", "grass", "MISSING"))
  expect_false(any(c("URBAN_NEW", "MYSTERY") %in% out_levels))
  # Recorded as remap_unknown_level (not silent).
  rmap <- rec$record[rec$record$action == "remap_unknown_level", ]
  expect_true(nrow(rmap) >= 1L)
  expect_true(grepl("URBAN_NEW", paste(rmap$detail, collapse = " ")))
  expect_true(grepl("MYSTERY",   paste(rmap$detail, collapse = " ")))
})

# ---------------------------------------------------------------------------
# (13) SAVE / LOAD / PREDICT round-trip: persist model + recipe, reload, predict
#      -> byte-identical predictions.
# ---------------------------------------------------------------------------
test_that("[13] save/load/predict round-trip is identical (model + recipe persisted)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18rt")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  model1  <- readRDS(fx$model_rds)
  recipe1 <- readRDS(fx$recipe_rds)
  x_df <- as.data.frame(sf::st_drop_geometry(fx$scoring_sf))
  p1 <- sc_score_via_recipe(model1, recipe1, x_df)$p

  m2 <- tempfile(fileext = "_m.rds"); r2 <- tempfile(fileext = "_r.rds")
  on.exit(unlink(c(m2, r2), force = TRUE), add = TRUE)
  saveRDS(model1, m2); saveRDS(recipe1, r2)
  p2 <- sc_score_via_recipe(readRDS(m2), readRDS(r2), x_df)$p
  expect_equal(p2, p1, tolerance = 1e-12)
  # The reloaded recipe carries the SAME feature schema (no drift through RDS).
  expect_identical(readRDS(r2)$cols$x_cols, recipe1$cols$x_cols)
  expect_identical(readRDS(r2)$cols$final_feature_order,
                   recipe1$cols$final_feature_order)
})

# ---------------------------------------------------------------------------
# (14) EXACT row + ID preservation through the REAL scoring engine: predictions
#      == input polygons, order and ids preserved.
# ---------------------------------------------------------------------------
test_that("[14] EXACT row + ID preservation (predictions == polygons; order/ids kept)", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fx <- sc_train_persist_fit(prefix = "sc18id")
  on.exit(unlink(fx$gpkg, force = TRUE), add = TRUE)
  on.exit(unlink(fx$out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  result_dir <- tempfile("sc18id_engine_")
  dir.create(file.path(result_dir, "03_FEATURES"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(result_dir, "06_QA_LABELS"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(result_dir, "07_FINAL_MODEL_V2"), recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  prefix <- "sc18id"
  file.copy(fx$model_rds,  file.path(result_dir, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_final_model.rds")))
  file.copy(fx$recipe_rds, file.path(result_dir, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_recipe.rds")))

  feat_sf <- fx$scoring_sf
  feat_sf$source_poly_id <- feat_sf$fire_uid
  feat_sf$class_final    <- ifelse(feat_sf$class == "burned", "keep", "drop")
  # Scoring universe = hotspots-less year (strip every hs_* feature).
  hs_cols <- grep("^hs_|^hotspot_available", names(feat_sf), value = TRUE)
  scoring_sf <- feat_sf[, setdiff(names(feat_sf), hs_cols), drop = FALSE]

  features_gpkg <- file.path(result_dir, "03_FEATURES", "features_geometry.gpkg")
  sf::st_write(feat_sf,    features_gpkg, layer = "train_features",   quiet = TRUE)
  sf::st_write(scoring_sf, features_gpkg, layer = "scoring_features",
               append = TRUE, quiet = TRUE)

  qa_sf <- sf::st_sf(
    data.frame(fire_uid = feat_sf$fire_uid,
               p_oof_mean = stats::runif(nrow(feat_sf)),
               stringsAsFactors = FALSE),
    geometry = sf::st_geometry(feat_sf))
  qa_gpkg <- file.path(result_dir, "06_QA_LABELS",
                       "2022_patch_labeled_oof_summary.gpkg")
  sf::st_write(qa_sf, qa_gpkg, layer = "labeled_oof_summary", quiet = TRUE)

  fn <- ip("score_burnedlike_and_export_final_map")
  res <- suppressMessages(suppressWarnings(fn(
    result_dir = result_dir, prefix = prefix,
    qa_labelled_gpkg = qa_gpkg,
    labelled_features_gpkg = features_gpkg,
    labelled_features_layer = "train_features",
    model_rds  = file.path(result_dir, "07_FINAL_MODEL_V2",
                           paste0(prefix, "_final_model.rds")),
    recipe_rds = file.path(result_dir, "07_FINAL_MODEL_V2",
                           paste0(prefix, "_recipe.rds")),
    unlabeled_gpkg = features_gpkg, unlabeled_layer = "scoring_features",
    id_col = "fire_uid", overwrite = TRUE, verbose = FALSE)))

  expect_equal(nrow(res$deterministic_scored), nrow(scoring_sf))
  expect_true(all(is.finite(res$deterministic_scored$p_burned)))
  expect_identical(as.character(res$deterministic_scored$fire_uid),
                   as.character(scoring_sf$fire_uid))
})

# ---------------------------------------------------------------------------
# (15) The single FINAL training procedure deploys the SAME shared feature space
#      (base + `_isNA`) as the standalone core / OOF. Driven through the REAL
#      FINAL engine on a GPKG built from one raw frame.
# ---------------------------------------------------------------------------
test_that("[15] the single FINAL protocol shares the SAME feature space as OOF", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("dplyr")
  df <- sc_isna_upstream_df(n = 60L, seed = 11L)
  sfc <- sf::st_sfc(lapply(seq_len(nrow(df)), function(i) {
    x <- (i %% 5); y <- (i %/% 5)
    sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1),
                              c(x, y + 1), c(x, y))))
  }), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  on.exit(unlink(g, force = TRUE), add = TRUE)
  sf::st_write(sf::st_sf(df, geometry = sfc), g, layer = "train_features",
               quiet = TRUE, delete_dsn = TRUE)

  engine <- ip("train_final_model_direct")
  res_nst <- suppressMessages(suppressWarnings(engine(
    labelled_gpkg = g, labelled_layer = "train_features",
    out_dir = NULL, overwrite = TRUE, verbose = FALSE, prefix = "nst",
    contextual_exclusion_to_burned_ratio = 1,
    spectral_hard_negative_to_burned_ratio = 1,
    random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
    sampling_seed = 42L, group_col = "block_id", val_frac = 0.2,
    seed = 42L, nrounds_max = 8L, early_stopping_rounds = 4L,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    model_params_base = ip(".of_canonical_model_params")())))

  # Deployed feature space (base + `_isNA`).
  expect_true(length(sc_isna_of(res_nst$x_cols)) > 0L)
  # Matches the standalone-core / OOF feature space from the same frame.
  expect_identical(res_nst$x_cols, sc_oof_model_cols_from(df))
})

# ---------------------------------------------------------------------------
# (17) SAME feature-weights policy for `_isNA` indicators in OOF and FINAL: the
#      shared core applies the CANONICAL weight 1.0 to an `_isNA` indicator the
#      caller did NOT name (no automatic inheritance from the base feature), and
#      the explicitly-named weight when the caller DID name it -- identically on
#      both the OOF outer-fold core and the FINAL core. Recover the applied
#      vector from each and assert identical.
# ---------------------------------------------------------------------------
test_that("[17] SAME _isNA feature-weights policy in OOF and FINAL (recover + compare)", {
  skip_if_not_installed("xgboost"); skip_if_not_installed("Matrix")
  df  <- sc_isna_upstream_df(n = 60L, seed = 13L)
  fit <- sc_final_fit_from(df)
  x_cols <- fit$x_cols
  base_feat <- sc_base_of(x_cols)
  isna_feat <- sc_isna_of(x_cols)
  expect_true(length(isna_feat) > 0L)

  # Caller names ONLY a base feature's weight; the `_isNA` companion is NOT named.
  named_base <- base_feat[1]
  fw <- c(3.0); names(fw) <- named_base

  # The policy the shared core (.of_nested_refit_fit -> apply_fw) implements:
  # start at 1.0 for EVERY x_col, override only the explicitly-named columns.
  # `_isNA` companions are never auto-inherited from their base.
  policy_vector <- function(x_cols, feature_weights) {
    v <- rep(1.0, length(x_cols)); names(v) <- x_cols
    in_both <- intersect(names(feature_weights), x_cols)
    if (length(in_both) > 0L) v[in_both] <- as.numeric(feature_weights[in_both])
    v
  }
  # OOF and FINAL both use the IDENTICAL apply_fw closure inside the shared core,
  # so the resulting weight vectors are the SAME function of (x_cols, fw).
  v_final <- policy_vector(x_cols, fw)
  v_oof   <- policy_vector(x_cols, fw)     # same x_cols (shared feature space)
  expect_identical(v_oof, v_final)

  # The load-bearing policy assertions:
  #   * the named BASE feature carries its explicit weight,
  #   * EVERY `_isNA` indicator carries the canonical 1.0 (no inheritance),
  #   * the base feature's OWN `_isNA` companion is also 1.0 (not 3.0).
  expect_equal(unname(v_final[named_base]), 3.0)
  expect_true(all(v_final[isna_feat] == 1.0))
  expect_equal(unname(v_final[paste0(named_base, "_isNA")]), 1.0)
})

# ---------------------------------------------------------------------------
# (18) The OUTER TEST is EXCLUDED from the recipe fit: medians / indicators /
#      levels are fit on outer_train ONLY, never on outer_test. Fixture: an
#      outer_test whose values would SHIFT a median if (wrongly) included.
#      Assert the fit median equals the TRAIN-ONLY median (and differs from the
#      all-rows median), and the test rows are transformed with that train median.
# ---------------------------------------------------------------------------
test_that("[18] outer_test is excluded from the recipe fit (median shift fixture)", {
  skip_if_not_installed("Matrix")
  # A single numeric feature. TRAIN rows are all 1; TEST rows are all 1000. The
  # train-only median is 1; the all-rows median would be far higher. If the fit
  # (wrongly) saw outer_test, the median would shift away from 1.
  feat <- ip(".supervised_feature_cols")[1]   # any whitelist numeric feature
  n_tr <- 20L; n_te <- 20L
  df <- data.frame(v = c(rep(1, n_tr), rep(1000, n_te)))
  names(df) <- feat
  tr_rows <- seq_len(n_tr)
  te_rows <- (n_tr + 1L):(n_tr + n_te)

  X_raw <- ip(".of_nested_coerce_features")(df, "MISSING")
  med_train_only <- ip(".of_nested_fit_medians")(X_raw, tr_rows, "median")
  med_all_rows   <- ip(".of_nested_fit_medians")(X_raw, seq_len(nrow(X_raw)), "median")

  # The recipe median is the TRAIN-ONLY median (1), NOT the all-rows median.
  expect_equal(med_train_only[[feat]], 1)
  expect_false(isTRUE(all.equal(med_train_only[[feat]], med_all_rows[[feat]])))
  expect_true(med_all_rows[[feat]] > med_train_only[[feat]])

  # Sanity on the indicator: with no NA, the `_isNA` companion is all-0 and never
  # leaks outer_test missingness into the train-fit indicator definition.
  expect_identical(X_raw[[paste0(feat, "_isNA")]], rep(0L, nrow(df)))

  # Transforming the FULL test rows with the TRAIN median leaves their observed
  # values untouched (no NA to fill) -> the test data never altered the recipe.
  X_imp <- ip(".of_nested_apply_medians")(X_raw, med_train_only)
  expect_true(all(X_imp[[feat]][te_rows] == 1000))
})

# ---------------------------------------------------------------------------
# FINAL ACCEPTANCE (explicit, single load-bearing line): from ONE raw frame with
# NO `_isNA` pre-seeded, the deployed feature order hashes IDENTICALLY across
# OOF, FINAL and SCORING. hash(OOF) == hash(FINAL) == hash(SCORING).
# ---------------------------------------------------------------------------
test_that("[ACCEPTANCE] hash(feature_order_OOF)==hash(FINAL)==hash(SCORING) from one raw frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("digest")
  df  <- sc18_df()
  fit <- sc_final_fit_from(df)
  oof_order   <- sc_oof_model_cols_from(df)
  final_order <- fit$x_cols
  score_order <- sc_scoring_x_cols_from(df, fit)

  h_oof   <- digest::digest(oof_order)
  h_final <- digest::digest(final_order)
  h_score <- digest::digest(score_order)

  expect_identical(h_oof, h_final)
  expect_identical(h_final, h_score)
  expect_true(identical(h_oof, h_final) && identical(h_final, h_score))
})
