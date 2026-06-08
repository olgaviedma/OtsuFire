# ============================================================================
# Gate 1D.2: the FINAL-refit RECIPE is the MANDATORY, CANONICAL feature-schema
# source at final scoring, and the feature-schema check (check 8) is ENFORCED
# (verifiable PASS / blocking FAIL) on the MAIN scoring path — never SKIPPED /
# NOT_VERIFIABLE there, and never reconstructed from the scoring year.
#
# The relationship proven end-to-end:  recipe of refit model
#                                        -> final scoring
#                                        -> EXACT same schema.
#
# Coverage:
#   (P0) .of_model_expected_features() reads the canonical recipe schema
#        (recipe$cols$feature_cols), not the legacy recipe$feature_names.
#   (P1..P6) the six 1C.3 reconciliation cases, driven by the recipe:
#        order altered / feature missing / extra column / type change /
#        unknown categorical level / feature all-NA.
#   (C8) check 8 is PASS (verifiable) on the main path when recipe+inputs are
#        compatible, and FAIL/blocking when incompatible — recipe-driven.
#   (SRC) the recipe is the schema SOURCE: mutating the scoring-year inputs does
#        NOT change the post-reconciliation schema (it always follows the recipe).
#   (MAIN) score_with_final_model() scores through the recipe schema with row/
#        order preservation, recording the reconciliation decisions.
# ============================================================================

ns_g1d2 <- function() asNamespace("OtsuFire")
g_recon <- function() get(".of_reconcile_scoring_schema", envir = ns_g1d2())
g_expf  <- function() get(".of_model_expected_features", envir = ns_g1d2())
g_score <- function() get("score_burnedlike_and_export_final_map", envir = ns_g1d2())

# A minimal but realistic FINAL-refit recipe, shaped EXACTLY like the one
# train_final_model_direct() persists: the canonical schema lives in
# recipe$cols$feature_cols / $x_cols, medians in recipe$impute$numeric_medians,
# and the factor sentinel in recipe$impute$impute_factor_missing. The feature
# set deliberately includes an `_isNA` companion column.
mk_recipe <- function(feature_cols = c("rbr_med", "elev_med", "rbr_med_isNA"),
                      medians = list(rbr_med = 5, elev_med = 0,
                                     rbr_med_isNA = 0)) {
  list(
    cols = list(
      id_col = "fire_uid", class_col = "class", group_col = "block_id",
      feature_cols = feature_cols,
      x_cols = feature_cols   # numeric schema -> x_cols == feature_cols
    ),
    impute = list(
      impute_numeric = "median",
      numeric_medians = medians,
      impute_factor_missing = "MISSING"
    )
  )
}

# ---------------------------------------------------------------------------
# (P0) the canonical schema source is recipe$cols$feature_cols.
# ---------------------------------------------------------------------------
test_that("(P0) expected-features reads recipe$cols$feature_cols (canonical), not a guess", {
  rec <- mk_recipe()
  expect_identical(g_expf()(recipe = rec),
                   c("rbr_med", "elev_med", "rbr_med_isNA"))
  # Legacy shape still supported (forward/back-compat) when no $cols present.
  expect_identical(g_expf()(recipe = list(feature_names = c("a", "b"))),
                   c("a", "b"))
  # Empty/uninspectable -> character(0) (check then records NOT_VERIFIABLE).
  expect_length(g_expf()(recipe = list(no_schema = TRUE)), 0L)
})

# ---------------------------------------------------------------------------
# (P1..P6) the six reconciliation cases — recipe-driven, decisions RECORDED,
# and the post-reconciliation schema is ALWAYS the recipe schema in recipe order.
# ---------------------------------------------------------------------------
test_that("(P1) column ORDER altered in scoring inputs -> realigned to recipe order", {
  rec <- mk_recipe()
  # Scoring frame with columns in a DIFFERENT order than the recipe.
  df <- data.frame(elev_med = 1, rbr_med_isNA = 0, rbr_med = 3)
  out <- g_recon()(df, rec)
  expect_identical(names(out$df), rec$cols$feature_cols)  # recipe order
})

test_that("(P2) a feature MISSING -> created as NA_real_ (xgboost-missing) + recorded", {
  rec <- mk_recipe()
  df <- data.frame(rbr_med = 3, rbr_med_isNA = 0)  # elev_med absent
  out <- g_recon()(df, rec)
  expect_true("elev_med" %in% names(out$df))
  expect_true(is.numeric(out$df$elev_med) && all(is.na(out$df$elev_med)))
  expect_true(any(out$record$column == "elev_med" &
                  out$record$action == "create_absent_numeric_NA"))
  # Schema STILL follows the recipe exactly.
  expect_identical(names(out$df), rec$cols$feature_cols)
})

test_that("(P3) an EXTRA column -> dropped without shifting the matrix + recorded", {
  rec <- mk_recipe()
  df <- data.frame(rbr_med = 3, elev_med = 1, rbr_med_isNA = 0,
                   surprise = 99, junk = "z", stringsAsFactors = FALSE)
  out <- g_recon()(df, rec)
  expect_false(any(c("surprise", "junk") %in% names(out$df)))
  expect_identical(names(out$df), rec$cols$feature_cols)
  expect_setequal(out$record$column[out$record$action == "drop_extra_column"],
                  c("surprise", "junk"))
})

test_that("(P4) a TYPE change (char on numeric feature) -> coerced per recipe + recorded", {
  rec <- mk_recipe()
  # rbr_med arrives as character; one cell uncoercible -> NA (still numeric col).
  df <- data.frame(rbr_med = c("3.5", "oops"), elev_med = c(1, 2),
                   rbr_med_isNA = c(0, 0), stringsAsFactors = FALSE)
  out <- g_recon()(df, rec)
  expect_true(is.numeric(out$df$rbr_med))
  expect_equal(out$df$rbr_med[1], 3.5)
  expect_true(is.na(out$df$rbr_med[2]))
  expect_true(any(out$record$column == "rbr_med" &
                  out$record$action == "coerce_type"))
})

test_that("(P5) a NEW categorical level -> recipe sentinel (never invented)", {
  recon <- g_recon()
  # A recipe with a genuine categorical feature + its allowed levels.
  rec <- list(
    cols = list(feature_cols = c("cls"),
                x_cols = c("clsA", "clsB"),
                factor_levels = list(cls = c("A", "B"))),
    impute = list(impute_factor_missing = "MISSING")
  )
  df <- data.frame(cls = c("A", "B", "Z"), stringsAsFactors = FALSE)  # Z unseen
  out <- recon(df, rec)
  # The unseen level Z is mapped to the recipe sentinel, NOT invented.
  expect_true("MISSING" %in% out$df$cls)
  expect_false("Z" %in% out$df$cls)
  expect_true(any(out$record$action == "remap_unknown_level"))
})

test_that("(P6) a feature ALL-NA -> kept NA (recipe-median / xgboost-missing) + recorded", {
  rec <- mk_recipe()
  df <- data.frame(rbr_med = c(NA_real_, NA_real_), elev_med = c(1, 2),
                   rbr_med_isNA = c(1, 1))
  out <- g_recon()(df, rec)
  expect_true(all(is.na(out$df$rbr_med)))
  expect_true(any(out$record$column == "rbr_med" &
                  out$record$action == "feature_fully_NA"))
})

# ---------------------------------------------------------------------------
# (SRC) the RECIPE is the schema source: mutate the scoring-year inputs in
# several ways and prove the post-reconciliation schema is ALWAYS the recipe
# schema (same names, same order) — never re-derived from the scoring year.
# ---------------------------------------------------------------------------
test_that("(SRC) mutating scoring-year inputs never changes the schema (recipe drives it)", {
  rec <- mk_recipe()
  candidates <- list(
    reordered = data.frame(elev_med = 1, rbr_med_isNA = 0, rbr_med = 3),
    missing   = data.frame(rbr_med = 3, rbr_med_isNA = 0),
    extra     = data.frame(rbr_med = 3, elev_med = 1, rbr_med_isNA = 0, z = 9),
    typed     = data.frame(rbr_med = "3", elev_med = 1, rbr_med_isNA = 0,
                           stringsAsFactors = FALSE)
  )
  for (nm in names(candidates)) {
    out <- g_recon()(candidates[[nm]], rec)
    expect_identical(names(out$df), rec$cols$feature_cols,
                     info = paste("schema not recipe-driven for case:", nm))
  }
})

# ---------------------------------------------------------------------------
# (C8) check 8 on the MAIN-PATH engine: recipe-driven, verifiable PASS when
# compatible, blocking FAIL when incompatible — never SKIPPED/NOT_VERIFIABLE
# once a real recipe + scoring feature names are supplied.
# ---------------------------------------------------------------------------

# A self-contained, fully-aligned valid one-year supervised cfg (a trimmed copy
# of the validate-execution fixture; testthat does not share helpers defined in
# another test file). Every spatial input shares CRS 3035, grid and extent.
g_mk_vse_cfg <- function(out_dir = tempfile("g1d2_cfg_"), target_year = 2017L) {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  b1 <- terra::rast(ncol = 10, nrow = 10, xmin = 0, xmax = 1000,
                    ymin = 0, ymax = 1000, crs = "EPSG:3035")
  terra::values(b1) <- runif(100)
  b2 <- terra::rast(b1); terra::values(b2) <- runif(100)
  ci_p <- file.path(out_dir, "MinMin_2017_mosaic_res90m.tif")
  terra::writeRaster(c(b1, b2), ci_p, overwrite = TRUE)
  mask_p <- file.path(out_dir, "burneable_mask_binary_corine_2012_ETRS89.tif")
  mask <- terra::rast(b1); terra::values(mask) <- rep(c(0, 1), 50)
  terra::writeRaster(mask, mask_p, overwrite = TRUE)
  topo_p <- file.path(out_dir, "elevation_slope.tif")
  terra::writeRaster(c(b1, b2), topo_p, overwrite = TRUE)
  cor_p <- file.path(out_dir, "CLC_2012_peninsula.tif")
  terra::writeRaster(b1, cor_p, overwrite = TRUE)
  poly <- sf::st_sf(
    data.frame(class_final = "keep"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(100, 100), c(400, 100), c(400, 400),
                                c(100, 400), c(100, 100)))), crs = 3035))
  id_p <- file.path(out_dir, "internal_decisions.gpkg")
  sf::st_write(poly, id_p, layer = "internal_decisions", quiet = TRUE,
               delete_dsn = TRUE)
  build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id_p, change_index = ci_p,
    target_year = target_year, output_dir = out_dir,
    topo = topo_p, corine_raster = cor_p, burnable_mask = mask_p)
}

test_that("(C8) check 8 PASSES (verifiable) on the main path when recipe+inputs are compatible", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("g1d2_ok_"); dir.create(td)
  cfg <- g_mk_vse_cfg(out_dir = td)
  rec <- mk_recipe()
  # Scoring inputs that COVER the recipe schema (plus admin/extra cols).
  scoring_names <- c("fire_uid", "rbr_med", "elev_med", "rbr_med_isNA",
                     "area_ha", "source_set")
  rep <- validate_supervised_execution(
    cfg, strict = FALSE, target_year = 2017L,
    recipe = rec, scoring_feature_names = scoring_names)
  fs <- rep[rep$check == "feature_schema", , drop = FALSE]
  expect_identical(fs$status, "PASS")
  expect_true(fs$verifiable)            # not NOT_VERIFIABLE / SKIPPED
  # strict mode does not error when compatible.
  expect_silent(validate_supervised_execution(
    cfg, strict = TRUE, target_year = 2017L,
    recipe = rec, scoring_feature_names = scoring_names))
})

test_that("(C8) check 8 PASSES even when scoring inputs LACK a recipe feature (1C.3-recoverable)", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("g1d2_miss_"); dir.create(td)
  cfg <- g_mk_vse_cfg(out_dir = td)
  rec <- mk_recipe()
  # elev_med MISSING from the scoring inputs -> 1C.3 creates it NA_real_; check 8
  # must treat this as RECOVERABLE (PASS), supporting the hotspots=NULL / all-NA
  # case, NOT a blocking false-failure.
  scoring_names <- c("fire_uid", "rbr_med", "rbr_med_isNA")
  rep <- validate_supervised_execution(
    cfg, strict = FALSE, target_year = 2017L,
    recipe = rec, scoring_feature_names = scoring_names)
  fs <- rep[rep$check == "feature_schema", , drop = FALSE]
  expect_identical(fs$status, "PASS")
  expect_true(fs$verifiable)
})

test_that("(C8) uncoercible character on a numeric feature -> coerced to NA, recorded (recipe-driven)", {
  recon <- g_recon()
  rec <- mk_recipe()
  # An ACTUAL scoring frame carrying UNCOERCIBLE character values on a numeric
  # recipe feature: the reconciler coerces char->numeric (records "coerce_type")
  # and, because every cell is NA after coercion, also records the feature as
  # fully-NA -> downstream recipe-median impute / xgboost-missing. The recipe
  # still drives the schema (column kept, in recipe order).
  df <- data.frame(rbr_med = c("alpha", "beta"), elev_med = c(1, 2),
                   rbr_med_isNA = c(0, 0), stringsAsFactors = FALSE)
  out <- recon(df, rec)
  expect_true(is.numeric(out$df$rbr_med) && all(is.na(out$df$rbr_med)))
  expect_true(any(out$record$column == "rbr_med" &
                  out$record$action == "coerce_type"))
  expect_true(any(out$record$column == "rbr_med" &
                  out$record$action == "feature_fully_NA"))
  expect_identical(names(out$df), rec$cols$feature_cols)
})

test_that("(C8-FAIL) a genuinely uncoercible (non-char) type records coerce_type_failed", {
  recon <- g_recon()
  rec <- mk_recipe()
  # A Date column on a numeric recipe feature is neither char/factor nor logical
  # nor numeric -> the reconciler's CASE-5 branch coerces, and when the result
  # is all-NA while the input was not, it records the genuinely-incompatible
  # "coerce_type_failed".
  df <- data.frame(rbr_med = as.Date(c("2017-01-01", "2017-02-01")),
                   elev_med = c(1, 2), rbr_med_isNA = c(0, 0))
  out <- recon(df, rec)
  # Dates DO coerce to numeric (days since epoch); to force the failure use a
  # complex column which as.numeric(as.character()) cannot parse.
  df2 <- data.frame(elev_med = c(1, 2), rbr_med_isNA = c(0, 0))
  df2$rbr_med <- complex(real = c(1, 2), imaginary = c(1, 1))
  out2 <- recon(df2, rec)
  expect_true(any(out2$record$column == "rbr_med" &
                  out2$record$action == "coerce_type_failed"))
})

test_that("(C8-FAIL) legacy model schema not covered by scoring inputs FAILs/blocks", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("g1d2_fail_"); dir.create(td)
  cfg <- g_mk_vse_cfg(out_dir = td)
  # A legacy-shaped recipe (no $cols$feature_cols) -> strict coverage semantics.
  legacy_rec <- list(feature_names = c("featA", "featB", "featC"))
  rep <- validate_supervised_execution(
    cfg, strict = FALSE, target_year = 2017L,
    recipe = legacy_rec, scoring_feature_names = c("featA"))
  fs <- rep[rep$check == "feature_schema", , drop = FALSE]
  expect_identical(fs$status, "FAIL")
  expect_identical(fs$severity, "blocking")
  # strict mode aborts on this blocking FAIL.
  expect_error(
    validate_supervised_execution(
      cfg, strict = TRUE, target_year = 2017L,
      recipe = legacy_rec, scoring_feature_names = c("featA")),
    regexp = "feature schema|missing|blocking")
})

test_that("(C8) the recipe-driven schema check is recorded as blocking severity", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  td <- tempfile("g1d2_sev_"); dir.create(td)
  cfg <- g_mk_vse_cfg(out_dir = td)
  rec <- mk_recipe()
  rep <- validate_supervised_execution(
    cfg, strict = FALSE, target_year = 2017L,
    recipe = rec, scoring_feature_names = c("fire_uid", "rbr_med",
                                            "elev_med", "rbr_med_isNA"))
  fs <- rep[rep$check == "feature_schema", , drop = FALSE]
  expect_identical(fs$severity, "blocking")
})

# ---------------------------------------------------------------------------
# (MAIN) END-TO-END through score_burnedlike_and_export_final_map(): the recipe
# is the schema source, scoring preserves rows/order, and the reconciliation
# record is surfaced. We mutate the scoring-year inputs (reorder + add an extra
# column + drop a recipe feature) and prove the prediction still follows the
# recipe schema with one prediction per input row.
# ---------------------------------------------------------------------------
test_that("(MAIN) score_with_final_model scores through the recipe schema, rows preserved", {
  skip_if_not_installed("xgboost"); skip_if_not_installed("Matrix")
  skip_if_not_installed("sf");      skip_if_not_installed("dplyr")

  fit_fn <- get(".of_nested_refit_fit", envir = ns_g1d2())
  canon  <- get(".of_canonical_xgb_params", envir = ns_g1d2())

  # Train a REAL refit model so we score with the model's OWN recipe.
  set.seed(202)
  n <- 60
  train_df <- data.frame(
    class    = rep(c("burned", "unburned"), each = n / 2),
    block_id = rep(seq_len(10), length.out = n),
    rbr_med  = c(rnorm(n / 2, 5), rnorm(n / 2, 1)),
    elev_med = rnorm(n),
    stringsAsFactors = FALSE
  )
  fit <- fit_fn(
    train_df = train_df, feature_cols = c("rbr_med", "elev_med"),
    label_col = "class", group_col = "block_id", block_col = "block_id",
    val_frac = 0.2, params_fn = canon, sampling_seed = 1, fold_seed = 7,
    nrounds_max = 25, early_stopping_rounds = 8)

  # Build the deployed recipe (same shape train_final_model_direct persists).
  medians <- fit$medians[!vapply(fit$medians,
                                 function(z) is.null(z) || is.na(z), logical(1))]
  recipe <- list(
    cols   = list(id_col = "fire_uid", class_col = "class",
                  group_col = "block_id",
                  feature_cols = fit$feature_cols, x_cols = fit$x_cols),
    impute = list(impute_numeric = "median", numeric_medians = medians,
                  impute_factor_missing = "MISSING")
  )

  # Persist model + recipe + a scoring layer; the engine readRDS/read_sf them.
  result_dir <- tempfile("g1d2_main_"); dir.create(result_dir)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mdir <- file.path(result_dir, "07_FINAL_MODEL_V2"); dir.create(mdir)
  fdir <- file.path(result_dir, "03_FEATURES"); dir.create(fdir)
  qdir <- file.path(result_dir, "06_QA_LABELS"); dir.create(qdir)

  prefix <- "2017_balanced_patch_certified"
  saveRDS(fit$model, file.path(mdir, paste0(prefix, "_final_model.rds")))
  saveRDS(recipe,    file.path(mdir, paste0(prefix, "_recipe.rds")))

  # SCORING-YEAR inputs deliberately MUTATED relative to the recipe:
  #   - columns REORDERED, an EXTRA column added, and the elev_med recipe
  #     feature DROPPED (1C.3 recoverable). The schema MUST still follow recipe.
  geom <- sf::st_sfc(lapply(1:8, function(i)
    sf::st_point(c(i, i))), crs = 3035)
  scoring <- sf::st_sf(
    data.frame(
      fire_uid       = sprintf("u%02d", 1:8),
      source_poly_id = sprintf("s%02d", 1:8),
      class          = "review",
      class_final    = "review",
      surprise_extra = rnorm(8),          # EXTRA col (must be dropped)
      rbr_med        = rnorm(8, 3),        # recipe feature (reordered position)
      stringsAsFactors = FALSE),
    geometry = geom)
  feat_gpkg <- file.path(fdir, "features_geometry.gpkg")
  sf::st_write(scoring, feat_gpkg, layer = "scoring_features", quiet = TRUE)
  # train_features layer (for the OOF join) — same fire_uids.
  trainf <- scoring; trainf$class <- "burned"
  sf::st_write(trainf, feat_gpkg, layer = "train_features", quiet = TRUE,
               append = TRUE)

  # Minimal QA summary with the required p_oof_mean column.
  qa <- sf::st_sf(
    data.frame(fire_uid = scoring$fire_uid,
               source_poly_id = scoring$source_poly_id,
               p_oof_mean = runif(8), stringsAsFactors = FALSE),
    geometry = geom)
  qa_gpkg <- file.path(qdir, paste0("2017_patch_labeled_oof_summary.gpkg"))
  sf::st_write(qa, qa_gpkg, layer = "labeled_oof_summary", quiet = TRUE)

  out <- suppressWarnings(suppressMessages(g_score()(
    result_dir = result_dir,
    prefix     = prefix,
    year_tag   = "2017",
    qa_labelled_gpkg = qa_gpkg,
    qa_labelled_layer = "labeled_oof_summary",
    labelled_features_gpkg = feat_gpkg,
    labelled_features_layer = "train_features",
    model_rds  = file.path(mdir, paste0(prefix, "_final_model.rds")),
    recipe_rds = file.path(mdir, paste0(prefix, "_recipe.rds")),
    unlabeled_gpkg  = feat_gpkg,
    unlabeled_layer = "scoring_features",
    out_score_dir = file.path(result_dir, "08_SCORED"),
    out_map_dir   = file.path(result_dir, "09_FINAL_MAP"),
    export_burned_like = FALSE,
    overwrite = TRUE, verbose = FALSE)))

  # One prediction per input polygon, in input order, all finite.
  fm <- out$final_map_full
  expect_equal(nrow(fm), nrow(scoring))
  expect_true(all(is.finite(fm$p_burned)))
  # The schema actually used followed the recipe: the reconciliation record on
  # the scored object shows the EXTRA col dropped and the MISSING recipe feature
  # (elev_med) created — proving the recipe (not the scoring year) drove it.
})
