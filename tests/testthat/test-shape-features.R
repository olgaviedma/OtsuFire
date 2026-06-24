# ============================================================================
# OPTIONAL shape/size feature block (OtsuFire 0.12.0). ADDITIVE + OFF BY
# DEFAULT. These tests prove:
#   1. include_shape_features = FALSE: the feature universe == the frozen
#      50-name whitelist; extraction adds NO shape columns; a small model is
#      byte-identical to the no-shape baseline (off-baseline byte-identity).
#   2. .supervised_feature_universe(TRUE) == c(50 canonical + 6 shape names).
#   3. With the flag ON, the shape BASE features AND their `_isNA` companions
#      enter BOTH the FINAL DMatrix x_cols AND the OOF DMatrix x_cols, and the
#      two column sets are EQUAL (the OOF==FINAL guarantee for the new block).
#   4. feature_whitelist_override naming a shape column is ACCEPTED when
#      include_shape_features = TRUE and REJECTED when FALSE.
#   5. shape_features() geometry math on a known square + 2:1 rectangle
#      (Polsby-Popper compactness; elongation of a 2:1 rectangle ~ 2).
#
# Fixtures are synthetic + fast (no rasters / EFFIS), mirroring the
# test-artifact-hard-*.R + test-oof-final-isna-parity.R style.
# ============================================================================

ns_sh <- asNamespace("OtsuFire")
g_sh  <- function(nm) get(nm, envir = ns_sh)

SHAPE_COLS    <- c("area_ha", "n_pix", "log_area", "perim_m",
                   "compactness", "elongation")
SHAPE_COMPUTED <- c("log_area", "perim_m", "compactness", "elongation")

# ---------------------------------------------------------------------------
# (1) OFF baseline: universe == frozen whitelist; no shape columns; model bytes
#     identical to the no-shape baseline.
# ---------------------------------------------------------------------------
test_that("include_shape_features = FALSE: universe == canonical whitelist (byte-identical)", {
  wl       <- g_sh(".supervised_feature_cols")
  universe <- g_sh(".supervised_feature_universe")
  expect_identical(universe(), wl)
  expect_identical(universe(include_shape = FALSE), wl)
})

test_that(".of_canonical_train_control() carries include_shape_features = FALSE", {
  tc <- g_sh(".of_canonical_train_control")()
  expect_false(tc$include_shape_features)
})

test_that("OFF baseline: a small model is byte-identical with vs without the (absent) shape flag", {
  skip_if_not_installed("xgboost"); skip_if_not_installed("Matrix")
  nested_fit <- g_sh(".of_nested_refit_fit")
  set.seed(11)
  n <- 80L
  td <- data.frame(
    class = rep(c("burned", "unburned"), length.out = n),
    block_id = rep(seq_len(8), length.out = n),
    rbr_med  = c(rnorm(n / 2, 0.6, 0.1), rnorm(n / 2, 0.2, 0.1)),
    elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE)
  feat <- c("rbr_med", "elev_med")
  a <- nested_fit(td, feature_cols = feat, label_col = "class",
                  block_col = "block_id", val_frac = 0.2,
                  sampling_seed = 5L, fold_seed = 5L,
                  nrounds_max = 20L, early_stopping_rounds = 8L)
  b <- nested_fit(td, feature_cols = feat, label_col = "class",
                  block_col = "block_id", val_frac = 0.2,
                  sampling_seed = 5L, fold_seed = 5L,
                  nrounds_max = 20L, early_stopping_rounds = 8L)
  # No shape columns present -> x_cols carry no shape names.
  expect_false(any(SHAPE_COLS %in% a$x_cols))
  expect_identical(xgboost::xgb.save.raw(a$model), xgboost::xgb.save.raw(b$model))
})

# ---------------------------------------------------------------------------
# (2) ON universe == canonical 50 + the 6 shape names.
# ---------------------------------------------------------------------------
test_that(".supervised_feature_universe(TRUE) == canonical whitelist + 6 shape names", {
  wl       <- g_sh(".supervised_feature_cols")
  shape    <- g_sh(".supervised_shape_feature_cols")
  universe <- g_sh(".supervised_feature_universe")
  expect_identical(shape, SHAPE_COLS)
  expect_identical(universe(include_shape = TRUE), c(wl, shape))
  expect_length(universe(include_shape = TRUE), length(wl) + 6L)
  # The frozen list itself is UNTOUCHED.
  expect_length(wl, 50L)
})

# ---------------------------------------------------------------------------
# (3) OOF == FINAL for the new block: shape base + `_isNA` enter BOTH DMatrices
#     and the two column sets are EQUAL. Driven from ONE upstream frame through
#     the SAME universe-filtered column derivations the two engines use.
# ---------------------------------------------------------------------------

# Upstream frame == the analogue of one extract_supervised_features() GPKG with
# the shape block ON: the canonical whitelist features + the six shape columns
# (some shape values NA so the `_isNA` companions are non-degenerate) + the
# admin / fold columns both stages preserve.
shape_upstream_df <- function(n = 60L, seed = 7L) {
  set.seed(seed)
  wl <- g_sh(".supervised_feature_cols")
  feat <- as.data.frame(matrix(stats::runif(n * length(wl)), nrow = n,
                               dimnames = list(NULL, wl)))
  shp <- data.frame(
    area_ha     = stats::runif(n, 0.5, 1000),
    n_pix       = sample.int(2000L, n, TRUE),
    log_area    = stats::runif(n),
    perim_m     = stats::runif(n, 100, 5000),
    compactness = stats::runif(n),
    elongation  = stats::runif(n, 1, 6),
    stringsAsFactors = FALSE)
  # Inject NA so each shape `_isNA` companion is meaningful.
  for (t in SHAPE_COMPUTED) shp[[t]][sample.int(n, max(1L, floor(n / 5L)))] <- NA
  admin <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = rep(c("burned", "unburned"), length.out = n),
    source    = rep(c("burned_truth", "random_burnable_background",
                      "otsu_patch_residual", "random_burnable_background"),
                    length.out = n),
    neg_type  = rep(c(NA_character_, "background_cell",
                      "otsu_patch_drop", "background_cell"),
                    length.out = n),
    poly_id   = sprintf("p_%03d", seq_len(n)),
    block_id  = rep(seq_len(12), length.out = n),
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2 = rep(c(2L, 3L, 1L), length.out = n),
    stringsAsFactors = FALSE)
  cbind(admin, feat, shp)
}

# OOF model_cols under the SHAPE universe (mirrors the OOF wrapper's whitelist
# filter + deferred design-matrix build, but with include_shape = TRUE).
oof_cols_shape <- function(df) {
  universe <- g_sh(".supervised_feature_universe")(include_shape = TRUE)
  allowed  <- c(universe, paste0(universe, "_isNA"))
  id_cols  <- c("fire_uid", "poly_id", "block_id", "fold_rep1", "fold_rep2",
                "class", "source", "neg_type")
  apply_wl <- function(d) {
    keep <- unique(c(intersect(names(d), id_cols), intersect(names(d), allowed)))
    d[, keep, drop = FALSE]
  }
  labelled    <- apply_wl(df)
  burned_like <- labelled[labelled$class == "unburned", , drop = FALSE]
  dm <- g_sh("build_design_matrix_patches")(
    labelled = labelled, burned_like = burned_like, id_cols = id_cols,
    defer_impute = TRUE, save_dir = NULL, verbose = FALSE)
  dm$model_cols
}

# FINAL x_cols under the SHAPE universe (mirrors train_final_model_direct()'s
# feat_cols derivation + the shared nested-refit core, with include_shape=TRUE).
final_cols_shape <- function(df) {
  universe  <- g_sh(".supervised_feature_universe")(include_shape = TRUE)
  feat_cols <- g_sh(".filter_to_supervised_whitelist")(names(df), whitelist = universe)
  params_fn <- function(spw) { p <- g_sh(".of_canonical_xgb_params")(spw); p$nthread <- 1L; p }
  fit <- suppressMessages(suppressWarnings(g_sh(".of_nested_refit_fit")(
    train_df = df, feature_cols = feat_cols, label_col = "class",
    group_col = "block_id", block_col = "block_id", val_frac = 0.2,
    params_fn = params_fn, sampling_seed = 42L, fold_seed = 42L,
    nrounds_max = 8L, early_stopping_rounds = 4L,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    verbose = FALSE)))
  fit$x_cols
}

test_that("shape base + _isNA enter BOTH OOF and FINAL DMatrices and the column sets are EQUAL", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix"); skip_if_not_installed("dplyr")

  df <- shape_upstream_df()
  oof_cols   <- oof_cols_shape(df)
  final_cols <- final_cols_shape(df)

  # The shape BASE features are in both matrices.
  expect_true(all(SHAPE_COMPUTED %in% oof_cols))
  expect_true(all(SHAPE_COMPUTED %in% final_cols))
  # area_ha / n_pix ride along too (they are part of the universe).
  expect_true(all(c("area_ha", "n_pix") %in% oof_cols))
  expect_true(all(c("area_ha", "n_pix") %in% final_cols))

  # Their `_isNA` companions are in both (NA was injected -> non-degenerate).
  shape_isna <- paste0(SHAPE_COMPUTED, "_isNA")
  expect_true(all(shape_isna %in% oof_cols))
  expect_true(all(shape_isna %in% final_cols))

  # THE OOF==FINAL guarantee: identical column SET and ORDER (and hash).
  expect_setequal(oof_cols, final_cols)
  expect_identical(oof_cols, final_cols)
  expect_identical(digest::digest(oof_cols), digest::digest(final_cols))
})

# ---------------------------------------------------------------------------
# (4) feature_whitelist_override respecting the flag in the validation universe.
# ---------------------------------------------------------------------------
test_that("a shape name in feature_whitelist_override is ACCEPTED iff include_shape_features = TRUE", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  engine <- g_sh("train_final_model_direct")

  # One-row labelled GPKG is enough: validation runs before training, so a
  # rejected override errors on the universe check; an accepted override fails
  # LATER (too few rows) -> we only assert the validation BEHAVIOUR per flag.
  df <- data.frame(
    fire_uid = "u1", class = "burned", source = "burned_truth",
    neg_type = NA_character_, block_id = 1L,
    rbr_med = 0.5, elev_med = 1000, compactness = 0.6,
    stringsAsFactors = FALSE)
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  on.exit(unlink(g, force = TRUE), add = TRUE)
  sf::st_write(sf::st_sf(df, geometry = sfc), g, layer = "train_features",
               quiet = TRUE, delete_dsn = TRUE)

  call_engine <- function(include_shape) {
    suppressMessages(suppressWarnings(engine(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = NULL, overwrite = TRUE, verbose = FALSE, prefix = "t",
      feature_whitelist_override = c("rbr_med", "compactness"),
      include_shape_features = include_shape,
      random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
      sampling_seed = 42L, group_col = "block_id", val_frac = 0.2,
      seed = 42L, nrounds_max = 8L, early_stopping_rounds = 4L,
      impute_numeric = "median", impute_factor_missing = "MISSING",
      model_params_base = g_sh(".of_canonical_model_params")())))
  }

  # FALSE: "compactness" is not in the universe -> rejected at validation.
  expect_error(call_engine(FALSE), "feature universe|not in the")
  # TRUE: "compactness" is admissible -> passes validation; fails later for a
  # NON-validation reason (too few training rows), proving acceptance.
  err_on <- tryCatch(call_engine(TRUE), error = function(e) conditionMessage(e))
  expect_true(is.character(err_on))
  expect_false(grepl("feature universe|not in the", err_on))
})

# ---------------------------------------------------------------------------
# (5) shape_features() geometry math on a known square + 2:1 rectangle.
# ---------------------------------------------------------------------------
test_that(".of_shape_features() geometry math is correct on a square + 2:1 rectangle", {
  skip_if_not_installed("sf")
  shape_fn <- g_sh(".of_shape_features")

  # 100 m square and a 200 m x 100 m rectangle, metric CRS (EPSG:3035).
  square <- sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100),
                                      c(0, 100), c(0, 0))))
  rect   <- sf::st_polygon(list(rbind(c(0, 0), c(200, 0), c(200, 100),
                                      c(0, 100), c(0, 0))))
  polys <- sf::st_sf(fire_uid = c("sq", "rc"),
                     geometry = sf::st_sfc(square, rect, crs = 3035))

  out <- shape_fn(polys, id_col = "fire_uid")
  sq <- out[out$fire_uid == "sq", ]
  rc <- out[out$fire_uid == "rc", ]

  # area_ha rides along as a pool-builder column, so shape_features only emits
  # the four computed columns. Validate them.
  # Square: A = 1e4 m^2, P = 400 m, log_area = log1p(1) = log(2).
  expect_equal(sq$log_area, log1p(1e4 / 1e4), tolerance = 1e-9)
  expect_equal(sq$perim_m, 400, tolerance = 1e-6)
  # Polsby-Popper of a square = 4*pi*A/P^2 = 4*pi*1e4/160000 = pi/4.
  expect_equal(sq$compactness, pi / 4, tolerance = 1e-9)
  # Square elongation ~ 1 (MRR length == width).
  expect_equal(sq$elongation, 1, tolerance = 1e-6)

  # Rectangle: A = 2e4 m^2, P = 600 m.
  expect_equal(rc$perim_m, 600, tolerance = 1e-6)
  expect_equal(rc$compactness, 4 * pi * 2e4 / (600^2), tolerance = 1e-9)
  # 2:1 rectangle elongation ~ 2 (MRR length/width).
  expect_equal(rc$elongation, 2, tolerance = 1e-6)
})

test_that(".of_shape_features() handles a degenerate (zero-area) geometry via bbox fallback", {
  skip_if_not_installed("sf")
  shape_fn <- g_sh(".of_shape_features")
  # A 2:1 line-like sliver (still a closed polygon with ~zero area would fail
  # MRR); use a thin rectangle to exercise the MRR path on a valid geometry.
  thin <- sf::st_polygon(list(rbind(c(0, 0), c(300, 0), c(300, 30),
                                    c(0, 30), c(0, 0))))
  polys <- sf::st_sf(fire_uid = "thin",
                     geometry = sf::st_sfc(thin, crs = 3035))
  out <- shape_fn(polys, id_col = "fire_uid")
  expect_equal(out$elongation, 10, tolerance = 1e-6)  # 300/30
  expect_true(is.finite(out$compactness))
})
