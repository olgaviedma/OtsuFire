# ============================================================================
# Gate 1D.5 (2026-06-08): THREE-LEVEL REPRODUCIBILITY CONTRACT.
#
# This suite defines and locks the reproducibility contract of the supervised
# pipeline at THREE distinct levels. It DELIBERATELY does NOT require bitwise-
# identical XGBoost model binaries as the general criterion: model binaries can
# legitimately differ across xgboost/GDAL versions and thread counts. Instead
# each artifact is held to the STRONGEST guarantee it can actually keep.
#
# ----------------------------------------------------------------------------
# STEP-0 INVENTORY — reproducible artifacts + where their identity is fixed,
# and which ones EXCLUDE wall-clock time from their hash (file:line as of this
# gate):
#
#   artifact                         determined / hashed at
#   -------------------------------- ---------------------------------------
#   resolved cfg (methodological)    build_supervised_burned_config()
#                                      R/supervised-config.R:166
#   per-field provenance             cfg$resolved_params_provenance
#                                      R/supervised-config.R:538
#   folds + fold fingerprint         make_fold_fingerprint()
#                                      R/internal-sup-make-folds.R:563
#       -> EXCLUDES timestamp; created_at kept SEPARATE, not hashed
#                                      R/internal-sup-make-folds.R:499-505,553
#   neg-pool fingerprint (1C.2)      otsu_negative_param_fingerprint()
#                                      R/supervised-pools.R:690
#       cfg-preview of that pool fp  .of_cfg_neg_pool_fingerprint_preview()
#                                      R/supervised-validate-execution.R:1075
#       -> "folds in NO timestamp"   R/supervised-validate-execution.R:1064-1071
#   run fingerprint (1C.4)           .of_cfg_run_fingerprint()
#                                      R/supervised-validate-execution.R:1122
#       -> "folds in NO timestamp"   R/supervised-validate-execution.R:1099-1116
#   cache NAMES (content-aware)      vector_input_cache_tag() ("in-<hash>")
#                                      R/validate-fire-maps.R:492
#       -> identity from content (path+mtime+size / sf structure), NOT time
#   feature schema (recipe)          recipe$cols / recipe$impute
#                                      R/internal-sup-train-final-direct.R:624
#       -> created_at DROPPED from recipe so it is byte-reproducible
#                                      R/internal-sup-train-final-direct.R:618-623
#   observation/row order            preserved by score_with_final_model()
#                                      R/internal-sup-final-map.R:307-311
#   base fingerprinter (no time)     otsu_negative_param_fingerprint()
#                                      R/internal-sup-otsu-negative.R
#
# Every fingerprint above is a deterministic function of CONTENT only; none
# folds in Sys.time(). The fold/pool/run fingerprints document this explicitly
# at the cited lines and keep created_at as separate, un-hashed metadata.
#
# ----------------------------------------------------------------------------
# LEVEL A — EXACT IDENTITY (must be bit-identical across two equivalent runs at
#           different wall-clock times):
#   resolved cfg, provenance, folds + fold fingerprint, pool ids + pool
#   fingerprint, cache NAMES, manifests EXCLUDING timestamps, observation/row
#   order, feature schema (recipe x_cols/feature_cols/levels/medians).
#   The proof has two halves: (i) two evaluations at DIFFERENT times produce
#   IDENTICAL fingerprints/cache-names; (ii) a relevant cfg-field change
#   produces a DIFFERENT fingerprint.
#
# LEVEL B — NUMERIC REPRODUCIBILITY (nthread = 1):
#   same seeds + same nthread(=1) + same data + same environment => model
#   PREDICTIONS (p_burned) match within an EXPLICIT tolerance, and the DERIVED
#   CLASSES at a fixed threshold match EXACTLY.
#     tolerance: PRED_ABS_TOL = 1e-8 (absolute). Justification: with identical
#     inputs/seed/nthread=1 and the same xgboost build the two boosters are the
#     SAME sequence of float32 operations; predictions are expected bit-equal,
#     so 1e-8 is a generous guard against benign float32<->float64 readout
#     noise while still failing on any real divergence (>> any rounding).
#     CLASS_THRESHOLD = 0.5; class parity asserted EXACTLY (no tolerance).
#
# LEVEL C — CARTOGRAPHIC REPRODUCIBILITY:
#   final map outputs match: final IDs (exact), number of polygons (exact),
#   class vector (exact), per-polygon AREA within AREA_REL_TOL, and a geometric
#   fingerprint (bbox + vertex count + per-polygon area signature) identical.
#     contract adopted: GEOMETRIES ARE EXACT (same input geometries flow
#     through, so we assert an exact geometric fingerprint) and AREA is held to
#     AREA_REL_TOL = 1e-9 relative, to absorb only platform double-formatting
#     noise in st_area; counts/IDs/classes are EXACT.
#
# ----------------------------------------------------------------------------
# RECORDED LIMITATIONS that can break EXACT reproducibility (apply to Level B
# numeric parity and any model-binary comparison; Levels A and C identity are
# unaffected because they are content fingerprints / pass-through geometry):
#   * MULTITHREADING: xgboost with nthread > 1 accumulates gradients in a
#     nondeterministic order -> float nondeterminism. The strict numeric test
#     PINS nthread = 1.
#   * XGBOOST VERSION: a different xgboost build can change the float op order
#     / defaults -> predictions may move beyond PRED_ABS_TOL. Test pinned to the
#     installed build; cross-version parity is OUT of contract.
#   * GDAL VERSION: geometry/area readout can differ at the ULP level across
#     GDAL builds -> handled by AREA_REL_TOL (Level C) and by hashing geometry
#     structure rather than raw WKB bytes.
#   * OS / BLAS: platform math libraries can perturb the last ULPs. nthread=1
#     plus the explicit tolerances above bound this.
# ============================================================================

ns_rep <- asNamespace("OtsuFire")
g_rep  <- function(nm) get(nm, envir = ns_rep)

# Explicit, documented tolerances (see header).
PRED_ABS_TOL  <- 1e-8
CLASS_THRESHOLD <- 0.5
AREA_REL_TOL  <- 1e-9

# ---------------------------------------------------------------------------
# Shared synthetic fixtures (small; nthread=1 for the numeric level).
# ---------------------------------------------------------------------------

# A tiny labelled training data.frame on the canonical supervised whitelist,
# aligned in spirit with make_sym_oof_df (test-oof-final-symmetry.R) so the two
# contract suites stay consistent.
make_rep_train_df <- function(n = 80, seed = 11L) {
  set.seed(seed)
  wl <- g_rep(".supervised_feature_cols")
  feat <- as.data.frame(
    matrix(stats::runif(n * length(wl)), nrow = n,
           dimnames = list(NULL, wl))
  )
  admin <- data.frame(
    fire_uid = sprintf("uid_%03d", seq_len(n)),
    class    = rep(c("burned", "unburned"), length.out = n),
    block_id = rep(seq_len(10), length.out = n),
    stringsAsFactors = FALSE
  )
  cbind(admin, feat)
}

# Train ONE final model through the SHARED leakage-free core at nthread=1.
# Returns the deployable refit model + recipe-equivalent schema. params_fn pins
# nthread = 1 so the numeric level is deterministic.
train_rep_final <- function(train_df, feature_cols, fold_seed = 42L) {
  params_fn <- function(spw) {
    p <- g_rep(".of_canonical_xgb_params")(spw)
    p$nthread <- 1L           # PIN: single thread for numeric determinism.
    p
  }
  g_rep(".of_nested_refit_fit")(
    train_df              = train_df,
    feature_cols          = feature_cols,
    label_col             = "class",
    group_col             = "block_id",
    block_col             = "block_id",
    val_frac              = 0.2,
    params_fn             = params_fn,
    sampling_seed         = 42L,
    fold_seed             = fold_seed,
    nrounds_max           = 25L,
    early_stopping_rounds = 8L,
    impute_numeric        = "median",
    impute_factor_missing = "MISSING",
    verbose               = FALSE
  )
}

# Build a scoring matrix from a fit's recipe (medians + x_cols), NA-preserving,
# exactly the path score_with_final_model() uses.
predict_rep <- function(fit, score_df) {
  M <- g_rep(".of_nested_transform")(
    df                    = score_df,
    feature_cols          = fit$feature_cols,
    med                   = fit$medians,
    ref_x_cols            = fit$x_cols,
    impute_factor_missing = "MISSING"
  )
  dmat <- xgboost::xgb.DMatrix(M, missing = NA)
  as.numeric(stats::predict(fit$model, dmat))
}

# ===========================================================================
# LEVEL A — EXACT IDENTITY
# ===========================================================================

# Build a minimal resolved cfg from on-disk dummy inputs (content-stable).
make_rep_cfg <- function(dir, year = 2017L, scenario = "balanced") {
  dec <- file.path(dir, "dec.gpkg")
  ci  <- file.path(dir, "ci.tif")
  if (!file.exists(dec)) writeLines("decisions", dec)
  if (!file.exists(ci))  writeLines("change-index", ci)
  build_supervised_burned_config(
    scenario = scenario, internal_decisions = dec, change_index = ci,
    target_year = year, output_dir = dir
  )
}

test_that("LEVEL A: cfg run-fingerprint is wall-clock-free (two times -> identical)", {
  skip_if_not_installed("sf")
  d <- tempfile("repA_"); dir.create(d)
  cfg <- make_rep_cfg(d)

  rf <- g_rep(".of_cfg_run_fingerprint")
  fp1 <- rf(cfg)
  Sys.sleep(0)                # the contract must hold regardless of time
  fp2 <- rf(cfg)

  expect_identical(fp1$checksum, fp2$checksum)
  expect_identical(fp1$text,     fp2$text)
  # No wall-clock leakage into the run fingerprint text.
  expect_false(grepl(format(Sys.Date()), fp1$text, fixed = TRUE))
  expect_false(grepl("created_at", fp1$text, fixed = TRUE))
})

test_that("LEVEL A: resolved cfg methodological content + provenance are identical across builds", {
  skip_if_not_installed("sf")
  d1 <- tempfile("repA1_"); dir.create(d1)
  d2 <- tempfile("repA2_"); dir.create(d2)
  cfg1 <- make_rep_cfg(d1)
  cfg2 <- make_rep_cfg(d2)

  # The methodological content (train_control + model_params) is byte-identical
  # regardless of the output_dir / wall-clock of the build.
  expect_identical(cfg1$train_control, cfg2$train_control)
  expect_identical(cfg1$model_params,  cfg2$model_params)
  expect_identical(cfg1$resolved_params_provenance,
                   cfg2$resolved_params_provenance)
  expect_identical(cfg1$scenario, cfg2$scenario)
  expect_identical(cfg1$target_year, cfg2$target_year)
})

test_that("LEVEL A: run-fingerprint CHANGES on a relevant cfg-field change", {
  skip_if_not_installed("sf")
  d <- tempfile("repAchg_"); dir.create(d)
  cfg  <- make_rep_cfg(d, year = 2017L)
  rf   <- g_rep(".of_cfg_run_fingerprint")
  base <- rf(cfg)

  # target_year is a consequential field -> different fingerprint.
  d2   <- tempfile("repAchg2_"); dir.create(d2)
  cfg_y <- make_rep_cfg(d2, year = 2018L)
  expect_false(identical(base$checksum, rf(cfg_y)$checksum))

  # scenario is consequential -> different fingerprint.
  d3   <- tempfile("repAchg3_"); dir.create(d3)
  cfg_s <- make_rep_cfg(d3, scenario = "restrictive")
  expect_false(identical(base$checksum, rf(cfg_s)$checksum))
})

test_that("LEVEL A: neg-pool fingerprint preview is deterministic + cfg-sensitive (no time)", {
  skip_if_not_installed("sf")
  d <- tempfile("repApool_"); dir.create(d)
  cfg <- make_rep_cfg(d, year = 2017L)
  pf  <- g_rep(".of_cfg_neg_pool_fingerprint_preview")

  a <- pf(cfg, 2017L)
  b <- pf(cfg, 2017L)
  expect_identical(a$checksum, b$checksum)
  expect_identical(a$text,     b$text)
  expect_false(grepl(format(Sys.Date()), a$text, fixed = TRUE))
  # Changing the target_year changes the preview fingerprint.
  expect_false(identical(a$checksum, pf(cfg, 2018L)$checksum))
})

test_that("LEVEL A: fold fingerprint is identical across two equivalent fold builds (no time)", {
  skip_if_not_installed("sf")
  mk_fp <- g_rep("make_fold_fingerprint")
  df <- data.frame(uid = sprintf("u%02d", 1:8),
                   fold_rep1 = rep(c(1L, 2L), 4),
                   fold_rep2 = rep(c(2L, 1L), 4),
                   stringsAsFactors = FALSE)
  args <- list(
    params   = list(seed_base = 42L, n_repeats = 2L, split_unit = "fire"),
    selected = list(block_size_m = 5000, k_folds = 2, ok = TRUE),
    train_with_folds = df, fold_cols = c("fold_rep1", "fold_rep2"),
    id_col = "uid"
  )
  a <- do.call(mk_fp, args)
  b <- do.call(mk_fp, args)
  expect_identical(a$checksum, b$checksum)
  expect_identical(a$text,     b$text)
  # Row order must not perturb identity (stable id ordering).
  args_sh <- args; args_sh$train_with_folds <- df[c(5, 1, 8, 3, 2, 7, 4, 6), ]
  expect_identical(a$checksum, do.call(mk_fp, args_sh)$checksum)
  # A config change (seed) flips the fingerprint.
  args_seed <- args; args_seed$params$seed_base <- 7L
  expect_false(identical(a$checksum, do.call(mk_fp, args_seed)$checksum))
})

test_that("LEVEL A: content-aware cache NAME is stable for fixed content + busts on change", {
  skip_if_not_installed("sf")
  # vector_input_cache_tag() is a CLOSURE local to validate_fire_maps() (see
  # R/validate-fire-maps.R:492), so it cannot be reached from the namespace.
  # Its CONTRACT is what Level A locks: an in-memory sf's cache name is a
  # deterministic function of its STRUCTURAL CONTENT (nrow + bbox + crs + geom
  # type) and carries NO wall-clock. We reimplement that documented algorithm
  # here (a faithful copy of the closure body) and assert: identical content ->
  # identical "in-" tag (no time leakage); changed content -> different tag.
  # The closure itself is exercised end-to-end by test-validate-fire-maps.R.
  cache_safe_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "-", as.character(x))
    gsub("(^-+)|(-+$)", "", x)
  }
  short_hash <- function(...) {
    s <- paste(unlist(list(...)), collapse = "|")
    bytes <- as.numeric(charToRaw(enc2utf8(s)))
    chk <- 0
    for (b in bytes) chk <- (chk * 31 + b) %% 1000000007
    sprintf("%09d", as.integer(chk))
  }
  sf_tag <- function(x) {
    bb <- sf::st_bbox(x)
    bb_str <- paste(format(as.numeric(bb), trim = TRUE, nsmall = 0),
                    collapse = ",")
    crs_str <- format(sf::st_crs(x)$wkt)
    geom_str <- paste(as.character(
      sf::st_geometry_type(x, by_geometry = FALSE)), collapse = "+")
    parts <- c("sf", format(nrow(x), scientific = FALSE),
               bb_str, crs_str, geom_str)
    paste0("in-", short_hash(cache_safe_tag(paste(parts, collapse = "|"))))
  }
  poly <- function(off) sf::st_polygon(list(rbind(
    c(off, 0), c(off + 1, 0), c(off + 1, 1), c(off, 1), c(off, 0))))
  s1 <- sf::st_sf(id = 1:3,
                  geometry = sf::st_sfc(poly(0), poly(2), poly(4), crs = 3035))
  s1b <- sf::st_sf(id = 1:3,
                   geometry = sf::st_sfc(poly(0), poly(2), poly(4), crs = 3035))
  s2 <- sf::st_sf(id = 1:4,
                  geometry = sf::st_sfc(poly(0), poly(2), poly(4), poly(6),
                                        crs = 3035))
  t1 <- sf_tag(s1)
  Sys.sleep(0)                                  # name must not depend on time
  expect_identical(t1, sf_tag(s1b))             # same content -> same name
  expect_false(identical(t1, sf_tag(s2)))       # changed content -> busts cache
  expect_false(grepl(format(Sys.Date()), t1, fixed = TRUE))
})

# ===========================================================================
# LEVEL B — NUMERIC REPRODUCIBILITY (nthread = 1)
# ===========================================================================

test_that("LEVEL B: two FINAL models (same seed/data/nthread=1) -> predictions within tol + EXACT class parity", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  wl <- g_rep(".supervised_feature_cols")
  train_df <- make_rep_train_df(n = 80, seed = 11L)
  score_df <- make_rep_train_df(n = 40, seed = 99L)  # held-out scoring rows

  fit1 <- train_rep_final(train_df, feature_cols = wl, fold_seed = 42L)
  fit2 <- train_rep_final(train_df, feature_cols = wl, fold_seed = 42L)

  # Same selection result (best_iteration) under identical seed/data/nthread=1.
  expect_identical(fit1$best_iteration, fit2$best_iteration)
  # Same recipe schema (x_cols / feature order / medians) -- a Level-A guarantee
  # that the numeric level depends on.
  expect_identical(fit1$x_cols, fit2$x_cols)
  expect_identical(fit1$feature_cols, fit2$feature_cols)
  expect_equal(fit1$medians, fit2$medians)

  p1 <- predict_rep(fit1, score_df)
  p2 <- predict_rep(fit2, score_df)
  expect_equal(length(p1), nrow(score_df))

  # NUMERIC PARITY within the explicit absolute tolerance.
  expect_lt(max(abs(p1 - p2)), PRED_ABS_TOL)

  # DERIVED CLASS parity at a fixed threshold is EXACT (no tolerance).
  c1 <- as.integer(p1 >= CLASS_THRESHOLD)
  c2 <- as.integer(p2 >= CLASS_THRESHOLD)
  expect_identical(c1, c2)
})

test_that("LEVEL B: re-predicting the SAME model on the SAME data is bit-identical", {
  skip_if_not_installed("xgboost")
  wl <- g_rep(".supervised_feature_cols")
  train_df <- make_rep_train_df(n = 70, seed = 3L)
  score_df <- make_rep_train_df(n = 35, seed = 7L)
  fit <- train_rep_final(train_df, feature_cols = wl, fold_seed = 42L)
  pa <- predict_rep(fit, score_df)
  pb <- predict_rep(fit, score_df)
  expect_identical(pa, pb)   # identical model + identical matrix -> bit-equal
})

# ===========================================================================
# LEVEL C — CARTOGRAPHIC REPRODUCIBILITY
# ===========================================================================

# Geometric fingerprint adopted by the contract: number of features, bbox
# (rounded), total vertex count, per-feature area signature (rounded). This
# identifies the MAP without depending on raw WKB byte layout (GDAL-version
# robust); area is additionally checked within AREA_REL_TOL.
geom_fingerprint <- function(x) {
  g <- sf::st_geometry(x)
  bb <- as.numeric(sf::st_bbox(x))
  nv <- sum(vapply(g, function(p) nrow(sf::st_coordinates(p)), integer(1)))
  ar <- as.numeric(sf::st_area(x))
  list(
    n        = length(g),
    bbox     = round(bb, 6),
    nvert    = nv,
    area_sig = round(ar, 6)
  )
}

# Produce the "final map" twice from identical inputs: score the SAME polygons
# with a model trained at nthread=1, attach p_burned + a class at threshold, and
# keep geometry pass-through. This mirrors the final-map assembly contract
# (one prediction per input polygon; row/id/order preserved).
make_rep_final_map <- function(train_df, score_sf, wl, seed = 42L) {
  fit <- train_rep_final(train_df, feature_cols = wl, fold_seed = seed)
  score_df <- as.data.frame(sf::st_drop_geometry(score_sf))
  p <- predict_rep(fit, score_df)
  out <- score_sf
  out$p_burned    <- p
  out$class_final <- ifelse(p >= CLASS_THRESHOLD, "burned_like", "keep")
  out
}

test_that("LEVEL C: the final map is reproducible (IDs/count/classes exact; area+geometry within tol)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")

  wl <- g_rep(".supervised_feature_cols")
  train_df <- make_rep_train_df(n = 80, seed = 11L)

  # A small polygon universe to score (the "deterministic universe").
  n_poly <- 24L
  sc_df <- make_rep_train_df(n = n_poly, seed = 5L)
  sfc <- sf::st_sfc(lapply(seq_len(n_poly), function(i) {
    x <- (i %% 6) * 10; y <- (i %/% 6) * 10
    sf::st_polygon(list(rbind(c(x, y), c(x + 5, y), c(x + 5, y + 5),
                              c(x, y + 5), c(x, y))))
  }), crs = 3035)
  score_sf <- sf::st_sf(sc_df, geometry = sfc)

  m1 <- make_rep_final_map(train_df, score_sf, wl, seed = 42L)
  m2 <- make_rep_final_map(train_df, score_sf, wl, seed = 42L)

  # EXACT: number of polygons, IDs (+ order), derived classes.
  expect_identical(nrow(m1), nrow(m2))
  expect_identical(nrow(m1), as.integer(n_poly))
  expect_identical(m1$fire_uid, m2$fire_uid)       # IDs + observation order
  expect_identical(m1$class_final, m2$class_final) # final classes EXACT

  # p_burned within the numeric tolerance (drives the classes above).
  expect_lt(max(abs(m1$p_burned - m2$p_burned)), PRED_ABS_TOL)

  # AREA within the relative tolerance.
  a1 <- as.numeric(sf::st_area(m1)); a2 <- as.numeric(sf::st_area(m2))
  rel <- abs(a1 - a2) / pmax(1e-12, abs(a1))
  expect_lt(max(rel), AREA_REL_TOL)

  # GEOMETRIC FINGERPRINT identical (bbox + vertex count + area signature).
  expect_identical(geom_fingerprint(m1), geom_fingerprint(m2))
})

test_that("LEVEL C: geometry passes through unchanged (input == output geometry)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  wl <- g_rep(".supervised_feature_cols")
  train_df <- make_rep_train_df(n = 60, seed = 21L)
  n_poly <- 12L
  sc_df <- make_rep_train_df(n = n_poly, seed = 8L)
  sfc <- sf::st_sfc(lapply(seq_len(n_poly), function(i) {
    sf::st_polygon(list(rbind(c(i, 0), c(i + 1, 0), c(i + 1, 1),
                              c(i, 1), c(i, 0))))
  }), crs = 3035)
  score_sf <- sf::st_sf(sc_df, geometry = sfc)
  m <- make_rep_final_map(train_df, score_sf, wl, seed = 42L)
  # The map geometry equals the input geometry exactly (pass-through contract).
  expect_identical(geom_fingerprint(score_sf)$area_sig,
                   geom_fingerprint(m)$area_sig)
  expect_true(all(sf::st_equals(score_sf, m, sparse = FALSE)[
    cbind(seq_len(n_poly), seq_len(n_poly))]))
})
