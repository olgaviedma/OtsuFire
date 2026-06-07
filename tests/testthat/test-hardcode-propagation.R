# 2026-06-05: SUPERVISED HARDCODE AUDIT — D2 (overwrite honesty) + D1
# (training-knob expose) + FRENTE 1 (OOF/FINAL hyperparameter UNIFICATION).
# These tests assert (a) the D1 training-control defaults now equal the
# UNIFIED CANONICAL values at every layer, (b) the OOF and FINAL `params = NULL`
# blocks resolve to the SAME canonical params (modulo scale_pos_weight) via the
# single source of truth, (c) an override actually reaches the engine, and
# (d) the D2 overwrite flag reaches the folds stage. They use formals()
# introspection and argument-capture stubs, so no ~60-min pipeline or raster
# stack is needed.
#
# FRENTE 1 (2026-06-05) updated these tests: the previous assertions pinned the
# defaults to the OLD asymmetric engine values (OOF nrounds 3000 / early_stop 75;
# FINAL seeds 999) and asserted the two were NOT unified. Natalia signed off on
# unifying them, so those assertions are now flipped to the canonical values
# (OOF nrounds 4000 / early_stop 80; all seeds 42) and to "OOF == FINAL".

ns <- asNamespace("OtsuFire")

# Canonical training-control values (FRENTE 1 single source of truth).
.CANON_NROUNDS    <- 4000
.CANON_EARLY_STOP <- 80
.CANON_SEED       <- 42

# --------------------------------------------------------------------------
# (a) Gate 1B (2026-06-07): cfg is now the SINGLE SOURCE OF TRUTH. The public
#     functions default their methodological knobs to NULL ("read from cfg")
#     and the internal-pure engines carry NO methodological defaults at all
#     (the args are REQUIRED). So the canonical values live ONLY in
#     build_supervised_burned_config() (.of_canonical_train_control /
#     .of_canonical_model_params), which is what these tests now assert.
# --------------------------------------------------------------------------
.canon_tc <- get(".of_canonical_train_control", envir = ns)()

test_that("Gate 1B: canonical train_control holds the UNIFIED canonical values", {
  expect_equal(.canon_tc$nrounds_max, .CANON_NROUNDS)
  expect_equal(.canon_tc$early_stop,  .CANON_EARLY_STOP)
  expect_equal(.canon_tc$seeds$oof_seed_base,       .CANON_SEED)
  expect_equal(.canon_tc$seeds$final_sampling_seed, .CANON_SEED)
  expect_equal(.canon_tc$seeds$final_seed,          .CANON_SEED)
  expect_equal(.canon_tc$val_frac, 0.15)
  expect_equal(.canon_tc$group_col, "block_id")
  expect_equal(.canon_tc$impute_numeric, "median")
  expect_equal(.canon_tc$impute_factor_missing, "MISSING")
  expect_equal(.canon_tc$caps$contextual, 0.25)
  expect_equal(.canon_tc$caps$spectral, 1.0)
  expect_equal(.canon_tc$caps$random, 1.0)
  expect_equal(.canon_tc$caps$otsu, 1.0)
  expect_equal(.canon_tc$training_protocol, "legacy")
  expect_equal(.canon_tc$oof_sampling, "capped")
})

test_that("Gate 1B: public methodological knobs default to NULL (= read from cfg)", {
  f_oof   <- formals(get("run_oof_diagnostics", envir = ns))
  f_final <- formals(get("train_final_burned_model", envir = ns))
  f_run   <- formals(get("run_oneyear_supervised_pipeline", envir = ns))
  for (nm in c("nrounds_max", "early_stop", "seed_base")) {
    expect_null(f_oof[[nm]])
  }
  for (nm in c("sampling_seed", "seed", "val_frac", "group_col", "nrounds_max",
               "early_stopping_rounds", "impute_numeric", "impute_factor_missing")) {
    expect_null(f_final[[nm]])
  }
  for (nm in c("oof_nrounds_max", "oof_early_stop", "oof_seed_base",
               "final_sampling_seed", "final_seed", "final_val_frac",
               "final_group_col", "final_nrounds_max", "final_early_stopping_rounds")) {
    expect_null(f_run[[nm]])
  }
})

test_that("Gate 1B: internal-pure engines carry NO methodological defaults (required args)", {
  f_eng_final <- formals(get("train_final_model_direct", envir = ns))
  f_eng_oof   <- formals(get("run_dm_oof_pipeline", envir = ns))
  # A required formal has an empty symbol as its "default".
  is_required <- function(x) is.symbol(x) && !nzchar(as.character(x))
  for (nm in c("nrounds_max", "early_stopping_rounds", "sampling_seed", "seed",
               "val_frac", "group_col", "impute_numeric", "impute_factor_missing",
               "model_params_base")) {
    expect_true(is_required(f_eng_final[[nm]]),
                info = paste("train_final_model_direct$", nm, "should be required"))
  }
  for (nm in c("nrounds_max", "early_stop", "seed_base")) {
    expect_true(is_required(f_eng_oof[[nm]]),
                info = paste("run_dm_oof_pipeline$", nm, "should be required"))
  }
})

test_that("FRENTE 1: OOF and FINAL params=NULL resolve to the SAME canonical block (modulo scale_pos_weight)", {
  # Both inline sites build from cfg$model_params (sourced from
  # .of_canonical_model_params()); assert the canonical builder is the single
  # source of truth and that the two sites produce identical lists apart from
  # scale_pos_weight, which each computes from its own labels.
  build <- get(".of_canonical_xgb_params", envir = ns)
  p_oof   <- build(scale_pos_weight = 3.0)   # pretend OOF training labels
  p_final <- build(scale_pos_weight = 7.0)   # pretend FINAL training labels
  # Everything except scale_pos_weight must be identical.
  drop_spw <- function(p) p[setdiff(names(p), "scale_pos_weight")]
  expect_identical(drop_spw(p_oof), drop_spw(p_final))
  # Canonical field values.
  expect_equal(p_oof$objective, "binary:logistic")
  expect_identical(p_oof$eval_metric, c("logloss", "aucpr")) # logloss FIRST
  expect_equal(p_oof$eval_metric[1], "logloss")              # drives early stop
  expect_equal(p_oof$eta, 0.05)
  expect_equal(p_oof$max_depth, 5)
  expect_equal(p_oof$min_child_weight, 5)
  expect_equal(p_oof$subsample, 0.8)
  expect_equal(p_oof$colsample_bytree, 0.75)
  expect_equal(p_oof$gamma, 0)
  expect_equal(p_oof$lambda, 1)
  expect_equal(p_oof$alpha, 0)
  # scale_pos_weight is the only per-site field.
  expect_equal(p_oof$scale_pos_weight, 3.0)
  expect_equal(p_final$scale_pos_weight, 7.0)
})

test_that("Gate 1B: top-level entry-point knobs default to NULL and a default cfg resolves to canonical", {
  # Gate 1B: the public entry no longer carries literal canonical defaults; the
  # knobs default to NULL ("read from cfg"). The canonical VALUES are resolved by
  # build_supervised_burned_config() into cfg$train_control, asserted here.
  f <- formals(get("run_oneyear_supervised_pipeline", envir = ns))
  for (nm in c("oof_nrounds_max", "oof_early_stop", "oof_seed_base",
               "final_sampling_seed", "final_seed", "final_val_frac",
               "final_group_col", "final_nrounds_max",
               "final_early_stopping_rounds", "final_impute_numeric",
               "final_impute_factor_missing")) {
    expect_null(f[[nm]])
  }
  tc <- .canon_tc
  expect_equal(tc$nrounds_max, .CANON_NROUNDS)
  expect_equal(tc$early_stop, .CANON_EARLY_STOP)
  expect_equal(tc$seeds$oof_seed_base, .CANON_SEED)
  expect_equal(tc$seeds$final_sampling_seed, .CANON_SEED)
  expect_equal(tc$seeds$final_seed, .CANON_SEED)
  expect_equal(tc$val_frac, 0.15)
  expect_equal(tc$group_col, "block_id")
  expect_equal(tc$impute_numeric, "median")
  expect_equal(tc$impute_factor_missing, "MISSING")
})

# --------------------------------------------------------------------------
# (b) An override REACHES the engine. We stub the engine functions with an
#     argument-capture closure via testthat::local_mocked_bindings and call the
#     public wrapper just far enough to hit the engine call. The wrappers do a
#     little I/O before the engine call, so we feed minimal on-disk inputs.
# --------------------------------------------------------------------------
test_that("D1 FINAL override reaches train_final_model_direct", {
  skip_if_not_installed("sf")
  testthat::skip_if(
    !exists("local_mocked_bindings", where = asNamespace("testthat")),
    "testthat::local_mocked_bindings not available"
  )

  captured <- new.env()
  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      args <- list(...)
      captured$nrounds_max <- args$nrounds_max
      captured$seed <- args$seed
      captured$val_frac <- args$val_frac
      captured$group_col <- args$group_col
      captured$early_stopping_rounds <- args$early_stopping_rounds
      captured$sampling_seed <- args$sampling_seed
      # Return a minimal shape so the wrapper's post-processing does not error.
      list(model = NULL, recipe = NULL, training_ok_sf = NULL, split = NULL,
           files = list())
    },
    .package = "OtsuFire"
  )

  # Build a minimal labelled GPKG + config so the wrapper reaches the engine.
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                              c(0,1), c(0,0)))), crs = 3035)
  tf <- sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc)
  tf_gpkg <- tempfile(fileext = ".gpkg")
  sf::st_write(tf, tf_gpkg, layer = "train_features", quiet = TRUE,
               delete_dsn = TRUE)
  ci <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), ci,
                     overwrite = TRUE)
  id <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), id, quiet = TRUE,
               delete_dsn = TRUE)
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )

  # Precision 1 (2026-06-07): these are DEPRECATED function-level shims; using
  # them is expected to emit "otsufire_deprecated_param" deprecation warnings
  # (one per overridden field) while STILL propagating. Muffle+count them so the
  # suite's incidental-warning total stays clean, assert at least one fired, then
  # confirm the values reached the engine.
  .ndep <- 0L
  withCallingHandlers(
    train_final_burned_model(
      train_features = tf_gpkg,
      config = cfg,
      out_dir = tempfile(),
      nrounds_max = 1234,
      seed = 7,
      val_frac = 0.42,
      group_col = "poly_id",
      early_stopping_rounds = 11,
      sampling_seed = 13
    ),
    otsufire_deprecated_param = function(w) {
      .ndep <<- .ndep + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(.ndep, 0L)

  expect_equal(captured$nrounds_max, 1234)
  expect_equal(captured$seed, 7)
  expect_equal(captured$val_frac, 0.42)
  expect_equal(captured$group_col, "poly_id")
  expect_equal(captured$early_stopping_rounds, 11)
  expect_equal(captured$sampling_seed, 13)
})

test_that("D1 OOF override reaches run_dm_oof_pipeline", {
  skip_if_not_installed("sf")
  testthat::skip_if(
    !exists("local_mocked_bindings", where = asNamespace("testthat")),
    "testthat::local_mocked_bindings not available"
  )

  captured <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      args <- list(...)
      captured$nrounds_max <- args$nrounds_max
      captured$early_stop <- args$early_stop
      captured$seed_base <- args$seed_base
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL),
           files = list())
    },
    .package = "OtsuFire"
  )

  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1),
                                              c(0,1), c(0,0)))), crs = 3035)
  tf <- sf::st_sf(fire_uid = "a", class = "burned",
                  fold_rep1 = 1L, fold_rep2 = 1L, geometry = sfc)
  ci <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), ci,
                     overwrite = TRUE)
  id <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), id, quiet = TRUE,
               delete_dsn = TRUE)
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = id,
    change_index = ci, target_year = 2025L
  )

  # Precision 1 (2026-06-07): DEPRECATED function-level shims -> muffle+count the
  # "otsufire_deprecated_param" warnings while asserting propagation.
  .ndep <- 0L
  withCallingHandlers(
    run_oof_diagnostics(
      train_features = tf,
      scoring_features = tf,
      config = cfg,
      out_dir = tempfile(), matrix_dir = tempfile(),
      nrounds_max = 321,
      early_stop = 9,
      seed_base = 1234
    ),
    otsufire_deprecated_param = function(w) {
      .ndep <<- .ndep + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(.ndep, 0L)

  expect_equal(captured$nrounds_max, 321)
  expect_equal(captured$early_stop, 9)
  expect_equal(captured$seed_base, 1234)
})

# --------------------------------------------------------------------------
# (c) D2: the orchestrator forwards the user's `overwrite` flag to the folds
#     stage make_spatial_folds() instead of a hardcoded TRUE. We assert it via
#     source inspection of the delegation call (the engine itself needs the
#     full raster stack to run, which is out of scope for a unit test).
# --------------------------------------------------------------------------
test_that("D2: orchestrator forwards overwrite (not a literal TRUE) to make_spatial_folds", {
  lines <- deparse(get("run_supervised_pipeline", envir = ns))
  expect_true(any(grepl("make_spatial_folds", lines, fixed = TRUE)))
  # Locate the make_spatial_folds(...) delegation and inspect the lines from the
  # call open-paren to the next close-paren at the call's indentation. We assert
  # the call forwards `overwrite = overwrite` and that no `overwrite = TRUE`
  # literal appears inside this delegation block. (The orchestrator's own
  # signature default `overwrite = TRUE` lives elsewhere in the deparse.)
  i_open <- which(grepl("make_spatial_folds\\(", lines))[1]
  expect_false(is.na(i_open))
  # The folds-call argument lines are between the open-paren line and the line
  # that closes it (the `time_step(...)` block). Scan forward until we see the
  # `overwrite` argument; assert the FIRST overwrite argument after the folds
  # call open is the forwarded variable, not a literal.
  tail_lines <- lines[i_open:min(length(lines), i_open + 25L)]
  ow_lines <- tail_lines[grepl("overwrite", tail_lines)]
  expect_true(length(ow_lines) >= 1L)
  expect_true(grepl("overwrite\\s*=\\s*overwrite", ow_lines[1]))
  expect_false(grepl("overwrite\\s*=\\s*TRUE", ow_lines[1]))
})
