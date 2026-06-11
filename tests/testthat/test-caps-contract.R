# =============================================================================
# Gate 1D.6 PART A -- NEGATIVE-POOL CAPS CONTRACT
# =============================================================================
# Contract for the FOUR negative-bucket caps (contextual / spectral / random /
# otsu). For each cap this file proves, end-to-end:
#
#   (1) the value STORED in cfg               -> cfg$train_control$caps$<name>;
#   (2) the value RECEIVED by the OOF engine  -> the resolved ratio handed to
#                                                run_dm_oof_pipeline();
#   (3) the value RECEIVED by the FINAL engine -> the resolved ratio handed to
#                                                 train_final_model_direct();
#   (4) the EFFECTIVE cap computed in the pool builder:
#           effective_cap = ceiling(n_burned * ratio)             [the FORMULA]
#       recovered from R/internal-sup-train-final-direct.R:222-226 (FINAL) and
#       R/internal-sup-oof-training.R:211-218,245-248 (OOF) -- BOTH identical;
#   (5) the NUMBER AVAILABLE in the bucket (rows of the bucket pool);
#   (6) the NUMBER SELECTED after the cap     -> n_selected = min(eff_cap, n_avail);
#   (7) NO residual default silently re-supplies a cap: a dropped required cap
#       ERRORS at the pool boundary (tie to Gate 1B single-source guarantees);
#   (8) ERROR on discrepancy: OOF and FINAL must receive the SAME cap as cfg;
#       the package-level parity guard `.of_assert_spectral_cap_parity`
#       (R/internal-sup-shim-resolver.R:169) fails fast on a silent reversion.
#
# Phase B is covered EXPLICITLY: a cfg built with cap_spectral = 2.0 reaches
# BOTH OOF and FINAL as 2.0 and does NOT silently revert to the 1.0 canonical
# default -- and the parity guard accepts 2.0/2.0/2.0 while rejecting any 1.0
# reversion.
#
# The cfg->OOF / cfg->FINAL CAPTURE parity (full single-source propagation) is
# established in test-cfg-single-source.R; this file does NOT duplicate it.
# Here we (a) bind the recovered effective-cap FORMULA to an executable
# assertion (n_selected <= effective_cap, selection respects it), (b) exercise
# the parity guard directly (incl. Phase B no-revert), and (c) prove the
# cfg-stored cap equals the resolved cap received by both engines through the
# single shim resolver, with no residual default.
# =============================================================================

ns <- asNamespace("OtsuFire")
`%||%` <- function(x, y) if (is.null(x)) y else x  # local mirror (internal-only)

# ---- small reusable cfg fixtures (Gate 1D.7: shared via
#      helper-supervised-contracts.R; kept as thin aliases so the assertions
#      below are untouched) -----------------------------------------------------
.caps_tif  <- sc_change_index_tif
.caps_gpkg <- sc_decisions_gpkg
.caps_cfg  <- sc_min_cfg

# The four canonical (default) cap ratios, in the bucket order used everywhere
# (contextual, spectral, random, otsu).
CANON_CAPS <- c(contextual = 0.25, spectral = 1.0, random = 1.0, otsu = 1.0)

# =============================================================================
# (4)+(5)+(6) EFFECTIVE-CAP FORMULA, recovered + bound to an assertion.
#     effective_cap = ceiling(n_burned * ratio); n_selected = min(eff, n_avail).
#     This faithfully reproduces the FINAL `target_* <- ceiling(n_burned*ratio)`
#     + `min(target, nrow(pool))` selection AND the OOF `pick()` closure. We
#     assert the invariant the production code must satisfy for every bucket.
# =============================================================================

# Faithful reproduction of the single selection rule used by BOTH engines.
effective_cap <- function(n_burned, ratio) ceiling(n_burned * ratio)
n_selected_after_cap <- function(n_burned, ratio, n_available) {
  if (!is.finite(ratio)) return(n_available)        # Inf disables the cap
  target <- effective_cap(n_burned, ratio)
  if (target <= 0L) return(0L)
  min(target, n_available)
}

test_that("PART A: effective-cap formula = ceiling(n_burned * ratio); n_selected <= eff_cap", {
  # A grid of (n_burned, ratio, n_available) spanning the regimes:
  #   under-cap (avail < eff), over-cap (avail > eff), exact, zero ratio, Inf.
  grid <- expand.grid(
    n_burned    = c(0L, 1L, 7L, 40L, 123L),
    ratio       = c(0, 0.25, 1.0, 2.0, Inf),
    n_available = c(0L, 3L, 50L, 1000L)
  )
  for (i in seq_len(nrow(grid))) {
    nb <- grid$n_burned[i]; rt <- grid$ratio[i]; av <- grid$n_available[i]
    eff <- effective_cap(nb, rt)
    sel <- n_selected_after_cap(nb, rt, av)

    # The formula itself: ceiling, never floor; 0 ratio -> 0 cap.
    if (is.finite(rt)) {
      expect_identical(eff, ceiling(nb * rt),
                       info = sprintf("nb=%d ratio=%g", nb, rt))
    }
    # INVARIANT 1: selection never exceeds availability.
    expect_lte(sel, av)
    # INVARIANT 2 (the cap is RESPECTED): when finite, selection never exceeds
    # the effective cap.
    if (is.finite(rt)) expect_lte(sel, eff)
    # INVARIANT 3: when availability is the binding constraint, take all of it.
    if (is.finite(rt) && av <= eff && eff > 0L) expect_equal(sel, av)
    # INVARIANT 4: when the cap is the binding constraint, take exactly the cap.
    if (is.finite(rt) && eff > 0L && av > eff) expect_equal(sel, eff)
    # INVARIANT 5: Inf disables -> take everything available.
    if (!is.finite(rt)) expect_equal(sel, av)
    # INVARIANT 6: zero ratio -> select nothing.
    if (rt == 0) expect_equal(sel, 0)
  }
})

test_that("PART A: the production FINAL pool builder uses ceiling(n_burned*ratio)", {
  # Source-level guard: the recovered formula is what the FINAL engine actually
  # computes (target_<bucket> <- ceiling(n_burned * <ratio>)), one per bucket,
  # and the selection is min(target, nrow(pool)).
  src <- paste(deparse(body(get("train_final_model_direct", envir = ns))),
               collapse = "\n")
  expect_true(grepl("ceiling(n_burned * contextual_exclusion_to_burned_ratio)",
                    src, fixed = TRUE))
  expect_true(grepl("ceiling(n_burned * spectral_hard_negative_to_burned_ratio)",
                    src, fixed = TRUE))
  expect_true(grepl("ceiling(n_burned * random_to_burned_ratio)",
                    src, fixed = TRUE))
  expect_true(grepl("ceiling(n_burned * otsu_unburned_to_burned_ratio)",
                    src, fixed = TRUE))
  # The selection is capped by min(target, nrow(pool)).
  expect_true(grepl("min(target_contextual, nrow(contextual_exclusion_pool))",
                    src, fixed = TRUE))
})

test_that("PART A: the production OOF cap closure uses the SAME ceiling(n_burned*ratio)", {
  src <- paste(deparse(body(get("run_oof_xgb", envir = ns))), collapse = "\n")
  # The pick() closure computes target <- ceiling(n_burned * ratio) and takes
  # n_take <- min(target, length(avail)) -- byte-identical rule to FINAL.
  expect_true(grepl("ceiling(n_burned * ratio)", src, fixed = TRUE))
  expect_true(grepl("min(target, length(avail))", src, fixed = TRUE))
  # The per-bucket audit records the SAME ceiling formula for each cap, so the
  # OOF audit's *_cap columns are the effective caps, not the raw ratios.
  expect_true(grepl("ceiling(n_burned * contextual_exclusion_to_burned_ratio)",
                    src, fixed = TRUE))
  expect_true(grepl("ceiling(n_burned * spectral_hard_negative_to_burned_ratio)",
                    src, fixed = TRUE))
})

# =============================================================================
# (1)->(2)/(3) cfg-STORED cap == cap RECEIVED by OOF and FINAL, via the single
#     shim resolver, with NO residual default. We resolve each cap exactly the
#     way the public boundary does (.of_resolve_methodological_shim with the
#     override NULL) and assert the resolved scalar is the cfg cap.
# =============================================================================

test_that("PART A: each cfg cap resolves (override=NULL) to itself -- no residual default", {
  resolve <- get(".of_resolve_methodological_shim", envir = ns)

  for (cfg in list(.caps_cfg(),                                   # all default
                   .caps_cfg(cap_contextual = 0.5, cap_spectral = 2.0,
                             cap_random = 0.75, cap_otsu = 1.5))) {  # all user
    caps <- cfg$train_control$caps
    prov <- cfg$resolved_params_provenance$train_control %||% list()
    pairs <- list(
      list("cap_contextual", caps$contextual),
      list("cap_spectral",   caps$spectral),
      list("cap_random",     caps$random),
      list("cap_otsu",       caps$otsu)
    )
    for (p in pairs) {
      key <- p[[1]]; cfg_val <- p[[2]]
      # override = NULL -> the resolver MUST return the cfg value unchanged.
      out <- suppressWarnings(resolve(
        override = NULL, cfg_value = cfg_val, param = key,
        cfg_provenance = prov[[key]] %||% "default", record = NULL))
      expect_equal(out, cfg_val, info = key)
    }
  }
})

test_that("PART A: a conflicting explicit cap (cfg-user vs function override) ERRORS, never silently picks", {
  resolve <- get(".of_resolve_methodological_shim", envir = ns)
  # cfg has an EXPLICIT user spectral cap of 2.0; a DIFFERENT function-level
  # override (1.0) is two explicit incompatible sources -> hard error.
  expect_error(
    resolve(override = 1.0, cfg_value = 2.0, param = "cap_spectral",
            cfg_provenance = "user", record = NULL),
    regexp = "Conflicting supervised parameter 'cap_spectral'"
  )
  # The SAME value is accepted silently (the two explicit sources agree).
  expect_silent(
    out <- resolve(override = 2.0, cfg_value = 2.0, param = "cap_spectral",
                   cfg_provenance = "user", record = NULL))
  expect_equal(out, 2.0)
})

# =============================================================================
# (8) PARITY GUARD: OOF and FINAL must receive the SAME spectral cap as cfg.
#     Direct exercise of `.of_assert_spectral_cap_parity` (the package-level
#     P2 guard the orchestrator runs before any model is fit).
# =============================================================================

test_that("PART A: spectral-cap parity guard PASSES when OOF == FINAL == cfg", {
  guard <- get(".of_assert_spectral_cap_parity", envir = ns)
  # Canonical (1.0) and Phase B (2.0) both pass when all three agree.
  expect_silent(guard(oof_spectral = 1.0, final_spectral = 1.0, cfg_spectral = 1.0))
  expect_equal(guard(oof_spectral = 2.0, final_spectral = 2.0, cfg_spectral = 2.0), 2.0)
})

test_that("PART A: spectral-cap parity guard FAILS on ANY divergence (OOF!=FINAL, OOF!=cfg, FINAL!=cfg)", {
  guard <- get(".of_assert_spectral_cap_parity", envir = ns)
  # OOF diverges from FINAL.
  expect_error(guard(oof_spectral = 2.0, final_spectral = 1.0, cfg_spectral = 2.0),
               regexp = "parity guard FAILED")
  # FINAL silently reverted to 1.0 while cfg/OOF are 2.0 (the classic regression).
  expect_error(guard(oof_spectral = 2.0, final_spectral = 1.0, cfg_spectral = 1.0),
               regexp = "parity guard FAILED")
  # cfg diverges from the (equal) engines.
  expect_error(guard(oof_spectral = 2.0, final_spectral = 2.0, cfg_spectral = 1.0),
               regexp = "parity guard FAILED")
  # A malformed (NULL/NA/non-scalar) cap is rejected, never coerced.
  expect_error(guard(oof_spectral = NULL, final_spectral = 1.0, cfg_spectral = 1.0),
               regexp = "single .*finite numeric")
  expect_error(guard(oof_spectral = NA_real_, final_spectral = 1.0, cfg_spectral = 1.0),
               regexp = "single .*finite numeric")
})

# =============================================================================
# PHASE B (EXPLICIT): cap_spectral = 2.0 reaches cfg, both engines, and the
#     parity guard -- and does NOT silently revert to 1.0.
# =============================================================================

test_that("PART A / Phase B: cap_spectral=2.0 is STORED in cfg as 2.0 (not the 1.0 default)", {
  cfg <- .caps_cfg(cap_spectral = 2.0)
  expect_equal(cfg$train_control$caps$spectral, 2.0)
  # Provenance records this as a USER value, not a default.
  prov <- cfg$resolved_params_provenance$train_control
  expect_equal(prov$cap_spectral, "user")
  # The OTHER three caps remain canonical (no collateral reversion/override).
  expect_equal(cfg$train_control$caps$contextual, CANON_CAPS[["contextual"]])
  expect_equal(cfg$train_control$caps$random,     CANON_CAPS[["random"]])
  expect_equal(cfg$train_control$caps$otsu,       CANON_CAPS[["otsu"]])
})

test_that("PART A / Phase B: cap_spectral=2.0 reaches BOTH OOF and FINAL and the guard accepts 2.0 (no revert to 1.0)", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))

  cfg <- .caps_cfg(cap_spectral = 2.0)
  cfg_spec <- cfg$train_control$caps$spectral
  expect_equal(cfg_spec, 2.0)

  # Capture the spectral cap each engine RECEIVES (mirror the single-source
  # capture pattern; here we only need the spectral slot). Gate 1D.7: the tiny
  # train fixtures are the shared sc_* fixtures (identical content).
  mk_oof_feats  <- sc_burned_train_sf
  mk_final_gpkg <- sc_burned_train_gpkg

  cap_env <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      cap_env$oof_spectral <- list(...)$spectral_hard_negative_to_burned_ratio
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL), files = list())
    },
    .package = "OtsuFire"
  )
  run_oof_diagnostics(train_features = mk_oof_feats(),
                      scoring_features = mk_oof_feats(), config = cfg,
                      out_dir = tempfile(), matrix_dir = tempfile())

  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      cap_env$final_spectral <- list(...)$spectral_hard_negative_to_burned_ratio
      list(model = NULL, recipe = NULL, training_ok_sf = NULL, split = NULL,
           files = list())
    },
    .package = "OtsuFire"
  )
  train_final_burned_model(train_features = mk_final_gpkg(), config = cfg,
                           out_dir = tempfile())

  # The Phase B cap reaches BOTH engines as 2.0 (NOT the 1.0 canonical default).
  expect_equal(cap_env$oof_spectral, 2.0)
  expect_equal(cap_env$final_spectral, 2.0)

  # The package-level parity guard ACCEPTS the Phase B triple (no revert): with
  # the captured engine caps + the cfg cap all 2.0, the guard passes silently.
  guard <- get(".of_assert_spectral_cap_parity", envir = ns)
  expect_equal(
    guard(oof_spectral = cap_env$oof_spectral,
          final_spectral = cap_env$final_spectral,
          cfg_spectral = cfg_spec),
    2.0)

  # And it would FAIL the moment any layer silently reverted FINAL to 1.0.
  expect_error(
    guard(oof_spectral = cap_env$oof_spectral, final_spectral = 1.0,
          cfg_spectral = cfg_spec),
    regexp = "parity guard FAILED")
})

# =============================================================================
# (7) NO RESIDUAL DEFAULT: a dropped required cap ERRORS at the pool boundary
#     (cannot silently revert a bucket to ratio 1.0). Tie to Gate 1B.
# =============================================================================

test_that("PART A: dropping a required cap ERRORS in the FINAL pool boundary (no silent 1.0 revert)", {
  fn <- get("train_final_model_direct", envir = ns)
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc),
               g, layer = "train_features", quiet = TRUE, delete_dsn = TRUE)

  # Supply every required methodological arg EXCEPT the spectral cap -> must
  # error (the required-arg guard refuses a silent default for the bucket).
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = tempfile(), prefix = "x", overwrite = TRUE, verbose = FALSE,
      contextual_exclusion_to_burned_ratio = 1,
      # spectral_hard_negative_to_burned_ratio DROPPED on purpose.
      random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
      sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
      nrounds_max = 80, early_stopping_rounds = 80, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      model_params_base = get(".of_canonical_model_params", envir = ns)()
    )),
    regexp = "spectral_hard_negative_to_burned_ratio"
  )
})

test_that("PART A: dropping a required cap ERRORS on the OOF path (no silent 1.0 revert)", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fn <- get("run_oof_xgb", envir = ns)
  # The four cap ratios are REQUIRED formals; a dropped one ERRORS rather than
  # reverting the bucket to 1.0.
  expect_error(
    suppressMessages(suppressWarnings(fn(
      XL_mat = Matrix::Matrix(matrix(0, 2, 1), sparse = TRUE), y = c(1L, 0L),
      labelled_df = data.frame(fire_uid = c("a", "b"),
                               class = c("burned", "unburned"),
                               block_id = c(1L, 2L), fold_rep1 = c(1L, 2L)),
      fold_cols = "fold_rep1",
      params = list(), nrounds_max = 4L, early_stop = 3L, seed_base = 42L,
      group_col = "block_id", val_frac = 0.15, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      prepared_labelled = data.frame(x = c(0, 1)), model_cols = "x",
      contextual_exclusion_to_burned_ratio = 1,
      # spectral_hard_negative_to_burned_ratio DROPPED on purpose.
      random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1
    ))),
    regexp = "spectral_hard_negative_to_burned_ratio"
  )
})
