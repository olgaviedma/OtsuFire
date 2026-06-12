# =============================================================================
# Gate 1D.6 PART A -- NEGATIVE-POOL CAPS CONTRACT
# =============================================================================
# Contract for the TWO negative-bucket caps (random / otsu).
# (GATE 6.1 removed the dead spectral bucket; GATE 6.5 removed the contextual
# (deterministic-drop) bucket end-to-end.) For each cap this file proves,
# end-to-end:
#
#   (1) the value STORED in cfg               -> cfg$train_control$caps$<name>;
#   (2) the value RECEIVED by the OOF engine  -> the resolved ratio handed to
#                                                run_dm_oof_pipeline();
#   (3) the value RECEIVED by the FINAL engine -> the resolved ratio handed to
#                                                 train_final_model_direct();
#   (4) the EFFECTIVE cap computed in the pool builder:
#           effective_cap = ceiling(n_burned * ratio)             [the FORMULA]
#   (5) the NUMBER AVAILABLE in the bucket (rows of the bucket pool);
#   (6) the NUMBER SELECTED after the cap     -> n_selected = min(eff_cap, n_avail);
#   (7) NO residual default silently re-supplies a cap: a dropped required cap
#       ERRORS at the pool boundary (tie to Gate 1B single-source guarantees).
#
# The cfg->OOF / cfg->FINAL CAPTURE parity (full single-source propagation) is
# established in test-cfg-single-source.R; this file does NOT duplicate it.
# =============================================================================

ns <- asNamespace("OtsuFire")
`%||%` <- function(x, y) if (is.null(x)) y else x  # local mirror (internal-only)

# ---- small reusable cfg fixtures (Gate 1D.7: shared via
#      helper-supervised-contracts.R; kept as thin aliases so the assertions
#      below are untouched) -----------------------------------------------------
.caps_tif  <- sc_change_index_tif
.caps_gpkg <- sc_decisions_gpkg
.caps_cfg  <- sc_min_cfg

# The two canonical (default) cap ratios, in the bucket order used everywhere
# (random, otsu). GATE 6.5: the contextual cap was removed.
CANON_CAPS <- c(random = 1.0, otsu = 1.0)

# =============================================================================
# (4)+(5)+(6) EFFECTIVE-CAP FORMULA, recovered + bound to an assertion.
#     effective_cap = ceiling(n_burned * ratio); n_selected = min(eff, n_avail).
# =============================================================================

effective_cap <- function(n_burned, ratio) ceiling(n_burned * ratio)
n_selected_after_cap <- function(n_burned, ratio, n_available) {
  if (!is.finite(ratio)) return(n_available)        # Inf disables the cap
  target <- effective_cap(n_burned, ratio)
  if (target <= 0L) return(0L)
  min(target, n_available)
}

test_that("PART A: effective-cap formula = ceiling(n_burned * ratio); n_selected <= eff_cap", {
  grid <- expand.grid(
    n_burned    = c(0L, 1L, 7L, 40L, 123L),
    ratio       = c(0, 0.25, 1.0, 2.0, Inf),
    n_available = c(0L, 3L, 50L, 1000L)
  )
  for (i in seq_len(nrow(grid))) {
    nb <- grid$n_burned[i]; rt <- grid$ratio[i]; av <- grid$n_available[i]
    eff <- effective_cap(nb, rt)
    sel <- n_selected_after_cap(nb, rt, av)

    if (is.finite(rt)) {
      expect_identical(eff, ceiling(nb * rt),
                       info = sprintf("nb=%d ratio=%g", nb, rt))
    }
    expect_lte(sel, av)
    if (is.finite(rt)) expect_lte(sel, eff)
    if (is.finite(rt) && av <= eff && eff > 0L) expect_equal(sel, av)
    if (is.finite(rt) && eff > 0L && av > eff) expect_equal(sel, eff)
    if (!is.finite(rt)) expect_equal(sel, av)
    if (rt == 0) expect_equal(sel, 0)
  }
})

test_that("PART A: the SHARED capping helper uses ceiling(n_burned*ratio); min(avail, cap)", {
  src <- paste(deparse(body(get(".of_cap_negative_buckets", envir = ns))),
               collapse = "\n")
  expect_true(grepl("ceiling(n_burned * cap_ratio)", src, fixed = TRUE))
  expect_true(grepl("n_cap_max >= length(avail)", src, fixed = TRUE))
  expect_true(grepl("sample.int(length(avail), n_cap_max)", src, fixed = TRUE))
})

test_that("PART A: BOTH the OOF and FINAL engines route capping through the SHARED resolver + helper", {
  oof_src   <- paste(deparse(body(get("run_oof_xgb", envir = ns))),
                     collapse = "\n")
  final_src <- paste(deparse(body(get("train_final_model_direct", envir = ns))),
                     collapse = "\n")
  for (s in list(oof_src, final_src)) {
    expect_true(grepl(".of_cap_negative_buckets", s, fixed = TRUE))
    expect_true(grepl(".of_resolve_supervised_eligibility", s, fixed = TRUE))
  }
  src <- paste(deparse(body(get(".of_cap_negative_buckets", envir = ns))),
               collapse = "\n")
  expect_true(grepl("ceiling(n_burned * cap_ratio)", src, fixed = TRUE))
  expect_true(grepl("ceiling(n_burned * random_to_burned_ratio)",
                    final_src, fixed = TRUE))
  # No spectral / contextual negative machinery survives in either engine.
  expect_false(grepl("spectral_hard_negative", oof_src, fixed = TRUE))
  expect_false(grepl("spectral_hard_negative", final_src, fixed = TRUE))
  expect_false(grepl("contextual_exclusion_to_burned_ratio", oof_src, fixed = TRUE))
  expect_false(grepl("contextual_exclusion_to_burned_ratio", final_src, fixed = TRUE))
})

# =============================================================================
# (1)->(2)/(3) cfg-STORED cap == cap RECEIVED by OOF and FINAL, via the single
#     shim resolver, with NO residual default.
# =============================================================================

test_that("PART A: each cfg cap resolves (override=NULL) to itself -- no residual default", {
  resolve <- get(".of_resolve_methodological_shim", envir = ns)

  for (cfg in list(.caps_cfg(),                                   # all default
                   # GATE 6.7: caps via the typed negative_pool_params block.
                   .caps_cfg(negative_pool_params =
                               list(caps = c(random = 0.75, otsu = 1.5))))) {  # all user
    caps <- cfg$train_control$caps
    prov <- cfg$resolved_params_provenance$train_control %||% list()
    pairs <- list(
      list("cap_random",     caps$random),
      list("cap_otsu",       caps$otsu)
    )
    for (p in pairs) {
      key <- p[[1]]; cfg_val <- p[[2]]
      out <- suppressWarnings(resolve(
        override = NULL, cfg_value = cfg_val, param = key,
        cfg_provenance = prov[[key]] %||% "default", record = NULL))
      expect_equal(out, cfg_val, info = key)
    }
  }
})

test_that("PART A: a conflicting explicit cap (cfg-user vs function override) ERRORS, never silently picks", {
  resolve <- get(".of_resolve_methodological_shim", envir = ns)
  expect_error(
    resolve(override = 1.0, cfg_value = 2.0, param = "cap_random",
            cfg_provenance = "user", record = NULL),
    regexp = "Conflicting supervised parameter 'cap_random'"
  )
  expect_silent(
    out <- resolve(override = 2.0, cfg_value = 2.0, param = "cap_random",
                   cfg_provenance = "user", record = NULL))
  expect_equal(out, 2.0)
})

# =============================================================================
# GATE 6.5: only two buckets exist; cap_spectral / cap_contextual are no longer
# builder args.
# =============================================================================

test_that("PART A / GATE 6.5: only two negative buckets + two cfg caps", {
  expect_equal(get(".of_valid_negative_buckets", envir = ns)(),
               c("random", "otsu"))
  cfg <- .caps_cfg()
  expect_setequal(names(cfg$train_control$caps), c("random", "otsu"))
})

test_that("PART A / GATE 6.5+6.7: passing cap_spectral / cap_contextual / cap_random / cap_otsu to the builder ERRORS (unused argument)", {
  expect_error(.caps_cfg(cap_spectral = 2.0),
               regexp = "cap_spectral|unused argument")
  expect_error(.caps_cfg(cap_contextual = 2.0),
               regexp = "cap_contextual|unused argument")
  # GATE 6.7: top-level cap_random / cap_otsu were REMOVED -> single-source caps
  # live ONLY in negative_pool_params$caps.
  expect_error(.caps_cfg(cap_random = 1.0),
               regexp = "cap_random|unused argument")
  expect_error(.caps_cfg(cap_otsu = 1.0),
               regexp = "cap_otsu|unused argument")
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

  # Supply every required methodological arg EXCEPT the random cap -> must
  # error (the required-arg guard refuses a silent default for the bucket).
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = tempfile(), prefix = "x", overwrite = TRUE, verbose = FALSE,
      # random_to_burned_ratio DROPPED on purpose.
      otsu_unburned_to_burned_ratio = 1,
      sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
      nrounds_max = 80, early_stopping_rounds = 80, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      model_params_base = get(".of_canonical_model_params", envir = ns)()
    )),
    regexp = "random_to_burned_ratio"
  )
})

test_that("PART A: dropping a required cap ERRORS on the OOF path (no silent 1.0 revert)", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fn <- get("run_oof_xgb", envir = ns)
  # The two cap ratios are REQUIRED formals; a dropped one ERRORS rather than
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
      # random_to_burned_ratio DROPPED on purpose.
      otsu_unburned_to_burned_ratio = 1
    ))),
    regexp = "random_to_burned_ratio"
  )
})
