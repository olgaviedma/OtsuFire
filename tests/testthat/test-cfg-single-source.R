# Gate 1B (2026-06-07): the cfg object is the SINGLE SOURCE OF TRUTH for the
# supervised methodological / training-control parameters. These tests assert:
#   (1) cfg$model_params + cfg$train_control carry the canonical schema/values,
#       sourced from the single canonical builders;
#   (2) the params that REACH the OOF and FINAL engines equal cfg
#       (recipe/param equality, not just non-error);
#   (3) a USER override in the builder (spectral cap = 2.0, nrounds_max = 10)
#       propagates to BOTH the OOF and FINAL engines IDENTICALLY;
#   (4) an internal-pure function ERRORS when a required methodological arg is
#       omitted (proves there is no silent default).

ns <- asNamespace("OtsuFire")

# Gate 1D.7: these tiny cfg fixtures are shared with the other supervised
# contract suites via helper-supervised-contracts.R (sc_* fixtures). Kept here
# as thin aliases so the assertions below are untouched.
mk_ss_tif  <- sc_change_index_tif
mk_ss_gpkg <- sc_decisions_gpkg
mk_ss_cfg  <- sc_min_cfg

# ---------------------------------------------------------------------------
# (1) Schema + canonical-equality of the two resolved cfg sections.
# ---------------------------------------------------------------------------
test_that("Gate 1B: cfg carries model_params + train_control with the canonical schema", {
  cfg <- mk_ss_cfg()

  # model_params = canonical xgb block MINUS scale_pos_weight.
  expect_type(cfg$model_params, "list")
  expect_false("scale_pos_weight" %in% names(cfg$model_params))
  expect_identical(cfg$model_params,
                   get(".of_canonical_model_params", envir = ns)())
  # Sanity: it equals the engine's canonical block with spw stripped.
  blk <- get(".of_canonical_xgb_params", envir = ns)(scale_pos_weight = 1)
  blk[["scale_pos_weight"]] <- NULL
  expect_identical(cfg$model_params, blk)

  # train_control schema.
  tc <- cfg$train_control
  expect_type(tc, "list")
  expect_setequal(
    names(tc),
    c("nrounds_max", "early_stop", "seeds", "val_frac", "group_col",
      "impute_numeric", "impute_factor_missing", "caps",
      "feature_whitelist_override", "feature_weights",
      "training_protocol", "oof_sampling")
  )
  expect_setequal(names(tc$seeds),
                  c("oof_seed_base", "final_sampling_seed", "final_seed"))
  expect_setequal(names(tc$caps),
                  c("contextual", "spectral", "random", "otsu"))
  # Canonical values (the INVARIANT: identical to the historical defaults).
  expect_identical(tc, get(".of_canonical_train_control", envir = ns)())
  expect_equal(tc$nrounds_max, 4000L)
  expect_equal(tc$early_stop, 80L)
  expect_equal(unlist(tc$seeds, use.names = FALSE), c(42L, 42L, 42L))
  expect_equal(tc$val_frac, 0.15)
  expect_equal(tc$group_col, "block_id")
  expect_equal(unname(unlist(tc$caps)), c(0.25, 1.0, 1.0, 1.0))
  # Fixed internal constants (traceability only; not user-settable).
  expect_equal(tc$training_protocol, "nested_refit")
  # oof_sampling is now a FIXED constant, not a user-settable toggle.
  expect_equal(tc$oof_sampling, "capped")
})

# ---------------------------------------------------------------------------
# (2) The params that REACH the OOF + FINAL engines equal cfg.
#     We mock the internal engines (run_dm_oof_pipeline / train_final_model_direct)
#     to CAPTURE the args the public wrappers hand them, and assert equality
#     against cfg$train_control / cfg$model_params.
# ---------------------------------------------------------------------------
capture_oof_args <- function(cfg, train_feats) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    run_dm_oof_pipeline = function(...) {
      a <- list(...)
      captured$nrounds_max <- a$nrounds_max
      captured$early_stop  <- a$early_stop
      captured$seed_base   <- a$seed_base
      captured$group_col   <- a$group_col
      captured$caps <- c(a$contextual_exclusion_to_burned_ratio,
                         a$spectral_hard_negative_to_burned_ratio,
                         a$random_to_burned_ratio,
                         a$otsu_unburned_to_burned_ratio)
      captured$params <- a$params
      list(dm = NULL, oof = list(oof_agg = NULL, oof_long = NULL),
           files = list())
    },
    .package = "OtsuFire"
  )
  run_oof_diagnostics(train_features = train_feats,
                      scoring_features = train_feats, config = cfg,
                      out_dir = tempfile(), matrix_dir = tempfile())
  captured
}

capture_final_args <- function(cfg, train_gpkg) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    train_final_model_direct = function(...) {
      a <- list(...)
      captured$nrounds_max           <- a$nrounds_max
      captured$early_stopping_rounds <- a$early_stopping_rounds
      captured$sampling_seed         <- a$sampling_seed
      captured$seed                  <- a$seed
      captured$val_frac              <- a$val_frac
      captured$group_col             <- a$group_col
      captured$impute_numeric        <- a$impute_numeric
      captured$impute_factor_missing <- a$impute_factor_missing
      captured$caps <- c(a$contextual_exclusion_to_burned_ratio,
                         a$spectral_hard_negative_to_burned_ratio,
                         a$random_to_burned_ratio,
                         a$otsu_unburned_to_burned_ratio)
      captured$model_params_base <- a$model_params_base
      list(model = NULL, recipe = NULL, training_ok_sf = NULL, split = NULL,
           files = list())
    },
    .package = "OtsuFire"
  )
  train_final_burned_model(train_features = train_gpkg, config = cfg,
                           out_dir = tempfile())
  captured
}

# Gate 1D.7: shared with the caps suite via helper-supervised-contracts.R.
mk_oof_train_feats  <- sc_burned_train_sf
mk_final_train_gpkg <- sc_burned_train_gpkg

test_that("Gate 1B: params reaching the OOF engine equal cfg$train_control / cfg$model_params", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_ss_cfg()
  cap <- capture_oof_args(cfg, mk_oof_train_feats())

  expect_equal(cap$nrounds_max, cfg$train_control$nrounds_max)
  expect_equal(cap$early_stop,  cfg$train_control$early_stop)
  expect_equal(cap$seed_base,   cfg$train_control$seeds$oof_seed_base)
  # Gate 1B residual: OOF's group_col is sourced from cfg$train_control (same
  # single source the FINAL stage uses), not a hardcoded literal.
  expect_equal(cap$group_col,   cfg$train_control$group_col)
  expect_equal(cap$caps, with(cfg$train_control$caps,
                              c(contextual, spectral, random, otsu)))
  # The xgb params built for OOF are cfg$model_params + a computed spw.
  expect_identical(cap$params[setdiff(names(cap$params), "scale_pos_weight")],
                   cfg$model_params)
})

test_that("Gate 1B: params reaching the FINAL engine equal cfg$train_control / cfg$model_params", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_ss_cfg()
  cap <- capture_final_args(cfg, mk_final_train_gpkg())

  expect_equal(cap$nrounds_max,           cfg$train_control$nrounds_max)
  expect_equal(cap$early_stopping_rounds, cfg$train_control$early_stop)
  expect_equal(cap$sampling_seed,         cfg$train_control$seeds$final_sampling_seed)
  expect_equal(cap$seed,                  cfg$train_control$seeds$final_seed)
  expect_equal(cap$val_frac,              cfg$train_control$val_frac)
  expect_equal(cap$group_col,             cfg$train_control$group_col)
  expect_equal(cap$impute_numeric,        cfg$train_control$impute_numeric)
  expect_equal(cap$impute_factor_missing, cfg$train_control$impute_factor_missing)
  expect_equal(cap$caps, with(cfg$train_control$caps,
                              c(contextual, spectral, random, otsu)))
  # The FINAL engine receives cfg$model_params as its model_params_base.
  expect_identical(cap$model_params_base, cfg$model_params)
})

# ---------------------------------------------------------------------------
# (3) A builder override (spectral cap = 2.0, nrounds_max = 10) propagates to
#     BOTH OOF and FINAL identically.
# ---------------------------------------------------------------------------
test_that("Gate 1B: a builder override propagates to BOTH OOF and FINAL identically", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_ss_cfg(cap_spectral = 2.0, nrounds_max = 10L)
  expect_equal(cfg$train_control$caps$spectral, 2.0)
  expect_equal(cfg$train_control$nrounds_max, 10L)

  cap_oof   <- capture_oof_args(cfg, mk_oof_train_feats())
  cap_final <- capture_final_args(cfg, mk_final_train_gpkg())

  # spectral cap == 2.0 in BOTH.
  expect_equal(cap_oof$caps[2],   2.0)
  expect_equal(cap_final$caps[2], 2.0)
  # nrounds_max == 10 in BOTH.
  expect_equal(cap_oof$nrounds_max,   10L)
  expect_equal(cap_final$nrounds_max, 10L)
  # The full cap vectors are identical across the two stages.
  expect_equal(cap_oof$caps, cap_final$caps)
})

# ---------------------------------------------------------------------------
# (3b) STRONG propagation test: values DELIBERATELY different from every
#      canonical default reach BOTH engines WITHOUT substitution.
# ---------------------------------------------------------------------------
test_that("Gate 1B: deliberately non-default values propagate to BOTH OOF and FINAL without substitution", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_ss_cfg(
    nrounds_max         = 137L,
    early_stop          = 9L,
    val_frac            = 0.222,
    cap_spectral        = 2.0,
    oof_seed_base       = 4242L,
    final_sampling_seed = 1234L,
    final_seed          = 9999L,
    # Gate 1B residual: a DELIBERATELY non-default group_col must reach BOTH the
    # OOF and FINAL engines from the single cfg$train_control (was hardcoded
    # "block_id" on the OOF path).
    group_col           = "custom_block_col"
  )
  # The cfg stores exactly these (no canonical fallback).
  expect_equal(cfg$train_control$nrounds_max, 137L)
  expect_equal(cfg$train_control$early_stop, 9L)
  expect_equal(cfg$train_control$val_frac, 0.222)
  expect_equal(cfg$train_control$caps$spectral, 2.0)
  expect_equal(cfg$train_control$seeds$oof_seed_base, 4242L)
  expect_equal(cfg$train_control$seeds$final_sampling_seed, 1234L)
  expect_equal(cfg$train_control$seeds$final_seed, 9999L)
  expect_equal(cfg$train_control$group_col, "custom_block_col")

  cap_oof   <- capture_oof_args(cfg, mk_oof_train_feats())
  cap_final <- capture_final_args(cfg, mk_final_train_gpkg())

  # OOF engine: the OOF-specific seed + shared nrounds/early-stop + spectral cap.
  expect_equal(cap_oof$nrounds_max, 137L)
  expect_equal(cap_oof$early_stop, 9L)
  expect_equal(cap_oof$seed_base, 4242L)
  expect_equal(cap_oof$caps[2], 2.0)
  # The non-default group_col reaches the OOF engine (no hardcoded "block_id").
  expect_equal(cap_oof$group_col, "custom_block_col")

  # FINAL engine: FINAL-specific seeds + val_frac + nrounds/early-stop + cap.
  expect_equal(cap_final$nrounds_max, 137L)
  expect_equal(cap_final$early_stopping_rounds, 9L)
  expect_equal(cap_final$val_frac, 0.222)
  expect_equal(cap_final$sampling_seed, 1234L)
  expect_equal(cap_final$seed, 9999L)
  expect_equal(cap_final$caps[2], 2.0)
  # The non-default group_col reaches the FINAL engine too.
  expect_equal(cap_final$group_col, "custom_block_col")

  # The shared knobs (nrounds, early-stop, the cap vector, group_col) are
  # IDENTICAL across the two stages -- proof of a single source of truth.
  expect_equal(cap_oof$nrounds_max, cap_final$nrounds_max)
  expect_equal(cap_oof$caps, cap_final$caps)
  expect_equal(cap_oof$group_col, cap_final$group_col)
})

# ---------------------------------------------------------------------------
# (3c) LEGACY-RECIPE EQUALITY: the CURRENT defaults resolve to EXACTLY the same
#      xgb params / training controls as the historical pre-refactor numbers.
#      This guards against the refactor silently changing a canonical value.
# ---------------------------------------------------------------------------
test_that("Gate 1B: current defaults equal the documented legacy recipe (params + controls)", {
  cfg <- mk_ss_cfg()

  # (a) xgb block: cfg$model_params + a spw == the historical
  #     .of_canonical_xgb_params(spw) block, field-for-field.
  spw <- 3.7  # any finite spw; the block must match except for this field.
  legacy_block <- get(".of_canonical_xgb_params", envir = ns)(scale_pos_weight = spw)
  cfg_block <- cfg$model_params
  cfg_block[["scale_pos_weight"]] <- spw
  # Same names, same order, same values.
  expect_identical(names(cfg_block), names(legacy_block))
  expect_identical(cfg_block, legacy_block)

  # (b) the documented legacy NUMBERS (FRENTE 1 / D1 unification + B1 defaults).
  expect_identical(legacy_block$booster, "gbtree")
  expect_identical(legacy_block$objective, "binary:logistic")
  expect_identical(legacy_block$eval_metric, c("logloss", "aucpr"))
  expect_equal(legacy_block$eta, 0.05)
  expect_equal(legacy_block$max_depth, 5)
  expect_equal(legacy_block$min_child_weight, 5)
  expect_equal(legacy_block$subsample, 0.8)
  expect_equal(legacy_block$colsample_bytree, 0.75)
  expect_equal(legacy_block$gamma, 0)
  expect_equal(legacy_block$lambda, 1)
  expect_equal(legacy_block$alpha, 0)

  # (c) training controls == the documented legacy numbers.
  tc <- cfg$train_control
  expect_equal(tc$nrounds_max, 4000L)
  expect_equal(tc$early_stop, 80L)
  expect_equal(tc$seeds$oof_seed_base, 42L)
  expect_equal(tc$seeds$final_sampling_seed, 42L)
  expect_equal(tc$seeds$final_seed, 42L)
  expect_equal(tc$val_frac, 0.15)
  expect_equal(tc$group_col, "block_id")
  expect_equal(tc$impute_numeric, "median")
  expect_equal(tc$impute_factor_missing, "MISSING")
  expect_equal(tc$caps$contextual, 0.25)
  expect_equal(tc$caps$spectral, 1.0)
  expect_equal(tc$caps$random, 1.0)
  expect_equal(tc$caps$otsu, 1.0)
  expect_equal(tc$training_protocol, "nested_refit")
  # oof_sampling is now a FIXED constant, not a user-settable toggle.
  expect_equal(tc$oof_sampling, "capped")
})

# ---------------------------------------------------------------------------
# (4) An internal-pure function ERRORS when a required methodological arg is
#     omitted (proves there is no silent revert to a hardcoded default).
# ---------------------------------------------------------------------------
test_that("Gate 1B: train_final_model_direct ERRORS when a required methodological arg is omitted", {
  fn <- get("train_final_model_direct", envir = ns)
  g <- mk_final_train_gpkg()
  # Supply everything EXCEPT nrounds_max -> must error (no silent default).
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = tempfile(), prefix = "x", overwrite = TRUE, verbose = FALSE,
      contextual_exclusion_to_burned_ratio = 1,
      spectral_hard_negative_to_burned_ratio = 1,
      random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
      sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
      early_stopping_rounds = 80, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      model_params_base = get(".of_canonical_model_params", envir = ns)()
      # nrounds_max DROPPED on purpose.
    )),
    regexp = "required resolved arg 'nrounds_max'"
  )
})

test_that("Gate 1B: train_final_model_direct ERRORS when model_params_base is omitted", {
  fn <- get("train_final_model_direct", envir = ns)
  g <- mk_final_train_gpkg()
  expect_error(
    suppressMessages(fn(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = tempfile(), prefix = "x", overwrite = TRUE, verbose = FALSE,
      contextual_exclusion_to_burned_ratio = 1,
      spectral_hard_negative_to_burned_ratio = 1,
      random_to_burned_ratio = 1, otsu_unburned_to_burned_ratio = 1,
      sampling_seed = 42, seed = 42, val_frac = 0.15, group_col = "block_id",
      nrounds_max = 80, early_stopping_rounds = 80, impute_numeric = "median",
      impute_factor_missing = "MISSING"
      # model_params_base DROPPED on purpose.
    )),
    regexp = "required resolved arg 'model_params_base'"
  )
})

test_that("Gate 1B: run_dm_oof_pipeline ERRORS when nrounds_max is omitted (no silent default)", {
  fn <- get("run_dm_oof_pipeline", envir = ns)
  expect_error(
    suppressMessages(fn(
      labelled = data.frame(), burned_like = data.frame(),
      labelled_df = data.frame(), params = list(), result_dir = tempdir(),
      early_stop = 80L, seed_base = 42L, group_col = "block_id"
      # nrounds_max DROPPED on purpose.
    )),
    regexp = "required resolved arg 'nrounds_max'"
  )
})

test_that("Gate 1B: run_dm_oof_pipeline ERRORS when group_col is omitted (no silent default)", {
  fn <- get("run_dm_oof_pipeline", envir = ns)
  expect_error(
    suppressMessages(fn(
      labelled = data.frame(), burned_like = data.frame(),
      labelled_df = data.frame(), params = list(), result_dir = tempdir(),
      nrounds_max = 4000L, early_stop = 80L, seed_base = 42L
      # group_col DROPPED on purpose -> no silent revert to a hardcoded literal.
    )),
    regexp = "required resolved arg 'group_col'"
  )
})

test_that("Gate 1B: run_oof_xgb ERRORS when group_col is omitted (no silent default)", {
  fn <- get("run_oof_xgb", envir = ns)
  expect_error(
    suppressMessages(fn(
      XL_mat = Matrix::Matrix(matrix(0, 1, 1), sparse = TRUE), y = 0L,
      labelled_df = data.frame(fire_uid = "a", class = "burned"),
      params = list(), nrounds_max = 4000L, early_stop = 80L, seed_base = 42L
      # group_col DROPPED on purpose.
    )),
    regexp = "required resolved arg 'group_col'"
  )
})

# ---------------------------------------------------------------------------
# (5) Builder validation of the new override arguments.
# ---------------------------------------------------------------------------
test_that("Gate 1B: builder validates the resolved-param overrides", {
  expect_error(mk_ss_cfg(val_frac = 0), regexp = "val_frac")
  expect_error(mk_ss_cfg(val_frac = 1), regexp = "val_frac")
  expect_error(mk_ss_cfg(nrounds_max = 0), regexp = "nrounds_max")
  expect_error(mk_ss_cfg(early_stop = -1), regexp = "early_stop")
  expect_error(mk_ss_cfg(cap_spectral = -0.1), regexp = "cap_spectral")
  # training_protocol is no longer a builder argument -> R "unused argument".
  expect_error(mk_ss_cfg(training_protocol = "nested_refit"),
               regexp = "unused argument")
  # oof_sampling is no longer a builder argument -> R "unused argument".
  expect_error(mk_ss_cfg(oof_sampling = "capped"), regexp = "unused argument")
  expect_error(mk_ss_cfg(model_params = list(scale_pos_weight = 3)),
               regexp = "scale_pos_weight")
})

# ---------------------------------------------------------------------------
# (6) Gate 1B (2026-06-07): PER-FIELD provenance for cfg$model_params.
#     The builder folds a PARTIAL model_params override onto the canonical block,
#     so provenance is tracked PER xgb field (canonical / requested / resolved /
#     provenance), mirroring train_control. A partial override of eta marks ONLY
#     eta "user"; every other field stays "default" and canonical.
# ---------------------------------------------------------------------------
test_that("Gate 1B: a PARTIAL model_params override yields per-field provenance (eta=user, rest default)", {
  cfg <- mk_ss_cfg(model_params = list(eta = 0.03))

  # The resolved block keeps eta = 0.03 and every other field canonical.
  canon <- get(".of_canonical_model_params", envir = ns)()
  expect_equal(cfg$model_params$eta, 0.03)
  for (nm in setdiff(names(canon), "eta")) {
    expect_identical(cfg$model_params[[nm]], canon[[nm]],
                     info = paste("field", nm, "must stay canonical"))
  }

  # Per-field provenance structure mirrors train_control: keyed by xgb field,
  # each a list(canonical, requested, resolved, provenance).
  mp <- cfg$resolved_params_provenance$model_params
  expect_type(mp, "list")
  expect_true(all(names(canon) %in% names(mp)))
  expect_false("scale_pos_weight" %in% names(mp))  # site-specific, excluded.

  # eta = user, with requested 0.03 and resolved 0.03; canonical recorded.
  expect_equal(mp$eta$provenance, "user")
  expect_equal(mp$eta$requested, "0.03")
  expect_equal(mp$eta$resolved, "0.03")
  expect_equal(mp$eta$canonical, "0.05")

  # Every OTHER field: default, requested = NA, resolved == canonical.
  for (nm in setdiff(names(canon), "eta")) {
    expect_equal(mp[[nm]]$provenance, "default",
                 info = paste("field", nm))
    expect_true(is.na(mp[[nm]]$requested), info = paste("field", nm))
    expect_equal(mp[[nm]]$resolved, mp[[nm]]$canonical,
                 info = paste("field", nm))
  }
})

test_that("Gate 1B: the all-default model_params case marks EVERY field 'default'", {
  cfg <- mk_ss_cfg()  # no model_params override.
  mp <- cfg$resolved_params_provenance$model_params
  canon <- get(".of_canonical_model_params", envir = ns)()
  expect_setequal(names(mp), names(canon))
  for (nm in names(mp)) {
    expect_equal(mp[[nm]]$provenance, "default", info = paste("field", nm))
    expect_true(is.na(mp[[nm]]$requested), info = paste("field", nm))
    expect_equal(mp[[nm]]$resolved, mp[[nm]]$canonical, info = paste("field", nm))
  }
})

test_that("Gate 1B: a multi-field partial override tags exactly the supplied fields 'user'", {
  cfg <- mk_ss_cfg(model_params = list(eta = 0.03, max_depth = 8))
  mp <- cfg$resolved_params_provenance$model_params
  canon <- get(".of_canonical_model_params", envir = ns)()

  expect_equal(mp$eta$provenance, "user")
  expect_equal(mp$max_depth$provenance, "user")
  expect_equal(mp$max_depth$requested, "8")
  expect_equal(cfg$model_params$max_depth, 8)
  # subsample and the rest stay default.
  for (nm in setdiff(names(canon), c("eta", "max_depth"))) {
    expect_equal(mp[[nm]]$provenance, "default", info = paste("field", nm))
  }
})

test_that("Gate 1B: a PARTIAL model_params override (eta=0.03) reaches BOTH OOF and FINAL with the rest canonical", {
  skip_if_not_installed("sf")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))
  cfg <- mk_ss_cfg(model_params = list(eta = 0.03))
  canon <- get(".of_canonical_model_params", envir = ns)()

  cap_oof   <- capture_oof_args(cfg, mk_oof_train_feats())
  cap_final <- capture_final_args(cfg, mk_final_train_gpkg())

  # OOF: cfg$model_params + computed spw. eta = 0.03, rest canonical.
  oof_block <- cap_oof$params
  expect_equal(oof_block$eta, 0.03)
  for (nm in setdiff(names(canon), "eta")) {
    expect_identical(oof_block[[nm]], canon[[nm]], info = paste("OOF field", nm))
  }
  # FINAL: model_params_base == cfg$model_params (no spw at this boundary).
  expect_identical(cap_final$model_params_base, cfg$model_params)
  expect_equal(cap_final$model_params_base$eta, 0.03)
  for (nm in setdiff(names(canon), "eta")) {
    expect_identical(cap_final$model_params_base[[nm]], canon[[nm]],
                     info = paste("FINAL field", nm))
  }
})

test_that("Gate 1B: the run-summary surfaces the per-field model_params provenance table", {
  cfg <- mk_ss_cfg(model_params = list(eta = 0.03))
  # The render helper flattens cfg$resolved_params_provenance$model_params into
  # the per-field table a downstream manifest documents.
  df <- get(".of_model_params_provenance_to_df", envir = ns)(
    cfg$resolved_params_provenance$model_params)
  expect_s3_class(df, "data.frame")
  expect_setequal(names(df),
                  c("param", "canonical", "requested", "resolved", "provenance"))
  eta_row <- df[df$param == "eta", , drop = FALSE]
  expect_equal(nrow(eta_row), 1L)
  expect_equal(eta_row$provenance, "user")
  expect_equal(eta_row$requested, "0.03")
  expect_equal(eta_row$resolved, "0.03")
  expect_equal(eta_row$canonical, "0.05")
  # max_depth is a default row with NA requested.
  md_row <- df[df$param == "max_depth", , drop = FALSE]
  expect_equal(md_row$provenance, "default")
  expect_true(is.na(md_row$requested))
})
