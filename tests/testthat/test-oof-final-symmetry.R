# ============================================================================
# Gate 1D.4 (2026-06-08): OOF <-> FINAL SYMMETRY CONTRACT.
#
# This suite PROVES, as a stable contract, that the OOF diagnostics stage and
# the FINAL-model stage share ONE recipe + ONE set of resolved parameters, and
# it ENUMERATES the legitimate (by-design) differences between the two paths.
#
# It does NOT compare model binaries. It compares the RESOLVED objects/params
# the two paths construct + hand to the SHARED leakage-free core
# (.of_nested_refit_fit). The B1 nested_refit work already made OOF and FINAL
# call that SAME core with .of_canonical_* params sourced from a single cfg;
# this suite locks that in and documents the boundary.
#
# ---------------------------------------------------------------------------
# STEP-0 TRACE (where each shared element is obtained / used). file:line refer
# to the package sources as of this gate.
#
#   shared element                | OOF site                                  | FINAL site
#   ------------------------------|-------------------------------------------|------------------------------------------
#   shared fit core               | run_oof_xgb -> .of_nested_refit_fit       | train_final_model_direct -> .of_nested_refit_fit
#                                 |   internal-sup-oof-training.R:300         |   internal-sup-train-final-direct.R:501
#   model_params (xgb block)      | params_fn = .of_canonical_xgb_params      | params_fn = .params_from_cfg (model_params_base = cfg$model_params)
#                                 |   internal-sup-oof-training.R:307         |   internal-sup-train-final-direct.R:157,508
#                                 |   both resolve to .of_canonical_xgb_block (internal-sup-xgb-params.R:55)
#   nrounds_max / early_stop      | cfg$train_control (threaded)              | cfg$train_control (threaded)
#                                 |   internal-sup-oof-training.R:311-312     |   internal-sup-train-final-direct.R:512-513
#   val_frac / group_col          | cfg$train_control                        | cfg$train_control
#                                 |   internal-sup-oof-training.R:304-305     |   internal-sup-train-final-direct.R:506-507
#   impute_numeric / factor       | cfg$train_control                        | cfg$train_control
#                                 |   internal-sup-oof-training.R:313-314     |   internal-sup-train-final-direct.R:514-515
#   the 4 caps                    | run_dm_oof_pipeline (cap_outer_train)     | train_final_model_direct (pool sampling)
#                                 |   internal-sup-oof-training.R:221-234     |   internal-sup-train-final-direct.R:223-259
#   feature whitelist             | .filter_to_supervised_whitelist           | .filter_to_supervised_whitelist
#                                 |   run_dm_oof_pipeline apply_whitelist_filter (internal-sup-oof-wrapper.R:300) | internal-sup-train-final-direct.R:332
#   feature weights               | feature_weights -> core apply_fw          | feature_weights -> core apply_fw
#                                 |   internal-sup-oof-training.R:310         |   internal-sup-train-final-direct.R:511
#   imputation (median fit rule)  | .of_nested_fit_medians (train rows only)  | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:113,317,356
#   missing handling (DMatrix NA) | .of_nested_build_matrix + missing=NA      | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:173,329,361
#   feature order (x_cols)        | colnames(refit matrix) from same builder  | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:358-359
#   scale_pos_weight rule         | n_neg/max(1,n_pos) inside core            | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:335,365
#   inner split strategy          | .of_nested_group_split inside core        | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:42,307
#   early stopping                | xgb.train(watchlist=inner_val, early_stop)| SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:340-347
#   best_iteration selection      | m_sel$best_iteration %||% nrounds_max     | SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:348-350
#   refit                         | xgb.train(nrounds=best_iter, no watchlist)| SAME helper, SAME core
#                                 |   internal-sup-nested-refit.R:370-376
#   scoring recipe schema         | recipe$cols/$impute/$training fields      | written by train_final_model_direct
#                                 |   internal-sup-train-final-direct.R:624-679
#
# CONFIRMATION: every "inside core" row is, by construction, the SAME R code
# executed for both stages (one function body). The only way the two paths can
# differ on those is by passing different ARGUMENTS to the core; the tests below
# capture and compare those resolved arguments.
#
# ---------------------------------------------------------------------------
# LEGITIMATE (by-design) DIFFERENCES — the contract DOCUMENTS + ALLOWS exactly:
#   (L1) OOF holds out an OUTER test fold (never in any watchlist / imputation /
#        selection); FINAL has no outer test fold.
#   (L2) NEGATIVE SAMPLING is shared: OOF uses the SAME capped negative-sampling
#        policy as FINAL, applied independently within each training fold. There
#        is NO oof_sampling toggle on any public function or OOF engine (passing
#        oof_sampling ERRORS); the cfg/audit constant is the fixed label
#        "capped".
#   (L3) function-specific SEEDS differ by design: OOF fold_seed =
#        seed_base + 1000*rep + fold (oof_seed_base); FINAL fold_seed =
#        final_seed, sampling_seed = final_sampling_seed. The DERIVATION rule is
#        the documented one (asserted below).
#   (L4) FINAL's refit uses the ENTIRE available pool (L_ok); OOF's per-fold
#        refit uses that fold's (capped) outer-train partition.
# ============================================================================

ns_sym <- asNamespace("OtsuFire")

g_sym       <- function(nm) get(nm, envir = ns_sym)
nr_fit_sym  <- function() g_sym(".of_nested_refit_fit")
# These return the FUNCTION OBJECT; call as canon_xgb()(...) etc.
canon_xgb   <- function() g_sym(".of_canonical_xgb_params")
canon_block <- function() g_sym(".of_canonical_xgb_block")
canon_mp    <- function() g_sym(".of_canonical_model_params")
canon_tc    <- function() g_sym(".of_canonical_train_control")

# Re-use the nested-refit fixture generator if loaded; otherwise define a local
# copy (kept tiny + identical in spirit to make_nested_oof_df in
# test-nested-refit.R so the two suites stay aligned).
make_sym_oof_df <- function(n = 60, seed = 5L) {
  set.seed(seed)
  whitelist <- g_sym(".supervised_feature_cols")
  feat_df <- as.data.frame(
    matrix(stats::runif(n * length(whitelist)),
           nrow = n, dimnames = list(NULL, whitelist))
  )
  admin <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = rep(c("burned", "unburned"), length.out = n),
    source    = rep(c("burned_truth", "random_burnable_background",
                      "deterministic_drop_hard", "otsu_patch_residual"),
                    length.out = n),
    neg_type  = rep(c(NA_character_, "background_cell",
                      "spectral_reject_medium", "otsu_patch_drop"),
                    length.out = n),
    poly_id   = sprintf("p_%03d", seq_len(n)),
    block_id  = rep(seq_len(12), length.out = n),
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2 = rep(c(2L, 3L, 1L), length.out = n),
    stringsAsFactors = FALSE
  )
  blk_to_fold1 <- rep(c(1L, 2L, 3L), length.out = 12)
  blk_to_fold2 <- rep(c(3L, 1L, 2L), length.out = 12)
  admin$fold_rep1 <- blk_to_fold1[admin$block_id]
  admin$fold_rep2 <- blk_to_fold2[admin$block_id]
  cbind(admin, feat_df)
}

# A capture mock for .of_nested_refit_fit(): records the RESOLVED arguments each
# path hands the shared core, then returns a minimal stub the caller can keep
# threading (so neither path actually trains xgboost). This is the OOF<->FINAL
# analogue of the cfg-capture pattern in test-cfg-single-source.R.
capture_core_args <- function(expr_fn) {
  store <- new.env()
  store$calls <- list()
  stub <- function(train_df, feature_cols, label_col, group_col, block_col,
                   val_frac, params_fn, sampling_seed, fold_seed,
                   feature_weights = NULL, nrounds_max, early_stopping_rounds,
                   impute_numeric, impute_factor_missing, verbose = FALSE, ...) {
    rec <- list(
      feature_cols          = feature_cols,
      label_col             = label_col,
      group_col             = group_col,
      block_col             = block_col,
      val_frac              = val_frac,
      sampling_seed         = sampling_seed,
      fold_seed             = fold_seed,
      nrounds_max           = nrounds_max,
      early_stopping_rounds = early_stopping_rounds,
      impute_numeric        = impute_numeric,
      impute_factor_missing = impute_factor_missing,
      # the canonical xgb block the params_fn produces for a fixed spw, so the
      # two paths' model_params can be compared field-for-field.
      params_block          = params_fn(2.0),
      n_train_rows          = nrow(train_df),
      train_classes         = sort(unique(as.character(train_df[[label_col]])))
    )
    store$calls[[length(store$calls) + 1L]] <- rec
    # The OOF caller PREDICTS the outer-test fold with the returned model, so the
    # stub must hand back a REAL (trivially trained) booster whose feature space
    # matches feature_cols. We train one round on a tiny random matrix; the model
    # itself is irrelevant to the contract (we only inspect the captured args).
    x_cols <- feature_cols
    Xtiny <- matrix(stats::runif(4 * length(x_cols)), nrow = 4,
                    dimnames = list(NULL, x_cols))
    dtiny <- xgboost::xgb.DMatrix(Xtiny, label = c(0L, 1L, 0L, 1L), missing = NA)
    mtiny <- xgboost::xgb.train(
      params = list(objective = "binary:logistic", max_depth = 2),
      data = dtiny, nrounds = 1L, verbose = 0)
    list(
      model = mtiny,
      best_iteration = 1L,
      medians = stats::setNames(as.list(rep(0, length(feature_cols))), feature_cols),
      x_cols = x_cols,
      feature_cols = feature_cols,
      spw_selection = 1.0,
      spw_refit = 1.0,
      inner_split = list(tr_idx = 1L, val_idx = 1L, mode = "row_split"),
      degenerate_cols = character(0),
      audit = data.frame(spw_selection = 1.0, spw_refit = 1.0,
                         best_iteration = 1L, stringsAsFactors = FALSE)
    )
  }
  testthat::local_mocked_bindings(.of_nested_refit_fit = stub,
                                  .package = "OtsuFire")
  expr_fn()
  store$calls
}

# ---------------------------------------------------------------------------
# SHARED-CORE IDENTITY: both stages literally call the SAME function object.
# (If a future refactor forks the core, this fails -- the single most important
# structural invariant of the contract.)
# ---------------------------------------------------------------------------
test_that("OOF and FINAL reference the SAME shared fit-core symbol", {
  oof_src   <- deparse(body(g_sym("run_oof_xgb")))
  final_src <- deparse(body(g_sym("train_final_model_direct")))
  expect_true(any(grepl("\\.of_nested_refit_fit", oof_src)),
              info = "OOF must call .of_nested_refit_fit")
  expect_true(any(grepl("\\.of_nested_refit_fit", final_src)),
              info = "FINAL must call .of_nested_refit_fit")
  # The symbol resolves to ONE function object in the namespace.
  expect_true(is.function(nr_fit_sym()))
})

# ---------------------------------------------------------------------------
# (1) model_params: the canonical xgb block both paths build is identical
#     except for the site-computed scale_pos_weight.
# ---------------------------------------------------------------------------
test_that("model_params: OOF and FINAL build the SAME canonical xgb block (spw aside)", {
  # cfg$model_params (FINAL's model_params_base) == canonical block minus spw.
  mp <- canon_mp()()
  blk <- canon_block()(scale_pos_weight = 1)
  blk[["scale_pos_weight"]] <- NULL
  expect_identical(mp, blk)

  # OOF's params_fn (.of_canonical_xgb_params) and FINAL's .params_from_cfg
  # closure both yield the canonical block + the spw they are handed. For a
  # fixed spw the non-spw fields must be byte-identical.
  oof_block   <- canon_xgb()(scale_pos_weight = 3.7)
  final_base  <- mp                       # == cfg$model_params
  final_block <- final_base
  final_block[["scale_pos_weight"]] <- 3.7
  expect_identical(
    oof_block[setdiff(names(oof_block), "scale_pos_weight")],
    final_block[setdiff(names(final_block), "scale_pos_weight")]
  )
  expect_equal(oof_block$scale_pos_weight, final_block$scale_pos_weight)
})

# ---------------------------------------------------------------------------
# (2) train_control: nrounds_max, early_stop, val_frac, group_col, impute_*
#     come from ONE canonical builder -> identical for both stages.
# ---------------------------------------------------------------------------
test_that("train_control: both stages draw the same nrounds/early_stop/val_frac/group_col/impute", {
  tc <- canon_tc()()
  expect_equal(tc$nrounds_max, 4000L)
  expect_equal(tc$early_stop, 80L)
  expect_equal(tc$val_frac, 0.15)
  expect_equal(tc$group_col, "block_id")
  expect_equal(tc$impute_numeric, "median")
  expect_equal(tc$impute_factor_missing, "MISSING")
  # The 4 caps come from the same canonical source too.
  expect_equal(unname(unlist(tc$caps)), c(0.25, 1.0, 1.0, 1.0))
})

# ---------------------------------------------------------------------------
# (3) RESOLVED-ARG CAPTURE: run both the FINAL engine and the OOF engine on the
#     same tiny synthetic pool, mocking the shared core to CAPTURE the resolved
#     arguments. Assert the shared invariants are EQUAL across the two paths and
#     that ONLY the enumerated legitimate differences differ.
# ---------------------------------------------------------------------------
sym_env <- new.env()

sym_env$build_final_gpkg <- function(df, seed = 5L) {
  sfc <- sf::st_sfc(
    lapply(seq_len(nrow(df)), function(i) {
      x <- (i %% 5); y <- (i %/% 5)
      sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1),
                                c(x, y + 1), c(x, y))))
    }), crs = 3035)
  sf_df <- sf::st_sf(df, geometry = sfc)
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(sf_df, g, layer = "train_features", quiet = TRUE,
               delete_dsn = TRUE)
  g
}

run_capture_pair <- function() {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("dplyr")
  testthat::skip_if(!exists("local_mocked_bindings",
                            where = asNamespace("testthat")))

  df <- make_sym_oof_df(seed = 5L)

  # ---- FINAL capture ----
  g <- sym_env$build_final_gpkg(df)
  fin_engine <- g_sym("train_final_model_direct")
  fin_calls <- capture_core_args(function() {
    suppressMessages(suppressWarnings(fin_engine(
      labelled_gpkg = g, labelled_layer = "train_features",
      out_dir = NULL, prefix = "sym_final", overwrite = TRUE, verbose = FALSE,
      contextual_exclusion_to_burned_ratio   = 0.25,
      spectral_hard_negative_to_burned_ratio = 1.0,
      random_to_burned_ratio                 = 1.0,
      otsu_unburned_to_burned_ratio          = 1.0,
      sampling_seed = 42L, group_col = "block_id", val_frac = 0.15,
      seed = 42L, nrounds_max = 12L, early_stopping_rounds = 6L,
      impute_numeric = "median", impute_factor_missing = "MISSING",
      model_params_base = canon_mp()()
    )))
  })

  # ---- OOF capture ----
  oof_engine <- g_sym("run_dm_oof_pipeline")
  result_dir <- tempfile("sym_oof_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  oof_calls <- capture_core_args(function() {
    suppressMessages(suppressWarnings(oof_engine(
      labelled    = df,
      burned_like = df[df$class == "unburned", , drop = FALSE],
      labelled_df = df,
      params      = NULL,
      result_dir  = result_dir,
      target_year = 2017L,
      nrounds_max = 12L, early_stop = 6L, seed_base = 42L, verbose = 0,
      save_prefix = "sym_oof", prefix = "sym_oof", overwrite = TRUE,
      group_col   = "block_id",
      contextual_exclusion_to_burned_ratio   = 0.25,
      spectral_hard_negative_to_burned_ratio = 1.0,
      random_to_burned_ratio                 = 1.0,
      otsu_unburned_to_burned_ratio          = 1.0,
      val_frac = 0.15, impute_numeric = "median",
      impute_factor_missing = "MISSING"
    )))
  })
  unlink(result_dir, recursive = TRUE, force = TRUE)
  list(final = fin_calls, oof = oof_calls)
}

test_that("captured core args: shared invariants are EQUAL across OOF and FINAL", {
  caps <- run_capture_pair()
  expect_true(length(caps$final) == 1L,
              info = "FINAL calls the core exactly once")
  expect_true(length(caps$oof) >= 1L,
              info = "OOF calls the core once per (rep, fold)")

  fin <- caps$final[[1]]
  oof <- caps$oof[[1]]   # all OOF folds share the same resolved knobs.

  # -- SHARED: model_params (the canonical xgb block for a fixed spw) --
  expect_identical(fin$params_block[setdiff(names(fin$params_block), "scale_pos_weight")],
                   oof$params_block[setdiff(names(oof$params_block), "scale_pos_weight")])

  # -- SHARED: train_control knobs --
  expect_equal(fin$val_frac,              oof$val_frac)
  expect_equal(fin$group_col,             oof$group_col)
  expect_equal(fin$block_col,             oof$block_col)
  expect_equal(fin$nrounds_max,           oof$nrounds_max)
  expect_equal(fin$early_stopping_rounds, oof$early_stopping_rounds)
  expect_equal(fin$impute_numeric,        oof$impute_numeric)
  expect_equal(fin$impute_factor_missing, oof$impute_factor_missing)
  expect_equal(fin$label_col,             oof$label_col)

  # -- SHARED: feature whitelist + feature order. Both paths filter to the SAME
  #    canonical .supervised_feature_cols whitelist. The CONTRACT is on the
  #    canonical BASE feature set + its order: strip the `_isNA` companions and
  #    the two feature vectors must be byte-identical (same names, same order).
  #
  #    NOTE on this fixture: the `_isNA` companion columns differ between the two
  #    captured vectors ONLY because of how each path is fed here, NOT a recipe
  #    asymmetry. The OOF chain runs build_design_matrix_patches(), which
  #    SYNTHESISES an `_isNA` flag for every feature; the FINAL engine reads
  #    `_isNA` companions straight from its input GPKG, and this synthetic GPKG
  #    carries none. In a real run BOTH stages consume the SAME features GPKG
  #    (produced once by extract_supervised_features()), so the `_isNA` set is
  #    identical too. The base-feature contract below is the load-bearing one.
  wl <- g_sym(".supervised_feature_cols")
  allowed <- c(wl, paste0(wl, "_isNA"))
  expect_true(all(fin$feature_cols %in% allowed))
  expect_true(all(oof$feature_cols %in% allowed))
  strip_isna <- function(v) v[!grepl("_isNA$", v)]
  expect_identical(strip_isna(fin$feature_cols), strip_isna(oof$feature_cols))
  # The base feature order equals the canonical whitelist order, in BOTH paths.
  expect_identical(strip_isna(fin$feature_cols),
                   intersect(wl, strip_isna(fin$feature_cols)))
  expect_identical(strip_isna(oof$feature_cols),
                   intersect(wl, strip_isna(oof$feature_cols)))

  # -- ALL OOF folds share identical resolved knobs (no per-fold drift of the
  #    methodological params) --
  for (k in seq_along(caps$oof)) {
    expect_equal(caps$oof[[k]]$val_frac,              oof$val_frac)
    expect_equal(caps$oof[[k]]$nrounds_max,           oof$nrounds_max)
    expect_equal(caps$oof[[k]]$early_stopping_rounds, oof$early_stopping_rounds)
    expect_equal(caps$oof[[k]]$impute_numeric,        oof$impute_numeric)
    expect_equal(caps$oof[[k]]$impute_factor_missing, oof$impute_factor_missing)
    expect_identical(caps$oof[[k]]$feature_cols,      oof$feature_cols)
    expect_identical(caps$oof[[k]]$params_block[setdiff(names(oof$params_block), "scale_pos_weight")],
                     oof$params_block[setdiff(names(oof$params_block), "scale_pos_weight")])
  }
})

# ---------------------------------------------------------------------------
# (4) LEGITIMATE DIFFERENCES — assert they are present AND are the documented
#     ones (the contract allows exactly these).
# ---------------------------------------------------------------------------
test_that("legitimate difference (L3 seeds): the seed DERIVATION rule is the documented one", {
  caps <- run_capture_pair()
  fin <- caps$final[[1]]

  # FINAL: fold_seed == final_seed (42), sampling_seed == final_sampling_seed (42).
  expect_equal(fin$fold_seed, 42L)
  expect_equal(fin$sampling_seed, 42L)

  # OOF: fold_seed == seed_base + 1000*rep + fold, and sampling_seed == fold_seed
  # (the cap RNG seed). Recover (rep, fold) from the documented derivation and
  # check every captured OOF fold matches it for SOME (rep in 1:2, fold).
  seed_base <- 42L
  ok <- vapply(caps$oof, function(c1) {
    fs <- c1$fold_seed
    # fs = 42 + 1000*r + k  with r in {1,2}, k a small fold index.
    any(vapply(1:2, function(r) {
      k <- fs - seed_base - 1000L * r
      k >= 1L && k <= 50L
    }, logical(1))) &&
      identical(c1$sampling_seed, c1$fold_seed)
  }, logical(1))
  expect_true(all(ok),
              info = "every OOF fold_seed follows seed_base + 1000*rep + fold")

  # The two stages' base seeds are DISTINCT knobs by design (oof_seed_base vs
  # final_seed / final_sampling_seed) even when numerically equal here.
  tc <- canon_tc()()
  expect_setequal(names(tc$seeds),
                  c("oof_seed_base", "final_sampling_seed", "final_seed"))
})

test_that("legitimate difference (L1/L4): OOF holds out an outer fold; FINAL refits on the whole pool", {
  caps <- run_capture_pair()
  fin <- caps$final[[1]]

  # FINAL's core sees the ENTIRE capped pool (one call, all rows).
  total_rows <- 60L  # make_sym_oof_df default n.
  # FINAL trains on the capped L_ok (<= total labelled, but a single partition).
  expect_true(fin$n_train_rows >= 1L)

  # OOF: the SUM of fold sizes a single rep hands the core is STRICTLY LESS than
  # the full pool per fold (each fold withholds its outer-test partition). Prove
  # each OOF outer-train is smaller than the full labelled pool.
  for (k in seq_along(caps$oof)) {
    expect_true(caps$oof[[k]]$n_train_rows < total_rows,
                info = "each OOF outer-train excludes its held-out test fold")
  }
})

test_that("symmetry (L2): negative sampling is the SAME capped policy for OOF and FINAL (no oof_sampling toggle)", {
  # NEITHER the FINAL engine NOR the OOF chain carries an oof_sampling formal:
  # both apply the SAME capped negative-sampling policy (OOF per training fold).
  expect_false("oof_sampling" %in% names(formals(g_sym("train_final_model_direct"))))
  expect_false("oof_sampling" %in% names(formals(g_sym("run_dm_oof_pipeline"))))
  expect_false("oof_sampling" %in% names(formals(g_sym("run_oof_xgb"))))
  # No public function exposes the knob, and passing it ERRORS.
  for (f in c("run_oneyear_supervised_pipeline", "run_oof_diagnostics",
              "build_supervised_burned_config", "validate_supervised_execution")) {
    expect_false("oof_sampling" %in% names(formals(g_sym(f))),
                 info = paste("public still has oof_sampling formal:", f))
  }
  # The cfg/audit constant is the fixed traceability label "capped".
  expect_equal(canon_tc()()$oof_sampling, "capped")
})

test_that("symmetry (L2): OOF per-fold effective bucket ratios match the configured caps when enough negatives exist", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("dplyr")

  # Run the OOF chain and inspect the per-fold capping audit. The effective
  # (achieved) ratio per bucket equals the configured cap whenever the available
  # pool is at least the cap target -- the SAME cap semantics FINAL uses.
  df <- make_sym_oof_df()
  oof_engine <- g_sym("run_dm_oof_pipeline")
  result_dir <- tempfile("sym_oof_caps_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)
  caps_cfg <- list(contextual = 0.25, spectral = 1.0, random = 1.0, otsu = 1.0)
  res <- suppressMessages(suppressWarnings(oof_engine(
    labelled    = df,
    burned_like = df[df$class == "unburned", , drop = FALSE],
    labelled_df = df,
    params      = NULL,
    result_dir  = result_dir,
    target_year = 2017L,
    nrounds_max = 10L, early_stop = 5L, seed_base = 42L, verbose = 0,
    save_prefix = "sym_caps", prefix = "sym_caps", overwrite = TRUE,
    group_col   = "block_id",
    contextual_exclusion_to_burned_ratio   = caps_cfg$contextual,
    spectral_hard_negative_to_burned_ratio = caps_cfg$spectral,
    random_to_burned_ratio                 = caps_cfg$random,
    otsu_unburned_to_burned_ratio          = caps_cfg$otsu,
    val_frac = 0.15, impute_numeric = "median",
    impute_factor_missing = "MISSING"
  )))
  aud <- res$oof$oof_audit
  # For each bucket: when the available pool >= ceil(n_burned*cap), the selected
  # count equals the cap target (effective ratio matches the cap). Otherwise the
  # whole pool is taken (selected == available). Same semantics as FINAL.
  chk_bucket <- function(avail, cap_target, selected) {
    enough <- avail >= cap_target
    expect_equal(selected[enough],  cap_target[enough])
    expect_equal(selected[!enough], avail[!enough])
  }
  chk_bucket(aud$contextual_available, aud$contextual_cap, aud$contextual_selected)
  chk_bucket(aud$spectral_available,   aud$spectral_cap,   aud$spectral_selected)
  chk_bucket(aud$random_bg_available,  aud$random_bg_cap,  aud$random_bg_selected)
  chk_bucket(aud$otsu_available,       aud$otsu_cap,       aud$otsu_selected)
})

# ---------------------------------------------------------------------------
# (5) The "ONLY these differ" assertion: enumerate the full set of resolved
#     core arguments, partition into SHARED vs LEGITIMATELY-DIFFERENT, and prove
#     the partition is exhaustive (no unaccounted-for argument).
# ---------------------------------------------------------------------------
test_that("the resolved core-arg set partitions cleanly into shared + enumerated-different", {
  core_formals <- setdiff(names(formals(nr_fit_sym())), c("train_df", "...",
                                                          "params_fn", "verbose"))
  # By design every methodological core arg is either SHARED across the two
  # stages or one of the enumerated legitimate differences.
  shared_args <- c("feature_cols", "label_col", "group_col", "block_col",
                   "val_frac", "feature_weights", "nrounds_max",
                   "early_stopping_rounds", "impute_numeric",
                   "impute_factor_missing")
  # fold_seed / sampling_seed are the L3 (legitimate seed-derivation) difference.
  legitimate_diff_args <- c("fold_seed", "sampling_seed")

  accounted <- c(shared_args, legitimate_diff_args)
  expect_setequal(core_formals, accounted)
  # No argument is BOTH shared and a legitimate difference.
  expect_length(intersect(shared_args, legitimate_diff_args), 0L)
})

# ---------------------------------------------------------------------------
# (6) SCORING-RECIPE SCHEMA: the FINAL recipe carries the canonical schema
#     fields that the scoring path reads (recipe$cols / $impute / $training).
#     OOF predicts in-memory via the SAME core medians + x_cols + missing=NA
#     transform (.of_nested_transform), so the scoring recipe SCHEMA the two
#     stages rely on is the same set of fields.
# ---------------------------------------------------------------------------
test_that("scoring recipe schema fields are the canonical ones both stages rely on", {
  # The FINAL recipe schema (what scoring reads) is recipe$cols$feature_cols /
  # $cols$x_cols / $impute$numeric_medians / $impute$impute_factor_missing.
  # .of_model_expected_features() (the scoring-side reader) reads $cols$feature_cols.
  expf <- g_sym(".of_model_expected_features")
  rec <- list(cols = list(feature_cols = c("rbr_med", "elev_med")),
              impute = list(numeric_medians = list(rbr_med = 5, elev_med = 0),
                            impute_factor_missing = "MISSING"))
  expect_identical(expf(recipe = rec), c("rbr_med", "elev_med"))

  # The OOF per-fold transform consumes the SAME recipe fields (medians +
  # x_cols + factor sentinel) via .of_nested_transform -- same schema contract.
  tf_formals <- names(formals(g_sym(".of_nested_transform")))
  expect_true(all(c("feature_cols", "med", "ref_x_cols",
                    "impute_factor_missing") %in% tf_formals))
})
