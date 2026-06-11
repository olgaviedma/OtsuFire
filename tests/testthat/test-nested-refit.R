# B1 (2026-06-07): unit tests for the nested_refit protocol. These exercise
# the shared core (.of_nested_refit_fit + median helpers) and the per-fold OOF
# logic (run_oof_xgb / run_dm_oof_pipeline) on small fixtures. The full pipeline
# is NOT run here.

ns <- asNamespace("OtsuFire")

nr_fit       <- function() get(".of_nested_refit_fit", envir = ns)
nr_fit_med   <- function() get(".of_nested_fit_medians", envir = ns)
nr_apply_med <- function() get(".of_nested_apply_medians", envir = ns)
nr_transform <- function() get(".of_nested_transform", envir = ns)
canon_params <- function() get(".of_canonical_xgb_params", envir = ns)
run_oof_xgb_fn <- function() get("run_oof_xgb", envir = ns)
run_dm_oof_fn  <- function() get("run_dm_oof_pipeline", envir = ns)

# ---------------------------------------------------------------------------
# (a) medians fit on inner_tr / refit-train only, never on val/test.
# ---------------------------------------------------------------------------
test_that("nested core fits medians on the fit rows only (train), never on test", {
  fn <- nr_fit_med()
  # Train rows (1:4) have values 1,2,3 (median 2); test rows have huge values
  # that would shift the median if leaked.
  X <- data.frame(a = c(1, 2, 3, NA, 1000, 2000, 3000, NA))
  med_tr <- fn(X, fit_rows = 1:4, impute_numeric = "median")
  expect_equal(med_tr$a, 2)

  # If (incorrectly) computed over all rows the median would be far larger.
  med_all <- fn(X, fit_rows = 1:8, impute_numeric = "median")
  expect_true(med_all$a > 2)
})

test_that("nested_refit_fit imputes held-out (inner-val) rows using train medians", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fn <- nr_fit()
  set.seed(11)
  n <- 60
  df <- data.frame(
    class    = rep(c("burned", "unburned"), each = n / 2),
    block_id = rep(seq_len(10), length.out = n),
    rbr_med  = c(rnorm(n / 2, 5), rnorm(n / 2, 1)),
    elev_med = rnorm(n),
    stringsAsFactors = FALSE
  )
  fit <- fn(
    train_df = df, feature_cols = c("rbr_med", "elev_med"),
    label_col = "class", group_col = "block_id", block_col = "block_id",
    val_frac = 0.2, params_fn = canon_params(),
    sampling_seed = 42, fold_seed = 7, nrounds_max = 20,
    early_stopping_rounds = 8
  )
  # Refit medians are fit on ALL train rows; both phases use train-derived
  # medians only. The recipe medians must be finite for non-degenerate cols.
  expect_true(is.finite(fit$medians$rbr_med))
  expect_true(is.finite(fit$medians$elev_med))
  expect_s3_class(fit$model, "xgb.Booster")
  expect_true(fit$best_iteration >= 1L)
  # spw computed from full-train labels.
  expect_equal(fit$spw_refit, sum(df$class == "unburned") / sum(df$class == "burned"))
})

# ---------------------------------------------------------------------------
# (d) degenerate all-NA-in-train column stays NA (not fabricated to 0).
# ---------------------------------------------------------------------------
test_that("degenerate all-NA-in-train column keeps NA median (no fabrication)", {
  fit_med <- nr_fit_med()
  apply_med <- nr_apply_med()
  X <- data.frame(b = c(NA, NA, NA, NA, 5, 6, 7, 8))
  med <- fit_med(X, fit_rows = 1:4, impute_numeric = "median")
  expect_true(is.na(med$b))
  out <- apply_med(X, med)
  # Train rows stay NA (xgboost treats as missing); no 0 fabricated.
  expect_true(all(is.na(out$b[1:4])))
})

test_that("nested_refit_fit reports degenerate columns in the audit", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fn <- nr_fit()
  set.seed(3)
  n <- 50
  df <- data.frame(
    class    = rep(c("burned", "unburned"), each = n / 2),
    block_id = rep(seq_len(10), length.out = n),
    rbr_med  = c(rnorm(n / 2, 5), rnorm(n / 2, 1)),
    dead_col = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  fit <- fn(
    train_df = df, feature_cols = c("rbr_med", "dead_col"),
    label_col = "class", block_col = "block_id", group_col = "block_id",
    val_frac = 0.2, params_fn = canon_params(),
    sampling_seed = 1, fold_seed = 1, nrounds_max = 15,
    early_stopping_rounds = 6
  )
  expect_true("dead_col" %in% fit$degenerate_cols)
  expect_true(is.na(fit$medians$dead_col))
})

# ---------------------------------------------------------------------------
# Build a synthetic OOF input pair for run_oof_xgb / run_dm_oof_pipeline.
# ---------------------------------------------------------------------------
make_nested_oof_df <- function(n = 60, seed = 5L) {
  set.seed(seed)
  whitelist <- get(".supervised_feature_cols", envir = ns)
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
  # Make folds respect block_id (no block spans two folds within a rep) so the
  # no-overlap assertions in the nested path hold.
  blk_to_fold1 <- rep(c(1L, 2L, 3L), length.out = 12)
  blk_to_fold2 <- rep(c(3L, 1L, 2L), length.out = 12)
  admin$fold_rep1 <- blk_to_fold1[admin$block_id]
  admin$fold_rep2 <- blk_to_fold2[admin$block_id]
  cbind(admin, feat_df)
}

# ---------------------------------------------------------------------------
# (c) caps applied to outer_train only; outer_test row count unchanged.
# (b) outer_test never enters any watchlist.
# (e) fold_seed keeps both repeat and fold.
# ---------------------------------------------------------------------------
test_that("nested OOF runs, predicts every labelled unit, and emits a per-fold audit", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("dplyr")

  df <- make_nested_oof_df()
  result_dir <- tempfile("nested_oof_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fn <- run_dm_oof_fn()
  res <- suppressMessages(suppressWarnings(fn(
    labelled    = df,
    burned_like = df[df$class == "unburned", , drop = FALSE],
    labelled_df = df,
    params      = NULL,
    result_dir  = result_dir,
    target_year = 2017L,
    nrounds_max = 12L,
    early_stop  = 6L,
    seed_base   = 100L,
    verbose     = 0,
    save_prefix = "nr_smoke",
    prefix      = "nr_smoke",
    overwrite   = TRUE,
    group_col   = "block_id",
    contextual_exclusion_to_burned_ratio   = 0.25,
    spectral_hard_negative_to_burned_ratio = 1.0,
    random_to_burned_ratio                 = 1.0,
    otsu_unburned_to_burned_ratio          = 1.0,
    # The per-fold core requires these resolved controls.
    val_frac = 0.15, impute_numeric = "median",
    impute_factor_missing = "MISSING"
  )))

  # Every labelled unit gets an OOF prediction (outer_test never sub-sampled).
  expect_equal(nrow(res$oof$oof_agg), nrow(df))
  # Audit present, one row per (rep, fold).
  aud <- res$oof$oof_audit
  expect_true(is.data.frame(aud))
  expect_true(all(c("rep", "fold", "spw_selection", "spw_refit",
                    "best_iteration", "fold_seed",
                    "outer_test_capped", "outer_test_used_for_fit",
                    "n_outer_test", "n_outer_train_pre",
                    "n_outer_train_post") %in% names(aud)))
  # OOF always uses the capped policy (constant stamped in the audit).
  expect_true(all(aud$oof_sampling == "capped"))
  # outer_test never capped / never used for fitting.
  expect_true(all(aud$outer_test_capped == FALSE))
  expect_true(all(aud$outer_test_used_for_fit == FALSE))
  # Caps shrink (or keep) the outer-train; never grow it.
  expect_true(all(aud$n_outer_train_post <= aud$n_outer_train_pre))
  # fold_seed encodes both rep and fold (seed_base + 1000*r + k).
  expect_equal(aud$fold_seed, 100L + 1000L * aud$rep + aud$fold)
})

# ---------------------------------------------------------------------------
# Per-fold capping LOG: the OOF audit records, PER FOLD: n_burned,
# available-per-bucket, max cap, selected, and effective (achieved) ratio.
# (Replaces the deleted "full"=no-capping test: capping ALWAYS runs now.)
# ---------------------------------------------------------------------------
test_that("nested OOF capping ALWAYS runs and records the per-fold capping audit", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("dplyr")

  df <- make_nested_oof_df(seed = 9L)
  result_dir <- tempfile("nested_oof_cap_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fn <- run_dm_oof_fn()
  res <- suppressMessages(suppressWarnings(fn(
    labelled    = df,
    burned_like = df[df$class == "unburned", , drop = FALSE],
    labelled_df = df,
    params      = NULL,
    result_dir  = result_dir,
    target_year = 2017L,
    nrounds_max = 10L, early_stop = 5L, seed_base = 50L, verbose = 0,
    save_prefix = "nr_cap", prefix = "nr_cap", overwrite = TRUE,
    group_col   = "block_id",
    contextual_exclusion_to_burned_ratio   = 0.25,
    spectral_hard_negative_to_burned_ratio = 1.0,
    random_to_burned_ratio                 = 1.0,
    otsu_unburned_to_burned_ratio          = 1.0,
    # Gate 1B (2026-06-07): nested path now requires these resolved controls.
    val_frac = 0.15, impute_numeric = "median",
    impute_factor_missing = "MISSING"
  )))
  aud <- res$oof$oof_audit
  # The capping always runs: the audit carries the per-fold, per-bucket schema.
  capping_fields <- c(
    "n_burned",
    "contextual_available", "contextual_cap", "contextual_selected",
    "contextual_effective_ratio",
    "spectral_available", "spectral_cap", "spectral_selected",
    "spectral_effective_ratio",
    "random_bg_available", "random_bg_cap", "random_bg_selected",
    "random_bg_effective_ratio",
    "otsu_available", "otsu_cap", "otsu_selected", "otsu_effective_ratio"
  )
  expect_true(all(capping_fields %in% names(aud)),
              info = paste("missing per-fold capping audit fields:",
                           paste(setdiff(capping_fields, names(aud)),
                                 collapse = ", ")))
  # Effective ratio = selected / n_burned (achieved ratio).
  ok <- aud$n_burned > 0
  expect_equal(aud$contextual_effective_ratio[ok],
               aud$contextual_selected[ok] / aud$n_burned[ok])
  expect_equal(aud$spectral_effective_ratio[ok],
               aud$spectral_selected[ok] / aud$n_burned[ok])
  # Selected never exceeds the cap (target ceiling) nor the available pool.
  expect_true(all(aud$contextual_selected <= aud$contextual_cap))
  expect_true(all(aud$contextual_selected <= aud$contextual_available))
})

# ---------------------------------------------------------------------------
# (f) a dropped cap arg in the OOF chain ERRORS (required formal).
# ---------------------------------------------------------------------------
test_that("dropped cap arg in run_dm_oof_pipeline(nested_refit) ERRORS", {
  df <- make_nested_oof_df(seed = 13L)
  result_dir <- tempfile("nested_oof_err_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fn <- run_dm_oof_fn()
  expect_error(
    suppressMessages(suppressWarnings(fn(
      labelled    = df,
      burned_like = df[df$class == "unburned", , drop = FALSE],
      labelled_df = df,
      params      = NULL,
      result_dir  = result_dir,
      target_year = 2017L,
      nrounds_max = 8L, early_stop = 4L, seed_base = 7L, verbose = 0,
      save_prefix = "nr_err", prefix = "nr_err", overwrite = TRUE,
      group_col   = "block_id",
      # Gate 1B: supply the always-needed nested controls so the spectral-cap
      # guard (not the val_frac/impute guard) is what surfaces.
      val_frac = 0.15, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      # spectral cap DROPPED on purpose.
      contextual_exclusion_to_burned_ratio   = 0.25,
      random_to_burned_ratio                 = 1.0,
      otsu_unburned_to_burned_ratio          = 1.0
    ))),
    regexp = "spectral_hard_negative_to_burned_ratio"
  )
})

# ---------------------------------------------------------------------------
# (Constraint 5) The nested OOF path ASSERTS no fire_uid / block_id overlap
# between outer_train and outer_test. Prove the guard is live by handing it a
# fold assignment where the SAME fire_uid appears in two folds (an outer
# train/test leak), and expect a hard stop.
# ---------------------------------------------------------------------------
test_that("nested OOF errors when outer_train and outer_test share a fire_uid (leak guard live)", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("dplyr")

  df <- make_nested_oof_df(seed = 31L)
  # Force a fire_uid to straddle two folds within rep1: copy uid of a fold-1 row
  # onto a fold-2 row, so that fire_uid is in both outer_train and outer_test.
  i1 <- which(df$fold_rep1 == 1L)[1]
  i2 <- which(df$fold_rep1 == 2L)[1]
  df$fire_uid[i2] <- df$fire_uid[i1]

  fn <- run_dm_oof_fn()
  expect_error(
    suppressMessages(suppressWarnings(fn(
      labelled    = df,
      burned_like = df[df$class == "unburned", , drop = FALSE],
      labelled_df = df,
      params      = NULL,
      result_dir  = tempfile("nested_oof_leak_"),
      target_year = 2017L,
      nrounds_max = 8L, early_stop = 4L, seed_base = 3L, verbose = 0,
      save_prefix = "nr_leak", prefix = "nr_leak", overwrite = TRUE,
      group_col   = "block_id",
      val_frac = 0.15, impute_numeric = "median",
      impute_factor_missing = "MISSING",
      contextual_exclusion_to_burned_ratio   = 0.25,
      spectral_hard_negative_to_burned_ratio = 1.0,
      random_to_burned_ratio                 = 1.0,
      otsu_unburned_to_burned_ratio          = 1.0
    ))),
    regexp = "fire_uid overlap"
  )
})

test_that("the 4 cap formals have NO default in run_oof_xgb / run_dm_oof_pipeline", {
  f_xgb <- formals(run_oof_xgb_fn())
  f_dm  <- formals(run_dm_oof_fn())
  for (nm in c("contextual_exclusion_to_burned_ratio",
               "spectral_hard_negative_to_burned_ratio",
               "random_to_burned_ratio",
               "otsu_unburned_to_burned_ratio")) {
    expect_true(identical(f_xgb[[nm]], quote(expr = )),
                info = paste("run_oof_xgb cap has a default:", nm))
    expect_true(identical(f_dm[[nm]], quote(expr = )),
                info = paste("run_dm_oof_pipeline cap has a default:", nm))
  }
})

# ---------------------------------------------------------------------------
# Single-protocol contract (Contract #9 + #10): no public function accepts
# training_protocol, no engine carries a training_protocol formal, and the cfg
# constant is the fixed traceability label.
# ---------------------------------------------------------------------------
test_that("Contract #10: no PUBLIC function accepts training_protocol", {
  for (f in c("run_oneyear_supervised_pipeline", "run_oof_diagnostics",
              "train_final_burned_model", "build_supervised_burned_config",
              "validate_supervised_execution")) {
    fm <- formals(get(f, envir = ns))
    expect_false("training_protocol" %in% names(fm),
                 info = paste("public still has training_protocol formal:", f))
  }
})

test_that("Contract #10: passing training_protocol to a public fn ERRORS", {
  # build_supervised_burned_config / train_final_burned_model / run_oof_diagnostics
  # have no `...`, so R raises "unused argument".
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                   internal_decisions = tempfile(fileext = ".gpkg"),
                                   change_index = tempfile(fileext = ".tif"),
                                   target_year = 2017L,
                                   training_protocol = "legacy"),
    regexp = "unused argument|training_protocol")
  # run_oneyear_supervised_pipeline has a `...` -> explicit guard message.
  expect_error(
    run_oneyear_supervised_pipeline(config = structure(list(),
      class = "otsufire_supervised_burned_config"),
      training_protocol = "nested_refit"),
    regexp = "training_protocol is no longer an argument")
})

test_that("Contract #9: no engine carries a training_protocol formal", {
  for (f in c("run_dm_oof_pipeline", "run_oof_xgb", "train_final_model_direct")) {
    fm <- formals(get(f, envir = ns))
    expect_false("training_protocol" %in% names(fm),
                 info = paste("engine still has training_protocol formal:", f))
  }
})

test_that("training_protocol survives only as a fixed internal cfg constant", {
  tc <- get(".of_canonical_train_control", envir = ns)()
  expect_equal(tc$training_protocol, "nested_refit")
})

# ---------------------------------------------------------------------------
# Single capped-policy contract: no public function accepts oof_sampling, no
# OOF engine carries an oof_sampling formal, passing oof_sampling ERRORS, and
# the cfg/audit constant is the fixed traceability label "capped".
# ---------------------------------------------------------------------------
test_that("no PUBLIC function accepts oof_sampling", {
  for (f in c("run_oneyear_supervised_pipeline", "run_oof_diagnostics",
              "build_supervised_burned_config", "validate_supervised_execution")) {
    fm <- formals(get(f, envir = ns))
    expect_false("oof_sampling" %in% names(fm),
                 info = paste("public still has oof_sampling formal:", f))
  }
})

test_that("no OOF engine carries an oof_sampling formal", {
  for (f in c("run_dm_oof_pipeline", "run_oof_xgb")) {
    fm <- formals(get(f, envir = ns))
    expect_false("oof_sampling" %in% names(fm),
                 info = paste("engine still has oof_sampling formal:", f))
  }
})

test_that("passing oof_sampling to a public fn ERRORS", {
  # build_supervised_burned_config / run_oof_diagnostics / validate_supervised_execution
  # have no `...`, so R raises "unused argument".
  expect_error(
    build_supervised_burned_config(scenario = "balanced",
                                   internal_decisions = tempfile(fileext = ".gpkg"),
                                   change_index = tempfile(fileext = ".tif"),
                                   target_year = 2017L,
                                   oof_sampling = "full"),
    regexp = "unused argument|oof_sampling")
  # run_oneyear_supervised_pipeline has a `...` -> explicit guard message.
  expect_error(
    run_oneyear_supervised_pipeline(config = structure(list(),
      class = "otsufire_supervised_burned_config"),
      oof_sampling = "full"),
    regexp = "oof_sampling is no longer an argument")
})

test_that("no source-level oof_sampling 'full' branch / validator survives in R/", {
  r_dir <- system.file("R", package = "OtsuFire")
  # Fall back to the source tree when not installed (load_all dev workflow).
  src_files <- list.files(
    file.path(testthat::test_path("..", ".."), "R"),
    pattern = "\\.R$", full.names = TRUE)
  if (length(src_files) == 0L && nzchar(r_dir)) {
    src_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  }
  skip_if(length(src_files) == 0L, "no R/ source files visible to grep")
  body_txt <- unlist(lapply(src_files, readLines, warn = FALSE))
  # No "full"-branch logic for oof sampling anywhere.
  expect_false(any(grepl('oof_sampling\\s*==\\s*"full"', body_txt)))
  expect_false(any(grepl('identical\\(\\s*oof_sampling\\s*,\\s*"full"', body_txt)))
  expect_false(any(grepl('match\\.arg\\([^)]*oof_sampling', body_txt)))
})

test_that("oof_sampling survives only as a fixed internal cfg/audit constant", {
  tc <- get(".of_canonical_train_control", envir = ns)()
  expect_equal(tc$oof_sampling, "capped")
})

test_that("Contract #9: NO legacy training branch/validator remains in R/", {
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if(!dir.exists(r_dir), "R/ source dir not available")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  src <- unlist(lapply(files, function(f) readLines(f, warn = FALSE)))
  # No legacy training-protocol gating: no `training_protocol == "legacy"`, no
  # match.arg(training_protocol, ...), no c("legacy","nested_refit") enum.
  expect_false(any(grepl('training_protocol\\s*==\\s*"legacy"', src)))
  expect_false(any(grepl('identical\\(training_protocol', src)))
  expect_false(any(grepl('match\\.arg\\(training_protocol', src)))
  expect_false(any(grepl('c\\("legacy",\\s*"nested_refit"\\)', src)))
})

# ---------------------------------------------------------------------------
# Contract #7 (source-level): OOF and FINAL share the SAME training core.
# ---------------------------------------------------------------------------
test_that("Contract #7: OOF and FINAL both call .of_nested_refit_fit (shared core)", {
  oof_src   <- deparse(get("run_oof_xgb", envir = ns))
  final_src <- deparse(get("train_final_model_direct", envir = ns))
  expect_true(any(grepl("\\.of_nested_refit_fit", oof_src)))
  expect_true(any(grepl("\\.of_nested_refit_fit", final_src)))
})

test_that("neg_type is in the OOF id_cols default (per-fold bucketing)", {
  idc <- eval(formals(run_dm_oof_fn())$id_cols)
  expect_true("neg_type" %in% idc)
})

# ---------------------------------------------------------------------------
# (b) structural: build_design_matrix_patches(defer_impute=TRUE) leaves NAs.
# ---------------------------------------------------------------------------
test_that("defer_impute returns prepared frame with un-imputed numeric NAs", {
  skip_if_not_installed("Matrix")
  bdm <- get("build_design_matrix_patches", envir = ns)
  df <- make_nested_oof_df(seed = 21L)
  # Inject an NA into a whitelist numeric feature.
  df$rbr_med[1] <- NA
  out <- suppressMessages(bdm(
    labelled = df, burned_like = df[df$class == "unburned", , drop = FALSE],
    defer_impute = TRUE, save_dir = NULL, verbose = FALSE
  ))
  expect_true(isTRUE(out$deferred))
  expect_true("prepared_labelled" %in% names(out))
  # The deferred frame must still carry the NA (not median-imputed).
  expect_true(is.na(out$prepared_labelled$rbr_med[1]))
  # And it must carry the _isNA companion flag for that row.
  expect_true("rbr_med_isNA" %in% names(out$prepared_labelled))
  expect_equal(out$prepared_labelled$rbr_med_isNA[1], 1L)
})

# ---------------------------------------------------------------------------
# (Deliverable 6) The OUTER_TEST medians do NOT affect the OOF model's
# imputation or prediction. Construct a fixture where a feature's median
# differs SHARPLY between the (outer-)train rows and the (outer-)test rows, and
# prove that:
#   (i)  the recipe median equals the TRAIN median (never the test median);
#   (ii) an NA in the held-out test frame is imputed with the TRAIN median, so
#        the transformed test design matrix is byte-identical regardless of what
#        the test rows' own (non-NA) values were;
#   (iii) consequently the model's prediction for a test row with a missing
#        feature is invariant to the test-set distribution of that feature.
# ---------------------------------------------------------------------------
test_that("outer_test medians never affect the OOF model (imputation/prediction depend on train only)", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  fit_fn       <- nr_fit()
  transform_fn <- nr_transform()
  fit_med      <- nr_fit_med()

  set.seed(101)
  n_tr <- 80
  # TRAIN frame: feature `rbr_med` is tightly centred near 5 (train median ~5);
  # `elev_med` is informative for the label so xgboost actually splits on it.
  train_df <- data.frame(
    class    = rep(c("burned", "unburned"), each = n_tr / 2),
    block_id = rep(seq_len(16), length.out = n_tr),
    rbr_med  = rnorm(n_tr, mean = 5, sd = 0.1),
    elev_med = c(rnorm(n_tr / 2, 3), rnorm(n_tr / 2, -3)),
    stringsAsFactors = FALSE
  )
  train_median_rbr <- stats::median(train_df$rbr_med)

  fit <- fit_fn(
    train_df = train_df, feature_cols = c("rbr_med", "elev_med"),
    label_col = "class", group_col = "block_id", block_col = "block_id",
    val_frac = 0.2, params_fn = canon_params(),
    sampling_seed = 1, fold_seed = 1, nrounds_max = 30,
    early_stopping_rounds = 10
  )

  # (i) recipe median equals the TRAIN median, not anything from a test set.
  expect_equal(fit$medians$rbr_med, train_median_rbr, tolerance = 1e-9)

  # Build TWO held-out test frames whose `rbr_med` distributions are wildly
  # different from each other AND from train (medians ~1000 vs ~-1000), but in
  # both the FIRST row's `rbr_med` is MISSING (NA) -> must be imputed from the
  # TRAIN median in both cases.
  mk_test <- function(center) {
    data.frame(
      rbr_med  = c(NA, rnorm(9, mean = center, sd = 0.1)),
      elev_med = rep(0, 10),
      stringsAsFactors = FALSE
    )
  }
  te_hi <- mk_test( 1000)
  te_lo <- mk_test(-1000)

  X_hi <- transform_fn(te_hi, feature_cols = c("rbr_med", "elev_med"),
                       med = fit$medians, ref_x_cols = fit$x_cols)
  X_lo <- transform_fn(te_lo, feature_cols = c("rbr_med", "elev_med"),
                       med = fit$medians, ref_x_cols = fit$x_cols)

  # (ii) the imputed (formerly NA) test cell carries the TRAIN median in BOTH
  # frames -- the test-set distribution is irrelevant to imputation.
  expect_equal(unname(X_hi[1, "rbr_med"]), train_median_rbr, tolerance = 1e-9)
  expect_equal(unname(X_lo[1, "rbr_med"]), train_median_rbr, tolerance = 1e-9)

  # (iii) prediction for the missing-feature test row is identical no matter
  # which (wildly different) test set it was drawn from: the model only ever
  # sees the train-derived imputation, never an outer-test statistic.
  p_hi <- predict(fit$model,
                  xgboost::xgb.DMatrix(X_hi[1, , drop = FALSE], missing = NA))
  p_lo <- predict(fit$model,
                  xgboost::xgb.DMatrix(X_lo[1, , drop = FALSE], missing = NA))
  expect_equal(p_hi, p_lo, tolerance = 1e-9)

  # Control: a median fit on the TEST rows would be ~1000 / ~-1000, i.e. nothing
  # like the train median -- confirming the fixture really is "sharply
  # different" and the recipe ignored it.
  med_from_test_hi <- fit_med(te_hi[-1, , drop = FALSE], fit_rows = seq_len(9),
                              impute_numeric = "median")
  expect_gt(abs(med_from_test_hi$rbr_med - train_median_rbr), 100)
})
