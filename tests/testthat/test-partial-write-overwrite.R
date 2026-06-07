# 2026-06-06 (partial-write overwrite fix): the supervised stages honored
# `overwrite` for GPKG datasets but the CSV/RDS/TXT sidecars used to clobber
# unconditionally. These tests assert the new contract: with overwrite=FALSE an
# existing sidecar is NOT clobbered, and with overwrite=TRUE the sidecar is
# (re)written exactly as before.

run_oof_xgb_internal <- function() {
  get("run_oof_xgb", envir = asNamespace("OtsuFire"))
}
oof_pipeline_internal2 <- function() {
  get("run_dm_oof_pipeline", envir = asNamespace("OtsuFire"))
}
whitelist_internal2 <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

# Minimal self-contained OOF input for run_oof_xgb(): a tiny sparse matrix, a
# binary label, a labelled_df carrying fire_uid/class + the two fold columns.
make_run_oof_xgb_fixture <- function(n = 24L, seed = 7L) {
  set.seed(seed)
  X <- Matrix::Matrix(
    matrix(stats::runif(n * 3L), nrow = n,
           dimnames = list(NULL, c("rbr_med", "elev_med", "slope_med"))),
    sparse = TRUE
  )
  y <- rep(c(1L, 0L), length.out = n)
  labelled_df <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = ifelse(y == 1L, "burned", "unburned"),
    fold_rep1 = rep(c(1L, 2L), length.out = n),
    fold_rep2 = rep(c(2L, 1L), length.out = n),
    stringsAsFactors = FALSE
  )
  list(X = X, y = y, labelled_df = labelled_df)
}

test_that("run_oof_xgb honors overwrite for its CSV sidecars", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- make_run_oof_xgb_fixture()
  out_dir <- tempfile("oofxgb_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fn <- run_oof_xgb_internal()

  # First run (overwrite=TRUE) writes the two sidecars.
  suppressWarnings(suppressMessages(fn(
    XL_mat = fx$X, y = fx$y, labelled_df = fx$labelled_df,
    fold_cols = c("fold_rep1", "fold_rep2"),
    params = list(booster = "gbtree", objective = "binary:logistic",
                  eval_metric = "logloss", eta = 0.1, max_depth = 2,
                  nthread = 1),
    nrounds_max = 4L, early_stop = 3L, seed_base = 1L,
    out_dir = out_dir, prefix = "smoke", verbose = 0, overwrite = TRUE
  )))

  agg_csv  <- file.path(out_dir, "smoke_oof_agg.csv")
  long_csv <- file.path(out_dir, "smoke_oof_long.csv")
  expect_true(file.exists(agg_csv))
  expect_true(file.exists(long_csv))

  # Replace both sidecars with sentinel content.
  sentinel <- "SENTINEL_DO_NOT_CLOBBER"
  writeLines(sentinel, agg_csv)
  writeLines(sentinel, long_csv)

  # Second run with overwrite=FALSE must NOT clobber the sidecars.
  suppressWarnings(suppressMessages(fn(
    XL_mat = fx$X, y = fx$y, labelled_df = fx$labelled_df,
    fold_cols = c("fold_rep1", "fold_rep2"),
    params = list(booster = "gbtree", objective = "binary:logistic",
                  eval_metric = "logloss", eta = 0.1, max_depth = 2,
                  nthread = 1),
    nrounds_max = 4L, early_stop = 3L, seed_base = 1L,
    out_dir = out_dir, prefix = "smoke", verbose = 0, overwrite = FALSE
  )))
  expect_identical(readLines(agg_csv),  sentinel)
  expect_identical(readLines(long_csv), sentinel)

  # Third run with overwrite=TRUE must replace the sidecars (back to real CSV).
  suppressWarnings(suppressMessages(fn(
    XL_mat = fx$X, y = fx$y, labelled_df = fx$labelled_df,
    fold_cols = c("fold_rep1", "fold_rep2"),
    params = list(booster = "gbtree", objective = "binary:logistic",
                  eval_metric = "logloss", eta = 0.1, max_depth = 2,
                  nthread = 1),
    nrounds_max = 4L, early_stop = 3L, seed_base = 1L,
    out_dir = out_dir, prefix = "smoke", verbose = 0, overwrite = TRUE
  )))
  expect_false(identical(readLines(agg_csv), sentinel))
})

# run_dm_oof_pipeline writes the OOF metric/summary CSV+TXT sidecars AFTER the
# design matrix is built. With overwrite=FALSE and a pre-existing sidecar (but
# no design-matrix bundle yet, so the stage proceeds), the sidecar must NOT be
# clobbered. This exercises the .write_if_allowed() guards added in the OOF
# wrapper.
make_oof_pipeline_fixture <- function(n = 30L, seed = 31L) {
  set.seed(seed)
  whitelist <- whitelist_internal2()
  feat_df <- as.data.frame(
    matrix(stats::runif(n * length(whitelist)),
           nrow = n, dimnames = list(NULL, whitelist))
  )
  for (nm in c("hotspot_available", "hs_in_poly", "hs_in_buffer",
               "hs_used_n", "hs_hiConf_n", "hs_support_present",
               "hs_no_support_when_available", "hs_only_buffer_support")) {
    if (nm %in% names(feat_df)) feat_df[[nm]] <- as.integer(round(feat_df[[nm]]))
  }
  admin_df <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = rep(c("burned", "unburned"), length.out = n),
    source    = rep(c("burned_truth", "random_burnable_background"),
                    length.out = n),
    poly_id   = sprintf("p_%03d", seq_len(n)),
    block_id  = rep(seq_len(5), length.out = n),
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2 = rep(c(2L, 1L, 3L), length.out = n),
    stringsAsFactors = FALSE
  )
  full_df <- cbind(admin_df, feat_df)
  bl_df <- full_df
  bl_df$class <- NULL
  bl_df$fire_uid <- sprintf("bl_%03d", seq_len(n))
  bl_df$poly_id  <- sprintf("bl_p_%03d", seq_len(n))
  list(labelled = full_df, burned_like = bl_df,
       labelled_df = full_df)
}

test_that("run_dm_oof_pipeline does not clobber CSV/TXT sidecars when overwrite=FALSE", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- make_oof_pipeline_fixture()
  result_dir <- tempfile("oofpipe_ow_")
  dir.create(file.path(result_dir, "05_OOF"), recursive = TRUE,
             showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Pre-create a sentinel OOF metrics-summary TXT (a sidecar the pipeline would
  # otherwise write). No design-matrix bundle exists yet, so the stage proceeds
  # past build_design_matrix_patches() and reaches the sidecar writes.
  prefix <- "2005_patch"
  summary_txt <- file.path(result_dir, "05_OOF",
                           paste0(prefix, "_oof_metrics_summary.txt"))
  sentinel <- "SENTINEL_SUMMARY_DO_NOT_CLOBBER"
  writeLines(sentinel, summary_txt)

  fn <- oof_pipeline_internal2()
  suppressWarnings(suppressMessages(fn(
    labelled    = fx$labelled,
    burned_like = fx$burned_like,
    labelled_df = fx$labelled_df,
    params      = list(booster = "gbtree", objective = "binary:logistic",
                       eval_metric = "logloss", eta = 0.1, max_depth = 3,
                       nthread = 1),
    result_dir  = result_dir,
    target_year = 2005L,
    nrounds_max = 6L, early_stop = 4L, seed_base = 11L,
    verbose     = 0, prefix = prefix, save_prefix = prefix,
    overwrite   = FALSE
  )))

  # The pre-existing summary TXT must NOT have been clobbered.
  expect_identical(readLines(summary_txt), sentinel)

  # ...but a fresh sidecar with no pre-existing file IS written (first run).
  agg_csv <- file.path(result_dir, "05_OOF", paste0(prefix, "_oof_agg.csv"))
  expect_true(file.exists(agg_csv))
  expect_false(identical(readLines(agg_csv), sentinel))
})
