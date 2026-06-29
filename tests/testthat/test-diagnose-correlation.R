# Tests for the correlation / redundancy module (H) of diagnose_training_pools().

# --- synthetic frame with a deliberately redundant pair -----------------------
# rbr_med and rbr_aw_med are near-perfectly correlated (aw = 2*med + small noise);
# persist_ratio and doy_iqr are independent. group separates on rbr_med.
make_corr_fixture <- function(n = 300, seed = 7) {
  set.seed(seed)
  grp  <- sample(c("a", "b", "c"), n, TRUE)
  base <- c(a = 200, b = 120, c = 260)[grp]
  rbr  <- stats::rnorm(n, base, 30)
  data.frame(
    pool_source   = grp,
    rbr_med       = rbr,
    rbr_aw_med    = 2 * rbr + stats::rnorm(n, 0, 1e-3),  # ~perfectly correlated
    persist_ratio = stats::runif(n),
    doy_iqr       = stats::rpois(n, 20),
    stringsAsFactors = FALSE
  )
}

test_that(".pool_correlation clusters a redundant pair and keeps the important one", {
  df  <- make_corr_fixture()
  fts <- c("rbr_med", "rbr_aw_med", "persist_ratio", "doy_iqr")
  # rbr_med has higher importance than rbr_aw_med -> should be the representative
  imp <- c(rbr_med = 0.50, rbr_aw_med = 0.40, persist_ratio = 0.05, doy_iqr = 0.02)
  cr  <- .pool_correlation(df, fts, importance = imp,
                           method = "spearman", thr = 0.9)
  expect_false(is.null(cr))
  # rbr_med + rbr_aw_med share one cluster
  cl_med <- cr$clusters$cluster[cr$clusters$feature == "rbr_med"]
  cl_aw  <- cr$clusters$cluster[cr$clusters$feature == "rbr_aw_med"]
  expect_identical(cl_med, cl_aw)
  # higher-importance feature is the representative
  expect_identical(unique(cr$clusters$representative[cr$clusters$cluster == cl_med]),
                   "rbr_med")
  # the other one is dropped from the reduced set
  expect_true("rbr_med" %in% cr$reduced_set)
  expect_false("rbr_aw_med" %in% cr$reduced_set)
  expect_true("rbr_aw_med" %in% cr$dropped$feature)
  # high_pairs is family-tagged
  expect_true(all(c("feature1", "feature2", "cor", "family1", "family2") %in%
                    names(cr$high_pairs)))
})

test_that("uncorrelated features remain singletons (all in reduced_set)", {
  df  <- make_corr_fixture()
  fts <- c("persist_ratio", "doy_iqr")          # independent
  cr  <- .pool_correlation(df, fts, method = "spearman", thr = 0.9)
  expect_false(is.null(cr))
  expect_setequal(cr$reduced_set, fts)          # nothing dropped
  expect_equal(nrow(cr$dropped), 0L)
  expect_equal(cr$n_reduced, cr$n_features)
})

test_that(".pool_correlation returns NULL on < 2 features (graceful)", {
  df <- make_corr_fixture()
  expect_null(.pool_correlation(df, "rbr_med"))
  expect_null(.pool_correlation(df, character(0)))
})

test_that("diagnose_training_pools() exposes a $correlation component", {
  out <- diagnose_training_pools(make_corr_fixture(), group_col = "pool_source",
                                 verbose = FALSE, make_plots = FALSE)
  expect_false(is.null(out$correlation))
  expect_true(all(c("matrix", "high_pairs", "clusters", "representative",
                    "reduced_set", "dropped") %in% names(out$correlation)))
  # redundancy present -> reduced_set strictly shorter than the full feature set
  expect_lt(length(out$correlation$reduced_set),
            length(out$feature_info$features_used))
  # representative chosen by epsilon^2: rbr_med separates groups, so it wins
  expect_true("rbr_med" %in% out$correlation$reduced_set)
  expect_false("rbr_aw_med" %in% out$correlation$reduced_set)
})

test_that("run_correlation = FALSE -> correlation is NULL", {
  out <- diagnose_training_pools(make_corr_fixture(), group_col = "pool_source",
                                 run_correlation = FALSE,
                                 verbose = FALSE, make_plots = FALSE)
  expect_null(out$correlation)
  # nothing else broke: documented structure intact
  for (nm in c("summary", "feature_ranking", "zscore_fingerprint", "pca"))
    expect_true(nm %in% names(out), info = nm)
})

test_that("correlation CSVs are written when output_dir is set", {
  od  <- withr::local_tempdir()
  out <- diagnose_training_pools(make_corr_fixture(), group_col = "pool_source",
                                 output_dir = od, make_plots = FALSE,
                                 verbose = FALSE)
  files <- list.files(od)
  expect_true(any(grepl("15_correlation_high_pairs.csv$", files)))
  expect_true(any(grepl("16_correlation_clusters.csv$", files)))
  expect_true(any(grepl("17_correlation_reduced_set.csv$", files)))
})

test_that("correlation heatmap is produced when ggplot2 is available", {
  skip_if_not_installed("ggplot2")
  od  <- withr::local_tempdir()
  out <- diagnose_training_pools(make_corr_fixture(), group_col = "pool_source",
                                 output_dir = od, make_plots = TRUE,
                                 verbose = FALSE)
  expect_true(!is.null(out$plot_paths$correlation_heatmap))
  expect_true(file.exists(out$plot_paths$correlation_heatmap))
})

test_that("custom corr_method / corr_threshold are honoured", {
  out <- diagnose_training_pools(make_corr_fixture(), group_col = "pool_source",
                                 corr_method = "pearson", corr_threshold = 0.95,
                                 verbose = FALSE, make_plots = FALSE)
  expect_identical(out$correlation$method, "pearson")
  expect_equal(out$correlation$threshold, 0.95)
})
