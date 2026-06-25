# --- shared synthetic fixture -------------------------------------------------
make_pool_fixture <- function(n = 300, seed = 1, add_geom = FALSE,
                              add_zerovar = FALSE, add_na = FALSE) {
  set.seed(seed)
  grp <- sample(c("otsu", "random", "high_confidence_keep"), n, TRUE)
  # make rbr_med separable by group so tests are meaningful
  base <- c(otsu = 200, random = 120, high_confidence_keep = 260)[grp]
  df <- data.frame(
    pool_source      = grp,
    fire_uid         = sample(sprintf("F%03d", 1:30), n, TRUE),
    training_label   = sample(c("burned", "unburned"), n, TRUE),
    used_in_training = sample(c(TRUE, FALSE), n, TRUE),
    VISUAL           = sample(c(0, 1), n, TRUE),
    year             = sample(2015:2017, n, TRUE),
    rbr_med          = stats::rnorm(n, base, 40),
    rbr_aw_med       = stats::rnorm(n, base - 20, 35),
    persist_ratio    = stats::runif(n),
    persist_delta    = stats::rnorm(n),
    area_ha          = stats::rexp(n, 0.1),
    doy_iqr          = stats::rpois(n, 20),
    stringsAsFactors = FALSE
  )
  if (add_zerovar) df$flat_feature <- 5
  if (add_na) df$rbr_med[sample(n, max(1, floor(n * 0.1)))] <- NA
  if (add_geom) {
    pts <- lapply(seq_len(n), function(i) sf::st_point(c(runif(1), runif(1))))
    df <- sf::st_sf(df, geometry = sf::st_sfc(pts, crs = 4326))
  }
  df
}

test_that("synthetic fixture helper builds a coherent pool", {
  expect_silent(make_pool_fixture(120))
})

test_that("(1) accepts a data.frame and returns the documented structure", {
  out <- diagnose_training_pools(make_pool_fixture(200), verbose = FALSE,
                                 make_plots = FALSE)
  expect_s3_class(out, "otsufire_pool_diagnostics")
  for (nm in c("summary", "settings", "feature_info", "missingness",
               "zero_variance", "kruskal", "pairwise", "feature_ranking",
               "zscore_fingerprint", "pca", "permanova", "rs_separability",
               "probe_classifier", "plot_paths", "output_dir"))
    expect_true(nm %in% names(out), info = nm)
})

test_that("(2) accepts sf input and drops geometry internally", {
  skip_if_not_installed("sf")
  out <- diagnose_training_pools(make_pool_fixture(200, add_geom = TRUE),
                                 verbose = FALSE, make_plots = FALSE)
  expect_s3_class(out, "otsufire_pool_diagnostics")
  expect_false(any(c("geom", "geometry") %in% out$feature_info$features_used))
})

test_that("(3) auto-detects numeric features", {
  out <- diagnose_training_pools(make_pool_fixture(200), verbose = FALSE,
                                 make_plots = FALSE)
  expect_true(all(c("rbr_med", "persist_ratio", "area_ha") %in%
                    out$feature_info$features_used))
})

test_that("(4) excludes geometry/ID/label/group/VISUAL/used columns", {
  out <- diagnose_training_pools(make_pool_fixture(200), verbose = FALSE,
                                 make_plots = FALSE)
  fu <- out$feature_info$features_used
  expect_false("pool_source" %in% fu)
  expect_false("training_label" %in% fu)
  expect_false("VISUAL" %in% fu)
  expect_false("used_in_training" %in% fu)
  expect_false("year" %in% fu)
  expect_false("fire_uid" %in% fu)
})

test_that("(5) removes zero-variance features and reports them", {
  out <- diagnose_training_pools(make_pool_fixture(200, add_zerovar = TRUE),
                                 verbose = FALSE, make_plots = FALSE)
  expect_true("flat_feature" %in% out$feature_info$zero_variance_removed)
  expect_false("flat_feature" %in% out$feature_info$features_used)
})

test_that("(6) handles NAs without erroring", {
  out <- diagnose_training_pools(make_pool_fixture(200, add_na = TRUE),
                                 verbose = FALSE, make_plots = FALSE)
  expect_s3_class(out, "otsufire_pool_diagnostics")
  expect_true(any(out$missingness$n_missing > 0))
})

test_that("(7) drops small groups from tests but reports them", {
  df <- make_pool_fixture(200)
  df$pool_source[1:5] <- "tiny_group"   # only 5 rows -> below min_n
  out <- diagnose_training_pools(df, min_n_per_group = 10, verbose = FALSE,
                                 make_plots = FALSE)
  expect_true("tiny_group" %in% out$summary$dropped_small_groups)
  expect_false("tiny_group" %in% out$kruskal$feature)  # not a feature
  expect_false("tiny_group" %in% unique(c(out$pairwise$group1, out$pairwise$group2)))
})

test_that("(8) Kruskal-Wallis results present with labelled effect sizes", {
  out <- diagnose_training_pools(make_pool_fixture(300), verbose = FALSE,
                                 make_plots = FALSE)
  expect_true(is.data.frame(out$kruskal))
  expect_true(all(c("feature", "statistic", "p", "p_adj",
                    "epsilon_squared", "eta_squared") %in% names(out$kruskal)))
  # rbr_med is built to separate groups -> should be significant
  rr <- out$kruskal[out$kruskal$feature == "rbr_med", ]
  expect_lt(rr$p_adj, 0.05)
})

test_that("(9) pairwise Wilcoxon returned when run_pairwise=TRUE", {
  out <- diagnose_training_pools(make_pool_fixture(300), run_pairwise = TRUE,
                                 verbose = FALSE, make_plots = FALSE)
  expect_true(is.data.frame(out$pairwise))
  expect_true(all(c("group1", "group2", "feature", "statistic", "p", "p_adj",
                    "median_diff", "direction") %in% names(out$pairwise)))
  out2 <- diagnose_training_pools(make_pool_fixture(300), run_pairwise = FALSE,
                                  verbose = FALSE, make_plots = FALSE)
  expect_null(out2$pairwise)
})

test_that("(10) PCA outputs returned when run_pca=TRUE", {
  out <- diagnose_training_pools(make_pool_fixture(300), run_pca = TRUE,
                                 verbose = FALSE, make_plots = FALSE)
  expect_true(!is.null(out$pca))
  expect_true(all(c("scores", "loadings", "explained", "centroid_distances",
                    "imputation") %in% names(out$pca)))
  expect_identical(out$pca$imputation, "median")
  out2 <- diagnose_training_pools(make_pool_fixture(300), run_pca = FALSE,
                                  verbose = FALSE, make_plots = FALSE)
  expect_null(out2$pca)
})

test_that("(11) writes CSVs when output_dir is set", {
  od <- withr::local_tempdir()
  out <- diagnose_training_pools(make_pool_fixture(300), output_dir = od,
                                 make_plots = FALSE, verbose = FALSE)
  files <- list.files(od)
  expect_true("00_console_summary.txt" %in% files)
  expect_true(any(grepl("01_pool_composition.csv$", files)))
  expect_true(any(grepl("04_kruskal_results.csv$", files)))
})

test_that("(12) writes plots when make_plots=TRUE or skips gracefully", {
  od <- withr::local_tempdir()
  out <- diagnose_training_pools(make_pool_fixture(300), output_dir = od,
                                 make_plots = TRUE, verbose = FALSE)
  # should not error regardless of ggplot2 availability
  expect_s3_class(out, "otsufire_pool_diagnostics")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_true(length(out$plot_paths) >= 1)
    expect_true(dir.exists(file.path(od, "figures")))
  }
})

test_that("(13) does NOT train or hold a final supervised model", {
  out <- diagnose_training_pools(make_pool_fixture(300), run_probe_classifier = FALSE,
                                 verbose = FALSE, make_plots = FALSE)
  # no fitted-model objects anywhere in the return
  expect_null(out$probe_classifier)
  flat <- unlist(out, recursive = TRUE)
  expect_false(any(vapply(out, function(z) inherits(z, c("xgb.Booster", "train")),
                          logical(1))))
  # runs without xgboost being attached
  expect_false("package:xgboost" %in% search())
})

test_that("(14) does NOT mutate the input object", {
  df <- make_pool_fixture(200)
  snap <- df
  invisible(diagnose_training_pools(df, verbose = FALSE, make_plots = FALSE))
  expect_identical(df, snap)
})

test_that("sf input is not mutated either", {
  skip_if_not_installed("sf")
  df <- make_pool_fixture(150, add_geom = TRUE)
  snap <- df
  invisible(diagnose_training_pools(df, verbose = FALSE, make_plots = FALSE))
  expect_identical(df, snap)
})

test_that("optional modules skip gracefully when dependency absent", {
  out <- diagnose_training_pools(make_pool_fixture(300),
                                 run_permanova = TRUE,
                                 run_rs_separability = TRUE,
                                 run_probe_classifier = TRUE,
                                 verbose = FALSE, make_plots = FALSE)
  # permanova: list with available flag whether or not vegan installed
  expect_true(is.list(out$permanova))
  expect_true("available" %in% names(out$permanova))
  expect_true(is.list(out$rs_separability))
  expect_true(is.list(out$probe_classifier))
})

test_that("print method emits the console summary", {
  out <- diagnose_training_pools(make_pool_fixture(200), verbose = FALSE,
                                 make_plots = FALSE)
  expect_output(print(out), "pool separability diagnostic")
})
