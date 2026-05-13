# 0.4.0 architectural refactor (Agent H): tests for the whitelist
# filter inside run_dm_oof_pipeline().
#
# Under 0.4.0, the OOF stage applies a whitelist filter to its
# `labelled` and `burned_like` inputs before they reach
# build_design_matrix_patches(). The OOF design bundle's
# `prep$matrix_colnames` must therefore be a subset of
# `.supervised_feature_cols ∪ _isNA companions`.

oof_pipeline_internal <- function() {
  get("run_dm_oof_pipeline", envir = asNamespace("OtsuFire"))
}

whitelist_internal <- function() {
  get(".supervised_feature_cols", envir = asNamespace("OtsuFire"))
}

# Build a synthetic OOF input pair (labelled + burned_like sf
# objects) carrying the 50 whitelist features plus a generous set of
# administrative + residual columns.
make_oof_fixture <- function(n = 30, seed = 23L) {
  set.seed(seed)
  whitelist <- whitelist_internal()
  rect <- function(i) sf::st_polygon(list(matrix(c(
    i, 0, i + 0.4, 0,
    i + 0.4, 0.4, i, 0.4, i, 0
  ), ncol = 2, byrow = TRUE)))

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
    fire_uid                 = sprintf("uid_%03d", seq_len(n)),
    class                    = rep(c("burned", "unburned"), length.out = n),
    source                   = rep(c("burned_truth",
                                       "random_burnable_background"),
                                    length.out = n),
    neg_type                 = rep(c(NA_character_, "background_cell"),
                                    length.out = n),
    poly_id                  = sprintf("p_%03d", seq_len(n)),
    block_id                 = rep(seq_len(5), length.out = n),
    fold_rep1                = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2                = rep(c(2L, 1L, 3L), length.out = n),
    cell_id                  = seq_len(n) + 9000L,
    intersects_deterministic = sample(c(TRUE, FALSE), n, replace = TRUE),
    legacy_decision          = sample(c("keep", "drop"), n, replace = TRUE),
    class_final              = rep(c("burned", "drop"), length.out = n),
    median_rbr               = stats::runif(n) * 1000,
    p_above_keep_q25         = stats::runif(n),
    qa_changed               = sample(0:1, n, replace = TRUE),
    p_oof_mean               = stats::runif(n),
    area_ha                  = stats::runif(n, 0.5, 5),
    n_pix                    = sample(10:200, n, replace = TRUE),
    stringsAsFactors         = FALSE
  )

  full_df <- cbind(admin_df, feat_df)
  geom <- sf::st_sfc(lapply(seq_len(n), rect), crs = 3035)
  labelled <- sf::st_sf(full_df, geometry = geom)

  # burned_like = "deterministic universe" without explicit class
  # column. Mirror the same feature set.
  bl_df <- full_df
  bl_df$class <- NULL
  bl_df$fire_uid <- sprintf("bl_%03d", seq_len(n))
  bl_df$poly_id <- sprintf("bl_p_%03d", seq_len(n))
  geom2 <- sf::st_sfc(lapply(seq_len(n), function(i) rect(i + 100)),
                       crs = 3035)
  burned_like <- sf::st_sf(bl_df, geometry = geom2)

  list(labelled = labelled, burned_like = burned_like)
}

test_that("OOF design bundle columns are within the whitelist", {
  skip_if_not_installed("sf")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  fx <- make_oof_fixture()
  whitelist <- whitelist_internal()
  allowed <- c(whitelist, paste0(whitelist, "_isNA"))

  result_dir <- tempfile("oof_whitelist_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Minimal xgb params; nrounds_max + early_stop kept tiny.
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "logloss",
    eta = 0.1,
    max_depth = 3,
    nthread = 1
  )

  # The OOF wrapper consumes labelled_df (with fold cols) separately.
  labelled_df <- sf::st_drop_geometry(fx$labelled)

  # Match orchestrator usage: pass plain data.frames (geometry
  # stripped) to run_dm_oof_pipeline().
  labelled_plain    <- sf::st_drop_geometry(fx$labelled)
  burned_like_plain <- sf::st_drop_geometry(fx$burned_like)

  fn <- oof_pipeline_internal()
  res <- suppressMessages(suppressWarnings(fn(
    labelled    = labelled_plain,
    burned_like = burned_like_plain,
    labelled_df = labelled_df,
    params      = params,
    result_dir  = result_dir,
    target_year = 2005L,
    nrounds_max = 6L,
    early_stop  = 4L,
    seed_base   = 11L,
    verbose     = 0,
    save_prefix = "oof_smoke",
    overwrite   = TRUE
  )))

  bundle_path <- list.files(file.path(result_dir, "04_MATRIX"),
                              pattern = "_design_bundle\\.rds$",
                              full.names = TRUE)
  expect_true(length(bundle_path) >= 1L)
  bundle <- readRDS(bundle_path[1])

  matrix_colnames <- bundle$prep$matrix_colnames
  expect_true(all(matrix_colnames %in% allowed),
              info = paste("non-whitelist OOF matrix cols:",
                            paste(setdiff(matrix_colnames, allowed),
                                  collapse = ", ")))

  # Specific assertions: known administrative / residual columns
  # do NOT appear in the OOF matrix.
  for (nm in c("neg_type", "fire_uid", "block_id", "fold_rep1",
                "fold_rep2", "median_rbr", "p_above_keep_q25",
                "qa_changed", "p_oof_mean", "area_ha", "n_pix",
                "legacy_decision", "intersects_deterministic")) {
    expect_false(nm %in% matrix_colnames,
                  info = paste("admin/residual leaked into OOF cols:", nm))
  }
})

test_that("drop_regex argument is deprecated in run_dm_oof_pipeline()", {
  fn <- oof_pipeline_internal()
  expect_equal(eval(formals(fn)$drop_regex), character(0))
})
