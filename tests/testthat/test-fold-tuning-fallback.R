# Bug 4 (OtsuFire 0.3.0): fold-tuning soft-fallback selection.
#
# The selection loop must prefer the highest-scoring OK candidate over
# any non-OK candidate. Only when no OK candidate exists may it fall
# back to the best non-OK, with a loud warning and a
# `_FOLD_FALLBACK.txt` audit file.

# Emulator of the 0.3.0 selection logic.
.select_best_fold_config <- function(candidates) {
  best_ok <- NULL
  best_nonok <- NULL
  for (ev in candidates) {
    if (isTRUE(ev$ok)) {
      if (is.null(best_ok) || ev$score > best_ok$score) best_ok <- ev
    } else {
      if (is.null(best_nonok) || ev$score > best_nonok$score) best_nonok <- ev
    }
  }
  if (!is.null(best_ok)) {
    list(picked = best_ok, fallback = FALSE)
  } else {
    list(picked = best_nonok, fallback = TRUE)
  }
}

# --- T12 ------------------------------------------------------------
test_that("T12: mixed OK/non-OK candidates -> OK wins even when its score is lower", {
  candidates <- list(
    list(ok = TRUE,  score = 50,   label = "ok_low"),
    list(ok = FALSE, score = 2000, label = "nonok_high"),
    list(ok = TRUE,  score = 60,   label = "ok_better")
  )
  pick <- .select_best_fold_config(candidates)
  expect_false(pick$fallback)
  expect_identical(pick$picked$label, "ok_better")
})

# --- T13 ------------------------------------------------------------
test_that("T13: no OK candidate -> fallback to best non-OK", {
  candidates <- list(
    list(ok = FALSE, score = 100, label = "nonok_low"),
    list(ok = FALSE, score = 250, label = "nonok_best"),
    list(ok = FALSE, score = 200, label = "nonok_mid")
  )
  pick <- .select_best_fold_config(candidates)
  expect_true(pick$fallback)
  expect_identical(pick$picked$label, "nonok_best")
})

test_that("T13b: make_block_folds writes _FOLD_FALLBACK.txt when fallback fires", {
  skip_if_not_installed("sf")

  # 2 burned + 6 unburned points across two well-separated clusters.
  # This is a degenerate setup: with k_candidates 5 and burned units = 2,
  # min_burned_units_per_fold >= 1 is not even possible per fold for k>=3.
  set.seed(42)
  pts <- list(
    sf::st_point(c(0, 0)),     sf::st_point(c(1000, 0)),
    sf::st_point(c(2000, 0)),  sf::st_point(c(0, 1000)),
    sf::st_point(c(50000, 0)), sf::st_point(c(51000, 0)),
    sf::st_point(c(50000, 1000)), sf::st_point(c(51000, 1000))
  )
  geom <- sf::st_sfc(pts, crs = 3035)
  sf_in <- sf::st_sf(
    fire_uid = sprintf("f%02d", seq_along(pts)),
    class    = c("burned", "unburned", "unburned", "unburned",
                 "burned", "unburned", "unburned", "unburned"),
    geometry = geom
  )
  out_dir <- tempfile("fold_fallback_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res <- tryCatch(
    suppressWarnings(suppressMessages(
      make_block_folds <- get("make_block_folds", envir = asNamespace("OtsuFire"))
    )),
    error = function(e) NULL
  )
  fn <- get("make_block_folds", envir = asNamespace("OtsuFire"))

  # Force impossibility: very high min_burned_units_per_fold.
  out <- tryCatch(
    suppressWarnings(suppressMessages(fn(
      train_labelled_sf = sf_in,
      class_col = "class",
      pos_lab = "burned",
      neg_lab = "unburned",
      fire_id_col = "fire_uid",
      split_unit = "polygon",
      block_sizes_m = c(5000, 3000),
      k_candidates = c(5, 4, 3),
      n_repeats = 1L,
      min_burned_units_per_fold = 50L,
      min_pos_blocks_per_fold = 50L,
      out_dir = out_dir,
      target_year = 2025L,
      out_prefix = "test",
      write_blocks_gpkg = FALSE,
      write_folds_csv = FALSE,
      write_train_with_folds_gpkg = FALSE,
      verbose = FALSE
    ))),
    error = function(e) e
  )

  fallback_path <- file.path(out_dir, "_FOLD_FALLBACK.txt")
  expect_true(file.exists(fallback_path),
              info = paste("expected", fallback_path, "to be written"))
})
