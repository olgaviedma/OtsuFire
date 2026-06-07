# Item 10 (Gate 1B, 2026-06-07): the REPRODUCIBLE fold fingerprint must NOT
# depend on wall-clock time. Two equivalent fold runs (same inputs, same seed,
# same config) must produce an IDENTICAL `fingerprint`, while the separate
# `created_at` metadata field may differ between the two runs.

ns <- asNamespace("OtsuFire")

mk_fold_sf <- function() {
  set.seed(42)
  pts <- list(
    sf::st_point(c(0, 0)),       sf::st_point(c(1000, 0)),
    sf::st_point(c(2000, 0)),    sf::st_point(c(0, 1000)),
    sf::st_point(c(1000, 1000)), sf::st_point(c(0, 2000)),
    sf::st_point(c(50000, 0)),   sf::st_point(c(51000, 0)),
    sf::st_point(c(50000, 1000)),sf::st_point(c(50000, 2000)),
    sf::st_point(c(51000, 2000)),sf::st_point(c(52000, 1000))
  )
  sf::st_sf(
    fire_uid = sprintf("f%02d", seq_along(pts)),
    class    = rep(c("burned", "unburned"), 6),
    geometry = sf::st_sfc(pts, crs = 3035)
  )
}

run_folds <- function(sf_in) {
  fn <- get("make_block_folds", envir = ns)
  suppressWarnings(suppressMessages(fn(
    train_labelled_sf = sf_in,
    class_col = "class", pos_lab = "burned", neg_lab = "unburned",
    fire_id_col = "fire_uid", split_unit = "fire",
    block_sizes_m = c(20000, 10000), k_candidates = c(3, 2),
    n_repeats = 1L, min_burned_units_per_fold = 1L, min_pos_blocks_per_fold = 1L,
    out_dir = NULL,
    write_blocks_gpkg = FALSE, write_folds_csv = FALSE,
    write_train_with_folds_gpkg = FALSE, verbose = FALSE
  )))
}

test_that("Item 10: two equivalent fold runs have IDENTICAL reproducible fingerprints", {
  skip_if_not_installed("sf")
  sf_in <- mk_fold_sf()

  out1 <- run_folds(sf_in)
  Sys.sleep(0)  # do not rely on time; equivalence must hold regardless
  out2 <- run_folds(sf_in)

  # The reproducible fingerprint is present and structured.
  expect_type(out1$fingerprint, "list")
  expect_true(all(c("text", "checksum") %in% names(out1$fingerprint)))

  # Equivalent runs -> identical reproducible fingerprint (text + checksum).
  expect_identical(out1$fingerprint$checksum, out2$fingerprint$checksum)
  expect_identical(out1$fingerprint$text,     out2$fingerprint$text)

  # The wall-clock time MUST NOT leak into the reproducible fingerprint.
  expect_false(grepl("created_at", out1$fingerprint$text, fixed = TRUE))
  expect_false(grepl(format(Sys.Date()), out1$fingerprint$text, fixed = TRUE))

  # created_at is SEPARATE metadata (not hashed); it exists and is a time.
  expect_s3_class(out1$created_at, "POSIXct")
  expect_s3_class(out2$created_at, "POSIXct")
})

test_that("Item 10: the fingerprint CHANGES when the fold configuration changes", {
  skip_if_not_installed("sf")
  sf_in <- mk_fold_sf()

  base <- run_folds(sf_in)
  # A different seed_base must change the reproducible fingerprint.
  fn <- get("make_block_folds", envir = ns)
  alt <- suppressWarnings(suppressMessages(fn(
    train_labelled_sf = sf_in,
    class_col = "class", pos_lab = "burned", neg_lab = "unburned",
    fire_id_col = "fire_uid", split_unit = "fire",
    block_sizes_m = c(20000, 10000), k_candidates = c(3, 2),
    n_repeats = 1L, min_burned_units_per_fold = 1L, min_pos_blocks_per_fold = 1L,
    seed_base = 7L,
    out_dir = NULL,
    write_blocks_gpkg = FALSE, write_folds_csv = FALSE,
    write_train_with_folds_gpkg = FALSE, verbose = FALSE
  )))

  expect_false(identical(base$fingerprint$text, alt$fingerprint$text))
})

test_that("Item 10: make_fold_fingerprint is wall-clock-free by construction", {
  fp_fn <- get("make_fold_fingerprint", envir = ns)
  df <- data.frame(uid = c("a", "b", "c"), fold_rep1 = c(1L, 2L, 1L),
                   stringsAsFactors = FALSE)
  args <- list(
    params   = list(seed_base = 42L, n_repeats = 1L, split_unit = "fire"),
    selected = list(block_size_m = 5000, k_folds = 3, ok = TRUE),
    train_with_folds = df, fold_cols = "fold_rep1", id_col = "uid"
  )
  a <- do.call(fp_fn, args)
  b <- do.call(fp_fn, args)
  expect_identical(a$checksum, b$checksum)
  expect_identical(a$text, b$text)
  # Row order of the input must not affect the fingerprint (stable id ordering).
  df_shuffled <- df[c(3, 1, 2), , drop = FALSE]
  args2 <- args; args2$train_with_folds <- df_shuffled
  c2 <- do.call(fp_fn, args2)
  expect_identical(a$checksum, c2$checksum)
})
