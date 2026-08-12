# Bug 6 (OtsuFire 0.3.0): `make_block_folds()` must accept sf inputs
# whose active geometry column is named "geom" (e.g., post GPKG
# round-trip).

# --- T15 ------------------------------------------------------------
test_that("T15: GPKG with active geom column 'geom' works in make_block_folds", {
  skip_if_not_installed("sf")

  set.seed(42)
  pts <- list(
    sf::st_point(c(0, 0)),     sf::st_point(c(1000, 0)),
    sf::st_point(c(2000, 0)),  sf::st_point(c(0, 1000)),
    sf::st_point(c(1000, 1000)), sf::st_point(c(0, 2000)),
    sf::st_point(c(50000, 0)), sf::st_point(c(51000, 0)),
    sf::st_point(c(50000, 1000)), sf::st_point(c(50000, 2000)),
    sf::st_point(c(51000, 2000)), sf::st_point(c(52000, 1000))
  )
  geom_sfc <- sf::st_sfc(pts, crs = 3035)
  sf_in <- sf::st_sf(
    fire_uid = sprintf("f%02d", seq_along(pts)),
    class    = c(rep(c("burned", "unburned"), 3),
                 rep(c("burned", "unburned"), 3)),
    geometry = geom_sfc
  )

  # Round-trip via GPKG to obtain an sf object whose active geom column
  # is named "geom" (typical of GeoPackage reads).
  f <- tempfile(fileext = ".gpkg")
  sf::st_write(sf_in, f, layer = "train", quiet = TRUE, delete_dsn = TRUE)
  sf_geom <- sf::read_sf(f, layer = "train")
  active_name <- attr(sf_geom, "sf_column")
  # If the driver kept it as "geom", great; otherwise force the name.
  if (!identical(active_name, "geom")) {
    names(sf_geom)[names(sf_geom) == active_name] <- "geom"
    attr(sf_geom, "sf_column") <- "geom"
  }
  expect_identical(attr(sf_geom, "sf_column"), "geom")

  fn <- get("make_block_folds", envir = asNamespace("OtsuFire"))
  out_dir <- tempfile("mbf_geom_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out <- tryCatch(
    suppressWarnings(suppressMessages(fn(
      train_labelled_sf = sf_geom,
      class_col = "class",
      pos_lab = "burned",
      neg_lab = "unburned",
      fire_id_col = "fire_uid",
      split_unit = "fire",
      block_sizes_m = c(20000, 10000),
      k_candidates = c(3, 2),
      n_repeats = 1L,
      min_burned_units_per_fold = 1L,
      min_pos_blocks_per_fold = 1L,
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
  expect_false(inherits(out, "error"),
               info = paste("make_block_folds errored on geom-named input:",
                            if (inherits(out, "error")) conditionMessage(out)))
})
