# Gate 1B PIECE 1 (2026-06-07): supervised input inventory in cfg$inputs.
#
# Every filesystem input the supervised pipeline consumes must be (a) a declared,
# VALIDATED field of cfg$inputs and (b) CONSUMED from cfg$inputs end-to-end — never
# reconstructed by filename convention behind the cfg's back. These tests assert:
#   * the single shared accessor .of_sup_input_path() honours the cfg path and
#     returns NULL for absent / in-memory / NULL-config inputs;
#   * the MAIN change_index and the OPTIONAL hotspots input are CONSUMED from
#     cfg$inputs by the orchestrator, the pool builder and the feature extractor
#     (source-inspection: each site routes through .of_sup_input_path / a cfg
#     lookup and the convention path is only a %||% fallback, never the primary);
#   * a pre-MODIS / no-hotspots config (hotspots = NULL) validates and proceeds.

ns <- asNamespace("OtsuFire")

mk_tmp_tif_inv <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}
mk_tmp_gpkg_inv <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}

# --------------------------------------------------------------------------
# 1) The shared cfg-input accessor honours the cfg path (sentinel) and is
#    NULL-safe for absent / in-memory / NULL-config inputs.
# --------------------------------------------------------------------------
test_that("PIECE 1: .of_sup_input_path returns the cfg path for a path spec", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  resolve <- get(".of_sup_input_path", envir = ns)

  sentinel_ci <- mk_tmp_tif_inv()
  sentinel_hs <- mk_tmp_gpkg_inv()
  id          <- mk_tmp_gpkg_inv()

  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id,
    change_index = sentinel_ci, hotspots = sentinel_hs, target_year = 2017L
  )

  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  # The accessor returns EXACTLY the sentinel cfg path, not any convention path.
  expect_identical(resolve(cfg, "change_index"), norm(sentinel_ci))
  expect_identical(resolve(cfg, "hotspots"),     norm(sentinel_hs))
})

test_that("PIECE 1: .of_sup_input_path is NULL-safe (absent / in-memory / NULL cfg)", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  resolve <- get(".of_sup_input_path", envir = ns)

  id <- mk_tmp_gpkg_inv(); ci <- mk_tmp_tif_inv()
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci,
    target_year = 2017L  # no hotspots -> optional input absent
  )
  expect_null(resolve(cfg, "hotspots"))          # absent optional input
  expect_null(resolve(NULL, "change_index"))     # NULL config
  expect_null(resolve(cfg, "not_an_input"))      # unknown key

  # In-memory change_index spec -> no on-disk path -> accessor returns NULL.
  ci_mem <- terra::rast(ncol = 4, nrow = 4, vals = 1:16)
  cfg2 <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci_mem,
    target_year = 2017L
  )
  expect_null(resolve(cfg2, "change_index"))
})

# --------------------------------------------------------------------------
# 2) Source-inspection: the three consumption sites read change_index and
#    hotspots from cfg$inputs (via .of_sup_input_path / a cfg lookup) and the
#    MinMin / Hotspots convention path is only a %||% fallback. This guards
#    against a regression that re-introduces an unconditional convention
#    reconstruction (the KNOWN defect: change_index validated-but-ignored).
# --------------------------------------------------------------------------
test_that("PIECE 1: orchestrator consumes change_index + hotspots from cfg$inputs", {
  lines <- deparse(get("run_supervised_pipeline", envir = ns))
  src <- paste(lines, collapse = "\n")

  # change_index is resolved from the cfg accessor, with the convention only as
  # the %||% fallback (cfg lookup appears BEFORE the MinMin convention).
  i_ci   <- regexpr('cfg_input_path\\("change_index"\\)|sup_input_path\\([^,]*,\\s*"change_index"\\)', src)
  i_conv <- regexpr("MinMin_", src, fixed = TRUE)
  expect_gt(i_ci, 0)
  expect_gt(i_conv, 0)
  expect_lt(i_ci, i_conv)   # cfg read precedes convention fallback

  # hotspots is resolved from the cfg accessor too (no unconditional
  # data_base/Hotspots reconstruction as the sole source).
  expect_match(src, 'cfg_input_path\\("hotspots"\\)|sup_input_path\\([^,]*,\\s*"hotspots"\\)')
})

test_that("PIECE 1: pool builder consumes change_index from cfg$inputs", {
  lines <- deparse(get("build_supervised_training_pools", envir = ns))
  src <- paste(lines, collapse = "\n")
  i_ci   <- regexpr('cfg_input_path\\("change_index"\\)|sup_input_path\\([^,]*,\\s*"change_index"\\)', src)
  i_conv <- regexpr("MinMin_", src, fixed = TRUE)
  expect_gt(i_ci, 0)
  expect_gt(i_conv, 0)
  expect_lt(i_ci, i_conv)
})

test_that("PIECE 1: feature extractor consumes change_index + hotspots from cfg$inputs", {
  lines <- deparse(get("extract_supervised_features", envir = ns))
  src <- paste(lines, collapse = "\n")
  expect_match(src, 'cfg_input_path\\("change_index"\\)|sup_input_path\\([^,]*,\\s*"change_index"\\)')
  expect_match(src, 'sup_input_path\\([^,]*,\\s*"hotspots"\\)|cfg_input_path\\("hotspots"\\)')
})

# --------------------------------------------------------------------------
# 3) Pre-MODIS / no-hotspots: hotspots is OPTIONAL (validated allow_null) so a
#    config with hotspots = NULL constructs and the optional input stays NULL.
#    This is the load-bearing "works for years WITHOUT hotspots" guarantee.
# --------------------------------------------------------------------------
test_that("PIECE 1: pre-MODIS config validates + proceeds with hotspots = NULL", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")

  id <- mk_tmp_gpkg_inv(); ci <- mk_tmp_tif_inv()
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci,
    target_year = 1985L  # pre-MODIS year: no hotspots layer exists
  )
  # hotspots is a declared cfg$inputs field, validated-as-optional -> NULL.
  expect_true("hotspots" %in% names(cfg$inputs))
  expect_null(cfg$inputs$hotspots)
  # The shared accessor reports NULL (absent), which every stage treats as
  # "no hotspots layer" (use_hotspots disabled) rather than erroring.
  resolve <- get(".of_sup_input_path", envir = ns)
  expect_null(resolve(cfg, "hotspots"))
})
