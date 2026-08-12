# Gate 1B PIECE 3 (2026-06-07): immediate/delayed change-index canonical cfg
# routes. Two canonical routes (immediate_change_index_path /
# delayed_change_index_path) with historical aliases (change_index /
# delayed_change_index). Eliminates the change_index-validated-but-ignored
# situation and the severity_raster_path duplication that shadowed the cfg route.

ns <- asNamespace("OtsuFire")

mk_tmp_tif_cr <- function() {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), f,
                     overwrite = TRUE)
  f
}
mk_tmp_gpkg_cr <- function() {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  f
}
norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

# --------------------------------------------------------------------------
# 1) The canonical route names are honoured and resolve onto the single
#    cfg$inputs field for the route.
# --------------------------------------------------------------------------
test_that("PIECE 3: immediate_change_index_path is honoured (canonical immediate route)", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_cr(); imm <- mk_tmp_tif_cr()
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id,
    immediate_change_index_path = imm, target_year = 2017L
  )
  expect_identical(cfg$inputs$change_index$path, norm(imm))
  resolve <- get(".of_sup_input_path", envir = ns)
  expect_identical(resolve(cfg, "change_index"), norm(imm))
})

test_that("PIECE 3: delayed_change_index_path is honoured (canonical delayed route)", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_cr(); imm <- mk_tmp_tif_cr(); del <- mk_tmp_tif_cr()
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id,
    immediate_change_index_path = imm, delayed_change_index_path = del,
    target_year = 2017L
  )
  expect_identical(cfg$inputs$delayed_change_index$path, norm(del))
})

# --------------------------------------------------------------------------
# 2) The historical aliases still work (NOT silently ignored) and map onto the
#    same cfg fields.
# --------------------------------------------------------------------------
test_that("PIECE 3: historical change_index / delayed_change_index aliases still work", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_cr(); ci <- mk_tmp_tif_cr(); del <- mk_tmp_tif_cr()
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id, change_index = ci,
    delayed_change_index = del, target_year = 2017L
  )
  expect_identical(cfg$inputs$change_index$path, norm(ci))
  expect_identical(cfg$inputs$delayed_change_index$path, norm(del))
})

# --------------------------------------------------------------------------
# 3) Supplying both forms of a route is fine if identical, errors on conflict.
# --------------------------------------------------------------------------
test_that("PIECE 3: same value via both names is accepted; a genuine conflict ERRORS", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  id <- mk_tmp_gpkg_cr(); ci <- mk_tmp_tif_cr(); other <- mk_tmp_tif_cr()

  # identical -> accepted
  cfg <- build_supervised_burned_config(
    run_label = "balanced", internal_decisions = id,
    change_index = ci, immediate_change_index_path = ci, target_year = 2017L
  )
  expect_identical(cfg$inputs$change_index$path, norm(ci))

  # conflicting immediate -> error
  expect_error(
    build_supervised_burned_config(
      run_label = "balanced", internal_decisions = id,
      change_index = ci, immediate_change_index_path = other,
      target_year = 2017L
    ),
    regexp = "Conflicting change-index routes"
  )
  # conflicting delayed -> error
  expect_error(
    build_supervised_burned_config(
      run_label = "balanced", internal_decisions = id, change_index = ci,
      delayed_change_index = ci, delayed_change_index_path = other,
      target_year = 2017L
    ),
    regexp = "Conflicting change-index routes"
  )
})

test_that("PIECE 3: the immediate route is REQUIRED (neither name supplied -> error)", {
  skip_if_not_installed("sf")
  id <- mk_tmp_gpkg_cr()
  expect_error(
    build_supervised_burned_config(
      run_label = "balanced", internal_decisions = id, target_year = 2017L
    ),
    regexp = "change_index.*required|required.*change_index|immediate_change_index_path"
  )
})

# --------------------------------------------------------------------------
# 4) Both unburned builders consume the IMMEDIATE change index from cfg via
#    severity_raster_path (no internal convention reconstruction shadowing the
#    cfg route). The pool builder threads the cfg change index to BOTH.
# --------------------------------------------------------------------------
test_that("PIECE 3: deterministic-decisions builder accepts + consumes severity_raster_path", {
  fn <- get("build_unburned_from_deterministic_decisions", envir = ns)
  expect_true("severity_raster_path" %in% names(formals(fn)))
  src <- paste(deparse(fn), collapse = "\n")
  # The cfg route (severity_raster_path) is consulted before the MinMin convention.
  i_cfg  <- regexpr("severity_raster_path", src)
  i_conv <- regexpr("MinMin_", src, fixed = TRUE)
  expect_gt(i_cfg, 0); expect_gt(i_conv, 0)
  expect_lt(i_cfg, i_conv)
})

test_that("PIECE 3: pool builder threads the cfg immediate change index to BOTH unburned builders", {
  src <- paste(deparse(get("build_supervised_training_pools", envir = ns)),
               collapse = "\n")
  # one_year_tif is the cfg change_index (PIECE 1) and is passed as
  # severity_raster_path to BOTH builders.
  n_sev <- length(gregexpr("severity_raster_path\\s*=\\s*one_year_tif", src)[[1]])
  expect_gte(n_sev, 2L)
})
