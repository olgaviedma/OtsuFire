# Gate 1B piece 5 / C3 (2026-06-08): fail-fast on missing REQUIRED supervised
# inputs / output routes. The legacy-unburned pipeline must STOP with a clear,
# actionable error when the change-index/severity raster or the output root is
# absent, instead of silently globbing the filesystem for a convention candidate
# raster or inventing an output directory under data_base/Results/...
#
# These tests prove:
#   (a) a missing REQUIRED input/route triggers a clear fail-fast error BEFORE
#       any heavy raster compute (the error fires at the top-of-function guards,
#       before process_otsu_rasters_/terra::rast are reached);
#   (b) the canonical, fully-supplied call path still passes the guards (no
#       regression in the fail-fast checks);
#   (c) no Sys.glob/list.files/first-existing candidate search remains on the
#       supervised required-input (severity raster) path, asserted by code
#       inspection of the function body;
#   (d) genuinely-optional inputs are unaffected by these guards.

ns_failfast_c3 <- function() asNamespace("OtsuFire")

mk_tmp_tif_c3 <- function() {
  skip_if_not_installed("terra")
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mk_tmp_gpkg_c3 <- function() {
  skip_if_not_installed("sf")
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
    crs = 3035
  )
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f,
               layer = "internal_decisions", quiet = TRUE, delete_dsn = TRUE)
  f
}

# ---------------------------------------------------------------------------
# (a) Missing REQUIRED severity / change-index raster -> fail-fast
# ---------------------------------------------------------------------------
test_that("C3: missing severity_raster_path fails fast (no candidate-search guess)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  fn <- get("build_otsu_negative_pipeline", envir = ns_failfast_c3())

  data_base <- tempfile("c3_db_")
  dir.create(data_base, recursive = TRUE, showWarnings = FALSE)

  expect_error(
    fn(
      target_year             = 2017L,
      scenario_name           = "balanced",
      data_base               = data_base,
      result_name             = "Min_Min",
      composite_base          = data_base,
      severity_raster_path    = NULL,           # REQUIRED input absent
      internal_decisions_path = mk_tmp_gpkg_c3(),
      out_root_dir            = tempfile("c3_out_"),
      verbose                 = FALSE
    ),
    regexp = "requires 'severity_raster_path'|change_index"
  )
})

# ---------------------------------------------------------------------------
# (a) Missing REQUIRED output root -> fail-fast (does NOT invent a directory)
# ---------------------------------------------------------------------------
test_that("C3: missing out_root_dir fails fast before heavy compute", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  fn <- get("build_otsu_negative_pipeline", envir = ns_failfast_c3())

  data_base <- tempfile("c3_db2_")
  dir.create(data_base, recursive = TRUE, showWarnings = FALSE)

  sev    <- mk_tmp_tif_c3()    # valid explicit change-index raster
  burn   <- mk_tmp_tif_c3()    # valid explicit burnable mask
  intern <- mk_tmp_gpkg_c3()   # valid internal decisions

  expect_error(
    fn(
      target_year             = 2017L,
      scenario_name           = "balanced",
      data_base               = data_base,
      result_name             = "Min_Min",
      composite_base          = data_base,
      severity_raster_path    = sev,
      internal_decisions_path = intern,
      burnable_mask_path      = burn,
      otsu_mode               = "burnable_only",  # skip corine/peninsula checks
      out_root_dir            = NULL,             # REQUIRED route absent
      verbose                 = FALSE
    ),
    regexp = "requires 'out_root_dir'|output_routes"
  )

  # The fail-fast must NOT have fabricated the historical convention directory.
  invented <- file.path(data_base, "Results", "2017", "Min_Min",
                        "SUPERVISED", "balanced", "_OTSU_NEGATIVE")
  expect_false(dir.exists(invented))
})

# ---------------------------------------------------------------------------
# (c) No candidate-search machinery remains on the required-input path
# ---------------------------------------------------------------------------
test_that("C3: legacy pipeline body has no severity candidate-search / glob", {
  fn <- get("build_otsu_negative_pipeline", envir = ns_failfast_c3())
  src <- paste(deparse(fn), collapse = "\n")

  # The first-existing convention search helper must not be invoked anymore.
  expect_no_match(src, "resolve_first_existing_path_otsu_negative", fixed = TRUE)
  # No filesystem discovery of an input on this path.
  expect_no_match(src, "Sys.glob", fixed = TRUE)
  expect_no_match(src, "list.files", fixed = TRUE)
  # The DOY convention candidates that the old glob walked are gone.
  expect_no_match(src, "DOY_", fixed = TRUE)

  # The helper itself was removed from the package namespace.
  expect_false(exists("resolve_first_existing_path_otsu_negative",
                      envir = ns_failfast_c3(), inherits = FALSE))
})

# ---------------------------------------------------------------------------
# (b) Canonical, fully-supplied call still passes the C3 guards (no regression)
# ---------------------------------------------------------------------------
test_that("C3: fully-supplied required inputs pass the fail-fast guards", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  fn <- get("build_otsu_negative_pipeline", envir = ns_failfast_c3())

  data_base <- tempfile("c3_db3_")
  dir.create(data_base, recursive = TRUE, showWarnings = FALSE)

  sev     <- mk_tmp_tif_c3()
  burn    <- mk_tmp_tif_c3()
  intern  <- mk_tmp_gpkg_c3()
  out_root <- tempfile("c3_out3_")

  # With every REQUIRED input/route supplied, the function must get PAST the C3
  # guards. It will still stop later inside the OTSU/polygonize stages (no real
  # tooling in the test env), so we only assert that the error is NOT one of the
  # C3 fail-fast messages -> the guards did not regress / fire spuriously.
  err <- tryCatch(
    {
      fn(
        target_year             = 2017L,
        scenario_name           = "balanced",
        data_base               = data_base,
        result_name             = "Min_Min",
        composite_base          = data_base,
        severity_raster_path    = sev,
        internal_decisions_path = intern,
        burnable_mask_path      = burn,
        otsu_mode               = "burnable_only",
        out_root_dir            = out_root,
        reuse_existing          = FALSE,
        verbose                 = FALSE
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  # The C3 guards must not be the cause of failure.
  if (!is.null(err)) {
    expect_no_match(err, "requires 'severity_raster_path'")
    expect_no_match(err, "requires 'out_root_dir'")
  }
  # The output root WAS created from the supplied route (not invented).
  expect_true(dir.exists(out_root))
})

# ---------------------------------------------------------------------------
# (d) Optional inputs stay optional: omitting the burnable mask falls back to
#     the documented standalone convention (NOT a fail-fast), so the C3 guards
#     do not over-reach into legitimately-optional behaviour.
# ---------------------------------------------------------------------------
test_that("C3: omitting optional burnable_mask_path is not a C3 fail-fast", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  fn <- get("build_otsu_negative_pipeline", envir = ns_failfast_c3())

  data_base <- tempfile("c3_db4_")
  dir.create(data_base, recursive = TRUE, showWarnings = FALSE)

  sev     <- mk_tmp_tif_c3()
  intern  <- mk_tmp_gpkg_c3()
  out_root <- tempfile("c3_out4_")

  # No burnable_mask_path -> the function builds the §N+25 convention path
  # (data_base/Corine_Masks/...). That file does not exist here, so the run
  # stops on the burnable-mask existence check -- NOT on a C3 required-input/
  # output-route guard.
  err <- tryCatch(
    {
      fn(
        target_year             = 2017L,
        scenario_name           = "balanced",
        data_base               = data_base,
        result_name             = "Min_Min",
        composite_base          = data_base,
        severity_raster_path    = sev,
        internal_decisions_path = intern,
        # burnable_mask_path omitted -> convention fallback (optional behaviour)
        otsu_mode               = "burnable_only",
        out_root_dir            = out_root,
        verbose                 = FALSE
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  expect_false(is.null(err))
  # It is the convention-path existence check that trips, not a C3 guard.
  expect_no_match(err, "requires 'severity_raster_path'")
  expect_no_match(err, "requires 'out_root_dir'")
})
