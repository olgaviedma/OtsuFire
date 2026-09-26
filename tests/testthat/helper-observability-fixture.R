# Shared synthetic fixture for the observability tests.
#
# Lives in a helper- file so BOTH test-observability-wholefire.R and
# test-observability-undetermined-policy.R can use it: testthat sources
# helper-*.R before every test file, but does not share top-level definitions
# across test files.
#
# Synthetic 10x10 grid of 30 m cells on EPSG:3035, everything burnable.
#   ref1 = top-left  quadrant (rows 1-5 x cols 1-5)  -> 25 cells
#   ref2 = top-right quadrant (rows 1-5 x cols 6-10) -> 25 cells
# The observability raster is built band-wise INSIDE ref2 so that its
# temporal_observable_fraction can be dialled to an exact value (k/25), which
# is what T2/T3/T4 need. ref1 is always fully observable.
#
# Rules under test
#   legacy_any         : observable <=> max(DOY) >= required DOY
#   wholefire_fraction : observable <=> TOF >= observability_min_fraction,
#                        plus removal of the excluded geometries from the
#                        evaluation domain.

# ---- fixture -----------------------------------------------------------------

# `ref2_obs_rows` = how many of the 5 rows of ref2 get a POST-fire DOY.
# ref2 spans raster rows 1-5, so k observable rows give TOF = 5*k/25 = k/5.
build_obs_fixture <- function(dir, ref2_obs_rows = 5L,
                              ref2_end_doy = 100,
                              ref2_start_doy = 90,
                              extra_ref_cols = NULL,
                              pred_over_ref2 = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  template <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 0, xmax = 300, ymin = 0, ymax = 300,
    crs = "EPSG:3035"
  )

  burnable <- template
  terra::values(burnable) <- 1
  burnable_path <- file.path(dir, "burnable.tif")
  terra::writeRaster(burnable, burnable_path, overwrite = TRUE, datatype = "FLT4S")

  mask_sf <- sf::st_sf(
    id = 1L,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0, 300, 0, 300, 300, 0, 300, 0, 0), ncol = 2, byrow = TRUE))), crs = 3035)
  )
  mask_path <- file.path(dir, "mask.shp")
  suppressWarnings(sf::st_write(mask_sf, mask_path, quiet = TRUE, delete_dsn = TRUE))

  ref1_poly <- sf::st_polygon(list(matrix(
    c(1, 151, 149, 151, 149, 299, 1, 299, 1, 151), ncol = 2, byrow = TRUE)))
  ref2_poly <- sf::st_polygon(list(matrix(
    c(151, 151, 299, 151, 299, 299, 151, 299, 151, 151), ncol = 2, byrow = TRUE)))

  ref_sf <- sf::st_sf(
    id        = 1:2,
    year      = c(2000L, 2000L),
    end_doy   = c(100, ref2_end_doy),
    start_doy = c(90, ref2_start_doy),
    geometry  = sf::st_sfc(ref1_poly, ref2_poly, crs = 3035)
  )
  if (!is.null(extra_ref_cols)) {
    for (nm in names(extra_ref_cols)) ref_sf[[nm]] <- extra_ref_cols[[nm]]
  }
  ref_path <- file.path(dir, "ref.gpkg")
  suppressWarnings(sf::st_write(ref_sf, ref_path, quiet = TRUE, delete_dsn = TRUE))

  # Observability DOY: 120 everywhere over ref1 (post-fire for required 100),
  # and over ref2 the top `ref2_obs_rows` raster rows get 200 (post-fire for
  # any required DOY used here) while the rest get 50 (pre-fire).
  obs <- template
  terra::values(obs) <- NA_real_
  m <- matrix(NA_real_, nrow = 10, ncol = 10)
  m[1:5, 1:5] <- 120
  if (ref2_obs_rows > 0L) m[seq_len(ref2_obs_rows), 6:10] <- 200
  if (ref2_obs_rows < 5L) m[(ref2_obs_rows + 1L):5L, 6:10] <- 50
  terra::values(obs) <- as.vector(t(m))
  obs_path <- file.path(dir, "obs.tif")
  terra::writeRaster(obs, obs_path, overwrite = TRUE, datatype = "FLT4S")

  # Prediction: the whole top-left quadrant, optionally extended over ref2.
  pred_xmax <- if (pred_over_ref2) 299 else 149
  pred_sf <- sf::st_sf(
    id = 1L,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(1, 151, pred_xmax, 151, pred_xmax, 299, 1, 299, 1, 151),
      ncol = 2, byrow = TRUE))), crs = 3035)
  )
  pred_path <- file.path(dir, "pred.gpkg")
  suppressWarnings(sf::st_write(pred_sf, pred_path, quiet = TRUE, delete_dsn = TRUE))

  # Prediction fully OUTSIDE any reference fire: bottom-left quadrant.
  fp_sf <- sf::st_sf(
    id = 1L,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(1, 1, 149, 1, 149, 149, 1, 149, 1, 1), ncol = 2, byrow = TRUE))), crs = 3035)
  )
  fp_path <- file.path(dir, "pred_outside.gpkg")
  suppressWarnings(sf::st_write(fp_sf, fp_path, quiet = TRUE, delete_dsn = TRUE))

  list(burnable = burnable_path, mask = mask_path, ref = ref_path,
       obs = obs_path, pred = pred_path, pred_outside = fp_path)
}

run_obs <- function(fx, dir, ...) {
  suppressWarnings(suppressMessages(
    validate_fire_maps(
      input_shapefile      = fx$pred,
      ref_shapefile        = fx$ref,
      mask_shapefile       = fx$mask,
      burnable_raster      = fx$burnable,
      year_target          = 2000L,
      validation_dir       = dir,
      observability_raster = fx$obs,
      metrics_type         = "all",
      force_reprocess_ref  = TRUE,
      force_reprocess_pred = TRUE,
      ...
    )
  ))
}
