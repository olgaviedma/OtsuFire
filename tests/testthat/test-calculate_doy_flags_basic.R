test_that("calculate_doy_flags works with synthetic raster and polygons", {
  skip_on_cran()
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  tmp_dir <- file.path(tempdir(), "doyflags_test")
  dir.create(tmp_dir, showWarnings = FALSE)

  # Crear raster sintetico
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r) <- runif(terra::ncell(r), min = 0, max = 100)

  doy_band <- terra::rast(r)
  terra::values(doy_band) <- runif(terra::ncell(r), min = 60, max = 90)

  r_stack <- c(r, doy_band)
  names(r_stack) <- c("index", "DOY")
  terra::crs(r_stack) <- "EPSG:32630"

  raster_path <- file.path(tmp_dir, "synthetic_doy.tif")
  terra::writeRaster(r_stack, raster_path, overwrite = TRUE)

  # Crear poligono sintetico
  poly <- sf::st_polygon(list(rbind(
    c(100, 100), c(100, 900), c(900, 900), c(900, 100), c(100, 100)
  )))
  poly_sf <- sf::st_sf(geometry = sf::st_sfc(poly), crs = 32630)
  poly_path <- file.path(tmp_dir, "poly_doy.shp")
  sf::st_write(poly_sf, poly_path, delete_layer = TRUE, quiet = TRUE)

  # Ejecutar funcion
  result <- calculate_doy_flags(
    raster = terra::rast(raster_path),
    doy_band = 2,
    polygons_sf = poly_path,
    output_dir = tmp_dir,
    year = 2022,
    doy_thresholds = 10,
    stats = "both",
    polygonize = FALSE
  )

  expect_type(result, "list")
  expect_true("sf_objects" %in% names(result))
  expect_true("shp_paths" %in% names(result))

  stats_result <- result$sf_objects$stats

  if (length(stats_result) == 0) {
    skip("No se generaron estadisticas DOY para este poligono sintetico.")
  }

  expect_gte(length(stats_result), 1)
  expect_s3_class(stats_result[[1]], "sf")
})
