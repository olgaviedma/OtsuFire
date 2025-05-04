test_that("mask_to_mosaic applies binary mask and optional clip correctly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  tmp_dir <- file.path(tempdir(), "mask_test")
  dir.create(tmp_dir, showWarnings = FALSE)

  # Crear un raster de mosaico sintético con 2 bandas (rbr, doy)
  r1 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = "EPSG:3035")
  terra::values(r1) <- runif(100, 0, 100)
  r2 <- terra::rast(r1); terra::values(r2) <- sample(1:365, 100, replace = TRUE)
  mosaic <- c(r1, r2)
  names(mosaic) <- c("rbr", "doy")
  mosaic_path <- file.path(tmp_dir, "synthetic_mosaic.tif")
  terra::writeRaster(mosaic, mosaic_path, overwrite = TRUE)

  # Crear un raster máscara binaria (mitad válida, mitad no válida)
  mask <- terra::rast(r1)
  terra::values(mask) <- c(rep(1, 50), rep(0, 50))
  mask_path <- file.path(tmp_dir, "binary_mask.tif")
  terra::writeRaster(mask, mask_path, overwrite = TRUE)

  # Crear shapefile de recorte (opcional)
  poly <- sf::st_polygon(list(rbind(
    c(100, 100), c(100, 800), c(800, 800), c(800, 100), c(100, 100)
  )))
  clip_sf <- sf::st_sf(geometry = sf::st_sfc(poly), crs = 3035)
  clip_path <- file.path(tmp_dir, "clip.shp")
  sf::st_write(clip_sf, clip_path, delete_layer = TRUE, quiet = TRUE)

  # Ejecutar la función
  result_path <- mask_to_mosaic(
    mosaic_path = mosaic_path,
    mask_raster_path = mask_path,
    shapefile_clip = clip_path
  )

  # Comprobaciones
  expect_type(result_path, "character")
  expect_true(file.exists(result_path))

  result_rast <- terra::rast(result_path)
  expect_equal(terra::nlyr(result_rast), 2)
  expect_equal(names(result_rast), c("rbr", "doy"))
  expect_true(any(is.na(terra::values(result_rast))))
})
