test_that("mosaic_reproject_resample works with synthetic rasters", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("gdalUtilities")

  # Crear carpeta temporal
  tmp_dir <- file.path(tempdir(), "mosaic_test")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # Crear dos raster tiles sintéticos con CRS proyectado (EPSG:3035)
  make_raster_tile <- function(name, offset_x = 0) {
    r1 <- terra::rast(nrows = 10, ncols = 10,
                      xmin = offset_x, xmax = offset_x + 1000,
                      ymin = 0, ymax = 1000,
                      crs = "EPSG:3035")
    r2 <- terra::rast(r1)

    terra::values(r1) <- runif(raster::ncell(r1), 0, 100)      # RBR
    terra::values(r2) <- sample(150:250, ncell(r2), replace = TRUE)  # DOY

    r_stack <- c(r1, r2)
    out_path <- file.path(tmp_dir, paste0(name, ".tif"))

    # ✅ Cambiar esta línea para que admita -9999 como nodata
    terra::writeRaster(r_stack, out_path, overwrite = TRUE, datatype = "FLT4S")

    return(out_path)
  }


  raster1 <- make_raster_tile("tile1", offset_x = 0)
  raster2 <- make_raster_tile("tile2", offset_x = 1000)

  # Crear shapefile de máscara que cubre ambos raster tiles
  mask_poly <- sf::st_polygon(list(rbind(
    c(0, 0), c(0, 1000), c(2000, 1000), c(2000, 0), c(0, 0)
  )))
  mask_sf <- sf::st_sf(geometry = sf::st_sfc(mask_poly), crs = 3035)
  mask_path <- file.path(tmp_dir, "mask.shp")
  sf::st_write(mask_sf, mask_path, quiet = TRUE)

  # Ejecutar la función con resolución adecuada (e.g., 100 m)
  result <- mosaic_reproject_resample(
    folder_path = tmp_dir,
    mask_path = mask_path,
    year = 2025,
    raster_pattern = "tile*.tif",
    crs_target = "EPSG:3035",
    res_target = 100,
    nodata_value = -9999
  )

  # Comprobaciones
  expect_type(result, "character")
  expect_true(file.exists(result))

  out_raster <- terra::rast(result)
  expect_equal(terra::nlyr(out_raster), 2)
  expect_false(any(is.infinite(terra::values(out_raster))))

  unlink(tmp_dir, recursive = TRUE)
})
