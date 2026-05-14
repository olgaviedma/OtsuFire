test_that("mask_mosaic_raster works with synthetic rasters", {
  # Configuracion de entorno
  Sys.setenv(
    PYTHON_EXE_PATH = "C:/ProgramData/anaconda3/python.exe",
    GDALWARP_PATH = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe"
  )

  python <- get_env_tool("PYTHON_EXE_PATH")
  gdalwarp_path <- get_env_tool("GDALWARP_PATH")

  # Saltar si herramientas clave no estan disponibles
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("gdalUtilities")
  skip_if(is.null(gdalwarp_path) || !file.exists(gdalwarp_path), "gdalwarp.exe no encontrado")
  skip_if(is.null(python) || !file.exists(python), "Python no disponible")

  tmp_dir <- file.path(tempdir(), "mosaic_test")
  dir.create(tmp_dir, showWarnings = FALSE)

  # Crear dos rasters sinteticos
  r1 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r1) <- runif(100, min = 0, max = 100)
  r2 <- terra::rast(r1)
  terra::values(r2) <- runif(100, min = 0, max = 100)
  terra::crs(r1) <- "EPSG:3035"
  terra::crs(r2) <- "EPSG:3035"

  r1_path <- file.path(tmp_dir, "masked_r1.tif")
  r2_path <- file.path(tmp_dir, "masked_r2.tif")
  terra::writeRaster(r1, r1_path, overwrite = TRUE)
  terra::writeRaster(r2, r2_path, overwrite = TRUE)

  # Poligono mascara que solapa parcialmente
  poly <- sf::st_polygon(list(rbind(
    c(100, 100), c(100, 800), c(800, 800), c(800, 100), c(100, 100)
  )))
  poly_sf <- sf::st_sf(geometry = sf::st_sfc(poly), crs = 3035)
  poly_path <- file.path(tmp_dir, "mask.shp")
  sf::st_write(poly_sf, poly_path, delete_layer = TRUE, quiet = TRUE)

  # Ejecutar funcion sin esperar warnings explicitos
  result <- suppressWarnings(
    mask_mosaic_raster(
      folder_path = tmp_dir,
      mask_shapefile = poly_path,
      gdalwarp_path = gdalwarp_path
    )
  )

  expect_type(result, "character")
  expect_true(file.exists(result) || is.na(result))
})

