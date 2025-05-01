# tests/testthat/test-calculate_doy_flags_polygonize.R

test_that("calculate_doy_flags runs with polygonize = TRUE and external scripts", {
  skip_on_cran()
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("data.table")

  python_path <- "C:/ProgramData/anaconda3/python.exe"  # AJUSTA
  gdal_polygonize_path <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"  # AJUSTA

  skip_if_not(file.exists(python_path), "Python no encontrado")
  skip_if_not(file.exists(gdal_polygonize_path), "gdal_polygonize.py no encontrado")

  library(terra)
  library(sf)

  tmp_dir <- file.path(tempdir(), "doyflags_polygon_test")
  dir.create(tmp_dir, showWarnings = FALSE)

  # Raster sintetico (2 bandas)
  r1 <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  values(r1) <- matrix(sample(1:200, 100, replace = TRUE), ncol = 1)
  r2 <- rast(r1)
  values(r2) <- matrix(sample(150:200, 100, replace = TRUE), ncol = 1)
  r_stack <- c(r1, r2)
  crs(r_stack) <- "EPSG:32630"
    names(r_stack) <- c("index", "DOY")
  raster_path <- file.path(tmp_dir, "synthetic_doy.tif")
  writeRaster(r_stack, raster_path, overwrite = TRUE)

  # Poligono de prueba
  poly <- st_polygon(list(rbind(
    c(100, 100), c(100, 500), c(500, 500), c(500, 100), c(100, 100)
  )))
  poly_sf <- st_sf(geometry = st_sfc(poly), crs = 32630)
  poly_path <- file.path(tmp_dir, "poly_doy.shp")
  st_write(poly_sf, poly_path, delete_layer = TRUE, quiet = TRUE)

  # Ejecutar con polygonize
  result <- calculate_doy_flags(
    raster = rast(raster_path),
    doy_band = 2,
    polygons_sf = poly_path,
    output_dir = tmp_dir,
    year = 2022,
    doy_thresholds = 10,
    stats = "both",
    polygonize = TRUE,
    python_exe = python_path,
    gdal_polygonize_script = gdal_polygonize_path
  )

  expect_type(result, "list")
  expect_true("shp_paths" %in% names(result))
  expect_true(length(result$shp_paths$stats) > 0)

  # Verificar existencia de shapefiles poligonizados
  shapefiles_created <- unlist(result$shp_paths$stats)
  expect_true(all(file.exists(shapefiles_created)))
})
