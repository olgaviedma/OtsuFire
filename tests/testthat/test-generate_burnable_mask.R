test_that("generate_burnable_mask funciona con datos sinteticos", {
  skip_on_cran()

  # Asegura que las rutas estan disponibles
  Sys.setenv(
    GDALWARP_PATH = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe",
    PYTHON_EXE_PATH = "C:/ProgramData/anaconda3/python.exe",
    GDAL_POLYGONIZE = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"
  )

  gdalwarp <- get_env_tool("GDALWARP_PATH")
  python <- get_env_tool("PYTHON_EXE_PATH")
  gdal_polygonize <- get_env_tool("GDAL_POLYGONIZE")

  skip_if(is.null(gdalwarp) || !file.exists(gdalwarp), "gdalwarp.exe no encontrado")
  skip_if(is.null(python) || !file.exists(python), "Python no disponible")
  skip_if(is.null(gdal_polygonize) || !file.exists(gdal_polygonize), "gdal_polygonize.py no disponible")

  # Crear raster sintetico tipo CORINE
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = "EPSG:3035")
  values(r) <- sample(c(311, 312, 313, 324), size = terra::ncell(r), replace = TRUE)  # codigos CORINE arbitrarios
  corine_path <- file.path(tempdir(), "corine_synthetic.tif")
  terra::writeRaster(r, corine_path, overwrite = TRUE)

  # Crear shapefile sintetico para enmascarar
  tmp_shp <- file.path(tempdir(), "mask_synthetic.shp")
  pts <- terra::vect(matrix(c(500, 500), ncol = 2), type = "points", crs = "EPSG:3035")
  poly <- terra::buffer(pts, width = 600)
  terra::writeVector(poly, tmp_shp, overwrite = TRUE)

  # Carpetas de salida
  out_raster <- file.path(tempdir(), "rasters_test")
  out_vector <- file.path(tempdir(), "shapes_test")
  dir.create(out_raster, showWarnings = FALSE)
  dir.create(out_vector, showWarnings = FALSE)

  # Ejecutar la funcion
  expect_error({
    generate_burnable_mask(
      corine_rasters = list(prueba = corine_path),
      peninsula_shapefile = tmp_shp,
      output_raster_dir = out_raster,
      output_vector_dir = out_vector,
      gdalwarp_path = gdalwarp,
      python_exe = python,
      gdal_polygonize_script = gdal_polygonize,
      reproject = TRUE,
      res = 90,
      to_wgs84 = FALSE,
      vectorize = FALSE
    )
  }, NA)

  # Verificar que se genera el raster de salida
  result_file <- file.path(out_raster, "burneable_mask_binary_prueba_ETRS89.tif")
  expect_true(file.exists(result_file))
})
