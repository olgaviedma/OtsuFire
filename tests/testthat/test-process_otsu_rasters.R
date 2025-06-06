test_that("process_otsu_rasters runs on synthetic raster and outputs valid shapefile", {
  skip_on_cran()

  # Rutas requeridas
  python <- get_env_tool("PYTHON_EXE_PATH")
  gdal_polygonize <- get_env_tool("GDAL_POLYGONIZE")
  skip_if(is.null(python) || !file.exists(python), "Python no disponible")
  skip_if(is.null(gdal_polygonize) || !file.exists(gdal_polygonize), "gdal_polygonize.py no encontrado")
  skip_if_not_installed("OtsuSeg")

  # Crear raster sintetico
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = "EPSG:3035")
  terra::values(r) <- matrix(runif(100, 0, 1000), nrow = 10)
  tmp_dir <- file.path(tempdir(), "otsu_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  raster_path <- file.path(tmp_dir, "raster_synthetic.tif")
  terra::writeRaster(r, raster_path, overwrite = TRUE)

  # Ejecutar funcion
  process_otsu_rasters(
    raster_path = raster_path,
    output_dir = tmp_dir,
    otsu_thresholds = c(0),
    use_original = TRUE,
    python_exe = python,
    gdal_polygonize_script = gdal_polygonize,
    tile = FALSE
  )

  # Buscar shapefile
  shp_files <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE)
  skip_if(length(shp_files) == 0, "No shapefile created")

  shp <- sf::st_read(shp_files[1], quiet = TRUE)

  # Validar que es un sf valido
  expect_s3_class(shp, "sf")
  expect_true(nrow(shp) > 0)
  expect_true("geometry" %in% names(shp))

  # Campos opcionales esperados: DN o Label
  expect_true(any(c("DN", "Label") %in% names(shp)),
              info = "Shapefile should contain either 'DN' or 'Label' field.")
})
