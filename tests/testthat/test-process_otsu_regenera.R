test_that("process_otsu_regenera runs on synthetic RBR raster", {
  skip_on_cran()

  # Set paths
  python <- get_env_tool("PYTHON_EXE_PATH")
  gdal_polygonize <- get_env_tool("GDAL_POLYGONIZE")
  skip_if(is.null(python) || !file.exists(python), "Python not available")
  skip_if(is.null(gdal_polygonize) || !file.exists(gdal_polygonize), "gdal_polygonize.py not found")

  tmpdir <- file.path(tempdir(), "regenera_test")
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

  # Create synthetic raster
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r) <- matrix(sample(-200:200, 100, replace = TRUE), nrow = 10)
  terra::crs(r) <- "EPSG:3035"
  rbr_path <- file.path(tmpdir, "RBR_P2.tif")
  terra::writeRaster(r, rbr_path, overwrite = TRUE)

  # Run function
  result <- process_otsu_regenera(
    rbr_post = list(P2 = rbr_path),
    output_dir = tmpdir,
    python_exe = python,
    gdal_polygonize_script = gdal_polygonize,
    fire_year = 2000,
    use_fixed_threshold = TRUE,
    fixed_threshold_value = -100
  )

  expect_true("P2" %in% names(result))
  p2_result <- result$P2[[1]]
  expect_true(all(c("raster", "shapefile") %in% names(p2_result)))
  expect_true(file.exists(p2_result$raster))
  expect_true(file.exists(p2_result$shapefile))

  shp <- sf::st_read(p2_result$shapefile, quiet = TRUE)
  expect_s3_class(shp, "sf")
  expect_gt(nrow(shp), 0)
})
