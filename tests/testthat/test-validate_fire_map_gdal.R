# tests/testthat/test-validate_fire_maps_gdal.R

test_that("validate_fire_maps runs with GDAL polygonization if available", {
  skip_on_cran()
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("data.table")

  python_path <- "C:/ProgramData/anaconda3/python.exe"  # AJUSTAR
  gdal_polygonize_path <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

  skip_if_not(file.exists(python_path), "Python no encontrado")
  skip_if_not(file.exists(gdal_polygonize_path), "gdal_polygonize.py no encontrado")

  tmp_dir <- file.path(tempdir(), "validation_gdal_test")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000, crs = "EPSG:32630")
  terra::values(r) <- 1
  burnable_path <- file.path(tmp_dir, "burnable.tif")
  terra::writeRaster(r, burnable_path, overwrite = TRUE)

  mask_poly <- sf::st_as_sfc(sf::st_bbox(r))
  mask_sf <- sf::st_sf(geometry = mask_poly, crs = 32630)
  mask_path <- file.path(tmp_dir, "mask.shp")
  sf::st_write(mask_sf, mask_path, delete_layer = TRUE, quiet = TRUE)

  ref_poly <- sf::st_polygon(list(rbind(
    c(300, 300), c(300, 1100), c(1100, 1100), c(1100, 300), c(300, 300)
  )))
  ref_sf <- sf::st_sf(year = 2022, geometry = sf::st_sfc(ref_poly), crs = 32630)
  ref_path <- file.path(tmp_dir, "ref_polygons.shp")
  sf::st_write(ref_sf, ref_path, delete_layer = TRUE, quiet = TRUE)

  det_poly <- sf::st_polygon(list(rbind(
    c(800, 800), c(800, 1500), c(1500, 1500), c(1500, 800), c(800, 800)
  )))
  det_sf <- sf::st_sf(geometry = sf::st_sfc(det_poly), crs = 32630)
  det_path <- file.path(tmp_dir, "detected_fire.shp")
  sf::st_write(det_sf, det_path, delete_layer = TRUE, quiet = TRUE)

  result <- validate_fire_maps(
    input_shapefile = det_path,
    ref_shapefile = ref_path,
    mask_shapefile = mask_path,
    burnable_raster = burnable_path,
    year_target = 2022,
    validation_dir = tmp_dir,
    use_gdal = TRUE,
    python_exe = python_path,
    gdal_polygonize_script = gdal_polygonize_path,
    metrics_type = "all",
    class_shape = NULL
  )


  expect_type(result, "list")
  expect_s3_class(result$metrics, "data.table")
  expect_s3_class(result$polygon_summary, "data.table")
  expect_true(nrow(result$polygon_summary) == 1)

  csv_metrics <- file.path(tmp_dir, "VALIDATION", "metrics_summary_2022.csv")
  csv_polygons <- file.path(tmp_dir, "VALIDATION", "polygon_summary_2022.csv")
  expect_true(file.exists(csv_metrics))
  expect_true(file.exists(csv_polygons))
})
