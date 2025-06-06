# tests/testthat/test-validate_fire_maps.R

test_that("validate_fire_maps runs on synthetic data (no GDAL)", {
  skip_on_cran()
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("data.table")

  tmp_dir <- file.path(tempdir(), "validation_test")
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
    c(200, 200), c(200, 1000), c(1000, 1000), c(1000, 200), c(200, 200)
  )))
  ref_sf <- sf::st_sf(year = 2022, geometry = sf::st_sfc(ref_poly), crs = 32630)
  ref_path <- file.path(tmp_dir, "ref_polygons.shp")
  sf::st_write(ref_sf, ref_path, delete_layer = TRUE, quiet = TRUE)

  det_poly <- sf::st_polygon(list(rbind(
    c(600, 600), c(600, 1400), c(1400, 1400), c(1400, 600), c(600, 600)
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
    use_gdal = FALSE,
    metrics_type = "all",
    class_shape = NULL
  )


  expect_type(result, "list")
  expect_true(all(c("metrics", "polygon_summary") %in% names(result)))
  expect_s3_class(result$metrics, "data.table")
  expect_s3_class(result$polygon_summary, "data.table")
  expect_true(nrow(result$polygon_summary) == 1)
  expect_true(nrow(result$metrics) == 1)

  expect_true(file.exists(file.path(tmp_dir, "VALIDATION", "metrics_summary_2022.csv")))
  expect_true(file.exists(file.path(tmp_dir, "VALIDATION", "polygon_summary_2022.csv")))
})

