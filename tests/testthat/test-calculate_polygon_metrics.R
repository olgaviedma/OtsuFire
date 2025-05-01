test_that("calculate_polygon_metrics works on synthetic shapefile", {
  skip_on_cran()
  skip_if_not_installed("sf")

  tmp_dir <- tempdir()
  shp_path <- file.path(tmp_dir, "fire_polygons_test.shp")

  # Crear poligono sintetico simple con CRS
  poly <- sf::st_polygon(list(rbind(
    c(0, 0), c(0, 100), c(100, 100), c(100, 0), c(0, 0)
  )))
  sf_poly <- sf::st_sf(id = 1, geometry = sf::st_sfc(poly), crs = 32630)
  sf::st_write(sf_poly, shp_path, delete_layer = TRUE, quiet = TRUE)

  # Llamar a la funcion
  result <- calculate_polygon_metrics(
    shapefile_paths = shp_path,
    output_dir = tmp_dir,
    area_min_ha = 0.01,
    bbox_h_min = 10,
    mnbbx_wd_min = 10,
    p_w_ratio_min = 1.0,
    h_w_ratio_min = 0.2
  )

  expect_type(result, "list")
  expect_named(result[[1]], c("metrics", "filtered", "polygons_all", "polygons_filtered"))
  expect_s3_class(result[[1]]$polygons_all, "sf")
})
