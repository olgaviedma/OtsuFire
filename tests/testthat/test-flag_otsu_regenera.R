test_that("flag_otsu_regenera() runs correctly with group_field = NULL (uses burned_id)", {
  skip_if_not_installed("sf")

  # Polígono quemado con columna burned_id
  burned_poly <- sf::st_sf(
    burned_id = 1L,
    CORINE_CLA = "311",
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
    )))),
    crs = 4326
  )

  # Polígono de regeneración que se solapa parcialmente
  regen_poly <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0.5, 0.5), c(1.5, 0.5), c(1.5, 1.5), c(0.5, 1.5), c(0.5, 0.5)
    )))),
    crs = 4326
  )

  # Guardar ambos en archivos temporales
  tmp_burned <- tempfile(fileext = ".geojson")
  tmp_regen <- tempfile(pattern = "regen_P1_thresh200_", fileext = ".geojson")
  sf::st_write(burned_poly, tmp_burned, quiet = TRUE)
  sf::st_write(regen_poly, tmp_regen, quiet = TRUE)

  # Ejecutar función y verificar que no lanza error
  expect_error({
    result <- OtsuFire::flag_otsu_regenera(
      burned_files = tmp_burned,
      regenera_files = tmp_regen,
      min_regen_ratio = c(P1 = 0.01),
      output_dir = tempdir(),
      replace_by_P1 = FALSE,
      save_no_regenera = "none",
      output_format = "geojson"
      # group_field = NULL por defecto, usa burned_id
    )
  }, NA)
})

