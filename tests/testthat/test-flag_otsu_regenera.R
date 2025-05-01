test_that("flag_otsu_regenera completes and writes outputs", {
  skip_on_cran()

  out_dir <- file.path(tempdir(), "flag_test")
  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE)

  # Crear geometria valida
  p <- sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))) |> sf::st_sfc(crs = 3035)
  burned_poly <- sf::st_sf(ID = 1, geometry = p)
  regenera_poly <- sf::st_sf(DN = 1, geometry = p)

  # Escribir archivos asegurando que se puedan sobrescribir
  burned_path <- file.path(out_dir, "burned.shp")
  regen_path <- file.path(out_dir, "regenera_P2_thresh100.shp")

  unlink(list.files(out_dir, pattern = "burned.*", full.names = TRUE))
  unlink(list.files(out_dir, pattern = "regenera.*", full.names = TRUE))
  sf::st_write(burned_poly, burned_path, delete_layer = TRUE, quiet = TRUE)
  sf::st_write(regenera_poly, regen_path, delete_layer = TRUE, quiet = TRUE)

  # Ejecutar funcion
  flag_otsu_regenera(
    burned_files = burned_path,
    regenera_files = regen_path,
    min_regen_ratio = 0.01,
    remove_no_regenera = FALSE,
    remove_condition = "any_year",
    output_dir = out_dir
  )

  # Verificar que al menos se haya generado un shapefile de salida
  result_files <- list.files(out_dir, pattern = "_rat.*\\.shp$", full.names = TRUE)
  expect_true(length(result_files) >= 1)

  # Verificar que el shapefile se puede leer correctamente
  shp <- sf::st_read(result_files[1], quiet = TRUE)
  expect_s3_class(shp, "sf")
})
