# Generar shapefile sintetico para tests

library(sf)

# Crear dos poligonos simples (cuadrados)
poly1 <- st_polygon(list(rbind(
  c(0, 0), c(0, 100), c(100, 100), c(100, 0), c(0, 0)
)))
poly2 <- st_polygon(list(rbind(
  c(200, 0), c(200, 50), c(300, 50), c(300, 0), c(200, 0)
)))

# Crear sf object
polys <- st_sf(
  id = 1:2,
  geometry = st_sfc(poly1, poly2),
  crs = 32630  # UTM Zone 30N (puedes ajustar a tu zona)
)

# Guardar en inst/extdata/
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
st_write(polys, "inst/extdata/fire_polygons_test.shp", delete_dsn = TRUE)
