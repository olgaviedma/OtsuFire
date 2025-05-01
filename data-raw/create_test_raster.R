# Cargar paquetes necesarios
library(terra)
library(sf)

# Crear carpeta si no existe
if (!dir.exists("inst/extdata")) dir.create("inst/extdata", recursive = TRUE)

# Raster 1: 100x100 pixeles, valores aleatorios
r1 <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin=0, ymax=100)
values(r1) <- runif(ncell(r1), 0, 100)
writeRaster(r1, "inst/extdata/r1.tif", overwrite=TRUE, gdal=TRUE, options=c("COMPRESS=LZW"))

# Raster 2: adyacente al este del anterior
r2 <- rast(nrows=100, ncols=100, xmin=100, xmax=200, ymin=0, ymax=100)
values(r2) <- runif(ncell(r2), 0, 100)
writeRaster(r2, "inst/extdata/r2.tif", overwrite=TRUE, gdal=TRUE, options=c("COMPRESS=LZW"))

# Shapefile de mascara: poligono rectangular que cubre parte de ambos
mask_poly <- st_as_sf(st_sfc(st_polygon(list(rbind(
  c(50, 20), c(150, 20), c(150, 80), c(50, 80), c(50, 20)
)))), crs = 4326)
st_write(mask_poly, "inst/extdata/mask.shp", delete_layer = TRUE)

# El shapefile generara automaticamente .shx, .dbf, etc., necesarios
