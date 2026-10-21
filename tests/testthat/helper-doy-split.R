# Synthetic fixture for the DOY split tests: 90 m grid (EPSG:3035), no real data.

# 60 x 60 cells of 90 m. Band 1 = RBR (dummy), band 2 = DOY: 168 for x < 2700,
# 240 elsewhere. Candidate 1 spans both DOY regions (mixed, 324 ha); candidate 2
# lies entirely in the DOY 240 region (homogeneous, 243 ha).
mk_doy_fixture <- function(doy_left = 168, doy_right = 240) {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 5400,
                   ymin = 0, ymax = 5400, crs = "EPSG:3035", nlyrs = 2)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- cbind(rep(500, nrow(xy)),
                            ifelse(xy[, 1] < 2700, doy_left, doy_right))
  names(r) <- c("rbr", "doy")
  ci <- tempfile(fileext = ".tif")
  terra::writeRaster(r, ci, overwrite = TRUE)
  box <- function(x0, y0, x1, y1) {
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  }
  cand <- sf::st_sf(
    CLUMP_ID = 1:2, YEAR = 2006L, RUN = "t",
    geometry = sf::st_sfc(box(1800, 3600, 3600, 5400),
                          box(3600, 0, 4950, 1800), crs = 3035))
  gp <- file.path(tempfile("run"), "02_REFINE", "BA_2006_REFINE_MERGED.gpkg")
  dir.create(dirname(gp), recursive = TRUE)
  sf::st_write(cand, gp, quiet = TRUE)
  list(ci = ci, cand = cand, gpkg = gp)
}
