# Polygon geometric metrics. Internal only (per D3 approval).
#
# Lightweight companion to the heavier engine-side metrics used by the
# scoring stage. The richer engine implementation is invoked through the
# dispatcher when it is actually needed; this helper is only used for
# sanity-level area/perimeter summaries inside input validation.

.calculate_polygon_metrics <- function(polys_sf) {
  if (is.null(polys_sf) || !inherits(polys_sf, "sf")) {
    stop(".calculate_polygon_metrics(): input must be an sf object.",
         call. = FALSE)
  }
  if (nrow(polys_sf) == 0L) {
    return(data.frame(area_ha = numeric(0),
                      perimeter_m = numeric(0),
                      compactness = numeric(0)))
  }
  area_m2 <- as.numeric(sf::st_area(polys_sf))
  perim_m <- as.numeric(sf::st_length(sf::st_cast(sf::st_geometry(polys_sf),
                                                   "MULTILINESTRING",
                                                   warn = FALSE)))
  compact <- ifelse(perim_m > 0,
                    4 * pi * area_m2 / (perim_m^2),
                    NA_real_)
  data.frame(
    area_ha     = area_m2 / 10000,
    perimeter_m = perim_m,
    compactness = compact
  )
}
