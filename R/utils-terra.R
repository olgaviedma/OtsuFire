# Shared terra helpers. Internal only.
#
# Scope (Phase A — stubs only):
# - grid_equal(): tolerant resolution/origin comparison.
# - ensure_raster_crs_and_grid(): project + resample a raster to a template.
# - clean_raster_extremes(): NA-out values outside a sanity window.

grid_equal <- function(r, template, tol = 1e-6) {
  .NotYetImplemented()
}

ensure_raster_crs_and_grid <- function(r, template, name,
                                       method_project = "bilinear",
                                       method_resample = "bilinear") {
  .NotYetImplemented()
}

clean_raster_extremes <- function(r, min_ok = -10000, max_ok = NULL, tag = "raster") {
  .NotYetImplemented()
}
