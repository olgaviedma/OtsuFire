# I/O helpers and output routing. Internal only.
#
# Scope (Phase A — stubs only):
# - normalize_path_safe(): portable path canonicalization.
# - resolve_crs_sf(): sf/terra CRS resolver from EPSG / WKT / crs object.
# - build_output_routes(): central implementation of OUTPUTS_ROUTES.csv,
#   returning a list of absolute paths keyed by output_name.

normalize_path_safe <- function(p) {
  .NotYetImplemented()
}

resolve_crs_sf <- function(x) {
  .NotYetImplemented()
}

.build_output_routes <- function(config, stage = c("mosaic", "deterministic", "supervised")) {
  .NotYetImplemented()
}
