# Shared sf helpers. Internal only.
#
# Scope (Phase A — stubs only):
# - safe_make_valid_polys(): validate, drop empties, cast to MULTIPOLYGON.
# - delete_shapefile_bundle(): remove .shp/.shx/.dbf/.prj/... siblings.
# - make_gpkg_safe(): sanitize column names/list-cols for GPKG writing.
# - assert_schema(): validate returned sf schemas against DATA_MODEL CSVs.

safe_make_valid_polys <- function(x) {
  .NotYetImplemented()
}

delete_shapefile_bundle <- function(shp_path) {
  .NotYetImplemented()
}

make_gpkg_safe <- function(x, sep = "|") {
  .NotYetImplemented()
}

.assert_schema <- function(x, expected_cols, name = "object") {
  .NotYetImplemented()
}
