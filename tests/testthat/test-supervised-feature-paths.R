# Gate 1B PIECE 2 (2026-06-07): every FEATURE/CANDIDATE input path is a
# cfg$inputs field (no hardcoded absolute paths). PIECE 2 ONLY makes the paths
# configurable + consumed from cfg.
#
# Gate 1B PIECE 4 (2026-06-08): the supervised ecoregion Otsu branch (otsu_mode
# "ecoregion" / "corine_ecoregion") and its `ecoregion_shapefile` cfg input were
# removed; the canonical negative-pool Otsu mode is `burnable_only` ("corine" is
# also supported, ecoregion-free). The former PIECE-2 ecoregion blocks (1-3)
# were deleted; the topo/corine block below is kept. The supervised burnable_only
# / no-ecoregion contract is asserted in test-supervised-config.R.

ns <- asNamespace("OtsuFire")

# --------------------------------------------------------------------------
# The feature extractor's standalone branch reads topo + corine from cfg.
# --------------------------------------------------------------------------
test_that("PIECE 2: feature extractor consumes topo + corine_raster from cfg", {
  src <- paste(deparse(get("extract_supervised_features", envir = ns)),
               collapse = "\n")
  expect_match(src, 'cfg_input_path\\("topo"\\)|sup_input_path\\([^,]*,\\s*"topo"\\)')
  expect_match(src, 'cfg_input_path\\("corine_raster"\\)|sup_input_path\\([^,]*,\\s*"corine_raster"\\)')
})
