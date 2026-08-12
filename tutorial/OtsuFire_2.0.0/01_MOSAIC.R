# =============================================================================
# OtsuFire 2.0.0 tutorial - STAGE 1: MOSAIC
# =============================================================================
#
# GOAL
#   Turn a folder of per-tile severity rasters into ONE change-index raster for
#   the year. Every later stage reads that single file.
#
# WHY IT MATTERS
#   The change index is the only spectral input to the whole workflow. If its
#   NoData handling or its extreme values are wrong, the Otsu thresholds in
#   stage 2 are wrong too, and nothing downstream can recover from it.
#
# WHAT YOU NEED
#   - a folder of tiles for one year, all in the same CRS
#   - a polygon of the study area (used to crop and mask)
#
# RUN 00_SETUP.R FIRST.
# =============================================================================

source("00_SETUP.R")


# -----------------------------------------------------------------------------
# 1) Where the tiles are
# -----------------------------------------------------------------------------
# One folder, one year. The pattern must match the tiles of THIS year only,
# and must not match an already-built mosaic (that would feed the output back
# into the input).
TILES_DIR <- file.path(COMPOSITE_DIR, RUN_NAME, "tiles")

TILE_PATTERN <- sprintf("MinMin_%d_tile*.tif", TARGET_YEAR)

MOSAIC_OUT_DIR <- file.path(COMPOSITE_DIR, RUN_NAME)


# -----------------------------------------------------------------------------
# 2) Build the mosaic
# -----------------------------------------------------------------------------
# change_index_mosaic() crops each tile to the mask, merges them, and writes a
# single GeoTIFF. It also returns a summary you should read before moving on.
#
# The arguments that matter, and why:
#
#   nodata_value = -9999
#       The value that means "no observation". Tiles often carry it as a plain
#       number; declaring it here stops it being averaged into real data.
#
#   cap_below = TRUE, lower_cap = -1000
#       Severity composites can carry huge negative outliers from clouds and
#       water. Capping them keeps the histogram usable: Otsu splits a
#       histogram, so a handful of extreme values can move the threshold for
#       the entire ecoregion.
#
#   method = "merge"
#       Writes a real merged raster. "vrt" builds a virtual mosaic instead,
#       which is faster and smaller but depends on the tiles staying put.

mosaic_res <- change_index_mosaic(
  folder_path    = TILES_DIR,
  mask           = STUDY_AREA_MASK,   # path, sf, or SpatVector all work
  year           = TARGET_YEAR,
  raster_pattern = TILE_PATTERN,
  output_dir     = MOSAIC_OUT_DIR,
  nodata_value   = -9999,
  cap_below      = TRUE,
  lower_cap      = -1000,
  method         = "merge",
  overwrite      = FALSE
)

# What you get back:
mosaic_res$mosaic_path     # the file that stage 2 will read
mosaic_res$tile_inventory  # which tiles went in
mosaic_res$method_used
mosaic_res$mask_summary


# -----------------------------------------------------------------------------
# 3) Look at the result before trusting it
# -----------------------------------------------------------------------------
# Two minutes here saves a whole wasted run.

r <- terra::rast(mosaic_res$mosaic_path)

terra::plot(r[[1]], main = sprintf("Change index %d", TARGET_YEAR))

# Range: a plausible RBR is roughly -500 to 1500. A max in the millions means
# NoData leaked into the data.
terra::global(r, c("min", "max"), na.rm = TRUE)

# Non-finite cells: Inf / -Inf / NaN. These must be zero. If they are not,
# clean them before continuing (see block 4).
terra::global(!is.finite(r), "sum")


# -----------------------------------------------------------------------------
# 4) Only if the check above found non-finite values
# -----------------------------------------------------------------------------
# Replace Inf/NaN with NA and write a clean copy. Then point CHANGE_INDEX in
# 00_SETUP.R at this cleaned file.

if (FALSE) {   # set to TRUE only if you need it
  r_clean <- terra::ifel(is.finite(r), r, NA)

  clean_path <- sub("\\.tif$", "_clean.tif", mosaic_res$mosaic_path)

  terra::writeRaster(
    r_clean, clean_path, overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      NAflag   = -9999,
      gdal     = c("TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES")
    )
  )
  cat("Cleaned mosaic written to:", clean_path, "\n")
}


# -----------------------------------------------------------------------------
# 5) The simpler alternative
# -----------------------------------------------------------------------------
# mosaic_from_tiles() does the same merging with fewer options. Use it when you
# just want the tiles joined and masked, with no year logic and no output
# naming rules.

if (FALSE) {
  mosaic_from_tiles(
    folder_path    = TILES_DIR,
    mask_path      = STUDY_AREA_MASK,
    raster_pattern = TILE_PATTERN,
    nodata_value   = -9999,
    cap_below      = TRUE,
    lower_cap      = -1000
  )
}


# -----------------------------------------------------------------------------
# WHAT YOU HAVE NOW
# -----------------------------------------------------------------------------
#   One change-index raster for the year, cropped to the study area, with
#   NoData declared and outliers capped.
#
#   Repeat this stage for the autumn-winter composite as well if you have it:
#   the supervised stage uses it as the "delayed change index".
#
# NEXT: 02_DETERMINISTIC.R
# =============================================================================
