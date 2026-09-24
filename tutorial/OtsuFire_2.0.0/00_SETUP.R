# =============================================================================
# OtsuFire 2.0.0 tutorial - STAGE 0: SETUP
# =============================================================================
#
# WHAT THIS FILE DOES
#   Loads the package and declares every path the other four files use.
#   It runs no analysis. Source it first, then run one stage file at a time.
#
# HOW TO USE THE TUTORIAL
#   source("00_SETUP.R")     <- edit the paths below first
#   then open 01_MOSAIC.R, 02_DETERMINISTIC.R, 03_SUPERVISED.R, 04_VALIDATION.R
#   and run them in that order, block by block.
#
# EVERY PATH BELOW IS A PLACEHOLDER. Replace the "/path/to/..." strings with
# real locations on your machine. Nothing runs until you do.
# =============================================================================


# -----------------------------------------------------------------------------
# 1) Load the package
# -----------------------------------------------------------------------------
# Once OtsuFire 2.0.0 is installed:
library(OtsuFire)

# While developing against a source checkout instead, use:
#   pkgload::load_all("/path/to/OtsuFire_v02_rebuild", quiet = TRUE)

suppressPackageStartupMessages({
  library(sf)
  library(terra)
})

# The workflow is planar (EPSG:3035). Spherical geometry off avoids s2 warnings.
sf::sf_use_s2(FALSE)


# -----------------------------------------------------------------------------
# 2) The year you are mapping
# -----------------------------------------------------------------------------
# Everything else is derived from this. Use an integer (the L suffix).
TARGET_YEAR <- 2020L

# A short, stable name for this product. It becomes a folder level in the
# output tree, so keep it constant across years you want to compare.
RUN_NAME <- "Min_Min"


# -----------------------------------------------------------------------------
# 3) Where your data lives
# -----------------------------------------------------------------------------
DATA_BASE <- "/path/to/1_DATA"          # root of the input data
OUTPUT_DIR <- "/path/to/1_DATA/Results" # root of everything the package writes

# Imagery -------------------------------------------------------------------
# The change index is the burn-severity composite for the fire summer (RBR or
# dNBR). One raster per year, already mosaicked - stage 01 builds it.
COMPOSITE_DIR <- file.path(DATA_BASE, "Imagery", "Composites_90m")

CHANGE_INDEX <- file.path(
  COMPOSITE_DIR, RUN_NAME,
  sprintf("MinMin_%d_mosaic_res90m.tif", TARGET_YEAR)
)

# The delayed change index is the autumn-winter composite. It gives the model
# the "rbr_aw" features, which separate a real scar from a summer artefact.
# Optional: set to NULL if you do not have it.
DELAYED_CHANGE_INDEX <- file.path(
  COMPOSITE_DIR, "Autumn",
  sprintf("mean_mean_%d_mosaic.tif", TARGET_YEAR)
)

# Land cover ------------------------------------------------------------------
# CORINE is five-yearly, so the CORINE year is not the fire year.
CORINE_YEAR <- if (TARGET_YEAR < 2000) 1990L else
               if (TARGET_YEAR < 2006) 2000L else
               if (TARGET_YEAR < 2012) 2006L else
               if (TARGET_YEAR < 2018) 2012L else 2018L

VEGETATION_MAP <- file.path(
  DATA_BASE, "Corine_Masks", sprintf("CLC_%d_peninsula.tif", CORINE_YEAR)
)

# Binary 0/1 raster: which pixels can burn at all. It bounds the search and,
# in the supervised stage, bounds where negatives may be drawn from.
BURNABLE_MASK <- file.path(
  DATA_BASE, "Corine_Masks",
  sprintf("burnable_mask_binary_corine_%d_ETRS89.tif", CORINE_YEAR)
)

# Study area ------------------------------------------------------------------
PENINSULA_SHP <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")
# STUDY_AREA_MASK is the PROCESSING extent for mosaicking. It is NOT the
# validation mask: since Reference V2 the two are different objects.
STUDY_AREA_MASK <- file.path(DATA_BASE, "Borders", "Iberian_Peninsula.shp")

# VALIDATION_MASK is the territory with valid independent reference coverage
# for TARGET_YEAR. Pre-2000 it is Portugal + the Spanish regions with official
# cartography that year; 2006+ it is the whole peninsula.
VALIDATION_MASK <- file.path(
  DATA_BASE, "Mask_StudyArea_FINAL",
  sprintf("mask_3035_%d_final.shp", TARGET_YEAR)
)

# Ecoregions: the deterministic stage runs one Otsu threshold per
# CORINE x ecoregion unit, not one global threshold.
ECOREGION_SHP <- file.path(DATA_BASE, "Ecoregions", "ecoregions_iberia.shp")

# Terrain: a two-band raster, elevation in band 1 and slope in band 2.
TOPO <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

# Active fire detections ------------------------------------------------------
# MODIS/VIIRS hotspots. These only exist from 2000 onward: for earlier years
# set HOTSPOTS <- NULL and the pipeline runs without the hotspot features.
HOTSPOTS <- if (TARGET_YEAR >= 2000L) {
  file.path(DATA_BASE, "Hotspots", sprintf("hotspots_iberia_%d.geojson", TARGET_YEAR))
} else {
  NULL
}

# Previous year's burned area -------------------------------------------------
# Used to penalise this year's patches that sit on last year's scar, which is
# usually the same scar seen twice rather than a new fire.
PREVIOUS_YEAR_BURNED <- file.path(
  DATA_BASE, "Fires", sprintf("burned_%d.shp", TARGET_YEAR - 1L)
)

# External reference (ground truth for validation) -----------------------------
# Reference V2 (authoritative). ERA1 1984-2005: Neves atlas for Portugal +
# official regional cartography for Spain. ERA2 2006-2025: EFFIS raw.
# Use the authoritative temporal layer (population_doy / obs_required_doy),
# NOT start_doy / end_doy, which are kept only for legacy compatibility.
# The burnable domain is NOT applied to this layer: validate_fire_maps()
# intersects VALIDATION_MASK with BURNABLE_MASK itself.
REFERENCE_BURNED_MAP <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano_FINAL",
  sprintf("Effis_CA_%d_maskKeep_summer.gpkg", TARGET_YEAR)
)


# -----------------------------------------------------------------------------
# 4) Check everything is really there
# -----------------------------------------------------------------------------
# Failing here takes one second. Failing inside the pipeline takes ten minutes.
check_inputs <- function() {
  needed <- c(
    CHANGE_INDEX = CHANGE_INDEX,
    VEGETATION_MAP = VEGETATION_MAP,
    BURNABLE_MASK = BURNABLE_MASK,
    PENINSULA_SHP = PENINSULA_SHP,
    STUDY_AREA_MASK = STUDY_AREA_MASK,
    VALIDATION_MASK = VALIDATION_MASK,
    TOPO = TOPO,
    REFERENCE_BURNED_MAP = REFERENCE_BURNED_MAP
  )
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    cat("MISSING INPUTS:\n")
    for (nm in names(missing)) cat(sprintf("  %-22s %s\n", nm, missing[[nm]]))
    stop("Fix the paths in 00_SETUP.R before continuing.", call. = FALSE)
  }
  cat("All required inputs found.\n")
  invisible(TRUE)
}

# Uncomment once you have filled in the real paths:
# check_inputs()

cat(sprintf(
  "Setup loaded. year=%d  corine=%d  run=%s\n", TARGET_YEAR, CORINE_YEAR, RUN_NAME
))
