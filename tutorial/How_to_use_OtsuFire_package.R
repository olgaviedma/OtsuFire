# =============================================================================
# How to use the OtsuFire package
# =============================================================================
#
# Contents
# --------
#   1. Introduction: installing and loading the package
#   2. Unsupervised stage
#       2.1 Year with active hotspot support
#
# =============================================================================


# =============================================================================
# 1. INTRODUCTION: INSTALLING AND LOADING THE PACKAGE
# =============================================================================
#
# OtsuFire is distributed through the Comprehensive R Archive Network
# (CRAN) and can be installed in the same way as any other R package.
# The recommended workflow is to install it once from the official
# source, attach it at the start of each session, and then call its
# functions directly.
#
# To install the package, run the following inside an R session
# connected to the internet:
#
#   install.packages("OtsuFire")
#
# This command pulls the latest released version from CRAN, together
# with its mandatory dependencies, and places the package in the
# default user library.
#
# Once installed, the package is loaded for use with:
#
#   library(OtsuFire)
#
# Two ancillary R packages are used throughout the tutorial alongside
# OtsuFire. The `sf` package provides vector data handling, while
# `terra` provides raster handling. Both are loaded explicitly at the
# top of each script:
#
#   library(sf)
#   library(terra)
#
# Some processing steps in the deterministic workflow delegate raster
# and vector operations to external GDAL utilities (gdal_polygonize,
# gdalwarp, ogr2ogr) and to a Python interpreter. These tools must be
# available on the host system; OtsuFire is then told where to find
# them through configuration paths set later in this tutorial.
#
# Once OtsuFire is installed and loaded, the package is ready to use.
# The next sections illustrate its workflow through worked examples,
# starting with the unsupervised stage.


# =============================================================================
# 2. UNSUPERVISED STAGE
# =============================================================================
#
# The unsupervised stage of OtsuFire converts an annual spectral-change
# composite into an auditable layer of candidate burned patches. It is
# composed of three operational steps: building a configuration object
# that bundles inputs and parameters, detecting candidate patches
# through stratified Otsu thresholding and seed-and-grow segmentation,
# and scoring those candidates through rule-based filtering into the
# canonical keep / review / drop decision layer.
#
# The behaviour of the stage depends on the data available for the
# target year. In particular, the second rule-based filter relies on
# active-fire hotspot detections, which are only available from 1995
# onwards. The tutorial therefore separates two cases: years with
# active hotspot support, where the full rule-based system can be
# applied, and earlier years, for which the workflow records the
# hotspot filter as not applied rather than as failed. The first case
# is illustrated below; earlier years are covered in a later
# subsection.


# -----------------------------------------------------------------------------
# 2.1 Year with active hotspot support
# -----------------------------------------------------------------------------
#
# This worked example runs the unsupervised stage for a target year
# in which active-fire hotspot detections are available. The example
# uses 2025 as the target year and the balanced detection preset, but
# the same script applies to any year from 1995 onwards by changing
# the two user settings at the top of the user-settings block.
#
# The script is presented here as a self-contained draft. Subsequent
# revisions of this tutorial will execute it block by block and
# discuss the intermediate outputs in detail.


# ------------------------------------------------------------
# Package load
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(OtsuFire)
  library(sf)
  library(terra)
})


# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------
# Choose one target year and one detection preset.
# Everything else below is derived from these two values.

base_dir <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"

target_year <- 2025L

selected_preset <- "balanced"

# Typical options:
#   "permissive"
#   "balanced"
#   "conservative"
#   "very_conservative"


# ------------------------------------------------------------
# TOOL PATHS
# ------------------------------------------------------------
# These paths are needed because the deterministic workflow delegates
# some raster and vector operations to external GDAL / Python tools.

python_exe <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/python.exe"
gdal_polygonize_script <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Scripts/gdal_polygonize.py"
gdalwarp_path <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/gdalwarp.exe"
ogr2ogr_exe <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/ogr2ogr.exe"


# ------------------------------------------------------------
# DERIVED YEAR SETTINGS
# ------------------------------------------------------------
# CORINE land-cover maps are not available every year.
# This block maps the target fire year to the corresponding CORINE
# reference year used by the workflow.
#
# Rule:
#   target_year < 2000     -> CORINE 1990
#   2000-2005              -> CORINE 2000
#   2006-2011              -> CORINE 2006
#   2012-2017              -> CORINE 2012
#   2018 onwards           -> CORINE 2018

corine_year <- if (target_year < 2000) {
  1990
} else if (target_year < 2006) {
  2000
} else if (target_year < 2012) {
  2006
} else if (target_year < 2018) {
  2012
} else {
  2018
}


# ------------------------------------------------------------
# INPUT PATHS
# ------------------------------------------------------------
# These are the spatial inputs consumed by the deterministic workflow.
#
# Some paths depend on target_year:
#   - change_index
#   - reference_burned_map
#   - previous_year_burned
#   - hotspots
#
# Others depend on corine_year:
#   - vegetation_map
#   - burnable_mask

change_index <- file.path(
  base_dir,
  "1_DATA/Imagery/Composites_90m/Min_Min",
  paste0("MinMin_", target_year, "_mosaic_res90m.tif")
)

vegetation_map <- file.path(
  base_dir,
  "1_DATA/Corine_Masks",
  paste0("CLC_", corine_year, "_peninsula.tif")
)

burnable_mask <- file.path(
  base_dir,
  "1_DATA/Corine_Masks",
  paste0(
    "burneable_mask_binary_corine",
    substr(corine_year, 3, 4),
    "_wgs84.tif"
  )
)

reference_burned_map <- file.path(
  base_dir,
  "1_DATA/Fires/Validation_fires_burneable_verano",
  paste0("Effis_CA_", target_year, "_maskKeep_summer.shp")
)

previous_year_burned <- if (target_year > 1985) {
  file.path(
    base_dir,
    "1_DATA/Fires/Validation_fires",
    paste0("Effis_CA_", target_year - 1, ".shp")
  )
} else {
  NULL
}

hotspots <- if (target_year < 1995) {
  NULL
} else {
  file.path(
    base_dir,
    "1_DATA/Hotspots",
    paste0("hotspots_iberia_", target_year, ".geojson")
  )
}

ecoregion_shapefile_path <- file.path(
  base_dir,
  "1_DATA/Ecoregion/ecoregiones_olson.shp"
)

peninsula_shapefile_path <- file.path(
  base_dir,
  "1_DATA/Borders/Iberian_Peninsula.shp"
)


# ------------------------------------------------------------
# QUICK INPUT CHECK
# ------------------------------------------------------------
# This check verifies that the expected files exist before the workflow
# starts.
#
# Optional inputs:
#   - previous_year_burned may be NULL for 1985.
#   - hotspots is NULL for pre-hotspot years.
#
# Required inputs should exist before detection is run.

paths_to_check <- list(
  change_index = change_index,
  vegetation_map = vegetation_map,
  burnable_mask = burnable_mask,
  reference_burned_map = reference_burned_map,
  previous_year_burned = previous_year_burned,
  hotspots = hotspots,
  ecoregion_shapefile_path = ecoregion_shapefile_path,
  peninsula_shapefile_path = peninsula_shapefile_path
)

input_check <- data.frame(
  input = names(paths_to_check),
  path = vapply(
    paths_to_check,
    function(x) if (is.null(x)) NA_character_ else x,
    character(1)
  ),
  exists = vapply(
    paths_to_check,
    function(x) if (is.null(x)) NA else file.exists(x),
    logical(1)
  )
)

input_check

# Stop early if any required input is missing.
# Keep optional inputs out of this list if NULL is expected.

required_inputs <- c(
  "change_index",
  "vegetation_map",
  "burnable_mask",
  "reference_burned_map",
  "ecoregion_shapefile_path",
  "peninsula_shapefile_path"
)

if (any(!input_check$exists[input_check$input %in% required_inputs])) {
  stop("Some required input files are missing. Check input_check.")
}


# ------------------------------------------------------------
# DETECTION PRESETS
# ------------------------------------------------------------
# These are the main detection knobs.
#
# Read them like this:
#
#   seed_threshold:
#     Minimum change value needed to create a burned seed pixel.
#     Higher = fewer seeds / more conservative detection.
#
#   growth_delta:
#     Amount by which the grow threshold is relaxed relative to the seed
#     threshold.
#     Higher = easier patch expansion.
#
#   minimum_growth_threshold:
#     Absolute floor below which region growing is not allowed to go.
#     Higher = tighter, more conservative patch growth.
#
#   minimum_seed_pixels:
#     Minimum number of seed pixels needed to keep a candidate component.
#     Higher = removes small isolated detections.
#
#   *_by_vegetation:
#     Per-class overrides. Use these when one global value is too coarse
#     and vegetation types need different behaviour.

preset_params <- list(

  permissive = list(
    expected = "Highest sensitivity. More candidate burned patches, but higher commission risk.",
    detect_params = list(
      seed_threshold = 300,
      growth_delta = 100,
      minimum_growth_threshold = 200,
      minimum_seed_pixels = 20L,
      seed_threshold_by_vegetation = c(
        "1" = 300, "2" = 350, "3" = 350, "4" = 350, "5" = 350,
        "6" = 350, "7" = 350, "8" = 500, "9" = 440, "10" = 650, "11" = 450
      ),
      growth_delta_by_vegetation = c(
        "1" = 100, "2" = 50, "3" = 100, "4" = 120, "5" = 100,
        "6" = 100, "7" = 100, "8" = 10, "9" = 15, "10" = 5, "11" = 15
      ),
      minimum_growth_threshold_by_vegetation = c(
        "1" = 200, "2" = 300, "3" = 250, "4" = 210, "5" = 220,
        "6" = 220, "7" = 220, "8" = 400, "9" = 380, "10" = 550, "11" = 380
      )
    )
  ),

  balanced = list(
    expected = "Default operational compromise between omission and commission.",
    detect_params = list(
      seed_threshold = 310,
      growth_delta = 90,
      minimum_growth_threshold = 240,
      minimum_seed_pixels = 30L,
      seed_threshold_by_vegetation = c(
        "1" = 350, "2" = 420, "3" = 380, "4" = 320, "5" = 320,
        "6" = 320, "7" = 320, "8" = 500, "9" = 470, "10" = 670, "11" = 480
      ),
      growth_delta_by_vegetation = c(
        "1" = 70, "2" = 20, "3" = 40, "4" = 100, "5" = 85,
        "6" = 85, "7" = 85, "8" = 10, "9" = 10, "10" = 5, "11" = 10
      ),
      minimum_growth_threshold_by_vegetation = c(
        "1" = 280, "2" = 360, "3" = 300, "4" = 230, "5" = 240,
        "6" = 240, "7" = 240, "8" = 430, "9" = 400, "10" = 560, "11" = 400
      )
    )
  ),

  conservative = list(
    expected = "More restrictive. Reduces false positives and over-expansion, with increased omission risk.",
    detect_params = list(
      seed_threshold = 400,
      growth_delta = 50,
      minimum_growth_threshold = 320,
      minimum_seed_pixels = 50L,
      seed_threshold_by_vegetation = c(
        "1" = 400, "2" = 470, "3" = 430, "4" = 400, "5" = 400,
        "6" = 400, "7" = 400, "8" = 520, "9" = 500, "10" = 700, "11" = 520
      ),
      growth_delta_by_vegetation = c(
        "1" = 50, "2" = 10, "3" = 20, "4" = 60, "5" = 50,
        "6" = 50, "7" = 50, "8" = 5, "9" = 5, "10" = 5, "11" = 5
      ),
      minimum_growth_threshold_by_vegetation = c(
        "1" = 330, "2" = 400, "3" = 370, "4" = 330, "5" = 330,
        "6" = 330, "7" = 330, "8" = 450, "9" = 430, "10" = 580, "11" = 430
      )
    )
  ),

  very_conservative = list(
    expected = paste(
      "Intermediate anti-banding preset for problematic historical years.",
      "Designed to suppress Landsat banding artefacts while preserving",
      "medium and large real burned patches."
    ),
    detect_params = list(
      seed_threshold = 430,
      growth_delta = 40,
      minimum_growth_threshold = 350,
      minimum_seed_pixels = 70L,
      seed_threshold_by_vegetation = c(
        "1" = 430, "2" = 480, "3" = 450, "4" = 430, "5" = 430,
        "6" = 430, "7" = 430, "8" = 540, "9" = 520, "10" = 720, "11" = 540
      ),
      growth_delta_by_vegetation = c(
        "1" = 40, "2" = 10, "3" = 20, "4" = 45, "5" = 40,
        "6" = 40, "7" = 40, "8" = 5, "9" = 5, "10" = 5, "11" = 5
      ),
      minimum_growth_threshold_by_vegetation = c(
        "1" = 360, "2" = 400, "3" = 370, "4" = 360, "5" = 360,
        "6" = 360, "7" = 360, "8" = 470, "9" = 450, "10" = 600, "11" = 450
      )
    )
  )
)

# Make sure the selected preset exists before building the config.
if (!selected_preset %in% names(preset_params)) {
  stop(
    "selected_preset must be one of: ",
    paste(names(preset_params), collapse = ", ")
  )
}

# Optional inspection:
preset_params[[selected_preset]]$expected
preset_params[[selected_preset]]$detect_params


# ------------------------------------------------------------
# REFINEMENT PARAMETERS
# ------------------------------------------------------------
# These are the current package defaults written explicitly so the user
# can see them and modify them if needed.

refine_params <- list(
  aoi_buffer_m = 5000,
  omission_buffer_m = 0,
  minimum_detected_area_m2 = 10,
  merge_overlaps = TRUE,
  merge_overlaps_buffer_m = 0
)

# Parameter meaning:
#
#   aoi_buffer_m:
#     Buffer used around AOI boundaries during refinement to reduce edge
#     artefacts.
#
#   omission_buffer_m:
#     Reserved cleanup buffer for omission-side adjustments.
#     0 means no extra omission buffer.
#
#   minimum_detected_area_m2:
#     Minimum polygon area kept after refinement.
#     Very permissive by default.
#
#   merge_overlaps:
#     If TRUE, overlapping refined polygons are merged.
#
#   merge_overlaps_buffer_m:
#     Optional buffer before overlap merging.
#     0 means merge only direct overlaps.


# ------------------------------------------------------------
# SCORING PARAMETERS
# ------------------------------------------------------------
# These parameters belong to the deterministic scoring stage itself.

scoring_params <- list(
  support_buffer_m = 0,
  previous_fire_exclusion_buffer_m = 90,
  previous_fire_cleanup_buffer_m = 90,
  minimum_remaining_area_m2 = 20000,
  reference_buffer_m = 90
)

# Parameter meaning:
#
#   support_buffer_m:
#     Buffer used when checking local contextual support around candidate
#     patches.
#
#   previous_fire_exclusion_buffer_m:
#     Buffer around previous-year burns used to detect temporal conflicts.
#
#   previous_fire_cleanup_buffer_m:
#     Cleanup buffer used after overlap removal with previous-year burns.
#
#   minimum_remaining_area_m2:
#     Minimum area a patch must retain after temporal cleanup.
#
#   reference_buffer_m:
#     Buffer used when comparing against support-reference areas in scoring.


# ------------------------------------------------------------
# OUTPUT PATHS
# ------------------------------------------------------------
# output_dir must be the root results folder only.
#
# The package itself will append:
#   <target_year>/DETERMINISTIC/<run_name>
#
# Therefore, do not manually add target_year or DETERMINISTIC here,
# otherwise they will be duplicated in the final path.

output_dir <- file.path(
  base_dir,
  "1_DATA/Results"
)

# This is the visible name of the run inside output_dir.
#
# Examples:
#   2025 + balanced      -> 2025_balanced
#   2025 + conservative  -> 2025_conservative
#   1985 + balanced      -> 1985_balanced

run_name <- paste0(target_year, "_", selected_preset)


# ------------------------------------------------------------
# STEP 1. Build the deterministic configuration
# ------------------------------------------------------------
# At this stage the workflow is only being prepared.
#
# Nothing is detected yet.
# Nothing is scored yet.
#
# This is the right moment to inspect:
#   - selected_preset
#   - detect_params
#   - refine_params
#   - scoring_params
#   - input paths
#
# Check the help if needed:
# ?build_burned_mapping_config

cfg <- build_burned_mapping_config(
  change_index = change_index,
  vegetation_map = vegetation_map,
  burnable_mask = burnable_mask,
  hotspots = hotspots,
  previous_year_burned = previous_year_burned,
  reference_burned_map = reference_burned_map,
  target_year = target_year,
  output_dir = output_dir,
  run_name = run_name,
  detect_params = preset_params[[selected_preset]]$detect_params,
  refine_params = refine_params,
  scoring_params = scoring_params,
  options = list(
    ecoregion_shapefile_path = ecoregion_shapefile_path,
    peninsula_shapefile_path = peninsula_shapefile_path,
    change_index_validation = list(
      expected_nodata = -9999,
      lower_cap = -1000,
      cap_below = TRUE,
      check_finite = TRUE
    )
  )
)

# Attach external tool paths to the validated config.

cfg$tool_paths$python_exe <- python_exe
cfg$tool_paths$gdal_polygonize_script <- gdal_polygonize_script
cfg$tool_paths$gdalwarp_path <- gdalwarp_path
cfg$tool_paths$ogr2ogr_exe <- ogr2ogr_exe

# Optional quick checks after STEP 1:
# cfg
# str(cfg, max.level = 2)
# cfg$detect_params
# cfg$refine_params
# cfg$scoring_params
# cfg$output_routes


# ------------------------------------------------------------
# STEP 2. Detect candidate burned patches
# ------------------------------------------------------------
# Detection uses the chosen preset.
#
# The preset controls:
#   - how easily burned seed pixels are created;
#   - how far candidate patches can grow;
#   - how much behaviour differs by vegetation / land-cover class.

det <- detect_burned_patches(
  config = cfg,
  write_outputs = TRUE,
  overwrite = TRUE
)

# Save the two most important detection output paths explicitly.

grown_candidates_path <- det$grown_patches_path
refined_candidates_path <- det$refined_patches_path

# Fail early if the expected detection outputs were not written.

stopifnot(file.exists(grown_candidates_path))
stopifnot(file.exists(refined_candidates_path))

# Optional quick checks after STEP 2:
grown_candidates_path
refined_candidates_path
file.exists(grown_candidates_path)
file.exists(refined_candidates_path)


# ------------------------------------------------------------
# STEP 3. Same-year scoring
# ------------------------------------------------------------
# This is the native deterministic scoring behaviour.
#
# Use this route when you want a standard run for one year and one
# preset, for example:
#   2025 balanced
#   2025 conservative
#
# Important:
#   keep_pool = NULL
#
# means:
#   use the package's built-in same-year local fallback.
#
# Check the help if needed:
# ?score_burned_patches

scored <- score_burned_patches(
  burned_candidates = refined_candidates_path,
  config = cfg,
  keep_pool = NULL,
  write_outputs = TRUE,
  overwrite = TRUE
)

# Optional quick checks after STEP 3:
scored
scored$internal_decisions_path
scored$keep_pool_summary
names(scored)

table(scored$internal_decisions$class_final, useNA = "ifany")
table(scored$internal_decisions$filter_3, useNA = "ifany")


# ------------------------------------------------------------
# VALIDATION PATHS
# ------------------------------------------------------------

validation_dir <- file.path(cfg$output_routes$base, "VALIDATION_KEEP_ONLY")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Keep only final polygons classified as "keep"
# ------------------------------------------------------------

keep_path <- file.path(validation_dir, "keep_only.gpkg")

keep_only <- scored$internal_decisions[
  scored$internal_decisions$class_final == "keep",
]

keep_only <- sf::st_zm(keep_only, drop = TRUE, what = "ZM")
keep_only <- sf::st_make_valid(keep_only)
keep_only <- keep_only[!sf::st_is_empty(keep_only), , drop = FALSE]

if (nrow(keep_only) == 0) {
  stop("No polygons with class_final == 'keep' were found. Validation cannot continue.")
}

if (file.exists(keep_path)) unlink(keep_path, force = TRUE)

sf::st_write(
  keep_only,
  keep_path,
  layer = "keep_only",
  quiet = TRUE
)

# ------------------------------------------------------------
# 2. Prepare validation mask
# ------------------------------------------------------------

validation_mask_path <- file.path(
  validation_dir,
  "support",
  "validation_mask_geometry_only.gpkg"
)

dir.create(dirname(validation_mask_path), recursive = TRUE, showWarnings = FALSE)

validation_mask_shapefile <- if (target_year < 2000) {
  file.path(
    base_dir,
    "1_DATA/Mask_StudyArea",
    paste0("mask_3035_", target_year, ".shp")
  )
} else {
  file.path(
    base_dir,
    "1_DATA/Mask_StudyArea",
    "mask_Peninsula_3035.shp"
  )
}

mask_sf <- sf::st_read(validation_mask_shapefile, quiet = TRUE) |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_make_valid()

mask_sf <- mask_sf[!sf::st_is_empty(mask_sf), 0, drop = FALSE]
mask_sf$mask_id <- seq_len(nrow(mask_sf))

if (nrow(mask_sf) == 0) {
  stop("The validation mask is empty after geometry cleaning.")
}

if (file.exists(validation_mask_path)) unlink(validation_mask_path, force = TRUE)

sf::st_write(
  mask_sf,
  validation_mask_path,
  layer = "mask",
  quiet = TRUE
)

# ------------------------------------------------------------
# 3. External reference burned map
# ------------------------------------------------------------

reference_burned_map <- file.path(
  base_dir,
  "1_DATA/Fires/Validation_fires_burneable_verano",
  paste0("Effis_CA_", target_year, "_maskKeep_summer.shp")
)

# ------------------------------------------------------------
# 4. Optional stratified validation by CORINE-derived strata
# ------------------------------------------------------------

strata_raster <- file.path(
  base_dir,
  "1_DATA/Corine_Masks/STRATA",
  paste0("strata_CLC_", corine_year, "_res30.tif")
)

strata_lut <- file.path(
  base_dir,
  "1_DATA/Corine_Masks/LUT/lut_full_strata8_v1.csv"
)

# ------------------------------------------------------------
# 5. Run validation
# ------------------------------------------------------------

val <- validate_fire_maps(
  input_shapefile = keep_path,
  ref_shapefile = reference_burned_map,
  mask_shapefile = validation_mask_path,
  burnable_raster = burnable_mask,
  year_target = target_year,
  validation_dir = validation_dir,
  force_reprocess_ref = FALSE,
  force_reprocess_pred = FALSE,
  metrics_type = "all",
  dissolve_ref_by = "id",
  dissolve_input_by = NULL,
  strata_raster = strata_raster,
  strata_lut = strata_lut,
  observability_raster = change_index,
  ref_end_doy_col = "end_doy",
  ref_start_doy_col = "start_doy"
)

# Optional quick checks after validation:
# val
# val$metrics
# val$polygon_summary
# val$reference_observability


# -------------------------------------------------------------------
# End of the year-with-active-hotspot-support draft.
# Later revisions will execute these steps interactively and discuss
# the intermediate outputs.
# -------------------------------------------------------------------
