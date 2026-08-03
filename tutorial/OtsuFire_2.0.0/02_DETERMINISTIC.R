# =============================================================================
# OtsuFire 2.0.0 tutorial - STAGE 2: DETERMINISTIC
# =============================================================================
#
# GOAL
#   Find burned patches with rules only - no training data, no model - and
#   label each one keep / review / drop.
#
# HOW IT WORKS, IN ONE PARAGRAPH
#   The change index is split into burned / unburned by an Otsu threshold. The
#   threshold is NOT global: it is computed per CORINE x ecoregion unit, so a
#   pine forest in the north and a shrubland in the south are cut at different
#   severities. Pixels above a strict "seed" threshold start a patch; the patch
#   then grows into neighbouring pixels above a looser threshold. Grown patches
#   pass through filters (size, previous-year overlap, spectral support) and
#   each ends up as keep, review or drop.
#
# THE OUTPUT THAT MATTERS
#   internal_decisions.gpkg - one polygon per candidate, with its class. This
#   file is the training input for stage 3. Everything else is diagnostics.
#
# RUN 00_SETUP.R FIRST.
# =============================================================================

source("00_SETUP.R")


# -----------------------------------------------------------------------------
# 1) The three parameter blocks
# -----------------------------------------------------------------------------
# In 2.0.0 there are no closed scenario presets any more. You pass the
# ecological controls explicitly, in three blocks.

# --- detect_params: where a patch starts and how far it grows ----------------
#
#   seed_threshold             severity a pixel needs to START a patch. Higher
#                              = fewer, more certain patches.
#   growth_delta               how far below the seed a neighbour may be and
#                              still be absorbed. Higher = patches grow more.
#   minimum_growth_threshold   hard floor: never grow below this severity,
#                              whatever growth_delta says.
#   minimum_seed_pixels        a patch needs this many seed pixels to exist.
#                              This is the main noise filter.
#
# The *_by_vegetation vectors override the global value per CORINE class. The
# names are CORINE class codes as character ("1".."11"). Use them when the same
# severity means different things in different vegetation: a low-severity fire
# in shrubland can be a real fire, while the same value in a dense forest is
# usually not.

detect_params <- list(
  seed_threshold           = 300,
  growth_delta             = 100,
  minimum_growth_threshold = 200,
  minimum_seed_pixels      = 20L,

  seed_threshold_by_vegetation = c(
    "1" = 300, "2" = 350, "3" = 350, "4" = 350, "5" = 350, "6" = 350,
    "7" = 350, "8" = 500, "9" = 440, "10" = 650, "11" = 450
  ),
  growth_delta_by_vegetation = c(
    "1" = 100, "2" = 50, "3" = 100, "4" = 120, "5" = 100, "6" = 100,
    "7" = 100, "8" = 10, "9" = 15, "10" = 5, "11" = 15
  ),
  minimum_growth_threshold_by_vegetation = c(
    "1" = 200, "2" = 300, "3" = 250, "4" = 210, "5" = 220, "6" = 220,
    "7" = 220, "8" = 400, "9" = 380, "10" = 550, "11" = 380
  )
)

# --- refine_params: geometry clean-up ----------------------------------------
#
#   aoi_buffer_m               working buffer around the area of interest
#   minimum_detected_area_m2   drop slivers below this area
#   merge_overlaps             merge patches that touch or overlap
#
refine_params <- list(
  aoi_buffer_m             = 5000,
  omission_buffer_m        = 0,
  minimum_detected_area_m2 = 10,
  merge_overlaps           = TRUE,
  merge_overlaps_buffer_m  = 0
)

# --- scoring_params: the filters that decide keep / review / drop ------------
#
#   previous_fire_exclusion_buffer_m   how close to last year's scar a patch
#                                      may sit before being penalised
#   minimum_remaining_area_m2          area a patch must still have AFTER the
#                                      exclusions are applied (20000 m2 = 2 ha)
#   reference_buffer_m                 tolerance when comparing to reference
#
scoring_params <- list(
  support_buffer_m                 = 0,
  previous_fire_exclusion_buffer_m = 90,
  previous_fire_cleanup_buffer_m   = 90,
  minimum_remaining_area_m2        = 20000,
  reference_buffer_m               = 90
)


# -----------------------------------------------------------------------------
# 2) Build the configuration
# -----------------------------------------------------------------------------
# The config is just a description of the run. It executes nothing. Building it
# validates your paths and types straight away, so a typo fails in a second
# instead of two hours in.

det_cfg <- build_burned_mapping_config(
  change_index         = CHANGE_INDEX,
  vegetation_map       = VEGETATION_MAP,
  burnable_mask        = BURNABLE_MASK,
  hotspots             = HOTSPOTS,               # NULL before 2000
  previous_year_burned = PREVIOUS_YEAR_BURNED,
  reference_burned_map = REFERENCE_BURNED_MAP,
  target_year          = TARGET_YEAR,
  output_dir           = OUTPUT_DIR,
  run_name             = RUN_NAME,

  detect_params  = detect_params,
  refine_params  = refine_params,
  scoring_params = scoring_params,

  options = list(
    ecoregion_shapefile_path = ECOREGION_SHP,
    peninsula_shapefile_path = PENINSULA_SHP,
    # Re-checks the mosaic assumptions from stage 1 at run time.
    change_index_validation = list(
      expected_nodata = -9999,
      lower_cap       = -1000,
      cap_below       = TRUE,
      check_finite    = TRUE
    )
  )
)


# -----------------------------------------------------------------------------
# 3) Run it
# -----------------------------------------------------------------------------
# One call does detection + scoring. Expect tens of minutes for a full
# peninsula year.

det <- run_deterministic_pipeline(
  det_cfg,
  write_outputs  = TRUE,
  overwrite      = FALSE,   # TRUE re-does work that already exists on disk
  run_validation = FALSE    # stage 4 handles validation properly
)

# The paths it wrote:
det$result_paths$grown_patches
det$result_paths$refined_patches
det$result_paths$internal_decisions   # <- the input to stage 3
det$result_paths$timing_csv


# -----------------------------------------------------------------------------
# 4) Look at what came out
# -----------------------------------------------------------------------------
dec <- sf::st_read(det$result_paths$internal_decisions, quiet = TRUE)

nrow(dec)
table(dec$class_final, useNA = "ifany")   # keep / review / drop counts

# A useful sanity check: total kept area, in hectares.
sum(as.numeric(sf::st_area(dec[dec$class_final == "keep", ])) / 10000)

# If "keep" is a handful of polygons, the thresholds are too strict for this
# year. If it is tens of thousands, they are too loose. Adjust seed_threshold
# and minimum_seed_pixels first - they move the result most.


# -----------------------------------------------------------------------------
# 5) Step by step, when you want to see the middle
# -----------------------------------------------------------------------------
# run_deterministic_pipeline() is a wrapper over two functions. Call them
# directly when you want to inspect the candidates before they are scored.

if (FALSE) {
  patches <- detect_burned_patches(det_cfg, write_outputs = TRUE)
  patches$grown_patches_path

  scored <- score_burned_patches(patches, det_cfg, write_outputs = TRUE)
  scored$internal_decisions_path
}


# -----------------------------------------------------------------------------
# 6) Optional: de-chaining
# -----------------------------------------------------------------------------
# Sometimes two real fires are joined by a thin neck of pixels and come out as
# one huge patch. apply_chain_cleaning() cuts necks narrower than a given
# width and splits such a patch back into its parts.
#
# It is opt-in and does nothing to well-shaped patches: only candidates whose
# inflation (1/sqrt(compactness)) exceeds inflation_thr are touched at all.

if (FALSE) {
  chained <- apply_chain_cleaning(
    refine_path        = det$result_paths$refined_patches,
    out_dir            = NULL,          # inferred from refine_path
    year               = TARGET_YEAR,
    dechain_distance_m = 180,           # neck width to cut: try 90 / 180 / 270
    pixel_size_m       = 90,
    inflation_thr      = 22.135,        # below this, a patch is left alone
    min_area_ha        = NULL,
    marker             = "connected",
    verbose            = TRUE
  )

  # Then score the de-chained patches instead of the original ones.
}


# -----------------------------------------------------------------------------
# WHAT YOU HAVE NOW
# -----------------------------------------------------------------------------
#   internal_decisions.gpkg: candidate polygons labelled keep / review / drop.
#
#   These labels are RULES, not ground truth. Stage 3 learns from them, and
#   stage 4 checks the result against an independent reference.
#
# NEXT: 03_SUPERVISED.R
# =============================================================================
