# =============================================================================
# BLOCK 0 - COMPLETE SETUP: THE WHOLE USER CONFIGURATION
# =============================================================================
#
# WHAT IT DOES
#   Loads the OtsuFire package in development mode and declares EXPLICITLY
#   every experiment variable the user has to configure:
#     - 2 combo parameters (year + scenario)
#     - 13 input paths (what the pipeline reads from disk)
#     - 18 methodological constants (policies + thresholds + toggles)
#
#   This block is the ONLY one the user has to understand and edit in
#   order to run the pipeline. Everything that follows (BLOCK 1+) simply
#   consumes these variables.
#
#
# TWO READING LEVELS
# ------------------
# Every section of the block carries narrative comments (what is done and
# why) plus DEBUG blocks (cross-references to the handoff with technical
# detail for diagnosing package problems).
#
#
# STATE OF THE API IN 0.5.0 (summary, for context)
# ------------------------------------------------
# OtsuFire 0.5.0 has an incomplete configuration API: of the 33 variables
# that conceptually make up the experiment, only some are accepted by the
# public config builder. The rest are built by the orchestrator by
# convention (paths) or read from hardcoded global constants (parameters).
# This is to be resolved in the post-paper refactor.
#
# Cross-references:
#   §N+25: 13 input paths - only 4 arrive via cfg$inputs in 0.5.0.
#   §N+26: negative_pool_policy - REMOVED from the package on
#          2026-06-05; all_sources is now the only (implicit) mode.
#   §N+27: burned_like_registry_path field - abandoned line of
#          investigation, REMOVED from the package on 2026-06-05.
#   §N+28: 21 constants hardcoded in the orchestrator - decide which
#          ones to expose after the refactor.
#
# So that this tutorial stays forward-compatible and doubles as a spec for
# the refactor, all 33 variables are declared explicitly HERE. Later, in
# BLOCK 1, they are split into:
#   - Layer B (commented out): the ideal post-refactor API, where the user
#                          passes all 33 variables to the config builder.
#   - Layer C (executable): the current 0.5.0 API, where only some can be
#                          passed; the rest are used in later blocks
#                          or simply documented.
#
# =============================================================================


# =============================================================================
# 1) LOADING THE PACKAGE
# =============================================================================

cat("\n========== BLOCK 0: SETUP ==========\n")

PKG_ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"

# Clean load: if it was already loaded, unload it first to make sure we
# pick up the current version of the code on disk (this matters during
# active development, where the package code keeps changing).
if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

cat("OtsuFire version:", as.character(utils::packageVersion("OtsuFire")), "\n")
stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

# DEBUG: pkgload::load_all exposes both the NAMESPACE exports and the
# internal functions. That is what allows OtsuFire:::xxx to be used to
# reach internals, which later blocks of the tutorial need in order to
# "open up the package" and call functions step by step instead of going
# through the run_oneyear_supervised_pipeline() wrapper.


# =============================================================================
# 2) COMBO PARAMETERS (year + scenario)
# =============================================================================
#
# To run the tutorial on a different combo, change ONLY these two lines.
# All the rest of BLOCK 0 (paths, constants) adapts automatically.

YEAR     <- 2005L
SCENARIO <- "balanced"

# DEBUG: YEAR must be an integer (L suffix). Passing a double (2005 with
# no L) can make the package functions that expect an integer target_year
# fail with coercion warnings. SCENARIO must be one of: "balanced", "lax",
# "original" (the three supported by the deterministic stage).

# Resolution of the five-yearly Corine year via an internal function.
# Mapping (see R/utils-corine.R::get_corine_year):
#   1985-2002 -> 2000 ; 2003-2008 -> 2006 ; 2009-2014 -> 2012 ;
#   2015+    -> 2018
# For YEAR=2005, CORINE_YEAR resolves to "2006" (returned as character).
CORINE_YEAR <- OtsuFire:::get_corine_year(YEAR)

cat(sprintf("\nCombo: year=%d, scenario=%s, corine=%s\n",
            YEAR, SCENARIO, CORINE_YEAR))

# DEBUG: get_corine_year returns character, not integer. That is why the
# sprintf calls for Corine paths use %s rather than %d. Should it ever be
# changed to integer, update those sprintf calls.


# =============================================================================
# 3) PROJECT PATHS
# =============================================================================

DATA_BASE <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA"
COMPOSITE <- file.path(DATA_BASE, "Imagery", "Composites_90m")
RESULTS   <- file.path(DATA_BASE, "Results")

# RESULT_NAME identifies the pipeline "run". Min_Min refers to the
# temporal compositing method behind the RBR mosaics (Min of Min).
# The package supports others (Mean_Mean, Min_Mean, etc.) but the
# canonical Baseline uses Min_Min.
RESULT_NAME <- "Min_Min"


# =============================================================================
# 4) THE 13 PIPELINE INPUTS
# =============================================================================
#
# These are the 13 files the supervised module orchestrator reads in order
# to process one combo. They are organised into 4 groups by the kind of
# dependency:
#
#   A) Year/scenario-dependent (4)  - change with (year, scenario)
#   B) External validation, year-dep (1) - changes with year only
#   C) Five-yearly Corine          (4) - change with CORINE_YEAR
#   D) Static                      (4) - never change
#
# DEBUG: HANDOFF §N+25. In 0.5.0 only 4 of these 13 reach the config
# builder via cfg$inputs (internal_decisions, change_index,
# delayed_change_index, hotspots). The other 9 are built by the
# orchestrator by convention at lines 906-934 of internal-sup-
# orchestrator.R, plus the burnable mask in internal helpers.
# The post-paper refactor will route ALL of them through cfg$inputs.


# --- A) Year/scenario-dependent (4) -----------------------------------------

# A1. Output of the deterministic module (main input of the supervised one).
INTERNAL_DECISIONS <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "05_DECISIONS", "internal_decisions.gpkg"
)

# A2. Summer mosaic (RBR + DOY of the fire summer).
CHANGE_INDEX <- file.path(COMPOSITE, RESULT_NAME,
                          sprintf("MinMin_%d_mosaic_res90m.tif", YEAR))

# A3. Autumn-winter mosaic (delayed change index, source G2_RBR_AW).
RBR_AUTUMN <- file.path(COMPOSITE, "Autumn",
                        sprintf("mean_mean_%d_mosaic.tif", YEAR))

# A4. MODIS hotspots. For pre-MODIS years (<2000) HOTSPOTS will not exist
# on disk; the pipeline detects this and proceeds without hotspot features.
HOTSPOTS <- file.path(DATA_BASE, "Hotspots",
                      sprintf("hotspots_iberia_%d.geojson", YEAR))


# --- B) Year-dependent external validation (1) ------------------------------

# B1. Effis-CA summer fires filtered by the burnable mask.
# This is the EXTERNAL ground truth used to validate the model predictions.
# DEBUG: HANDOFF §N+25. In 0.5.0 it corresponds to the field
# cfg$inputs$reference_burned_map (cosmetic) plus ref_tif_path/ref_shp_path
# (built by convention). The refactor will consolidate both under the
# public reference_burned_map field.
EFFIS_CA_TIF <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.tif", YEAR)
)

# The companion shp is derivable from the tif (swap .tif for .shp).
# It is declared explicitly because the 0.5.0 orchestrator reads them as
# separate paths. After the refactor it will be derivable internally from
# reference_burned_map.
EFFIS_CA_SHP <- sub("\\.tif$", ".shp", EFFIS_CA_TIF)


# --- C) Five-yearly Corine (4) ----------------------------------------------

# C1. Corine Land Cover raster (land cover).
# Source of the G3_CORINE features.
CORINE_RASTER <- file.path(DATA_BASE, "Corine_Masks",
                           sprintf("CLC_%s_peninsula.tif", CORINE_YEAR))

# C2. Strata raster (stratification derived from Corine for the stratified
# spatial sampling in STEP B1 - block folds).
CORINE_STRATA <- file.path(DATA_BASE, "Corine_Masks", "STRATA",
                           sprintf("strata_CLC_%s_res30.tif", CORINE_YEAR))

# C3. Look-up table for the strata (does not depend on the year).
CORINE_LUT <- file.path(DATA_BASE, "Corine_Masks", "LUT",
                        "lut_full_strata8_v1.csv")

# C4. Binary burnable-surface mask derived from Corine.
# DEBUG: HANDOFF §N+25 §1.3. This is the only input the orchestrator does
# NOT validate in its initial stopifnot block. It is read inside
# internal-sup-unburned-deterministic.R at line 170. If it is missing, the
# pipeline fails 5-10 min after starting rather than up front. That is why
# it is verified by hand here in BLOCK 0.
BURNEABLE_MASK <- file.path(
  DATA_BASE, "Corine_Masks",
  sprintf("burneable_mask_binary_corine_%s_ETRS89.tif", CORINE_YEAR)
)


# --- D) Static (4) ----------------------------------------------------------

# D1. Border of the Iberian Peninsula.
PENINSULA_SHP <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")

# D2. Topography: stack with DEM (band 1) and slope (band 2).
# Source of the G4_TOPO features.
TOPO <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

# D3, D4. Study-area mask in EPSG:3035 (raster + shp).
MASK_TIF <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")
MASK_SHP <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.shp")


# =============================================================================
# 5) METHODOLOGICAL CONSTANTS (18 variables)
# =============================================================================
#
# These are the methodological decisions of the canonical Baseline, grouped
# by function. Most users do NOT need to change them: the values here are
# the ones from the published Baseline.
#
# DEBUG: HANDOFF §N+28. In 0.5.0:
#   - 8 of these are accepted via cfg$options, plus 1 at cfg top level.
#   - 9 are hardcoded in internal-sup-orchestrator.R and are NOT
#     configurable without editing the package code.
# The tutorial declares them ALL here so that Layer B of BLOCK 1 can show
# the ideal post-refactor API.


# --- 5.1) Negative sampling policy (1) --------------------------------------

# Sources of negatives that enter the training pool.
# NOTE (2026-06-05): the negatives policy is no longer configurable; the
# package ALWAYS uses all_sources (the only mode, implicit). The old
# "deterministic_direct" was removed. This local constant is kept purely
# as documentation; it is not passed to the config builder.
#   all_sources - 4 sources: internal_keep_qc (burned) + deterministic
#                 drops + random background + Otsu unburned patches.
NEGATIVE_POOL_POLICY <- "all_sources"


# --- 5.2) Minimum-size guard for the burned pool (1) ------------------------

# If the deterministic stage yields fewer than N burned polygons for the
# combo, the supervised stage aborts cleanly with an informative error (a
# model cannot be trained on so few positives).
# This is the guard that aborted the three 1997 combos (n=0/1, see §N+24).
MIN_BURNED_POOL_N <- 5L


# --- 5.3) Reproducibility: random seeds (2) ---------------------------------

# Seed of the random_burnable_background sampling (all_sources Part 1).
RANDOM_SEED <- 42L

# Seed of the Otsu unburned sampling (all_sources branch, Part 2).
# DEBUG: HANDOFF §N+28. The package has TWO independent seeds, both
# hardcoded to 42. The refactor could unify them into one if the
# independence turns out not to be worth keeping.
LEGACY_RANDOM_SEED <- 42L


# --- 5.4) Between-year temporal filter (3) ----------------------------------
#
# These penalise polygons of the current year that overlap too much with
# polygons of the previous year (proxy: "persistent fires are probably not
# real fires"). DEBUG: HANDOFF §N+28 - currently hardcoded.

# Spatial-overlap threshold for declaring a temporal conflict.
CURRENTYEAR_PREYEAR_OVERLAP_THR <- 0.70

# Minimum hotspot density (per ha) for a polygon to pass the filter.
CURRENTYEAR_HOTSPOT_DENSITY_THR <- 0.001

# Floor of the penalty score (it never drops below this value).
CURRENTYEAR_TEMPORAL_PENALTY_FLOOR <- 0.10


# --- 5.5) Otsu pipeline - main parameters (4) -------------------------------
#
# These control the Otsu branch of the all_sources mode (legacy_otsu_*).

# How Otsu is applied. "burnable_only" = only over the pixels the Corine
# mask considers burnable (this keeps Otsu from seeing water, urban, etc.
# and distorting the thresholds).
LEGACY_OTSU_MODE <- "burnable_only"

# Minimum RBR threshold for an Otsu patch to be considered at all.
# A guard against negative spectral noise.
LEGACY_OTSU_THRESHOLD <- 0

# Reference threshold used to scale the other relative Otsu thresholds.
LEGACY_REFERENCE_OTSU_THRESHOLD <- 100

# Cap on the otsu_patch_residual pool. Without it, the ~13,000 candidate
# patches would swamp the training set.
LEGACY_SAMPLE_N <- 2000L


# --- 5.6) Otsu pipeline - confidence thresholds (2) -------------------------
#
# DEBUG: HANDOFF §N+28 Tier 1. Hardcoded in 0.5.0, exposed after the refactor.

# Threshold above which an Otsu patch counts as "high confidence".
OTSU_KEEP_HI <- 0.45

# Threshold below which an Otsu patch counts as "low confidence" (drop).
OTSU_DROP_LO <- 0.15


# --- 5.7) Pipeline toggles (4) ----------------------------------------------
#
# These allow pipeline stages to be skipped (handy for partial reruns).
# DEBUG: HANDOFF §N+28 Tier 2. Hardcoded to TRUE in 0.5.0.

DO_POOLS    <- TRUE  # STEP A: pools (burned + unburned)
DO_FOLDS    <- TRUE  # STEP B1-B2: spatial block folds
DO_FEATURES <- TRUE  # STEP B3: extract features per polygon
DO_MODEL    <- TRUE  # STEP C: OOF + final model + scoring + map


# --- 5.8) I/O control (3) ---------------------------------------------------

# If TRUE and Otsu outputs from a previous run exist, they are reused
# instead of regenerated. Only affects the Otsu branch of all_sources.
LEGACY_REUSE_EXISTING <- TRUE

# If TRUE, intermediate Otsu outputs are written to disk.
LEGACY_WRITE_OUTPUT <- TRUE

# If TRUE, enables the detailed STEP A4 log (the "[unb] ..." messages).
UNB_VERBOSE <- TRUE


# =============================================================================
# 6) TUTORIAL OUTPUT DIRECTORY
# =============================================================================
#
# So as not to trample the real run, the tutorial outputs are redirected to
# a parallel folder, SUPERVISED_TUTORIAL/ (not SUPERVISED/).

TUTORIAL_RESULT_DIR <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "SUPERVISED_TUTORIAL", SCENARIO
)
dir.create(TUTORIAL_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 7) GLOBAL VERIFICATION
# =============================================================================
#
# The 13 inputs are validated BEFORE the pipeline starts. If something is
# missing we fail in under a second instead of 5-10 minutes later.

stopifnot(
  # A) Year-dependent
  file.exists(INTERNAL_DECISIONS),
  file.exists(CHANGE_INDEX),
  file.exists(RBR_AUTUMN),
  file.exists(HOTSPOTS),
  # B) External validation
  file.exists(EFFIS_CA_TIF),
  file.exists(EFFIS_CA_SHP),
  # C) Corine
  file.exists(CORINE_RASTER),
  file.exists(CORINE_STRATA),
  file.exists(CORINE_LUT),
  file.exists(BURNEABLE_MASK),
  # D) Static
  file.exists(PENINSULA_SHP),
  file.exists(TOPO),
  file.exists(MASK_TIF),
  file.exists(MASK_SHP)
)


# =============================================================================
# 8) VISUAL SUMMARY OF THE CONFIGURATION
# =============================================================================

cat("\n========== EXPERIMENT CONFIGURATION ==========\n")
cat("\n[1] COMBO PARAMETERS\n")
cat("  YEAR          :", YEAR, "\n")
cat("  SCENARIO      :", SCENARIO, "\n")
cat("  CORINE_YEAR   :", CORINE_YEAR, "\n")
cat("  RESULT_NAME   :", RESULT_NAME, "\n")

cat("\n[2] THE 13 INPUTS (verified)\n")

cat("\n  A) Year/scenario-dependent (4):\n")
cat("    INTERNAL_DECISIONS :", basename(INTERNAL_DECISIONS), "\n")
cat("    CHANGE_INDEX       :", basename(CHANGE_INDEX), "\n")
cat("    RBR_AUTUMN         :", basename(RBR_AUTUMN), "\n")
cat("    HOTSPOTS           :", basename(HOTSPOTS), "\n")

cat("\n  B) External validation (1, +shp):\n")
cat("    EFFIS_CA_TIF       :", basename(EFFIS_CA_TIF), "\n")
cat("    EFFIS_CA_SHP       :", basename(EFFIS_CA_SHP), "\n")

cat("\n  C) Five-yearly Corine (4):\n")
cat("    CORINE_RASTER      :", basename(CORINE_RASTER), "\n")
cat("    CORINE_STRATA      :", basename(CORINE_STRATA), "\n")
cat("    CORINE_LUT         :", basename(CORINE_LUT), "\n")
cat("    BURNEABLE_MASK     :", basename(BURNEABLE_MASK), "\n")

cat("\n  D) Static (4):\n")
cat("    PENINSULA_SHP      :", basename(PENINSULA_SHP), "\n")
cat("    TOPO               :", basename(TOPO), "\n")
cat("    MASK_TIF           :", basename(MASK_TIF), "\n")
cat("    MASK_SHP           :", basename(MASK_SHP), "\n")

cat("\n[3] METHODOLOGICAL CONSTANTS\n")

cat("\n  3.1 Negative sampling policy\n")
cat("    NEGATIVE_POOL_POLICY              :", NEGATIVE_POOL_POLICY, "\n")
cat("    : always all_sources (the only mode, implicit since 2026-06-05).\n")
cat("      Combines 4 sources: internal_keep_qc (burned) + deterministic\n")
cat("      drops + random_background + Otsu_patches.\n")

cat("\n  3.2 Minimum-size guard for the burned pool\n")
cat("    MIN_BURNED_POOL_N                 :", MIN_BURNED_POOL_N, "\n")
cat("    : combos with fewer than N burned polygons abort cleanly.\n")
cat("      Lowering it to 3 admits years with very few fires (risk of an\n")
cat("      unreliable model). Raising it to 10 discards more pre-MODIS years.\n")

cat("\n  3.3 Reproducibility: random seeds\n")
cat("    RANDOM_SEED / LEGACY_RANDOM_SEED  :", RANDOM_SEED, "/", LEGACY_RANDOM_SEED, "\n")
cat("    : seeds of the random sampling (background + Otsu).\n")
cat("      Changing them produces different training pools. The final\n")
cat("      metrics vary little (~ +/- 0.01 in F1). Useful for robustness tests.\n")

cat("\n  3.4 Between-year temporal filter\n")
cat("    CURRENTYEAR_PREYEAR_OVERLAP_THR   :", CURRENTYEAR_PREYEAR_OVERLAP_THR, "\n")
cat("    : spatial-overlap threshold against previous-year polygons.\n")
cat("      Above 70% the polygon is treated as a 'persistent fire' (suspicious).\n")
cat("      Raising it to 0.85 = stricter filter, more drops.\n")

cat("    CURRENTYEAR_HOTSPOT_DENSITY_THR   :", CURRENTYEAR_HOTSPOT_DENSITY_THR, "\n")
cat("    : minimum hotspot density (per ha) to validate a polygon.\n")
cat("      Raising it demands more thermal evidence. Only applies from 1995.\n")

cat("    CURRENTYEAR_TEMPORAL_PENALTY_FLOOR:", CURRENTYEAR_TEMPORAL_PENALTY_FLOOR, "\n")
cat("    : floor of the temporal penalty score (it never goes below this).\n")
cat("      Prevents penalty=0 from wiping the polygon out entirely.\n")

cat("\n  3.5 Otsu pipeline: main parameters\n")
cat("    LEGACY_OTSU_MODE                  :", LEGACY_OTSU_MODE, "\n")
cat("    : applies Otsu only over the pixels Corine deems burnable.\n")
cat("      Without this, Otsu would see water/urban and the thresholds\n")
cat("      would come out distorted.\n")

cat("    LEGACY_OTSU_THRESHOLD             :", LEGACY_OTSU_THRESHOLD, "\n")
cat("    : minimum RBR threshold for an Otsu patch to be considered.\n")
cat("      A guard against negative spectral noise.\n")

cat("    LEGACY_REFERENCE_OTSU_THRESHOLD   :", LEGACY_REFERENCE_OTSU_THRESHOLD, "\n")
cat("    : reference threshold used to scale the other Otsu thresholds.\n")

cat("    LEGACY_SAMPLE_N                   :", LEGACY_SAMPLE_N, "\n")
cat("    : cap on the otsu_patch_residual pool. Without it, the ~13,000\n")
cat("      candidate patches would swamp the training set.\n")

cat("\n  3.6 Otsu pipeline: confidence thresholds\n")
cat("    OTSU_KEEP_HI / OTSU_DROP_LO       :", OTSU_KEEP_HI, "/", OTSU_DROP_LO, "\n")
cat("    : probability thresholds for classifying Otsu patches.\n")
cat("      KEEP_HI=0.45 -> high-confidence burned. DROP_LO=0.15 -> low.\n")
cat("      Raising KEEP_HI = stricter about accepting something as burned.\n")

cat("\n  3.7 Pipeline toggles\n")
cat("    DO_POOLS / FOLDS / FEATURES / MODEL:",
    DO_POOLS, "/", DO_FOLDS, "/", DO_FEATURES, "/", DO_MODEL, "\n")
cat("    : run STEP A (pools) / B1-B2 (folds) / B3 (features) /\n")
cat("      C-D (model). In 0.5.0 all hardcoded to TRUE; the refactor\n")
cat("      will allow partial reruns (e.g. regenerate the model only).\n")

cat("\n  3.8 I/O control\n")
cat("    LEGACY_REUSE_EXISTING             :", LEGACY_REUSE_EXISTING, "\n")
cat("    : reuses Otsu outputs from previous runs when they exist.\n")
cat("      FALSE = always regenerate from scratch (slower, more robust).\n")

cat("    LEGACY_WRITE_OUTPUT               :", LEGACY_WRITE_OUTPUT, "\n")
cat("    : writes intermediate Otsu outputs to disk. FALSE = in memory\n")
cat("      only (not inspectable after the fact).\n")

cat("    UNB_VERBOSE                       :", UNB_VERBOSE, "\n")
cat("    : detailed STEP A4 log (the '[unb] ...' messages).\n")


cat("\n[4] TUTORIAL OUTPUT DIRECTORY\n")
cat("    ", TUTORIAL_RESULT_DIR, "\n")

cat("\nBLOCK 0 complete. ",
    "13 inputs + 18 methodological constants declared and verified.\n")

