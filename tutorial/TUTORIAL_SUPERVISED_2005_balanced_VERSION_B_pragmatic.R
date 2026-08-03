# =============================================================================
# TUTORIAL_SUPERVISED_2005_balanced_VERSION_B_pragmatic.R
# =============================================================================
#
# OBJETIVO
#   Reproduce, step by step, the OtsuFire 0.5.0 supervised pipeline for
#   year 2005, scenario balanced, in the CANONICAL BASELINE configuration,
#   calling the internal functions (OtsuFire:::internal_function) instead
#   de a la unica funcion publica `run_oneyear_supervised_pipeline()`.
#
# CANONICAL BASELINE CONFIGURATION (confirmed by Natalia, see HANDOFF
# §N+23):
#   - feature_weights              = NULL  (uniform weight = 1.0)
#   - feature_whitelist_override   = NULL  (the package's natural canon: 51)
#   - contextual_exclusion_to_burned_ratio   = 0.25
#   - spectral_hard_negative_to_burned_ratio = 1.0
#   - random_to_burned_ratio                 = 1.0
#   - otsu_unburned_to_burned_ratio          = 1.0
#   - reuse_upstream                         = FALSE
#   In 2005 (post-MODIS) the model uses the 51 nominal features
#   (including the 13 hs_* ones).
#
# STRATEGY OF THIS TUTORIAL (VERSION B - pragmatic)
#   1) Build cfg with the public function build_supervised_burned_config().
#      This fills in the ~30 internal Otsu parameters automatically
#      legacy y unburned, evitando reproducir el dispatcher a mano.
#   2) Asignar tool_paths a cfg.
#   3) Call the internal functions (OtsuFire:::xxx) by hand, in the
#      exact order internal-sup-orchestrator.R follows, stopping
#      after each STEP to inspect the outputs.
#
# OUTPUT DIRECTORY
#   So as not to trample the nightly agent's run, this tutorial writes to:
#     Results/2005/Min_Min/SUPERVISED_TUTORIAL/balanced/
#   en vez de en Results/2005/Min_Min/SUPERVISED/balanced/.
#
# HOW TO USE
#   - Clean restart of R (Ctrl+Shift+F10).
#   - Run it block by block (select and Ctrl+Enter).
#   - Inspect the variables each block leaves in memory.
#   - Compare against the wrapper's output at the end.
#
# REFERENCES INTO THE PACKAGE CODE
#   - Wrapper publico:           R/supervised-run.R
#   - Orchestrator (corazon):    R/internal-sup-orchestrator.R
#                                run_supervised_pipeline() linea 719
#                                STEP A1-A6, B1-B3, C0-C4
#   - Funciones internas:        R/internal-sup-*.R y R/supervised-*.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
})

# =============================================================================
# BLOCK 0 - SETUP: LOADING THE PACKAGE AND PATHS
# =============================================================================
#
# WHAT IT DOES
#   Loads the OtsuFire package from the source directory (without installing
#   it) via pkgload::load_all(). Defines the global project paths.
#
# WHY
#   load_all() is the standard way to work on a package under development:
#   it sees the NAMESPACE exports but also allows access to internal
#   functions via OtsuFire:::xxx.
#
# OUTPUTS EN MEMORIA
#   - PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, etc. (paths)
#   - YEAR=2005L, SCENARIO="balanced"
#   - INTERNAL_DECISIONS, CHANGE_INDEX, HOTSPOTS (paths a inputs criticos)
#
# =============================================================================
# =============================================================================
# BLOCK 0 - SETUP: LOADING THE PACKAGE AND PATHS
# =============================================================================
#
# WHAT IT DOES
#   Loads the OtsuFire package from the source directory (without installing
#   it) via pkgload::load_all(). Defines the global project paths and the
#   4 critical input paths of the experiment (Y/E = year/scenario):
#     - internal_decisions.gpkg          (deterministic stage)
#     - mosaico summer (RBR + DOY)       (composite ano)
#     - mosaico autumn-winter            (composite delayed)
#     - hotspots geojson                 (post-MODIS, NULL en pre-MODIS)
#
# WHY
#   load_all() is the standard way to work on a package under development:
#   it sees the NAMESPACE exports but also allows access to internal
#   functions via OtsuFire:::xxx.
#
# OUTPUTS EN MEMORIA
#   - PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, TUTORIAL_RESULT_DIR (paths)
#   - YEAR=2005L, SCENARIO="balanced", RESULT_NAME="Min_Min"
#   - INTERNAL_DECISIONS, CHANGE_INDEX, RBR_AUTUMN, HOTSPOTS
#     (paths to the 4 critical inputs of the combo)
#
# =============================================================================

# =============================================================================
# BLOCK 0 - SETUP: PACKAGE AND THE 14 SUPERVISED PIPELINE INPUTS
# =============================================================================
#
# WHAT IT DOES
#   Loads the OtsuFire package from source (without installing it) via
#   pkgload::load_all() and declares the 14 critical input paths of the
#   supervised module, organised by category.
#
# WHY 14 INPUTS Y NO 4
# ------------------------
# The config builder build_supervised_burned_config() accepts only 4 paths
# of inputs (internal_decisions, change_index, hotspots, and the cosmetic
# delayed_change_index). But the real pipeline reads 14 different files:
#
#   - 4 are passed to the config builder (year/scenario-dependent + delayed)
#   - 10 are built by the orchestrator by convention from data_base,
#     composite_base, target_year y corine_year
#
# This inconsistency is documented in HANDOFF §N+25 as technical debt
# in the package. The post-paper refactor will wire it up: all 14 will
# come in through the config builder.
#
# For the tutorial, ALL 14 are declared explicitly up front. That way:
#   1. Every dataset the pipeline will touch is visible at a glance.
#   2. Validation happens up front (in <1 second) instead of failing 5-10
#      min later if some file is missing.
#   3. When the refactor wires the 10 implicit ones into cfg, this tutorial
#      will keep working: only the config builder call changes,
#      back in BLOCK 1.
#
#
# OUTPUTS EN MEMORIA
# ------------------
#   Paths de proyecto:
#     PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, TUTORIAL_RESULT_DIR
#
#   Combo objetivo:
#     YEAR=2005L, SCENARIO="balanced", RESULT_NAME="Min_Min", CORINE_YEAR
#
#   Inputs (14):
#     A) Year/scenario-dependent (4):
#        INTERNAL_DECISIONS, CHANGE_INDEX, RBR_AUTUMN, HOTSPOTS
#     B) Validacion externa year-dependent (2):
#        EFFIS_CA_TIF, EFFIS_CA_SHP
#     C) Corine quinquenal (4):
#        CORINE_RASTER, CORINE_STRATA, CORINE_LUT, BURNEABLE_MASK
#     D) Estaticos (4):
#        PENINSULA_SHP, TOPO, MASK_TIF, MASK_SHP
#
#
# REFERENCIA EN EL PAQUETE
# ------------------------
#   internal-sup-orchestrator.R lineas 906-934 - construccion de paths
#   internal-sup-orchestrator.R lineas 939-950 - validacion stopifnot
#   internal-sup-unburned-deterministic.R linea 170 - burneable mask
#   internal-sup-unburned-legacy.R linea 655 - burneable mask
#
# =============================================================================

cat("\n========== BLOCK 0: SETUP ==========\n")

# --- Loading the package -------------------------------------------------------
PKG_ROOT  <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"

if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

cat("OtsuFire version:", as.character(utils::packageVersion("OtsuFire")), "\n")
stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

# --- Project paths ------------------------------------------------------
DATA_BASE <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA"
COMPOSITE <- file.path(DATA_BASE, "Imagery", "Composites_90m")
RESULTS   <- file.path(DATA_BASE, "Results")

# --- Combo objetivo ----------------------------------------------------------
YEAR        <- 2005L
SCENARIO    <- "balanced"
RESULT_NAME <- "Min_Min"

# Resolution of the five-yearly Corine year via the package's internal function.
# Mapeo (segun el codigo de OtsuFire):
#   1985-2002 -> CLC2000 ; 2003-2008 -> CLC2006 ; 2009-2014 -> CLC2012 ;
#   2015+    -> CLC2018
# For 2005, CORINE_YEAR resolves to 2006.
CORINE_YEAR <- OtsuFire:::get_corine_year(YEAR)
cat(sprintf("Combo: year=%d, scenario=%s, corine=%s\n",
            YEAR, SCENARIO, CORINE_YEAR))

# === A. INPUTS YEAR/SCENARIO-DEPENDENT (4) ===================================

# A1. Output of the deterministic module, main input of the supervised one.
INTERNAL_DECISIONS <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "05_DECISIONS", "internal_decisions.gpkg"
)

# A2. Summer mosaic (RBR + DOY of the fire summer).
CHANGE_INDEX <- file.path(COMPOSITE, RESULT_NAME,
                          sprintf("MinMin_%d_mosaic_res90m.tif", YEAR))

# A3. Mosaico autumn-winter (delayed change index, fuente de features G2_RBR_AW).
#     Cosmetic in cfg under 0.5.0 (HANDOFF §N+25), but declared here anyway.
RBR_AUTUMN <- file.path(COMPOSITE, "Autumn",
                        sprintf("mean_mean_%d_mosaic.tif", YEAR))

# A4. Hotspots MODIS. Pre-MODIS (<2000): NULL.
HOTSPOTS <- file.path(DATA_BASE, "Hotspots",
                      sprintf("hotspots_iberia_%d.geojson", YEAR))


# === B. VALIDACION EXTERNA YEAR-DEPENDENT (2) =================================

# B1. Effis-CA summer fires filtered by the burnable mask (raster).
EFFIS_CA_TIF <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.tif", YEAR)
)

# B2. The same thing in shp format.
EFFIS_CA_SHP <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.shp", YEAR)
)


# === C. CORINE QUINQUENAL (4) =================================================

CORINE_RASTER <- file.path(DATA_BASE, "Corine_Masks",
                           sprintf("CLC_%s_peninsula.tif", CORINE_YEAR))

CORINE_STRATA <- file.path(DATA_BASE, "Corine_Masks", "STRATA",
                           sprintf("strata_CLC_%s_res30.tif", CORINE_YEAR))


# C3. Look-up table for the strata.
CORINE_LUT <- file.path(DATA_BASE, "Corine_Masks", "LUT",
                        "lut_full_strata8_v1.csv")

# C4. Mascara binaria de superficie quemable derivada de Corine.
#     NOTE: this input is NOT validated by the orchestrator up front.
#     It is only validated inside internal-sup-unburned-deterministic.R line 170
#     (when STEP A4 runs). If it is missing, the pipeline fails 5-10 min later
#     after starting. Documented in HANDOFF §N+25.
BURNEABLE_MASK <- file.path(
  DATA_BASE, "Corine_Masks",
  sprintf("burneable_mask_binary_corine_%s_ETRS89.tif", CORINE_YEAR)
)


# === D. ESTATICOS (4) =========================================================

# D1. Frontera de la Peninsula Iberica (polygon shapefile).
PENINSULA_SHP <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")

# D2. Topography: stack with DEM and slope.
TOPO <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

# D3. Study-area mask in EPSG:3035 (raster).
MASK_TIF <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")

# D4. Idem en shp.
MASK_SHP <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.shp")


# === GLOBAL VERIFICATION (all 14 files must exist) ===========================

stopifnot(
  # A) Year-dependent
  file.exists(INTERNAL_DECISIONS),
  file.exists(CHANGE_INDEX),
  file.exists(RBR_AUTUMN),
  file.exists(HOTSPOTS),
  # B) Validacion externa
  file.exists(EFFIS_CA_TIF),
  file.exists(EFFIS_CA_SHP),
  # C) Corine
  file.exists(CORINE_RASTER),
  file.exists(CORINE_STRATA),
  file.exists(CORINE_LUT),
  file.exists(BURNEABLE_MASK),
  # D) Estaticos
  file.exists(PENINSULA_SHP),
  file.exists(TOPO),
  file.exists(MASK_TIF),
  file.exists(MASK_SHP)
)

cat("\nLos 14 inputs criticos verificados:\n")
cat("\n  A) Year/scenario-dependent (4):\n")
cat("    internal_decisions   :", basename(INTERNAL_DECISIONS), "\n")
cat("    change_index (summer):", basename(CHANGE_INDEX), "\n")
cat("    delayed (autumn)     :", basename(RBR_AUTUMN), "\n")
cat("    hotspots             :", basename(HOTSPOTS), "\n")

cat("\n  B) Validacion externa (2):\n")
cat("    Effis-CA tif         :", basename(EFFIS_CA_TIF), "\n")
cat("    Effis-CA shp         :", basename(EFFIS_CA_SHP), "\n")

cat("\n  C) Corine quinquenal (4) [CORINE_YEAR=", CORINE_YEAR, "]:\n", sep="")
cat("    Corine raster        :", basename(CORINE_RASTER), "\n")
cat("    Corine strata        :", basename(CORINE_STRATA), "\n")
cat("    Corine LUT           :", basename(CORINE_LUT), "\n")
cat("    Burneable mask       :", basename(BURNEABLE_MASK), "\n")

cat("\n  D) Estaticos (4):\n")
cat("    Peninsula shp        :", basename(PENINSULA_SHP), "\n")
cat("    Topografia           :", basename(TOPO), "\n")
cat("    Mask studyarea tif   :", basename(MASK_TIF), "\n")
cat("    Mask studyarea shp   :", basename(MASK_SHP), "\n")


# === TUTORIAL OUTPUT DIRECTORY ===============================================

# So as not to trample the real run, the tutorial outputs are redirected to
# a parallel folder, SUPERVISED_TUTORIAL/.
TUTORIAL_RESULT_DIR <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "SUPERVISED_TUTORIAL", SCENARIO
)
dir.create(TUTORIAL_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("\nTutorial output directory:\n  ", TUTORIAL_RESULT_DIR, "\n")

# =============================================================================
# BLOCK 1 - CONFIG: build_supervised_burned_config()
# =============================================================================
#
# WHAT IT DOES
#   Builds a cfg object of class 'otsufire_supervised_burned_config' holding
#   ALL the paths, parameters and options the pipeline needs.
#   The function runs NOTHING: it only configures. Execution is done by
#   run_oneyear_supervised_pipeline(cfg).
#
#
# PACKAGE DESIGN: THE CONFIG + RUN PATTERN
# ---------------------------------------
# OtsuFire sigue el patron clasico de software cientifico (sklearn, MLflow,
# targets...): SEPARATE configuration from execution. There are two functions
# publicas centrales en el modulo supervised:
#
#   1) build_supervised_burned_config()  -> "what I want to do"
#   2) run_oneyear_supervised_pipeline() -> "hazlo"
#
# Advantages of this pattern:
#   - EARLY VALIDATION. The builder checks paths and types. If something is
#     wrong the error appears in <1s rather than 2 hours into the pipeline.
#   - REPRODUCIBILITY. cfg is an immutable R object that can be serialised
#     with saveRDS(cfg, "baseline_experiment.rds"). Anyone holding that
#     rds plus the data reproduces the exact experiment.
#   - PROGRAMACION DE EXPERIMENTOS. configs <- list(baseline=..., e1=...);
#     lapply(configs, run_oneyear_supervised_pipeline). Clean.
#
#
# ANATOMY OF THE cfg OBJECT (14 components in this case)
# -----------------------------------------------------
# str(cfg, max.level = 1) revela:
#
#   $ scenario                  - chr "balanced"     (operational scenario)
#   $ target_year               - int 2005           (ano objetivo)
#   $ inputs                    - List of 5          (paths a inputs)
#   $ run_name                  - chr "Min_Min"      (name of the run)
#   $ output_dir                - chr ".../Results"  (raiz de outputs)
#   $ output_routes             - List of 19         (paths derivados)
#   (NOTA: burned_like_registry_path y negative_pool_policy fueron
#    removed from the package on 2026-06-05; the burned-like registry no
#    longer exists and the negatives policy is always all_sources, implicit.)
#   $ min_burned_pool_n         - int 5              (guard minimo positivos)
#   $ engine_root               - NULL               (legacy, unused)
#   $ supervised_engine_root    - NULL               (legacy, unused)
#   $ scripts_root              - NULL               (legacy, unused)
#   $ tool_paths                - List of 4          (Anaconda3 binaries)
#   $ options                   - List of 11         (resto de opciones)
#
#
# THE 5 FIELDS OF cfg$inputs AND THEIR STATUS IN 0.5.0
# ------------------------------------------------
# cfg stores 5 input paths, but NOT all of them are consumed by the
# pipeline en la version actual (0.5.0). Conviene saberlo:
#
#   $ internal_decisions   -> IS consumed (STEP A1)
#   $ change_index         -> IS consumed (STEP A4 + feature extraction)
#   $ delayed_change_index -> COSMETIC in 0.5.0; el orchestrator construye
#                             the path by convention from composite_base
#                             see HANDOFF §N+25 for details
#   $ hotspots             -> IS consumed (STEP B3)
#   $ reference_burned_map -> RESERVED for future external validation,
#                             not wired up yet in 0.5.0
#
# IMPORTANT: the config builder accepts only 4 input paths; the pipeline
# actually reads 14 (HANDOFF §N+25). The other 10 are built by the
# orchestrator by convention. BLOCK 0 declares all 14 explicitly so they
# can be validated up front. Here in BLOCK 1 only the 4 the public API
# accepts are passed to the builder.
#
# The tutorial passes delayed_change_index to the builder even though it is
# cosmetic, for two reasons:
#   1. When you inspect cfg$inputs you see every real input of the
#      experiment. If someone serialises cfg with saveRDS() to
#      reproducibilidad, el path autumn quedara documentado.
#   2. When the post-paper refactor happens (HANDOFF §N+25) the field
#      will be consumed. Passing the path now makes the script
#      forward-compatible: same call, different package behaviour,
#      with no change in the user's script.
#
#
# CLASE S3
# --------
# cfg has class = c("otsufire_supervised_burned_config", "list").
# It is a lightweight S3 class: a labelled list. When
# run_oneyear_supervised_pipeline() receives cfg, the first thing it does is
# check inherits(config, "otsufire_supervised_burned_config"). Passing a
# plain list fails with a clear message. It is a type contract.
#
#
# METHODOLOGICAL DECISIONS FOR THE BASELINE
# --------------------------------------
#   - Negative pool (always all_sources, implicit; no longer an option)
#       Uses the 4 sub-pools: internal_keep_qc (burned), deterministic_drop_hard,
#       random_burnable_background, otsu_patch_residual. (La antigua
#       alternativa "deterministic_direct" fue eliminada el 2026-06-05.)
#
#   - legacy_otsu_mode = "burnable_only"
#       Otsu is applied only over the pixels the Corine mask deems
#       burnable. Without this, Otsu would see water, urban, etc. and the
#       saldrian distorsionados.
#
#   - legacy_sample_n = 2000L
#       Cap on the otsu_patch_residual pool. Without it the ~13,000 candidate
#       candidatos saturarian el training.
#
#   - legacy_otsu_threshold = 0
#       Discards Otsu candidates with RBR < 0. A minimal guard against noise.
#
#
# REFERENCIA EN EL PAQUETE
# ------------------------
#   R/supervised-config.R (15 KB) - definicion de build_supervised_burned_config
#   NAMESPACE linea 5             - export(build_supervised_burned_config)
#   NAMESPACE linea 1             - S3method(print, otsufire_supervised_burned_config)
#
# =============================================================================

cat("\n========== BLOCK 1: CONFIG ==========\n")

# RBR_AUTUMN, INTERNAL_DECISIONS, CHANGE_INDEX, HOTSPOTS are already defined
# and verified back in BLOCK 0. They are used directly here.


# --- Construir el cfg -------------------------------------------------------

cfg <- build_supervised_burned_config(
  run_label                 = SCENARIO,
  internal_decisions        = INTERNAL_DECISIONS,
  change_index              = CHANGE_INDEX,        # mosaico summer (RBR + DOY)
  delayed_change_index      = RBR_AUTUMN,          # autumn-winter (cosmetic, §N+25)
  hotspots                  = HOTSPOTS,
  target_year               = YEAR,
  output_dir                = RESULTS,
  run_name                  = RESULT_NAME,
  options = list(
    data_base                       = DATA_BASE,
    composite_base                  = COMPOSITE,
    result_name                     = RESULT_NAME,
    # negative-pool policy is always all_sources now (implicit; only mode)
    legacy_otsu_mode                = "burnable_only",
    legacy_otsu_threshold           = 0,
    legacy_reference_otsu_threshold = 100,
    legacy_sample_n                 = 2000L,
    legacy_reuse_existing           = TRUE,
    legacy_write_output             = TRUE,
    unb_verbose                     = TRUE
  )
)

# Tool paths (Anaconda on this machine). They are assigned AFTER the builder
# because they are machine-specific and not part of the experiment
# configuration. That keeps cfg portable: for a collaborator on another
# machine only tool_paths change.
cfg$tool_paths$python_exe              <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/python.exe"
cfg$tool_paths$gdal_polygonize_script  <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Scripts/gdal_polygonize.py"
cfg$tool_paths$gdalwarp_path           <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/gdalwarp.exe"
cfg$tool_paths$ogr2ogr_exe             <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/ogr2ogr.exe"

stopifnot(all(sapply(cfg$tool_paths, file.exists)))
cat("Tool paths OK.\n")


# --- INSPECTING cfg -----------------------------------------------------

cat("\ncfg class:\n")
print(class(cfg))
# Resultado esperado:
#   [1] "otsufire_supervised_burned_config" "list"
#
# The S3 class ("otsufire_supervised_burned_config") comes first. It lets
# run_oneyear_supervised_pipeline() validate that this is a legitimate cfg.
# The second ("list") is R's base class, which lets cfg behave like a
# list: cfg$inputs, cfg$options$legacy_sample_n, etc.

cat("\nTop-level structure of cfg:\n")
str(cfg, max.level = 1)

# --- cfg$inputs: paths a INPUTS criticos -----------------------------------
#
# These are the paths to the files the pipeline reads. You passed them as
# arguments to the builder; here is where they are stored.
#
# Reminder (HANDOFF §N+25): the config builder accepts only 4 paths,
# but the pipeline reads 14 different files. The other 10 (Corine, topo,
# mask, peninsula, EFFIS validation, burnable mask) are built by the
# orchestrator by convention. Those do NOT appear here.

cat("\ncfg$inputs:\n")
str(cfg$inputs, max.level = 1)
# Resultado:
#   $ internal_decisions  : List of 3   <- internal_decisions.gpkg + metadata
#                                          CONSUMIDO en STEP A1
#   $ change_index        : List of 3   <- mosaico summer + metadata
#                                          CONSUMIDO en STEP A4 + extraction
#   $ delayed_change_index: List of 3   <- autumn-winter mosaic + metadata
#                                          COSMETIC in 0.5.0 (HANDOFF §N+25)
#                                          the autumn path is built by the
#                                          orchestrator by convention in
#                                          composite_base/Autumn/
#   $ hotspots            : List of 3   <- hotspots.geojson + metadata
#                                          CONSUMED in STEP B3
#   $ reference_burned_map: NULL        <- reserved for the future, not wired up


# --- cfg$options: the operational knobs -------------------------------
#
# These are the parameters you control to define Baseline versus
# experiments. 11 are passed to the builder; the rest get defaults.

cat("\ncfg$options (11 passed explicitly; the rest are defaults):\n")
str(cfg$options, max.level = 1)
# Resultado:
#   $ data_base                      : "C:/.../1_DATA"
#                                       USED by the orchestrator to
#                                       build 9 implicit paths (§N+25)
#   $ composite_base                 : "C:/.../Composites_90m"
#                                       USED to build the autumn path
#                                       and the summer mosaic (§N+25)
#   $ result_name                    : "Min_Min"
#   (negative_pool_policy was removed; all_sources is always implicit)
#   $ legacy_otsu_mode               : "burnable_only"
#   $ legacy_otsu_threshold          : 0
#   $ legacy_reference_otsu_threshold: 100
#   $ legacy_sample_n                : 2000
#   $ legacy_reuse_existing          : TRUE
#   $ legacy_write_output             : TRUE
#   $ unb_verbose                    : TRUE


# --- cfg$output_routes: paths DERIVADOS ------------------------------------
#
# El builder calcula estos solos a partir de output_dir + target_year +
# result_name + scenario. That is why there are 19 paths without you having
# passed any. The SUPERVISED/<scenario>/<NN_FOLDER>/ layout is canonical
# to the package.

cat("\ncfg$output_routes (19 paths derivados automaticamente):\n")
str(cfg$output_routes, max.level = 1)
# Subcarpetas: 01_POOLS, 02_FOLDS, 03_FEATURES, 04_MATRIX, 05_OOF,
# 07_FINAL_MODEL_V2, 08_SCORED, 09_FINAL_MAP, 11_CONSISTENCY_CHECKS, 99_LOGS.


# --- REDIRECTING OUTPUTS FOR THE TUTORIAL --------------------------------
#
# IMPORTANT: the orchestrator writes ALL outputs under
# cfg$output_routes$base. Left as it is, it would trample the real
# SUPERVISED/balanced/ run. For the tutorial it is redirected to a folder
# alternativa SUPERVISED_TUTORIAL/.

cfg$output_routes$base       <- TUTORIAL_RESULT_DIR
cfg$output_routes$result_dir <- TUTORIAL_RESULT_DIR

cat("\nBLOCK 1 complete. cfg is ready to feed the pipeline.\n")



cat("\n========== BLOCK 2: DIRECTORIES AND RASTERS ==========\n")

# Directory structure (same as the orchestrator, line 779)
dirs <- list(
  `01_POOLS`          = file.path(TUTORIAL_RESULT_DIR, "01_POOLS"),
  `02_FOLDS`          = file.path(TUTORIAL_RESULT_DIR, "02_FOLDS"),
  `03_FEATURES`       = file.path(TUTORIAL_RESULT_DIR, "03_FEATURES"),
  `04_MATRIX`         = file.path(TUTORIAL_RESULT_DIR, "04_MATRIX"),
  `05_OOF`            = file.path(TUTORIAL_RESULT_DIR, "05_OOF"),
  `07_FINAL_MODEL_V2` = file.path(TUTORIAL_RESULT_DIR, "07_FINAL_MODEL_V2"),
  `08_SCORED`         = file.path(TUTORIAL_RESULT_DIR, "08_SCORED"),
  `09_FINAL_MAP`      = file.path(TUTORIAL_RESULT_DIR, "09_FINAL_MAP"),
  `99_LOGS_EMPTY`     = file.path(TUTORIAL_RESULT_DIR, "99_LOGS_EMPTY")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- Paths a rasters comunes (orchestrator linea 906) -----------------------
peninsula_shp <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")
topo_path     <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

corine_year   <- OtsuFire:::get_corine_year(YEAR)
corine_path   <- file.path(DATA_BASE, "Corine_Masks",
                           paste0("CLC_", corine_year, "_peninsula.tif"))

mask_tif_path <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")

rbr_aw_tif <- file.path(COMPOSITE, "Autumn",
                        paste0("mean_mean_", YEAR, "_mosaic.tif"))

stopifnot(file.exists(peninsula_shp), file.exists(topo_path),
          file.exists(corine_path), file.exists(mask_tif_path),
          file.exists(rbr_aw_tif))

cat("Inputs auxiliares verificados.\n")
cat("  Corine year:", corine_year, "\n")

# --- Load rasters and align them to the template (orchestrator line 952) ---
cat("\nCargando y alineando rasters...\n")

topo       <- terra::rast(topo_path)
rbr_stack  <- terra::rast(CHANGE_INDEX)
rbr_summer <- rbr_stack[[1]]
doy_post   <- rbr_stack[[2]]

template_r <- rbr_summer

doy_post <- OtsuFire:::align_to_template(
  r = doy_post, template = template_r, method = "near",
  name = "doy_post", verbose = FALSE
)

dem_r <- OtsuFire:::align_to_template(
  r = topo[[1]], template = template_r, method = "bilinear",
  name = "dem", verbose = FALSE
)

slope_r <- OtsuFire:::align_to_template(
  r = topo[[2]], template = template_r, method = "bilinear",
  name = "slope", verbose = FALSE
)

corine_r <- OtsuFire:::align_to_template(
  r = terra::rast(corine_path), template = template_r, method = "near",
  name = "corine_r", verbose = FALSE
)

rbr_aw <- OtsuFire:::align_to_template(
  r = terra::rast(rbr_aw_tif)[[1]], template = template_r, method = "bilinear",
  name = "rbr_aw", verbose = FALSE
)

# --- Hotspots ---------------------------------------------------------------
hotspots_sf_base <- sf::read_sf(HOTSPOTS)
cat(sprintf("Hotspots cargados: %d puntos\n", nrow(hotspots_sf_base)))

# --- INSPECCION ---
cat("\nResumen rasters alineados:\n")
cat(sprintf("  rbr_summer: %d x %d, res %.0f m\n",
            nrow(rbr_summer), ncol(rbr_summer), terra::res(rbr_summer)[1]))
cat(sprintf("  doy_post:   %d x %d\n", nrow(doy_post), ncol(doy_post)))
cat(sprintf("  dem:        %d x %d\n", nrow(dem_r), ncol(dem_r)))
cat(sprintf("  slope:      %d x %d\n", nrow(slope_r), ncol(slope_r)))
cat(sprintf("  corine:     %d x %d\n", nrow(corine_r), ncol(corine_r)))
cat(sprintf("  rbr_aw:     %d x %d\n", nrow(rbr_aw), ncol(rbr_aw)))

# =============================================================================
# BLOCK 3 - STEP A1: READ internal_decisions.gpkg
# =============================================================================
#
# WHAT IT DOES
#   Reads the GPKG produced by the deterministic phase, carrying the
#   per-polygon classification: keep / drop / review.
#
# INPUTS  (de disco)
#   internal_decisions.gpkg, layer "internal_decisions"
#
# OUTPUTS (en memoria)
#   internal_sf - sf with every classified polygon of the year
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R linea 1026 (STEP A1)
#
# =============================================================================

cat("\n========== BLOCK 3: STEP A1 - READ internal_decisions ==========\n")

internal_sf <- sf::read_sf(INTERNAL_DECISIONS, layer = "internal_decisions") |>
  dplyr::mutate(class = as.character(class_final))

# The package functions sanitize_polygons and ensure_area_ha live in the
# internal safety kit. They are called from there.
sfkit <- OtsuFire:::make_sf_safety_kit(
  dirs = dirs, result_dir = TUTORIAL_RESULT_DIR, target_year = YEAR
)
sanitize_polygons    <- sfkit$sanitize_polygons
ensure_area_ha       <- sfkit$ensure_area_ha
drop_empty_sf        <- sfkit$drop_empty
check_sf             <- sfkit$check_sf
check_sf_if_nonempty <- sfkit$check_sf_if_nonempty
to_crs_safe          <- sfkit$to_crs_safe
safe_read_gpkg       <- sfkit$safe_read_gpkg
safe_write_gpkg      <- sfkit$safe_write_gpkg

internal_sf <- internal_sf |>
  sanitize_polygons() |>
  ensure_area_ha()

internal_sf <- drop_empty_sf(internal_sf, tag = "internal_sf",
                             dump_dir = dirs$`99_LOGS_EMPTY`)
check_sf(internal_sf, "internal_sf")

crs_master <- sf::st_crs(internal_sf)

cat(sprintf("internal_sf: %d poligonos\n", nrow(internal_sf)))
cat("Distribucion de clases:\n")
print(table(internal_sf$class))

# =============================================================================
# BLOCK 4 - STEP A2 + A3: AUDIT and POOLS
# =============================================================================
#
# WHAT IT DOES
#   A2: marks every polygon with source = "internal".
#   A3: applies audit_deterministic_pools(), which validates and reclassifies
#       polygones marginales, y construye 3 pools:
#         burned_pool   (class == "keep")    -> positivos
#         review_pool   (class == "review")  -> doubtful, not used in training
#         scoring_pool  (all)                -> the universe to score at the end
#
# INPUTS  (memoria)
#   internal_clean (= internal_sf con source="internal")
#
# OUTPUTS (memoria + disco)
#   poolsA$internal_qc, $burned_pool, $review_pool, $scoring_pool
#   GPKG de QA en dirs$`01_POOLS`
#
# REFERENCIA EN EL PAQUETE
#   audit_deterministic_pools() en R/internal-sup-qa-pools.R
#   internal-sup-orchestrator.R linea 1065 (STEP A3)
#
# =============================================================================

cat("\n========== BLOCK 4: STEP A2 + A3 - AUDIT + POOLS ==========\n")

# A2: marcar source
internal_clean <- internal_sf |>
  dplyr::mutate(source = "internal") |>
  sanitize_polygons() |>
  ensure_area_ha()

# A3: auditoria + construir pools
qa_det <- OtsuFire:::audit_deterministic_pools(
  internal_sf  = internal_clean,
  out_dir      = dirs$`01_POOLS`,
  prefix       = sprintf("%d_%s_deterministic_pool_qa", YEAR, SCENARIO),
  save_outputs = TRUE,
  verbose      = TRUE,
  raw_class_col = "class"
)

internal_qc <- qa_det$audited_sf |>
  dplyr::mutate(
    raw_class = as.character(raw_class),
    class     = as.character(class_audited)
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

# Construir burned_pool (positivos)
burned_pool <- internal_qc |>
  dplyr::filter(class == "keep") |>
  dplyr::mutate(
    source   = "internal_keep_qc",
    class    = "burned",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_B_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

burned_pool <- drop_empty_sf(burned_pool, tag = "burned_pool",
                             dump_dir = dirs$`99_LOGS_EMPTY`)

# Construir review_pool
review_pool <- internal_qc |>
  dplyr::filter(class == "review") |>
  dplyr::mutate(
    source   = "internal_review_qc",
    class    = "review",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_R_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

review_pool <- drop_empty_sf(review_pool, tag = "review_pool",
                             dump_dir = dirs$`99_LOGS_EMPTY`)

# Build scoring_pool (every polygon to be scored at the end)
scoring_pool <- internal_qc |>
  dplyr::mutate(
    class    = as.character(raw_class),
    source   = "deterministic_all_qc",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_D_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

scoring_pool <- drop_empty_sf(scoring_pool, tag = "scoring_pool",
                              dump_dir = dirs$`99_LOGS_EMPTY`)

# --- INSPECCION ---
cat("\n--- Resultado A3 ---\n")
cat(sprintf("internal_qc:   %d poligonos\n", nrow(internal_qc)))
cat(sprintf("burned_pool:   %d (positivos)\n", nrow(burned_pool)))
cat(sprintf("review_pool:   %d (dudosos)\n", nrow(review_pool)))
cat(sprintf("scoring_pool:  %d (universo)\n", nrow(scoring_pool)))

# Sanity check: burned_pool cannot be too small
n_burned <- nrow(burned_pool)
if (n_burned < 5L) {
  stop(sprintf("burned_pool insuficiente (n=%d, requerido>=5)", n_burned))
}
cat(sprintf("Sanity check: n_burned=%d >= 5 OK\n", n_burned))

# =============================================================================
# BLOCK 5 - STEP A4: GENERATE UNBURNED (3 sub-pools)
# =============================================================================
#
# WHAT IT DOES
#   Builds the negative pools in two steps:
#     Parte 1 (build_unburned_from_deterministic_decisions):
#       - unburned_hard:   polygons classified as drop by the deterministic stage
#       - unburned_random: celdas de fondo quemable fuera de buffers
#     Parte 2 (build_unburned_from_legacy_pipeline):
#       - otsu_sampled:    parches Otsu unburned residual (cap 2000)
#   Then it combines both parts into unburned_final_raw.
#
# DECISIONES METODOLOGICAS
#   negative pool: always all_sources, implicit (the combined sub-pools)
#   exclude_buffer_m: separa negativos de positivos
#   n_random_cells: cuantos puntos random extraer
#   random_rbr_q: maximum RBR percentile for something to count as "background"
#   legacy_otsu_mode: "burnable_only" -> Otsu only where the mask allows
#   legacy_sample_n: 2000L (cap on the otsu_patch_residual pool)
#
# INPUTS
#   internal_decisions.gpkg, RBR raster, mascara, hotspots (informativo)
#
# OUTPUTS (memoria + disco GPKG)
#   unb$unburned_hard
#   unb$unburned_random
#   unb$unburned_final_raw  (all negatives combined, before sampling)
#
# REFERENCIA EN EL PAQUETE
#   build_unburned_from_deterministic_decisions() en R/internal-sup-unburned-deterministic.R
#   build_unburned_from_legacy_pipeline()         en R/internal-sup-unburned-legacy.R
#   internal-sup-orchestrator.R lineas 1180-1372 (STEP A4)
#
# IMPORTANTE - DURACION
#   This is the slowest stage of the pipeline when reuse_existing=FALSE, as
#   it has to run Otsu + polygonize in Python (~10-20 min).
#   For 2005/balanced an older run already left outputs on disk; with
#   legacy_reuse_existing = TRUE the polygonize step is reused and only
#   recalculan unburned_hard + unburned_random (1-2 min).
#
# =============================================================================

cat("\n========== BLOCK 5: STEP A4 - GENERATE UNBURNED ==========\n")
cat("(this stage can take 1-20 min depending on caches)\n")

# Output of the unburned-deterministic step
unb_out_gpkg <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "UNBURNED", sprintf("%d_%s_unburned.gpkg", YEAR, SCENARIO)
)
dir.create(dirname(unb_out_gpkg), recursive = TRUE, showWarnings = FALSE)

# Folder for legacy Otsu outputs (shared across scenarios)
legacy_unb_root_dir <- file.path(TUTORIAL_RESULT_DIR, "_LEGACY_UNBURNED")
dir.create(legacy_unb_root_dir, recursive = TRUE, showWarnings = FALSE)

t0_a4 <- Sys.time()

# --- A4 Parte 1: deterministic drops + random background -------------------
res_det <- OtsuFire:::build_unburned_from_deterministic_decisions(
  target_year             = YEAR,
  scenario_name           = SCENARIO,
  data_base               = DATA_BASE,
  result_name             = RESULT_NAME,
  composite_base          = COMPOSITE,
  exclude_buffer_m        = cfg$options$unb_excl_buffer_m %||% 500,
  n_random_cells          = cfg$options$unb_n_random_cells %||% 15000L,
  random_rbr_q            = cfg$options$unb_random_rbr_q %||% 0.5,
  random_seed             = cfg$options$unb_random_seed %||% 42L,
  random_patch_size_cells = cfg$options$unb_random_patch_size_cells %||% 11L,
  overwrite_output        = TRUE,
  out_gpkg                = unb_out_gpkg,
  verbose                 = TRUE
)

unburned_hard   <- to_crs_safe(res_det$unburned_hard,   crs_master) |>
  dplyr::mutate(source = "deterministic_drop_hard")
unburned_random <- to_crs_safe(res_det$unburned_random, crs_master) |>
  dplyr::mutate(source = "random_burnable_background")
det_final_raw   <- to_crs_safe(res_det$unburned_final,  crs_master) |>
  dplyr::mutate(source = as.character(source))
exclusion_buffer <- to_crs_safe(res_det$exclusion_buffer, crs_master)

cat(sprintf("\nA4-Parte1 OK: hard=%d, random=%d\n",
            nrow(unburned_hard), nrow(unburned_random)))

# --- A4 Parte 2: Otsu unburned residual ------------------------------------
res_otsu <- OtsuFire:::build_unburned_from_legacy_pipeline(
  target_year              = YEAR,
  scenario_name            = SCENARIO,
  data_base                = DATA_BASE,
  result_name              = RESULT_NAME,
  composite_base           = COMPOSITE,
  severity_raster_path     = CHANGE_INDEX,
  legacy_code_dir          = NULL,
  python_exe               = cfg$tool_paths$python_exe,
  gdal_polygonize_script   = cfg$tool_paths$gdal_polygonize_script,
  gdalwarp_path            = cfg$tool_paths$gdalwarp_path,
  ogr2ogr_exe              = cfg$tool_paths$ogr2ogr_exe,
  otsu_mode                = cfg$options$legacy_otsu_mode %||% "burnable_only",
  otsu_threshold           = cfg$options$legacy_otsu_threshold %||% 0,
  reference_otsu_threshold = cfg$options$legacy_reference_otsu_threshold %||% 100,
  sample_n                 = cfg$options$legacy_sample_n %||% 2000L,
  reuse_existing           = cfg$options$legacy_reuse_existing %||% TRUE,
  write_unburned           = cfg$options$legacy_write_output %||% TRUE,
  out_root_dir             = legacy_unb_root_dir,
  verbose                  = TRUE
  # Nota: el resto de argumentos (min_pixels, buffers_m, core_thr, alpha_boost,
  # etc.) take internal defaults. To control them all, see
  # version A of the tutorial.
)

otsu_pool    <- to_crs_safe(res_otsu$unburned$legacy_unburned_pool,    crs_master)
otsu_sampled <- to_crs_safe(res_otsu$unburned$legacy_unburned_sampled, crs_master)
if (nrow(otsu_sampled) == 0L && nrow(otsu_pool) > 0L) otsu_sampled <- otsu_pool

otsu_raw <- otsu_sampled |>
  dplyr::mutate(
    source = "otsu_patch_residual",
    neg_type = dplyr::case_when(
      as.character(legacy_decision) == "drop"   ~ "otsu_patch_drop",
      as.character(legacy_decision) == "review" ~ "otsu_patch_review",
      as.character(legacy_decision) == "keep"   ~ "otsu_patch_keep",
      TRUE ~ as.character(neg_type)
    )
  )

cat(sprintf("A4-Parte2 OK: otsu_pool=%d, otsu_sampled=%d\n",
            nrow(otsu_pool), nrow(otsu_sampled)))

# --- A4 Parte 3: combinar fuentes ------------------------------------------
.bind_sf_safe <- function(a, b) {
  df_a <- sf::st_drop_geometry(a)
  df_b <- sf::st_drop_geometry(b)
  geom <- c(sf::st_geometry(a), sf::st_geometry(b))
  combined <- dplyr::bind_rows(df_a, df_b)
  combined[["geometry"]] <- geom
  sf::st_as_sf(combined, sf_column_name = "geometry", crs = sf::st_crs(a))
}
unburned_final_raw <- .bind_sf_safe(det_final_raw, otsu_raw)

# unburned_pool with its final labels
unburned_pool <- unburned_final_raw |>
  sanitize_polygons() |>
  ensure_area_ha() |>
  dplyr::mutate(
    class    = "unburned",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_U_", dplyr::row_number())
  )

unburned_pool <- drop_empty_sf(unburned_pool, tag = "unburned_pool",
                               dump_dir = dirs$`99_LOGS_EMPTY`)

elapsed_a4 <- as.numeric(difftime(Sys.time(), t0_a4, units = "mins"))
cat(sprintf("\nSTEP A4 completo en %.1f min\n", elapsed_a4))
cat(sprintf("unburned_pool total: %d\n", nrow(unburned_pool)))
cat("Composition by source:\n")
print(table(unburned_pool$source))

# =============================================================================
# BLOCK 6 - STEP A5 + A6: train_labeled and writing the pools
# =============================================================================
#
# WHAT IT DOES
#   A5: Concatena burned_pool + unburned_pool en train_labeled (positivos +
#       all negatives, with no caps applied yet).
#   A6: Writes every pool to the GPKG 01_POOLS/<year>_<scenario>_pools.gpkg
#       as separate layers, for later visual inspection.
#
# OUTPUTS
#   train_labeled (memoria)
#   gpkg_pools_out (on disk) with layers: burned_pool, unburned_pool, review_pool,
#                                     scoring_pool, train_labeled,
#                                     unburned_hard, unburned_random,
#                                     unburned_final_raw, exclusion_buffer
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R lineas 1384-1428 (STEP A5, A6)
#
# =============================================================================

cat("\n========== BLOCK 6: STEP A5 + A6 - train_labeled + writing ==========\n")

# A5: train_labeled
train_labeled <- dplyr::bind_rows(burned_pool, unburned_pool) |>
  sanitize_polygons() |>
  ensure_area_ha()

train_labeled <- drop_empty_sf(train_labeled, tag = "train_labeled",
                               dump_dir = dirs$`99_LOGS_EMPTY`)
check_sf(train_labeled, "train_labeled")
stopifnot(!anyDuplicated(train_labeled$fire_uid))

cat(sprintf("train_labeled: burned=%d | unburned=%d\n",
            sum(train_labeled$class == "burned"),
            sum(train_labeled$class == "unburned")))

# A6: escritura
gpkg_pools_out <- file.path(dirs$`01_POOLS`,
                            sprintf("%d_%s_pools.gpkg", YEAR, SCENARIO))

safe_write_gpkg(burned_pool,    gpkg_pools_out, layer = "burned_pool",    tag = "burned_pool_w")
safe_write_gpkg(unburned_pool,  gpkg_pools_out, layer = "unburned_pool",  tag = "unburned_pool_w")
safe_write_gpkg(review_pool,    gpkg_pools_out, layer = "review_pool",    tag = "review_pool_w")
safe_write_gpkg(scoring_pool,   gpkg_pools_out, layer = "scoring_pool",   tag = "scoring_pool_w")
safe_write_gpkg(train_labeled,  gpkg_pools_out, layer = "train_labeled",  tag = "train_labeled_w")
if (nrow(unburned_hard) > 0)
  safe_write_gpkg(unburned_hard, gpkg_pools_out, layer = "unburned_hard", tag = "unburned_hard_w")
if (nrow(unburned_random) > 0)
  safe_write_gpkg(unburned_random, gpkg_pools_out, layer = "unburned_random", tag = "unburned_random_w")
if (nrow(unburned_final_raw) > 0)
  safe_write_gpkg(unburned_final_raw, gpkg_pools_out, layer = "unburned_final_raw", tag = "unburned_final_raw_w")
if (nrow(exclusion_buffer) > 0)
  safe_write_gpkg(exclusion_buffer, gpkg_pools_out, layer = "exclusion_buffer", tag = "exclusion_buffer_w")

cat(sprintf("Pools escritos en: %s\n", gpkg_pools_out))
cat("Layers written:\n")
print(sf::st_layers(gpkg_pools_out)$name)

# =============================================================================
# BLOCK 7 - STEP B2: SPATIAL FOLDS (make_block_folds)
# =============================================================================
#
# WHAT IT DOES
#   Builds spatial folds from contiguous blocks to avoid leakage
#   between train and test. It does 2 repetitions (fold_rep1, fold_rep2) with
#   different seeds, to average out the OOF variance.
#
# DECISIONES METODOLOGICAS
#   block_sizes_m = c(5000, 3000, 2000) - prueba 3 tamanos, elige el mejor
#   k_candidates  = c(5, 4, 3)          - prueba 3 numeros de folds
#   min_burned_units_per_fold = 3       - each fold must hold >=3 burned
#   min_pos_blocks_per_fold   = 10      - each fold must hold >=10 blocks
#                                         with at least 1 burned
#   n_repeats = 2                        - dos repeticiones
#   seed_base = 42                       - reproducibilidad
#
# OUTPUTS (memoria + disco)
#   res_folds$selected            - which block_size was chosen
#   res_folds$saved$train_gpkg    - GPKG with polygons + columns
#                                   block_id, fold_rep1, fold_rep2
#
# REFERENCIA EN EL PAQUETE
#   make_block_folds() en R/internal-sup-make-folds.R
#   internal-sup-orchestrator.R linea 1450 (STEP B2)
#
# =============================================================================

cat("\n========== BLOCK 7: STEP B2 - SPATIAL FOLDS ==========\n")

train_labeled_sf <- safe_read_gpkg(gpkg_pools_out, "train_labeled", "train_labeled")

res_folds <- OtsuFire:::make_block_folds(
  train_labelled_sf         = train_labeled_sf,
  fire_id_col               = "fire_uid",
  split_unit                = "fire",
  block_sizes_m             = c(5000, 3000, 2000),
  k_candidates              = c(5, 4, 3),
  min_burned_units_per_fold = 3,
  min_pos_blocks_per_fold   = 10,
  n_repeats                 = 2,
  seed_base                 = 42,
  lonlat_action             = "transform",
  out_dir                   = dirs$`02_FOLDS`,
  target_year               = YEAR,
  verbose                   = TRUE,
  write_train_with_folds_gpkg = TRUE,
  write_blocks_gpkg           = TRUE,
  write_folds_csv             = TRUE
)

train_with_folds_gpkg <- res_folds$saved$train_gpkg
blocks_gpkg           <- res_folds$saved$blocks_gpkg
bs_sel                <- res_folds$selected$block_size_m

cat(sprintf("\nFolds OK. block_size_m elegido: %d\n", bs_sel))
cat(sprintf("train_with_folds_gpkg: %s\n", train_with_folds_gpkg))

# Inspeccion: ver distribucion de fold_rep1
twf <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "twf")
cat("\nDistribucion fold_rep1 x class:\n")
print(table(twf$fold_rep1, twf$class))
cat("\nDistribucion fold_rep2 x class:\n")
print(table(twf$fold_rep2, twf$class))

# =============================================================================
# BLOCK 8 - STEP B3: EXTRACT FEATURES
# =============================================================================
#
# WHAT IT DOES
#   For each polygon (train + scoring), extracts the 51 canonical features
#   of the whitelist, by family:
#     RBR summer:   rbr_p10, rbr_med, rbr_p90, rbr_iqr (4)
#     RBR autumn:   rbr_aw_p10, rbr_aw_med, rbr_aw_p90, rbr_aw_iqr,
#                   rbr_aw_valid_frac (5)
#     DOY:          doy_p10, doy_med, doy_p90, doy_iqr (4)
#     Topografia:   elev_*, slope_* (12)
#     Corine:       cor_agri_frac, cor_forest_frac, ... (7)
#     Persistencia: persist_delta, persist_ratio (2)
#     Hotspots:     hs_in_poly, hs_in_buffer, hs_used_n, hs_min_dist_m,
#                   hs_frp_sum, hs_frp_max, hs_conf_mean, hs_hiConf_n,
#                   hs_support_present, hs_no_support_when_available,
#                   hs_only_buffer_support, hotspot_available (12)
#                   [0.6.2: hs_any removed from the canonical whitelist]
#     Validez:      otras (4)
#
# DECISIONES METODOLOGICAS
#   use_hotspots = TRUE en 2005 (post-MODIS)
#   hs_buffer_m = 1000   - radius for nearby hotspots
#   hs_start_month/end_month = 6,10  - estacion fuego
#   hi_conf_thr = 0.8    - threshold for hs_hiConf_n
#   cor_groups = Corine groups, already in cfg$options$cor_groups
#
# OUTPUTS (disco)
#   features_geometry.gpkg with two layers:
#     train_features    - training polygons with all the features
#     scoring_features  - polygons of the scoring universe, with features
#
# REFERENCIA EN EL PAQUETE
#   extract_features() en R/internal-sup-extract-features.R
#   internal-sup-orchestrator.R lineas 1681-1741 (STEP B3.3)
#
# =============================================================================

cat("\n========== BLOCK 8: STEP B3 - EXTRACT FEATURES ==========\n")
cat("(this stage takes 5-15 min depending on the size of the area)\n")

# Re-read from disk to make sure extract_features sees the data exactly
# as the real orchestrator will see it
train_folds        <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "tf")
scoring_candidates <- safe_read_gpkg(gpkg_pools_out, "scoring_pool", "sp")

# Reproject hotspots if the CRS does not match
hotspots_sf <- hotspots_sf_base
if (nrow(hotspots_sf) > 0) {
  if (sf::st_crs(hotspots_sf) != sf::st_crs(train_folds)) {
    hotspots_sf <- sf::st_transform(hotspots_sf, sf::st_crs(train_folds))
  }
}
use_hotspots_flag <- nrow(hotspots_sf) > 0
cat(sprintf("use_hotspots_flag: %s (n=%d)\n", use_hotspots_flag, nrow(hotspots_sf)))

# Corine groups - taken from cfg
cor_groups <- cfg$options$cor_groups

t0_b3 <- Sys.time()

feat <- OtsuFire:::extract_features(
  train_folds = train_folds,
  unlabeled   = scoring_candidates,

  id_col    = "fire_uid",
  class_col = "class",
  pos_lab   = "burned",
  neg_lab   = "unburned",

  build_features = TRUE,

  rbr_summer = rbr_summer,
  doy_post   = doy_post,
  rbr_aw     = rbr_aw,
  nbr_pre    = NULL,
  nbr_post   = NULL,
  dnbr       = NULL,
  dem        = dem_r,
  slope      = slope_r,
  corine_r   = corine_r,

  hotspots   = hotspots_sf,
  frp_col    = "frp",
  conf_col   = "confidence",
  year_col   = "year",
  month_col  = "month",
  hs_start_month = 6,
  hs_end_month   = 10,
  hi_conf_thr    = 0.8,

  year_target           = YEAR,
  hs_year_min_available = 2000,
  hs_buffer_m           = 1000,
  hs_missing_value      = -9999,
  hs_use_season_filter  = TRUE,

  cor_groups = cor_groups,

  use_doy        = TRUE,
  use_aw         = TRUE,
  use_nbr        = FALSE,
  use_hotspots   = use_hotspots_flag,

  max_cells_in_memory       = NULL,
  return_features           = TRUE,
  return_features_geometry  = TRUE,
  save_features_dir         = dirs$`03_FEATURES`,
  save_features_format      = "rds",
  save_features_gpkg        = TRUE,
  train_features_basename   = "train_features",
  scoring_features_basename = "scoring_features",
  train_features_layer      = "train_features",
  scoring_features_layer    = "scoring_features",
  verbose                   = TRUE
)

elapsed_b3 <- as.numeric(difftime(Sys.time(), t0_b3, units = "mins"))
cat(sprintf("\nSTEP B3 completo en %.1f min\n", elapsed_b3))

features_gpkg <- file.path(dirs$`03_FEATURES`, "features_geometry.gpkg")
cat("Layers generated:\n")
print(sf::st_layers(features_gpkg)$name)

# Inspeccion de columnas
train_features_check <- safe_read_gpkg(features_gpkg, "train_features", "tfc")
feature_cols <- setdiff(names(train_features_check),
                        c("fire_uid", "class", "source", "poly_id",
                          "block_id", "fold_rep1", "fold_rep2", "geometry"))
cat(sprintf("\nFeatures extraidas: %d\n", length(feature_cols)))
cat("Familias detectadas:\n")
fam_counts <- list(
  RBR_SW    = sum(grepl("^rbr_(p10|med|p90|iqr|sd|mean)$", feature_cols)),
  RBR_AW    = sum(grepl("^rbr_aw_", feature_cols)),
  DOY       = sum(grepl("^doy_", feature_cols)),
  TOPO      = sum(grepl("^(elev|slope)_", feature_cols)),
  CORINE    = sum(grepl("^cor_", feature_cols)),
  PERSIST   = sum(grepl("^persist_", feature_cols)),
  HOTSPOTS  = sum(grepl("^(hs_|hotspot_)", feature_cols))
)
print(fam_counts)
cat(sprintf("Total checksum: %d\n", sum(unlist(fam_counts))))

# =============================================================================
# BLOCK 9 - STEP C0 + C1: READ FEATURES + XGBOOST PARAMS
# =============================================================================
#
# WHAT IT DOES
#   C0: reads train_features and scoring_features from the GPKG, joining the
#       folds if needed.
#   C1: builds the XGBoost hyperparameters. It computes scale_pos_weight
#       as n_neg/n_pos to balance the classes in the loss.
#
# DECISIONES METODOLOGICAS (XGBOOST)
#   booster          = "gbtree"
#   objective        = "binary:logistic"
#   eval_metric      = c("logloss", "aucpr")  - logloss first because it is
#                                              quien drivea early stopping
#   eta              = 0.06   (learning rate moderado)
#   max_depth        = 5      (shallow trees, for regularisation)
#   subsample        = 0.85   (row sub-sampling per tree)
#   colsample_bytree = 0.75   (column sub-sampling per tree)
#   min_child_weight = 5      (regularizacion)
#   gamma            = 0      (splits with no extra penalty)
#   scale_pos_weight = n_neg/n_pos  (balanceo automatico de clases)
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R lineas 1751-1832 (STEP C0, C1)
#
# =============================================================================

cat("\n========== BLOCK 9: STEP C0 + C1 - FEATURES + XGB PARAMS ==========\n")

train_features   <- safe_read_gpkg(features_gpkg, "train_features",   "tf2")
scoring_features <- safe_read_gpkg(features_gpkg, "scoring_features", "sf2")

# Make sure train_features carries the folds
needed_folds <- c("block_id", "fold_rep1", "fold_rep2")
if (!all(needed_folds %in% names(train_features))) {
  cat("Joining folds into train_features...\n")
  twf2 <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "twf2")
  folds_df <- sf::st_drop_geometry(twf2) |>
    dplyr::select(fire_uid, dplyr::any_of(needed_folds))
  train_features <- train_features |> dplyr::left_join(folds_df, by = "fire_uid")
}

# Drop geometry for modelling
labelled    <- sf::st_drop_geometry(train_features)
burned_like <- sf::st_drop_geometry(scoring_features)
labelled_df <- labelled

# Construir params XGBoost
n_pos <- sum(labelled$class == "burned",   na.rm = TRUE)
n_neg <- sum(labelled$class == "unburned", na.rm = TRUE)
spw   <- if (n_pos > 0) n_neg / n_pos else 1

params <- list(
  booster          = "gbtree",
  objective        = "binary:logistic",
  eval_metric      = c("logloss", "aucpr"),
  eta              = 0.06,
  max_depth        = 5,
  subsample        = 0.85,
  colsample_bytree = 0.75,
  min_child_weight = 5,
  gamma            = 0,
  scale_pos_weight = spw
)

cat(sprintf("n_pos=%d, n_neg=%d, scale_pos_weight=%.3f\n", n_pos, n_neg, spw))
cat("XGB params:\n")
print(params)

# =============================================================================
# BLOCK 10 - STEP C2: OOF PIPELINE (per-fold training + diagnostics)
# =============================================================================
#
# WHAT IT DOES
#   For each repetition (fold_rep1, fold_rep2) and each fold k=1..5:
#     1. Marks fold k as test and the rest as train.
#     2. Aplica caps de sampling:
#          contextual_exclusion (deterministic_drop_hard, neg_type != spectral) cap 0.25 x n_burned
#          spectral_hard_negative (deterministic_drop_hard con neg_type == "spectral_reject_medium") cap 1.0
#          random_burnable_background cap 1.0
#          otsu_patch_residual cap 1.0
#     3. Entrena XGBoost en train, predice en test.
#     4. Acumula predicciones.
#   In the end each row has n_preds=2 predictions (one per repetition).
#   Calcula metricas (accuracy, precision, recall, etc.) en multiples
#   thresholds, identifica el threshold optimo segun varios criterios.
#
# DECISIONES METODOLOGICAS
#   feature_whitelist_override = NULL  (canon 51)
#   feature_weights            = NULL  (uniforme 1.0)
#   median_from = "labelled"   - NA imputation by each feature's median
#                                en el set labelled
#
# OUTPUTS (disco)
#   05_OOF/<prefix>_oof_metrics_summary.txt
#   05_OOF/<prefix>_labeled_oof_summary.gpkg     - each polygon with p_oof_mean
#   05_OOF/<prefix>_oof_agg.csv
#   05_OOF/<prefix>_oof_long.csv
#   05_OOF/<prefix>_oof_metrics_by_threshold.csv
#   05_OOF/<prefix>_oof_best_thresholds.csv
#   04_MATRIX/<prefix>_dm.rds  (DesignMatrix)
#
# REFERENCIA EN EL PAQUETE
#   run_dm_oof_pipeline() en R/internal-sup-oof-wrapper.R
#   internal-sup-orchestrator.R lineas 1834-1870 (STEP C2)
#
# IMPORTANTE
#   This block takes 10-20 min for 2005 (the slowest modelling stage).
#
# =============================================================================

cat("\n========== BLOCK 10: STEP C2 - OOF PIPELINE ==========\n")
cat("(this stage takes 10-20 min)\n")

prefix_oof <- sprintf("%d_%s_patch", YEAR, SCENARIO)
prefix     <- sprintf("%d_%s_patch_certified", YEAR, SCENARIO)

t0_c2 <- Sys.time()

pipe1 <- OtsuFire:::run_dm_oof_pipeline(
  labelled    = labelled,
  burned_like = burned_like,
  labelled_df = labelled_df,
  params      = params,

  result_dir  = TUTORIAL_RESULT_DIR,
  target_year = YEAR,
  prefix      = prefix_oof,

  id_cols     = c("fire_uid", "class", "source", "poly_id",
                  "block_id", "fold_rep1", "fold_rep2"),
  cat_cols    = character(0),
  hs_n_col    = "hs_used_n",
  hs_conf_col = "hs_conf_mean",
  hs_frp_col  = "hs_frp_max",
  median_from = "labelled",
  save_dir_dm = dirs$`04_MATRIX`,
  save_prefix = paste0(YEAR, "_", SCENARIO, "_patch"),
  overwrite   = TRUE,

  # ----- BASELINE: NULL en ambos -----
  feature_whitelist_override = NULL,
  feature_weights            = NULL,

  labelled_gpkg  = train_with_folds_gpkg,
  labelled_layer = "train_with_folds"
)

elapsed_c2 <- as.numeric(difftime(Sys.time(), t0_c2, units = "mins"))
cat(sprintf("\nSTEP C2 completo en %.1f min\n", elapsed_c2))

cat("Files generated by the OOF pipeline:\n")
print(pipe1$files)

# Inspecting the summary
oof_summary_path <- file.path(dirs$`05_OOF`,
                              paste0(prefix_oof, "_oof_metrics_summary.txt"))
cat("\n--- oof_metrics_summary.txt ---\n")
cat(readLines(oof_summary_path), sep = "\n")

# =============================================================================
# BLOCK 11 - STEP C3: FINAL MODEL + SCORING + FINAL MAP
# =============================================================================
#
# WHAT IT DOES
#   1. Trains the FINAL model on ALL the labelled data (no OOF).
#   2. Applies the model to the scoring_features universe to score it.
#   3. Builds the final map, with the class assigned by the optimal threshold
#      inherited from the OOF stage.
#   4. Identifies burned_like (unlabelled polygons that the model
#      considera quemados).
#
# DECISIONES METODOLOGICAS
#   contextual_exclusion_to_burned_ratio = 0.25
#   spectral_hard_negative_to_burned_ratio = 1.0
#   random_to_burned_ratio = 1.0
#   otsu_unburned_to_burned_ratio = 1.0
#   feature_whitelist_override = NULL
#   feature_weights = NULL
#   preyear_overlap_threshold = cfg$options$currentyear_preyear_overlap_thr
#   hotspot_density_threshold = cfg$options$currentyear_hotspot_density_thr
#   temporal_penalty_floor    = cfg$options$currentyear_temporal_penalty_floor
#
# OUTPUTS (disco)
#   07_FINAL_MODEL_V2/<prefix>_certified_final_model.rds
#   07_FINAL_MODEL_V2/<prefix>_certified_recipe.rds
#   07_FINAL_MODEL_V2/<prefix>_certified_meta.txt
#   07_FINAL_MODEL_V2/<prefix>_certified_feature_importance.csv
#   07_FINAL_MODEL_V2/<prefix>_certified_model_summary.txt
#   07_FINAL_MODEL_V2/<prefix>_certified_training_ok.csv
#   07_FINAL_MODEL_V2/<prefix>_certified_training_ok.gpkg
#   08_SCORED/<prefix>_certified_scored.csv
#   08_SCORED/<prefix>_certified_scored.gpkg
#   09_FINAL_MAP/<prefix>_certified_final_map.gpkg
#   09_FINAL_MAP/<prefix>_certified_burned_like_scored.gpkg
#
# REFERENCIA EN EL PAQUETE
#   run_train_final_model_and_export_final_map() en
#     R/internal-sup-train-final-wrapper.R
#   Funcion interna real: train_final_model_direct() en
#     R/internal-sup-train-final-direct.R
#   internal-sup-orchestrator.R lineas 1875-1925 (STEP C3)
#
# =============================================================================

cat("\n========== BLOCK 11: STEP C3 - FINAL MODEL + SCORING ==========\n")
cat("(this stage takes 5-15 min)\n")

oof_agg_csv      <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_oof_agg.csv"))
oof_summary_gpkg <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_labeled_oof_summary.gpkg"))
stopifnot(file.exists(oof_agg_csv), file.exists(oof_summary_gpkg))

t0_c3 <- Sys.time()

pipe2 <- OtsuFire:::run_train_final_model_and_export_final_map(
  result_dir     = TUTORIAL_RESULT_DIR,
  qa             = oof_agg_csv,
  labelled_gpkg  = features_gpkg,
  labelled_layer = "train_features",
  out_dir        = dirs$`07_FINAL_MODEL_V2`,
  prefix         = prefix,
  overwrite      = TRUE,
  verbose        = TRUE,

  qa_labelled_gpkg  = oof_summary_gpkg,
  qa_labelled_layer = "labeled_oof_summary",

  labelled_features_gpkg  = features_gpkg,
  labelled_features_layer = "train_features",

  model_rds  = file.path(dirs$`07_FINAL_MODEL_V2`,
                         paste0(prefix, "_final_model.rds")),
  recipe_rds = file.path(dirs$`07_FINAL_MODEL_V2`,
                         paste0(prefix, "_recipe.rds")),

  unlabeled_gpkg  = features_gpkg,
  unlabeled_layer = "scoring_features",

  out_score_dir = dirs$`08_SCORED`,
  out_map_dir   = dirs$`09_FINAL_MAP`,

  export_burned_like        = TRUE,
  preyear_overlap_threshold = cfg$options$currentyear_preyear_overlap_thr %||% 0.5,
  hotspot_density_threshold = cfg$options$currentyear_hotspot_density_thr %||% 1.0,
  temporal_penalty_floor    = cfg$options$currentyear_temporal_penalty_floor %||% 0.1,

  # ----- BASELINE: canonical caps + no overrides + no weights -----
  contextual_exclusion_to_burned_ratio   = 0.25,
  spectral_hard_negative_to_burned_ratio = 1.0,
  random_to_burned_ratio                 = 1.0,
  otsu_unburned_to_burned_ratio          = 1.0,
  feature_whitelist_override             = NULL,
  feature_weights                        = NULL
)

elapsed_c3 <- as.numeric(difftime(Sys.time(), t0_c3, units = "mins"))
cat(sprintf("\nSTEP C3 completo en %.1f min\n", elapsed_c3))

# Inspecting the final model's meta
meta_path <- file.path(dirs$`07_FINAL_MODEL_V2`,
                       paste0(prefix, "_meta.txt"))
if (file.exists(meta_path)) {
  cat("\n--- meta.txt of the final model ---\n")
  cat(readLines(meta_path), sep = "\n")
}

# =============================================================================
# BLOCK 12 - FINAL VERIFICATION: canonical Baseline
# =============================================================================
#
# WHAT IT DOES
#   Checks that meta.txt confirms the canonical Baseline configuration:
#     feature_whitelist_override_applied = FALSE
#     feature_weights_applied            = FALSE
#   And that the metrics make sense (drops low probability, keeps high).
#
# =============================================================================

cat("\n========== BLOCK 12: BASELINE VERIFICATION ==========\n")

if (file.exists(meta_path)) {
  meta_lines <- readLines(meta_path)
  override_line <- grep("feature_whitelist_override_applied", meta_lines, value = TRUE)
  weights_line  <- grep("feature_weights_applied", meta_lines, value = TRUE)
  n_x_line      <- grep("^n_x_cols:", meta_lines, value = TRUE)
  n_feat_line   <- grep("^n_feature_cols:", meta_lines, value = TRUE)

  cat("Relevant lines of meta.txt:\n")
  cat(" ", override_line, "\n")
  cat(" ", weights_line, "\n")
  cat(" ", n_x_line, "\n")
  cat(" ", n_feat_line, "\n")

  is_baseline <- grepl("FALSE", override_line) && grepl("FALSE", weights_line)
  if (is_baseline) {
    cat("\nBaseline canonical CONFIRMADO.\n")
  } else {
    cat("\nWARNING: meta.txt does NOT reflect the canonical Baseline.\n")
  }
}

# Verifying the final map
final_map_path <- file.path(dirs$`09_FINAL_MAP`,
                            paste0(prefix, "_final_map.gpkg"))
if (file.exists(final_map_path)) {
  fm <- sf::read_sf(final_map_path, layer = "final_map_full")

  drops <- fm$p_burned_model[fm$class_final == "drop"]
  keeps <- fm$p_burned_model[fm$class_final == "keep"]

  drops_med <- median(drops, na.rm = TRUE)
  keeps_med <- median(keeps, na.rm = TRUE)

  cat(sprintf("\nDrops mediana: %.4f (esperado <0.1)\n", drops_med))
  cat(sprintf("Keeps mediana: %.4f (esperado >0.9)\n", keeps_med))

  if (drops_med < 0.1 && keeps_med > 0.9) {
    cat("\nBASELINE 0.5.0 PASSED.\n")
  } else {
    cat("\nBaseline does NOT meet the expected thresholds.\n")
  }
}

cat("\n========== TUTORIAL COMPLETADO ==========\n")
cat(sprintf("Outputs en: %s\n", TUTORIAL_RESULT_DIR))
