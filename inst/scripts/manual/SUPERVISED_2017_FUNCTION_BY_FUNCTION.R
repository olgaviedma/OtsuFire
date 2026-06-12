# =============================================================================
# OtsuFire — SUPERVISED 2017, FUNCTION BY FUNCTION (manual interactive script)
# CANONICAL VERSIONED TEMPLATE (inst/scripts/manual/).
#
# PURPOSE: run the supervised burned-area phase ONE FUNCTION AT A TIME, with a
# readable CHECK after every block, so the whole pipeline can be inspected and
# debugged interactively. Each of the 13 blocks is INDEPENDENT: run a block,
# read its CHECK, STOP. A block never calls the next one and never overwrites a
# prior block's result object.
#
# IMPORTANT — FROZEN CODE CONTRACT:
#   * Uses the INSTALLED OtsuFire from a configurable ISOLATED library.
#   * NEVER pkgload::load_all(). NEVER library(OtsuFire) without lib.loc.
#   * Asserts find.package() resolves INSIDE the isolated lib (so an older
#     install elsewhere cannot shadow it).
#
# MODEL_PROFILE switch (00_SETUP): "FULL" or "NO_HOTSPOT". The ONLY difference is
# explicit feature configuration:
#   FULL       = canonical 50-base-feature whitelist, use_hotspots = TRUE.
#   NO_HOTSPOT = feature_whitelist_override = the 38 non-hotspot base features
#                (50 minus the 12 hotspot-lineage features), use_hotspots = FALSE.
#
# EFFIS_MODE switch (00_SETUP), block 13:
#   "STANDARD_90M"        = RECOMMENDED default. EPSG:3035, 90 m, terra todisk,
#                           modest memfrac, dedicated tempdir, memory logging,
#                           safe tempdir cleanup. Low-memory.
#   "NATIVE_HIGH_MEMORY"  = high-consumption; guarded behind ALLOW_HIGH_MEMORY.
#   Both modes call the SAME validate_fire_maps() with the SAME metric
#   definitions; only resolution / memory strategy differ.
#
# NOTE: p_burned is a MODEL SCORE, not a calibrated probability. Use it for
# ranking / thresholding; pick the operating threshold from EFFIS validation.
# =============================================================================


# #############################################################################
# 00_SETUP — libpaths / version pin / find.package check / profile / paths /
#            EFFIS mode / output folder. Run this FIRST; every other block
#            assumes the objects it defines are in the session.
# #############################################################################

## ---- (a) ISOLATED LIB + VERSION PIN (REQUIRED) ----------------------------
# SNAPSHOT_LIB: an isolated library containing ONLY the frozen OtsuFire to use.
# Fill this once the 0.6.2 isolated lib exists. TODO(Natalia): set this path.
SNAPSHOT_LIB <- "<PATH_TO_0.6.2_ISOLATED_LIB>"   # e.g. C:/.../00_SNAPSHOTS/<snap>/Rlib

.libPaths(c(SNAPSHOT_LIB, .libPaths()))
library(OtsuFire, lib.loc = SNAPSHOT_LIB)

# Version pin: require >= 0.6.2 once the 0.6.2 snapshot exists. Until then 0.6.1
# is tolerated with a printed TODO (hs_any already removed -> whitelist 50).
.MIN_VERSION <- "0.6.2"
.pkg_ver <- packageVersion("OtsuFire")
if (.pkg_ver < .MIN_VERSION) {
  if (.pkg_ver >= "0.6.1") {
    message(sprintf(
      "[00_SETUP][TODO] OtsuFire %s < %s. Tolerated for now (>=0.6.1). Re-point SNAPSHOT_LIB at the 0.6.2 isolated lib when available, then this becomes a hard requirement.",
      as.character(.pkg_ver), .MIN_VERSION))
  } else {
    stop(sprintf("[00_SETUP] OtsuFire %s < 0.6.1 — too old. Point SNAPSHOT_LIB at the frozen 0.6.x isolated lib.",
                 as.character(.pkg_ver)))
  }
}

# Assert find.package resolves INSIDE the isolated lib (fails if an old install
# would shadow). This is the load-bearing guard against picking up a stale copy.
.pkg_path <- normalizePath(find.package("OtsuFire"), winslash = "/")
.iso_path <- normalizePath(SNAPSHOT_LIB, winslash = "/")
if (!startsWith(.pkg_path, .iso_path)) {
  stop(sprintf(paste0(
    "[00_SETUP] OtsuFire resolved OUTSIDE the isolated lib:\n  find.package = %s\n  SNAPSHOT_LIB = %s\n",
    "An old install would shadow the frozen code. Fix .libPaths / SNAPSHOT_LIB."),
    .pkg_path, .iso_path))
}

suppressPackageStartupMessages({ library(sf); library(terra) })
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- (b) MODEL PROFILE SWITCH ---------------------------------------------
# The ONLY methodological difference between the two profiles.
MODEL_PROFILE <- "FULL"          # "FULL" or "NO_HOTSPOT"
stopifnot(MODEL_PROFILE %in% c("FULL", "NO_HOTSPOT"))

# The 12 hotspot-lineage base features (hs_any is GONE from the whitelist;
# do NOT reference it). NO_HOTSPOT removes exactly these 12.
HOTSPOT_LINEAGE_BASE <- c(
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support")
stopifnot(length(HOTSPOT_LINEAGE_BASE) == 12L)

# Canonical 50 base features (FULL). NO_HOTSPOT override = these 50 minus the 12.
FULL_BASE_50 <- c(
  "rbr_valid_frac", "rbr_p10", "rbr_med", "rbr_p90", "rbr_iqr",
  "cor_open_frac", "cor_wetlands_frac", "cor_water_frac", "cor_herbaceous_frac",
  "cor_urban_frac", "cor_agri_frac", "cor_forest_frac",
  "elev_valid_frac", "elev_p10", "elev_med", "elev_p90", "elev_iqr", "elev_mean", "elev_sd",
  "slope_valid_frac", "slope_p10", "slope_med", "slope_p90", "slope_iqr", "slope_mean", "slope_sd",
  "doy_valid_frac", "doy_p10", "doy_med", "doy_p90", "doy_iqr",
  "rbr_aw_valid_frac", "rbr_aw_p10", "rbr_aw_med", "rbr_aw_p90", "rbr_aw_iqr",
  "persist_delta", "persist_ratio",
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support")
stopifnot(length(FULL_BASE_50) == 50L, all(HOTSPOT_LINEAGE_BASE %in% FULL_BASE_50))
NOHS_BASE_38 <- setdiff(FULL_BASE_50, HOTSPOT_LINEAGE_BASE)
stopifnot(length(NOHS_BASE_38) == 38L)

# Profile-resolved knobs.
FEATURE_WHITELIST_OVERRIDE <- if (MODEL_PROFILE == "NO_HOTSPOT") NOHS_BASE_38 else NULL
USE_HOTSPOTS               <- if (MODEL_PROFILE == "NO_HOTSPOT") FALSE else TRUE
EXPECTED_BASE_N            <- if (MODEL_PROFILE == "NO_HOTSPOT") 38L else 50L

## ---- (c) PATHS (REQUIRED) -------------------------------------------------
# TODO(Natalia): fill every placeholder. These mirror the long-run runners.
base_dir       <- "<PATH_TO_FIRE_MAPPING_BASE>"            # e.g. C:/.../00_FIRE_MAPPING
data_base      <- file.path(base_dir, "1_DATA")
results_base   <- file.path(data_base, "Results")
composite_base <- file.path(data_base, "Imagery", "Composites_90m")

target_year <- 2017L
scenario    <- "balanced"
result_name <- "Min_Min"

internal_decisions <- file.path(results_base, as.character(target_year), "DETERMINISTIC",
                                sprintf("%d_%s", target_year, scenario),
                                "05_DECISIONS", "internal_decisions.gpkg")
change_index       <- file.path(composite_base, result_name,
                                sprintf("MinMin_%d_mosaic_res90m.tif", target_year))
delayed_change_index <- file.path(composite_base, "Autumn",
                                  sprintf("mean_mean_%d_mosaic.tif", target_year))
hotspots           <- file.path(data_base, "Hotspots",
                                sprintf("hotspots_iberia_%d.geojson", target_year))
reference_burned_map <- file.path(data_base, "Fires/Validation_fires_burneable_verano",
                                  sprintf("Effis_CA_%d_maskKeep_summer.shp", target_year))
burnable_mask_run  <- file.path(data_base, "Corine_Masks",
                                sprintf("burneable_mask_binary_corine12_3035_alignedMinMin%d.tif",
                                        target_year))
validation_mask_shapefile <- file.path(data_base, "Mask_StudyArea", "mask_Peninsula_3035.shp")
strata_raster <- file.path(data_base, "Corine_Masks/STRATA", "strata_CLC_2012_res30.tif")
strata_lut    <- file.path(data_base, "Corine_Masks/LUT/lut_full_strata8_v1.csv")

# External tools for the all_sources Otsu residual negative builder. TODO(Natalia): fill.
python_exe             <- "<PATH_TO_python.exe>"
gdal_polygonize_script <- "<PATH_TO_gdal_polygonize.py>"
gdalwarp_path          <- "<PATH_TO_gdalwarp.exe>"
ogr2ogr_exe            <- "<PATH_TO_ogr2ogr.exe>"

## ---- (d) CANONICAL CONTROLS -----------------------------------------------
# GATE 6.7 (2026-06-12): the two operative negative-bucket caps (random + otsu).
CAP_RANDOM <- 1.0; CAP_OTSU <- 1.0
SEED_BASE      <- 42L
FOLD_COLS      <- c("fold_rep1", "fold_rep2")
OPERATING_THRESHOLD <- 0.50
NROUNDS_MAX    <- 4000L
EARLY_STOP     <- 80L
XGB_NTHREAD    <- 6L
Sys.setenv(OMP_NUM_THREADS = as.character(XGB_NTHREAD))

## ---- (e) EFFIS MODE -------------------------------------------------------
EFFIS_MODE <- "STANDARD_90M"     # RECOMMENDED default. "STANDARD_90M" or "NATIVE_HIGH_MEMORY".
stopifnot(EFFIS_MODE %in% c("STANDARD_90M", "NATIVE_HIGH_MEMORY"))
# NATIVE_HIGH_MEMORY must NOT run accidentally: it requires this explicit opt-in.
ALLOW_HIGH_MEMORY <- FALSE
EFFIS_MEMFRAC     <- 0.4         # modest terra memory fraction for STANDARD_90M

## ---- (f) EXCLUSIVE MANUAL OUTPUT FOLDER -----------------------------------
# Every block writes ONLY under here (never the read-only data tree, never the
# long-run tree). TODO(Natalia): set this to a fresh, manual-only folder.
MANUAL_OUT <- "<PATH_TO_MANUAL_OUTPUT_FOLDER>"   # e.g. C:/.../00_MANUAL/SUPERVISED_2017_<STAMP>
dir.create(MANUAL_OUT, recursive = TRUE, showWarnings = FALSE)
MANUAL_ENGINE <- file.path(MANUAL_OUT, "ENGINE_ROUTES")
dir.create(MANUAL_ENGINE, recursive = TRUE, showWarnings = FALSE)

## ---- (g) PRINT THE ENVIRONMENT --------------------------------------------
cat("\n================= 00_SETUP CHECK =================\n")
cat(sprintf("OtsuFire version : %s\n", as.character(.pkg_ver)))
cat(sprintf("find.package     : %s\n", .pkg_path))
cat(sprintf("isolated lib OK  : %s\n", startsWith(.pkg_path, .iso_path)))
cat(sprintf("MODEL_PROFILE    : %s  (expected base features = %d)\n", MODEL_PROFILE, EXPECTED_BASE_N))
cat(sprintf("use_hotspots     : %s\n", USE_HOTSPOTS))
cat(sprintf("whitelist ovrd n : %s\n", if (is.null(FEATURE_WHITELIST_OVERRIDE)) "NULL (FULL 50)" else length(FEATURE_WHITELIST_OVERRIDE)))
cat(sprintf("caps             : rnd=%.2f otsu=%.2f\n",
            CAP_RANDOM, CAP_OTSU))
cat(sprintf("seed / folds     : %d / %s\n", SEED_BASE, paste(FOLD_COLS, collapse = "+")))
cat(sprintf("EFFIS_MODE       : %s  (ALLOW_HIGH_MEMORY=%s)\n", EFFIS_MODE, ALLOW_HIGH_MEMORY))
cat(sprintf("MANUAL_OUT       : %s\n", MANUAL_OUT))
cat(sprintf("target year      : %d  scenario=%s  run_name=%s\n", target_year, scenario, result_name))
cat("00_SETUP done. STOP. Next: run 01_BUILD_CFG.\n")
cat("=================================================\n")


# #############################################################################
# 01_BUILD_CFG — build_supervised_burned_config (single source of truth for all
#                methodological params). Keeps B01_cfg. Prints resolved caps,
#                provenance, and feature override count.
# #############################################################################
B01_cfg <- build_supervised_burned_config(
  scenario             = scenario,
  internal_decisions   = internal_decisions,
  change_index         = change_index,
  hotspots             = hotspots,
  reference_burned_map = reference_burned_map,
  target_year          = target_year,
  output_dir           = MANUAL_ENGINE,
  run_name             = result_name,
  burnable_mask        = burnable_mask_run,
  nrounds_max          = NROUNDS_MAX,
  early_stop           = EARLY_STOP,
  oof_seed_base = SEED_BASE, final_sampling_seed = SEED_BASE, final_seed = SEED_BASE,
  # GATE 6.7 (2026-06-12): typed PUBLIC negative_pool_params (single caps source)
  # + technical runtime_options. OtsuFire always uses inner-early-stopping
  # selection + full-data refit; OOF uses the SAME capped negative-sampling
  # policy as the FINAL model.
  negative_pool_params = list(
    random = list(n_cells = 1500L, rbr_quantile = 0.50),
    otsu   = list(candidate_threshold = 0, reference_threshold = 100),
    caps   = c(random = CAP_RANDOM, otsu = CAP_OTSU)
  ),
  runtime_options = list(reuse_existing = TRUE, write_outputs = TRUE,
                         verbose = TRUE),
  feature_whitelist_override = FEATURE_WHITELIST_OVERRIDE,
  options = list(
    data_base = data_base, composite_base = composite_base, result_name = result_name))
B01_cfg$tool_paths$python_exe             <- python_exe
B01_cfg$tool_paths$gdal_polygonize_script <- gdal_polygonize_script
B01_cfg$tool_paths$gdalwarp_path          <- gdalwarp_path
B01_cfg$tool_paths$ogr2ogr_exe            <- ogr2ogr_exe

cat("\n================= 01_BUILD_CFG CHECK =================\n")
.caps <- B01_cfg$train_control$caps
cat("resolved caps:\n"); print(unlist(.caps))
cat("training          : early-stopping selection + full-data refit (fixed)\n")
cat(sprintf("oof_sampling (fixed): %s  (same capped policy as FINAL; not a knob)\n",
            B01_cfg$train_control$oof_sampling %||% "capped"))
.ovr <- B01_cfg$train_control$feature_whitelist_override
cat(sprintf("feature_whitelist_override n : %s\n", if (is.null(.ovr)) "NULL (FULL 50)" else length(.ovr)))
.prov <- B01_cfg$resolved_params_provenance
if (!is.null(.prov)) {
  cat(sprintf("provenance(caps)      : %s\n", .prov$train_control$caps %||% "?"))
  cat(sprintf("provenance(whitelist) : %s\n", .prov$train_control$feature_whitelist_override %||% "(default)"))
}
if (MODEL_PROFILE == "NO_HOTSPOT")
  stopifnot(identical(sort(.ovr), sort(NOHS_BASE_38)))
cat("01_BUILD_CFG done. STOP. Next: run 02_VALIDATE_INPUTS.\n")
cat("=====================================================\n")


# #############################################################################
# 02_VALIDATE_INPUTS — validate_supervised_execution(cfg, strict = TRUE).
#                      Pre-run structural / schema / inputs / year check.
#                      Keeps B02_validation. Prints the structured report.
# #############################################################################
B02_validation <- validate_supervised_execution(
  config = B01_cfg, strict = TRUE,
  feature_whitelist_override = FEATURE_WHITELIST_OVERRIDE,
  data_base = data_base, composite_base = composite_base, result_name = result_name)

cat("\n================= 02_VALIDATE_INPUTS CHECK =================\n")
cat("structured report:\n")
print(B02_validation)
cat("\nExpected highlights to confirm:\n")
cat("  * sanitize_decisions: PASS_WITH_SANITIZATION 1278 -> 1274 (4 empty geoms dropped)\n")
cat("  * year checks: target_year present + consistent\n")
cat("  * feature schema check: PASS (recipe schema covered by inputs)\n")
cat("02_VALIDATE_INPUTS done. STOP. Next: run 03_SANITIZE_GEOMETRIES.\n")
cat("===========================================================\n")


# #############################################################################
# 03_SANITIZE_GEOMETRIES — inspect the geometry sanitization audit.
#   The sanitize step is exercised inside the pools build; here we surface the
#   audit counts explicitly. Keeps B03_sanitized (the sanitize audit).
#   Expected: input = 1278, empty = 4, usable = 1274.
# #############################################################################
# The pools builder runs sanitize internally and writes a QA audit; we read the
# decisions layer and the sanitize audit to confirm the counts BEFORE the heavy
# pools build (this block is read-only on the inputs).
.dec <- sf::st_read(internal_decisions, quiet = TRUE)
.n_in    <- nrow(.dec)
.empty   <- sum(sf::st_is_empty(.dec))
.usable  <- .n_in - .empty
B03_sanitized <- list(input = .n_in, empty = .empty, usable = .usable,
                      decisions_path = internal_decisions)

cat("\n================= 03_SANITIZE_GEOMETRIES CHECK =================\n")
cat(sprintf("decisions input rows : %d\n", B03_sanitized$input))
cat(sprintf("empty geometries     : %d\n", B03_sanitized$empty))
cat(sprintf("usable rows          : %d\n", B03_sanitized$usable))
cat("Expected (2017 balanced): input=1278, empty=4, usable=1274.\n")
if (!(.n_in == 1278L && .empty == 4L && .usable == 1274L))
  message("[03] NOTE: counts differ from the canonical 1278/4/1274 — confirm the decisions input is the 2017 balanced one.")
cat("03_SANITIZE_GEOMETRIES done. STOP. Next: run 04_BUILD_POOLS.\n")
cat("===============================================================\n")


# #############################################################################
# 04_BUILD_POOLS — build_supervised_training_pools. Keeps B04_pools.
#   Prints positives, negatives per bucket, available, caps, selected, fingerprint.
# #############################################################################
B04_pools <- build_supervised_training_pools(
  config = B01_cfg, write_outputs = TRUE, overwrite = TRUE)

cat("\n================= 04_BUILD_POOLS CHECK =================\n")
.tl <- B04_pools$train_labeled
.n_burned <- sum(as.character(.tl$class) == "burned", na.rm = TRUE)
.n_unburned <- sum(as.character(.tl$class) == "unburned", na.rm = TRUE)
cat(sprintf("train_labeled rows : %d  (burned=%d, unburned=%d)\n",
            nrow(.tl), .n_burned, .n_unburned))
cat(sprintf("scoring_pool rows  : %d\n", nrow(B04_pools$scoring_pool)))
# Negative-pool buckets (all_sources policy). Source column name may vary by
# build; print whatever bucket/source breakdown is available.
.src_col <- intersect(c("neg_source", "source", "pool_source", "bucket"), names(.tl))
if (length(.src_col)) {
  cat("negative buckets (selected):\n")
  print(table(sf::st_drop_geometry(.tl)[.tl$class == "unburned", .src_col[1]]))
} else {
  cat("(negative bucket/source column not found on train_labeled; see pools QA CSVs in 01_POOLS)\n")
}
cat(sprintf("caps applied : rnd=%.2f otsu=%.2f\n",
            CAP_RANDOM, CAP_OTSU))
cat(sprintf("pools fingerprint : %s\n",
            B04_pools$fingerprint %||% B04_pools$pool_fingerprint %||% "(see pools sidecar)"))
stopifnot(.n_burned > 0)
cat("04_BUILD_POOLS done. STOP. Next: run 05_EXTRACT_FEATURES.\n")
cat("=======================================================\n")


# #############################################################################
# 05_EXTRACT_FEATURES — extract_supervised_features. Keeps B05_features.
#   Prints the feature columns (train + scoring).
#   NOTE: extraction always builds the FULL feature frame; the profile's feature
#   SUBSET is applied later via the recipe/whitelist (block 06), so NO_HOTSPOT
#   re-uses the same extracted frame (a subset selection, not a re-extraction).
# #############################################################################
# Folds are needed by extract_supervised_features; build them once here for the
# extraction call. (Block 07 re-derives folds independently for the leakage check.)
.folds_for_feats <- make_spatial_folds(
  train_labelled = B04_pools$train_labeled, config = B01_cfg,
  split_unit = "fire", out_dir = NULL, write_outputs = TRUE, overwrite = TRUE)
B05_features <- extract_supervised_features(
  train_with_folds = .folds_for_feats$train_with_folds,
  scoring_pool = B04_pools$scoring_pool, config = B01_cfg,
  use_hotspots = USE_HOTSPOTS, out_dir = NULL, write_outputs = TRUE)

cat("\n================= 05_EXTRACT_FEATURES CHECK =================\n")
cat(sprintf("train_features   : %d rows x %d cols\n",
            nrow(B05_features$train_features), ncol(B05_features$train_features)))
cat(sprintf("scoring_features : %d rows x %d cols\n",
            nrow(B05_features$scoring_features), ncol(B05_features$scoring_features)))
cat("train_features columns:\n")
print(names(B05_features$train_features))
.hs_in_frame <- intersect(HOTSPOT_LINEAGE_BASE, names(B05_features$train_features))
cat(sprintf("hotspot-lineage cols present in extracted frame: %d of 12\n", length(.hs_in_frame)))
cat("05_EXTRACT_FEATURES done. STOP. Next: run 06_FIT_AND_APPLY_RECIPE.\n")
cat("============================================================\n")


# #############################################################################
# 06_FIT_AND_APPLY_RECIPE — fit the shared recipe (medians/levels/weights/_isNA)
#   on the training features and inspect it. Keeps B06_recipe.
#   The recipe is the SINGLE shared object used by OOF, FINAL and scoring.
#   Here we obtain it by training the FINAL model (which fits+returns the recipe)
#   on a HOLD — instead we surface the recipe from a lightweight FINAL fit later.
#   To keep this block self-contained and cheap, we fit ONLY the recipe via the
#   FINAL trainer's recipe and inspect schema sizes.
#   Expected: FULL base=50 _isNA=50 total=100 ; NO_HOTSPOT base=38 _isNA=38 total=76.
# #############################################################################
# The recipe is produced as a side product of the FINAL model fit; we run a
# minimal FINAL fit here purely to materialize + inspect the recipe object.
# (Block 09 runs the canonical FINAL with the OOF aggregate; this fit is for
# recipe inspection and is kept separate so blocks stay independent.)
.recipe_dir <- file.path(MANUAL_OUT, "B06_recipe"); dir.create(.recipe_dir, showWarnings = FALSE, recursive = TRUE)
.b06_model <- train_final_burned_model(
  train_features = B05_features$train_features, config = B01_cfg,
  oof_agg = NULL, out_dir = .recipe_dir, overwrite = TRUE, verbose = FALSE)
B06_recipe <- .b06_model$recipe

cat("\n================= 06_FIT_AND_APPLY_RECIPE CHECK =================\n")
.fn   <- .b06_model$model$feature_names
.base <- .fn[!grepl("_isNA$", .fn)]; .isna <- .fn[grepl("_isNA$", .fn)]
cat(sprintf("recipe feature schema : base=%d  _isNA=%d  total=%d\n",
            length(.base), length(.isna), length(.fn)))
cat(sprintf("expected for %s        : base=%d _isNA=%d total=%d\n",
            MODEL_PROFILE, EXPECTED_BASE_N, EXPECTED_BASE_N, 2L * EXPECTED_BASE_N))
.fp <- digest::digest(paste(sort(.fn), collapse = "|"), algo = "md5")
cat(sprintf("feature_schema_fingerprint : %s\n", .fp))
.meds <- B06_recipe$medians %||% B06_recipe$imputation_medians
if (!is.null(.meds)) { cat("recipe medians (head):\n"); print(utils::head(.meds)) }
.lvls <- B06_recipe$levels %||% B06_recipe$factor_levels
if (!is.null(.lvls)) { cat("recipe factor levels:\n"); print(.lvls) }
.wts <- B06_recipe$feature_weights %||% B06_recipe$weights
if (!is.null(.wts)) { cat("feature weights (head):\n"); print(utils::head(.wts)) }
stopifnot(length(.base) == EXPECTED_BASE_N, length(.fn) == 2L * EXPECTED_BASE_N)
cat("06_FIT_AND_APPLY_RECIPE done. STOP. Next: run 07_BUILD_FOLDS.\n")
cat("================================================================\n")


# #############################################################################
# 07_BUILD_FOLDS — make_spatial_folds (independent derivation for the leakage
#   check). Keeps B07_folds. Prints groups, train/test sizes, leakage check.
# #############################################################################
B07_folds <- make_spatial_folds(
  train_labelled = B04_pools$train_labeled, config = B01_cfg,
  split_unit = "fire", out_dir = NULL, write_outputs = TRUE, overwrite = TRUE)

cat("\n================= 07_BUILD_FOLDS CHECK =================\n")
.twf <- B07_folds$train_with_folds
cat(sprintf("train_with_folds rows : %d\n", nrow(.twf)))
for (fc in FOLD_COLS) {
  if (fc %in% names(.twf)) {
    cat(sprintf("fold column %-10s groups:\n", fc))
    print(table(sf::st_drop_geometry(.twf)[[fc]], useNA = "ifany"))
  }
}
# Leakage check: a fire's polygons must not straddle train/test within a rep.
.grp_col <- intersect(c("fire_id", "fire", "group_id"), names(.twf))
if (length(.grp_col) && all(FOLD_COLS %in% names(.twf))) {
  .leak <- sapply(FOLD_COLS, function(fc) {
    tt <- tapply(sf::st_drop_geometry(.twf)[[fc]], sf::st_drop_geometry(.twf)[[.grp_col[1]]],
                 function(v) length(unique(v)))
    sum(tt > 1L)
  })
  cat("groups split across folds (should be 0 for each rep):\n"); print(.leak)
  stopifnot(all(.leak == 0L))
} else {
  cat("(group column not found; leakage check skipped — confirm split_unit='fire')\n")
}
cat("07_BUILD_FOLDS done. STOP. Next: run 08_RUN_OOF.\n")
cat("=======================================================\n")


# #############################################################################
# 08_RUN_OOF — run_oof_diagnostics. Keeps B08_oof.
#   Prints caps, scale_pos_weight, best_iteration, schema fingerprint,
#   refit, predicted rows, OOF metrics (with the hotspot-label-circularity caveat).
#
#   ELIGIBILITY (internal, no functional change to drive here): training rows are
#   chosen by EXPLICIT class, never by negation. class=="burned" -> positive;
#   class=="unburned" + a valid bucket (contextual/spectral/random/otsu) ->
#   negative; otsu review/keep -> excluded+logged; NA/unknown/unbucketed -> error.
#   Review/keep/NA/unknown can NEVER become negatives. OOF and FINAL (block 09)
#   share the SAME internal eligibility resolver + the SAME capping helper, so the
#   two stages cannot diverge on eligibility or caps (they differ only in the
#   retained outer-test fold and the seed).
# #############################################################################
.oof_dir <- file.path(MANUAL_OUT, "B08_oof/05_OOF"); dir.create(.oof_dir, showWarnings = FALSE, recursive = TRUE)
.mat_dir <- file.path(MANUAL_OUT, "B08_oof/04_MATRIX"); dir.create(.mat_dir, showWarnings = FALSE, recursive = TRUE)
B08_oof <- run_oof_diagnostics(
  train_features = B05_features$train_features, scoring_features = B05_features$scoring_features,
  config = B01_cfg, fold_cols = FOLD_COLS, out_dir = .oof_dir, matrix_dir = .mat_dir,
  labelled_gpkg = B07_folds$train_with_folds_gpkg, labelled_layer = "train_with_folds")

cat("\n================= 08_RUN_OOF CHECK =================\n")
cat("training          : early-stopping selection + full-data refit (fixed)\n")
cat(sprintf("oof_sampling (fixed): %s  (same capped policy as FINAL; not a knob)\n",
            B01_cfg$train_control$oof_sampling %||% "capped"))
cat(sprintf("scale_pos_weight  : %s\n", B08_oof$scale_pos_weight %||% B08_oof$spw %||% "(see audit)"))
cat(sprintf("best_iteration    : %s\n", B08_oof$best_iteration %||% "(see audit)"))
cat(sprintf("schema fingerprint: %s\n", B08_oof$feature_schema_fingerprint %||% "(see fingerprints CSV)"))
cat(sprintf("refit             : %s\n", B08_oof$refit %||% "(nested_refit)"))
cat(sprintf("predicted rows    : %s\n", B08_oof$n_predicted %||% nrow(B08_oof$oof_agg %||% data.frame())))
cat("OOF metrics:\n")
print(B08_oof$metrics %||% B08_oof$oof_metrics %||% "(see oof_metrics_by_threshold.csv)")
cat("\nCAVEAT: hotspot-derived labels can make hotspot features partly circular;\n")
cat("        OOF metrics for the FULL profile must be read with that caveat. The\n")
cat("        external EFFIS validation (block 13) is the unbiased accuracy check.\n")
cat("08_RUN_OOF done. STOP. Next: run 09_TRAIN_FINAL_AND_REFIT.\n")
cat("===================================================\n")


# #############################################################################
# 09_TRAIN_FINAL_AND_REFIT — train_final_burned_model with the OOF aggregate.
#   Keeps B09_final. Prints best_iteration, fingerprint, refit, saved paths.
#
#   ELIGIBILITY (internal): the FINAL pool builder selects its negatives through
#   the SAME shared eligibility resolver + capping helper as OOF (block 08) —
#   explicit burned positives + explicit unburned-with-valid-bucket negatives;
#   review/keep/NA/unknown never become negatives. No knob to set here.
# #############################################################################
.final_dir <- file.path(MANUAL_OUT, "B09_final/07_FINAL_MODEL"); dir.create(.final_dir, showWarnings = FALSE, recursive = TRUE)
B09_final <- train_final_burned_model(
  train_features = B05_features$train_features, config = B01_cfg,
  oof_agg = B08_oof$oof_agg_csv, out_dir = .final_dir, overwrite = TRUE, verbose = TRUE)

cat("\n================= 09_TRAIN_FINAL_AND_REFIT CHECK =================\n")
.fn9 <- B09_final$model$feature_names
.fp9 <- digest::digest(paste(sort(.fn9), collapse = "|"), algo = "md5")
cat(sprintf("best_iteration : %s\n", B09_final$recipe$training$best_iteration %||% "?"))
cat(sprintf("feature schema : base=%d total=%d  fingerprint=%s\n",
            sum(!grepl("_isNA$", .fn9)), length(.fn9), .fp9))
cat(sprintf("refit          : %s\n", B09_final$recipe$training$refit %||% "(nested_refit)"))
cat(sprintf("saved model    : %s\n", B09_final$model_path %||% B09_final$saved_model_path %||% "(in out_dir)"))
cat(sprintf("saved recipe   : %s\n", B09_final$recipe_path %||% B09_final$saved_recipe_path %||% "(in out_dir)"))
cat("09_TRAIN_FINAL_AND_REFIT done. STOP. Next: run 10_SAVE_LOAD_MODEL_RECIPE.\n")
cat("=================================================================\n")


# #############################################################################
# 10_SAVE_LOAD_MODEL_RECIPE — save + reload the model/recipe and assert the
#   reloaded objects produce IDENTICAL predictions. Keeps B10_saveload.
# #############################################################################
.sl_dir <- file.path(MANUAL_OUT, "B10_saveload"); dir.create(.sl_dir, showWarnings = FALSE, recursive = TRUE)
.model_rds  <- file.path(.sl_dir, "final_model.rds")
.recipe_rds <- file.path(.sl_dir, "final_recipe.rds")
saveRDS(B09_final$model,  .model_rds)
saveRDS(B09_final$recipe, .recipe_rds)
.model2  <- readRDS(.model_rds)
.recipe2 <- readRDS(.recipe_rds)
# Score a small slice through both the in-memory and reloaded model to compare.
.score_one <- function(model, recipe, feats) {
  # apply_supervised_recipe(df, recipe) is the internal recipe applier (fit once,
  # applied identically to OOF/FINAL/scoring). Returns the model-ready matrix.
  mat <- tryCatch(
    OtsuFire:::apply_supervised_recipe(feats, recipe),
    error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  as.numeric(predict(model, newdata = as.matrix(mat)))
}
.slice <- utils::head(B05_features$scoring_features, 200L)
.p1 <- tryCatch(.score_one(B09_final$model, B09_final$recipe, .slice), error = function(e) NULL)
.p2 <- tryCatch(.score_one(.model2, .recipe2, .slice), error = function(e) NULL)
B10_saveload <- list(model_rds = .model_rds, recipe_rds = .recipe_rds,
                     identical_preds = !is.null(.p1) && !is.null(.p2) && isTRUE(all.equal(.p1, .p2)))

cat("\n================= 10_SAVE_LOAD_MODEL_RECIPE CHECK =================\n")
cat(sprintf("saved model  : %s\n", .model_rds))
cat(sprintf("saved recipe : %s\n", .recipe_rds))
if (is.null(.p1) || is.null(.p2)) {
  cat("NOTE: could not score the slice via the internal recipe applier in this build;\n")
  cat("      fall back to comparing xgb.Booster raw bytes instead.\n")
  .raw1 <- xgboost::xgb.save.raw(B09_final$model)
  .raw2 <- xgboost::xgb.save.raw(.model2)
  B10_saveload$identical_preds <- identical(.raw1, .raw2)
  cat(sprintf("reloaded model raw-bytes identical : %s\n", B10_saveload$identical_preds))
} else {
  cat(sprintf("reloaded predictions identical (n=%d) : %s\n", length(.p1), B10_saveload$identical_preds))
}
stopifnot(isTRUE(B10_saveload$identical_preds))
cat("10_SAVE_LOAD_MODEL_RECIPE done. STOP. Next: run 11_SCORE_POLYGONS.\n")
cat("==================================================================\n")


# #############################################################################
# 11_SCORE_POLYGONS — score_supervised_burned_map. Keeps B11_scored.
#   Pass labelled_features = the EXPLICIT features GPKG from cfg's route (so the
#   score step never silently reconstructs a wrong path). Prints input rows,
#   score rows, id/order preservation, NO silent loss.
# #############################################################################
# The labelled-features GPKG written by block 05/07 lives at cfg's route.
.labelled_features_path <- B01_cfg$output_routes$features_geometry_gpkg
stopifnot(file.exists(.labelled_features_path))
.sc_dir  <- file.path(MANUAL_OUT, "B11_scored/08_SCORED"); dir.create(.sc_dir, showWarnings = FALSE, recursive = TRUE)
.map_dir <- file.path(MANUAL_OUT, "B11_scored/09_FINAL_MAP"); dir.create(.map_dir, showWarnings = FALSE, recursive = TRUE)
B11_scored <- score_supervised_burned_map(
  scoring_features = B05_features$scoring_features, model = B09_final$model, recipe = B09_final$recipe,
  config = B01_cfg, oof_summary = B08_oof$labeled_oof_summary_gpkg,
  labelled_features = .labelled_features_path, export_burned_like = TRUE,
  out_map_dir = .map_dir, out_score_dir = .sc_dir, overwrite = TRUE, verbose = TRUE)

cat("\n================= 11_SCORE_POLYGONS CHECK =================\n")
.fmf <- B11_scored$final_map_full
cat(sprintf("scoring input rows : %d\n", nrow(B05_features$scoring_features)))
cat(sprintf("scored map rows    : %d\n", nrow(.fmf)))
cat(sprintf("labelled_features  : %s\n", .labelled_features_path))
.id_col <- intersect(c("poly_id", "id", "geom_id"), names(.fmf))
if (length(.id_col)) {
  cat(sprintf("ids preserved (n unique) : %d\n", length(unique(.fmf[[.id_col[1]]]))))
}
cat("Expected (2017 balanced FULL): 1274 input -> 1274 scored (no silent loss).\n")
if (nrow(.fmf) != nrow(B05_features$scoring_features))
  message("[11] NOTE: scored rows != scoring rows — investigate before trusting the map.")
cat("11_SCORE_POLYGONS done. STOP. Next: run 12_EXPORT_MAP.\n")
cat("==========================================================\n")


# #############################################################################
# 12_EXPORT_MAP — inspect the exported map layers. Keeps B12_map.
#   Layers: deterministic_scored (all scored), final_map_full (all),
#   final_map (current-year public / thresholded). Prints counts, the filter
#   applied, and the ids removed + reason (the current_year_public_drop temporal
#   contract — e.g. the 1274 vs 1273 finding).
# #############################################################################
B12_map <- list(
  deterministic_scored = B11_scored$deterministic_scored %||% B11_scored$scored,
  final_map_full       = B11_scored$final_map_full,
  final_map            = B11_scored$final_map)

cat("\n================= 12_EXPORT_MAP CHECK =================\n")
.nz <- function(x) if (is.null(x)) NA_integer_ else nrow(x)
cat(sprintf("deterministic_scored rows : %s (all scored polygons)\n", .nz(B12_map$deterministic_scored)))
cat(sprintf("final_map_full rows       : %s (all polygons, unfiltered)\n", .nz(B12_map$final_map_full)))
cat(sprintf("final_map rows            : %s (current-year public / thresholded)\n", .nz(B12_map$final_map)))
.full_n <- .nz(B12_map$final_map_full); .pub_n <- .nz(B12_map$final_map)
if (!is.na(.full_n) && !is.na(.pub_n) && .full_n != .pub_n) {
  .removed <- .full_n - .pub_n
  cat(sprintf("filter applied : current_year_public_drop -> removed %d polygon(s)\n", .removed))
  cat("reason : temporal contract (current-year public layer drops polygons not\n")
  cat("         meeting the public/current-year rule; e.g. the 1274 vs 1273 finding).\n")
  # Try to surface which id(s) were dropped.
  .idc <- intersect(c("poly_id", "id", "geom_id"), names(B12_map$final_map_full))
  if (length(.idc) && length(intersect(.idc, names(B12_map$final_map)))) {
    .dropped <- setdiff(B12_map$final_map_full[[.idc[1]]], B12_map$final_map[[.idc[1]]])
    cat(sprintf("dropped id(s)  : %s\n", paste(utils::head(.dropped, 20), collapse = ", ")))
  }
} else {
  cat("filter applied : none (final_map == final_map_full this run)\n")
}
cat("12_EXPORT_MAP done. STOP. Next: run 13_VALIDATE_EFFIS.\n")
cat("======================================================\n")


# #############################################################################
# 13_VALIDATE_EFFIS — validate_fire_maps in the selected EFFIS_MODE.
#   Keeps B13_effis. Prints resolution, CRS, observability, omission,
#   commission, F1, IoU.
#
#   STANDARD_90M (default)      : EPSG:3035, 90 m, terra todisk + modest memfrac +
#                                 dedicated tempdir + memory logging + safe cleanup.
#   NATIVE_HIGH_MEMORY          : high-consumption; guarded behind ALLOW_HIGH_MEMORY.
#   Both modes call the SAME validate_fire_maps() with the SAME metric definitions;
#   only the resolution / memory strategy differ.
# #############################################################################
.eff_dir <- file.path(MANUAL_OUT, "B13_effis"); dir.create(.eff_dir, showWarnings = FALSE, recursive = TRUE)
.ttmp    <- file.path(.eff_dir, "terra_tmp"); dir.create(.ttmp, showWarnings = FALSE, recursive = TRUE)

.log_mem <- function(tag) {
  free <- tryCatch(round(as.numeric(suppressWarnings(system2("powershell",
    c("-NoProfile","-Command","(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"),
    stdout = TRUE))[1]) / 1e6, 1), error = function(e) NA)
  cat(sprintf("[effis-mem][%s] free physical mem GB = %s\n", tag, free))
}

if (EFFIS_MODE == "NATIVE_HIGH_MEMORY") {
  if (!isTRUE(ALLOW_HIGH_MEMORY))
    stop("[13] EFFIS_MODE=NATIVE_HIGH_MEMORY requires ALLOW_HIGH_MEMORY <- TRUE (explicit opt-in). It is HIGH-CONSUMPTION; STANDARD_90M is recommended.")
  warning("[13] EFFIS NATIVE_HIGH_MEMORY: HIGH memory consumption. Native-resolution rasterization may OOM on large extents.")
  .effis_burnable <- burnable_mask_run   # caller may swap to a native-res mask if intended
  .effis_memfrac  <- 0.8
} else {
  # STANDARD_90M (recommended): 3035 90 m aligned burnable, modest memfrac.
  .effis_burnable <- burnable_mask_run
  .effis_memfrac  <- EFFIS_MEMFRAC
}

terra::terraOptions(todisk = TRUE, memfrac = .effis_memfrac, tempdir = .ttmp, progress = 0)
.log_mem("before")
cat(sprintf("[13] EFFIS_MODE=%s  memfrac=%.2f  tempdir=%s\n", EFFIS_MODE, .effis_memfrac, .ttmp))

# Threshold the public/full map at the operating threshold, then validate.
.fm <- B12_map$final_map_full
.p  <- suppressWarnings(as.numeric(sf::st_drop_geometry(.fm)$p_burned_model))
.keep <- .fm[is.finite(.p) & .p >= OPERATING_THRESHOLD, , drop = FALSE]
.keep <- sf::st_make_valid(sf::st_zm(.keep, drop = TRUE, what = "ZM"))
.keep <- .keep[!sf::st_is_empty(.keep), , drop = FALSE]
.keep_path <- file.path(.eff_dir, "thresholded_burned.gpkg")
if (file.exists(.keep_path)) unlink(.keep_path, force = TRUE)
sf::st_write(.keep, .keep_path, layer = "thresholded_burned", quiet = TRUE)

.mask_sf <- sf::st_make_valid(sf::st_zm(sf::st_read(validation_mask_shapefile, quiet = TRUE),
                                        drop = TRUE, what = "ZM"))
.mask_sf <- .mask_sf[!sf::st_is_empty(.mask_sf), 0, drop = FALSE]; .mask_sf$mask_id <- seq_len(nrow(.mask_sf))
.mask_path <- file.path(.eff_dir, "validation_mask_geometry_only.gpkg")
if (file.exists(.mask_path)) unlink(.mask_path, force = TRUE)
sf::st_write(.mask_sf, .mask_path, layer = "mask", quiet = TRUE)

.has_strata <- file.exists(strata_raster) && file.exists(strata_lut)
B13_effis <- validate_fire_maps(
  input_shapefile = .keep_path, ref_shapefile = reference_burned_map,
  mask_shapefile = .mask_path,
  burnable_raster = .effis_burnable,
  year_target = target_year, validation_dir = .eff_dir,
  force_reprocess_ref = TRUE, force_reprocess_pred = TRUE,
  metrics_type = "all", dissolve_ref_by = "id", dissolve_input_by = NULL,
  strata_raster = if (.has_strata) strata_raster else NULL,
  strata_lut    = if (.has_strata) strata_lut    else NULL,
  observability_raster = change_index,
  ref_end_doy_col = "end_doy", ref_start_doy_col = "start_doy")
.log_mem("after")
# safe cleanup of the terra tempdir
tryCatch(unlink(.ttmp, recursive = TRUE, force = TRUE), error = function(e) NULL)

cat("\n================= 13_VALIDATE_EFFIS CHECK =================\n")
.m <- B13_effis$metrics %||% B13_effis
.get <- function(nm) {
  v <- tryCatch(.m[[nm]], error = function(e) NULL)
  if (is.null(v)) tryCatch(.m[[which(tolower(names(.m)) == tolower(nm))[1]]], error = function(e) NA) else v
}
cat(sprintf("EFFIS_MODE     : %s\n", EFFIS_MODE))
cat(sprintf("resolution (m) : 90 (3035 aligned)\n"))
cat(sprintf("CRS            : EPSG:3035\n"))
cat(sprintf("observability  : %s\n", .get("observability") %||% "(see 02_OBSERVABILITY)"))
cat(sprintf("omission       : %s\n", .get("omission") %||% .get("omission_error")))
cat(sprintf("commission     : %s\n", .get("commission") %||% .get("commission_error")))
cat(sprintf("F1             : %s\n", .get("F1") %||% .get("f1")))
cat(sprintf("IoU            : %s\n", .get("IoU") %||% .get("iou")))
cat("Full metric table:\n"); print(.m)
cat("13_VALIDATE_EFFIS done. STOP. Pipeline complete.\n")
cat("==========================================================\n")
