# =============================================================================
# LONG <YEAR> BALANCED — COMMON HEADER (sourced by every stage script)
# CANONICAL VERSIONED TEMPLATE. The single source of truth lives here
# (inst/scripts/long_run/00_common.R). Any copy under 00_LONG_RUNS/<run>/
# _ORCHESTRATION/ is a WORKING COPY generated from this template.
#
# Paired experiment: FULL (thermal-enhanced, all features incl hotspots) vs
# NO-HOTSPOT (spectral-historical ablation).
# FROZEN CODE: an INSTALLED OtsuFire in an ISOLATED library (NOT load_all).
# This header pins the version, resolves the official paths/caps/seeds, defines
# the feature lineage (FULL 50 base vs NO-HOTSPOT 38 base), provides logging
# helpers, and exposes the SHARED upstream products as EXPLICIT absolute routes
# (the G2.3 route-propagation fix: no junctions). It performs NO heavy compute.
#
# ===========================================================================
# CONFIG BLOCK — EDIT EVERY <PATH_TO_...> / <...> PLACEHOLDER BEFORE RUNNING.
# These are TEMPLATE placeholders on purpose: the canonical copy must NOT bake
# any one specific run's absolute paths. Fill them for your run, then source.
# ===========================================================================
options(warn = 1)

# ---- ROOT + ISOLATED LIB (REQUIRED) -----------------------------------------
# LONG_ROOT: the run root that will hold SHARED/ + the two profile folders.
LONG_ROOT <- "<PATH_TO_LONG_RUN_ROOT>"      # e.g. C:/.../00_LONG_RUNS/LONG_<YEAR>_BALANCED_<STAMP>
# RLIB: an ISOLATED library that contains ONLY the frozen OtsuFire build to use.
RLIB      <- "<PATH_TO_ISOLATED_RLIB>"      # e.g. C:/.../00_SMOKE/<smoke>/Rlib  (or the 0.6.x snapshot lib)
.libPaths(c(RLIB, .libPaths()))

suppressPackageStartupMessages({
  library(OtsuFire, lib.loc = RLIB)
  library(sf); library(terra)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- VERSION PIN (REQUIRED) -------------------------------------------------
# The frozen OtsuFire version this run is pinned to (must be the version in RLIB).
PIN_VERSION <- "0.6.1"   # bump to "0.6.2" when the 0.6.2 snapshot lib is used.
.assert_version_pin <- function() {
  ver <- as.character(packageVersion("OtsuFire"))
  pkgpath <- normalizePath(find.package("OtsuFire"), winslash = "/")
  in_iso  <- startsWith(pkgpath, normalizePath(RLIB, winslash = "/"))
  cat(sprintf("[pin] packageVersion=%s | find.package=%s | isolated=%s\n",
              ver, pkgpath, in_iso))
  stopifnot(ver == PIN_VERSION, in_iso)
  invisible(TRUE)
}

# ---- SNAPSHOT PROVENANCE (REQUIRED for the env snapshot) --------------------
SNAPSHOT_ID     <- "<SNAPSHOT_ID>"                 # e.g. GATE1_FINAL_v0.6.1_<STAMP>
SNAPSHOT_DIR    <- "<PATH_TO_SNAPSHOT_DIR>"
TARBALL_PATH    <- file.path(SNAPSHOT_DIR, "PACKAGE", sprintf("OtsuFire_%s.tar.gz", PIN_VERSION))
TARBALL_SHA256  <- "<TARBALL_SHA256>"

# ---- OFFICIAL DATA PATHS (REQUIRED) -----------------------------------------
base_dir       <- "<PATH_TO_FIRE_MAPPING_BASE>"    # e.g. C:/.../00_FIRE_MAPPING
data_base      <- file.path(base_dir, "1_DATA")
results_base   <- file.path(data_base, "Results")
composite_base <- file.path(data_base, "Imagery", "Composites_90m")

target_year <- 2017L          # <YEAR>
scenario    <- "balanced"
result_name <- "Min_Min"
corine_year <- 2012
use_hotspots_FULL <- TRUE      # FULL profile: hotspots ON (only valid for year > 2000)

internal_decisions <- file.path(results_base, as.character(target_year), "DETERMINISTIC",
                                sprintf("%d_%s", target_year, scenario),
                                "05_DECISIONS", "internal_decisions.gpkg")
change_index       <- file.path(composite_base, result_name,
                                sprintf("MinMin_%d_mosaic_res90m.tif", target_year))
hotspots           <- file.path(data_base, "Hotspots",
                                sprintf("hotspots_iberia_%d.geojson", target_year))
reference_burned_map <- file.path(data_base, "Fires/Validation_fires_burneable_verano",
                                  sprintf("Effis_CA_%d_maskKeep_summer.shp", target_year))
# EFFIS native-res WGS84 burnable (potential OOM source) -- provenance only.
burnable_mask_wgs84 <- file.path(data_base, "Corine_Masks",
                                 "burneable_mask_binary_corine12_wgs84.tif")
# 3035 90m composite-aligned burnable: used for BOTH the supervised RUN and the
# EFFIS validation (memory-safe).
burnable_mask_run  <- file.path(data_base, "Corine_Masks",
                                sprintf("burneable_mask_binary_corine12_3035_alignedMinMin%d.tif",
                                        target_year))
validation_mask_shapefile <- file.path(data_base, "Mask_StudyArea", "mask_Peninsula_3035.shp")
strata_raster <- file.path(data_base, "Corine_Masks/STRATA", "strata_CLC_2012_res30.tif")
strata_lut    <- file.path(data_base, "Corine_Masks/LUT/lut_full_strata8_v1.csv")
# topo / corine / delayed_change_index are resolved CANONICALLY by the cfg builder
# (cfg$inputs); these mirror that resolution for the env-snapshot hash.
topo_path   <- file.path(data_base, "Topography", "elevation_slope.tif")
corine_raster <- file.path(data_base, "Corine_Masks", "CLC_2012_peninsula.tif")
delayed_change_index <- file.path(composite_base, "Autumn",
                                  sprintf("mean_mean_%d_mosaic.tif", target_year))

# ---- EXTERNAL TOOL PATHS (REQUIRED; used by the legacy all_sources Otsu) -----
python_exe             <- "<PATH_TO_python.exe>"
gdal_polygonize_script <- "<PATH_TO_gdal_polygonize.py>"
gdalwarp_path          <- "<PATH_TO_gdalwarp.exe>"
ogr2ogr_exe            <- "<PATH_TO_ogr2ogr.exe>"

# ---- CANONICAL CONTROLS (FULL training; the canonical, NOT smoke values) -----
CAP_CONTEXTUAL <- 0.25; CAP_SPECTRAL <- 2.0; CAP_RANDOM <- 1.0; CAP_OTSU <- 1.0
SEED_BASE      <- 42L
FOLD_COLS      <- c("fold_rep1", "fold_rep2")
OPERATING_THRESHOLD <- 0.50
NROUNDS_MAX    <- 4000L   # CANONICAL full training (smoke used 60)
EARLY_STOP     <- 80L     # CANONICAL (smoke used 20)
DROP_RESCUE_THRESHOLD <- 0.50

# FIXED nthread for BOTH profiles (fair comparison). Adjust to leave headroom for
# terra/sf during EFFIS (8 logical CPUs -> 6 is a safe default).
XGB_NTHREAD <- 6L
Sys.setenv(OMP_NUM_THREADS = as.character(XGB_NTHREAD))

# ===========================================================================
# END CONFIG BLOCK. Nothing below needs per-run editing.
# ===========================================================================

# ---- FEATURE LINEAGE (derived from FROZEN code, not just prefix) ------------
# Source of truth: .supervised_feature_cols (R/internal-sup-residual-cols.R),
# the canonical 50-name whitelist (0.6.2: hs_any REMOVED; resolved space
# UNCHANGED, still 50 base + 50 _isNA = 100). doy_* features come from doy_post
# (change-index band 2); they are SPECTRAL/temporal and STAY in BOTH profiles.
# The 12 hotspot-lineage base features are produced EXCLUSIVELY from the
# active-fire/MODIS hotspot data path (extract-features hotspot_features()).
HOTSPOT_LINEAGE_BASE <- c(
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support"
)  # 12 genuine hotspot features (hs_any is GONE from the whitelist; removed-count == 12).

# The exact 50 base features the FULL model sees (50 base + 50 _isNA).
FULL_BASE_50 <- c(
  # RBR same-window (5)
  "rbr_valid_frac", "rbr_p10", "rbr_med", "rbr_p90", "rbr_iqr",
  # CORINE (7)
  "cor_open_frac", "cor_wetlands_frac", "cor_water_frac", "cor_herbaceous_frac",
  "cor_urban_frac", "cor_agri_frac", "cor_forest_frac",
  # elevation (7)
  "elev_valid_frac", "elev_p10", "elev_med", "elev_p90", "elev_iqr", "elev_mean", "elev_sd",
  # slope (7)
  "slope_valid_frac", "slope_p10", "slope_med", "slope_p90", "slope_iqr", "slope_mean", "slope_sd",
  # DOY (5) -- SPECTRAL/temporal (change-index band 2); STAY
  "doy_valid_frac", "doy_p10", "doy_med", "doy_p90", "doy_iqr",
  # RBR all-window + persistence (7)
  "rbr_aw_valid_frac", "rbr_aw_p10", "rbr_aw_med", "rbr_aw_p90", "rbr_aw_iqr",
  "persist_delta", "persist_ratio",
  # Hotspots (12) -- REMOVED in NO-HOTSPOT
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support"
)
stopifnot(length(FULL_BASE_50) == 50L,
          all(HOTSPOT_LINEAGE_BASE %in% FULL_BASE_50))

# NO-HOTSPOT base = 50 minus the 12 hotspot-lineage features = 38. Passed as
# feature_whitelist_override (provenance "user"; a subset of the 50-name canonical
# whitelist). Its _isNA companions follow automatically.
NOHS_BASE_38 <- setdiff(FULL_BASE_50, HOTSPOT_LINEAGE_BASE)
stopifnot(length(NOHS_BASE_38) == 38L)

# Lineage table (Feature | Family | FULL | NO-HOTSPOT | Source).
.feature_family <- function(f) {
  if (f %in% HOTSPOT_LINEAGE_BASE) return("hotspot")
  if (grepl("^rbr_aw_", f)) return("rbr_allwindow")
  if (grepl("^rbr_", f))    return("rbr_samewindow")
  if (grepl("^persist", f)) return("persistence")
  if (grepl("^cor_", f))    return("corine")
  if (grepl("^elev_", f))   return("topo_elevation")
  if (grepl("^slope_", f))  return("topo_slope")
  if (grepl("^doy_", f))    return("doy_spectral")
  "other"
}
.feature_source <- function(f) {
  fam <- .feature_family(f)
  switch(fam,
    hotspot        = "MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT",
    rbr_allwindow  = "delayed_change_index RBR all-window composite (spectral)",
    rbr_samewindow = "change_index RBR summer band 1 (spectral)",
    persistence    = "RBR same-window vs all-window persistence (spectral)",
    corine         = "CORINE land-cover fractions (zonal)",
    topo_elevation = "DEM elevation zonal stats (topography)",
    topo_slope     = "Slope zonal stats (topography)",
    doy_spectral   = "change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS)",
    "other")
}
build_feature_lineage <- function() {
  data.frame(
    Feature    = FULL_BASE_50,
    Family     = vapply(FULL_BASE_50, .feature_family, character(1)),
    FULL       = "in",
    NO_HOTSPOT = ifelse(FULL_BASE_50 %in% HOTSPOT_LINEAGE_BASE, "out", "in"),
    Source     = vapply(FULL_BASE_50, .feature_source, character(1)),
    row.names  = NULL, stringsAsFactors = FALSE
  )
}

# ---- LOGGING HELPER (per-script log + console) ------------------------------
make_logger <- function(log_file) {
  function(...) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
    cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
  }
}

# ---- LONG-RUN engine output base (keeps ALL writes under the long-run tree) --
# The cfg's canonical output_routes derive from output_dir. We point output_dir
# at a long-run results base so the engine writes pools/folds/features/sidecars
# UNDER the long-run folders (NOT the read-only 1_DATA/Results tree). A fresh
# pools_dir also means NO stale negative-pool fingerprint sidecar.
LONG_RESULTS_BASE <- file.path(LONG_ROOT, "SHARED", "ENGINE_ROUTES")
dir.create(LONG_RESULTS_BASE, recursive = TRUE, showWarnings = FALSE)

# ---- SHARED UPSTREAM CANONICAL ROUTES (G2.3 route-propagation fix) -----------
# The shared upstream (pools / folds / features) is built ONCE under SHARED with
# out_base = LONG_RESULTS_BASE. Each profile builds its OWN cfg with a DIFFERENT
# out_base (its own ENGINE_ROUTES), so cfg$output_routes$features_geometry_gpkg
# reconstructs to the PROFILE folder, NOT to SHARED. The scoring step reads the
# labelled-features GPKG from config$output_routes$features_geometry_gpkg when
# labelled_features = NULL -- that was the SHARED-features-not-visible failure.
#
# CANONICAL FIX (no junctions, no convention reconstruction): expose the SHARED
# upstream products as EXPLICIT absolute paths here; each profile passes the
# SHARED features GPKG to score_supervised_burned_map(labelled_features = ...)
# and validates+hashes every shared input BEFORE any heavy compute.
SHARED_SUP_BASE <- file.path(LONG_RESULTS_BASE, as.character(target_year),
                             result_name, "SUPERVISED", scenario)
SHARED_FEATURES_GPKG  <- file.path(SHARED_SUP_BASE, "03_FEATURES",
                                   "features_geometry.gpkg")
SHARED_TRAIN_FEATURES_RDS   <- file.path(SHARED_SUP_BASE, "03_FEATURES",
                                         "train_features.rds")
SHARED_SCORING_FEATURES_RDS <- file.path(SHARED_SUP_BASE, "03_FEATURES",
                                         "scoring_features.rds")
SHARED_FOLDS_GPKG     <- file.path(SHARED_SUP_BASE, "02_FOLDS",
                                   sprintf("%d_train_with_folds_5000m.gpkg",
                                           target_year))
SHARED_POOLS_GPKG     <- file.path(SHARED_SUP_BASE, "01_POOLS",
                                   sprintf("%d_%s_pools.gpkg",
                                           target_year, scenario))
# RDS handoffs of the in-memory products (written by 10_shared_upstream.R).
SHARED_POOLS_RDS <- file.path(LONG_ROOT, "SHARED", "pools.rds")
SHARED_FOLDS_RDS <- file.path(LONG_ROOT, "SHARED", "folds.rds")
SHARED_FEATS_RDS <- file.path(LONG_ROOT, "SHARED", "feats.rds")

# Map of the shared inputs each profile MUST resolve from SHARED (label -> path).
shared_inputs_map <- function() {
  list(
    features_geometry_gpkg = SHARED_FEATURES_GPKG,
    train_features_rds     = SHARED_TRAIN_FEATURES_RDS,
    scoring_features_rds   = SHARED_SCORING_FEATURES_RDS,
    folds_gpkg             = SHARED_FOLDS_GPKG,
    pools_gpkg             = SHARED_POOLS_GPKG,
    pools_rds              = SHARED_POOLS_RDS,
    folds_rds              = SHARED_FOLDS_RDS,
    feats_rds              = SHARED_FEATS_RDS
  )
}

# Fail-fast existence + SHA256 validation of every shared input. Aborts with a
# clear, actionable error BEFORE any heavy compute if a path is unresolvable or
# unhashable. Returns a data.frame and (optionally) writes a provenance CSV.
validate_shared_inputs <- function(log = NULL, out_csv = NULL,
                                    map = shared_inputs_map()) {
  emit <- if (is.function(log)) log else function(...) cat(..., "\n")
  rows <- vector("list", length(map))
  for (i in seq_along(map)) {
    label <- names(map)[i]; path <- map[[i]]
    exists_ok <- is.character(path) && length(path) == 1L &&
                 nzchar(path) && file.exists(path)
    sha <- NA_character_; sz <- NA_real_
    if (exists_ok) {
      sz  <- tryCatch(as.numeric(file.info(path)$size), error = function(e) NA_real_)
      sha <- tryCatch(
        as.character(digest::digest(file = path, algo = "sha256")),
        error = function(e) NA_character_)
    }
    rows[[i]] <- data.frame(label = label, path = path, exists = exists_ok,
                            sha256 = sha, size_bytes = sz,
                            stringsAsFactors = FALSE)
    emit(sprintf("[shared-input] %-22s exists=%-5s sha256=%s",
                 label, exists_ok,
                 if (is.na(sha)) "<MISSING>" else substr(sha, 1, 16)))
  }
  df <- do.call(rbind, rows)

  missing <- df[!df$exists, , drop = FALSE]
  if (nrow(missing) > 0L) {
    stop(sprintf(paste0(
      "FAIL-FAST: %d SHARED upstream input(s) unresolvable BEFORE heavy compute.\n",
      "The canonical architecture requires every shared product to resolve from ",
      "its SHARED absolute path (no junctions). Missing:\n%s\n",
      "Re-run 10_shared_upstream.R or fix SHARED_* in 00_common.R."),
      nrow(missing),
      paste(sprintf("  - %s: %s", missing$label, missing$path),
            collapse = "\n")), call. = FALSE)
  }
  no_hash <- df[df$exists & is.na(df$sha256), , drop = FALSE]
  if (nrow(no_hash) > 0L) {
    stop(sprintf(paste0(
      "FAIL-FAST: %d SHARED input(s) exist but their SHA256 could not be ",
      "computed (unreadable / locked):\n%s"),
      nrow(no_hash),
      paste(sprintf("  - %s: %s", no_hash$label, no_hash$path),
            collapse = "\n")), call. = FALSE)
  }
  if (!is.null(out_csv)) {
    utils::write.csv(df, out_csv, row.names = FALSE)
    emit(sprintf("[shared-input] provenance written: %s", out_csv))
  }
  invisible(df)
}

# ---- CFG BUILDERS (FULL vs NO-HOTSPOT differ ONLY via explicit config) -------
# Both use the single canonical training procedure (inner-early-stopping
# selection + full-data refit), identical caps/seeds. The ONLY difference is
# feature_whitelist_override (NO-HOTSPOT) + use_hotspots flag.
.mk_long_cfg <- function(oof_sampling,
                         feature_whitelist_override = NULL,
                         out_base = LONG_RESULTS_BASE) {
  cfg <- build_supervised_burned_config(
    scenario             = scenario,
    internal_decisions   = internal_decisions,
    change_index         = change_index,
    hotspots             = hotspots,
    reference_burned_map = reference_burned_map,
    target_year          = target_year,
    output_dir           = out_base,
    run_name             = result_name,
    burnable_mask        = burnable_mask_run,
    nrounds_max          = NROUNDS_MAX,
    early_stop           = EARLY_STOP,
    cap_contextual = CAP_CONTEXTUAL, cap_spectral = CAP_SPECTRAL,
    cap_random = CAP_RANDOM, cap_otsu = CAP_OTSU,
    oof_seed_base = SEED_BASE, final_sampling_seed = SEED_BASE, final_seed = SEED_BASE,
    oof_sampling = oof_sampling,
    feature_whitelist_override = feature_whitelist_override,
    options = list(
      data_base = data_base, composite_base = composite_base, result_name = result_name,
      legacy_otsu_mode = "burnable_only", legacy_otsu_threshold = 0,
      legacy_reference_otsu_threshold = 100, legacy_sample_n = 2000L,
      legacy_reuse_existing = TRUE, legacy_write_output = TRUE, unb_verbose = TRUE)
  )
  cfg$tool_paths$python_exe             <- python_exe
  cfg$tool_paths$gdal_polygonize_script <- gdal_polygonize_script
  cfg$tool_paths$gdalwarp_path          <- gdalwarp_path
  cfg$tool_paths$ogr2ogr_exe            <- ogr2ogr_exe
  cfg
}

# ---- EFFIS (dedicated step; memory-safe; STANDARD memfrac+todisk strategy) ---
# This is the long-run EFFIS job logic. It mirrors the manual script's
# STANDARD_90M EFFIS mode: terraOptions(todisk=TRUE, modest memfrac, dedicated
# tempdir) + the SAME validate_fire_maps contract. The metric definitions are
# identical regardless of resolution/memory strategy.
run_long_effis <- function(final_map_full, out_val_dir, log,
                           memfrac = 0.4) {
  dir.create(out_val_dir, recursive = TRUE, showWarnings = FALSE)
  TTMP <- file.path(out_val_dir, "terra_tmp"); dir.create(TTMP, recursive = TRUE, showWarnings = FALSE)
  terra::terraOptions(todisk = TRUE, memfrac = memfrac, tempdir = TTMP, progress = 0)
  log("[effis] terraOptions todisk=TRUE memfrac=", memfrac, " tempdir=", TTMP)
  fm <- final_map_full
  p  <- suppressWarnings(as.numeric(sf::st_drop_geometry(fm)$p_burned_model))
  keep <- fm[is.finite(p) & p >= OPERATING_THRESHOLD, , drop = FALSE]
  keep <- sf::st_zm(keep, drop = TRUE, what = "ZM"); keep <- sf::st_make_valid(keep)
  keep <- keep[!sf::st_is_empty(keep), , drop = FALSE]
  log("[effis] polygons kept at thr ", OPERATING_THRESHOLD, " = ", nrow(keep))
  if (nrow(keep) == 0L) { log("[effis] 0 polygons >= thr; EFFIS skipped."); return(NULL) }
  keep_path <- file.path(out_val_dir, "thresholded_burned.gpkg")
  if (file.exists(keep_path)) unlink(keep_path, force = TRUE)
  sf::st_write(keep, keep_path, layer = "thresholded_burned", quiet = TRUE)
  mask_sf <- sf::st_make_valid(sf::st_zm(sf::st_read(validation_mask_shapefile, quiet = TRUE),
                                         drop = TRUE, what = "ZM"))
  mask_sf <- mask_sf[!sf::st_is_empty(mask_sf), 0, drop = FALSE]; mask_sf$mask_id <- seq_len(nrow(mask_sf))
  mask_path <- file.path(out_val_dir, "validation_mask_geometry_only.gpkg")
  if (file.exists(mask_path)) unlink(mask_path, force = TRUE)
  sf::st_write(mask_sf, mask_path, layer = "mask", quiet = TRUE)
  has_strata <- file.exists(strata_raster) && file.exists(strata_lut)
  effis <- validate_fire_maps(
    input_shapefile = keep_path, ref_shapefile = reference_burned_map,
    mask_shapefile = mask_path,
    burnable_raster = burnable_mask_run,   # 3035 90m aligned (memory-safe)
    year_target = target_year, validation_dir = out_val_dir,
    force_reprocess_ref = TRUE, force_reprocess_pred = TRUE,
    metrics_type = "all", dissolve_ref_by = "id", dissolve_input_by = NULL,
    strata_raster = if (has_strata) strata_raster else NULL,
    strata_lut    = if (has_strata) strata_lut    else NULL,
    observability_raster = change_index,
    ref_end_doy_col = "end_doy", ref_start_doy_col = "start_doy")
  # safe cleanup of the terra tempdir
  tryCatch(unlink(TTMP, recursive = TRUE, force = TRUE), error = function(e) NULL)
  effis
}

invisible(TRUE)
