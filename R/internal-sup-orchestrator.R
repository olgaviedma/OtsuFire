# ============================================================
# SUPERVISED MASTER PIPELINE - MULTIYEAR + MULTISCENARIO
# Reads deterministic decisions from 05_DECISIONS/internal_decisions.gpkg
# keep   -> burned candidates for supervised training
# review -> candidate polygons for supervised scoring
# drop   -> excluded from positives, used for hard negatives, and also scored later
# Folds, features, OOF and final mapping are run per scenario
# ============================================================

# -------------------------------------------------------------------------
# Current supervised design
# - Reads deterministic decisions from:
#   DETERMINISTIC/05_DECISIONS/internal_decisions.gpkg
# - Interprets deterministic labels as:
#   * keep   -> burned candidates for supervised training
#   * review -> candidate polygons for later supervised scoring
#   * drop   -> excluded from positives; source of hard negatives and also scored
#
# Main stages
# - A1-A3: read deterministic decisions and build burned/review pools
# - A4:    build unburned pool for training negatives
# - B1-B3: spatial folds + feature extraction
# - C1-C2: OOF diagnostics
# - D1-D2: final model + final map
#
# Note
# This script already uses the new deterministic decisions as input.
# Negative-class generation uses a single policy (2026-06-05):
# - all_sources: deterministic drops + random background + Otsu current-year patches
# -------------------------------------------------------------------------

# ============================== 1) PACKAGES ==============================
# Block 7: top-level library() calls removed. All these packages are
# now declared in DESCRIPTION Imports and resolved via the package
# namespace at load time.
#
# Gate 1C.4 (2026-06-08) cfg-isolation: the former top-level
# `sf::sf_use_s2(FALSE)` was REMOVED. At namespace load it mutated the user's
# GLOBAL sf S2 setting permanently with no restore (a session-state leak), and
# in installed-package mode it was DEAD anyway (top-level R/*.R only re-runs
# under pkgload::load_all()). The S2 toggle is now scoped + restored INSIDE
# run_supervised_pipeline() via on.exit (mirroring the AS07 fix already applied
# in internal-sup-otsu-negative.R / validate-fire-maps.R), so the supervised
# run does not contaminate the caller's session and two runs are isolated.

# ============================== 2) GLOBAL CONFIG =========================
# Block 4B: project-root resolution is dispatcher-injected at runtime.
find_project_paths_file <- function(start = getwd()) NULL
project_paths <- list()
data_base <- NULL
# BUG 3 Phase 1b (2026-06-05): removed the `script_base` and
# `supervised_functions_dir` placeholders. Both were always NULL: `script_base`
# was dead plumbing into verify_otsu_negative_helpers() (which ignores it),
# and `supervised_functions_dir` was never read. The matching dispatcher
# bindings were dropped too.
composite_base <- NULL
result_name <- NULL
# 2026-06-05: SUPERVISED output base is injected by the dispatcher from the
# config's output_routes (output_dir + run_name). When NULL (e.g. the
# standalone-script path), result_dir falls back to the data_base/result_name
# reconstruction below.
supervised_output_base <- NULL

# ---- Years / scenarios to run
years_to_run     <- c(2012)   # 
scenarios_to_run <- c("balanced")

# ---- Prefixes base
prefix_oof_base <- "patch"
prefix_base     <- "patch_certified"

# ---- Toggles
DO_POOLS    <- TRUE
DO_FOLDS    <- TRUE
DO_FEATURES <- TRUE
DO_MODEL    <- TRUE
# §N+27 RESOLVED (2026-06-05): the abandoned supervised burned-like
# "RBR_TRAINING_REGISTRY" research line was removed. The constants
# LEGACY_RBR_TRAINING_REGISTRY_PATH / RBR_TRAINING_REGISTRY_PATH /
# REGISTRY_PROB_THR / REGISTRY_TOP_FRACTION, the append_supervised_registry()
# helper, and STEP C4 are gone. The deterministic phase keeps its own
# registry READER (build_rbr_keep_pool_from_registry() in
# internal-det-canonical-decisions.R) — unrelated to this supervised writer.

# ------------------------------ CURRENT-YEAR TEMPORAL FILTER -------------
CURRENTYEAR_PREYEAR_OVERLAP_THR <- 0.70
CURRENTYEAR_HOTSPOT_DENSITY_THR <- 0.001
CURRENTYEAR_TEMPORAL_PENALTY_FLOOR <- 0.10

# ------------------------------ UNBURNED OPTIONS -------------------------
# Negative-pool strategy - FIXED to a single policy (§N+26, 2026-06-05):
#
#   "all_sources"  - Full operational mode (used in MASTER pipeline).
#       Combines the TWO current-year unburned sources into the training pool
#       (GATE 6.5, 2026-06-12: the deterministic-drop "contextual" bucket was
#       removed — a deterministic drop does NOT automatically become an unburned
#       label):
#         (1) Random burnable-background cells -> random_background_sampled
#         (2) Otsu current-year unburned patches -> otsu_unburned_sampled
#       Each source routes to its own training_group in train_final_model_direct()
#       with an independent cap. No historical/cross-year content is used.
#       See HANDOFF sections 23-27 for taxonomy and composition rules.
#
# §N+26 (2026-06-05): the `UNB_STRATEGY` binding was REMOVED. The pipeline is
#   always `all_sources`; the user-selectable `negative_pool_policy` knob, the
#   deprecated `otsu_unburned_generation` alias and the retired
#   `deterministic_direct` mode are gone, the dispatcher no longer supplies a
#   strategy binding, and the only former reader (the `common_unb_dir` output
#   selection) is now hardwired to the all_sources path. (Unrelated to the D4a
#   empty-Otsu-pool robustness fallback in internal-sup-otsu-negative.R, which
#   degrades all_sources to deterministic_direct *semantics* and is kept.)
# B7 (2026-06-06): the `MIN_BURNED_POOL_N <- 5L` placeholder was removed. The
# live F7 sparse-year guard (years with fewer keep patches abort cleanly) reads
# `config$min_burned_pool_n` directly in build_supervised_training_pools()
# (supervised-pools.R); this uppercase orchestrator binding was never read, and
# the dispatcher no longer injects it.

UNB_VERBOSE  <- TRUE
# GATE 6.4 (2026-06-11): `UNB_OTSU_NEG_CODE_DIR` placeholder removed (Otsu
# negative helpers are in-package; no live consumer).

# Shared deterministic-decision negative parameters (det drops + random background;
# used by all_sources' build_unburned_from_deterministic_decisions() call).
UNB_EXCL_BUFFER_M   <- 500
UNB_RANDOM_RBR_Q    <- 0.50
UNB_N_RANDOM_CELLS  <- 1500
UNB_RANDOM_SEED     <- 42
UNB_RANDOM_PATCH_SIZE_CELLS <- 3

# Otsu pipeline parameters (used by all_sources mode)
UNB_OTSU_NEG_MODE <- "burnable_only"
UNB_OTSU_NEG_THRESHOLD <- 0
UNB_OTSU_NEG_REFERENCE_THRESHOLD <- 100
UNB_OTSU_NEG_MIN_THRESHOLD_VALUE <- 0
UNB_OTSU_NEG_MIN_PIXELS <- 8
UNB_OTSU_NEG_BUFFERS_M <- 90
UNB_OTSU_NEG_CORE_THR <- 0.60
UNB_OTSU_NEG_ALPHA_BOOST <- 0.25
UNB_OTSU_NEG_MIN_BASE_BOOST <- 0.35
UNB_OTSU_NEG_DIST_POWER <- 1
UNB_OTSU_NEG_KEEP_HI <- 0.45
UNB_OTSU_NEG_DROP_LO <- 0.15
UNB_OTSU_NEG_USE_DROP <- TRUE
# GATE 6.4 (2026-06-11): the dead `UNB_OTSU_NEG_USE_REVIEW` /
# `UNB_OTSU_NEG_USE_KEEP` / `UNB_OTSU_NEG_REVIEW_MAX_S_PATCH` /
# `UNB_OTSU_NEG_KEEP_MAX_S_PATCH` placeholders were removed (Otsu review/keep are
# never negatives — they were already excluded by train_final_model_direct()).
# Only the `drop` path remains.
UNB_OTSU_NEG_DROP_MAX_S_PATCH <- 0.15
# GATE 6.2 (2026-06-11): `UNB_OTSU_NEG_SAMPLE_N` / `UNB_OTSU_NEG_SAMPLE_PROPS` /
# `UNB_OTSU_NEG_RANDOM_SEED` placeholders removed with the Otsu generation-side
# pre-thinning. cap_otsu (otsu_unburned_to_burned_ratio) is the sole Otsu
# selector; the full valid Otsu drop pool flows on.
UNB_OTSU_NEG_EXCL_BUFFER_M <- 0
UNB_OTSU_NEG_MIN_AREA_HA <- 0
UNB_OTSU_NEG_REUSE_EXISTING <- TRUE
UNB_OTSU_NEG_WRITE_OUTPUT <- TRUE

# ======================================================================
# USER CONFIGURATION - fill in these paths for your machine before running.
# These are external tool binaries that are not bundled with the package.
# ======================================================================

# Python executable (used by GDAL scripts)
# Example: "C:/ProgramData/Anaconda3/python.exe"
python_exe <- NULL  # injected by dispatcher

# gdalwarp executable
# Example: "C:/ProgramData/Anaconda3/Library/bin/gdalwarp.exe"
gdalwarp_path <- NULL  # injected by dispatcher

# gdal_polygonize.py script (distributed with GDAL)
# Example: "C:/ProgramData/Anaconda3/Scripts/gdal_polygonize.py"
gdal_polygonize_script <- NULL  # injected by dispatcher

# ogr2ogr executable
# Example: "C:/ProgramData/Anaconda3/Library/bin/ogr2ogr.exe"
ogr2ogr_exe <- NULL  # injected by dispatcher

# ArcGIS Python (optional - leave NULL if not using ArcGIS)
arcgis_python <- NULL  # injected by dispatcher
# arcgis_python <- "C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/python.exe"

# ======================================================================
# PATH VALIDATION - do not edit below this line.
# ======================================================================



























# ============================== 3) SOURCE FUNCTIONS ======================
# Block 4B: every supervised function file has been migrated into the
# package as internal-sup-*.R, so no external sourcing is needed here.

# active functions
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)

# active wrappers
# (migrated: function now in-package)
# (migrated: function now in-package)
# (migrated: function now in-package)

# ============================== 4) COMMON HELPERS ========================
get_corine_year <- function(y) {
  if (y >= 1984 && y <= 1999) "1990"
  else if (y <= 2005) "2000"
  else if (y <= 2011) "2006"
  else if (y <= 2017) "2012"
  else "2018"
}

reclass_matrix <- matrix(c(
  1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,
  12,2,13,2,14,2,15,2,16,2,17,2,18,2,19,2,20,2,21,2,
  22,3,23,4,24,5,25,6,
  26,7,31,7,32,7,
  27,8,28,8,29,8,
  30,9,33,10,34,11,
  35,12,36,12,37,12,38,12,39,12,40,12,41,12,42,12,43,12,44,12
), ncol = 2, byrow = TRUE)

cor_groups <- list(
  open       = c(30, 31, 32, 33, 34),
  wetlands   = c(35, 36, 37, 38, 39),
  water      = c(40, 41, 42, 43, 44),
  herbaceous = c(26, 18),
  urban      = c(1,2,3,4,5,6,7,8,9,10,11),
  agri       = c(12,13,14,15,16,17,19,20,21,22),
  forest     = c(23,24,25,27,28,29)
)

align_to_template <- function(r, template, method = c("near", "bilinear"), name = deparse(substitute(r)), verbose = TRUE) {
  method <- match.arg(method)
  
  if (!inherits(r, "SpatRaster")) stop(sprintf("%s is not a SpatRaster", name))
  if (!inherits(template, "SpatRaster")) stop("template is not a SpatRaster")
  
  msg0 <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # 1) CRS
  same_crs <- FALSE
  try({
    same_crs <- terra::same.crs(r, template)
  }, silent = TRUE)
  
  if (!isTRUE(same_crs)) {
    msg0("Reprojecting raster '%s' to template CRS (method=%s)...", name, method)
    r <- terra::project(r, template, method = method)
  }
  
  # 2) Grid geometry
  same_geom <- FALSE
  try({
    same_geom <- terra::compareGeom(r, template, stopOnError = FALSE)
  }, silent = TRUE)
  
  if (!isTRUE(same_geom)) {
    msg0("Grid mismatch detected for '%s' -> resampling to template (method=%s).", name, method)
    msg0("  '%s' res=(%s,%s) origin=(%s,%s)", name,
         terra::res(r)[1], terra::res(r)[2],
         terra::origin(r)[1], terra::origin(r)[2])
    msg0("  template res=(%s,%s) origin=(%s,%s)",
         terra::res(template)[1], terra::res(template)[2],
         terra::origin(template)[1], terra::origin(template)[2])
    
    r <- terra::resample(r, template, method = method)
  }
  
  # 3) Final hard check
  ok_final <- terra::compareGeom(r, template, stopOnError = FALSE)
  if (!isTRUE(ok_final)) {
    stop(sprintf("Raster '%s' still does not match template after alignment.", name))
  }
  
  r
}

# ============================== 5) MAIN FUNCTION =========================
#
# OtsuFire 0.5.0 (2026-05-09): the supervised orchestrator now exposes
# the four sampling caps, the new `feature_whitelist_override` and
# `feature_weights` hooks for OOF/FINAL parity (replacing the removed
# `additional_drop_cols` deny-list), and a `reuse_upstream` execution
# mode for downstream-only reruns. Defaults reproduce the historical
# behaviour byte-for-byte. See `_RELEASE_0_5_0_REPORT.md` for details.
run_supervised_pipeline <- function(target_year, scenario,
                                    min_burned_pool_n = 5L,
                                    overwrite                              = TRUE,
                                    # Gate 1B (2026-06-07): these methodological
                                    # params are REQUIRED resolved args with NO
                                    # defaults. The single source of truth is
                                    # cfg$train_control / cfg$model_params,
                                    # resolved by build_supervised_burned_config()
                                    # and threaded by the dispatcher
                                    # .of_run_supervised_oneyear(). A dropped arg
                                    # ERRORS in the guard block below.
                                    random_to_burned_ratio,
                                    otsu_unburned_to_burned_ratio,
                                    feature_whitelist_override,
                                    feature_weights,
                                    reuse_upstream                         = FALSE,
                                    oof_nrounds_max,
                                    oof_early_stop,
                                    oof_seed_base,
                                    final_sampling_seed,
                                    final_seed,
                                    final_val_frac,
                                    final_group_col,
                                    final_nrounds_max,
                                    final_early_stopping_rounds,
                                    final_impute_numeric,
                                    final_impute_factor_missing,
                                    internal_decisions_path                = NULL,
                                    engine_bindings                        = NULL,
                                    config                                 = NULL) {
  # Gate 1B: required-arg guard (no silent methodological defaults).
  .req <- c("random_to_burned_ratio", "otsu_unburned_to_burned_ratio",
            "feature_whitelist_override", "feature_weights",
            "oof_nrounds_max", "oof_early_stop", "oof_seed_base",
            "final_sampling_seed", "final_seed", "final_val_frac",
            "final_group_col", "final_nrounds_max", "final_early_stopping_rounds",
            "final_impute_numeric", "final_impute_factor_missing")
  for (.nm in .req) {
    if (eval(call("missing", as.name(.nm)))) {
      stop("run_supervised_pipeline(): required resolved arg '", .nm,
           "' is missing (no methodological default; threaded from ",
           "cfg$train_control / cfg$model_params via the dispatcher).",
           call. = FALSE)
    }
  }

  # BUG 3 Phase 1b (2026-06-05): explicit engine bindings.
  # The dispatcher previously re-environmented this function so its ~40
  # config-derived bare-name lookups (data_base, composite_base, UNB_*,
  # REGISTRY_*, tool paths, supervised_output_base, prefix_*, ...) resolved
  # against an injected env. That `environment(pipeline) <- eng_env` magic is
  # gone; the dispatcher now passes the SAME named list as `engine_bindings`
  # and we bind every entry as a LOCAL in this execution frame. This MUST be
  # the first executable statement, before ANY use of an injected name, so
  # all downstream lexical lookups see the injected values. Binding the whole
  # list (rather than enumerating names) guarantees completeness. The explicit
  # params above are NOT among the binding keys, so list2env cannot clobber
  # them. When NULL (standalone-script path) the top-of-file NULL placeholders
  # apply, exactly as before.
  if (!is.null(engine_bindings)) list2env(engine_bindings, envir = environment())

  # 2026-06-05: `overwrite` is now honored end-to-end (previously the
  # orchestrator hardcoded TRUE at every write site, ignoring the user's
  # flag). overwrite=TRUE reproduces the historical behaviour (regenerate /
  # clobber stage outputs). overwrite=FALSE forwards to the unburned builder
  # and the OOF / final-model wrappers, which skip/reuse existing on-disk
  # stage outputs instead of clobbering them. Default TRUE preserves the
  # standalone-script behaviour.
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be a single TRUE or FALSE.", call. = FALSE)
  }

  # Gate 1C.4 cfg-isolation: scope the sf S2 toggle to THIS run and restore the
  # caller's prior value on exit (replaces the removed top-level
  # sf::sf_use_s2(FALSE) leak). The supervised geometry ops assume planar S2-off
  # semantics; we set it here and put it back exactly as we found it, so a run
  # leaves no global session residue for the next run / the user's session.
  .prev_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(.prev_s2)), add = TRUE)

  cat("\n============================================================\n")
  cat(sprintf("RUNNING SUPERVISED PIPELINE | year=%s | scenario=%s\n", target_year, scenario))
  cat("============================================================\n")

  # ------------------------------ ROOT PATHS -----------------------------
  # OUTPUT base. The config's output_routes (output_dir + run_name) is the
  # single source of truth for WHERE supervised outputs are written, injected
  # by the dispatcher as `supervised_output_base`. Fall back to the
  # data_base/result_name reconstruction only when it is not supplied (the
  # standalone-script path). data_base + result_name remain the authority for
  # INPUT imagery/composite/decision paths below.
  result_dir <- if (!is.null(supervised_output_base) &&
                    is.character(supervised_output_base) &&
                    length(supervised_output_base) == 1L &&
                    nzchar(supervised_output_base)) {
    supervised_output_base
  } else {
    file.path(
      data_base, "Results", target_year, result_name, "SUPERVISED", scenario
    )
  }

  deterministic_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", scenario
  )

  # The deterministic decisions GPKG comes ONLY from the user-supplied
  # internal_decisions_path threaded down from
  # build_supervised_burned_config(). It is the single source of truth; there
  # is NO convention-based reconstruction. Fail fast if it was not provided.
  if (is.null(internal_decisions_path) || !nzchar(internal_decisions_path)) {
    stop(
      "run_supervised_pipeline() requires 'internal_decisions_path'. ",
      "Pass the deterministic decisions .gpkg path via ",
      "build_supervised_burned_config(internal_decisions=...). ",
      "There is no convention-based fallback.",
      call. = FALSE
    )
  }
  internal_decisions_gpkg <- internal_decisions_path
  # The supervised run still writes its own 05_DECISIONS outputs under the
  # SUPERVISED result tree; this directory is created/used for outputs only,
  # NOT to locate the deterministic decisions input.
  deterministic_decisions_dir <- file.path(deterministic_dir, "05_DECISIONS")

  # Carpeta comon por a+/-o para reutilizar OTSU + polygonize
  deterministic_common_unb_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", "_COMMON_UNBURNED"
  )
  otsu_negative_root_dir <- file.path(result_dir, "_OTSU_NEGATIVE")

  # Salida final unburned por escenario
  unb_out_gpkg <- file.path(
    deterministic_dir, "UNBURNED",
    sprintf("%d_%s_unburned.gpkg", target_year, scenario)
  )
  unb_out_path_actual <- unb_out_gpkg
  
  dir.create(deterministic_common_unb_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(unb_out_gpkg), recursive = TRUE, showWarnings = FALSE)
  dir.create(otsu_negative_root_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Prefixes
  prefix_oof <- sprintf("%d_%s_%s", target_year, scenario, prefix_oof_base)
  prefix     <- sprintf("%d_%s_%s", target_year, scenario, prefix_base)
  
  # ------------------------------ DIRS -----------------------------------
  dirs <- list(
    `01_POOLS`          = file.path(result_dir, "01_POOLS"),
    `02_FOLDS`          = file.path(result_dir, "02_FOLDS"),
    `03_FEATURES`       = file.path(result_dir, "03_FEATURES"),
    `04_MATRIX`         = file.path(result_dir, "04_MATRIX"),
    `05_OOF`            = file.path(result_dir, "05_OOF"),
    `07_FINAL_MODEL_V2` = file.path(result_dir, "07_FINAL_MODEL_V2"),
    `08_SCORED`         = file.path(result_dir, "08_SCORED"),
    `09_FINAL_MAP`      = file.path(result_dir, "09_FINAL_MAP"),
    `0_DETERMINISTIC`   = deterministic_decisions_dir,
    `99_LOGS_EMPTY`     = file.path(result_dir, "99_LOGS_EMPTY")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  
  # ------------------------------ SAFETY KIT -----------------------------
  sfkit <- make_sf_safety_kit(
    dirs        = dirs,
    result_dir  = result_dir,
    target_year = target_year
  )
  
  # helpers del kit
  msg                  <- sfkit$msg
  `%||%`               <- sfkit$`%||%`
  announce_start       <- sfkit$announce_start
  announce_end         <- sfkit$announce_end
  
  sanitize_polygons    <- sfkit$sanitize_polygons
  ensure_area_ha       <- sfkit$ensure_area_ha
  drop_empty_sf        <- sfkit$drop_empty
  check_sf             <- sfkit$check_sf
  check_sf_if_nonempty <- sfkit$check_sf_if_nonempty
  empty_sf             <- sfkit$empty_sf
  to_crs_safe          <- sfkit$to_crs_safe
  
  safe_read_gpkg       <- sfkit$safe_read_gpkg
  safe_write_gpkg      <- sfkit$safe_write_gpkg
  read_layer_if_exists <- sfkit$read_layer_if_exists
  
  time_step            <- sfkit$time_step
  save_timelog         <- sfkit$save_timelog
  timing_csv           <- sfkit$timing_csv
  
  announce_start(sprintf("SUPERVISED MASTER pipeline | year=%s | scenario=%s", target_year, scenario))
  on.exit({
    save_timelog()
    announce_end(sprintf("SUPERVISED MASTER pipeline | year=%s | scenario=%s", target_year, scenario))
  }, add = TRUE)

  if (isTRUE(reuse_upstream)) {
    DO_POOLS <- FALSE
    DO_FOLDS <- FALSE
    DO_FEATURES <- FALSE
    DO_MODEL <- TRUE
    msg("[reuse_upstream=TRUE] Skipping STEP A, B1-B2, B3. Consuming existing pools/folds/features from disk.")
  }

  gpkg_pools_out <- file.path(
    dirs$`01_POOLS`,
    sprintf("%d_%s_pools.gpkg", target_year, scenario)
  )
  features_gpkg <- file.path(dirs$`03_FEATURES`, "features_geometry.gpkg")
  train_with_folds_reuse_gpkg <- NULL
  blocks_reuse_gpkg <- NULL

  if (isTRUE(reuse_upstream)) {
    train_with_folds_matches <- list.files(
      dirs$`02_FOLDS`,
      pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$", target_year),
      full.names = TRUE
    )
    blocks_matches <- list.files(
      dirs$`02_FOLDS`,
      pattern = sprintf("^%d_blocks_\\d+m\\.gpkg$", target_year),
      full.names = TRUE
    )
    train_with_folds_reuse_gpkg <- if (length(train_with_folds_matches) > 0) {
      train_with_folds_matches[1]
    } else {
      NULL
    }
    blocks_reuse_gpkg <- if (length(blocks_matches) > 0) {
      blocks_matches[1]
    } else {
      NULL
    }

    missing_inputs <- character(0)
    if (!file.exists(gpkg_pools_out)) {
      missing_inputs <- c(missing_inputs, gpkg_pools_out)
    }
    if (is.null(train_with_folds_reuse_gpkg) || !file.exists(train_with_folds_reuse_gpkg)) {
      missing_inputs <- c(
        missing_inputs,
        file.path(dirs$`02_FOLDS`, sprintf("%d_train_with_folds_<size>m.gpkg", target_year))
      )
    }
    if (!file.exists(features_gpkg)) {
      missing_inputs <- c(missing_inputs, features_gpkg)
    }
    if (length(missing_inputs) > 0) {
      stop(
        paste0(
          "[reuse_upstream=TRUE] Missing upstream artefact(s). Run a Baseline first. Missing: ",
          paste(missing_inputs, collapse = "; ")
        ),
        call. = FALSE
      )
    }
  }

  # §N+26 (2026-06-05): the deprecated-strategy normalization block was removed.
  # The pipeline is always `all_sources`; there is no longer an
  # `otsu_unburned_generation` alias (or any `UNB_STRATEGY` binding) to redirect.

  # ============================== 5b) COMMON INPUTS ======================
  # §N+25 (2026-06-05): the supervised-RUN inputs are now wired to the config
  # (build_supervised_burned_config). The orchestrator reads each RUN input's
  # path from `config$inputs` when supplied, else falls back to EXACTLY the
  # historical convention path. The cfg defaults are built from the SAME
  # centralized helper (.of_supervised_convention_paths) the config builder
  # uses, so cfg value and convention fallback are byte-identical. A local
  # accessor returns the on-disk path of a cfg input spec, or NULL when the
  # input is NULL / in-memory (in-memory specs are not produced for these
  # raster/vector inputs in the run path, but NULL-safety is kept).
  # Gate 1B: route every cfg-input lookup through the single shared accessor
  # `.of_sup_input_path()` (supervised-config.R) so the orchestrator, the pool
  # builder and the feature extractor resolve inputs identically and cfg$inputs
  # is the one source of truth.
  .cfg_input_path <- function(name) .of_sup_input_path(config, name)

  corine_year         <- get_corine_year(target_year)

  # ---- RUN inputs: cfg value %||% convention ----------------------------
  # peninsula_shapefile is a RUN input only for the CORINE Otsu modes inside
  # the Otsu negative builder; the orchestrator itself never reads it. It is
  # threaded to the unburned builders via supervised-pools.R (config$inputs),
  # so it is intentionally NOT reconstructed or validated here.
  topo_path <- .cfg_input_path("topo") %||%
    file.path(data_base, "Topography", "elevation_slope.tif")
  # hotspots is an OPTIONAL input wired to cfg$inputs$hotspots (cfg validates
  # it allow_null = TRUE so pre-MODIS years can run with hotspots = NULL). The
  # cfg path is the single source of truth; the convention path is only the
  # fallback for a cfg that did not carry an explicit hotspots path (e.g. the
  # standalone-script path with a data_base set). When neither the cfg nor a
  # data_base supplies one, hotspots_path is NULL and the layer is treated as
  # absent below (use_hotspots disabled), exactly as for a missing file.
  hotspots_path <- .cfg_input_path("hotspots") %||%
    (if (!is.null(data_base) && nzchar(data_base))
       file.path(data_base, "Hotspots",
                 paste0("hotspots_iberia_", target_year, ".geojson"))
     else NULL)
  corine_raster_path  <- .cfg_input_path("corine_raster") %||%
    file.path(data_base, "Corine_Masks", paste0("CLC_", corine_year, "_peninsula.tif"))

  # change_index: the MAIN annual change-index (summer RBR + DOY_post) raster.
  # It is a REQUIRED, validated field of cfg$inputs (build_supervised_burned_config
  # rejects a NULL change_index). It is CONSUMED from cfg$inputs here — NOT
  # reconstructed by the historical MinMin_<year>_mosaic_res90m.tif filename
  # convention. The convention is retained ONLY as a fallback for the
  # standalone-script path (config == NULL), so behaviour is byte-identical for
  # existing scripts while the packaged pipeline is single-source-of-truth.
  one_year_tif <- .cfg_input_path("change_index") %||% {
    cand <- file.path(
      composite_base, result_name,
      paste0("MinMin_", target_year, "_mosaic_res90m.tif")
    )
    if (!file.exists(cand)) {
      cand <- file.path(
        composite_base, "Min_Min",
        paste0("MinMin_", target_year, "_mosaic_res90m.tif")
      )
    }
    cand
  }

  rbr_aw_tif <- .cfg_input_path("delayed_change_index") %||%
    file.path(
      composite_base, "Autumn",
      paste0("mean_mean_", target_year, "_mosaic.tif")
    )

  # ------------------------------ CHECKS ---------------------------------
  # §N+25: the VALIDATION-only inputs (strata_tif / strata_lut / mask_tif /
  # mask_shp / ref_tif / ref_shp EFFIS) were previously CONSTRUCTED and
  # stopifnot-validated here but never CONSUMED by the supervised run. They
  # belong to the separate validate_fire_maps() call (validation is invoked
  # separately, exactly like the deterministic phase), so they are removed
  # from the run start: not built, not validated. The deterministic-vs-
  # supervised CONSISTENCY check (dispatcher) uses internal_decisions +
  # final_map, NOT EFFIS, and is unaffected.
  #
  # The deterministic decisions input is the user-supplied path; assert it
  # exists directly. We no longer assert dir.exists(deterministic_dir) /
  # deterministic_decisions_dir as input gates: the decisions may live
  # anywhere the user wrote them, not necessarily under the Min_Min tree.
  if (!file.exists(internal_decisions_gpkg)) {
    stop(sprintf(
      "internal_decisions GPKG not found: %s (user-supplied path).",
      internal_decisions_gpkg
    ), call. = FALSE)
  }
  stopifnot(file.exists(one_year_tif))
  stopifnot(file.exists(corine_raster_path))
  stopifnot(file.exists(topo_path))
  stopifnot(file.exists(rbr_aw_tif))

  # ------------------------------ GATE 1C.4 SEMANTIC/SPATIAL FAIL-FAST ----
  # Single canonical semantic + spatial validator. Runs ONCE here, AFTER the
  # cfg has been resolved and the REQUIRED inputs confirmed to exist (the Gate
  # 1B / C3 path layer above), but BEFORE any heavy compute (the STATIC OBJECTS
  # raster reads / alignment below and the entire pools->folds->features->model
  # chain). It fails fast (no artefact written) on CRS/overlap/empty-raster/
  # zero-burnable-mask/wrong-year/missing-column/incompatible-type/incompatible-
  # schema/contradictory-cfg/foreign-cache problems, naming the offending input
  # + cfg field. The checks are metadata/geometry-level (rasters opened
  # header-only via terra::rast); only the burnable-mask alignment touches cells,
  # reusing the Gate 1C.1 helper. cfg-isolated: reads ONLY config + the on-disk
  # inputs it points at, writes nothing. Skipped on the standalone-script path
  # (config == NULL), which has no cfg to validate against.
  if (!is.null(config)) {
    # strict = TRUE: the SAME engine that builds the structured report, run in
    # fail-fast mode here. A blocking-severity check FAIL aborts the run with a
    # single aggregated error before any heavy compute (no duplicated checks).
    validate_supervised_execution(
      config                     = config,
      strict                     = TRUE,
      target_year                = target_year,
      reuse_upstream             = reuse_upstream,
      feature_whitelist_override = feature_whitelist_override,
      data_base                  = data_base,
      composite_base             = composite_base,
      result_name                = result_name
    )
  }

  # ------------------------------ STATIC OBJECTS -------------------------
  topo       <- terra::rast(topo_path)
  
  rbr_stack  <- terra::rast(one_year_tif)   # band1 = RBR summer, band2 = DOY_post
  rbr_summer <- rbr_stack[[1]]
  doy_post   <- rbr_stack[[2]]
  
  # template maestro del a+/-o
  template_r <- rbr_summer
  
  # Alinear doy_post al propio template por seguridad
  doy_post <- align_to_template(
    r = doy_post,
    template = template_r,
    method = "near",
    name = "doy_post",
    verbose = TRUE
  )
  
  dem_r <- align_to_template(
    r = topo[[1]],
    template = template_r,
    method = "bilinear",
    name = "dem",
    verbose = TRUE
  )
  
  slope_r <- align_to_template(
    r = topo[[2]],
    template = template_r,
    method = "bilinear",
    name = "slope",
    verbose = TRUE
  )
  
  corine_r <- align_to_template(
    r = terra::rast(corine_raster_path),
    template = template_r,
    method = "near",
    name = "corine_r",
    verbose = TRUE
  )
  
  rbr_aw <- align_to_template(
    r = terra::rast(rbr_aw_tif)[[1]],
    template = template_r,
    method = "bilinear",
    name = "rbr_aw",
    verbose = TRUE
  )
  
  # ------------------------------ VECTORS --------------------------------
  # No ecoregion I/O in the supervised phase: ecoregions belong to the
  # deterministic delineation stage and were removed from the supervised
  # feature engine on 2026-06-05.

  hotspots_sf_base <- if (!is.null(hotspots_path) && file.exists(hotspots_path)) {
    sf::read_sf(hotspots_path)
  } else {
    msg("WARNING: no existe hotspots: %s (se desactiva use_hotspots)",
        hotspots_path %||% "<none: cfg$inputs$hotspots = NULL>")
    sf::st_sf(
      frp = numeric(),
      confidence = numeric(),
      year = integer(),
      month = integer(),
      geometry = sf::st_sfc(crs = sf::NA_crs_)
    )[0, ]
  }
  
  # ============================== 6) POOLS ===============================
  # BUG 3 Phase 2 (2026-06-05): STEP A (training-pool construction, A1-A6) is
  # now a single implementation living in the exported
  # build_supervised_training_pools(); the orchestrator DELEGATES to it instead
  # of inlining the A1-A6 sequence (read internal_decisions -> QA relabel ->
  # burned/review/scoring pools -> both all_sources unburned builders -> merge
  # train_labeled -> write 01_POOLS/<year>_<scenario>_pools.gpkg). Behaviour is
  # byte-identical:
  #   * config is the supervised config S3 object (threaded explicitly by the
  #     dispatcher). The stage resolves the deterministic decisions path, every
  #     UNB_* / UNB_OTSU_NEG_* parameter and the four tool paths straight from
  #     config (same defaults as .of_supervised_engine_bindings()), so the two
  #     unburned-builder calls, the random-background seed (UNB_RANDOM_SEED = 42;
  #     GATE 6.2 removed the Otsu-negative sampling seed) and the
  #     det-then-Otsu-negative call ordering are unchanged -> identical RNG
  #     stream + random negatives.
  #   * deterministic_decisions = NULL so the stage uses the canonical
  #     config$inputs$internal_decisions path (== the orchestrator's
  #     internal_decisions_gpkg, the single source of truth).
  #   * write_outputs = TRUE writes the SAME pools.gpkg layers + QA files;
  #     overwrite is forwarded to the deterministic-decisions builder.
  # The returned pools_gpkg / train_labeled rewire the downstream folds + feature
  # consumers (which previously read the in-frame gpkg_pools_out / train_labeled).
  if (isTRUE(DO_POOLS)) {
    stopifnot(inherits(config, "otsufire_supervised_burned_config"))
    res_pools <- time_step("A build_supervised_training_pools", {
      build_supervised_training_pools(
        config                  = config,
        deterministic_decisions = NULL,
        write_outputs           = TRUE,
        overwrite               = overwrite
      )
    })

    if (!is.null(res_pools$pools_gpkg)) gpkg_pools_out <- res_pools$pools_gpkg

    msg("DONE - pools saved: %s", gpkg_pools_out)
  }
  
  # ============================== 7) SPATIAL FOLDS =======================
  train_with_folds_gpkg <- NULL
  blocks_gpkg <- NULL
  if (isTRUE(reuse_upstream)) {
    train_with_folds_gpkg <- train_with_folds_reuse_gpkg
    blocks_gpkg <- blocks_reuse_gpkg
  }
  
  if (isTRUE(DO_FOLDS)) {
    
    stopifnot(exists("make_block_folds"))
    
    msg("STEP B1 - Read train_labeled from 01_POOLS")
    train_labeled_sf <- time_step("B0 Read train_labeled (for folds)", {
      safe_read_gpkg(gpkg_pools_out, "train_labeled", "train_labeled")
    })
    
    # BUG 3 Phase 2 PILOT (2026-06-05): the folds stage is now a single
    # implementation living in the exported make_spatial_folds(); the
    # orchestrator DELEGATES to it instead of inlining the make_block_folds()
    # call. Behaviour is byte-identical: the same train_labeled_sf object is
    # forwarded, the seed_base/n_repeats/acceptance-gate values are pinned to
    # the historical orchestrator values inside make_spatial_folds(), out_dir
    # is dirs$`02_FOLDS` (== config$output_routes$folds_dir), and the year /
    # filenames are unchanged. `config` is the supervised config S3 object,
    # threaded explicitly into run_supervised_pipeline() by the dispatcher
    # (see internal-supervised-dispatch.R); the engine's list2env bindings do
    # not carry it, so it is passed as a dedicated parameter.
    msg("STEP B2 - Make spatial folds (delegated to make_spatial_folds())")
    res_folds <- time_step("B1 make_block_folds", {
      stopifnot(inherits(config, "otsufire_supervised_burned_config"))
      make_spatial_folds(
        train_labelled = train_labeled_sf,
        config         = config,
        split_unit     = "fire",
        block_sizes_m  = c(5000, 3000, 2000),
        k_candidates   = c(5L, 4L, 3L),
        out_dir        = dirs$`02_FOLDS`,
        write_outputs  = TRUE,
        # 2026-06-05 (D2 BAD-hardcode fix): forward the user's overwrite flag
        # instead of a literal TRUE. This was the ONLY remaining overwrite
        # override; the other four stage delegations already pass
        # `overwrite = overwrite`. Default overwrite=TRUE -> byte-identical.
        overwrite      = overwrite
      )
    })

    if (!is.null(res_folds$train_with_folds_gpkg)) train_with_folds_gpkg <- res_folds$train_with_folds_gpkg
    if (!is.null(res_folds$blocks_gpkg))           blocks_gpkg          <- res_folds$blocks_gpkg

    bs_sel <- tryCatch(res_folds$selected$block_size_m, error = function(e) NULL)

    if (is.null(train_with_folds_gpkg)) {
      cand <- if (!is.null(bs_sel)) {
        file.path(dirs$`02_FOLDS`, sprintf("%d_train_with_folds_%dm.gpkg", target_year, bs_sel))
      } else {
        lf <- list.files(dirs$`02_FOLDS`, pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$", target_year), full.names = TRUE)
        if (length(lf)) lf[1] else NA_character_
      }
      if (is.character(cand) && length(cand) == 1 && !is.na(cand) && file.exists(cand)) train_with_folds_gpkg <- cand
    }
    
    if (is.null(blocks_gpkg)) {
      cand <- if (!is.null(bs_sel)) {
        file.path(dirs$`02_FOLDS`, sprintf("%d_blocks_%dm.gpkg", target_year, bs_sel))
      } else {
        lf <- list.files(dirs$`02_FOLDS`, pattern = sprintf("^%d_blocks_\\d+m\\.gpkg$", target_year), full.names = TRUE)
        if (length(lf)) lf[1] else NA_character_
      }
      if (is.character(cand) && length(cand) == 1 && !is.na(cand) && file.exists(cand)) blocks_gpkg <- cand
    }
    
    msg("FOLDS outputs:")
    msg(" - train_with_folds_gpkg: %s", train_with_folds_gpkg %||% "NULL (no encontrado)")
    msg(" - blocks_gpkg:           %s", blocks_gpkg %||% "NULL (no encontrado)")
  }
  
  # ============================== 8) FEATURE EXTRACTION ==================
  need_features_for_model <- isTRUE(DO_MODEL) && !file.exists(features_gpkg)
  
  if (isTRUE(DO_FEATURES) || isTRUE(need_features_for_model)) {
    if (!isTRUE(DO_FEATURES) && isTRUE(need_features_for_model)) {
      msg("STEP B3 - features_geometry.gpkg no existe pero DO_MODEL=TRUE. Se fuerzan features.")
    }

    if (is.null(train_with_folds_gpkg) || !file.exists(train_with_folds_gpkg)) {
      lf <- list.files(
        dirs$`02_FOLDS`,
        pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$", target_year),
        full.names = TRUE
      )
      if (length(lf) > 0) train_with_folds_gpkg <- lf[1]
    }

    stopifnot(!is.null(train_with_folds_gpkg), file.exists(train_with_folds_gpkg))
    stopifnot(file.exists(gpkg_pools_out))

    # BUG 3 Phase 2 (2026-06-05): STEP B3 (feature extraction) is now a single
    # implementation living in the exported extract_supervised_features(); the
    # orchestrator DELEGATES to it instead of inlining the staleness-rebuild
    # block + the extract_features() engine call. Behaviour is byte-identical:
    #   * train_with_folds / scoring_pool are passed as the resolved on-disk
    #     GPKG paths (train_with_folds_gpkg / gpkg_pools_out), read back inside
    #     the stage with the SAME safe_read_gpkg closures.
    #   * the aligned raster stack (rbr_summer/doy_post/rbr_aw/dem_r/slope_r/
    #     corine_r), the hotspot layer (hotspots_sf_base) and the CORINE group
    #     LUT (cor_groups) are passed through the dot-prefixed internal args so
    #     the stage does NOT re-load or re-align them -> identical raster inputs.
    #   * use_hotspots is passed TRUE so the stage's runtime flag reduces to the
    #     historical use_hotspots_flag = nrow(hotspots_sf) > 0 (the stage's
    #     public default FALSE would force the block off; production must keep
    #     the runtime behaviour).
    #   * out_dir is dirs$`03_FEATURES` (== config$output_routes$features_dir);
    #     the engine writes features_geometry.gpkg + the two *_features.rds with
    #     the same basenames/layer names. The staleness-rebuild logic (reuse vs
    #     rebuild) is replicated verbatim inside the stage.
    res_features <- time_step("B3 extract_supervised_features", {
      stopifnot(inherits(config, "otsufire_supervised_burned_config"))
      extract_supervised_features(
        train_with_folds = train_with_folds_gpkg,
        scoring_pool     = gpkg_pools_out,
        config           = config,
        use_hotspots     = TRUE,
        # OPTIONAL shape/size block (OtsuFire 0.12.0): compute the shape columns
        # when the cfg flag is ON so the OOF + FINAL stages can use them. SAME
        # cfg field the modeling block threads into OOF/FINAL. OFF (default) ->
        # no shape join -> byte-identical features.
        use_shape        = isTRUE(config$train_control$include_shape_features),
        out_dir          = dirs$`03_FEATURES`,
        write_outputs    = TRUE,
        .aligned_rasters = list(
          rbr_summer = rbr_summer,
          doy_post   = doy_post,
          rbr_aw     = rbr_aw,
          dem_r      = dem_r,
          slope_r    = slope_r,
          corine_r   = corine_r
        ),
        .hotspots_sf = hotspots_sf_base,
        .cor_groups  = cor_groups
      )
    })

    if (!is.null(res_features$features_geometry_gpkg)) {
      features_gpkg <- res_features$features_geometry_gpkg
    }
  }
  
  # ============================== 9) MODELING ============================
  # Gate 1E (2026-06-09): the feature-schema parity manifest is assembled inside
  # the modeling block (OOF -> FINAL -> scoring); pre-initialise it to NULL so the
  # run return can surface it whether or not modeling ran this invocation.
  schema_parity_manifest <- NULL
  # PHASE 2 (artifact_hard): consolidated training-pool layer path, set inside
  # the modeling block; pre-initialised so the run return surfaces it (or NULL).
  training_pool_layer_path <- NULL
  if (isTRUE(DO_MODEL)) {

    msg("STEP C0 - Read features_geometry.gpkg (03_FEATURES)")
    gpkg_features <- file.path(dirs$`03_FEATURES`, "features_geometry.gpkg")
    stopifnot(file.exists(gpkg_features))
    
    if (is.null(train_with_folds_gpkg) || !file.exists(train_with_folds_gpkg)) {
      lf <- list.files(
        dirs$`02_FOLDS`,
        pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$", target_year),
        full.names = TRUE
      )
      if (length(lf) > 0) train_with_folds_gpkg <- lf[1]
    }
    
    if (is.null(blocks_gpkg) || !file.exists(blocks_gpkg)) {
      lf <- list.files(
        dirs$`02_FOLDS`,
        pattern = sprintf("^%d_blocks_\\d+m\\.gpkg$", target_year),
        full.names = TRUE
      )
      if (length(lf) > 0) blocks_gpkg <- lf[1]
    }
    
    train_features <- time_step("C0 Read train_features", {
      safe_read_gpkg(gpkg_features, "train_features", "train_features")
    })
    scoring_features <- time_step("C0b Read scoring_features", {
      safe_read_gpkg(gpkg_features, "scoring_features", "scoring_features")
    })
    
    needed_folds <- c("block_id","fold_rep1","fold_rep2")
    if (!all(needed_folds %in% names(train_features))) {
      msg("STEP C0c - train_features no tiene folds; haciendo join desde train_with_folds")
      
      stopifnot(!is.null(train_with_folds_gpkg) && file.exists(train_with_folds_gpkg))
      train_with_folds <- time_step("C0c Read train_with_folds", {
        safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "train_with_folds")
      })
      
      folds_df <- sf::st_drop_geometry(train_with_folds) |>
        select(fire_uid, any_of(needed_folds))
      
      train_features <- time_step("C0d Join folds into train_features", {
        train_features |> left_join(folds_df, by = "fire_uid")
      })
    }

    # ============================================================
    # PHASE 2 (artifact_hard hard-negative mining): ADDITIVE + OFF BY DEFAULT.
    # Promote artifact_hard negatives from the deterministic scoring universe
    # into the training features AFTER feature extraction (+ fold-join) and
    # BEFORE the OOF/FINAL stages, so the promoted rows get fold assignments and
    # appear in the OOF predictions, then are trained on by the FINAL model. The
    # scoring_features layer is returned UNCHANGED (the promoted rows STAY in the
    # scoring universe too). When the feature is disabled
    # (negative_pool_params$artifact_hard$enabled = FALSE, the default), this is
    # a STRICT NO-OP and the entire modeling chain is byte-identical to today.
    # ------------------------------------------------------------
    .ah_cfg <- config$negative_pool_params$artifact_hard
    .ah_on  <- isTRUE(.ah_cfg$enabled)
    # Resolve the per-row weighting inputs ONCE here (single source = cfg). With
    # the feature OFF these are NULL / character(0) so the OOF + FINAL stages
    # resolve NO weight vector -> byte-identical.
    .ah_source_weights <- config$negative_pool_params$source_weights
    .ah_source_tag     <- if (.ah_on) "artifact_hard" else character(0)
    .ah_total_weight_ratio <-
      config$negative_pool_params$artifact_hard$total_weight_ratio %||% 0.10
    # PHASE 2 (artifact_hard): the promotion result is hoisted here (NULL when
    # the feature is OFF) so the consolidated training-pool layer can record the
    # per-row eligibility / branch / threshold after the FINAL stage.
    promo              <- NULL
    .gpkg_train_aug    <- gpkg_features
    if (.ah_on) {
      msg("STEP C0e - PHASE 2 artifact_hard promotion (enabled)")
      promo <- time_step("C0e promote_artifact_hard_negatives", {
        promote_artifact_hard_negatives(
          train_features   = train_features,
          scoring_features = scoring_features,
          config           = config,
          enable_derived   = FALSE
        )
      })
      train_aug <- promo$train_features
      msg(" - artifact_hard promoted rows: %d", promo$n_promoted)

      if (promo$n_promoted > 0L) {
        # Assign block_id / fire_uid / fold_rep to the promoted rows (mirror the
        # verified runner run_one()): unique block_id "AH_*", a unique fire_uid
        # disjoint from the existing rows, and a distributed fold per fold_rep so
        # every promoted row gets a held-out OOF prediction.
        train_aug <- sf::st_as_sf(train_aug)
        train_aug$fire_uid <- as.character(train_aug$fire_uid)
        ah_mask <- !is.na(train_aug$source) & train_aug$source == "artifact_hard"
        n_ah <- sum(ah_mask)
        train_aug$block_id[ah_mask] <- paste0("AH_", seq_len(n_ah))
        uid <- as.character(train_aug$fire_uid[ah_mask])
        bad <- is.na(uid) | uid == "" | duplicated(uid) |
          (uid %in% as.character(train_aug$fire_uid[!ah_mask]))
        uid[bad] <- paste0("AH_uid_", which(ah_mask)[bad])
        train_aug$fire_uid[ah_mask] <- uid
        # Detect the fold-rep columns DYNAMICALLY (do NOT hard-code fold_rep1/2):
        # identical when only fold_rep1/fold_rep2 exist, robust if make_spatial_folds
        # ever produces more repeats. (make_spatial_folds is NOT re-run; promoted
        # rows get synthetic folds across whatever fold_rep* columns are present.)
        ah_fold_cols <- grep("^fold_rep\\d+$", names(train_aug), value = TRUE)
        kvals <- sort(unique(unlist(lapply(
          ah_fold_cols, function(fc) as.integer(train_aug[[fc]])))))
        kvals <- kvals[is.finite(kvals)]
        nk <- length(kvals)
        if (nk > 0L) {
          set.seed(20260615)
          for (fc in ah_fold_cols) {
            assign_k <- kvals[(seq_len(n_ah) - 1L) %% nk + 1L]
            assign_k <- sample(assign_k)
            train_aug[[fc]][ah_mask] <- assign_k
          }
        }

        # Persist the augmented train_features into a SEPARATE features GPKG (the
        # scoring_features layer is copied UNCHANGED). The OOF + FINAL + scoring
        # stages below read from this augmented GPKG instead of gpkg_features, so
        # the promoted rows flow through the whole modeling chain. The original
        # gpkg_features is left untouched.
        .gpkg_train_aug <- file.path(dirs$`03_FEATURES`,
                                     "features_geometry_artifact_hard.gpkg")
        if (file.exists(.gpkg_train_aug)) unlink(.gpkg_train_aug, force = TRUE)
        sf::st_write(train_aug, .gpkg_train_aug, layer = "train_features",
                     quiet = TRUE)
        sf::st_write(sf::st_as_sf(scoring_features), .gpkg_train_aug,
                     layer = "scoring_features", quiet = TRUE, append = TRUE)
        # The in-memory train_features handed to run_oof_diagnostics below is the
        # augmented frame.
        train_features <- train_aug

        # Per-source traceability: write the promotion audit alongside the model
        # outputs.
        .ah_audit_csv <- file.path(dirs$`07_FINAL_MODEL_V2`,
                                   paste0(prefix, "_artifact_hard_promotion_audit.csv"))
        if (isTRUE(overwrite) || !file.exists(.ah_audit_csv)) {
          utils::write.csv(promo$audit, .ah_audit_csv, row.names = FALSE)
        }
      } else {
        msg(" - no rows promoted; modeling chain unchanged")
      }
    }

    # BUG 3 Phase 2 (2026-06-05): the OOF-diagnostics stage (STEP C1-C2:
    # build XGB params -> design matrix -> OOF) is now a single exported
    # function, `run_oof_diagnostics()` (R/supervised-oof.R). The
    # orchestrator DELEGATES to it. Behaviour is byte-identical: the
    # function rebuilds the SAME XGB params block inline (params = NULL),
    # forwards the SAME run_dm_oof_pipeline(...) args/values, the SAME
    # whitelist override / feature weights, and writes the SAME 05_OOF /
    # 04_MATRIX files under the SAME prefix. `config` is threaded in from
    # the dispatcher (Phase 1b). The Bug-5 note (character is canonical for
    # `class`; factor conversion happens only inside
    # build_design_matrix_patches) still holds inside run_oof_diagnostics.
    stopifnot(inherits(config, "otsufire_supervised_burned_config"))
    # OPTIONAL shape/size block (OtsuFire 0.12.0): ONE config-derived local fed
    # to BOTH run_oof_diagnostics() and train_final_burned_model() (so OOF and
    # FINAL receive the value from the SAME cfg field and cannot diverge), and
    # already fed to extract_supervised_features() above so the columns are
    # actually computed when the flag is on. OFF (FALSE) -> byte-identical.
    include_shape_features <- isTRUE(config$train_control$include_shape_features)
    # OOF cross-validates over the fold-repeat columns. Detect them DYNAMICALLY
    # from the (possibly augmented) training frame that actually enters OOF, ordered
    # by repeat number, instead of assuming exactly fold_rep1/fold_rep2. Identical
    # to the old hard-coded c("fold_rep1","fold_rep2") when only those two exist
    # (the fold-join above guarantees train_features carries the fold columns).
    fold_cols <- .of_fold_rep_cols(train_features)
    if (length(fold_cols) == 0L) {
      stop("run_supervised_pipeline(): no fold_rep* columns found in the training ",
           "frame fed to OOF (expected at least fold_rep1 -- did make_spatial_folds() ",
           "run and the fold-join succeed?).", call. = FALSE)
    }
    .fold_na <- vapply(fold_cols,
                       function(fc) anyNA(train_features[[fc]]), logical(1))
    if (any(.fold_na)) {
      warning("run_supervised_pipeline(): fold column(s) with NA in training rows: ",
              paste(fold_cols[.fold_na], collapse = ", "),
              " -- those rows get no held-out OOF prediction.", call. = FALSE)
    }
    msg("STEP C1-C2 - DM -> OOF (diagnostic) via run_oof_diagnostics()")
    pipe1 <- time_step("C2 run_dm_oof_pipeline", {
      stopifnot(!is.null(train_with_folds_gpkg) && file.exists(train_with_folds_gpkg))
      stopifnot(!is.null(blocks_gpkg) && file.exists(blocks_gpkg))

      oof_res <- run_oof_diagnostics(
        train_features   = train_features,
        scoring_features = scoring_features,
        config           = config,
        fold_cols        = fold_cols,   # detected dynamically just above
        # params = NULL: build from the canonical source of truth
        # (.of_canonical_xgb_params): logloss-first eval_metric, eta 0.05,
        # depth 5, min_child_weight 5, subsample 0.8, colsample 0.75, gamma 0,
        # lambda 1, alpha 0, scale_pos_weight = neg/pos. FRENTE 1: SAME block
        # as the FINAL stage.
        params           = NULL,
        feature_whitelist_override = feature_whitelist_override,
        feature_weights            = feature_weights,
        # OPTIONAL shape/size block (OtsuFire 0.12.0): the SAME cfg-derived local
        # fed to FINAL below (single cfg field -> OOF/FINAL cannot diverge).
        include_shape_features     = include_shape_features,
        # 2026-06-05 (D1 expose): forward OOF training knobs.
        nrounds_max      = oof_nrounds_max,
        early_stop       = oof_early_stop,
        seed_base        = oof_seed_base,
        out_dir          = dirs$`05_OOF`,
        matrix_dir       = dirs$`04_MATRIX`,
        # PHASE 2 (artifact_hard): when promotion ran, the labelled-summary
        # geometry must come from the AUGMENTED train_features layer (which
        # carries the promoted rows), not the pre-promotion train_with_folds
        # GPKG. Off / no-promotion -> the canonical train_with_folds GPKG, so the
        # labeled_oof_summary join is byte-identical to today.
        labelled_gpkg    = if (identical(.gpkg_train_aug, gpkg_features)) {
                             train_with_folds_gpkg
                           } else {
                             .gpkg_train_aug
                           },
        labelled_layer   = if (identical(.gpkg_train_aug, gpkg_features)) {
                             "train_with_folds"
                           } else {
                             "train_features"
                           },
        overwrite        = overwrite,
        # PHASE 2 (artifact_hard): per-row source weighting + the artifact_hard
        # source tag (OFF -> NULL / character(0) -> byte-identical).
        source_weights       = .ah_source_weights,
        artifact_hard_source = .ah_source_tag,
        total_weight_ratio   = .ah_total_weight_ratio,
        # Precision 1 (2026-06-07): shims already resolved at the public
        # boundary; suppress a second deprecation warning here.
        .internal_resolved = TRUE,
        # The SAME cap ratios forwarded to FINAL (so OOF and FINAL see identical
        # caps). OOF always uses the capped negative-sampling policy. GATE 6.5:
        # contextual cap removed (random + otsu only).
        random_to_burned_ratio                 = random_to_burned_ratio,
        otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
        val_frac          = final_val_frac,
        impute_numeric    = final_impute_numeric,
        impute_factor_missing = final_impute_factor_missing
      )
      # Preserve the historical `pipe1` shape consumed below
      # (pipe1$files), and keep the inner wrapper object available.
      oof_res$pipe
    })

    print(pipe1$files)
    
    msg("STEP C3 - Train final model + export final map (modular)")
    stopifnot(inherits(config, "otsufire_supervised_burned_config"))
    pipe2 <- time_step("C3 train_final_model_and_export_final_map", {

      oof_agg_csv <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_oof_agg.csv"))
      oof_summary_gpkg <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_labeled_oof_summary.gpkg"))
      stopifnot(file.exists(oof_agg_csv))
      stopifnot(file.exists(oof_summary_gpkg))

      # ---- D) TRAIN final model (07_FINAL_MODEL_V2) -------------------------
      # Faithful split of the former fused wrapper: same engine args/values
      # (the 3 caps, whitelist/weights, the OOF `qa` aggregate, overwrite).
      tm <- train_final_burned_model(
        # PHASE 2 (artifact_hard): read the AUGMENTED train_features (with the
        # promoted rows) when promotion ran; otherwise the canonical features
        # GPKG -> byte-identical.
        train_features = .gpkg_train_aug,
        config         = config,
        oof_agg        = oof_agg_csv,
        random_to_burned_ratio                 = random_to_burned_ratio,
        otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
        # 0.5.0: same feature space + weights as the OOF stage (KB1/KB2).
        feature_whitelist_override = feature_whitelist_override,
        feature_weights            = feature_weights,
        # OPTIONAL shape/size block (OtsuFire 0.12.0): the SAME cfg-derived local
        # fed to OOF above (single cfg field -> OOF/FINAL cannot diverge).
        include_shape_features     = include_shape_features,
        # PHASE 2 (artifact_hard): per-row source weighting + the artifact_hard
        # source tag (OFF -> NULL / character(0) -> byte-identical).
        source_weights             = .ah_source_weights,
        artifact_hard_source       = .ah_source_tag,
        total_weight_ratio         = .ah_total_weight_ratio,
        # 2026-06-05 (D1 expose): forward FINAL training knobs.
        sampling_seed              = final_sampling_seed,
        seed                       = final_seed,
        val_frac                   = final_val_frac,
        group_col                  = final_group_col,
        nrounds_max                = final_nrounds_max,
        early_stopping_rounds      = final_early_stopping_rounds,
        impute_numeric             = final_impute_numeric,
        impute_factor_missing      = final_impute_factor_missing,
        # Gate 1E (2026-06-09): thread the CANONICAL OOF structural feature-schema
        # fingerprint into the FINAL refit so it ASSERTS its own structural
        # contract == the OOF contract BEFORE training/saving (ERROR on mismatch).
        canonical_oof_fingerprint  =
          if (!is.null(pipe1$schema_guard)) pipe1$schema_guard$canonical else NULL,
        labelled_layer = "train_features",
        out_dir        = dirs$`07_FINAL_MODEL_V2`,
        overwrite      = overwrite,
        verbose        = TRUE
      )

      # ---- E) SCORE deterministic universe + export final map --------------
      # Same engine args/values: the trained model/recipe handed off in
      # memory, the labelled-features + OOF-summary join inputs, the
      # CURRENTYEAR_* temporal thresholds, the 08_SCORED/09_FINAL_MAP dirs.
      sm <- score_supervised_burned_map(
        # PHASE 2 (artifact_hard): the scoring universe is the UNCHANGED
        # scoring_features layer (same row count whether or not promotion ran).
        # labelled_features is the AUGMENTED train_features (so the promoted
        # rows' held-out p_oof can join through to p_burned_eval). Off -> the
        # canonical features GPKG -> byte-identical.
        scoring_features = .gpkg_train_aug,
        model            = tm$model,
        recipe           = tm$recipe,
        config           = config,
        oof_summary      = oof_summary_gpkg,
        labelled_features = .gpkg_train_aug,
        export_burned_like = TRUE,
        preyear_overlap_threshold = CURRENTYEAR_PREYEAR_OVERLAP_THR,
        hotspot_density_threshold = CURRENTYEAR_HOTSPOT_DENSITY_THR,
        temporal_penalty_floor    = CURRENTYEAR_TEMPORAL_PENALTY_FLOOR,
        qa_labelled_layer       = "labeled_oof_summary",
        labelled_features_layer = "train_features",
        scoring_layer           = "scoring_features",
        out_map_dir   = dirs$`09_FINAL_MAP`,
        out_score_dir = dirs$`08_SCORED`,
        overwrite     = overwrite,
        verbose       = TRUE
      )

      list(train = tm, score = sm)
    })

    print(pipe2$train$direct$files)
    print(pipe2$score$scored$files)

    # Gate 1E (2026-06-09): runtime feature-schema PARITY MANIFEST + Phase B
    # abort. Assemble the OOF (per-fold + canonical) / FINAL / scoring structural
    # fingerprints, the n_base/n_indicators/n_total counts, the contract version
    # and the guard RESULT into the run summary the package writes. The run
    # ABORTS here if the three legs are not compatible. The OOF (per-fold
    # equality + canonical), FINAL (vs canonical)
    # and scoring (vs saved FINAL) assertions already fired upstream; this is the
    # consolidated record + the cross-leg Phase B gate.
    schema_parity_manifest <- time_step("C3 feature-schema parity manifest", {
      mani <- .of_build_schema_parity_manifest(
        oof_guard  = pipe1$schema_guard,
        final_fp   = pipe2$train$schema_fingerprint,
        scoring_fp = pipe2$score$scoring_schema_fingerprint,
        # OtsuFire always runs the leakage-free protocol (per-fold refit recipe
        # fingerprints), so the cross-leg schema-parity gate is always enforced.
        phase_b    = TRUE
      )
      manifest_path <- file.path(dirs$`07_FINAL_MODEL_V2`,
                                 paste0(prefix_oof, "_feature_schema_parity.txt"))
      if (isTRUE(overwrite) || !file.exists(manifest_path)) {
        writeLines(c(
          "[feature_schema_parity_guard]",
          paste0("contract_version: ", mani$contract_version),
          paste0("phase_b: ", mani$phase_b),
          paste0("guard_result: ", mani$guard_result),
          paste0("n_base: ", mani$n_base),
          paste0("n_indicators: ", mani$n_indicators),
          paste0("n_total: ", mani$n_total),
          paste0("oof_canonical_fingerprint: ", mani$oof_canonical_fingerprint),
          paste0("final_fingerprint: ", mani$final_fingerprint),
          paste0("scoring_fingerprint: ", mani$scoring_fingerprint),
          "[oof_per_fold_fingerprints]",
          if (length(mani$oof_per_fold_fingerprints)) {
            paste0(names(mani$oof_per_fold_fingerprints), ": ",
                   mani$oof_per_fold_fingerprints)
          } else "none"
        ), con = manifest_path)
      }
      mani$manifest_path <- manifest_path
      mani
    })
    msg("[Gate 1E] feature-schema parity guard: %s (contract %s) | OOF=%s FINAL=%s SCORING=%s",
        schema_parity_manifest$guard_result,
        schema_parity_manifest$contract_version,
        schema_parity_manifest$oof_canonical_fingerprint,
        schema_parity_manifest$final_fingerprint,
        schema_parity_manifest$scoring_fingerprint)

    # §N+27 (2026-06-05): STEP C4 (append_supervised_registry) removed —
    # abandoned research line. `pipe2$score$burned_like_scored` remains a
    # normal public output (written by C3 as `_burned_like_scored.gpkg`);
    # it is simply no longer routed into a registry.

    # ============================================================
    # PHASE 2 (artifact_hard): CONSOLIDATED supervised training-pool layer.
    # A single, one-row-per-example view of the EXACT training rows + weights the
    # engine received (used_in_training <=> FINAL capped frame; sample_weight
    # re-resolved with the engine helper), plus the deterministic origin and the
    # artifact_hard provenance. Written ALONGSIDE the existing separate
    # burned/unburned pool outputs (additive). Emitted in BOTH modes: with
    # artifact_hard OFF it describes the baseline pool (no artifact_hard rows,
    # weight 1). It is a NEW file and does not change the model or any existing
    # artifact, so OFF-baseline parity of the existing outputs is preserved.
    # ------------------------------------------------------------
    training_pool_layer_path <- time_step("C3b supervised_training_pool layer", {
      tryCatch(
        .of_write_supervised_training_pool(
          train_features   = train_features,
          scoring_features = scoring_features,
          training_ok      = pipe2$train$training_ok,
          oof_agg          = oof_agg_csv,
          promo            = promo,
          config           = config,
          target_year      = target_year,
          out_dir          = dirs$`07_FINAL_MODEL_V2`,
          overwrite        = overwrite,
          id_col           = "fire_uid"),
        error = function(e) {
          msg("WARN consolidated training-pool layer not written: %s",
              conditionMessage(e))
          NULL
        })
    })
    if (!is.null(training_pool_layer_path)) {
      msg("Consolidated training-pool layer: %s", training_pool_layer_path)
    }
  }

  save_timelog()
  msg("Timing CSV saved: %s", timing_csv)
  cat(sprintf("SUPERVISED MASTER pipeline finished at: %s | year=%s | scenario=%s\n",
              format(Sys.time()), target_year, scenario))
  
  invisible(list(
    year        = target_year,
    scenario    = scenario,
    result_dir  = result_dir,
    timing_csv  = timing_csv,
    unb_out_gpkg = unb_out_path_actual,
    # §N+26: always all_sources, so common_unb_dir is always the Otsu residual
    # negative root. (The deterministic_common_unb_dir is still created above
    # as a byte-identical filesystem side effect, but is no longer selected.)
    common_unb_dir = otsu_negative_root_dir,
    # Gate 1E (2026-06-09): the runtime feature-schema parity guard manifest
    # (OOF/FINAL/scoring fingerprints, counts, contract version, pass/abort).
    # NULL when modeling did not run this invocation.
    feature_schema_parity = schema_parity_manifest,
    # PHASE 2 (artifact_hard): path to the consolidated supervised training-pool
    # layer (NULL when modeling did not run / it could not be written).
    training_pool_layer =
      if (exists("training_pool_layer_path", inherits = FALSE))
        training_pool_layer_path else NULL
  ))
}
