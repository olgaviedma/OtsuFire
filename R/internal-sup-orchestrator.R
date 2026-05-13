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
# Negative-class generation can use either:
# - deterministic drops + random background  (deterministic_direct mode)
# - all_sources: deterministic drops + random background + Otsu current-year patches
# -------------------------------------------------------------------------

# ============================== 1) PACKAGES ==============================
# Block 7: top-level library() calls removed. All these packages are
# now declared in DESCRIPTION Imports and resolved via the package
# namespace at load time.
sf::sf_use_s2(FALSE)

# ============================== 2) GLOBAL CONFIG =========================
# Block 4B: project-root resolution is dispatcher-injected at runtime.
find_project_paths_file <- function(start = getwd()) NULL
project_paths <- list()
data_base <- NULL
script_base <- NULL
composite_base <- NULL
result_name <- NULL
supervised_functions_dir <- NULL

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
# ------------------------------ REGISTRY ---------------------------------
LEGACY_RBR_TRAINING_REGISTRY_PATH <- project_paths$rbr_training_registry_path
RBR_TRAINING_REGISTRY_PATH <- NULL
REGISTRY_PROB_THR <- 0.90
REGISTRY_TOP_FRACTION <- 0.01

# ------------------------------ CURRENT-YEAR TEMPORAL FILTER -------------
CURRENTYEAR_PREYEAR_OVERLAP_THR <- 0.70
CURRENTYEAR_HOTSPOT_DENSITY_THR <- 0.001
CURRENTYEAR_TEMPORAL_PENALTY_FLOOR <- 0.10

# ------------------------------ UNBURNED OPTIONS -------------------------
# Negative-pool strategy - two modes are supported:
#
#   "deterministic_direct"  - Mode A / canonical baseline.
#       Negatives come from the current year's deterministic output only:
#       drop polygons (hard negatives) plus random burnable-background cells.
#       No Otsu patches. Use for benchmarking, reproducibility tests, and
#       package-vs-no-package equivalence validation.
#       *** This standalone script uses Mode A as the canonical baseline. ***
#
#   "all_sources"  - Full operational mode (used in MASTER pipeline).
#       Combines all three current-year unburned sources into the training pool:
#         (1) Deterministic drop polygons  -> contextual_exclusion / spectral_hard_negative
#         (2) Random burnable-background cells -> random_background_sampled
#         (3) Otsu current-year unburned patches -> otsu_unburned_sampled
#       Each source routes to its own training_group in train_final_model_direct()
#       with an independent cap. No historical/cross-year content is used.
#       See HANDOFF sections 23-27 for taxonomy and composition rules.
#
# DEPRECATED: "otsu_unburned_generation" is redirected to "all_sources" with a
#   warning. That old mode incorrectly excluded deterministic drops and background
#   cells; it is no longer a valid standalone mode.
#
# Do NOT silently unify all modes under "deterministic_direct" without an
# explicit methodological decision.
# ---- Minimum burned-pool size (F7 guard: years with fewer keep patches abort cleanly)
MIN_BURNED_POOL_N <- 5L

UNB_STRATEGY <- "deterministic_direct"
UNB_VERBOSE  <- TRUE
UNB_LEGACY_CODE_DIR <- NULL  # injected by dispatcher

# deterministic_direct strategy
UNB_EXCL_BUFFER_M   <- 500
UNB_RANDOM_RBR_Q    <- 0.50
UNB_N_RANDOM_CELLS  <- 1500
UNB_RANDOM_SEED     <- 42
UNB_RANDOM_PATCH_SIZE_CELLS <- 3

# Otsu pipeline parameters (used by all_sources mode)
UNB_LEGACY_OTSU_MODE <- "burnable_only"
UNB_LEGACY_OTSU_THRESHOLD <- 0
UNB_LEGACY_REFERENCE_OTSU_THRESHOLD <- 100
UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE <- 0
UNB_LEGACY_MIN_PIXELS <- 8
UNB_LEGACY_BUFFERS_M <- 90
UNB_LEGACY_CORE_THR <- 0.60
UNB_LEGACY_ALPHA_BOOST <- 0.25
UNB_LEGACY_MIN_BASE_BOOST <- 0.35
UNB_LEGACY_DIST_POWER <- 1
UNB_LEGACY_KEEP_HI <- 0.45
UNB_LEGACY_DROP_LO <- 0.15
UNB_LEGACY_USE_DROP <- TRUE
# AS12 (0.3.0): defaults flipped to FALSE so review/keep rows do not
# enter the legacy pool. They were unconditionally excluded by
# `train_final_model_direct()` anyway, so the previous TRUE defaults
# wasted ~30% of UNB_LEGACY_SAMPLE_N on rows the trainer dropped.
UNB_LEGACY_USE_REVIEW <- FALSE
UNB_LEGACY_USE_KEEP <- FALSE
UNB_LEGACY_DROP_MAX_S_PATCH <- 0.15
UNB_LEGACY_REVIEW_MAX_S_PATCH <- 0.45
UNB_LEGACY_KEEP_MAX_S_PATCH <- 0.70
UNB_LEGACY_SAMPLE_N <- 2000
UNB_LEGACY_SAMPLE_PROPS <- c(drop = 0.70, review = 0.25, keep = 0.05)
UNB_LEGACY_EXCL_BUFFER_M <- 0
UNB_LEGACY_MIN_AREA_HA <- 0
UNB_LEGACY_RANDOM_SEED <- 42
UNB_LEGACY_REUSE_EXISTING <- TRUE
UNB_LEGACY_WRITE_OUTPUT <- TRUE

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

append_supervised_registry <- function(
  burned_like_sf,
  registry_path,
  target_year,
  scenario,
  prob_thr = 0.95,
  top_fraction = 0.01,
  registry_layer = "burned_high_conf_registry",
  verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  build_registry_support_paths <- function(path) {
    list(
      registry_path = path,
      registry_dir = dirname(path),
      ranked_csv = file.path(dirname(path), paste0(registry_layer, "_ranked.csv")),
      summary_txt = file.path(dirname(path), paste0(registry_layer, "_summary.txt")),
      candidate_layer = paste0(registry_layer, "_candidates"),
      legacy_registry_path = LEGACY_RBR_TRAINING_REGISTRY_PATH
    )
  }

  resolve_registry_spec <- function() {
    if (!is.null(registry_path) && nzchar(registry_path)) {
      return(build_registry_support_paths(registry_path))
    }

    get_fire_mapping_registry_paths(
      data_base = data_base,
      scenario = scenario,
      registry_basename = "RBR_TRAINING_REGISTRY.gpkg",
      registry_layer = registry_layer
    )
  }

  compute_registry_scores <- function(x) {
    n_rows <- if (is.null(x)) 0L else nrow(x)

    score_current <- if ("p_burned_current_year" %in% names(x)) {
      suppressWarnings(as.numeric(x[["p_burned_current_year"]]))
    } else {
      rep(NA_real_, n_rows)
    }

    score_fallback <- if ("p_burned" %in% names(x)) {
      suppressWarnings(as.numeric(x[["p_burned"]]))
    } else {
      rep(NA_real_, n_rows)
    }

    score_value <- ifelse(is.finite(score_current), score_current, score_fallback)
    score_used <- ifelse(
      is.finite(score_current),
      "p_burned_current_year",
      ifelse(is.finite(score_fallback), "p_burned", NA_character_)
    )

    list(
      score_value = score_value,
      score_used = score_used
    )
  }

  registry_spec <- resolve_registry_spec()
  registry_path <- registry_spec$registry_path

  if (is.null(registry_path) || !nzchar(registry_path)) {
    return(list(path = registry_path, n_added = 0L, n_total = NA_integer_))
  }

  top_fraction <- suppressWarnings(as.numeric(top_fraction)[1])
  if (!is.finite(top_fraction) || top_fraction <= 0 || top_fraction > 1) {
    stop("'top_fraction' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }

  bind_sf_registry <- function(x, y) {
    crs_x <- sf::st_crs(x)
    crs_y <- sf::st_crs(y)
    if (!is.na(crs_x) && !is.na(crs_y) && crs_x != crs_y) {
      y <- sf::st_transform(y, crs_x)
    } else if (is.na(crs_x) && !is.na(crs_y)) {
      sf::st_crs(x) <- crs_y
    } else if (!is.na(crs_x) && is.na(crs_y)) {
      sf::st_crs(y) <- crs_x
    }

    geom_x <- attr(x, "sf_column")
    geom_y <- attr(y, "sf_column")
    if (!identical(geom_x, geom_y)) {
      names(y)[names(y) == geom_y] <- geom_x
      attr(y, "sf_column") <- geom_x
    }

    cols_all <- union(names(x), names(y))
    for (nm in setdiff(cols_all, names(x))) x[[nm]] <- NA
    for (nm in setdiff(cols_all, names(y))) y[[nm]] <- NA

    is_numeric_like <- function(v) {
      is.numeric(v) || is.integer(v) || is.logical(v)
    }

    for (nm in setdiff(cols_all, geom_x)) {
      if (!nm %in% names(x) || !nm %in% names(y)) next

      vx <- x[[nm]]
      vy <- y[[nm]]

      if (inherits(vx, "factor")) vx <- as.character(vx)
      if (inherits(vy, "factor")) vy <- as.character(vy)

      if (is_numeric_like(vx) && is_numeric_like(vy)) {
        x[[nm]] <- as.numeric(vx)
        y[[nm]] <- as.numeric(vy)
      } else if (!identical(class(vx)[1], class(vy)[1])) {
        x[[nm]] <- as.character(vx)
        y[[nm]] <- as.character(vy)
      } else {
        x[[nm]] <- vx
        y[[nm]] <- vy
      }
    }

    x <- x[, cols_all, drop = FALSE]
    y <- y[, cols_all, drop = FALSE]
    dplyr::bind_rows(x, y)
  }

  maybe_seed_registry_from_legacy <- function() {
    if (file.exists(registry_path)) {
      return(0L)
    }

    legacy_path <- registry_spec$legacy_registry_path
    if (is.null(legacy_path) || !nzchar(legacy_path) || !file.exists(legacy_path)) {
      return(0L)
    }

    legacy_layers <- tryCatch(sf::st_layers(legacy_path), error = function(e) NULL)
    if (is.null(legacy_layers)) {
      return(0L)
    }

    legacy_source_layer <- if (registry_spec$candidate_layer %in% legacy_layers$name) {
      registry_spec$candidate_layer
    } else if (registry_layer %in% legacy_layers$name) {
      registry_layer
    } else {
      return(0L)
    }

    legacy_sf <- tryCatch(
      sf::read_sf(legacy_path, layer = legacy_source_layer, quiet = TRUE),
      error = function(e) NULL
    )

    if (!inherits(legacy_sf, "sf") || !nrow(legacy_sf)) {
      return(0L)
    }

    if (!"scenario" %in% names(legacy_sf)) {
      msg(
        "Legacy registry migration skipped for scenario '%s': missing 'scenario' column in %s",
        scenario,
        legacy_path
      )
      return(0L)
    }

    legacy_sf <- legacy_sf |>
      dplyr::filter(
        tolower(trimws(as.character(.data$scenario))) ==
          tolower(trimws(as.character(scenario)))
      )

    if (!nrow(legacy_sf)) {
      return(0L)
    }

    dir.create(registry_spec$registry_dir, recursive = TRUE, showWarnings = FALSE)
    sf::st_write(
      legacy_sf,
      registry_path,
      layer = registry_spec$candidate_layer,
      delete_layer = TRUE,
      quiet = TRUE
    )

    msg(
      "Seeded scenario registry from legacy shared registry: %s | migrated=%d",
      registry_path,
      nrow(legacy_sf)
    )

    nrow(legacy_sf)
  }

  migrated_n <- maybe_seed_registry_from_legacy()
  dir.create(registry_spec$registry_dir, recursive = TRUE, showWarnings = FALSE)

  keep <- NULL
  if (inherits(burned_like_sf, "sf") && nrow(burned_like_sf)) {
    keep_scores <- compute_registry_scores(burned_like_sf)
    score_col <- if ("p_burned_current_year" %in% names(burned_like_sf)) {
      "p_burned_current_year"
    } else if ("p_burned" %in% names(burned_like_sf)) {
      "p_burned"
    } else {
      NA_character_
    }

    keep <- burned_like_sf |>
      dplyr::mutate(registry_score_value = keep_scores$score_value) |>
      dplyr::filter(
        .data$source_set == "deterministic",
        as.character(.data$class_final) != "keep",
        is.finite(.data$registry_score_value),
        .data$registry_score_value >= prob_thr
      ) |>
      dplyr::mutate(
        year = target_year,
        scenario = scenario,
        source_type = "supervised_high_conf_burned",
        registry_score_col = if (is.na(score_col)) "p_burned" else score_col,
        registry_prob_threshold = prob_thr,
        registry_added_at = as.character(Sys.time())
      )
  }

  existing <- NULL
  lyr_info <- tryCatch(sf::st_layers(registry_path), error = function(e) NULL)
  if (!is.null(lyr_info) && registry_spec$candidate_layer %in% lyr_info$name) {
    existing <- sf::read_sf(
      registry_path,
      layer = registry_spec$candidate_layer,
      quiet = TRUE
    )
  } else if (!is.null(lyr_info) && registry_layer %in% lyr_info$name) {
    existing <- sf::read_sf(
      registry_path,
      layer = registry_layer,
      quiet = TRUE
    )
  }
  # Guard: keep only rows that belong to this scenario.
  # Pre-existing files may contain entries from other scenarios if they were
  # created or seeded before the registry-separation fix (2026-04-20, HANDOFF 43).
  if (inherits(existing, "sf") && nrow(existing) && "scenario" %in% names(existing)) {
    n_before <- nrow(existing)
    existing <- existing |>
      dplyr::filter(
        tolower(trimws(as.character(.data$scenario))) ==
          tolower(trimws(as.character(scenario)))
      )
    n_dropped <- n_before - nrow(existing)
    if (n_dropped > 0L) {
      msg(
        "Registry scenario guard: dropped %d cross-scenario candidates from '%s' (scenario filter: '%s')",
        n_dropped,
        registry_path,
        scenario
      )
    }
  }

  combined <- if (inherits(existing, "sf") && nrow(existing) &&
                  inherits(keep, "sf") && nrow(keep)) {
    bind_sf_registry(existing, keep)
  } else if (inherits(existing, "sf") && nrow(existing)) {
    existing
  } else if (inherits(keep, "sf") && nrow(keep)) {
    keep
  } else {
    NULL
  }

  if (!inherits(combined, "sf") || !nrow(combined)) {
    msg(
      "Registry append skipped: no scenario candidates available for '%s' after migration and current-year filtering.",
      scenario
    )
    return(list(
      path = registry_path,
      n_added = 0L,
      n_total = 0L,
      n_candidates_total = 0L,
      ranked_csv = registry_spec$ranked_csv,
      summary_txt = registry_spec$summary_txt
    ))
  }

  dedupe_keys <- intersect(c("fire_uid", "year", "scenario", "source_type"), names(combined))
  if (length(dedupe_keys)) {
    combined <- combined |>
      dplyr::distinct(dplyr::across(dplyr::all_of(dedupe_keys)), .keep_all = TRUE)
  }

  combined_scores <- compute_registry_scores(combined)
  combined <- combined |>
    dplyr::mutate(
      registry_score_value = combined_scores$score_value,
      registry_score_used = combined_scores$score_used
    ) |>
    dplyr::filter(
      is.finite(.data$registry_score_value),
      .data$registry_score_value >= prob_thr
    )

  n_candidates_total <- nrow(combined)
  n_selected <- if (n_candidates_total > 0L) {
    max(1L, as.integer(ceiling(n_candidates_total * top_fraction)))
  } else {
    0L
  }

  if (!n_candidates_total) {
    msg(
      "Registry append skipped: no scenario candidates remained above %.2f for '%s'.",
      prob_thr,
      scenario
    )
    return(list(
      path = registry_path,
      n_added = if (inherits(keep, "sf")) nrow(keep) else 0L,
      n_total = 0L,
      n_candidates_total = 0L,
      ranked_csv = registry_spec$ranked_csv,
      summary_txt = registry_spec$summary_txt
    ))
  }

  if ("year" %in% names(combined) && "fire_uid" %in% names(combined)) {
    combined <- combined |>
      dplyr::arrange(
        dplyr::desc(.data$registry_score_value),
        dplyr::desc(.data$year),
        .data$fire_uid
      )
  } else if ("year" %in% names(combined)) {
    combined <- combined |>
      dplyr::arrange(
        dplyr::desc(.data$registry_score_value),
        dplyr::desc(.data$year)
      )
  } else {
    combined <- combined |>
      dplyr::arrange(dplyr::desc(.data$registry_score_value))
  }

  combined <- combined |>
    dplyr::mutate(
      registry_rank = dplyr::row_number(),
      registry_top_fraction_used = top_fraction
    )

  selected <- combined |>
    dplyr::slice_head(n = n_selected) |>
    dplyr::mutate(
      selected_registry_high_conf = TRUE,
      registry_threshold_used = prob_thr
    )

  sf::st_write(
    combined,
    registry_path,
    layer = registry_spec$candidate_layer,
    delete_layer = TRUE,
    quiet = TRUE
  )

  sf::st_write(
    selected,
    registry_path,
    layer = registry_layer,
    delete_layer = TRUE,
    quiet = TRUE
  )

  ranked_df <- selected |>
    sf::st_drop_geometry() |>
    dplyr::arrange(dplyr::desc(.data$registry_score_value))

  utils::write.csv(ranked_df, registry_spec$ranked_csv, row.names = FALSE)
  writeLines(c(
    paste0("registry_path: ", registry_path),
    paste0("registry_layer: ", registry_layer),
    paste0("registry_candidate_layer: ", registry_spec$candidate_layer),
    paste0("registry_scenario_dir: ", registry_spec$registry_dir),
    paste0("updated_at: ", as.character(Sys.time())),
    paste0("target_year_appended: ", target_year),
    paste0("scenario_appended: ", scenario),
    paste0("registry_score_rule: descending p_burned_current_year then p_burned"),
    paste0("registry_threshold_used: ", prob_thr),
    paste0("registry_top_fraction_used: ", top_fraction),
    paste0("n_migrated_from_legacy: ", migrated_n),
    paste0("n_added_current_run: ", if (inherits(keep, "sf")) nrow(keep) else 0L),
    paste0("n_candidates_total: ", n_candidates_total),
    paste0("n_selected_registry: ", nrow(selected))
  ), registry_spec$summary_txt)

  msg(
    "Registry updated: %s | added=%d | candidates=%d | selected_top=%.2f%% -> %d",
    registry_path,
    if (inherits(keep, "sf")) nrow(keep) else 0L,
    n_candidates_total,
    top_fraction * 100,
    nrow(selected)
  )
  list(
    path = registry_path,
    n_added = if (inherits(keep, "sf")) nrow(keep) else 0L,
    n_total = nrow(selected),
    n_candidates_total = n_candidates_total,
    ranked_csv = registry_spec$ranked_csv,
    summary_txt = registry_spec$summary_txt,
    candidate_layer = registry_spec$candidate_layer
  )
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
                                    contextual_exclusion_to_burned_ratio   = 0.25,
                                    spectral_hard_negative_to_burned_ratio = 1.0,
                                    random_to_burned_ratio                 = 1.0,
                                    otsu_unburned_to_burned_ratio          = 1.0,
                                    feature_whitelist_override             = NULL,
                                    feature_weights                        = NULL,
                                    reuse_upstream                         = FALSE) {
  
  cat("\n============================================================\n")
  cat(sprintf("RUNNING SUPERVISED PIPELINE | year=%s | scenario=%s\n", target_year, scenario))
  cat("============================================================\n")

  resolved_registry_path <- if (!is.null(RBR_TRAINING_REGISTRY_PATH) &&
                                nzchar(RBR_TRAINING_REGISTRY_PATH)) {
    RBR_TRAINING_REGISTRY_PATH
  } else {
    get_fire_mapping_registry_paths(
      data_base = data_base,
      scenario = scenario,
      registry_basename = "RBR_TRAINING_REGISTRY.gpkg",
      registry_layer = "burned_high_conf_registry"
    )$registry_path
  }
  
  # ------------------------------ ROOT PATHS -----------------------------
  result_dir <- file.path(
    data_base, "Results", target_year, result_name, "SUPERVISED", scenario
  )
  
  deterministic_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", scenario
  )
  
  deterministic_decisions_dir <- file.path(deterministic_dir, "05_DECISIONS")
  internal_decisions_gpkg <- file.path(deterministic_decisions_dir, "internal_decisions.gpkg")
  
  # Carpeta comon por a+/-o para reutilizar OTSU + polygonize
  deterministic_common_unb_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", "_COMMON_UNBURNED"
  )
  legacy_unb_root_dir <- file.path(result_dir, "_LEGACY_UNBURNED")
  
  # Salida final unburned por escenario
  unb_out_gpkg <- file.path(
    deterministic_dir, "UNBURNED",
    sprintf("%d_%s_unburned.gpkg", target_year, scenario)
  )
  unb_out_path_actual <- unb_out_gpkg
  
  dir.create(deterministic_common_unb_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(unb_out_gpkg), recursive = TRUE, showWarnings = FALSE)
  dir.create(legacy_unb_root_dir, recursive = TRUE, showWarnings = FALSE)
  
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

  # Normalize deprecated strategy names to their canonical replacements.
  # "otsu_unburned_generation" was incorrect: it excluded deterministic drops
  # and random background cells, leaving only Otsu patches in the pool.
  # "all_sources" is the correct operational mode that includes all three sources.
  if (identical(UNB_STRATEGY, "otsu_unburned_generation")) {
    warning(
      paste0(
        "'otsu_unburned_generation' is deprecated as a standalone unburned strategy. ",
        "That mode excluded deterministic drops and random background cells. ",
        "Redirecting to 'all_sources', which combines all three current-year sources."
      ),
      call. = FALSE
    )
    UNB_STRATEGY <- "all_sources"
  }

  # ============================== 5b) COMMON INPUTS ======================
  peninsula_shapefile <- file.path(data_base, "Borders", "Iberian_peninsula.shp")
  topo_path           <- file.path(data_base, "Topography", "elevation_slope.tif")
  hotspots_path       <- file.path(data_base, "Hotspots", paste0("hotspots_iberia_", target_year, ".geojson"))
  
  corine_year         <- get_corine_year(target_year)
  corine_raster_path  <- file.path(data_base, "Corine_Masks", paste0("CLC_", corine_year, "_peninsula.tif"))
  strata_tif_path     <- file.path(data_base, "Corine_Masks", "STRATA", paste0("strata_CLC_", corine_year, "_res30.tif"))
  strata_lut_csv      <- file.path(data_base, "Corine_Masks", "LUT", "lut_full_strata8_v1.csv")
  mask_tif_path       <- file.path(data_base, "Mask_StudyArea", "mask_Peninsula_3035.tif")
  mask_shp_path       <- file.path(data_base, "Mask_StudyArea", "mask_Peninsula_3035.shp")
  ref_tif_path        <- file.path(data_base, "Fires", "Validation_fires_burneable_verano", paste0("Effis_CA_", target_year, "_maskKeep_summer.tif"))
  ref_shp_path        <- file.path(data_base, "Fires", "Validation_fires_burneable_verano", paste0("Effis_CA_", target_year, "_maskKeep_summer.shp"))
  
  one_year_tif <- file.path(
    composite_base, result_name,
    paste0("MinMin_", target_year, "_mosaic_res90m.tif")
  )
  if (!file.exists(one_year_tif)) {
    one_year_tif <- file.path(
      composite_base, "Min_Min",
      paste0("MinMin_", target_year, "_mosaic_res90m.tif")
    )
  }
  
  rbr_aw_tif <- file.path(
    composite_base, "Autumn",
    paste0("mean_mean_", target_year, "_mosaic.tif")
  )
  
  # ------------------------------ CHECKS ---------------------------------
  stopifnot(dir.exists(deterministic_dir))
  stopifnot(dir.exists(deterministic_decisions_dir))
  stopifnot(file.exists(internal_decisions_gpkg))
  stopifnot(file.exists(one_year_tif))
  stopifnot(file.exists(peninsula_shapefile))
  stopifnot(file.exists(corine_raster_path))
  stopifnot(file.exists(strata_tif_path))
  stopifnot(file.exists(strata_lut_csv))
  stopifnot(file.exists(mask_tif_path))
  stopifnot(file.exists(mask_shp_path))
  stopifnot(file.exists(ref_tif_path))
  stopifnot(file.exists(ref_shp_path))
  stopifnot(file.exists(topo_path))
  stopifnot(file.exists(rbr_aw_tif))
  
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
  # Phase B (0.2.2): skip ecoregion I/O here because extract_features()
  # no longer computes eco_major unless use_ecoregions = TRUE.
  
  hotspots_sf_base <- if (file.exists(hotspots_path)) {
    sf::read_sf(hotspots_path)
  } else {
    msg("WARNING: no existe hotspots: %s (se desactiva use_hotspots)", hotspots_path)
    sf::st_sf(
      frp = numeric(),
      confidence = numeric(),
      year = integer(),
      month = integer(),
      geometry = sf::st_sfc(crs = sf::NA_crs_)
    )[0, ]
  }
  
  # ============================== 6) POOLS ===============================
  if (isTRUE(DO_POOLS)) {
    
    # ---------------------------------------------------------------------
    # A1. Read internal_decisions from deterministic pipeline (scenario)
    # ---------------------------------------------------------------------
    msg("STEP A1 - Read internal_decisions")
    internal_sf <- time_step("A1 Read internal_decisions", {
      msg("Using internal_decisions_gpkg: %s", internal_decisions_gpkg)
      
      x <- sf::read_sf(internal_decisions_gpkg, layer = "internal_decisions") |>
        mutate(
          class = as.character(class_final)
        ) |>
        sanitize_polygons() |>
        ensure_area_ha()
      
      x <- drop_empty_sf(x, tag = "internal_sf", dump_dir = dirs$`99_LOGS_EMPTY`)
      check_sf(x, "internal_sf")
      x
    })
    
    crs_master <- sf::st_crs(internal_sf)
    stopifnot(!is.na(crs_master))
    
    # ---------------------------------------------------------------------
    # A2. Build clean internal set
    # ---------------------------------------------------------------------
    msg("STEP A2 - Build internal_clean")
    internal_clean <- time_step("A2 Build internal_clean", {
      x <- internal_sf |>
        mutate(source = "internal") |>
        sanitize_polygons() |>
        ensure_area_ha()
      
      x <- drop_empty_sf(x, tag = "internal_clean", dump_dir = dirs$`99_LOGS_EMPTY`)
      check_sf(x, "internal_clean")
      x
    })
    
    # ---------------------------------------------------------------------
    # A3. QA burned/review pools from deterministic decisions
    # keep   -> candidate burned positives
    # review -> candidate unlabeled polygons
    # ---------------------------------------------------------------------
    msg("STEP A3 - QA deterministic pools and build burned_pool + review_pool")
    poolsA <- time_step("A3 QA deterministic pools + build supervised pools", {
      qa_det <- audit_deterministic_pools(
        internal_sf = internal_clean,
        out_dir = dirs$`01_POOLS`,
        prefix = sprintf("%d_%s_deterministic_pool_qa", target_year, scenario),
        save_outputs = TRUE,
        verbose = TRUE,
        raw_class_col = "class",
      )

      internal_qc <- qa_det$audited_sf |>
        mutate(
          raw_class = as.character(raw_class),
          class = as.character(class_audited)
        ) |>
        sanitize_polygons() |>
        ensure_area_ha()

      burned_pool <- internal_qc |>
        filter(class == "keep") |>
        mutate(source = "internal_keep_qc") |>
        sanitize_polygons() |>
        ensure_area_ha() |>
        mutate(
          class = "burned",
          fire_uid = paste0(target_year, "_", scenario, "_B_", row_number())
        )

      burned_pool <- drop_empty_sf(
        burned_pool,
        tag = "burned_pool",
        dump_dir = dirs$`99_LOGS_EMPTY`
      )
      check_sf_if_nonempty(burned_pool, "burned_pool")

      review_pool <- internal_qc |>
        filter(class == "review") |>
        mutate(source = "internal_review_qc") |>
        sanitize_polygons() |>
        ensure_area_ha() |>
        mutate(
          class = "review",
          fire_uid = paste0(target_year, "_", scenario, "_R_", row_number())
        )

      review_pool <- drop_empty_sf(
        review_pool,
        tag = "review_pool",
        dump_dir = dirs$`99_LOGS_EMPTY`
      )
      check_sf_if_nonempty(review_pool, "review_pool")

      scoring_pool <- internal_qc |>
        mutate(
          class = as.character(raw_class),
          source = "deterministic_all_qc",
          fire_uid = paste0(target_year, "_", scenario, "_D_", row_number())
        ) |>
        sanitize_polygons() |>
        ensure_area_ha()

      scoring_pool <- drop_empty_sf(
        scoring_pool,
        tag = "scoring_pool",
        dump_dir = dirs$`99_LOGS_EMPTY`
      )
      check_sf_if_nonempty(scoring_pool, "scoring_pool")

      list(
        internal_qc = internal_qc,
        burned_pool = burned_pool,
        review_pool = review_pool,
        scoring_pool = scoring_pool,
        qa_det = qa_det
      )
    })
    
    internal_qc  <- poolsA$internal_qc
    burned_pool  <- poolsA$burned_pool
    review_pool  <- poolsA$review_pool
    scoring_pool <- poolsA$scoring_pool
    
    # ---------------------------------------------------------------------
    # A3b. Pool sanity checks
    # ---------------------------------------------------------------------
    msg("STEP A3b - Sanity checks on pools")
    time_step("A3b Pool sanity checks", {
      
      n_burned <- nrow(burned_pool)
      n_rev    <- nrow(review_pool)
      
      msg("Pools: burned=%d | review=%d", n_burned, n_rev)
      
      if (n_burned < min_burned_pool_n) {
        stop(sprintf(
          "Supervised run aborted: insufficient burned pool (n=%d, required>=%d). Year likely has too few detected fires for supervised modeling.",
          n_burned, min_burned_pool_n
        ))
      }
      
      if (n_rev == 0) {
        stop("Supervised run aborted: review pool is empty.")
      }
      
      NULL
    })
    
    # ---------------------------------------------------------------------
    # A4. Generate unburned datasets
    # Strategy can be either:
    # - deterministic_direct  (Mode A: det drops + random background; no Otsu)
    # - all_sources           (operational: det drops + random background + Otsu)
    # Final negatives are then merged with burned_pool to create train_labeled.
    # ---------------------------------------------------------------------
    msg("STEP A4 - Generate unburned datasets [%s]", UNB_STRATEGY)
    unb <- time_step("A4 Generate unburned datasets", {
      if (identical(UNB_STRATEGY, "deterministic_direct")) {
        # Mode A / validation baseline: deterministic drops + random background only.
        # No Otsu patches. Use for benchmarking and equivalence validation.
        res_unb <- build_unburned_from_deterministic_decisions(
          target_year = target_year,
          scenario_name = scenario,
          data_base = data_base,
          result_name = result_name,
          composite_base = composite_base,
          exclude_buffer_m = UNB_EXCL_BUFFER_M,
          n_random_cells = UNB_N_RANDOM_CELLS,
          random_rbr_q = UNB_RANDOM_RBR_Q,
          random_seed = UNB_RANDOM_SEED,
          random_patch_size_cells = UNB_RANDOM_PATCH_SIZE_CELLS,
          overwrite_output = TRUE,
          out_gpkg = unb_out_gpkg,
          verbose = UNB_VERBOSE
        )

        out_path <- res_unb$out_gpkg
        stopifnot(file.exists(out_path))

        unburned_hard <- to_crs_safe(res_unb$unburned_hard, crs_master) |>
          mutate(source = "deterministic_drop_hard")

        unburned_random <- to_crs_safe(res_unb$unburned_random, crs_master) |>
          mutate(source = "random_burnable_background")

        unburned_final_raw <- to_crs_safe(res_unb$unburned_final, crs_master) |>
          mutate(source = as.character(source))

        exclusion_buffer <- to_crs_safe(res_unb$exclusion_buffer, crs_master)

      } else if (identical(UNB_STRATEGY, "all_sources")) {
        # Full operational mode: all three current-year unburned sources combined.
        # Each source gets its own training_group in train_final_model_direct():
        #   deterministic_drop_hard  -> contextual_exclusion / spectral_hard_negative
        #   random_burnable_background -> random_background_sampled
        #   otsu_patch_residual       -> otsu_unburned_sampled
        # unburned_hard / unburned_random below reflect the deterministic sub-pools
        # only (written as auxiliary GPKG layers in step A6; not used in modeling).

        # ---- Part 1: deterministic drops + random burnable background ----
        res_det <- build_unburned_from_deterministic_decisions(
          target_year             = target_year,
          scenario_name           = scenario,
          data_base               = data_base,
          result_name             = result_name,
          composite_base          = composite_base,
          exclude_buffer_m        = UNB_EXCL_BUFFER_M,
          n_random_cells          = UNB_N_RANDOM_CELLS,
          random_rbr_q            = UNB_RANDOM_RBR_Q,
          random_seed             = UNB_RANDOM_SEED,
          random_patch_size_cells = UNB_RANDOM_PATCH_SIZE_CELLS,
          overwrite_output        = TRUE,
          out_gpkg                = unb_out_gpkg,
          verbose                 = UNB_VERBOSE
        )
        out_path <- res_det$out_gpkg
        stopifnot(file.exists(out_path))

        unburned_hard   <- to_crs_safe(res_det$unburned_hard, crs_master) |>
          mutate(source = "deterministic_drop_hard")
        unburned_random <- to_crs_safe(res_det$unburned_random, crs_master) |>
          mutate(source = "random_burnable_background")
        det_final_raw   <- to_crs_safe(res_det$unburned_final, crs_master) |>
          mutate(source = as.character(source))
        exclusion_buffer <- to_crs_safe(res_det$exclusion_buffer, crs_master)

        # ---- Part 2: Otsu current-year unburned patches ----
        res_otsu <- build_unburned_from_legacy_pipeline(
          target_year              = target_year,
          scenario_name            = scenario,
          data_base                = data_base,
          result_name              = result_name,
          composite_base           = composite_base,
          severity_raster_path     = one_year_tif,
          legacy_code_dir          = get0("UNB_LEGACY_CODE_DIR", ifnotfound = NULL),
          script_base              = script_base,  # F2 fix: explicit path avoids getwd() dependency
          # AS01 (0.3.0): tool paths injected by dispatcher from config$tool_paths.
          python_exe               = python_exe,
          gdal_polygonize_script   = gdal_polygonize_script,
          gdalwarp_path            = gdalwarp_path,
          ogr2ogr_exe              = ogr2ogr_exe,
          otsu_mode                = UNB_LEGACY_OTSU_MODE,
          otsu_threshold           = UNB_LEGACY_OTSU_THRESHOLD,
          reference_otsu_threshold = UNB_LEGACY_REFERENCE_OTSU_THRESHOLD,
          min_otsu_threshold_value = UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE,
          min_pixels               = UNB_LEGACY_MIN_PIXELS,
          buffers_m                = UNB_LEGACY_BUFFERS_M,
          core_thr                 = UNB_LEGACY_CORE_THR,
          alpha_boost              = UNB_LEGACY_ALPHA_BOOST,
          min_base_boost           = UNB_LEGACY_MIN_BASE_BOOST,
          dist_power               = UNB_LEGACY_DIST_POWER,
          keep_hi                  = UNB_LEGACY_KEEP_HI,
          drop_lo                  = UNB_LEGACY_DROP_LO,
          use_drop                 = UNB_LEGACY_USE_DROP,
          use_review               = UNB_LEGACY_USE_REVIEW,
          use_keep                 = UNB_LEGACY_USE_KEEP,
          drop_max_s_patch         = UNB_LEGACY_DROP_MAX_S_PATCH,
          review_max_s_patch       = UNB_LEGACY_REVIEW_MAX_S_PATCH,
          keep_max_s_patch         = UNB_LEGACY_KEEP_MAX_S_PATCH,
          sample_n                 = UNB_LEGACY_SAMPLE_N,
          sample_props             = UNB_LEGACY_SAMPLE_PROPS,
          exclude_buffer_m         = UNB_LEGACY_EXCL_BUFFER_M,
          min_area_ha              = UNB_LEGACY_MIN_AREA_HA,
          random_seed              = UNB_LEGACY_RANDOM_SEED,
          reuse_existing           = UNB_LEGACY_REUSE_EXISTING,
          write_unburned           = UNB_LEGACY_WRITE_OUTPUT,
          out_root_dir             = legacy_unb_root_dir,
          verbose                  = UNB_VERBOSE
        )
        otsu_pool    <- to_crs_safe(res_otsu$unburned$legacy_unburned_pool, crs_master)
        otsu_sampled <- to_crs_safe(res_otsu$unburned$legacy_unburned_sampled, crs_master)
        if (nrow(otsu_sampled) == 0L && nrow(otsu_pool) > 0L) otsu_sampled <- otsu_pool

        otsu_raw <- otsu_sampled |>
          mutate(
            source = "otsu_patch_residual",
            neg_type = dplyr::case_when(
              as.character(.data$legacy_decision) == "drop"   ~ "otsu_patch_drop",
              as.character(.data$legacy_decision) == "review" ~ "otsu_patch_review",
              as.character(.data$legacy_decision) == "keep"   ~ "otsu_patch_keep",
              TRUE ~ as.character(.data$neg_type)
            )
          )

        # ---- Part 3: combine all sources ----
        # All three sources are present; train_final_model_direct() routes them
        # to independent training_groups via source/neg_type fields.
        #
        # Use geometry-safe binding: drop geometry, bind data frames, re-attach.
        # dplyr::bind_rows on sf objects with mismatched column counts can
        # produce duplicate "geometry" column names (dplyr treats the sfc column
        # as a regular column when column sets differ).
        .bind_sf_safe <- function(a, b) {
          df_a <- sf::st_drop_geometry(a)
          df_b <- sf::st_drop_geometry(b)
          geom <- c(sf::st_geometry(a), sf::st_geometry(b))
          combined <- dplyr::bind_rows(df_a, df_b)
          combined[["geometry"]] <- geom
          sf::st_as_sf(combined, sf_column_name = "geometry",
                       crs = sf::st_crs(a))
        }
        unburned_final_raw <- .bind_sf_safe(det_final_raw, otsu_raw)
        # res_unb is only assigned in the deterministic_direct branch;
        # set to NULL here so the shared return list below resolves correctly.
        res_unb <- NULL

      } else {
        stop(
          sprintf(
            "Unknown UNB_STRATEGY='%s'. Expected 'deterministic_direct' or 'all_sources'.",
            UNB_STRATEGY
          )
        )
      }

      check_sf(unburned_final_raw, "unburned_final_raw")
      check_sf_if_nonempty(unburned_hard, "unburned_hard")
      check_sf_if_nonempty(unburned_random, "unburned_random")
      check_sf_if_nonempty(exclusion_buffer, "exclusion_buffer")

      unburned_pool <- unburned_final_raw |>
        sanitize_polygons() |>
        ensure_area_ha() |>
        mutate(
          class = "unburned",
          fire_uid = paste0(target_year, "_", scenario, "_U_", row_number())
        )

      unburned_pool <- drop_empty_sf(
        unburned_pool,
        tag = "unburned_pool",
        dump_dir = dirs$`99_LOGS_EMPTY`
      )
      check_sf(unburned_pool, "unburned_pool")

      msg("UNB out path: %s", out_path %||% NA_character_)
      msg("Unburned final count: %d", nrow(unburned_pool))

      list(
        wrapper_res        = res_unb,
        wrapper_out_gpkg   = out_path,
        unburned_hard      = unburned_hard,
        unburned_random    = unburned_random,
        unburned_final_raw = unburned_final_raw,
        exclusion_buffer   = exclusion_buffer,
        unburned_pool      = unburned_pool
      )
    })
    
    unburned_hard      <- unb$unburned_hard
    unburned_random    <- unb$unburned_random
    unburned_final_raw <- unb$unburned_final_raw
    exclusion_buffer   <- unb$exclusion_buffer
    unburned_pool      <- unb$unburned_pool
    unb_out_path_actual <- unb$wrapper_out_gpkg %||% unb_out_gpkg
    
    # ---------------------------------------------------------------------
    # A5. Final labelled training set = burned + unburned
    # ---------------------------------------------------------------------
    msg("STEP A5 - Build train_labeled = burned + unburned")
    train_labeled <- time_step("A5 Build train_labeled", {
      x <- bind_rows(
        burned_pool,
        unburned_pool
      ) |>
        sanitize_polygons() |>
        ensure_area_ha()
      
      x <- drop_empty_sf(x, tag = "train_labeled", dump_dir = dirs$`99_LOGS_EMPTY`)
      check_sf(x, "train_labeled")
      stopifnot(!anyDuplicated(x$fire_uid))
      
      msg("train_labeled counts: burned=%d | unburned=%d",
          sum(as.character(x$class) == "burned"),
          sum(as.character(x$class) == "unburned"))
      x
    })
    
    # ---------------------------------------------------------------------
    # A6. Write pools + auxiliary unburned datasets
    # ---------------------------------------------------------------------
    msg("STEP A6 - Write pools to 01_POOLS")
    time_step("A6 Write pools GPKG", {
      safe_write_gpkg(burned_pool,       gpkg_pools_out, layer = "burned_pool",       tag = "burned_pool_write")
      safe_write_gpkg(unburned_pool,     gpkg_pools_out, layer = "unburned_pool",     tag = "unburned_pool_write")
      safe_write_gpkg(review_pool,       gpkg_pools_out, layer = "review_pool",       tag = "review_pool_write")
      safe_write_gpkg(scoring_pool,      gpkg_pools_out, layer = "scoring_pool",      tag = "scoring_pool_write")
      safe_write_gpkg(train_labeled,     gpkg_pools_out, layer = "train_labeled",     tag = "train_labeled_write")
      
      if (nrow(unburned_hard) > 0) {
        safe_write_gpkg(unburned_hard, gpkg_pools_out, layer = "unburned_hard", tag = "unburned_hard_write")
      }
      if (nrow(unburned_random) > 0) {
        safe_write_gpkg(unburned_random, gpkg_pools_out, layer = "unburned_random", tag = "unburned_random_write")
      }
      if (nrow(unburned_final_raw) > 0) {
        safe_write_gpkg(unburned_final_raw, gpkg_pools_out, layer = "unburned_final_raw", tag = "unburned_final_raw_write")
      }
      if (nrow(exclusion_buffer) > 0) {
        safe_write_gpkg(exclusion_buffer, gpkg_pools_out, layer = "exclusion_buffer", tag = "exclusion_buffer_write")
      }
      
      NULL
    })
    
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
    
    msg("STEP B2 - Make spatial folds")
    res_folds <- time_step("B1 make_block_folds", {
      make_block_folds(
        train_labelled_sf = train_labeled_sf,
        fire_id_col       = "fire_uid",
        split_unit        = "fire",
        block_sizes_m = c(5000, 3000, 2000),
        k_candidates  = c(5, 4, 3),
        min_burned_units_per_fold = 3,
        min_pos_blocks_per_fold   = 10,
        n_repeats   = 2,
        seed_base   = 42,
        lonlat_action = "transform",
        out_dir     = dirs$`02_FOLDS`,
        target_year = target_year,
        verbose     = TRUE,
        write_train_with_folds_gpkg = TRUE,
        write_blocks_gpkg = TRUE,
        write_folds_csv   = TRUE
      )
    })
    
    if (!is.null(res_folds$saved)) {
      if (!is.null(res_folds$saved$train_gpkg))  train_with_folds_gpkg <- res_folds$saved$train_gpkg
      if (!is.null(res_folds$saved$blocks_gpkg)) blocks_gpkg          <- res_folds$saved$blocks_gpkg
    }
    
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
    feature_layers <- function(dsn) {
      if (!file.exists(dsn)) return(character(0))
      tryCatch(as.character(sf::st_layers(dsn)$name), error = function(e) character(0))
    }

    feature_layer_n <- function(dsn, layer) {
      if (!file.exists(dsn)) return(NA_integer_)
      if (!layer %in% feature_layers(dsn)) return(NA_integer_)
      out <- tryCatch(
        nrow(sf::read_sf(dsn, layer = layer, quiet = TRUE)),
        error = function(e) NA_integer_
      )
      as.integer(out)
    }

    safe_remove_features_gpkg <- function(path) {
      if (!file.exists(path)) return(invisible(TRUE))
      unlink(path, force = TRUE)
      if (file.exists(path)) file.remove(path)
      if (file.exists(path)) {
        stop("No se pudo borrar features_geometry.gpkg antes de reconstruir: ", path, call. = FALSE)
      }
      invisible(TRUE)
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
    
    should_rebuild_features <- FALSE
    rebuild_reasons <- character(0)
    
    if (file.exists(features_gpkg)) {
      feat_time  <- file.info(features_gpkg)$mtime
      pools_time <- file.info(gpkg_pools_out)$mtime
      folds_time <- file.info(train_with_folds_gpkg)$mtime
      feat_layers <- feature_layers(features_gpkg)
      
      if (isTRUE(feat_time < pools_time)) {
        should_rebuild_features <- TRUE
        rebuild_reasons <- c(rebuild_reasons, "features_geometry.gpkg es mas antiguo que 01_POOLS")
      }
      if (isTRUE(feat_time < folds_time)) {
        should_rebuild_features <- TRUE
        rebuild_reasons <- c(rebuild_reasons, "features_geometry.gpkg es mas antiguo que 02_FOLDS")
      }
      
      n_train_expected <- tryCatch(
        nrow(safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "train_with_folds_count_check")),
        error = function(e) NA_integer_
      )
      n_unl_expected <- tryCatch(
        nrow(safe_read_gpkg(gpkg_pools_out, "scoring_pool", "scoring_pool_count_check")),
        error = function(e) NA_integer_
      )
      n_train_feat <- feature_layer_n(features_gpkg, "train_features")
      n_unl_feat   <- feature_layer_n(features_gpkg, "scoring_features")

      if (!all(c("train_features", "scoring_features") %in% feat_layers)) {
        should_rebuild_features <- TRUE
        rebuild_reasons <- c(
          rebuild_reasons,
          sprintf(
            "features_geometry.gpkg no contiene todas las capas requeridas (layers actuales: %s)",
            paste(feat_layers, collapse = ", ")
          )
        )
      }
      
      if (!is.na(n_train_expected) && !is.na(n_train_feat) && n_train_expected != n_train_feat) {
        should_rebuild_features <- TRUE
        rebuild_reasons <- c(rebuild_reasons, sprintf("train_features tiene %d filas y train_with_folds %d", n_train_feat, n_train_expected))
      }
      if (!is.na(n_unl_expected) && !is.na(n_unl_feat) && n_unl_expected != n_unl_feat) {
        should_rebuild_features <- TRUE
        rebuild_reasons <- c(rebuild_reasons, sprintf("scoring_features tiene %d filas y scoring_pool %d", n_unl_feat, n_unl_expected))
      }
      
      if (isTRUE(should_rebuild_features)) {
        msg("STEP B3 - features_geometry.gpkg existe pero esta obsoleto. Se reconstruye.")
        for (rr in unique(rebuild_reasons)) msg("  * %s", rr)
        safe_remove_features_gpkg(features_gpkg)
      } else {
      msg("STEP B3 - features_geometry.gpkg ya existe, SKIP: %s", features_gpkg)
      
      }
    }
    
    if (!file.exists(features_gpkg)) {
      
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
      
      stopifnot(!is.null(train_with_folds_gpkg), file.exists(train_with_folds_gpkg))
      stopifnot(file.exists(gpkg_pools_out))
      
      msg("STEP B3.1 - Read train_with_folds + scoring_pool")
      train_folds <- time_step("B3.1 Read train_with_folds", {
        safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "train_with_folds")
      })
      
      scoring_candidates <- time_step("B3.2 Read scoring_pool", {
        safe_read_gpkg(gpkg_pools_out, "scoring_pool", "scoring_pool")
      })
      
      hotspots_sf <- hotspots_sf_base
      if (nrow(hotspots_sf) > 0) {
        if (is.na(sf::st_crs(hotspots_sf))) {
          msg("WARNING: hotspots_sf has NA CRS; disabling hotspots.")
          hotspots_sf <- hotspots_sf[0, ]
        } else if (sf::st_crs(hotspots_sf) != sf::st_crs(train_folds)) {
          hotspots_sf <- sf::st_transform(hotspots_sf, sf::st_crs(train_folds))
        }
      }
      
      use_hotspots_flag <- nrow(hotspots_sf) > 0
      
      msg("STEP B3.3 - extract_features (train + unlabeled)")
      
      msg("STEP B3.2b - Check raster geometry before extract_features")
      time_step("B3.2b Check raster geometry", {
        raster_list <- list(
          rbr_summer = rbr_summer,
          doy_post   = doy_post,
          rbr_aw     = rbr_aw,
          dem_r      = dem_r,
          slope_r    = slope_r,
          corine_r   = corine_r
        )
        
        for (nm in names(raster_list)) {
          rr <- raster_list[[nm]]
          ok <- terra::compareGeom(rr, rbr_summer, stopOnError = FALSE)
          
          msg("Raster check: %s | ok=%s | nrow=%s ncol=%s | res=(%s,%s) | origin=(%s,%s)",
              nm, ok, nrow(rr), ncol(rr),
              terra::res(rr)[1], terra::res(rr)[2],
              terra::origin(rr)[1], terra::origin(rr)[2])
          
          if (!isTRUE(ok)) {
            stop(sprintf("Raster geometry mismatch before extract_features: %s does not match rbr_summer", nm))
          }
        }
        
        NULL
      })
      
      feat <- time_step("B3.3 extract_features", {
        extract_features(
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
          
          ecoregions = NULL,
          eco_id_col = "EnZ_name",
          
          hotspots   = hotspots_sf,
          frp_col    = "frp",
          conf_col   = "confidence",
          year_col   = "year",
          month_col  = "month",
          hs_start_month = 6,
          hs_end_month   = 10,
          hi_conf_thr    = 0.8,
          
          year_target = target_year,
          hs_year_min_available = 2000,
          hs_buffer_m = 1000,
          hs_missing_value = -9999,
          hs_use_season_filter = TRUE,
          
          cor_groups = cor_groups,
          
          use_doy      = TRUE,
          use_aw       = TRUE,
          use_nbr      = FALSE,
          use_hotspots = use_hotspots_flag,
          use_ecoregions = FALSE,
          
          max_cells_in_memory      = NULL,
          return_features          = TRUE,
          return_features_geometry = TRUE,
          save_features_dir        = dirs$`03_FEATURES`,
          save_features_format     = "rds",
          save_features_gpkg       = TRUE,
          train_features_basename  = "train_features",
          scoring_features_basename = "scoring_features",
          train_features_layer     = "train_features",
          scoring_features_layer   = "scoring_features",
          verbose                  = TRUE
        )
      })
      
      stopifnot(file.exists(features_gpkg))
      msg("DONE - features saved: %s", features_gpkg)
    }
  }
  
  # ============================== 9) MODELING ============================
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
    
    # Bug 5 (0.3.0): character is the canonical type for `class` and
    # other label-like columns. Conversion to factor happens only
    # inside `build_design_matrix_patches`, immediately before the
    # `model.matrix`/`sparse.model.matrix` call.
    labelled <- sf::st_drop_geometry(train_features)
    burned_like <- sf::st_drop_geometry(scoring_features)
    
    labelled_df <- labelled
    
    msg("STEP C1 - Build XGB params")
    params <- time_step("C1 Build params", {
      n_pos <- sum(labelled$class == "burned",   na.rm = TRUE)
      n_neg <- sum(labelled$class == "unburned", na.rm = TRUE)
      spw   <- if (n_pos > 0) n_neg / n_pos else 1
      
      list(
        booster = "gbtree",
        objective = "binary:logistic",
        # OOF repair (HANDOFF Section 87): eval_metric is now a vector.
        # In the installed xgboost build the FIRST metric in eval_metric is
        # used for early stopping (verified empirically on the 1985 fold:
        # "aucpr" alone -> best_iter=1; "logloss" alone -> best_iter=137;
        # c("logloss","aucpr") -> best_iter=137; c("aucpr","logloss") ->
        # best_iter=1).  We therefore put logloss FIRST so early stopping
        # is driven by logloss (which does not saturate) while AUCPR
        # remains visible in the per-iteration eval_log for diagnostics.
        eval_metric = c("logloss", "aucpr"),
        eta = 0.06,
        max_depth = 5,
        subsample = 0.85,
        colsample_bytree = 0.75,
        min_child_weight = 5,
        gamma = 0,
        scale_pos_weight = spw
      )
    })
    
    msg("STEP C2 - DM -> OOF (diagnostic)")
    pipe1 <- time_step("C2 run_dm_oof_pipeline", {
      stopifnot(!is.null(train_with_folds_gpkg) && file.exists(train_with_folds_gpkg))
      stopifnot(!is.null(blocks_gpkg) && file.exists(blocks_gpkg))
      
      run_dm_oof_pipeline(
        labelled    = labelled,
        burned_like = burned_like,
        labelled_df = labelled_df,
        params      = params,
        
        result_dir  = result_dir,
        target_year = target_year,
        prefix      = prefix_oof,
        
        id_cols     = c("fire_uid","class","source","poly_id","block_id","fold_rep1","fold_rep2"),
        # 0.4.0 (Agent H): the orchestrator does not pass drop_regex.
        # The OOF wrapper applies the canonical whitelist filter
        # (`.supervised_feature_cols`) directly; drop_regex is a no-op
        # since 0.4.0.
        cat_cols    = c("eco_major"),
        hs_n_col    = "hs_used_n",
        hs_conf_col = "hs_conf_mean",
        hs_frp_col  = "hs_frp_max",
        median_from = "labelled",
        save_dir_dm = dirs$`04_MATRIX`,
        save_prefix = paste0(target_year, "_", scenario, "_patch"),
        overwrite   = TRUE,
        # 0.5.0: pass the active whitelist override and per-feature
        # weights down to the OOF stage so it trains under the same
        # conditions as the final model (KB1/KB2 symmetry).
        feature_whitelist_override = feature_whitelist_override,
        feature_weights            = feature_weights,

        labelled_gpkg = train_with_folds_gpkg,
        labelled_layer = "train_with_folds"
      )
    })
    
    print(pipe1$files)
    
    msg("STEP C3 - Train final model + export final map (wrapper)")
    pipe2 <- time_step("C3 run_train_final_model_and_export_final_map", {
      
      oof_agg_csv <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_oof_agg.csv"))
      oof_summary_gpkg <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_labeled_oof_summary.gpkg"))
      stopifnot(file.exists(oof_agg_csv))
      stopifnot(file.exists(oof_summary_gpkg))
      
      run_train_final_model_and_export_final_map(
        result_dir     = result_dir,
        qa             = oof_agg_csv,
        labelled_gpkg  = gpkg_features,
        labelled_layer = "train_features",
        out_dir        = dirs$`07_FINAL_MODEL_V2`,
        prefix         = prefix,
        overwrite      = TRUE,
        verbose        = TRUE,
        
        qa_labelled_gpkg  = oof_summary_gpkg,
        qa_labelled_layer = "labeled_oof_summary",
        
        labelled_features_gpkg  = gpkg_features,
        labelled_features_layer = "train_features",
        
        model_rds  = file.path(dirs$`07_FINAL_MODEL_V2`, paste0(prefix, "_final_model.rds")),
        recipe_rds = file.path(dirs$`07_FINAL_MODEL_V2`, paste0(prefix, "_recipe.rds")),
        
        unlabeled_gpkg  = gpkg_features,
        unlabeled_layer = "scoring_features",
        
        out_score_dir = dirs$`08_SCORED`,
        out_map_dir   = dirs$`09_FINAL_MAP`,
        
        export_burned_like = TRUE,
        preyear_overlap_threshold = CURRENTYEAR_PREYEAR_OVERLAP_THR,
        hotspot_density_threshold = CURRENTYEAR_HOTSPOT_DENSITY_THR,
        temporal_penalty_floor = CURRENTYEAR_TEMPORAL_PENALTY_FLOOR,

        # Phase B pass-through. Defaults reproduce historical behaviour.
        contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
        spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
        random_to_burned_ratio                 = random_to_burned_ratio,
        otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
        # 0.5.0: pass the active whitelist override and per-feature
        # weights down to the FINAL model. Combined with the OOF
        # call above, this guarantees both stages share the same
        # feature space and weight vector.
        feature_whitelist_override             = feature_whitelist_override,
        feature_weights                        = feature_weights
      )
    })

    print(pipe2$m2$files)
    print(pipe2$out$files)
    
    msg("STEP C4 - Append supervised high-confidence burned to multiyear registry")
    registry_res <- time_step("C4 append_supervised_registry", {
      append_supervised_registry(
        burned_like_sf = pipe2$out$burned_like_scored,
        registry_path = resolved_registry_path,
        target_year = target_year,
        scenario = scenario,
        prob_thr = REGISTRY_PROB_THR,
        top_fraction = REGISTRY_TOP_FRACTION,
        verbose = TRUE
      )
    })
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
    common_unb_dir = if (identical(UNB_STRATEGY, "all_sources")) legacy_unb_root_dir else deterministic_common_unb_dir,
    registry_path = resolved_registry_path
  ))
}

if (FALSE) {
  # ============================== 6) RUN ALL ===============================
  batch_results <- list()
  
  for (yy in years_to_run) {
    for (scn in scenarios_to_run) {
      
      run_id <- paste(yy, scn, sep = "_")
      
      res <- tryCatch({
        out <- run_supervised_pipeline(target_year = yy, scenario = scn,
                                       min_burned_pool_n = MIN_BURNED_POOL_N)
        data.frame(
          year        = yy,
          scenario    = scn,
          status      = "OK",
          error       = NA_character_,
          result_dir  = if (!is.null(out$result_dir)) out$result_dir else NA_character_,
          timing_csv  = if (!is.null(out$timing_csv)) out$timing_csv else NA_character_,
          unb_out_gpkg = if (!is.null(out$unb_out_gpkg)) out$unb_out_gpkg else NA_character_,
          common_unb_dir = if (!is.null(out$common_unb_dir)) out$common_unb_dir else NA_character_,
          registry_path = if (!is.null(out$registry_path)) out$registry_path else NA_character_,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message(sprintf("ERROR in year=%s scenario=%s -> %s", yy, scn, e$message))
        err_unb_path <- if (identical(UNB_STRATEGY, "all_sources")) {
          file.path(
            data_base, "Results", yy, result_name, "SUPERVISED", scn,
            "_LEGACY_UNBURNED", "05_UNBURNED",
            sprintf("%d_%s_legacy_unburned.gpkg", yy, scn)
          )
        } else {
          file.path(
            data_base, "Results", yy, result_name, "DETERMINISTIC", scn,
            "UNBURNED", sprintf("%d_%s_unburned.gpkg", yy, scn)
          )
        }
        err_common_unb_dir <- if (identical(UNB_STRATEGY, "all_sources")) {
          file.path(data_base, "Results", yy, result_name, "SUPERVISED", scn, "_LEGACY_UNBURNED")
        } else {
          file.path(data_base, "Results", yy, result_name, "DETERMINISTIC", "_COMMON_UNBURNED")
        }
        data.frame(
          year        = yy,
          scenario    = scn,
          status      = "ERROR",
          error       = e$message,
          result_dir  = file.path(data_base, "Results", yy, result_name, "SUPERVISED", scn),
          timing_csv  = NA_character_,
          unb_out_gpkg = err_unb_path,
          common_unb_dir = err_common_unb_dir,
          registry_path = if (!is.null(RBR_TRAINING_REGISTRY_PATH) &&
                              nzchar(RBR_TRAINING_REGISTRY_PATH)) {
            RBR_TRAINING_REGISTRY_PATH
          } else {
            get_fire_mapping_registry_paths(
              data_base = data_base,
              scenario = scn,
              registry_basename = "RBR_TRAINING_REGISTRY.gpkg",
              registry_layer = "burned_high_conf_registry"
            )$registry_path
          },
          stringsAsFactors = FALSE
        )
      })
      
      batch_results[[run_id]] <- res
    }
  }
  
  batch_results <- dplyr::bind_rows(batch_results)
  print(batch_results)
  
  batch_out_dir <- file.path(data_base, "Results", result_name)
  dir.create(batch_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  batch_csv <- file.path(batch_out_dir, "supervised_batch_results.csv")
  write.csv(batch_results, batch_csv, row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("ALL RUNS FINISHED\n")
  cat(sprintf("Batch results saved: %s\n", batch_csv))
  cat("============================================================\n")
}
