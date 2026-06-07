#' Build a configuration object for the supervised burned-area workflow
#'
#' @description
#' Public entry point of the one-year supervised workflow. Collects the
#' deterministic decision layer, raster support layers, optional hotspot and
#' reference inputs, output location, and advanced options into a single
#' normalized configuration object of
#' class `otsufire_supervised_burned_config`.
#'
#' Contract source: `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `SUPERVISED_ONEYEAR_INPUTS_FINAL.csv`, `SUPERVISED_ONEYEAR_OUTPUTS_ROUTES.csv`.
#'
#' @param scenario Character scalar. Methodological preset. One of
#'   `"balanced"` (default), `"original"`, `"lax"`, `"restrictive"`.
#' @param internal_decisions sf POLYGON layer or a GPKG path. Deterministic
#'   decision layer produced for the same year and scenario. Required.
#' @param change_index SpatRaster or raster path. Main annual change-index
#'   raster. Required.
#' @param delayed_change_index SpatRaster or raster path. Optional delayed
#'   (autumn-winter) change-index raster for persistence-style features.
#'   RUN input (§N+25): when supplied it is CONSUMED by the orchestrator as
#'   the `rbr_aw` feature layer; when `NULL` it defaults to the convention
#'   path `<composite_base>/Autumn/mean_mean_<target_year>_mosaic.tif`.
#' @param hotspots sf POINT layer or path. Optional hotspot layer for the
#'   target year.
#' @param reference_burned_map sf POLYGON / raster / path. Optional external
#'   burned reference used by `validate_fire_maps()` on thresholded
#'   supervised outputs.
#' @param peninsula_shapefile sf / SpatVector / path. Optional Iberian
#'   peninsula border polygon. RUN input (§N+25): consumed by the legacy
#'   unburned Otsu builder only when `options$legacy_otsu_mode` crops CORINE
#'   to the peninsula (modes `"corine"` / `"corine_ecoregion"`). When `NULL`
#'   it defaults to the convention path
#'   `<data_base>/Borders/Iberian_peninsula.shp`.
#' @param topo SpatRaster / path. Optional two-band topography raster
#'   (band1 = DEM, band2 = slope). RUN input (§N+25). When `NULL` it defaults
#'   to the convention path `<data_base>/Topography/elevation_slope.tif`. The
#'   DEM and slope live in a single file by convention, so they are exposed as
#'   one `topo` input rather than two separate paths.
#' @param corine_raster SpatRaster / path. Optional CORINE land-cover raster
#'   for the target year's CORINE epoch. RUN input (§N+25). When `NULL` it
#'   defaults to the convention path
#'   `<data_base>/Corine_Masks/CLC_<corine_year>_peninsula.tif`.
#' @param burnable_mask SpatRaster / path. Optional binary burnable-area mask
#'   for the target year's CORINE epoch. RUN input (§N+25): consumed by both
#'   all_sources unburned builders. When `NULL` it defaults to the convention
#'   path
#'   `<data_base>/Corine_Masks/burneable_mask_binary_corine_<corine_year>_ETRS89.tif`.
#' @param target_year Integer scalar. Target year. Required.
#' @param output_dir Character. Root output directory. Defaults to
#'   `tempdir()`.
#' @param run_name Character. Stable identifier for the supervised run.
#'   Defaults to `"supervised_burned_map"`.
#' @param min_burned_pool_n Integer. Sparse-year guard. The pipeline aborts
#'   before training when fewer than this many keep-class polygons remain
#'   after QA. Default `5L` (HANDOFF Section 51).
#' @param options Named list of advanced options. Recognized keys include
#'   `legacy_*`, `currentyear_*`, `unb_verbose`, `result_name`,
#'   `data_base`, plus overrides forwarded to the engine. The
#'   `negative_pool_policy` / `otsu_unburned_generation` keys were removed in
#'   §N+26 (2026-06-05): the negative-pool policy is always `"all_sources"` and
#'   is no longer user-settable. Passing `negative_pool_policy = "all_sources"`
#'   is accepted as a no-op; any other stale value (including
#'   `"deterministic_direct"`) raises a clear error (guarded, never silently
#'   ignored). The returned cfg still carries `negative_pool_policy =
#'   "all_sources"` as an informational constant. The
#'   `engine_root`, `supervised_engine_root` and `scripts_root` keys were
#'   removed on 2026-06-05 (the engine is now fully in-package).
#'   B2 (2026-06-05) additionally exposes the shared deterministic-decision
#'   negative-pool knobs `unb_excl_buffer_m` (500), `unb_n_random_cells`
#'   (1500), `unb_random_rbr_q` (0.50), `unb_random_seed` (42),
#'   `unb_random_patch_size_cells` (3); the supervised run-prefix names
#'   `prefix_oof_base` (`"patch"`) and `prefix_base` (`"patch_certified"`);
#'   and the remaining legacy Otsu knobs `legacy_min_otsu_threshold_value`
#'   (0), `legacy_min_pixels` (8), `legacy_buffers_m` (90), `legacy_core_thr`
#'   (0.60), `legacy_alpha_boost` (0.25), `legacy_min_base_boost` (0.35),
#'   `legacy_dist_power` (1), `legacy_keep_hi` (0.45), `legacy_drop_lo`
#'   (0.15), `legacy_excl_buffer_m` (0), `legacy_min_area_ha` (0). Every
#'   default matches the value previously hardcoded in the orchestrator, so
#'   passing nothing reproduces the historical run byte-for-byte.
#'
#'   Methodologically significant keys (HANDOFF §N+28 "Tier 1" — the
#'   validity-affecting knobs worth recording for reproducibility and the
#'   natural ablation candidates): the current-year temporal-adjustment
#'   thresholds `currentyear_preyear_overlap_thr` (0.70, previous-year
#'   overlap gate), `currentyear_hotspot_density_thr` (0.001, hotspot-density
#'   gate) and `currentyear_temporal_penalty_floor` (0.10, temporal-penalty
#'   floor); the negative-pool seeds `unb_random_seed` (42) and
#'   `legacy_random_seed` (42, the legacy Otsu sampler seed — set BOTH for
#'   full reproducibility from the cfg alone); and the Otsu confidence
#'   thresholds `legacy_keep_hi` (0.45) / `legacy_drop_lo` (0.15). The
#'   seventh Tier-1 knob, the sparse-year guard `min_burned_pool_n`, is the
#'   top-level `min_burned_pool_n` argument (documented above), NOT an
#'   `options` key — kept top-level by design for discoverability.
#'
#' @return S3 object of class `otsufire_supervised_burned_config`.
#'
#' @family workflow
#' @export
build_supervised_burned_config <- function(
    scenario = c("balanced", "original", "lax", "restrictive"),
    internal_decisions,
    change_index,
    delayed_change_index = NULL,
    hotspots = NULL,
    reference_burned_map = NULL,
    target_year,
    output_dir = tempdir(),
    run_name = "supervised_burned_map",
    min_burned_pool_n = 5L,
    # §N+25 (2026-06-05): supervised-RUN inputs wired to the config. All
    # default NULL; when NULL they fall back to the EXACT convention path
    # previously built by the orchestrator / unburned helpers, so passing
    # nothing reproduces the historical run byte-for-byte. The
    # VALIDATION-only inputs (strata / study-area mask / EFFIS reference)
    # are NOT part of the supervised run; they belong to the separate
    # validate_fire_maps() call, exactly like the deterministic phase.
    peninsula_shapefile = NULL,
    topo = NULL,
    corine_raster = NULL,
    burnable_mask = NULL,
    options = list()
) {
  scenario <- match.arg(scenario)

  if (missing(target_year) || is.null(target_year)) {
    stop("'target_year' is required.", call. = FALSE)
  }
  target_year <- suppressWarnings(as.integer(target_year)[1L])
  if (is.na(target_year) || target_year < 1900L || target_year > 2100L) {
    stop("'target_year' must be an integer in [1900, 2100].", call. = FALSE)
  }

  if (missing(internal_decisions) || is.null(internal_decisions)) {
    stop("'internal_decisions' is required.", call. = FALSE)
  }
  if (missing(change_index) || is.null(change_index)) {
    stop("'change_index' is required.", call. = FALSE)
  }

  if (!is.character(run_name) || length(run_name) != 1L || !nzchar(run_name)) {
    stop("'run_name' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("'output_dir' must be a single non-empty character string.", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  min_burned_pool_n <- suppressWarnings(as.integer(min_burned_pool_n)[1L])
  if (is.na(min_burned_pool_n) || min_burned_pool_n < 0L) {
    stop("'min_burned_pool_n' must be a non-negative integer.", call. = FALSE)
  }

  if (!is.list(options)) {
    stop("'options' must be a list().", call. = FALSE)
  }
  if (length(options) > 0L && (is.null(names(options)) || any(!nzchar(names(options))))) {
    stop("'options' must be a named list (all elements must have names).", call. = FALSE)
  }

  # ---- §N+25: resolve the supervised-RUN input paths --------------------
  # Centralized convention construction (`.of_supervised_convention_paths`)
  # so the cfg, the orchestrator and the unburned helpers agree on EXACTLY
  # the same default paths. The convention requires `data_base` (and, for
  # delayed_change_index, `composite_base` + `result_name`), which live in
  # `options`. When `data_base` is absent (e.g. an ahead-of-batch config or a
  # unit test that does not exercise the run), the convention path cannot be
  # built; in that case the user-supplied value is normalized as-is and NULL
  # stays NULL, deferring resolution to runtime exactly as before.
  conv <- .of_supervised_convention_paths(
    data_base      = options$data_base,
    composite_base = options$composite_base,
    result_name    = options$result_name %||% "Min_Min",
    target_year    = target_year
  )

  delayed_change_index <- delayed_change_index %||% conv$delayed_change_index
  peninsula_shapefile  <- peninsula_shapefile  %||% conv$peninsula_shapefile
  topo                 <- topo                 %||% conv$topo
  corine_raster        <- corine_raster        %||% conv$corine_raster
  burnable_mask        <- burnable_mask        %||% conv$burnable_mask

  inputs <- list(
    internal_decisions   = .of_normalize_input_spec(internal_decisions,
                                                     "internal_decisions",
                                                     allow_null = FALSE),
    change_index         = .of_normalize_input_spec(change_index,
                                                     "change_index",
                                                     allow_null = FALSE),
    delayed_change_index = .of_normalize_input_spec(delayed_change_index,
                                                     "delayed_change_index",
                                                     allow_null = TRUE),
    hotspots             = .of_normalize_input_spec(hotspots, "hotspots",
                                                     allow_null = TRUE),
    reference_burned_map = .of_normalize_input_spec(reference_burned_map,
                                                     "reference_burned_map",
                                                     allow_null = TRUE),
    # §N+25 RUN inputs wired to cfg (convention defaults applied above).
    peninsula_shapefile  = .of_normalize_input_spec(peninsula_shapefile,
                                                     "peninsula_shapefile",
                                                     allow_null = TRUE),
    topo                 = .of_normalize_input_spec(topo, "topo",
                                                     allow_null = TRUE),
    corine_raster        = .of_normalize_input_spec(corine_raster,
                                                     "corine_raster",
                                                     allow_null = TRUE),
    burnable_mask        = .of_normalize_input_spec(burnable_mask,
                                                     "burnable_mask",
                                                     allow_null = TRUE)
  )

  # ---- §N+25: fail-fast validation of the newly-wired RUN inputs --------
  # When a wired RUN input resolves to an on-disk PATH spec, assert the file
  # exists NOW (config-construction is <1s) instead of 5-10 min into the run.
  # This only fires for path specs; in-memory objects and NULL (deferred) are
  # skipped and re-validated downstream. The convention default for a
  # `data_base`-less config is NULL, so such configs are not blocked here.
  # `internal_decisions` / `change_index` / `hotspots` keep their historical
  # config-time behaviour (existence enforced downstream, not here). The
  # VALIDATION-only inputs (strata / mask / EFFIS reference) are intentionally
  # absent here: they belong to the separate validate_fire_maps() call.
  for (.nm in c("delayed_change_index", "peninsula_shapefile", "topo",
                "corine_raster", "burnable_mask")) {
    .of_check_input_file(inputs[[.nm]], .nm)
  }

  # BUG 3 Phase 1a (2026-06-05): the supervised engine is fully in-package, so
  # the engine-root / scripts-root discovery here was inert metadata that only
  # fed the now-removed external PROJECT_PATHS.R fallback. Removed the
  # `engine_root`, `supervised_engine_root` and `scripts_root` config fields,
  # their getwd()-walk finders (`.of_find_supervised_engine_root()`,
  # `.of_find_supervised_scripts_root()`) and the `.of_find_engine_root()`
  # call from the supervised path. `.of_find_engine_root()` itself is retained
  # in deterministic-config.R for the deterministic stage.

  # External tool paths (AS01, OtsuFire 0.3.0). Mirrors the deterministic
  # stage convention from `deterministic-config.R`. All NULL by default;
  # the legacy unburned pipeline (`all_sources` mode) errors loudly if
  # the path it needs is missing. Callers may pass either flat options
  # (e.g. options$python_exe) or a nested options$tool_paths list.
  tool_paths <- list(
    python_exe             = options$python_exe             %||%
                              options$tool_paths$python_exe,
    gdal_polygonize_script = options$gdal_polygonize_script %||%
                              options$tool_paths$gdal_polygonize_script,
    gdalwarp_path          = options$gdalwarp_path          %||%
                              options$tool_paths$gdalwarp_path,
    ogr2ogr_exe            = options$ogr2ogr_exe            %||%
                              options$tool_paths$ogr2ogr_exe
  )

  # Negative-pool policy (§N+26, 2026-06-05): FULLY REMOVED as a user knob.
  # The supervised pipeline is ALWAYS `all_sources` (deterministic drops +
  # random burnable background + Otsu current-year patches). The resolver
  # `.of_resolve_negative_pool_policy()`, the deprecated
  # `otsu_unburned_generation` alias, and the `deterministic_direct` operational
  # mode are gone. HONEST GUARD: a config that carries a stale
  # `negative_pool_policy` / `otsu_unburned_generation` option errors clearly
  # rather than being silently ignored. The historical canonical value
  # `"all_sources"` is still accepted as a no-op so existing scripts keep
  # working. `neg_pool` is then pinned to the fixed constant below.
  for (.stale_key in c("negative_pool_policy", "negative_pool_strategy",
                       "otsu_unburned_generation")) {
    if (!is.null(options[[.stale_key]]) &&
        !identical(options[[.stale_key]], "all_sources")) {
      stop(
        "`negative_pool_policy`/`otsu_unburned_generation` was removed in ",
        "2026-06; the supervised pipeline is always 'all_sources'. ",
        "Remove this option (you passed `", .stale_key, " = ",
        as.character(options[[.stale_key]])[1L], "`).",
        call. = FALSE
      )
    }
  }
  neg_pool <- "all_sources"

  output_routes <- .of_build_supervised_output_routes(
    output_dir = output_dir, target_year = target_year,
    run_name = run_name, scenario = scenario
  )

  # §N+27 RESOLVED (2026-06-05): the supervised burned-like registry was an
  # abandoned research line. The `burned_like_registry_path` parameter, its
  # default-path resolution (via .resolve_registry_path()), and the cfg field
  # are removed. The deterministic phase keeps its own registry READER
  # (build_rbr_keep_pool_from_registry); .resolve_registry_path() in
  # R/utils-registry.R is preserved for that side.

  cfg <- list(
    scenario                 = scenario,
    target_year              = target_year,
    inputs                   = inputs,
    run_name                 = run_name,
    output_dir               = output_dir,
    output_routes            = output_routes,
    # §N+26: informational constant only. No longer user-settable; always
    # "all_sources". Kept on the cfg so saveRDS(cfg) still records the policy
    # and downstream introspection / print() stays self-documenting.
    negative_pool_policy     = neg_pool,
    min_burned_pool_n        = min_burned_pool_n,
    tool_paths               = tool_paths,
    options                  = options
  )
  class(cfg) <- c("otsufire_supervised_burned_config", "list")
  cfg
}

#' @export
print.otsufire_supervised_burned_config <- function(x, ...) {
  cat("<otsufire_supervised_burned_config>\n")
  cat("  scenario        :", x$scenario, "\n")
  cat("  target_year     :", x$target_year, "\n")
  cat("  run_name        :", x$run_name, "\n")
  cat("  output_dir      :", x$output_dir, "\n")
  cat("  negative_pool   :", x$negative_pool_policy, "\n")
  cat("  min_burned_pool :", x$min_burned_pool_n, "\n")
  cat("  inputs          :\n")
  for (nm in names(x$inputs)) {
    v <- x$inputs[[nm]]
    cat(sprintf("    - %-22s: %s\n", nm,
                if (is.null(v)) "<NULL>" else
                  sprintf("%s [%s]",
                          if (!is.null(v$path)) v$path else "<in-memory>",
                          v$type)))
  }
  invisible(x)
}

# --- internal helpers --------------------------------------------------

#' Centralized supervised-RUN input-path conventions (§N+25)
#'
#' Single source of truth for the by-convention default paths of the
#' supervised-RUN inputs. `build_supervised_burned_config()` calls this to
#' fill the cfg defaults; the orchestrator and the unburned helpers fall back
#' to the SAME constructions when the cfg does not carry an explicit path, so
#' the cfg-default and the runtime-fallback always agree (byte-identical).
#'
#' Returns a named list whose entries are character paths, or NULL when
#' `data_base` (resp. `composite_base`) is not available to build them. The
#' CORINE epoch is resolved with `.of_corine_year_for()` (the same mapping the
#' orchestrator's `get_corine_year()` uses).
#'
#' @keywords internal
#' @noRd
.of_supervised_convention_paths <- function(data_base, composite_base,
                                            result_name, target_year) {
  has_db <- is.character(data_base) && length(data_base) == 1L &&
            nzchar(data_base)
  has_cb <- is.character(composite_base) && length(composite_base) == 1L &&
            nzchar(composite_base)
  corine_year <- .of_corine_year_for(target_year)

  list(
    # AW (delayed) change-index mosaic lives under composite_base/Autumn.
    delayed_change_index = if (has_cb) {
      file.path(composite_base, "Autumn",
                paste0("mean_mean_", target_year, "_mosaic.tif"))
    } else NULL,
    peninsula_shapefile = if (has_db) {
      file.path(data_base, "Borders", "Iberian_peninsula.shp")
    } else NULL,
    topo = if (has_db) {
      file.path(data_base, "Topography", "elevation_slope.tif")
    } else NULL,
    corine_raster = if (has_db) {
      file.path(data_base, "Corine_Masks",
                paste0("CLC_", corine_year, "_peninsula.tif"))
    } else NULL,
    burnable_mask = if (has_db) {
      file.path(data_base, "Corine_Masks",
                paste0("burneable_mask_binary_corine_", corine_year,
                       "_ETRS89.tif"))
    } else NULL
  )
}

#' CORINE epoch for a target year (§N+25 convention helper).
#'
#' Mirrors the orchestrator's `get_corine_year()` and the unburned helpers'
#' `get_corine_year_unb*()` so every site resolves the same CORINE raster /
#' burnable-mask filename. Returns a character year string.
#'
#' @keywords internal
#' @noRd
.of_corine_year_for <- function(y) {
  y <- suppressWarnings(as.integer(y)[1L])
  if (!is.na(y) && y >= 1984L && y <= 1999L) "1990"
  else if (!is.na(y) && y <= 2005L) "2000"
  else if (!is.na(y) && y <= 2011L) "2006"
  else if (!is.na(y) && y <= 2017L) "2012"
  else "2018"
}

# BUG 3 Phase 1a (2026-06-05): removed the dead supervised getwd()-walk
# finders `.of_find_supervised_engine_root()` and
# `.of_find_supervised_scripts_root()`. They only populated inert
# `supervised_engine_root` / `scripts_root` config metadata that fed the
# now-removed external PROJECT_PATHS.R fallback. The supervised engine is
# in-package; nothing reads these roots anymore.

# §N+26 (2026-06-05): `.of_resolve_negative_pool_policy()` REMOVED. The
# supervised pipeline no longer resolves a negative-pool policy — it is always
# `all_sources`. The user-facing knob, the `otsu_unburned_generation` alias and
# the `deterministic_direct` mode are gone; a stale option now triggers the
# honest guard in `build_supervised_burned_config()` instead of a silent
# redirect. (This is unrelated to the D4a empty-Otsu-pool robustness fallback
# in internal-sup-unburned-legacy.R, which is retained.)

#' @keywords internal
#' @noRd
.of_build_supervised_output_routes <- function(output_dir, target_year,
                                                run_name, scenario) {
  base <- file.path(output_dir, as.character(target_year), run_name,
                    "SUPERVISED", scenario)
  prefix_oof <- sprintf("%d_%s_patch", target_year, scenario)
  prefix     <- sprintf("%d_%s_patch_certified", target_year, scenario)
  # Bug 1 (0.3.0): folder names match the operational runtime
  # (`07_FINAL_MODEL_V2`, `08_SCORED`, `09_FINAL_MAP`,
  # `11_CONSISTENCY_CHECKS`). Previously the config advertised
  # `06_FINAL_MODEL`/`07_SCORED`/`08_FINAL_MAP`/`09_CONSISTENCY`,
  # which did not exist on disk — public callers reading
  # `config$output_routes` directly were getting nonexistent paths.
  list(
    base                = base,
    pools_dir           = file.path(base, "01_POOLS"),
    folds_dir           = file.path(base, "02_FOLDS"),
    features_dir        = file.path(base, "03_FEATURES"),
    matrix_dir          = file.path(base, "04_MATRIX"),
    oof_dir             = file.path(base, "05_OOF"),
    final_model_dir     = file.path(base, "07_FINAL_MODEL_V2"),
    scored_dir          = file.path(base, "08_SCORED"),
    final_map_dir       = file.path(base, "09_FINAL_MAP"),
    consistency_dir     = file.path(base, "11_CONSISTENCY_CHECKS"),
    logs_dir            = file.path(base, "99_LOGS"),
    pools_gpkg          = file.path(base, "01_POOLS",
                                     sprintf("%d_%s_pools.gpkg",
                                             target_year, scenario)),
    train_with_folds_gpkg = file.path(base, "02_FOLDS",
                                       sprintf("%d_train_with_folds_5000m.gpkg",
                                               target_year)),
    features_geometry_gpkg = file.path(base, "03_FEATURES",
                                        "features_geometry.gpkg"),
    oof_agg_csv         = file.path(base, "05_OOF",
                                     paste0(prefix_oof, "_oof_agg.csv")),
    final_model_rds     = file.path(base, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_final_model.rds")),
    final_map_gpkg      = file.path(base, "09_FINAL_MAP",
                                     paste0(prefix, "_final_map.gpkg")),
    burned_like_gpkg    = file.path(base, "09_FINAL_MAP",
                                     paste0(prefix, "_burned_like_scored.gpkg")),
    timing_csv          = file.path(base, "99_LOGS",
                                     paste0(target_year, "_timing_steps.csv"))
  )
}
