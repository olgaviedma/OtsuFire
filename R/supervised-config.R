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
#'   raster (band1 = summer RBR, band2 = DOY_post). Required. Gate 1B
#'   (2026-06-07): CONSUMED from `cfg$inputs$change_index` end-to-end by the
#'   orchestrator, the pool builder (legacy-Otsu severity raster) and the
#'   feature extractor — it is NO LONGER reconstructed by the historical
#'   `MinMin_<year>_mosaic_res90m.tif` filename convention. The convention path
#'   survives only as a fallback for the standalone-script path (no cfg) or an
#'   in-memory spec, so existing scripts stay byte-identical.
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
#' @param ecoregion_shapefile sf / SpatVector / path. Optional ecoregion border
#'   polygon (Olson ecoregions). RUN input. Gate 1B PIECE 2 (2026-06-07):
#'   consumed by the legacy unburned Otsu builder ONLY in the `ecoregion` /
#'   `corine_ecoregion` modes (the default mode is `burnable_only`, which never
#'   reads it). Previously a `<data_base>/Ecoregion/ecoregiones_olson.shp`
#'   hardcode inside the builder; now a configurable cfg input. When `NULL` it
#'   defaults to that convention path. NOTE: PIECE 2 only makes the path
#'   configurable; the ecoregion Otsu branch is removed later in PIECE 4.
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
#' @param nrounds_max,early_stop Integer or `NULL`. XGBoost max boosting rounds
#'   and early-stopping patience. `NULL` uses the canonical defaults
#'   (`4000` / `80`). Stored in `cfg$train_control`.
#' @param oof_seed_base,final_sampling_seed,final_seed Integer or `NULL`. The
#'   three distinct RNG seeds (OOF block-CV; FINAL negative-pool sampling; FINAL
#'   train/val split + xgboost). `NULL` uses the canonical `42` for each. Stored
#'   in `cfg$train_control$seeds`.
#' @param val_frac Numeric in (0, 1) or `NULL`. FINAL / nested-refit inner
#'   validation fraction. `NULL` uses `0.15`.
#' @param group_col Character or `NULL`. Grouping column for the grouped
#'   train/val split. `NULL` uses `"block_id"`.
#' @param impute_numeric,impute_factor_missing Character or `NULL`. Numeric
#'   imputation rule (`"median"` / `"zero"`) and the missing-factor sentinel.
#'   `NULL` uses `"median"` / `"MISSING"`.
#' @param cap_contextual,cap_spectral,cap_random,cap_otsu Numeric `>= 0`
#'   (`Inf` disables) or `NULL`. The four negative-bucket caps as a multiple of
#'   `n_burned`. `NULL` uses `0.25` / `1.0` / `1.0` / `1.0`. Stored in
#'   `cfg$train_control$caps`.
#' @param feature_whitelist_override,feature_weights Optional subset of
#'   `.supervised_feature_cols` and a named numeric per-feature weight vector,
#'   or `NULL`. Stored in `cfg$train_control` and consumed by BOTH OOF and FINAL.
#' @param training_protocol,oof_sampling Character or `NULL`. Training protocol
#'   (`"legacy"` / `"nested_refit"`) and OOF per-fold negative sampling mode
#'   (`"capped"` / `"full"`). `NULL` uses `"legacy"` / `"capped"`.
#' @param model_params Named list or `NULL`. Optional override of the canonical
#'   XGBoost hyperparameter block stored in `cfg$model_params` (the canonical
#'   block MINUS `scale_pos_weight`, which is computed at train time). Supplying
#'   it merges your fields onto the canonical block; `scale_pos_weight` is
#'   rejected here.
#'
#' @section Gate 1B (2026-06-07) cfg single source of truth:
#' `cfg$model_params` (the methodological XGBoost block, without
#' `scale_pos_weight`) and `cfg$train_control` (caps, seeds, nrounds/early-stop,
#' val_frac, group_col, imputation rules, whitelist/weights, protocol toggles)
#' are the SINGLE SOURCE OF TRUTH for the supervised methodological parameters.
#' The canonical defaults live ONLY in `.of_canonical_model_params()` /
#' `.of_canonical_train_control()`; the builder arguments above let a user
#' override them. Every downstream stage ([run_oneyear_supervised_pipeline()],
#' [run_oof_diagnostics()], [train_final_burned_model()] and the internal
#' engines) reads these resolved cfg sections rather than carrying its own
#' defaults.
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
    # Gate 1B PIECE 2 (2026-06-07): the ecoregion border used by the legacy Otsu
    # builder's ecoregion / corine_ecoregion modes is now a configurable cfg
    # input (was a data_base/Ecoregion hardcode). NULL -> convention default.
    ecoregion_shapefile = NULL,
    # ---- Gate 1B (2026-06-07): resolved methodological / training params ----
    # The cfg is the SINGLE SOURCE OF TRUTH for these. Every argument defaults
    # to NULL = "use the canonical default" (the canonical defaults live ONLY
    # in .of_canonical_train_control() / .of_canonical_model_params(), nowhere
    # else). A non-NULL value OVERRIDES the canonical default and is validated +
    # stored in cfg$train_control / cfg$model_params, which is what every
    # downstream stage then reads.
    nrounds_max                = NULL,
    early_stop                 = NULL,
    oof_seed_base              = NULL,
    final_sampling_seed        = NULL,
    final_seed                 = NULL,
    val_frac                   = NULL,
    group_col                  = NULL,
    impute_numeric             = NULL,
    impute_factor_missing      = NULL,
    cap_contextual             = NULL,
    cap_spectral               = NULL,
    cap_random                 = NULL,
    cap_otsu                   = NULL,
    feature_whitelist_override = NULL,
    feature_weights            = NULL,
    training_protocol          = NULL,
    oof_sampling               = NULL,
    model_params               = NULL,
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
  ecoregion_shapefile  <- ecoregion_shapefile  %||% conv$ecoregion_shapefile

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
                                                     allow_null = TRUE),
    # Gate 1B PIECE 2: ecoregion border (legacy Otsu ecoregion modes).
    ecoregion_shapefile  = .of_normalize_input_spec(ecoregion_shapefile,
                                                     "ecoregion_shapefile",
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
  # Gate 1B PIECE 2: `ecoregion_shapefile` is deliberately NOT in this
  # config-time existence loop. It is consumed ONLY by the non-default legacy
  # Otsu `ecoregion` / `corine_ecoregion` modes (the default is
  # `burnable_only`), and its convention default (data_base/Ecoregion) is
  # absent in most setups; the legacy builder asserts its existence inside the
  # ecoregion branch, exactly when (and only when) it is actually needed.
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

  # ---- Gate 1B (2026-06-07): resolve cfg$model_params + cfg$train_control ----
  # The cfg becomes the SINGLE SOURCE OF TRUTH for all supervised
  # methodological / training-control parameters. Canonical defaults come from
  # .of_canonical_model_params() / .of_canonical_train_control() (the ONLY
  # places they are defined); user-supplied builder args override them and are
  # validated here. Every downstream stage reads cfg$model_params /
  # cfg$train_control, so the 5-layer default duplication is eliminated.
  # Capture whether model_params was explicitly supplied BEFORE it is folded
  # onto the canonical block (provenance "user" vs "default").
  .model_params_user_set <- !is.null(model_params)
  model_params <- .of_resolve_supervised_model_params(model_params)
  train_control <- .of_resolve_supervised_train_control(
    nrounds_max                = nrounds_max,
    early_stop                 = early_stop,
    oof_seed_base              = oof_seed_base,
    final_sampling_seed        = final_sampling_seed,
    final_seed                 = final_seed,
    val_frac                   = val_frac,
    group_col                  = group_col,
    impute_numeric             = impute_numeric,
    impute_factor_missing      = impute_factor_missing,
    cap_contextual             = cap_contextual,
    cap_spectral               = cap_spectral,
    cap_random                 = cap_random,
    cap_otsu                   = cap_otsu,
    feature_whitelist_override = feature_whitelist_override,
    feature_weights            = feature_weights,
    training_protocol          = training_protocol,
    oof_sampling               = oof_sampling
  )

  # ---- Gate 1B / Precision 1 (2026-06-07): provenance of each methodological
  # train_control field. Tagged "user" when the caller explicitly supplied the
  # corresponding builder argument (non-NULL), "default" when the canonical
  # default was used. The public-boundary shim resolver
  # (.of_resolve_methodological_shim) consults this to (a) ERROR when a
  # function-level override conflicts with an EXPLICIT user value, and (b) emit a
  # deprecation warning when a function-level override is applied over a
  # canonical default. The provenance record is exposed on the cfg AND written
  # into the run-level summary the pipeline emits.
  tc_provenance <- list(
    nrounds_max           = if (is.null(nrounds_max))           "default" else "user",
    early_stop            = if (is.null(early_stop))            "default" else "user",
    oof_seed_base         = if (is.null(oof_seed_base))         "default" else "user",
    final_sampling_seed   = if (is.null(final_sampling_seed))   "default" else "user",
    final_seed            = if (is.null(final_seed))            "default" else "user",
    val_frac              = if (is.null(val_frac))              "default" else "user",
    group_col             = if (is.null(group_col))             "default" else "user",
    impute_numeric        = if (is.null(impute_numeric))        "default" else "user",
    impute_factor_missing = if (is.null(impute_factor_missing)) "default" else "user",
    cap_contextual        = if (is.null(cap_contextual))        "default" else "user",
    cap_spectral          = if (is.null(cap_spectral))          "default" else "user",
    cap_random            = if (is.null(cap_random))            "default" else "user",
    cap_otsu              = if (is.null(cap_otsu))              "default" else "user",
    feature_whitelist_override =
      if (is.null(feature_whitelist_override)) "default" else "user",
    feature_weights       = if (is.null(feature_weights))       "default" else "user",
    training_protocol     = if (is.null(training_protocol))     "default" else "user",
    oof_sampling          = if (is.null(oof_sampling))          "default" else "user",
    model_params          = if (.model_params_user_set)         "user"    else "default"
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
    # Gate 1B (2026-06-07): resolved methodological / training params. THE
    # single source of truth read by every supervised stage.
    model_params             = model_params,
    train_control            = train_control,
    # Gate 1B / Precision 1: per-field provenance ("default" canonical vs "user"
    # explicitly set in the builder). The public-boundary shim resolver enriches
    # this into a full requested/resolved record at run time; the static builder
    # record below is the authoritative "was this field user-set?" source the
    # conflict-error logic relies on.
    resolved_params_provenance = list(train_control = tc_provenance),
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
  if (!is.null(x$train_control)) {
    tc <- x$train_control
    cat("  train_control   : protocol=", tc$training_protocol,
        " nrounds=", tc$nrounds_max, " early_stop=", tc$early_stop,
        " val_frac=", tc$val_frac, "\n", sep = "")
    cat("    caps          : contextual=", tc$caps$contextual,
        " spectral=", tc$caps$spectral, " random=", tc$caps$random,
        " otsu=", tc$caps$otsu, "\n", sep = "")
  }
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

#' Resolve + validate cfg$model_params (Gate 1B).
#'
#' Canonical default = .of_canonical_model_params() (the canonical xgb block
#' minus scale_pos_weight). A user override (named list) is merged on top of
#' the canonical block and validated. scale_pos_weight is rejected here: it is
#' site-specific (computed at train time), never a cfg field.
#'
#' @keywords internal
#' @noRd
.of_resolve_supervised_model_params <- function(model_params) {
  base <- .of_canonical_model_params()
  if (is.null(model_params)) return(base)
  if (!is.list(model_params) || is.null(names(model_params)) ||
      any(!nzchar(names(model_params)))) {
    stop("'model_params' must be NULL or a named list of XGBoost params.",
         call. = FALSE)
  }
  if ("scale_pos_weight" %in% names(model_params)) {
    stop("'model_params' must not set 'scale_pos_weight': it is computed at ",
         "train time from the training labels, not a cfg field.", call. = FALSE)
  }
  for (nm in names(model_params)) base[[nm]] <- model_params[[nm]]
  base
}

#' Resolve + validate cfg$train_control (Gate 1B).
#'
#' Each argument is NULL (= canonical default from .of_canonical_train_control())
#' or a user override. Validates numeric ranges, match.arg-style enums, and
#' scalar types, then returns the resolved list stored as cfg$train_control.
#'
#' @keywords internal
#' @noRd
.of_resolve_supervised_train_control <- function(
    nrounds_max, early_stop, oof_seed_base, final_sampling_seed, final_seed,
    val_frac, group_col, impute_numeric, impute_factor_missing,
    cap_contextual, cap_spectral, cap_random, cap_otsu,
    feature_whitelist_override, feature_weights,
    training_protocol, oof_sampling) {

  tc <- .of_canonical_train_control()

  # --- positive-integer / positive-numeric controls --------------------------
  chk_pos_int <- function(v, nm) {
    vi <- suppressWarnings(as.integer(v)[1L])
    if (is.na(vi) || vi <= 0L) {
      stop(sprintf("'%s' must be a single positive integer.", nm),
           call. = FALSE)
    }
    vi
  }
  if (!is.null(nrounds_max)) tc$nrounds_max <- chk_pos_int(nrounds_max, "nrounds_max")
  if (!is.null(early_stop))  tc$early_stop  <- chk_pos_int(early_stop, "early_stop")

  if (!is.null(oof_seed_base))
    tc$seeds$oof_seed_base <- chk_pos_int(oof_seed_base, "oof_seed_base")
  if (!is.null(final_sampling_seed))
    tc$seeds$final_sampling_seed <- chk_pos_int(final_sampling_seed,
                                                "final_sampling_seed")
  if (!is.null(final_seed))
    tc$seeds$final_seed <- chk_pos_int(final_seed, "final_seed")

  # --- val_frac in (0, 1) ----------------------------------------------------
  if (!is.null(val_frac)) {
    if (!is.numeric(val_frac) || length(val_frac) != 1L || is.na(val_frac) ||
        val_frac <= 0 || val_frac >= 1) {
      stop("'val_frac' must be a single number in (0, 1).", call. = FALSE)
    }
    tc$val_frac <- as.numeric(val_frac)
  }

  # --- group_col: single non-empty string or NULL ----------------------------
  if (!is.null(group_col)) {
    if (!is.character(group_col) || length(group_col) != 1L ||
        !nzchar(group_col)) {
      stop("'group_col' must be a single non-empty character string or NULL.",
           call. = FALSE)
    }
    tc$group_col <- group_col
  }

  # --- imputation rules ------------------------------------------------------
  if (!is.null(impute_numeric)) {
    tc$impute_numeric <- match.arg(impute_numeric, c("median", "zero"))
  }
  if (!is.null(impute_factor_missing)) {
    if (!is.character(impute_factor_missing) ||
        length(impute_factor_missing) != 1L ||
        is.na(impute_factor_missing) || !nzchar(impute_factor_missing)) {
      stop("'impute_factor_missing' must be a single non-empty character ",
           "string.", call. = FALSE)
    }
    tc$impute_factor_missing <- impute_factor_missing
  }

  # --- negative-bucket caps: numeric >= 0 (Inf allowed = disable) -------------
  chk_cap <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop(sprintf("'%s' must be a single number >= 0 (Inf disables the cap).",
                   nm), call. = FALSE)
    }
    as.numeric(v)
  }
  if (!is.null(cap_contextual)) tc$caps$contextual <- chk_cap(cap_contextual, "cap_contextual")
  if (!is.null(cap_spectral))   tc$caps$spectral   <- chk_cap(cap_spectral, "cap_spectral")
  if (!is.null(cap_random))     tc$caps$random     <- chk_cap(cap_random, "cap_random")
  if (!is.null(cap_otsu))       tc$caps$otsu       <- chk_cap(cap_otsu, "cap_otsu")

  # --- feature whitelist / weights (pass-through, light validation) ----------
  if (!is.null(feature_whitelist_override)) {
    if (!is.character(feature_whitelist_override) ||
        length(feature_whitelist_override) < 1L) {
      stop("'feature_whitelist_override' must be NULL or a non-empty ",
           "character vector.", call. = FALSE)
    }
    tc$feature_whitelist_override <- feature_whitelist_override
  }
  if (!is.null(feature_weights)) {
    if (!is.numeric(feature_weights) || is.null(names(feature_weights)) ||
        any(!nzchar(names(feature_weights)))) {
      stop("'feature_weights' must be NULL or a NAMED numeric vector.",
           call. = FALSE)
    }
    tc$feature_weights <- feature_weights
  }

  # --- protocol toggles ------------------------------------------------------
  if (!is.null(training_protocol)) {
    tc$training_protocol <- match.arg(training_protocol,
                                      c("legacy", "nested_refit"))
  }
  if (!is.null(oof_sampling)) {
    tc$oof_sampling <- match.arg(oof_sampling, c("capped", "full"))
  }

  tc
}

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
    # Gate 1B PIECE 2: ecoregion border convention (legacy Otsu ecoregion modes).
    ecoregion_shapefile = if (has_db) {
      file.path(data_base, "Ecoregion", "ecoregiones_olson.shp")
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

#' Resolve a cfg$inputs on-disk path (Gate 1B input inventory).
#'
#' Single shared accessor for the supervised-input subsystem: returns the
#' normalized on-disk PATH of `config$inputs[[name]]` when that input is a
#' `type = "path"` spec, else `NULL` (the input is absent, in-memory, or the
#' config itself is NULL). Every supervised stage that needs an input's file
#' path (orchestrator, pool builder, feature extractor) consumes the input
#' THROUGH this accessor, so `cfg$inputs` is the single source of truth and no
#' stage reconstructs an input path by filename convention behind the cfg's
#' back. A `NULL` return lets the caller apply its documented optional-input
#' behaviour (skip / convention-fallback for the standalone-script path).
#'
#' @param config supervised config S3 object, or `NULL`.
#' @param name character scalar input name (a key of `cfg$inputs`).
#' @return character path or `NULL`.
#'
#' @keywords internal
#' @noRd
.of_sup_input_path <- function(config, name) {
  if (is.null(config) || is.null(config$inputs)) return(NULL)
  sp <- config$inputs[[name]]
  if (is.null(sp)) return(NULL)
  if (!is.null(sp$path) && length(sp$path) == 1L && !is.na(sp$path) &&
      nzchar(sp$path) && identical(sp$type, "path")) {
    sp$path
  } else {
    NULL
  }
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
