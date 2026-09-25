# Internal supervised dispatcher. NOT exported.
#
# Block 8: previously the dispatcher injected per-call configuration by
# unlocking + reassigning top-level placeholder bindings inside the
# `OtsuFire` namespace, then restoring the previous values via on.exit.
# That works but trips R CMD check ("possibly unsafe calls:
# unlockBinding(nm, ns)") and is poor citizenship for a CRAN-bound
# package. An interim refactor replaced it with a per-call execution
# environment whose parent was the package namespace, switching the
# orchestrator function's environment so its lexical lookups (data_base,
# composite_base, UNB_OTSU_NEG_*, REGISTRY_*, ...) resolved there first.
#
# BUG 3 Phase 1b (2026-06-05): that env-injection bridge is now removed too.
# `.of_supervised_engine_bindings(config)` still builds the named list of
# config -> orchestrator-global values, but the dispatcher now passes it
# EXPLICITLY to `run_supervised_pipeline(engine_bindings = ...)`, which
# list2env()s it into its own execution frame as its first statement. No
# `environment(pipeline) <- env` magic, no namespace mutation, and no
# per-call env construction (`.of_make_supervised_engine_env()` deleted).

# BUG 3 Phase 1a (2026-06-05): removed `.of_require_supervised_engine()`.
# It had no callers and only existed to locate/source an external supervised
# engine file (returning `pipeline_path = NULL`). The supervised engine is now
# in-package, so the function was dead.

# Build the list of names + values that the orchestrator expects to find
# at top level (data_base, UNB_*, REGISTRY_*, ...). Pure data-only.
#' @keywords internal
#' @noRd
.of_supervised_engine_bindings <- function(config) {
  # GATE 6.7 (2026-06-12): the negative-pool methodological knobs, caps, random
  # seed and runtime toggles are sourced from the TYPED cfg blocks
  # (negative_pool_params / negative_pool_runtime / negative_pool_internal +
  # train_control$seeds$random_seed), NOT from free-form config$options$unb_* /
  # otsu_negative_*. There is ONE source of truth; the old free-form options are
  # no longer read.
  np  <- config$negative_pool_params   %||% .of_negative_pool_params_defaults()
  rt  <- config$negative_pool_runtime  %||% .of_runtime_options_defaults()
  ni  <- config$negative_pool_internal %||%
           .of_supervised_negative_pool_internal_defaults()
  sds <- (config$train_control %||% list())$seeds %||% list()
  list(
    data_base                = config$options$data_base,
    # BUG 3 Phase 1b (2026-06-05): dropped the `script_base = NULL` binding (and
    # its orchestrator placeholder). `script_base` was dead plumbing into
    # verify_otsu_negative_helpers(), which ignores it; the arg has been
    # removed end-to-end.
    composite_base           = config$options$composite_base,
    result_name              = config$options$result_name %||% "Min_Min",
    # 2026-06-05: the SUPERVISED output base is the SINGLE SOURCE OF TRUTH
    # for WHERE outputs are written, driven by output_dir + run_name via
    # `.of_build_supervised_output_routes()`. The orchestrator consumes this
    # instead of reconstructing the path from data_base + result_name (which
    # only resolve INPUT imagery/composite/decision paths). For Natalia's
    # default config (output_dir=<data_base>/Results, run_name=result_name)
    # this resolves to the identical path, so outputs do not move.
    supervised_output_base   = config$output_routes$base,
    # BUG 3 Phase 1b (2026-06-05): dropped the `supervised_functions_dir = NULL`
    # binding (and its orchestrator placeholder). It was inert metadata the
    # orchestrator never read.
    years_to_run             = c(config$target_year),
    scenarios_to_run         = c(config$scenario),
    DO_POOLS                 = TRUE,
    DO_FOLDS                 = TRUE,
    DO_FEATURES              = TRUE,
    DO_MODEL                 = TRUE,
    # §N+26 (2026-06-05): the `UNB_STRATEGY = config$negative_pool_policy`
    # binding was removed. The negative-pool policy is always `all_sources`;
    # the orchestrator no longer carries a `UNB_STRATEGY` binding at all (its
    # only former reader, the common_unb_dir selection, is hardwired to the
    # all_sources path), so nothing needs to be threaded here.
    # GATE 6.7 (2026-06-12): runtime verbose comes from the typed runtime block.
    UNB_VERBOSE              = rt$verbose,
    # GATE 6.7: random-background negative params come from the typed PUBLIC block
    # (n_cells / rbr_quantile), the internal-defaults block (exclusion buffer,
    # patch size) and the canonical seeds block (random_seed). The free-form
    # config$options$unb_* reads are GONE (single source of truth).
    UNB_EXCL_BUFFER_M        = ni$random_exclusion_buffer_m,
    UNB_N_RANDOM_CELLS       = np$random$n_cells,
    UNB_RANDOM_RBR_Q         = np$random$rbr_quantile,
    UNB_RANDOM_SEED          = sds$random_seed %||% 42L,
    UNB_RANDOM_PATCH_SIZE_CELLS = ni$random_patch_size_cells,
    # B2 (2026-06-05): supervised output-naming prefixes. Defaults reproduce
    # the historical "patch" / "patch_certified" run prefixes exactly.
    prefix_oof_base          = config$options$prefix_oof_base %||% "patch",
    prefix_base              = config$options$prefix_base %||% "patch_certified",
    # GATE 6.7: Otsu residual negative params. The two PUBLIC knobs
    # (candidate_threshold / reference_threshold) come from negative_pool_params;
    # the rest are FIXED internal defaults (negative_pool_internal).
    UNB_OTSU_NEG_MODE        = ni$otsu_mode,
    UNB_OTSU_NEG_THRESHOLD   = np$otsu$candidate_threshold,
    UNB_OTSU_NEG_REFERENCE_THRESHOLD = np$otsu$reference_threshold,
    UNB_OTSU_NEG_MIN_THRESHOLD_VALUE = ni$otsu_min_threshold_value,
    UNB_OTSU_NEG_MIN_PIXELS  = ni$otsu_min_pixels,
    UNB_OTSU_NEG_BUFFERS_M   = ni$otsu_buffers_m,
    UNB_OTSU_NEG_CORE_THR    = ni$otsu_core_thr,
    UNB_OTSU_NEG_ALPHA_BOOST = ni$otsu_alpha_boost,
    UNB_OTSU_NEG_MIN_BASE_BOOST = ni$otsu_min_base_boost,
    UNB_OTSU_NEG_DIST_POWER  = ni$otsu_dist_power,
    UNB_OTSU_NEG_KEEP_HI     = ni$otsu_keep_hi,
    UNB_OTSU_NEG_DROP_LO     = ni$otsu_drop_lo,
    UNB_OTSU_NEG_EXCL_BUFFER_M = ni$otsu_excl_buffer_m,
    UNB_OTSU_NEG_MIN_AREA_HA = ni$otsu_min_area_ha,
    # GATE 6.7: reuse_existing / write_outputs are TECHNICAL runtime toggles.
    UNB_OTSU_NEG_REUSE_EXISTING = rt$reuse_existing,
    UNB_OTSU_NEG_WRITE_OUTPUT   = rt$write_outputs,
    UNB_OTSU_NEG_USE_DROP    = ni$otsu_use_drop,
    UNB_OTSU_NEG_DROP_MAX_S_PATCH   = ni$otsu_drop_max_s_patch,
    # AS01 (0.3.0): plumb external-tool paths from config$tool_paths.
    # Default NULL — the Otsu negative pipeline stops with a clear message if a
    # path it needs is missing, instead of falling through to a hard-coded
    # Olga.Viedma default.
    python_exe               = config$tool_paths$python_exe,
    gdal_polygonize_script   = config$tool_paths$gdal_polygonize_script,
    gdalwarp_path            = config$tool_paths$gdalwarp_path,
    ogr2ogr_exe              = config$tool_paths$ogr2ogr_exe,
    CURRENTYEAR_PREYEAR_OVERLAP_THR =
      config$options$currentyear_preyear_overlap_thr %||% 0.70,
    CURRENTYEAR_HOTSPOT_DENSITY_THR =
      config$options$currentyear_hotspot_density_thr %||% 0.001,
    CURRENTYEAR_TEMPORAL_PENALTY_FLOOR =
      config$options$currentyear_temporal_penalty_floor %||% 0.10
    # §N+27 (2026-06-05): RBR_TRAINING_REGISTRY_PATH / REGISTRY_PROB_THR /
    # REGISTRY_TOP_FRACTION bindings removed (abandoned supervised registry).
    # B7 (2026-06-06): the `MIN_BURNED_POOL_N = config$min_burned_pool_n`
    # binding was removed. It was injected into the orchestrator frame but the
    # live sparse-year guard reads `config$min_burned_pool_n` directly in
    # build_supervised_training_pools() (supervised-pools.R); the uppercase
    # binding + its orchestrator placeholder were never read.
  )
}

#' Internal adapter: runs the full validated probabilistic refinement chain for one year.
#'
#' OtsuFire 0.5.0 (2026-05-09): the dispatcher forwards the supervised
#' final-model sampling caps, the new `feature_whitelist_override` and
#' `feature_weights` hooks (replacing the removed `additional_drop_cols`
#' deny-list), and the `reuse_upstream` execution toggle. Defaults
#' reproduce historical behaviour byte-for-byte.
#'
#' @keywords internal
#' @noRd
.of_run_supervised_oneyear <- function(config, run_consistency = TRUE,
                                       overwrite = FALSE,
                                       # Gate 1B (2026-06-07): these methodological
                                       # params are REQUIRED resolved args with NO
                                       # defaults. The single source of truth is
                                       # cfg$train_control / cfg$model_params; the
                                       # public run_oneyear_supervised_pipeline()
                                       # resolves them (cfg + explicit overrides)
                                       # and passes them here. A dropped arg ERRORS
                                       # in the guard block below.
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
                                       final_impute_factor_missing) {
  # Gate 1B: required-arg guard (no silent methodological defaults).
  .req <- c("random_to_burned_ratio", "otsu_unburned_to_burned_ratio",
            "feature_whitelist_override", "feature_weights",
            "oof_nrounds_max", "oof_early_stop", "oof_seed_base",
            "final_sampling_seed", "final_seed", "final_val_frac",
            "final_group_col", "final_nrounds_max", "final_early_stopping_rounds",
            "final_impute_numeric", "final_impute_factor_missing")
  for (.nm in .req) {
    if (eval(call("missing", as.name(.nm)))) {
      stop(".of_run_supervised_oneyear(): required resolved arg '", .nm,
           "' is missing (no methodological default; resolved from ",
           "cfg$train_control / cfg$model_params by the public entry).",
           call. = FALSE)
    }
  }
  ns <- asNamespace("OtsuFire")

  # Resolve the internal_decisions path ONCE, here at the top of the chain.
  # `build_supervised_burned_config()` requires `internal_decisions`, and
  # `.of_normalize_input_spec()` stores a character path as a
  # `type = "path"` spec whose `$path` is the normalized filesystem path.
  # That user-supplied path is THE single source of truth for the
  # deterministic decisions GPKG at every downstream read site
  # (orchestrator, both unburned builders, consistency check). There is NO
  # fallback to the historical Min_Min convention path: if the user did not
  # supply a usable on-disk GPKG path we fail fast here with a clear
  # message, rather than silently reconstructing a path by convention.
  id_spec <- config$inputs$internal_decisions
  eff_internal_decisions <- if (!is.null(id_spec) &&
                                 is.character(id_spec$path) &&
                                 length(id_spec$path) == 1L &&
                                 !is.na(id_spec$path) &&
                                 nzchar(id_spec$path)) {
    id_spec$path
  } else {
    NULL
  }
  if (is.null(eff_internal_decisions)) {
    stop(
      "Supervised one-year pipeline requires an on-disk 'internal_decisions' ",
      "GPKG path. The config either has no internal_decisions or it was ",
      "supplied as an in-memory object (sf/SpatVector). Write the ",
      "deterministic decisions to a .gpkg (layer 'internal_decisions') and ",
      "pass that file path to build_supervised_burned_config(internal_decisions=...). ",
      "There is no convention-based fallback.",
      call. = FALSE
    )
  }

  # BUG 3 Phase 1b (2026-06-05): explicit config passing replaces the
  # env-injection bridge. We build the SAME named list the orchestrator's
  # top-level lexical lookups expect (data_base, composite_base, UNB_*,
  # REGISTRY_*, tool paths, supervised_output_base, prefix_*, ...) and hand
  # it to `run_supervised_pipeline()` via the new `engine_bindings` param,
  # which list2env()s it into the pipeline's own execution frame as the very
  # first statement. The previous `eng_env <- new.env(parent = ns)` +
  # `environment(pipeline) <- eng_env` magic (and the
  # `.of_make_supervised_engine_env()` helper) are gone. No namespace
  # mutation; the orchestrator function value is used unchanged.
  bindings <- .of_supervised_engine_bindings(config)
  pipeline <- get("run_supervised_pipeline", envir = ns, inherits = FALSE)

  globals_budget <- config$options$future_globals_maxsize %||% (8 * 1024^3)
  old_globals <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = globals_budget)
  on.exit(options(future.globals.maxSize = old_globals), add = TRUE)

  out <- pipeline(
    target_year       = config$target_year,
    scenario          = config$scenario,
    min_burned_pool_n = config$min_burned_pool_n,
    # 2026-06-05: thread the user's overwrite flag into the engine
    # (previously dropped here; the engine hardcoded overwrite=TRUE).
    overwrite         = overwrite,
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    feature_whitelist_override             = feature_whitelist_override,
    feature_weights                        = feature_weights,
    reuse_upstream                         = reuse_upstream,
    # 2026-06-05 (D1 expose): forward OOF + FINAL training knobs.
    oof_nrounds_max                        = oof_nrounds_max,
    oof_early_stop                         = oof_early_stop,
    oof_seed_base                          = oof_seed_base,
    final_sampling_seed                    = final_sampling_seed,
    final_seed                             = final_seed,
    final_val_frac                         = final_val_frac,
    final_group_col                        = final_group_col,
    final_nrounds_max                      = final_nrounds_max,
    final_early_stopping_rounds            = final_early_stopping_rounds,
    final_impute_numeric                   = final_impute_numeric,
    final_impute_factor_missing            = final_impute_factor_missing,
    internal_decisions_path                = eff_internal_decisions,
    engine_bindings                        = bindings,
    # BUG 3 Phase 2 PILOT (2026-06-05): thread the supervised config S3 object
    # into the orchestrator so its modular stages (make_spatial_folds(), ...)
    # can resolve output_routes/target_year directly. The engine_bindings list
    # carries only unpacked scalars, not the config object itself.
    config                                 = config
  )

  # Block 6: consistency check is package-internal and runs entirely
  # from the config; no namespace mutation needed.
  consistency_out <- NULL
  if (isTRUE(run_consistency)) {
    data_base   <- config$options$data_base
    if (is.null(data_base) || !nzchar(data_base) || !dir.exists(data_base)) {
      warning("run_consistency=TRUE requires options$data_base. ",
              "Skipping consistency check.", call. = FALSE)
    } else {
      # The user-supplied internal_decisions path is the single source of
      # truth. No convention reconstruction.
      det_gpkg <- eff_internal_decisions
      # 2026-06-05: the consistency check reads the FINAL MAP, a supervised
      # OUTPUT. Its location is the config's output_routes (output_dir +
      # run_name), NOT a data_base/result_name reconstruction, so it always
      # matches where the engine actually wrote.
      fmap_gpkg <- config$output_routes$final_map_gpkg
      cdir      <- config$output_routes$consistency_dir
      consistency_out <- tryCatch(
        get(".of_run_consistency_check", envir = ns, inherits = FALSE)(
          target_year = config$target_year,
          scenario    = config$scenario,
          internal_decisions_gpkg = det_gpkg,
          final_map_gpkg = fmap_gpkg,
          out_dir = cdir
        ),
        error = function(e) {
          warning("Consistency check failed: ", conditionMessage(e),
                  call. = FALSE)
          NULL
        }
      )
    }
  }

  list(legacy_run_summary = out, consistency = consistency_out)
}
