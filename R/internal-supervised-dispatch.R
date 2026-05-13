# Internal supervised dispatcher. NOT exported.
#
# Block 8: previously the dispatcher injected per-call configuration by
# unlocking + reassigning top-level placeholder bindings inside the
# `OtsuFire` namespace, then restoring the previous values via on.exit.
# That works but trips R CMD check ("possibly unsafe calls:
# unlockBinding(nm, ns)") and is poor citizenship for a CRAN-bound
# package. The refactor replaces it with a per-call execution
# environment whose parent is the package namespace. We then switch
# the orchestrator function's environment to that env, so its lexical
# lookups (data_base, composite_base, UNB_LEGACY_*, REGISTRY_*, ...)
# resolve there first and fall through to the namespace for everything
# else. No namespace bindings are unlocked or modified.

#' @keywords internal
#' @noRd
.of_require_supervised_engine <- function(config) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("Internal error: dispatcher expects otsufire_supervised_burned_config.",
         call. = FALSE)
  }
  list(scripts_root = config$scripts_root, pipeline_path = NULL)
}

# Build the list of names + values that the orchestrator expects to find
# at top level (data_base, UNB_*, REGISTRY_*, ...). Pure data-only.
#' @keywords internal
#' @noRd
.of_supervised_engine_bindings <- function(config) {
  list(
    data_base                = config$options$data_base,
    script_base              = config$scripts_root,
    composite_base           = config$options$composite_base,
    result_name              = config$options$result_name %||% "Min_Min",
    supervised_functions_dir = config$supervised_engine_root,
    years_to_run             = c(config$target_year),
    scenarios_to_run         = c(config$scenario),
    DO_POOLS                 = TRUE,
    DO_FOLDS                 = TRUE,
    DO_FEATURES              = TRUE,
    DO_MODEL                 = TRUE,
    UNB_STRATEGY             = config$negative_pool_policy,
    UNB_VERBOSE              = config$options$unb_verbose %||% TRUE,
    UNB_LEGACY_CODE_DIR      = config$options$legacy_code_dir,
    UNB_LEGACY_OTSU_MODE     = config$options$legacy_otsu_mode %||% "burnable_only",
    UNB_LEGACY_OTSU_THRESHOLD = config$options$legacy_otsu_threshold %||% 0,
    UNB_LEGACY_REFERENCE_OTSU_THRESHOLD =
      config$options$legacy_reference_otsu_threshold %||% 100,
    UNB_LEGACY_SAMPLE_N      = config$options$legacy_sample_n %||% 2000,
    UNB_LEGACY_SAMPLE_PROPS  = config$options$legacy_sample_props %||%
      c(drop = 0.70, review = 0.25, keep = 0.05),
    UNB_LEGACY_REUSE_EXISTING = config$options$legacy_reuse_existing %||% TRUE,
    UNB_LEGACY_WRITE_OUTPUT   = config$options$legacy_write_output %||% TRUE,
    UNB_LEGACY_USE_DROP      = config$options$legacy_use_drop   %||% TRUE,
    UNB_LEGACY_USE_REVIEW    = config$options$legacy_use_review %||% FALSE,
    UNB_LEGACY_USE_KEEP      = config$options$legacy_use_keep   %||% FALSE,
    # AS03 (0.3.0): expose RANDOM_SEED and the three S_PATCH caps. Other
    # UNB_LEGACY_* parameters remain pinned at orchestrator defaults.
    UNB_LEGACY_RANDOM_SEED   = config$options$legacy_random_seed %||% 42L,
    UNB_LEGACY_DROP_MAX_S_PATCH   = config$options$legacy_drop_max_s_patch   %||% 0.15,
    UNB_LEGACY_REVIEW_MAX_S_PATCH = config$options$legacy_review_max_s_patch %||% 0.45,
    UNB_LEGACY_KEEP_MAX_S_PATCH   = config$options$legacy_keep_max_s_patch   %||% 0.70,
    # AS01 (0.3.0): plumb external-tool paths from config$tool_paths.
    # Default NULL — the legacy pipeline stops with a clear message if a
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
      config$options$currentyear_temporal_penalty_floor %||% 0.10,
    RBR_TRAINING_REGISTRY_PATH = config$burned_like_registry_path,
    REGISTRY_PROB_THR        = config$options$registry_prob_threshold %||% 0.90,
    REGISTRY_TOP_FRACTION    = config$options$registry_top_fraction %||% 0.01,
    MIN_BURNED_POOL_N        = config$min_burned_pool_n
  )
}

# Build a per-call execution env that overrides the package namespace
# for the orchestrator's top-level lexical lookups, without touching
# any namespace bindings.
#' @keywords internal
#' @noRd
.of_make_supervised_engine_env <- function(config) {
  ns <- asNamespace("OtsuFire")
  env <- new.env(parent = ns)
  bindings <- .of_supervised_engine_bindings(config)
  for (nm in names(bindings)) {
    # NULL is a legitimate value for some of these names (e.g.
    # UNB_LEGACY_CODE_DIR when no legacy_code_dir was configured);
    # assigning NULL into the local env is correct and overrides the
    # namespace-level placeholder while we run.
    assign(nm, bindings[[nm]], envir = env)
  }
  env
}

#' Internal adapter: runs the full validated supervised chain for one year.
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
                                       contextual_exclusion_to_burned_ratio   = 0.25,
                                       spectral_hard_negative_to_burned_ratio = 1.0,
                                       random_to_burned_ratio                 = 1.0,
                                       otsu_unburned_to_burned_ratio          = 1.0,
                                       feature_whitelist_override             = NULL,
                                       feature_weights                        = NULL,
                                       reuse_upstream                         = FALSE) {
  ns <- asNamespace("OtsuFire")

  # Per-call env carrying all configuration the orchestrator looks up
  # via lexical scoping. Parent is `ns`, so any name not defined here
  # still resolves through the package namespace.
  eng_env <- .of_make_supervised_engine_env(config)

  # Re-environment the orchestrator function so its body sees eng_env
  # first. We copy the function value out of the namespace; the
  # original `ns$run_supervised_pipeline` is unchanged. No namespace
  # mutation.
  pipeline <- get("run_supervised_pipeline", envir = ns, inherits = FALSE)
  environment(pipeline) <- eng_env

  globals_budget <- config$options$future_globals_maxsize %||% (8 * 1024^3)
  old_globals <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = globals_budget)
  on.exit(options(future.globals.maxSize = old_globals), add = TRUE)

  out <- pipeline(
    target_year       = config$target_year,
    scenario          = config$scenario,
    min_burned_pool_n = config$min_burned_pool_n,
    contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
    spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    feature_whitelist_override             = feature_whitelist_override,
    feature_weights                        = feature_weights,
    reuse_upstream                         = reuse_upstream
  )

  # Block 6: consistency check is package-internal and runs entirely
  # from the config; no namespace mutation needed.
  consistency_out <- NULL
  if (isTRUE(run_consistency)) {
    data_base   <- config$options$data_base
    result_name <- config$options$result_name %||% "Min_Min"
    if (is.null(data_base) || !nzchar(data_base) || !dir.exists(data_base)) {
      warning("run_consistency=TRUE requires options$data_base. ",
              "Skipping consistency check.", call. = FALSE)
    } else {
      det_gpkg <- file.path(data_base, "Results", config$target_year,
                            result_name, "DETERMINISTIC", config$scenario,
                            "05_DECISIONS", "internal_decisions.gpkg")
      prefix   <- sprintf("%d_%s_patch_certified",
                          config$target_year, config$scenario)
      fmap_gpkg <- file.path(data_base, "Results", config$target_year,
                             result_name, "SUPERVISED", config$scenario,
                             "09_FINAL_MAP",
                             paste0(prefix, "_final_map.gpkg"))
      cdir <- file.path(data_base, "Results", config$target_year,
                        result_name, "SUPERVISED", config$scenario,
                        "11_CONSISTENCY_CHECKS")
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
