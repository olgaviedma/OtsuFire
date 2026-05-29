#' Build a configuration object for the supervised burned-area workflow
#'
#' @description
#' Public entry point of the one-year supervised workflow. Collects the
#' deterministic decision layer, raster support layers, optional hotspot and
#' reference inputs, scenario-driven negative-pool policy, output location,
#' and advanced options into a single normalized configuration object of
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
#'   change-index raster for persistence-style features.
#' @param hotspots sf POINT layer or path. Optional hotspot layer for the
#'   target year.
#' @param reference_burned_map sf POLYGON / raster / path. Optional external
#'   burned reference used by `validate_fire_maps()` on thresholded
#'   supervised outputs.
#' @param burned_like_registry_path Character. Path to the scenario-specific
#'   multiyear burned-like training registry GPKG. If `NULL`, defaults to
#'   `<output_dir>/Results/<SCENARIO_UPPER>/RBR_TRAINING_REGISTRY.gpkg`.
#' @param target_year Integer scalar. Target year. Required.
#' @param output_dir Character. Root output directory. Defaults to
#'   `tempdir()`.
#' @param run_name Character. Stable identifier for the supervised run.
#'   Defaults to `"supervised_burned_map"`.
#' @param min_burned_pool_n Integer. Sparse-year guard. The pipeline aborts
#'   before training when fewer than this many keep-class polygons remain
#'   after QA. Default `5L` (HANDOFF Section 51).
#' @param options Named list of advanced options. Recognized keys include
#'   `engine_root`, `supervised_engine_root`, `negative_pool_policy`
#'   (`"all_sources"` | `"deterministic_direct"`),
#'   `registry_prob_threshold`, `registry_top_fraction`,
#'   `legacy_*`, `currentyear_*`, `unb_verbose`, `result_name`,
#'   `scripts_root`, `data_base`, plus overrides forwarded to the engine.
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
    burned_like_registry_path = NULL,
    target_year,
    output_dir = tempdir(),
    run_name = "supervised_burned_map",
    min_burned_pool_n = 5L,
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
                                                     allow_null = TRUE)
  )

  # Locate both engine roots (deterministic helper reused; supervised dir
  # under 00_FUNCTIONS/04_SUPERVISED_FUNCTIONS).
  engine_root <- options$engine_root %||% .of_find_engine_root()
  if (!is.null(engine_root)) {
    engine_root <- normalizePath(engine_root, winslash = "/", mustWork = FALSE)
    if (!dir.exists(engine_root)) {
      stop("options$engine_root does not exist: ", engine_root, call. = FALSE)
    }
  }
  sup_engine_root <- options$supervised_engine_root %||%
                      .of_find_supervised_engine_root()
  if (!is.null(sup_engine_root)) {
    sup_engine_root <- normalizePath(sup_engine_root, winslash = "/",
                                      mustWork = FALSE)
    if (!dir.exists(sup_engine_root)) {
      stop("options$supervised_engine_root does not exist: ",
           sup_engine_root, call. = FALSE)
    }
  }

  # Block 4B: engine is now in-package. scripts_root is only used for the
  # optional consistency check (which still lives under 00_USAGE/03_ANALYSIS).
  # Resolution is best-effort; a missing scripts_root is NOT an error.
  scripts_root <- options$scripts_root %||%
                   .of_find_supervised_scripts_root()
  if (!is.null(scripts_root)) {
    scripts_root <- normalizePath(scripts_root, winslash = "/",
                                   mustWork = FALSE)
    if (!dir.exists(scripts_root)) scripts_root <- NULL
  }

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

  # Scenario-driven negative-pool default (Phase 2A resolved 2026-04-17).
  neg_pool <- .of_resolve_negative_pool_policy(scenario, options)

  output_routes <- .of_build_supervised_output_routes(
    output_dir = output_dir, target_year = target_year,
    run_name = run_name, scenario = scenario
  )

  # Registry path: scenario-specific under the user's output_dir/Results.
  # `.of_resolve_registry_path()` was removed in handoff §N+5.1 when the
  # deterministic builder dropped registry_path; the supervised side kept a
  # stale call. Redirect to the canonical helper in utils-registry.R.
  if (is.null(burned_like_registry_path) ||
      !nzchar(burned_like_registry_path)) {
    burned_like_registry_path <- .resolve_registry_path(
      results_root = file.path(output_dir, "Results"),
      scenario     = scenario
    )
  }
  burned_like_registry_path <- normalizePath(burned_like_registry_path,
                                              winslash = "/", mustWork = FALSE)

  cfg <- list(
    scenario                 = scenario,
    target_year              = target_year,
    inputs                   = inputs,
    run_name                 = run_name,
    output_dir               = output_dir,
    output_routes            = output_routes,
    burned_like_registry_path = burned_like_registry_path,
    negative_pool_policy     = neg_pool,
    min_burned_pool_n        = min_burned_pool_n,
    engine_root              = engine_root,
    supervised_engine_root   = sup_engine_root,
    scripts_root             = scripts_root,
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
  cat("  registry_path   :", x$burned_like_registry_path %||% "<unresolved>", "\n")
  cat("  engine_root     :", x$engine_root %||% "<unresolved>", "\n")
  cat("  sup_engine_root :", x$supervised_engine_root %||% "<unresolved>", "\n")
  cat("  scripts_root    :", x$scripts_root %||% "<unresolved>", "\n")
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

#' @keywords internal
#' @noRd
.of_find_supervised_engine_root <- function() {
  candidates <- unique(c(getwd(), dirname(getwd()), dirname(dirname(getwd()))))
  for (cand in candidates) {
    cur <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    repeat {
      p1 <- file.path(cur, "2_SCRIPTS", "00_FUNCTIONS",
                      "04_SUPERVISED_FUNCTIONS")
      if (dir.exists(p1)) return(p1)
      p2 <- file.path(cur, "00_FUNCTIONS", "04_SUPERVISED_FUNCTIONS")
      if (dir.exists(p2)) return(p2)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

#' @keywords internal
#' @noRd
.of_find_supervised_scripts_root <- function() {
  candidates <- unique(c(getwd(), dirname(getwd()), dirname(dirname(getwd()))))
  for (cand in candidates) {
    cur <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "PROJECT_PATHS.R")) &&
          dir.exists(file.path(cur, "00_FUNCTIONS")) &&
          dir.exists(file.path(cur, "00_USAGE"))) {
        return(cur)
      }
      if (file.exists(file.path(cur, "2_SCRIPTS", "PROJECT_PATHS.R")) &&
          dir.exists(file.path(cur, "2_SCRIPTS", "00_FUNCTIONS")) &&
          dir.exists(file.path(cur, "2_SCRIPTS", "00_USAGE"))) {
        return(file.path(cur, "2_SCRIPTS"))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

#' @keywords internal
#' @noRd
.of_resolve_negative_pool_policy <- function(scenario, options = list()) {
  override <- options$negative_pool_policy %||%
              options$negative_pool_strategy %||% NULL
  if (!is.null(override) && identical(override, "otsu_unburned_generation")) {
    warning("'otsu_unburned_generation' is deprecated. Redirecting to 'all_sources'.",
            call. = FALSE)
    override <- "all_sources"
  }
  policy <- if (!is.null(override)) {
    override
  } else if (tolower(as.character(scenario)) %in%
             c("balanced", "strict", "restrictive")) {
    "all_sources"
  } else {
    "deterministic_direct"
  }
  if (!policy %in% c("deterministic_direct", "all_sources")) {
    stop("Unsupported negative-pool policy. Use 'deterministic_direct' or 'all_sources'.",
         call. = FALSE)
  }
  policy
}

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
