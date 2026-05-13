#' Build a deterministic burned-mapping configuration object
#'
#' @description
#' Public entry point for the deterministic workflow. Collects the spatial
#' inputs, scenario preset, output location, and advanced options needed by
#' the rest of the deterministic stages into a single normalized object of
#' class `otsufire_burned_mapping_config`. The function does not run any
#' heavy processing: it validates arguments, stores them, and resolves
#' registry and output routes.
#'
#' Pass the returned object to [detect_burned_patches()],
#' [score_burned_patches()], or [run_deterministic_pipeline()].
#'
#' Contract source: `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `DETERMINISTIC_INPUTS_FINAL.csv`, `OUTPUTS_ROUTES.csv`.
#'
#' @param scenario Character scalar. Methodological preset. One of
#'   `"balanced"` (default), `"original"`, `"lax"`, or `"restrictive"`.
#' @param change_index Raster path, `terra::SpatRaster`, or `NULL`.
#'   Annual change-index raster (RBR, dNBR, RdNBR, ...). Required before
#'   detection runs; may be `NULL` at config-build time for configs that
#'   will be used only with pre-built candidates (e.g. direct calls to
#'   [score_burned_patches()]).
#' @param vegetation_map Raster, sf object, or path. Optional vegetation or
#'   land-cover stratification layer (CORINE-compatible in 0.2.x).
#' @param burnable_mask Raster, sf object, or path. Required: binary mask
#'   of burnable area (1 burnable, 0 non-burnable).
#' @param hotspots sf POINT layer or path. Optional hotspot layer for the
#'   target year.
#' @param previous_year_burned sf POLYGON layer or path. Optional
#'   previous-year burned map for temporal conflict assessment.
#' @param reference_burned_map sf POLYGON, raster, or path. Optional burned
#'   reference used by [validate_fire_maps()] when validation is requested.
#' @param target_year Integer scalar. Target year. Required: used for
#'   naming outputs and filtering hotspots/previous-year layers.
#' @param output_dir Character scalar. Root output directory.
#'   Defaults to `tempdir()`.
#' @param run_name Character scalar. Stable run identifier used to name
#'   output folders and files. Defaults to `"deterministic_burned_map"`.
#' @param options Named list. Container for advanced options. Recognized
#'   keys include `engine_root`, `deterministic_seed`, `whitebox_exe`,
#'   `gdalwarp_path`, `tool_paths`, `vegetation_reclass`, `strata_rules`,
#'   `hotspot_fields`, `target_crs`, `tolerances`.
#'
#' @return An S3 object of class `otsufire_burned_mapping_config` (a named
#'   list). Stable fields: `scenario`, `target_year`, `inputs`,
#'   `output_routes`, `options`, `engine_root`, `deterministic_seed`,
#'   `tool_paths`, `registry_path`.
#'
#' @family workflow
#' @export
build_burned_mapping_config <- function(
    scenario = c("balanced", "original", "lax", "restrictive"),
    change_index = NULL,
    vegetation_map = NULL,
    burnable_mask,
    hotspots = NULL,
    previous_year_burned = NULL,
    reference_burned_map = NULL,
    target_year,
    output_dir = tempdir(),
    run_name = "deterministic_burned_map",
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

  if (missing(burnable_mask) || is.null(burnable_mask)) {
    stop("'burnable_mask' is required.", call. = FALSE)
  }

  if (!is.character(run_name) || length(run_name) != 1L || !nzchar(run_name)) {
    stop("'run_name' must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("'output_dir' must be a single non-empty character string.", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  if (!is.list(options)) {
    stop("'options' must be a list().", call. = FALSE)
  }
  if (length(options) > 0L && (is.null(names(options)) || any(!nzchar(names(options))))) {
    stop("'options' must be a named list (all elements must have names).", call. = FALSE)
  }

  inputs <- list(
    change_index = .of_normalize_input_spec(change_index, "change_index",
                                            allow_null = TRUE),
    vegetation_map = .of_normalize_input_spec(vegetation_map, "vegetation_map",
                                              allow_null = TRUE),
    burnable_mask = .of_normalize_input_spec(burnable_mask, "burnable_mask",
                                              allow_null = FALSE),
    hotspots = .of_normalize_input_spec(hotspots, "hotspots",
                                         allow_null = TRUE),
    previous_year_burned = .of_normalize_input_spec(previous_year_burned,
                                                     "previous_year_burned",
                                                     allow_null = TRUE),
    reference_burned_map = .of_normalize_input_spec(reference_burned_map,
                                                     "reference_burned_map",
                                                     allow_null = TRUE)
  )

  deterministic_seed <- .of_coerce_seed(options$deterministic_seed,
                                        default = 12345L)

  tool_paths <- list(
    whitebox_exe  = options$whitebox_exe  %||% options$tool_paths$whitebox_exe,
    gdalwarp_path = options$gdalwarp_path %||% options$tool_paths$gdalwarp_path
  )

  # Block 4A: deterministic engine is now in the package itself; engine_root
  # is optional and only used for legacy diagnostics. Resolution is best-effort.
  engine_root <- options$engine_root %||% .of_find_engine_root()
  if (!is.null(engine_root)) {
    engine_root <- normalizePath(engine_root, winslash = "/", mustWork = FALSE)
    if (!dir.exists(engine_root)) engine_root <- NULL
  }

  output_routes <- .of_build_output_routes(
    output_dir = output_dir,
    target_year = target_year,
    run_name = run_name,
    scenario = scenario
  )

  registry_path <- options$registry_path %||% .of_resolve_registry_path(
    output_dir = output_dir,
    scenario = scenario,
    explicit = options$registry_path
  )

  cfg <- list(
    scenario = scenario,
    target_year = target_year,
    inputs = inputs,
    run_name = run_name,
    output_dir = output_dir,
    output_routes = output_routes,
    registry_path = registry_path,
    engine_root = engine_root,
    deterministic_seed = deterministic_seed,
    tool_paths = tool_paths,
    options = options
  )

  class(cfg) <- c("otsufire_burned_mapping_config", "list")
  cfg
}

#' @export
print.otsufire_burned_mapping_config <- function(x, ...) {
  cat("<otsufire_burned_mapping_config>\n")
  cat("  scenario     :", x$scenario, "\n")
  cat("  target_year  :", x$target_year, "\n")
  cat("  run_name     :", x$run_name, "\n")
  cat("  output_dir   :", x$output_dir, "\n")
  cat("  engine_root  :", x$engine_root %||% "<unresolved>", "\n")
  cat("  deterministic_seed:", x$deterministic_seed, "\n")
  cat("  registry_path:", x$registry_path %||% "<auto>", "\n")
  cat("  inputs       :\n")
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

.of_normalize_input_spec <- function(x, name, allow_null) {
  if (is.null(x)) {
    if (!allow_null) {
      stop(sprintf("'%s' is required and cannot be NULL.", name), call. = FALSE)
    }
    return(NULL)
  }

  if (inherits(x, "SpatRaster")) {
    return(list(type = "SpatRaster", path = NA_character_, value = x))
  }
  if (inherits(x, "SpatVector")) {
    return(list(type = "SpatVector", path = NA_character_, value = x))
  }
  if (inherits(x, "sf")) {
    return(list(type = "sf", path = NA_character_, value = x))
  }
  if (is.character(x) && length(x) == 1L && nzchar(x)) {
    # Existence is not enforced at config-build time: users may build a config
    # ahead of a batch that materializes inputs later. Downstream stages
    # re-validate file existence before use.
    return(list(
      type = "path",
      path = normalizePath(x, winslash = "/", mustWork = FALSE),
      value = NULL
    ))
  }

  stop(sprintf(
    "'%s' must be a SpatRaster, SpatVector, sf object, or a single character path.",
    name
  ), call. = FALSE)
}

.of_coerce_seed <- function(x, default = 12345L) {
  if (is.null(x)) return(as.integer(default))
  y <- suppressWarnings(as.integer(x)[1L])
  if (is.na(y)) {
    stop("options$deterministic_seed must be coercible to an integer.", call. = FALSE)
  }
  y
}

.of_find_engine_root <- function() {
  candidates <- unique(c(
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd()))
  ))
  for (cand in candidates) {
    cur <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    repeat {
      probe <- file.path(cur, "2_SCRIPTS", "00_FUNCTIONS",
                         "02_DETERMINISTIC_FUNCTIONS")
      if (dir.exists(probe)) return(probe)
      # Also accept being launched from inside 2_SCRIPTS/
      probe2 <- file.path(cur, "00_FUNCTIONS", "02_DETERMINISTIC_FUNCTIONS")
      if (dir.exists(probe2)) return(probe2)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

.of_build_output_routes <- function(output_dir, target_year, run_name, scenario) {
  base <- file.path(output_dir, as.character(target_year), run_name,
                    "DETERMINISTIC", scenario)
  list(
    base          = base,
    grow_dir      = file.path(base, "01_GROW"),
    refine_dir    = file.path(base, "02_REFINE"),
    scoring_dir   = file.path(base, "05_DECISIONS"),
    validation_dir = file.path(base, "06_VALIDATION"),
    timing_dir    = file.path(base, "99_LOGS_TIMING"),
    # Canonical file names per DETERMINISTIC_OUTPUTS_FINAL.csv.
    grow_vector           = file.path(base, "01_GROW",
                                       paste0("BA_", target_year,
                                              "_OTSUGROW_CORI_ECOREG.shp")),
    refine_merged_gpkg    = file.path(base, "02_REFINE",
                                       paste0("BA_", target_year,
                                              "_REFINE_MERGED.gpkg")),
    internal_decisions    = file.path(base, "05_DECISIONS",
                                       "internal_decisions.gpkg"),
    reference_decisions   = file.path(base, "05_DECISIONS",
                                       "reference_decisions.gpkg"),
    validation_workbook   = file.path(base, "06_VALIDATION",
                                       paste0("validation_ALL_", target_year,
                                              "_", scenario, "_res30.xlsx")),
    timing_csv            = file.path(base, "99_LOGS_TIMING",
                                       paste0(target_year, "_timing_steps.csv"))
  )
}

.of_resolve_registry_path <- function(output_dir, scenario, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = FALSE))
  }
  scenario_dir <- toupper(scenario)

  # Scenario-specific registry convention:
  #   <results_root>/<SCENARIO_UPPER>/RBR_TRAINING_REGISTRY.gpkg
  # If the caller's output_dir is already under a "Results" folder (either
  # ending in .../Results or already containing a .../Results/... segment),
  # we reuse that segment instead of appending a second "Results/".
  out_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  parts <- strsplit(out_norm, "/", fixed = TRUE)[[1L]]

  results_root <- NULL
  res_idx <- which(parts == "Results")
  if (length(res_idx) > 0L) {
    # Use the deepest "Results" segment as the anchor.
    anchor <- tail(res_idx, 1L)
    results_root <- paste(parts[seq_len(anchor)], collapse = "/")
  } else {
    # No "Results" segment — append one under output_dir.
    results_root <- file.path(out_norm, "Results")
  }

  file.path(results_root, scenario_dir, "RBR_TRAINING_REGISTRY.gpkg")
}
