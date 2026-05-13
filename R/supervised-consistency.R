#' Compare deterministic and supervised outputs for one run
#'
#' @description
#' Public consistency-check entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. As of Block 6 this
#' wrapper executes the consistency check directly via the package-internal
#' engine; the same engine is also invoked from
#' [run_oneyear_supervised_pipeline()] when `run_consistency = TRUE`.
#'
#' @param deterministic_decisions GPKG path to the deterministic
#'   `internal_decisions.gpkg` (must contain layer `internal_decisions`).
#'   sf objects are not supported in 0.2.x — pass a path instead.
#' @param final_map GPKG path to the supervised final-map file
#'   (must contain layers `final_map_full` and `final_map`).
#' @param out_dir Output folder; defaults to
#'   `config$output_routes$consistency_dir`.
#' @param config Optional `otsufire_supervised_burned_config`. When
#'   supplied, `target_year` and `scenario` are taken from it. If NULL,
#'   the wrapper accepts a numeric `target_year` (column 1 of the GPKG
#'   path's parent year folder) — see Examples in the package vignette.
#' @param target_year Optional integer; ignored when `config` is supplied.
#' @param scenario Optional character; ignored when `config` is supplied.
#'
#' @family modular
#' @export
check_supervised_consistency <- function(deterministic_decisions, final_map,
                                         out_dir = NULL, config = NULL,
                                         target_year = NULL,
                                         scenario = NULL) {
  if (missing(deterministic_decisions) || is.null(deterministic_decisions)) {
    stop("'deterministic_decisions' is required.", call. = FALSE)
  }
  if (missing(final_map) || is.null(final_map)) {
    stop("'final_map' is required.", call. = FALSE)
  }
  if (!is.character(deterministic_decisions) ||
      length(deterministic_decisions) != 1L) {
    stop("'deterministic_decisions' must be a GPKG path (length-1 character).",
         call. = FALSE)
  }
  if (!is.character(final_map) || length(final_map) != 1L) {
    stop("'final_map' must be a GPKG path (length-1 character).",
         call. = FALSE)
  }
  if (!is.null(config) &&
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be an otsufire_supervised_burned_config object.",
         call. = FALSE)
  }
  if (is.null(out_dir)) {
    if (is.null(config)) {
      stop("Supply either 'out_dir' or a valid 'config'.", call. = FALSE)
    }
    out_dir <- config$output_routes$consistency_dir
  }
  if (!is.null(config)) {
    target_year <- config$target_year
    scenario    <- config$scenario
  }
  if (is.null(target_year) || is.null(scenario)) {
    stop("'target_year' and 'scenario' must be supplied either directly ",
         "or via 'config'.", call. = FALSE)
  }
  ns <- asNamespace("OtsuFire")
  res <- tryCatch(
    ns$.of_run_consistency_check(
      target_year = as.integer(target_year),
      scenario    = as.character(scenario),
      internal_decisions_gpkg = deterministic_decisions,
      final_map_gpkg = final_map,
      out_dir = out_dir
    ),
    error = function(e) NULL
  )
  prefix <- sprintf("%d_%s", as.integer(target_year), as.character(scenario))
  c(
    list(
      consistency_issues       = NULL,
      consistency_issues_gpkg  = file.path(out_dir,
                                  paste0(prefix, "_consistency_issues.gpkg")),
      consistency_summary_csv  = file.path(out_dir,
                                  paste0(prefix, "_consistency_summary.csv")),
      consistency_summary_txt  = file.path(out_dir,
                                  paste0(prefix, "_consistency_summary.txt")),
      pburned_overall_summary  = file.path(out_dir,
                                  paste0(prefix, "_pburned_overall_summary.csv")),
      pburned_by_class_input   = file.path(out_dir,
                                  paste0(prefix, "_pburned_by_class_input.csv")),
      pburned_by_source_set    = file.path(out_dir,
                                  paste0(prefix, "_pburned_by_source_set.csv"))
    ),
    list(summary_row = res)
  )
}
