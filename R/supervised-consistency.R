#' Check consistency between Otsu-guided patch decisions and supervised outputs
#'
#' @description
#' Run the consistency check between the decision layer for Otsu-guided
#' patches and the supervised map outputs from the same run.
#'
#' Use this function after [score_supervised_burned_map()]. The same check is
#' invoked by [run_oneyear_supervised_pipeline()] when
#' `run_consistency = TRUE`.
#'
#' This check compares workflow outputs. To evaluate mapping accuracy against
#' external reference fire perimeters, use [validate_fire_maps()].
#'
#' @param deterministic_decisions Character scalar. Path to a GeoPackage
#'   containing the `internal_decisions` layer for the Otsu-guided patches.
#'   Supply a file path; in-memory `sf` objects are not supported by this
#'   interface.
#' @param final_map Character scalar. Path to the supervised map GeoPackage.
#'   Must contain both `final_map_full` and `final_map` layers.
#' @param out_dir Character scalar or `NULL`. Output directory for
#'   consistency-check products. When `NULL`, uses
#'   `config$output_routes$consistency_dir`. Supply explicitly when no
#'   configuration is provided.
#' @param config Optional `otsufire_supervised_burned_config` object created by
#'   [build_supervised_burned_config()]. Supplies the target year, scenario,
#'   and default output directory.
#' @param target_year Optional integer scalar. Target year used when
#'   `config = NULL`. Ignored when `config` is supplied.
#' @param scenario Optional character scalar. Scenario identifier used when
#'   `config = NULL`. Ignored when `config` is supplied.
#'
#' @section Required input layers:
#' Supply the decision and supervised output files corresponding to the same
#' target year and scenario.
#'
#' | Input | Required layers |
#' |---|---|
#' | `deterministic_decisions` | `internal_decisions` |
#' | `final_map` | `final_map_full`, `final_map` |
#'
#' The argument name `deterministic_decisions` is retained for API
#' compatibility. It refers to the decision layer for Otsu-guided patches.
#'
#' A GeoPackage containing only thresholded burned polygons is not a
#' substitute for the supervised map file required here.
#'
#' @section Complete and public map layers:
#' The two supervised layers have different purposes:
#' * `final_map_full` contains the complete scored candidate set.
#' * `final_map` contains the public output after temporal exclusions and
#'   public-column selection.
#'
#' These layers can legitimately have different row counts. A difference in
#' counts alone does not establish an inconsistency.
#'
#' Likewise, agreement between the initial patch decisions and supervised
#' results is not a measure of accuracy against independent observations.
#'
#' @section Run identifiers and output location:
#' When `config` is supplied, the function takes the target year and scenario
#' from that object. Explicit `target_year` and `scenario` arguments are
#' ignored.
#'
#' Without a configuration, supply `target_year`, `scenario`, and `out_dir`
#' explicitly.
#'
#' @return A list with the paths of the written artefacts
#'   (`consistency_issues_gpkg`, `consistency_summary_csv`,
#'   `consistency_summary_txt`, `pburned_overall_summary`,
#'   `pburned_by_class_input`, `pburned_by_source_set`) and `summary_row`,
#'   the one-row summary of the check (`NULL` if the check could not run).
#'
#' @seealso [build_supervised_burned_config()], [score_supervised_burned_map()],
#'   [run_oneyear_supervised_pipeline()], [validate_supervised_execution()],
#'   [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' # Check outputs using the configuration from the supervised run
#' check_supervised_consistency(
#'   deterministic_decisions = "data/internal_decisions_2022.gpkg",
#'   final_map = "results/supervised_final_map_2022.gpkg",
#'   config = cfg
#' )
#'
#' # Alternatively, supply run identifiers and the output directory
#' check_supervised_consistency(
#'   deterministic_decisions = "data/internal_decisions_2022.gpkg",
#'   final_map = "results/supervised_final_map_2022.gpkg",
#'   out_dir = "results/consistency_2022",
#'   target_year = 2022L,
#'   scenario = "balanced"
#' )
#' }
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
