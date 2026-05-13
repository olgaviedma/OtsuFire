#' Update the scenario-specific burned-like multiyear registry
#'
#' @description
#' Public registry-update entry point defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. In 0.2.x this is a
#' contract-shaped wrapper; real execution happens inside
#' [run_oneyear_supervised_pipeline()]. `prob_threshold` is the
#' eligibility floor; the actual retained set is the top 1% by descending
#' `p_burned` across all appended years.
#'
#' @param burned_like_scored sf or GPKG path.
#' @param burned_like_registry_path Character path to the scenario-specific
#'   registry GPKG. Created automatically under
#'   `Results/<SCENARIO_UPPER>/` when absent.
#' @param prob_threshold Numeric scalar in (0,1). Default `0.90`.
#' @param config Optional `otsufire_supervised_burned_config`.
#' @param write_outputs Logical.
#'
#' @family modular
#' @export
update_burned_like_registry <- function(burned_like_scored,
                                        burned_like_registry_path,
                                        prob_threshold = 0.90,
                                        config = NULL,
                                        write_outputs = TRUE) {
  if (missing(burned_like_scored) || is.null(burned_like_scored)) {
    stop("'burned_like_scored' is required.", call. = FALSE)
  }
  if (missing(burned_like_registry_path) ||
      is.null(burned_like_registry_path)) {
    if (is.null(config) ||
        !inherits(config, "otsufire_supervised_burned_config")) {
      stop("Supply either 'burned_like_registry_path' or a valid 'config'.",
           call. = FALSE)
    }
    burned_like_registry_path <- config$burned_like_registry_path
  }
  if (!is.character(burned_like_registry_path) ||
      length(burned_like_registry_path) != 1L ||
      !nzchar(burned_like_registry_path)) {
    stop("'burned_like_registry_path' must be a single non-empty path.",
         call. = FALSE)
  }
  if (!is.numeric(prob_threshold) || length(prob_threshold) != 1L ||
      prob_threshold < 0 || prob_threshold > 1) {
    stop("'prob_threshold' must be a single numeric in [0, 1].", call. = FALSE)
  }
  registry_dir <- dirname(burned_like_registry_path)
  list(
    registry_candidates_layer = "burned_high_conf_registry_candidates",
    registry_layer            = "burned_high_conf_registry",
    registry_path             = burned_like_registry_path,
    registry_ranked_csv       = file.path(registry_dir,
                                   "burned_high_conf_registry_ranked.csv"),
    registry_summary_txt      = file.path(registry_dir,
                                   "burned_high_conf_registry_summary.txt"),
    prob_threshold            = prob_threshold,
    note = paste0(
      "In 0.2.x the registry is maintained by run_oneyear_supervised_pipeline(); ",
      "calling update_burned_like_registry() standalone validates inputs."
    )
  )
}
