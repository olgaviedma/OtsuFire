#' Build burned and negative training pools for supervised learning
#'
#' @description
#' DEFERRED: full implementation planned for 0.4.0. Use
#' `run_oneyear_supervised_pipeline()` for operational runs.
#'
#' Pool-construction function defined by
#' `SUPERVISED_ONEYEAR_PUBLIC_FUNCTION_CONTRACTS.csv`. Reads deterministic
#' decisions, performs conservative QA relabelling, assembles the burned
#' pool and the scenario-driven negative pool (Phase 2A four-bucket
#' taxonomy), and returns the labelled training set that feeds
#' [make_spatial_folds()].
#'
#' In 0.3.x this is a contract-shaped entry point only. The heavy engine
#' execution is performed through the internal supervised dispatcher when
#' invoked from [run_oneyear_supervised_pipeline()]. Calling this
#' function standalone validates inputs and returns the contract-shaped
#' output layer paths resolved from `config$output_routes`; the files
#' themselves are written by the full orchestrator.
#'
#' @param config `otsufire_supervised_burned_config` object.
#' @param deterministic_decisions Optional sf or GPKG path. When omitted,
#'   the canonical `config$inputs$internal_decisions` is used.
#' @param write_outputs Logical. Whether to write pool GPKG + QA files.
#' @param overwrite Logical. Whether to replace existing pool outputs.
#'
#' @return A named list with contract-shaped fields: `burned_pool`,
#'   `unburned_pool`, `review_pool`, `scoring_pool`, `train_labeled`,
#'   `pools_gpkg`, `deterministic_pool_qa_*`, and
#'   `negative_pool_support_outputs`.
#'
#' @keywords internal
#' @noRd
build_supervised_training_pools <- function(config,
                                            deterministic_decisions = NULL,
                                            write_outputs = TRUE,
                                            overwrite = FALSE) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(deterministic_decisions)) {
    if (!(inherits(deterministic_decisions, c("sf", "SpatVector")) ||
          (is.character(deterministic_decisions) &&
           length(deterministic_decisions) == 1L &&
           file.exists(deterministic_decisions)))) {
      stop("'deterministic_decisions' must be sf, SpatVector, or an existing path.",
           call. = FALSE)
    }
  }
  .of_check_input_file(config$inputs$internal_decisions, "internal_decisions")

  pools_gpkg <- config$output_routes$pools_gpkg
  prefix     <- sprintf("%d_%s", config$target_year, config$scenario)
  pools_dir  <- config$output_routes$pools_dir

  list(
    burned_pool       = NULL,
    unburned_pool     = NULL,
    review_pool       = NULL,
    scoring_pool      = NULL,
    train_labeled     = NULL,
    pools_gpkg        = pools_gpkg,
    deterministic_pool_qa_audited    = file.path(
      pools_dir, paste0(prefix, "_deterministic_pool_qa_audited_internal.gpkg")
    ),
    deterministic_pool_qa_summary    = file.path(
      pools_dir, paste0(prefix, "_deterministic_pool_qa_summary.csv")
    ),
    deterministic_pool_qa_transitions = file.path(
      pools_dir, paste0(prefix, "_deterministic_pool_qa_transitions.csv")
    ),
    deterministic_pool_qa_reasons    = file.path(
      pools_dir, paste0(prefix, "_deterministic_pool_qa_reasons.csv")
    ),
    negative_pool_support_outputs    = NULL,
    note = paste0(
      "In 0.2.x pools are materialized by run_oneyear_supervised_pipeline(); ",
      "calling build_supervised_training_pools() standalone validates inputs ",
      "and returns the canonical output paths."
    )
  )
}
