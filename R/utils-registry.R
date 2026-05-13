# Burned-like registry path resolution. Internal only.
#
# Scenario-specific path per handoff Section 8. No legacy-shared fallback:
# scenario separation is enforced (HANDOFF Section 43).

.resolve_registry_path <- function(results_root, scenario,
                                   basename = "RBR_TRAINING_REGISTRY.gpkg") {
  if (missing(results_root) || is.null(results_root) || !nzchar(results_root)) {
    stop(".resolve_registry_path(): 'results_root' is required.", call. = FALSE)
  }
  if (missing(scenario) || is.null(scenario) || !nzchar(scenario)) {
    stop(".resolve_registry_path(): 'scenario' is required.", call. = FALSE)
  }
  file.path(results_root, toupper(scenario), basename)
}

.registry_summary_paths <- function(registry_path) {
  if (is.null(registry_path) || !nzchar(registry_path)) return(NULL)
  dir <- dirname(registry_path)
  list(
    ranked_csv  = file.path(dir, "burned_high_conf_registry_ranked.csv"),
    summary_txt = file.path(dir, "burned_high_conf_registry_summary.txt")
  )
}
