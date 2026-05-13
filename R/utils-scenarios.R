# Scenario presets. Internal only.
#
# Scope (Phase A — stubs only):
# - .scenario_preset(): resolve one of "balanced", "original", "lax",
#   "restrictive" into its preset parameter bundle, loaded from
#   inst/extdata/scenario_presets/. "balanced" is the primary manuscript
#   scenario per handoff Section 45.

.scenario_preset <- function(scenario = c("balanced", "original", "lax", "restrictive")) {
  .NotYetImplemented()
}

.list_available_scenarios <- function() {
  .NotYetImplemented()
}
