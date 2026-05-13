# Phase 2A four-bucket unburned taxonomy helpers. Internal only.
#
# Closed 2026-04-17 (handoff Sections 32-38). Buckets:
# - contextual_exclusion         cap = 0.25x
# - spectral_hard_negative       cap = 1x, supply=0 in modern years
# - random_background_sampled    cap = 1x
# - otsu_unburned_sampled        cap = 1x
#
# Scope (Phase A — stubs only):
# - .build_all_sources_pool(): assemble unburned pool from the four buckets.
# - .bucket_cap(): return the cap multiplier for a given bucket name.
# - .resolve_unburned_strategy(): map a strategy name to its implementation
#   (all_sources is the operational default; deterministic_direct is the
#   canonical equivalence baseline; historical_registry_augment was rejected).

.build_all_sources_pool <- function(inputs, caps, year, scenario) {
  .NotYetImplemented()
}

.bucket_cap <- function(bucket = c("contextual_exclusion",
                                   "spectral_hard_negative",
                                   "random_background_sampled",
                                   "otsu_unburned_sampled")) {
  .NotYetImplemented()
}

.resolve_unburned_strategy <- function(strategy = c("all_sources", "deterministic_direct")) {
  .NotYetImplemented()
}
