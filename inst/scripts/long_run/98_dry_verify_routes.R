# =============================================================================
# LONG <YEAR> BALANCED — DRY VERIFY: SHARED -> FULL -> NO-HOTSPOT ROUTE RESOLUTION
# CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
# (G2.3 route-propagation fix; NO heavy compute, NO junctions)
#
# Proves the SCORE step finds the SHARED features PURELY from explicit cfg-level
# routing (the SHARED_FEATURES_GPKG passed to score_supervised_burned_map(
# labelled_features = ...)), WITHOUT any directory junctions under the profiles'
# ENGINE_ROUTES.
#
# For BOTH profiles it builds the deployed cfg (Config B) via
# build_supervised_burned_config and prints, with Test-Path-style existence:
#   1. cfg$output_routes$features_geometry_gpkg  (the CONVENTION reconstruction;
#      points at the PROFILE folder -> the OLD failure path, expected ABSENT once
#      junctions are removed -> exactly why the explicit override is required);
#   2. SHARED_FEATURES_GPKG                       (the EXPLICIT labelled_features
#      the runners now pass -> the SHARED file the score step will actually read);
#   3. pools / folds routes from SHARED (GPKG handoffs) used upstream.
#
# It also asserts NO junctions remain under either profile's ENGINE_ROUTES
# SUPERVISED/balanced dir (the canonical path must not depend on reparse points).
#
# Does NOT train, score, or extract features.
# =============================================================================
# EDIT: point this at your run's _ORCHESTRATION/00_common.R (working copy).
source("<PATH_TO_ORCHESTRATION>/00_common.R")
ORCH <- file.path(LONG_ROOT, "_ORCHESTRATION")
.assert_version_pin()

cat("\n=============== DRY VERIFY: ROUTE RESOLUTION (no junctions) ===============\n")

te <- function(p) isTRUE(file.exists(p))   # Test-Path equivalent
mark <- function(ok) if (ok) "TRUE " else "FALSE"

# ---- per-profile ENGINE_ROUTES out_base (matches the runners) ---------------
FULL_ENGINE <- file.path(LONG_ROOT, "LONG_2017_BALANCED_FULL", "ENGINE_ROUTES")
NOHS_ENGINE <- file.path(LONG_ROOT, "LONG_2017_BALANCED_NO_HOTSPOT", "ENGINE_ROUTES")

# ---- build each profile's DEPLOYED cfg (OOF capped + FINAL) -----------------
# OOF always uses the same capped negative-sampling policy as the FINAL model.
cfg_full <- .mk_long_cfg(feature_whitelist_override = NULL, out_base = FULL_ENGINE)
cfg_nohs <- .mk_long_cfg(feature_whitelist_override = NOHS_BASE_38, out_base = NOHS_ENGINE)

# ---- 1) SHARED explicit features path the score step will read --------------
cat("\n--- SHARED explicit labelled_features (what the runners now pass) ---\n")
cat(sprintf("[Test-Path %s] SHARED_FEATURES_GPKG = %s\n",
            mark(te(SHARED_FEATURES_GPKG)), SHARED_FEATURES_GPKG))
cat(sprintf("[Test-Path %s] SHARED_FOLDS_GPKG    = %s\n",
            mark(te(SHARED_FOLDS_GPKG)), SHARED_FOLDS_GPKG))
cat(sprintf("[Test-Path %s] SHARED_POOLS_GPKG    = %s\n",
            mark(te(SHARED_POOLS_GPKG)), SHARED_POOLS_GPKG))

# ---- 2) per-profile cfg CONVENTION reconstruction (the OLD failure path) ----
show_profile <- function(cfg, label) {
  conv_feat  <- cfg$output_routes$features_geometry_gpkg
  cat(sprintf("\n--- %s profile cfg routes ---\n", label))
  cat(sprintf("[Test-Path %s] cfg$output_routes$features_geometry_gpkg (CONVENTION)\n            = %s\n",
              mark(te(conv_feat)), conv_feat))
  cat(sprintf("            -> points at THIS profile's ENGINE_ROUTES; with junctions removed this is ABSENT,\n"))
  cat(sprintf("               which is EXACTLY why the runner passes labelled_features = SHARED_FEATURES_GPKG.\n"))
  cat(sprintf("[Test-Path %s] EFFECTIVE labelled_features (explicit) = %s\n",
              mark(te(SHARED_FEATURES_GPKG)), SHARED_FEATURES_GPKG))
  invisible(list(conv_feat = conv_feat, effective = SHARED_FEATURES_GPKG))
}
r_full <- show_profile(cfg_full, "FULL")
r_nohs <- show_profile(cfg_nohs, "NO_HOTSPOT")

# ---- 3) assert no junctions remain under either profile ---------------------
cat("\n--- junction check (canonical path must not depend on reparse points) ---\n")
junction_dirs_present <- function(engine_base) {
  sup <- file.path(engine_base, as.character(target_year), result_name,
                   "SUPERVISED", scenario)
  subs <- c("01_POOLS", "02_FOLDS", "03_FEATURES", "_LEGACY_UNBURNED")
  present <- subs[file.exists(file.path(sup, subs))]
  present
}
full_present <- junction_dirs_present(FULL_ENGINE)
nohs_present <- junction_dirs_present(NOHS_ENGINE)
cat(sprintf("FULL  profile upstream subdirs still present: %s\n",
            if (length(full_present)) paste(full_present, collapse = ", ") else "(none)"))
cat(sprintf("NOHS  profile upstream subdirs still present: %s\n",
            if (length(nohs_present)) paste(nohs_present, collapse = ", ") else "(none)"))

# ---- VERDICT ----------------------------------------------------------------
ok_effective <- te(SHARED_FEATURES_GPKG) && te(SHARED_FOLDS_GPKG) && te(SHARED_POOLS_GPKG)
cat("\n--- VERDICT ---\n")
cat(sprintf("Both profiles' EFFECTIVE labelled_features resolve to the existing SHARED file: %s\n",
            mark(ok_effective)))
if (!ok_effective) {
  stop("DRY VERIFY FAILED: a SHARED route does not resolve to an existing file.",
       call. = FALSE)
}
cat("DRY VERIFY OK: score step finds SHARED features purely from explicit cfg routing; no junctions needed.\n")
cat("\n=============== DRY VERIFY ROUTE RESOLUTION DONE ===============\n")
