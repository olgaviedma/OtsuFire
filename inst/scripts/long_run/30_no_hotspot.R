# =============================================================================
# LONG <YEAR> BALANCED — NO-HOTSPOT PROFILE (spectral-historical ablation)
# CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
# Only runs if FULL closed OK (the master gates on FULL_SUCCESS.marker).
# OOF nested_refit capped + FINAL nested_refit + score + map + EFFIS.
# Optional EXTRA diagnostic: OOF nested_refit FULL, reusing the SAME FINAL
# (no second FINAL). NO legacy config for no-hotspot.
# Difference vs FULL is ONLY explicit configuration:
#   feature_whitelist_override = NOHS_BASE_38 (38 base = 50 minus 12 hotspot)
#   use_hotspots = FALSE (belt-and-suspenders; does NOT change pools/labels --
#   pools/folds/features are the SHARED upstream; only the model feature set
#   differs). The SHARED feature GPKG is REUSED; nothing upstream is regenerated.
# =============================================================================
# EDIT: point this at your run's _ORCHESTRATION/00_common.R (working copy).
source("<PATH_TO_ORCHESTRATION>/00_common.R")

PROF <- file.path(LONG_ROOT, "LONG_2017_BALANCED_NO_HOTSPOT")
SHARED <- file.path(LONG_ROOT, "SHARED")
FULL_DIR <- file.path(LONG_ROOT, "LONG_2017_BALANCED_FULL")
dir.create(PROF, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(PROF, "no_hotspot.log"); log <- make_logger(LOG)
t0 <- Sys.time()
log("==== NO-HOTSPOT PROFILE START ====")
.assert_version_pin()

# Gate: SHARED + FULL must have succeeded.
stopifnot(file.exists(file.path(SHARED, "SHARED_SUCCESS.marker")),
          file.exists(file.path(FULL_DIR, "FULL_SUCCESS.marker")))

# G2.3 route-propagation fix: FAIL-FAST validate + SHA256 every SHARED upstream
# input BEFORE any heavy compute (same contract as the FULL profile). All shared
# products must resolve from their SHARED absolute paths -- no junctions.
shared_provenance <- validate_shared_inputs(
  log = log, out_csv = file.path(PROF, "SHARED_INPUTS_PROVENANCE.csv"))
log("[nohs] SHARED inputs validated + hashed (", nrow(shared_provenance),
    " products); features GPKG = ", SHARED_FEATURES_GPKG)

pools <- readRDS(SHARED_POOLS_RDS)
folds <- readRDS(SHARED_FOLDS_RDS)
feats <- readRDS(SHARED_FEATS_RDS)
log("[nohs] loaded SHARED upstream (REUSED; nothing regenerated).")

.warns <- character(0)
wh <- function(w) { .warns[[length(.warns) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }

# NO-HOTSPOT cfgs: explicit feature_whitelist_override (38 non-hotspot base).
# build_supervised_burned_config validates the override as a SUBSET of the
# canonical .supervised_feature_cols and records provenance "user". We do NOT
# touch the package whitelist.
stopifnot(length(NOHS_BASE_38) == 38L)
NOHS_ENGINE <- file.path(PROF, "ENGINE_ROUTES")
dir.create(NOHS_ENGINE, recursive = TRUE, showWarnings = FALSE)
cfg_B <- .mk_long_cfg("nested_refit", "capped", feature_whitelist_override = NOHS_BASE_38, out_base = NOHS_ENGINE)
cfg_C <- .mk_long_cfg("nested_refit", "full",   feature_whitelist_override = NOHS_BASE_38, out_base = NOHS_ENGINE)
saveRDS(list(B = cfg_B, C = cfg_C), file.path(PROF, "cfgs_no_hotspot.rds"))
stopifnot(identical(sort(cfg_B$train_control$feature_whitelist_override), sort(NOHS_BASE_38)))
log("[nohs] cfgs built; feature_whitelist_override n=",
    length(cfg_B$train_control$feature_whitelist_override),
    " (provenance=",
    cfg_B$resolved_params_provenance$train_control$feature_whitelist_override %||% "?", ")")

dir_of <- function(sub) { d <- file.path(PROF, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

# ---- OOF nested_refit capped + FINAL nested_refit ---------------------------
log("[nohs] OOF nested_refit capped + FINAL nested_refit (use_hotspots=FALSE) ...")
B_oof_dir <- dir_of("CONFIG_B/05_OOF"); B_mat_dir <- dir_of("CONFIG_B/04_MATRIX")
B_fm_dir  <- dir_of("CONFIG_B/07_FINAL_MODEL"); B_sc_dir <- dir_of("CONFIG_B/08_SCORED"); B_map_dir <- dir_of("CONFIG_B/09_FINAL_MAP")
B_oof <- withCallingHandlers(run_oof_diagnostics(
  train_features = feats$train_features, scoring_features = feats$scoring_features,
  config = cfg_B, fold_cols = FOLD_COLS, out_dir = B_oof_dir, matrix_dir = B_mat_dir,
  labelled_gpkg = folds$train_with_folds_gpkg, labelled_layer = "train_with_folds"), warning = wh)
B_model <- withCallingHandlers(train_final_burned_model(
  train_features = feats$train_features, config = cfg_B, oof_agg = B_oof$oof_agg_csv,
  out_dir = B_fm_dir, overwrite = TRUE, verbose = TRUE), warning = wh)
saveRDS(B_oof, file.path(PROF, "B_oof.rds")); saveRDS(B_model, file.path(PROF, "B_model.rds"))
log("[nohs] OOF+FINAL done; FINAL best_iteration=", B_model$recipe$training$best_iteration %||% NA)

# ---- OPTIONAL EXTRA DIAGNOSTIC: OOF nested FULL (reuse same FINAL) -----------
log("[nohs] EXTRA diagnostic: OOF nested_refit FULL (reuses the SAME FINAL) ...")
C_oof_dir <- dir_of("CONFIG_C_oof_full_diagnostic/05_OOF")
C_mat_dir <- dir_of("CONFIG_C_oof_full_diagnostic/04_MATRIX")
C_oof <- withCallingHandlers(run_oof_diagnostics(
  train_features = feats$train_features, scoring_features = feats$scoring_features,
  config = cfg_C, fold_cols = FOLD_COLS, out_dir = C_oof_dir, matrix_dir = C_mat_dir,
  labelled_gpkg = folds$train_with_folds_gpkg, labelled_layer = "train_with_folds"), warning = wh)
saveRDS(C_oof, file.path(PROF, "C_oof.rds"))
log("[nohs] EXTRA OOF-full diagnostic done (no FINAL retrain).")

# ---- SCORE + MAP (deploy the nested_refit FINAL) ----------------------------
# G2.3 route-propagation fix: pass the SHARED labelled-features GPKG EXPLICITLY
# (same as the FULL profile) so the score step does NOT reconstruct the features
# path under THIS profile's ENGINE_ROUTES. No junctions involved.
log("[nohs] SCORE + MAP ...")
stopifnot(file.exists(SHARED_FEATURES_GPKG))
log("[nohs] scoring labelled_features <- SHARED ", SHARED_FEATURES_GPKG)
B_scored <- withCallingHandlers(score_supervised_burned_map(
  scoring_features = feats$scoring_features, model = B_model$model, recipe = B_model$recipe,
  config = cfg_B, oof_summary = B_oof$labeled_oof_summary_gpkg,
  labelled_features = SHARED_FEATURES_GPKG, export_burned_like = TRUE,
  out_map_dir = B_map_dir, out_score_dir = B_sc_dir, overwrite = TRUE, verbose = TRUE), warning = wh)
saveRDS(B_scored, file.path(PROF, "B_scored.rds"))
fm <- B_scored$final_map_full
log("[nohs] scored map rows=", nrow(fm), " (scoring rows=", nrow(feats$scoring_features), ")")

# ---- FEATURE FINGERPRINT (this profile's own schema; differs from FULL) -----
fn_nohs <- B_model$model$feature_names
base_nohs <- fn_nohs[!grepl("_isNA$", fn_nohs)]; isna_nohs <- fn_nohs[grepl("_isNA$", fn_nohs)]
fp_nohs <- digest::digest(paste(sort(fn_nohs), collapse = "|"), algo = "md5")
removed <- setdiff(FULL_BASE_50, base_nohs)
writeLines(c(
  "profile: NO_HOTSPOT",
  paste0("n_base: ", length(base_nohs)),
  paste0("n_isNA: ", length(isna_nohs)),
  paste0("final_feature_order_n: ", length(fn_nohs)),
  paste0("feature_schema_fingerprint: ", fp_nohs),
  paste0("removed_vs_FULL_n: ", length(removed)),
  "", "=== removed (hotspot-lineage) ===", sort(removed),
  "", "=== base features ===", sort(base_nohs),
  "", "=== isNA companions ===", sort(isna_nohs)),
  file.path(PROF, "NO_HOTSPOT_feature_list.txt"))
log("[nohs] feature schema: base=", length(base_nohs), " isNA=", length(isna_nohs),
    " removed=", length(removed), " fingerprint=", fp_nohs)

# ---- EFFIS (dedicated, memory-safe) -----------------------------------------
log("[nohs] EFFIS validation ...")
effis <- run_long_effis(fm, file.path(PROF, "VALIDATION_EFFIS"), log)
if (!is.null(effis)) {
  saveRDS(effis, file.path(PROF, "effis.rds"))
  m <- effis$metrics %||% effis
  df <- tryCatch(as.data.frame(m), error = function(e) NULL)
  if (!is.null(df)) {
    utils::write.csv(df, file.path(PROF, "EFFIS_metrics.csv"), row.names = FALSE)
    log("[nohs] EFFIS metrics written.")
  }
}

# ---- SUMMARY + manifest -----------------------------------------------------
writeLines(if (length(.warns)) .warns else "(none)", file.path(PROF, "WARNINGS.txt"))
md <- c("# LONG BALANCED — NO-HOTSPOT profile (spectral-historical ablation)",
  "", sprintf("- Snapshot: %s (tarball SHA256 %s)", SNAPSHOT_ID, TARBALL_SHA256),
  sprintf("- Training: nrounds_max=%d early_stop=%d seed=%d nthread=%d (canonical full).",
          NROUNDS_MAX, EARLY_STOP, SEED_BASE, XGB_NTHREAD),
  sprintf("- Caps: contextual %.2f / spectral %.2f / random %.2f / otsu %.2f.",
          CAP_CONTEXTUAL, CAP_SPECTRAL, CAP_RANDOM, CAP_OTSU),
  sprintf("- feature_whitelist_override: %d base (50 minus 12 hotspot-lineage).", length(NOHS_BASE_38)),
  sprintf("- Resolved model features: %d base + %d _isNA.", length(base_nohs), length(isna_nohs)),
  sprintf("- feature_schema_fingerprint: %s", fp_nohs),
  sprintf("- Removed (hotspot-lineage): %s", paste(sort(removed), collapse = ", ")),
  "- OOF nested_refit capped + FINAL nested_refit (DEPLOYED) + OOF-full diagnostic (reuses FINAL). No legacy.",
  "", "## EFFIS metrics (validate_fire_maps, metrics_type='all')",
  if (exists("df") && !is.null(df)) paste(utils::capture.output(print(t(df))), collapse = "\n") else "(EFFIS not produced)")
writeLines(md, file.path(PROF, "SUMMARY.md"))

dur <- as.numeric(Sys.time() - t0, units = "mins")
log(sprintf("==== NO-HOTSPOT PROFILE DONE (%.2f min) ====", dur))
writeLines(c(sprintf("NO_HOTSPOT_OK %s", format(Sys.time())),
             sprintf("n_base=%d removed=%d scored_rows=%d fingerprint=%s",
                     length(base_nohs), length(removed), nrow(fm), fp_nohs)),
           file.path(PROF, "NO_HOTSPOT_SUCCESS.marker"))
