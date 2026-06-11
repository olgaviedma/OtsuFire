# =============================================================================
# LONG <YEAR> BALANCED — FULL PROFILE (thermal-enhanced; all features incl hotspots)
# CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
# OtsuFire always uses ONE training procedure (inner-early-stopping selection +
# full-data refit); there is no training-protocol choice. The OOF diagnostic uses
# the SAME capped negative-sampling policy as the FINAL model, applied
# independently within each training fold (there is no oof_sampling choice):
#   DEPLOYED = OOF (capped) + FINAL  (the deployed FINAL)
# Then: score + map + EFFIS using the DEPLOYED FINAL.
# Consumes the SHARED pools/folds/features. Builds its OWN matrices/recipes/
# models/predictions/fingerprints/scoring/map (NOT shared).
# =============================================================================
# EDIT: point this at your run's _ORCHESTRATION/00_common.R (working copy).
source("<PATH_TO_ORCHESTRATION>/00_common.R")

PROF <- file.path(LONG_ROOT, "LONG_2017_BALANCED_FULL")
SHARED <- file.path(LONG_ROOT, "SHARED")
dir.create(PROF, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(PROF, "full.log"); log <- make_logger(LOG)
t0 <- Sys.time()
log("==== FULL PROFILE START ====")
.assert_version_pin()

# Gate: SHARED must have succeeded.
stopifnot(file.exists(file.path(SHARED, "SHARED_SUCCESS.marker")))

# G2.3 route-propagation fix: FAIL-FAST validate + SHA256 every SHARED upstream
# input BEFORE any heavy compute. The shared features GPKG, pools/folds GPKGs and
# RDS handoffs must ALL resolve from their SHARED absolute paths (NO junctions,
# NO convention reconstruction). Aborts here with a clear error if anything is
# unresolvable or unhashable, long before training/scoring.
shared_provenance <- validate_shared_inputs(
  log = log, out_csv = file.path(PROF, "SHARED_INPUTS_PROVENANCE.csv"))
log("[full] SHARED inputs validated + hashed (", nrow(shared_provenance),
    " products); features GPKG = ", SHARED_FEATURES_GPKG)

pools <- readRDS(SHARED_POOLS_RDS)
folds <- readRDS(SHARED_FOLDS_RDS)
feats <- readRDS(SHARED_FEATS_RDS)
n_burned <- sum(as.character(pools$train_labeled$class) == "burned", na.rm = TRUE)
log("[full] loaded SHARED upstream; n_burned=", n_burned)

.warns <- character(0)
wh <- function(w) { .warns[[length(.warns) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }

# FULL profile cfgs (all features; NO whitelist override; hotspots ON).
# Per-profile engine route base (isolated from the other profile + from SHARED).
FULL_ENGINE <- file.path(PROF, "ENGINE_ROUTES")
dir.create(FULL_ENGINE, recursive = TRUE, showWarnings = FALSE)
cfg_B <- .mk_long_cfg(feature_whitelist_override = NULL, out_base = FULL_ENGINE)
saveRDS(list(B = cfg_B), file.path(PROF, "cfgs_full.rds"))
log("[full] cfg built (DEPLOYED: OOF capped + FINAL).")

dir_of <- function(sub) { d <- file.path(PROF, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

# ---- DEPLOYED: OOF capped + FINAL (the deployed model) ----------------------
log("[full] DEPLOYED (OOF capped + FINAL) ...")
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
log("[full] DEPLOYED done; FINAL best_iteration=", B_model$recipe$training$best_iteration %||% NA)

# ---- DEPLOY THE FINAL: SCORE + MAP ------------------------------------------
# G2.3 route-propagation fix: pass the SHARED labelled-features GPKG EXPLICITLY
# via `labelled_features` so the score step does NOT fall back to
# cfg_B$output_routes$features_geometry_gpkg (which reconstructs to THIS profile's
# ENGINE_ROUTES, where no features GPKG exists -- the original LONG-run failure).
# The scoring universe is the in-memory feats$scoring_features (already explicit)
# and oof_summary is this profile's own OOF (explicit). No junctions involved.
stopifnot(file.exists(SHARED_FEATURES_GPKG))
log("[full] scoring labelled_features <- SHARED ", SHARED_FEATURES_GPKG)
B_scored <- withCallingHandlers(score_supervised_burned_map(
  scoring_features = feats$scoring_features, model = B_model$model, recipe = B_model$recipe,
  config = cfg_B, oof_summary = B_oof$labeled_oof_summary_gpkg,
  labelled_features = SHARED_FEATURES_GPKG, export_burned_like = TRUE,
  out_map_dir = B_map_dir, out_score_dir = B_sc_dir, overwrite = TRUE, verbose = TRUE), warning = wh)
saveRDS(B_scored, file.path(PROF, "B_scored.rds"))
fm <- B_scored$final_map_full
log("[full] scored map rows=", nrow(fm), " (scoring rows=", nrow(feats$scoring_features), ")")

# ---- FEATURE FINGERPRINT (within-profile parity must hold OOF=FINAL=scoring) -
fn_full <- B_model$model$feature_names
base_full <- fn_full[!grepl("_isNA$", fn_full)]; isna_full <- fn_full[grepl("_isNA$", fn_full)]
fp_full <- digest::digest(paste(sort(fn_full), collapse = "|"), algo = "md5")
writeLines(c(
  paste0("profile: FULL"),
  paste0("n_base: ", length(base_full)),
  paste0("n_isNA: ", length(isna_full)),
  paste0("final_feature_order_n: ", length(fn_full)),
  paste0("feature_schema_fingerprint: ", fp_full),
  "", "=== base features ===", sort(base_full),
  "", "=== isNA companions ===", sort(isna_full)),
  file.path(PROF, "FULL_feature_list.txt"))
log("[full] feature schema: base=", length(base_full), " isNA=", length(isna_full),
    " fingerprint=", fp_full)

# ---- EFFIS (dedicated, memory-safe) -----------------------------------------
log("[full] EFFIS validation ...")
effis <- run_long_effis(fm, file.path(PROF, "VALIDATION_EFFIS"), log)
if (!is.null(effis)) {
  saveRDS(effis, file.path(PROF, "effis.rds"))
  m <- effis$metrics %||% effis
  df <- tryCatch(as.data.frame(m), error = function(e) NULL)
  if (!is.null(df)) {
    utils::write.csv(df, file.path(PROF, "EFFIS_metrics.csv"), row.names = FALSE)
    log("[full] EFFIS metrics written.")
  }
}

# ---- SUMMARY + manifest -----------------------------------------------------
writeLines(if (length(.warns)) .warns else "(none)", file.path(PROF, "WARNINGS.txt"))
md <- c("# LONG BALANCED — FULL profile (thermal-enhanced; all features incl hotspots)",
  "", sprintf("- Snapshot: %s (tarball SHA256 %s)", SNAPSHOT_ID, TARBALL_SHA256),
  sprintf("- Training: inner-early-stopping selection + full-data refit; nrounds_max=%d early_stop=%d seed=%d nthread=%d.",
          NROUNDS_MAX, EARLY_STOP, SEED_BASE, XGB_NTHREAD),
  sprintf("- Caps: contextual %.2f / spectral %.2f / random %.2f / otsu %.2f.",
          CAP_CONTEXTUAL, CAP_SPECTRAL, CAP_RANDOM, CAP_OTSU),
  sprintf("- Features: %d base + %d _isNA (FULL whitelist; hotspots included).",
          length(base_full), length(isna_full)),
  sprintf("- feature_schema_fingerprint: %s", fp_full),
  "- DEPLOYED: OOF (capped) + FINAL; OOF uses the same capped policy as FINAL.",
  "- Deployed FINAL scored + mapped + EFFIS validated.",
  "", "## EFFIS metrics (validate_fire_maps, metrics_type='all')",
  if (exists("df") && !is.null(df)) paste(utils::capture.output(print(t(df))), collapse = "\n") else "(EFFIS not produced)")
writeLines(md, file.path(PROF, "SUMMARY.md"))

dur <- as.numeric(Sys.time() - t0, units = "mins")
log(sprintf("==== FULL PROFILE DONE (%.2f min) ====", dur))
writeLines(c(sprintf("FULL_OK %s", format(Sys.time())),
             sprintf("B_best_iteration=%s scored_rows=%d fingerprint=%s",
                     B_model$recipe$training$best_iteration %||% NA, nrow(fm), fp_full)),
           file.path(PROF, "FULL_SUCCESS.marker"))
