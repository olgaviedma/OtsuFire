# =============================================================================
# LONG <YEAR> BALANCED — SHARED UPSTREAM (built ONCE; consumed by BOTH profiles)
# CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
# candidates -> internal labels -> positive+negative pools -> spatial folds ->
# extract_supervised_features (FULL feature set incl hotspots).
# Both profiles reuse THIS feature GPKG; they only select different feature
# subsets via feature_whitelist_override. Do NOT regenerate candidates / labels /
# pools / folds for the no-hotspot profile.
# =============================================================================
# EDIT: point this at your run's _ORCHESTRATION/00_common.R (working copy).
source("<PATH_TO_ORCHESTRATION>/00_common.R")

SHARED <- file.path(LONG_ROOT, "SHARED")
dir.create(SHARED, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(SHARED, "shared_upstream.log")
log <- make_logger(LOG)
t0 <- Sys.time()

log("==== SHARED UPSTREAM START ====")
.assert_version_pin()

# Build the SHARED/upstream cfg. Pools/folds/features are hotspot-INCLUSIVE: use
# the FULL cfg (no whitelist override, hotspots ON). Training is always
# inner-early-stopping selection + full-data refit, and OOF uses the same capped
# negative-sampling policy as the FINAL model. This is the experimental control
# reused by both profiles.
cfg_shared <- .mk_long_cfg(feature_whitelist_override = NULL)
stopifnot(all(vapply(cfg_shared$tool_paths, file.exists, logical(1))))
saveRDS(cfg_shared, file.path(SHARED, "cfg_shared_upstream.rds"))
log("[shared] cfg built; caps=",
    paste(names(unlist(cfg_shared$train_control$caps)), unlist(cfg_shared$train_control$caps),
          sep = "=", collapse = " "))

.warns <- character(0)
wh <- function(w) { .warns[[length(.warns) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }

# ---- POOLS ------------------------------------------------------------------
log("[shared] STAGE pools ...")
pools <- withCallingHandlers(
  build_supervised_training_pools(config = cfg_shared, write_outputs = TRUE, overwrite = TRUE),
  warning = wh)
saveRDS(pools, file.path(SHARED, "pools.rds"))
n_burned <- sum(as.character(pools$train_labeled$class) == "burned", na.rm = TRUE)
log("[shared] pools: train_labeled=", nrow(pools$train_labeled), " n_burned=", n_burned,
    " scoring_pool=", nrow(pools$scoring_pool))
stopifnot(n_burned > 0)

# ---- FOLDS ------------------------------------------------------------------
log("[shared] STAGE folds ...")
folds <- withCallingHandlers(
  make_spatial_folds(train_labelled = pools$train_labeled, config = cfg_shared,
    split_unit = "fire", out_dir = NULL, write_outputs = TRUE, overwrite = TRUE),
  warning = wh)
saveRDS(folds, file.path(SHARED, "folds.rds"))
log("[shared] folds: train_with_folds=", nrow(folds$train_with_folds),
    " gpkg=", folds$train_with_folds_gpkg %||% "(in-memory)")

# ---- FEATURES (FULL feature set incl hotspots) ------------------------------
log("[shared] STAGE features (use_hotspots=TRUE; FULL feature set) ...")
feats <- withCallingHandlers(
  extract_supervised_features(train_with_folds = folds$train_with_folds,
    scoring_pool = pools$scoring_pool, config = cfg_shared,
    use_hotspots = TRUE, out_dir = NULL, write_outputs = TRUE),
  warning = wh)
saveRDS(feats, file.path(SHARED, "feats.rds"))
log("[shared] features: train=", nrow(feats$train_features), " x ", ncol(feats$train_features),
    " ; scoring=", nrow(feats$scoring_features), " x ", ncol(feats$scoring_features))

# Confirm the SHARED feature frame carries the hotspot columns (so NO-HOTSPOT is
# a true subset selection, not a re-extraction).
hs_present <- intersect(HOTSPOT_LINEAGE_BASE, names(feats$train_features))
log("[shared] hotspot-lineage cols present in feature frame: ", length(hs_present),
    " of ", length(HOTSPOT_LINEAGE_BASE))

# ---- feature lineage table (shared provenance) ------------------------------
lineage <- build_feature_lineage()
utils::write.csv(lineage, file.path(SHARED, "feature_lineage.csv"), row.names = FALSE)
log("[shared] feature lineage table written (", nrow(lineage), " base features).")

# ---- SHARED manifest --------------------------------------------------------
manifest <- data.frame(
  key = c("n_burned", "train_features_rows", "train_features_cols",
          "scoring_features_rows", "scoring_features_cols",
          "fold_cols", "seed_base", "cap_contextual", "cap_spectral",
          "cap_random", "cap_otsu", "nrounds_max", "early_stop", "xgb_nthread",
          "hotspot_cols_present"),
  value = c(n_burned, nrow(feats$train_features), ncol(feats$train_features),
            nrow(feats$scoring_features), ncol(feats$scoring_features),
            paste(FOLD_COLS, collapse = "+"), SEED_BASE, CAP_CONTEXTUAL, CAP_SPECTRAL,
            CAP_RANDOM, CAP_OTSU, NROUNDS_MAX, EARLY_STOP, XGB_NTHREAD,
            length(hs_present)),
  stringsAsFactors = FALSE)
utils::write.csv(manifest, file.path(SHARED, "shared_manifest.csv"), row.names = FALSE)

writeLines(if (length(.warns)) .warns else "(none)", file.path(SHARED, "shared_WARNINGS.txt"))
dur <- as.numeric(Sys.time() - t0, units = "mins")
log(sprintf("==== SHARED UPSTREAM DONE (%.2f min) ====", dur))

# SUCCESS marker (the master gates on this).
writeLines(c(sprintf("SHARED_UPSTREAM_OK %s", format(Sys.time())),
             sprintf("n_burned=%d train_features=%dx%d scoring=%dx%d",
                     n_burned, nrow(feats$train_features), ncol(feats$train_features),
                     nrow(feats$scoring_features), ncol(feats$scoring_features))),
           file.path(SHARED, "SHARED_SUCCESS.marker"))
