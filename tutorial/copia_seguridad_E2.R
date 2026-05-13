# ============================================================
# Snapshot E2 antes de lanzar E4
# ============================================================

YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
RESULTS_BASE <- file.path(ROOT, "1_DATA", "Results")

phase_b_root <- file.path(RESULTS_BASE, "_PHASE_B_RESULTS", "2005_balanced")
production_dir <- file.path(RESULTS_BASE, as.character(YEAR), "Min_Min",
                            "SUPERVISED", SCENARIO)

e2_snap <- file.path(phase_b_root, "E2_no_hotspots_persist5x")
if (dir.exists(e2_snap)) {
  cat("E2 snapshot ya existe, lo borramos para overwrite limpio.\n")
  unlink(e2_snap, recursive = TRUE, force = TRUE)
}

cat("Copiando E2:", production_dir, "->", e2_snap, "\n")
dir.create(e2_snap, recursive = TRUE)
R.utils::copyDirectory(production_dir, e2_snap,
                       recursive = TRUE, overwrite = TRUE,
                       private = TRUE, copy.mode = TRUE, copy.date = TRUE)

cat("\n¿Meta confirma E2?:\n")
print(grep("whitelist_override|weights_applied|weights_nondefault|n_x_cols|best_iteration",
           readLines(file.path(e2_snap, "07_FINAL_MODEL_V2",
                               "2005_balanced_patch_certified_meta.txt")),
           value = TRUE))