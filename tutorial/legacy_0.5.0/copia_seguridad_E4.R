# ============================================================
# E4 snapshot, before closing the session
# ============================================================

YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
RESULTS_BASE <- file.path(ROOT, "1_DATA", "Results")

phase_b_root <- file.path(RESULTS_BASE, "_PHASE_B_RESULTS", "2005_balanced")
production_dir <- file.path(RESULTS_BASE, as.character(YEAR), "Min_Min",
                            "SUPERVISED", SCENARIO)

e4_snap <- file.path(phase_b_root, "E4_no_hotspots_corine5x")
if (dir.exists(e4_snap)) {
  cat("E4 snapshot already exists; removing it for a clean overwrite.\n")
  unlink(e4_snap, recursive = TRUE, force = TRUE)
}

cat("Copiando E4:", production_dir, "->", e4_snap, "\n")
dir.create(e4_snap, recursive = TRUE)
R.utils::copyDirectory(production_dir, e4_snap,
                       recursive = TRUE, overwrite = TRUE,
                       private = TRUE, copy.mode = TRUE, copy.date = TRUE)

cat("\nDoes meta confirm E4?:\n")
print(grep("whitelist_override|weights_applied|weights_n|weights_nondefault|n_x_cols|best_iteration",
           readLines(file.path(e4_snap, "07_FINAL_MODEL_V2",
                               "2005_balanced_patch_certified_meta.txt")),
           value = TRUE))

cat("\nSubdirs en _PHASE_B_RESULTS/2005_balanced:\n")
print(list.dirs(phase_b_root, recursive = FALSE, full.names = FALSE))
