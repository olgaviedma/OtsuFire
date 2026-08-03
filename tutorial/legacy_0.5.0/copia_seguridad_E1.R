# ============================================================
# E1 snapshot, before moving on to more experiments
# ============================================================

YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
RESULTS_BASE <- file.path(ROOT, "1_DATA", "Results")

phase_b_root <- file.path(RESULTS_BASE, "_PHASE_B_RESULTS", "2005_balanced")
production_dir <- file.path(RESULTS_BASE, as.character(YEAR), "Min_Min",
                            "SUPERVISED", SCENARIO)

# Snapshot E1
e1_snap <- file.path(phase_b_root, "E1_no_hotspots")
if (dir.exists(e1_snap)) {
  cat("E1 snapshot already exists; removing it for a clean overwrite.\n")
  unlink(e1_snap, recursive = TRUE, force = TRUE)
}

cat("Copiando E1:", production_dir, "->", e1_snap, "\n")
dir.create(e1_snap, recursive = TRUE)
R.utils::copyDirectory(production_dir, e1_snap,
                       recursive = TRUE, overwrite = TRUE,
                       private = TRUE, copy.mode = TRUE, copy.date = TRUE)

# Verificar
cat("\nSubdirs en E1 snapshot:\n")
print(list.dirs(e1_snap, recursive = FALSE, full.names = FALSE))

cat("\nDoes meta confirm E1?:\n")
meta_path <- file.path(e1_snap, "07_FINAL_MODEL_V2",
                       "2005_balanced_patch_certified_meta.txt")
print(grep("whitelist_override|n_x_cols|best_iteration|n_training",
           readLines(meta_path), value = TRUE))



# ============================================================
# Snapshot E1 antes de lanzar E2
# ============================================================

YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
RESULTS_BASE <- file.path(ROOT, "1_DATA", "Results")

phase_b_root <- file.path(RESULTS_BASE, "_PHASE_B_RESULTS", "2005_balanced")
production_dir <- file.path(RESULTS_BASE, as.character(YEAR), "Min_Min",
                            "SUPERVISED", SCENARIO)

e1_snap <- file.path(phase_b_root, "E1_no_hotspots")
if (dir.exists(e1_snap)) {
  cat("E1 snapshot already exists; removing it for a clean overwrite.\n")
  unlink(e1_snap, recursive = TRUE, force = TRUE)
}

cat("Copiando E1:", production_dir, "->", e1_snap, "\n")
dir.create(e1_snap, recursive = TRUE)
R.utils::copyDirectory(production_dir, e1_snap,
                       recursive = TRUE, overwrite = TRUE,
                       private = TRUE, copy.mode = TRUE, copy.date = TRUE)

cat("\nSubdirs en E1 snapshot:\n")
print(list.dirs(e1_snap, recursive = FALSE, full.names = FALSE))

cat("\nDoes meta confirm E1?:\n")
print(grep("whitelist_override|n_x_cols|best_iteration",
           readLines(file.path(e1_snap, "07_FINAL_MODEL_V2",
                               "2005_balanced_patch_certified_meta.txt")),
           value = TRUE))

