# ============================================================
# Snapshot Baseline antes de lanzar E1
# ============================================================

YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
RESULTS_BASE <- file.path(ROOT, "1_DATA", "Results")

phase_b_root <- file.path(RESULTS_BASE, "_PHASE_B_RESULTS", "2005_balanced")
production_dir <- file.path(RESULTS_BASE, as.character(YEAR), "Min_Min",
                            "SUPERVISED", SCENARIO)

# Crear directorio raíz si no existe
if (!dir.exists(phase_b_root)) dir.create(phase_b_root, recursive = TRUE)

# Snapshot Baseline
baseline_snap <- file.path(phase_b_root, "Baseline")
if (dir.exists(baseline_snap)) {
  cat("Baseline snapshot ya existe, lo borramos para overwrite limpio.\n")
  unlink(baseline_snap, recursive = TRUE, force = TRUE)
}

cat("Copiando Baseline:", production_dir, "->", baseline_snap, "\n")
dir.create(baseline_snap, recursive = TRUE)
R.utils::copyDirectory(production_dir, baseline_snap,
                       recursive = TRUE, overwrite = TRUE,
                       private = TRUE, copy.mode = TRUE, copy.date = TRUE)

# Verificar
cat("\nSubdirs en Baseline snapshot:\n")
print(list.dirs(baseline_snap, recursive = FALSE, full.names = FALSE))
cat("\nMeta confirma Baseline?:\n")
meta_path <- file.path(baseline_snap, "07_FINAL_MODEL_V2",
                       "2005_balanced_patch_certified_meta.txt")
meta <- readLines(meta_path)
print(grep("ratio|whitelist_override|weights_applied|n_x_cols",
           meta, value = TRUE))