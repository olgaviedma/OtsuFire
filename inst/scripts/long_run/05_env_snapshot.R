# =============================================================================
# LONG <YEAR> BALANCED — ENV SNAPSHOT (recorded BEFORE any heavy step)
# CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
# Writes <ORCH>/env_snapshot.txt: disk/mem/tempdir, fixed nthread, R version,
# OtsuFire pin + find.package, snapshot id + tarball SHA256, and SHA256 of every
# resolved input. NO heavy compute.
# =============================================================================
# EDIT: point this at your run's _ORCHESTRATION/00_common.R (working copy).
source("<PATH_TO_ORCHESTRATION>/00_common.R")
ORCH <- file.path(LONG_ROOT, "_ORCHESTRATION")
OUT  <- file.path(ORCH, "env_snapshot.txt")
.assert_version_pin()

sha256 <- function(p) {
  if (is.null(p) || is.na(p) || !file.exists(p)) return("MISSING")
  tryCatch(as.character(tools::sha256sum(p)), error = function(e) {
    tryCatch(digest::digest(p, algo = "sha256", file = TRUE), error = function(e2) "ERR")
  })
}

free_gb <- tryCatch({
  fi <- suppressWarnings(system2("powershell", c("-NoProfile","-Command",
    "(Get-PSDrive C).Free"), stdout = TRUE)); round(as.numeric(fi[1]) / 1e9, 1)
}, error = function(e) NA)
mem_gb <- tryCatch({
  mi <- suppressWarnings(system2("powershell", c("-NoProfile","-Command",
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE))
  round(as.numeric(mi[1]) / 1e6, 1)
}, error = function(e) NA)

inputs <- c(
  internal_decisions = internal_decisions,
  change_index = change_index,
  delayed_change_index = delayed_change_index,
  hotspots = hotspots,
  topo = topo_path,
  corine_raster = corine_raster,
  corine_strata_raster = strata_raster,
  corine_strata_lut = strata_lut,
  burnable_mask_run_3035 = burnable_mask_run,
  burnable_mask_wgs84 = burnable_mask_wgs84,
  peninsula_mask = validation_mask_shapefile,
  effis_reference = reference_burned_map)

lines <- c(
  "==== LONG BALANCED — ENV SNAPSHOT ====",
  paste0("timestamp: ", format(Sys.time())),
  paste0("R version: ", R.version.string),
  paste0("OtsuFire packageVersion: ", as.character(packageVersion("OtsuFire"))),
  paste0("find.package(OtsuFire): ", normalizePath(find.package("OtsuFire"), winslash = "/")),
  paste0("isolated Rlib: ", RLIB),
  paste0("snapshot id: ", SNAPSHOT_ID),
  paste0("tarball: ", TARBALL_PATH),
  paste0("tarball SHA256 (declared): ", TARBALL_SHA256),
  paste0("tarball SHA256 (computed): ", sha256(TARBALL_PATH)),
  paste0("xgboost nthread (FIXED, both profiles): ", XGB_NTHREAD),
  paste0("OMP_NUM_THREADS: ", Sys.getenv("OMP_NUM_THREADS")),
  paste0("free disk C: GB: ", free_gb),
  paste0("free physical mem GB: ", mem_gb),
  paste0("tempdir(): ", tempdir()),
  paste0("nrounds_max: ", NROUNDS_MAX, " ; early_stop: ", EARLY_STOP, " ; seed: ", SEED_BASE),
  paste0("caps: contextual=", CAP_CONTEXTUAL, " spectral=", CAP_SPECTRAL,
         " random=", CAP_RANDOM, " otsu=", CAP_OTSU),
  "", "==== RESOLVED INPUT SHA256 ====")
for (nm in names(inputs)) {
  lines <- c(lines, sprintf("%-26s %s | %s", nm, sha256(inputs[[nm]]), inputs[[nm]]))
}
writeLines(lines, OUT)
cat(paste(lines, collapse = "\n"), "\n")
cat("\n[env] snapshot written to ", OUT, "\n")
