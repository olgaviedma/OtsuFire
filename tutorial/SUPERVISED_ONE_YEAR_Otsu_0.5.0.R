# ============================================================
# Baseline 2005/balanced — OtsuFire 0.5.0 — desde cero
# ============================================================

suppressPackageStartupMessages({
  library(sf)
})

?build_burned_mapping_config()

# Cargar paquete
PKG_ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"
if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

ver <- utils::packageVersion("OtsuFire")
cat("OtsuFire version:", as.character(ver), "\n")
stopifnot(ver >= "0.5.0")

# Paths
YEAR <- 2005L
SCENARIO <- "balanced"
ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING"
DATA_BASE <- file.path(ROOT, "1_DATA")
RESULTS_BASE <- file.path(DATA_BASE, "Results")
COMPOSITE_BASE <- file.path(DATA_BASE, "Imagery", "Composites_90m")
RESULT_NAME <- "Min_Min"

INTERNAL_DECISIONS <- file.path(
  RESULTS_BASE, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "05_DECISIONS", "internal_decisions.gpkg"
)
CHANGE_INDEX <- file.path(COMPOSITE_BASE, RESULT_NAME,
                          sprintf("MinMin_%d_mosaic_res90m.tif", YEAR))
HOTSPOTS <- file.path(DATA_BASE, "Hotspots",
                      sprintf("hotspots_iberia_%d.geojson", YEAR))

stopifnot(file.exists(INTERNAL_DECISIONS),
          file.exists(CHANGE_INDEX),
          file.exists(HOTSPOTS))

# Construir config
cfg <- OtsuFire::build_supervised_burned_config(
  scenario = SCENARIO,
  internal_decisions = INTERNAL_DECISIONS,
  change_index = CHANGE_INDEX,
  hotspots = HOTSPOTS,
  target_year = YEAR,
  output_dir = RESULTS_BASE,
  run_name = RESULT_NAME,
  options = list(
    data_base = DATA_BASE,
    composite_base = COMPOSITE_BASE,
    result_name = RESULT_NAME,
    # negative-pool policy is always all_sources now (implicit; only mode)
    legacy_otsu_mode = "burnable_only",
    legacy_otsu_threshold = 0,
    legacy_reference_otsu_threshold = 100,
    legacy_sample_n = 2000L,
    legacy_reuse_existing = TRUE,
    legacy_write_output = TRUE,
    unb_verbose = TRUE
  )
)

# Tool paths Anaconda3
cfg$tool_paths$python_exe              <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/python.exe"
cfg$tool_paths$gdal_polygonize_script  <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Scripts/gdal_polygonize.py"
cfg$tool_paths$gdalwarp_path           <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/gdalwarp.exe"
cfg$tool_paths$ogr2ogr_exe             <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/ogr2ogr.exe"

stopifnot(all(sapply(cfg$tool_paths, file.exists)))
cat("✓ Tool paths OK.\n\n")

# Lanzar pipeline (Baseline = todos defaults, sin override, sin pesos)
cat("=== LANZANDO BASELINE 0.5.0 (esperar ~60-70 min) ===\n")
t0 <- Sys.time()

result <- OtsuFire::run_oneyear_supervised_pipeline(
  config = cfg,
  run_consistency = FALSE,
  overwrite = TRUE,
  contextual_exclusion_to_burned_ratio = 0.25,
  spectral_hard_negative_to_burned_ratio = 1.0,
  random_to_burned_ratio = 1.0,
  otsu_unburned_to_burned_ratio = 1.0,
  feature_whitelist_override = NULL,    # Baseline = canonical 51
  feature_weights = NULL,                # Baseline = todo peso 1.0
  reuse_upstream = FALSE
)

elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nPipeline finished in %.1f min.\n", elapsed_min))

# Verificaciones
final_map_path <- result$final_map_gpkg
fm <- sf::read_sf(final_map_path, layer = "final_map_full")

drops <- fm$p_burned_model[fm$class_final == "drop"]
keeps <- fm$p_burned_model[fm$class_final == "keep"]

cat("\n=== DROPS ===\n"); print(summary(drops))
cat(sprintf("Mediana: %.4f\n", median(drops, na.rm = TRUE)))

cat("\n=== KEEPS ===\n"); print(summary(keeps))
cat(sprintf("Mediana: %.4f\n", keeps_median <- median(keeps, na.rm = TRUE)))

# Verificar que NO se aplicó override ni weights
meta_path <- file.path(dirname(dirname(final_map_path)),
                       "07_FINAL_MODEL_V2",
                       "2005_balanced_patch_certified_meta.txt")
meta_lines <- readLines(meta_path)
override_line <- grep("feature_whitelist_override_applied", meta_lines, value = TRUE)
weights_line <- grep("feature_weights_applied", meta_lines, value = TRUE)
n_x_line <- grep("^n_x_cols:", meta_lines, value = TRUE)
cat("\n=== Verificación API 0.5.0 ===\n")
cat(override_line, "\n")
cat(weights_line, "\n")
cat(n_x_line, "\n")
# Esperado:
# feature_whitelist_override_applied: FALSE
# feature_weights_applied: FALSE
# n_x_cols: 50 o 51 (según dónde quede el conteo en 0.5.0)

# Veredicto
drops_median <- median(drops, na.rm = TRUE)
cat("\n========================================\n")
if (drops_median < 0.1 && keeps_median > 0.9) {
  cat("✅ BASELINE 0.5.0 PASSED.\n")
  cat(sprintf("   Drops mediana: %.4f | Keeps mediana: %.4f\n",
              drops_median, keeps_median))
} else {
  cat("❌ BASELINE 0.5.0 FAILED.\n")
  cat(sprintf("   Drops: %.4f, Keeps: %.4f\n", drops_median, keeps_median))
}
cat("========================================\n")