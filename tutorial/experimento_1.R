# ============================================================
# E1 — Hotspots fuera (las 13 features G6)
# Reusa pools/folds/features de Baseline en Min_Min
# Sobreescribe 04_MATRIX/ → 09_FINAL_MAP/ con resultados E1
# ============================================================

suppressPackageStartupMessages({
  library(sf)
})

# ---- Setup ----
PKG_ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"
if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

ver <- utils::packageVersion("OtsuFire")
cat("OtsuFire version:", as.character(ver), "\n")
stopifnot(ver >= "0.5.0")

# ---- Definir el override ----
canonical <- OtsuFire:::.supervised_feature_cols

HOTSPOT_BLOCK <- c(
  "hotspot_available",
  "hs_in_poly", "hs_in_buffer", "hs_used_n", "hs_min_dist_m",
  "hs_frp_sum", "hs_frp_max", "hs_conf_mean", "hs_hiConf_n",
  "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support"
)

stopifnot(all(HOTSPOT_BLOCK %in% canonical))

E1_whitelist <- setdiff(canonical, HOTSPOT_BLOCK)
cat("E1 whitelist size:", length(E1_whitelist), "(expected 38)\n")

# ---- Reconstruir config (mismo que Baseline) ----
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

cfg <- OtsuFire::build_supervised_burned_config(
  run_label = SCENARIO,
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

cfg$tool_paths$python_exe              <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/python.exe"
cfg$tool_paths$gdal_polygonize_script  <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Scripts/gdal_polygonize.py"
cfg$tool_paths$gdalwarp_path           <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/gdalwarp.exe"
cfg$tool_paths$ogr2ogr_exe             <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/ogr2ogr.exe"

stopifnot(all(sapply(cfg$tool_paths, file.exists)))

# ---- Lanzar E1 con reuse_upstream ----
cat("\n=== LANZANDO E1 — Hotspots fuera (~10-20 min) ===\n")
t0 <- Sys.time()

result_e1 <- OtsuFire::run_oneyear_supervised_pipeline(
  config = cfg,
  run_consistency = FALSE,
  overwrite = TRUE,
  contextual_exclusion_to_burned_ratio = 0.25,
  spectral_hard_negative_to_burned_ratio = 1.0,
  random_to_burned_ratio = 1.0,
  otsu_unburned_to_burned_ratio = 1.0,
  feature_whitelist_override = E1_whitelist,
  feature_weights = NULL,
  reuse_upstream = TRUE
)

elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nE1 finished in %.1f min.\n", elapsed_min))

# ---- Verificación rápida ----
final_map_path <- result_e1$final_map_gpkg
cat("\nfinal_map_path:", final_map_path, "\n")
fm <- sf::read_sf(final_map_path, layer = "final_map_full")
cat("Rows:", nrow(fm), "\n\n")

drops <- fm$p_burned_model[fm$class_final == "drop"]
keeps <- fm$p_burned_model[fm$class_final == "keep"]
review <- fm$p_burned_model[fm$class_final == "review"]

cat("=== DROPS ===\n"); print(summary(drops))
cat(sprintf("Mediana: %.4f\n\n", median(drops, na.rm = TRUE)))

cat("=== KEEPS ===\n"); print(summary(keeps))
cat(sprintf("Mediana: %.4f\n\n", median(keeps, na.rm = TRUE)))

cat("=== REVIEW ===\n"); print(summary(review))
cat(sprintf("Mediana: %.4f\n", median(review, na.rm = TRUE)))

# ¿Cuántos reviews rescata E1?
n_rescued_e1 <- sum(review > 0.5, na.rm = TRUE)
cat(sprintf("\nReviews rescatados por E1 (p>0.5): %d / %d\n",
            n_rescued_e1, length(review)))
cat("(Comparar con Baseline: 30 / 1114)\n")

# Verificar override aplicado
meta_e1_path <- file.path(dirname(dirname(final_map_path)),
                          "07_FINAL_MODEL_V2",
                          "2005_balanced_patch_certified_meta.txt")
meta_e1 <- readLines(meta_e1_path)
cat("\n=== Verificación E1 meta ===\n")
print(grep("whitelist_override|n_x_cols|n_feature|best_iteration|n_training",
           meta_e1, value = TRUE))









# Cargar OOF metrics de E1
oof_dir_e1 <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results",
  "2005/Min_Min/SUPERVISED/balanced/05_OOF"
)

cat("=== E1 OOF metrics summary ===\n\n")
cat(readLines(file.path(oof_dir_e1,
                        "2005_balanced_patch_oof_metrics_summary.txt")),
    sep = "\n")

# Y feature importance de E1
fi_path_e1 <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results",
  "2005/Min_Min/SUPERVISED/balanced/07_FINAL_MODEL_V2",
  "2005_balanced_patch_certified_feature_importance.csv"
)
fi_e1 <- read.csv(fi_path_e1, stringsAsFactors = FALSE)

cat("\n\n=== E1 Top 15 features ===\n")
print(head(fi_e1[, c("Feature", "Gain")], 15))

cat("\n=== E1 ¿alguna hotspot aparece? (no debería) ===\n")
hs_in_e1 <- fi_e1[grepl("^(hs_|hotspot_)", fi_e1$Feature), ]
if (nrow(hs_in_e1) > 0) {
  cat("ERROR: Hotspot features in E1!\n")
  print(hs_in_e1)
} else {
  cat("OK: ninguna hotspot en E1.\n")
}

cat("\n=== E1 número total de features con Gain > 0 ===\n")
cat(nrow(fi_e1), "/ 38 disponibles\n")




# Cargar E1 final_map_full
fm_e1 <- sf::read_sf(final_map_path, layer = "final_map_full")

# Cargar Baseline final_map_full (del snapshot)
baseline_fm_path <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results",
  "_PHASE_B_RESULTS/2005_balanced/Baseline/09_FINAL_MAP",
  "2005_balanced_patch_certified_final_map.gpkg"
)
fm_baseline <- sf::read_sf(baseline_fm_path, layer = "final_map_full")

# Subset reviews
reviews_baseline <- fm_baseline[fm_baseline$class_final == "review", ]
reviews_e1 <- fm_e1[fm_e1$class_final == "review", ]

# Identificar rescatados (p > 0.5) en cada uno
rescued_b_ids <- reviews_baseline$poly_id[reviews_baseline$p_burned_model > 0.5]
rescued_e1_ids <- reviews_e1$poly_id[reviews_e1$p_burned_model > 0.5]

cat("Rescatados Baseline:", length(rescued_b_ids), "\n")
cat("Rescatados E1:", length(rescued_e1_ids), "\n")

# Intersección
common <- intersect(rescued_b_ids, rescued_e1_ids)
only_baseline <- setdiff(rescued_b_ids, rescued_e1_ids)
only_e1 <- setdiff(rescued_e1_ids, rescued_b_ids)

cat("\n=== Comparación de rescatados ===\n")
cat("En AMBOS (sólidos):", length(common), "\n")
cat("Solo Baseline (perdidos por E1):", length(only_baseline), "\n")
cat("Solo E1 (nuevos sin hotspot):", length(only_e1), "\n")

# Para los "nuevos solo E1": tienen hotspot dentro?
new_e1 <- reviews_baseline[reviews_baseline$poly_id %in% only_e1, ]
cat("\n=== De los", length(only_e1), "nuevos rescatados solo E1 ===\n")
cat("Tienen hotspot dentro (hs_in_poly > 0)?:\n")
print(table(new_e1$hs_in_poly > 0, useNA = "always"))

cat("\nDistribución de RBR med:\n")
print(summary(new_e1$rbr_med))

cat("\nDistribución de área:\n")
print(summary(new_e1$area_ha))





# Ver perfil de los 14 perdidos
lost_baseline <- reviews_baseline[reviews_baseline$poly_id %in% only_baseline, ]
cat("=== 14 burned perdidos por E1 ===\n")
cat("\np_burned_model en Baseline:\n")
print(summary(lost_baseline$p_burned_model))

cat("\nhs_in_poly:\n")
print(table(lost_baseline$hs_in_poly, useNA = "always"))

cat("\nrbr_med:\n")
print(summary(lost_baseline$rbr_med))

cat("\narea_ha:\n")
print(summary(lost_baseline$area_ha))

# Y su p_burned_model en E1
lost_e1 <- reviews_e1[reviews_e1$poly_id %in% only_baseline, ]
cat("\np_burned_model en E1 (los mismos polígonos):\n")
print(summary(lost_e1$p_burned_model))

# Y los 16 sólidos: ¿qué probabilidad les da E1?
common_e1 <- reviews_e1[reviews_e1$poly_id %in% common, ]
cat("\n=== 16 sólidos (en ambos) ===\n")
cat("p_burned_model en E1:\n")
print(summary(common_e1$p_burned_model))

cat("\nhs_in_poly de los 16 sólidos (de Baseline):\n")
common_baseline <- reviews_baseline[reviews_baseline$poly_id %in% common, ]
print(summary(common_baseline$hs_in_poly))





