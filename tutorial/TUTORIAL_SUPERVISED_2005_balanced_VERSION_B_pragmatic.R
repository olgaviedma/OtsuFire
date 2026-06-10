# =============================================================================
# TUTORIAL_SUPERVISED_2005_balanced_VERSION_B_pragmatic.R
# =============================================================================
#
# OBJETIVO
#   Reproducir paso a paso el pipeline supervisado de OtsuFire 0.5.0 para el
#   ano 2005 escenario balanced en configuracion CANONICAL BASELINE,
#   llamando a las funciones internas (OtsuFire:::funcion_interna) en lugar
#   de a la unica funcion publica `run_oneyear_supervised_pipeline()`.
#
# CONFIGURACION CANONICAL BASELINE (confirmada por Natalia, ver HANDOFF
# §N+23):
#   - feature_weights              = NULL  (peso uniforme = 1.0)
#   - feature_whitelist_override   = NULL  (canon natural del paquete: 51)
#   - contextual_exclusion_to_burned_ratio   = 0.25
#   - spectral_hard_negative_to_burned_ratio = 1.0
#   - random_to_burned_ratio                 = 1.0
#   - otsu_unburned_to_burned_ratio          = 1.0
#   - reuse_upstream                         = FALSE
#   En 2005 (post-MODIS) el modelo usara las 51 features nominales
#   (incluyendo las 13 hs_*).
#
# ESTRATEGIA DE ESTE TUTORIAL (VERSION B - pragmatica)
#   1) Construir cfg con la funcion publica build_supervised_burned_config().
#      Esto rellena automaticamente los ~30 parametros internos de Otsu
#      legacy y unburned, evitando reproducir el dispatcher a mano.
#   2) Asignar tool_paths a cfg.
#   3) Llamar manualmente a las funciones internas (OtsuFire:::xxx) en el
#      orden exacto que sigue internal-sup-orchestrator.R, deteniendonos
#      tras cada STEP para inspeccionar los outputs.
#
# DIRECTORIO DE SALIDA
#   Para no machacar el run del agente nocturno, este tutorial escribe en:
#     Results/2005/Min_Min/SUPERVISED_TUTORIAL/balanced/
#   en vez de en Results/2005/Min_Min/SUPERVISED/balanced/.
#
# COMO USAR
#   - Restart R limpio (Ctrl+Shift+F10).
#   - Ejecutar bloque a bloque (seleccionar y Ctrl+Enter).
#   - Inspeccionar las variables que cada bloque deja en memoria.
#   - Comparar al final con el output del wrapper.
#
# REFERENCIAS EN EL CODIGO DEL PAQUETE
#   - Wrapper publico:           R/supervised-run.R
#   - Orchestrator (corazon):    R/internal-sup-orchestrator.R
#                                run_supervised_pipeline() linea 719
#                                STEP A1-A6, B1-B3, C0-C4
#   - Funciones internas:        R/internal-sup-*.R y R/supervised-*.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
})

# =============================================================================
# BLOQUE 0 - SETUP: CARGA DEL PAQUETE Y PATHS
# =============================================================================
#
# QUE HACE
#   Carga el paquete OtsuFire desde el directorio de fuentes (sin instalarlo)
#   con pkgload::load_all(). Define los paths globales del proyecto.
#
# POR QUE
#   load_all() es la forma estandar de trabajar con un paquete en desarrollo:
#   ve los exports del NAMESPACE pero tambien permite acceso a funciones
#   internas con OtsuFire:::xxx.
#
# OUTPUTS EN MEMORIA
#   - PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, etc. (paths)
#   - YEAR=2005L, SCENARIO="balanced"
#   - INTERNAL_DECISIONS, CHANGE_INDEX, HOTSPOTS (paths a inputs criticos)
#
# =============================================================================
# =============================================================================
# BLOQUE 0 - SETUP: CARGA DEL PAQUETE Y PATHS
# =============================================================================
#
# QUE HACE
#   Carga el paquete OtsuFire desde el directorio de fuentes (sin instalarlo)
#   con pkgload::load_all(). Define los paths globales del proyecto y los
#   4 paths de inputs criticos del experimento (Y/E = year/scenario):
#     - internal_decisions.gpkg          (deterministic stage)
#     - mosaico summer (RBR + DOY)       (composite ano)
#     - mosaico autumn-winter            (composite delayed)
#     - hotspots geojson                 (post-MODIS, NULL en pre-MODIS)
#
# POR QUE
#   load_all() es la forma estandar de trabajar con un paquete en desarrollo:
#   ve los exports del NAMESPACE pero tambien permite acceso a funciones
#   internas con OtsuFire:::xxx.
#
# OUTPUTS EN MEMORIA
#   - PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, TUTORIAL_RESULT_DIR (paths)
#   - YEAR=2005L, SCENARIO="balanced", RESULT_NAME="Min_Min"
#   - INTERNAL_DECISIONS, CHANGE_INDEX, RBR_AUTUMN, HOTSPOTS
#     (paths a los 4 inputs criticos del combo)
#
# =============================================================================

# =============================================================================
# BLOQUE 0 - SETUP: PAQUETE Y LOS 14 INPUTS DEL PIPELINE SUPERVISED
# =============================================================================
#
# QUE HACE
#   Carga el paquete OtsuFire desde fuentes (sin instalarlo) con
#   pkgload::load_all() y declara los 14 paths de inputs criticos del
#   modulo supervised, organizados por categoria.
#
# POR QUE 14 INPUTS Y NO 4
# ------------------------
# El config builder build_supervised_burned_config() solo acepta 4 paths
# de inputs (internal_decisions, change_index, hotspots, y el cosmetico
# delayed_change_index). Pero el pipeline real lee 14 ficheros distintos:
#
#   - 4 los pasamos al config builder (year/scenario-dependent + delayed)
#   - 10 los construye el orchestrator por convencion desde data_base,
#     composite_base, target_year y corine_year
#
# Esta inconsistencia esta documentada en HANDOFF §N+25 como deuda tecnica
# del paquete. El refactor post-paper la cableara: todos los 14 entraran
# por el config builder.
#
# Para el tutorial, declaramos LOS 14 explicitamente al inicio. Asi:
#   1. Vemos de un vistazo todos los datos que el pipeline va a tocar.
#   2. Validamos al inicio (en <1 segundo) en lugar de fallar 5-10 min
#      despues si algun fichero falta.
#   3. Cuando el refactor cablee los 10 implicitos al cfg, este tutorial
#      seguira funcionando: solo cambiara la llamada al config builder
#      en BLOQUE 1.
#
#
# OUTPUTS EN MEMORIA
# ------------------
#   Paths de proyecto:
#     PKG_ROOT, DATA_BASE, COMPOSITE, RESULTS, TUTORIAL_RESULT_DIR
#
#   Combo objetivo:
#     YEAR=2005L, SCENARIO="balanced", RESULT_NAME="Min_Min", CORINE_YEAR
#
#   Inputs (14):
#     A) Year/scenario-dependent (4):
#        INTERNAL_DECISIONS, CHANGE_INDEX, RBR_AUTUMN, HOTSPOTS
#     B) Validacion externa year-dependent (2):
#        EFFIS_CA_TIF, EFFIS_CA_SHP
#     C) Corine quinquenal (4):
#        CORINE_RASTER, CORINE_STRATA, CORINE_LUT, BURNEABLE_MASK
#     D) Estaticos (4):
#        PENINSULA_SHP, TOPO, MASK_TIF, MASK_SHP
#
#
# REFERENCIA EN EL PAQUETE
# ------------------------
#   internal-sup-orchestrator.R lineas 906-934 - construccion de paths
#   internal-sup-orchestrator.R lineas 939-950 - validacion stopifnot
#   internal-sup-unburned-deterministic.R linea 170 - burneable mask
#   internal-sup-unburned-legacy.R linea 655 - burneable mask
#
# =============================================================================

cat("\n========== BLOQUE 0: SETUP ==========\n")

# --- Carga del paquete -------------------------------------------------------
PKG_ROOT  <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"

if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

cat("OtsuFire version:", as.character(utils::packageVersion("OtsuFire")), "\n")
stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

# --- Paths del proyecto ------------------------------------------------------
DATA_BASE <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA"
COMPOSITE <- file.path(DATA_BASE, "Imagery", "Composites_90m")
RESULTS   <- file.path(DATA_BASE, "Results")

# --- Combo objetivo ----------------------------------------------------------
YEAR        <- 2005L
SCENARIO    <- "balanced"
RESULT_NAME <- "Min_Min"

# Resolucion del ano Corine quinquenal usando la funcion interna del paquete.
# Mapeo (segun el codigo de OtsuFire):
#   1985-2002 -> CLC2000 ; 2003-2008 -> CLC2006 ; 2009-2014 -> CLC2012 ;
#   2015+    -> CLC2018
# Para 2005, CORINE_YEAR resuelve a 2006.
CORINE_YEAR <- OtsuFire:::get_corine_year(YEAR)
cat(sprintf("Combo: year=%d, scenario=%s, corine=%s\n",
            YEAR, SCENARIO, CORINE_YEAR))

# === A. INPUTS YEAR/SCENARIO-DEPENDENT (4) ===================================

# A1. Salida del modulo deterministic, principal entrada del supervised.
INTERNAL_DECISIONS <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "05_DECISIONS", "internal_decisions.gpkg"
)

# A2. Mosaico summer (RBR + DOY del verano del fuego).
CHANGE_INDEX <- file.path(COMPOSITE, RESULT_NAME,
                          sprintf("MinMin_%d_mosaic_res90m.tif", YEAR))

# A3. Mosaico autumn-winter (delayed change index, fuente de features G2_RBR_AW).
#     Cosmetico en cfg en 0.5.0 (HANDOFF §N+25), pero lo declaramos aqui.
RBR_AUTUMN <- file.path(COMPOSITE, "Autumn",
                        sprintf("mean_mean_%d_mosaic.tif", YEAR))

# A4. Hotspots MODIS. Pre-MODIS (<2000): NULL.
HOTSPOTS <- file.path(DATA_BASE, "Hotspots",
                      sprintf("hotspots_iberia_%d.geojson", YEAR))


# === B. VALIDACION EXTERNA YEAR-DEPENDENT (2) =================================

# B1. Effis-CA fuegos verano filtrados por mascara quemable (raster).
EFFIS_CA_TIF <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.tif", YEAR)
)

# B2. Lo mismo en formato shp.
EFFIS_CA_SHP <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.shp", YEAR)
)


# === C. CORINE QUINQUENAL (4) =================================================

CORINE_RASTER <- file.path(DATA_BASE, "Corine_Masks",
                           sprintf("CLC_%s_peninsula.tif", CORINE_YEAR))

CORINE_STRATA <- file.path(DATA_BASE, "Corine_Masks", "STRATA",
                           sprintf("strata_CLC_%s_res30.tif", CORINE_YEAR))


# C3. Look-up table de los strata.
CORINE_LUT <- file.path(DATA_BASE, "Corine_Masks", "LUT",
                        "lut_full_strata8_v1.csv")

# C4. Mascara binaria de superficie quemable derivada de Corine.
#     OJO: este input NO esta validado por el orchestrator al inicio.
#     Se valida solo dentro de internal-sup-unburned-deterministic.R linea 170
#     (cuando STEP A4 ejecuta). Si falta, el pipeline falla 5-10 min despues
#     de iniciado. Documentado en HANDOFF §N+25.
BURNEABLE_MASK <- file.path(
  DATA_BASE, "Corine_Masks",
  sprintf("burneable_mask_binary_corine_%s_ETRS89.tif", CORINE_YEAR)
)


# === D. ESTATICOS (4) =========================================================

# D1. Frontera de la Peninsula Iberica (polygon shapefile).
PENINSULA_SHP <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")

# D2. Topografia: stack con DEM y slope.
TOPO <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

# D3. Mascara del area de estudio en EPSG:3035 (raster).
MASK_TIF <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")

# D4. Idem en shp.
MASK_SHP <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.shp")


# === VERIFICACION GLOBAL (los 14 ficheros deben existir) ======================

stopifnot(
  # A) Year-dependent
  file.exists(INTERNAL_DECISIONS),
  file.exists(CHANGE_INDEX),
  file.exists(RBR_AUTUMN),
  file.exists(HOTSPOTS),
  # B) Validacion externa
  file.exists(EFFIS_CA_TIF),
  file.exists(EFFIS_CA_SHP),
  # C) Corine
  file.exists(CORINE_RASTER),
  file.exists(CORINE_STRATA),
  file.exists(CORINE_LUT),
  file.exists(BURNEABLE_MASK),
  # D) Estaticos
  file.exists(PENINSULA_SHP),
  file.exists(TOPO),
  file.exists(MASK_TIF),
  file.exists(MASK_SHP)
)

cat("\nLos 14 inputs criticos verificados:\n")
cat("\n  A) Year/scenario-dependent (4):\n")
cat("    internal_decisions   :", basename(INTERNAL_DECISIONS), "\n")
cat("    change_index (summer):", basename(CHANGE_INDEX), "\n")
cat("    delayed (autumn)     :", basename(RBR_AUTUMN), "\n")
cat("    hotspots             :", basename(HOTSPOTS), "\n")

cat("\n  B) Validacion externa (2):\n")
cat("    Effis-CA tif         :", basename(EFFIS_CA_TIF), "\n")
cat("    Effis-CA shp         :", basename(EFFIS_CA_SHP), "\n")

cat("\n  C) Corine quinquenal (4) [CORINE_YEAR=", CORINE_YEAR, "]:\n", sep="")
cat("    Corine raster        :", basename(CORINE_RASTER), "\n")
cat("    Corine strata        :", basename(CORINE_STRATA), "\n")
cat("    Corine LUT           :", basename(CORINE_LUT), "\n")
cat("    Burneable mask       :", basename(BURNEABLE_MASK), "\n")

cat("\n  D) Estaticos (4):\n")
cat("    Peninsula shp        :", basename(PENINSULA_SHP), "\n")
cat("    Topografia           :", basename(TOPO), "\n")
cat("    Mask studyarea tif   :", basename(MASK_TIF), "\n")
cat("    Mask studyarea shp   :", basename(MASK_SHP), "\n")


# === DIRECTORIO DE SALIDA DEL TUTORIAL =======================================

# Para no machacar el run real, redirigimos los outputs del tutorial a
# una carpeta paralela SUPERVISED_TUTORIAL/.
TUTORIAL_RESULT_DIR <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "SUPERVISED_TUTORIAL", SCENARIO
)
dir.create(TUTORIAL_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("\nDirectorio de outputs del tutorial:\n  ", TUTORIAL_RESULT_DIR, "\n")

# =============================================================================
# BLOQUE 1 - CONFIG: build_supervised_burned_config()
# =============================================================================
#
# QUE HACE
#   Construye un objeto cfg de clase 'otsufire_supervised_burned_config' que
#   contiene TODOS los paths, parametros y opciones que el pipeline necesita.
#   La funcion NO ejecuta nada: solo configura. La ejecucion la hace
#   run_oneyear_supervised_pipeline(cfg).
#
#
# DISENO DEL PAQUETE: PATRON CONFIG + RUN
# ---------------------------------------
# OtsuFire sigue el patron clasico de software cientifico (sklearn, MLflow,
# targets...): SEPARAR la configuracion de la ejecucion. Hay dos funciones
# publicas centrales en el modulo supervised:
#
#   1) build_supervised_burned_config()  -> "que quiero hacer"
#   2) run_oneyear_supervised_pipeline() -> "hazlo"
#
# Ventajas de este patron:
#   - VALIDACION TEMPRANA. El builder comprueba paths y tipos. Si algo falla,
#     el error sale en <1s en vez de a las 2 horas de pipeline.
#   - REPRODUCIBILIDAD. cfg es un objeto R inmutable que se puede serializar
#     con saveRDS(cfg, "experimento_baseline.rds"). Cualquier persona con
#     ese rds + los datos reproduce el experimento exacto.
#   - PROGRAMACION DE EXPERIMENTOS. configs <- list(baseline=..., e1=...);
#     lapply(configs, run_oneyear_supervised_pipeline). Limpio.
#
#
# ANATOMIA DEL OBJETO cfg (14 componentes en este caso)
# -----------------------------------------------------
# str(cfg, max.level = 1) revela:
#
#   $ scenario                  - chr "balanced"     (escenario operacional)
#   $ target_year               - int 2005           (ano objetivo)
#   $ inputs                    - List of 5          (paths a inputs)
#   $ run_name                  - chr "Min_Min"      (nombre del run)
#   $ output_dir                - chr ".../Results"  (raiz de outputs)
#   $ output_routes             - List of 19         (paths derivados)
#   (NOTA: burned_like_registry_path y negative_pool_policy fueron
#    eliminados del paquete el 2026-06-05; el registro burned-like ya no
#    existe y la policy de negativos es siempre all_sources, implicita.)
#   $ min_burned_pool_n         - int 5              (guard minimo positivos)
#   $ engine_root               - NULL               (legacy, no usado)
#   $ supervised_engine_root    - NULL               (legacy, no usado)
#   $ scripts_root              - NULL               (legacy, no usado)
#   $ tool_paths                - List of 4          (Anaconda3 binaries)
#   $ options                   - List of 11         (resto de opciones)
#
#
# LOS 5 CAMPOS DE cfg$inputs Y SU ESTADO EN 0.5.0
# ------------------------------------------------
# El cfg almacena 5 paths a inputs, pero NO todos son consumidos por el
# pipeline en la version actual (0.5.0). Conviene saberlo:
#
#   $ internal_decisions   -> SI consumido (STEP A1)
#   $ change_index         -> SI consumido (STEP A4 + features extraction)
#   $ delayed_change_index -> COSMETICO en 0.5.0; el orchestrator construye
#                             el path por convencion desde composite_base
#                             ver HANDOFF §N+25 para detalles
#   $ hotspots             -> SI consumido (STEP B3)
#   $ reference_burned_map -> RESERVADO para futura validacion externa,
#                             todavia no cableado en 0.5.0
#
# IMPORTANTE: el config builder solo acepta 4 paths de inputs; el pipeline
# en realidad lee 14 (HANDOFF §N+25). Los otros 10 los construye el
# orchestrator por convencion. En BLOQUE 0 declaramos los 14 explicitamente
# para validar al inicio. Aqui en BLOQUE 1 solo pasamos al builder los 4
# que la API publica acepta.
#
# Para el tutorial pasamos delayed_change_index al builder aunque sea
# cosmetico, por dos razones:
#   1. Cuando inspecciones cfg$inputs, veras todos los inputs reales del
#      experimento. Si alguien serializa el cfg con saveRDS() para
#      reproducibilidad, el path autumn quedara documentado.
#   2. Cuando se aborde el refactor post-paper (HANDOFF §N+25), el campo
#      sera consumido. Pasar el path ahora hace que el script sea
#      forward-compatible: misma llamada, distinto comportamiento del
#      paquete sin cambios en el script del usuario.
#
#
# CLASE S3
# --------
# El cfg tiene class = c("otsufire_supervised_burned_config", "list").
# Es una clase S3 ligera: una lista con etiqueta. Cuando
# run_oneyear_supervised_pipeline() recibe el cfg, lo primero que hace es
# verificar inherits(config, "otsufire_supervised_burned_config"). Si pasas
# una lista normal, falla con un mensaje claro. Es un contrato de tipos.
#
#
# DECISIONES METODOLOGICAS PARA BASELINE
# --------------------------------------
#   - Pool de negativos (siempre all_sources, implicito; ya no es opcion)
#       Usa los 4 sub-pools: internal_keep_qc (burned), deterministic_drop_hard,
#       random_burnable_background, otsu_patch_residual. (La antigua
#       alternativa "deterministic_direct" fue eliminada el 2026-06-05.)
#
#   - legacy_otsu_mode = "burnable_only"
#       El Otsu solo se aplica sobre pixeles que la mascara Corine considera
#       quemables. Sin esto, el Otsu veria agua, urbano, etc. y los umbrales
#       saldrian distorsionados.
#
#   - legacy_sample_n = 2000L
#       Cap del pool otsu_patch_residual. Sin esto, los ~13.000 parches
#       candidatos saturarian el training.
#
#   - legacy_otsu_threshold = 0
#       Descarta candidatos Otsu con RBR < 0. Guard minimo contra ruido.
#
#
# REFERENCIA EN EL PAQUETE
# ------------------------
#   R/supervised-config.R (15 KB) - definicion de build_supervised_burned_config
#   NAMESPACE linea 5             - export(build_supervised_burned_config)
#   NAMESPACE linea 1             - S3method(print, otsufire_supervised_burned_config)
#
# =============================================================================

cat("\n========== BLOQUE 1: CONFIG ==========\n")

# RBR_AUTUMN, INTERNAL_DECISIONS, CHANGE_INDEX, HOTSPOTS ya estan definidos
# y verificados en el BLOQUE 0. Los usamos directamente aqui.


# --- Construir el cfg -------------------------------------------------------

cfg <- build_supervised_burned_config(
  scenario                  = SCENARIO,
  internal_decisions        = INTERNAL_DECISIONS,
  change_index              = CHANGE_INDEX,        # mosaico summer (RBR + DOY)
  delayed_change_index      = RBR_AUTUMN,          # autumn-winter (cosmetico, §N+25)
  hotspots                  = HOTSPOTS,
  target_year               = YEAR,
  output_dir                = RESULTS,
  run_name                  = RESULT_NAME,
  options = list(
    data_base                       = DATA_BASE,
    composite_base                  = COMPOSITE,
    result_name                     = RESULT_NAME,
    # negative-pool policy is always all_sources now (implicit; only mode)
    legacy_otsu_mode                = "burnable_only",
    legacy_otsu_threshold           = 0,
    legacy_reference_otsu_threshold = 100,
    legacy_sample_n                 = 2000L,
    legacy_reuse_existing           = TRUE,
    legacy_write_output             = TRUE,
    unb_verbose                     = TRUE
  )
)

# Tool paths (Anaconda en este equipo). Se asignan DESPUES del builder
# porque son especificos del equipo y no parte de la configuracion del
# experimento. Asi el cfg sigue siendo portable: un colaborador en otro
# equipo solo cambia tool_paths.
cfg$tool_paths$python_exe              <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/python.exe"
cfg$tool_paths$gdal_polygonize_script  <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Scripts/gdal_polygonize.py"
cfg$tool_paths$gdalwarp_path           <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/gdalwarp.exe"
cfg$tool_paths$ogr2ogr_exe             <- "C:/Users/Olga.Viedma/AppData/Local/Anaconda3/Library/bin/ogr2ogr.exe"

stopifnot(all(sapply(cfg$tool_paths, file.exists)))
cat("Tool paths OK.\n")


# --- INSPECCION DEL cfg -----------------------------------------------------

cat("\nClase del cfg:\n")
print(class(cfg))
# Resultado esperado:
#   [1] "otsufire_supervised_burned_config" "list"
#
# La clase S3 ("otsufire_supervised_burned_config") es la primera. Permite
# que run_oneyear_supervised_pipeline() valide que es un cfg legitimo. La
# segunda ("list") es la clase base de R, que permite que cfg se comporte
# como una lista: cfg$inputs, cfg$options$legacy_sample_n, etc.

cat("\nEstructura de primer nivel del cfg:\n")
str(cfg, max.level = 1)

# --- cfg$inputs: paths a INPUTS criticos -----------------------------------
#
# Son las rutas a los ficheros que el pipeline lee. Tu las pasaste como
# argumentos al builder; aqui estan almacenadas.
#
# Recordatorio (HANDOFF §N+25): el config builder solo acepta 4 paths,
# pero el pipeline lee 14 ficheros distintos. Los otros 10 (Corine, topo,
# mask, peninsula, EFFIS validacion, burneable mask) los construye el
# orchestrator por convencion. Esos NO aparecen aqui.

cat("\ncfg$inputs:\n")
str(cfg$inputs, max.level = 1)
# Resultado:
#   $ internal_decisions  : List of 3   <- internal_decisions.gpkg + metadatos
#                                          CONSUMIDO en STEP A1
#   $ change_index        : List of 3   <- mosaico summer + metadatos
#                                          CONSUMIDO en STEP A4 + extraction
#   $ delayed_change_index: List of 3   <- mosaico autumn-winter + metadatos
#                                          COSMETICO en 0.5.0 (HANDOFF §N+25)
#                                          el path autumn lo construye el
#                                          orchestrator por convencion en
#                                          composite_base/Autumn/
#   $ hotspots            : List of 3   <- hotspots.geojson + metadatos
#                                          CONSUMIDO en STEP B3
#   $ reference_burned_map: NULL        <- reservado para futuro, no cableado


# --- cfg$options: las perillas operacionales -------------------------------
#
# Aqui estan los parametros que tu controlas para definir Baseline vs
# experimentos. Pasamos 11 al builder; el resto los rellena con defaults.

cat("\ncfg$options (11 que pasamos explicitos; el resto defaults):\n")
str(cfg$options, max.level = 1)
# Resultado:
#   $ data_base                      : "C:/.../1_DATA"
#                                       USADO por el orchestrator para
#                                       construir 9 paths implicitos (§N+25)
#   $ composite_base                 : "C:/.../Composites_90m"
#                                       USADO para construir el path autumn
#                                       y el mosaico summer (§N+25)
#   $ result_name                    : "Min_Min"
#   (negative_pool_policy fue eliminado; all_sources es siempre implicito)
#   $ legacy_otsu_mode               : "burnable_only"
#   $ legacy_otsu_threshold          : 0
#   $ legacy_reference_otsu_threshold: 100
#   $ legacy_sample_n                : 2000
#   $ legacy_reuse_existing          : TRUE
#   $ legacy_write_output             : TRUE
#   $ unb_verbose                    : TRUE


# --- cfg$output_routes: paths DERIVADOS ------------------------------------
#
# El builder calcula estos solos a partir de output_dir + target_year +
# result_name + scenario. Por eso son 19 paths sin que tu hayas pasado
# ninguno. La estructura SUPERVISED/<scenario>/<NN_FOLDER>/ es canonica
# del paquete.

cat("\ncfg$output_routes (19 paths derivados automaticamente):\n")
str(cfg$output_routes, max.level = 1)
# Subcarpetas: 01_POOLS, 02_FOLDS, 03_FEATURES, 04_MATRIX, 05_OOF,
# 07_FINAL_MODEL_V2, 08_SCORED, 09_FINAL_MAP, 11_CONSISTENCY_CHECKS, 99_LOGS.


# --- REDIRECCION DE OUTPUTS PARA EL TUTORIAL --------------------------------
#
# IMPORTANTE: el orchestrator escribe TODOS los outputs bajo
# cfg$output_routes$base. Si lo dejasemos como esta, machacariamos el run
# real de SUPERVISED/balanced/. Para el tutorial redirigimos a una carpeta
# alternativa SUPERVISED_TUTORIAL/.

cfg$output_routes$base       <- TUTORIAL_RESULT_DIR
cfg$output_routes$result_dir <- TUTORIAL_RESULT_DIR

cat("\nBLOQUE 1 completado. cfg listo para alimentar el pipeline.\n")



cat("\n========== BLOQUE 2: DIRECTORIOS Y RASTERS ==========\n")

# Estructura de directorios (igual que el orchestrator linea 779)
dirs <- list(
  `01_POOLS`          = file.path(TUTORIAL_RESULT_DIR, "01_POOLS"),
  `02_FOLDS`          = file.path(TUTORIAL_RESULT_DIR, "02_FOLDS"),
  `03_FEATURES`       = file.path(TUTORIAL_RESULT_DIR, "03_FEATURES"),
  `04_MATRIX`         = file.path(TUTORIAL_RESULT_DIR, "04_MATRIX"),
  `05_OOF`            = file.path(TUTORIAL_RESULT_DIR, "05_OOF"),
  `07_FINAL_MODEL_V2` = file.path(TUTORIAL_RESULT_DIR, "07_FINAL_MODEL_V2"),
  `08_SCORED`         = file.path(TUTORIAL_RESULT_DIR, "08_SCORED"),
  `09_FINAL_MAP`      = file.path(TUTORIAL_RESULT_DIR, "09_FINAL_MAP"),
  `99_LOGS_EMPTY`     = file.path(TUTORIAL_RESULT_DIR, "99_LOGS_EMPTY")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- Paths a rasters comunes (orchestrator linea 906) -----------------------
peninsula_shp <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")
topo_path     <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

corine_year   <- OtsuFire:::get_corine_year(YEAR)
corine_path   <- file.path(DATA_BASE, "Corine_Masks",
                           paste0("CLC_", corine_year, "_peninsula.tif"))

mask_tif_path <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")

rbr_aw_tif <- file.path(COMPOSITE, "Autumn",
                        paste0("mean_mean_", YEAR, "_mosaic.tif"))

stopifnot(file.exists(peninsula_shp), file.exists(topo_path),
          file.exists(corine_path), file.exists(mask_tif_path),
          file.exists(rbr_aw_tif))

cat("Inputs auxiliares verificados.\n")
cat("  Corine year:", corine_year, "\n")

# --- Cargar rasters y alinearlos al template (orchestrator linea 952) -------
cat("\nCargando y alineando rasters...\n")

topo       <- terra::rast(topo_path)
rbr_stack  <- terra::rast(CHANGE_INDEX)
rbr_summer <- rbr_stack[[1]]
doy_post   <- rbr_stack[[2]]

template_r <- rbr_summer

doy_post <- OtsuFire:::align_to_template(
  r = doy_post, template = template_r, method = "near",
  name = "doy_post", verbose = FALSE
)

dem_r <- OtsuFire:::align_to_template(
  r = topo[[1]], template = template_r, method = "bilinear",
  name = "dem", verbose = FALSE
)

slope_r <- OtsuFire:::align_to_template(
  r = topo[[2]], template = template_r, method = "bilinear",
  name = "slope", verbose = FALSE
)

corine_r <- OtsuFire:::align_to_template(
  r = terra::rast(corine_path), template = template_r, method = "near",
  name = "corine_r", verbose = FALSE
)

rbr_aw <- OtsuFire:::align_to_template(
  r = terra::rast(rbr_aw_tif)[[1]], template = template_r, method = "bilinear",
  name = "rbr_aw", verbose = FALSE
)

# --- Hotspots ---------------------------------------------------------------
hotspots_sf_base <- sf::read_sf(HOTSPOTS)
cat(sprintf("Hotspots cargados: %d puntos\n", nrow(hotspots_sf_base)))

# --- INSPECCION ---
cat("\nResumen rasters alineados:\n")
cat(sprintf("  rbr_summer: %d x %d, res %.0f m\n",
            nrow(rbr_summer), ncol(rbr_summer), terra::res(rbr_summer)[1]))
cat(sprintf("  doy_post:   %d x %d\n", nrow(doy_post), ncol(doy_post)))
cat(sprintf("  dem:        %d x %d\n", nrow(dem_r), ncol(dem_r)))
cat(sprintf("  slope:      %d x %d\n", nrow(slope_r), ncol(slope_r)))
cat(sprintf("  corine:     %d x %d\n", nrow(corine_r), ncol(corine_r)))
cat(sprintf("  rbr_aw:     %d x %d\n", nrow(rbr_aw), ncol(rbr_aw)))

# =============================================================================
# BLOQUE 3 - STEP A1: LEER internal_decisions.gpkg
# =============================================================================
#
# QUE HACE
#   Lee el GPKG producido por la fase deterministic con la clasificacion
#   por poligono: keep / drop / review.
#
# INPUTS  (de disco)
#   internal_decisions.gpkg, layer "internal_decisions"
#
# OUTPUTS (en memoria)
#   internal_sf - sf con todos los poligonos del ano clasificados
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R linea 1026 (STEP A1)
#
# =============================================================================

cat("\n========== BLOQUE 3: STEP A1 - LEER internal_decisions ==========\n")

internal_sf <- sf::read_sf(INTERNAL_DECISIONS, layer = "internal_decisions") |>
  dplyr::mutate(class = as.character(class_final))

# Las funciones del paquete sanitize_polygons y ensure_area_ha estan en el
# safety kit interno. Las llamamos desde alli.
sfkit <- OtsuFire:::make_sf_safety_kit(
  dirs = dirs, result_dir = TUTORIAL_RESULT_DIR, target_year = YEAR
)
sanitize_polygons    <- sfkit$sanitize_polygons
ensure_area_ha       <- sfkit$ensure_area_ha
drop_empty_sf        <- sfkit$drop_empty
check_sf             <- sfkit$check_sf
check_sf_if_nonempty <- sfkit$check_sf_if_nonempty
to_crs_safe          <- sfkit$to_crs_safe
safe_read_gpkg       <- sfkit$safe_read_gpkg
safe_write_gpkg      <- sfkit$safe_write_gpkg

internal_sf <- internal_sf |>
  sanitize_polygons() |>
  ensure_area_ha()

internal_sf <- drop_empty_sf(internal_sf, tag = "internal_sf",
                             dump_dir = dirs$`99_LOGS_EMPTY`)
check_sf(internal_sf, "internal_sf")

crs_master <- sf::st_crs(internal_sf)

cat(sprintf("internal_sf: %d poligonos\n", nrow(internal_sf)))
cat("Distribucion de clases:\n")
print(table(internal_sf$class))

# =============================================================================
# BLOQUE 4 - STEP A2 + A3: AUDIT y POOLS
# =============================================================================
#
# QUE HACE
#   A2: marca todos los poligonos con source = "internal".
#   A3: aplica audit_deterministic_pools() que valida y reclasifica
#       polygones marginales, y construye 3 pools:
#         burned_pool   (class == "keep")    -> positivos
#         review_pool   (class == "review")  -> dudosos, no se usan en train
#         scoring_pool  (todos)              -> universo a puntuar al final
#
# INPUTS  (memoria)
#   internal_clean (= internal_sf con source="internal")
#
# OUTPUTS (memoria + disco)
#   poolsA$internal_qc, $burned_pool, $review_pool, $scoring_pool
#   GPKG de QA en dirs$`01_POOLS`
#
# REFERENCIA EN EL PAQUETE
#   audit_deterministic_pools() en R/internal-sup-qa-pools.R
#   internal-sup-orchestrator.R linea 1065 (STEP A3)
#
# =============================================================================

cat("\n========== BLOQUE 4: STEP A2 + A3 - AUDIT + POOLS ==========\n")

# A2: marcar source
internal_clean <- internal_sf |>
  dplyr::mutate(source = "internal") |>
  sanitize_polygons() |>
  ensure_area_ha()

# A3: auditoria + construir pools
qa_det <- OtsuFire:::audit_deterministic_pools(
  internal_sf  = internal_clean,
  out_dir      = dirs$`01_POOLS`,
  prefix       = sprintf("%d_%s_deterministic_pool_qa", YEAR, SCENARIO),
  save_outputs = TRUE,
  verbose      = TRUE,
  raw_class_col = "class"
)

internal_qc <- qa_det$audited_sf |>
  dplyr::mutate(
    raw_class = as.character(raw_class),
    class     = as.character(class_audited)
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

# Construir burned_pool (positivos)
burned_pool <- internal_qc |>
  dplyr::filter(class == "keep") |>
  dplyr::mutate(
    source   = "internal_keep_qc",
    class    = "burned",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_B_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

burned_pool <- drop_empty_sf(burned_pool, tag = "burned_pool",
                             dump_dir = dirs$`99_LOGS_EMPTY`)

# Construir review_pool
review_pool <- internal_qc |>
  dplyr::filter(class == "review") |>
  dplyr::mutate(
    source   = "internal_review_qc",
    class    = "review",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_R_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

review_pool <- drop_empty_sf(review_pool, tag = "review_pool",
                             dump_dir = dirs$`99_LOGS_EMPTY`)

# Construir scoring_pool (todos los poligonos para puntuar al final)
scoring_pool <- internal_qc |>
  dplyr::mutate(
    class    = as.character(raw_class),
    source   = "deterministic_all_qc",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_D_", dplyr::row_number())
  ) |>
  sanitize_polygons() |>
  ensure_area_ha()

scoring_pool <- drop_empty_sf(scoring_pool, tag = "scoring_pool",
                              dump_dir = dirs$`99_LOGS_EMPTY`)

# --- INSPECCION ---
cat("\n--- Resultado A3 ---\n")
cat(sprintf("internal_qc:   %d poligonos\n", nrow(internal_qc)))
cat(sprintf("burned_pool:   %d (positivos)\n", nrow(burned_pool)))
cat(sprintf("review_pool:   %d (dudosos)\n", nrow(review_pool)))
cat(sprintf("scoring_pool:  %d (universo)\n", nrow(scoring_pool)))

# Sanity check: burned_pool no puede ser muy pequeno
n_burned <- nrow(burned_pool)
if (n_burned < 5L) {
  stop(sprintf("burned_pool insuficiente (n=%d, requerido>=5)", n_burned))
}
cat(sprintf("Sanity check: n_burned=%d >= 5 OK\n", n_burned))

# =============================================================================
# BLOQUE 5 - STEP A4: GENERAR UNBURNED (3 sub-pools)
# =============================================================================
#
# QUE HACE
#   Construye los pools de negativos en dos pasos:
#     Parte 1 (build_unburned_from_deterministic_decisions):
#       - unburned_hard:   poligonos clasificados drop por la deterministic
#       - unburned_random: celdas de fondo quemable fuera de buffers
#     Parte 2 (build_unburned_from_legacy_pipeline):
#       - otsu_sampled:    parches Otsu unburned residual (cap 2000)
#   Luego combina las dos partes en unburned_final_raw.
#
# DECISIONES METODOLOGICAS
#   pool de negativos: siempre all_sources, implicito (los sub-pools combinados)
#   exclude_buffer_m: separa negativos de positivos
#   n_random_cells: cuantos puntos random extraer
#   random_rbr_q: percentil maximo de RBR para considerar "fondo"
#   legacy_otsu_mode: "burnable_only" -> Otsu solo donde la mascara permite
#   legacy_sample_n: 2000L (cap del pool otsu_patch_residual)
#
# INPUTS
#   internal_decisions.gpkg, RBR raster, mascara, hotspots (informativo)
#
# OUTPUTS (memoria + disco GPKG)
#   unb$unburned_hard
#   unb$unburned_random
#   unb$unburned_final_raw  (todos los negativos combinados, sin sampleo)
#
# REFERENCIA EN EL PAQUETE
#   build_unburned_from_deterministic_decisions() en R/internal-sup-unburned-deterministic.R
#   build_unburned_from_legacy_pipeline()         en R/internal-sup-unburned-legacy.R
#   internal-sup-orchestrator.R lineas 1180-1372 (STEP A4)
#
# IMPORTANTE - DURACION
#   Esta etapa es la mas lenta del pipeline si reuse_existing=FALSE porque
#   tiene que correr Otsu + polygonize en Python (~10-20 min).
#   En 2005/balanced ya hay una corrida vieja con outputs en disco; con
#   legacy_reuse_existing = TRUE reutilizamos el polygonize y solo se
#   recalculan unburned_hard + unburned_random (1-2 min).
#
# =============================================================================

cat("\n========== BLOQUE 5: STEP A4 - GENERAR UNBURNED ==========\n")
cat("(esta etapa puede tardar 1-20 min segun caches)\n")

# Output del unburned-deterministic
unb_out_gpkg <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "UNBURNED", sprintf("%d_%s_unburned.gpkg", YEAR, SCENARIO)
)
dir.create(dirname(unb_out_gpkg), recursive = TRUE, showWarnings = FALSE)

# Carpeta de outputs Otsu legacy (compartida entre escenarios)
legacy_unb_root_dir <- file.path(TUTORIAL_RESULT_DIR, "_LEGACY_UNBURNED")
dir.create(legacy_unb_root_dir, recursive = TRUE, showWarnings = FALSE)

t0_a4 <- Sys.time()

# --- A4 Parte 1: deterministic drops + random background -------------------
res_det <- OtsuFire:::build_unburned_from_deterministic_decisions(
  target_year             = YEAR,
  scenario_name           = SCENARIO,
  data_base               = DATA_BASE,
  result_name             = RESULT_NAME,
  composite_base          = COMPOSITE,
  exclude_buffer_m        = cfg$options$unb_excl_buffer_m %||% 500,
  n_random_cells          = cfg$options$unb_n_random_cells %||% 15000L,
  random_rbr_q            = cfg$options$unb_random_rbr_q %||% 0.5,
  random_seed             = cfg$options$unb_random_seed %||% 42L,
  random_patch_size_cells = cfg$options$unb_random_patch_size_cells %||% 11L,
  overwrite_output        = TRUE,
  out_gpkg                = unb_out_gpkg,
  verbose                 = TRUE
)

unburned_hard   <- to_crs_safe(res_det$unburned_hard,   crs_master) |>
  dplyr::mutate(source = "deterministic_drop_hard")
unburned_random <- to_crs_safe(res_det$unburned_random, crs_master) |>
  dplyr::mutate(source = "random_burnable_background")
det_final_raw   <- to_crs_safe(res_det$unburned_final,  crs_master) |>
  dplyr::mutate(source = as.character(source))
exclusion_buffer <- to_crs_safe(res_det$exclusion_buffer, crs_master)

cat(sprintf("\nA4-Parte1 OK: hard=%d, random=%d\n",
            nrow(unburned_hard), nrow(unburned_random)))

# --- A4 Parte 2: Otsu unburned residual ------------------------------------
res_otsu <- OtsuFire:::build_unburned_from_legacy_pipeline(
  target_year              = YEAR,
  scenario_name            = SCENARIO,
  data_base                = DATA_BASE,
  result_name              = RESULT_NAME,
  composite_base           = COMPOSITE,
  severity_raster_path     = CHANGE_INDEX,
  legacy_code_dir          = NULL,
  python_exe               = cfg$tool_paths$python_exe,
  gdal_polygonize_script   = cfg$tool_paths$gdal_polygonize_script,
  gdalwarp_path            = cfg$tool_paths$gdalwarp_path,
  ogr2ogr_exe              = cfg$tool_paths$ogr2ogr_exe,
  otsu_mode                = cfg$options$legacy_otsu_mode %||% "burnable_only",
  otsu_threshold           = cfg$options$legacy_otsu_threshold %||% 0,
  reference_otsu_threshold = cfg$options$legacy_reference_otsu_threshold %||% 100,
  sample_n                 = cfg$options$legacy_sample_n %||% 2000L,
  reuse_existing           = cfg$options$legacy_reuse_existing %||% TRUE,
  write_unburned           = cfg$options$legacy_write_output %||% TRUE,
  out_root_dir             = legacy_unb_root_dir,
  verbose                  = TRUE
  # Nota: el resto de argumentos (min_pixels, buffers_m, core_thr, alpha_boost,
  # etc.) toman defaults internos. Si quisieras controlarlos todos, ver
  # version A del tutorial.
)

otsu_pool    <- to_crs_safe(res_otsu$unburned$legacy_unburned_pool,    crs_master)
otsu_sampled <- to_crs_safe(res_otsu$unburned$legacy_unburned_sampled, crs_master)
if (nrow(otsu_sampled) == 0L && nrow(otsu_pool) > 0L) otsu_sampled <- otsu_pool

otsu_raw <- otsu_sampled |>
  dplyr::mutate(
    source = "otsu_patch_residual",
    neg_type = dplyr::case_when(
      as.character(legacy_decision) == "drop"   ~ "otsu_patch_drop",
      as.character(legacy_decision) == "review" ~ "otsu_patch_review",
      as.character(legacy_decision) == "keep"   ~ "otsu_patch_keep",
      TRUE ~ as.character(neg_type)
    )
  )

cat(sprintf("A4-Parte2 OK: otsu_pool=%d, otsu_sampled=%d\n",
            nrow(otsu_pool), nrow(otsu_sampled)))

# --- A4 Parte 3: combinar fuentes ------------------------------------------
.bind_sf_safe <- function(a, b) {
  df_a <- sf::st_drop_geometry(a)
  df_b <- sf::st_drop_geometry(b)
  geom <- c(sf::st_geometry(a), sf::st_geometry(b))
  combined <- dplyr::bind_rows(df_a, df_b)
  combined[["geometry"]] <- geom
  sf::st_as_sf(combined, sf_column_name = "geometry", crs = sf::st_crs(a))
}
unburned_final_raw <- .bind_sf_safe(det_final_raw, otsu_raw)

# unburned_pool con etiquetas finales
unburned_pool <- unburned_final_raw |>
  sanitize_polygons() |>
  ensure_area_ha() |>
  dplyr::mutate(
    class    = "unburned",
    fire_uid = paste0(YEAR, "_", SCENARIO, "_U_", dplyr::row_number())
  )

unburned_pool <- drop_empty_sf(unburned_pool, tag = "unburned_pool",
                               dump_dir = dirs$`99_LOGS_EMPTY`)

elapsed_a4 <- as.numeric(difftime(Sys.time(), t0_a4, units = "mins"))
cat(sprintf("\nSTEP A4 completo en %.1f min\n", elapsed_a4))
cat(sprintf("unburned_pool total: %d\n", nrow(unburned_pool)))
cat("Composicion por source:\n")
print(table(unburned_pool$source))

# =============================================================================
# BLOQUE 6 - STEP A5 + A6: train_labeled y escritura de pools
# =============================================================================
#
# QUE HACE
#   A5: Concatena burned_pool + unburned_pool en train_labeled (positivos +
#       todos los negativos, sin aplicar caps todavia).
#   A6: Escribe todos los pools al GPKG 01_POOLS/<year>_<scenario>_pools.gpkg
#       como capas separadas para inspeccion visual posterior.
#
# OUTPUTS
#   train_labeled (memoria)
#   gpkg_pools_out (disco) con capas: burned_pool, unburned_pool, review_pool,
#                                     scoring_pool, train_labeled,
#                                     unburned_hard, unburned_random,
#                                     unburned_final_raw, exclusion_buffer
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R lineas 1384-1428 (STEP A5, A6)
#
# =============================================================================

cat("\n========== BLOQUE 6: STEP A5 + A6 - train_labeled + escritura ==========\n")

# A5: train_labeled
train_labeled <- dplyr::bind_rows(burned_pool, unburned_pool) |>
  sanitize_polygons() |>
  ensure_area_ha()

train_labeled <- drop_empty_sf(train_labeled, tag = "train_labeled",
                               dump_dir = dirs$`99_LOGS_EMPTY`)
check_sf(train_labeled, "train_labeled")
stopifnot(!anyDuplicated(train_labeled$fire_uid))

cat(sprintf("train_labeled: burned=%d | unburned=%d\n",
            sum(train_labeled$class == "burned"),
            sum(train_labeled$class == "unburned")))

# A6: escritura
gpkg_pools_out <- file.path(dirs$`01_POOLS`,
                            sprintf("%d_%s_pools.gpkg", YEAR, SCENARIO))

safe_write_gpkg(burned_pool,    gpkg_pools_out, layer = "burned_pool",    tag = "burned_pool_w")
safe_write_gpkg(unburned_pool,  gpkg_pools_out, layer = "unburned_pool",  tag = "unburned_pool_w")
safe_write_gpkg(review_pool,    gpkg_pools_out, layer = "review_pool",    tag = "review_pool_w")
safe_write_gpkg(scoring_pool,   gpkg_pools_out, layer = "scoring_pool",   tag = "scoring_pool_w")
safe_write_gpkg(train_labeled,  gpkg_pools_out, layer = "train_labeled",  tag = "train_labeled_w")
if (nrow(unburned_hard) > 0)
  safe_write_gpkg(unburned_hard, gpkg_pools_out, layer = "unburned_hard", tag = "unburned_hard_w")
if (nrow(unburned_random) > 0)
  safe_write_gpkg(unburned_random, gpkg_pools_out, layer = "unburned_random", tag = "unburned_random_w")
if (nrow(unburned_final_raw) > 0)
  safe_write_gpkg(unburned_final_raw, gpkg_pools_out, layer = "unburned_final_raw", tag = "unburned_final_raw_w")
if (nrow(exclusion_buffer) > 0)
  safe_write_gpkg(exclusion_buffer, gpkg_pools_out, layer = "exclusion_buffer", tag = "exclusion_buffer_w")

cat(sprintf("Pools escritos en: %s\n", gpkg_pools_out))
cat("Capas escritas:\n")
print(sf::st_layers(gpkg_pools_out)$name)

# =============================================================================
# BLOQUE 7 - STEP B2: SPATIAL FOLDS (make_block_folds)
# =============================================================================
#
# QUE HACE
#   Construye folds espaciales por bloques contiguos para evitar leakage
#   entre train/test. Hace 2 repeticiones (fold_rep1, fold_rep2) con seed
#   distinta para promediar la varianza de OOF.
#
# DECISIONES METODOLOGICAS
#   block_sizes_m = c(5000, 3000, 2000) - prueba 3 tamanos, elige el mejor
#   k_candidates  = c(5, 4, 3)          - prueba 3 numeros de folds
#   min_burned_units_per_fold = 3       - cada fold debe tener >=3 burned
#   min_pos_blocks_per_fold   = 10      - cada fold debe tener >=10 bloques
#                                         con al menos 1 burned
#   n_repeats = 2                        - dos repeticiones
#   seed_base = 42                       - reproducibilidad
#
# OUTPUTS (memoria + disco)
#   res_folds$selected            - cual block_size se eligio
#   res_folds$saved$train_gpkg    - GPKG con poligonos + columnas
#                                   block_id, fold_rep1, fold_rep2
#
# REFERENCIA EN EL PAQUETE
#   make_block_folds() en R/internal-sup-make-folds.R
#   internal-sup-orchestrator.R linea 1450 (STEP B2)
#
# =============================================================================

cat("\n========== BLOQUE 7: STEP B2 - SPATIAL FOLDS ==========\n")

train_labeled_sf <- safe_read_gpkg(gpkg_pools_out, "train_labeled", "train_labeled")

res_folds <- OtsuFire:::make_block_folds(
  train_labelled_sf         = train_labeled_sf,
  fire_id_col               = "fire_uid",
  split_unit                = "fire",
  block_sizes_m             = c(5000, 3000, 2000),
  k_candidates              = c(5, 4, 3),
  min_burned_units_per_fold = 3,
  min_pos_blocks_per_fold   = 10,
  n_repeats                 = 2,
  seed_base                 = 42,
  lonlat_action             = "transform",
  out_dir                   = dirs$`02_FOLDS`,
  target_year               = YEAR,
  verbose                   = TRUE,
  write_train_with_folds_gpkg = TRUE,
  write_blocks_gpkg           = TRUE,
  write_folds_csv             = TRUE
)

train_with_folds_gpkg <- res_folds$saved$train_gpkg
blocks_gpkg           <- res_folds$saved$blocks_gpkg
bs_sel                <- res_folds$selected$block_size_m

cat(sprintf("\nFolds OK. block_size_m elegido: %d\n", bs_sel))
cat(sprintf("train_with_folds_gpkg: %s\n", train_with_folds_gpkg))

# Inspeccion: ver distribucion de fold_rep1
twf <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "twf")
cat("\nDistribucion fold_rep1 x class:\n")
print(table(twf$fold_rep1, twf$class))
cat("\nDistribucion fold_rep2 x class:\n")
print(table(twf$fold_rep2, twf$class))

# =============================================================================
# BLOQUE 8 - STEP B3: EXTRACT FEATURES
# =============================================================================
#
# QUE HACE
#   Para cada poligono (train + scoring), extrae las 51 features canonicas
#   del whitelist por familia:
#     RBR summer:   rbr_p10, rbr_med, rbr_p90, rbr_iqr (4)
#     RBR autumn:   rbr_aw_p10, rbr_aw_med, rbr_aw_p90, rbr_aw_iqr,
#                   rbr_aw_valid_frac (5)
#     DOY:          doy_p10, doy_med, doy_p90, doy_iqr (4)
#     Topografia:   elev_*, slope_* (12)
#     Corine:       cor_agri_frac, cor_forest_frac, ... (7)
#     Persistencia: persist_delta, persist_ratio (2)
#     Hotspots:     hs_in_poly, hs_in_buffer, hs_used_n, hs_min_dist_m,
#                   hs_frp_sum, hs_frp_max, hs_conf_mean, hs_hiConf_n,
#                   hs_support_present, hs_no_support_when_available,
#                   hs_only_buffer_support, hotspot_available (12)
#                   [0.6.2: hs_any removed from the canonical whitelist]
#     Validez:      otras (4)
#
# DECISIONES METODOLOGICAS
#   use_hotspots = TRUE en 2005 (post-MODIS)
#   hs_buffer_m = 1000   - radio para hotspots cercanos
#   hs_start_month/end_month = 6,10  - estacion fuego
#   hi_conf_thr = 0.8    - umbral para hs_hiConf_n
#   cor_groups = grupos Corine ya en cfg$options$cor_groups
#
# OUTPUTS (disco)
#   features_geometry.gpkg con dos capas:
#     train_features    - poligonos de train con todas las features
#     scoring_features  - poligonos del universo a puntuar con features
#
# REFERENCIA EN EL PAQUETE
#   extract_features() en R/internal-sup-extract-features.R
#   internal-sup-orchestrator.R lineas 1681-1741 (STEP B3.3)
#
# =============================================================================

cat("\n========== BLOQUE 8: STEP B3 - EXTRACT FEATURES ==========\n")
cat("(esta etapa tarda 5-15 min segun tamano del area)\n")

# Re-leer del disco para asegurarnos de que extract_features ve los datos
# tal como los vera el orchestrator real
train_folds        <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "tf")
scoring_candidates <- safe_read_gpkg(gpkg_pools_out, "scoring_pool", "sp")

# Reproyectar hotspots si CRS no coincide
hotspots_sf <- hotspots_sf_base
if (nrow(hotspots_sf) > 0) {
  if (sf::st_crs(hotspots_sf) != sf::st_crs(train_folds)) {
    hotspots_sf <- sf::st_transform(hotspots_sf, sf::st_crs(train_folds))
  }
}
use_hotspots_flag <- nrow(hotspots_sf) > 0
cat(sprintf("use_hotspots_flag: %s (n=%d)\n", use_hotspots_flag, nrow(hotspots_sf)))

# Grupos Corine - tomados del cfg
cor_groups <- cfg$options$cor_groups

t0_b3 <- Sys.time()

feat <- OtsuFire:::extract_features(
  train_folds = train_folds,
  unlabeled   = scoring_candidates,

  id_col    = "fire_uid",
  class_col = "class",
  pos_lab   = "burned",
  neg_lab   = "unburned",

  build_features = TRUE,

  rbr_summer = rbr_summer,
  doy_post   = doy_post,
  rbr_aw     = rbr_aw,
  nbr_pre    = NULL,
  nbr_post   = NULL,
  dnbr       = NULL,
  dem        = dem_r,
  slope      = slope_r,
  corine_r   = corine_r,

  hotspots   = hotspots_sf,
  frp_col    = "frp",
  conf_col   = "confidence",
  year_col   = "year",
  month_col  = "month",
  hs_start_month = 6,
  hs_end_month   = 10,
  hi_conf_thr    = 0.8,

  year_target           = YEAR,
  hs_year_min_available = 2000,
  hs_buffer_m           = 1000,
  hs_missing_value      = -9999,
  hs_use_season_filter  = TRUE,

  cor_groups = cor_groups,

  use_doy        = TRUE,
  use_aw         = TRUE,
  use_nbr        = FALSE,
  use_hotspots   = use_hotspots_flag,

  max_cells_in_memory       = NULL,
  return_features           = TRUE,
  return_features_geometry  = TRUE,
  save_features_dir         = dirs$`03_FEATURES`,
  save_features_format      = "rds",
  save_features_gpkg        = TRUE,
  train_features_basename   = "train_features",
  scoring_features_basename = "scoring_features",
  train_features_layer      = "train_features",
  scoring_features_layer    = "scoring_features",
  verbose                   = TRUE
)

elapsed_b3 <- as.numeric(difftime(Sys.time(), t0_b3, units = "mins"))
cat(sprintf("\nSTEP B3 completo en %.1f min\n", elapsed_b3))

features_gpkg <- file.path(dirs$`03_FEATURES`, "features_geometry.gpkg")
cat("Capas generadas:\n")
print(sf::st_layers(features_gpkg)$name)

# Inspeccion de columnas
train_features_check <- safe_read_gpkg(features_gpkg, "train_features", "tfc")
feature_cols <- setdiff(names(train_features_check),
                        c("fire_uid", "class", "source", "poly_id",
                          "block_id", "fold_rep1", "fold_rep2", "geometry"))
cat(sprintf("\nFeatures extraidas: %d\n", length(feature_cols)))
cat("Familias detectadas:\n")
fam_counts <- list(
  RBR_SW    = sum(grepl("^rbr_(p10|med|p90|iqr|sd|mean)$", feature_cols)),
  RBR_AW    = sum(grepl("^rbr_aw_", feature_cols)),
  DOY       = sum(grepl("^doy_", feature_cols)),
  TOPO      = sum(grepl("^(elev|slope)_", feature_cols)),
  CORINE    = sum(grepl("^cor_", feature_cols)),
  PERSIST   = sum(grepl("^persist_", feature_cols)),
  HOTSPOTS  = sum(grepl("^(hs_|hotspot_)", feature_cols))
)
print(fam_counts)
cat(sprintf("Total checksum: %d\n", sum(unlist(fam_counts))))

# =============================================================================
# BLOQUE 9 - STEP C0 + C1: LEER FEATURES + PARAMS XGBOOST
# =============================================================================
#
# QUE HACE
#   C0: lee train_features y scoring_features del GPKG, hace join con folds
#       si hace falta.
#   C1: construye los hyperparametros de XGBoost. Calcula scale_pos_weight
#       como n_neg/n_pos para balancear clases en el loss.
#
# DECISIONES METODOLOGICAS (XGBOOST)
#   booster          = "gbtree"
#   objective        = "binary:logistic"
#   eval_metric      = c("logloss", "aucpr")  - logloss primero porque es
#                                              quien drivea early stopping
#   eta              = 0.06   (learning rate moderado)
#   max_depth        = 5      (arboles poco profundos para regularizar)
#   subsample        = 0.85   (sub-muestreo de filas por arbol)
#   colsample_bytree = 0.75   (sub-muestreo de cols por arbol)
#   min_child_weight = 5      (regularizacion)
#   gamma            = 0      (split sin penalizacion adicional)
#   scale_pos_weight = n_neg/n_pos  (balanceo automatico de clases)
#
# REFERENCIA EN EL PAQUETE
#   internal-sup-orchestrator.R lineas 1751-1832 (STEP C0, C1)
#
# =============================================================================

cat("\n========== BLOQUE 9: STEP C0 + C1 - FEATURES + XGB PARAMS ==========\n")

train_features   <- safe_read_gpkg(features_gpkg, "train_features",   "tf2")
scoring_features <- safe_read_gpkg(features_gpkg, "scoring_features", "sf2")

# Asegurar que train_features tiene los folds
needed_folds <- c("block_id", "fold_rep1", "fold_rep2")
if (!all(needed_folds %in% names(train_features))) {
  cat("Joining folds into train_features...\n")
  twf2 <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds", "twf2")
  folds_df <- sf::st_drop_geometry(twf2) |>
    dplyr::select(fire_uid, dplyr::any_of(needed_folds))
  train_features <- train_features |> dplyr::left_join(folds_df, by = "fire_uid")
}

# Drop geometry para modelado
labelled    <- sf::st_drop_geometry(train_features)
burned_like <- sf::st_drop_geometry(scoring_features)
labelled_df <- labelled

# Construir params XGBoost
n_pos <- sum(labelled$class == "burned",   na.rm = TRUE)
n_neg <- sum(labelled$class == "unburned", na.rm = TRUE)
spw   <- if (n_pos > 0) n_neg / n_pos else 1

params <- list(
  booster          = "gbtree",
  objective        = "binary:logistic",
  eval_metric      = c("logloss", "aucpr"),
  eta              = 0.06,
  max_depth        = 5,
  subsample        = 0.85,
  colsample_bytree = 0.75,
  min_child_weight = 5,
  gamma            = 0,
  scale_pos_weight = spw
)

cat(sprintf("n_pos=%d, n_neg=%d, scale_pos_weight=%.3f\n", n_pos, n_neg, spw))
cat("XGB params:\n")
print(params)

# =============================================================================
# BLOQUE 10 - STEP C2: OOF PIPELINE (entrenamiento por folds + diagnostico)
# =============================================================================
#
# QUE HACE
#   Para cada repeticion (fold_rep1, fold_rep2) y para cada fold k=1..5:
#     1. Marca fold k como test, resto como train.
#     2. Aplica caps de sampling:
#          contextual_exclusion (deterministic_drop_hard, neg_type != spectral) cap 0.25 x n_burned
#          spectral_hard_negative (deterministic_drop_hard con neg_type == "spectral_reject_medium") cap 1.0
#          random_burnable_background cap 1.0
#          otsu_patch_residual cap 1.0
#     3. Entrena XGBoost en train, predice en test.
#     4. Acumula predicciones.
#   Al final, cada fila tiene n_preds=2 predicciones (una por cada rep).
#   Calcula metricas (accuracy, precision, recall, etc.) en multiples
#   thresholds, identifica el threshold optimo segun varios criterios.
#
# DECISIONES METODOLOGICAS
#   feature_whitelist_override = NULL  (canon 51)
#   feature_weights            = NULL  (uniforme 1.0)
#   median_from = "labelled"   - imputacion de NA por mediana de cada feature
#                                en el set labelled
#
# OUTPUTS (disco)
#   05_OOF/<prefix>_oof_metrics_summary.txt
#   05_OOF/<prefix>_labeled_oof_summary.gpkg     - cada poligono con p_oof_mean
#   05_OOF/<prefix>_oof_agg.csv
#   05_OOF/<prefix>_oof_long.csv
#   05_OOF/<prefix>_oof_metrics_by_threshold.csv
#   05_OOF/<prefix>_oof_best_thresholds.csv
#   04_MATRIX/<prefix>_dm.rds  (DesignMatrix)
#
# REFERENCIA EN EL PAQUETE
#   run_dm_oof_pipeline() en R/internal-sup-oof-wrapper.R
#   internal-sup-orchestrator.R lineas 1834-1870 (STEP C2)
#
# IMPORTANTE
#   Este bloque tarda 10-20 min en 2005 (es la etapa mas lenta del modelado).
#
# =============================================================================

cat("\n========== BLOQUE 10: STEP C2 - OOF PIPELINE ==========\n")
cat("(esta etapa tarda 10-20 min)\n")

prefix_oof <- sprintf("%d_%s_patch", YEAR, SCENARIO)
prefix     <- sprintf("%d_%s_patch_certified", YEAR, SCENARIO)

t0_c2 <- Sys.time()

pipe1 <- OtsuFire:::run_dm_oof_pipeline(
  labelled    = labelled,
  burned_like = burned_like,
  labelled_df = labelled_df,
  params      = params,

  result_dir  = TUTORIAL_RESULT_DIR,
  target_year = YEAR,
  prefix      = prefix_oof,

  id_cols     = c("fire_uid", "class", "source", "poly_id",
                  "block_id", "fold_rep1", "fold_rep2"),
  cat_cols    = character(0),
  hs_n_col    = "hs_used_n",
  hs_conf_col = "hs_conf_mean",
  hs_frp_col  = "hs_frp_max",
  median_from = "labelled",
  save_dir_dm = dirs$`04_MATRIX`,
  save_prefix = paste0(YEAR, "_", SCENARIO, "_patch"),
  overwrite   = TRUE,

  # ----- BASELINE: NULL en ambos -----
  feature_whitelist_override = NULL,
  feature_weights            = NULL,

  labelled_gpkg  = train_with_folds_gpkg,
  labelled_layer = "train_with_folds"
)

elapsed_c2 <- as.numeric(difftime(Sys.time(), t0_c2, units = "mins"))
cat(sprintf("\nSTEP C2 completo en %.1f min\n", elapsed_c2))

cat("Files generados por OOF pipeline:\n")
print(pipe1$files)

# Inspeccion del summary
oof_summary_path <- file.path(dirs$`05_OOF`,
                              paste0(prefix_oof, "_oof_metrics_summary.txt"))
cat("\n--- oof_metrics_summary.txt ---\n")
cat(readLines(oof_summary_path), sep = "\n")

# =============================================================================
# BLOQUE 11 - STEP C3: FINAL MODEL + SCORING + FINAL MAP
# =============================================================================
#
# QUE HACE
#   1. Entrena el modelo FINAL con TODOS los datos labelled (sin OOF).
#   2. Aplica el modelo al universo scoring_features para puntuar.
#   3. Construye el mapa final con clase asignada por threshold optimo
#      heredado del OOF.
#   4. Identifica burned_like (poligonos no etiquetados que el modelo
#      considera quemados).
#
# DECISIONES METODOLOGICAS
#   contextual_exclusion_to_burned_ratio = 0.25
#   spectral_hard_negative_to_burned_ratio = 1.0
#   random_to_burned_ratio = 1.0
#   otsu_unburned_to_burned_ratio = 1.0
#   feature_whitelist_override = NULL
#   feature_weights = NULL
#   preyear_overlap_threshold = cfg$options$currentyear_preyear_overlap_thr
#   hotspot_density_threshold = cfg$options$currentyear_hotspot_density_thr
#   temporal_penalty_floor    = cfg$options$currentyear_temporal_penalty_floor
#
# OUTPUTS (disco)
#   07_FINAL_MODEL_V2/<prefix>_certified_final_model.rds
#   07_FINAL_MODEL_V2/<prefix>_certified_recipe.rds
#   07_FINAL_MODEL_V2/<prefix>_certified_meta.txt
#   07_FINAL_MODEL_V2/<prefix>_certified_feature_importance.csv
#   07_FINAL_MODEL_V2/<prefix>_certified_model_summary.txt
#   07_FINAL_MODEL_V2/<prefix>_certified_training_ok.csv
#   07_FINAL_MODEL_V2/<prefix>_certified_training_ok.gpkg
#   08_SCORED/<prefix>_certified_scored.csv
#   08_SCORED/<prefix>_certified_scored.gpkg
#   09_FINAL_MAP/<prefix>_certified_final_map.gpkg
#   09_FINAL_MAP/<prefix>_certified_burned_like_scored.gpkg
#
# REFERENCIA EN EL PAQUETE
#   run_train_final_model_and_export_final_map() en
#     R/internal-sup-train-final-wrapper.R
#   Funcion interna real: train_final_model_direct() en
#     R/internal-sup-train-final-direct.R
#   internal-sup-orchestrator.R lineas 1875-1925 (STEP C3)
#
# =============================================================================

cat("\n========== BLOQUE 11: STEP C3 - FINAL MODEL + SCORING ==========\n")
cat("(esta etapa tarda 5-15 min)\n")

oof_agg_csv      <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_oof_agg.csv"))
oof_summary_gpkg <- file.path(dirs$`05_OOF`, paste0(prefix_oof, "_labeled_oof_summary.gpkg"))
stopifnot(file.exists(oof_agg_csv), file.exists(oof_summary_gpkg))

t0_c3 <- Sys.time()

pipe2 <- OtsuFire:::run_train_final_model_and_export_final_map(
  result_dir     = TUTORIAL_RESULT_DIR,
  qa             = oof_agg_csv,
  labelled_gpkg  = features_gpkg,
  labelled_layer = "train_features",
  out_dir        = dirs$`07_FINAL_MODEL_V2`,
  prefix         = prefix,
  overwrite      = TRUE,
  verbose        = TRUE,

  qa_labelled_gpkg  = oof_summary_gpkg,
  qa_labelled_layer = "labeled_oof_summary",

  labelled_features_gpkg  = features_gpkg,
  labelled_features_layer = "train_features",

  model_rds  = file.path(dirs$`07_FINAL_MODEL_V2`,
                         paste0(prefix, "_final_model.rds")),
  recipe_rds = file.path(dirs$`07_FINAL_MODEL_V2`,
                         paste0(prefix, "_recipe.rds")),

  unlabeled_gpkg  = features_gpkg,
  unlabeled_layer = "scoring_features",

  out_score_dir = dirs$`08_SCORED`,
  out_map_dir   = dirs$`09_FINAL_MAP`,

  export_burned_like        = TRUE,
  preyear_overlap_threshold = cfg$options$currentyear_preyear_overlap_thr %||% 0.5,
  hotspot_density_threshold = cfg$options$currentyear_hotspot_density_thr %||% 1.0,
  temporal_penalty_floor    = cfg$options$currentyear_temporal_penalty_floor %||% 0.1,

  # ----- BASELINE: caps canonicos + sin overrides + sin weights -----
  contextual_exclusion_to_burned_ratio   = 0.25,
  spectral_hard_negative_to_burned_ratio = 1.0,
  random_to_burned_ratio                 = 1.0,
  otsu_unburned_to_burned_ratio          = 1.0,
  feature_whitelist_override             = NULL,
  feature_weights                        = NULL
)

elapsed_c3 <- as.numeric(difftime(Sys.time(), t0_c3, units = "mins"))
cat(sprintf("\nSTEP C3 completo en %.1f min\n", elapsed_c3))

# Inspeccion del meta del modelo final
meta_path <- file.path(dirs$`07_FINAL_MODEL_V2`,
                       paste0(prefix, "_meta.txt"))
if (file.exists(meta_path)) {
  cat("\n--- meta.txt del modelo final ---\n")
  cat(readLines(meta_path), sep = "\n")
}

# =============================================================================
# BLOQUE 12 - VERIFICACION FINAL: Baseline canonico
# =============================================================================
#
# QUE HACE
#   Verifica que el meta.txt confirme la configuracion Baseline canonical:
#     feature_whitelist_override_applied = FALSE
#     feature_weights_applied            = FALSE
#   Y que las metricas tengan sentido (drops baja probabilidad, keeps alta).
#
# =============================================================================

cat("\n========== BLOQUE 12: VERIFICACION BASELINE ==========\n")

if (file.exists(meta_path)) {
  meta_lines <- readLines(meta_path)
  override_line <- grep("feature_whitelist_override_applied", meta_lines, value = TRUE)
  weights_line  <- grep("feature_weights_applied", meta_lines, value = TRUE)
  n_x_line      <- grep("^n_x_cols:", meta_lines, value = TRUE)
  n_feat_line   <- grep("^n_feature_cols:", meta_lines, value = TRUE)

  cat("Lineas relevantes del meta.txt:\n")
  cat(" ", override_line, "\n")
  cat(" ", weights_line, "\n")
  cat(" ", n_x_line, "\n")
  cat(" ", n_feat_line, "\n")

  is_baseline <- grepl("FALSE", override_line) && grepl("FALSE", weights_line)
  if (is_baseline) {
    cat("\nBaseline canonical CONFIRMADO.\n")
  } else {
    cat("\nATENCION: el meta.txt NO refleja Baseline canonical.\n")
  }
}

# Verificacion del mapa final
final_map_path <- file.path(dirs$`09_FINAL_MAP`,
                            paste0(prefix, "_final_map.gpkg"))
if (file.exists(final_map_path)) {
  fm <- sf::read_sf(final_map_path, layer = "final_map_full")

  drops <- fm$p_burned_model[fm$class_final == "drop"]
  keeps <- fm$p_burned_model[fm$class_final == "keep"]

  drops_med <- median(drops, na.rm = TRUE)
  keeps_med <- median(keeps, na.rm = TRUE)

  cat(sprintf("\nDrops mediana: %.4f (esperado <0.1)\n", drops_med))
  cat(sprintf("Keeps mediana: %.4f (esperado >0.9)\n", keeps_med))

  if (drops_med < 0.1 && keeps_med > 0.9) {
    cat("\nBASELINE 0.5.0 PASSED.\n")
  } else {
    cat("\nBaseline NO pasa los umbrales esperados.\n")
  }
}

cat("\n========== TUTORIAL COMPLETADO ==========\n")
cat(sprintf("Outputs en: %s\n", TUTORIAL_RESULT_DIR))
