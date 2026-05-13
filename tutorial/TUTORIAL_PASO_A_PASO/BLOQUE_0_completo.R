# =============================================================================
# BLOQUE 0 - SETUP COMPLETO: TODA LA CONFIGURACION DEL USUARIO
# =============================================================================
#
# QUE HACE
#   Carga el paquete OtsuFire en modo desarrollo y declara EXPLICITAMENTE
#   todas las variables del experimento que el usuario debe configurar:
#     - 2 parametros del combo (year + scenario)
#     - 13 paths de inputs (lo que el pipeline lee de disco)
#     - 18 constantes metodologicas (politicas + thresholds + toggles)
#
#   Este bloque es el UNICO que el usuario tiene que entender y modificar
#   para correr el pipeline. Todo lo que viene despues (BLOQUE 1+) es
#   consumo de estas variables.
#
#
# DOS NIVELES DE LECTURA
# ----------------------
# Cada seccion del bloque incluye comentarios narrativos (que se hace y
# por que) + bloques de DEBUG (cross-references al handoff con detalles
# tecnicos para diagnosticar problemas del paquete).
#
#
# ESTADO DE LA API EN 0.5.0 (resumen para contextualizar)
# -------------------------------------------------------
# El paquete OtsuFire 0.5.0 tiene una API de configuracion incompleta:
# de las 33 variables que conceptualmente forman el experimento, solo
# parte se aceptan via el config builder publico. El resto las construye
# el orchestrator por convencion (paths) o las lee de constantes globales
# hardcoded (parametros). Esto se resolvera en el refactor post-paper.
#
# Cross-references:
#   §N+25: 13 paths de inputs - solo 4 entran via cfg$inputs en 0.5.0.
#   §N+26: politica negative_pool_policy - redundancia a simplificar.
#   §N+27: campo burned_like_registry_path - linea de investigacion
#          abortada, eliminar en refactor.
#   §N+28: 21 constantes hardcoded en orchestrator - decidir cuales
#          exponer post-refactor.
#
# Para que este tutorial sea forward-compatible y sirva como spec del
# refactor, declaramos AQUI las 33 variables explicitamente. Despues
# en BLOQUE 1 las separamos en:
#   - Capa B (comentada): API ideal post-refactor donde el usuario pasa
#                          las 33 variables al config builder.
#   - Capa C (ejecutable): API actual 0.5.0 donde solo se pueden pasar
#                          algunas; el resto se usan en bloques posteriores
#                          o quedan documentadas.
#
# =============================================================================


# =============================================================================
# 1) CARGA DEL PAQUETE
# =============================================================================

cat("\n========== BLOQUE 0: SETUP ==========\n")

PKG_ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"

# Carga limpia: si ya estaba cargado, descargar primero para asegurar
# que tomamos la version actual del codigo en disco (importante en
# desarrollo activo donde el codigo del paquete cambia).
if ("OtsuFire" %in% loadedNamespaces()) {
  try(unloadNamespace("OtsuFire"), silent = TRUE)
}
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

cat("OtsuFire version:", as.character(utils::packageVersion("OtsuFire")), "\n")
stopifnot(utils::packageVersion("OtsuFire") >= "0.5.0")

# DEBUG: pkgload::load_all expone tanto los exports del NAMESPACE como
# las funciones internas. Esto permite usar OtsuFire:::xxx para acceder
# a internas, lo que el tutorial necesita en bloques posteriores para
# "abrir la nevera" del paquete y llamar a funciones step-by-step en
# vez de pasar por el wrapper run_oneyear_supervised_pipeline().


# =============================================================================
# 2) PARAMETROS DEL COMBO (year + scenario)
# =============================================================================
#
# Para correr el tutorial con otro combo, cambiar SOLO estas dos lineas.
# Todo el resto del BLOQUE 0 (paths, constantes) se adapta automaticamente.

YEAR     <- 2005L
SCENARIO <- "balanced"

# DEBUG: YEAR debe ser integer (sufijo L). Si pasas un double (2005 sin L),
# las funciones del paquete que esperan target_year integer pueden fallar
# con coercion warnings. SCENARIO debe ser uno de: "balanced", "lax",
# "original" (los tres soportados por el deterministic stage).

# Resolucion del año Corine quinquenal usando funcion interna.
# Mapeo (ver R/utils-corine.R::get_corine_year):
#   1985-2002 -> 2000 ; 2003-2008 -> 2006 ; 2009-2014 -> 2012 ;
#   2015+    -> 2018
# Para YEAR=2005, CORINE_YEAR resuelve a "2006" (devuelto como character).
CORINE_YEAR <- OtsuFire:::get_corine_year(YEAR)

cat(sprintf("\nCombo: year=%d, scenario=%s, corine=%s\n",
            YEAR, SCENARIO, CORINE_YEAR))

# DEBUG: get_corine_year devuelve character, no integer. Por eso los
# sprintf de paths Corine usan %s en lugar de %d. Si en el futuro
# se cambia a integer, actualizar los sprintf.


# =============================================================================
# 3) PATHS DEL PROYECTO
# =============================================================================

DATA_BASE <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA"
COMPOSITE <- file.path(DATA_BASE, "Imagery", "Composites_90m")
RESULTS   <- file.path(DATA_BASE, "Results")

# RESULT_NAME identifica el "run" del pipeline. Min_Min hace referencia
# al metodo de composicion temporal de los mosaicos RBR (Min de Min).
# El paquete soporta otros (Mean_Mean, Min_Mean, etc.) pero el Baseline
# canonico usa Min_Min.
RESULT_NAME <- "Min_Min"


# =============================================================================
# 4) LOS 13 INPUTS DEL PIPELINE
# =============================================================================
#
# Estos son los 13 ficheros que el orchestrator del modulo supervised lee
# para procesar un combo. Estan organizados en 4 grupos por tipo de
# dependencia:
#
#   A) Year/scenario-dependent (4)  - cambian con (year, scenario)
#   B) Validacion externa year-dep  (1) - cambia solo con year
#   C) Corine quinquenal           (4) - cambian con CORINE_YEAR
#   D) Estaticos                   (4) - no cambian nunca
#
# DEBUG: HANDOFF §N+25. En 0.5.0, solo 4 de estos 13 entran al config
# builder via cfg$inputs (internal_decisions, change_index,
# delayed_change_index, hotspots). Los otros 9 los construye el
# orchestrator por convencion en lineas 906-934 de internal-sup-
# orchestrator.R + el burneable mask en helpers internos.
# El refactor post-paper hara que TODOS entren via cfg$inputs.


# --- A) Year/scenario-dependent (4) -----------------------------------------

# A1. Salida del modulo deterministic (entrada principal del supervised).
INTERNAL_DECISIONS <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "DETERMINISTIC", SCENARIO,
  "05_DECISIONS", "internal_decisions.gpkg"
)

# A2. Mosaico summer (RBR + DOY del verano del fuego).
CHANGE_INDEX <- file.path(COMPOSITE, RESULT_NAME,
                          sprintf("MinMin_%d_mosaic_res90m.tif", YEAR))

# A3. Mosaico autumn-winter (delayed change index, fuente G2_RBR_AW).
RBR_AUTUMN <- file.path(COMPOSITE, "Autumn",
                        sprintf("mean_mean_%d_mosaic.tif", YEAR))

# A4. Hotspots MODIS. Para anos pre-MODIS (<2000), HOTSPOTS no existira
# en disco; el pipeline detecta y procesa sin features de hotspots.
HOTSPOTS <- file.path(DATA_BASE, "Hotspots",
                      sprintf("hotspots_iberia_%d.geojson", YEAR))


# --- B) Validacion externa year-dependent (1) -------------------------------

# B1. Effis-CA fuegos verano filtrados por mascara quemable.
# Es el ground truth EXTERNO para validar las predicciones del modelo.
# DEBUG: HANDOFF §N+25. En 0.5.0 corresponde al campo
# cfg$inputs$reference_burned_map (cosmetico) + ref_tif_path/ref_shp_path
# (construidos por convencion). El refactor consolidara ambos bajo el
# campo publico reference_burned_map.
EFFIS_CA_TIF <- file.path(
  DATA_BASE, "Fires", "Validation_fires_burneable_verano",
  sprintf("Effis_CA_%d_maskKeep_summer.tif", YEAR)
)

# El shp companion es derivable desde el tif (cambia .tif por .shp).
# Lo declaramos explicitamente porque el orchestrator 0.5.0 los lee como
# rutas separadas. Post-refactor sera derivable internamente desde
# reference_burned_map.
EFFIS_CA_SHP <- sub("\\.tif$", ".shp", EFFIS_CA_TIF)


# --- C) Corine quinquenal (4) -----------------------------------------------

# C1. Raster Corine Land Cover (cobertura del suelo).
# Fuente de features G3_CORINE.
CORINE_RASTER <- file.path(DATA_BASE, "Corine_Masks",
                           sprintf("CLC_%s_peninsula.tif", CORINE_YEAR))

# C2. Strata raster (estratificacion derivada de Corine para el muestreo
# espacial estratificado en STEP B1 - block folds).
CORINE_STRATA <- file.path(DATA_BASE, "Corine_Masks", "STRATA",
                           sprintf("strata_CLC_%s_res30.tif", CORINE_YEAR))

# C3. Look-up table de los strata (no depende del año).
CORINE_LUT <- file.path(DATA_BASE, "Corine_Masks", "LUT",
                        "lut_full_strata8_v1.csv")

# C4. Mascara binaria de superficie quemable derivada de Corine.
# DEBUG: HANDOFF §N+25 §1.3. Este input es el unico que el orchestrator
# NO valida en su bloque inicial de stopifnot. Lo lee dentro de
# internal-sup-unburned-deterministic.R linea 170. Si falta, el pipeline
# falla 5-10 min despues de iniciado, no al inicio. Por eso lo
# verificamos manualmente en BLOQUE 0.
BURNEABLE_MASK <- file.path(
  DATA_BASE, "Corine_Masks",
  sprintf("burneable_mask_binary_corine_%s_ETRS89.tif", CORINE_YEAR)
)


# --- D) Estaticos (4) -------------------------------------------------------

# D1. Frontera de la Peninsula Iberica.
PENINSULA_SHP <- file.path(DATA_BASE, "Borders", "Iberian_peninsula.shp")

# D2. Topografia: stack con DEM (banda 1) y slope (banda 2).
# Fuente de features G4_TOPO.
TOPO <- file.path(DATA_BASE, "Topography", "elevation_slope.tif")

# D3, D4. Mascara del area de estudio en EPSG:3035 (raster + shp).
MASK_TIF <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.tif")
MASK_SHP <- file.path(DATA_BASE, "Mask_StudyArea", "mask_Peninsula_3035.shp")


# =============================================================================
# 5) CONSTANTES METODOLOGICAS (18 variables)
# =============================================================================
#
# Estas son las decisiones metodologicas del Baseline canonico. Estan
# agrupadas por funcion. La mayoria de usuarios NO necesita cambiarlas:
# los valores aqui son los del Baseline publicado.
#
# DEBUG: HANDOFF §N+28. En 0.5.0:
#   - 8 de estas se aceptan via cfg$options + 1 via cfg top-level.
#   - 9 estan hardcoded en internal-sup-orchestrator.R y NO son
#     configurables sin modificar el codigo del paquete.
# El tutorial las declara TODAS aqui para que la Capa B del BLOQUE 1
# muestre la API ideal post-refactor.


# --- 5.1) Politica de muestreo de negativos (1) -----------------------------

# Define que fuentes de negativos entran al training pool.
# Valores soportados:
#   "all_sources"          - 3 fuentes: deterministic drops + random
#                            background + Otsu unburned patches.
#                            Es el modo OPERACIONAL del Baseline canonico.
#   "deterministic_direct" - 2 fuentes: solo drops + random. Sin Otsu.
#                            Modo de validacion interna del paquete.
#                            DEBUG: HANDOFF §N+26 - eliminar en refactor.
NEGATIVE_POOL_POLICY <- "all_sources"


# --- 5.2) Guard de tamaño minimo del pool burned (1) ------------------------

# Si el deterministic stage produce menos de N polygonos burned para el
# combo, el supervised aborta limpiamente con error informativo (no se
# puede entrenar un modelo con tan pocos positivos).
# Este guard fue el que aborto los 3 combos de 1997 (n=0/1, ver §N+24).
MIN_BURNED_POOL_N <- 5L


# --- 5.3) Reproducibilidad: seeds aleatorios (2) ----------------------------

# Seed del muestreo random_burnable_background (rama deterministic_direct
# y rama all_sources Part 1).
RANDOM_SEED <- 42L

# Seed del muestreo Otsu unburned (rama all_sources Part 2).
# DEBUG: HANDOFF §N+28. El paquete tiene DOS seeds independientes,
# ambos hardcoded a 42. El refactor podria unificarlos a uno solo si
# se decide que la independencia no aporta valor.
LEGACY_RANDOM_SEED <- 42L


# --- 5.4) Filtro temporal entre años (3) ------------------------------------
#
# Penalizan polygonos del año actual que solapan demasiado con polygonos
# del año anterior (proxy: "fuegos persistentes probablemente no son
# fuegos reales"). DEBUG: HANDOFF §N+28 - actualmente hardcoded.

# Umbral de overlap espacial para considerar conflicto temporal.
CURRENTYEAR_PREYEAR_OVERLAP_THR <- 0.70

# Densidad minima de hotspots/ha para que un polygono pase el filtro.
CURRENTYEAR_HOTSPOT_DENSITY_THR <- 0.001

# Floor de la penalty score (no baja de este valor).
CURRENTYEAR_TEMPORAL_PENALTY_FLOOR <- 0.10


# --- 5.5) Otsu pipeline - parametros principales (4) ------------------------
#
# Controlan la rama Otsu del all_sources mode (legacy_otsu_*).

# Modo de aplicacion del Otsu. "burnable_only" = solo sobre pixeles
# que la mascara Corine considera quemables (evita que el Otsu vea
# agua, urbano, etc. y distorsione los umbrales).
LEGACY_OTSU_MODE <- "burnable_only"

# Umbral minimo de RBR para que un parche Otsu sea considerado.
# Guard contra ruido espectral negativo.
LEGACY_OTSU_THRESHOLD <- 0

# Umbral de referencia para escalar otros umbrales relativos del Otsu.
LEGACY_REFERENCE_OTSU_THRESHOLD <- 100

# Cap del pool otsu_patch_residual. Sin esto, los ~13.000 parches
# candidatos saturarian el training set.
LEGACY_SAMPLE_N <- 2000L


# --- 5.6) Otsu pipeline - thresholds de confianza (2) -----------------------
#
# DEBUG: HANDOFF §N+28 Tier 1. Hardcoded en 0.5.0, expuesto post-refactor.

# Threshold para considerar un parche Otsu como "alta confianza".
OTSU_KEEP_HI <- 0.45

# Threshold para considerar un parche Otsu como "baja confianza" (drop).
OTSU_DROP_LO <- 0.15


# --- 5.7) Toggles del pipeline (4) ------------------------------------------
#
# Permiten saltar etapas del pipeline (util para reruns parciales).
# DEBUG: HANDOFF §N+28 Tier 2. Hardcoded a TRUE en 0.5.0.

DO_POOLS    <- TRUE  # STEP A: pools (burned + unburned)
DO_FOLDS    <- TRUE  # STEP B1-B2: spatial block folds
DO_FEATURES <- TRUE  # STEP B3: extract features per polygon
DO_MODEL    <- TRUE  # STEP C: OOF + final model + scoring + map


# --- 5.8) Control de I/O (3) ------------------------------------------------

# Si TRUE y existen outputs Otsu de un run previo, los reutiliza en
# lugar de regenerarlos. Solo afecta a la rama Otsu del all_sources.
LEGACY_REUSE_EXISTING <- TRUE

# Si TRUE, escribe outputs intermedios del Otsu a disco.
LEGACY_WRITE_OUTPUT <- TRUE

# Si TRUE, activa el log detallado del STEP A4 (mensajes "[unb] ...").
UNB_VERBOSE <- TRUE


# =============================================================================
# 6) DIRECTORIO DE SALIDA DEL TUTORIAL
# =============================================================================
#
# Para no machacar el run real, redirigimos los outputs del tutorial a
# una carpeta paralela SUPERVISED_TUTORIAL/ (no SUPERVISED/).

TUTORIAL_RESULT_DIR <- file.path(
  RESULTS, as.character(YEAR), RESULT_NAME, "SUPERVISED_TUTORIAL", SCENARIO
)
dir.create(TUTORIAL_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 7) VERIFICACION GLOBAL
# =============================================================================
#
# Validamos que los 13 inputs existen ANTES de empezar el pipeline.
# Si algo falta, fallamos en <1 segundo en lugar de 5-10 minutos despues.

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


# =============================================================================
# 8) RESUMEN VISUAL DE LA CONFIGURACION
# =============================================================================

cat("\n========== CONFIGURACION DEL EXPERIMENTO ==========\n")
cat("\n[1] PARAMETROS DEL COMBO\n")
cat("  YEAR          :", YEAR, "\n")
cat("  SCENARIO      :", SCENARIO, "\n")
cat("  CORINE_YEAR   :", CORINE_YEAR, "\n")
cat("  RESULT_NAME   :", RESULT_NAME, "\n")

cat("\n[2] LOS 13 INPUTS (verificados)\n")

cat("\n  A) Year/scenario-dependent (4):\n")
cat("    INTERNAL_DECISIONS :", basename(INTERNAL_DECISIONS), "\n")
cat("    CHANGE_INDEX       :", basename(CHANGE_INDEX), "\n")
cat("    RBR_AUTUMN         :", basename(RBR_AUTUMN), "\n")
cat("    HOTSPOTS           :", basename(HOTSPOTS), "\n")

cat("\n  B) Validacion externa (1, +shp):\n")
cat("    EFFIS_CA_TIF       :", basename(EFFIS_CA_TIF), "\n")
cat("    EFFIS_CA_SHP       :", basename(EFFIS_CA_SHP), "\n")

cat("\n  C) Corine quinquenal (4):\n")
cat("    CORINE_RASTER      :", basename(CORINE_RASTER), "\n")
cat("    CORINE_STRATA      :", basename(CORINE_STRATA), "\n")
cat("    CORINE_LUT         :", basename(CORINE_LUT), "\n")
cat("    BURNEABLE_MASK     :", basename(BURNEABLE_MASK), "\n")

cat("\n  D) Estaticos (4):\n")
cat("    PENINSULA_SHP      :", basename(PENINSULA_SHP), "\n")
cat("    TOPO               :", basename(TOPO), "\n")
cat("    MASK_TIF           :", basename(MASK_TIF), "\n")
cat("    MASK_SHP           :", basename(MASK_SHP), "\n")

cat("\n[3] CONSTANTES METODOLOGICAS\n")

cat("\n  3.1 Politica de muestreo de negativos\n")
cat("    NEGATIVE_POOL_POLICY              :", NEGATIVE_POOL_POLICY, "\n")
cat("    : combina 3 fuentes de negativos en el training pool:\n")
cat("      deterministic_drops + random_background + Otsu_patches.\n")
cat("      La alternativa 'deterministic_direct' usa solo las 2 primeras.\n")

cat("\n  3.2 Guard de tamaño minimo del pool burned\n")
cat("    MIN_BURNED_POOL_N                 :", MIN_BURNED_POOL_N, "\n")
cat("    : combos con menos de N polygonos burned abortan limpiamente.\n")
cat("      Bajar a 3 admite anos con muy pocos fuegos (riesgo de\n")
cat("      modelo poco fiable). Subir a 10 descarta mas anos pre-MODIS.\n")

cat("\n  3.3 Reproducibilidad: seeds aleatorios\n")
cat("    RANDOM_SEED / LEGACY_RANDOM_SEED  :", RANDOM_SEED, "/", LEGACY_RANDOM_SEED, "\n")
cat("    : seeds del muestreo random (background + Otsu).\n")
cat("      Cambiar genera training pools distintos. Las metricas finales\n")
cat("      varian poco (~ +/- 0.01 en F1). Util para tests de robustez.\n")

cat("\n  3.4 Filtro temporal entre anos\n")
cat("    CURRENTYEAR_PREYEAR_OVERLAP_THR   :", CURRENTYEAR_PREYEAR_OVERLAP_THR, "\n")
cat("    : umbral de solape espacial con polygonos del ano anterior.\n")
cat("      Si supera 70%, se considera 'fuego persistente' (sospechoso).\n")
cat("      Subir a 0.85 = filtro mas estricto, mas drops.\n")

cat("    CURRENTYEAR_HOTSPOT_DENSITY_THR   :", CURRENTYEAR_HOTSPOT_DENSITY_THR, "\n")
cat("    : densidad minima de hotspots/ha para validar un polygono.\n")
cat("      Subir = exigir mas evidencia termica. Solo aplica >= 1995.\n")

cat("    CURRENTYEAR_TEMPORAL_PENALTY_FLOOR:", CURRENTYEAR_TEMPORAL_PENALTY_FLOOR, "\n")
cat("    : floor de la penalty score temporal (no baja de aqui).\n")
cat("      Evita que penalty=0 elimine completamente al polygono.\n")

cat("\n  3.5 Otsu pipeline: parametros principales\n")
cat("    LEGACY_OTSU_MODE                  :", LEGACY_OTSU_MODE, "\n")
cat("    : aplica Otsu solo sobre pixeles quemables segun Corine.\n")
cat("      Sin esto el Otsu veria agua/urbano y los umbrales saldrian\n")
cat("      distorsionados.\n")

cat("    LEGACY_OTSU_THRESHOLD             :", LEGACY_OTSU_THRESHOLD, "\n")
cat("    : umbral minimo de RBR para que un parche Otsu se considere.\n")
cat("      Guard contra ruido espectral negativo.\n")

cat("    LEGACY_REFERENCE_OTSU_THRESHOLD   :", LEGACY_REFERENCE_OTSU_THRESHOLD, "\n")
cat("    : umbral de referencia para escalar otros umbrales del Otsu.\n")

cat("    LEGACY_SAMPLE_N                   :", LEGACY_SAMPLE_N, "\n")
cat("    : cap del pool otsu_patch_residual. Sin esto, los ~13.000\n")
cat("      parches candidatos saturarian el training set.\n")

cat("\n  3.6 Otsu pipeline: thresholds de confianza\n")
cat("    OTSU_KEEP_HI / OTSU_DROP_LO       :", OTSU_KEEP_HI, "/", OTSU_DROP_LO, "\n")
cat("    : umbrales de probabilidad para clasificar parches Otsu.\n")
cat("      KEEP_HI=0.45 -> alta confianza burned. DROP_LO=0.15 -> baja.\n")
cat("      Subir KEEP_HI = mas estricto en aceptar como burned.\n")

cat("\n  3.7 Toggles del pipeline\n")
cat("    DO_POOLS / FOLDS / FEATURES / MODEL:",
    DO_POOLS, "/", DO_FOLDS, "/", DO_FEATURES, "/", DO_MODEL, "\n")
cat("    : ejecutar STEP A (pools) / B1-B2 (folds) / B3 (features) /\n")
cat("      C-D (modelo). En 0.5.0 todos hardcoded a TRUE; el refactor\n")
cat("      permitira reruns parciales (ej: solo regenerar el modelo).\n")

cat("\n  3.8 Control de I/O\n")
cat("    LEGACY_REUSE_EXISTING             :", LEGACY_REUSE_EXISTING, "\n")
cat("    : reutiliza outputs Otsu de runs previos si existen.\n")
cat("      FALSE = regenerar siempre desde cero (mas lento, mas robusto).\n")

cat("    LEGACY_WRITE_OUTPUT               :", LEGACY_WRITE_OUTPUT, "\n")
cat("    : escribe outputs intermedios del Otsu a disco. FALSE = solo\n")
cat("      en memoria (no inspeccionable post-hoc).\n")

cat("    UNB_VERBOSE                       :", UNB_VERBOSE, "\n")
cat("    : log detallado del STEP A4 (mensajes '[unb] ...').\n")


cat("\n[4] DIRECTORIO DE OUTPUTS DEL TUTORIAL\n")
cat("    ", TUTORIAL_RESULT_DIR, "\n")

cat("\nBLOQUE 0 completado. ",
    "13 inputs + 18 constantes metodologicas declaradas y verificadas.\n")

