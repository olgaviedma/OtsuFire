# =============================================================================
# OtsuFire 2.0.0 tutorial - STAGE 4: VALIDATION
# =============================================================================
#
# GOAL
#   Turn the continuous scores from stage 3 into a burned-area map, and measure
#   how good it is against an INDEPENDENT reference (EFFIS).
#
# THE ONE THING TO GET RIGHT
#   There are two different validations in this workflow and they answer
#   different questions. Do not compare their numbers.
#
#     INTERNAL (out-of-fold, stage 3, 05_OOF)
#       Question: does the model reproduce the DETERMINISTIC DECISIONS?
#       Reference: stage 2's keep/drop labels - which are rules, not truth.
#       Use it to detect leakage and overfitting.
#
#     EXTERNAL (this file, validate_fire_maps)
#       Question: does the map match REAL FIRES?
#       Reference: EFFIS perimeters, independent of the whole workflow.
#       This is the number that goes in a paper.
#
#   An OOF AUC of 0.95 and an external F1 of 0.40 are not a contradiction. The
#   first says the model learned stage 2's rules; the second says those rules
#   miss real fires. Both can be true at once.
#
# RUN THE PREVIOUS FILES FIRST.
# =============================================================================

source("00_SETUP.R")


# -----------------------------------------------------------------------------
# 1) Inputs
# -----------------------------------------------------------------------------
# The scored map from stage 3:
FINAL_MAP <- file.path(OUTPUT_DIR, "09_FINAL_MAP",
                       sprintf("%d_final_map.gpkg", TARGET_YEAR))

# Where validation products go:
VALIDATION_DIR <- file.path(OUTPUT_DIR, "VALIDATION")
dir.create(VALIDATION_DIR, recursive = TRUE, showWarnings = FALSE)

# Optional stratification: reports metrics per land-cover stratum as well as
# overall. Set both to NULL to skip.
STRATA_RASTER <- file.path(DATA_BASE, "Corine_Masks", "STRATA",
                           sprintf("strata_CLC_%d_res30.tif", CORINE_YEAR))
STRATA_LUT    <- file.path(DATA_BASE, "Corine_Masks", "LUT",
                           "lut_full_strata8_v1.csv")

fm <- sf::st_read(FINAL_MAP, layer = "final_map_full", quiet = TRUE)


# -----------------------------------------------------------------------------
# 2) Choose the threshold from the data, not by habit
# -----------------------------------------------------------------------------
# p_burned_model is a score, so the best cut differs per year. Picking 0.5
# because it looks round is the most common mistake at this stage.
#
# The idea below: rasterise the CONTINUOUS score once onto the reference grid,
# then count TP/FP/FN for every candidate threshold at the same time. Running
# validate_fire_maps 40 times instead would take hours.

GRID <- seq(0.30, 0.70, by = 0.01)

# The reference, rasterised onto a common grid. Simple version: rasterise the
# EFFIS polygons directly. (validate_fire_maps caches its own reference grid;
# if you have already run it once, reuse that cached mask instead.)
ref_sf <- sf::st_read(REFERENCE_BURNED_MAP, quiet = TRUE)
ref_sf <- sf::st_make_valid(ref_sf)
if (sf::st_crs(ref_sf)$epsg != 3035) ref_sf <- sf::st_transform(ref_sf, 3035)

# ~90 m is fine for CHOOSING a threshold; the official metrics below run at the
# reference's own resolution.
template <- terra::rast(terra::vect(ref_sf), resolution = 90)
ref_r <- terra::rasterize(terra::vect(ref_sf), template, field = 1, background = 0)

pred_sf <- sf::st_make_valid(fm[is.finite(fm$p_burned_model), ])
if (sf::st_crs(pred_sf)$epsg != 3035) pred_sf <- sf::st_transform(pred_sf, 3035)

score_r <- terra::rasterize(terra::vect(pred_sf), template,
                            field = "p_burned_model", fun = "max", background = NA)

# Discretise to hundredths and cross-tabulate once.
score_bin <- terra::app(score_r, function(x) round(x * 100))
ct  <- terra::crosstab(c(score_bin, ref_r), useNA = FALSE)
ctm <- as.matrix(ct)
bins <- as.numeric(rownames(ctm))
col_burned   <- if ("1" %in% colnames(ctm)) ctm[, "1"] else rep(0, nrow(ctm))
col_unburned <- if ("0" %in% colnames(ctm)) ctm[, "0"] else rep(0, nrow(ctm))
n_ref_burned <- sum(terra::freq(ref_r)$count[terra::freq(ref_r)$value == 1], na.rm = TRUE)

sweep <- do.call(rbind, lapply(GRID, function(thr) {
  sel <- bins >= (thr * 100 - 1e-9)
  TP <- sum(col_burned[sel]); FP <- sum(col_unburned[sel]); FN <- n_ref_burned - TP
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA_real_
  data.frame(threshold = thr, TP = TP, FP = FP, FN = FN,
             precision = round(precision, 4), recall = round(recall, 4),
             F1 = round(f1, 4))
}))

best_threshold <- sweep$threshold[which.max(sweep$F1)]

print(sweep[which.max(sweep$F1), ])
cat(sprintf("Chosen threshold for %d: %.2f\n", TARGET_YEAR, best_threshold))

plot(sweep$threshold, sweep$F1, type = "l",
     xlab = "threshold", ylab = "pixel F1",
     main = sprintf("Threshold sweep %d", TARGET_YEAR))
abline(v = best_threshold, lty = 2)

write.csv(sweep, file.path(VALIDATION_DIR,
          sprintf("threshold_sweep_%d.csv", TARGET_YEAR)), row.names = FALSE)


# -----------------------------------------------------------------------------
# 3) Cut the map at that threshold
# -----------------------------------------------------------------------------
burned <- fm[is.finite(fm$p_burned_model) & fm$p_burned_model >= best_threshold, ]
burned <- sf::st_make_valid(burned)
burned <- burned[!sf::st_is_empty(burned), ]

thresholded_path <- file.path(VALIDATION_DIR, "thresholded_burned.gpkg")
if (file.exists(thresholded_path)) unlink(thresholded_path, force = TRUE)
sf::st_write(burned, thresholded_path, layer = "thresholded_burned", quiet = TRUE)

cat(sprintf("%d patches kept, %.0f ha\n",
            nrow(burned), sum(as.numeric(sf::st_area(burned))) / 10000))


# -----------------------------------------------------------------------------
# 4) The official validation
# -----------------------------------------------------------------------------
# Note the argument names: input_shapefile / ref_shapefile, not
# predicted / reference.
#
#   mask_shapefile   the domain metrics are computed inside. For years before
#                    2000 this is usually a smaller mask, because reference
#                    data does not cover the whole peninsula yet.
#   burnable_raster  restricts the comparison to pixels that can burn, so bare
#                    rock is not counted as a correct rejection.
#   metrics_type     "all" gives both pixel-based and polygon-based metrics.
#
#   observability_raster + observability_mode + ref_obs_doy_col: TEMPORAL
#   observability. You cannot miss a fire your imagery never saw. With
#   "wholefire_fraction" a reference fire is evaluated only if at least 75 %
#   of its burnable pixels carry a composite DOY at or after its
#   obs_required_doy (the authoritative Reference V2 date; start_doy/end_doy
#   are not consulted). A fire that fails is removed from the reference AND
#   from the evaluation domain, so it counts neither as omission nor, through
#   a correct detection over it, as commission. A fire with no usable date is
#   reported as UNDETERMINED_OBS_DATE and still evaluated. The decision is per
#   fire: nothing is trimmed at pixel level, and this is not a cloud mask.

val <- validate_fire_maps(
  input_shapefile = thresholded_path,
  ref_shapefile   = REFERENCE_BURNED_MAP,
  mask_shapefile  = VALIDATION_MASK,
  burnable_raster = BURNABLE_MASK,
  year_target     = TARGET_YEAR,
  validation_dir  = VALIDATION_DIR,

  metrics_type      = "all",
  dissolve_ref_by   = "id",
  dissolve_input_by = NULL,

  strata_raster = STRATA_RASTER,
  strata_lut    = STRATA_LUT,

  observability_raster = CHANGE_INDEX,          # its "doy" band is used
  observability_mode   = "wholefire_fraction",  # >= 75 % of the fire observed
  ref_obs_doy_col      = "obs_required_doy",    # Reference V2 authoritative date

  # TRUE forces the cached reference to be rebuilt. The cache is content-aware
  # and normally rebuilds itself when anything relevant changes.
  force_reprocess_ref  = FALSE,
  force_reprocess_pred = TRUE
)

print(val$metrics)

write.csv(val$metrics,
          file.path(VALIDATION_DIR, sprintf("metrics_%d.csv", TARGET_YEAR)),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 5) How to read those numbers
# -----------------------------------------------------------------------------
#   omission   real fire the map missed.        Lower is better.
#   commission map says burned, reference says no. Lower is better.
#   F1         the balance of the two.          Higher is better.
#   IoU        overlap area / union area. Stricter than F1: it punishes a
#              perimeter that is roughly right but too fat or too thin.
#
# Commission is not automatically an error. The reference has its own
# omissions - EFFIS misses small fires - so a patch counted as commission may
# be a real fire the reference never recorded. Look at a few before concluding
# the model is over-detecting.


# -----------------------------------------------------------------------------
# 6) Consistency between stages 2 and 3
# -----------------------------------------------------------------------------
# A different question again: do the deterministic and supervised outputs of
# the same run agree with each other? Useful when a run looks odd and you need
# to know which stage introduced the problem.

if (FALSE) {
  cons <- check_supervised_consistency(
    deterministic_decisions = INTERNAL_DECISIONS,
    final_map               = FINAL_MAP,
    config                  = sup_cfg,
    target_year             = TARGET_YEAR
  )
}
# This also runs automatically inside
# run_oneyear_supervised_pipeline(..., run_consistency = TRUE).


# -----------------------------------------------------------------------------
# THE WHOLE WORKFLOW, IN FOUR LINES
# -----------------------------------------------------------------------------
#   1 MOSAIC         tiles            -> one change-index raster per year
#   2 DETERMINISTIC  change index     -> candidate patches, keep/review/drop
#   3 SUPERVISED     those decisions  -> a model, and p_burned for every patch
#   4 VALIDATION     p_burned         -> a burned map, and honest metrics
#
# Repeating this for a new year means changing TARGET_YEAR in 00_SETUP.R and
# running the four files again. Choose the threshold per year: it is not a
# constant.
# =============================================================================
