# =============================================================================
# OtsuFire 2.0.0 tutorial - STAGE 3: SUPERVISED
# =============================================================================
#
# GOAL
#   Learn from the deterministic decisions, then score every candidate patch
#   with a burned probability-like value, p_burned.
#
# WHY DO THIS AT ALL, IF STAGE 2 ALREADY DECIDED
#   Stage 2 applies fixed rules. It is transparent but rigid: a patch either
#   clears a threshold or it does not. Stage 3 learns which COMBINATIONS of
#   severity, persistence, terrain, land cover and thermal evidence actually
#   mark a fire, and it returns a continuous score instead of a hard label.
#   That score is what stage 4 thresholds, per year, against the reference.
#
# WHERE THE TRAINING DATA COMES FROM
#   Positives: the "keep" polygons from stage 2.
#   Negatives: NOT "everything else". Exactly two buckets, both inside the
#     burnable mask:
#       random - random background cells
#       otsu   - Otsu residual patches from this year
#     Review and keep-Otsu patches are excluded and logged. A row that is
#     neither a positive nor one of those two buckets is an error, never a
#     silent negative.
#
# RUN 00_SETUP.R AND 02_DETERMINISTIC.R FIRST.
# =============================================================================

source("00_SETUP.R")


# -----------------------------------------------------------------------------
# 1) Point at the deterministic output
# -----------------------------------------------------------------------------
# Either reuse the object from stage 2:
#   INTERNAL_DECISIONS <- det$result_paths$internal_decisions
# or name the file directly:

INTERNAL_DECISIONS <- file.path(
  OUTPUT_DIR, as.character(TARGET_YEAR), RUN_NAME,
  "DETERMINISTIC", "05_DECISIONS", "internal_decisions.gpkg"
)

stopifnot(file.exists(INTERNAL_DECISIONS))


# -----------------------------------------------------------------------------
# 2) The negative pool
# -----------------------------------------------------------------------------
# This is the block that most changes the model, so it is worth understanding.
#
#   random$n_cells       how many background cells to draw
#   random$rbr_quantile  draw them at this severity percentile. 0.50 means the
#                        background looks "average", not suspiciously dark: a
#                        model trained against only dark background learns
#                        "bright = fire", which is exactly the mistake we want
#                        to avoid.
#   otsu$candidate_threshold / reference_threshold
#                        severity limits for the Otsu residual bucket.
#   caps                 how many negatives per positive, per bucket. c(1, 1)
#                        means at most one random and one Otsu negative per
#                        burned polygon. The cap is
#                        ceiling(n_burned * ratio); Inf disables a cap.
#
# The caps are the single source of truth and are mirrored into
# cfg$train_control$caps. Both the OOF folds and the final model route through
# the same capping helper, so they cannot diverge.

negative_pool_params <- list(
  random = list(n_cells = 1500L, rbr_quantile = 0.50),
  otsu   = list(candidate_threshold = 0, reference_threshold = 100),
  caps   = c(random = 1.0, otsu = 1.0)
)


# -----------------------------------------------------------------------------
# 3) Build the configuration
# -----------------------------------------------------------------------------
# run_label only NAMES the run. It has no methodological effect: it decides
# which folder the outputs land in, nothing else.

sup_cfg <- build_supervised_burned_config(
  run_label            = "tutorial",
  internal_decisions   = INTERNAL_DECISIONS,
  change_index         = CHANGE_INDEX,           # summer  -> rbr features
  delayed_change_index = DELAYED_CHANGE_INDEX,   # autumn  -> rbr_aw features
  hotspots             = HOTSPOTS,               # NULL before 2000
  reference_burned_map = REFERENCE_BURNED_MAP,
  target_year          = TARGET_YEAR,
  output_dir           = OUTPUT_DIR,
  run_name             = RUN_NAME,

  # TRUE: outputs hang directly off output_dir, instead of the deeper
  # <year>/<run>/SUPERVISED/<label>/ tree. Easier to find while learning.
  flat_output_routes   = TRUE,

  peninsula_shapefile  = PENINSULA_SHP,
  topo                 = TOPO,
  corine_raster        = VEGETATION_MAP,
  burnable_mask        = BURNABLE_MASK,

  # Training control -----------------------------------------------------
  nrounds_max          = 2000L,   # ceiling on boosting rounds
  early_stop           = 50L,     # stop after this many rounds with no gain
  val_frac             = 0.2,     # inner validation split, for early stopping

  # Seeds. Same seeds + same inputs = same model.
  oof_seed_base        = 42L,
  final_sampling_seed  = 42L,
  final_seed           = 42L,
  random_seed          = 42L,

  impute_numeric        = "median",
  impute_factor_missing = "MISSING",

  negative_pool_params = negative_pool_params,

  runtime_options = list(
    reuse_existing = TRUE,   # reuse heavy intermediate products if present
    write_outputs  = TRUE,
    verbose        = TRUE
  ),

  options = list(
    data_base      = DATA_BASE,
    composite_base = COMPOSITE_DIR,
    result_name    = RUN_NAME
  )
)


# -----------------------------------------------------------------------------
# 4) Run the whole supervised stage
# -----------------------------------------------------------------------------
# This is the long one: pools, spatial folds, feature extraction, out-of-fold
# diagnostics, final model, scoring. Budget an hour or more for a full year.

sup <- run_oneyear_supervised_pipeline(
  config          = sup_cfg,
  run_consistency = TRUE,    # also cross-check against the deterministic output
  overwrite       = FALSE
)

# The products, in the order they are made:
sup$pools_gpkg              # 01 burned + unburned training pools
sup$train_with_folds_gpkg   # 02 spatial folds
sup$features_geometry_gpkg  # 03 the feature table
sup$oof_agg_csv             # 05 out-of-fold diagnostics
sup$final_model_rds         # 07 the trained model
sup$final_map_gpkg          # 09 the scored map  <- stage 4 reads this
sup$burned_like_gpkg        # unlabelled patches the model scores as burned
sup$timing_csv


# -----------------------------------------------------------------------------
# 5) Read the scored map
# -----------------------------------------------------------------------------
fm <- sf::st_read(sup$final_map_gpkg, layer = "final_map_full", quiet = TRUE)

# The three score columns, and the difference between them:
#
#   p_burned_model  the model's score for every patch.
#   p_burned_oof    for TRAINING patches only: the score the model gave when
#                   that patch was held out. This is the honest one for them.
#   p_burned_eval   p_burned_oof where it exists, p_burned_model elsewhere.
#                   Use this when you want one column for everything.
summary(fm$p_burned_model)

hist(fm$p_burned_model, breaks = 50,
     main = "p_burned_model", xlab = "score")


# -----------------------------------------------------------------------------
# 6) READ THIS BEFORE USING THE SCORES
# -----------------------------------------------------------------------------
# p_burned is a SCORE, not a calibrated probability. A value of 0.8 does not
# mean "80% chance of being burned". It orders patches sensibly under the
# training distribution, and that is all it promises.
#
# Two consequences:
#   - 0.5 is not a meaningful default cut. Stage 4 picks the cut from the data.
#   - Scores are not comparable across years unless you calibrate them.
#
# The out-of-fold metrics in sup$oof_agg_csv are also NOT ground truth: the
# labels they score against are stage 2's decisions, which are rules. Only
# stage 4 compares against an independent reference.


# -----------------------------------------------------------------------------
# 7) Step by step, when you want to see the middle
# -----------------------------------------------------------------------------
# The wrapper above chains these six functions. Call them one at a time when
# you want to inspect or change something between stages - most often the
# training pool, before any model is fitted.

if (FALSE) {
  pools <- build_supervised_training_pools(
    config                  = sup_cfg,
    deterministic_decisions = INTERNAL_DECISIONS,
    write_outputs           = TRUE
  )

  # Are positives and negatives actually separable? Run this BEFORE fitting.
  # It never trains the final model; it just measures separability and writes
  # plots and tables. If the two groups overlap almost completely, no amount
  # of tuning will save the run.
  diag <- diagnose_training_pools(
    pools      = pools,
    output_dir = file.path(OUTPUT_DIR, "POOL_DIAGNOSTICS"),
    prefix     = sprintf("tutorial_%d", TARGET_YEAR),
    make_plots = TRUE
  )

  # Careful with the spelling: the ARGUMENT is train_labelled (two l's), the
  # FIELD the pools object returns is train_labeled (one l).
  folds <- make_spatial_folds(
    train_labelled = pools$train_labeled,
    config         = sup_cfg,
    write_outputs  = TRUE
  )

  feats <- extract_supervised_features(
    train_with_folds = folds$train_with_folds,
    scoring_pool     = pools$scoring_pool,
    config           = sup_cfg,
    write_outputs    = TRUE
  )

  oof <- run_oof_diagnostics(
    train_features   = feats$train_features,
    scoring_features = feats$scoring_features,
    config           = sup_cfg
  )

  final <- train_final_burned_model(
    train_features = feats$train_features,
    config         = sup_cfg,
    oof_agg        = oof$oof_agg
  )

  scored <- score_supervised_burned_map(
    scoring_features = feats$scoring_features,
    model            = final$model,
    recipe           = final$recipe,
    config           = sup_cfg,
    oof_summary      = oof$labeled_oof_summary   # note: not oof$oof_summary
  )
}


# -----------------------------------------------------------------------------
# WHY THE FOLDS ARE SPATIAL
# -----------------------------------------------------------------------------
# Fires are spatially clustered. With random folds, a patch and its neighbour
# land in train and test respectively, the model recognises the neighbourhood
# rather than the fire, and the metrics come out far too good. Block folds put
# whole spatial blocks on one side of the split, which removes that leak.
#
# NEXT: 04_VALIDATION.R
# =============================================================================
