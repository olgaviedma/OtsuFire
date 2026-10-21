# OtsuFire 2.1.0

Update of 2.0.0. Nothing that worked in 2.0.0 changes by default: every new
behaviour is opt-in, and runs that do not use it give identical results.

## New

* `run_split_rescore()` updates a year that has already been processed with
  the date-aware split of mixed candidates. It reads the existing scored
  candidates, splits those that mix post-event acquisition dates, extracts
  the predictors of the sub-patches, scores them with a given model and
  recipe (no retraining), applies the threshold and writes the final map, an
  audit table and a reproducibility recipe. The deterministic stage is not
  repeated. This is the recommended flow for updating existing years.
* `split_patches_by_doy()` separates candidate patches whose pixels were
  observed on dates more than `doy_gap_days` (default 30) apart, so that a
  real fire joined to haze, cloud residue or another event is scored in
  parts. Validated against official perimeters in 2006 (F1 0.44 to 0.67),
  1986, 1992 and 2017 (a year without haze: unchanged within 0.003).
  Sub-patches smaller than about 100 pixels at the edge of real fires can
  lower recall; see the vignette "Updating years with the DOY split".
* `refine_params = list(split_by_doy = TRUE, doy_gap_days = 30)` runs the
  split inside `run_deterministic_pipeline()`, between detection and scoring.
  It is off by default and has not been validated on real data; use
  `run_split_rescore()` to update years.
* Every `run_split_rescore()` run writes `RECIPE.yml` (package commit,
  parameters including `split_by_doy`, `doy_gap_days` and the
  `artifact_hard` weight when given, input paths, model training
  parameters) and `sessionInfo.txt` next to its results.
* The public API is now 32 exported functions.

## Bug fixes

* Out-of-fold diagnostics now fit every fold with the same XGBoost parameters
  as the final model. Previously `run_oof_xgb()` ignored the `params` it was
  given and always used the package defaults, so a `model_params` override
  (`build_supervised_burned_config()`) changed only the final model. Runs with
  the default `model_params` are unaffected and give identical results.

# OtsuFire 2.0.0

## A complete new release

OtsuFire 2.0.0 is a complete new release, rebuilt from scratch. It is not
an update of the 0.1.x line and is not backward compatible with it: code
written for 0.1.x will not run unchanged. OtsuFire 0.1.4 remains available
on CRAN.

The package maps burned areas from annual change-index rasters (such as
RBR or dNBR) through four stages. Its public API is 31 exported functions,
and each stage is driven by a configuration object, so a run can be
reproduced from its resolved configuration.

## The four stages

* **Mosaic.** `change_index_mosaic()` and `mosaic_from_tiles()` build an
  annual change-index mosaic from raster tiles; `clean_raster_file()`
  checks and optionally repairs raster files.
* **Otsu-guided segmentation.** `build_burned_mapping_config()` and
  `run_deterministic_pipeline()` estimate Otsu thresholds per land-cover
  and ecoregion unit, grow seed-supported candidate patches, apply
  rule-based filters and classify every candidate as `keep`, `review` or
  `drop`. The stage can also be run step by step with
  `detect_burned_patches()` and `score_burned_patches()`, and the optional
  `apply_chain_cleaning()` (alias `apply_chain()`) separates candidates
  joined by narrow connections.
  `collect_keep_reference_samples()` and
  `build_keep_pool_from_samples()` build spectral-support references from
  other years.
* **Probabilistic refinement.** `build_supervised_burned_config()` and
  `run_oneyear_supervised_pipeline()` train a supervised XGBoost model on
  labels derived from the Otsu-guided segmentation and score every
  candidate patch with `p_burned`. Each step is also available on its own:
  `validate_supervised_execution()`, `build_supervised_training_pools()`,
  `make_spatial_folds()`, `extract_supervised_features()`,
  `run_oof_diagnostics()`, `train_final_burned_model()`,
  `score_supervised_burned_map()` and `write_thresholded_burned()`.
* **Validation.** `validate_fire_maps()` compares any burned-area map with
  reference fire perimeters, with pixel-based and polygon/area-based
  metrics, optional stratification and a whole-fire temporal-observability
  rule. `check_supervised_consistency()` compares the outputs of the two
  mapping stages of the same run.

## Training pools and their curation

* Burned labels are the `keep` patches. Unburned labels come from up to
  three pools: burnable background, moderately burned-like drop patches
  and, optionally, strongly burned-like drop patches. The first two are
  capped relative to the burned pool; the third is balanced by weight.
  Review patches are never used as training labels.
* `diagnose_training_pools()` assesses how separable the pools are before
  a model is fitted. `export_visual_validation_pools()` and
  `apply_visual_validation()` support a visual review of the training pool.
  `promote_artifact_hard_negatives()`, `subsample_random_background()`,
  `add_pool_shape_features()` and `assemble_era_training_pool()` build and
  combine curated pools, including pools spanning several years.

## Good to know

* `p_burned` is a model score, not a calibrated probability. Choose the
  threshold for the binary map according to your mapping objective.
* Out-of-fold diagnostics measure agreement with labels derived from the
  Otsu-guided segmentation; they are not an independent accuracy
  assessment. Use `validate_fire_maps()` against independent reference
  data for that.
* The Otsu-guided segmentation stage requires the `whitebox` package
  (WhiteboxTools). The training-pool construction of the probabilistic
  refinement stage requires Python with the GDAL bindings.
