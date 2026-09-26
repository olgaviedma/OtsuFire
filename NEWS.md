# OtsuFire 2.0.0

## A complete new release

OtsuFire 2.0.0 is a complete new release, rebuilt from scratch. It is not
an update of the 0.1.x line and is not backward compatible with it: code
written for 0.1.x will not run unchanged. OtsuFire 0.1.4 remains available
on CRAN.

The package maps burned areas from annual change-index rasters (such as
RBR or dNBR) through four stages. Its public API is 30 exported functions,
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
  joined by narrow connections. `collect_keep_reference_samples()` and
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
