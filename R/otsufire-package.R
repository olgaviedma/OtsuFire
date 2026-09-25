#' OtsuFire: a self-labelling framework for wall-to-wall burned-area mapping across sensors and decades
#'
#' @description
#' Reproducible multi-year burned-area mapping from change-index rasters.
#' Provides four layers that can be used together or independently:
#' a mosaic stage, an Otsu-guided segmentation workflow, a
#' one-year probabilistic refinement workflow, and a shared
#' workflow-independent validation utility.
#'
#' The stage names follow the accompanying paper. In function and folder
#' names, the Otsu-guided segmentation stage appears as "deterministic"
#' (`run_deterministic_pipeline()`, `DETERMINISTIC/`) and the probabilistic
#' refinement stage, which trains a supervised XGBoost model, as "supervised"
#' (`run_oneyear_supervised_pipeline()`, `SUPERVISED/`).
#'
#' @section Workflow entry points:
#' - [change_index_mosaic()] — annual mosaic of a change index.
#' - [build_burned_mapping_config()] + [run_deterministic_pipeline()] —
#'   Otsu-guided segmentation workflow.
#' - [build_supervised_burned_config()] +
#'   [run_oneyear_supervised_pipeline()] — one-year probabilistic refinement workflow.
#' - [validate_fire_maps()] — shared workflow-independent validation
#'   utility; usable inside the Otsu-guided segmentation pipeline and directly on
#'   any compatible burned-area prediction layer (including probabilistic refinement
#'   outputs after thresholding `p_burned`).
#'
#' @section Modular utilities:
#' [detect_burned_patches()], [score_burned_patches()],
#' [build_supervised_training_pools()], [make_spatial_folds()],
#' [extract_supervised_features()], [run_oof_diagnostics()],
#' [train_final_burned_model()], [score_supervised_burned_map()],
#' [check_supervised_consistency()].
#'
#' @docType package
#' @name OtsuFire-package
#' @aliases OtsuFire
#' @keywords internal
"_PACKAGE"
