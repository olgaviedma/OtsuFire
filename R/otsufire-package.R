#' OtsuFire: a self-labelling framework for wall-to-wall burned-area mapping across sensors and decades
#'
#' @description
#' \subsection{Package, workflow book and data}{
#' The OtsuFire package and its source code are available on
#' [GitHub](https://github.com/olgaviedma/OtsuFire).
#'
#' For step-by-step guidance on implementing and using OtsuFire, see the
#' [OtsuFire workflow book](https://github.com/nataliaquintero/OtsuFire-book).
#' The companion repository also includes the Google Earth Engine scripts
#' used for image selection, masking and compositing.
#'
#' Burned-area maps for the full time series are available on
#' [Zenodo](https://doi.org/10.5281/zenodo.22893763).
#'
#' Tables S13–S15 of the accompanying manuscript's supplementary material
#' describe the functions included in the package.
#' }
#'
#' Reproducible multi-year burned-area mapping from change-index rasters.
#' Provides four layers that can be used together or independently:
#' a mosaic stage, an Otsu-guided segmentation workflow, a
#' one-year probabilistic refinement workflow, and a shared
#' workflow-independent validation utility.
#'
#' Stages are named by what they do; function and folder names describe how
#' they work. Otsu-guided segmentation is \emph{deterministic} (explicit rules,
#' with data-adaptive Otsu thresholds; no classifier is fitted):
#' [run_deterministic_pipeline()], `DETERMINISTIC/`. Probabilistic
#' refinement is \emph{supervised} (a trained XGBoost model):
#' [run_oneyear_supervised_pipeline()], `SUPERVISED/`.
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
