#' OtsuFire: a self-labelling framework for wall-to-wall burned-area mapping across sensors and decades
#'
#' @description
#' Reproducible multi-year burned-area mapping from change-index rasters.
#' Provides four layers that can be used together or independently:
#' a mosaic stage, a deterministic Otsu/grow segmentation workflow, a
#' supervised one-year probabilistic workflow, and a shared
#' workflow-independent validation utility.
#'
#' @section Workflow entry points:
#' - [change_index_mosaic()] — annual mosaic of a change index.
#' - [build_burned_mapping_config()] + [run_deterministic_pipeline()] —
#'   deterministic workflow.
#' - [build_supervised_burned_config()] +
#'   [run_oneyear_supervised_pipeline()] — one-year supervised workflow.
#' - [validate_fire_maps()] — shared workflow-independent validation
#'   utility; usable inside the deterministic pipeline and directly on
#'   any compatible burned-area prediction layer (including supervised
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
