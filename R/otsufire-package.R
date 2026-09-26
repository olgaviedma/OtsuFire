#' OtsuFire: a self-labelling framework for wall-to-wall burned-area mapping across sensors and decades
#'
#' @description
#' `OtsuFire` is an open-source R package for long-term, wall-to-wall
#' burned-area reconstruction. It implements a self-labelling, weakly
#' supervised framework that generates training labels without requiring
#' prior knowledge of fire locations or reliable target-year reference data.
#'
#' The framework separates **supervision generation** from **final
#' classification**: Otsu-guided segmentation and physically interpretable
#' rules provide the evidence used to construct training pools, which
#' then train a probabilistic patch-level classifier.
#'
#' The package provides four composable stages:
#'
#' 1. **Mosaic preparation:** assemble image tiles into an annual
#'    change-index raster.
#'
#' 2. **Otsu-guided segmentation and rule-based filtering:** delineate
#'    candidate patches using Otsu thresholds and seed-and-grow
#'    segmentation, then apply spectral, spatial and contextual filters
#'    to assign `keep`, `review` or `drop` decisions.
#'
#' 3. **Probabilistic refinement:** construct burned and unburned training
#'    pools, train an XGBoost classifier, and estimate `p_burned` for
#'    candidate patches to produce the final burned-area map.
#'
#' 4. **Independent validation:** compare mapped burned areas with an
#'    external reference. This utility can also evaluate maps generated
#'    outside OtsuFire.
#'
#' Wall-to-wall coverage describes the mapping domain: the whole study area
#' is searched for candidate patches. Probabilistic scoring then operates on
#' those candidate patches.
#'
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
#' Zenodo (\doi{10.5281/zenodo.22893763}).
#' }
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
