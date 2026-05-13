#' Train the final supervised model and export the final burned-area map
#'
#' Internal wrapper that chains [train_final_model_direct()] (final-model
#' training on the labelled supervised pool) with
#' `score_burnedlike_and_export_final_map()` (scoring of the deterministic
#' universe and export of the final-map GeoPackage).
#'
#' @section Sampling caps and feature space (Phase B experiments):
#' The four `*_to_burned_ratio` arguments plus `feature_whitelist_override`
#' and `feature_weights` are pass-through hooks to
#' [train_final_model_direct()] used by the Phase B experiment matrix.
#' Their defaults reproduce the historical behaviour byte-for-byte: a
#' caller that does not pass them obtains the same final model and final
#' map as before.
#'
#' @param contextual_exclusion_to_burned_ratio Numeric. Cap on the
#'   contextual-exclusion negative pool. Default `0.25`.
#' @param spectral_hard_negative_to_burned_ratio Numeric. Cap on the
#'   spectral-hard-negative pool. Default `1.0`.
#' @param random_to_burned_ratio Numeric. Cap on random-burnable-background
#'   negatives. Default `1.0`.
#' @param otsu_unburned_to_burned_ratio Numeric. Cap on Otsu current-year
#'   unburned-patch negatives. Default `1.0`.
#' @param feature_whitelist_override Character vector or `NULL`. Subset of
#'   `.supervised_feature_cols` for both the OOF and final-model stages.
#'   See [run_oneyear_supervised_pipeline()] for details.
#' @param feature_weights Named numeric vector or `NULL`. Per-feature
#'   weights forwarded to xgboost. See [run_oneyear_supervised_pipeline()]
#'   for details.
#'
#' @keywords internal
#' @noRd
run_train_final_model_and_export_final_map <- function(
    # ------------------- TRAIN (mismos nombres que tú) -------------------
    run_dir,
    qa,
    labelled_gpkg,
    labelled_layer,
    out_dir,
    prefix,
    overwrite = TRUE,
    verbose   = TRUE,

    # ------------------- SCORE/MAP (mismos nombres que tú) ----------------
    result_dir = run_dir,

    qa_labelled_gpkg,
    qa_labelled_layer,

    labelled_features_gpkg,
    labelled_features_layer,

    model_rds  = file.path(result_dir, "07_FINAL_MODEL_V2", paste0(prefix, "_final_model.rds")),
    recipe_rds = file.path(result_dir, "07_FINAL_MODEL_V2", paste0(prefix, "_recipe.rds")),

    unlabeled_gpkg,
    unlabeled_layer,

    out_score_dir = file.path(result_dir, "08_SCORED"),
    out_map_dir   = file.path(result_dir, "09_FINAL_MAP"),

    export_burned_like = TRUE,
    preyear_overlap_threshold = 0.70,
    hotspot_density_threshold = 0.001,
    temporal_penalty_floor = 0.10,

    # ------------------- Phase B sampling caps (pass-through) -------------
    # Defaults reproduce historical behaviour byte-for-byte.
    contextual_exclusion_to_burned_ratio   = 0.25,
    spectral_hard_negative_to_burned_ratio = 1.0,
    random_to_burned_ratio                 = 1.0,
    otsu_unburned_to_burned_ratio          = 1.0,

    # ------------------- Feature space (0.5.0) ----------------------------
    feature_whitelist_override = NULL,
    feature_weights            = NULL
) {

  if (!exists("train_final_model_direct")) stop("No encuentro train_final_model_direct() cargada (source 06D_FUNCTION_TRAIN_FINAL_MODEL_DIRECT.R).")
  if (!exists("score_burnedlike_and_export_final_map")) stop("No encuentro score_burnedlike_and_export_final_map() cargada (source 07C_FUNCTION_FINAL_MAP_FIRE_LEVEL_PROB_ONLY.R).")

  dir.create(out_dir,        recursive = TRUE, showWarnings = FALSE)
  dir.create(out_score_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_map_dir,    recursive = TRUE, showWarnings = FALSE)

  # 1) Entrena modelo final
  m2 <- train_final_model_direct(
    qa = qa,
    labelled_gpkg  = labelled_gpkg,
    labelled_layer = labelled_layer,
    out_dir        = out_dir,
    prefix         = prefix,
    overwrite      = overwrite,
    verbose        = verbose,
    contextual_exclusion_to_burned_ratio   = contextual_exclusion_to_burned_ratio,
    spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
    random_to_burned_ratio                 = random_to_burned_ratio,
    otsu_unburned_to_burned_ratio          = otsu_unburned_to_burned_ratio,
    feature_whitelist_override             = feature_whitelist_override,
    feature_weights                        = feature_weights
  )

  # 2) Si el train devolvió rutas, las priorizo (sin cambiar tu lógica)
  if (is.list(m2) && !is.null(m2$files$model_rds))  model_rds  <- m2$files$model_rds
  if (is.list(m2) && !is.null(m2$files$recipe_rds)) recipe_rds <- m2$files$recipe_rds

  # 3) Score + export mapa final
  out <- score_burnedlike_and_export_final_map(
    result_dir = result_dir,
    prefix     = prefix,

    qa_labelled_gpkg  = qa_labelled_gpkg,
    qa_labelled_layer = qa_labelled_layer,

    labelled_features_gpkg  = labelled_features_gpkg,
    labelled_features_layer = labelled_features_layer,

    model_rds  = model_rds,
    recipe_rds = recipe_rds,

    unlabeled_gpkg  = unlabeled_gpkg,
    unlabeled_layer = unlabeled_layer,

    out_score_dir = out_score_dir,
    out_map_dir   = out_map_dir,

    export_burned_like = export_burned_like,
    preyear_overlap_threshold = preyear_overlap_threshold,
    hotspot_density_threshold = hotspot_density_threshold,
    temporal_penalty_floor = temporal_penalty_floor,
    overwrite = overwrite,
    verbose   = verbose
  )

  invisible(list(
    m2   = m2,
    out  = out,
    files = list(
      model_rds  = model_rds,
      recipe_rds = recipe_rds
    )
  ))
}
