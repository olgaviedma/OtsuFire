# =============================================================================
# supervised-threshold.R  (Phase 3A)
#
# Public threshold + write step for the supervised burned map. Packages the
# decision/write half of the manual champion closure `write_thr()` (from
# 00_USAGE/02_SUPERVISED_USAGE/exp_champion_dropsummerlevel.R), WITHOUT the
# model-scoring half (`dscore()`): it takes an ALREADY-scored candidate set
# (`scored_all`, carrying a per-row model probability column) and writes the
# thresholded burned polygons. No model, no retraining, no re-run.
#
# Champion semantics reproduced exactly:
#   keep <- scored[is.finite(p) & p >= threshold, ]
#   keep <- st_make_valid(st_zm(keep, drop=TRUE, what="ZM"))
#   keep <- keep[!st_is_empty(keep), ]
#   -> thresholded_burned.gpkg (layer "thresholded_burned"), WAL/SHM removed.
# Champion threshold = 0.50, score column = "p_burned_model".
# =============================================================================

#' Apply a score threshold and export burned-area polygons
#'
#' @description
#' Select already-scored candidate polygons using an explicit score threshold
#' and write the selected polygons to a GeoPackage.
#'
#' A candidate is retained when its score is finite and greater than or equal
#' to `threshold`. Optional geometry cleaning removes Z/M dimensions, repairs
#' invalid geometries, and removes empty geometries before export.
#'
#' Use this function after supervised scoring to create a thresholded
#' burned-area product. It does not train a model, calculate predictions, or
#' apply temporal filtering.
#'
#' @param scored An `sf` object containing scored candidate polygons, or a
#'   GeoPackage path. When a path is supplied, the function reads the
#'   `scored_all` layer if present; otherwise, it reads the first layer.
#' @param out_dir Character scalar. Output directory, created if needed.
#' @param threshold Numeric scalar. Minimum score required to retain a
#'   candidate. Values equal to the threshold are included. Default: `0.50`.
#' @param score_col Character scalar. Name of the numeric score column used
#'   for selection. Default: `"p_burned_model"`. Set this explicitly when the
#'   input uses another score field, such as `"p_burned"`.
#' @param layer Character scalar. Name of the output GeoPackage layer.
#'   Default: `"thresholded_burned"`. This argument does not select the input
#'   layer.
#' @param filename Character scalar. Output filename within `out_dir`.
#'   Default: `"thresholded_burned.gpkg"`.
#' @param overwrite Logical scalar. When `FALSE`, an existing output file
#'   causes an error. When `TRUE`, the entire output file and any associated
#'   `-wal` and `-shm` sidecar files are removed before writing. Default:
#'   `FALSE`.
#' @param clean_geometry Logical scalar. Whether to remove Z/M dimensions,
#'   repair invalid geometries with `sf::st_make_valid()`, and remove empty
#'   geometries before writing. Default: `TRUE`.
#' @param verbose Logical scalar. Whether to print a brief processing summary.
#'   Default: `TRUE`.
#'
#' @section Threshold selection:
#' The selection rule is `is.finite(score) & score >= threshold`, where
#' `score` is the column identified by `score_col`.
#'
#' Missing and non-finite scores, including `NA`, `NaN`, `Inf`, and `-Inf`,
#' are excluded.
#'
#' The default threshold of `0.50` is an operating value, not an
#' automatically selected optimum. Model scores are not necessarily
#' calibrated probabilities. Choose the threshold using the intended mapping
#' objective and appropriate evaluation results.
#'
#' @section Input layers and score fields:
#' GeoPackages may contain several scored layers with different candidate
#' sets. When `scored` is a path, automatic selection uses `scored_all` if
#' available and otherwise the first layer.
#'
#' To threshold a specific layer, read it with `sf::st_read()` and pass the
#' resulting `sf` object.
#'
#' For outputs from [score_supervised_burned_map()]:
#' * use `final_map` to threshold candidates retained after the current-year
#'   temporal filter;
#' * use `final_map_full` to threshold the complete scored candidate set;
#' * specify the score column present in the selected layer.
#'
#' This function applies only the score threshold and optional geometry
#' cleaning. It does not reproduce temporal exclusions when given the complete
#' candidate set.
#'
#' @section Geometry cleaning:
#' When `clean_geometry = TRUE`, selected candidates are processed in this
#' order:
#' 1. remove Z/M dimensions;
#' 2. repair invalid geometries;
#' 3. remove empty geometries.
#'
#' The CRS and attribute columns, including the selected score column, are
#' retained. Geometry repair may alter geometry structure, and removing empty
#' geometries may reduce the number of exported features.
#'
#' When `clean_geometry = FALSE`, these cleaning operations are skipped.
#'
#' @section Existing outputs:
#' With `overwrite = TRUE`, replacement applies to the whole output
#' GeoPackage, including any other layers it contains. Use a dedicated output
#' filename for the thresholded product.
#'
#' @return An invisibly returned named list containing:
#'
#' | Field | Contents |
#' |---|---|
#' | `path` | Path to the output GeoPackage. |
#' | `layer` | Output layer name. |
#' | `threshold` | Applied score threshold. |
#' | `score_col` | Score column used for selection. |
#' | `n_in` | Number of input candidates. |
#' | `n_kept` | Number of retained output features. |
#' | `n_nonfinite` | Number of input candidates with missing or non-finite scores. |
#' | `n_dropped_empty` | Number of empty geometries removed during cleaning. |
#' | `crs_epsg` | EPSG identifier reported for the output CRS, when available. |
#' | `overwrite` | Applied overwrite setting. |
#'
#' Assign the result to an object to inspect the summary and output path.
#'
#' @seealso [score_supervised_burned_map()], [train_final_burned_model()],
#'   [run_oof_diagnostics()], [validate_fire_maps()],
#'   [run_oneyear_supervised_pipeline()].
#'
#' @examples
#' \dontrun{
#' # Threshold a scored_all layer containing p_burned_model
#' res <- write_thresholded_burned(
#'   scored = "results/scored_all.gpkg",
#'   out_dir = "results/thresholded",
#'   threshold = 0.50,
#'   score_col = "p_burned_model"
#' )
#'
#' res$n_in
#' res$n_kept
#' res$path
#'
#' # Select a specific layer from a probabilistic refinement final-map GeoPackage
#' public_map <- sf::st_read(
#'   "path/to/final_map.gpkg",
#'   layer = "final_map",
#'   quiet = TRUE
#' )
#'
#' # Threshold the public layer using its p_burned field
#' res_public <- write_thresholded_burned(
#'   scored = public_map,
#'   out_dir = "results/thresholded",
#'   threshold = 0.50,
#'   score_col = "p_burned",
#'   filename = "public_thresholded_burned.gpkg"
#' )
#'
#' # Inspect the exported polygons
#' burned <- sf::st_read(
#'   res_public$path,
#'   layer = res_public$layer,
#'   quiet = TRUE
#' )
#'
#' plot(sf::st_geometry(burned))
#' }
#'
#' @export
write_thresholded_burned <- function(scored,
                                     out_dir,
                                     threshold      = 0.50,
                                     score_col      = "p_burned_model",
                                     layer          = "thresholded_burned",
                                     filename       = "thresholded_burned.gpkg",
                                     overwrite      = FALSE,
                                     clean_geometry = TRUE,
                                     verbose        = TRUE) {
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("write_thresholded_burned(): 'threshold' must be a single finite number.", call. = FALSE)
  }
  if (missing(out_dir) || !is.character(out_dir) || length(out_dir) != 1L) {
    stop("write_thresholded_burned(): 'out_dir' must be a single directory path.", call. = FALSE)
  }

  # ---- read scored set (NEVER scores; input must already be scored) ----
  if (is.character(scored) && length(scored) == 1L) {
    if (!file.exists(scored)) stop("write_thresholded_burned(): scored file not found: ", scored, call. = FALSE)
    layers <- sf::st_layers(scored)$name
    rl <- if ("scored_all" %in% layers) "scored_all" else layers[1]
    scored <- sf::st_read(scored, layer = rl, quiet = TRUE)
  }
  if (!inherits(scored, "sf")) {
    stop("write_thresholded_burned(): 'scored' must be a path or an sf object.", call. = FALSE)
  }
  if (!score_col %in% names(scored)) {
    stop("write_thresholded_burned(): score column '", score_col, "' not found. ",
         "Pass a scored set (the supervised model output) or set 'score_col'.", call. = FALSE)
  }

  n_in <- nrow(scored)
  pm   <- suppressWarnings(as.numeric(scored[[score_col]]))
  n_nonfinite <- sum(!is.finite(pm))

  keep <- scored[is.finite(pm) & pm >= threshold, , drop = FALSE]
  n_dropped_empty <- 0L
  if (isTRUE(clean_geometry)) {
    keep <- sf::st_make_valid(sf::st_zm(keep, drop = TRUE, what = "ZM"))
    n_pre <- nrow(keep)
    keep  <- keep[!sf::st_is_empty(keep), , drop = FALSE]
    n_dropped_empty <- n_pre - nrow(keep)
  }

  # ---- write (robust overwrite: clears WAL/SHM sidecars) ----
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, filename)
  if (file.exists(out_path) && !isTRUE(overwrite)) {
    stop("write_thresholded_burned(): output exists and overwrite=FALSE: ", out_path, call. = FALSE)
  }
  for (s in c("", "-wal", "-shm")) if (file.exists(paste0(out_path, s))) unlink(paste0(out_path, s), force = TRUE)
  sf::st_write(keep, out_path, layer = layer, quiet = TRUE)

  crs_epsg <- sf::st_crs(keep)$epsg
  if (isTRUE(verbose)) {
    message(sprintf(
      "[write_thresholded_burned] %s >= %.4f : kept %d / %d (non-finite %d, empty dropped %d) -> %s (EPSG:%s)",
      score_col, threshold, nrow(keep), n_in, n_nonfinite, n_dropped_empty, out_path,
      if (is.na(crs_epsg)) "NA" else crs_epsg))
  }

  invisible(list(
    path            = out_path,
    layer           = layer,
    threshold       = threshold,
    score_col       = score_col,
    n_in            = n_in,
    n_kept          = nrow(keep),
    n_nonfinite     = n_nonfinite,
    n_dropped_empty = n_dropped_empty,
    crs_epsg        = crs_epsg,
    overwrite       = isTRUE(overwrite)
  ))
}
