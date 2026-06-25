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

#' Threshold a scored supervised candidate set and write thresholded_burned.gpkg
#'
#' @description
#' Takes candidates already scored by the supervised model (a \code{scored_all}
#' set carrying a per-row probability column) and writes the burned polygons
#' that pass an explicit probability threshold. Use it as the final, inspectable
#' threshold-and-write step of the supervised pipeline; it does only that step
#' (no scoring, no model, no retraining).
#'
#' The keep rule is \code{is.finite(score) & score >= threshold}, followed by
#' dropping Z/M dimensions, \code{\link[sf]{st_make_valid}}, and removal of empty
#' geometries. CRS and all attribute columns (including the score column) are
#' preserved.
#'
#' @param scored A scored candidate set: either a path to a GeoPackage (the
#'   \code{scored_all} layer is used when present, otherwise the first layer),
#'   or an in-memory \code{sf}.
#' @param out_dir Directory to write into (created if needed).
#' @param threshold Numeric operating threshold on \code{score_col}. Default
#'   \code{0.50}.
#' @param score_col Name of the model-probability column. Default
#'   \code{"p_burned_model"} (written by the supervised scorer).
#' @param layer Output layer name. Default \code{"thresholded_burned"}.
#' @param filename Output file name. Default \code{"thresholded_burned.gpkg"}.
#' @param overwrite Logical. If \code{FALSE} (default) and the output file
#'   already exists, the function stops; if \code{TRUE} the file (and any
#'   \code{-wal}/\code{-shm} sidecars) is removed first.
#' @param clean_geometry Logical (default \code{TRUE}). Drop Z/M dimensions,
#'   make geometries valid, and drop empties before writing.
#' @param verbose Logical (default \code{TRUE}). Emit a one-line summary message.
#' @return Invisibly, a list with \code{path}, \code{layer}, \code{threshold},
#'   \code{score_col}, \code{n_in}, \code{n_kept}, \code{n_nonfinite},
#'   \code{n_dropped_empty}, \code{crs_epsg} and \code{overwrite}.
#'
#' @seealso
#' \code{\link{score_supervised_burned_map}},
#' \code{\link{train_final_burned_model}},
#' \code{\link{run_oneyear_supervised_pipeline}}
#'
#' @examples
#' \dontrun{
#' # Take a scored candidate set and write the thresholded burned polygons.
#' res <- write_thresholded_burned(
#'   ".../balanced/scored_all.gpkg",
#'   out_dir   = ".../my_run/balanced",
#'   threshold = 0.50,
#'   overwrite = TRUE)
#' res$n_in; res$n_kept; res$path
#' }
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
