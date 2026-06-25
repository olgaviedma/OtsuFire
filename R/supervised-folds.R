#' Build spatial block cross-validation folds for the supervised pipeline
#'
#' @description
#' Assigns spatially blocked cross-validation folds to the labelled training
#' pool, so that out-of-fold diagnostics measure spatial transferability rather
#' than memorisation of nearby polygons. This is the folds stage of the
#' supervised pipeline; its `train_with_folds` output then feeds
#' [extract_supervised_features()] and [run_oof_diagnostics()].
#'
#' Run it after the pools stage ([build_supervised_training_pools()]) and before
#' feature extraction. Pass the labelled training pool and the configuration
#' from [build_supervised_burned_config()]. Calling this function directly
#' produces the same fold outputs as the equivalent stage of the full pipeline
#' [run_oneyear_supervised_pipeline()].
#'
#' @section What it does:
#' \enumerate{
#'   \item Builds candidate square block grids over the polygons, trying each
#'     `block_sizes_m` from larger (stricter) to smaller (looser).
#'   \item Auto-tunes the `(block_size, k)` combination against per-fold
#'     acceptance gates (minimum burned units and positive blocks per fold).
#'   \item Assigns repeated block-CV folds, keeping each `split_unit` (a whole
#'     fire or a single polygon) together within one fold.
#'   \item Propagates the unit-level folds back onto the polygon table.
#'   \item Writes the `02_FOLDS` outputs and returns the in-memory objects.
#' }
#' When no candidate passes the acceptance gate, the best non-passing candidate
#' is used, a `_FOLD_FALLBACK.txt` audit file is written, and a loud
#' `warning()` is raised.
#'
#' @section Outputs:
#' Written under the `02_FOLDS` folder:
#' \itemize{
#'   \item `<year>_blocks_<bs>m.gpkg` — the selected block grid.
#'   \item `<year>_folds_<bs>m_k<k>_r<reps>_unit-<unit>.csv` — fold assignments.
#'   \item `<year>_train_with_folds_<bs>m.gpkg` — polygons with `fold_rep*`
#'     columns.
#' }
#'
#' @param train_labelled sf POLYGON layer OR a single GPKG path. The labelled
#'   training pool produced by the pools stage (the `train_labeled` object /
#'   the `train_labeled` layer of `01_POOLS/<year>_<scenario>_pools.gpkg`).
#'   Must carry the `class` column (with `burned`/`unburned` labels) and, for
#'   `split_unit = "fire"`, the `fire_uid` column. When a path is supplied the
#'   `train_labeled` layer is read.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Used to derive `out_dir` (the
#'   `02_FOLDS` folder under `config$output_routes$base`) and the
#'   output-naming year (`config$target_year`).
#' @param split_unit Character scalar. Unit kept together across folds:
#'   `"fire"` (default, recommended for transferability) or `"poly"`.
#' @param block_sizes_m Numeric vector. Candidate block sizes in meters, tried
#'   from stricter (larger) to looser. Default `c(5000, 3000, 2000)`.
#' @param k_candidates Integer vector. Candidate fold counts, each `>= 2`.
#'   Default `c(5L, 4L, 3L)`.
#' @param out_dir Character or `NULL`. Output folder. Defaults to
#'   `config$output_routes$folds_dir` (the `02_FOLDS` folder).
#' @param write_outputs Logical. Whether to write the blocks GPKG, folds CSV
#'   and train-with-folds GPKG. Default `TRUE`.
#' @param overwrite Logical. Whether to clobber existing fold outputs. Default
#'   `FALSE`. When `FALSE` and all three canonical outputs already exist on
#'   disk, the stage reuses them instead of recomputing (the recomputation is
#'   skipped and the existing objects/paths are returned).
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `blocks` — sf grid of the selected block size.
#'     \item `folds_table` — data.frame of polygon-level fold assignments
#'       (the contents written to the folds CSV).
#'     \item `train_with_folds` — sf polygons + `fold_rep*` columns.
#'     \item `blocks_gpkg`, `folds_csv`, `train_with_folds_gpkg` — written
#'       paths (or the canonical target paths when `write_outputs = FALSE`).
#'     \item `selected` — the chosen `(block_size_m, k_folds, ok, ...)`.
#'   }
#'
#' @seealso
#' [build_supervised_burned_config()], [build_supervised_training_pools()],
#' [extract_supervised_features()], [run_oof_diagnostics()],
#' [validate_supervised_execution()], [run_oneyear_supervised_pipeline()]
#'
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L
#' )
#' pools <- build_supervised_training_pools(cfg)
#' folds <- make_spatial_folds(pools$pools_gpkg, config = cfg)
#' folds$train_with_folds_gpkg
#' }
make_spatial_folds <- function(train_labelled,
                               config,
                               split_unit = "fire",
                               block_sizes_m = c(5000, 3000, 2000),
                               k_candidates = c(5L, 4L, 3L),
                               out_dir = NULL,
                               write_outputs = TRUE,
                               overwrite = FALSE) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(train_labelled) || is.null(train_labelled)) {
    stop("'train_labelled' is required.", call. = FALSE)
  }
  if (!(inherits(train_labelled, c("sf", "SpatVector")) ||
        (is.character(train_labelled) && length(train_labelled) == 1L &&
         file.exists(train_labelled)))) {
    stop("'train_labelled' must be sf, SpatVector, or an existing path.",
         call. = FALSE)
  }
  if (missing(config) || is.null(config) ||
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  # split_unit: accept the public "fire"/"poly" spelling and translate to the
  # internal make_block_folds() spelling ("fire"/"polygon").
  if (!is.character(split_unit) || length(split_unit) != 1L ||
      !split_unit %in% c("fire", "poly", "polygon")) {
    stop("'split_unit' must be one of 'fire' or 'poly'.", call. = FALSE)
  }
  split_unit_internal <- if (identical(split_unit, "poly")) "polygon" else split_unit
  if (!is.numeric(block_sizes_m) || length(block_sizes_m) < 1L ||
      any(block_sizes_m <= 0)) {
    stop("'block_sizes_m' must be a positive numeric vector.", call. = FALSE)
  }
  if (!is.numeric(k_candidates) || length(k_candidates) < 1L ||
      any(k_candidates < 2L)) {
    stop("'k_candidates' must be an integer-like vector with values >= 2.",
         call. = FALSE)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L ||
      is.na(write_outputs)) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve out_dir + year from config
  # ---------------------------------------------------------------------------
  if (is.null(out_dir)) {
    out_dir <- config$output_routes$folds_dir
  }
  if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L ||
      !nzchar(out_dir)) {
    stop("Could not resolve 'out_dir' (02_FOLDS) from config.", call. = FALSE)
  }
  target_year <- config$target_year

  # ---------------------------------------------------------------------------
  # 2) Resolve input sf (read the train_labeled layer if a path was supplied)
  # ---------------------------------------------------------------------------
  if (inherits(train_labelled, "SpatVector")) {
    train_labelled_sf <- sf::st_as_sf(train_labelled)
  } else if (is.character(train_labelled)) {
    train_labelled_sf <- sf::read_sf(train_labelled, layer = "train_labeled")
  } else {
    train_labelled_sf <- train_labelled
  }

  # ---------------------------------------------------------------------------
  # 3) overwrite=FALSE short-circuit: reuse existing canonical outputs.
  #    make_block_folds() recomputes + clobbers via delete_layer=TRUE, so the
  #    historical (overwrite=TRUE) path is the recompute. We only skip when all
  #    three canonical outputs already exist on disk.
  # ---------------------------------------------------------------------------
  if (isTRUE(write_outputs) && !isTRUE(overwrite)) {
    existing <- .of_find_existing_fold_outputs(out_dir, target_year)
    if (!is.null(existing)) {
      twf_sf <- sf::read_sf(existing$train_with_folds_gpkg,
                            layer = "train_with_folds")
      blk_sf <- sf::read_sf(existing$blocks_gpkg, layer = "blocks")
      folds_tbl <- utils::read.csv(existing$folds_csv,
                                   stringsAsFactors = FALSE)
      return(list(
        blocks                = blk_sf,
        folds_table           = folds_tbl,
        train_with_folds      = twf_sf,
        blocks_gpkg           = existing$blocks_gpkg,
        folds_csv             = existing$folds_csv,
        train_with_folds_gpkg = existing$train_with_folds_gpkg,
        selected              = NULL,
        reused                = TRUE
      ))
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Build folds (faithful MOVE of the orchestrator STEP B2 call). The
  #    seed_base / n_repeats / acceptance-gate values are pinned here at the
  #    exact historical orchestrator values so folds are byte-identical.
  # ---------------------------------------------------------------------------
  res_folds <- make_block_folds(
    train_labelled_sf = train_labelled_sf,
    fire_id_col       = "fire_uid",
    split_unit        = split_unit_internal,
    block_sizes_m     = block_sizes_m,
    k_candidates      = k_candidates,
    min_burned_units_per_fold = 3,
    min_pos_blocks_per_fold   = 10,
    n_repeats   = 2,
    seed_base   = 42,
    lonlat_action = "transform",
    out_dir     = out_dir,
    target_year = target_year,
    verbose     = TRUE,
    write_train_with_folds_gpkg = isTRUE(write_outputs),
    write_blocks_gpkg = isTRUE(write_outputs),
    write_folds_csv   = isTRUE(write_outputs)
  )

  # ---------------------------------------------------------------------------
  # 5) Assemble the documented return (objects + paths). Resolve any path the
  #    writer did not populate from the selected block size, mirroring the
  #    orchestrator's defensive path recovery.
  # ---------------------------------------------------------------------------
  bs_sel <- tryCatch(res_folds$selected$block_size_m, error = function(e) NULL)

  blocks_gpkg <- res_folds$saved$blocks_gpkg
  folds_csv   <- res_folds$saved$folds_csv
  train_with_folds_gpkg <- res_folds$saved$train_gpkg

  if (is.null(train_with_folds_gpkg) && !is.null(bs_sel)) {
    cand <- file.path(out_dir,
                      sprintf("%d_train_with_folds_%dm.gpkg", target_year, bs_sel))
    if (file.exists(cand)) train_with_folds_gpkg <- cand
  }
  if (is.null(blocks_gpkg) && !is.null(bs_sel)) {
    cand <- file.path(out_dir,
                      sprintf("%d_blocks_%dm.gpkg", target_year, bs_sel))
    if (file.exists(cand)) blocks_gpkg <- cand
  }

  # Polygon-level folds table (the contents written to the folds CSV).
  folds_table <- if (!is.null(folds_csv) && file.exists(folds_csv)) {
    utils::read.csv(folds_csv, stringsAsFactors = FALSE)
  } else {
    sf::st_drop_geometry(res_folds$train_with_folds)
  }

  list(
    blocks                = res_folds$blocks,
    folds_table           = folds_table,
    train_with_folds      = res_folds$train_with_folds,
    blocks_gpkg           = blocks_gpkg,
    folds_csv             = folds_csv,
    train_with_folds_gpkg = train_with_folds_gpkg,
    selected              = res_folds$selected,
    diagnostics           = res_folds$diagnostics,
    reused                = FALSE
  )
}

#' Locate an existing complete set of canonical fold outputs in 02_FOLDS.
#'
#' Returns a list with `blocks_gpkg`, `folds_csv`, `train_with_folds_gpkg`
#' when all three exist for `target_year`, else `NULL`.
#'
#' @keywords internal
#' @noRd
.of_find_existing_fold_outputs <- function(out_dir, target_year) {
  if (is.null(out_dir) || !dir.exists(out_dir)) return(NULL)
  twf <- list.files(
    out_dir,
    pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$", target_year),
    full.names = TRUE
  )
  blk <- list.files(
    out_dir,
    pattern = sprintf("^%d_blocks_\\d+m\\.gpkg$", target_year),
    full.names = TRUE
  )
  csv <- list.files(
    out_dir,
    pattern = sprintf("^%d_folds_\\d+m_k\\d+_r\\d+_unit-.+\\.csv$", target_year),
    full.names = TRUE
  )
  if (length(twf) >= 1L && length(blk) >= 1L && length(csv) >= 1L) {
    return(list(
      blocks_gpkg           = blk[[1L]],
      folds_csv             = csv[[1L]],
      train_with_folds_gpkg = twf[[1L]]
    ))
  }
  NULL
}
