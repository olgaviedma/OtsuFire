#' Assign spatial cross-validation folds to training polygons
#'
#' @description
#' Assign repeated spatial block cross-validation folds to the labelled
#' training pool.
#'
#' The function evaluates candidate block sizes and fold counts, selects a
#' partition, and assigns folds while keeping observations from the same
#' splitting unit together. These assignments support the evaluation of model
#' performance on spatially held-out observations.
#'
#' Run this stage after [build_supervised_training_pools()] and before
#' feature extraction and out-of-fold diagnostics.
#'
#' @param train_labelled `sf` polygon object or GeoPackage path containing the
#'   labelled training pool. Must contain a `class` column with `"burned"` and
#'   `"unburned"` labels. With `split_unit = "fire"`, a `fire_uid` column is
#'   also required. When a path is supplied, the function reads the
#'   `train_labeled` layer.
#' @param config Required object of class `otsufire_supervised_burned_config`,
#'   created with [build_supervised_burned_config()]. Supplies the target year
#'   and default output location.
#' @param split_unit Character scalar. Unit kept together within each
#'   cross-validation repetition: `"fire"` groups polygons by `fire_uid`;
#'   `"poly"` treats each polygon as a separate unit. Default: `"fire"`.
#' @param block_sizes_m Numeric vector of candidate square-block sizes in
#'   metres. Larger sizes are tried before smaller sizes. Default:
#'   `c(5000, 3000, 2000)`.
#' @param k_candidates Integer vector of candidate fold counts. Each value
#'   must be at least 2. Default: `c(5L, 4L, 3L)`.
#' @param out_dir Character scalar or `NULL`. Output directory. When `NULL`,
#'   uses `config$output_routes$folds_dir`, normally the run's `02_FOLDS`
#'   folder.
#' @param write_outputs Logical scalar. Whether to write the block grid,
#'   fold-assignment table, and training layer with folds. Default: `TRUE`.
#' @param overwrite Logical scalar. Whether existing fold outputs may be
#'   replaced. When `FALSE` and all three expected output files exist, they
#'   are reused without recomputing the partition. Default: `FALSE`.
#'
#' @section Input labels and identifiers:
#' The input is the labelled training set returned as `train_labeled` by
#' [build_supervised_training_pools()] or stored in the `train_labeled` layer
#' of its output GeoPackage.
#'
#' The argument name is `train_labelled`, while the pool object and
#' GeoPackage layer use `train_labeled`.
#'
#' With `split_unit = "fire"`, polygons sharing a `fire_uid` are kept
#' together within each repetition. Ensure that these identifiers represent
#' the intended grouping units.
#'
#' @section Partition selection:
#' The function:
#' 1. Builds candidate square-block grids over the training polygons.
#' 2. Evaluates combinations of block size and fold count against acceptance
#'    criteria, including minimum burned-unit and positive-block counts per
#'    fold.
#' 3. Selects a partition.
#' 4. Assigns repeated spatial folds while preserving the selected splitting
#'    units.
#' 5. Transfers fold assignments back to the polygon table.
#' 6. Writes outputs when requested.
#'
#' The selected settings and acceptance result are returned in `selected`.
#'
#' Larger blocks generally impose broader spatial grouping. Block size alone
#' does not guarantee a minimum distance between every training and test
#' polygon, particularly near block boundaries.
#'
#' @section Fire-level and polygon-level splitting:
#' | Setting | Behaviour |
#' |---|---|
#' | `split_unit = "fire"` | Keep all polygons sharing a `fire_uid` in the same fold within a repetition. |
#' | `split_unit = "poly"` | Use individual polygons as splitting units, without an additional shared-fire grouping constraint. |
#'
#' Use fire-level splitting when multiple polygons belong to the same event
#' and should not be divided between training and evaluation sets.
#'
#' @section Repeated cross-validation:
#' The returned training layer contains `fold_rep*` columns, with one
#' fold-assignment column per repetition.
#'
#' Grouping constraints apply separately within each repetition. A unit may
#' receive a different fold assignment in another repetition.
#'
#' @section Fallback partitions:
#' If no candidate partition passes the acceptance criteria, the function
#' selects the best available non-passing candidate, raises a warning, and
#' records the fallback in a `_FOLD_FALLBACK.txt` audit file.
#'
#' A returned result therefore does not necessarily mean that an acceptable
#' partition was found. Inspect `selected`, including its `ok` value, before
#' continuing to out-of-fold evaluation.
#'
#' A fallback may indicate that the available burned examples are too sparse
#' or spatially concentrated for the requested block sizes and fold counts.
#'
#' @section Reusing outputs:
#' With `overwrite = FALSE`, the function reuses existing outputs when all
#' three expected files are present.
#'
#' Before reusing them, ensure that they correspond to the current training
#' polygons, labels, grouping identifiers, and fold settings. Use
#' `overwrite = TRUE` to recompute assignments after changing those inputs.
#'
#' @section Output files:
#' Outputs are written under `out_dir`, normally the run's `02_FOLDS` folder.
#'
#' | Filename | Contents |
#' |---|---|
#' | `<year>_blocks_<bs>m.gpkg` | Selected spatial block grid. |
#' | `<year>_folds_<bs>m_k<k>_r<reps>_unit-<unit>.csv` | Polygon-level fold assignments. |
#' | `<year>_train_with_folds_<bs>m.gpkg` | Training polygons with `fold_rep*` columns. |
#'
#' Here, `<bs>` is the selected block size, `<k>` is the fold count, `<reps>`
#' is the number of repetitions, and `<unit>` is the splitting unit.
#'
#' @return A named list containing spatial objects, fold assignments,
#'   selected settings, and output paths.
#'
#' | Field | Contents |
#' |---|---|
#' | `blocks` | Selected block grid as an `sf` object. |
#' | `folds_table` | `data.frame` containing polygon-level fold assignments. |
#' | `train_with_folds` | Training polygons as an `sf` object with `fold_rep*` columns. |
#' | `blocks_gpkg` | Path to the block-grid GeoPackage. |
#' | `folds_csv` | Path to the fold-assignment CSV. |
#' | `train_with_folds_gpkg` | Path to the training GeoPackage with fold assignments. |
#' | `selected` | Selected partition settings and diagnostics, including `block_size_m`, `k_folds`, and `ok`. |
#'
#' When `write_outputs = FALSE`, path fields contain the expected target
#' locations. Their presence in the returned list does not mean that files
#' were written.
#'
#' @seealso [build_supervised_burned_config()],
#'   [build_supervised_training_pools()], [extract_supervised_features()],
#'   [run_oof_diagnostics()], [validate_supervised_execution()],
#'   [run_oneyear_supervised_pipeline()].
#'
#' @examples
#' \dontrun{
#' # Configure the probabilistic refinement workflow
#' config <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022"
#' )
#'
#' # Build the labelled training pool
#' pools <- build_supervised_training_pools(
#'   config = config,
#'   write_outputs = TRUE
#' )
#'
#' # Assign spatial folds from the written pool
#' folds <- make_spatial_folds(
#'   train_labelled = pools$pools_gpkg,
#'   config = config,
#'   split_unit = "fire",
#'   block_sizes_m = c(5000, 3000, 2000),
#'   k_candidates = c(5L, 4L, 3L)
#' )
#'
#' # Inspect the selected partition and acceptance result
#' folds$selected
#'
#' # Stop here if the selected partition failed the acceptance criteria
#' if (!isTRUE(folds$selected$ok)) {
#'   stop("Review the fallback partition before continuing.")
#' }
#'
#' # Locate the training layer with fold assignments
#' folds$train_with_folds_gpkg
#'
#' # Inspect the fold columns
#' fold_columns <- grep(
#'   "^fold_rep",
#'   names(folds$train_with_folds),
#'   value = TRUE
#' )
#'
#' head(
#'   sf::st_drop_geometry(folds$train_with_folds)[
#'     , c("class", fold_columns), drop = FALSE
#'   ]
#' )
#'
#' # Alternatively, supply the in-memory labelled pool
#' folds_in_memory <- make_spatial_folds(
#'   train_labelled = pools$train_labeled,
#'   config = config,
#'   split_unit = "fire",
#'   write_outputs = FALSE,
#'   overwrite = TRUE
#' )
#' }
#'
#' @family workflow
#' @export
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
