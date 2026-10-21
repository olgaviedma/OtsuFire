#' Split mixed candidates by DOY and rescore them (update an already processed year)
#'
#' @description
#' Recommended way to **update a year that has already been processed**
#' with the date-aware split of [split_patches_by_doy()]. For one year the
#' function:
#' \enumerate{
#'   \item reads the scored candidates that already exist (the output of the
#'     supervised stage, `deterministic_scored` / `final_map_full`);
#'   \item splits the candidates that mix post-event acquisition dates with
#'     [split_patches_by_doy()];
#'   \item extracts the predictors of the new sub-patches from the same
#'     rasters as the original run;
#'   \item scores them with a given model and recipe (**no retraining**);
#'   \item replaces each split candidate by its sub-patches, applies the
#'     score threshold and writes the final map, an audit table and a
#'     reproducibility recipe.
#' }
#' The deterministic stage is **not** repeated: candidates, their scores and
#' every non-split candidate are those of the original run.
#'
#' Enabling `refine_params = list(split_by_doy = TRUE)` in
#' [build_burned_mapping_config()] is another route (the split happens inside
#' [run_deterministic_pipeline()], before scoring). That route is optional and
#' has **not** been validated on real data; this function is the validated
#' flow (see the vignette "Updating years with the DOY split").
#'
#' @param year Integer scalar. Target year.
#' @param scored_candidates Character path to a scored candidate GeoPackage
#'   (for instance `.../SCORE/<year>_..._deterministic_scored.gpkg`, or the
#'   `final_map_full` layer of `.../MAP/..._final_map.gpkg`) or the same
#'   layer as an `sf` object. It needs the columns `poly_id` and
#'   `p_burned_model`.
#' @param change_index Character path to the summer change-index raster
#'   (band 1 RBR, band 2 DOY of the post-event image).
#' @param delayed_change_index Character path to the autumn change-index
#'   raster (band 1 is used).
#' @param topo Character path to the topography raster (band 1 elevation,
#'   band 2 slope).
#' @param corine_raster Character path to the land-cover raster of the
#'   CORINE epoch of `year`.
#' @param model An xgboost model or the path to its `.rds`.
#' @param recipe The recipe of that model or the path to its `.rds`.
#' @param threshold Numeric scalar in `[0, 1]`. Score threshold of the final
#'   map (for the published results, the data-driven threshold of
#'   `THRESHOLD_FINAL`).
#' @param out_dir Character scalar. Run directory. Must be new or empty of
#'   outputs unless `overwrite = TRUE`.
#' @param doy_gap_days,min_group_share,min_subpatch_px,min_area_ha Passed to
#'   [split_patches_by_doy()]. Defaults: `30`, `0.05`, `9L`, `50`.
#' @param window_margin_m Numeric scalar. Margin (metres) added around the
#'   extent of the sub-patches when cropping the rasters for feature
#'   extraction. Default `3000`.
#' @param extra_config Optional list (or object coercible to a list) stored
#'   verbatim in the recipe, for instance the full supervised configuration
#'   of the original run.
#' @param artifact_hard_weight_ratio Optional numeric. Weight ratio of the
#'   `artifact_hard` pool used to train the model, stored in the recipe (the
#'   model recipe itself does not record it).
#' @param overwrite Logical scalar. Default `FALSE`.
#' @param verbose Logical scalar. Default `TRUE`.
#'
#' @details
#' Predictors are the 14 of the final models (summer/autumn RBR statistics,
#' persistence, CORINE fractions, elevation and slope dispersion), extracted
#' with the same internals as [extract_supervised_features()] with
#' `use_doy = TRUE`, `use_aw = TRUE`, no hotspots and no shape block. A model
#' trained on a different feature set needs the matching recipe; columns the
#' recipe expects but the extraction does not produce are filled as missing
#' by the recipe schema reconciliation.
#'
#' Topography dispersion depends slightly on the cropped window, so metrics
#' can move in the fifth decimal when the set of split candidates changes.
#'
#' @section Output:
#' Under `out_dir`:
#' \itemize{
#'   \item `DOYSPLIT/dsplit_<G>d/...`: output of [split_patches_by_doy()].
#'   \item `subpatches_scored.gpkg`: scored sub-patches with their parent id.
#'   \item `final_map_split.gpkg` (layer `final_map_full`): all candidates,
#'     with split parents replaced by their sub-patches.
#'   \item `final_map_split_thr<NN>.gpkg`: candidates with
#'     `p_burned_model >= threshold`.
#'   \item `run_split_rescore_audit.csv`: counts and parameters.
#'   \item `RECIPE.yml` and `sessionInfo.txt`: see below.
#' }
#' The recipe records the package commit, the parameters (including
#' `split_by_doy`, `doy_gap_days` and the `artifact_hard` weight when given),
#' the input paths, the model training parameters and the threshold.
#'
#' @return Invisibly, a list with `final_map_path`, `thresholded_path`,
#'   `subpatches_path`, `audit` (one-row data frame), `recipe_path`,
#'   `split` (the result of [split_patches_by_doy()]) and `n_kept`.
#'
#' @seealso [split_patches_by_doy()], [score_supervised_burned_map()],
#'   [write_thresholded_burned()], [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' res <- run_split_rescore(
#'   year = 2006L,
#'   scored_candidates = "Results/2006/.../SCORE/2006_deterministic_scored.gpkg",
#'   change_index = "MinMin_2006_mosaic_res90m.tif",
#'   delayed_change_index = "mean_mean_2006_mosaic.tif",
#'   topo = "elevation_slope.tif",
#'   corine_raster = "CLC_2006_peninsula.tif",
#'   model = "final_model.rds", recipe = "recipe.rds",
#'   threshold = 0.51, out_dir = "Results_v2.1/2006"
#' )
#' res$thresholded_path
#' }
#'
#' @family workflow
#' @export
run_split_rescore <- function(year,
                              scored_candidates,
                              change_index,
                              delayed_change_index,
                              topo,
                              corine_raster,
                              model,
                              recipe,
                              threshold,
                              out_dir,
                              doy_gap_days = 30,
                              min_group_share = 0.05,
                              min_subpatch_px = 9L,
                              min_area_ha = 50,
                              window_margin_m = 3000,
                              extra_config = NULL,
                              artifact_hard_weight_ratio = NULL,
                              overwrite = FALSE,
                              verbose = TRUE) {
  # ---- validation (parameters first, files after) ---------------------------
  if (!is.numeric(year) || length(year) != 1L || !is.finite(year) || year != round(year)) {
    stop("'year' must be a single whole number.", call. = FALSE)
  }
  year <- as.integer(year)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("'threshold' must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(window_margin_m) || length(window_margin_m) != 1L ||
      !is.finite(window_margin_m) || window_margin_m < 0) {
    stop("'window_margin_m' must be a single non-negative number.", call. = FALSE)
  }
  if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    stop("'out_dir' must be a single directory path.", call. = FALSE)
  }
  for (nm in c("overwrite", "verbose")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", nm), call. = FALSE)
    }
  }
  if (!is.null(artifact_hard_weight_ratio) &&
      (!is.numeric(artifact_hard_weight_ratio) || length(artifact_hard_weight_ratio) != 1L ||
       !is.finite(artifact_hard_weight_ratio) || artifact_hard_weight_ratio < 0)) {
    stop("'artifact_hard_weight_ratio' must be NULL or a single non-negative number.", call. = FALSE)
  }
  # split parameters are validated by the same checker that split_patches_by_doy() uses
  .of_check_doy_split_params(doy_gap_days, min_group_share, min_subpatch_px, min_area_ha)
  for (nm in c("change_index", "delayed_change_index", "topo", "corine_raster")) {
    p <- get(nm)
    if (!is.character(p) || length(p) != 1L || !file.exists(p)) {
      stop(sprintf("'%s' must be the path to an existing file.", nm), call. = FALSE)
    }
  }
  if (is.character(model) && length(model) == 1L) {
    if (!file.exists(model)) stop("model file not found: ", model, call. = FALSE)
    model_path <- model; model <- readRDS(model)
  } else model_path <- NA_character_
  if (is.character(recipe) && length(recipe) == 1L) {
    if (!file.exists(recipe)) stop("recipe file not found: ", recipe, call. = FALSE)
    recipe_in <- recipe; recipe <- readRDS(recipe)
  } else recipe_in <- NA_character_
  if (!is.list(recipe) || is.null(recipe$cols$x_cols)) {
    stop("'recipe' does not look like a supervised model recipe (no cols$x_cols).", call. = FALSE)
  }
  if (is.character(scored_candidates) && length(scored_candidates) == 1L) {
    if (!file.exists(scored_candidates)) stop("scored_candidates not found: ", scored_candidates, call. = FALSE)
    cand_path <- scored_candidates
    lyr <- sf::st_layers(cand_path)$name
    scored <- sf::st_read(cand_path, layer = if ("final_map_full" %in% lyr) "final_map_full" else lyr[1], quiet = TRUE)
  } else if (inherits(scored_candidates, "sf")) {
    scored <- scored_candidates; cand_path <- NA_character_
  } else {
    stop("'scored_candidates' must be a path or an sf object.", call. = FALSE)
  }
  miss <- setdiff(c("poly_id", "p_burned_model"), names(scored))
  if (length(miss)) stop("scored_candidates lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  final_path <- file.path(out_dir, "final_map_split.gpkg")
  if (file.exists(final_path) && !isTRUE(overwrite)) {
    stop("Output exists and overwrite = FALSE: ", final_path, call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  say <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  t0 <- Sys.time()

  # ---- 1) split --------------------------------------------------------------
  say("[run_split_rescore] %d: splitting %d candidates (G = %g d)", year, nrow(scored), doy_gap_days)
  sp <- split_patches_by_doy(refine_path = scored, change_index = change_index, year = year,
                             doy_gap_days = doy_gap_days, min_group_share = min_group_share,
                             min_subpatch_px = min_subpatch_px, min_area_ha = min_area_ha,
                             id_col = "poly_id", out_dir = NULL, verbose = FALSE)
  # persist the split next to the other outputs (split_patches_by_doy() needs a path-style out_dir)
  sdir <- file.path(out_dir, "DOYSPLIT", sprintf("dsplit_%gd", doy_gap_days))
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  split_gpkg <- file.path(sdir, sprintf("BA_%d_REFINE_MERGED_dsplit%gd.gpkg", year, doy_gap_days))
  if (file.exists(split_gpkg)) unlink(split_gpkg, force = TRUE)
  sf::st_write(sp$patches, split_gpkg, quiet = TRUE)
  pa <- sp$patches_audit
  if (!is.null(pa) && nrow(pa)) utils::write.csv(pa, file.path(sdir, "doysplit_patches.csv"), row.names = FALSE)
  if (!is.null(sp$parents_audit)) utils::write.csv(sp$parents_audit, file.path(sdir, "doysplit_parents.csv"), row.names = FALSE)

  n_sub <- if (is.null(pa)) 0L else nrow(pa)
  scored_sub <- NULL
  if (n_sub > 0L) {
    # ---- 2) features of the sub-patches -------------------------------------
    g <- sp$patches[match(pa$CLUMP_ID, sp$patches$CLUMP_ID), "CLUMP_ID"]
    g$sub_id <- as.character(pa$sub_id); g$parent_id <- as.character(pa$parent_id)
    g$area_ha <- pa$area_ha; g$n_px <- pa$n_px; g$doy_med <- pa$doy_med
    ft <- .of_rescore_extract(g, year, change_index, delayed_change_index, topo, corine_raster,
                              window_margin_m)
    # ---- 3) score ------------------------------------------------------------
    g$p_burned <- .of_rescore_predict(sf::st_drop_geometry(ft), recipe, model)
    for (k in setdiff(names(ft), c("fire_uid", "class", attr(ft, "sf_column")))) {
      if (!k %in% names(g)) g[[k]] <- ft[[k]]
    }
    scored_sub <- g
    sub_path <- file.path(out_dir, "subpatches_scored.gpkg")
    if (file.exists(sub_path)) unlink(sub_path, force = TRUE)
    sf::st_write(sf::st_transform(scored_sub, sf::st_crs(scored)), sub_path, quiet = TRUE)
  } else sub_path <- NA_character_

  # ---- 4) assemble the split map ------------------------------------------------
  if (is.null(scored_sub)) {
    m_split <- scored
  } else {
    base <- scored[!(as.character(scored$poly_id) %in% unique(scored_sub$parent_id)), ]
    gcol <- attr(base, "sf_column")
    add <- sf::st_transform(scored_sub, sf::st_crs(base))[, 0]
    add$poly_id <- scored_sub$sub_id; add$p_burned_model <- scored_sub$p_burned
    for (cn in setdiff(names(base), c(names(add), gcol))) add[[cn]] <- NA
    sf::st_geometry(add) <- gcol
    m_split <- rbind(base, add[, names(base)])
  }
  if (file.exists(final_path)) unlink(final_path, force = TRUE)
  sf::st_write(m_split, final_path, layer = "final_map_full", quiet = TRUE)
  w <- write_thresholded_burned(scored = m_split, out_dir = out_dir, threshold = threshold,
                                score_col = "p_burned_model",
                                filename = sprintf("final_map_split_thr%02d.gpkg", round(threshold * 100)),
                                overwrite = TRUE, verbose = FALSE)

  n_thr_orig <- sum(is.finite(scored$p_burned_model) & scored$p_burned_model >= threshold)
  audit <- data.frame(year = year, threshold = threshold, doy_gap_days = doy_gap_days,
                      n_candidates_original = nrow(scored), n_examined = sp$n_examined,
                      n_split_parents = sp$n_mixed, n_subpatches = sp$n_subpatches,
                      n_split_errors = sp$n_errors, n_candidates_split_map = nrow(m_split),
                      n_kept_original = n_thr_orig, n_kept_split = w$n_kept,
                      area_kept_original_ha = round(sum(as.numeric(sf::st_area(scored[is.finite(scored$p_burned_model) & scored$p_burned_model >= threshold, ]))) / 1e4, 1),
                      area_kept_split_ha = round(sum(as.numeric(sf::st_area(sf::st_read(w$path, quiet = TRUE)))) / 1e4, 1),
                      minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  utils::write.csv(audit, file.path(out_dir, "run_split_rescore_audit.csv"), row.names = FALSE)

  # ---- 5) recipe -----------------------------------------------------------------
  rp <- .of_write_run_recipe(
    out_dir = out_dir, function_name = "run_split_rescore",
    params = list(year = year, threshold = threshold, split_by_doy = TRUE,
                  doy_gap_days = doy_gap_days, min_group_share = min_group_share,
                  min_subpatch_px = min_subpatch_px, min_area_ha = min_area_ha,
                  window_margin_m = window_margin_m,
                  artifact_hard_weight_ratio = if (is.null(artifact_hard_weight_ratio)) "not recorded" else artifact_hard_weight_ratio),
    inputs = list(scored_candidates = cand_path, change_index = change_index,
                  delayed_change_index = delayed_change_index, topo = topo,
                  corine_raster = corine_raster, model = model_path, recipe = recipe_in),
    model_info = list(training = recipe$training, params = recipe$params,
                      x_cols = recipe$cols$x_cols,
                      schema_hash = recipe$schema_fingerprint$hash),
    extra_config = extra_config)
  say("[run_split_rescore] %d done: %d parents split into %d sub-patches; kept %d -> %d (thr %.2f)",
      year, sp$n_mixed, sp$n_subpatches, n_thr_orig, w$n_kept, threshold)

  invisible(list(final_map_path = final_path, thresholded_path = w$path, subpatches_path = sub_path,
                 audit = audit, recipe_path = rp, split = sp, n_kept = w$n_kept))
}

# ---- internals ---------------------------------------------------------------

# Predictors of sub-patches, from the same rasters as the original run.
.of_rescore_extract <- function(g, year, change_index, delayed_change_index, topo,
                                corine_raster, window_margin_m) {
  r <- terra::rast(change_index)
  bb <- sf::st_bbox(sf::st_transform(g, terra::crs(r)))
  e <- terra::ext(bb["xmin"] - window_margin_m, bb["xmax"] + window_margin_m,
                  bb["ymin"] - window_margin_m, bb["ymax"] + window_margin_m)
  e <- terra::intersect(e, terra::ext(r))
  tmpl <- terra::crop(r[[1]], e)
  tp <- terra::rast(topo)
  tc <- terra::crop(tp[[1:2]], .of_snap_extent(terra::project(e, terra::crs(r), terra::crs(tp)), tp))
  dem <- align_to_template(tc[[1]], tmpl, "bilinear", "dem", verbose = FALSE)
  slp <- align_to_template(tc[[2]], tmpl, "bilinear", "slope", verbose = FALSE)
  clc <- align_to_template(terra::rast(corine_raster), tmpl, "near", "corine", verbose = FALSE)
  aw  <- align_to_template(terra::rast(delayed_change_index)[[1]], tmpl, "bilinear", "rbr_aw", verbose = FALSE)
  dp  <- align_to_template(terra::crop(r[[2]], e), tmpl, "near", "doy", verbose = FALSE)
  u <- g
  u$fire_uid <- as.character(g$sub_id); u$class <- NA_character_
  u <- u[, c("fire_uid", "class")]
  tr <- u[1:2, ]; tr$class <- c("burned", "unburned"); tr$block_id <- c("a", "b")
  tr$fold_rep1 <- 1:2; tr$fold_rep2 <- 1:2
  ft <- extract_features(train_folds = tr, unlabeled = u, id_col = "fire_uid", class_col = "class",
                         build_features = TRUE, rbr_summer = tmpl, doy_post = dp, rbr_aw = aw,
                         dem = dem, slope = slp, corine_r = clc, hotspots = NULL, use_hotspots = FALSE,
                         cor_groups = cor_groups, use_doy = TRUE, use_aw = TRUE, use_nbr = FALSE,
                         use_shape = FALSE, year_target = year, return_features = TRUE,
                         return_features_geometry = TRUE, verbose = FALSE)$unl_feat
  ft[match(as.character(g$sub_id), ft$fire_uid), ]
}

# Score a feature table with a recipe and a model (no refit).
.of_rescore_predict <- function(df, recipe, model) {
  rec <- .of_reconcile_scoring_schema(as.data.frame(df), recipe)
  Xr <- .of_nested_coerce_features(rec$df, recipe$impute$impute_factor_missing %||% "MISSING")
  Xi <- .of_nested_apply_medians(Xr, recipe$impute$numeric_medians %||% list())
  M <- .of_nested_build_matrix(Xi, ref_cols = recipe$cols$x_cols)
  stats::predict(model, xgboost::xgb.DMatrix(M, missing = NA))
}

# Reproducibility recipe: package commit, parameters, inputs, model info, sessionInfo().
.of_write_run_recipe <- function(out_dir, function_name, params, inputs, model_info = NULL,
                                 extra_config = NULL) {
  pkg_dir <- tryCatch(system.file(package = "OtsuFire"), error = function(e) "")
  commit <- branch <- NA_character_
  if (nzchar(pkg_dir) && dir.exists(file.path(pkg_dir, ".git"))) {
    gi <- function(a) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(pkg_dir), a), stdout = TRUE, stderr = FALSE))[1],
                               error = function(e) NA_character_)
    commit <- gi(c("rev-parse", "HEAD")); branch <- gi(c("rev-parse", "--abbrev-ref", "HEAD"))
    dirty <- tryCatch(length(suppressWarnings(system2("git", c("-C", shQuote(pkg_dir), "status", "--porcelain", "--untracked-files=no"), stdout = TRUE, stderr = FALSE))) > 0L,
                      error = function(e) NA)
  } else dirty <- NA
  rec <- list(
    function_name = function_name,
    written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    package = list(name = "OtsuFire", version = as.character(utils::packageVersion("OtsuFire")),
                   git_commit = commit, git_branch = branch, git_tracked_changes = dirty),
    params = params, inputs = inputs, model = model_info,
    extra_config = if (is.null(extra_config)) NULL else as.list(extra_config))
  path <- file.path(out_dir, "RECIPE.yml")
  writeLines(c("# OtsuFire run recipe (reproducibility record)", .of_to_yaml(rec)), path)
  utils::capture.output(print(utils::sessionInfo()), file = file.path(out_dir, "sessionInfo.txt"))
  path
}

# Minimal YAML-style serialiser (no extra dependency): lists nest, vectors become inline sequences.
.of_to_yaml <- function(x, indent = 0L) {
  pad <- strrep("  ", indent)
  if (is.null(x)) return(paste0(pad, "~"))
  fmt <- function(v) if (is.character(v)) paste0("\"", gsub("\"", "'", v), "\"") else as.character(v)
  if (!is.list(x)) {
    v <- vapply(x, function(z) if (is.na(z)) "null" else fmt(z), character(1))
    return(paste0(pad, if (length(v) == 1L) v else paste0("[", paste(v, collapse = ", "), "]")))
  }
  if (!length(x)) return(paste0(pad, "{}"))
  nm <- names(x); if (is.null(nm)) nm <- paste0("item", seq_along(x)); nm[!nzchar(nm)] <- "item"
  out <- character(0)
  for (k in seq_along(x)) {
    el <- x[[k]]
    if (is.list(el) && length(el)) out <- c(out, paste0(pad, nm[k], ":"), .of_to_yaml(el, indent + 1L))
    else out <- c(out, paste0(pad, nm[k], ": ", trimws(.of_to_yaml(el, 0L))))
  }
  out
}

# terra::project() aggregates a fine source before resampling to a coarser template, and the aggregation
# blocks start at the first cell of the (cropped) source. Cropping on block boundaries of the full raster
# makes the derived values independent of the window, so a run does not change when the set of split
# candidates (hence the window) changes. 12 is a common multiple of every plausible block size.
.of_snap_extent <- function(e, full, block = 12L) {
  c0 <- terra::colFromX(full, e[1]); c1 <- terra::colFromX(full, e[2])
  r0 <- terra::rowFromY(full, e[4]); r1 <- terra::rowFromY(full, e[3])
  c0 <- if (is.na(c0)) 1L else c0; r0 <- if (is.na(r0)) 1L else r0
  c1 <- if (is.na(c1)) terra::ncol(full) else c1; r1 <- if (is.na(r1)) terra::nrow(full) else r1
  c0 <- ((c0 - 1L) %/% block) * block + 1L; r0 <- ((r0 - 1L) %/% block) * block + 1L
  c1 <- min(terra::ncol(full), ceiling(c1 / block) * block + block)
  r1 <- min(terra::nrow(full), ceiling(r1 / block) * block + block)
  hx <- terra::res(full)[1] / 2; hy <- terra::res(full)[2] / 2
  terra::ext(terra::xFromCol(full, c0) - hx, terra::xFromCol(full, c1) + hx,
             terra::yFromRow(full, r1) - hy, terra::yFromRow(full, r0) + hy)
}
