# Historical keep-like reference helpers for the deterministic stage.
#
# These two functions provide an EXPLICIT (V1) way to build a reusable
# spectral-support reference (`keep_pool`) from previously written
# deterministic decision layers:
#
#   1. collect_keep_reference_samples()  -> select reliable historical
#      keep-like polygons using a fixed strict_hotspot rule.
#   2. build_keep_pool_from_samples()    -> turn those samples + per-year
#      change-index rasters into a reusable keep_pool object.
#
# They remain explicit helpers: the package does not auto-discover or
# auto-build historical keep pools on its own. However, the keep_pool
# produced by build_keep_pool_from_samples() is now accepted directly by
# score_burned_patches() and run_deterministic_pipeline(keep_pool = ...).

# ---------------------------------------------------------------------------
# internal helper: row-bind a list of sf objects with column alignment.
# Different decision files may carry slightly different optional columns;
# we keep the union (no aggressive reduction) and fill missing ones with NA.
# ---------------------------------------------------------------------------
.of_rbind_sf_aligned <- function(sf_list) {
  sf_list <- Filter(Negate(is.null), sf_list)
  if (length(sf_list) == 0L) stop("Internal error: empty sf list to bind.")
  crs0 <- sf::st_crs(sf_list[[1L]])

  geoms <- vector("list", length(sf_list))
  dfs   <- vector("list", length(sf_list))
  for (i in seq_along(sf_list)) {
    x <- sf_list[[i]]
    if (!identical(sf::st_crs(x), crs0)) {
      x <- sf::st_transform(x, crs0)
    }
    geoms[[i]] <- sf::st_geometry(x)
    dfs[[i]]   <- sf::st_drop_geometry(x)
  }

  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })

  combined_df   <- do.call(rbind, dfs)
  combined_geom <- do.call(c, geoms)
  out <- sf::st_sf(combined_df, geometry = combined_geom)
  sf::st_crs(out) <- crs0
  out
}

#' Collect reliable historical keep-like reference samples
#'
#' @description
#' Collect high-confidence historical \emph{keep}-like polygons from one
#' or more previously written deterministic decision layers, using a
#' single fixed selection rule (\code{"strict_hotspot"}).
#'
#' In the public deterministic workflow, this helper is meant for
#' building explicit \emph{external} reference pools from trusted years
#' other than the year currently being scored. The selected samples are
#' intended to feed
#' \code{\link[=build_keep_pool_from_samples]{build_keep_pool_from_samples()}},
#' which turns them into a reusable spectral-support \code{keep_pool}.
#'
#' This is the first of a two-step explicit workflow:
#' \enumerate{
#'   \item \strong{Select} reliable historical samples (this function).
#'   \item \strong{Build} a reusable \code{keep_pool} from those samples
#'     (\code{build_keep_pool_from_samples()}).
#' }
#'
#' \strong{Explicit helper workflow.} Neither this function nor
#' \code{build_keep_pool_from_samples()} is called automatically by the
#' deterministic pipeline. They are standalone helpers. If you choose to
#' use them, the resulting \code{keep_pool} can then be passed
#' explicitly to
#' \code{score_burned_patches(..., keep_pool = ...)} or
#' \code{run_deterministic_pipeline(..., keep_pool = ...)}.
#'
#' @details
#' \strong{Current helper design.} The caller supplies the exact decision
#' layer paths in \code{decision_paths}. There is no autodiscovery, and
#' the selection rule described below is fixed in the public helper.
#'
#' \strong{Fixed \code{strict_hotspot} rule.} A polygon is selected only
#' when \emph{all} of the following hold:
#' \itemize{
#'   \item \code{class_final == "keep"};
#'   \item \code{filter_1 == "keep"};
#'   \item \code{filter_2 == "keep"};
#'   \item \code{filter_3 == "keep"};
#'   \item \code{area_ha >= min_area_ha};
#'   \item \code{conf_area == "ok"};
#'   \item \code{conf_pix == "ok"};
#'   \item if a \code{preyear_action} column exists,
#'     \code{preyear_action != "drop"}.
#' }
#'
#' \code{run_label} is derived per polygon as the first available of
#' \code{run_name}, then \code{scenario}, otherwise \code{NA_character_}.
#'
#' @param decision_paths Character vector of existing paths to
#'   deterministic decision files. This argument is required; the helper
#'   does not search for files automatically.
#' @param target_year Optional scalar year. When supplied, polygons with
#'   \code{year == target_year} are excluded. This is the recommended
#'   setting when the helper is being used to build an explicit
#'   historical reference for scoring a specific target year, because it
#'   avoids leaking the scored year into its own explicit reference.
#'   Helper-generated explicit pools that still include the target year
#'   are rejected later by
#'   \code{\link[=score_burned_patches]{score_burned_patches()}} and
#'   \code{\link[=run_deterministic_pipeline]{run_deterministic_pipeline()}}
#'   when their metadata reveal that overlap.
#' @param years_include Optional numeric/integer vector. When supplied,
#'   only these years are kept.
#' @param years_exclude Optional numeric/integer vector. When supplied,
#'   these years are excluded.
#' @param run_names_include Optional character vector. When supplied, only
#'   polygons whose \code{run_label} is in this set are kept.
#' @param layer Character scalar. Layer to read from each path. Default
#'   \code{"internal_decisions"}.
#' @param min_area_ha Numeric scalar > 0. Minimum polygon area (ha) for
#'   the strict rule. Default \code{100}.
#' @param quiet Logical scalar. Passed to the layer reader. Default
#'   \code{FALSE}.
#'
#' @return An object of class
#'   \code{c("otsufire_keep_reference_samples", "list")} with fields:
#' \describe{
#'   \item{\code{samples}}{\code{sf} of selected polygons (original useful
#'     columns preserved, plus \code{run_label} and \code{decision_path}).}
#'   \item{\code{summary}}{list with \code{n_files_input},
#'     \code{n_polygons_input}, \code{n_polygons_selected},
#'     \code{years_used}, \code{run_labels_used},
#'     \code{target_year_excluded}, \code{min_area_ha}, \code{rule}.}
#'   \item{\code{exclusion_counts}}{auditable \code{data.frame} with
#'     columns \code{rule} and \code{n} (stable schema; see Details).}
#'   \item{\code{decision_paths}}{the input paths.}
#' }
#'
#' The \code{exclusion_counts} schema is stable and always contains, in
#' order: \code{input_total}, \code{excluded_target_year},
#' \code{excluded_years_include}, \code{excluded_years_exclude},
#' \code{excluded_run_names_include}, \code{failed_class_final_keep},
#' \code{failed_filter_1_keep}, \code{failed_filter_2_keep},
#' \code{failed_filter_3_keep}, \code{failed_min_area_ha},
#' \code{failed_conf_area}, \code{failed_conf_pix},
#' \code{failed_preyear_drop}, \code{selected}. The \code{excluded_*}
#' rows count polygons removed sequentially by the optional pre-filters;
#' each \code{failed_*} row counts polygons (among those entering the
#' strict rule) that independently fail that single predicate (counts may
#' overlap and are for auditing, not a partition); \code{selected} counts
#' polygons passing every predicate.
#'
#' @seealso
#' \code{\link[=build_keep_pool_from_samples]{build_keep_pool_from_samples()}},
#' \code{\link[=score_burned_patches]{score_burned_patches()}}
#' @export
collect_keep_reference_samples <- function(decision_paths,
                                           target_year = NULL,
                                           years_include = NULL,
                                           years_exclude = NULL,
                                           run_names_include = NULL,
                                           layer = "internal_decisions",
                                           min_area_ha = 100,
                                           quiet = FALSE) {

  # ---- initial validation ------------------------------------------------
  if (missing(decision_paths) || is.null(decision_paths) ||
      !is.character(decision_paths) || length(decision_paths) == 0L) {
    stop("'decision_paths' is required and must be a non-empty character vector.",
         call. = FALSE)
  }
  miss_paths <- decision_paths[!file.exists(decision_paths)]
  if (length(miss_paths) > 0L) {
    stop("'decision_paths': path(s) do not exist:\n  ",
         paste(miss_paths, collapse = "\n  "), call. = FALSE)
  }
  if (!is.character(layer) || length(layer) != 1L || is.na(layer) ||
      !nzchar(layer)) {
    stop("'layer' must be a non-empty scalar string.", call. = FALSE)
  }
  if (!is.numeric(min_area_ha) || length(min_area_ha) != 1L ||
      !is.finite(min_area_ha) || min_area_ha <= 0) {
    stop("'min_area_ha' must be a single positive numeric value.",
         call. = FALSE)
  }
  if (!is.null(target_year) &&
      (length(target_year) != 1L || !is.numeric(target_year) ||
       !is.finite(target_year))) {
    stop("'target_year' must be a single numeric/integer value when supplied.",
         call. = FALSE)
  }
  if (!is.null(years_include) &&
      (!is.numeric(years_include) || any(!is.finite(years_include)))) {
    stop("'years_include' must be a numeric/integer vector when supplied.",
         call. = FALSE)
  }
  if (!is.null(years_exclude) &&
      (!is.numeric(years_exclude) || any(!is.finite(years_exclude)))) {
    stop("'years_exclude' must be a numeric/integer vector when supplied.",
         call. = FALSE)
  }
  if (!is.null(run_names_include) && !is.character(run_names_include)) {
    stop("'run_names_include' must be a character vector when supplied.",
         call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L) {
    stop("'quiet' must be TRUE or FALSE.", call. = FALSE)
  }

  required_cols <- c("year", "class_final", "filter_1", "filter_2",
                     "filter_3", "area_ha", "conf_area", "conf_pix")

  # ---- read + normalize each file ---------------------------------------
  per_file <- vector("list", length(decision_paths))
  for (i in seq_along(decision_paths)) {
    p <- decision_paths[[i]]

    lyr_info <- tryCatch(sf::st_layers(p), error = function(e) NULL)
    if (is.null(lyr_info) || !layer %in% lyr_info$name) {
      stop("Cannot find layer '", layer, "' in: ", p, call. = FALSE)
    }
    x <- tryCatch(
      sf::st_read(p, layer = layer, quiet = isTRUE(quiet)),
      error = function(e)
        stop("Failed to read layer '", layer, "' from: ", p,
             "\n  ", conditionMessage(e), call. = FALSE)
    )
    if (!inherits(x, "sf")) {
      stop("Layer '", layer, "' in ", p, " is not an sf object.",
           call. = FALSE)
    }

    miss_cols <- setdiff(required_cols, names(x))
    if (length(miss_cols) > 0L) {
      stop("Missing required column(s) in ", p, " [", layer, "]: ",
           paste(miss_cols, collapse = ", "), call. = FALSE)
    }

    run_name_v <- if ("run_name" %in% names(x))
      as.character(x[["run_name"]]) else rep(NA_character_, nrow(x))
    scenario_v <- if ("scenario" %in% names(x))
      as.character(x[["scenario"]]) else rep(NA_character_, nrow(x))
    run_label <- run_name_v
    use_scn <- (is.na(run_label) | !nzchar(run_label)) &
      !is.na(scenario_v) & nzchar(scenario_v)
    run_label[use_scn] <- scenario_v[use_scn]
    run_label[is.na(run_label) | !nzchar(run_label)] <- NA_character_

    x[["run_label"]] <- run_label
    x[["decision_path"]] <- p
    per_file[[i]] <- x
  }

  all_sf <- .of_rbind_sf_aligned(per_file)
  input_total <- nrow(all_sf)

  # ---- optional pre-filters (sequential, counted) -----------------------
  excluded_target_year <- 0L
  excluded_years_include <- 0L
  excluded_years_exclude <- 0L
  excluded_run_names_include <- 0L

  yr <- suppressWarnings(as.integer(all_sf[["year"]]))

  if (!is.null(target_year)) {
    drop_idx <- which(yr == as.integer(target_year))
    excluded_target_year <- length(drop_idx)
    if (length(drop_idx) > 0L) {
      all_sf <- all_sf[-drop_idx, , drop = FALSE]
      yr <- yr[-drop_idx]
    }
  }
  if (!is.null(years_include)) {
    keep_idx <- which(yr %in% as.integer(years_include))
    excluded_years_include <- nrow(all_sf) - length(keep_idx)
    all_sf <- all_sf[keep_idx, , drop = FALSE]
    yr <- yr[keep_idx]
  }
  if (!is.null(years_exclude)) {
    drop_idx <- which(yr %in% as.integer(years_exclude))
    excluded_years_exclude <- length(drop_idx)
    if (length(drop_idx) > 0L) {
      all_sf <- all_sf[-drop_idx, , drop = FALSE]
      yr <- yr[-drop_idx]
    }
  }
  if (!is.null(run_names_include)) {
    keep_idx <- which(all_sf[["run_label"]] %in% run_names_include)
    excluded_run_names_include <- nrow(all_sf) - length(keep_idx)
    all_sf <- all_sf[keep_idx, , drop = FALSE]
  }

  # ---- fixed strict_hotspot rule (independent failure audit) ------------
  n_in <- nrow(all_sf)
  if (n_in == 0L) {
    stop("No polygons remain after the optional pre-filters; ",
         "cannot apply the strict_hotspot rule.", call. = FALSE)
  }

  area_ha <- suppressWarnings(as.numeric(all_sf[["area_ha"]]))

  ok_class   <- !is.na(all_sf[["class_final"]]) & all_sf[["class_final"]] == "keep"
  ok_f1      <- !is.na(all_sf[["filter_1"]])    & all_sf[["filter_1"]] == "keep"
  ok_f2      <- !is.na(all_sf[["filter_2"]])    & all_sf[["filter_2"]] == "keep"
  ok_f3      <- !is.na(all_sf[["filter_3"]])    & all_sf[["filter_3"]] == "keep"
  ok_area    <- is.finite(area_ha) & area_ha >= min_area_ha
  ok_cfarea  <- !is.na(all_sf[["conf_area"]])   & all_sf[["conf_area"]] == "ok"
  ok_cfpix   <- !is.na(all_sf[["conf_pix"]])    & all_sf[["conf_pix"]] == "ok"
  if ("preyear_action" %in% names(all_sf)) {
    ok_preyear <- is.na(all_sf[["preyear_action"]]) |
      all_sf[["preyear_action"]] != "drop"
  } else {
    ok_preyear <- rep(TRUE, n_in)
  }

  selected_mask <- ok_class & ok_f1 & ok_f2 & ok_f3 &
    ok_area & ok_cfarea & ok_cfpix & ok_preyear

  samples <- all_sf[selected_mask, , drop = FALSE]
  if (nrow(samples) == 0L) {
    stop("No polygons satisfy the strict_hotspot rule; ",
         "no reference samples could be collected.", call. = FALSE)
  }

  years_used <- sort(unique(suppressWarnings(as.integer(samples[["year"]]))))
  run_labels_used <- unique(as.character(samples[["run_label"]]))

  exclusion_counts <- data.frame(
    rule = c("input_total",
             "excluded_target_year",
             "excluded_years_include",
             "excluded_years_exclude",
             "excluded_run_names_include",
             "failed_class_final_keep",
             "failed_filter_1_keep",
             "failed_filter_2_keep",
             "failed_filter_3_keep",
             "failed_min_area_ha",
             "failed_conf_area",
             "failed_conf_pix",
             "failed_preyear_drop",
             "selected"),
    n = c(input_total,
          excluded_target_year,
          excluded_years_include,
          excluded_years_exclude,
          excluded_run_names_include,
          sum(!ok_class),
          sum(!ok_f1),
          sum(!ok_f2),
          sum(!ok_f3),
          sum(!ok_area),
          sum(!ok_cfarea),
          sum(!ok_cfpix),
          sum(!ok_preyear),
          nrow(samples)),
    stringsAsFactors = FALSE
  )
  exclusion_counts$n <- as.integer(exclusion_counts$n)

  out <- list(
    samples = samples,
    summary = list(
      n_files_input        = length(decision_paths),
      n_polygons_input     = as.integer(input_total),
      n_polygons_selected  = as.integer(nrow(samples)),
      years_used           = years_used,
      run_labels_used      = run_labels_used,
      target_year_excluded = if (is.null(target_year)) NULL
                             else as.integer(target_year),
      min_area_ha          = as.numeric(min_area_ha),
      rule                 = "strict_hotspot"
    ),
    exclusion_counts = exclusion_counts,
    decision_paths   = as.character(decision_paths)
  )
  class(out) <- c("otsufire_keep_reference_samples", "list")
  out
}

#' Build a reusable keep_pool from historical reference samples
#'
#' @description
#' Build a reusable spectral-support \code{keep_pool} object from clean
#' historical samples (typically produced by
#' \code{\link[=collect_keep_reference_samples]{collect_keep_reference_samples()}})
#' and a set of per-year change-index rasters. The returned object is
#' structurally compatible with the explicit \code{keep_pool} argument
#' already accepted by
#' \code{\link[=score_burned_patches]{score_burned_patches()}}.
#'
#' This is the second of a two-step explicit workflow:
#' \enumerate{
#'   \item \strong{Select} reliable historical samples
#'     (\code{collect_keep_reference_samples()}).
#'   \item \strong{Build} a reusable \code{keep_pool} from those samples
#'     (this function).
#' }
#'
#' \strong{Explicit helper workflow.} This function is not called
#' automatically by
#' \code{\link[=score_burned_patches]{score_burned_patches()}} or
#' \code{\link[=run_deterministic_pipeline]{run_deterministic_pipeline()}}.
#' It is a standalone helper; the resulting \code{keep_pool} can be
#' passed explicitly via
#' \code{score_burned_patches(..., keep_pool = ...)} or
#' \code{run_deterministic_pipeline(..., keep_pool = ...)}.
#'
#' In the supported public workflow, these explicit pools should come
#' from trusted \emph{external} years. Same-year explicit pools are not a
#' supported public mode; for same-year scoring, use the native local
#' fallback (\code{keep_pool = NULL}).
#'
#' @details
#' \strong{Current helper design.} The caller supplies exact raster paths
#' in \code{change_index_paths}, which must be a \emph{named} list or
#' vector keyed by year (for example
#' \code{list(`2010` = "...", `2011` = "...")}). There is no
#' autodiscovery. This function performs \strong{no methodological
#' filtering}: it assumes \code{samples} are already clean, as returned
#' by \code{collect_keep_reference_samples()}.
#'
#' For each year present in \code{samples}, the corresponding raster is
#' opened, sample geometries are reprojected to the raster CRS if needed,
#' pixel values are extracted per polygon, finite values are accumulated
#' to derive reference quantiles, and a per-polygon median is accumulated
#' into \code{keep_medians}. Non-finite values are ignored.
#'
#' \strong{Fixed internal parameters} (not exposed): \code{qref_prob =
#' 0.05}, \code{promote_p_above_ref = 0.50}, \code{promote_percentile =
#' 0.10}, \code{min_area_ha = 10}. \code{min_pix} is also not exposed but
#' is \emph{derived} (not \code{NULL}) using the exact same logic as the
#' scoring engine, \code{max(ceiling(min_area_ha / pixel_area_ha), 13)},
#' where \code{pixel_area_ha} comes from the change-index raster
#' resolution. This guarantees functional compatibility: when a
#' \code{keep_pool} is supplied, the engine does not re-derive
#' \code{min_pix} and uses \code{keep_pool$min_pix} directly for the
#' \code{conf_pix} check. An internal cap and chunking are applied to the
#' accumulated pixel values to bound memory; the full finite-pixel count
#' is reported in \code{metadata$n_pixels_used}.
#'
#' @param samples \code{sf} of clean keep-like reference polygons (e.g.
#'   \code{collect_keep_reference_samples()$samples}). Must be non-empty.
#' @param change_index_paths Named list or named character vector mapping
#'   year (as the name) to an existing change-index raster path. Every
#'   year present in \code{samples[[year_col]]} must have an entry.
#' @param year_col Character scalar. Column in \code{samples} holding the
#'   year. Default \code{"year"}.
#' @param sample_id_col Character scalar. Column in \code{samples} holding
#'   a per-polygon identifier. Default \code{"source_poly_id"}.
#' @param min_samples Integer/numeric scalar \code{>= 1}. Minimum number
#'   of sample polygons required. Default \code{25L}.
#' @param quiet Logical scalar. Suppress progress messages. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{c("otsufire_keep_pool", "list")} with
#'   fields \code{min_area_ha}, \code{min_pix}, \code{promote_percentile},
#'   \code{promote_p_above_ref}, \code{qref_prob}, \code{qref_keep},
#'   \code{q25_keep}, \code{pixel_area_ha}, \code{keep_medians}, and
#'   \code{metadata} (\code{n_samples}, \code{n_years}, \code{years_used},
#'   \code{n_pixels_used}, \code{year_col}, \code{sample_id_col},
#'   \code{change_index_paths}). The first nine fields match the
#'   \code{keep_pool} contract consumed by
#'   \code{\link[=score_burned_patches]{score_burned_patches()}}.
#'
#' @seealso
#' \code{\link[=collect_keep_reference_samples]{collect_keep_reference_samples()}},
#' \code{\link[=score_burned_patches]{score_burned_patches()}}
#' @export
build_keep_pool_from_samples <- function(samples,
                                         change_index_paths,
                                         year_col = "year",
                                         sample_id_col = "source_poly_id",
                                         min_samples = 25L,
                                         quiet = FALSE) {

  # fixed V1 internal parameters (not exposed to the caller).
  # `min_pix` is NOT exposed either, but it is NOT left NULL: it is
  # derived below with the exact same logic the scoring engine uses
  # (score_rbr_keep_classes()), so the returned keep_pool is functionally
  # compatible (the engine does not re-derive min_pix when keep_pool is
  # supplied; it uses keep_pool$min_pix directly for the conf_pix check).
  qref_prob           <- 0.05
  promote_p_above_ref <- 0.50
  promote_percentile  <- 0.10
  min_area_ha         <- 10
  max_keep_samples    <- 5e6
  chunk_size          <- 200L

  msg <- function(...) if (!isTRUE(quiet)) message(...)

  # ---- initial validation ------------------------------------------------
  if (missing(samples) || is.null(samples) || !inherits(samples, "sf") ||
      nrow(samples) == 0L) {
    stop("'samples' must be a non-empty sf object.", call. = FALSE)
  }
  if (missing(change_index_paths) || is.null(change_index_paths) ||
      !(is.list(change_index_paths) || is.character(change_index_paths))) {
    stop("'change_index_paths' must be a named list or named vector.",
         call. = FALSE)
  }
  cip_names <- names(change_index_paths)
  if (length(change_index_paths) == 0L || is.null(cip_names) ||
      any(is.na(cip_names)) || any(!nzchar(cip_names))) {
    stop("'change_index_paths' must be named by year (all names non-empty).",
         call. = FALSE)
  }
  cip_paths <- vapply(seq_along(change_index_paths),
                      function(k) as.character(change_index_paths[[k]])[1L],
                      character(1))
  miss_paths <- cip_paths[!file.exists(cip_paths)]
  if (length(miss_paths) > 0L) {
    stop("'change_index_paths': raster path(s) do not exist:\n  ",
         paste(miss_paths, collapse = "\n  "), call. = FALSE)
  }
  if (!is.character(year_col) || length(year_col) != 1L ||
      is.na(year_col) || !nzchar(year_col)) {
    stop("'year_col' must be a non-empty scalar string.", call. = FALSE)
  }
  if (!is.character(sample_id_col) || length(sample_id_col) != 1L ||
      is.na(sample_id_col) || !nzchar(sample_id_col)) {
    stop("'sample_id_col' must be a non-empty scalar string.", call. = FALSE)
  }
  if (!is.numeric(min_samples) || length(min_samples) != 1L ||
      !is.finite(min_samples) || min_samples < 1) {
    stop("'min_samples' must be a single number >= 1.", call. = FALSE)
  }

  miss_cols <- setdiff(c(year_col, sample_id_col), names(samples))
  if (length(miss_cols) > 0L) {
    stop("Missing required column(s) in 'samples': ",
         paste(miss_cols, collapse = ", "), call. = FALSE)
  }

  sample_years <- suppressWarnings(as.integer(samples[[year_col]]))
  if (any(is.na(sample_years))) {
    stop("'samples[[year_col]]' contains non-numeric/NA years.",
         call. = FALSE)
  }
  years_used <- sort(unique(sample_years))
  miss_years <- setdiff(as.character(years_used), as.character(cip_names))
  if (length(miss_years) > 0L) {
    stop("'change_index_paths' is missing entries for year(s): ",
         paste(miss_years, collapse = ", "), call. = FALSE)
  }

  if (nrow(samples) < min_samples) {
    stop("Not enough samples: nrow(samples) = ", nrow(samples),
         " < min_samples = ", min_samples, ".", call. = FALSE)
  }

  path_for_year <- function(y) {
    as.character(change_index_paths[[as.character(y)]])[1L]
  }

  # ---- per-year extraction ----------------------------------------------
  keep_vals      <- numeric(0)
  keep_medians   <- numeric(0)
  n_pixels_used  <- 0L
  pixel_area_ha  <- NA_real_

  for (y in years_used) {
    yr_idx  <- which(sample_years == y)
    polys_y <- samples[yr_idx, , drop = FALSE]
    r <- terra::rast(path_for_year(y))

    if (!is.finite(pixel_area_ha)) {
      res_xy <- terra::res(r)
      pixel_area_ha <- as.numeric(res_xy[1L] * res_xy[2L]) / 10000
    }

    crs_r <- sf::st_crs(terra::crs(r))
    if (is.na(sf::st_crs(polys_y))) {
      stop("Sample geometries for year ", y, " have an NA CRS.",
           call. = FALSE)
    }
    if (!identical(sf::st_crs(polys_y), crs_r)) {
      polys_y <- sf::st_transform(polys_y, crs_r)
    }

    msg("Year ", y, ": extracting ", nrow(polys_y), " sample polygons.")

    for (k0 in seq(1L, nrow(polys_y), by = chunk_size)) {
      k1    <- min(nrow(polys_y), k0 + chunk_size - 1L)
      chunk <- polys_y[k0:k1, , drop = FALSE]

      v <- terra::vect(chunk)
      if (!terra::same.crs(v, r)) v <- terra::project(v, terra::crs(r))

      ex <- as.data.frame(terra::extract(r, v))
      if (nrow(ex) == 0L) next

      nms    <- names(ex)
      id_col <- if ("ID" %in% nms) "ID" else nms[1L]
      val_col <- setdiff(nms, c(id_col, "coverage_fraction",
                                "fraction", "weight", "weights"))[1L]
      if (is.na(val_col)) next

      id  <- suppressWarnings(as.integer(ex[[id_col]]))
      val <- suppressWarnings(as.numeric(ex[[val_col]]))

      fin <- is.finite(val)
      finite_vals <- val[fin]
      if (length(finite_vals) > 0L) {
        n_pixels_used <- n_pixels_used + length(finite_vals)
        keep_vals <- c(keep_vals, finite_vals)
        if (length(keep_vals) > max_keep_samples) {
          keep_vals <- sample(keep_vals, max_keep_samples)
        }
      }

      ok <- is.finite(id) & is.finite(val) & id >= 1L & id <= nrow(chunk)
      if (any(ok)) {
        meds <- tapply(val[ok], id[ok], stats::median)
        keep_medians <- c(keep_medians, as.numeric(meds))
      }
    }
  }

  if (length(keep_vals) == 0L) {
    stop("Sample pixel extraction returned 0 finite values.",
         call. = FALSE)
  }
  keep_medians <- keep_medians[is.finite(keep_medians)]
  if (length(keep_medians) == 0L) {
    stop("All per-polygon medians are non-finite (samples may fall ",
         "outside the rasters).", call. = FALSE)
  }

  qref_keep <- as.numeric(stats::quantile(keep_vals, probs = qref_prob,
                                          na.rm = TRUE, names = FALSE))
  q25_keep  <- as.numeric(stats::quantile(keep_vals, probs = 0.25,
                                          na.rm = TRUE, names = FALSE))

  # Derive min_pix with the SAME logic as the scoring engine
  # (score_rbr_keep_classes(): max(ceiling(min_area_ha / pixel_area_ha), 13)).
  # When keep_pool is supplied, the engine does NOT re-derive min_pix; it
  # reads keep_pool$min_pix directly for the conf_pix check, so this must
  # be a finite value, not NULL.
  min_pix <- as.integer(max(ceiling(min_area_ha / pixel_area_ha), 13))

  msg(sprintf("qref_keep (P%.0f) = %.4f; q25_keep = %.4f; n_pix = %d; min_pix = %d",
              qref_prob * 100, qref_keep, q25_keep, n_pixels_used, min_pix))

  out <- list(
    min_area_ha         = min_area_ha,
    min_pix             = min_pix,
    promote_percentile  = promote_percentile,
    promote_p_above_ref = promote_p_above_ref,
    qref_prob           = qref_prob,
    qref_keep           = qref_keep,
    q25_keep            = q25_keep,
    pixel_area_ha       = pixel_area_ha,
    keep_medians        = keep_medians,
    metadata = list(
      n_samples          = as.integer(nrow(samples)),
      n_years            = length(years_used),
      years_used         = years_used,
      n_pixels_used      = as.integer(n_pixels_used),
      year_col           = year_col,
      sample_id_col      = sample_id_col,
      change_index_paths = change_index_paths
    )
  )
  class(out) <- c("otsufire_keep_pool", "list")
  out
}
