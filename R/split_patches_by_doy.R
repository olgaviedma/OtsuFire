#' Split candidate patches that mix several acquisition dates
#'
#' @description
#' Optionally split refined candidate burned patches whose pixels carry
#' clearly different post-event acquisition dates (day of year, DOY). A patch
#' that joins a real fire with haze, cloud residue or another event observed
#' weeks earlier or later is scored as a single object, so the whole mixture
#' can be accepted. Splitting it by date lets [score_burned_patches()] judge
#' each homogeneous part on its own.
#'
#' The function reads the DOY band (band 2) of the change-index raster,
#' groups the DOY values of each candidate, and replaces candidates with two
#' or more date groups by one sub-patch per connected component of each group.
#' Candidates with a single date group are returned unchanged. The original
#' refined layer is preserved.
#'
#' Run this function after [detect_burned_patches()] and pass the output to
#' [score_burned_patches()]. It is **off by default** in
#' [run_deterministic_pipeline()]; enable it with
#' `refine_params = list(split_by_doy = TRUE)` in
#' [build_burned_mapping_config()].
#'
#' @param refine_path Character scalar or `sf` polygon object. Path to the
#'   refined candidate GeoPackage, typically
#'   `02_REFINE/BA_<year>_REFINE_MERGED.gpkg` (use the
#'   `refined_patches_path` returned by [detect_burned_patches()]), or the
#'   same layer already read into memory.
#' @param change_index Character scalar path or `SpatRaster`. Change-index
#'   raster with at least two bands. Band 1 is the summer RBR and band 2 the
#'   DOY of the post-event image (`DOY_post`). Only band 2 is read.
#' @param year Integer scalar or `NULL`. Target year used in output names.
#'   When `NULL`, it is extracted from the input filename.
#' @param doy_gap_days Numeric scalar. Two date groups are kept apart only
#'   when their median DOY differ by more than this number of days. Larger
#'   values split less. Default: `30`.
#' @param min_group_share Numeric scalar in `[0, 1)`. A date group with fewer
#'   than `max(10, min_group_share * n_pixels)` pixels is absorbed into the
#'   group with the closest median DOY before the comparison. Default:
#'   `0.05`.
#' @param min_subpatch_px Integer scalar. Connected components of a date
#'   group with fewer pixels than this are dropped. Default: `9L`.
#' @param min_area_ha Numeric scalar. Only candidates with at least this area
#'   (hectares) are examined; smaller ones are returned unchanged. Default:
#'   `50`.
#' @param id_col Character scalar or `NULL`. Column with the candidate
#'   identifier, recorded in the audit table. When `NULL`, `"CLUMP_ID"` or
#'   `"poly_id"` is used if present, otherwise the row number.
#' @param out_dir Character scalar or `NULL`. Run directory under which
#'   `DOYSPLIT/dsplit_<G>d/` is created. When `NULL` and `refine_path` is a
#'   path, uses `dirname(dirname(refine_path))`. When `NULL` and `refine_path`
#'   is an `sf` object, nothing is written and the result is returned in
#'   memory.
#' @param overwrite Logical scalar. Whether existing outputs may be replaced.
#'   The original refined layer is never modified. Default: `FALSE`.
#' @param verbose Logical scalar. Whether to display progress messages.
#'   Default: `TRUE`.
#'
#' @section Grouping rule:
#' For each examined candidate, the DOY of every pixel whose centre falls
#' inside the polygon is read. Pixels with the same DOY form an initial group.
#' Groups smaller than the minimum share are absorbed into the nearest group
#' (by median DOY). Then the two groups with the closest medians are merged
#' repeatedly while their medians differ by `doy_gap_days` or less. A candidate
#' is *mixed* when two or more groups remain.
#'
#' Grouping uses the **medians** of the groups, not the gap between
#' consecutive DOY values: single linkage on consecutive values was rejected
#' because rare intermediate dates chained distinct groups together.
#'
#' @section Sub-patches:
#' The group labels are smoothed with a 3 x 3 modal filter. Each group is then
#' split into 8-connected components, and components smaller than
#' `min_subpatch_px` pixels are dropped. Sub-patch polygons follow the pixel
#' footprint on the grid of the DOY raster, derived from cell centres. Their
#' outline therefore differs slightly from the smoothed outline of the
#' original candidate, and the area of a split candidate is not exactly
#' conserved. Candidates that are not split keep their original geometry.
#'
#' @section Limitations:
#' * Only DOY separates the parts. A mixture observed on a single date is not
#'   split.
#' * Sub-patches are scored with models trained on whole patches. Their
#'   predictors (area, dispersion) change, so the scoring model is applied
#'   outside its exact training distribution.
#' * Splitting can lower recall on the edges of a real fire whose pixels were
#'   observed on dates more than `doy_gap_days` apart. In 2017 (a year without
#'   haze) a gap of 30 days left the F1 against EFFIS unchanged within 0.003
#'   and a gap of 15 days lowered it by 0.007.
#' * Output polygons carry only `CLUMP_ID` and the run metadata columns. The
#'   link between each output polygon and its parent is stored in the audit
#'   table, not in the layer.
#'
#' @section Output location:
#' Outputs are written under `<out_dir>/DOYSPLIT/dsplit_<G>d/`:
#' \itemize{
#'   \item `BA_<year>_REFINE_MERGED_dsplit<G>d.gpkg`: candidate layer for
#'     scoring.
#'   \item `doysplit_audit.csv`: one-row summary (counts, areas, parameters).
#'   \item `doysplit_patches.csv`: one row per output polygon that comes from
#'     a split, with parent identifier, date group, pixel count and DOY
#'     statistics.
#'   \item `doysplit_parents.csv`: one row per examined candidate.
#'   \item `doysplit_params.yml`: parameters used.
#' }
#'
#' @return A named list with `patches` (the output `sf` layer), `output_path`,
#'   `audit_path`, `patches_audit_path`, `parents_audit_path`, `params_path`
#'   (all `NA` when nothing is written), and the counters
#'   `n_candidates_before`, `n_candidates_after`, `n_examined`, `n_mixed`,
#'   `n_subpatches`, `n_errors` (candidates that raised an error while being
#'   split; they are kept whole and a warning is issued), `area_before_ha`, `area_after_ha`, plus the parameters
#'   `doy_gap_days`, `min_group_share`, `min_subpatch_px` and `min_area_ha`.
#'
#' @seealso [build_burned_mapping_config()], [detect_burned_patches()],
#'   [apply_chain_cleaning()], [score_burned_patches()].
#'
#' @examples
#' \dontrun{
#' detection <- detect_burned_patches(config, write_outputs = TRUE)
#'
#' split <- split_patches_by_doy(
#'   refine_path  = detection$refined_patches_path,
#'   change_index = "data/MinMin_2006_mosaic_res90m.tif",
#'   year         = 2006L,
#'   doy_gap_days = 30
#' )
#' split$n_mixed
#' split$n_subpatches
#'
#' scored <- score_burned_patches(
#'   burned_candidates = split$output_path,
#'   config = config
#' )
#' }
#'
#' @export
split_patches_by_doy <- function(refine_path,
                                 change_index,
                                 year = NULL,
                                 doy_gap_days = 30,
                                 min_group_share = 0.05,
                                 min_subpatch_px = 9L,
                                 min_area_ha = 50,
                                 id_col = NULL,
                                 out_dir = NULL,
                                 overwrite = FALSE,
                                 verbose = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE))
    stop("split_patches_by_doy() requires the 'sf' and 'terra' packages.", call. = FALSE)
  say <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## --- parameter validation ---
  .of_check_doy_split_params(doy_gap_days, min_group_share, min_subpatch_px, min_area_ha)
  min_subpatch_px <- as.integer(min_subpatch_px)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(id_col) && (!is.character(id_col) || length(id_col) != 1L))
    stop("id_col must be NULL or a single column name.", call. = FALSE)

  in_memory <- inherits(refine_path, "sf")
  if (!in_memory) {
    if (!is.character(refine_path) || length(refine_path) != 1L || !file.exists(refine_path))
      stop("refine_path must be an existing file path or an sf object.", call. = FALSE)
  }

  ## --- DOY raster (band 2) ---
  if (inherits(change_index, "SpatRaster")) {
    ci <- change_index
  } else if (is.character(change_index) && length(change_index) == 1L && file.exists(change_index)) {
    ci <- terra::rast(change_index)
  } else {
    stop("change_index must be an existing raster path or a SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(ci) < 2L)
    stop("change_index must have at least 2 bands (band 2 = DOY of the post-event image).", call. = FALSE)
  doy <- ci[[2L]]

  ## --- candidates ---
  ref <- if (in_memory) refine_path else sf::st_read(refine_path, quiet = TRUE)
  ref <- sf::st_make_valid(ref)
  ref <- ref[!sf::st_is_empty(ref), , drop = FALSE]
  if (nrow(ref) == 0L) stop("No candidates in refine_path.", call. = FALSE)
  if (is.na(sf::st_crs(ref))) stop("refine layer has no CRS.", call. = FALSE)

  if (is.null(year)) {
    nm <- if (in_memory) "" else basename(refine_path)
    m <- regmatches(nm, regexpr("(?<=BA_)[0-9]{4}", nm, perl = TRUE))
    year <- if (length(m)) as.integer(m) else NA_integer_
  }
  if (is.null(id_col)) id_col <- intersect(c("CLUMP_ID", "poly_id"), names(ref))[1L]
  ids <- if (!is.na(id_col)) as.character(ref[[id_col]]) else as.character(seq_len(nrow(ref)))

  gtag <- format(doy_gap_days, trim = TRUE)
  write_out <- !(in_memory && is.null(out_dir))
  if (write_out) {
    if (is.null(out_dir)) out_dir <- dirname(dirname(refine_path))
    sdir <- file.path(out_dir, "DOYSPLIT", sprintf("dsplit_%sd", gtag))
    yr <- ifelse(is.na(year), "NA", year)
    out_path <- file.path(sdir, sprintf("BA_%s_REFINE_MERGED_dsplit%sd.gpkg", yr, gtag))
    if (file.exists(out_path) && !isTRUE(overwrite))
      stop("Output exists and overwrite = FALSE: ", out_path, call. = FALSE)
  }

  area_ha <- as.numeric(sf::st_area(ref)) / 1e4
  examine <- which(area_ha >= min_area_ha)
  say("split_patches_by_doy: %d candidates | %d examined (>= %g ha) | gap = %s d", nrow(ref), length(examine), min_area_ha, gtag)

  ## --- DOY pixels of the examined candidates ---
  ref_r <- if (!sf::st_crs(ref) == sf::st_crs(terra::crs(doy))) sf::st_transform(ref, terra::crs(doy)) else ref
  ex_by <- list()
  if (length(examine)) {
    ex <- terra::extract(doy, terra::vect(ref_r[examine, ]), ID = TRUE, cells = TRUE)
    names(ex)[2L] <- "doy"
    ex <- ex[is.finite(ex$doy), , drop = FALSE]
    ex_by <- split(ex, ex$ID)
  }

  subs <- list(); stat <- vector("list", length(examine)); mixed_idx <- integer(0)
  n_errors <- 0L; first_error <- NULL
  for (k in seq_along(examine)) {
    i <- examine[k]
    d <- ex_by[[as.character(k)]]
    n_px <- if (is.null(d)) 0L else nrow(d)
    stat[[k]] <- data.frame(parent_id = ids[i], n_px = n_px, mixed = FALSE, n_groups = NA_integer_,
                            n_sub = 0L, kept_px_frac = NA_real_, stringsAsFactors = FALSE)
    if (n_px < 20L) next
    res <- tryCatch(.of_split_patch_by_doy(d, doy, doy_gap_days, min_group_share, min_subpatch_px),
                    error = function(e) e)
    if (inherits(res, "error")) {
      n_errors <- n_errors + 1L; first_error <- first_error %||% conditionMessage(res)
      next
    }
    if (is.null(res)) { stat[[k]]$n_groups <- 1L; next }
    pol <- res$pol
    pol$parent_id <- ids[i]
    pol$sub_id <- sprintf("%s_G%s_%d", ids[i], gtag, seq_len(nrow(pol)))
    pol$area_ha <- as.numeric(sf::st_area(pol)) / 1e4
    subs[[length(subs) + 1L]] <- pol; mixed_idx <- c(mixed_idx, i)
    stat[[k]]$mixed <- TRUE; stat[[k]]$n_groups <- res$n_groups
    stat[[k]]$n_sub <- nrow(pol); stat[[k]]$kept_px_frac <- res$n_px_kept / res$n_px_total
  }
  if (n_errors > 0L)
    warning(sprintf("split_patches_by_doy: %d candidate(s) could not be split and were kept whole. First error: %s",
                    n_errors, first_error), call. = FALSE)
  parents_audit <- if (length(stat)) do.call(rbind, stat) else
    data.frame(parent_id = character(), n_px = integer(), mixed = logical(), n_groups = integer(),
               n_sub = integer(), kept_px_frac = numeric())

  ## --- assemble: untouched candidates + sub-patches, renumber CLUMP_ID ---
  meta_cols <- intersect(c("YEAR", "RUN", "LB", "DELTA", "MINSEED", "INDEX", "WORKFLOW", "ENGINE", "SOURCE_AOI_FILE"), names(ref))
  build <- function(geom, src, ix) {
    s <- sf::st_sf(geometry = sf::st_geometry(geom))
    for (cc in meta_cols) s[[cc]] <- if (cc == "SOURCE_AOI_FILE") src else ref[[cc]][ix]
    s
  }
  keep_idx <- setdiff(seq_len(nrow(ref)), mixed_idx)
  parts <- list()
  if (length(keep_idx)) parts[["untouched"]] <- build(ref[keep_idx, ], "untouched", keep_idx)
  sub_tbl <- NULL
  if (length(subs)) {
    sub_all <- do.call(rbind, lapply(subs, function(z) z[, c("sub_id", "parent_id", "group", "n_px", "doy_med", "doy_min", "doy_max", "area_ha")]))
    sf::st_crs(sub_all) <- sf::st_crs(ref_r)
    if (!sf::st_crs(ref) == sf::st_crs(ref_r)) sub_all <- sf::st_transform(sub_all, sf::st_crs(ref))
    parent_row <- match(sub_all$parent_id, ids)
    pieces <- build(sub_all, "doy_split", parent_row)
    parts[["split"]] <- pieces
    sub_tbl <- sf::st_drop_geometry(sub_all)
  }
  n_sub_in <- if (is.null(sub_tbl)) 0L else nrow(sub_tbl)
  src <- c(rep("untouched", length(keep_idx)), rep("split", n_sub_in))
  out <- do.call(rbind, parts)
  out <- sf::st_make_valid(out)
  nonempty <- !sf::st_is_empty(out)
  out <- out[nonempty, , drop = FALSE]; src <- src[nonempty]
  out <- sf::st_cast(out, "MULTIPOLYGON", warn = FALSE)
  out[["CLUMP_ID"]] <- seq_len(nrow(out))
  out <- out[, c("CLUMP_ID", meta_cols)]
  sf::st_crs(out) <- sf::st_crs(ref)

  patches_audit <- NULL
  if (!is.null(sub_tbl)) {
    keep_sub <- nonempty[length(keep_idx) + seq_len(n_sub_in)]
    patches_audit <- cbind(CLUMP_ID = out$CLUMP_ID[src == "split"], sub_tbl[keep_sub, , drop = FALSE])
  }
  n_split_rows <- if (is.null(patches_audit)) 0L else nrow(patches_audit)

  n_before <- nrow(ref); n_after <- nrow(out)
  area_before_ha <- sum(area_ha, na.rm = TRUE)
  area_after_ha <- sum(as.numeric(sf::st_area(out)), na.rm = TRUE) / 1e4
  n_mixed <- length(mixed_idx)
  say("split_patches_by_doy DONE: %d -> %d candidates | %d mixed | %d sub-patches", n_before, n_after, n_mixed, n_split_rows)

  audit <- data.frame(
    year = year, n_candidates_before = n_before, n_candidates_after = n_after,
    n_examined = length(examine), n_mixed = n_mixed, n_subpatches = n_split_rows, n_errors = n_errors,
    area_before_ha = round(area_before_ha, 2), area_after_ha = round(area_after_ha, 2),
    doy_gap_days = doy_gap_days, min_group_share = min_group_share,
    min_subpatch_px = min_subpatch_px, min_area_ha = min_area_ha)

  paths <- list(output_path = NA_character_, audit_path = NA_character_, patches_audit_path = NA_character_,
                parents_audit_path = NA_character_, params_path = NA_character_)
  if (write_out) {
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    for (s in c("", "-wal", "-shm")) if (file.exists(paste0(out_path, s))) unlink(paste0(out_path, s), force = TRUE)
    sf::st_write(out, out_path, quiet = TRUE)
    paths$output_path <- out_path
    paths$audit_path <- file.path(sdir, "doysplit_audit.csv")
    paths$patches_audit_path <- file.path(sdir, "doysplit_patches.csv")
    paths$parents_audit_path <- file.path(sdir, "doysplit_parents.csv")
    paths$params_path <- file.path(sdir, "doysplit_params.yml")
    utils::write.csv(audit, paths$audit_path, row.names = FALSE)
    utils::write.csv(if (is.null(patches_audit)) data.frame() else patches_audit, paths$patches_audit_path, row.names = FALSE)
    utils::write.csv(parents_audit, paths$parents_audit_path, row.names = FALSE)
    writeLines(c(
      "# split_patches_by_doy parameters",
      sprintf("year: %s", year),
      sprintf("refine_path: %s", if (in_memory) "<in-memory>" else refine_path),
      sprintf("output_path: %s", out_path),
      sprintf("doy_gap_days: %g", doy_gap_days),
      sprintf("min_group_share: %g", min_group_share),
      sprintf("min_subpatch_px: %d", min_subpatch_px),
      sprintf("min_area_ha: %g", min_area_ha),
      sprintf("n_candidates_before: %d", n_before),
      sprintf("n_candidates_after: %d", n_after),
      sprintf("n_mixed: %d", n_mixed),
      sprintf("n_subpatches: %d", n_split_rows)
    ), paths$params_path)
  }

  c(list(patches = out), paths,
    list(patches_audit = patches_audit, parents_audit = parents_audit,
         n_candidates_before = n_before, n_candidates_after = n_after,
         n_examined = length(examine), n_mixed = n_mixed, n_subpatches = n_split_rows, n_errors = n_errors,
         area_before_ha = area_before_ha, area_after_ha = area_after_ha,
         doy_gap_days = doy_gap_days, min_group_share = min_group_share,
         min_subpatch_px = min_subpatch_px, min_area_ha = min_area_ha))
}

#' Validate the numeric parameters of the DOY split
#' @keywords internal
#' @noRd
.of_check_doy_split_params <- function(doy_gap_days, min_group_share, min_subpatch_px, min_area_ha) {
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!num1(doy_gap_days) || doy_gap_days <= 0)
    stop("doy_gap_days must be a single positive number.", call. = FALSE)
  if (!num1(min_group_share) || min_group_share < 0 || min_group_share >= 1)
    stop("min_group_share must be a single number in [0, 1).", call. = FALSE)
  if (!num1(min_subpatch_px) || min_subpatch_px < 1 || min_subpatch_px != round(min_subpatch_px))
    stop("min_subpatch_px must be a single whole number >= 1.", call. = FALSE)
  if (!num1(min_area_ha) || min_area_ha < 0)
    stop("min_area_ha must be a single number >= 0.", call. = FALSE)
  invisible(TRUE)
}

#' Group the DOY values of one patch (medians separated by more than `gap`)
#'
#' Groups with fewer than `max(10, min_share * n)` pixels are absorbed into
#' the group with the closest median; then the two groups with the closest
#' medians are merged while they differ by `gap` days or less.
#' @param v Numeric vector of DOY values.
#' @return Integer group labels, same length as `v`.
#' @keywords internal
#' @noRd
.of_doy_group_labels <- function(v, gap, min_share = 0.05) {
  g <- as.integer(factor(v))
  repeat {
    tab <- table(g); if (length(tab) < 2L) break
    small <- names(tab)[tab < max(10, min_share * length(g))]; if (!length(small)) break
    s1 <- small[which.min(tab[small])]; meds <- tapply(v, g, stats::median); oth <- setdiff(names(meds), s1)
    tgt <- oth[which.min(abs(meds[oth] - meds[s1]))]; g[g == as.integer(s1)] <- as.integer(tgt)
  }
  repeat {
    ids <- sort(unique(g)); if (length(ids) < 2L) break
    meds <- vapply(ids, function(k) stats::median(v[g == k]), numeric(1)); o <- order(meds)
    d <- diff(meds[o]); j <- which.min(d)
    if (d[j] > gap) break
    g[g == ids[o][j + 1L]] <- ids[o][j]
  }
  as.integer(factor(g))
}

#' Split one patch into sub-patches by DOY group
#'
#' The template grid is built from cell centres (extent = centres +/- one
#' cell), as in the validated proof of concept, so sub-patch outlines are
#' offset half a cell from the grid of the DOY raster.
#' @param d Data frame with `doy` and `cell` (cells of `doy_rast` inside the patch).
#' @return `NULL` when the patch is not mixed (or nothing survives), otherwise a
#'   list with `pol` (sf), `n_groups`, `n_px_total`, `n_px_kept`.
#' @keywords internal
#' @noRd
.of_split_patch_by_doy <- function(d, doy_rast, gap, min_share, min_px) {
  g <- .of_doy_group_labels(d$doy, gap, min_share)
  if (length(unique(g)) < 2L) return(NULL)
  xy <- terra::xyFromCell(doy_rast, d$cell); rs <- terra::res(doy_rast)[1L]
  e <- terra::ext(min(xy[, 1]) - rs, max(xy[, 1]) + rs, min(xy[, 2]) - rs, max(xy[, 2]) + rs)
  tmpl <- terra::rast(e, resolution = rs, crs = terra::crs(doy_rast)); lab <- tmpl; terra::values(lab) <- NA
  lab[terra::cellFromXY(tmpl, xy)] <- g
  lab <- terra::focal(lab, w = 3, fun = "modal", na.rm = TRUE, na.policy = "omit")
  comp <- tmpl; terra::values(comp) <- NA; off <- 0L
  for (k in sort(unique(stats::na.omit(terra::values(lab)[, 1])))) {
    pk <- terra::patches(terra::ifel(lab == k, 1, NA), directions = 8)
    m <- terra::global(pk, "max", na.rm = TRUE)[1, 1]
    if (!is.finite(m)) next
    sel <- !is.na(terra::values(pk)[, 1]); comp[sel] <- terra::values(pk)[sel, 1] + off; off <- off + m
  }
  cv <- terra::values(comp)[, 1]; if (!length(stats::na.omit(unique(cv)))) return(NULL)
  cnt <- table(cv); keep_ids <- as.integer(names(cnt)[cnt >= min_px])
  if (!length(keep_ids)) return(NULL)
  # (values are filtered in memory: %in% has no SpatRaster method unless terra is attached)
  comp2 <- comp; cv2 <- cv; cv2[!(cv2 %in% keep_ids)] <- NA; terra::values(comp2) <- cv2
  pol <- sf::st_as_sf(terra::as.polygons(comp2, dissolve = TRUE)); names(pol)[1L] <- "comp"
  cellc <- terra::cellFromXY(tmpl, xy); compc <- cv[cellc]
  st <- do.call(rbind, lapply(split(seq_along(compc), compc), function(ix) {
    data.frame(comp = as.integer(compc[ix[1L]]), group = as.integer(names(which.max(table(g[ix])))),
               n_px = length(ix), doy_med = stats::median(d$doy[ix]),
               doy_min = min(d$doy[ix]), doy_max = max(d$doy[ix]))
  }))
  st <- st[!is.na(st$comp), , drop = FALSE]
  pol <- merge(pol, st, by = "comp")
  list(pol = pol, n_groups = length(unique(g)), n_px_total = nrow(d), n_px_kept = sum(pol$n_px))
}
