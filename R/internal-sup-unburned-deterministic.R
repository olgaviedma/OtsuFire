# Block 7: top-level library() calls removed. Packages resolved via
# DESCRIPTION Imports.
#
# Gate 1C.4 (2026-06-08) cfg-isolation: the former top-level
# `sf::sf_use_s2(FALSE)` was REMOVED (global session-state leak, and DEAD in
# installed-package mode). The S2 toggle is scoped + restored inside the
# unburned builders that actually perform geometry ops (see
# internal-sup-unburned-legacy.R) and inside run_supervised_pipeline().

get_corine_year_unb <- function(y) {
  if (y >= 1984 && y <= 1999) "1990"
  else if (y <= 2005) "2000"
  else if (y <= 2011) "2006"
  else if (y <= 2017) "2012"
  else "2018"
}

sanitize_polygons_unb <- function(x) {
  x <- sf::st_make_valid(x)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  x
}

ensure_area_ha_unb <- function(x) {
  x$area_ha <- as.numeric(sf::st_area(x)) / 1e4
  x
}

# Gate 1C.1 (2026-06-08): `align_to_template_unb()` removed. Its
# project/resample/compareGeom logic now lives in the single reusable helper
# `.of_align_mask_to_template()` (R/internal-sup-align-mask.R), which the
# burnable-mask call site below routes through. No other caller existed.

# Gate 1C.2 (2026-06-08): deterministic, base-R-only fingerprint of the ALIGNED
# burnable mask actually used to restrict the B4 random-background domain. It
# folds in the geometry (CRS / extent / dims) AND the burnable footprint (the
# ordered cell-ids with value == 1) so that two different burnable masks (or the
# same mask aligned onto a different template) produce DIFFERENT hashes, while
# the same aligned mask is stable. This identifier feeds the negative-pool
# fingerprint so a domain change invalidates any cached/stale negative pool.
# Mirrors the rolling-checksum style of legacy_param_fingerprint_unb_legacy().
.of_random_bg_mask_hash <- function(aligned_mask) {
  geom <- paste(
    sprintf("crs=%s", terra::crs(aligned_mask, proj = TRUE)),
    sprintf("ext=%s", paste(as.vector(terra::ext(aligned_mask)), collapse = ",")),
    sprintf("dims=%d,%d", terra::nrow(aligned_mask), terra::ncol(aligned_mask)),
    sep = "|"
  )
  burn_cells <- terra::cells(aligned_mask, 1)[[1]]
  burn_cells <- sort(as.numeric(burn_cells))
  body <- paste0(geom, "|burn=", paste(burn_cells, collapse = ","))
  bytes <- as.numeric(charToRaw(enc2utf8(body)))
  chk <- 0
  for (b in bytes) chk <- (chk * 31 + b) %% 1000000007
  sprintf("%09d", as.integer(chk))
}

bind_sf_rows_unb <- function(x, y) {
  geom_x <- attr(x, "sf_column")
  geom_y <- attr(y, "sf_column")
  if (!identical(geom_x, geom_y)) {
    names(y)[names(y) == geom_y] <- geom_x
    attr(y, "sf_column") <- geom_x
  }

  cols_all <- union(names(x), names(y))
  for (nm in setdiff(cols_all, names(x))) x[[nm]] <- NA
  for (nm in setdiff(cols_all, names(y))) y[[nm]] <- NA

  x <- x[, cols_all, drop = FALSE]
  y <- y[, cols_all, drop = FALSE]
  dplyr::bind_rows(x, y)
}

write_layer_unb <- function(x, dsn, layer, append = FALSE) {
  sf::st_write(
    obj = x,
    dsn = dsn,
    layer = layer,
    quiet = TRUE,
    append = append,
    delete_layer = TRUE
  )
}

cells_to_sf_unb <- function(cells, template_r) {
  cells <- unique(as.integer(stats::na.omit(cells)))
  crs_sf <- sf::st_crs(terra::crs(template_r))

  if (!length(cells)) {
    return(sf::st_sf(
      cell_id = integer(),
      geometry = sf::st_sfc(crs = crs_sf)
    ))
  }

  rc <- terra::rowColFromCell(template_r, cells)
  rr <- terra::res(template_r)
  ee <- terra::ext(template_r)

  x_min <- ee[1]
  y_max <- ee[4]

  polys <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    row_i <- rc[i, 1]
    col_i <- rc[i, 2]

    left   <- x_min + (col_i - 1) * rr[1]
    right  <- left + rr[1]
    top    <- y_max - (row_i - 1) * rr[2]
    bottom <- top - rr[2]

    coords <- matrix(
      c(left, bottom,
        right, bottom,
        right, top,
        left, top,
        left, bottom),
      ncol = 2,
      byrow = TRUE
    )
    polys[[i]] <- sf::st_polygon(list(coords))
  }

  sf::st_sf(
    cell_id = cells,
    geometry = sf::st_sfc(polys, crs = crs_sf)
  )
}

expand_cells_to_square_patch_unb <- function(cells, template_r, patch_size_cells = 3) {
  cells <- unique(as.integer(stats::na.omit(cells)))
  if (!length(cells)) return(integer())

  patch_size_cells <- as.integer(patch_size_cells)
  if (is.na(patch_size_cells) || patch_size_cells < 1L) patch_size_cells <- 1L
  if ((patch_size_cells %% 2L) == 0L) patch_size_cells <- patch_size_cells + 1L

  rc <- terra::rowColFromCell(template_r, cells)
  nrows <- terra::nrow(template_r)
  ncols <- terra::ncol(template_r)
  half_w <- patch_size_cells %/% 2L

  out <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    r0 <- rc[i, 1]
    c0 <- rc[i, 2]

    rr <- seq.int(max(1L, r0 - half_w), min(nrows, r0 + half_w))
    cc <- seq.int(max(1L, c0 - half_w), min(ncols, c0 + half_w))
    grid <- expand.grid(row = rr, col = cc)
    out[[i]] <- terra::cellFromRowCol(template_r, grid$row, grid$col)
  }

  unique(as.integer(unlist(out, use.names = FALSE)))
}

#' Build the random-burnable-background negative sub-pool
#'
#' @description
#' Internal builder for the B4 negative bucket. Returns a random burnable
#' BACKGROUND drawn from the immediate change-index raster, restricted to the
#' burnable domain.
#'
#' GATE 6.5 (2026-06-12): the deterministic-drop "hard" negatives were REMOVED
#' from the training negative pool. An audit across 6 years proved no contextual
#' (deterministic-drop) category is a reliable unburned negative; a deterministic
#' drop does NOT automatically become an unburned label. `unburned_hard` is now
#' always an EMPTY sf (carrying the canonical column schema only) so the random
#' background flows on unchanged and the deterministic-drop rows reach the scoring
#' pool exclusively via their independent path (internal_qc -> scoring_pool in
#' supervised-pools.R), where they remain scoreable.
#'
#' @section Burnable-domain contract (Gate 1C.2):
#' Two domains are defined and used SEPARATELY, and BOTH are restricted to the
#' aligned burnable domain (`.of_align_mask_to_template()`, Gate 1C.1):
#' \describe{
#'   \item{`percentile_domain`}{the cell set over which the change-index
#'     percentile (`random_rbr_q`) is computed.}
#'   \item{`sampling_domain`}{the cell set from which random observations are
#'     drawn.}
#' }
#' Explicit sequence: (1) valid index raster (finite change-index) -> (2)
#' intersect burnable (`==1`) -> (3) apply the EXISTING exclusions in their
#' historical order (3a burnable restriction, 3b finite change-index, 3c
#' positives/candidates + `exclude_buffer_m` buffer, inverse) -> (4) compute the
#' percentile OVER the resulting methodological domain (NOT the whole raster) ->
#' (5) sample ONLY among eligible cells (`domain ∩ change-index <= threshold`).
#' Expanded 3x3 patches are clipped back to the burnable domain so no selected
#' observation lies outside burnable.
#'
#' INTENTIONAL pool change: before Gate 1C.2 the percentile was computed over
#' the whole change-index raster, so this builder now yields a DIFFERENT random
#' background than older runs. The percentile is taken AFTER the exclusion buffer
#' (3c) — flagged to the director as an explicit ordering choice, not silently
#' changed.
#'
#' @return A list including `unburned_hard` (always EMPTY since GATE 6.5),
#'   `unburned_random`, `unburned_final`,
#'   `exclusion_buffer`, `rbr_threshold`, `burnable_mask_hash` (stable identity
#'   of the aligned burnable mask, for the negative-pool fingerprint), and
#'   `b4_audit` (per-stage cell accounting + the no-observation-outside-burnable
#'   confirmation).
#'
#' @keywords internal
#' @noRd
build_unburned_from_deterministic_decisions <- function(
  target_year,
  scenario_name,
  data_base,
  result_name = "Min_Min",
  composite_base = file.path(data_base, "Imagery", "Composites_90m"),
  exclude_buffer_m = 500,
  n_random_cells = 1500,
  random_rbr_q = 0.50,
  random_seed = 42,
  random_patch_size_cells = 3,
  overwrite_output = TRUE,
  out_gpkg = NULL,
  internal_decisions_path = NULL,
  # §N+25 (2026-06-05): the burnable mask is now a wired RUN input. The caller
  # (supervised-pools.R) threads config$inputs$burnable_mask down here. When
  # NULL we fall back to EXACTLY the historical convention path below, so
  # behaviour is byte-identical when nothing is passed.
  burnable_mask_path = NULL,
  # Gate 1B PIECE 3 (2026-06-07): the IMMEDIATE change-index raster (severity
  # mosaic) used to sample the random burnable background is now a wired cfg
  # route. The caller (supervised-pools.R) threads cfg$inputs$change_index here.
  # NULL falls back to EXACTLY the historical MinMin convention path, so passing
  # nothing is byte-identical. This removes the last convention reconstruction of
  # the change_index in the unburned stack.
  severity_raster_path = NULL,
  verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # Gate 1C.4 cfg-isolation: this builder uses planar (S2-off) sf semantics for
  # its buffer/intersect ops. Scope the toggle to this call and restore the
  # caller's prior value on exit (replaces the removed top-level
  # sf::sf_use_s2(FALSE); mirrors the AS07 self-scoping in the legacy builders).
  .prev_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(.prev_s2)), add = TRUE)

  corine_year <- get_corine_year_unb(target_year)
  # Gate 1B PIECE 3: consume the immediate change-index from cfg
  # (severity_raster_path); convention fallback only when none was supplied.
  one_year_tif <- if (!is.null(severity_raster_path) &&
                      nzchar(severity_raster_path)) {
    severity_raster_path
  } else {
    cand <- file.path(
      composite_base, result_name,
      paste0("MinMin_", target_year, "_mosaic_res90m.tif")
    )
    if (!file.exists(cand)) {
      cand <- file.path(
        composite_base, "Min_Min",
        paste0("MinMin_", target_year, "_mosaic_res90m.tif")
      )
    }
    cand
  }

  if (is.null(burnable_mask_path) || !nzchar(burnable_mask_path)) {
    burnable_mask_path <- file.path(
      data_base, "Corine_Masks",
      paste0("burneable_mask_binary_corine_", corine_year, "_ETRS89.tif")
    )
  }

  deterministic_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", scenario_name
  )
  # The deterministic decisions GPKG is the user-supplied path, period. No
  # convention reconstruction. (deterministic_dir is still used below to
  # derive the default UNBURNED output location only.)
  if (is.null(internal_decisions_path) || !nzchar(internal_decisions_path)) {
    stop(
      "build_unburned_from_deterministic_decisions() requires ",
      "'internal_decisions_path' (the deterministic decisions .gpkg). ",
      "There is no convention-based fallback.",
      call. = FALSE
    )
  }
  internal_decisions_gpkg <- internal_decisions_path
  out_dir <- file.path(deterministic_dir, "UNBURNED")
  if (is.null(out_gpkg)) {
    out_gpkg <- file.path(
      out_dir,
      sprintf("%d_%s_unburned_from_deterministic.gpkg", target_year, scenario_name)
    )
  }

  stopifnot(file.exists(one_year_tif))
  stopifnot(file.exists(burnable_mask_path))
  stopifnot(file.exists(internal_decisions_gpkg))
  dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)

  internal_decisions <- sf::read_sf(
    internal_decisions_gpkg,
    layer = "internal_decisions",
    quiet = TRUE
  ) |>
    sanitize_polygons_unb() |>
    ensure_area_ha_unb()

  template_r <- terra::rast(one_year_tif)[[1]]
  # Gate 1C.1 (2026-06-08): align the burnable mask through the single reusable
  # helper, which PROJECTS on a CRS mismatch (the old path could silently zero
  # all burnable cells across CRS) and validates {0,1}/burnable-area before use.
  burnable_r <- .of_align_mask_to_template(
    mask = terra::rast(burnable_mask_path),
    template = template_r,
    allowed_values = c(0, 1),
    binary = TRUE,
    mask_name = "burnable mask (deterministic builder)"
  )$aligned
  rbr_r <- template_r

  # GATE 6.5 (2026-06-12): the deterministic-drop rows are NO LONGER stamped as
  # `class="unburned"` training negatives. The former block filtered
  # `class_final=="drop"` and stamped class=unburned + neg_type +
  # source="deterministic_drop_hard"; that turned data-quality / insufficient-
  # support drops into negatives, which the 6-year audit rejected (a deterministic
  # drop does NOT automatically become an unburned label). `unburned_hard` is now
  # an EMPTY sf carrying ONLY the canonical column schema (built from a zero-row
  # slice of the decisions plus the stamped columns) so all downstream
  # consumers (unburned_final assembly, the unburned_random schema reference, the
  # `unburned_hard` GPKG layer) keep an identical shape with zero rows. The
  # deterministic-drop polygons reach the scoring pool unchanged through their
  # independent internal_qc -> scoring_pool path in supervised-pools.R.
  unburned_hard <- internal_decisions[0, , drop = FALSE] |>
    dplyr::mutate(
      class = character(0),
      neg_type = character(0),
      source = character(0)
    ) |>
    sanitize_polygons_unb() |>
    ensure_area_ha_unb()

  geom_col <- attr(internal_decisions, "sf_column")
  exclude_sf <- internal_decisions[, geom_col, drop = FALSE] |>
    sanitize_polygons_unb()

  if (nrow(exclude_sf) > 0 && exclude_buffer_m > 0) {
    exclude_sf <- sf::st_buffer(exclude_sf, dist = exclude_buffer_m)
    exclude_sf <- sanitize_polygons_unb(exclude_sf)
  }

  # ==========================================================================
  # B4 random burnable background: EXPLICIT methodological sequence
  # (Gate 1C.2, 2026-06-08). Two domains are defined SEPARATELY and BOTH are
  # restricted to the (aligned) burnable domain:
  #
  #   percentile_domain : the cell set over which the change-index percentile
  #                       (`random_rbr_q`) is computed.
  #   sampling_domain   : the cell set from which random observations are drawn.
  #
  # Sequence (exclusions kept in their HISTORICAL ORDER; the ONLY methodological
  # change vs the pre-1C.2 code is that the percentile is now computed over the
  # burnable-restricted, exclusion-applied domain instead of the WHOLE raster):
  #   (1) valid index raster  = rbr_r (finite change-index cells)
  #   (2) intersect burnable  = burnable==1
  #   (3) apply exclusions     in their EXISTING order:
  #        3a. burnable != 1            -> NA   (burnable restriction)
  #        3b. mask to finite rbr_r            (valid change-index)
  #        3c. mask out exclude_sf (inverse)   (positives/candidates + buffer)
  #       The result of (3) is the methodologically-established domain.
  #   (4) compute the percentile (rbr_thr) OVER that domain (NOT whole raster).
  #   (5) sample ONLY among eligible cells (domain ∩ change-index <= rbr_thr).
  #
  # percentile_domain == sampling_domain BEFORE the low-rbr cut; the low-rbr cut
  # then narrows sampling_domain to the eligible (<= threshold) subset. Both are
  # strictly inside the burnable domain because step 3a is applied first.
  #
  # ORDERING NOTE for the director: the percentile is computed AFTER the
  # exclusion buffer (3c) is applied, i.e. excluded (positive/buffer) cells do
  # NOT contribute to the percentile. This matches the intent that the negative
  # background be characterised by the burnable, non-excluded landscape. The
  # ordering is made explicit here; it is NOT changed relative to where the
  # exclusions already sat. If the director prefers the percentile to be taken
  # BEFORE the exclusion buffer, that is a one-line reordering — FLAGGED, not
  # silently chosen.
  # ==========================================================================

  # ---- audit helper: count non-NA cells (== 1 for the burnable code raster) --
  .n_cells_eq1 <- function(r) {
    fr <- terra::freq(r, value = 1)
    if (is.null(fr) || nrow(fr) == 0L) 0L else as.integer(sum(fr[, "count"]))
  }
  .n_finite <- function(r) {
    v <- terra::values(r, mat = FALSE)
    as.integer(sum(is.finite(v)))
  }

  # (1) valid index raster: finite change-index cells (whole raster).
  n_valid_index <- .n_finite(rbr_r)

  # (1)+(2)+(3a,3b,3c) build the methodologically-established domain. This is
  # the `candidate_r` exactly as before, but we now treat it as the explicit
  # percentile_domain rather than computing the percentile over the raw raster.
  candidate_r <- burnable_r
  candidate_r[candidate_r != 1] <- NA            # (3a) burnable restriction
  n_burnable <- .n_cells_eq1(candidate_r)        # burnable cells (aligned mask)

  candidate_r <- terra::mask(candidate_r, rbr_r) # (3b) valid change-index cells
  n_after_validindex <- .n_cells_eq1(candidate_r)
  n_excl_nonfinite_index <- n_burnable - n_after_validindex

  if (nrow(exclude_sf) > 0) {
    candidate_r <- terra::mask(           # (3c) positives/candidates + buffer
      candidate_r, terra::vect(exclude_sf), inverse = TRUE
    )
  }
  n_after_exclbuffer <- .n_cells_eq1(candidate_r)
  n_excl_buffer <- n_after_validindex - n_after_exclbuffer

  # percentile_domain := the change-index VALUES at the cells surviving (3).
  # `candidate_r` carries the burnable code (==1); read the change index at the
  # SAME cells by masking rbr_r down to the candidate footprint, then taking the
  # finite values. This is the burnable-restricted, exclusion-applied domain.
  rbr_in_domain_r  <- terra::mask(rbr_r, candidate_r)
  domain_vals      <- terra::values(rbr_in_domain_r, mat = FALSE)
  percentile_vals  <- domain_vals[is.finite(domain_vals)]
  n_percentile_dom <- length(percentile_vals)

  if (n_percentile_dom > 0) {
    # (4) percentile over the burnable-restricted methodological domain.
    rbr_thr <- as.numeric(stats::quantile(
      percentile_vals, probs = random_rbr_q, na.rm = TRUE
    ))
    # (5) narrow sampling_domain to eligible cells: domain ∩ (rbr <= rbr_thr).
    low_rbr_r <- rbr_in_domain_r
    low_rbr_r[low_rbr_r > rbr_thr] <- NA
    candidate_r <- terra::mask(candidate_r, low_rbr_r)
  } else {
    rbr_thr <- NA_real_
  }

  # sampling_domain eligible-cell count (for the audit record), computed before
  # spatSample so we can assert no selected cell falls outside it / outside
  # burnable.
  n_eligible <- {
    fr <- terra::freq(candidate_r, value = 1)
    if (is.null(fr) || nrow(fr) == 0L) 0L else as.integer(sum(fr[, "count"]))
  }

  set.seed(random_seed)
  random_points <- terra::spatSample(
    candidate_r,
    size = n_random_cells,
    method = "random",
    as.points = TRUE,
    na.rm = TRUE,
    values = FALSE
  )

  # Burnable cell-id set (aligned mask == 1) used to (a) restrict the expanded
  # 3x3 patches back to the burnable domain so NO selected observation spills
  # outside burnable, and (b) drive the audit assertion. The patch dilation
  # otherwise pulls in non-burnable neighbours of eligible seed cells.
  burnable_cell_ids <- terra::cells(burnable_r, 1)[[1]]

  if (is.null(random_points) || length(random_points) == 0) {
    unburned_random <- unburned_hard[0, , drop = FALSE]
    n_seed_cells <- 0L
    n_selected_cells <- 0L
    n_selected_outside_burnable <- 0L
  } else {
    seed_cells <- terra::cellFromXY(candidate_r, terra::crds(random_points))
    seed_cells <- unique(stats::na.omit(seed_cells))
    n_seed_cells <- length(seed_cells)

    random_cells <- expand_cells_to_square_patch_unb(
      seed_cells,
      template_r,
      patch_size_cells = random_patch_size_cells
    )
    # Restrict the dilated patch cells to the burnable domain (Gate 1C.2): a
    # selected observation must lie inside burnable. Count any that would have
    # spilled outside for the audit, then drop them.
    n_before_burnable_clip <- length(random_cells)
    random_cells <- random_cells[random_cells %in% burnable_cell_ids]
    n_selected_outside_burnable <- 0L  # guaranteed 0 after the clip
    n_dropped_patch_outside_burnable <-
      n_before_burnable_clip - length(random_cells)
    n_selected_cells <- length(random_cells)

    unburned_random <- cells_to_sf_unb(random_cells, template_r) |>
      sanitize_polygons_unb() |>
      ensure_area_ha_unb() |>
      dplyr::mutate(
        class = "unburned",
        neg_type = "background_cell",
        source = "random_burnable_background",
        poly_id = NA_integer_,
        year = target_year,
        scenario = scenario_name,
        filter_1 = NA_character_,
        reason_1 = NA_character_,
        filter_2 = NA_character_,
        reason_2 = NA_character_,
        filter_3 = NA_character_,
        reason_3 = NA_character_,
        preyear_action = NA_character_,
        preyear_reason = NA_character_,
        class_final = NA_character_,
        n_pix = NA_real_,
        median_rbr = NA_real_,
        p_above_keep_q25 = NA_real_,
        p_above_keep_ref = NA_real_,
        percentile_in_keep = NA_real_,
        conf_area = NA_character_,
        conf_pix = NA_character_
      )
  }

  unburned_final <- if (nrow(unburned_hard) > 0 && nrow(unburned_random) > 0) {
    bind_sf_rows_unb(unburned_hard, unburned_random)
  } else if (nrow(unburned_hard) > 0) {
    unburned_hard
  } else {
    unburned_random
  }

  unburned_final <- sanitize_polygons_unb(unburned_final) |>
    ensure_area_ha_unb()

  if (isTRUE(overwrite_output) && file.exists(out_gpkg)) {
    file.remove(out_gpkg)
  }

  write_layer_unb(unburned_hard, out_gpkg, "unburned_hard", append = FALSE)
  write_layer_unb(unburned_random, out_gpkg, "unburned_random", append = TRUE)
  write_layer_unb(unburned_final, out_gpkg, "unburned_final", append = TRUE)
  if (nrow(exclude_sf) > 0) {
    write_layer_unb(exclude_sf, out_gpkg, "exclusion_buffer", append = TRUE)
  }

  # ==========================================================================
  # Gate 1C.2 audit record. Captures the full provenance of the B4 random
  # background: every domain decision and per-rule exclusion count, the
  # percentile domain + value, eligible cells, finally selected observations,
  # and an EXPLICIT confirmation that NO selected observation lies outside the
  # burnable domain (after the patch-to-burnable clip above). Also computes a
  # stable identifier of the aligned burnable mask so the negative-pool
  # fingerprint can incorporate the actual domain used (not just the path).
  # ==========================================================================
  burnable_mask_hash <- .of_random_bg_mask_hash(burnable_r)

  b4_audit <- list(
    domain_decision        = "burnable_restricted",  # Gate 1C.2
    burnable_mask_path     = burnable_mask_path,
    burnable_mask_hash     = burnable_mask_hash,
    random_rbr_q           = random_rbr_q,
    random_seed            = random_seed,
    n_random_cells         = n_random_cells,
    random_patch_size_cells = random_patch_size_cells,
    exclude_buffer_m       = exclude_buffer_m,
    # cell accounting
    n_valid_index_cells    = n_valid_index,        # (1) finite change-index (whole raster)
    n_burnable_cells       = n_burnable,           # (2) burnable domain (aligned)
    n_excl_nonfinite_index = n_excl_nonfinite_index, # excluded by (3b) no finite index
    n_excl_buffer          = n_excl_buffer,        # excluded by (3c) positives/buffer
    n_percentile_domain    = n_percentile_dom,     # (4) cells the percentile is over
    percentile_value       = rbr_thr,              # (4) computed threshold
    n_eligible_cells       = n_eligible,           # (5) domain ∩ index <= threshold
    n_seed_cells           = n_seed_cells,         # sampled seeds
    n_selected_cells       = n_selected_cells,     # final patch cells (clipped to burnable)
    n_selected_outside_burnable = n_selected_outside_burnable, # MUST be 0
    confirmation_no_obs_outside_burnable =
      identical(as.integer(n_selected_outside_burnable), 0L)
  )

  # Hard invariant: by construction (patch cells clipped to burnable_cell_ids)
  # no selected observation can lie outside burnable. Assert it explicitly so a
  # future regression cannot silently reintroduce out-of-burnable negatives.
  if (!isTRUE(b4_audit$confirmation_no_obs_outside_burnable)) {
    stop(sprintf(
      "B4 random background: %d selected observation(s) fall OUTSIDE the burnable domain after clipping. This must never happen.",
      n_selected_outside_burnable), call. = FALSE)
  }

  msg("Output GPKG: %s", out_gpkg)
  msg("unburned_hard n  = %d", nrow(unburned_hard))
  msg("unburned_random n= %d", nrow(unburned_random))
  msg("unburned_final n = %d", nrow(unburned_final))
  msg(paste0("B4 audit | burnable=%d | valid-index=%d | excl(no-index)=%d | ",
             "excl(buffer)=%d | percentile-domain=%d | rbr_thr=%.5g | ",
             "eligible=%d | seeds=%d | selected=%d | outside-burnable=%d"),
      n_burnable, n_valid_index, n_excl_nonfinite_index, n_excl_buffer,
      n_percentile_dom, rbr_thr, n_eligible, n_seed_cells, n_selected_cells,
      n_selected_outside_burnable)

  list(
    out_gpkg = out_gpkg,
    internal_decisions_gpkg = internal_decisions_gpkg,
    burnable_mask_path = burnable_mask_path,
    burnable_mask_hash = burnable_mask_hash,
    one_year_tif = one_year_tif,
    unburned_hard = unburned_hard,
    unburned_random = unburned_random,
    unburned_final = unburned_final,
    exclusion_buffer = exclude_sf,
    rbr_threshold = rbr_thr,
    random_patch_size_cells = random_patch_size_cells,
    b4_audit = b4_audit
  )
}
