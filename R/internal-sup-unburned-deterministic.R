# Block 7: top-level library() calls removed. Packages resolved via
# DESCRIPTION Imports.
sf::sf_use_s2(FALSE)

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

align_to_template_unb <- function(r, template, method = c("near", "bilinear")) {
  method <- match.arg(method)
  if (!terra::same.crs(r, template)) {
    r <- terra::project(r, template, method = method)
  }
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    r <- terra::resample(r, template, method = method)
  }
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    stop("Raster still does not match template after alignment.", call. = FALSE)
  }
  r
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
  burnable_r <- align_to_template_unb(
    terra::rast(burnable_mask_path),
    template_r,
    method = "near"
  )
  rbr_r <- template_r

  unburned_hard <- internal_decisions |>
    dplyr::filter(
      class_final == "drop"
    ) |>
    dplyr::mutate(
      class = "unburned",
      neg_type = dplyr::case_when(
        filter_1 == "drop" & reason_1 == "outside_burnable" ~ "drop_outside_burnable",
        filter_1 == "drop" & reason_1 == "corine_na_high"   ~ "geo_excluded_hot",
        filter_1 == "drop" & !is.na(reason_1) ~ paste0("drop_", reason_1),
        filter_2 == "drop" & !is.na(reason_2) ~ paste0("drop_", reason_2),
        filter_3 == "drop" & reason_3 == "rbr_far_from_model" ~ "spectral_reject_medium",
        filter_3 == "drop" & !is.na(reason_3) ~ paste0("drop_", reason_3),
        TRUE ~ "drop_hard"
      ),
      source = "deterministic_drop_hard"
    ) |>
    sanitize_polygons_unb() |>
    ensure_area_ha_unb()

  if ("poly_id" %in% names(unburned_hard)) {
    unburned_hard <- unburned_hard |>
      dplyr::distinct(poly_id, .keep_all = TRUE)
  }

  geom_col <- attr(internal_decisions, "sf_column")
  exclude_sf <- internal_decisions[, geom_col, drop = FALSE] |>
    sanitize_polygons_unb()

  if (nrow(exclude_sf) > 0 && exclude_buffer_m > 0) {
    exclude_sf <- sf::st_buffer(exclude_sf, dist = exclude_buffer_m)
    exclude_sf <- sanitize_polygons_unb(exclude_sf)
  }

  candidate_r <- burnable_r
  candidate_r[candidate_r != 1] <- NA
  candidate_r <- terra::mask(candidate_r, rbr_r)

  if (nrow(exclude_sf) > 0) {
    candidate_r <- terra::mask(candidate_r, terra::vect(exclude_sf), inverse = TRUE)
  }

  r_vals <- terra::values(rbr_r, mat = FALSE)
  ok_vals <- r_vals[is.finite(r_vals)]
  if (length(ok_vals) > 0) {
    rbr_thr <- as.numeric(stats::quantile(ok_vals, probs = random_rbr_q, na.rm = TRUE))
    low_rbr_r <- rbr_r
    low_rbr_r[low_rbr_r > rbr_thr] <- NA
    candidate_r <- terra::mask(candidate_r, low_rbr_r)
  } else {
    rbr_thr <- NA_real_
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

  if (is.null(random_points) || length(random_points) == 0) {
    unburned_random <- unburned_hard[0, , drop = FALSE]
  } else {
    random_cells <- terra::cellFromXY(candidate_r, terra::crds(random_points))
    random_cells <- unique(stats::na.omit(random_cells))
    random_cells <- expand_cells_to_square_patch_unb(
      random_cells,
      template_r,
      patch_size_cells = random_patch_size_cells
    )

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

  msg("Output GPKG: %s", out_gpkg)
  msg("unburned_hard n  = %d", nrow(unburned_hard))
  msg("unburned_random n= %d", nrow(unburned_random))
  msg("unburned_final n = %d", nrow(unburned_final))

  list(
    out_gpkg = out_gpkg,
    internal_decisions_gpkg = internal_decisions_gpkg,
    burnable_mask_path = burnable_mask_path,
    one_year_tif = one_year_tif,
    unburned_hard = unburned_hard,
    unburned_random = unburned_random,
    unburned_final = unburned_final,
    exclusion_buffer = exclude_sf,
    rbr_threshold = rbr_thr,
    random_patch_size_cells = random_patch_size_cells
  )
}
