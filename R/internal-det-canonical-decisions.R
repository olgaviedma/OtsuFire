# Canonical decision outputs for the redesigned deterministic workflow.

.ensure_source_poly_id <- function(x, id_col = "source_poly_id") {
  stopifnot(inherits(x, "sf"))
  if (!id_col %in% names(x)) {
    x[[id_col]] <- seq_len(nrow(x))
  }
  x
}

.map_internal_filter_1 <- function(flag_internal) {
  x <- as.character(flag_internal)
  out <- rep(NA_character_, length(x))
  out[x == "keep"] <- "keep"
  out[x == "review"] <- "review"
  out[grepl("^drop_", x)] <- "drop"
  out
}

.map_internal_reason_1 <- function(why_flag) {
  x <- as.character(why_flag)
  out <- x
  out[grepl("^intersects_stage1_seed$", x)] <- "intersects_seed"
  out[grepl("burnable_corine_ok|no_intersection_stage1_seed", x)] <- "no_seed_support_but_burnable"
  out[grepl("outside_burnable_corine|burnable_corine_empty", x)] <- "outside_burnable"
  out[grepl("^corine_na_frac_gt_", x)] <- "corine_na_high"
  out
}

.map_reference_filter_1 <- function(flag_ref_internal) {
  x <- as.character(flag_ref_internal)
  out <- rep(NA_character_, length(x))
  out[x == "ref"] <- "keep"
  out[x == "review"] <- "review"
  out[grepl("^drop_", x)] <- "drop"
  out
}

.map_reference_reason_1 <- function(why_flag_ref) {
  x <- as.character(why_flag_ref)
  out <- x
  out[grepl("burnable_corine_ok", x)] <- "inside_burnable"
  out[grepl("outside_burnable_corine|burnable_corine_empty", x)] <- "outside_burnable"
  out[grepl("^corine_na_frac_gt_", x)] <- "corine_na_high"
  out
}

.map_reference_filter_2 <- function(flag_ref) {
  x <- as.character(flag_ref)
  out <- rep(NA_character_, length(x))
  out[x == "keep_ref"] <- "keep"
  out[x == "review_ref"] <- "review"
  out
}

.map_reference_reason_2 <- function(why_ref) {
  x <- as.character(why_ref)
  out <- x
  out[grepl("^intersects_internal_keep$", x)] <- "supported_by_internal_keep"
  out[grepl("^no_intersection_internal_keep$", x)] <- "no_internal_keep_support"
  out
}

.derive_filter_3 <- function(flag_rbr, conf_area, conf_pix, median_rbr) {
  n <- length(flag_rbr)
  filter_3 <- rep(NA_character_, n)
  reason_3 <- rep("not_applied", n)

  has_rbr <- !(is.na(flag_rbr) & is.na(median_rbr))
  if (!any(has_rbr)) {
    return(list(filter_3 = filter_3, reason_3 = reason_3))
  }

  filter_3[has_rbr] <- "review"
  reason_3[has_rbr] <- "model_unavailable"

  missing_rbr <- has_rbr & !is.finite(as.numeric(median_rbr))
  filter_3[missing_rbr] <- "review"
  reason_3[missing_rbr] <- "missing_rbr"

  too_small <- has_rbr & is.character(conf_area) & conf_area == "too_small"
  filter_3[too_small] <- "drop"
  reason_3[too_small] <- "too_small"

  too_few <- has_rbr & is.character(conf_pix) & conf_pix == "too_few_pix"
  filter_3[too_few] <- "drop"
  reason_3[too_few] <- "too_few_pixels"

  keep_like <- has_rbr & is.character(flag_rbr) & flag_rbr == "keep"
  filter_3[keep_like] <- "keep"
  reason_3[keep_like] <- "rbr_matches_model"

  borderline <- has_rbr & !too_small & !too_few & !keep_like
  filter_3[borderline] <- "review"
  reason_3[borderline] <- "rbr_borderline"

  list(filter_3 = filter_3, reason_3 = reason_3)
}

.derive_class_final <- function(filter_1, filter_2, filter_3) {
  filters <- data.frame(
    filter_1 = as.character(filter_1),
    filter_2 = as.character(filter_2),
    filter_3 = as.character(filter_3),
    stringsAsFactors = FALSE
  )

  apply(filters, 1, function(row) {
    vals <- row[!is.na(row) & nzchar(row)]
    if (!length(vals)) return(NA_character_)

    # If exactly one active filter drops the polygon but all other active filters keep it,
    # keep the polygon for review instead of discarding it outright.
    n_drop <- sum(vals == "drop")
    n_keep <- sum(vals == "keep")
    if (n_drop == 1 && (n_drop + n_keep) == length(vals) && n_keep >= 1) {
      return("review")
    }

    if (any(vals == "drop")) return("drop")
    if (all(vals == "keep")) return("keep")
    "review"
  })
}

.aggregate_area_by_id <- function(x, id_col = "source_poly_id") {
  stopifnot(inherits(x, "sf"))
  if (!nrow(x)) {
    return(data.frame(source_poly_id = numeric(0), area_m2 = numeric(0)))
  }
  data.frame(
    source_poly_id = x[[id_col]],
    area_m2 = as.numeric(sf::st_area(x)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$source_poly_id) |>
    dplyr::summarise(area_m2 = sum(.data$area_m2, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
}

.collapse_geometry_by_id <- function(x, id_col = "source_poly_id") {
  stopifnot(inherits(x, "sf"))
  x <- .ensure_source_poly_id(x, id_col = id_col)
  if (!nrow(x)) {
    return(x[, c(id_col, attr(x, "sf_column")), drop = FALSE])
  }

  x[, c(id_col, attr(x, "sf_column")), drop = FALSE] |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(.groups = "drop")
}

.empty_multipolygon_sfc <- function(n, crs) {
  sf::st_sfc(
    lapply(seq_len(n), function(i) sf::st_multipolygon()),
    crs = crs
  )
}

.prepare_reference_support <- function(internal_flagged) {
  if (is.null(internal_flagged) || !inherits(internal_flagged, "sf") || !nrow(internal_flagged)) {
    return(NULL)
  }

  df <- sf::st_drop_geometry(internal_flagged)
  if (!"source_poly_id" %in% names(df)) return(NULL)

  hit <- if ("hit_ref_any" %in% names(df)) df$hit_ref_any else NA

  data.frame(
    source_poly_id = df$source_poly_id,
    hit_ref_any = as.logical(hit),
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$source_poly_id) |>
    dplyr::summarise(hit_ref_any = any(.data$hit_ref_any %in% TRUE, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
}

.prepare_hotspot_support <- function(base_sf, hotspots_sf, id_col = "source_poly_id") {
  stopifnot(inherits(base_sf, "sf"))
  if (is.null(hotspots_sf) || !inherits(hotspots_sf, "sf") || !nrow(hotspots_sf) || !id_col %in% names(base_sf)) {
    return(NULL)
  }

  base_use <- .ensure_source_poly_id(base_sf, id_col = id_col)
  hs_use <- hotspots_sf

  if (sf::st_crs(base_use) != sf::st_crs(hs_use)) {
    hs_use <- sf::st_transform(hs_use, sf::st_crs(base_use))
  }

  hit_hotspot <- lengths(sf::st_intersects(base_use, hs_use)) > 0

  data.frame(
    source_poly_id = base_use[[id_col]],
    hit_hotspot = as.logical(hit_hotspot),
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$source_poly_id) |>
    dplyr::summarise(hit_hotspot = any(.data$hit_hotspot %in% TRUE, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
}

.join_nongeom_by_id <- function(base_sf, add_sf, id_col = "source_poly_id", prefer_cols = NULL) {
  stopifnot(inherits(base_sf, "sf"))
  if (is.null(add_sf) || !inherits(add_sf, "sf") || !nrow(add_sf) || !id_col %in% names(add_sf)) {
    return(base_sf)
  }

  add_df <- sf::st_drop_geometry(add_sf)
  keep_cols <- setdiff(names(add_df), setdiff(names(base_sf), id_col))
  if (!is.null(prefer_cols)) {
    keep_cols <- unique(c(id_col, intersect(prefer_cols, names(add_df)), setdiff(keep_cols, c(id_col, prefer_cols))))
  }
  add_df <- add_df[, unique(keep_cols), drop = FALSE]
  dup <- duplicated(add_df[[id_col]])
  if (any(dup)) {
    add_df <- add_df[!dup, , drop = FALSE]
  }
  dplyr::left_join(base_sf, add_df, by = id_col)
}

.make_gpkg_safe_local <- function(x, sep = "|") {
  stopifnot(inherits(x, "sf"))

  df <- sf::st_drop_geometry(x)
  is_list <- vapply(df, is.list, logical(1))
  if (any(is_list)) {
    for (col in names(is_list)[is_list]) {
      x[[col]] <- vapply(x[[col]], function(v) {
        if (is.null(v) || length(v) == 0) return(NA_character_)
        if (is.atomic(v)) return(paste(as.character(v), collapse = sep))
        paste(utils::capture.output(str(v)), collapse = " ")
      }, character(1))
    }
  }

  geom_col <- attr(x, "sf_column")
  nm <- names(x)
  idx_attr <- which(nm != geom_col)

  nm2 <- nm
  nm2[idx_attr] <- tolower(nm2[idx_attr])
  nm2[idx_attr] <- gsub("[^a-z0-9_]+", "_", nm2[idx_attr])
  nm2[idx_attr] <- gsub("^_+|_+$", "", nm2[idx_attr])
  nm2[idx_attr] <- ifelse(nm2[idx_attr] == "", paste0("v", idx_attr), nm2[idx_attr])
  nm2[idx_attr] <- make.unique(nm2[idx_attr], sep = "_")
  names(x) <- nm2

  for (col in names(sf::st_drop_geometry(x))) {
    if (inherits(x[[col]], c("integer64", "POSIXct", "POSIXlt", "Date"))) {
      x[[col]] <- as.character(x[[col]])
    }
  }

  x
}

.select_internal_public_columns <- function(x) {
  stopifnot(inherits(x, "sf"))
  keep <- c(
    "poly_id", "source_poly_id", "year", "scenario",
    "filter_1", "reason_1",
    "filter_2", "reason_2",
    "filter_3", "reason_3",
    "preyear_action", "preyear_reason", "preyear_overlap_frac",
    "class_final",
    "area_ha", "n_pix", "median_rbr",
    "p_above_keep_q25", "p_above_keep_ref", "percentile_in_keep",
    "conf_area", "conf_pix"
  )
  keep <- intersect(keep, names(x))
  x[, c(keep, attr(x, "sf_column")), drop = FALSE]
}

.select_reference_public_columns <- function(x) {
  stopifnot(inherits(x, "sf"))
  keep <- c(
    "poly_id", "source_poly_id", "year", "scenario",
    "filter_1", "reason_1",
    "filter_2", "reason_2",
    "filter_3", "reason_3",
    "class_final",
    "area_ha", "n_pix", "median_rbr",
    "p_above_keep_q25", "p_above_keep_ref", "percentile_in_keep",
    "conf_area", "conf_pix"
  )
  keep <- intersect(keep, names(x))
  x[, c(keep, attr(x, "sf_column")), drop = FALSE]
}

build_rbr_keep_pool_from_registry <- function(
  registry_path,
  rbr_rast,
  registry_layer = "burned_high_conf_registry",
  target_year = NULL,
  min_samples = 25L,
  verbose = TRUE,
  out_dir = file.path(tempdir(), "registry_keep_pool")
) {
  if (is.null(registry_path) || !nzchar(registry_path) || !file.exists(registry_path)) {
    return(NULL)
  }

  lyr_info <- tryCatch(sf::st_layers(registry_path), error = function(e) NULL)
  if (is.null(lyr_info) || !registry_layer %in% lyr_info$name) {
    return(NULL)
  }

  samples <- sf::st_read(registry_path, layer = registry_layer, quiet = !isTRUE(verbose))
  if (!inherits(samples, "sf") || !nrow(samples)) {
    return(NULL)
  }

  if (!is.null(target_year) && "year" %in% names(samples)) {
    samples <- samples[as.character(samples$year) != as.character(target_year), , drop = FALSE]
  }
  if (!nrow(samples) || nrow(samples) < min_samples) {
    return(NULL)
  }

  samples <- .ensure_source_poly_id(samples)
  samples$pool_flag <- "keep"

  res <- score_rbr_keep_classes(
    polys_sf = samples,
    rbr_rast = rbr_rast,
    out_dir = out_dir,
    keep_pool = NULL,
    use_keep_common_pool = FALSE,
    keep_fallback_col = "pool_flag",
    keep_fallback_value = "keep",
    save_outputs = FALSE,
    verbose = isTRUE(verbose),
    quiet = !isTRUE(verbose)
  )

  list(
    keep_pool = res$keep_pool,
    n_samples = nrow(samples),
    source = "registry",
    registry_path = registry_path,
    registry_layer = registry_layer
  )
}

build_rbr_keep_pool_from_local_support <- function(
  polys_sf,
  rbr_rast,
  hotspots_sf = NULL,
  id_col = "source_poly_id",
  min_samples = 10L,
  verbose = TRUE,
  out_dir = file.path(tempdir(), "local_keep_pool")
) {
  if (is.null(polys_sf) || !inherits(polys_sf, "sf") || !nrow(polys_sf)) {
    return(NULL)
  }

  samples <- .ensure_source_poly_id(polys_sf, id_col = id_col)
  if (!"flag_internal" %in% names(samples)) {
    return(NULL)
  }

  samples <- samples[as.character(samples$flag_internal) == "keep", , drop = FALSE]
  if (!nrow(samples)) {
    return(NULL)
  }

  hotspot_df <- .prepare_hotspot_support(samples, hotspots_sf, id_col = id_col)
  if (!is.null(hotspot_df)) {
    samples <- dplyr::left_join(samples, hotspot_df, by = id_col)
    samples_hs <- samples[samples$hit_hotspot %in% TRUE, , drop = FALSE]
    if (nrow(samples_hs)) {
      samples <- samples_hs
    }
  }

  if (!nrow(samples) || nrow(samples) < min_samples) {
    return(NULL)
  }

  samples$pool_flag <- "keep"

  res <- score_rbr_keep_classes(
    polys_sf = samples,
    rbr_rast = rbr_rast,
    out_dir = out_dir,
    keep_pool = NULL,
    use_keep_common_pool = FALSE,
    keep_fallback_col = "pool_flag",
    keep_fallback_value = "keep",
    save_outputs = FALSE,
    verbose = isTRUE(verbose),
    quiet = !isTRUE(verbose)
  )

  list(
    keep_pool = res$keep_pool,
    n_samples = nrow(samples),
    source = if (is.null(hotspot_df)) "local_seed_keep_no_hotspots" else "local_seed_keep_hotspot"
  )
}

build_internal_decisions <- function(
  internal_flagged,
  internal_clean = NULL,
  internal_rbr = NULL,
  hotspots_sf = NULL,
  target_year,
  scenario_name,
  preyear_available = FALSE
) {
  stopifnot(inherits(internal_flagged, "sf"))

  base <- .ensure_source_poly_id(internal_flagged)
  base$year <- target_year
  base$scenario <- scenario_name

  base$filter_1 <- .map_internal_filter_1(base$flag_internal)
  base$reason_1 <- .map_internal_reason_1(base$why_flag)

  base$filter_2 <- NA_character_
  base$reason_2 <- "not_applied"

  hotspot_df <- .prepare_hotspot_support(base, hotspots_sf)
  if (!is.null(hotspot_df)) {
    base <- dplyr::left_join(base, hotspot_df, by = "source_poly_id")
    has_hotspot <- !is.na(base$hit_hotspot)
    base$filter_2[has_hotspot & base$hit_hotspot] <- "keep"
    base$reason_2[has_hotspot & base$hit_hotspot] <- "intersects_hotspot"
    base$filter_2[has_hotspot & !base$hit_hotspot] <- "review"
    base$reason_2[has_hotspot & !base$hit_hotspot] <- "no_hotspot_support"
  }

  base$preyear_action <- "not_applied"
  base$preyear_reason <- "not_applied"
  if (isTRUE(preyear_available) && !is.null(internal_clean) && inherits(internal_clean, "sf") && nrow(internal_clean)) {
    raw_keep_review <- base[base$flag_internal %in% c("keep", "review"), , drop = FALSE]
    raw_area <- .aggregate_area_by_id(raw_keep_review)
    clean_area <- .aggregate_area_by_id(.ensure_source_poly_id(internal_clean))

    area_df <- dplyr::left_join(raw_area, clean_area, by = "source_poly_id", suffix = c("_raw", "_clean"))
    area_df$area_m2_clean[is.na(area_df$area_m2_clean)] <- 0
    area_df$preyear_overlap_frac <- NA_real_
    ok_area <- is.finite(area_df$area_m2_raw) & area_df$area_m2_raw > 0
    area_df$preyear_overlap_frac[ok_area] <- 1 - (area_df$area_m2_clean[ok_area] / area_df$area_m2_raw[ok_area])
    area_df$preyear_overlap_frac <- pmax(0, pmin(1, area_df$preyear_overlap_frac))

    base <- dplyr::left_join(base, area_df, by = "source_poly_id")

    if ("preyear_overlap_frac.x" %in% names(base) || "preyear_overlap_frac.y" %in% names(base)) {
      px <- if ("preyear_overlap_frac.x" %in% names(base)) suppressWarnings(as.numeric(base$preyear_overlap_frac.x)) else rep(NA_real_, nrow(base))
      py <- if ("preyear_overlap_frac.y" %in% names(base)) suppressWarnings(as.numeric(base$preyear_overlap_frac.y)) else rep(NA_real_, nrow(base))
      base$preyear_overlap_frac <- dplyr::coalesce(py, px)
      base$preyear_overlap_frac.x <- NULL
      base$preyear_overlap_frac.y <- NULL
    } else if (!"preyear_overlap_frac" %in% names(base)) {
      base$preyear_overlap_frac <- NA_real_
    }

    idx <- which(base$filter_1 %in% c("keep", "review"))
    if (length(idx)) {
      unchanged <- is.finite(base$area_m2_raw[idx]) &
        is.finite(base$area_m2_clean[idx]) &
        abs(base$area_m2_raw[idx] - base$area_m2_clean[idx]) < 1e-6
      removed <- is.finite(base$area_m2_clean[idx]) & base$area_m2_clean[idx] <= 1e-6
      reduced <- !unchanged & !removed

      base$preyear_action[idx[unchanged]] <- "keep"
      base$preyear_reason[idx[unchanged]] <- "no_overlap_previous_year"

      base$preyear_action[idx[reduced]] <- "review"
      base$preyear_reason[idx[reduced]] <- "overlap_previous_year_removed"

      base$preyear_action[idx[removed]] <- "drop"
      base$preyear_reason[idx[removed]] <- "previous_year_conflict"
    }
  } else {
    base$preyear_overlap_frac <- NA_real_
  }

  if (!is.null(internal_clean) && inherits(internal_clean, "sf")) {
    clean_geom <- .collapse_geometry_by_id(internal_clean)
    geom_idx <- match(base$source_poly_id, clean_geom$source_poly_id)
    has_clean_geom <- !is.na(geom_idx)
    geom_out <- sf::st_geometry(base)
    if (any(has_clean_geom)) {
      geom_out[has_clean_geom] <- sf::st_geometry(clean_geom)[geom_idx[has_clean_geom]]
    }
    removed_idx <- which(base$flag_internal %in% c("keep", "review") & !has_clean_geom)
    if (length(removed_idx)) {
      geom_out[removed_idx] <- .empty_multipolygon_sfc(
        length(removed_idx),
        sf::st_crs(base)
      )
    }
    sf::st_geometry(base) <- geom_out
  }

  if (!is.null(internal_rbr) && inherits(internal_rbr, "sf")) {
    base <- .join_nongeom_by_id(base, internal_rbr, prefer_cols = c(
      "n_pix", "median_rbr", "p_above_keep_q25", "p_above_keep_ref",
      "area_ha", "conf_area", "conf_pix", "percentile_in_keep",
      "keep_like", "flag_rbr"
    ))
  }

  rbr_map <- .derive_filter_3(
    flag_rbr = base$flag_rbr %||% NA_character_,
    conf_area = base$conf_area %||% NA_character_,
    conf_pix = base$conf_pix %||% NA_character_,
    median_rbr = base$median_rbr %||% NA_real_
  )

  base$filter_3 <- rbr_map$filter_3
  base$reason_3 <- rbr_map$reason_3

  idx_drop_filter1 <- which(base$filter_1 == "drop")
  if (length(idx_drop_filter1)) {
    base$filter_3[idx_drop_filter1] <- NA_character_
    base$reason_3[idx_drop_filter1] <- "not_applied"
  }

  base$class_final <- .derive_class_final(base$filter_1, base$filter_2, base$filter_3)

  base$poly_id <- seq_len(nrow(base))
  .select_internal_public_columns(base)
}

build_reference_decisions <- function(
  ref_flagged,
  ref_validation = NULL,
  ref_rbr = NULL,
  target_year,
  scenario_name
) {
  stopifnot(inherits(ref_flagged, "sf"))

  base <- .ensure_source_poly_id(ref_flagged)
  base$year <- target_year
  base$scenario <- scenario_name

  base$filter_1 <- .map_reference_filter_1(base$flag_ref_internal)
  base$reason_1 <- .map_reference_reason_1(base$why_flag_ref)

  base$filter_2 <- NA_character_
  base$reason_2 <- "not_applied"
  if (!is.null(ref_validation) && inherits(ref_validation, "sf") && nrow(ref_validation)) {
    base <- .join_nongeom_by_id(base, ref_validation, prefer_cols = c("flag_ref", "why_ref"))
    base$filter_2 <- .map_reference_filter_2(base$flag_ref)
    base$reason_2 <- .map_reference_reason_2(base$why_ref)
    base$reason_2[is.na(base$reason_2)] <- "not_applied"
  }

  if (!is.null(ref_rbr) && inherits(ref_rbr, "sf")) {
    base <- .join_nongeom_by_id(base, ref_rbr, prefer_cols = c(
      "n_pix", "median_rbr", "p_above_keep_q25", "p_above_keep_ref",
      "area_ha", "conf_area", "conf_pix", "percentile_in_keep",
      "keep_like", "flag_rbr"
    ))
  }

  rbr_map <- .derive_filter_3(
    flag_rbr = base$flag_rbr %||% NA_character_,
    conf_area = base$conf_area %||% NA_character_,
    conf_pix = base$conf_pix %||% NA_character_,
    median_rbr = base$median_rbr %||% NA_real_
  )
  base$filter_3 <- rbr_map$filter_3
  base$reason_3 <- rbr_map$reason_3

  idx_drop_filter1 <- which(base$filter_1 == "drop")
  if (length(idx_drop_filter1)) {
    base$filter_3[idx_drop_filter1] <- NA_character_
    base$reason_3[idx_drop_filter1] <- "not_applied"
  }

  base$class_final <- .derive_class_final(base$filter_1, base$filter_2, base$filter_3)
  base$poly_id <- seq_len(nrow(base))
  .select_reference_public_columns(base)
}

write_canonical_decision_output <- function(
  sf_obj,
  out_dir,
  gpkg_name,
  layer_name,
  overwrite = TRUE,
  quiet = TRUE,
  arcgis_fix = TRUE,
  output_epsg = 3035
) {
  stopifnot(inherits(sf_obj, "sf"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sf_out <- sf_obj
  if (isTRUE(arcgis_fix)) {
    sf_out <- sf::st_make_valid(sf_out)
    sf_out <- sf::st_zm(sf_out, drop = TRUE, what = "ZM")
    sf_out <- sf::st_cast(sf_out, "MULTIPOLYGON", warn = FALSE)
    if (!is.null(output_epsg) && is.finite(output_epsg)) {
      sf_out <- tryCatch(
        sf::st_transform(sf_out, output_epsg),
        error = function(e) {
          sf::st_crs(sf_out) <- output_epsg
          sf_out
        }
      )
    }
  }

  sf_out <- .make_gpkg_safe_local(sf_out)
  gpkg_path <- file.path(out_dir, gpkg_name)
  if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

  sf::st_write(sf_out, gpkg_path, layer = layer_name, driver = "GPKG", delete_dsn = TRUE, quiet = quiet)
  gpkg_path
}
