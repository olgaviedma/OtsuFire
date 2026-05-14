# ============================================================
# Phase-4 diagnostics: helpers + low-level functions + wrapper
# (histograms + two-sample tests, optionally by class)
# ============================================================

# ----------------------------
# Internal helpers
# ----------------------------

# small helper: "x %||% y"
`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) == 1) {
    xx <- as.character(x)
    if (!is.na(xx) && nzchar(xx)) return(xx)
  }
  y
}

# Extract raster values for a chunk of polygons, returning columns: ID, val
# (ID is local to the input vector chunk: 1..nrow(chunk))
#' @keywords internal
extract_vals_chunk <- function(rbr_rast, v_chunk) {
  stopifnot(inherits(rbr_rast, "SpatRaster"), inherits(v_chunk, "SpatVector"))

  # Use first layer if multi-band
  rr <- rbr_rast[[1]]

  ex <- terra::extract(rr, v_chunk)
  if (is.null(ex) || nrow(ex) == 0) {
    return(data.frame(ID = integer(0), val = numeric(0)))
  }

  val_col <- setdiff(names(ex), "ID")[1]
  if (is.na(val_col) || !nzchar(val_col)) {
    return(data.frame(ID = integer(0), val = numeric(0)))
  }

  data.frame(
    ID  = as.integer(ex$ID),
    val = suppressWarnings(as.numeric(ex[[val_col]])),
    stringsAsFactors = FALSE
  )
}

# Add polygon median of RBR (per polygon)
#' @keywords internal
add_polygon_median_rbr <- function(polys_sf, rbr_rast, chunk_size = 500, colname = "median_rbr") {
  stopifnot(inherits(polys_sf, "sf"), inherits(rbr_rast, "SpatRaster"))

  if (colname %in% names(polys_sf) && any(is.finite(suppressWarnings(as.numeric(polys_sf[[colname]]))))) {
    return(polys_sf)
  }

  v <- terra::vect(polys_sf)
  if (!terra::same.crs(v, rbr_rast)) v <- terra::project(v, terra::crs(rbr_rast))

  polys_sf[[colname]] <- NA_real_

  idx_all <- seq_len(nrow(polys_sf))
  groups <- split(idx_all, ceiling(seq_along(idx_all) / chunk_size))

  for (g in groups) {
    ex <- extract_vals_chunk(rbr_rast, v[g])
    if (nrow(ex) == 0) next

    spl <- split(ex$val, ex$ID)  # ID = 1..length(g)
    for (id_local in names(spl)) {
      i_global <- g[as.integer(id_local)]
      vals <- spl[[id_local]]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next
      polys_sf[[colname]][i_global] <- stats::median(vals)
    }
  }

  polys_sf
}

# Quantile vertical lines for histogram overlays
#' @keywords internal
quantile_lines <- function(x_ref) {
  x_ref <- x_ref[is.finite(x_ref)]
  if (length(x_ref) < 10) return(NULL)
  as.numeric(stats::quantile(x_ref, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
}

# Overlaid histogram with Y as percent of polygons
#' @keywords internal
save_overlaid_hist_percent <- function(x1, x2, name1, name2, png_path,
                                       breaks_n = 40,
                                       xlab = "Median RBR",
                                       ylab = "% of polygons",
                                       main = NULL,
                                       vlines = NULL) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 5 || length(x2) < 5) return(invisible(FALSE))

  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  rng <- range(c(x1, x2), na.rm = TRUE)
  brks <- pretty(rng, n = breaks_n)

  h1 <- graphics::hist(x1, breaks = brks, plot = FALSE)
  h2 <- graphics::hist(x2, breaks = brks, plot = FALSE)

  pct1 <- 100 * h1$counts / max(1, sum(h1$counts))
  pct2 <- 100 * h2$counts / max(1, sum(h2$counts))
  ymax <- max(c(pct1, pct2), na.rm = TRUE)

  grDevices::png(png_path, width = 1400, height = 900, res = 140)
  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(op); grDevices::dev.off() }, add = TRUE)

  if (is.null(main)) main <- paste("Overlaid histogram:", name1, "vs", name2)

  col1 <- grDevices::adjustcolor("black", alpha.f = 0.25)
  col2 <- grDevices::adjustcolor("gray40", alpha.f = 0.25)

  graphics::plot(NA, xlim = rng, ylim = c(0, ymax * 1.05),
                 xlab = xlab, ylab = ylab, main = main)

  graphics::rect(h1$breaks[-length(h1$breaks)], 0, h1$breaks[-1], pct1, col = col1, border = "black")
  graphics::rect(h2$breaks[-length(h2$breaks)], 0, h2$breaks[-1], pct2, col = col2, border = "gray40")

  if (!is.null(vlines) && length(vlines) > 0) {
    for (vl in vlines) graphics::abline(v = vl, lty = 2)
  }

  graphics::legend("topright",
                   legend = c(paste0(name1, " (n=", length(x1), ")"),
                              paste0(name2, " (n=", length(x2), ")")),
                   fill = c(col1, col2),
                   border = c("black", "gray40"),
                   bty = "n")

  invisible(TRUE)
}

# ----------------------------
# Low-level: histograms
# ----------------------------

#' Save Phase-4 overlaid histograms (single variable)
#'
#' Low-level function to generate overlaid percent histograms comparing key groups:
#' `keep_common`, `review_internal`, and `review_external` (reference review).
#' Optionally creates class-based histograms (same label in internal vs ref) and
#' pairwise class histograms (with a manifest CSV).
#'
#' @param internal_sf,ref_sf `sf` polygon layers (or lists with `scored_sf`).
#' @param rbr_rast Optional `terra::SpatRaster`. Required only if `value_col="median_rbr"`
#'   is missing/empty and must be computed.
#' @param out_dir Output folder for PNGs (created if needed).
#' @param buffer_m Numeric. Buffer (meters) for `ref_keep` in geometry fallback.
#' @param chunk_size Integer. Chunk size for `add_polygon_median_rbr()`.
#' @param value_col Character. Numeric column to compare.
#' @param internal_keep_col,internal_review_col Character. Columns defining keep/review.
#' @param internal_keep_value,internal_review_value Character. Keep/review values.
#' @param ref_keep_col,ref_review_col Character. Columns defining keep/review in reference.
#' @param ref_keep_value,ref_review_value Character. Keep/review values in reference.
#' @param use_internal_effis_flag Logical. Prefer `keep_common_value` in `internal_effis_col`.
#' @param internal_effis_col Character. Cross-flag column in internal.
#' @param keep_common_value Character. Label for common keep polygons.
#' @param review_int_no_effis_value Character. Label for internal review with no reference.
#' @param make_extra_no_effis Logical. Make extra overlays using review_int_no_effis.
#' @param compare_class3 Logical. Run class comparisons.
#' @param class_col_internal,class_col_ref Class column names in internal/ref.
#' @param class_levels_internal,class_levels_ref Optional level filters.
#' @param class_pairwise Logical. Pairwise class histograms.
#' @param class_pairwise_mode `"cross_dataset"` or `"all_pairs"`.
#' @param class_pairwise_min_n Minimum finite values per group for pairwise plots (NULL uses `min_n`).
#' @param class_pairwise_max_pairs Safety cap on number of pairwise plots.
#' @param class_pairwise_include_within_dataset Logical. Used when `mode="all_pairs"`.
#' @param min_n Minimum finite values per group to plot.
#' @param verbose Logical.
#'
#' @return Named list of output paths (PNGs + optional manifest CSV).
#' @keywords internal
  phase4_histograms_internal_external <- function(
    internal_sf,
    ref_sf,
    rbr_rast = NULL,
    out_dir,
    buffer_m = 0,
    chunk_size = 500,

    value_col = "median_rbr",

    internal_keep_col   = "flag_internal",
    internal_keep_value = "keep",
    internal_review_col = "flag_internal",
    internal_review_value = "review",

    ref_keep_col   = "flag_ref",
    ref_keep_value = "keep_ref",
    ref_review_col = "flag_ref",
    ref_review_value = "review_ref",

    use_internal_effis_flag = TRUE,
    internal_effis_col = "flag_internal_effis",
    keep_common_value = "keep_common",
    review_int_no_effis_value = "review_int_no_effis",
    make_extra_no_effis = FALSE,

    compare_class3 = TRUE,
    class_col_internal = "class3_id",
    class_col_ref      = "class3_id",
    class_levels_internal = NULL,
    class_levels_ref      = NULL,

    class_pairwise = TRUE,
    class_pairwise_mode = c("cross_dataset", "all_pairs"),
    class_pairwise_min_n = NULL,
    class_pairwise_max_pairs = Inf,
    class_pairwise_include_within_dataset = FALSE,

    min_n = 5,
    verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.list(internal_sf) && "scored_sf" %in% names(internal_sf)) internal_sf <- internal_sf$scored_sf
  if (is.list(ref_sf)      && "scored_sf" %in% names(ref_sf))      ref_sf      <- ref_sf$scored_sf

  stopifnot(inherits(internal_sf, "sf"), inherits(ref_sf, "sf"))

  class_pairwise_mode <- match.arg(class_pairwise_mode)

  # ---- Ensure numeric column exists (or compute median_rbr if needed) --------
  need_compute_median <- identical(value_col, "median_rbr") && (
    (!"median_rbr" %in% names(internal_sf) || sum(is.finite(suppressWarnings(as.numeric(internal_sf$median_rbr)))) == 0) ||
      (!"median_rbr" %in% names(ref_sf)    || sum(is.finite(suppressWarnings(as.numeric(ref_sf$median_rbr)))) == 0)
  )

  if (isTRUE(need_compute_median)) {
    if (is.null(rbr_rast)) stop("value_col='median_rbr' missing/empty and rbr_rast is NULL.")
    stopifnot(inherits(rbr_rast, "SpatRaster"))
    internal_sf <- add_polygon_median_rbr(internal_sf, rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
    ref_sf      <- add_polygon_median_rbr(ref_sf,      rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
  } else {
    if (!value_col %in% names(internal_sf) || !value_col %in% names(ref_sf)) {
      stop("Missing value_col in data. value_col='", value_col, "' must exist in both internal_sf and ref_sf.")
    }
  }

  # Cast to numeric safely
  internal_sf[[value_col]] <- suppressWarnings(as.numeric(internal_sf[[value_col]]))
  ref_sf[[value_col]]      <- suppressWarnings(as.numeric(ref_sf[[value_col]]))

  # CRS align for geometry fallback (keep_common)
  if (sf::st_crs(ref_sf) != sf::st_crs(internal_sf)) {
    ref_sf <- sf::st_transform(ref_sf, sf::st_crs(internal_sf))
  }

  # ---- Group checks ----------------------------------------------------------
  if (!internal_keep_col %in% names(internal_sf)) stop("internal_sf missing internal_keep_col: ", internal_keep_col)
  if (!internal_review_col %in% names(internal_sf)) stop("internal_sf missing internal_review_col: ", internal_review_col)
  if (!ref_keep_col %in% names(ref_sf)) stop("ref_sf missing ref_keep_col: ", ref_keep_col)
  if (!ref_review_col %in% names(ref_sf)) stop("ref_sf missing ref_review_col: ", ref_review_col)

  internal_keep   <- internal_sf[as.character(internal_sf[[internal_keep_col]]) == internal_keep_value, ]
  internal_review <- internal_sf[as.character(internal_sf[[internal_review_col]]) == internal_review_value, ]

  ref_keep   <- ref_sf[as.character(ref_sf[[ref_keep_col]]) == ref_keep_value, ]
  ref_review <- ref_sf[as.character(ref_sf[[ref_review_col]]) == ref_review_value, ]

  if (nrow(internal_keep) == 0) stop("No internal keep polygons (internal_keep_col/value).")
  if (nrow(ref_keep) == 0) msg("Warning: No ref keep polygons; keep_common fallback-by-geometry may be empty.")
  if (nrow(ref_review) == 0) msg("Warning: No ref review polygons; external overlays may be empty.")

  # ---- keep_common -----------------------------------------------------------
  keep_common <- internal_keep

  if (isTRUE(use_internal_effis_flag) &&
      internal_effis_col %in% names(internal_sf) &&
      any(as.character(internal_sf[[internal_effis_col]]) == keep_common_value, na.rm = TRUE)) {

    keep_common <- internal_sf[as.character(internal_sf[[internal_effis_col]]) == keep_common_value, ]

  } else {
    # Geometry-based keep_common (fallback): internal_keep intersecting ref_keep
    if (nrow(ref_keep) > 0) {
      ik <- sf::st_make_valid(internal_keep)
      rk <- sf::st_make_valid(ref_keep)
      if (is.numeric(buffer_m) && buffer_m > 0) rk <- sf::st_buffer(rk, dist = buffer_m)
      hit_ref_keep <- lengths(sf::st_intersects(ik, rk)) > 0
      keep_common <- internal_keep[hit_ref_keep, ]
      if (nrow(keep_common) == 0) {
        msg("Warning: keep_common by intersects is empty; fallback to internal_keep.")
        keep_common <- internal_keep
      }
    } else {
      keep_common <- internal_keep
    }
  }

  # ---- review_internal -------------------------------------------------------
  review_internal_label <- "review_internal"
  if (nrow(internal_review) == 0 &&
      internal_effis_col %in% names(internal_sf) &&
      any(as.character(internal_sf[[internal_effis_col]]) == review_int_no_effis_value, na.rm = TRUE)) {

    internal_review <- internal_sf[as.character(internal_sf[[internal_effis_col]]) == review_int_no_effis_value, ]
    review_internal_label <- "review_int_no_effis"
    msg("Note: internal review empty; using '%s' as review_internal.", review_int_no_effis_value)
  }

  keep_common_vals     <- keep_common[[value_col]]
  review_internal_vals <- internal_review[[value_col]]
  review_external_vals <- ref_review[[value_col]]

  ok_vec <- function(x, n = min_n) sum(is.finite(x)) >= n

  # vertical lines: Q25/Q50/Q75 of reference group (keep_common)
  v_keep <- if (ok_vec(keep_common_vals)) quantile_lines(keep_common_vals) else NULL
  v_revI <- if (ok_vec(review_internal_vals)) quantile_lines(review_internal_vals) else NULL

  out_paths <- list()

  safe_tag <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9_]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x
  }

  short_tag <- function(x, max_len = 24) {
    x <- safe_tag(x)
    if (nchar(x) > max_len) substr(x, 1, max_len) else x
  }

  vtag <- short_tag(value_col, max_len = 24)

  # 1) keep_common vs review_internal
  if (ok_vec(keep_common_vals) && ok_vec(review_internal_vals)) {
    p1 <- file.path(out_dir, sprintf("hist_%s_keep_common_vs_%s.png", vtag, short_tag(review_internal_label, 20)))
    ok1 <- save_overlaid_hist_percent(
      x1 = keep_common_vals, x2 = review_internal_vals,
      name1 = "keep_common", name2 = review_internal_label,
      png_path = p1,
      main = sprintf("%s: keep_common vs %s", value_col, review_internal_label),
      vlines = v_keep
    )
    if (isTRUE(ok1)) out_paths$p1 <- p1
  } else {
    msg("Skipping p1: insufficient data in keep_common or review_internal.")
  }

  # 2) keep_common vs review_external
  if (ok_vec(keep_common_vals) && ok_vec(review_external_vals)) {
    p2 <- file.path(out_dir, sprintf("hist_%s_keep_common_vs_review_external.png", vtag))
    ok2 <- save_overlaid_hist_percent(
      x1 = keep_common_vals, x2 = review_external_vals,
      name1 = "keep_common", name2 = "review_external",
      png_path = p2,
      main = sprintf("%s: keep_common vs review_external (ref)", value_col),
      vlines = v_keep
    )
    if (isTRUE(ok2)) out_paths$p2 <- p2
  } else {
    msg("Skipping p2: insufficient data in keep_common or review_external.")
  }

  # 3) review_internal vs review_external
  if (ok_vec(review_internal_vals) && ok_vec(review_external_vals)) {
    p3 <- file.path(out_dir, sprintf("hist_%s_%s_vs_review_external.png", vtag, short_tag(review_internal_label, 20)))
    ok3 <- save_overlaid_hist_percent(
      x1 = review_internal_vals, x2 = review_external_vals,
      name1 = review_internal_label, name2 = "review_external",
      png_path = p3,
      main = sprintf("%s: %s vs review_external (ref)", value_col, review_internal_label),
      vlines = v_revI
    )
    if (isTRUE(ok3)) out_paths$p3 <- p3
  } else {
    msg("Skipping p3: insufficient data in review_internal or review_external.")
  }

  # ---- Optional extras: review_int_no_effis overlays -------------------------
  if (isTRUE(make_extra_no_effis) &&
      internal_effis_col %in% names(internal_sf) &&
      any(as.character(internal_sf[[internal_effis_col]]) == review_int_no_effis_value, na.rm = TRUE) &&
      review_internal_label != "review_int_no_effis") {

    review_no_effis <- internal_sf[as.character(internal_sf[[internal_effis_col]]) == review_int_no_effis_value, ]
    review_no_effis_vals <- review_no_effis[[value_col]]

    if (ok_vec(keep_common_vals) && ok_vec(review_no_effis_vals)) {
      p4 <- file.path(out_dir, sprintf("hist_%s_keep_common_vs_review_int_no_effis.png", vtag))
      ok4 <- save_overlaid_hist_percent(
        x1 = keep_common_vals, x2 = review_no_effis_vals,
        name1 = "keep_common", name2 = "review_int_no_effis",
        png_path = p4,
        main = sprintf("%s: keep_common vs review_int_no_effis", value_col),
        vlines = v_keep
      )
      if (isTRUE(ok4)) out_paths$p4 <- p4
    }

    if (ok_vec(review_no_effis_vals) && ok_vec(review_external_vals)) {
      p5 <- file.path(out_dir, sprintf("hist_%s_review_int_no_effis_vs_review_external.png", vtag))
      ok5 <- save_overlaid_hist_percent(
        x1 = review_no_effis_vals, x2 = review_external_vals,
        name1 = "review_int_no_effis", name2 = "review_external",
        png_path = p5,
        main = sprintf("%s: review_int_no_effis vs review_external (ref)", value_col),
        vlines = if (ok_vec(review_no_effis_vals)) quantile_lines(review_no_effis_vals) else NULL
      )
      if (isTRUE(ok5)) out_paths$p5 <- p5
    }
  }

  # ============================================================
  # CLASS (basic): same label in both datasets
  # ============================================================
  if (isTRUE(compare_class3)) {

    if (!class_col_internal %in% names(internal_sf)) {
      msg("Skipping class comparisons: internal missing class_col_internal='%s'.", class_col_internal)
    } else if (!class_col_ref %in% names(ref_sf)) {
      msg("Skipping class comparisons: ref missing class_col_ref='%s'.", class_col_ref)
    } else {

      internal_cls <- as.character(internal_sf[[class_col_internal]])
      ref_cls      <- as.character(ref_sf[[class_col_ref]])

      lv_int <- sort(unique(stats::na.omit(internal_cls)))
      lv_ref <- sort(unique(stats::na.omit(ref_cls)))

      if (!is.null(class_levels_internal)) lv_int <- intersect(lv_int, class_levels_internal)
      if (!is.null(class_levels_ref))      lv_ref <- intersect(lv_ref, class_levels_ref)

      lv_common <- intersect(lv_int, lv_ref)

      if (length(lv_common) == 0) {
        msg("Skipping same-class internal vs ref: no common labels between internal('%s') and ref('%s').",
            class_col_internal, class_col_ref)
      } else {

        coltag <- short_tag(paste0(class_col_internal, "_", class_col_ref), 18)

        for (lv in lv_common) {
          ii <- internal_sf[as.character(internal_sf[[class_col_internal]]) == lv, , drop = FALSE]
          rr <- ref_sf[as.character(ref_sf[[class_col_ref]]) == lv, , drop = FALSE]

          xi <- ii[[value_col]]
          xr <- rr[[value_col]]

          if (ok_vec(xi) && ok_vec(xr)) {
            p <- file.path(out_dir, sprintf("hist_%s_same_%s_%s_int_vs_ref.png", vtag, coltag, short_tag(lv, 18)))
            ok <- save_overlaid_hist_percent(
              x1 = xi, x2 = xr,
              name1 = paste0("internal:", lv), name2 = paste0("ref:", lv),
              png_path = p,
              main = sprintf("%s: internal vs ref (label=%s)", value_col, lv),
              vlines = quantile_lines(xi)
            )
            if (isTRUE(ok)) out_paths[[paste0("class_same_", safe_tag(lv))]] <- p
          } else {
            msg("Skipping same-class '%s': insufficient data (min_n=%s).", lv, min_n)
          }
        }
      }
    }
  }

  # ============================================================
  # CLASS (pairwise): pairs of (dataset x class) combinations
  # ============================================================
  if (isTRUE(compare_class3) && isTRUE(class_pairwise)) {

    n_min_pair <- if (is.null(class_pairwise_min_n)) min_n else class_pairwise_min_n

    if (!class_col_internal %in% names(internal_sf) || !class_col_ref %in% names(ref_sf)) {
      msg("Skipping class pairwise: missing class column in internal or ref (internal='%s', ref='%s').",
          class_col_internal, class_col_ref)
    } else {

      mk_groups <- function(sfobj, class_col, dataset_tag, lvl_filter = NULL) {
        cc <- as.character(sfobj[[class_col]])
        lvls <- sort(unique(stats::na.omit(cc)))
        if (!is.null(lvl_filter)) lvls <- intersect(lvls, lvl_filter)

        out <- list()
        for (lv in lvls) {
          sub <- sfobj[as.character(sfobj[[class_col]]) == lv, , drop = FALSE]
          out[[paste0(dataset_tag, ":", lv)]] <- sub[[value_col]]
        }
        out
      }

      groups_internal <- mk_groups(internal_sf, class_col_internal, "internal", class_levels_internal)
      groups_ref      <- mk_groups(ref_sf,      class_col_ref,      "ref",      class_levels_ref)

      groups_all <- c(groups_internal, groups_ref)

      keep_group <- vapply(groups_all, function(x) sum(is.finite(x)) >= n_min_pair, logical(1))
      groups_all <- groups_all[keep_group]

      if (length(groups_all) < 2) {
        msg("Skipping class pairwise: <2 groups with >=%s finite values.", n_min_pair)
      } else {

        keys <- names(groups_all)

        parse_key <- function(k) {
          sp <- strsplit(k, ":", fixed = TRUE)[[1]]
          list(dataset = sp[1], label = paste(sp[-1], collapse = ":"))
        }

        pairs <- list()

        if (class_pairwise_mode == "cross_dataset") {
          keys_i <- keys[grepl("^internal:", keys)]
          keys_r <- keys[grepl("^ref:", keys)]
          for (ki in keys_i) {
            for (kr in keys_r) {
              pairs[[length(pairs) + 1L]] <- c(ki, kr)
            }
          }
        } else {
          # all_pairs
          nkeys <- length(keys)
          if (nkeys >= 2) {
            for (a in seq_len(nkeys - 1L)) {
              for (b in (a + 1L):nkeys) {
                ka <- keys[a]; kb <- keys[b]
                pa <- parse_key(ka); pb <- parse_key(kb)
                if (!isTRUE(class_pairwise_include_within_dataset) && pa$dataset == pb$dataset) next
                pairs[[length(pairs) + 1L]] <- c(ka, kb)
              }
            }
          }
        }

        if (length(pairs) == 0) {
          msg("Skipping class pairwise: no valid pairs after filtering.")
        } else {

          if (is.finite(class_pairwise_max_pairs) && length(pairs) > class_pairwise_max_pairs) {
            msg("Capping class pairwise plots: %s -> %s pairs (class_pairwise_max_pairs).",
                length(pairs), class_pairwise_max_pairs)
            pairs <- pairs[seq_len(class_pairwise_max_pairs)]
          }

          coltag <- short_tag(paste0(class_col_internal, "_", class_col_ref), 18)
          manifest_path <- file.path(out_dir, sprintf("manifest_pairs_%s_%s_%s.csv", vtag, coltag, class_pairwise_mode))

          manifest_rows <- list()

          for (i in seq_along(pairs)) {
            k1 <- pairs[[i]][1]
            k2 <- pairs[[i]][2]
            x1 <- groups_all[[k1]]
            x2 <- groups_all[[k2]]

            if (!ok_vec(x1, n = n_min_pair) || !ok_vec(x2, n = n_min_pair)) next

            png_path <- file.path(out_dir, sprintf("hist_%s_pair_%s_%03d.png", vtag, coltag, i))

            ok <- save_overlaid_hist_percent(
              x1 = x1, x2 = x2,
              name1 = k1, name2 = k2,
              png_path = png_path,
              main = sprintf("%s: %s vs %s", value_col, k1, k2),
              vlines = quantile_lines(x1)
            )

            if (isTRUE(ok)) {
              out_paths[[sprintf("pair_%03d", i)]] <- png_path
              manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
                pair_id = i,
                group1 = k1,
                group2 = k2,
                n1 = sum(is.finite(x1)),
                n2 = sum(is.finite(x2)),
                png = png_path,
                stringsAsFactors = FALSE
              )
            }
          }

          if (length(manifest_rows) > 0) {
            manifest <- do.call(rbind, manifest_rows)
            utils::write.csv(manifest, manifest_path, row.names = FALSE)
            out_paths$pairwise_manifest <- manifest_path
            msg("Wrote pairwise manifest: %s", manifest_path)
          } else {
            msg("Pairwise manifest not written (no plots created after filters).")
          }
        }
      }
    }
  }

  out_paths
}

#' Save Phase-4 overlaid histograms (multiple variables)
#'
#' Convenience wrapper around \code{phase4_histograms_internal_external()}
#' to iterate over `value_cols`. If `median_rbr` is needed and missing/empty, it is
#' computed once and reused.
#'
#' @inheritParams phase4_histograms_internal_external
#' @param value_cols Character vector of numeric columns to analyze.
#'
#' @return Named list (by `value_col`) of output-path lists.
#' @keywords internal
save_phase4_histograms_multi_vars <- function(
    internal_sf,
    ref_sf,
    rbr_rast = NULL,
    out_dir,
    value_cols = c("median_rbr"),

    buffer_m = 0,
    chunk_size = 500,

    internal_keep_col   = "flag_internal",
    internal_keep_value = "keep",
    internal_review_col = "flag_internal",
    internal_review_value = "review",

    ref_keep_col   = "flag_ref",
    ref_keep_value = "keep_ref",
    ref_review_col = "flag_ref",
    ref_review_value = "review_ref",

    use_internal_effis_flag = TRUE,
    internal_effis_col = "flag_internal_effis",
    keep_common_value = "keep_common",
    review_int_no_effis_value = "review_int_no_effis",
    make_extra_no_effis = FALSE,

    compare_class3 = TRUE,

    class_col_internal = "class3_id",
    class_col_ref      = "class3_id",
    class_levels_internal = NULL,
    class_levels_ref      = NULL,

    class_pairwise = TRUE,
    class_pairwise_mode = c("cross_dataset", "all_pairs"),
    class_pairwise_min_n = 20,
    class_pairwise_max_pairs = Inf,
    class_pairwise_include_within_dataset = FALSE,

    min_n = 5,
    verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  if (is.list(internal_sf) && "scored_sf" %in% names(internal_sf)) internal_sf <- internal_sf$scored_sf
  if (is.list(ref_sf)      && "scored_sf" %in% names(ref_sf))      ref_sf      <- ref_sf$scored_sf

  stopifnot(inherits(internal_sf, "sf"), inherits(ref_sf, "sf"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  class_pairwise_mode <- match.arg(class_pairwise_mode)

  need_median <- "median_rbr" %in% value_cols && (
    (!"median_rbr" %in% names(internal_sf) || sum(is.finite(suppressWarnings(as.numeric(internal_sf$median_rbr)))) == 0) ||
      (!"median_rbr" %in% names(ref_sf)    || sum(is.finite(suppressWarnings(as.numeric(ref_sf$median_rbr)))) == 0)
  )

  if (isTRUE(need_median)) {
    if (is.null(rbr_rast)) stop("Need median_rbr but rbr_rast is NULL.")
    stopifnot(inherits(rbr_rast, "SpatRaster"))
    msg("Computing median_rbr once for internal/ref (reused across variables).")
    internal_sf <- add_polygon_median_rbr(internal_sf, rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
    ref_sf      <- add_polygon_median_rbr(ref_sf,      rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
  }

  out_all <- list()

  for (vcol in value_cols) {
    msg("---- Histograms for value_col = '%s' -> %s", vcol, out_dir)

    out_all[[vcol]] <- phase4_histograms_internal_external2(
      internal_sf = internal_sf,
      ref_sf      = ref_sf,
      rbr_rast    = rbr_rast,
      out_dir     = out_dir,
      buffer_m    = buffer_m,
      chunk_size  = chunk_size,

      value_col   = vcol,

      internal_keep_col   = internal_keep_col,
      internal_keep_value = internal_keep_value,
      internal_review_col = internal_review_col,
      internal_review_value = internal_review_value,

      ref_keep_col   = ref_keep_col,
      ref_keep_value = ref_keep_value,
      ref_review_col = ref_review_col,
      ref_review_value = ref_review_value,

      use_internal_effis_flag = use_internal_effis_flag,
      internal_effis_col = internal_effis_col,
      keep_common_value = keep_common_value,
      review_int_no_effis_value = review_int_no_effis_value,
      make_extra_no_effis = make_extra_no_effis,

      compare_class3 = compare_class3,

      class_col_internal = class_col_internal,
      class_col_ref      = class_col_ref,
      class_levels_internal = class_levels_internal,
      class_levels_ref      = class_levels_ref,

      class_pairwise = class_pairwise,
      class_pairwise_mode = class_pairwise_mode,
      class_pairwise_min_n = class_pairwise_min_n,
      class_pairwise_max_pairs = class_pairwise_max_pairs,
      class_pairwise_include_within_dataset = class_pairwise_include_within_dataset,

      min_n = min_n,
      verbose = verbose
    )
  }

  out_all
}




# ----------------------------
# Low-level: tests
# ----------------------------

#' Two-sample distribution test helper (Wilcoxon, optional KS, Cliff's delta)
#'
#' Computes a two-sample comparison between numeric vectors:
#' Wilcoxon p-value, optional KS p-value, and Cliff's delta (via Mann-Whitney U).
#'
#' @param v1,v2 Numeric vectors.
#' @param name1,name2 Character labels used in the output `pair`.
#' @param min_n Integer. Minimum sample size per vector to run tests.
#' @param do_ks Logical. If TRUE, attempts a KS test (may fail for ties).
#'
#' @return A list with counts, medians, and test statistics/p-values.
#' @keywords internal
test_two_samples <- function(v1, v2, name1, name2, min_n = 10, do_ks = TRUE) {
  v1 <- v1[is.finite(v1)]
  v2 <- v2[is.finite(v2)]

  out <- list(
    pair = paste(name1, "vs", name2),
    n1 = length(v1),
    n2 = length(v2),
    median1 = if (length(v1) > 0) stats::median(v1) else NA_real_,
    median2 = if (length(v2) > 0) stats::median(v2) else NA_real_,
    wilcox_p = NA_real_,
    ks_p = NA_real_,
    cliff_delta = NA_real_
  )

  if (length(v1) < min_n || length(v2) < min_n) return(out)

  wil <- suppressWarnings(stats::wilcox.test(v1, v2, exact = FALSE))
  out$wilcox_p <- unname(wil$p.value)

  n1 <- length(v1); n2 <- length(v2)
  r  <- rank(c(v1, v2), ties.method = "average")
  W1 <- sum(r[seq_len(n1)])
  U1 <- W1 - n1 * (n1 + 1) / 2
  out$cliff_delta <- as.numeric((2 * U1) / (n1 * n2) - 1)

  if (isTRUE(do_ks)) {
    ks <- suppressWarnings(tryCatch(stats::ks.test(v1, v2), error = function(e) NULL))
    if (!is.null(ks)) out$ks_p <- unname(ks$p.value)
  }

  out
}

#' Phase-4 distribution tests for internal vs reference polygons
#'
#' Computes two-sample comparisons between `keep_common`, `review_internal`,
#' `review_external` (reference review), and optional `review_int_no_effis`.
#' Optionally performs class-based comparisons.
#'
#' @param internal_sf,ref_sf `sf` polygon layers (or lists with `scored_sf`).
#' @param rbr_rast Optional `terra::SpatRaster`. Only required if `"median_rbr"` must be computed.
#' @param buffer_m Numeric. Buffer (meters) for `ref_keep` in geometry fallback.
#' @param chunk_size Integer. Chunk size for `add_polygon_median_rbr()`.
#' @param value_col Character. Numeric column to test.
#' @param ensure_median_rbr Logical. If TRUE and `value_col=="median_rbr"`, computes
#'   `median_rbr` if missing/empty (requires `rbr_rast`).
#' @param use_internal_effis_flag,internal_effis_col,keep_common_value,review_int_no_effis_value
#'   Logic for defining `keep_common` and `review_int_no_effis`.
#' @param internal_flag_col Character or NULL. If NULL, auto-detects `flag_internal` then `flag_rbr`.
#' @param internal_keep_value,internal_review_value Character. Internal keep/review values.
#' @param ref_flag_col Character. Reference flag column (default `"flag_ref"`).
#' @param ref_keep_value,ref_review_value Character. Reference keep/review values.
#' @param compare_class3 Logical. If TRUE, run class-based tests.
#' @param class_col_internal,class_col_ref Character. Class columns in internal/ref.
#' @param class_levels_internal,class_levels_ref Optional character vectors restricting class levels.
#' @param class_pairwise Which class comparisons to run. See function definition.
#' @param drop_na_class Logical. Drop NA/empty class labels before testing.
#' @param min_n Integer. Minimum finite values per group to run tests.
#' @param do_ks Logical. If TRUE, attempts a KS test (may fail for ties).
#' @param save_outputs Logical. If TRUE, writes CSV/TXT outputs to `out_dir`.
#' @param out_dir Output folder when saving.
#' @param prefix Filename prefix for saved outputs.
#' @param verbose Logical.
#'
#' @return A list with data.frames for base and class comparisons plus output paths.
#' @keywords internal
test_phase4_distributions2 <- function(
    internal_sf,
    ref_sf,
    rbr_rast = NULL,
    buffer_m = 0,
    chunk_size = 500,

    value_col = "median_rbr",
    ensure_median_rbr = TRUE,

    use_internal_effis_flag = TRUE,
    internal_effis_col = "flag_internal_effis",
    keep_common_value = "keep_common",
    review_int_no_effis_value = "review_int_no_effis",

    internal_flag_col = NULL,
    internal_keep_value = "keep",
    internal_review_value = "review",
    ref_flag_col = "flag_ref",
    ref_keep_value = "keep_ref",
    ref_review_value = "review_ref",

    compare_class3 = TRUE,
    class_col_internal = "class3_id",
    class_col_ref      = "class3_id",
    class_levels_internal = NULL,
    class_levels_ref      = NULL,
    class_pairwise = c("all_pairs_internal_vs_ref",
                       "same_class_internal_vs_ref",
                       "within_internal",
                       "within_ref",
                       "all"),
    drop_na_class = TRUE,

    min_n = 10,
    do_ks = TRUE,

    save_outputs = FALSE,
    out_dir = NULL,
    prefix = "phase4",
    verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  class_pairwise <- match.arg(class_pairwise)

  if (is.list(internal_sf) && "scored_sf" %in% names(internal_sf)) internal_sf <- internal_sf$scored_sf
  if (is.list(ref_sf)      && "scored_sf" %in% names(ref_sf))      ref_sf      <- ref_sf$scored_sf

  stopifnot(inherits(internal_sf, "sf"), inherits(ref_sf, "sf"))

  if (isTRUE(save_outputs)) {
    if (is.null(out_dir) || !nzchar(out_dir)) stop("save_outputs=TRUE requires out_dir.")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # auto-detect internal_flag_col if missing
  if (is.null(internal_flag_col)) {
    if ("flag_internal" %in% names(internal_sf)) internal_flag_col <- "flag_internal"
    else if ("flag_rbr" %in% names(internal_sf)) internal_flag_col <- "flag_rbr"
    else stop("internal_sf must have flag_internal or flag_rbr (or pass internal_flag_col).")
  }

  if (!internal_flag_col %in% names(internal_sf)) stop("internal_flag_col not found in internal_sf.")
  if (!ref_flag_col %in% names(ref_sf)) stop("ref_flag_col not found in ref_sf (expected flag_ref by default).")

  # optional: ensure median_rbr exists
  if (isTRUE(ensure_median_rbr) && identical(value_col, "median_rbr")) {
    if (is.null(rbr_rast)) stop("ensure_median_rbr=TRUE requires rbr_rast when value_col='median_rbr'.")
    stopifnot(inherits(rbr_rast, "SpatRaster"))

    if (!"median_rbr" %in% names(internal_sf) || sum(is.finite(suppressWarnings(as.numeric(internal_sf$median_rbr)))) == 0) {
      internal_sf <- add_polygon_median_rbr(internal_sf, rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
    }
    if (!"median_rbr" %in% names(ref_sf) || sum(is.finite(suppressWarnings(as.numeric(ref_sf$median_rbr)))) == 0) {
      ref_sf <- add_polygon_median_rbr(ref_sf, rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
    }
  }

  if (!value_col %in% names(internal_sf)) stop("value_col not found in internal_sf: ", value_col)
  if (!value_col %in% names(ref_sf)) stop("value_col not found in ref_sf: ", value_col)

  internal_sf <- sf::st_make_valid(internal_sf)
  ref_sf      <- sf::st_make_valid(ref_sf)

  # CRS align (only needed for geometry fallback keep_common)
  if (sf::st_crs(ref_sf) != sf::st_crs(internal_sf)) {
    ref_sf <- sf::st_transform(ref_sf, sf::st_crs(internal_sf))
  }

  # -------------------------
  # BASE GROUPS
  # -------------------------
  internal_keep   <- internal_sf[as.character(internal_sf[[internal_flag_col]]) == internal_keep_value, ]
  internal_review <- internal_sf[as.character(internal_sf[[internal_flag_col]]) == internal_review_value, ]

  ref_keep   <- ref_sf[as.character(ref_sf[[ref_flag_col]]) == ref_keep_value, ]
  ref_review <- ref_sf[as.character(ref_sf[[ref_flag_col]]) == ref_review_value, ]

  if (nrow(internal_keep) == 0) stop("No internal keep polygons found.")
  if (nrow(ref_review) == 0) warning("No ref review polygons found (ref_review is empty).")

  # keep_common (prefer Phase1b label if present)
  keep_common <- internal_keep
  review_int_no_effis <- internal_sf[0, , drop = FALSE]

  if (isTRUE(use_internal_effis_flag) && internal_effis_col %in% names(internal_sf)) {
    if (any(internal_sf[[internal_effis_col]] == keep_common_value, na.rm = TRUE)) {
      keep_common <- internal_sf[internal_sf[[internal_effis_col]] == keep_common_value, ]
    } else {
      # geometry fallback: internal_keep intersecting ref_keep
      if (nrow(ref_keep) > 0 && nrow(internal_keep) > 0) {
        ik <- sf::st_make_valid(internal_keep)
        rk <- sf::st_make_valid(ref_keep)
        if (is.numeric(buffer_m) && buffer_m > 0) rk <- sf::st_buffer(rk, dist = buffer_m)
        hit <- lengths(sf::st_intersects(ik, rk)) > 0
        keep_common <- internal_keep[hit, ]
        if (nrow(keep_common) == 0) keep_common <- internal_keep
      }
    }

    if (any(internal_sf[[internal_effis_col]] == review_int_no_effis_value, na.rm = TRUE)) {
      review_int_no_effis <- internal_sf[internal_sf[[internal_effis_col]] == review_int_no_effis_value, ]
    }
  }

  v_keep_common         <- suppressWarnings(as.numeric(keep_common[[value_col]]))
  v_review_internal     <- suppressWarnings(as.numeric(internal_review[[value_col]]))
  v_review_int_no_effis <- suppressWarnings(as.numeric(review_int_no_effis[[value_col]]))
  v_review_external     <- suppressWarnings(as.numeric(ref_review[[value_col]]))

  # run base comparisons
  base_results <- list(
    keep_common_vs_review_internal =
      test_two_samples(v_keep_common, v_review_internal,
                       "keep_common", "review_internal", min_n = min_n, do_ks = do_ks),

    keep_common_vs_review_int_no_effis =
      test_two_samples(v_keep_common, v_review_int_no_effis,
                       "keep_common", "review_int_no_effis", min_n = min_n, do_ks = do_ks),

    keep_common_vs_review_external =
      test_two_samples(v_keep_common, v_review_external,
                       "keep_common", "review_external", min_n = min_n, do_ks = do_ks),

    review_internal_vs_review_external =
      test_two_samples(v_review_internal, v_review_external,
                       "review_internal", "review_external", min_n = min_n, do_ks = do_ks),

    review_int_no_effis_vs_review_external =
      test_two_samples(v_review_int_no_effis, v_review_external,
                       "review_int_no_effis", "review_external", min_n = min_n, do_ks = do_ks)
  )

  base_df <- do.call(rbind, lapply(base_results, function(r) {
    data.frame(
      var = value_col,
      comparison = r$pair,
      n1 = r$n1, n2 = r$n2,
      median1 = r$median1, median2 = r$median2,
      cliff_delta = r$cliff_delta,
      wilcox_p = r$wilcox_p,
      ks_p = r$ks_p,
      stringsAsFactors = FALSE
    )
  }))

  # Holm on base family
  base_df$wilcox_p_adj_holm <- NA_real_
  base_df$ks_p_adj_holm     <- NA_real_
  w_ok <- is.finite(base_df$wilcox_p)
  k_ok <- is.finite(base_df$ks_p)
  if (any(w_ok)) base_df$wilcox_p_adj_holm[w_ok] <- stats::p.adjust(base_df$wilcox_p[w_ok], method = "holm")
  if (any(k_ok)) base_df$ks_p_adj_holm[k_ok]     <- stats::p.adjust(base_df$ks_p[k_ok], method = "holm")

  # -------------------------
  # CLASS COMPARISONS
  # -------------------------
  results_by_class_df <- data.frame(stringsAsFactors = FALSE)
  results_by_class_list <- list()

  if (isTRUE(compare_class3)) {
    if (!class_col_internal %in% names(internal_sf)) stop("class_col_internal not found in internal_sf: ", class_col_internal)
    if (!class_col_ref %in% names(ref_sf)) stop("class_col_ref not found in ref_sf: ", class_col_ref)

    int_cls <- as.character(internal_sf[[class_col_internal]])
    ref_cls <- as.character(ref_sf[[class_col_ref]])

    if (isTRUE(drop_na_class)) {
      internal_sf2 <- internal_sf[!is.na(int_cls) & nzchar(int_cls), , drop = FALSE]
      ref_sf2      <- ref_sf[!is.na(ref_cls) & nzchar(ref_cls), , drop = FALSE]
    } else {
      internal_sf2 <- internal_sf
      ref_sf2      <- ref_sf
    }

    int_levels <- sort(unique(as.character(internal_sf2[[class_col_internal]])))
    ref_levels <- sort(unique(as.character(ref_sf2[[class_col_ref]])))

    if (!is.null(class_levels_internal)) int_levels <- intersect(int_levels, class_levels_internal)
    if (!is.null(class_levels_ref))      ref_levels <- intersect(ref_levels, class_levels_ref)

    run_one <- function(v1, v2, nm1, nm2, family, c1 = NA_character_, c2 = NA_character_) {
      r <- test_two_samples(v1, v2, nm1, nm2, min_n = min_n, do_ks = do_ks)
      data.frame(
        var = value_col,
        family = family,
        internal_class = c1,
        ref_class = c2,
        comparison = r$pair,
        n1 = r$n1, n2 = r$n2,
        median1 = r$median1, median2 = r$median2,
        cliff_delta = r$cliff_delta,
        wilcox_p = r$wilcox_p,
        ks_p = r$ks_p,
        stringsAsFactors = FALSE
      )
    }

    dfs <- list()

    if (class_pairwise %in% c("same_class_internal_vs_ref", "all")) {
      common <- intersect(int_levels, ref_levels)
      if (length(common) > 0) {
        tmp <- do.call(rbind, lapply(common, function(cc) {
          v1 <- suppressWarnings(as.numeric(internal_sf2[as.character(internal_sf2[[class_col_internal]]) == cc, , drop = FALSE][[value_col]]))
          v2 <- suppressWarnings(as.numeric(ref_sf2[as.character(ref_sf2[[class_col_ref]]) == cc, , drop = FALSE][[value_col]]))
          run_one(v1, v2,
                  paste0("internal_", cc), paste0("ref_", cc),
                  family = "same_class_internal_vs_ref",
                  c1 = cc, c2 = cc)
        }))
        dfs[["same_class_internal_vs_ref"]] <- tmp
      }
    }

    if (class_pairwise %in% c("all_pairs_internal_vs_ref", "all")) {
      tmp <- do.call(rbind, lapply(int_levels, function(ci) {
        do.call(rbind, lapply(ref_levels, function(cr) {
          v1 <- suppressWarnings(as.numeric(internal_sf2[as.character(internal_sf2[[class_col_internal]]) == ci, , drop = FALSE][[value_col]]))
          v2 <- suppressWarnings(as.numeric(ref_sf2[as.character(ref_sf2[[class_col_ref]]) == cr, , drop = FALSE][[value_col]]))
          run_one(v1, v2,
                  paste0("internal_", ci), paste0("ref_", cr),
                  family = "all_pairs_internal_vs_ref",
                  c1 = ci, c2 = cr)
        }))
      }))
      dfs[["all_pairs_internal_vs_ref"]] <- tmp
    }

    if (class_pairwise %in% c("within_internal", "all")) {
      if (length(int_levels) >= 2) {
        comb <- utils::combn(int_levels, 2, simplify = FALSE)
        tmp <- do.call(rbind, lapply(comb, function(cc) {
          a <- cc[1]; b <- cc[2]
          v1 <- suppressWarnings(as.numeric(internal_sf2[as.character(internal_sf2[[class_col_internal]]) == a, , drop = FALSE][[value_col]]))
          v2 <- suppressWarnings(as.numeric(internal_sf2[as.character(internal_sf2[[class_col_internal]]) == b, , drop = FALSE][[value_col]]))
          run_one(v1, v2,
                  paste0("internal_", a), paste0("internal_", b),
                  family = "within_internal",
                  c1 = a, c2 = b)
        }))
        dfs[["within_internal"]] <- tmp
      }
    }

    if (class_pairwise %in% c("within_ref", "all")) {
      if (length(ref_levels) >= 2) {
        comb <- utils::combn(ref_levels, 2, simplify = FALSE)
        tmp <- do.call(rbind, lapply(comb, function(cc) {
          a <- cc[1]; b <- cc[2]
          v1 <- suppressWarnings(as.numeric(ref_sf2[as.character(ref_sf2[[class_col_ref]]) == a, , drop = FALSE][[value_col]]))
          v2 <- suppressWarnings(as.numeric(ref_sf2[as.character(ref_sf2[[class_col_ref]]) == b, , drop = FALSE][[value_col]]))
          run_one(v1, v2,
                  paste0("ref_", a), paste0("ref_", b),
                  family = "within_ref",
                  c1 = a, c2 = b)
        }))
        dfs[["within_ref"]] <- tmp
      }
    }

    if (length(dfs) > 0) {
      results_by_class_df <- do.call(rbind, dfs)

      # Holm correction per family (separately)
      results_by_class_df$wilcox_p_adj_holm <- NA_real_
      results_by_class_df$ks_p_adj_holm     <- NA_real_

      for (fam in unique(results_by_class_df$family)) {
        idx <- which(results_by_class_df$family == fam)
        w_ok <- is.finite(results_by_class_df$wilcox_p[idx])
        k_ok <- is.finite(results_by_class_df$ks_p[idx])

        if (any(w_ok)) {
          results_by_class_df$wilcox_p_adj_holm[idx[w_ok]] <-
            stats::p.adjust(results_by_class_df$wilcox_p[idx[w_ok]], method = "holm")
        }
        if (any(k_ok)) {
          results_by_class_df$ks_p_adj_holm[idx[k_ok]] <-
            stats::p.adjust(results_by_class_df$ks_p[idx[k_ok]], method = "holm")
        }
      }

      results_by_class_list <- dfs
    }
  }

  # -------------------------
  # SAVE
  # -------------------------
  out_paths <- list(csv_base = NA_character_, txt_base = NA_character_,
                    csv_by_class = NA_character_, txt_by_class = NA_character_)

  if (isTRUE(save_outputs)) {
    csv_base <- file.path(out_dir, sprintf("%s_tests_base.csv", prefix))
    utils::write.csv(base_df, csv_base, row.names = FALSE)
    out_paths$csv_base <- csv_base

    if (nrow(results_by_class_df) > 0) {
      csv_cls <- file.path(out_dir, sprintf("%s_tests_by_class.csv", prefix))
      utils::write.csv(results_by_class_df, csv_cls, row.names = FALSE)
      out_paths$csv_by_class <- csv_cls
    }

    txt_base <- file.path(out_dir, sprintf("%s_tests_base.txt", prefix))
    writeLines(c(
      "PHASE 4 TESTS (base)",
      sprintf("value_col: %s", value_col),
      sprintf("min_n: %s", min_n),
      "",
      paste(capture.output(print(base_df, row.names = FALSE)), collapse = "\n")
    ), txt_base)
    out_paths$txt_base <- txt_base

    if (nrow(results_by_class_df) > 0) {
      txt_cls <- file.path(out_dir, sprintf("%s_tests_by_class.txt", prefix))
      writeLines(c(
        "PHASE 4 TESTS (by class)",
        sprintf("value_col: %s", value_col),
        sprintf("class_col_internal: %s", class_col_internal),
        sprintf("class_col_ref: %s", class_col_ref),
        sprintf("class_pairwise: %s", class_pairwise),
        sprintf("min_n: %s", min_n),
        "",
        paste(capture.output(print(results_by_class_df, row.names = FALSE)), collapse = "\n")
      ), txt_cls)
      out_paths$txt_by_class <- txt_cls
    }

    msg("Wrote: %s", out_paths$csv_base)
    if (!is.na(out_paths$csv_by_class)) msg("Wrote: %s", out_paths$csv_by_class)
  }

  list(
    results_df = base_df,
    results_list = base_results,
    results_by_class_df = results_by_class_df,
    results_by_class_list = results_by_class_list,
    out_paths = out_paths,
    sizes = list(
      internal_keep = nrow(internal_keep),
      internal_review = nrow(internal_review),
      keep_common = nrow(keep_common),
      review_int_no_effis = nrow(review_int_no_effis),
      ref_keep = nrow(ref_keep),
      ref_review = nrow(ref_review)
    )
  )
}

# ----------------------------
# Exported wrapper
# ----------------------------

#' Phase-4 diagnostics for internal vs reference burned-area polygons
#'
#' Runs a combined diagnostic workflow to compare numeric distributions between an
#' internal polygon layer (your map) and a reference polygon layer (e.g., EFFIS):
#' \itemize{
#'   \item Overlaid percentage histograms for key group comparisons (optional).
#'   \item Two-sample distribution tests (Wilcoxon; optional KS; Cliff's delta),
#'         including optional class-stratified comparisons (optional).
#' }
#'
#' The function accepts `sf` objects (polygons) or list outputs that contain a
#' `scored_sf` element (it will automatically extract it).
#'
#' If `value_cols` includes `"median_rbr"` and that column is missing or has no
#' finite values in either dataset, the function computes it once using
#' `add_polygon_median_rbr()` (requires `rbr_rast`).
#'
#' @param internal_sf `sf` polygon layer (internal), or a list with `scored_sf`.
#' @param ref_sf `sf` polygon layer (reference), or a list with `scored_sf`.
#' @param rbr_rast Optional `terra::SpatRaster`. Required only if `"median_rbr"` must be computed.
#' @param out_dir Output folder (base). Subfolders are created when needed.
#' @param value_cols Character vector of numeric columns to analyze.
#' @param do_histograms Logical. If TRUE, histogram PNGs are generated.
#' @param do_tests Logical. If TRUE, two-sample tests are computed.
#'
#' @param hist_dir Folder for histogram outputs. Default: `file.path(out_dir, "phase4_histograms")`.
#' @param tests_dir Folder for test outputs. Default: `file.path(out_dir, "phase4_tests")`.
#'
#' @param buffer_m Numeric. Buffer (meters) applied to `ref_keep` only for the geometry
#'   fallback that derives `keep_common` from `internal_keep` intersecting `ref_keep`.
#' @param chunk_size Integer. Chunk size passed to `add_polygon_median_rbr()`.
#'
#' @param internal_flag_col Character. Flag column used to define keep/review groups in the
#'   internal dataset. If NULL, it is auto-detected (`flag_internal`, then `flag_rbr`).
#' @param internal_keep_value,internal_review_value Character. Keep/review labels in internal flags.
#' @param ref_flag_col Character. Flag column in the reference dataset.
#' @param ref_keep_value,ref_review_value Character. Keep/review labels in reference flags.
#'
#' @param use_internal_effis_flag Logical. If TRUE and `internal_effis_col` exists,
#'   `keep_common` prefers `keep_common_value` from that column.
#' @param internal_effis_col Character. Column holding Phase-1b / cross-flag labels.
#' @param keep_common_value Character. Label for common keep polygons.
#' @param review_int_no_effis_value Character. Label for "review internal, no reference".
#' @param make_extra_no_effis Logical. If TRUE, extra histogram overlays are generated using
#'   `review_int_no_effis_value` when present.
#'
#' @param compare_class3 Logical. If TRUE, run class-based comparisons.
#' @param class_col_internal,class_col_ref Character. Class columns in internal/ref.
#' @param class_levels_internal,class_levels_ref Optional character vectors restricting class levels (NULL = all).
#'
#' @param hist_class_pairwise Logical. If TRUE, generate pairwise class histograms.
#' @param hist_class_pairwise_mode One of `"cross_dataset"` or `"all_pairs"`.
#' @param hist_class_pairwise_min_n Integer. Minimum finite values per group for pairwise class histograms.
#' @param hist_class_pairwise_max_pairs Integer. Safety cap for number of pairwise plots.
#' @param hist_class_pairwise_include_within_dataset Logical. Only used when mode is `"all_pairs"`.
#'
#' @param test_class_pairwise One of: `"all_pairs_internal_vs_ref"`, `"same_class_internal_vs_ref"`,
#'   `"within_internal"`, `"within_ref"`, or `"all"`.
#' @param drop_na_class Logical. If TRUE, drops NA/empty class labels before testing.
#'
#' @param min_n_hist Integer. Minimum finite values per group to draw a histogram overlay.
#' @param min_n_test Integer. Minimum finite values per group to run tests.
#' @param do_ks Logical. If TRUE, compute KS p-values (in addition to Wilcoxon).
#'
#' @param save_test_outputs Logical. If TRUE, writes CSV/TXT summaries to `tests_dir`.
#' @param test_prefix Character. Prefix for saved test filenames.
#' @param return_sf Logical. If TRUE, returns updated `internal_sf` and `ref_sf` (useful when `median_rbr` was computed).
#' @param verbose Logical.
#'
#' @return A list with:
#' \itemize{
#'   \item `hist_paths`: named list by `value_col` (each an output-path list), or empty.
#'   \item `tests`: named list by `value_col` (each a `test_phase4_distributions2()` result), or empty.
#'   \item `out_dirs`: list with `hist_dir` and `tests_dir`.
#'   \item `internal_sf` / `ref_sf`: only when `return_sf=TRUE`.
#' }
#'
#' @export
phase4_diagnostics_internal_external <- function(
    internal_sf,
    ref_sf,
    rbr_rast = NULL,
    out_dir,
    value_cols = c("median_rbr"),

    do_histograms = TRUE,
    do_tests = TRUE,

    hist_dir = file.path(out_dir, "phase4_histograms"),
    tests_dir = file.path(out_dir, "phase4_tests"),

    buffer_m = 0,
    chunk_size = 500,

    internal_flag_col = NULL,
    internal_keep_value = "keep",
    internal_review_value = "review",

    ref_flag_col = "flag_ref",
    ref_keep_value = "keep_ref",
    ref_review_value = "review_ref",

    use_internal_effis_flag = TRUE,
    internal_effis_col = "flag_internal_effis",
    keep_common_value = "keep_common",
    review_int_no_effis_value = "review_int_no_effis",
    make_extra_no_effis = FALSE,

    compare_class3 = TRUE,
    class_col_internal = "class3_id",
    class_col_ref = "class3_id",
    class_levels_internal = NULL,
    class_levels_ref = NULL,

    hist_class_pairwise = TRUE,
    hist_class_pairwise_mode = c("cross_dataset", "all_pairs"),
    hist_class_pairwise_min_n = 20,
    hist_class_pairwise_max_pairs = Inf,
    hist_class_pairwise_include_within_dataset = FALSE,

    test_class_pairwise = c("all_pairs_internal_vs_ref",
                            "same_class_internal_vs_ref",
                            "within_internal",
                            "within_ref",
                            "all"),
    drop_na_class = TRUE,

    min_n_hist = 5,
    min_n_test = 10,
    do_ks = TRUE,

    save_test_outputs = FALSE,
    test_prefix = "phase4",
    return_sf = FALSE,
    verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # allow list outputs (scored_sf)
  if (is.list(internal_sf) && "scored_sf" %in% names(internal_sf)) internal_sf <- internal_sf$scored_sf
  if (is.list(ref_sf)      && "scored_sf" %in% names(ref_sf))      ref_sf      <- ref_sf$scored_sf
  stopifnot(inherits(internal_sf, "sf"), inherits(ref_sf, "sf"))

  # auto-detect internal_flag_col
  if (is.null(internal_flag_col)) {
    if ("flag_internal" %in% names(internal_sf)) internal_flag_col <- "flag_internal"
    else if ("flag_rbr" %in% names(internal_sf)) internal_flag_col <- "flag_rbr"
  }
  if (is.null(internal_flag_col) && (isTRUE(do_histograms) || isTRUE(do_tests))) {
    stop("internal_flag_col could not be auto-detected (expected 'flag_internal' or 'flag_rbr'). Please pass internal_flag_col.")
  }

  # compute median_rbr once if needed (missing OR empty)
  need_median <- "median_rbr" %in% value_cols && (
    (!"median_rbr" %in% names(internal_sf) || sum(is.finite(suppressWarnings(as.numeric(internal_sf$median_rbr)))) == 0) ||
      (!"median_rbr" %in% names(ref_sf)    || sum(is.finite(suppressWarnings(as.numeric(ref_sf$median_rbr)))) == 0)
  )
  if (isTRUE(need_median)) {
    if (is.null(rbr_rast)) stop("median_rbr is needed but rbr_rast is NULL.")
    stopifnot(inherits(rbr_rast, "SpatRaster"))
    msg("Computing median_rbr once for internal/ref (reused across value_cols).")
    internal_sf <- add_polygon_median_rbr(internal_sf, rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
    ref_sf      <- add_polygon_median_rbr(ref_sf,      rbr_rast, chunk_size = chunk_size, colname = "median_rbr")
  }

  out <- list(
    hist_paths = list(),
    tests = list(),
    out_dirs = list(hist_dir = hist_dir, tests_dir = tests_dir)
  )

  # histograms
  if (isTRUE(do_histograms)) {
    dir.create(hist_dir, recursive = TRUE, showWarnings = FALSE)
    hist_class_pairwise_mode <- match.arg(hist_class_pairwise_mode)

    out$hist_paths <- save_phase4_histograms_multi_vars(
      internal_sf = internal_sf,
      ref_sf      = ref_sf,
      rbr_rast    = rbr_rast,
      out_dir     = hist_dir,
      value_cols  = value_cols,

      buffer_m    = buffer_m,
      chunk_size  = chunk_size,

      internal_keep_col     = internal_flag_col,
      internal_keep_value   = internal_keep_value,
      internal_review_col   = internal_flag_col,
      internal_review_value = internal_review_value,

      ref_keep_col     = ref_flag_col,
      ref_keep_value   = ref_keep_value,
      ref_review_col   = ref_flag_col,
      ref_review_value = ref_review_value,

      use_internal_effis_flag = use_internal_effis_flag,
      internal_effis_col = internal_effis_col,
      keep_common_value = keep_common_value,
      review_int_no_effis_value = review_int_no_effis_value,
      make_extra_no_effis = make_extra_no_effis,

      compare_class3 = compare_class3,

      class_col_internal = class_col_internal,
      class_col_ref      = class_col_ref,
      class_levels_internal = class_levels_internal,
      class_levels_ref      = class_levels_ref,

      class_pairwise = hist_class_pairwise,
      class_pairwise_mode = hist_class_pairwise_mode,
      class_pairwise_min_n = hist_class_pairwise_min_n,
      class_pairwise_max_pairs = hist_class_pairwise_max_pairs,
      class_pairwise_include_within_dataset = hist_class_pairwise_include_within_dataset,

      min_n = min_n_hist,
      verbose = verbose
    )
  }

  # tests
  if (isTRUE(do_tests)) {
    dir.create(tests_dir, recursive = TRUE, showWarnings = FALSE)
    test_class_pairwise <- match.arg(test_class_pairwise)

    for (vcol in value_cols) {
      msg("---- Tests for value_col = '%s' -> %s", vcol, tests_dir)

      out$tests[[vcol]] <- test_phase4_distributions2(
        internal_sf = internal_sf,
        ref_sf      = ref_sf,
        rbr_rast    = rbr_rast,

        buffer_m    = buffer_m,
        chunk_size  = chunk_size,

        value_col   = vcol,
        ensure_median_rbr = FALSE,  # wrapper already computed if needed

        use_internal_effis_flag = use_internal_effis_flag,
        internal_effis_col = internal_effis_col,
        keep_common_value = keep_common_value,
        review_int_no_effis_value = review_int_no_effis_value,

        internal_flag_col = internal_flag_col,
        internal_keep_value = internal_keep_value,
        internal_review_value = internal_review_value,

        ref_flag_col = ref_flag_col,
        ref_keep_value = ref_keep_value,
        ref_review_value = ref_review_value,

        compare_class3 = compare_class3,
        class_col_internal = class_col_internal,
        class_col_ref      = class_col_ref,
        class_levels_internal = class_levels_internal,
        class_levels_ref      = class_levels_ref,
        class_pairwise = test_class_pairwise,
        drop_na_class = drop_na_class,

        min_n = min_n_test,
        do_ks = do_ks,

        save_outputs = save_test_outputs,
        out_dir = tests_dir,
        prefix = test_prefix,
        verbose = verbose
      )
    }
  }

  if (isTRUE(return_sf)) {
    out$internal_sf <- internal_sf
    out$ref_sf <- ref_sf
  }

  out
}
