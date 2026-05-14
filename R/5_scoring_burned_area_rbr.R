#' Phase 2 scoring: RBR-based keep-like flagging and optional class-3 labeling
#'
#' Computes per-polygon RBR summary statistics from an RBR raster and scores each polygon
#' relative to a "keep" reference pool (typically produced in Phase 1). The function:
#' (i) builds or reuses a keep-pool distribution, (ii) extracts per-polygon pixel values
#' from \code{rbr_rast}, (iii) computes robust features (median, pixel count, exceedance
#' probabilities), (iv) assigns a binary \code{keep_like} decision, and (optionally)
#' (v) assigns a 3-level quality class (\code{class3_id}) based on percentile thresholds
#' within the keep distribution.
#'
#' The keep pool is selected from \code{polys_sf} using either a common pool
#' (\code{keep_common_col == keep_common_value}) or a fallback pool
#' (\code{keep_fallback_col == keep_fallback_value}). You can also supply \code{keep_pool}
#' to reuse a previously computed pool and avoid recomputation.
#'
#' Outputs can be written to a GeoPackage (\code{driver = "GPKG"}) or Shapefile
#' (\code{driver = "ESRI Shapefile"}). When \code{arcgis_fix = TRUE}, the object written
#' is normalized to improve ArcGIS compatibility (valid geometries, MULTIPOLYGON, Z/M dropped,
#' optional reprojection to \code{output_epsg}).
#'
#' @param polys_sf \code{sf} polygon layer to score (Phase-2 target polygons).
#' @param rbr_rast \code{terra::SpatRaster}. Single-band RBR raster (first layer is used).
#' @param out_dir Character. Output folder (created if missing).
#' @param prefix Character. Prefix used for output filenames (default: \code{"phase2"}).
#'
#' @param keep_pool Optional list. If provided, the keep pool is reused and not recomputed.
#' Must contain: \code{min_area_ha}, \code{min_pix}, \code{promote_percentile},
#' \code{promote_p_above_ref}, \code{qref_prob}, \code{qref_keep}, \code{q25_keep},
#' \code{pixel_area_ha}, \code{keep_medians}.
#' @param use_keep_common_pool Logical. If TRUE, try to build the keep pool from
#' \code{keep_common_col == keep_common_value} first.
#' @param keep_common_col Character. Column used to define the common keep pool.
#' @param keep_common_value Character. Value indicating membership in the common keep pool.
#' @param keep_fallback_col Character. Fallback column used if the common pool is empty.
#' @param keep_fallback_value Character. Fallback value indicating keep membership.
#'
#' @param max_keep_samples Numeric. Max number of keep-pool pixels to retain for quantiles.
#' @param above_keep_prob Numeric in (0,1). Quantile probability used to define \code{qref_keep}.
#'
#' @param min_area_ha Numeric. Minimum polygon area (ha) required for confident scoring.
#' @param min_pix Integer or NULL. Minimum number of pixels required; if NULL it is derived
#' from \code{min_area_ha} and raster resolution (minimum 13 pixels).
#'
#' @param promote_percentile Numeric in (0,1). Minimum percentile (within keep distribution)
#' required to classify a polygon as \code{keep_like}.
#' @param promote_p_above_ref Numeric in (0,1). Minimum proportion of pixels above \code{qref_keep}
#' required to classify a polygon as \code{keep_like}.
#'
#' @param do_class3 Logical. If TRUE, create a 3-level class column (\code{out_col}).
#' @param out_col Character. Output column name for class-3 labeling (default: \code{"class3_id"}).
#' @param good_percentile Numeric in (0,1). Percentile threshold above which polygons are labeled
#' \code{"good_identified"}.
#' @param bad_percentile Numeric in (0,1). Percentile threshold below which polygons are labeled
#' \code{"bad_identified"}.
#' @param require_conf_ok Logical. If TRUE, polygons failing area/pixel confidence are forced to
#' \code{"ambiguous"}.
#'
#' @param chunk_size Integer. Chunk size for iterating polygons during raster extraction.
#' @param extract_exact Logical. Kept for API compatibility (not used; extraction uses
#' \code{terra::extract()} long-form output).
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @param save_outputs Logical. If TRUE, writes outputs to disk.
#' @param driver Character. \code{"GPKG"} or \code{"ESRI Shapefile"}.
#' @param overwrite Logical. If TRUE, overwrite existing outputs.
#' @param quiet Logical. Passed to \code{sf::st_write()}.
#' @param save_splits Logical. If TRUE, also saves \code{keep_only} and \code{review_only} layers.
#'
#' @param arcgis_fix Logical. If TRUE, normalizes geometries and CRS for ArcGIS-friendly outputs.
#' @param output_epsg Integer. EPSG code used when \code{arcgis_fix = TRUE} (default 3035).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{scored_sf}: input polygons with scoring fields appended.
#'   \item \code{keep_pool}: list describing the keep pool (can be reused).
#'   \item \code{out_paths}: list with output paths (if \code{save_outputs = TRUE}).
#' }
#'
#' @details
#' Per-polygon fields added include:
#' \describe{
#'   \item{\code{n_pix}}{Number of extracted pixels.}
#'   \item{\code{median_rbr}}{Median RBR within the polygon.}
#'   \item{\code{p_above_keep_q25}}{Proportion of pixels above keep-pool Q25.}
#'   \item{\code{p_above_keep_ref}}{Proportion of pixels above \code{qref_keep}.}
#'   \item{\code{percentile_in_keep}}{ECDF percentile of \code{median_rbr} within keep medians.}
#'   \item{\code{keep_like}}{0/1 decision.}
#'   \item{\code{flag_rbr}}{\code{"keep"} or \code{"review"} based on \code{keep_like}.}
#'   \item{\code{class3_id}}{Optional: \code{"good_identified"}, \code{"bad_identified"}, \code{"ambiguous"}.}
#' }
#'
#' @seealso \code{\link{scoring_burned_area_stage2}}
#'
#' @examples
#' \dontrun{
#' # After Phase 1:
#' p1 <- scoring_burned_area_stage2(..., save_outputs = FALSE)
#'
#' out_growth <- file.path(tempdir(), "phase2_demo")
#' dir.create(out_growth, recursive = TRUE, showWarnings = FALSE)
#'
#' p2 <- score_rbr_keep_like_class3(
#'   polys_sf  = p1$internal$stage2_keep_review_clean,
#'   rbr_rast  = rbr_rast,
#'   out_dir   = out_growth,
#'   prefix    = "BA_2012_phase2",
#'   save_outputs = TRUE,
#'   driver = "GPKG"
#' )
#' }
#'
#' @export
score_rbr_keep_classes <- function(
    polys_sf,
    rbr_rast,
    out_dir,
    prefix = "phase2",

    keep_pool = NULL,
    use_keep_common_pool = TRUE,
    keep_common_col   = "flag_internal_effis",
    keep_common_value = "keep_common",
    keep_fallback_col = "flag_internal",
    keep_fallback_value = "keep",

    max_keep_samples = 5e6,
    above_keep_prob  = 0.05,

    # QUALITY thresholds
    min_area_ha = 10,
    min_pix = NULL,

    # DECISION thresholds (keep_like)
    promote_percentile  = 0.10,
    promote_p_above_ref = 0.70,

    # CLASS3 thresholds
    do_class3 = TRUE,
    out_col = "class3_id",
    good_percentile = 0.05,
    bad_percentile  = 0.01,
    require_conf_ok = TRUE,

    # PERFORMANCE
    chunk_size = 500,
    extract_exact = FALSE,   # kept for API compatibility; we use terra::extract() long form
    verbose = TRUE,

    # SAVING (inside)
    save_outputs = TRUE,
    driver = "GPKG",
    overwrite = TRUE,
    quiet = TRUE,
    save_splits = FALSE,

    # ArcGIS compatibility
    arcgis_fix = TRUE,
    output_epsg = 3035,      # ETRS89 / LAEA Europe

    # --- NEW: CORINE NA handling for export/viewing
    normalize_corine_field = TRUE,
    corine_na_col = "corine_na_frac",
    corine_alt_cols = c("NA_CORINE_fraction", "na_corine_fraction", "corine_na_fraction"),
    fill_na_corine = FALSE,
    fill_na_corine_value = 0
) {
  stopifnot(inherits(polys_sf, "sf"))
  stopifnot(inherits(rbr_rast, "SpatRaster"))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Could not create out_dir: ", out_dir)

  if (!driver %in% c("GPKG", "ESRI Shapefile")) {
    stop("driver must be 'GPKG' or 'ESRI Shapefile'.")
  }

  msg <- function(...) if (isTRUE(verbose)) message(...)

  # -------------------------
  # Helpers
  # -------------------------
  safe_make_valid_sf <- function(x) {
    x <- sf::st_make_valid(x)
    empty <- sf::st_is_empty(x)
    if (any(empty)) x <- x[!empty, , drop = FALSE]

    gt <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
    if (any(gt %in% "GEOMETRYCOLLECTION")) {
      x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
      empty <- sf::st_is_empty(x)
      if (any(empty)) x <- x[!empty, , drop = FALSE]
    }
    x
  }

  delete_shapefile_bundle <- function(path_shp) {
    base <- tools::file_path_sans_ext(path_shp)
    exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx", ".xml")
    for (e in exts) {
      f <- paste0(base, e)
      if (file.exists(f)) try(file.remove(f), silent = TRUE)
    }
  }

  make_gpkg_safe <- function(x, sep = "|") {
    stopifnot(inherits(x, "sf"))

    # list-cols -> character
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

    # gpkg-safe names
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
    x
  }

  pick_extract_cols <- function(ex_df) {
    nms <- names(ex_df)
    id_col <- if ("ID" %in% nms) "ID" else nms[1]

    cov_candidates <- c("coverage_fraction", "fraction", "weight", "weights")
    cov_col <- cov_candidates[cov_candidates %in% nms][1]
    if (is.na(cov_col)) cov_col <- NULL

    value_cols <- setdiff(nms, c(id_col, cov_candidates))
    if (length(value_cols) < 1) stop("Could not identify value column from terra::extract output.")
    val_col <- value_cols[1]

    list(id_col = id_col, val_col = val_col, cov_col = cov_col)
  }

  # -------------------------
  # Validate and align CRS
  # -------------------------
  polys <- safe_make_valid_sf(polys_sf)
  if (nrow(polys) == 0) stop("polys_sf is empty after validation.")

  # Normalize CORINE NA field name (optional)
  if (isTRUE(normalize_corine_field)) {
    nms <- names(polys)
    if (!corine_na_col %in% nms) {
      alt_hit <- corine_alt_cols[corine_alt_cols %in% nms][1]
      if (!is.na(alt_hit)) {
        polys[[corine_na_col]] <- polys[[alt_hit]]
      }
    }
  }

  r <- rbr_rast[[1]]

  if (is.na(sf::st_crs(polys))) stop("polys_sf has NA CRS.")
  crs_r_wkt <- terra::crs(r, proj = TRUE)
  if (is.na(crs_r_wkt) || !nzchar(crs_r_wkt)) stop("rbr_rast has NA/empty CRS.")

  crs_poly_wkt <- sf::st_crs(polys)$wkt
  crs_r_sf <- sf::st_crs(crs_r_wkt)
  if (!identical(crs_poly_wkt, crs_r_sf$wkt)) {
    polys <- sf::st_transform(polys, crs_r_sf)
  }

  res_xy <- terra::res(r)
  pixel_area_m2 <- as.numeric(res_xy[1] * res_xy[2])
  pixel_area_ha <- pixel_area_m2 / 10000

  if (is.null(min_pix) || !is.finite(min_pix)) {
    min_pix <- max(ceiling(min_area_ha / pixel_area_ha), 13)
  }

  polys$._row_id <- seq_len(nrow(polys))

  # -------------------------
  # Keep pool computation (if keep_pool == NULL)
  # -------------------------
  if (is.null(keep_pool)) {

    keep_idx <- integer(0)

    if (isTRUE(use_keep_common_pool) && keep_common_col %in% names(polys)) {
      keep_idx <- which(as.character(polys[[keep_common_col]]) == keep_common_value)
    }

    if (length(keep_idx) == 0) {
      if (!keep_fallback_col %in% names(polys)) {
        stop("No keep pool found. Missing keep_fallback_col in polys: ", keep_fallback_col)
      }
      keep_idx <- which(as.character(polys[[keep_fallback_col]]) == keep_fallback_value)
    }

    if (length(keep_idx) == 0) stop("Keep pool is empty (no polygons matched keep pool rules).")

    keep_polys <- polys[keep_idx, , drop = FALSE]
    msg("Keep pool polys: ", nrow(keep_polys))

    # (A) sample keep-pool pixels for qref_keep and q25_keep
    keep_vals <- numeric(0)

    for (k0 in seq(1, nrow(keep_polys), by = chunk_size)) {
      k1 <- min(nrow(keep_polys), k0 + chunk_size - 1)
      chunk <- keep_polys[k0:k1, , drop = FALSE]

      v <- terra::vect(chunk)
      if (!terra::same.crs(v, r)) v <- terra::project(v, terra::crs(r))

      ex <- terra::extract(r, v)
      ex <- as.data.frame(ex)
      if (nrow(ex) == 0) next

      cols <- pick_extract_cols(ex)
      vals <- suppressWarnings(as.numeric(ex[[cols$val_col]]))
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next

      keep_vals <- c(keep_vals, vals)
      if (length(keep_vals) > max_keep_samples) {
        keep_vals <- sample(keep_vals, max_keep_samples)
      }
    }

    if (length(keep_vals) == 0) stop("Keep pool pixel extraction returned 0 finite values.")

    qref_keep <- as.numeric(stats::quantile(keep_vals, probs = above_keep_prob, na.rm = TRUE, names = FALSE))
    q25_keep  <- as.numeric(stats::quantile(keep_vals, probs = 0.25,          na.rm = TRUE, names = FALSE))

    msg(sprintf("qref_keep (P%.0f) = %.3f; q25_keep = %.3f",
                above_keep_prob * 100, qref_keep, q25_keep))

    # (B) keep-pool medians for ECDF(percentile_in_keep)
    keep_medians <- rep(NA_real_, nrow(keep_polys))

    for (k0 in seq(1, nrow(keep_polys), by = chunk_size)) {
      k1 <- min(nrow(keep_polys), k0 + chunk_size - 1)
      chunk <- keep_polys[k0:k1, , drop = FALSE]

      v <- terra::vect(chunk)
      if (!terra::same.crs(v, r)) v <- terra::project(v, terra::crs(r))

      ex <- terra::extract(r, v)
      ex <- as.data.frame(ex)
      mm <- rep(NA_real_, nrow(chunk))

      if (nrow(ex) > 0) {
        cols <- pick_extract_cols(ex)
        id <- suppressWarnings(as.integer(ex[[cols$id_col]]))
        val <- suppressWarnings(as.numeric(ex[[cols$val_col]]))

        ok <- is.finite(id) & is.finite(val) & id >= 1 & id <= nrow(chunk)
        if (any(ok)) {
          meds <- tapply(val[ok], id[ok], stats::median)
          mm[as.integer(names(meds))] <- as.numeric(meds)
        }
      }

      keep_medians[k0:k1] <- mm
    }

    keep_medians <- keep_medians[is.finite(keep_medians)]
    if (length(keep_medians) == 0) stop("Keep pool medians are all NA/inf (keep polys may be outside raster).")

    keep_pool <- list(
      min_area_ha = min_area_ha,
      min_pix     = min_pix,
      promote_percentile  = promote_percentile,
      promote_p_above_ref = promote_p_above_ref,
      qref_prob   = above_keep_prob,
      qref_keep   = qref_keep,
      q25_keep    = q25_keep,
      pixel_area_ha = pixel_area_ha,
      keep_medians = keep_medians
    )

  } else {
    req <- c("min_area_ha","min_pix","promote_percentile","promote_p_above_ref",
             "qref_prob","qref_keep","q25_keep","pixel_area_ha","keep_medians")
    miss <- setdiff(req, names(keep_pool))
    if (length(miss) > 0) stop("keep_pool is missing fields: ", paste(miss, collapse = ", "))

    min_area_ha <- keep_pool$min_area_ha
    min_pix     <- keep_pool$min_pix
    promote_percentile  <- keep_pool$promote_percentile
    promote_p_above_ref <- keep_pool$promote_p_above_ref
    above_keep_prob <- keep_pool$qref_prob

    msg("Using provided keep_pool (no recomputation).")
    qref_keep <- keep_pool$qref_keep
    q25_keep  <- keep_pool$q25_keep
  }

  keep_ecdf <- stats::ecdf(keep_pool$keep_medians)
  qref_keep <- keep_pool$qref_keep
  q25_keep  <- keep_pool$q25_keep

  # -------------------------
  # Extract per-polygon stats
  # -------------------------
  n <- nrow(polys)
  out_npix <- integer(n)
  out_med  <- rep(NA_real_, n)
  out_pq25 <- rep(NA_real_, n)
  out_pref <- rep(NA_real_, n)

  for (i0 in seq(1, n, by = chunk_size)) {
    i1 <- min(n, i0 + chunk_size - 1)
    chunk <- polys[i0:i1, , drop = FALSE]

    v <- terra::vect(chunk)
    if (!terra::same.crs(v, r)) v <- terra::project(v, terra::crs(r))

    ex <- terra::extract(r, v)
    ex <- as.data.frame(ex)

    npix <- rep(0L, nrow(chunk))
    med  <- rep(NA_real_, nrow(chunk))
    pq25 <- rep(NA_real_, nrow(chunk))
    pref <- rep(NA_real_, nrow(chunk))

    if (nrow(ex) > 0) {
      cols <- pick_extract_cols(ex)
      id <- suppressWarnings(as.integer(ex[[cols$id_col]]))
      val <- suppressWarnings(as.numeric(ex[[cols$val_col]]))

      ok <- is.finite(id) & is.finite(val) & id >= 1 & id <= nrow(chunk)
      if (any(ok)) {
        id_ok <- id[ok]
        val_ok <- val[ok]

        n_by <- tapply(val_ok, id_ok, length)
        m_by <- tapply(val_ok, id_ok, stats::median)
        q_by <- tapply(val_ok > q25_keep, id_ok, mean)
        r_by <- tapply(val_ok > qref_keep, id_ok, mean)

        idx <- as.integer(names(n_by))
        npix[idx] <- as.integer(n_by)
        med[idx]  <- as.numeric(m_by)
        pq25[idx] <- as.numeric(q_by)
        pref[idx] <- as.numeric(r_by)
      }
    }

    out_npix[i0:i1] <- npix
    out_med[i0:i1]  <- med
    out_pq25[i0:i1] <- pq25
    out_pref[i0:i1] <- pref
  }

  # -------------------------
  # Build scored sf
  # -------------------------
  scored <- polys
  scored$n_pix <- out_npix
  scored$median_rbr <- out_med
  scored$p_above_keep_q25 <- out_pq25
  scored$p_above_keep_ref <- out_pref

  scored$area_ha <- as.numeric(sf::st_area(scored)) / 10000

  scored$conf_area <- ifelse(scored$area_ha >= min_area_ha, "ok", "too_small")
  scored$conf_pix  <- ifelse(scored$n_pix   >= min_pix,     "ok", "too_few_pix")

  scored$percentile_in_keep <- ifelse(
    is.finite(scored$median_rbr),
    keep_ecdf(scored$median_rbr),
    NA_real_
  )

  scored$keep_like <- with(scored,
                           as.integer(
                             conf_area == "ok" &
                               conf_pix  == "ok" &
                               is.finite(percentile_in_keep) & (percentile_in_keep >= promote_percentile) &
                               is.finite(p_above_keep_ref)   & (p_above_keep_ref   >= promote_p_above_ref)
                           )
  )

  scored$flag_rbr <- ifelse(scored$keep_like == 1, "keep", "review")

  # -------------------------
  # CLASS3 classification
  # -------------------------
  if (isTRUE(do_class3)) {
    cls <- rep("ambiguous", nrow(scored))

    ok_conf <- (scored$conf_area == "ok" & scored$conf_pix == "ok")
    pct <- scored$percentile_in_keep

    good <- is.finite(pct) & (pct >= good_percentile)
    bad  <- is.finite(pct) & (pct <= bad_percentile)

    cls[good] <- "good_identified"
    cls[bad]  <- "bad_identified"

    if (isTRUE(require_conf_ok)) cls[!ok_conf] <- "ambiguous"

    scored[[out_col]] <- cls

    attr(scored, "idclass_percentiles") <- pct
    attr(scored, "idclass_pool_quartiles_median_rbr") <- stats::quantile(
      keep_pool$keep_medians,
      probs = c(bad_percentile, good_percentile, 0.25, 0.50, 0.75),
      na.rm = TRUE
    )
  }

  scored$._row_id <- NULL

  # -------------------------
  # Export-only: make CORINE field numeric + optionally fill NA
  # -------------------------
  if (corine_na_col %in% names(scored)) {
    scored[[corine_na_col]] <- suppressWarnings(as.numeric(scored[[corine_na_col]]))
    if (isTRUE(fill_na_corine)) {
      scored[[corine_na_col]][is.na(scored[[corine_na_col]])] <- fill_na_corine_value
    }
  }

  # -------------------------
  # Saving outputs
  # -------------------------
  out_paths <- list(
    container = NA_character_,
    scored = NA_character_,
    keep_only = NA_character_,
    review_only = NA_character_,
    keep_pool_txt = NA_character_
  )

  if (isTRUE(save_outputs)) {

    scored_out <- scored

    if (isTRUE(arcgis_fix)) {
      scored_out <- sf::st_make_valid(scored_out)
      scored_out <- sf::st_zm(scored_out, drop = TRUE, what = "ZM")
      scored_out <- sf::st_cast(scored_out, "MULTIPOLYGON", warn = FALSE)

      if (!is.null(output_epsg) && is.finite(output_epsg)) {
        scored_out <- tryCatch(
          sf::st_transform(scored_out, output_epsg),
          error = function(e) {
            sf::st_crs(scored_out) <- output_epsg
            scored_out
          }
        )
      }
    }

    scored_w <- make_gpkg_safe(scored_out)

    if (driver == "GPKG") {
      gpkg_path <- file.path(out_dir, sprintf("%s_scored.gpkg", prefix))
      if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

      sf::st_write(scored_w, gpkg_path, layer = "scored", driver = "GPKG",
                   delete_dsn = TRUE, quiet = quiet)

      out_paths$container <- gpkg_path
      out_paths$scored <- paste0(gpkg_path, " | layer=scored")

      if (isTRUE(save_splits)) {
        keep_only   <- scored_w[scored_w$keep_like == 1, , drop = FALSE]
        review_only <- scored_w[scored_w$keep_like == 0, , drop = FALSE]
        sf::st_write(keep_only,   gpkg_path, layer = "keep_only",   append = TRUE, quiet = quiet)
        sf::st_write(review_only, gpkg_path, layer = "review_only", append = TRUE, quiet = quiet)
        out_paths$keep_only   <- paste0(gpkg_path, " | layer=keep_only")
        out_paths$review_only <- paste0(gpkg_path, " | layer=review_only")
      }

    } else {
      p_scored <- file.path(out_dir, sprintf("%s_scored.shp", prefix))
      if (isTRUE(overwrite)) delete_shapefile_bundle(p_scored)
      sf::st_write(scored_w, p_scored, driver = "ESRI Shapefile", quiet = quiet)
      out_paths$container <- out_dir
      out_paths$scored <- p_scored
    }

    txt_path <- file.path(out_dir, sprintf("%s_keep_pool.txt", prefix))
    lines <- c(
      "KEEP POOL SUMMARY",
      sprintf("min_area_ha: %s", keep_pool$min_area_ha),
      sprintf("min_pix: %s", keep_pool$min_pix),
      sprintf("promote_percentile: %s", keep_pool$promote_percentile),
      sprintf("promote_p_above_ref: %s", keep_pool$promote_p_above_ref),
      sprintf("qref_prob: %s", keep_pool$qref_prob),
      sprintf("qref_keep: %s", keep_pool$qref_keep),
      sprintf("q25_keep: %s", keep_pool$q25_keep),
      sprintf("pixel_area_ha: %s", keep_pool$pixel_area_ha),
      sprintf("keep_medians_n: %s", length(keep_pool$keep_medians))
    )
    writeLines(lines, txt_path)
    out_paths$keep_pool_txt <- txt_path
  }

  return(list(
    scored_sf = scored,
    keep_pool = keep_pool,
    out_paths = out_paths
  ))
}
