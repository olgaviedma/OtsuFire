#' Internal patch-level feature extraction
#'
#' @description
#' Internal supervised feature-extraction engine used by the package
#' orchestrator.
#'
#' @details
#'   Ecoregions belong exclusively to the Otsu-guided segmentation stage
#'   and are not part of the supervised feature set (2026-06-05).
#' @keywords internal
#' @noRd
extract_features <- function(
    train_folds,
    unlabeled,
    
    # --- ID + class ---
    id_col    = "fire_uid",
    class_col = "class",
    pos_lab   = "burned",
    neg_lab   = "unburned",
    
    # --- build features from layers? ---
    build_features = TRUE,
    
    # Rasters (terra SpatRaster)
    rbr_summer = NULL,     # REQUIRED if build_features=TRUE
    doy_post   = NULL,     # optional
    rbr_aw     = NULL,     # optional
    nbr_pre    = NULL,     # optional
    nbr_post   = NULL,     # optional
    dnbr       = NULL,     # optional
    dem        = NULL,     # REQUIRED if build_features=TRUE
    slope      = NULL,     # REQUIRED if build_features=TRUE
    corine_r   = NULL,     # REQUIRED if build_features=TRUE (categorical)

    # Hotspots (sf points)
    hotspots   = NULL,     # optional sf POINTS
    frp_col    = "FRP",
    conf_col   = "confidence",
    year_col   = "year",
    month_col  = "month",
    hs_start_month = 6,
    hs_end_month   = 10,
    hi_conf_thr    = 0.8,
    
    # CORINE groups: named list -> columns cor_<name>_frac
    cor_groups = NULL,     # REQUIRED if build_features=TRUE
    
    # Feature switches
    use_doy = TRUE,
    use_aw  = TRUE,
    use_nbr = FALSE,
    use_hotspots = TRUE,
    # OPTIONAL shape/size block (OtsuFire 0.12.0, OFF by default). When
    # TRUE the per-polygon shape_features() helper computes log_area,
    # perim_m, compactness and elongation and they are left-joined onto
    # BOTH the train and scoring feature tables. With use_shape = FALSE
    # no join runs => extraction output is byte-identical to today.
    use_shape = FALSE,

    # --- Hotspots by year logic ---
    year_target = NULL,            # e.g., 2022
    hs_year_min_available = 2000,
    hs_buffer_m = 2000,
    hs_missing_value = -9999,
    hs_use_season_filter = FALSE,
    
    # exactextractr memory
    max_cells_in_memory = NULL,
    
    # --- CRS + grid handling ---
    crs_target = 3035,                 # default EPSG:3035
    reproject_rasters = TRUE,
    align_rasters_to_template = TRUE,  # if TRUE, resample to template grid when mismatch
    template_raster = NULL,            # if NULL uses rbr_summer (after CRS ops)
    grid_check_strict = FALSE,         # if TRUE and align=FALSE -> stop on mismatch
    grid_tol = 1e-6,                   # tolerance for res/origin comparisons
    
    # --- raster cleaning (extreme nodata) ---
    rbr_aw_min_ok = -10000,            # values < this set to NA
    
    # --- return/save ---
    return_features = TRUE,
    return_features_geometry = TRUE,
    save_features_dir = NULL,
    save_features_format = c("rds","csv"),
    save_features_gpkg = TRUE,
    train_features_basename = "train_features",
    scoring_features_basename = "scoring_features",
    train_features_layer = "train_features",
    scoring_features_layer = "scoring_features",

    # 0.4.0 architectural refactor (Agent H): `extract_features()` is
    # ADD-only. It receives `train_folds` and `unlabeled` and returns
    # the same sf objects with the 50 generated features ADDED. It
    # does NOT drop any input column. Filtering is done at the LAST
    # mile, inside `train_final_model_direct()`, against the
    # canonical whitelist `.supervised_feature_cols`.
    #
    # The `strict_input` and `allow_passthrough_cols` arguments are
    # retained for backward-compatibility with 0.3.x callers but are
    # no-ops since 0.4.0. They emit a one-line message when
    # `strict_input = TRUE` and otherwise have no effect on the data.
    strict_input = TRUE,
    allow_passthrough_cols = character(0),

    verbose = TRUE
) {

  # Block 7: library() calls removed; deps in DESCRIPTION Imports.
  save_features_format <- match.arg(save_features_format)

  # ---------------------------------------------------------------
  # 0.4.0 architectural refactor: input-contract enforcement moved.
  # `extract_features()` no longer drops anything. The whitelist
  # filter at the model-matrix step (inside
  # `train_final_model_direct()`) is the single source of truth for
  # which columns reach the supervised model.
  # ---------------------------------------------------------------
  if (isTRUE(strict_input)) {
    message("extract_features(): input-contract enforcement moved to ",
            "train_final_model_direct() since OtsuFire 0.4.0; ",
            "all input columns are preserved here.")
  }
  if (length(allow_passthrough_cols) > 0L) {
    message("extract_features(): allow_passthrough_cols is a no-op ",
            "since OtsuFire 0.4.0; all input columns are preserved here.")
  }
  # End 0.4.0 input contract no-op.
  if (!is.null(max_cells_in_memory)) {
    options(exactextractr.max_cells_in_memory = max_cells_in_memory)
  }
  
  # ---------------------------
  # Helpers
  # ---------------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  assert_cols <- function(df, cols) {
    miss <- setdiff(cols, names(df))
    if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  
  drop_geom <- function(x) if (inherits(x, "sf")) sf::st_drop_geometry(x) else x
  
  safe_make_valid <- function(x) {
    if (!inherits(x, "sf")) return(x)
    ok <- suppressWarnings(sf::st_is_valid(x))
    if (all(ok, na.rm = TRUE)) return(x)
    sf::st_make_valid(x)
  }
  
  drop_empty_sf <- function(x) {
    if (!inherits(x, "sf")) return(x)
    g <- sf::st_geometry(x)
    ok <- !sf::st_is_empty(g)
    ok[is.na(ok)] <- FALSE
    x[ok, , drop = FALSE]
  }
  
  as_num_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))
  as_int_safe <- function(x) suppressWarnings(as.integer(as.character(x)))
  
  # --- FIX: weighted quantiles (needed by zonal_numeric_stats) ---
  wtd_quantile <- function(x, w, probs = c(0.1, 0.5, 0.9)) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    x <- x[ok]; w <- w[ok]
    if (!length(x)) return(setNames(rep(NA_real_, length(probs)), as.character(probs)))
    o <- order(x)
    x <- x[o]; w <- w[o]
    cw <- cumsum(w) / sum(w)
    q <- sapply(probs, function(p) {
      if (p <= cw[1]) return(x[1])
      if (p >= cw[length(cw)]) return(x[length(x)])
      j <- which(cw >= p)[1]
      j0 <- j - 1
      x0 <- x[j0]; x1 <- x[j]
      c0 <- cw[j0]; c1 <- cw[j]
      if (!is.finite(c0) || !is.finite(c1) || c1 == c0) return(x1)
      x0 + (p - c0) * (x1 - x0) / (c1 - c0)
    })
    names(q) <- as.character(probs)
    q
  }
  
  # --- clean extreme nodata inside rasters ---
  clean_raster_extremes <- function(r, min_ok = -10000, max_ok = NULL, tag = "raster") {
    if (is.null(r)) return(r)
    if (!inherits(r, "SpatRaster")) return(r)
    
    if (!is.null(max_ok)) {
      r2 <- terra::ifel(r < min_ok | r > max_ok, NA, r)
    } else {
      r2 <- terra::ifel(r < min_ok, NA, r)
    }
    
    if (isTRUE(verbose)) {
      n_bad <- tryCatch(terra::global(r < min_ok, "sum", na.rm = TRUE)[1,1], error = function(e) NA)
      msg("clean_raster_extremes(%s): set NA for values < %s (bad pixels=%s)", tag, min_ok, n_bad)
    }
    r2
  }
  
  # --- CRS resolver ---
  resolve_crs_sf <- function(x) {
    if (is.null(x)) x <- 3035
    crs <- sf::st_crs(x)
    if (is.na(crs)) stop("crs_target could not be resolved. Provide EPSG (e.g., 3035) or valid CRS string/WKT.")
    crs
  }
  crs_sf <- resolve_crs_sf(crs_target)
  
  crs_to_terra <- function(crs_sf, raw) {
    if (is.numeric(raw) && length(raw) == 1 && is.finite(raw)) return(paste0("EPSG:", raw))
    if (is.character(raw) && length(raw) == 1) return(raw)
    if (!is.null(crs_sf$wkt) && !is.na(crs_sf$wkt)) return(crs_sf$wkt)
    if (!is.null(crs_sf$epsg) && !is.na(crs_sf$epsg)) return(paste0("EPSG:", crs_sf$epsg))
    stop("Could not build a CRS string for terra projection.")
  }
  crs_terra <- crs_to_terra(crs_sf, crs_target)
  
  to_crs_safe_sf <- function(x, crs_target_sf, name = "sf") {
    if (!inherits(x, "sf")) return(x)
    if (nrow(x) == 0) return(x)
    if (is.na(sf::st_crs(x))) stop(sprintf("Object '%s' has NA CRS; cannot transform safely.", name))
    if (sf::st_crs(x) != crs_target_sf) {
      if (isTRUE(verbose)) msg("Transforming %s to target CRS (EPSG:%s)...", name, crs_target_sf$epsg %||% "custom")
      x <- sf::st_transform(x, crs_target_sf)
    }
    x
  }
  
  # --- grid check (robust, tolerant) ---
  grid_equal <- function(r, template, tol = 1e-6) {
    if (is.null(r) || is.null(template)) return(TRUE)
    if (!inherits(r, "SpatRaster") || !inherits(template, "SpatRaster")) return(FALSE)
    
    rr <- terra::res(r); rt <- terra::res(template)
    or <- terra::origin(r); ot <- terra::origin(template)
    
    if (any(!is.finite(rr)) || any(!is.finite(rt)) || any(!is.finite(or)) || any(!is.finite(ot))) return(FALSE)
    
    same_res <- all(abs(rr - rt) <= tol)
    same_ori <- all(abs(or - ot) <= tol)
    
    same_res && same_ori
  }
  
  ensure_raster_crs_and_grid <- function(r, template, name,
                                         method_project = "bilinear",
                                         method_resample = "bilinear") {
    if (is.null(r)) return(r)
    if (!inherits(r, "SpatRaster")) stop(sprintf("Raster '%s' must be a terra SpatRaster.", name))
    
    # CRS present?
    r_crs <- terra::crs(r)
    if (is.na(r_crs) || r_crs == "") stop(sprintf("Raster '%s' has NA/empty CRS; cannot reproject safely.", name))
    
    # Target CRS mismatch?
    raster_in_target <- tryCatch(sf::st_crs(r_crs) == crs_sf, error = function(e) FALSE)
    if (!isTRUE(raster_in_target)) {
      if (!isTRUE(reproject_rasters)) {
        stop(sprintf("CRS mismatch: raster '%s' not in target CRS. Set reproject_rasters=TRUE.", name))
      }
      if (isTRUE(verbose)) msg("Reprojecting raster '%s' to target CRS (project method=%s)...", name, method_project)
      
      # Best: project directly to template grid if available
      if (!is.null(template) && inherits(template, "SpatRaster")) {
        r <- terra::project(r, template, method = method_project)
      } else {
        r <- terra::project(r, crs_terra, method = method_project)
      }
    }
    
    # Grid alignment vs template
    if (!is.null(template) && inherits(template, "SpatRaster")) {
      if (isTRUE(align_rasters_to_template)) {
        if (!grid_equal(r, template, tol = grid_tol)) {
          if (isTRUE(verbose)) {
            rr <- terra::res(r); rt <- terra::res(template)
            or <- terra::origin(r); ot <- terra::origin(template)
            msg("Grid mismatch detected for '%s' -> resampling to template (method=%s).", name, method_resample)
            msg("  '%s' res=(%s,%s) origin=(%s,%s)", name, rr[1], rr[2], or[1], or[2])
            msg("  template res=(%s,%s) origin=(%s,%s)", rt[1], rt[2], ot[1], ot[2])
          }
          r <- terra::resample(r, template, method = method_resample)
        }
      } else if (isTRUE(grid_check_strict)) {
        if (!grid_equal(r, template, tol = grid_tol)) {
          stop(sprintf(
            "Grid mismatch: raster '%s' does not match template res/origin.\nSet align_rasters_to_template=TRUE to auto-resample.",
            name
          ))
        }
      }
    }
    
    r
  }
  
  parse_confidence01 <- function(x) {
    xn <- as_num_safe(x)
    ok_num <- sum(!is.na(xn)) / max(1, length(xn)) > 0.5
    if (ok_num) {
      out <- ifelse(xn > 1, xn/100, xn)
      out <- pmin(pmax(out, 0), 1)
      out[is.na(out)] <- 0
      return(out)
    }
    xs <- tolower(trimws(as.character(x)))
    out <- rep(0, length(xs))
    out[xs %in% c("h","high")]     <- 0.9
    out[xs %in% c("n","nominal")]  <- 0.6
    out[xs %in% c("l","low")]      <- 0.3
    out[is.na(xs) | xs == ""] <- 0
    out
  }
  
  zonal_numeric_stats <- function(polys, r, prefix,
                                  probs = c(0.10,0.50,0.90), compute_sd = FALSE) {
    polys <- safe_make_valid(polys)
    polys2 <- drop_empty_sf(polys)
    
    all_ids <- unique(as.character(polys[[id_col]]))
    
    if (nrow(polys2) == 0) {
      out0 <- tibble(!!id_col := all_ids)
      out0[[paste0(prefix, "_valid_frac")]] <- NA_real_
      out0[[paste0(prefix, "_p10")]] <- NA_real_
      out0[[paste0(prefix, "_med")]] <- NA_real_
      out0[[paste0(prefix, "_p90")]] <- NA_real_
      out0[[paste0(prefix, "_iqr")]] <- NA_real_
      if (compute_sd) {
        out0[[paste0(prefix, "_mean")]] <- NA_real_
        out0[[paste0(prefix, "_sd")]] <- NA_real_
      }
      return(out0)
    }
    
    ex <- exactextractr::exact_extract(r, polys2, include_cols = id_col)
    
    out <- lapply(ex, function(df) {
      v <- df$value
      w <- df$coverage_fraction
      
      denom_all <- sum(df$coverage_fraction, na.rm = TRUE)
      ok <- is.finite(v) & is.finite(w) & w > 0
      v2 <- v[ok]; w2 <- w[ok]
      valid_frac <- if (denom_all > 0) sum(w2, na.rm = TRUE) / denom_all else 0
      
      if (length(v2) == 0) {
        res <- list(valid_frac=0, p10=NA_real_, p50=NA_real_, p90=NA_real_, iqr=NA_real_)
        if (compute_sd) res <- c(res, list(mean=NA_real_, sd=NA_real_))
        return(res)
      }
      
      qs  <- wtd_quantile(v2, w2, probs = probs)
      q25 <- wtd_quantile(v2, w2, probs = 0.25)[["0.25"]]
      q75 <- wtd_quantile(v2, w2, probs = 0.75)[["0.75"]]
      
      res <- list(
        valid_frac = valid_frac,
        p10 = qs[[as.character(probs[1])]],
        p50 = qs[[as.character(probs[2])]],
        p90 = qs[[as.character(probs[3])]],
        iqr = q75 - q25
      )
      
      if (compute_sd) {
        m <- weighted.mean(v2, w2, na.rm = TRUE)
        s <- sqrt(weighted.mean((v2 - m)^2, w2, na.rm = TRUE))
        res$mean <- m
        res$sd   <- s
      }
      res
    })
    
    df_out <- dplyr::bind_rows(out)
    df_out <- dplyr::bind_cols(tibble(!!id_col := polys2[[id_col]]), df_out)
    
    nm_map <- c(
      valid_frac = paste0(prefix, "_valid_frac"),
      p10 = paste0(prefix, "_p10"),
      p50 = paste0(prefix, "_med"),
      p90 = paste0(prefix, "_p90"),
      iqr = paste0(prefix, "_iqr"),
      mean = paste0(prefix, "_mean"),
      sd   = paste0(prefix, "_sd")
    )
    names(df_out) <- ifelse(names(df_out) %in% names(nm_map), nm_map[names(df_out)], names(df_out))
    
    tibble(!!id_col := all_ids) %>% dplyr::left_join(df_out, by = id_col)
  }
  
  corine_group_fractions <- function(polys, groups) {
    polys <- safe_make_valid(polys)
    polys2 <- drop_empty_sf(polys)
    all_ids <- unique(as.character(polys[[id_col]]))
    
    if (nrow(polys2) == 0) {
      base <- tibble(!!id_col := all_ids)
      for (nm in names(groups)) base[[paste0("cor_", nm, "_frac")]] <- NA_real_
      return(base)
    }
    
    ex <- exactextractr::exact_extract(corine_r, polys2, include_cols = id_col)
    
    frac_group <- function(df, codes) {
      v <- df$value; w <- df$coverage_fraction
      ok <- is.finite(v) & is.finite(w) & w > 0
      v <- v[ok]; w <- w[ok]
      denom <- sum(w, na.rm = TRUE)
      if (denom <= 0) return(0)
      sum(w[v %in% codes], na.rm = TRUE) / denom
    }
    
    out <- lapply(ex, function(df) lapply(groups, function(codes) frac_group(df, codes)))
    df_out <- dplyr::bind_rows(lapply(out, as.data.frame))
    names(df_out) <- paste0("cor_", names(groups), "_frac")
    df_out <- dplyr::bind_cols(tibble(!!id_col := polys2[[id_col]]), df_out)
    
    tibble(!!id_col := all_ids) %>% dplyr::left_join(df_out, by = id_col)
  }
  
  hotspot_features <- function(polys) {
    polys <- safe_make_valid(polys)
    polys2 <- drop_empty_sf(polys)
    all_ids <- unique(as.character(polys[[id_col]]))
    base_ids <- tibble(!!id_col := all_ids)
    
    yy <- year_target
    if (is.null(yy) && year_col %in% names(polys)) yy <- as_int_safe(polys[[year_col]][1])
    
    out_nodata <- function() {
      # B5 (2026-06-06): for years/polygons with NO hotspot data (use_hotspots
      # FALSE, year < hs_year_min_available, or hotspots NULL), the MEASURED
      # hotspot quantities are genuinely UNKNOWN, not zero and not -9999. We
      # emit NA_real_ so the design-matrix builder (internal-sup-create-matrix.R)
      # flags them via `<col>_isNA = 1` and imputes them to the column median.
      # This is the documented "impute + flag" path and lets a model trained on
      # hotspot-years follow the learned missing/default direction when applied
      # to a hotspot-less year (historical / LOYO application) instead of
      # ingesting -9999 as a real extreme number.
      #
      # Per-column decision (guiding principle: a MEASURED hotspot quantity with
      # no data -> NA so it is flagged + imputed; a STRUCTURAL availability flag
      # that is legitimately 0 when no data -> kept 0L):
      #   * hotspot_available = 0L          KEPT  : real flag "no hotspot data
      #                                             this year"; the model uses it.
      #   * hs_in_poly / hs_in_buffer / hs_used_n / hs_min_dist_m /
      #     hs_frp_sum / hs_frp_max / hs_conf_mean / hs_hiConf_n  -> NA_real_
      #                                     : measured counts/distances/FRP/conf,
      #                                       unknown without data.
      #   * hs_no_support_when_available = 0L  KEPT : its semantics are
      #     conditioned on availability ("...when available"); with
      #     hotspot_available == 0 it is correctly 0 (not "no support when
      #     available"), consistent with the per-polygon computation below.
      #   * hs_support_present / hs_only_buffer_support -> NA_real_ : these encode
      #     whether hotspot support exists / is buffer-only, which is undefined
      #     (unknown) when there is no hotspot data at all.
      base_ids %>%
        dplyr::transmute(
          !!id_col,
          hotspot_available = 0L,
          hs_in_poly   = NA_real_,
          hs_in_buffer = NA_real_,
          hs_used_n    = NA_real_,
          hs_min_dist_m = NA_real_,
          hs_frp_sum   = NA_real_,
          hs_frp_max   = NA_real_,
          hs_conf_mean = NA_real_,
          hs_hiConf_n  = NA_real_,
          hs_support_present = NA_real_,
          hs_no_support_when_available = 0L,
          hs_only_buffer_support = NA_real_
        )
    }
    out_available_empty <- function() {
      base_ids %>%
        dplyr::transmute(
          !!id_col,
          hotspot_available = 1L,
          hs_in_poly   = 0L,
          hs_in_buffer = 0L,
          hs_used_n    = 0L,
          hs_min_dist_m = hs_buffer_m,
          hs_frp_sum   = 0,
          hs_frp_max   = 0,
          hs_conf_mean = 0,
          hs_hiConf_n  = 0L,
          hs_support_present = 0L,
          hs_no_support_when_available = 1L,
          hs_only_buffer_support = 0L
        )
    }
    
    if (!isTRUE(use_hotspots)) return(out_nodata())
    if (is.null(yy) || !is.finite(yy) || yy < hs_year_min_available) return(out_nodata())
    if (is.null(hotspots)) return(out_nodata())
    if (nrow(polys2) == 0) return(out_available_empty())
    
    hs <- safe_make_valid(hotspots)
    if (sf::st_crs(polys2) != sf::st_crs(hs)) hs <- sf::st_transform(hs, sf::st_crs(polys2))
    
    if (year_col %in% names(hs)) {
      hs <- hs %>% dplyr::mutate(.year = as_int_safe(.data[[year_col]])) %>% dplyr::filter(.year == yy)
    }
    if (isTRUE(hs_use_season_filter) && month_col %in% names(hs)) {
      hs <- hs %>% dplyr::mutate(.month = as_int_safe(.data[[month_col]]))
      m <- hs$.month
      in_season <- if (hs_start_month <= hs_end_month) (m >= hs_start_month & m <= hs_end_month) else (m >= hs_start_month | m <= hs_end_month)
      hs <- hs[in_season, , drop = FALSE]
    }
    
    if (nrow(hs) == 0) return(out_available_empty())
    
    polys_id <- polys2 %>% dplyr::select(dplyr::all_of(id_col))
    
    inside_list <- sf::st_intersects(polys_id, hs)
    n_inside <- lengths(inside_list)
    
    buffer_list <- sf::st_is_within_distance(polys_id, hs, dist = hs_buffer_m)
    n_buffer <- lengths(buffer_list)
    
    used_idx <- vector("list", length = nrow(polys_id))
    used_where <- rep("none", nrow(polys_id))
    for (i in seq_len(nrow(polys_id))) {
      if (n_inside[i] > 0) {
        used_idx[[i]] <- inside_list[[i]]
        used_where[i] <- "inside"
      } else if (n_buffer[i] > 0) {
        used_idx[[i]] <- buffer_list[[i]]
        used_where[i] <- "buffer"
      } else {
        used_idx[[i]] <- integer(0)
        used_where[i] <- "none"
      }
    }
    
    hs_min <- numeric(nrow(polys_id))
    for (i in seq_len(nrow(polys_id))) {
      if (used_where[i] == "inside") {
        hs_min[i] <- 0
      } else if (used_where[i] == "buffer") {
        idx <- used_idx[[i]]
        d <- sf::st_distance(sf::st_geometry(polys_id[i, ]), sf::st_geometry(hs[idx, , drop = FALSE]))
        hs_min[i] <- as.numeric(min(d))
      } else {
        hs_min[i] <- hs_buffer_m
      }
    }
    
    hs$.frp <- if (frp_col %in% names(hs)) as_num_safe(hs[[frp_col]]) else NA_real_
    hs$.conf01 <- if (conf_col %in% names(hs)) parse_confidence01(hs[[conf_col]]) else NA_real_
    
    frp_sum <- numeric(nrow(polys_id))
    frp_max <- numeric(nrow(polys_id))
    conf_mean <- numeric(nrow(polys_id))
    hi_n <- integer(nrow(polys_id))
    
    for (i in seq_len(nrow(polys_id))) {
      idx <- used_idx[[i]]
      if (length(idx) == 0) {
        frp_sum[i] <- 0; frp_max[i] <- 0; conf_mean[i] <- 0; hi_n[i] <- 0L
      } else {
        frp_vals <- hs$.frp[idx]
        conf_vals <- hs$.conf01[idx]
        frp_sum[i] <- sum(frp_vals, na.rm = TRUE)
        frp_max[i] <- suppressWarnings(max(frp_vals, na.rm = TRUE)); if (!is.finite(frp_max[i])) frp_max[i] <- 0
        cm <- mean(conf_vals, na.rm = TRUE); conf_mean[i] <- if (is.finite(cm)) cm else 0
        hi_n[i] <- sum(conf_vals >= hi_conf_thr, na.rm = TRUE)
      }
    }
    
    df_polys2 <- tibble(
      !!id_col := as.character(polys2[[id_col]]),
      hotspot_available = 1L,
      hs_in_poly   = as.integer(n_inside),
      hs_in_buffer = as.integer(n_buffer),
      hs_used_n    = as.integer(ifelse(n_inside > 0, n_inside, ifelse(n_buffer > 0, n_buffer, 0L))),
      hs_min_dist_m = as.numeric(hs_min),
      hs_frp_sum   = as.numeric(frp_sum),
      hs_frp_max   = as.numeric(frp_max),
      hs_conf_mean = as.numeric(conf_mean),
      hs_hiConf_n  = as.integer(hi_n)
    )
    
    base_ids %>%
      dplyr::left_join(df_polys2, by = id_col) %>%
      dplyr::mutate(
        hotspot_available = ifelse(is.na(hotspot_available), 1L, hotspot_available),
        hs_in_poly   = ifelse(is.na(hs_in_poly), 0L, hs_in_poly),
        hs_in_buffer = ifelse(is.na(hs_in_buffer), 0L, hs_in_buffer),
        hs_used_n    = ifelse(is.na(hs_used_n), 0L, hs_used_n),
        hs_min_dist_m = ifelse(is.na(hs_min_dist_m), hs_buffer_m, hs_min_dist_m),
        hs_frp_sum   = ifelse(is.na(hs_frp_sum), 0, hs_frp_sum),
        hs_frp_max   = ifelse(is.na(hs_frp_max), 0, hs_frp_max),
        hs_conf_mean = ifelse(is.na(hs_conf_mean), 0, hs_conf_mean),
        hs_hiConf_n  = ifelse(is.na(hs_hiConf_n), 0L, hs_hiConf_n),
        hs_support_present = ifelse(hotspot_available == 1L & hs_used_n > 0L, 1L, 0L),
        hs_no_support_when_available = ifelse(hotspot_available == 1L & hs_used_n == 0L, 1L, 0L),
        hs_only_buffer_support = ifelse(hotspot_available == 1L & hs_in_poly == 0L & hs_in_buffer > 0L, 1L, 0L)
      )
  }

  # OPTIONAL shape/size feature helper (OtsuFire 0.12.0). Guards empties
  # with the engine's existing safe_make_valid / drop_empty_sf helpers,
  # then delegates the geometry math to the package-level, unit-testable
  # `.of_shape_features()` (defined in internal-sup-residual-cols.R). The
  # right_join to base_ids preserves coverage + order (NA for any id
  # dropped as empty/invalid). CRS is the metric target (EPSG:3035).
  shape_features <- function(polys) {
    all_ids  <- unique(as.character(polys[[id_col]]))
    base_ids <- tibble(!!id_col := all_ids)

    polys  <- safe_make_valid(polys)
    polys2 <- drop_empty_sf(polys)

    if (nrow(polys2) == 0) {
      return(base_ids %>%
               dplyr::mutate(
                 log_area    = NA_real_,
                 perim_m     = NA_real_,
                 compactness = NA_real_,
                 elongation  = NA_real_
               ))
    }

    shp <- .of_shape_features(polys2, id_col = id_col)
    base_ids %>% dplyr::left_join(shp, by = id_col)
  }

  # ---------------------------
  # Checks (inputs)
  # ---------------------------
  if (!inherits(train_folds, "sf") || !inherits(unlabeled, "sf")) {
    stop("train_folds and unlabeled must be sf objects.")
  }
  
  assert_cols(train_folds, c(id_col, class_col))
  assert_cols(train_folds, c("block_id", "fold_rep1"))
  
  if (!id_col %in% names(unlabeled)) {
    msg("unlabeled missing '%s' -> creating sequential patch id", id_col)
    unlabeled[[id_col]] <- seq_len(nrow(unlabeled))
  }
  
  train_folds <- train_folds %>% dplyr::filter(.data[[class_col]] %in% c(pos_lab, neg_lab))
  
  if (anyDuplicated(as.character(train_folds[[id_col]])) > 0) {
    stop("PATCH MODE: id_col has duplicates in train_folds. Use a unique patch id (e.g., id_col='poly_id').")
  }
  if (anyDuplicated(as.character(unlabeled[[id_col]])) > 0) {
    stop("PATCH MODE: id_col has duplicates in unlabeled.")
  }

  if (build_features) {
    if (is.null(rbr_summer) || is.null(dem) || is.null(slope) || is.null(corine_r)) {
      stop("build_features=TRUE requires rbr_summer, dem, slope, corine_r.")
    }
    if (is.null(cor_groups) || !is.list(cor_groups) || is.null(names(cor_groups))) {
      stop("Please provide cor_groups as a NAMED list.")
    }
  }
  
  # ---------------------------
  # CRS handling (sf)
  # ---------------------------
  train_folds <- to_crs_safe_sf(train_folds, crs_sf, "train_folds")
  unlabeled   <- to_crs_safe_sf(unlabeled,   crs_sf, "unlabeled")
  if (!is.null(hotspots)) hotspots <- to_crs_safe_sf(hotspots, crs_sf, "hotspots")
  
  # infer year_target once
  if (is.null(year_target) && year_col %in% names(train_folds)) {
    year_target <- as_int_safe(train_folds[[year_col]][1])
    if (isTRUE(verbose)) msg("year_target inferred from train_folds[['%s']] = %s", year_col, year_target)
  }
  
  # ---------------------------
  # CRS + grid handling (rasters)
  # ---------------------------
  if (build_features) {
    # Choose template:
    # - user template_raster if provided
    # - else rbr_summer
    template <- template_raster %||% rbr_summer
    if (is.null(template) || !inherits(template, "SpatRaster")) stop("template_raster/rbr_summer must be a SpatRaster.")
    
    # Ensure template is in target CRS (it defines grid afterwards)
    template <- ensure_raster_crs_and_grid(
      r = template,
      template = NULL,
      name = "template_raster",
      method_project = "bilinear",
      method_resample = "bilinear"
    )
    
    # If template came from rbr_summer, sync it
    rbr_summer <- template
    
    # Now make all rasters match CRS + grid (project + resample)
    dem      <- ensure_raster_crs_and_grid(dem,      template, "dem",      "bilinear", "bilinear")
    slope    <- ensure_raster_crs_and_grid(slope,    template, "slope",    "bilinear", "bilinear")
    corine_r <- ensure_raster_crs_and_grid(corine_r, template, "corine_r", "near",     "near")
    
    if (!is.null(doy_post)) doy_post <- ensure_raster_crs_and_grid(doy_post, template, "doy_post", "near",     "near")
    if (!is.null(rbr_aw))   rbr_aw   <- ensure_raster_crs_and_grid(rbr_aw,   template, "rbr_aw",   "bilinear", "bilinear")
    if (!is.null(nbr_pre))  nbr_pre  <- ensure_raster_crs_and_grid(nbr_pre,  template, "nbr_pre",  "bilinear", "bilinear")
    if (!is.null(nbr_post)) nbr_post <- ensure_raster_crs_and_grid(nbr_post, template, "nbr_post", "bilinear", "bilinear")
    if (!is.null(dnbr))     dnbr     <- ensure_raster_crs_and_grid(dnbr,     template, "dnbr",     "bilinear", "bilinear")
  }
  
  # ---------------------------
  # PATCH-level meta
  # ---------------------------
  # 0.4.0 architectural refactor: NO column drops here. All input
  # columns flow through unchanged; the 50 generated features are
  # ADDED via the left_joins below. Filtering is done at the
  # last mile in `train_final_model_direct()`.
  train_meta <- drop_geom(train_folds)
  unl_meta   <- drop_geom(unlabeled)
  
  # ---------------------------
  # Build features
  # ---------------------------
  if (!build_features) {
    train_feat <- train_meta
    unl_feat   <- unl_meta
    
    train_feat_tbl <- tibble(!!id_col := as.character(train_folds[[id_col]]))
    unl_feat_tbl   <- tibble(!!id_col := as.character(unlabeled[[id_col]]))
    
  } else {
    msg("Extracting features per PATCH (no dissolve) ...")
    
    # Clean rbr_aw extremes once (if used)
    rbr_aw_clean <- NULL
    if (use_aw && !is.null(rbr_aw)) {
      rbr_aw_clean <- clean_raster_extremes(rbr_aw, min_ok = rbr_aw_min_ok, tag = "rbr_aw")
    }
    
    train_feat_tbl <- tibble(!!id_col := as.character(train_folds[[id_col]])) %>%
      dplyr::left_join(zonal_numeric_stats(train_folds, rbr_summer, prefix = "rbr"), by = id_col) %>%
      dplyr::left_join(corine_group_fractions(train_folds, cor_groups), by = id_col)

    train_feat_tbl <- train_feat_tbl %>%
      dplyr::left_join(zonal_numeric_stats(train_folds, dem,   prefix = "elev",  compute_sd = TRUE), by = id_col) %>%
      dplyr::left_join(zonal_numeric_stats(train_folds, slope, prefix = "slope", compute_sd = TRUE), by = id_col)
    
    if (use_doy && !is.null(doy_post)) {
      train_feat_tbl <- train_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(train_folds, doy_post, prefix = "doy"), by = id_col)
    }
    
    if (use_aw && !is.null(rbr_aw_clean)) {
      train_feat_tbl <- train_feat_tbl %>%
        dplyr::left_join(zonal_numeric_stats(train_folds, rbr_aw_clean, prefix = "rbr_aw"), by = id_col) %>%
        dplyr::mutate(
          persist_delta = rbr_aw_med - rbr_med,
          persist_ratio = (rbr_aw_med + 1e-6) / (rbr_med + 1e-6)
        )
    }
    
    if (use_nbr) {
      if (!is.null(nbr_pre))  train_feat_tbl <- train_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(train_folds, nbr_pre,  prefix = "nbr_pre"),  by = id_col)
      if (!is.null(nbr_post)) train_feat_tbl <- train_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(train_folds, nbr_post, prefix = "nbr_post"), by = id_col)
      if (!is.null(dnbr))     train_feat_tbl <- train_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(train_folds, dnbr,     prefix = "dnbr"),     by = id_col)
    }
    
    if (use_hotspots) {
      train_feat_tbl <- train_feat_tbl %>% dplyr::left_join(hotspot_features(train_folds), by = id_col)
    }

    # OPTIONAL shape/size block (TRAIN). OFF by default -> no join -> the
    # train feature table is byte-identical to today.
    if (use_shape) {
      train_feat_tbl <- dplyr::left_join(train_feat_tbl, shape_features(train_folds), by = id_col)
    }

    unl_feat_tbl <- tibble(!!id_col := as.character(unlabeled[[id_col]])) %>%
      dplyr::left_join(zonal_numeric_stats(unlabeled, rbr_summer, prefix = "rbr"), by = id_col) %>%
      dplyr::left_join(corine_group_fractions(unlabeled, cor_groups), by = id_col)

    unl_feat_tbl <- unl_feat_tbl %>%
      dplyr::left_join(zonal_numeric_stats(unlabeled, dem,   prefix = "elev",  compute_sd = TRUE), by = id_col) %>%
      dplyr::left_join(zonal_numeric_stats(unlabeled, slope, prefix = "slope", compute_sd = TRUE), by = id_col)
    
    if (use_doy && !is.null(doy_post)) {
      unl_feat_tbl <- unl_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(unlabeled, doy_post, prefix = "doy"), by = id_col)
    }
    
    if (use_aw && !is.null(rbr_aw_clean)) {
      unl_feat_tbl <- unl_feat_tbl %>%
        dplyr::left_join(zonal_numeric_stats(unlabeled, rbr_aw_clean, prefix = "rbr_aw"), by = id_col) %>%
        dplyr::mutate(
          persist_delta = rbr_aw_med - rbr_med,
          persist_ratio = (rbr_aw_med + 1e-6) / (rbr_med + 1e-6)
        )
    }
    
    if (use_nbr) {
      if (!is.null(nbr_pre))  unl_feat_tbl <- unl_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(unlabeled, nbr_pre,  prefix = "nbr_pre"),  by = id_col)
      if (!is.null(nbr_post)) unl_feat_tbl <- unl_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(unlabeled, nbr_post, prefix = "nbr_post"), by = id_col)
      if (!is.null(dnbr))     unl_feat_tbl <- unl_feat_tbl %>% dplyr::left_join(zonal_numeric_stats(unlabeled, dnbr,     prefix = "dnbr"),     by = id_col)
    }
    
    if (use_hotspots) {
      unl_feat_tbl <- unl_feat_tbl %>% dplyr::left_join(hotspot_features(unlabeled), by = id_col)
    }

    # OPTIONAL shape/size block (SCORING). Same single extract_features()
    # invocation feeds both train + scoring pools, so both layers get the
    # columns. OFF by default -> no join -> byte-identical scoring table.
    if (use_shape) {
      unl_feat_tbl <- dplyr::left_join(unl_feat_tbl, shape_features(unlabeled), by = id_col)
    }

    # When the shape block is ON, the four computed columns (log_area,
    # perim_m, compactness, elongation) may already exist on the input
    # pools as legacy deterministic-residual columns. Dropping the input
    # copies before the join makes the freshly-computed shape values
    # authoritative and avoids dplyr `.x` / `.y` name collisions. With
    # use_shape = FALSE these columns are never in the feature table, so
    # this prunes nothing and the join is byte-identical to today.
    if (use_shape) {
      .shape_computed <- c("log_area", "perim_m", "compactness", "elongation")
      train_meta  <- train_meta[, setdiff(names(train_meta), .shape_computed), drop = FALSE]
      unl_meta    <- unl_meta[,   setdiff(names(unl_meta),   .shape_computed), drop = FALSE]
      train_folds <- train_folds[, setdiff(names(train_folds), .shape_computed), drop = FALSE]
      unlabeled   <- unlabeled[,   setdiff(names(unlabeled),   .shape_computed), drop = FALSE]
    }

    train_feat <- train_meta %>% dplyr::left_join(train_feat_tbl, by = id_col)
    unl_feat   <- unl_meta   %>% dplyr::left_join(unl_feat_tbl,   by = id_col)
  }

  # ---------------------------
  # Optionally add geometry back
  # IMPORTANT: join ONLY pure feature tables to avoid duplicate columns
  # ---------------------------
  if (isTRUE(return_features_geometry) || isTRUE(save_features_gpkg)) {
    # 0.4.0: ADD-only — preserve every input column verbatim and
    # left-join the generated feature table.
    train_feat_geom <- train_folds %>% dplyr::left_join(train_feat_tbl, by = id_col)
    unl_feat_geom   <- unlabeled   %>% dplyr::left_join(unl_feat_tbl,   by = id_col)
  } else {
    train_feat_geom <- NULL
    unl_feat_geom   <- NULL
  }
  
  # ---------------------------
  # Save (optional)
  # ---------------------------
  if (!is.null(save_features_dir)) {
    dir.create(save_features_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (save_features_format == "rds") {
      saveRDS(train_feat, file.path(save_features_dir, paste0(train_features_basename, ".rds")))
      saveRDS(unl_feat,   file.path(save_features_dir, paste0(scoring_features_basename, ".rds")))
    } else {
      utils::write.csv(train_feat, file.path(save_features_dir, paste0(train_features_basename, ".csv")), row.names = FALSE)
      utils::write.csv(unl_feat,   file.path(save_features_dir, paste0(scoring_features_basename, ".csv")), row.names = FALSE)
    }
    
    if (isTRUE(save_features_gpkg)) {
      if (is.null(train_feat_geom) || is.null(unl_feat_geom)) {
        stop("save_features_gpkg=TRUE requires return_features_geometry=TRUE (sf outputs).")
      }
      gpkg <- file.path(save_features_dir, "features_geometry.gpkg")
      if (file.exists(gpkg)) {
        unlink(gpkg, force = TRUE)
        if (file.exists(gpkg)) file.remove(gpkg)
        if (file.exists(gpkg)) {
          stop("Could not overwrite existing features_geometry.gpkg: ", gpkg)
        }
      }
      sf::st_write(train_feat_geom, gpkg, layer = train_features_layer, delete_dsn = TRUE, quiet = TRUE)
      sf::st_write(unl_feat_geom,   gpkg, layer = scoring_features_layer, append = TRUE, quiet = TRUE)
      gpkg_layers <- tryCatch(as.character(sf::st_layers(gpkg)$name), error = function(e) character(0))
      if (!all(c(train_features_layer, scoring_features_layer) %in% gpkg_layers)) {
        stop("features_geometry.gpkg was written without the required train/scoring layers.")
      }
    }
    
    msg("Saved features in: %s", save_features_dir)
  }
  
  # ---------------------------
  # Return
  # ---------------------------
  if (!isTRUE(return_features)) {
    return(invisible(NULL))
  }
  
  list(
    train_feat = if (isTRUE(return_features_geometry)) train_feat_geom else train_feat,
    unl_feat   = if (isTRUE(return_features_geometry)) unl_feat_geom   else unl_feat
  )
}
