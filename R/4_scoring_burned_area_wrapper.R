# -------------------------------------------------------------------------
# Burned area scoring + optional validation (Stage-2 refinement)
# -------------------------------------------------------------------------
# This file is written to be package-safe:
# - No library() calls
# - Uses explicit namespaces (sf::, terra::, stats::)
#
# Exported functions:
# - erase_overlap_with_preyear()
# - scoring_internal_burned_area()
# - scoring_reference_burned_area()
# - scoring_validation_burned_area()
# - scoring_burned_area_stage2()  (wrapper: internal + reference scoring + optional validation)
# -------------------------------------------------------------------------

# ============================
# Internal helper: suppress only common sf warnings we want to ignore
# ============================
.sf_quiet <- function(expr, enable = TRUE) {
  if (!isTRUE(enable)) return(eval.parent(substitute(expr)))

  patterns <- c(
    "x is already of type POLYGON",
    "attribute variables are assumed to be spatially constant throughout all geometries"
  )

  withCallingHandlers(
    eval.parent(substitute(expr)),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (any(vapply(patterns, function(p) grepl(p, msg, fixed = TRUE), logical(1)))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# ============================
# Internal helper: robust polygon sf (valid, polygonal, multipolygon, non-empty)
# ============================
.safe_make_valid_polys <- function(x, suppress_sf_warnings = TRUE) {
  if (is.null(x)) return(NULL)

  x <- sf::st_as_sf(x)
  if (is.na(sf::st_crs(x))) stop("Input polygons have NA CRS. Assign CRS before processing.")

  # drop empty
  emp <- sf::st_is_empty(sf::st_geometry(x))
  if (any(emp)) x <- x[!emp, , drop = FALSE]
  if (nrow(x) == 0) return(NULL)

  # make valid only where needed
  ok <- .sf_quiet(sf::st_is_valid(x), enable = suppress_sf_warnings)
  ok[is.na(ok)] <- FALSE
  if (!all(ok)) {
    x_bad <- x[!ok, , drop = FALSE]
    x[!ok, ] <- .sf_quiet(sf::st_make_valid(x_bad), enable = suppress_sf_warnings)
  }

  # keep polygonal only (avoid "already POLYGON" warnings)
  gt <- as.character(sf::st_geometry_type(sf::st_geometry(x), by_geometry = TRUE))
  if (!all(gt %in% c("POLYGON", "MULTIPOLYGON"))) {
    x <- .sf_quiet(sf::st_collection_extract(x, "POLYGON", warn = FALSE), enable = suppress_sf_warnings)
    if (is.null(x) || nrow(x) == 0) return(NULL)
  }

  # cast geometry only (avoid attribute warnings)
  geom <- .sf_quiet(sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE), enable = suppress_sf_warnings)
  sf::st_geometry(x) <- geom

  emp2 <- sf::st_is_empty(sf::st_geometry(x))
  if (any(emp2)) x <- x[!emp2, , drop = FALSE]
  if (nrow(x) == 0) return(NULL)

  x
}

# ============================
# Internal helper: detect binary-like raster values (0/1)
# ============================
#' @keywords internal
.is_binary_like <- function(x) {
  if (is.null(x)) return(FALSE)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(FALSE)
  ux <- unique(x)
  all(ux %in% c(0, 1))
}

# ============================
# Internal helper: add per-polygon median from a single-layer raster (chunked)
# ============================
#' @keywords internal
.add_polygon_median_rbr <- function(polys_sf, rbr_rast, chunk_size = 500L, colname = "median_rbr") {
  if (is.null(polys_sf) || nrow(polys_sf) == 0) {
    if (!is.null(polys_sf) && !colname %in% names(polys_sf)) polys_sf[[colname]] <- numeric(0)
    return(polys_sf)
  }
  if (!inherits(rbr_rast, "SpatRaster")) stop("rbr_rast must be a terra::SpatRaster.")
  if (terra::nlyr(rbr_rast) > 1) rbr_rast <- rbr_rast[[1]]

  out <- polys_sf
  out[[colname]] <- NA_real_

  v_all <- terra::vect(out)
  if (!terra::same.crs(v_all, rbr_rast)) v_all <- terra::project(v_all, terra::crs(rbr_rast))

  n <- nrow(out)
  idx_all <- seq_len(n)
  chunk_size <- as.integer(chunk_size)
  if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- 500L

  med_fun <- function(x, ...) {
    if (length(x) == 0) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  }

  starts <- seq(1L, n, by = chunk_size)
  for (s in starts) {
    e <- min(n, s + chunk_size - 1L)
    ii <- idx_all[s:e]
    v  <- v_all[ii]

    ex <- terra::extract(rbr_rast, v, fun = med_fun, exact = FALSE)
    ex <- as.data.frame(ex)

    # ex: ID (1..len(ii)), value
    if (nrow(ex) > 0) {
      ids  <- as.integer(ex[[1]])
      vals <- suppressWarnings(as.numeric(ex[[2]]))
      vals[!is.finite(vals)] <- NA_real_
      out[[colname]][ii[ids]] <- vals
    }
  }

  out
}

#' Erase overlaps with a pre-year fire mask (robust, warning-safe)
#'
#' Removes (differences out) the areas of `stage2_sf` overlapping `preyear_sf`,
#' optionally buffering the mask outward and shaving thin strips after erase.
#'
#' Notes:
#' - This function performs `st_difference()` only on polygons that actually intersect the mask,
#'   and keeps non-intersecting polygons unchanged (fast + avoids unnecessary geometry churn).
#' - If `post_shave_m > 0`, the result is buffered inward then outward to remove thin slivers.
#'
#' @param stage2_sf sf. Polygons to erase.
#' @param preyear_sf sf. Mask polygons (e.g., year-1 fires).
#' @param mask_buffer_m numeric. Outward buffer on mask (meters).
#' @param post_shave_m numeric. Shave width (meters) after erase (negative then positive buffer).
#' @param min_area_m2 numeric. Drop fragments smaller than this area after erase.
#' @param quiet logical. If `FALSE`, prints a message when no overlaps are found.
#' @param suppress_sf_warnings logical. If `TRUE`, suppress common sf warnings.
#'
#' @return sf. Same attributes as `stage2_sf`, possibly with rows expanded or removed after erase.
#' @export
erase_overlap_with_preyear <- function(
    stage2_sf,
    preyear_sf,
    mask_buffer_m = 0,
    post_shave_m  = 0,
    min_area_m2   = 0,
    quiet = TRUE,
    suppress_sf_warnings = TRUE
) {
  stopifnot(inherits(stage2_sf, "sf"), inherits(preyear_sf, "sf"))
  if (nrow(stage2_sf) == 0) return(stage2_sf)

  # CRS align (preyear -> stage2)
  if (sf::st_crs(stage2_sf) != sf::st_crs(preyear_sf)) {
    preyear_sf <- sf::st_transform(preyear_sf, sf::st_crs(stage2_sf))
  }

  # ------------------------------------------------------------
  # Auto-project if CRS is longlat and distances/areas are in meters
  # ------------------------------------------------------------
  did_project <- FALSE
  crs_in <- sf::st_crs(stage2_sf)

  needs_metric <- (is.numeric(mask_buffer_m) && mask_buffer_m > 0) ||
    (is.numeric(post_shave_m)  && post_shave_m  > 0) ||
    (is.numeric(min_area_m2)   && min_area_m2   > 0)

  if (isTRUE(needs_metric) && isTRUE(sf::st_is_longlat(stage2_sf))) {

    # Pick a local UTM zone from bbox center (works globally)
    bb <- sf::st_bbox(stage2_sf)
    lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
    lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)

    # Guard against weird bbox/NA
    if (!is.finite(lon) || !is.finite(lat)) {
      # fallback: EPSG:3857 (meters-ish) if bbox is unusable
      target_crs <- sf::st_crs(3857)
    } else {
      zone <- floor((lon + 180) / 6) + 1
      zone <- max(1, min(60, zone))
      epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
      target_crs <- sf::st_crs(epsg)
    }

    if (!isTRUE(quiet)) {
      msg_epsg <- if (!is.na(target_crs$epsg)) target_crs$epsg else "custom"
      message("Input CRS is longlat; projecting to metric CRS (EPSG:", msg_epsg, ") for erase/buffer/area ops.")
    }

    stage2_sf  <- sf::st_transform(stage2_sf,  target_crs)
    preyear_sf <- sf::st_transform(preyear_sf, target_crs)
    did_project <- TRUE
  }

  # Robustify geometries
  stage2_sf  <- .safe_make_valid_polys(stage2_sf,  suppress_sf_warnings = suppress_sf_warnings)
  preyear_sf <- .safe_make_valid_polys(preyear_sf, suppress_sf_warnings = suppress_sf_warnings)
  if (is.null(stage2_sf)  || nrow(stage2_sf) == 0) return(stage2_sf)
  if (is.null(preyear_sf) || nrow(preyear_sf) == 0) return(stage2_sf)

  # Union mask geometry
  mask_union <- .sf_quiet(sf::st_union(sf::st_geometry(preyear_sf)), enable = suppress_sf_warnings)
  mask_union <- .sf_quiet(sf::st_make_valid(mask_union), enable = suppress_sf_warnings)
  if (length(mask_union) > 1) mask_union <- .sf_quiet(sf::st_union(mask_union), enable = suppress_sf_warnings)

  # Optional outward buffer on mask
  if (is.numeric(mask_buffer_m) && mask_buffer_m > 0) {
    mask_union <- .sf_quiet(sf::st_buffer(mask_union, dist = mask_buffer_m), enable = suppress_sf_warnings)
    mask_union <- .sf_quiet(sf::st_make_valid(mask_union), enable = suppress_sf_warnings)
  }

  hits <- lengths(sf::st_intersects(stage2_sf, mask_union)) > 0
  if (!any(hits)) {
    if (!isTRUE(quiet)) message("No overlaps found. Returning unchanged.")
    # return in original CRS if we projected
    out0 <- stage2_sf
    if (isTRUE(did_project) && !is.na(crs_in)) out0 <- sf::st_transform(out0, crs_in)
    return(out0)
  }

  keep_nohit <- stage2_sf[!hits, , drop = FALSE]
  to_erase   <- stage2_sf[hits,  , drop = FALSE]

  # IMPORTANT: st_difference on sf (can expand rows)
  erased <- .sf_quiet(sf::st_difference(to_erase, mask_union), enable = suppress_sf_warnings)
  erased <- erased[!sf::st_is_empty(sf::st_geometry(erased)), , drop = FALSE]
  erased <- .safe_make_valid_polys(erased, suppress_sf_warnings = suppress_sf_warnings)

  # Optional shave (meters; only meaningful in metric CRS)
  if (!is.null(erased) && nrow(erased) > 0 && is.numeric(post_shave_m) && post_shave_m > 0) {
    er <- .sf_quiet(sf::st_buffer(erased, dist = -post_shave_m), enable = suppress_sf_warnings)
    er <- er[!sf::st_is_empty(sf::st_geometry(er)), , drop = FALSE]
    if (nrow(er) > 0) {
      er <- .sf_quiet(sf::st_buffer(er, dist = post_shave_m), enable = suppress_sf_warnings)
      er <- .safe_make_valid_polys(er, suppress_sf_warnings = suppress_sf_warnings)
      erased <- er
    } else {
      erased <- erased[0, , drop = FALSE]
    }
  }

  # Area filter (m2; meaningful in metric CRS)
  if (!is.null(erased) && nrow(erased) > 0 && is.numeric(min_area_m2) && min_area_m2 > 0) {
    a <- as.numeric(sf::st_area(erased))
    erased <- erased[a >= min_area_m2, , drop = FALSE]
  }

  out <- if (is.null(erased) || nrow(erased) == 0) keep_nohit else rbind(keep_nohit, erased)
  out <- .safe_make_valid_polys(out, suppress_sf_warnings = suppress_sf_warnings)

  # Back-transform to original CRS if we projected
  if (isTRUE(did_project) && !is.na(crs_in)) {
    out <- sf::st_transform(out, crs_in)
    out <- .safe_make_valid_polys(out, suppress_sf_warnings = suppress_sf_warnings)
  }

  out
}
#' Score and flag Stage-2 burned area polygons using Stage-1 support and a burnable CORINE mask
#'
#' Flags Stage-2 polygons as `keep`, `review`, or `drop_*` using:
#' 1) intersection with Stage-1 support (seeds): keep vs nohit_stage1
#' 2) burnable mask check for no-hit polygons: review vs drop_unburnable
#'
#' Optionally:
#' - Computes a CORINE NA fraction diagnostic (`corine_na_frac`) and can re-flag polygons above a threshold.
#' - Computes per-polygon median values from an RBR raster (default: `median_rbr`) for keep/review outputs.
#' - Erases overlaps with a pre-year mask for keep/review outputs.
#' - Saves outputs to disk (GPKG or Shapefile).
#'
#' @param polys_stage2 sf. Stage-2 polygons to flag.
#' @param polys_stage1 sf. Stage-1 polygons used as support.
#' @param burnable_corine Burnable mask. Either a terra::SpatRaster (CORINE-like) or an sf layer.
#' @param burnable_classes Numeric/integer vector of burnable classes when `burnable_corine` is categorical.
#'   If `NULL`, the raster must be binary-like (values in {0,1}).
#' @param rbr_rast Optional terra::SpatRaster. If provided and `compute_median=TRUE`, median RBR is computed.
#' @param support_buffer_m Numeric. Buffer (meters) applied to Stage-1 support before intersects.
#' @param max_corine_na_frac Numeric or `NULL`. If set, polygons with NA fraction greater than this
#'   are re-flagged as `"drop_corine_na"`.
#' @param corine_na_scope Character. Where to compute CORINE NA fraction. One of:
#'   `"nohit_stage1"`, `"keep"`, `"keep_nohit"`, `"all"`.
#' @param support_method Character. `"raster"` (fast if burnable_corine is raster) or `"vector"`.
#' @param touches Logical. If `TRUE`, terra::rasterize uses touches=TRUE to reduce boundary gaps.
#' @param make_valid Character. `"auto"`, `"always"`, or `"none"` for validity fixing during intermediate steps.
#' @param na_exact Logical. If `TRUE`, compute CORINE NA fraction using exact extraction (slower).
#' @param burnable_exact Logical. If `TRUE`, compute burnable check using exact extraction (slower).
#' @param clip_to_burnable Logical. If `TRUE`, clip Stage-2 polygons to burnable mask (raster required).
#' @param clip_mask_buffer_m Numeric. Optional dilation buffer (meters) on burnable mask during clipping.
#' @param clip_min_piece_area_m2 Numeric. Drop clipped slivers smaller than this area.
#' @param clip_keep_largest_piece Logical. If `TRUE`, keep only largest clipped piece per original polygon.
#' @param preyear_polys Optional sf. Pre-year polygons to erase from outputs.
#' @param erase_mask_buffer_m Numeric. Buffer (meters) applied outward to erase mask before difference.
#' @param erase_post_shave_m Numeric. Optional shave width (meters) after erase to remove thin strips.
#' @param erase_min_area_m2 Numeric. Drop erased fragments smaller than this area.
#' @param compute_median Logical. If `TRUE` and `rbr_rast` is provided, compute median RBR.
#' @param median_colname Character. Name of the median column to create/fill.
#' @param median_chunk_size Integer. Chunk size for polygon extraction.
#' @param save_outputs Logical. If `TRUE`, write outputs and a summary file.
#' @param out_dir Character. Output directory used when `save_outputs=TRUE`.
#' @param prefix Character. Prefix for output filenames when saving.
#' @param driver Character. `"GPKG"` or `"ESRI Shapefile"`.
#' @param overwrite Logical. Overwrite outputs when saving.
#' @param quiet Logical. If `TRUE`, suppress messages in sf::st_write.
#' @param save_splits Logical. If `TRUE`, also save keep-only and review-only layers.
#' @param suppress_sf_warnings Logical. If `TRUE`, suppress common sf warnings (recommended).
#'
#' @return A list with:
#' - stage2_flagged, stage2_flagged_clean, stage2_keep_review, stage2_keep_review_clean, stage2_dropped, stage2_keep_only, stage2_keep_only_clean, stage2_review_only, stage2_review_only_clean
#' - stage2_keep_review_erased, stage2_keep_only_erased, stage2_review_only_erased (if preyear_polys provided)
#' - out_paths (if save_outputs=TRUE)
#'
#' @export
scoring_internal_burned_area <- function(
  polys_stage2,
  polys_stage1,
  burnable_corine,
  burnable_classes = NULL,

  rbr_rast = NULL,

  support_buffer_m = 0,

  max_corine_na_frac = 0.7,
  corine_na_scope = c("nohit_stage1", "keep", "keep_nohit", "all"),

  support_method = c("raster", "vector"),
  touches = TRUE,
  make_valid = c("auto", "always", "none"),
  na_exact = FALSE,
  burnable_exact = FALSE,

  clip_to_burnable = FALSE,
  clip_mask_buffer_m = 0,
  clip_min_piece_area_m2 = 0,
  clip_keep_largest_piece = FALSE,

  preyear_polys = NULL,
  erase_mask_buffer_m = 0,
  erase_post_shave_m  = 0,
  erase_min_area_m2   = 0,

  compute_median = TRUE,
  median_colname = "median_rbr",
  median_chunk_size = 500L,

  save_outputs = FALSE,
  out_dir = NULL,
  prefix = "BA",
  driver = "GPKG",
  overwrite = TRUE,
  quiet = TRUE,
  save_splits = FALSE,

  suppress_sf_warnings = TRUE
){

  pkgs <- c("sf", "terra", "stats")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing required packages: ", paste(miss, collapse = ", "))

  corine_na_scope <- match.arg(corine_na_scope)
  support_method  <- match.arg(support_method)
  make_valid      <- match.arg(make_valid)

  make_valid_selective <- function(x) {
    if (make_valid == "none") return(x)
    if (make_valid == "always") return(.sf_quiet(sf::st_make_valid(x), enable = suppress_sf_warnings))
    ok <- .sf_quiet(sf::st_is_valid(x), enable = suppress_sf_warnings)
    ok[is.na(ok)] <- FALSE
    if (all(ok)) return(x)
    x_bad <- x[!ok, , drop = FALSE]
    x[!ok, ] <- .sf_quiet(sf::st_make_valid(x_bad), enable = suppress_sf_warnings)
    x
  }

  # Validate inputs
  stopifnot(inherits(polys_stage2, "sf"), inherits(polys_stage1, "sf"))
  polys_stage2 <- .safe_make_valid_polys(polys_stage2, suppress_sf_warnings = suppress_sf_warnings)
  polys_stage1 <- .safe_make_valid_polys(polys_stage1, suppress_sf_warnings = suppress_sf_warnings)
  if (is.null(polys_stage2) || nrow(polys_stage2) == 0) stop("polys_stage2 is empty after validation.")
  if (is.null(polys_stage1) || nrow(polys_stage1) == 0) stop("polys_stage1 is empty after validation.")

  # Optional: clip stage2 to burnable first (fast raster method)
  if (isTRUE(clip_to_burnable)) {
    if (!inherits(burnable_corine, "SpatRaster")) {
      stop("clip_to_burnable=TRUE requires burnable_corine as a terra::SpatRaster.")
    }

    clip_fast <- function(polys_sf, burnable_corine_rast, burnable_classes,
                          mask_buffer_m, min_piece_area_m2, keep_largest_piece, touches) {

      if (nrow(polys_sf) == 0) return(polys_sf)
      crs0 <- sf::st_crs(polys_sf)
      if (is.na(crs0)) stop("polys_sf CRS is NA.")

      polys_sf <- .safe_make_valid_polys(polys_sf, suppress_sf_warnings = suppress_sf_warnings)
      if (is.null(polys_sf) || nrow(polys_sf) == 0) return(polys_sf)

      if (!"._row_id" %in% names(polys_sf)) polys_sf$._row_id <- seq_len(nrow(polys_sf))

      r <- burnable_corine_rast[[1]]

      v <- terra::vect(polys_sf)
      if (!terra::same.crs(v, r)) v <- terra::project(v, r)

      e <- terra::ext(v)

      if (is.numeric(mask_buffer_m) && mask_buffer_m > 0) {
        res_xy <- terra::res(r)
        px <- min(res_xy, na.rm = TRUE)
        buf <- max(px, mask_buffer_m)
        e <- terra::ext(e[1] - buf, e[2] + buf, e[3] - buf, e[4] + buf)
      }

      r_sub <- terra::crop(r, e)
      vals <- terra::values(r_sub, mat = FALSE)
      if (all(is.na(vals))) return(polys_sf[0, , drop = FALSE])

      if (is.null(burnable_classes)) {
        if (!.is_binary_like(vals[!is.na(vals)])) {
          stop("burnable_classes is NULL but burnable_corine is not binary-like.")
        }
        burn <- terra::ifel(r_sub == 1, 1, NA)
      } else {
        burn <- terra::ifel(r_sub %in% burnable_classes, 1, NA)
      }

      if (all(is.na(terra::values(burn)))) return(polys_sf[0, , drop = FALSE])

      if (is.numeric(mask_buffer_m) && mask_buffer_m > 0) {
        res_xy <- terra::res(burn)
        px <- min(res_xy, na.rm = TRUE)
        npx <- max(1L, as.integer(ceiling(mask_buffer_m / px)))
        w <- matrix(1, nrow = 2 * npx + 1, ncol = 2 * npx + 1)
        burn <- terra::focal(burn, w = w, fun = "max", na.policy = "omit", fillvalue = NA)
        burn <- terra::ifel(burn >= 1, 1, NA)
      }

      id_r <- terra::rasterize(v, burn, field = "._row_id", touches = touches)
      keep_id <- terra::ifel(!is.na(burn) & !is.na(id_r), id_r, NA)
      if (all(is.na(terra::values(keep_id)))) return(polys_sf[0, , drop = FALSE])

      pv <- terra::as.polygons(keep_id, dissolve = TRUE, na.rm = TRUE)
      if (is.null(pv) || nrow(pv) == 0) return(polys_sf[0, , drop = FALSE])

      names(pv)[1] <- "._row_id"
      pv_sf <- sf::st_as_sf(pv)

      attrs <- sf::st_drop_geometry(polys_sf)
      attrs <- attrs[match(pv_sf$._row_id, attrs$._row_id), , drop = FALSE]
      pv_sf <- cbind(pv_sf, attrs[setdiff(names(attrs), "._row_id")])

      if (is.numeric(min_piece_area_m2) && min_piece_area_m2 > 0 && nrow(pv_sf) > 0) {
        a <- as.numeric(sf::st_area(pv_sf))
        pv_sf <- pv_sf[a >= min_piece_area_m2, , drop = FALSE]
      }
      if (nrow(pv_sf) == 0) return(polys_sf[0, , drop = FALSE])

      if (isTRUE(keep_largest_piece)) {
        pv_sf$._area <- as.numeric(sf::st_area(pv_sf))
        pv_sf <- pv_sf[order(pv_sf$._row_id, -pv_sf$._area), ]
        pv_sf <- pv_sf[!duplicated(pv_sf$._row_id), , drop = FALSE]
        pv_sf$._area <- NULL
      }

      if (sf::st_crs(pv_sf) != crs0) pv_sf <- sf::st_transform(pv_sf, crs0)
      pv_sf
    }

    polys_stage2 <- clip_fast(
      polys_sf = polys_stage2,
      burnable_corine_rast = burnable_corine,
      burnable_classes = burnable_classes,
      mask_buffer_m = clip_mask_buffer_m,
      min_piece_area_m2 = clip_min_piece_area_m2,
      keep_largest_piece = clip_keep_largest_piece,
      touches = touches
    )

    polys_stage2 <- .safe_make_valid_polys(polys_stage2, suppress_sf_warnings = suppress_sf_warnings)
    if (is.null(polys_stage2) || nrow(polys_stage2) == 0) stop("After clip_to_burnable, polys_stage2 is empty.")
  }

  # CRS align stage1 -> stage2
  if (sf::st_crs(polys_stage2) != sf::st_crs(polys_stage1)) {
    polys_stage1 <- sf::st_transform(polys_stage1, sf::st_crs(polys_stage2))
  }

  stage2_valid <- make_valid_selective(polys_stage2)
  stage1_valid <- make_valid_selective(polys_stage1)

  # Optional buffer on stage1 support
  stage1_support <- stage1_valid
  if (is.numeric(support_buffer_m) && support_buffer_m > 0) {
    g <- .sf_quiet(sf::st_buffer(sf::st_geometry(stage1_support), dist = support_buffer_m), enable = suppress_sf_warnings)
    stage1_support <- sf::st_sf(geometry = g)
    sf::st_crs(stage1_support) <- sf::st_crs(stage2_valid)
    stage1_support <- make_valid_selective(stage1_support)
  }

  # Prepare burnable raster logic (if raster)
  burnable_r_na   <- NULL
  burnable_r_mask <- NULL
  if (inherits(burnable_corine, "SpatRaster")) {
    burnable_r_na <- burnable_corine[[1]]
    vals0 <- terra::values(burnable_r_na, mat = FALSE)
    vals0 <- vals0[!is.na(vals0)]

    if (is.null(burnable_classes)) {
      if (!.is_binary_like(vals0)) {
        stop("burnable_corine is categorical but burnable_classes is NULL. Provide burnable_classes or use a binary mask.")
      }
      burnable_r_mask <- terra::ifel(burnable_r_na == 1, 1, NA)
    } else {
      burnable_r_mask <- terra::ifel(burnable_r_na %in% burnable_classes, 1, NA)
    }
  }

  # 1) hit_stage1 (keep vs nohit_stage1)
  hit_stage1 <- rep(FALSE, nrow(stage2_valid))

  if (support_method == "raster" && inherits(burnable_corine, "SpatRaster")) {
    r_template <- burnable_r_na

    v_stage1 <- terra::vect(stage1_support)
    if (!terra::same.crs(v_stage1, r_template)) v_stage1 <- terra::project(v_stage1, terra::crs(r_template))

    r_sup <- terra::rasterize(v_stage1, r_template, field = 1, background = 0, touches = touches)

    v_stage2 <- terra::vect(stage2_valid)
    if (!terra::same.crs(v_stage2, r_template)) v_stage2 <- terra::project(v_stage2, terra::crs(r_template))

    sup_fun <- function(x, ...) as.integer(length(x) > 0 && any(x == 1, na.rm = TRUE))
    ex_sup <- terra::extract(r_sup, v_stage2, fun = sup_fun, exact = FALSE)
    ex_sup <- as.data.frame(ex_sup)

    hit_stage1 <- rep(FALSE, nrow(stage2_valid))
    if (nrow(ex_sup) > 0) {
      ids  <- as.integer(ex_sup[[1]])
      vals <- as.integer(ex_sup[[2]])
      hit_stage1[ids] <- as.logical(vals)
      hit_stage1[is.na(hit_stage1)] <- FALSE
    }
  } else {
    hit_stage1 <- lengths(sf::st_intersects(stage2_valid, stage1_support)) > 0
  }

  out <- stage2_valid
  out$flag_internal  <- ifelse(hit_stage1, "keep", "nohit_stage1")
  out$why_flag       <- ifelse(hit_stage1, "intersects_stage1_seed", "no_intersection_stage1_seed")
  out$corine_na_frac <- NA_real_

  # OPTIONAL: corine_na_frac diagnostic (+ optional drop)
  if (inherits(burnable_corine, "SpatRaster")) {
    idx_na <- switch(
      corine_na_scope,
      "nohit_stage1" = which(out$flag_internal == "nohit_stage1"),
      "keep"         = which(out$flag_internal == "keep"),
      "keep_nohit"   = which(out$flag_internal %in% c("keep", "nohit_stage1")),
      "all"          = seq_len(nrow(out))
    )

    if (length(idx_na) > 0) {
      r_na <- burnable_r_na
      v_na <- terra::vect(out[idx_na, , drop = FALSE])
      if (!terra::same.crs(v_na, r_na)) v_na <- terra::project(v_na, terra::crs(r_na))

      if (!isTRUE(na_exact)) {
        na_fun <- function(x, ...) if (length(x) == 0) 1 else mean(is.na(x))
        ex <- terra::extract(r_na, v_na, fun = na_fun, exact = FALSE)
        ex <- as.data.frame(ex)

        na_frac <- rep(1, length(idx_na))
        if (nrow(ex) > 0) {
          ids  <- as.integer(ex[[1]])
          vals <- as.numeric(ex[[2]])
          vals[!is.finite(vals)] <- 1
          na_frac[ids] <- vals
        }
      } else {
        ex <- terra::extract(r_na, v_na, exact = TRUE)
        ex <- as.data.frame(ex)
        if (!"ID" %in% names(ex)) stop("terra::extract output does not contain 'ID' column.")

        cov_candidates <- c("coverage_fraction", "fraction", "weight", "weights")
        cov_col <- cov_candidates[cov_candidates %in% names(ex)][1]
        if (is.na(cov_col)) cov_col <- NULL

        value_cols <- setdiff(names(ex), c("ID", cov_candidates))
        if (length(value_cols) < 1) stop("Could not identify value column from terra::extract output.")
        val_col <- value_cols[1]

        cov <- if (!is.null(cov_col)) as.numeric(ex[[cov_col]]) else rep(1, nrow(ex))
        val <- ex[[val_col]]
        ids <- as.integer(ex[["ID"]])

        tmp <- tapply(seq_len(nrow(ex)), ids, function(ii) {
          cc <- cov[ii]
          vv <- val[ii]
          tot <- sum(cc, na.rm = TRUE)
          if (!is.finite(tot) || tot <= 0) return(1)
          na_cov <- sum(cc[is.na(vv)], na.rm = TRUE)
          as.numeric(na_cov / tot)
        })

        na_frac <- rep(1, length(idx_na))
        na_frac[as.integer(names(tmp))] <- as.numeric(tmp)
        na_frac[!is.finite(na_frac)] <- 1
      }

      out$corine_na_frac[idx_na] <- na_frac

      if (!is.null(max_corine_na_frac) && is.finite(max_corine_na_frac)) {
        drop_na <- na_frac > max_corine_na_frac
        if (any(drop_na)) {
          out$flag_internal[idx_na[drop_na]] <- "drop_corine_na"
          out$why_flag[idx_na[drop_na]] <- sprintf("corine_na_frac_gt_%.2f", max_corine_na_frac)
        }
      }
    }
  }

  # 2) burnable mask check for no-hit
  idx <- which(out$flag_internal == "nohit_stage1")
  if (length(idx) > 0) {

    if (inherits(burnable_corine, "sf")) {

      burn_sf <- .safe_make_valid_polys(burnable_corine, suppress_sf_warnings = suppress_sf_warnings)
      if (is.null(burn_sf) || nrow(burn_sf) == 0) {
        out$flag_internal[idx] <- "drop_unburnable"
        out$why_flag[idx]      <- "burnable_corine_empty"
      } else {
        if (sf::st_crs(out) != sf::st_crs(burn_sf)) burn_sf <- sf::st_transform(burn_sf, sf::st_crs(out))
        burn_sf <- make_valid_selective(burn_sf)
        hit_burnable <- lengths(sf::st_intersects(out[idx, , drop = FALSE], burn_sf)) > 0
        out$flag_internal[idx] <- ifelse(hit_burnable, "review", "drop_unburnable")
        out$why_flag[idx]      <- ifelse(hit_burnable, "burnable_corine_ok", "outside_burnable_corine")
      }

    } else if (inherits(burnable_corine, "SpatRaster")) {

      r_burn <- burnable_r_mask
      v_burn <- terra::vect(out[idx, , drop = FALSE])
      if (!terra::same.crs(v_burn, r_burn)) v_burn <- terra::project(v_burn, terra::crs(r_burn))

      if (!isTRUE(burnable_exact)) {
        burn_fun <- function(x, ...) as.integer(length(x) > 0 && any(x == 1, na.rm = TRUE))
        ex <- terra::extract(r_burn, v_burn, fun = burn_fun, exact = FALSE)
        ex <- as.data.frame(ex)

        hit_burnable <- rep(FALSE, length(idx))
        if (nrow(ex) > 0) {
          ids  <- as.integer(ex[[1]])
          vals <- as.integer(ex[[2]])
          hit_burnable[ids] <- as.logical(vals)
          hit_burnable[is.na(hit_burnable)] <- FALSE
        }
      } else {
        ex <- terra::extract(r_burn, v_burn, exact = TRUE)
        ex <- as.data.frame(ex)
        if (!"ID" %in% names(ex)) stop("terra::extract output does not contain 'ID' column.")

        value_cols <- setdiff(names(ex), c("ID", "coverage_fraction", "fraction", "weight", "weights"))
        if (length(value_cols) < 1) stop("Could not identify value column from terra::extract output.")
        val <- suppressWarnings(as.integer(ex[[value_cols[1]]]))
        ids <- as.integer(ex[["ID"]])

        tmp <- tapply(seq_len(nrow(ex)), ids, function(ii) any(val[ii] == 1, na.rm = TRUE))
        hit_burnable <- rep(FALSE, length(idx))
        hit_burnable[as.integer(names(tmp))] <- as.logical(tmp)
        hit_burnable[is.na(hit_burnable)] <- FALSE
      }

      out$flag_internal[idx] <- ifelse(hit_burnable, "review", "drop_unburnable")
      out$why_flag[idx]      <- ifelse(hit_burnable, "burnable_corine_ok", "outside_burnable_corine")

    } else {
      stop("burnable_corine must be an sf or a terra::SpatRaster.")
    }
  }

  # Split outputs
  stage2_keep_review <- out[out$flag_internal %in% c("keep", "review"), , drop = FALSE]
  stage2_dropped     <- out[out$flag_internal %in% c("drop_unburnable", "drop_corine_na"), , drop = FALSE]
  stage2_keep_only   <- out[out$flag_internal == "keep", , drop = FALSE]
  stage2_review_only <- out[out$flag_internal == "review", , drop = FALSE]

  # Optional: compute median RBR
  if (isTRUE(compute_median) && !is.null(rbr_rast)) {
    if (!inherits(rbr_rast, "SpatRaster")) stop("rbr_rast must be a terra::SpatRaster.")
    if (terra::nlyr(rbr_rast) > 1) rbr_rast <- rbr_rast[[1]]

    stage2_keep_review <- .add_polygon_median_rbr(stage2_keep_review, rbr_rast, chunk_size = median_chunk_size, colname = median_colname)
    stage2_keep_only   <- .add_polygon_median_rbr(stage2_keep_only,   rbr_rast, chunk_size = median_chunk_size, colname = median_colname)
    stage2_review_only <- .add_polygon_median_rbr(stage2_review_only, rbr_rast, chunk_size = median_chunk_size, colname = median_colname)
  }

  # Optional: erase overlaps with preyear
  stage2_keep_review_erased <- NULL
  stage2_keep_only_erased   <- NULL
  stage2_review_only_erased <- NULL

  if (!is.null(preyear_polys)) {
    preyear_polys <- .safe_make_valid_polys(preyear_polys, suppress_sf_warnings = suppress_sf_warnings)
    if (!is.null(preyear_polys) && nrow(preyear_polys) > 0) {
      stage2_keep_review_erased <- erase_overlap_with_preyear(
        stage2_sf = stage2_keep_review,
        preyear_sf = preyear_polys,
        mask_buffer_m = erase_mask_buffer_m,
        post_shave_m  = erase_post_shave_m,
        min_area_m2   = erase_min_area_m2,
        quiet = quiet,
        suppress_sf_warnings = suppress_sf_warnings
      )
      stage2_keep_only_erased <- erase_overlap_with_preyear(
        stage2_sf = stage2_keep_only,
        preyear_sf = preyear_polys,
        mask_buffer_m = erase_mask_buffer_m,
        post_shave_m  = erase_post_shave_m,
        min_area_m2   = erase_min_area_m2,
        quiet = quiet,
        suppress_sf_warnings = suppress_sf_warnings
      )
      stage2_review_only_erased <- erase_overlap_with_preyear(
        stage2_sf = stage2_review_only,
        preyear_sf = preyear_polys,
        mask_buffer_m = erase_mask_buffer_m,
        post_shave_m  = erase_post_shave_m,
        min_area_m2   = erase_min_area_m2,
        quiet = quiet,
        suppress_sf_warnings = suppress_sf_warnings
      )
    }
  }

  # Optional: save outputs
  out_paths <- list(
    container = NA_character_,
    flagged = NA_character_,
    keep_review = NA_character_,
    dropped = NA_character_,
    keep_only = NA_character_,
    review_only = NA_character_,
    keep_review_erased = NA_character_,
    keep_only_erased = NA_character_,
    review_only_erased = NA_character_,
    summary = NA_character_
  )

  if (isTRUE(save_outputs)) {
    if (is.null(out_dir) || !nzchar(out_dir)) stop("out_dir must be provided when save_outputs=TRUE")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    n_all   <- nrow(out)
    n_keep  <- sum(out$flag_internal == "keep", na.rm = TRUE)
    n_rev   <- sum(out$flag_internal == "review", na.rm = TRUE)
    n_drop_unburn <- sum(out$flag_internal == "drop_unburnable", na.rm = TRUE)
    n_drop_na     <- sum(out$flag_internal == "drop_corine_na", na.rm = TRUE)
    n_nohit <- sum(out$flag_internal == "nohit_stage1", na.rm = TRUE)

    if (driver == "GPKG") {
      gpkg_path <- file.path(out_dir, sprintf("%s_phase1.gpkg", prefix))
      if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

      sf::st_write(out, gpkg_path, layer = "stage2_flagged", quiet = quiet)
      sf::st_write(stage2_keep_review, gpkg_path, layer = "stage2_keep_review", append = TRUE, quiet = quiet)
      sf::st_write(stage2_dropped, gpkg_path, layer = "stage2_dropped", append = TRUE, quiet = quiet)

      if (isTRUE(save_splits)) {
        sf::st_write(stage2_keep_only, gpkg_path, layer = "stage2_keep_only", append = TRUE, quiet = quiet)
        sf::st_write(stage2_review_only, gpkg_path, layer = "stage2_review_only", append = TRUE, quiet = quiet)
      }

      if (!is.null(stage2_keep_review_erased)) {
        sf::st_write(stage2_keep_review_erased, gpkg_path, layer = "stage2_keep_review_erased", append = TRUE, quiet = quiet)
      }
      if (!is.null(stage2_keep_only_erased) && isTRUE(save_splits)) {
        sf::st_write(stage2_keep_only_erased, gpkg_path, layer = "stage2_keep_only_erased", append = TRUE, quiet = quiet)
      }
      if (!is.null(stage2_review_only_erased) && isTRUE(save_splits)) {
        sf::st_write(stage2_review_only_erased, gpkg_path, layer = "stage2_review_only_erased", append = TRUE, quiet = quiet)
      }

      out_paths$container   <- gpkg_path
      out_paths$flagged     <- paste0(gpkg_path, " | layer=stage2_flagged")
      out_paths$keep_review <- paste0(gpkg_path, " | layer=stage2_keep_review")
      out_paths$dropped     <- paste0(gpkg_path, " | layer=stage2_dropped")
      if (isTRUE(save_splits)) {
        out_paths$keep_only   <- paste0(gpkg_path, " | layer=stage2_keep_only")
        out_paths$review_only <- paste0(gpkg_path, " | layer=stage2_review_only")
      }
      if (!is.null(stage2_keep_review_erased)) out_paths$keep_review_erased <- paste0(gpkg_path, " | layer=stage2_keep_review_erased")
      if (!is.null(stage2_keep_only_erased) && isTRUE(save_splits)) out_paths$keep_only_erased <- paste0(gpkg_path, " | layer=stage2_keep_only_erased")
      if (!is.null(stage2_review_only_erased) && isTRUE(save_splits)) out_paths$review_only_erased <- paste0(gpkg_path, " | layer=stage2_review_only_erased")

    } else if (driver == "ESRI Shapefile") {

      p_flagged <- file.path(out_dir, sprintf("%s_stage2_flagged.shp", prefix))
      p_keeprev <- file.path(out_dir, sprintf("%s_stage2_keep_review.shp", prefix))
      p_drop    <- file.path(out_dir, sprintf("%s_stage2_dropped.shp", prefix))

      sf::st_write(out, p_flagged, delete_dsn = overwrite, quiet = quiet)
      sf::st_write(stage2_keep_review, p_keeprev, delete_dsn = overwrite, quiet = quiet)
      sf::st_write(stage2_dropped, p_drop, delete_dsn = overwrite, quiet = quiet)

      out_paths$container   <- out_dir
      out_paths$flagged     <- p_flagged
      out_paths$keep_review <- p_keeprev
      out_paths$dropped     <- p_drop

      if (isTRUE(save_splits)) {
        p_keep <- file.path(out_dir, sprintf("%s_stage2_keep_only.shp", prefix))
        p_rev  <- file.path(out_dir, sprintf("%s_stage2_review_only.shp", prefix))
        sf::st_write(stage2_keep_only, p_keep, delete_dsn = overwrite, quiet = quiet)
        sf::st_write(stage2_review_only, p_rev, delete_dsn = overwrite, quiet = quiet)
        out_paths$keep_only   <- p_keep
        out_paths$review_only <- p_rev
      }

      if (!is.null(stage2_keep_review_erased)) {
        p_er <- file.path(out_dir, sprintf("%s_stage2_keep_review_erased.shp", prefix))
        sf::st_write(stage2_keep_review_erased, p_er, delete_dsn = overwrite, quiet = quiet)
        out_paths$keep_review_erased <- p_er
      }
      if (!is.null(stage2_keep_only_erased) && isTRUE(save_splits)) {
        p_er <- file.path(out_dir, sprintf("%s_stage2_keep_only_erased.shp", prefix))
        sf::st_write(stage2_keep_only_erased, p_er, delete_dsn = overwrite, quiet = quiet)
        out_paths$keep_only_erased <- p_er
      }
      if (!is.null(stage2_review_only_erased) && isTRUE(save_splits)) {
        p_er <- file.path(out_dir, sprintf("%s_stage2_review_only_erased.shp", prefix))
        sf::st_write(stage2_review_only_erased, p_er, delete_dsn = overwrite, quiet = quiet)
        out_paths$review_only_erased <- p_er
      }

    } else {
      stop("driver must be 'GPKG' or 'ESRI Shapefile'.")
    }

    summary_path <- file.path(out_dir, sprintf("%s_phase1_summary.txt", prefix))
    summary_lines <- c(
      "PHASE 1 SUMMARY (stage2 flagged by stage1 support + burnable mask + optional NA diagnostic/drop)",
      sprintf("support_buffer_m: %s", support_buffer_m),
      sprintf("support_method: %s", support_method),
      sprintf("touches: %s", touches),
      sprintf("corine_na_scope: %s", corine_na_scope),
      sprintf("max_corine_na_frac: %s", ifelse(is.null(max_corine_na_frac), "NULL", max_corine_na_frac)),
      sprintf("na_exact (slow): %s", na_exact),
      sprintf("burnable_exact (slow): %s", burnable_exact),
      sprintf("clip_to_burnable: %s", clip_to_burnable),
      "",
      sprintf("N stage2 total: %s", n_all),
      sprintf("N keep: %s", n_keep),
      sprintf("N review: %s", n_rev),
      sprintf("N dropped (unburnable): %s", n_drop_unburn),
      sprintf("N dropped (corine NA): %s", n_drop_na),
      sprintf("N nohit_stage1 (pre-burnable split): %s", n_nohit),
      "",
      sprintf("Saved container: %s", out_paths$container)
    )
    writeLines(summary_lines, summary_path)
    out_paths$summary <- summary_path
  }


  # --------------------------------------------------
  # Convenience "clean" aliases (prefer erased outputs when present)
  # --------------------------------------------------
  stage2_keep_review_clean <- if (!is.null(stage2_keep_review_erased)) stage2_keep_review_erased else stage2_keep_review
  stage2_keep_only_clean   <- if (!is.null(stage2_keep_only_erased))   stage2_keep_only_erased   else stage2_keep_only
  stage2_review_only_clean <- if (!is.null(stage2_review_only_erased)) stage2_review_only_erased else stage2_review_only

  # Also provide a cleaned flagged layer for downstream validation/workflows.
  # This keeps dropped polygons unchanged and replaces keep/review geometries with the erased versions.
  stage2_flagged_clean <- out
  if (!is.null(stage2_keep_review_erased)) {
    dropped2 <- stage2_dropped
    keep2    <- stage2_keep_review_erased

    # Align columns to avoid rbind issues (e.g., median columns exist only in keep/review)
    miss_d <- setdiff(names(keep2), names(dropped2))
    if (length(miss_d) > 0) for (nm in miss_d) dropped2[[nm]] <- NA
    miss_k <- setdiff(names(dropped2), names(keep2))
    if (length(miss_k) > 0) for (nm in miss_k) keep2[[nm]] <- NA

    keep2 <- keep2[, names(dropped2), drop = FALSE]
    stage2_flagged_clean <- rbind(dropped2, keep2)
  }

  list(
    stage2_flagged     = out,
    stage2_flagged_clean = stage2_flagged_clean,
    stage2_keep_review = stage2_keep_review,
    stage2_keep_review_clean = stage2_keep_review_clean,
    stage2_dropped     = stage2_dropped,
    stage2_keep_only   = stage2_keep_only,
    stage2_keep_only_clean   = stage2_keep_only_clean,
    stage2_review_only = stage2_review_only,
    stage2_review_only_clean = stage2_review_only_clean,

    stage2_keep_review_erased = stage2_keep_review_erased,
    stage2_keep_only_erased   = stage2_keep_only_erased,
    stage2_review_only_erased = stage2_review_only_erased,

    out_paths = out_paths
  )
}

#' Score/prepare reference burned area polygons (e.g., EFFIS) with the same
#' per-polygon diagnostics used for internal polygons
#'
#' This helper is designed to make reference polygons directly comparable to the
#' internally detected polygons by producing:
#' - a burnable-domain quality flag (`flag_ref_internal`) with optional drops for
#'   high NA coverage in the burnable raster and/or non-burnable locations
#' - optional per-polygon median RBR (`median_colname`)
#' - optional pre-year erase (same geometry cleaning logic as internal)
#'
#' Unlike `scoring_internal_burned_area()`, this function does not use Stage-1 support
#' to define keep/review. Instead, it focuses on harmonizing per-polygon diagnostics
#' and optional domain filtering.
#'
#' @param ref_polys sf. Reference polygons (e.g., EFFIS).
#' @param burnable_corine sf or terra::SpatRaster. Burnable mask used to check domain
#'   and compute NA coverage when raster-based.
#' @param burnable_classes Integer vector or NULL. If `burnable_corine` is categorical,
#'   these values are treated as burnable. If NULL, `burnable_corine` must be binary-like (0/1).
#' @param rbr_rast Optional terra::SpatRaster. If provided and `compute_median=TRUE`,
#'   adds per-polygon median to `median_colname`.
#' @param max_corine_na_frac Numeric or NULL. If finite, polygons with NA fraction in
#'   `burnable_corine` greater than this value are flagged as `drop_corine_na`.
#' @param na_exact Logical. If TRUE, uses exact extraction weights when computing NA fraction.
#' @param burnable_exact Logical. If TRUE, uses exact extraction weights for burnable hit test.
#' @param clip_to_burnable Logical. If TRUE, clip polygons to burnable pixels before scoring.
#'   Requires `burnable_corine` as a raster.
#' @param clip_mask_buffer_m Numeric. Buffer (meters) for expanding burnable mask in raster clip.
#' @param clip_min_piece_area_m2 Numeric. Minimum area (m2) for clipped pieces to keep.
#' @param clip_keep_largest_piece Logical. If TRUE, keep only the largest clipped piece per polygon.
#' @param preyear_polys Optional sf. If provided, overlaps are erased from kept polygons.
#' @param erase_mask_buffer_m Numeric. Buffer on erase mask (meters).
#' @param erase_post_shave_m Numeric. Post-erase shave (meters).
#' @param erase_min_area_m2 Numeric. Minimum area (m2) to keep after erase.
#' @param compute_median Logical. If TRUE, compute median RBR for kept polygons.
#' @param median_colname Character. Column name for median RBR.
#' @param median_chunk_size Integer. Chunk size for median extraction.
#' @param save_outputs Logical. Save outputs to disk.
#' @param out_dir Character. Output directory when saving.
#' @param prefix Character. Prefix for filenames.
#' @param driver Character. "GPKG" or "ESRI Shapefile".
#' @param overwrite Logical.
#' @param quiet Logical.
#' @param save_splits Logical. If TRUE, save kept/dropped splits too.
#' @param touches Logical. Passed to terra::rasterize in clip-to-burnable.
#' @param suppress_sf_warnings Logical. Suppress common sf warnings.
#'
#' @return list with:
#' - ref_flagged: reference polygons with diagnostics and `flag_ref_internal`
#' - ref_flagged_clean: flagged layer where kept geometries are replaced by erased geometries (when available)
#' - ref_kept, ref_kept_clean, ref_dropped
#' - ref_kept_erased (NULL if no erase was performed)
#' - out_paths (paths when saved)
#' @export
scoring_reference_burned_area <- function(
  ref_polys,
  burnable_corine,
  burnable_classes = NULL,
  rbr_rast = NULL,

  max_corine_na_frac = 0.7,
  na_exact = FALSE,
  burnable_exact = FALSE,

  clip_to_burnable = FALSE,
  clip_mask_buffer_m = 0,
  clip_min_piece_area_m2 = 0,
  clip_keep_largest_piece = FALSE,

  preyear_polys = NULL,
  erase_mask_buffer_m = 0,
  erase_post_shave_m  = 0,
  erase_min_area_m2   = 0,

  compute_median = TRUE,
  median_colname = "median_rbr",
  median_chunk_size = 500L,

  save_outputs = FALSE,
  out_dir = NULL,
  prefix = "BA",
  driver = "GPKG",
  overwrite = TRUE,
  quiet = TRUE,
  save_splits = FALSE,

  touches = TRUE,
  suppress_sf_warnings = TRUE
){
  pkgs <- c("sf", "terra", "stats")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing required packages: ", paste(miss, collapse = ", "))

  stopifnot(inherits(ref_polys, "sf"))
  ref_polys <- .safe_make_valid_polys(ref_polys, suppress_sf_warnings = suppress_sf_warnings)
  if (is.null(ref_polys) || nrow(ref_polys) == 0) stop("ref_polys is empty after validation.")

  # Optional: clip reference polygons to burnable first (fast raster method)
  if (isTRUE(clip_to_burnable)) {
    if (!inherits(burnable_corine, "SpatRaster")) {
      stop("clip_to_burnable=TRUE requires burnable_corine as a terra::SpatRaster.")
    }

    # NOTE: duplicated from scoring_internal_burned_area() to keep behavior identical.
    clip_fast <- function(polys_sf, burnable_corine_rast, burnable_classes,
                          mask_buffer_m, min_piece_area_m2, keep_largest_piece, touches) {

      if (nrow(polys_sf) == 0) return(polys_sf)
      crs0 <- sf::st_crs(polys_sf)
      if (is.na(crs0)) stop("polys_sf CRS is NA.")

      polys_sf <- .safe_make_valid_polys(polys_sf, suppress_sf_warnings = suppress_sf_warnings)
      if (is.null(polys_sf) || nrow(polys_sf) == 0) return(polys_sf)

      if (!"._row_id" %in% names(polys_sf)) polys_sf$._row_id <- seq_len(nrow(polys_sf))

      r <- burnable_corine_rast[[1]]

      v <- terra::vect(polys_sf)
      if (!terra::same.crs(v, r)) v <- terra::project(v, r)

      e <- terra::ext(v)

      if (is.numeric(mask_buffer_m) && mask_buffer_m > 0) {
        res_xy <- terra::res(r)
        px <- min(res_xy, na.rm = TRUE)
        buf <- max(px, mask_buffer_m)
        e <- terra::ext(e[1] - buf, e[2] + buf, e[3] - buf, e[4] + buf)
      }

      r_sub <- terra::crop(r, e)
      vals <- terra::values(r_sub, mat = FALSE)
      if (all(is.na(vals))) return(polys_sf[0, , drop = FALSE])

      if (is.null(burnable_classes)) {
        if (!.is_binary_like(vals[!is.na(vals)])) {
          stop("burnable_classes is NULL but burnable_corine is not binary-like.")
        }
        burn <- terra::ifel(r_sub == 1, 1, NA)
      } else {
        burn <- terra::ifel(r_sub %in% burnable_classes, 1, NA)
      }

      if (all(is.na(terra::values(burn)))) return(polys_sf[0, , drop = FALSE])

      if (is.numeric(mask_buffer_m) && mask_buffer_m > 0) {
        res_xy <- terra::res(burn)
        px <- min(res_xy, na.rm = TRUE)
        npx <- max(1L, as.integer(ceiling(mask_buffer_m / px)))
        w <- matrix(1, nrow = 2 * npx + 1, ncol = 2 * npx + 1)
        burn <- terra::focal(burn, w = w, fun = "max", na.policy = "omit", fillvalue = NA)
        burn <- terra::ifel(burn >= 1, 1, NA)
      }

      id_r <- terra::rasterize(v, burn, field = "._row_id", touches = touches)
      keep_id <- terra::ifel(!is.na(burn) & !is.na(id_r), id_r, NA)
      if (all(is.na(terra::values(keep_id)))) return(polys_sf[0, , drop = FALSE])

      pv <- terra::as.polygons(keep_id, dissolve = TRUE, na.rm = TRUE)
      if (is.null(pv) || nrow(pv) == 0) return(polys_sf[0, , drop = FALSE])

      names(pv)[1] <- "._row_id"
      pv_sf <- sf::st_as_sf(pv)

      attrs <- sf::st_drop_geometry(polys_sf)
      attrs <- attrs[match(pv_sf$._row_id, attrs$._row_id), , drop = FALSE]
      pv_sf <- cbind(pv_sf, attrs[setdiff(names(attrs), "._row_id")])

      if (is.numeric(min_piece_area_m2) && min_piece_area_m2 > 0 && nrow(pv_sf) > 0) {
        a <- as.numeric(sf::st_area(pv_sf))
        pv_sf <- pv_sf[a >= min_piece_area_m2, , drop = FALSE]
      }
      if (nrow(pv_sf) == 0) return(polys_sf[0, , drop = FALSE])

      if (isTRUE(keep_largest_piece)) {
        pv_sf$._area <- as.numeric(sf::st_area(pv_sf))
        pv_sf <- pv_sf[order(pv_sf$._row_id, -pv_sf$._area), ]
        pv_sf <- pv_sf[!duplicated(pv_sf$._row_id), , drop = FALSE]
        pv_sf$._area <- NULL
      }

      if (sf::st_crs(pv_sf) != crs0) pv_sf <- sf::st_transform(pv_sf, crs0)
      pv_sf
    }

    ref_polys <- clip_fast(
      polys_sf = ref_polys,
      burnable_corine_rast = burnable_corine,
      burnable_classes = burnable_classes,
      mask_buffer_m = clip_mask_buffer_m,
      min_piece_area_m2 = clip_min_piece_area_m2,
      keep_largest_piece = clip_keep_largest_piece,
      touches = touches
    )

    ref_polys <- .safe_make_valid_polys(ref_polys, suppress_sf_warnings = suppress_sf_warnings)
    if (is.null(ref_polys) || nrow(ref_polys) == 0) stop("After clip_to_burnable, ref_polys is empty.")
  }

  out <- ref_polys
  out$flag_ref_internal <- "ref"
  out$why_flag_ref      <- "reference_input"
  out$corine_na_frac    <- NA_real_

  # --- burnable raster prep (if raster) ---
  burnable_r_na   <- NULL
  burnable_r_mask <- NULL
  if (inherits(burnable_corine, "SpatRaster")) {
    burnable_r_na <- burnable_corine[[1]]
    vals0 <- terra::values(burnable_r_na, mat = FALSE)
    vals0 <- vals0[!is.na(vals0)]

    if (is.null(burnable_classes)) {
      if (!.is_binary_like(vals0)) {
        stop("burnable_corine is categorical but burnable_classes is NULL. Provide burnable_classes or use a binary mask.")
      }
      burnable_r_mask <- terra::ifel(burnable_r_na == 1, 1, NA)
    } else {
      burnable_r_mask <- terra::ifel(burnable_r_na %in% burnable_classes, 1, NA)
    }
  }

  # --- corine_na_frac diagnostic (+ optional drop), raster mode only ---
  if (inherits(burnable_corine, "SpatRaster")) {
    idx_na <- seq_len(nrow(out))
    if (length(idx_na) > 0) {
      r_na <- burnable_r_na
      v_na <- terra::vect(out[idx_na, , drop = FALSE])
      if (!terra::same.crs(v_na, r_na)) v_na <- terra::project(v_na, terra::crs(r_na))

      if (!isTRUE(na_exact)) {
        na_fun <- function(x, ...) if (length(x) == 0) 1 else mean(is.na(x))
        ex <- terra::extract(r_na, v_na, fun = na_fun, exact = FALSE)
        ex <- as.data.frame(ex)

        na_frac <- rep(1, length(idx_na))
        if (nrow(ex) > 0) {
          ids  <- as.integer(ex[[1]])
          vals <- as.numeric(ex[[2]])
          vals[!is.finite(vals)] <- 1
          na_frac[ids] <- vals
        }
      } else {
        ex <- terra::extract(r_na, v_na, exact = TRUE)
        ex <- as.data.frame(ex)
        if (!"ID" %in% names(ex)) stop("terra::extract output does not contain 'ID' column.")

        cov_candidates <- c("coverage_fraction", "fraction", "weight", "weights")
        cov_col <- cov_candidates[cov_candidates %in% names(ex)][1]
        if (is.na(cov_col)) cov_col <- NULL

        value_cols <- setdiff(names(ex), c("ID", cov_candidates))
        if (length(value_cols) < 1) stop("Could not identify value column from terra::extract output.")
        val_col <- value_cols[1]

        cov <- if (!is.null(cov_col)) as.numeric(ex[[cov_col]]) else rep(1, nrow(ex))
        val <- ex[[val_col]]
        ids <- as.integer(ex[["ID"]])

        tmp <- tapply(seq_len(nrow(ex)), ids, function(ii) {
          cc <- cov[ii]
          vv <- val[ii]
          tot <- sum(cc, na.rm = TRUE)
          if (!is.finite(tot) || tot <= 0) return(1)
          na_cov <- sum(cc[is.na(vv)], na.rm = TRUE)
          as.numeric(na_cov / tot)
        })

        na_frac <- rep(1, length(idx_na))
        na_frac[as.integer(names(tmp))] <- as.numeric(tmp)
        na_frac[!is.finite(na_frac)] <- 1
      }

      out$corine_na_frac[idx_na] <- na_frac

      if (!is.null(max_corine_na_frac) && is.finite(max_corine_na_frac)) {
        drop_na <- na_frac > max_corine_na_frac
        if (any(drop_na)) {
          out$flag_ref_internal[idx_na[drop_na]] <- "drop_corine_na"
          out$why_flag_ref[idx_na[drop_na]] <- sprintf("corine_na_frac_gt_%.2f", max_corine_na_frac)
        }
      }
    }
  }

  # --- burnable mask check (all non-dropped) ---
  idx <- which(out$flag_ref_internal == "ref")
  if (length(idx) > 0) {

    if (inherits(burnable_corine, "sf")) {

      burn_sf <- .safe_make_valid_polys(burnable_corine, suppress_sf_warnings = suppress_sf_warnings)
      if (is.null(burn_sf) || nrow(burn_sf) == 0) {
        out$flag_ref_internal[idx] <- "drop_unburnable"
        out$why_flag_ref[idx]      <- "burnable_corine_empty"
      } else {
        if (sf::st_crs(out) != sf::st_crs(burn_sf)) burn_sf <- sf::st_transform(burn_sf, sf::st_crs(out))
        hit_burnable <- lengths(sf::st_intersects(out[idx, , drop = FALSE], burn_sf)) > 0
        out$flag_ref_internal[idx] <- ifelse(hit_burnable, "ref", "drop_unburnable")
        out$why_flag_ref[idx]      <- ifelse(hit_burnable, "burnable_corine_ok", "outside_burnable_corine")
      }

    } else if (inherits(burnable_corine, "SpatRaster")) {

      r_burn <- burnable_r_mask
      v_burn <- terra::vect(out[idx, , drop = FALSE])
      if (!terra::same.crs(v_burn, r_burn)) v_burn <- terra::project(v_burn, terra::crs(r_burn))

      if (!isTRUE(burnable_exact)) {
        burn_fun <- function(x, ...) as.integer(length(x) > 0 && any(x == 1, na.rm = TRUE))
        ex <- terra::extract(r_burn, v_burn, fun = burn_fun, exact = FALSE)
        ex <- as.data.frame(ex)

        hit_burnable <- rep(FALSE, length(idx))
        if (nrow(ex) > 0) {
          ids  <- as.integer(ex[[1]])
          vals <- as.integer(ex[[2]])
          hit_burnable[ids] <- as.logical(vals)
          hit_burnable[is.na(hit_burnable)] <- FALSE
        }
      } else {
        ex <- terra::extract(r_burn, v_burn, exact = TRUE)
        ex <- as.data.frame(ex)
        if (!"ID" %in% names(ex)) stop("terra::extract output does not contain 'ID' column.")

        value_cols <- setdiff(names(ex), c("ID", "coverage_fraction", "fraction", "weight", "weights"))
        if (length(value_cols) < 1) stop("Could not identify value column from terra::extract output.")
        val <- suppressWarnings(as.integer(ex[[value_cols[1]]]))
        ids <- as.integer(ex[["ID"]])

        tmp <- tapply(seq_len(nrow(ex)), ids, function(ii) any(val[ii] == 1, na.rm = TRUE))
        hit_burnable <- rep(FALSE, length(idx))
        hit_burnable[as.integer(names(tmp))] <- as.logical(tmp)
        hit_burnable[is.na(hit_burnable)] <- FALSE
      }

      out$flag_ref_internal[idx] <- ifelse(hit_burnable, "ref", "drop_unburnable")
      out$why_flag_ref[idx]      <- ifelse(hit_burnable, "burnable_corine_ok", "outside_burnable_corine")

    } else {
      stop("burnable_corine must be an sf or a terra::SpatRaster.")
    }
  }

  # Split outputs
  ref_kept    <- out[out$flag_ref_internal == "ref", , drop = FALSE]
  ref_dropped <- out[out$flag_ref_internal %in% c("drop_unburnable", "drop_corine_na"), , drop = FALSE]

  # Optional: compute median RBR for kept polygons
  if (isTRUE(compute_median) && !is.null(rbr_rast) && nrow(ref_kept) > 0) {
    if (!inherits(rbr_rast, "SpatRaster")) stop("rbr_rast must be a terra::SpatRaster.")
    if (terra::nlyr(rbr_rast) > 1) rbr_rast <- rbr_rast[[1]]
    ref_kept <- .add_polygon_median_rbr(ref_kept, rbr_rast, chunk_size = median_chunk_size, colname = median_colname)
  }

  # Optional: erase overlaps with preyear (kept polygons only)
  ref_kept_erased <- NULL
  if (!is.null(preyear_polys) && nrow(ref_kept) > 0) {
    preyear_polys <- .safe_make_valid_polys(preyear_polys, suppress_sf_warnings = suppress_sf_warnings)
    if (!is.null(preyear_polys) && nrow(preyear_polys) > 0) {
      ref_kept_erased <- erase_overlap_with_preyear(
        stage2_sf = ref_kept,
        preyear_sf = preyear_polys,
        mask_buffer_m = erase_mask_buffer_m,
        post_shave_m  = erase_post_shave_m,
        min_area_m2   = erase_min_area_m2,
        quiet = quiet,
        suppress_sf_warnings = suppress_sf_warnings
      )
    }
  }

  # Clean aliases
  ref_kept_clean <- if (!is.null(ref_kept_erased)) ref_kept_erased else ref_kept

  # Build flagged_clean (replace kept geometries by erased when available)
  ref_flagged_clean <- out
  if (!is.null(ref_kept_erased)) {
    dropped2 <- ref_dropped
    kept2    <- ref_kept_erased

    miss_d <- setdiff(names(kept2), names(dropped2))
    if (length(miss_d) > 0) for (nm in miss_d) dropped2[[nm]] <- NA
    miss_k <- setdiff(names(dropped2), names(kept2))
    if (length(miss_k) > 0) for (nm in miss_k) kept2[[nm]] <- NA

    kept2 <- kept2[, names(dropped2), drop = FALSE]
    ref_flagged_clean <- rbind(dropped2, kept2)
  }


  out_paths <- list(
    container = NA_character_,
    flagged = NA_character_,
    kept = NA_character_,
    dropped = NA_character_,
    kept_erased = NA_character_,
    summary = NA_character_
  )

  if (isTRUE(save_outputs)) {
    if (is.null(out_dir) || !nzchar(out_dir)) stop("out_dir must be provided when save_outputs=TRUE")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (driver == "GPKG") {
      gpkg_path <- file.path(out_dir, sprintf("%s_reference.gpkg", prefix))
      if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

      sf::st_write(out, gpkg_path, layer = "ref_flagged", quiet = quiet)
      sf::st_write(ref_flagged_clean, gpkg_path, layer = "ref_flagged_clean", append = TRUE, quiet = quiet)
      sf::st_write(ref_kept, gpkg_path, layer = "ref_kept", append = TRUE, quiet = quiet)
      sf::st_write(ref_dropped, gpkg_path, layer = "ref_dropped", append = TRUE, quiet = quiet)

      if (!is.null(ref_kept_erased)) {
        sf::st_write(ref_kept_erased, gpkg_path, layer = "ref_kept_erased", append = TRUE, quiet = quiet)
      }

      out_paths$container <- gpkg_path
      out_paths$flagged   <- paste0(gpkg_path, " | layer=ref_flagged")
      out_paths$kept      <- paste0(gpkg_path, " | layer=ref_kept")
      out_paths$dropped   <- paste0(gpkg_path, " | layer=ref_dropped")
      if (!is.null(ref_kept_erased)) out_paths$kept_erased <- paste0(gpkg_path, " | layer=ref_kept_erased")

    } else if (driver == "ESRI Shapefile") {

      p_flagged <- file.path(out_dir, sprintf("%s_ref_flagged.shp", prefix))
      p_kept    <- file.path(out_dir, sprintf("%s_ref_kept.shp", prefix))
      p_drop    <- file.path(out_dir, sprintf("%s_ref_dropped.shp", prefix))

      sf::st_write(out, p_flagged, delete_dsn = overwrite, quiet = quiet)
      if (isTRUE(save_splits)) {
        sf::st_write(ref_kept, p_kept, delete_dsn = overwrite, quiet = quiet)
        sf::st_write(ref_dropped, p_drop, delete_dsn = overwrite, quiet = quiet)
      }

      if (!is.null(ref_kept_erased) && isTRUE(save_splits)) {
        p_kept_e <- file.path(out_dir, sprintf("%s_ref_kept_erased.shp", prefix))
        sf::st_write(ref_kept_erased, p_kept_e, delete_dsn = overwrite, quiet = quiet)
        out_paths$kept_erased <- p_kept_e
      }

      out_paths$container <- out_dir
      out_paths$flagged   <- p_flagged
      if (isTRUE(save_splits)) {
        out_paths$kept    <- p_kept
        out_paths$dropped <- p_drop
      }

    } else {
      stop("driver must be 'GPKG' or 'ESRI Shapefile'.")
    }

    summary_path <- file.path(out_dir, sprintf("%s_reference_summary.txt", prefix))
    summary_lines <- c(
      "REFERENCE SCORING SUMMARY",
      sprintf("N ref total: %s", nrow(out)),
      sprintf("N ref kept: %s", nrow(ref_kept)),
      sprintf("N ref dropped_unburnable: %s", sum(out$flag_ref_internal == "drop_unburnable", na.rm = TRUE)),
      sprintf("N ref dropped_corine_na: %s", sum(out$flag_ref_internal == "drop_corine_na", na.rm = TRUE)),
      sprintf("compute_median: %s", compute_median),
      sprintf("median_colname: %s", median_colname),
      sprintf("Saved container: %s", out_paths$container)
    )
    writeLines(summary_lines, summary_path)
    out_paths$summary <- summary_path
  }

  list(
    ref_flagged = out,
    ref_flagged_clean = ref_flagged_clean,
    ref_kept = ref_kept,
    ref_kept_clean = ref_kept_clean,
    ref_dropped = ref_dropped,
    ref_kept_erased = ref_kept_erased,
    out_paths = out_paths
  )
}
# ============================
# Validation helpers (internal)
# ============================

#' @keywords internal
.flag_reference_against_internal_keep <- function(
  ref_polys,
  internal_polys,
  internal_keep_value = "keep",
  buffer_m = 0,
  suppress_sf_warnings = TRUE
){
  stopifnot(inherits(ref_polys, "sf"), inherits(internal_polys, "sf"))
  if (!"flag_internal" %in% names(internal_polys)) {
    stop("internal_polys must contain column 'flag_internal'.")
  }

  ref_polys      <- .safe_make_valid_polys(ref_polys, suppress_sf_warnings = suppress_sf_warnings)
  internal_polys <- .safe_make_valid_polys(internal_polys, suppress_sf_warnings = suppress_sf_warnings)

  internal_keep <- internal_polys[internal_polys$flag_internal %in% internal_keep_value, , drop = FALSE]
  if (nrow(internal_keep) == 0) stop("No internal keep polygons found.")

  if (sf::st_crs(ref_polys) != sf::st_crs(internal_keep)) {
    ref_polys <- sf::st_transform(ref_polys, sf::st_crs(internal_keep))
  }

  internal_keep_use <- internal_keep
  if (is.numeric(buffer_m) && buffer_m > 0) {
    internal_keep_use <- .sf_quiet(sf::st_buffer(internal_keep_use, dist = buffer_m), enable = suppress_sf_warnings)
    internal_keep_use <- .safe_make_valid_polys(internal_keep_use, suppress_sf_warnings = suppress_sf_warnings)
  }

  hit_internal_keep <- lengths(sf::st_intersects(ref_polys, internal_keep_use)) > 0

  out <- ref_polys
  out$flag_ref <- ifelse(hit_internal_keep, "keep_ref", "review_ref")
  out$why_ref  <- ifelse(hit_internal_keep, "intersects_internal_keep", "no_intersection_internal_keep")

  list(
    ref_flagged = out,
    keep_ref_only   = out[out$flag_ref == "keep_ref", , drop = FALSE],
    review_ref_only = out[out$flag_ref == "review_ref", , drop = FALSE]
  )
}

#' @keywords internal
.flag_internal_keep_against_reference <- function(
  internal_polys,
  ref_polys,
  internal_keep_value = "keep",
  buffer_m = 0,
  suppress_sf_warnings = TRUE
){
  stopifnot(inherits(internal_polys, "sf"), inherits(ref_polys, "sf"))
  if (!"flag_internal" %in% names(internal_polys)) stop("flag_internal not found in internal_polys.")

  internal_polys <- .safe_make_valid_polys(internal_polys, suppress_sf_warnings = suppress_sf_warnings)
  ref_polys      <- .safe_make_valid_polys(ref_polys, suppress_sf_warnings = suppress_sf_warnings)

  if (sf::st_crs(ref_polys) != sf::st_crs(internal_polys)) {
    ref_polys <- sf::st_transform(ref_polys, sf::st_crs(internal_polys))
  }

  ref_use <- ref_polys
  if (is.numeric(buffer_m) && buffer_m > 0) {
    ref_use <- .sf_quiet(sf::st_buffer(ref_use, dist = buffer_m), enable = suppress_sf_warnings)
    ref_use <- .safe_make_valid_polys(ref_use, suppress_sf_warnings = suppress_sf_warnings)
  }

  is_keep <- internal_polys$flag_internal %in% internal_keep_value
  hit_ref <- rep(NA, nrow(internal_polys))
  if (any(is_keep)) {
    hit_ref[is_keep] <- lengths(sf::st_intersects(internal_polys[is_keep, , drop = FALSE], ref_use)) > 0
  }

  out <- internal_polys
  out$hit_ref_any <- hit_ref

  out$flag_internal_effis <- as.character(out$flag_internal)
  out$flag_internal_effis[is_keep & !is.na(hit_ref) & hit_ref]  <- "keep_common"
  out$flag_internal_effis[is_keep & !is.na(hit_ref) & !hit_ref] <- "review_int_no_effis"

  list(internal_flagged = out)
}

#' Score validation between internal polygons and reference polygons
#'
#' Produces two complementary products:
#' 1) reference polygons flagged as keep_ref/review_ref against internal keep
#' 2) internal keep polygons flagged as keep_common/review_int_no_effis against reference
#'
#' @param internal_polys sf. Internal polygons (must contain column `flag_internal`).
#' @param ref_polys sf. Reference polygons (e.g., EFFIS).
#' @param internal_keep_value Character vector. Values in `flag_internal` treated as "keep".
#' @param buffer_m numeric. Tolerance buffer (meters) used for both directions.
#' @param save_outputs logical. Save outputs to disk.
#' @param out_dir character. Output directory when saving.
#' @param prefix character. Prefix for filenames.
#' @param driver character. "GPKG" or "ESRI Shapefile".
#' @param overwrite logical.
#' @param quiet logical.
#' @param save_splits logical. If TRUE, save keep/ref splits too.
#' @param suppress_sf_warnings logical. Suppress sf warnings.
#'
#' @return list with `ref_flagged`, `internal_flagged`, and `out_paths` (if saved).
#' @export
scoring_validation_burned_area <- function(
  internal_polys,
  ref_polys,
  internal_keep_value = "keep",
  buffer_m = 0,

  save_outputs = FALSE,
  out_dir = NULL,
  prefix = "BA",
  driver = "GPKG",
  overwrite = TRUE,
  quiet = TRUE,
  save_splits = FALSE,

  suppress_sf_warnings = TRUE
){
  stopifnot(inherits(internal_polys, "sf"), inherits(ref_polys, "sf"))
  if (!"flag_internal" %in% names(internal_polys)) {
    stop("internal_polys must contain column 'flag_internal'.")
  }

  tmp_ref <- .flag_reference_against_internal_keep(
    ref_polys = ref_polys,
    internal_polys = internal_polys,
    internal_keep_value = internal_keep_value,
    buffer_m = buffer_m,
    suppress_sf_warnings = suppress_sf_warnings
  )

  tmp_int <- .flag_internal_keep_against_reference(
    internal_polys = internal_polys,
    ref_polys = ref_polys,
    internal_keep_value = internal_keep_value,
    buffer_m = buffer_m,
    suppress_sf_warnings = suppress_sf_warnings
  )

  out_paths <- list(container = NA_character_, summary = NA_character_)

  if (isTRUE(save_outputs)) {
    if (is.null(out_dir) || !nzchar(out_dir)) stop("out_dir must be provided when save_outputs=TRUE")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (driver == "GPKG") {
      gpkg_path <- file.path(out_dir, sprintf("%s_validation.gpkg", prefix))
      if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

      sf::st_write(tmp_ref$ref_flagged, gpkg_path, layer = "ref_flagged", quiet = quiet)
      sf::st_write(tmp_int$internal_flagged, gpkg_path, layer = "internal_flagged", append = TRUE, quiet = quiet)

      if (isTRUE(save_splits)) {
        sf::st_write(tmp_ref$keep_ref_only, gpkg_path, layer = "ref_keep_ref", append = TRUE, quiet = quiet)
        sf::st_write(tmp_ref$review_ref_only, gpkg_path, layer = "ref_review_ref", append = TRUE, quiet = quiet)
      }

      out_paths$container <- gpkg_path

    } else if (driver == "ESRI Shapefile") {

      p_ref <- file.path(out_dir, sprintf("%s_ref_flagged.shp", prefix))
      p_int <- file.path(out_dir, sprintf("%s_internal_flagged.shp", prefix))
      sf::st_write(tmp_ref$ref_flagged, p_ref, delete_dsn = overwrite, quiet = quiet)
      sf::st_write(tmp_int$internal_flagged, p_int, delete_dsn = overwrite, quiet = quiet)

      out_paths$container <- out_dir

    } else {
      stop("driver must be 'GPKG' or 'ESRI Shapefile'.")
    }

    summary_path <- file.path(out_dir, sprintf("%s_validation_summary.txt", prefix))
    summary_lines <- c(
      "VALIDATION SUMMARY (internal vs reference)",
      sprintf("internal_keep_value: %s", paste(internal_keep_value, collapse = ",")),
      sprintf("buffer_m: %s", buffer_m),
      "",
      sprintf("N ref total: %s", nrow(tmp_ref$ref_flagged)),
      sprintf("N ref keep_ref: %s", sum(tmp_ref$ref_flagged$flag_ref == "keep_ref", na.rm = TRUE)),
      sprintf("N ref review_ref: %s", sum(tmp_ref$ref_flagged$flag_ref == "review_ref", na.rm = TRUE)),
      "",
      sprintf("Saved container: %s", out_paths$container)
    )
    writeLines(summary_lines, summary_path)
    out_paths$summary <- summary_path
  }

  list(
    ref_flagged = tmp_ref$ref_flagged,
    keep_ref_only = tmp_ref$keep_ref_only,
    review_ref_only = tmp_ref$review_ref_only,
    internal_flagged = tmp_int$internal_flagged,
    out_paths = out_paths
  )
}

#' Stage-2 burned area scoring with optional preyear erase and optional validation
#'
#' Convenience wrapper that:
#' 1) runs `scoring_internal_burned_area()` (Stage-2 flagged using Stage-1 support + burnable mask),
#' 2) optionally erases overlaps with a pre-year mask (already included in internal scoring when `preyear_polys` is set),
#' 3) optionally runs `scoring_validation_burned_area()` if `ref_polys` is provided.
#'
#' @inheritParams scoring_internal_burned_area
#' @param ref_polys Optional sf. Reference polygons (e.g., EFFIS). If provided, validation is run.
#' @param ref_buffer_m Numeric. Tolerance buffer (meters) used for both directions in validation.
#' @param internal_keep_value Character vector. Internal keep label(s) used in validation (default "keep").
#' @param validate_use_clean Logical. If TRUE (default), validation is performed on the final keep/review geometry
#' delivered for analysis by using `stage2_keep_review_clean` (i.e., `stage2_keep_review_erased` when available).
#' If FALSE, validation uses the raw keep/review layer (`stage2_keep_review`).

#' @param score_ref_polys Logical. If TRUE (default) and `ref_polys` is provided, the wrapper
#' also runs `scoring_reference_burned_area()` to harmonize reference outputs (burnable/NA diagnostics,
#' optional median RBR, and optional preyear erase).
#' @param validate_ref_use_clean Logical. Only used when `score_ref_polys=TRUE`. If TRUE (default),
#' validation uses the cleaned reference geometry (`ref_kept_clean` or `ref_flagged_clean`) when available.
#' @param validate_ref_use_kept Logical. Only used when `score_ref_polys=TRUE`. If TRUE (default),
#' validation uses the kept reference layer (`ref_kept*`) rather than the full flagged reference layer.

#'
#' @return A list with:
#' - internal: output of `scoring_internal_burned_area()`
#' - reference: output of `scoring_reference_burned_area()` (or NULL if `ref_polys` is NULL or `score_ref_polys=FALSE`)
#' - validation: output of `scoring_validation_burned_area()` (or NULL if `ref_polys` is NULL)
#'
#' @export
scoring_burned_area_stage2 <- function(
  polys_stage2,
  polys_stage1,
  burnable_corine,
  burnable_classes = NULL,

  rbr_rast = NULL,

  support_buffer_m = 0,

  max_corine_na_frac = 0.7,
  corine_na_scope = c("nohit_stage1", "keep", "keep_nohit", "all"),

  support_method = c("raster", "vector"),
  touches = TRUE,
  make_valid = c("auto", "always", "none"),
  na_exact = FALSE,
  burnable_exact = FALSE,

  clip_to_burnable = FALSE,
  clip_mask_buffer_m = 0,
  clip_min_piece_area_m2 = 0,
  clip_keep_largest_piece = FALSE,

  preyear_polys = NULL,
  erase_mask_buffer_m = 0,
  erase_post_shave_m  = 0,
  erase_min_area_m2   = 0,

  compute_median = TRUE,
  median_colname = "median_rbr",
  median_chunk_size = 500L,

  ref_polys = NULL,
  ref_buffer_m = 0,
  internal_keep_value = "keep",
  validate_use_clean = TRUE,

  # Reference scoring (for EFFIS or other ref polygons)
  score_ref_polys = TRUE,
  validate_ref_use_clean = TRUE,
  validate_ref_use_kept  = TRUE,

  save_outputs = FALSE,
  out_dir = NULL,
  prefix = "BA",
  driver = "GPKG",
  overwrite = TRUE,
  quiet = TRUE,
  save_splits = FALSE,

  suppress_sf_warnings = TRUE
){
  internal <- scoring_internal_burned_area(
    polys_stage2 = polys_stage2,
    polys_stage1 = polys_stage1,
    burnable_corine = burnable_corine,
    burnable_classes = burnable_classes,
    rbr_rast = rbr_rast,
    support_buffer_m = support_buffer_m,
    max_corine_na_frac = max_corine_na_frac,
    corine_na_scope = corine_na_scope,
    support_method = support_method,
    touches = touches,
    make_valid = make_valid,
    na_exact = na_exact,
    burnable_exact = burnable_exact,
    clip_to_burnable = clip_to_burnable,
    clip_mask_buffer_m = clip_mask_buffer_m,
    clip_min_piece_area_m2 = clip_min_piece_area_m2,
    clip_keep_largest_piece = clip_keep_largest_piece,
    preyear_polys = preyear_polys,
    erase_mask_buffer_m = erase_mask_buffer_m,
    erase_post_shave_m = erase_post_shave_m,
    erase_min_area_m2 = erase_min_area_m2,
    compute_median = compute_median,
    median_colname = median_colname,
    median_chunk_size = median_chunk_size,
    save_outputs = save_outputs,
    out_dir = out_dir,
    prefix = prefix,
    driver = driver,
    overwrite = overwrite,
    quiet = quiet,
    save_splits = save_splits,
    suppress_sf_warnings = suppress_sf_warnings
  )

  reference <- NULL
  validation <- NULL

  if (!is.null(ref_polys)) {

    if (isTRUE(score_ref_polys)) {
      reference <- scoring_reference_burned_area(
        ref_polys = ref_polys,
        burnable_corine = burnable_corine,
        burnable_classes = burnable_classes,
        rbr_rast = rbr_rast,

        max_corine_na_frac = max_corine_na_frac,
        na_exact = na_exact,
        burnable_exact = burnable_exact,

        clip_to_burnable = clip_to_burnable,
        clip_mask_buffer_m = clip_mask_buffer_m,
        clip_min_piece_area_m2 = clip_min_piece_area_m2,
        clip_keep_largest_piece = clip_keep_largest_piece,

        preyear_polys = preyear_polys,
        erase_mask_buffer_m = erase_mask_buffer_m,
        erase_post_shave_m  = erase_post_shave_m,
        erase_min_area_m2   = erase_min_area_m2,

        compute_median = compute_median,
        median_colname = median_colname,
        median_chunk_size = median_chunk_size,

        save_outputs = save_outputs,
        out_dir = out_dir,
        prefix = prefix,
        driver = driver,
        overwrite = overwrite,
        quiet = quiet,
        save_splits = save_splits,

        touches = touches,
        suppress_sf_warnings = suppress_sf_warnings
      )
    }

    ref_for_validation <- ref_polys
    if (isTRUE(score_ref_polys) && !is.null(reference)) {
      if (isTRUE(validate_ref_use_kept)) {
        ref_for_validation <- if (isTRUE(validate_ref_use_clean)) reference$ref_kept_clean else reference$ref_kept
      } else {
        ref_for_validation <- if (isTRUE(validate_ref_use_clean)) reference$ref_flagged_clean else reference$ref_flagged
      }
    }

    validation <- scoring_validation_burned_area(
      internal_polys = if (isTRUE(validate_use_clean)) internal$stage2_keep_review_clean else internal$stage2_keep_review,
      ref_polys = ref_for_validation,
      internal_keep_value = internal_keep_value,
      buffer_m = ref_buffer_m,
      save_outputs = save_outputs,
      out_dir = out_dir,
      prefix = prefix,
      driver = driver,
      overwrite = overwrite,
      quiet = quiet,
      save_splits = save_splits,
      suppress_sf_warnings = suppress_sf_warnings
    )
  }

  list(internal = internal, reference = reference, validation = validation)
}
