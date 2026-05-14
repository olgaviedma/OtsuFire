#' Segmentation refinement around detected polygons using AOI bounding boxes
#'
#' @description
#' Builds axis-aligned AOIs (bounding boxes) around selected detected polygons,
#' optionally merges overlapping AOIs, and runs a user-supplied processing
#' function (`process_fun`) independently for each AOI (typically in parallel).
#' Each AOI run is restartable via per-folder DONE/RUNNING/ERROR markers.
#'
#' @details
#' This helper supports “targeted rescue/refinement” workflows: instead of
#' re-running a full-scene segmentation, it re-processes only neighborhoods that
#' are likely to contain errors (e.g., around large detections or around omission
#' areas inferred from reference perimeters). It can operate with or without
#' external reference polygons.
#'
#' The workflow is:
#' \enumerate{
#'   \item \strong{Load detections}. Provide `detected_polys` (sf or path), or
#'   derive detections from `binary_final_raster` (0/1).
#'   \item \strong{CRS harmonization}. Detected (and reference, if provided)
#'   polygons are transformed to `target_crs` (must be projected; meters).
#'   \item \strong{Select detections to refine}.
#'   \itemize{
#'     \item If `reference_polys` is \code{NULL}: select the \code{top_n_if_no_ref}
#'     detections by area (after filtering by \code{min_detected_area_m2}).
#'     \item If `reference_polys` is provided: compute an omission layer as
#'     \code{reference - union(detected)}, optionally buffer it by
#'     \code{omission_buffer_m}, then select detections intersecting the omission.
#'     If omission is empty or no intersections occur, fall back to
#'     \code{top_n_if_no_ref} by area.
#'   }
#'   \item \strong{Build AOIs}. For each selected detection, build a bbox expanded
#'   by \code{aoi_buffer_m}. If \code{clip_aois_to_raster_extent=TRUE}, AOIs are
#'   numerically clamped to the raster bbox (in \code{target_crs}).
#'   \item \strong{Merge overlaps (optional)}. If \code{merge_overlaps=TRUE}, AOIs
#'   intersecting (optionally after a small \code{merge_overlaps_buffer_m} applied
#'   only for grouping) are dissolved to reduce the number of AOI tasks.
#'   \item \strong{Per-AOI execution}. For each AOI:
#'   crop/mask \code{raster_path} in source CRS, project the crop to
#'   \code{target_crs} (optionally at \code{rescue_args$resolution}), align an
#'   optional CORINE raster to the projected template, clip optional ecoregions to
#'   the AOI, then call \code{process_fun(...)}. Progress is reported via
#'   \pkg{progressr}.
#' }
#'
#' \strong{Restartability and markers}
#' \itemize{
#'   \item Each AOI uses a dedicated folder: \code{base_output_dir/AOI_###/}.
#'   \item \code{RUNNING.ok} is created at start; removed on success or error.
#'   \item \code{DONE.ok} is created on success; completed AOIs are skipped on re-runs.
#'   \item \code{ERROR.txt} stores the error message when an AOI fails.
#' }
#'
#' \strong{How \code{process_fun} is called}
#' \code{process_fun} must be a function. It is called per AOI as:
#' \preformatted{
#' process_fun(
#'   raster_path = <projected AOI raster>,
#'   output_dir  = <AOI folder>,
#'   year        = <year>,
#'   corine_raster_path       = <aligned CORINE raster or NULL>,
#'   reclassify_corine        = FALSE,
#'   corine_classes           = <corine_classes>,
#'   ecoregion_shapefile_path = <clipped ecoregions file or NULL>,
#'   ecoregion_field          = <ecoregion_field>,
#'   ... additional args from rescue_args ...
#' )
#' }
#'
#' \strong{Argument compatibility and “unused argument” failures}
#' If \code{process_fun} does not include \code{...} in its signature, arguments
#' that are not part of \code{formals(process_fun)} can trigger an \dQuote{unused
#' arguments} error and prevent AOI outputs. This helper can pre-check and/or
#' drop such arguments (depending on \code{stop_on_unused_args}).
#'
#' \strong{Ecoregion attribute field}
#' Provide \code{ecoregion_field} as a valid attribute name in the ecoregion
#' layer (recommended for reproducibility). If \code{ecoregion_field} is
#' \code{NULL}, the implementation may attempt a backward-compatible detection
#' (e.g., \code{"EnS_name"} then \code{"EnZ_name"}), otherwise it errors.
#'
#' \strong{Resolution control}
#' If \code{rescue_args$resolution} is provided (single numeric), it is used as
#' the target projected resolution when projecting the AOI raster.
#'
#' @param n_workers Integer. Number of parallel workers (capped to available cores minus one).
#' @param raster_path Character. Path to the source raster to be refined (e.g., RBR).
#' @param base_output_dir Character. Output directory where AOI subfolders are created.
#' @param year Integer or character. Passed through to \code{process_fun} for naming/logging.
#' @param detected_polys sf object or character path. Detected polygons used to build AOIs.
#' If \code{NULL}, \code{binary_final_raster} must be provided.
#' @param binary_final_raster Character path. Binary raster (0/1) polygonized to derive detections
#' when \code{detected_polys} is \code{NULL}.
#' @param reference_polys sf object or character path. Optional reference perimeters used only
#' to prioritize detections to refine (omission-based selection).
#' @param target_crs Projected CRS in meters. Accepts numeric EPSG (e.g., 3035), \code{"EPSG:3035"},
#' WKT, or \code{sf::st_crs()}. Geographic CRS (lon/lat) are rejected.
#' @param aoi_buffer_m Numeric. Buffer (meters) applied to each detection bbox.
#' @param omission_buffer_m Numeric. Buffer (meters) applied to omission polygons before intersection
#' with detections (only used when \code{reference_polys} is provided).
#' @param top_n_if_no_ref Numeric (integer-like), allowing \code{Inf}. Fallback when \code{reference_polys}
#' is \code{NULL}/empty (or omission cannot be derived): select the top \code{N} detections by area.
#' Use this to cap runtime (e.g., 50). Use \code{Inf} to refine all eligible detections.
#' @param min_detected_area_m2 Numeric. Minimum detected polygon area (m^2) to be eligible for refinement.
#' @param clip_aois_to_raster_extent Logical. If \code{TRUE}, clamp AOIs to the raster bbox in \code{target_crs}.
#' @param merge_overlaps Logical. If \code{TRUE}, dissolve overlapping AOIs to reduce tasks.
#' @param merge_overlaps_buffer_m Numeric. Optional buffering (meters) used only for grouping before merging AOIs.
#' @param verbose Logical. If \code{TRUE}, print progress messages during AOI construction/merging.
#' @param corine_raster_path Character. Optional CORINE raster path; if provided, it is cropped per AOI and aligned
#' to the projected AOI raster template.
#' @param ecoregion_shapefile_path Character. Optional ecoregion polygon layer path; if provided, it is clipped per AOI
#' and written inside the AOI folder.
#' @param ecoregion_field Character. Name of the ecoregion attribute field used by \code{process_fun}.
#' @param corine_classes Integer/numeric vector. Optional CORINE classes passed to \code{process_fun}.
#' Interpretation depends on \code{process_fun}.
#' @param process_fun Function. Per-AOI processing function to run (required).
#' @param rescue_args Named list. Extra arguments forwarded to \code{process_fun} via \code{do.call()}.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{selected_detected}: sf. Detections selected for refinement (in \code{target_crs}).
#'   \item \code{aois}: sf. AOI polygons (in \code{target_crs}), possibly merged.
#'   \item \code{results}: list. Per-AOI outputs from \code{process_fun} (same order as \code{aois}).
#'   Elements are \code{NULL} when an AOI errors, and also \code{NULL} when skipped due to an existing \code{DONE.ok}.
#' }
#'
#' @note
#' The \pkg{whitebox} package is only required if your \code{process_fun} uses
#' \code{grow_engine = "whitebox"} (WhiteboxTools-based region growing). If you
#' use \code{grow_engine = "terra"}, \pkg{whitebox} is not needed.
#'
#' @examples
#' \dontrun{
#' # Run a segmentation function only inside AOIs built from detected polygons.
#' process_fun <- process_otsu_rasters_grow
#'
#' res <- segmentation_refinement(
#'   n_workers = 4,
#'   raster_path = "RBR_2012.tif",
#'   base_output_dir = "OUTPUT/REFINE_2012",
#'   year = 2012,
#'   detected_polys = "OUTPUT/STAGE1_2012_detected.gpkg",
#'   reference_polys = NULL,
#'   target_crs = "EPSG:3035",
#'   aoi_buffer_m = 5000,
#'   top_n_if_no_ref = Inf,
#'   min_detected_area_m2 = 10,
#'   merge_overlaps = TRUE,
#'   corine_raster_path = "AUX/CORINE_2018.tif",
#'   corine_classes = 1:12,
#'   ecoregion_shapefile_path = "AUX/ECOREGIONS.gpkg",
#'   ecoregion_field = "ECO_ID",
#'   process_fun = process_fun,
#'   rescue_args = list(
#'     trim_percentiles = list(min = 0.05, max = 0.99),
#'     otsu_value_range = c(0, 1500),
#'     otsu_thresholds = c(300),
#'     grow_delta = 100,
#'     min_grow_threshold_value = 200,
#'     min_seed_pixels_per_component = 5L,
#'     grow_engine = "whitebox",
#'     wbt_diag = FALSE,
#'     segment_by_intersection = TRUE,
#'     crop_to_units = TRUE,
#'     tile = FALSE,
#'     resolution = 90
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link[future]{plan}}, \code{\link[future.apply]{future_lapply}},
#' \code{\link[progressr]{with_progress}}
#'
#' @importFrom sf st_as_sf st_sf st_sfc st_crs st_transform st_make_valid st_is_valid
#' @importFrom sf st_geometry st_is_empty st_cast st_collection_extract st_union st_difference
#' @importFrom sf st_intersection st_bbox st_as_sfc st_buffer st_intersects st_area st_is_longlat
#' @importFrom sf st_polygon st_read st_write
#' @importFrom terra rast crs ext as.polygons vect project same.crs crop mask ifel values res writeRaster
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#' @importFrom parallelly availableCores
#' @importFrom tools file_path_sans_ext
#' @export
utils::globalVariables(c("AREA_M2", "REFINE_REASON", "AOI_ID", "AOI_GROUP"))

segmentation_refinement <- function(
    n_workers = 6,
    raster_path,
    base_output_dir,
    year,
    detected_polys = NULL,
    binary_final_raster = NULL,
    reference_polys = NULL,

    target_crs,

    aoi_buffer_m = 5000,
    omission_buffer_m = 0,
    top_n_if_no_ref = 50,
    min_detected_area_m2 = 0,
    clip_aois_to_raster_extent = TRUE,
    merge_overlaps = TRUE,
    merge_overlaps_buffer_m = 0,
    verbose = TRUE,

    corine_raster_path = NULL,
    ecoregion_shapefile_path = NULL,
    ecoregion_field = NULL,
    corine_classes = NULL,

    process_fun,
    rescue_args = list(),

    write_aoi_per_folder_shp = TRUE,
    standardize_ba_output = TRUE,
    keep_original_ba = FALSE,
    stop_on_unused_args = TRUE
) {

  msg <- function(...) if (isTRUE(verbose)) message(...)

  # ----------------------------
  # Package checks
  # ----------------------------
  pkgs <- c("sf", "terra", "future", "future.apply", "progressr", "parallelly", "tools", "utils")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing required packages: ", paste(miss, collapse = ", "))

  if (missing(process_fun) || is.null(process_fun) || !is.function(process_fun)) {
    stop("process_fun must be a valid function.")
  }
  if (missing(target_crs) || is.null(target_crs)) {
    stop("target_crs must be provided (projected CRS in meters).")
  }
  if (missing(raster_path) || is.null(raster_path) || !file.exists(raster_path)) {
    stop("raster_path must exist: ", raster_path)
  }
  if (missing(base_output_dir) || is.null(base_output_dir) || !nzchar(base_output_dir)) {
    stop("base_output_dir must be a valid directory path.")
  }

  if (!is.list(rescue_args)) stop("rescue_args must be a list().")
  if (length(rescue_args) > 0 && (is.null(names(rescue_args)) || any(!nzchar(names(rescue_args))))) {
    stop("rescue_args must be a named list (all elements must have names).")
  }

  # ----------------------------
  # Internal helpers
  # ----------------------------
  delete_shapefile_bundle <- function(shp_path) {
    base <- tools::file_path_sans_ext(shp_path)
    exts <- c(".shp",".shx",".dbf",".prj",".cpg",".qix",".fix",".sbn",".sbx",".xml")
    for (e in exts) {
      f <- paste0(base, e)
      if (file.exists(f)) try(file.remove(f), silent = TRUE)
    }
    invisible(TRUE)
  }

  write_shp_clean <- function(x_sf, shp_path, quiet = TRUE) {
    delete_shapefile_bundle(shp_path)
    sf::st_write(x_sf, shp_path, driver = "ESRI Shapefile", quiet = quiet)
    shp_path
  }

  copy_shapefile_bundle <- function(src_shp, dst_shp) {
    delete_shapefile_bundle(dst_shp)
    src_base <- tools::file_path_sans_ext(src_shp)
    dst_base <- tools::file_path_sans_ext(dst_shp)
    exts <- c(".shp",".shx",".dbf",".prj",".cpg",".qix",".fix",".sbn",".sbx",".xml")
    for (e in exts) {
      fs <- paste0(src_base, e)
      if (file.exists(fs)) file.copy(fs, paste0(dst_base, e), overwrite = TRUE)
    }
    dst_shp
  }

  move_shapefile_bundle <- function(src_shp, dst_shp) {
    delete_shapefile_bundle(dst_shp)
    src_base <- tools::file_path_sans_ext(src_shp)
    dst_base <- tools::file_path_sans_ext(dst_shp)
    exts <- c(".shp",".shx",".dbf",".prj",".cpg",".qix",".fix",".sbn",".sbx",".xml")

    for (e in exts) {
      fs <- paste0(src_base, e)
      fd <- paste0(dst_base, e)
      if (file.exists(fs)) {
        ok <- suppressWarnings(file.rename(fs, fd))
        if (!isTRUE(ok)) {
          file.copy(fs, fd, overwrite = TRUE)
          try(file.remove(fs), silent = TRUE)
        }
      }
    }
    dst_shp
  }

  safe_make_valid_polys <- function(x) {
    if (is.null(x)) return(NULL)

    if (inherits(x, "SpatVector")) {
      x <- sf::st_as_sf(x)
    } else if (inherits(x, "sfc")) {
      x <- sf::st_sf(geometry = x)
    } else {
      x <- sf::st_as_sf(x)
    }

    if (is.na(sf::st_crs(x))) stop("Input polygons have NA CRS. Assign CRS before processing.")

    emp <- sf::st_is_empty(sf::st_geometry(x))
    if (any(emp)) x <- x[!emp, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)

    inv <- !sf::st_is_valid(x)
    if (any(inv)) x <- sf::st_make_valid(x)

    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
    if (is.null(x) || nrow(x) == 0) return(NULL)

    sf::st_geometry(x) <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)

    emp2 <- sf::st_is_empty(sf::st_geometry(x))
    if (any(emp2)) x <- x[!emp2, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)

    x
  }

  read_sf_polys <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "sf") || inherits(x, "sfc") || inherits(x, "SpatVector")) {
      return(safe_make_valid_polys(x))
    }
    if (is.character(x) && length(x) == 1 && file.exists(x)) {
      out <- sf::st_read(x, quiet = TRUE)
      return(safe_make_valid_polys(out))
    }
    stop("Input polygons must be an sf/sfc/SpatVector object or a valid file path.")
  }

  polys_from_binary_raster <- function(binary_raster_path, dissolve = TRUE) {
    if (is.null(binary_raster_path) || !file.exists(binary_raster_path)) {
      stop("binary_final_raster does not exist: ", binary_raster_path)
    }
    r <- terra::rast(binary_raster_path)[[1]]
    r1 <- terra::ifel(r == 1, 1, NA)
    vv <- terra::values(r1, mat = FALSE)
    if (is.null(vv) || all(is.na(vv))) return(NULL)
    v <- terra::as.polygons(r1, dissolve = dissolve, na.rm = TRUE)
    if (is.null(v) || nrow(v) == 0) return(NULL)
    safe_make_valid_polys(sf::st_as_sf(v))
  }

  select_detected_for_refine <- function(det_polys,
                                         ref_polys = NULL,
                                         omission_buffer_m = 0,
                                         top_n_if_no_ref = 50,
                                         min_area_m2 = 0) {

    det_polys <- safe_make_valid_polys(det_polys)
    if (is.null(det_polys) || nrow(det_polys) == 0) return(NULL)

    det_area <- as.numeric(sf::st_area(det_polys))
    det_polys <- det_polys[det_area >= min_area_m2, , drop = FALSE]
    if (nrow(det_polys) == 0) return(NULL)

    if (is.null(ref_polys)) {
      det_polys$AREA_M2 <- as.numeric(sf::st_area(det_polys))
      det_polys <- det_polys[order(det_polys$AREA_M2, decreasing = TRUE), ]
      if (is.finite(top_n_if_no_ref) && nrow(det_polys) > top_n_if_no_ref) det_polys <- det_polys[1:top_n_if_no_ref, ]
      det_polys$REFINE_REASON <- "no_reference_top_area"
      return(det_polys)
    }

    ref_polys <- safe_make_valid_polys(ref_polys)
    if (is.null(ref_polys) || nrow(ref_polys) == 0) {
      det_polys$AREA_M2 <- as.numeric(sf::st_area(det_polys))
      det_polys <- det_polys[order(det_polys$AREA_M2, decreasing = TRUE), ]
      if (is.finite(top_n_if_no_ref) && nrow(det_polys) > top_n_if_no_ref) det_polys <- det_polys[1:top_n_if_no_ref, ]
      det_polys$REFINE_REASON <- "reference_empty_fallback_top_area"
      return(det_polys)
    }

    det_u <- sf::st_union(det_polys)
    omis <- suppressWarnings(sf::st_difference(sf::st_make_valid(ref_polys), sf::st_make_valid(det_u)))
    omis <- sf::st_collection_extract(omis, "POLYGON", warn = FALSE)

    if (is.null(omis) || length(omis) == 0) {
      det_polys$AREA_M2 <- as.numeric(sf::st_area(det_polys))
      det_polys <- det_polys[order(det_polys$AREA_M2, decreasing = TRUE), ]
      if (is.finite(top_n_if_no_ref) && nrow(det_polys) > top_n_if_no_ref) det_polys <- det_polys[1:top_n_if_no_ref, ]
      det_polys$REFINE_REASON <- "no_omission_fallback_top_area"
      return(det_polys)
    }

    omis_sf <- safe_make_valid_polys(sf::st_as_sf(omis))
    if (!is.null(omis_sf) && nrow(omis_sf) > 0 && omission_buffer_m > 0) {
      omis_sf <- sf::st_buffer(omis_sf, omission_buffer_m)
    }

    hit <- sf::st_intersects(det_polys, omis_sf, sparse = TRUE)
    det_sel <- det_polys[lengths(hit) > 0, , drop = FALSE]

    if (nrow(det_sel) == 0) {
      det_polys$AREA_M2 <- as.numeric(sf::st_area(det_polys))
      det_polys <- det_polys[order(det_polys$AREA_M2, decreasing = TRUE), ]
      if (is.finite(top_n_if_no_ref) && nrow(det_polys) > top_n_if_no_ref) det_polys <- det_polys[1:top_n_if_no_ref, ]
      det_polys$REFINE_REASON <- "omission_no_intersect_fallback_top_area"
      return(det_polys)
    }

    det_sel$REFINE_REASON <- "intersects_omission"
    det_sel
  }

  normalize_target_crs <- function(target_crs) {
    if (inherits(target_crs, "crs")) {
      crs_sf <- target_crs
    } else if (is.numeric(target_crs) && length(target_crs) == 1) {
      crs_sf <- sf::st_crs(as.integer(target_crs))
    } else if (is.character(target_crs) && length(target_crs) == 1) {
      if (grepl("^EPSG:\\d+$", target_crs, ignore.case = TRUE)) {
        epsg <- as.integer(sub("EPSG:", "", toupper(target_crs)))
        crs_sf <- sf::st_crs(epsg)
      } else {
        crs_sf <- sf::st_crs(target_crs)
      }
    } else {
      stop("target_crs must be numeric EPSG, 'EPSG:xxxx', WKT string, or sf::st_crs().")
    }

    if (is.na(crs_sf)) stop("target_crs could not be parsed to a valid sf CRS.")
    if (isTRUE(sf::st_is_longlat(crs_sf))) stop("target_crs must be projected (meters), not lon/lat.")

    wkt <- crs_sf$wkt
    if (is.null(wkt) || !nzchar(wkt)) stop("target_crs has empty WKT (cannot pass to terra).")

    list(crs_sf = crs_sf, wkt = wkt)
  }

  make_bbox_aois_projected_fast <- function(det_sel,
                                            buffer_m = 5000,
                                            raster_path,
                                            target_crs,
                                            clip_to_raster_bbox = TRUE,
                                            verbose = TRUE) {

    if (is.null(det_sel) || nrow(det_sel) == 0) return(NULL)

    crs_target <- sf::st_crs(target_crs)
    if (is.na(crs_target)) stop("[AOI] target_crs is NA/invalid.")

    det_p <- if (sf::st_crs(det_sel) != crs_target) sf::st_transform(det_sel, crs_target) else det_sel

    rb <- NULL
    if (isTRUE(clip_to_raster_bbox)) {
      r <- terra::rast(raster_path)
      ext_poly <- terra::as.polygons(terra::ext(r))
      terra::crs(ext_poly) <- terra::crs(r, proj = TRUE)
      ext_sf <- sf::st_as_sf(ext_poly)
      sf::st_crs(ext_sf) <- sf::st_crs(terra::crs(r, proj = TRUE))
      ext_sf <- sf::st_transform(ext_sf, crs_target)
      bb <- sf::st_bbox(ext_sf)
      rb <- c(xmin = as.numeric(bb["xmin"]),
              ymin = as.numeric(bb["ymin"]),
              xmax = as.numeric(bb["xmax"]),
              ymax = as.numeric(bb["ymax"]))
    }

    geoms <- sf::st_geometry(det_p)
    n <- length(geoms)

    xmin <- numeric(n); xmax <- numeric(n); ymin <- numeric(n); ymax <- numeric(n)
    for (i in seq_len(n)) {
      b <- sf::st_bbox(geoms[[i]])
      xmin[i] <- as.numeric(b["xmin"]) - buffer_m
      xmax[i] <- as.numeric(b["xmax"]) + buffer_m
      ymin[i] <- as.numeric(b["ymin"]) - buffer_m
      ymax[i] <- as.numeric(b["ymax"]) + buffer_m
    }

    if (!is.null(rb)) {
      xmin <- pmax(xmin, rb["xmin"])
      xmax <- pmin(xmax, rb["xmax"])
      ymin <- pmax(ymin, rb["ymin"])
      ymax <- pmin(ymax, rb["ymax"])
    }

    ok <- is.finite(xmin) & is.finite(xmax) & is.finite(ymin) & is.finite(ymax) & (xmin < xmax) & (ymin < ymax)
    if (!any(ok)) stop("[AOI] All AOIs dropped (invalid bboxes).")

    idx <- which(ok)

    polys <- lapply(idx, function(k) {
      m <- matrix(c(
        xmin[k], ymin[k],
        xmax[k], ymin[k],
        xmax[k], ymax[k],
        xmin[k], ymax[k],
        xmin[k], ymin[k]
      ), ncol = 2, byrow = TRUE)
      sf::st_polygon(list(m))
    })

    sfc <- sf::st_sfc(polys, crs = crs_target)
    aoi_sf <- sf::st_sf(AOI_ID = seq_along(idx), geometry = sfc)

    if (isTRUE(verbose)) {
      message(sprintf("[AOI] built=%d (from %d input) | target_crs=%s",
                      nrow(aoi_sf), n, crs_target$input))
    }

    aoi_sf
  }

  merge_overlapping_aois <- function(aoi_sf, merge_buffer_m = 0, verbose = TRUE) {
    aoi_sf <- safe_make_valid_polys(aoi_sf)
    if (is.null(aoi_sf) || nrow(aoi_sf) == 0) return(aoi_sf)

    crs0 <- sf::st_crs(aoi_sf)
    if (is.na(crs0)) stop("[AOI] AOI CRS is NA inside merge_overlapping_aois().")

    g_group <- sf::st_geometry(aoi_sf)
    if (is.numeric(merge_buffer_m) && length(merge_buffer_m) == 1 && merge_buffer_m > 0) {
      g_group <- sf::st_buffer(g_group, merge_buffer_m)
    }

    ix <- sf::st_intersects(g_group, sparse = TRUE)
    n <- length(ix)
    if (n == 0) return(aoi_sf)

    comp <- integer(n)
    cid <- 0L

    for (i in seq_len(n)) {
      if (comp[i] != 0L) next
      cid <- cid + 1L
      queue <- i
      comp[i] <- cid

      while (length(queue) > 0) {
        v <- queue[[1]]
        queue <- queue[-1]
        neigh <- ix[[v]]
        if (length(neigh) == 0) next
        for (u in neigh) {
          if (comp[u] == 0L) {
            comp[u] <- cid
            queue <- c(queue, u)
          }
        }
      }
    }

    aoi_sf$AOI_GROUP <- comp
    groups <- split(aoi_sf, aoi_sf$AOI_GROUP)

    merged_list <- lapply(groups, function(x) {
      g <- sf::st_union(sf::st_make_valid(sf::st_geometry(x)))
      g_poly <- sf::st_collection_extract(sf::st_sfc(g, crs = crs0), "POLYGON", warn = FALSE)
      if (length(g_poly) == 0) return(NULL)

      g2 <- sf::st_union(g_poly)
      g2_sfc <- sf::st_sfc(g2, crs = crs0)
      g2_sfc <- sf::st_cast(g2_sfc, "MULTIPOLYGON", warn = FALSE)
      if (length(g2_sfc) == 0 || sf::st_is_empty(g2_sfc)[1]) return(NULL)

      sf::st_sf(AOI_GROUP = unique(x$AOI_GROUP)[1], geometry = g2_sfc)
    })

    merged_list <- Filter(Negate(is.null), merged_list)
    if (length(merged_list) == 0) stop("[AOI] merge_overlapping_aois produced no valid merged AOIs.")

    merged <- do.call(rbind, merged_list)
    sf::st_geometry(merged) <- sf::st_make_valid(sf::st_geometry(merged))
    emp <- sf::st_is_empty(sf::st_geometry(merged))
    if (any(emp)) merged <- merged[!emp, , drop = FALSE]
    if (nrow(merged) == 0) stop("[AOI] merged AOIs are empty after union/make_valid.")

    merged$AOI_ID <- seq_len(nrow(merged))

    if (isTRUE(verbose)) message(sprintf("[AOI] Merged overlaps: %d -> %d AOIs", nrow(aoi_sf), nrow(merged)))
    merged
  }

  crop_project_raster_to_aoi <- function(raster_path,
                                         aoi_geom,
                                         out_path_src,
                                         out_path_proj,
                                         target_crs,
                                         method_project = "bilinear",
                                         project_res = NULL) {

    crs_info <- normalize_target_crs(target_crs)
    target_wkt <- crs_info$wkt

    r <- terra::rast(raster_path)

    aoi_sf <- safe_make_valid_polys(aoi_geom)
    if (is.null(aoi_sf) || nrow(aoi_sf) == 0) stop("[CROP] AOI geometry is empty/invalid.")
    v <- terra::vect(aoi_sf)

    if (!isTRUE(terra::same.crs(v, r))) v <- terra::project(v, r)

    rr <- terra::crop(r, v)
    rr <- terra::mask(rr, v)
    terra::writeRaster(rr, out_path_src, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

    same_crs <- isTRUE(terra::same.crs(rr, target_wkt))
    same_res <- TRUE
    if (!is.null(project_res)) same_res <- isTRUE(all(terra::res(rr) == project_res))

    if (same_crs && same_res) {
      rr_p <- rr
    } else {
      rr_p <- if (!is.null(project_res)) {
        terra::project(rr, target_wkt, method = method_project, res = project_res)
      } else {
        terra::project(rr, target_wkt, method = method_project)
      }
    }

    terra::writeRaster(rr_p, out_path_proj, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    out_path_proj
  }

  crop_align_raster_to_aoi <- function(aux_raster_path, aoi_geom, template_raster_path,
                                       out_path, method = "near") {

    aux <- terra::rast(aux_raster_path)
    tmpl <- terra::rast(template_raster_path)

    aoi_sf <- safe_make_valid_polys(aoi_geom)
    if (is.null(aoi_sf) || nrow(aoi_sf) == 0) stop("[AUX] AOI geometry is empty/invalid.")
    v <- terra::vect(aoi_sf)

    v_aux <- if (!isTRUE(terra::same.crs(v, aux))) terra::project(v, aux) else v
    aux_crop <- terra::crop(aux, v_aux)
    aux_crop <- terra::mask(aux_crop, v_aux)

    aux_aligned <- terra::project(aux_crop, tmpl, method = method)
    terra::writeRaster(aux_aligned, out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    out_path
  }

  clip_ecoregions_to_aoi <- function(ecoregion_shapefile_path,
                                     aoi_geom,
                                     out_path,
                                     quiet = TRUE) {
    if (is.null(ecoregion_shapefile_path) || !file.exists(ecoregion_shapefile_path)) return(NULL)

    aoi_sf <- safe_make_valid_polys(aoi_geom)
    if (is.null(aoi_sf) || nrow(aoi_sf) == 0) return(NULL)

    aoi_u <- sf::st_union(sf::st_make_valid(sf::st_geometry(aoi_sf)))
    aoi_u <- sf::st_collection_extract(sf::st_sfc(aoi_u, crs = sf::st_crs(aoi_sf)), "POLYGON", warn = FALSE)
    if (length(aoi_u) == 0 || sf::st_is_empty(aoi_u)[1]) return(NULL)
    aoi_u <- sf::st_union(aoi_u)
    aoi_u <- sf::st_sfc(aoi_u, crs = sf::st_crs(aoi_sf))

    eco <- sf::st_read(ecoregion_shapefile_path, quiet = quiet)
    eco <- safe_make_valid_polys(eco)
    if (is.null(eco) || nrow(eco) == 0) return(NULL)

    eco <- sf::st_transform(eco, sf::st_crs(aoi_u))

    bb <- sf::st_bbox(aoi_u)
    eco_bb <- suppressWarnings(sf::st_intersects(eco, sf::st_as_sfc(bb), sparse = FALSE))
    if (is.matrix(eco_bb)) eco_bb <- eco_bb[, 1]
    eco <- eco[which(eco_bb), , drop = FALSE]
    if (nrow(eco) == 0) return(NULL)

    eco_clip <- suppressWarnings(sf::st_intersection(sf::st_make_valid(eco), sf::st_make_valid(sf::st_sf(geometry = aoi_u))))
    eco_clip <- safe_make_valid_polys(eco_clip)
    if (is.null(eco_clip) || nrow(eco_clip) == 0) return(NULL)

    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(out_path)) try(unlink(out_path, force = TRUE), silent = TRUE)
    sf::st_write(eco_clip, out_path, driver = "GPKG", delete_dsn = TRUE, quiet = quiet)
    out_path
  }

  # Standardize per-AOI BA output WITHOUT duplicates by default (move/rename bundle).
  ensure_ba_aoi_shp <- function(aoi_dir, i, keep_original = FALSE) {

    dst <- file.path(aoi_dir, sprintf("BA_AOI_%03d.shp", i))

    # If already standardized, return
    if (file.exists(dst)) return(dst)

    shp_all <- list.files(aoi_dir, pattern = "\\.shp$", full.names = TRUE)
    bn <- basename(shp_all)

    # Exclude AOI geometry and already standardized outputs
    keep <- !grepl("^AOI_", bn) &
      !grepl("^AOIs_", bn) &
      !grepl("^BA_AOI_\\d+\\.shp$", bn)

    shp_all <- shp_all[keep]
    bn <- basename(shp_all)

    if (length(shp_all) > 0) {

      # Prefer BA-like names; else use all candidates
      ba_cand <- shp_all[grepl("^BA", bn, ignore.case = TRUE) | grepl("BA_", bn, ignore.case = TRUE)]
      cand <- if (length(ba_cand) > 0) ba_cand else shp_all

      # Pick largest shapefile (.shp size only; good proxy)
      sz <- file.info(cand)$size
      src <- cand[which.max(sz)]

      # If keep_original=FALSE, move/rename bundle to avoid duplicates
      if (isTRUE(keep_original)) {
        copy_shapefile_bundle(src, dst)
      } else {
        move_shapefile_bundle(src, dst)
      }

      return(dst)
    }

    # Fallback: look for GPKG and export the most likely layer
    gpkg_all <- list.files(aoi_dir, pattern = "\\.gpkg$", full.names = TRUE)
    gpkg_all <- gpkg_all[!grepl("AOIs_debug\\.gpkg$", basename(gpkg_all))]

    if (length(gpkg_all) > 0) {

      pick <- gpkg_all[grepl("BA", basename(gpkg_all), ignore.case = TRUE)]
      srcg <- if (length(pick) > 0) pick[1] else gpkg_all[which.max(file.info(gpkg_all)$size)]

      lyr_info <- sf::st_layers(srcg)
      if (!is.null(lyr_info) && nrow(lyr_info) > 0) {

        prefer_names <- c("polys","POLYS","ba","BA","scored","SCORED")
        lyr_names <- lyr_info$name

        idx_pref <- which(lyr_names %in% prefer_names)
        if (length(idx_pref) > 0) {
          lyr <- lyr_names[idx_pref[1]]
        } else if ("features" %in% names(lyr_info)) {
          lyr <- lyr_info$name[which.max(lyr_info$features)]
        } else {
          lyr <- lyr_names[1]
        }

        x <- sf::st_read(srcg, layer = lyr, quiet = TRUE)
        x <- safe_make_valid_polys(x)

        if (!is.null(x) && nrow(x) > 0) {
          write_shp_clean(x, dst, quiet = TRUE)
          return(dst)
        }
      }
    }

    NULL
  }

  strip_ecoregion_cols <- function(x_sf,
                                   ecoregion_field = NULL,
                                   eco_layer_path = NULL,
                                   verbose = FALSE) {
    stopifnot(inherits(x_sf, "sf"))
    nms <- names(x_sf)
    geom_col <- attr(x_sf, "sf_column")
    if (is.null(geom_col) || !nzchar(geom_col)) geom_col <- "geometry"

    msg <- function(...) if (isTRUE(verbose)) message(...)

    # ---- 1) Pattern-based drop (common eco/intersection fields) ----
    pat <- paste0(
      "(",
      "^ens($|[._])|^enz($|[._])|",
      "ecoreg|ecoregion|ecorreg|",
      "(^|_)eco(_|$)|",
      "cor_eco|coreco|corine_eco|",
      "intersection_eco|eco_inter",
      ")"
    )
    drop_pat <- grep(pat, nms, ignore.case = TRUE, value = TRUE)

    # ---- 2) Exact eco field names (read from clipped eco layer) ----
    eco_fields <- character(0)
    if (!is.null(eco_layer_path) && is.character(eco_layer_path) &&
        length(eco_layer_path) == 1 && file.exists(eco_layer_path)) {

      # Read only schema names cheaply via st_layers + st_read(quiet=TRUE)
      eco_tmp <- try(sf::st_read(eco_layer_path, quiet = TRUE), silent = TRUE)
      if (!inherits(eco_tmp, "try-error")) {
        eco_fields <- setdiff(names(eco_tmp), attr(eco_tmp, "sf_column"))
      }
    }

    # Include user-provided ecoregion_field too
    base_fields <- unique(c(ecoregion_field, eco_fields))
    base_fields <- base_fields[!is.na(base_fields) & nzchar(base_fields)]

    # Build variants that commonly appear after joins/intersections and after Shapefile truncation
    variants <- unique(unlist(lapply(base_fields, function(f) {
      f0 <- f
      f1 <- tolower(f0)
      f2 <- toupper(f0)
      f3 <- gsub("[^a-zA-Z0-9_]+", "_", f0)

      # Shapefile truncation to 10 chars
      trunc10 <- unique(c(substr(f0, 1, 10), substr(f1, 1, 10), substr(f2, 1, 10), substr(f3, 1, 10)))

      unique(c(f0, f1, f2, f3, trunc10))
    })))

    # Match columns that are exactly the variant OR start with it plus a join suffix
    drop_by_fields <- unique(unlist(lapply(variants, function(v) {
      if (!nzchar(v)) return(character(0))
      rx <- paste0("^", gsub("([\\.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1", v),
                   "($|[._].+|_[0-9]+$|_x$|_y$|\\.x$|\\.y$)")
      grep(rx, nms, ignore.case = TRUE, value = TRUE)
    })))

    drop <- unique(c(drop_pat, drop_by_fields))

    # Never drop geometry
    drop <- setdiff(drop, geom_col)

    if (length(drop) == 0) return(x_sf)

    msg("[strip_ecoregion_cols] Dropping: ", paste(drop, collapse = ", "))
    x_sf[, setdiff(nms, drop), drop = FALSE]
  }

  # ----------------------------
  # Setup
  # ----------------------------
  crs_info <- normalize_target_crs(target_crs)
  target_crs_sf <- crs_info$crs_sf

  dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

  old_opts <- options(progressr.enable = TRUE)
  on.exit(options(old_opts), add = TRUE)

  # ----------------------------
  # Load detected polygons
  # ----------------------------
  det_polys <- NULL
  if (!is.null(detected_polys)) {
    det_polys <- read_sf_polys(detected_polys)
  } else if (!is.null(binary_final_raster)) {
    det_polys <- polys_from_binary_raster(binary_final_raster, dissolve = FALSE)
  }
  if (is.null(det_polys) || nrow(det_polys) == 0) {
    stop("No detected polygons available. Provide detected_polys or binary_final_raster.")
  }

  det_polys <- sf::st_transform(det_polys, target_crs_sf)

  # Optional reference
  ref_polys <- NULL
  if (!is.null(reference_polys)) {
    ref_polys <- read_sf_polys(reference_polys)
    if (!is.null(ref_polys) && nrow(ref_polys) > 0) {
      ref_polys <- sf::st_transform(ref_polys, target_crs_sf)
    } else {
      ref_polys <- NULL
    }
  }

  det_sel <- select_detected_for_refine(
    det_polys = det_polys,
    ref_polys = ref_polys,
    omission_buffer_m = omission_buffer_m,
    top_n_if_no_ref = top_n_if_no_ref,
    min_area_m2 = min_detected_area_m2
  )
  if (is.null(det_sel) || nrow(det_sel) == 0) stop("No detected polygons selected for refinement.")

  # Build AOIs
  aoi_sf <- make_bbox_aois_projected_fast(
    det_sel = det_sel,
    buffer_m = aoi_buffer_m,
    raster_path = raster_path,
    target_crs = target_crs_sf,
    clip_to_raster_bbox = isTRUE(clip_aois_to_raster_extent),
    verbose = verbose
  )
  if (is.null(aoi_sf) || nrow(aoi_sf) == 0) stop("Failed to build AOIs.")

  if (isTRUE(merge_overlaps)) {
    aoi_sf <- merge_overlapping_aois(aoi_sf, merge_buffer_m = merge_overlaps_buffer_m, verbose = verbose)
    if (is.null(aoi_sf) || nrow(aoi_sf) == 0) stop("Failed to merge AOIs.")
    suppressWarnings(sf::st_write(aoi_sf, file.path(base_output_dir, "AOIs_debug.gpkg"), delete_dsn = TRUE, quiet = TRUE))
  }

  n_aoi <- nrow(aoi_sf)

  # Projected resolution (optional)
  project_res <- NULL
  if (!is.null(rescue_args$resolution) && is.numeric(rescue_args$resolution) && length(rescue_args$resolution) == 1) {
    project_res <- rescue_args$resolution
  }

  # ----------------------------
  # PRE-FLIGHT: argument compatibility
  # ----------------------------
  base_args <- list(
    raster_path = "DUMMY",
    output_dir = "DUMMY",
    year = year,
    corine_raster_path = "DUMMY",
    reclassify_corine = FALSE,
    corine_classes = corine_classes,
    ecoregion_shapefile_path = "DUMMY",
    ecoregion_field = ecoregion_field
  )
  call_args_preview <- c(base_args, rescue_args)

  fn_formals <- names(formals(process_fun))
  has_dots <- "..." %in% fn_formals

  if (!has_dots) {
    extras <- setdiff(names(call_args_preview), fn_formals)
    if (length(extras) > 0 && isTRUE(stop_on_unused_args)) {
      stop(
        "process_fun does not accept some arguments you are passing (this can trigger 'unused arguments' and prevent BA outputs).\n",
        "Unexpected args:\n- ", paste(extras, collapse = "\n- "), "\n\n",
        "Fix: verify you are calling the correct process_fun (e.g., OtsuFire::process_otsu_rasters) and that its signature includes these parameters.\n",
        "Tip: check `names(formals(process_fun))` and `environmentName(environment(process_fun))`."
      )
    }
  }

  # ----------------------------
  # Parallel plan
  # ----------------------------
  n_workers <- min(n_workers, max(1, parallelly::availableCores() - 1))
  old_plan <- future::plan()
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(old_plan), add = TRUE)

  # ----------------------------
  # Run per AOI
  # ----------------------------
  out_list <- progressr::with_progress({
    p <- progressr::progressor(along = seq_len(n_aoi))

    future.apply::future_lapply(seq_len(n_aoi), function(i) {

      aoi_i <- aoi_sf[i, , drop = FALSE]
      aoi_dir <- file.path(base_output_dir, sprintf("AOI_%03d", i))
      dir.create(aoi_dir, recursive = TRUE, showWarnings = FALSE)

      running_flag <- file.path(aoi_dir, "RUNNING.ok")
      done_flag <- file.path(aoi_dir, "DONE.ok")
      err_flag <- file.path(aoi_dir, "ERROR.txt")

      if (file.exists(done_flag)) {
        p(sprintf("AOI %d/%d (skipped: DONE)", i, n_aoi))
        return(invisible(NULL))
      }

      file.create(running_flag)

      res <- tryCatch({

        # Save AOI geometry per folder
        if (isTRUE(write_aoi_per_folder_shp)) {
          aoi_shp <- file.path(aoi_dir, sprintf("AOI_%03d.shp", i))
          write_shp_clean(aoi_i, aoi_shp, quiet = TRUE)
        }

        # Crop + project raster
        r_crop_src <- file.path(aoi_dir, "rbr_crop_src.tif")
        r_crop_proj <- file.path(aoi_dir, "rbr_crop_proj.tif")

        rbr_proj_path <- crop_project_raster_to_aoi(
          raster_path = raster_path,
          aoi_geom = aoi_i,
          out_path_src = r_crop_src,
          out_path_proj = r_crop_proj,
          target_crs = target_crs_sf,
          method_project = "bilinear",
          project_res = project_res
        )

        # Align CORINE per AOI
        cor_crop_path <- NULL
        if (!is.null(corine_raster_path) && file.exists(corine_raster_path)) {
          cor_crop_path <- file.path(aoi_dir, "corine_crop_aligned.tif")
          crop_align_raster_to_aoi(
            aux_raster_path = corine_raster_path,
            aoi_geom = aoi_i,
            template_raster_path = rbr_proj_path,
            out_path = cor_crop_path,
            method = "near"
          )
        }

        # Clip ecoregions per AOI
        eco_crop_path <- NULL
        if (!is.null(ecoregion_shapefile_path) && file.exists(ecoregion_shapefile_path)) {
          eco_crop_path <- file.path(aoi_dir, "ecoregions_crop.gpkg")
          eco_crop_path <- clip_ecoregions_to_aoi(
            ecoregion_shapefile_path = ecoregion_shapefile_path,
            aoi_geom = aoi_i,
            out_path = eco_crop_path,
            quiet = TRUE
          )
        }

        base_call_args <- list(
          raster_path = rbr_proj_path,
          output_dir = aoi_dir,
          year = year,
          corine_raster_path = cor_crop_path,
          reclassify_corine = FALSE,
          corine_classes = corine_classes,
          ecoregion_shapefile_path = eco_crop_path,
          ecoregion_field = ecoregion_field
        )

        call_args <- c(base_call_args, rescue_args)

        # Filter args if no "..." (avoid unused arguments crash)
        fn_formals_local <- names(formals(process_fun))
        has_dots_local <- "..." %in% fn_formals_local

        dropped <- character(0)
        if (!has_dots_local) {
          dropped <- setdiff(names(call_args), fn_formals_local)
          call_args <- call_args[names(call_args) %in% fn_formals_local]
        }

        # Debug files per AOI
        writeLines(capture.output(str(call_args)), file.path(aoi_dir, "CALL_ARGS_USED.txt"))
        if (length(dropped) > 0) writeLines(dropped, file.path(aoi_dir, "DROPPED_ARGS.txt"))

        out <- do.call(process_fun, call_args)

        # Standardize BA output name (and avoid duplicates by default)
        if (isTRUE(standardize_ba_output)) {
          ba_path <- ensure_ba_aoi_shp(aoi_dir, i, keep_original = isTRUE(keep_original_ba))

          if (!is.null(ba_path) && file.exists(ba_path)) {
            ba_sf <- sf::st_read(ba_path, quiet = TRUE)
            ba_sf <- safe_make_valid_polys(ba_sf)

            # Drop ecoregion/intersection columns robustly
            ba_sf <- strip_ecoregion_cols(
              ba_sf,
              ecoregion_field = ecoregion_field,
              eco_layer_path = eco_crop_path,
              verbose = FALSE
            )

            write_shp_clean(ba_sf, ba_path, quiet = TRUE)
          }
        }

        if (file.exists(running_flag)) file.remove(running_flag)
        file.create(done_flag)

        out

      }, error = function(e) {
        if (file.exists(running_flag)) file.remove(running_flag)
        writeLines(conditionMessage(e), con = err_flag)
        NULL
      })

      p(sprintf("AOI %d/%d", i, n_aoi))
      res

    }, future.seed = TRUE, future.packages = c("sf", "terra"))
  })

  list(
    selected_detected = det_sel,
    aois = aoi_sf,
    results = out_list
  )
}
