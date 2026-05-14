#' Stratified Otsu + seed-supported growth segmentation (CORINE x ecoregions)
#'
#' \code{process_otsu_rasters_grow()} implements a two-threshold burned-area segmentation:
#' (1) unit-specific Otsu thresholds are estimated per CORINE x ecoregion unit, and
#' (2) candidate components are grown from conservative seeds and retained only if
#' they contain enough seed support.
#'
#' @description
#' Workflow \code{"otsu_grow"} (implemented here) uses:
#' \itemize{
#'   \item \strong{Seed threshold (THR_SEED)}: per-unit Otsu estimate (\code{OTSU_RAW})
#'   constrained by a run-level global lower bound (\code{LB_GLOBAL}) and, optionally,
#'   a classwise lower bound (\code{LB_CLASS}).
#'   \item \strong{Growth threshold (THR_GROW)}: expands candidates down from seeds
#'   using a delta (\code{DELTA_USED}) and an optional minimum growth floor.
#' }
#' Candidate components are labeled on the growth mask and kept only if they contain
#' at least \code{min_seed_pixels_per_component} seed pixels (seed-supported growth).
#'
#' @details
#' \section{Required inputs (for \code{workflow="otsu_grow"})}{
#' \itemize{
#'   \item Provide either \code{raster_path} (single-band severity index) OR both
#'   \code{nbr_pre_path} and \code{nbr_post_path} (to compute \code{index_type}).
#'   \item \code{"otsu_grow"} requires stratification into units, therefore it needs
#'   \code{corine_raster_path}, \code{ecoregion_shapefile_path}, and
#'   \code{segment_by_intersection=TRUE}.
#' }
#' }
#'
#' \section{Units (CORINE x ecoregion)}{
#' Units are built as an integer unit raster (one value per CORINE x ecoregion combination).
#' Otsu is estimated independently per unit. Use \code{ecoregion_touches=TRUE} to reduce
#' boundary gaps when rasterizing ecoregions.
#' }
#'
#' \section{Otsu estimation (does NOT change raster values)}{
#' Otsu is estimated on real severity values using a smoothed-histogram approach:
#' values are linearly scaled to 0--255 for histogram-based Otsu and then mapped back
#' to real units. Two optional filters affect only the values used for estimation:
#' \itemize{
#'   \item \code{trim_percentiles}: trims tails (e.g., 0.05--0.99) before estimating Otsu.
#'   \item \code{otsu_value_range}: hard filter \code{c(lo, hi)} before estimating Otsu.
#' }
#' A single global Otsu threshold is also estimated from a large random sample and used
#' \emph{only} as a fallback when unit-level estimation is unreliable (small units, low
#' samples, or failed Otsu), if \code{fallback_to_global=TRUE}.
#' }
#'
#' \section{Seed threshold (recommended, non-redundant)}{
#' \strong{Important:} despite its name, \code{otsu_thresholds} does \emph{not} change how Otsu
#' is estimated. It defines run-specific \emph{global seed lower bounds} (LB_GLOBAL) and
#' produces one run per value.
#'
#' Recommended seed construction:
#' \preformatted{
#' THR_SEED = max(OTSU_RAW, LB_GLOBAL, LB_CLASS)
#' }
#' where:
#' \itemize{
#'   \item \code{OTSU_RAW}: unit-specific Otsu estimate (or global fallback if needed),
#'   \item \code{LB_GLOBAL}: run-level global seed lower bound from \code{otsu_thresholds},
#'   \item \code{LB_CLASS}: optional classwise seed lower bound from \code{otsu_min_by_class}.
#' }
#' \code{otsu_min_by_class} supports \code{$corine}, \code{$ecoregion}, and \code{$intersection}.
#' Each entry should be named by the class identifier (stored as character keys).
#' }
#'
#' \section{Growth threshold (expansion from seeds)}{
#' Growth expands candidate pixels down from THR_SEED:
#' \preformatted{
#' THR_GROW = max(THR_SEED - DELTA_USED, GROW_FLOOR_USED)
#' }
#' where:
#' \itemize{
#'   \item \code{DELTA_USED} comes from \code{grow_delta} (default) and can be overridden by
#'   \code{grow_delta_by_class}.
#'   \item \code{GROW_FLOOR_USED} comes from \code{min_grow_threshold_value} and/or
#'   \code{min_grow_threshold_by_class} (recommended when you want to prevent overly permissive growth).
#' }
#' Candidate components are labeled on the growth mask and retained only if they contain at
#' least \code{min_seed_pixels_per_component} seed pixels.
#' }
#'
#' \section{Advanced/legacy seed floors (usually redundant)}{
#' \code{min_otsu_threshold_value} and \code{min_otsu_threshold_by_class} provide an additional
#' \emph{seed floor} applied after lower bounds. In practice, they are often redundant if you
#' already use \code{otsu_thresholds} (LB_GLOBAL) and \code{otsu_min_by_class} (LB_CLASS).
#' They only change THR_SEED when they are more restrictive than both lower bounds.
#' For most users, leaving them as \code{NULL} avoids overlapping controls.
#' }
#'
#' \section{Outputs}{
#' For each run, the function writes:
#' \itemize{
#'   \item a vector file (\code{.shp} or \code{.geojson}) with one polygon per retained component, and
#'   \item a tab-delimited \code{*_log.txt} with unit-level thresholds, bounds applied, and Otsu sampling diagnostics.
#' }
#' The return value is a named list keyed by the run label.
#' }
#'
#' @name process_otsu_rasters_grow
#' @rdname process_otsu_rasters_grow
#'
#' @param raster_path Character. Path to a single-band severity raster (e.g., RBR or dNBR).
#' Provide this OR both \code{nbr_pre_path} and \code{nbr_post_path}.
#' @param nbr_pre_path Character. Optional path to pre-fire NBR raster (used only if \code{nbr_post_path} is also provided).
#' @param nbr_post_path Character. Optional path to post-fire NBR raster (used only if \code{nbr_pre_path} is also provided).
#' @param output_dir Character. Output directory. Created if it does not exist.
#' @param year Integer/character. Used for naming outputs and metadata fields. If NULL, inferred from \code{raster_path} or \code{output_dir}.
#'
#' @param otsu_thresholds Numeric vector. Run-level \strong{global seed lower bounds} (LB_GLOBAL). One run per value.
#' Larger values force more conservative seeds. Does not affect Otsu estimation (\code{OTSU_RAW}).
#' @param otsu_min_by_class List or NULL. Optional \strong{classwise seed lower bounds} (LB_CLASS) applied with LB_GLOBAL.
#' Expected structure: \code{list(corine=..., ecoregion=..., intersection=...)}. Values may be named lists or named vectors.
#'
#' @param trim_percentiles List or NULL. Optional trimming used only for Otsu estimation:
#' \code{list(min=..., max=...)} in [0,1] with \code{min < max}.
#' @param otsu_value_range Numeric length-2 or NULL. Optional hard filter \code{c(lo, hi)} used only for Otsu estimation.
#'
#' @param corine_raster_path Character or NULL. CORINE raster path (required for \code{workflow="otsu_grow"}).
#' @param reclassify_corine Logical. If TRUE, CORINE is cropped to \code{peninsula_shapefile}, reclassified with \code{reclass_matrix},
#' and optionally reprojected via \code{gdalwarp_path}.
#' @param reclass_matrix Matrix. Reclassification matrix (2 columns: value,new_value; or 3 columns: from,to,becomes). Used only if \code{reclassify_corine=TRUE}.
#' @param peninsula_shapefile Character. Polygon mask used to crop CORINE (and optionally ecoregions). Required if \code{reclassify_corine=TRUE}.
#' @param output_corine_raster_dir Character. Output directory for the reclassified CORINE raster (when \code{reclassify_corine=TRUE}).
#' @param reproject Logical. If TRUE and \code{reclassify_corine=TRUE}, reprojects CORINE output to EPSG:3035 via GDAL.
#' @param resolution Numeric. Target pixel size (meters) used when \code{reproject=TRUE}.
#'
#' @param ecoregion_shapefile_path Character or NULL. Ecoregion polygons path (required for \code{workflow="otsu_grow"}).
#' @param ecoregion_field Character or NULL. Attribute field identifying ecoregion classes. If NULL, tries \code{"EnS_name"} then \code{"EnZ_name"}.
#' @param ecoregion_classes Vector or NULL. Optional subset of ecoregion classes to keep.
#' @param ecoregion_touches Logical. Passed to \code{terra::rasterize(touches=...)} to reduce boundary gaps.
#'
#' @param grow_delta Numeric. Default delta for growth: \code{THR_GROW = THR_SEED - DELTA}. Length 1 or \code{length(otsu_thresholds)}.
#' @param grow_delta_by_class List or NULL. Optional classwise overrides for \code{DELTA_USED} (same structure as \code{otsu_min_by_class}).
#' @param min_grow_threshold_value Numeric or NULL. Optional global minimum for \code{THR_GROW}. Length 1 or \code{length(otsu_thresholds)}.
#' @param min_grow_threshold_by_class List or NULL. Optional classwise minimum for \code{THR_GROW} (same structure as \code{otsu_min_by_class}).
#'
#' @param min_unit_pixels Integer. Minimum pixels in a unit to attempt unit-level Otsu; otherwise uses global fallback when \code{fallback_to_global=TRUE}.
#' @param min_unit_sample Integer. Minimum sample size for unit-level Otsu; otherwise uses global fallback when \code{fallback_to_global=TRUE}.
#' @param max_global_sample Integer. Maximum sample size for global Otsu estimation.
#' @param max_unit_sample Integer. Maximum per-unit sample size for unit-level Otsu estimation.
#' @param fallback_to_global Logical. If TRUE, global Otsu is used as fallback for unreliable unit estimates. If FALSE, such units are skipped (NA thresholds).
#'
#' @param mask_edge_n_pixels Integer >= 0. Erode the unit mask inward by N pixels before thresholding (reduces border artifacts).
#' @param min_seed_pixels_per_component Integer >= 1. Minimum seed pixels required to keep a candidate component (seed-supported components only).
#'
#' @param grow_engine Character. \code{"whitebox"} or \code{"terra"} for labeling candidate components.
#' @param wbt_diag Logical. WhiteboxTools clump parameter \code{diag} (TRUE=8-neighbor; FALSE=4-neighbor, more conservative).
#' @param wbt_zero_back Logical. WhiteboxTools clump parameter \code{zero_back}.
#'
#' @param tile Logical. If TRUE, polygonization is performed in tiles to reduce memory use.
#' @param n_rows Integer. Number of tile rows (when \code{tile=TRUE}).
#' @param n_cols Integer. Number of tile cols (when \code{tile=TRUE}).
#' @param tile_overlap Numeric. Tile overlap distance in map units (same units as raster CRS).
#'
#' @param index_type Character. Either \code{"RBR"} or \code{"dNBR"} when computing from NBR pre/post. Ignored when \code{raster_path} is provided.
#' @param segment_by_intersection Logical. If TRUE, builds units as CORINE x ecoregion. Required for \code{workflow="otsu_grow"} here.
#' @param corine_classes Numeric/integer vector or NULL. Optional subset of CORINE classes to keep.
#' @param output_format Character. Output vector format: \code{"shp"} or \code{"geojson"}.
#' @param workflow Character. \code{"otsu_grow"} (implemented) or \code{"classic_otsu"} (not implemented here).
#'
#' @param min_otsu_threshold_value Numeric or NULL. Advanced/legacy seed floor (usually redundant). Only affects THR_SEED if more restrictive than LB_GLOBAL and LB_CLASS.
#' @param min_otsu_threshold_by_class List or NULL. Advanced/legacy classwise seed floors (usually redundant). Same structure as \code{otsu_min_by_class}.
#'
#' @param save_debug_rasters Logical. If TRUE, writes intermediate rasters to \code{output_dir/DEBUG}.
#' @param crop_to_units Logical. If TRUE, crops rasters to the extent of valid units (speed-up).
#' @param gdalwarp_path Character. GDAL \code{gdalwarp} executable (used only when \code{reclassify_corine=TRUE} and \code{reproject=TRUE}).
#'
#' @return
#' A named list (one element per run). Each element is a list with:
#' \itemize{
#'   \item \code{out_file}: path to the saved vector output (or \code{NA} if empty),
#'   \item \code{log_file}: path to the per-unit log table (\code{.txt}).
#' }
#'
#' @note
#' Packages required: \pkg{terra}, \pkg{sf}, \pkg{stringr}, \pkg{glue}.
#' Otsu helpers are provided by \pkg{otsuSeg} (or \pkg{OtsuSeg}) in \code{Suggests:}.
#' \pkg{whitebox} is required only when \code{grow_engine="whitebox"}.
#'
#' @examples
#' \dontrun{
#' # Typical run (classwise seeds and classwise growth control)
#' otsu_min_by_class <- list(
#'   corine = list("6"=350,"7"=350,"2"=350,"3"=300,"4"=350,"8"=350)
#' )
#' grow_delta_by_class <- list(
#'   corine = list("6"=50,"7"=50,"2"=50,"3"=50,"4"=50,"8"=50)
#' )
#' min_grow_threshold_by_class <- list(
#'   corine = c("6"=300,"7"=300,"2"=300,"3"=250,"4"=300,"8"=300)
#' )
#'
#' res <- process_otsu_rasters_grow(
#'   raster_path = "RBR_2025_EPSG3035.tif",
#'   output_dir  = "OUT/2025_growth",
#'   year = 2025,
#'   corine_raster_path = "AUX/CORINE_2018.tif",
#'   ecoregion_shapefile_path = "AUX/ECOREGIONS.gpkg",
#'   ecoregion_field = "EnZ_name",
#'   segment_by_intersection = TRUE,
#'   ecoregion_touches = TRUE,
#'   otsu_value_range = c(0, 1500),
#'   trim_percentiles = list(min = 0.05, max = 0.99),
#'   otsu_thresholds = c(300),          # LB_GLOBAL
#'   otsu_min_by_class = otsu_min_by_class,  # LB_CLASS
#'   grow_delta = 100,
#'   grow_delta_by_class = grow_delta_by_class,
#'   min_grow_threshold_value = 200,
#'   min_grow_threshold_by_class = min_grow_threshold_by_class,
#'   mask_edge_n_pixels = 1L,
#'   min_seed_pixels_per_component = 30L,
#'   grow_engine = "whitebox",
#'   wbt_diag = FALSE,
#'   crop_to_units = TRUE,
#'   save_debug_rasters = TRUE
#' )
#' }
#'
#' @importFrom terra rast crs compareGeom nlyr ncell res values global spatSample
#' @importFrom terra crop mask trim resample rasterize vect ifel classify patches focal
#' @importFrom terra as.polygons ext writeRaster as.int freq
#' @importFrom sf st_read st_write st_transform st_crs st_make_valid st_intersection
#' @importFrom sf st_as_sf st_geometry_type st_is_valid st_union st_sf st_sfc
#' @importFrom stringr str_extract
#' @importFrom glue glue
#' @importFrom stats quantile
#' @importFrom graphics hist
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table getExportedValue
#' @export
#'
utils::globalVariables(c(
  "CLUMP_ID",
  "ECO_ID_INTERNAL",
  "UNIT_ID"
))

process_otsu_rasters_grow <- function(
    raster_path = NULL,
    nbr_pre_path = NULL,
    nbr_post_path = NULL,
    output_dir,
    year = NULL,
    otsu_thresholds = c(0, 50, 100),
    trim_percentiles = NULL,
    otsu_value_range = NULL,      # OPTIONAL: range filter (lo,hi) used ONLY for Otsu estimation
    use_original = FALSE,

    corine_raster_path = NULL,
    reclassify_corine = FALSE,
    reclass_matrix = NULL,
    peninsula_shapefile = NULL,
    output_corine_raster_dir = NULL,
    output_corine_vector_dir = NULL,
    reproject = TRUE,
    resolution = 100,

    ecoregion_shapefile_path = NULL,
    ecoregion_field = NULL,
    ecoregion_classes = NULL,
    ecoregion_touches = TRUE,   # NEW: rasterize with touches=TRUE to reduce boundary gaps

    min_otsu_threshold_value = NULL,
    min_otsu_threshold_by_class = NULL,
    otsu_min_by_class = NULL,

    gdalwarp_path = "gdalwarp",

    n_rows = 2,
    n_cols = 3,
    tile_overlap = 1000,
    tile = TRUE,

    index_type = "RBR",
    segment_by_intersection = TRUE,
    corine_classes = NULL,

    output_format = c("shp", "geojson"),

    # -------------------- Otsu + growth controls --------------------
    workflow = c("otsu_grow", "classic_otsu"),
    grow_delta = 50,
    grow_delta_by_class = NULL,            # list(corine=..., ecoregion=..., intersection=...)
    min_grow_threshold_value = NULL,       # optional clamp for thr_grow
    min_grow_threshold_by_class = NULL,    # NEW: list(corine=..., ecoregion=..., intersection=...)

    min_unit_pixels = 10000,
    min_unit_sample = 1000,
    max_global_sample = 2000000,
    max_unit_sample = 200000,
    fallback_to_global = TRUE,

    # -------------------- NEW: growth engine selection --------------------
    grow_engine = c("whitebox", "terra"),
    wbt_diag = TRUE,
    wbt_zero_back = TRUE,

    # -------------------- Border control + seed support (NEW) --------------------
    mask_edge_n_pixels = 1L,               # integer >= 0; erode unit mask inward by N pixels (reduces border artifacts)
    min_seed_pixels_per_component = 15L,   # integer >= 1; minimum seed pixels needed to keep a candidate component

    # Debug / performance
    save_debug_rasters = FALSE,
    crop_to_units = TRUE                   # NEW: crop rasters to unit extent to speed up growth/polygonization
) {

  # -------------------------------------------------------------------------
  # 0) Basic validation / setup
  # -------------------------------------------------------------------------
  output_format <- match.arg(output_format, choices = c("shp", "geojson"))
  workflow <- match.arg(workflow, choices = c("otsu_grow", "classic_otsu"))
  grow_engine <- match.arg(grow_engine, choices = c("whitebox", "terra"))

  tile <- as.logical(tile)[1]
  if (is.na(tile)) stop("Argument 'tile' must be a single logical value (TRUE or FALSE).")
  if (missing(output_dir) || is.null(output_dir)) stop("'output_dir' must be provided.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # NOTE: use_original is not used in this workflow
  if (isTRUE(use_original)) {
    message("[NOTE] use_original=TRUE is not used in this workflow; raster values are always preserved.")
  }

  # ---- Canonical otsu thresholds vector -------------------------------------
  otsu_thresholds_vec <- suppressWarnings(as.numeric(otsu_thresholds))
  if (length(otsu_thresholds_vec) < 1) stop("'otsu_thresholds' must have at least 1 value.")
  if (any(!is.finite(otsu_thresholds_vec))) stop("'otsu_thresholds' must be finite numeric values (no NA/Inf).")
  n_runs <- length(otsu_thresholds_vec)

  # Validate trim_percentiles (optional)
  if (!is.null(trim_percentiles)) {
    if (!is.list(trim_percentiles) || !all(c("min", "max") %in% names(trim_percentiles))) {
      stop("trim_percentiles must be a list with named elements 'min' and 'max' (each in [0,1]).")
    }
    pmin_t <- suppressWarnings(as.numeric(trim_percentiles$min[1]))
    pmax_t <- suppressWarnings(as.numeric(trim_percentiles$max[1]))
    if (!is.finite(pmin_t) || !is.finite(pmax_t) || pmin_t < 0 || pmax_t > 1 || pmin_t >= pmax_t) {
      stop("trim_percentiles$min and $max must be finite in [0,1] with min < max.")
    }
  }

  # Validate optional Otsu value range (applies ONLY to Otsu estimation)
  if (!is.null(otsu_value_range)) {
    if (length(otsu_value_range) != 2 || any(!is.finite(otsu_value_range))) {
      stop("otsu_value_range must be a numeric vector of length 2 (lo, hi).")
    }
    otsu_value_range <- as.numeric(otsu_value_range)
  }

  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
  if (!requireNamespace("glue", quietly = TRUE)) stop("Package 'glue' is required.")

  # Otsu smooth histogram package (Suggests)
  pkg_name <- NULL
  if (requireNamespace("otsuSeg", quietly = TRUE)) {
    pkg_name <- "otsuSeg"
  } else if (requireNamespace("OtsuSeg", quietly = TRUE)) {
    pkg_name <- "OtsuSeg"
  } else {
    stop("Package 'otsuSeg' (or 'OtsuSeg') is required for smooth_histogram() and otsu_threshold_smoothed().")
  }
  smooth_histogram_fun <- getExportedValue(pkg_name, "smooth_histogram")
  otsu_threshold_smoothed_fun <- getExportedValue(pkg_name, "otsu_threshold_smoothed")

  # Whitebox (recommended)
  if (grow_engine == "whitebox") {
    if (!requireNamespace("whitebox", quietly = TRUE)) {
      stop("grow_engine='whitebox' requires package 'whitebox'. Install it (e.g., install.packages('whitebox')).")
    }
    if (!("wbt_clump" %in% getNamespaceExports("whitebox"))) {
      stop("Package 'whitebox' is installed but does not export wbt_clump(). Please update the package.")
    }
  }

  # -------------------------------------------------------------------------
  # Helpers (normalization)
  # -------------------------------------------------------------------------
  normalize_per_run_numeric <- function(x, n, name, allow_null = FALSE) {
    if (is.null(x)) {
      if (isTRUE(allow_null)) return(rep(NA_real_, n))
      stop(sprintf("'%s' must be provided (cannot be NULL).", name))
    }
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 1) return(rep(x, n))
    if (length(x) == n) return(x)
    stop(sprintf("'%s' must be length 1 or %d (same as otsu_thresholds).", name, n))
  }

  normalize_per_run_int <- function(x, n, name, default = NULL) {
    if (is.null(x)) {
      if (!is.null(default)) return(rep(as.integer(default), n))
      stop(sprintf("'%s' must be provided (cannot be NULL).", name))
    }
    x <- suppressWarnings(as.integer(x))
    if (length(x) == 1) return(rep(x, n))
    if (length(x) == n) return(x)
    stop(sprintf("'%s' must be length 1 or %d (same as otsu_thresholds).", name, n))
  }

  grow_delta_vec <- normalize_per_run_numeric(grow_delta, n_runs, "grow_delta", allow_null = FALSE)
  if (any(!is.finite(grow_delta_vec))) stop("'grow_delta' must be finite numeric values (no NA/Inf).")

  min_grow_floor_vec <- normalize_per_run_numeric(min_grow_threshold_value, n_runs, "min_grow_threshold_value", allow_null = TRUE)

  min_seed_vec <- normalize_per_run_int(min_seed_pixels_per_component, n_runs, "min_seed_pixels_per_component", default = 15L)
  min_seed_vec[!is.finite(min_seed_vec)] <- 1L
  min_seed_vec[min_seed_vec < 1L] <- 1L

  mask_edge_n_pixels <- as.integer(mask_edge_n_pixels)[1]
  if (is.na(mask_edge_n_pixels) || mask_edge_n_pixels < 0) stop("mask_edge_n_pixels must be an integer >= 0.")

  # Infer year if not provided
  if (is.null(year)) {
    if (!is.null(raster_path)) {
      raster_base <- tools::file_path_sans_ext(basename(raster_path))
      year_match <- stringr::str_extract(raster_base, "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    } else if (!is.null(output_dir)) {
      year_match <- stringr::str_extract(basename(normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "(19|20)[0-9]{2}")
      year <- if (!is.na(year_match)) year_match else "unknown"
    } else {
      year <- "unknown"
    }
  }

  message(glue::glue(
    "[START] process_otsu_rasters_grow | workflow={workflow} | index_type={index_type} | grow_engine={grow_engine} | tile={tile} | grid={n_rows}x{n_cols} | overlap={tile_overlap} | format={output_format}"
  ))

  # -------------------------------------------------------------------------
  # Helpers
  # -------------------------------------------------------------------------
  sanitize_label <- function(x, max_len = 80) {
    x <- as.character(x)
    if (length(x) > 1) x <- paste(x, collapse = "_")
    x <- x[1]
    x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (nchar(x) > max_len) x <- substr(x, 1, max_len)
    x
  }

  fmt_num_label <- function(x) {
    x <- suppressWarnings(as.numeric(x)[1])
    if (!is.finite(x)) return("NA")

    # string sin notación científica, sin espacios
    s <- formatC(x, format = "fg", digits = 12, flag = "#")
    s <- gsub("\\s+", "", s)

    # quita ceros finales si es decimal (300.000 -> 300)
    s <- sub("\\.?0+$", "", s)

    # seguro para filename: "-" y "." no dan problemas, pero lo estandarizamos:
    s <- gsub("-", "m", s)     # -50 -> m50
    s <- gsub("\\.", "p", s)   # 12.5 -> 12p5
    s
  }

  label_with_suffix <- function(base, suffix, max_len = 80) {
    base <- sanitize_label(base, max_len = max_len)
    suffix <- sanitize_label(suffix, max_len = max_len)
    if (nchar(suffix) >= max_len) return(substr(suffix, 1, max_len))
    avail <- max_len - nchar(suffix)
    if (avail < 1) return(substr(suffix, 1, max_len))
    base2 <- sanitize_label(base, max_len = avail)
    paste0(base2, suffix)
  }

  safe_make_valid_polygon_sf <- function(x) {
    if (is.null(x)) return(NULL)
    x <- sf::st_as_sf(x)
    gt <- sf::st_geometry_type(x)
    x <- x[gt %in% c("POLYGON", "MULTIPOLYGON"), , drop = FALSE]
    if (nrow(x) == 0) return(NULL)
    if (any(!sf::st_is_valid(x))) x <- sf::st_make_valid(x)
    x
  }

  fix_names_for_shp <- function(x) {
    if (is.null(x)) return(NULL)
    names(x) <- make.unique(substr(names(x), 1, 10))
    x
  }

  delete_existing_vector <- function(out_file, output_format) {
    output_format <- match.arg(output_format, choices = c("shp", "geojson"))
    if (output_format == "shp") {
      shp_base <- tools::file_path_sans_ext(out_file)
      for (ext_i in c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".sbn", ".sbx", ".shp.xml")) {
        f <- paste0(shp_base, ext_i)
        if (file.exists(f)) try(file.remove(f), silent = TRUE)
      }
    } else {
      if (file.exists(out_file)) try(file.remove(out_file), silent = TRUE)
    }
  }

  bind_sf_list_safe <- function(lst) {
    lst <- lst[!vapply(lst, is.null, logical(1))]
    lst <- lst[vapply(lst, function(x) inherits(x, "sf") && nrow(x) > 0, logical(1))]
    if (length(lst) == 0) return(NULL)

    all_names <- Reduce(union, lapply(lst, names))
    lst2 <- lapply(lst, function(x) {
      miss <- setdiff(all_names, names(x))
      for (m in miss) x[[m]] <- NA
      x <- x[, all_names, drop = FALSE]
      sf::st_as_sf(x)
    })
    out <- do.call(rbind, lst2)
    sf::st_as_sf(out)
  }

  # Polygonize clump-id raster (one polygon per CLUMP_ID, with tile-safe union by ID)
  polygonize_id_tiles <- function(id_raster, tile = TRUE, n_rows = 2, n_cols = 3, tile_overlap = 0) {
    if (!inherits(id_raster, "SpatRaster")) stop("id_raster must be a terra SpatRaster")

    poly_one <- function(rr) {
      v <- try(terra::as.polygons(rr, dissolve = TRUE, na.rm = TRUE), silent = TRUE)
      if (inherits(v, "try-error") || is.null(v) || nrow(v) == 0) return(NULL)
      sf::st_as_sf(v)
    }

    if (!isTRUE(tile)) {
      out <- poly_one(id_raster)
      out <- safe_make_valid_polygon_sf(out)
      if (is.null(out)) return(sf::st_sf(geometry = sf::st_sfc(crs = 3035))[0, ])
      out <- sf::st_transform(out, 3035)
      id_field <- setdiff(names(out), attr(out, "sf_column"))
      if (length(id_field) > 0) names(out)[names(out) == id_field[1]] <- "CLUMP_ID"
      return(out)
    }

    ext_r <- terra::ext(id_raster)
    tile_width  <- (ext_r[2] - ext_r[1]) / n_cols
    tile_height <- (ext_r[4] - ext_r[3]) / n_rows

    tiles_sf <- list()
    idx <- 1L

    for (i in 0:(n_rows - 1)) {
      for (j in 0:(n_cols - 1)) {

        xmin_tile <- max(ext_r[1] + j * tile_width  - tile_overlap, ext_r[1])
        xmax_tile <- min(ext_r[1] + (j + 1) * tile_width  + tile_overlap, ext_r[2])
        ymin_tile <- max(ext_r[3] + i * tile_height - tile_overlap, ext_r[3])
        ymax_tile <- min(ext_r[3] + (i + 1) * tile_height + tile_overlap, ext_r[4])

        rr <- terra::crop(id_raster, terra::ext(xmin_tile, xmax_tile, ymin_tile, ymax_tile))
        mx <- try(terra::global(rr, "max", na.rm = TRUE)[1, 1], silent = TRUE)
        if (inherits(mx, "try-error") || is.na(mx)) { rm(rr); gc(FALSE); next }

        vv <- poly_one(rr)
        vv <- safe_make_valid_polygon_sf(vv)
        if (!is.null(vv) && nrow(vv) > 0) tiles_sf[[idx]] <- vv

        idx <- idx + 1L
        rm(rr, vv); gc(FALSE)
      }
    }

    if (length(tiles_sf) == 0) return(sf::st_sf(geometry = sf::st_sfc(crs = 3035))[0, ])

    out <- bind_sf_list_safe(tiles_sf)
    out <- safe_make_valid_polygon_sf(out)
    if (is.null(out) || nrow(out) == 0) return(sf::st_sf(geometry = sf::st_sfc(crs = 3035))[0, ])

    out <- sf::st_transform(out, 3035)

    id_field <- setdiff(names(out), attr(out, "sf_column"))
    if (length(id_field) == 0) return(out)
    id_field <- id_field[1]
    names(out)[names(out) == id_field] <- "CLUMP_ID"

    spl <- split(out, out$CLUMP_ID)
    out2 <- lapply(spl, function(g) {
      g1 <- g[1, , drop = FALSE]
      g1$geometry <- sf::st_union(g$geometry)
      g1
    })
    out2 <- do.call(rbind, out2)
    sf::st_as_sf(out2)
  }

  # Erode the unit mask inward by N pixels (border artifact control)
  erode_unit_mask <- function(unit_r, n = 1L) {
    if (!inherits(unit_r, "SpatRaster")) stop("unit_r must be a terra SpatRaster")
    n <- as.integer(n)[1]
    if (is.na(n) || n <= 0) return(unit_r)

    m <- terra::ifel(!is.na(unit_r), 1, 0)
    w <- matrix(1, 3, 3)

    inner <- m
    for (i in seq_len(n)) {
      inner <- terra::focal(inner, w = w, fun = "min", na.rm = FALSE, fillvalue = 0)
    }

    unit_r2 <- terra::mask(unit_r, inner, maskvalues = 0)
    rm(m, inner); gc(FALSE)
    unit_r2
  }

  # terra::freq version-safe (drop NA)
  freq_no_na <- function(x) {
    f <- terra::freq(x)
    if (is.null(f) || length(f) == 0) return(NULL)

    f <- as.data.frame(f)

    if (all(c("value", "count") %in% names(f))) {
      out <- f[, c("value", "count"), drop = FALSE]
    } else if (ncol(f) >= 3) {
      out <- f[, 2:3, drop = FALSE]
      names(out) <- c("value", "count")
    } else if (ncol(f) >= 2) {
      out <- f[, 1:2, drop = FALSE]
      names(out) <- c("value", "count")
    } else {
      return(NULL)
    }

    out <- out[!is.na(out$value), , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    out
  }

  # Otsu threshold on real values, using smoothed histogram on 0..255 scaled values
  compute_otsu_threshold <- function(vals_real, pmin = NULL, pmax = NULL, value_range = NULL) {

    n_input <- length(vals_real)
    vals <- vals_real
    vals <- vals[is.finite(vals)]
    n_finite <- length(vals)

    # Optional hard value range filter (only for Otsu estimation)
    n_range_drop <- 0L
    if (!is.null(value_range) && length(value_range) == 2 && all(is.finite(value_range))) {
      lo <- as.numeric(value_range[1])
      hi <- as.numeric(value_range[2])
      if (lo > hi) { tmp <- lo; lo <- hi; hi <- tmp }
      vals2 <- vals[vals >= lo & vals <= hi]
      n_range_drop <- as.integer(length(vals) - length(vals2))
      vals <- vals2
    }

    if (length(vals) == 0) {
      return(list(ok = FALSE, thr_real = NA_real_, thr_scaled = NA_real_, n_used = 0L,
                  n_input = n_input, n_finite = n_finite, n_range_drop = n_range_drop, n_trim_drop = 0L,
                  min_val = NA_real_, max_val = NA_real_, note = "no_values"))
    }

    # Percentile trimming
    n_trim_drop <- 0L
    if (!is.null(pmin) && !is.null(pmax)) {
      q1 <- stats::quantile(vals, probs = pmin, na.rm = TRUE)
      q2 <- stats::quantile(vals, probs = pmax, na.rm = TRUE)
      vals2 <- vals[vals >= q1 & vals <= q2]
      n_trim_drop <- as.integer(length(vals) - length(vals2))
      vals <- vals2
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) {
        return(list(ok = FALSE, thr_real = NA_real_, thr_scaled = NA_real_, n_used = 0L,
                    n_input = n_input, n_finite = n_finite, n_range_drop = n_range_drop, n_trim_drop = n_trim_drop,
                    min_val = NA_real_, max_val = NA_real_, note = "no_values_after_trim"))
      }
    }

    min_val <- min(vals)
    max_val <- max(vals)
    rng <- max_val - min_val
    if (!is.finite(rng) || rng <= 0) {
      return(list(ok = FALSE, thr_real = NA_real_, thr_scaled = NA_real_, n_used = length(vals),
                  n_input = n_input, n_finite = n_finite, n_range_drop = n_range_drop, n_trim_drop = n_trim_drop,
                  min_val = min_val, max_val = max_val, note = "degenerate_range"))
    }

    vals_scaled <- (vals - min_val) / rng * 255
    vals_scaled[vals_scaled < 0] <- 0
    vals_scaled[vals_scaled > 255] <- 255

    h <- graphics::hist(vals_scaled, breaks = 256, plot = FALSE)
    sm <- smooth_histogram_fun(h$counts)
    sm[is.na(sm)] <- 0

    thr_scaled <- otsu_threshold_smoothed_fun(sm, h$mids)
    thr_real <- (thr_scaled / 255) * rng + min_val

    list(ok = TRUE, thr_real = as.numeric(thr_real), thr_scaled = as.numeric(thr_scaled),
         n_used = length(vals),
         n_input = n_input, n_finite = n_finite, n_range_drop = n_range_drop, n_trim_drop = n_trim_drop,
         min_val = min_val, max_val = max_val, note = "ok")
  }

  # Resolve classwise lower bound (otsu_min)
  resolve_otsu_min_local <- function(corine_class = NULL, ecoregion_name = NULL, cor_eco_name = NULL,
                                     otsu_min_by_class = NULL) {
    if (is.null(otsu_min_by_class)) return(NULL)

    if (!is.null(cor_eco_name) && !is.null(otsu_min_by_class$intersection)) {
      v <- otsu_min_by_class$intersection[[as.character(cor_eco_name)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    if (!is.null(corine_class) && !is.null(otsu_min_by_class$corine)) {
      v <- otsu_min_by_class$corine[[as.character(corine_class)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    if (!is.null(ecoregion_name) && !is.null(otsu_min_by_class$ecoregion)) {
      v <- otsu_min_by_class$ecoregion[[as.character(ecoregion_name)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    NULL
  }

  # Helper: safely pick a scalar numeric from either a named vector OR a named list
  pick_class_value <- function(map, key) {
    if (is.null(map) || is.null(key)) return(NULL)
    key <- as.character(key)[1]

    v <- NULL
    if (is.list(map)) {
      v <- map[[key]]
    } else {
      vv <- map[key]
      if (length(vv) >= 1) v <- vv[[1]]
    }

    if (is.null(v) || length(v) < 1) return(NULL)
    v <- suppressWarnings(as.numeric(v)[1])
    if (!is.finite(v)) return(NULL)
    v
  }

  # Resolve classwise floor (min_otsu_threshold_value) [robust: NULL/list/named-vector]
  resolve_min_floor_local <- function(corine_class = NULL, ecoregion_name = NULL, cor_eco_name = NULL,
                                      min_otsu_threshold_value_default = NULL,
                                      min_otsu_threshold_by_class = NULL) {

    base <- min_otsu_threshold_value_default
    if (!is.null(base)) {
      base <- suppressWarnings(as.numeric(base)[1])
      if (!is.finite(base)) base <- NULL
    }

    if (is.null(min_otsu_threshold_by_class)) return(base)

    take_max <- function(a, b) {
      if (is.null(b) || !is.finite(b)) return(a)
      if (is.null(a) || !is.finite(a)) return(b)
      max(a, b)
    }

    if (!is.null(cor_eco_name) && !is.null(min_otsu_threshold_by_class$intersection)) {
      v <- pick_class_value(min_otsu_threshold_by_class$intersection, cor_eco_name)
      base <- take_max(base, v)
    }
    if (!is.null(corine_class) && !is.null(min_otsu_threshold_by_class$corine)) {
      v <- pick_class_value(min_otsu_threshold_by_class$corine, corine_class)
      base <- take_max(base, v)
    }
    if (!is.null(ecoregion_name) && !is.null(min_otsu_threshold_by_class$ecoregion)) {
      v <- pick_class_value(min_otsu_threshold_by_class$ecoregion, ecoregion_name)
      base <- take_max(base, v)
    }

    base
  }

  # Resolve classwise floor for THR_GROW (optional) [robust: NULL/list/named-vector]
  resolve_min_grow_floor_local <- function(corine_class = NULL, ecoregion_name = NULL, cor_eco_name = NULL,
                                           min_grow_threshold_value_default = NULL,
                                           min_grow_threshold_by_class = NULL) {

    out <- min_grow_threshold_value_default
    if (!is.null(out)) {
      out <- suppressWarnings(as.numeric(out)[1])
      if (!is.finite(out)) out <- NULL
    }

    if (is.null(min_grow_threshold_by_class)) return(out)

    take_max <- function(a, b) {
      if (is.null(b) || !is.finite(b)) return(a)
      if (is.null(a) || !is.finite(a)) return(b)
      max(a, b)
    }

    if (!is.null(cor_eco_name) && !is.null(min_grow_threshold_by_class$intersection)) {
      v <- pick_class_value(min_grow_threshold_by_class$intersection, cor_eco_name)
      out <- take_max(out, v)
    }
    if (!is.null(corine_class) && !is.null(min_grow_threshold_by_class$corine)) {
      v <- pick_class_value(min_grow_threshold_by_class$corine, corine_class)
      out <- take_max(out, v)
    }
    if (!is.null(ecoregion_name) && !is.null(min_grow_threshold_by_class$ecoregion)) {
      v <- pick_class_value(min_grow_threshold_by_class$ecoregion, ecoregion_name)
      out <- take_max(out, v)
    }

    out
  }

  # Resolve delta by class
  resolve_delta_local <- function(corine_class = NULL, ecoregion_name = NULL, cor_eco_name = NULL,
                                  grow_delta_default = 0, grow_delta_by_class = NULL) {
    out <- grow_delta_default
    if (is.null(grow_delta_by_class)) return(out)

    if (!is.null(cor_eco_name) && !is.null(grow_delta_by_class$intersection)) {
      v <- grow_delta_by_class$intersection[[as.character(cor_eco_name)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    if (!is.null(corine_class) && !is.null(grow_delta_by_class$corine)) {
      v <- grow_delta_by_class$corine[[as.character(corine_class)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    if (!is.null(ecoregion_name) && !is.null(grow_delta_by_class$ecoregion)) {
      v <- grow_delta_by_class$ecoregion[[as.character(ecoregion_name)]]
      if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
    }
    out
  }

  # Whitebox clump wrapper (file-based)
  wbt_clump_safe <- function(in_tif, out_tif, diag = TRUE, zero_back = TRUE) {
    out_dir <- dirname(out_tif)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    whitebox::wbt_clump(
      i = in_tif,
      output = out_tif,
      diag = diag,
      zero_back = zero_back
    )
    if (!file.exists(out_tif)) stop("WhiteboxTools clump did not create output: ", out_tif)
    out_tif
  }

  # -------------------------------------------------------------------------
  # 1) Load / compute severity raster (single-band)
  # -------------------------------------------------------------------------
  if (!is.null(nbr_pre_path) && !is.null(nbr_post_path)) {

    if (!file.exists(nbr_pre_path)) stop("Pre-fire NBR raster not found: ", nbr_pre_path)
    if (!file.exists(nbr_post_path)) stop("Post-fire NBR raster not found: ", nbr_post_path)

    nbr_pre  <- terra::rast(nbr_pre_path)
    nbr_post <- terra::rast(nbr_post_path)

    if (!terra::compareGeom(nbr_pre, nbr_post, stopOnError = FALSE)) {
      stop("NBR pre and post rasters must have the same extent, resolution, and CRS.")
    }

    if (tolower(index_type) == "dnbr") {
      r <- (nbr_pre - nbr_post) * 1000
    } else if (tolower(index_type) == "rbr") {
      r <- (nbr_pre - nbr_post) * 1000 / (nbr_pre + 1.001)
    } else {
      stop("`index_type` must be either 'RBR' or 'dNBR'")
    }
    r <- r[[1]]

  } else if (!is.null(raster_path)) {

    if (!file.exists(raster_path)) stop("Raster not found: ", raster_path)
    r <- terra::rast(raster_path)[[1]]

  } else {
    stop("You must provide either 'raster_path' or both 'nbr_pre_path' and 'nbr_post_path'.")
  }

  message(sprintf(
    "[RASTER] loaded | nlyr=%d | ncell=%s | res=%s",
    terra::nlyr(r),
    format(terra::ncell(r), big.mark = ","),
    paste(terra::res(r), collapse = "x")
  ))

  # -------------------------------------------------------------------------
  # 2) Output directories
  # -------------------------------------------------------------------------
  figures_dir <- file.path(output_dir, "FIGURES")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  debug_dir <- file.path(output_dir, "DEBUG")
  if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)

  # -------------------------------------------------------------------------
  # 3) CORINE preparation (optional)
  # -------------------------------------------------------------------------
  corine_masked <- NULL

  if (!is.null(corine_raster_path)) {

    if (!file.exists(corine_raster_path)) stop("CORINE raster not found: ", corine_raster_path)

    corine_raster <- terra::rast(corine_raster_path)
    corine_raster_reclassed <- corine_raster

    message(sprintf("[CORINE] loaded | path=%s", corine_raster_path))

    if (reclassify_corine) {
      if (is.null(peninsula_shapefile)) stop("reclassify_corine=TRUE but peninsula_shapefile is NULL.")
      if (is.null(output_corine_raster_dir)) stop("output_corine_raster_dir must be provided when reclassify_corine=TRUE.")
      if (is.null(reclass_matrix)) stop("reclass_matrix must be provided when reclassify_corine=TRUE.")
      if (!dir.exists(output_corine_raster_dir)) dir.create(output_corine_raster_dir, recursive = TRUE)

      peninsula <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula <- sf::st_transform(peninsula, sf::st_crs(terra::crs(corine_raster, proj = TRUE)))
      peninsula <- sf::st_make_valid(peninsula)

      cropped <- terra::crop(corine_raster, terra::vect(peninsula))
      masked  <- terra::mask(cropped, terra::vect(peninsula))

      reclass_matrix <- as.matrix(reclass_matrix)

      if (ncol(reclass_matrix) == 2) {
        class_r <- terra::classify(masked, rcl = reclass_matrix, others = NA)
      } else if (ncol(reclass_matrix) == 3) {
        class_r <- terra::classify(masked, rcl = reclass_matrix, include.lowest = TRUE, right = FALSE)
      } else {
        stop("reclass_matrix must have 2 columns (value,new_value) or 3 columns (from,to,becomes).")
      }

      input_name <- tools::file_path_sans_ext(basename(corine_raster_path))
      out_filename <- paste0("reclassified_", input_name, "_", format(Sys.Date(), "%Y%j"), ".tif")
      out_reproj <- file.path(output_corine_raster_dir, out_filename)

      if (reproject) {
        tmp_unproj <- tempfile(fileext = ".tif")
        terra::writeRaster(class_r, tmp_unproj, overwrite = TRUE)
        if (file.exists(out_reproj)) file.remove(out_reproj)

        cmd <- glue::glue('"{gdalwarp_path}" -t_srs "EPSG:3035" -tr {resolution} {resolution} -r near -co "COMPRESS=LZW" "{tmp_unproj}" "{out_reproj}"')
        system(cmd)
        if (!file.exists(out_reproj)) stop("gdalwarp did not create: ", out_reproj)

      } else {
        terra::writeRaster(class_r, out_reproj, overwrite = TRUE, datatype = "INT2U", NAflag = 0, gdal = c("COMPRESS=LZW"))
        if (!file.exists(out_reproj)) stop("writeRaster failed: ", out_reproj)
      }

      corine_raster_reclassed <- terra::rast(out_reproj)
    }

    # Align CORINE to r
    corine_resampled <- terra::resample(corine_raster_reclassed, r, method = "near")
    corine_masked <- terra::mask(corine_resampled, r)

    if (!is.null(corine_classes)) {
      corine_masked <- terra::ifel(corine_masked %in% corine_classes, corine_masked, NA)
    }

    cor_n <- try(terra::global(terra::ifel(!is.na(corine_masked), 1, 0), "sum", na.rm = TRUE)[1, 1], silent = TRUE)
    message(sprintf("[CORINE] valid cells after mask/filter=%s",
                    if (!inherits(cor_n, "try-error") && is.finite(cor_n)) format(cor_n, big.mark = ",") else "NA"))

    if (all(is.na(terra::values(corine_masked)))) stop("All CORINE values removed after filtering/masking.")
  }

  # -------------------------------------------------------------------------
  # 4) Ecoregions preparation (optional)
  # -------------------------------------------------------------------------
  ecoregions_sf <- NULL
  ecoregion_classes_internal <- NULL

  if (!is.null(ecoregion_shapefile_path)) {

    if (!file.exists(ecoregion_shapefile_path)) stop("Ecoregion shapefile not found: ", ecoregion_shapefile_path)

    ecoregions_sf <- sf::st_read(ecoregion_shapefile_path, quiet = TRUE)

    if (!is.null(peninsula_shapefile) && file.exists(peninsula_shapefile)) {
      peninsula_sf <- sf::st_read(peninsula_shapefile, quiet = TRUE)
      peninsula_sf <- sf::st_transform(peninsula_sf, sf::st_crs(ecoregions_sf))
      peninsula_sf <- sf::st_make_valid(peninsula_sf)
      ecoregions_sf <- suppressWarnings(sf::st_intersection(ecoregions_sf, peninsula_sf))
    }

    if (nrow(ecoregions_sf) == 0) stop("Ecoregions shapefile has no geometries after clipping.")
    ecoregions_sf <- sf::st_make_valid(ecoregions_sf)

    if (is.null(ecoregion_field)) {
      if ("EnS_name" %in% names(ecoregions_sf)) {
        ecoregion_field <- "EnS_name"
      } else if ("EnZ_name" %in% names(ecoregions_sf)) {
        ecoregion_field <- "EnZ_name"
      } else {
        stop("You must specify ecoregion_field. No default field found.")
      }
    }
    if (!(ecoregion_field %in% names(ecoregions_sf))) {
      stop("Field '", ecoregion_field, "' not found in ecoregion shapefile.")
    }

    r_crs <- sf::st_crs(terra::crs(r, proj = TRUE))
    if (!is.na(r_crs)) ecoregions_sf <- sf::st_transform(ecoregions_sf, r_crs)

    if (is.null(ecoregion_classes)) {
      ecoregion_classes_internal <- unique(ecoregions_sf[[ecoregion_field]])
    } else {
      ecoregion_classes_internal <- ecoregion_classes
      ecoregions_sf <- ecoregions_sf[ecoregions_sf[[ecoregion_field]] %in% ecoregion_classes_internal, , drop = FALSE]
    }

    message(sprintf("[ECOREG] classes_n=%d", length(ecoregion_classes_internal)))
  }

  # -------------------------------------------------------------------------
  # 5) STRATIFIED OTSU + GROWTH
  # -------------------------------------------------------------------------
  if (workflow != "otsu_grow") stop("workflow='classic_otsu' is not implemented here. Use workflow='otsu_grow'.")

  if (!isTRUE(segment_by_intersection)) {
    stop("workflow='otsu_grow' requires segment_by_intersection=TRUE (CORINE x ECOREGION).")
  }
  if (is.null(corine_masked) || is.null(ecoregions_sf)) {
    stop("workflow='otsu_grow' requires BOTH corine_raster_path and ecoregion_shapefile_path.")
  }

  message("[WORKFLOW] otsu_grow: building CORINE x ECOREGION unit raster...")

  eco_factor <- as.factor(ecoregions_sf[[ecoregion_field]])
  eco_levels <- levels(eco_factor)
  ecoregions_sf$ECO_ID_INTERNAL <- as.integer(eco_factor)

  eco_r <- terra::rasterize(
    terra::vect(ecoregions_sf),
    terra::rast(corine_masked),
    field = "ECO_ID_INTERNAL",
    touches = isTRUE(ecoregion_touches)
  )
  eco_r <- terra::mask(eco_r, corine_masked)

  eco_n <- try(terra::global(terra::ifel(!is.na(eco_r), 1, 0), "sum", na.rm = TRUE)[1, 1], silent = TRUE)
  cor_n <- try(terra::global(terra::ifel(!is.na(corine_masked), 1, 0), "sum", na.rm = TRUE)[1, 1], silent = TRUE)
  if (!inherits(eco_n, "try-error") && !inherits(cor_n, "try-error") && is.finite(eco_n) && is.finite(cor_n) && cor_n > 0) {
    cov <- 100 * (eco_n / cor_n)
    message(sprintf("[ECOREG] raster coverage over CORINE = %.2f%%", cov))
    if (cov < 95) {
      warning(sprintf(
        "Low ecoregion coverage (%.2f%%). If full coverage is expected, check CRS, clipping, filtering, geometry validity. Consider ecoregion_touches=TRUE and avoid unintended clipping.",
        cov
      ))
    }
  }

  unit_r <- corine_masked * 10000 + eco_r
  unit_r <- terra::ifel(!is.na(unit_r), unit_r, NA)

  if (isTRUE(crop_to_units)) {
    unit_trim <- terra::trim(unit_r)
    if (!is.null(unit_trim)) {
      r <- terra::crop(r, unit_trim)
      corine_masked <- terra::crop(corine_masked, unit_trim)
      eco_r <- terra::crop(eco_r, unit_trim)
      unit_r <- unit_trim
      rm(unit_trim); gc(FALSE)
      message("[CROP] cropped rasters to unit extent (speed-up enabled).")
    }
  }

  # Remove N-pixel border inside unit mask
  if (mask_edge_n_pixels > 0) {
    message(sprintf("[BORDER] eroding unit mask inward by %d pixel(s) (each pixel = %.0f m)",
                    mask_edge_n_pixels, terra::res(unit_r)[1]))

    unit_r <- erode_unit_mask(unit_r, n = mask_edge_n_pixels)

    r            <- terra::mask(r, unit_r)
    corine_masked <- terra::mask(corine_masked, unit_r)
    eco_r        <- terra::mask(eco_r, unit_r)
  }

  unit_freq <- freq_no_na(unit_r)
  if (is.null(unit_freq) || nrow(unit_freq) == 0) stop("No CORINE x ECOREGION units could be built (unit_freq empty).")
  names(unit_freq) <- c("UNIT_ID", "N_PIX")
  unit_freq$UNIT_ID <- as.integer(unit_freq$UNIT_ID)
  unit_freq$N_PIX <- as.numeric(unit_freq$N_PIX)

  message(sprintf("[UNITS] n_units=%d | total_valid_pixels=%s", nrow(unit_freq), format(sum(unit_freq$N_PIX), big.mark = ",")))

  message(sprintf("[SAMPLE] drawing global sample up to %s points...", format(max_global_sample, big.mark = ",")))
  s <- c(r, unit_r)
  names(s) <- c("RBR", "UNIT_ID")

  sample_df <- try(terra::spatSample(
    s,
    size = min(max_global_sample, sum(unit_freq$N_PIX)),
    method = "random",
    na.rm = TRUE,
    as.df = TRUE
  ), silent = TRUE)

  sample_ok <- !(inherits(sample_df, "try-error") || is.null(sample_df) || nrow(sample_df) == 0)
  if (!sample_ok) message("[SAMPLE][WARN] terra::spatSample failed -> per-unit fallback sampling will be used.")

  if (sample_ok) {
    if (!all(c("RBR", "UNIT_ID") %in% names(sample_df))) names(sample_df)[1:2] <- c("RBR", "UNIT_ID")
    sample_df$UNIT_ID <- as.integer(sample_df$UNIT_ID)
    sample_df$RBR <- as.numeric(sample_df$RBR)
    sample_df <- sample_df[is.finite(sample_df$RBR) & is.finite(sample_df$UNIT_ID), , drop = FALSE]
  }

  global_vals <- if (sample_ok) sample_df$RBR else as.numeric(terra::values(r))
  if (!is.null(trim_percentiles) && all(c("min", "max") %in% names(trim_percentiles))) {
    pmin_g <- trim_percentiles$min[1]; pmax_g <- trim_percentiles$max[1]
  } else { pmin_g <- NULL; pmax_g <- NULL }

  otsu_global <- compute_otsu_threshold(global_vals, pmin = pmin_g, pmax = pmax_g, value_range = otsu_value_range)
  if (!isTRUE(otsu_global$ok) || !is.finite(otsu_global$thr_real)) stop("Failed to compute global Otsu threshold.")
  thr_global_otsu <- otsu_global$thr_real
  message(sprintf("[OTSU_GLOBAL] thr=%.3f | n_used=%s | n_input=%s | drop_range=%s | drop_trim=%s",
                  thr_global_otsu,
                  format(otsu_global$n_used, big.mark = ","),
                  format(otsu_global$n_input, big.mark = ","),
                  format(otsu_global$n_range_drop, big.mark = ","),
                  format(otsu_global$n_trim_drop, big.mark = ",")))

  samples_by_unit <- if (sample_ok) split(sample_df$RBR, sample_df$UNIT_ID) else NULL

  unit_table <- unit_freq
  unit_table$CORINE_CLASS <- as.integer(floor(unit_table$UNIT_ID / 10000))
  unit_table$ECO_ID_INTERNAL <- as.integer(unit_table$UNIT_ID - unit_table$CORINE_CLASS * 10000)
  unit_table$ECO_CLASS <- ifelse(
    unit_table$ECO_ID_INTERNAL >= 1 & unit_table$ECO_ID_INTERNAL <= length(eco_levels),
    eco_levels[unit_table$ECO_ID_INTERNAL],
    NA_character_
  )
  unit_table$COR_ECO_LABEL <- paste0(unit_table$CORINE_CLASS, "_", gsub("[^A-Za-z0-9]", "", unit_table$ECO_CLASS))

  unit_table$OTSU_RAW <- NA_real_
  unit_table$OTSU_NOTE <- NA_character_
  unit_table$N_SAMPLE <- 0L

  unit_table$OTSU_N_INPUT <- 0L
  unit_table$OTSU_N_FINITE <- 0L
  unit_table$OTSU_N_RANGE_DROP <- 0L
  unit_table$OTSU_N_TRIM_DROP <- 0L

  if (!is.null(trim_percentiles) && all(c("min", "max") %in% names(trim_percentiles))) {
    pmin_u <- trim_percentiles$min[1]; pmax_u <- trim_percentiles$max[1]
  } else { pmin_u <- NULL; pmax_u <- NULL }

  message("[OTSU_UNITS] computing per-unit Otsu thresholds...")

  for (k in seq_len(nrow(unit_table))) {

    uid <- unit_table$UNIT_ID[k]
    n_pix <- unit_table$N_PIX[k]

    if (!is.finite(n_pix) || n_pix < min_unit_pixels) {
      if (isTRUE(fallback_to_global)) {
        unit_table$OTSU_RAW[k] <- thr_global_otsu
        unit_table$OTSU_NOTE[k] <- "fallback_small_unit"
      } else {
        unit_table$OTSU_RAW[k] <- NA_real_
        unit_table$OTSU_NOTE[k] <- "skip_small_unit"
      }
      next
    }

    vals_u <- NULL
    if (sample_ok) {
      vals_u <- samples_by_unit[[as.character(uid)]]
      if (!is.null(vals_u)) {
        vals_u <- vals_u[is.finite(vals_u)]
        if (length(vals_u) > max_unit_sample) vals_u <- sample(vals_u, size = max_unit_sample)
        unit_table$N_SAMPLE[k] <- length(vals_u)
      }
    }

    if (is.null(vals_u) || length(vals_u) < min_unit_sample) {
      r_u <- terra::ifel(unit_r == uid, r, NA)
      vals_all <- terra::values(r_u)
      vals_all <- vals_all[is.finite(vals_all)]
      if (length(vals_all) > max_unit_sample) vals_all <- sample(vals_all, size = max_unit_sample)
      vals_u <- vals_all
      unit_table$N_SAMPLE[k] <- length(vals_u)
    }

    if (length(vals_u) < min_unit_sample) {
      if (isTRUE(fallback_to_global)) {
        unit_table$OTSU_RAW[k] <- thr_global_otsu
        unit_table$OTSU_NOTE[k] <- "fallback_low_sample"
      } else {
        unit_table$OTSU_RAW[k] <- NA_real_
        unit_table$OTSU_NOTE[k] <- "skip_low_sample"
      }
      next
    }

    o <- compute_otsu_threshold(vals_u, pmin = pmin_u, pmax = pmax_u, value_range = otsu_value_range)

    unit_table$OTSU_N_INPUT[k] <- as.integer(o$n_input)
    unit_table$OTSU_N_FINITE[k] <- as.integer(o$n_finite)
    unit_table$OTSU_N_RANGE_DROP[k] <- as.integer(o$n_range_drop)
    unit_table$OTSU_N_TRIM_DROP[k] <- as.integer(o$n_trim_drop)

    if (!isTRUE(o$ok) || !is.finite(o$thr_real)) {
      if (isTRUE(fallback_to_global)) {
        unit_table$OTSU_RAW[k] <- thr_global_otsu
        unit_table$OTSU_NOTE[k] <- paste0("fallback_otsu_fail_", o$note)
      } else {
        unit_table$OTSU_RAW[k] <- NA_real_
        unit_table$OTSU_NOTE[k] <- paste0("skip_otsu_fail_", o$note)
      }
    } else {
      unit_table$OTSU_RAW[k] <- o$thr_real
      unit_table$OTSU_NOTE[k] <- "ok"
    }
  }

  message(sprintf("[OTSU_UNITS] done | ok_units=%d | fallback/skip_units=%d",
                  sum(unit_table$OTSU_NOTE == "ok", na.rm = TRUE),
                  sum(unit_table$OTSU_NOTE != "ok", na.rm = TRUE)))

  # -------------------------------------------------------------------------
  # PRECOMPUTE BASE RUN LABELS + DUPLICATE FLAGS (prevents overwriting outputs)
  # Output label MUST include: otsu_threshold (run LB), grow_delta, min_seed_pixels
  # -------------------------------------------------------------------------
  base_labels_raw <- vapply(seq_len(n_runs), function(i) {
    paste0(
      "otsu",  fmt_num_label(otsu_thresholds_vec[i]),
      "_d",    fmt_num_label(grow_delta_vec[i]),
      "_seed", as.integer(min_seed_vec[i]),
      "_",     grow_engine
    )
  }, character(1))

  base_labels <- vapply(base_labels_raw, sanitize_label, character(1), max_len = 80)

  dup_flag <- duplicated(base_labels) | duplicated(base_labels, fromLast = TRUE)
  if (any(dup_flag)) {
    message(sprintf(
      "[RUNLABEL] detected %d duplicated run label(s); adding _runiXX suffix to avoid overwriting.",
      sum(dup_flag)
    ))
  }



  results <- list()

  for (run_i in seq_len(n_runs)) {

    run_lb <- as.numeric(otsu_thresholds_vec[run_i])
    run_delta    <- as.numeric(grow_delta_vec[run_i])
    run_min_grow <- as.numeric(min_grow_floor_vec[run_i])
    run_min_seed <- as.integer(min_seed_vec[run_i])

    run_label_base <- base_labels[run_i]
    if (isTRUE(dup_flag[run_i])) {
      run_label <- label_with_suffix(
        base = run_label_base,
        suffix = paste0("_runi", sprintf("%02d", run_i)),
        max_len = 80
      )
    } else {
      run_label <- run_label_base
    }

    message(sprintf(
      "[RUN] lower_bound=%s | delta_default=%s | min_seed_pixels=%s | min_grow_floor=%s",
      as.character(run_lb),
      as.character(run_delta),
      as.character(run_min_seed),
      ifelse(is.finite(run_min_grow), as.character(run_min_grow), "NA")
    ))

    unit_run <- unit_table
    unit_run$LB_GLOBAL <- run_lb
    unit_run$LB_CLASS <- NA_real_
    unit_run$FLOOR_USED <- NA_real_
    unit_run$DELTA_USED <- NA_real_
    unit_run$THR_SEED <- NA_real_
    unit_run$THR_GROW <- NA_real_
    unit_run$LB_APPLIED <- FALSE
    unit_run$FLOOR_APPLIED <- FALSE
    unit_run$GROW_CLAMPED <- FALSE
    unit_run$EXCLUDED_NO_THR <- FALSE
    unit_run$GROW_FLOOR_USED <- NA_real_
    unit_run$GROW_FLOOR_APPLIED <- FALSE

    for (k in seq_len(nrow(unit_run))) {

      cls <- unit_run$CORINE_CLASS[k]
      eco <- unit_run$ECO_CLASS[k]
      ce  <- unit_run$COR_ECO_LABEL[k]

      thr <- suppressWarnings(as.numeric(unit_run$OTSU_RAW[k])[1])

      lb_class <- resolve_otsu_min_local(cls, eco, ce, otsu_min_by_class)
      lb_class <- suppressWarnings(as.numeric(lb_class)[1])
      if (!is.finite(lb_class)) lb_class <- NULL
      unit_run$LB_CLASS[k] <- if (!is.null(lb_class)) lb_class else NA_real_

      floor_local <- resolve_min_floor_local(cls, eco, ce, min_otsu_threshold_value, min_otsu_threshold_by_class)
      floor_local <- suppressWarnings(as.numeric(floor_local)[1])
      if (!is.finite(floor_local)) floor_local <- NULL
      unit_run$FLOOR_USED[k] <- if (!is.null(floor_local)) floor_local else NA_real_

      delta_local <- resolve_delta_local(cls, eco, ce, run_delta, grow_delta_by_class)
      delta_local <- suppressWarnings(as.numeric(delta_local)[1])
      if (!is.finite(delta_local)) delta_local <- run_delta
      unit_run$DELTA_USED[k] <- delta_local

      if (!is.finite(thr)) {
        if (isTRUE(fallback_to_global)) {
          thr <- thr_global_otsu
        } else {
          unit_run$THR_SEED[k] <- NA_real_
          unit_run$THR_GROW[k] <- NA_real_
          unit_run$EXCLUDED_NO_THR[k] <- TRUE
          next
        }
      }

      lb_final <- run_lb
      if (!is.null(lb_class) && is.finite(lb_class)) lb_final <- max(lb_final, lb_class)

      thr_seed <- thr
      lb_applied <- FALSE
      floor_applied <- FALSE

      if (is.finite(lb_final) && thr_seed < lb_final) {
        thr_seed <- lb_final
        lb_applied <- TRUE
      }

      if (!is.null(floor_local) && is.finite(floor_local) && thr_seed < floor_local) {
        thr_seed <- floor_local
        floor_applied <- TRUE
      }

      thr_grow <- thr_seed - delta_local

      grow_floor_local <- resolve_min_grow_floor_local(
        corine_class = cls,
        ecoregion_name = eco,
        cor_eco_name = ce,
        min_grow_threshold_value_default = ifelse(is.finite(run_min_grow), run_min_grow, NULL),
        min_grow_threshold_by_class = min_grow_threshold_by_class
      )
      grow_floor_local <- suppressWarnings(as.numeric(grow_floor_local)[1])
      if (!is.finite(grow_floor_local)) grow_floor_local <- NULL

      unit_run$GROW_FLOOR_USED[k] <- if (!is.null(grow_floor_local)) grow_floor_local else NA_real_

      grow_clamped <- FALSE
      grow_floor_applied <- FALSE

      if (!is.null(grow_floor_local) && thr_grow < grow_floor_local) {
        thr_grow <- grow_floor_local
        grow_clamped <- TRUE
        grow_floor_applied <- TRUE
      }

      if (is.finite(thr_seed) && is.finite(thr_grow) && thr_grow > thr_seed) {
        thr_grow <- thr_seed
        grow_clamped <- TRUE
      }

      unit_run$THR_SEED[k] <- thr_seed
      unit_run$THR_GROW[k] <- thr_grow

      unit_run$LB_APPLIED[k] <- lb_applied
      unit_run$FLOOR_APPLIED[k] <- floor_applied
      unit_run$GROW_CLAMPED[k] <- grow_clamped
      unit_run$GROW_FLOOR_APPLIED[k] <- grow_floor_applied
    }

    # ------------------------------------------------------------------
    # Per-unit thresholds as rasters (filter out NA thresholds)
    # ------------------------------------------------------------------
    rcl_seed <- unit_run[is.finite(unit_run$THR_SEED), c("UNIT_ID", "THR_SEED"), drop = FALSE]
    if (nrow(rcl_seed) == 0) {
      thr_seed_r <- terra::ifel(!is.na(unit_r), NA_real_, NA_real_)
    } else {
      thr_seed_r <- terra::classify(unit_r, rcl = as.matrix(rcl_seed), others = NA)
    }

    rcl_grow <- unit_run[is.finite(unit_run$THR_GROW), c("UNIT_ID", "THR_GROW"), drop = FALSE]
    if (nrow(rcl_grow) == 0) {
      thr_grow_r <- terra::ifel(!is.na(unit_r), NA_real_, NA_real_)
    } else {
      thr_grow_r <- terra::classify(unit_r, rcl = as.matrix(rcl_grow), others = NA)
    }

    # ------------------------------------------------------------------
    # Seed mask (per-unit THR_SEED)
    # ------------------------------------------------------------------
    seed_mask <- terra::ifel(r >= thr_seed_r, 1, NA)

    # ------------------------------------------------------------------
    # Candidate mask (per-unit THR_GROW), clump on binary
    # ------------------------------------------------------------------
    cand_uid <- terra::ifel(r >= thr_grow_r, unit_r, NA)
    cand_bin  <- terra::ifel(!is.na(cand_uid), 1, NA)   # 1/NA (used by terra patches)
    cand_bin0 <- terra::ifel(is.na(cand_bin), 0, 1)     # 0/1 (used by WBT clump)
    cand_bin0 <- terra::as.int(cand_bin0)
    cand_bin0 <- terra::ifel(is.na(unit_r), 0, cand_bin0)

    message("[GROW] labeling candidate components on binary mask...")

    if (grow_engine == "terra") {

      # IMPORTANT FIX:
      # terra::patches() groups non-NA cells. Using 0/1 would clump the entire 0 background.
      cl_in <- cand_bin  # 1/NA
      cl <- terra::patches(cl_in, directions = 8)
      cl <- terra::ifel(cl == 0, NA, cl)
      rm(cl_in); gc(FALSE)

    } else {

      # IMPORTANT FIX:
      # Do not always write debug rasters. Use temp files unless save_debug_rasters=TRUE.
      if (isTRUE(save_debug_rasters)) {
        cand_file <- file.path(debug_dir, paste0("cand_bin0_", year, "_", run_label, ".tif"))
        cl_file   <- file.path(debug_dir, paste0("clumps_", year, "_", run_label, ".tif"))
      } else {
        cand_file <- tempfile(fileext = ".tif")
        cl_file   <- tempfile(fileext = ".tif")
      }

      terra::writeRaster(
        cand_bin0, cand_file,
        overwrite = TRUE,
        datatype  = "INT1U",
        gdal      = c("COMPRESS=LZW")
      )

      wbt_clump_safe(
        cand_file, cl_file,
        diag      = isTRUE(wbt_diag),
        zero_back = isTRUE(wbt_zero_back)
      )

      cl <- terra::rast(cl_file)
      cl <- terra::ifel(cl == 0, NA, cl)

      if (!isTRUE(save_debug_rasters)) {
        try(file.remove(cand_file), silent = TRUE)
        try(file.remove(cl_file), silent = TRUE)
      }
    }

    # ------------------------------------------------------------------
    # Keep only components that contain seed pixels (and enough seed support)
    # ------------------------------------------------------------------
    run_min_seed <- as.integer(run_min_seed)[1]
    if (is.na(run_min_seed) || run_min_seed < 1L) run_min_seed <- 1L

    cl_seed <- terra::mask(cl, seed_mask)
    fq <- freq_no_na(cl_seed)

    if (!is.null(fq) && nrow(fq) > 0) {
      fq_before <- nrow(fq)
      if (!all(c("value", "count") %in% names(fq))) {
        message("[GROW] freq table missing expected columns -> empty output.")
        fq <- NULL
      } else {
        fq <- fq[is.finite(fq$value) & is.finite(fq$count), , drop = FALSE]
        fq <- fq[fq$count >= run_min_seed, , drop = FALSE]

        message(sprintf(
          "[GROW] components with seeds: %d | kept after min_seed_pixels (%d): %d",
          fq_before, run_min_seed, ifelse(is.null(fq), 0, nrow(fq))
        ))
      }
    }

    if (is.null(fq) || nrow(fq) == 0) {
      message("[GROW] no seed-supported components found -> empty output.")
      keep_id_r <- terra::ifel(!is.na(seed_mask), NA, NA)
      final_bin <- keep_id_r
    } else {

      ids <- as.integer(fq$value)
      ids <- ids[is.finite(ids)]
      ids <- unique(ids)

      if (length(ids) == 0) {
        message("[GROW] seed-supported ids empty after filtering -> empty output.")
        keep_id_r <- terra::ifel(!is.na(seed_mask), NA, NA)
        final_bin <- keep_id_r
      } else {
        keep_id_r <- terra::classify(cl, rcl = cbind(ids, ids), others = NA)
        final_bin <- terra::ifel(!is.na(keep_id_r), 1, NA)
      }
    }

    # ------------------------------------------------------------------
    # Debug outputs (optional)
    # ------------------------------------------------------------------
    if (isTRUE(save_debug_rasters)) {

      seed_pix <- try(
        terra::global(terra::ifel(!is.na(seed_mask), 1, 0), "sum", na.rm = TRUE)[1, 1],
        silent = TRUE
      )
      cand_pix <- try(
        terra::global(cand_bin0, "sum", na.rm = TRUE)[1, 1],
        silent = TRUE
      )
      cl_max <- try(
        terra::global(cl, "max", na.rm = TRUE)[1, 1],
        silent = TRUE
      )

      message(sprintf(
        "[DBG] seed_pix=%s | cand_pix=%s | clumps_max=%s",
        if (!inherits(seed_pix, "try-error") && is.finite(seed_pix)) format(seed_pix, big.mark = ",") else "NA",
        if (!inherits(cand_pix, "try-error") && is.finite(cand_pix)) format(cand_pix, big.mark = ",") else "NA",
        if (!inherits(cl_max, "try-error") && is.finite(cl_max)) as.character(cl_max) else "NA"
      ))

      terra::writeRaster(seed_mask, file.path(debug_dir, paste0("seed_mask_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      terra::writeRaster(cand_uid,  file.path(debug_dir, paste0("cand_uid_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      terra::writeRaster(cand_bin0, file.path(debug_dir, paste0("cand_bin0_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      terra::writeRaster(cl,        file.path(debug_dir, paste0("clumps_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      terra::writeRaster(keep_id_r, file.path(debug_dir, paste0("keep_id_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      terra::writeRaster(final_bin, file.path(debug_dir, paste0("binary_final_", year, "_", run_label, ".tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    }

    # ------------------------------------------------------------------
    # Polygonize kept components BY clump ID (one polygon per component)
    # ------------------------------------------------------------------
    message("[VECTOR] polygonizing kept components (by clump ID)...")
    polys <- polygonize_id_tiles(
      keep_id_r,
      tile = tile, n_rows = n_rows, n_cols = n_cols, tile_overlap = tile_overlap
    )

    if (!is.null(polys) && inherits(polys, "sf") && nrow(polys) > 0) {
      polys$YEAR <- as.character(year)
      polys$RUN <- as.character(run_label)
      polys$LB <- as.numeric(run_lb)
      polys$DELTA <- as.numeric(run_delta)
      polys$MINSEED <- as.integer(run_min_seed)
      polys$INDEX <- as.character(index_type)
      polys$WORKFLOW <- "otsu_grow"
      polys$ENGINE <- as.character(grow_engine)
    }

    ext <- if (output_format == "geojson") ".geojson" else ".shp"
    out_file <- file.path(output_dir, paste0("BA_", year, "_OTSUGROW_CORI_ECOREG_", run_label, ext))
    log_file <- file.path(output_dir, paste0("BA_", year, "_OTSUGROW_CORI_ECOREG_", run_label, "_log.txt"))

    if (!is.null(polys) && inherits(polys, "sf") && nrow(polys) > 0) {
      polys_out <- if (output_format == "shp") fix_names_for_shp(polys) else polys
      delete_existing_vector(out_file, output_format)
      sf::st_write(polys_out, out_file, append = FALSE, quiet = TRUE)
      message(sprintf("[SAVE] %s", out_file))
    } else {
      message("[SAVE] empty polygons -> no vector written.")
      out_file <- NA_character_
    }

    log_df <- unit_run
    log_df$YEAR <- as.character(year)
    log_df$RUN <- as.character(run_label)
    log_df$MIN_SEED_PIX <- as.integer(run_min_seed)
    log_df$INDEX <- as.character(index_type)
    log_df$ENGINE <- as.character(grow_engine)

    keep_cols <- c(
      "YEAR","RUN","INDEX","ENGINE","MIN_SEED_PIX","UNIT_ID","CORINE_CLASS","ECO_ID_INTERNAL","ECO_CLASS","COR_ECO_LABEL",
      "N_PIX","N_SAMPLE","OTSU_RAW","OTSU_NOTE","OTSU_N_INPUT","OTSU_N_FINITE","OTSU_N_RANGE_DROP","OTSU_N_TRIM_DROP",
      "LB_GLOBAL","LB_CLASS","LB_APPLIED",
      "FLOOR_USED","FLOOR_APPLIED",
      "DELTA_USED","THR_SEED","THR_GROW","GROW_CLAMPED","GROW_FLOOR_USED","GROW_FLOOR_APPLIED","EXCLUDED_NO_THR"
    )
    keep_cols <- keep_cols[keep_cols %in% names(log_df)]
    log_df <- log_df[, keep_cols, drop = FALSE]
    utils::write.table(log_df, file = log_file, row.names = FALSE, sep = "\t", quote = FALSE)
    message(sprintf("[LOG] %s", log_file))

    results[[run_label]] <- list(out_file = out_file, log_file = log_file)

    rm(thr_seed_r, thr_grow_r, seed_mask, cand_uid, cand_bin, cand_bin0, cl, cl_seed, fq, keep_id_r, final_bin, polys)
    gc(FALSE)
  }

  return(results)
}
