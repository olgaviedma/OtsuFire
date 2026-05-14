# ============================================================
# prepare_fire_polys_for_validation()
# Super-complete + robust geometry handling + customizable output names
# ============================================================

#' Prepare burned-area polygons (reference and detected) for validation
#'
#' Standardizes two polygon datasets (reference and detected) for validation:
#' optionally filters by DOY/attributes, clips to AOI and (optionally) to a
#' burnable mask raster, and optionally removes detections that intersect any
#' reference polygons dropped by ANY stage (SAFE erase rule).
#'
#' Key outputs are controlled by:
#' - output_prefix: prepended to filenames (sanitized)
#' - ref_output_tag / det_output_tag: tag inside filenames (sanitized)
#' - out_scheme + out_patterns: fully customizable filename templates
#'
#' SAFE erase rule (when remove_detected_if_intersects_dropped_ref = TRUE):
#'   remove detected if intersects(DROPPED_ANY) AND NOT intersects(KEEP_FINAL)
#'
#' Prepare burned-area polygons (reference and detected) for validation
#'
#' Standardizes two polygon datasets (reference and detected) for validation.
#' Works on single files or batch mode (multi-year) when inputs are directories
#' or vectors of files.
#'
#' Main steps (per year):
#' \enumerate{
#'   \item Read polygons (SHP/GPKG), auto-pick a layer for GPKG when needed.
#'   \item Reproject to \code{target_epsg} and enforce valid MULTIPOLYGON geometries.
#'   \item Optionally filter reference by DOY (e.g., summer) and/or attributes.
#'   \item Optionally filter detected polygons by attributes.
#'   \item Clip polygons to an AOI mask.
#'   \item Optionally clip to a burnable mask raster (raster-based clip) and compute
#'         per-polygon burnable area ratios.
#'   \item Write ArcGIS-friendly GeoPackages (WKT1_ESRI) plus optional CSV ratios and
#'         optional rasterized outputs.
#' }
#'
#' Dropped-reference masking (key behavior):
#' \itemize{
#'   \item If \code{remove_detected_if_intersects_dropped_ref = TRUE}, the function
#'         builds a "dropped reference" mask as any reference polygon that does NOT
#'         reach the \emph{final} reference output (i.e., dropped by DOY/attribute filters,
#'         AOI clip, and/or burnable clip).
#'   \item The dropped mask is built to be comparable to detections (which are AOI-clipped),
#'         so reference polygons dropped by filters are first clipped to the AOI before
#'         contributing geometry. Polygons dropped because they never intersect the AOI are
#'         tracked in counters but typically do not contribute geometry (they cannot intersect
#'         AOI-clipped detections).
#'   \item \strong{SAFE erase rule}: detected polygons are removed only if they intersect
#'         the dropped-reference mask \emph{and do not intersect} the \emph{final kept reference}
#'         geometry:
#'         \deqn{remove = intersects(dropped\_any) \;\&\; !intersects(keep\_final)}
#'         This prevents deleting detections that still overlap reference polygons retained after
#'         all filtering/clipping steps (common when dropped and kept reference polygons overlap).
#'   \item The SAFE removal is applied after AOI clip (and again after burnable clip if enabled).
#' }
#'
#' AOI selection precedence (per year):
#' \itemize{
#'   \item If \code{aoi_shapefile} is provided: it is always used.
#'   \item Else, if \code{use_origin_aoi = TRUE}: derive AOI from reference \code{origin}
#'         using \code{borders_shapefile}.
#'   \item Else, if \code{aoi_mask_dir} is provided: try year-based AOI masks from disk.
#' }
#'
#' Burnable mask selection:
#' \itemize{
#'   \item If \code{burnable_raster} is provided: it is used for all years.
#'   \item Else, if \code{burnable_mask_dir} is provided: a year-specific CORINE mask is
#'         auto-picked using \code{corine_year_fun} (or the internal default).
#' }
#'
#' @param ref_path Character. Reference polygons: a file, a vector of files, or a directory.
#' @param det_path Character. Detected polygons: a file, a vector of files, or a directory.
#' @param out_dir Character. Output folder.
#'
#' @param aoi_shapefile Character or NULL. Fixed AOI polygon. If NULL, AOI can be derived
#'   by \code{origin} or picked from \code{aoi_mask_dir}.
#' @param burnable_raster Character or NULL. Fixed burnable mask raster. If NULL, can be
#'   auto-picked from \code{burnable_mask_dir}.
#'
#' @param target_epsg Integer. Target EPSG for processing and outputs.
#' @param ref_layer Character or NULL. Layer name if \code{ref_path} is a GeoPackage. If NULL,
#'   the function auto-picks a suitable layer.
#' @param det_layer Character or NULL. Layer name if \code{det_path} is a GeoPackage. If NULL,
#'   the function auto-picks a suitable layer.
#' @param year_target Integer or NULL. If not NULL, forces processing for this year.
#'
#' @param apply_ref_doy_filter Logical. Apply DOY filter to reference polygons.
#' @param ref_doy_col Character. DOY field in reference polygons.
#' @param ref_doy_min,ref_doy_max Integers. Inclusive DOY window.
#'
#' @param apply_ref_attr_filter Logical. Apply attribute filter to reference polygons.
#' @param ref_filter_col Character vector. Candidate columns for reference attribute filter
#'   (first existing column is used).
#' @param ref_filter_values Character vector. Values to keep for reference polygons.
#'
#' @param apply_det_filter Logical. Apply attribute filter to detected polygons.
#' @param det_filter_col Character vector. Candidate columns for detected attribute filter.
#' @param det_filter_values Character vector. Values to keep for detected polygons.
#'
#' @param remove_detected_if_intersects_dropped_ref Logical. If TRUE, remove detected polygons
#'   using the SAFE erase rule based on reference polygons dropped by ANY stage (filters, AOI, burnable).
#' @param dropped_ref_buffer_m Numeric. Optional buffer (meters; in \code{target_epsg}) applied to
#'   the dropped-reference union before intersection tests.
#' @param dropped_ref_stage Character. Geometry stage used to build the dropped-reference mask:
#'   \code{"aoi"} (recommended) or \code{"burnable"}. In practice, \code{"aoi"} is usually required
#'   if you want polygons dropped by the burnable mask to still contribute geometry for removal.
#'
#' @param binary_burnable Logical. If TRUE, burnable raster is treated as binary (>=0.5 kept).
#'   If FALSE, \code{burnable_classes} may be used.
#' @param burnable_classes Numeric vector or NULL. Classes to keep if \code{binary_burnable = FALSE}.
#'
#' @param do_burnable_clip Logical. If TRUE, clip polygons by burnable raster.
#' @param dissolve_by_id Logical. Dissolve clipped polygons back to their input polygon id.
#'
#' @param save_aoi_only Logical. Write AOI-clipped polygons (GPKG).
#' @param save_burnable Logical. Write AOI+burnable-clipped polygons (GPKG).
#' @param save_ratio_csv Logical. Write per-polygon burnable ratio CSV.
#'
#' @param det_output_tag,ref_output_tag Character or NULL. Optional custom tags for output names.
#' @param output_prefix Character or NULL. Optional prefix for output filenames.
#'
#' @param crop_burnable_to_aoi Logical. Crop burnable raster to AOI bbox (performance).
#' @param verbose Logical. Print progress messages.
#'
#' @param file_pattern Character. Regex used when expanding directories.
#' @param year_regex Character. Regex used to extract year from filenames.
#' @param multi_year Logical or NULL. If NULL, auto-detected. If TRUE, process each year found in filenames.
#' @param use_year_subdirs Logical. If TRUE and \code{multi_year=TRUE}, write outputs under \code{out_dir/YYYY/}.
#'
#' @param aoi_mask_dir Character or NULL. Folder containing year-specific AOI masks.
#' @param aoi_latest_from_year Integer. Years >= this use the "latest" AOI mask.
#' @param aoi_latest_shp,aoi_latest_tif Character. Filenames in \code{aoi_mask_dir}.
#' @param aoi_year_shp_template,aoi_year_tif_template Character templates (sprintf) for older years.
#'
#' @param use_origin_aoi Logical. If TRUE and \code{aoi_shapefile} is NULL, derive AOI from \code{origin}.
#' @param borders_shapefile Character or NULL. Required for origin-based AOI.
#' @param borders_nameunit_col Character. Administrative unit field in borders shapefile.
#' @param origin_col Character. Origin field in reference polygons.
#' @param origin_full_regex Character. Regex: if any normalized origin matches it, AOI becomes full peninsula.
#' @param origin_to_nameunit Named character vector or NULL. Mapping \code{origin_norm -> nameunit_norm}.
#' @param origin_fallback Character. Behavior when no origin can be mapped: \code{"full"} or \code{"stop"}.
#'
#' @param burnable_mask_dir Character or NULL. Folder containing CORINE burnable masks by year.
#' @param burnable_mask_template Character. Template filename in \code{burnable_mask_dir} (sprintf).
#' @param corine_year_fun Function or NULL. Optional function(year)->"1990"/"2000"/...
#'
#' @param save_rasterized Logical. If TRUE, rasterize final polygons to a template grid.
#' @param raster_template Character. One of \code{"aoi_tif"} or \code{"burnable_raster"}.
#' @param raster_layer_name Character. Reserved (not used; kept for symmetry).
#' @param raster_datatype Character. GDAL datatype for output raster.
#' @param raster_gdal Character vector. GDAL creation options.
#'
#' @param disable_s2 Logical. Disable s2 for safer polygon operations.
#' @param on_empty Character. Behavior when a dataset becomes empty: \code{"stop"} or \code{"skip"}.
#'
#' @return A list. In single-year mode, returns a flat list containing output paths plus a \code{$years} entry.
#'   In multi-year mode, returns a list with \code{$years} keyed by year strings. Per year:
#'   \itemize{
#'     \item \code{$reference} includes paths and counters \code{dropped_by_filter_n}, \code{dropped_by_aoi_n},
#'           \code{dropped_by_burn_n} describing how many reference polygons were dropped at each stage.
#'     \item \code{$detected} includes paths and counters \code{erased_by_droppedref_aoi_n},
#'           \code{erased_by_droppedref_burn_n} describing SAFE removals at AOI and burnable stages.
#'   }
#'
#'#' @examples
#' \dontrun{
#' # Example 1: Single-year validation preparation using a fixed AOI
#' # and a fixed burnable raster.
#' #
#' # This example is suitable for recent years where a complete reference
#' # AOI is already available, for example 2025.
#'
#' base_dir <- "C:/OtsuFire/ZENODO/exdata/2025/OTSU_MIN"
#'
#' scoring_dir <- file.path(
#'   base_dir,
#'   "CORINE_ECOREG_seed_30pix"
#' )
#'
#' out_dir <- file.path(
#'   scoring_dir,
#'   "PREP_VALIDATION"
#' )
#'
#' dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#'
#' ref_path <- file.path(
#'   scoring_dir,
#'   "phase2_external_scoring_classes_30",
#'   "poly_metrics",
#'   "external_scored_only_orig_metrics.gpkg"
#' )
#'
#' det_path <- file.path(
#'   scoring_dir,
#'   "phase2_internal_scoring_classes_30",
#'   "poly_metrics",
#'   "internal_scored_only_orig_metrics.gpkg"
#' )
#'
#' res_2025 <- prepare_fire_polys_for_validation(
#'   ref_path = ref_path,
#'   det_path = det_path,
#'   out_dir  = out_dir,
#'
#'   year_target = 2025,
#'   target_epsg = 3035,
#'
#'   ref_layer = NULL,
#'   det_layer = NULL,
#'
#'   aoi_shapefile = "C:/OtsuFire/ZENODO/exdata/iberian_peninsula_proj_final.shp",
#'   use_origin_aoi = FALSE,
#'
#'   burnable_raster = "C:/OtsuFire/ZENODO/exdata/corine18_burnable_mask_binary.tif",
#'   do_burnable_clip = TRUE,
#'
#'   apply_ref_doy_filter = TRUE,
#'   ref_doy_col = "start_doy",
#'   ref_doy_min = 152,
#'   ref_doy_max = 274,
#'
#'   apply_ref_attr_filter = FALSE,
#'   apply_det_filter = FALSE,
#'
#'   remove_detected_if_intersects_dropped_ref = TRUE,
#'   dropped_ref_stage = "aoi",
#'   dropped_ref_buffer_m = 0,
#'
#'   output_prefix = "PREP_",
#'   verbose = TRUE
#' )
#'
#'
#' # Example 2: Multi-year validation preparation using origin-based AOI masks
#' # for years before 2000 and standard AOI masks from 2000 onwards.
#' #
#' # This is useful when pre-2000 reference fire polygons do not come from
#' # complete EFFIS annual series and the validation AOI must be restricted to
#' # the CCAA/Portugal regions indicated by the fire polygon field 'origin'.
#' #
#' # Spanish origin values are matched to NAMEUNIT. Portuguese fires are
#' # assigned to Portugal when origin is 'Landsat' or 'EFFIS+Landsat'.
#' # Ambiguous values such as 'EFFIS+CA' are resolved internally by largest
#' # spatial overlap with the administrative borders, if that logic is enabled
#' # in the function.
#'
#' ref_dir <- "C:/OtsuFire/ZENODO/exdata/reference_fires_1984_2025"
#' det_dir <- "C:/OtsuFire/ZENODO/exdata/detected_fires_1984_2025"
#'
#' out_dir_multi <- "C:/OtsuFire/ZENODO/exdata/PREP_VALIDATION_1984_2025"
#'
#' dir.create(out_dir_multi, recursive = TRUE, showWarnings = FALSE)
#'
#' origin_to_nameunit <- c(
#'   "Andalucia"     = "Andalucía",
#'   "Cataluña"      = "Cataluña/Catalunya",
#'   "CLM"           = "Castilla-La Mancha",
#'   "Madrid"        = "Comunidad de Madrid",
#'   "Navarra"       = "Comunidad Foral de Navarra",
#'   "Valencia"      = "Comunitat Valenciana",
#'   "Landsat"       = "Portugal",
#'   "EFFIS+Landsat" = "Portugal"
#' )
#'
#' res_multi <- prepare_fire_polys_for_validation(
#'   ref_path = ref_dir,
#'   det_path = det_dir,
#'   out_dir  = out_dir_multi,
#'
#'   multi_year = TRUE,
#'   use_year_subdirs = TRUE,
#'   year_regex = "\\\\d{4}",
#'
#'   target_epsg = 3035,
#'   ref_layer = NULL,
#'   det_layer = "scored",
#'
#'   # For years < 2000, the AOI is built from the fire 'origin' field.
#'   # For years >= 2000, the standard AOI picker is used.
#'   aoi_shapefile = NULL,
#'   use_origin_aoi = TRUE,
#'   aoi_latest_from_year = 2000L,
#'
#'   borders_shapefile = "C:/OtsuFire/ZENODO/exdata/iberian_peninsula_comunidades.shp",
#'   borders_nameunit_col = "NAMEUNIT",
#'   origin_col = "origin",
#'   origin_to_nameunit = origin_to_nameunit,
#'   origin_fallback = "stop",
#'
#'   # Standard AOI masks used from year 2000 onwards.
#'   aoi_mask_dir = "C:/OtsuFire/ZENODO/exdata/AOI_masks",
#'   aoi_latest_shp = "Iberian_Peninsula_mask_3035.shp",
#'   aoi_latest_tif = "Iberian_Peninsula_mask_3035.tif",
#'   aoi_year_shp_template = "mask_3035_%s.shp",
#'   aoi_year_tif_template = "mask_3035_%s.tif",
#'
#'   # Year-dependent CORINE burnable masks.
#'   burnable_raster = NULL,
#'   burnable_mask_dir = "C:/OtsuFire/ZENODO/exdata/CORINE_BURNABLE_MASKS",
#'   burnable_mask_template = "burneable_mask_binary_corine_%s_ETRS89.tif",
#'   do_burnable_clip = TRUE,
#'   binary_burnable = TRUE,
#'
#'   apply_ref_doy_filter = TRUE,
#'   ref_doy_col = "start_doy",
#'   ref_doy_min = 152,
#'   ref_doy_max = 274,
#'
#'   apply_ref_attr_filter = FALSE,
#'   apply_det_filter = TRUE,
#'   det_filter_col = c("merged_flag_value", "class3_id"),
#'   det_filter_values = c("good_identified"),
#'
#'   remove_detected_if_intersects_dropped_ref = TRUE,
#'   dropped_ref_stage = "aoi",
#'   dropped_ref_buffer_m = 0,
#'
#'   output_prefix = "PREP_",
#'   verbose = TRUE
#' )
#' }
#'
#' @importFrom sf st_as_sf st_read st_write st_transform st_make_valid st_crs
#' @importFrom sf st_zm st_union st_geometry st_drop_geometry
#' @importFrom sf st_collection_extract st_cast st_is_valid st_is_empty
#' @importFrom sf st_filter st_intersection st_sfc st_sf st_layers
#' @importFrom sf sf_use_s2 st_intersects st_buffer st_combine
#' @importFrom terra rast project crop vect same.crs crs rasterize as.polygons
#' @importFrom terra ext res ncell nlyr values mask writeRaster ifel
#' @importFrom dplyr bind_rows filter group_by summarise across first transmute
#' @importFrom dplyr left_join mutate
#' @importFrom rlang .data
#' @importFrom tools file_ext
#' @importFrom utils write.csv
#' @export
prepare_fire_polys_for_validation <- function(
    # --- inputs (file, vector of files, or directory)
  ref_path,
  det_path,
  out_dir,

  # --- AOI / burnable inputs
  aoi_shapefile = NULL,
  burnable_raster = NULL,

  # --- CRS
  target_epsg = 3035,

  # --- reading layers
  ref_layer = NULL,
  det_layer = "scored",
  year_target = NULL,

  # --- Reference DOY filter
  apply_ref_doy_filter = TRUE,
  ref_doy_col = "start_doy",
  ref_doy_min = 152,
  ref_doy_max = 274,

  # --- Reference attribute filter
  apply_ref_attr_filter = FALSE,
  ref_filter_col = c("merged_flag_value", "class3_id"),
  ref_filter_values = c("good_identified"),

  # --- Remove detected polys that intersect reference polys dropped by any stage
  remove_detected_if_intersects_dropped_ref = TRUE,
  dropped_ref_buffer_m = 0,
  dropped_ref_stage = c("aoi", "burnable"),

  # --- Detected filter
  apply_det_filter = TRUE,
  det_filter_col = c("merged_flag_value", "class3_id"),
  det_filter_values = c("good_identified"),

  # --- Burnable raster definition
  binary_burnable = TRUE,
  burnable_classes = NULL,

  # --- Burnable clip behavior
  do_burnable_clip = TRUE,
  dissolve_by_id = TRUE,

  # --- outputs
  save_aoi_mask  = TRUE,
  save_aoi_only  = TRUE,
  save_burnable  = TRUE,
  save_ratio_csv = TRUE,

  # --- output tags
  det_output_tag = NULL,
  ref_output_tag = NULL,
  output_prefix  = NULL,

  # --- output naming / writing control
  out_scheme = c("default", "compact"),
  out_patterns = NULL,
  out_layer_polys = "polys",
  out_layer_aoi   = "aoi",
  overwrite = TRUE,

  # --- performance
  crop_burnable_to_aoi = TRUE,
  verbose = TRUE,

  # =============================================================================
  # Batch + year-based processing
  # =============================================================================
  file_pattern = "\\.(shp|gpkg)$",
  year_regex = "\\d{4}",
  multi_year = NULL,
  use_year_subdirs = TRUE,

  # =============================================================================
  # Standard AOI picker
  # Used for years >= aoi_latest_from_year unless a fixed aoi_shapefile is supplied.
  # =============================================================================
  aoi_mask_dir = NULL,
  aoi_latest_from_year = 2000L,
  aoi_latest_shp = "Iberian_Peninsula_mask_3035.shp",
  aoi_latest_tif = "Iberian_Peninsula_mask_3035.tif",
  aoi_year_shp_template = "mask_3035_%s.shp",
  aoi_year_tif_template = "mask_3035_%s.tif",

  # =============================================================================
  # Origin-based AOI
  # Used only for years < origin_aoi_before_year.
  # Intended for pre-2000 reference fire polygons, before complete annual
  # EFFIS-based validation masks are available.
  # For years >= origin_aoi_before_year, the standard AOI picker is used.
  # =============================================================================
  use_origin_aoi = TRUE,
  origin_aoi_before_year = 2000L,
  borders_shapefile = NULL,
  borders_nameunit_col = "NAMEUNIT",
  origin_col = "origin",
  origin_to_nameunit = NULL,
  resolve_unmatched_origin_by_overlap = TRUE,
  origin_fallback = c("stop", "full"),

  # =============================================================================
  # Burnable CORINE picker
  # =============================================================================
  burnable_mask_dir = NULL,
  burnable_mask_template = "burneable_mask_binary_corine_%s_ETRS89.tif",
  corine_year_fun = NULL,

  # --- optional rasterization output
  save_rasterized = FALSE,
  raster_template = c("aoi_tif", "burnable_raster"),
  raster_datatype = "INT1U",
  raster_gdal = c("COMPRESS=LZW", "TILED=YES"),

  # --- sf robustness
  disable_s2 = TRUE,
  arcgis_wkt1 = TRUE,

  # --- behavior if empty after filtering/clipping
  on_empty = c("stop", "skip")
) {

  suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(dplyr)
  })

  out_scheme <- match.arg(out_scheme)
  dropped_ref_stage <- match.arg(dropped_ref_stage)
  origin_fallback <- match.arg(origin_fallback)
  on_empty <- match.arg(on_empty)
  raster_template <- match.arg(raster_template)

  msg <- function(...) if (isTRUE(verbose)) message(...)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------------------------
  # Small utilities
  # ------------------------------------------------------------------
  `%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

  .safe_slug <- function(x, max_len = 90) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    x <- gsub("_+", "_", x)
    if (nchar(x) > max_len) x <- substr(x, 1, max_len)
    if (!nzchar(x)) x <- "X"
    x
  }

  .safe_overwrite <- function(path) {
    if (!file.exists(path)) return(invisible(TRUE))
    if (!isTRUE(overwrite)) return(invisible(FALSE))
    try(unlink(path), silent = TRUE)
    if (file.exists(path)) stop("Cannot overwrite (file locked?): ", path)
    invisible(TRUE)
  }

  .extract_year_from_path <- function(p, year_regex = "\\d{4}") {
    b <- basename(p)
    m <- regexpr(year_regex, b, perl = TRUE)
    if (m[1] == -1) return(NA_integer_)
    as.integer(substr(b, m[1], m[1] + attr(m, "match.length") - 1))
  }

  .expand_paths <- function(x) {
    if (is.null(x)) return(character(0))
    if (length(x) == 1 && dir.exists(x)) {
      return(list.files(x, pattern = file_pattern, full.names = TRUE))
    }
    x
  }

  # ------------------------------------------------------------------
  # Output naming (patterns)
  # ------------------------------------------------------------------
  .render_tpl <- function(tpl, pref, tag = "", epsg = "", stage = "") {
    out <- as.character(tpl)
    out <- gsub("\\{pref\\}",  pref,  out)
    out <- gsub("\\{tag\\}",   tag,   out)
    out <- gsub("\\{epsg\\}",  as.character(epsg),  out)
    out <- gsub("\\{stage\\}", stage, out)
    out
  }

  .default_patterns <- function(scheme) {
    if (identical(scheme, "compact")) {
      return(list(
        aoi_mask       = "{pref}AOI_epsg{epsg}.gpkg",
        polys_aoi      = "{pref}{tag}.gpkg",
        polys_burnable = "{pref}{tag}_burn.gpkg",
        ratio          = "{pref}{tag}_ratio.csv",
        raster         = "{pref}{tag}_{stage}.tif"
      ))
    }
    list(
      aoi_mask       = "{pref}AOI_mask_epsg{epsg}.gpkg",
      polys_aoi      = "{pref}{tag}__AOI.gpkg",
      polys_burnable = "{pref}{tag}__AOI_BURNABLE.gpkg",
      ratio          = "{pref}{tag}__burnable_ratio.csv",
      raster         = "{pref}{tag}__{stage}.tif"
    )
  }

  .get_patterns <- function() {
    p <- .default_patterns(out_scheme)
    if (!is.null(out_patterns)) {
      if (!is.list(out_patterns)) stop("out_patterns must be a named list.")
      p[names(out_patterns)] <- out_patterns
    }
    needed <- c("aoi_mask","polys_aoi","polys_burnable","ratio","raster")
    miss <- setdiff(needed, names(p))
    if (length(miss) > 0) stop("out_patterns missing: ", paste(miss, collapse = ", "))
    p
  }

  .patterns <- .get_patterns()

  .outfile <- function(kind, year_out_dir, pref, tag = "", epsg = "", stage = "") {
    if (!kind %in% names(.patterns)) stop("Unknown output kind: ", kind)
    fn <- .render_tpl(.patterns[[kind]], pref = pref, tag = tag, epsg = epsg, stage = stage)
    file.path(year_out_dir, fn)
  }

  # ------------------------------------------------------------------
  # ArcGIS-friendly GPKG writer (WKT1_ESRI)
  # ------------------------------------------------------------------
  st_write_gpkg_arcgis <- function(x, dsn, layer, delete_dsn = FALSE, quiet = TRUE) {
    x <- sf::st_as_sf(x)
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")

    old <- Sys.getenv("OGR_WKT_FORMAT", unset = NA_character_)
    if (isTRUE(arcgis_wkt1)) {
      Sys.setenv(OGR_WKT_FORMAT = "WKT1_ESRI")
      on.exit({
        if (!is.na(old)) Sys.setenv(OGR_WKT_FORMAT = old) else Sys.unsetenv("OGR_WKT_FORMAT")
      }, add = TRUE)
    }

    sf::st_write(
      x, dsn,
      layer = layer,
      driver = "GPKG",
      delete_dsn = delete_dsn,
      quiet = quiet
    )
  }

  # ------------------------------------------------------------------
  # Robust sf utilities (prevents "st_bbox applied to list" failures)
  # IMPORTANT: Never rename columns on an sf object directly.
  #            We rename on st_drop_geometry() then reattach geometry.
  # ------------------------------------------------------------------
  fix_sf <- function(x, epsg = NULL) {
    x <- sf::st_as_sf(x)
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")

    if (nrow(x) == 0) {
      if (!is.null(epsg)) x <- sf::st_set_crs(x, epsg)
      return(x)
    }

    if (any(!sf::st_is_valid(x))) x <- sf::st_make_valid(x)
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    if (nrow(x) == 0) {
      if (!is.null(epsg)) x <- sf::st_set_crs(x, epsg)
      return(x)
    }

    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    if (nrow(x) == 0) {
      if (!is.null(epsg)) x <- sf::st_set_crs(x, epsg)
      return(x)
    }

    geom <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)

    df <- sf::st_drop_geometry(x)
    names(df) <- make.unique(names(df), sep = "_")

    out <- sf::st_sf(df, geometry = geom)
    sf::st_crs(out) <- sf::st_crs(x)

    if (!is.null(epsg)) {
      if (is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs(epsg)
      out <- sf::st_transform(out, epsg)
    }

    out
  }

  as_sfc_any <- function(g, crs = NA) {
    if (inherits(g, "sf")) g <- sf::st_geometry(g)
    if (inherits(g, "sfc")) return(g)
    if (inherits(g, "sfg")) return(sf::st_sfc(g, crs = crs))
    if (is.list(g)) return(sf::st_sfc(g, crs = crs))
    stop("as_sfc_any: cannot coerce to sfc from class: ", paste(class(g), collapse = ", "))
  }

  normalize_sfc1 <- function(g, crs) {
    g <- as_sfc_any(g, crs = sf::st_crs(crs))
    if (is.na(sf::st_crs(g))) sf::st_crs(g) <- sf::st_crs(crs)
    g <- sf::st_make_valid(g)
    if (length(g) > 1) {
      uu <- sf::st_union(g)
      if (inherits(uu, "sfg")) uu <- sf::st_sfc(uu, crs = sf::st_crs(g))
      g <- uu
    }
    g
  }

  rbind_sf_safe <- function(lst, epsg = NULL) {
    lst <- Filter(function(z) inherits(z, "sf"), lst)
    if (length(lst) == 0) return(fix_sf(sf::st_sf(geometry = sf::st_sfc(crs = epsg %||% target_epsg))[0, ], epsg = epsg %||% target_epsg))

    lst <- lapply(lst, function(z) fix_sf(z, epsg = epsg))

    # align attributes
    gcol <- attr(lst[[1]], "sf_column")
    if (is.null(gcol) || !nzchar(gcol) || !(gcol %in% names(lst[[1]]))) gcol <- "geometry"

    all_cols <- unique(unlist(lapply(lst, function(z) setdiff(names(z), gcol))))
    lst2 <- lapply(lst, function(z) {
      gcolz <- attr(z, "sf_column")
      if (is.null(gcolz) || !nzchar(gcolz) || !(gcolz %in% names(z))) gcolz <- "geometry"
      for (cc in setdiff(all_cols, names(z))) z[[cc]] <- NA
      z <- z[, c(all_cols, gcolz), drop = FALSE]
      z
    })

    out <- suppressWarnings(do.call(rbind, lst2))
    fix_sf(out, epsg = epsg)
  }

  read_sf_any <- function(path, layer = NULL,
                          prefer_layers = c("polys", "scored_metrics", "scored", "ref", "internal_flagged", "scored_metrics")) {
    if (!file.exists(path)) stop("File not found: ", path)
    ext <- tolower(tools::file_ext(path))

    if (ext == "gpkg") {
      lyr_names <- sf::st_layers(path)$name
      if (!is.null(layer) && nzchar(layer) && !(layer %in% lyr_names)) {
        stop("Layer '", layer, "' not found in: ", basename(path),
             "\nAvailable: ", paste(lyr_names, collapse = ", "))
      }
      if (is.null(layer) || !nzchar(layer)) {
        hit <- prefer_layers[prefer_layers %in% lyr_names]
        layer <- if (length(hit) > 0) hit[1] else lyr_names[1]
      }
      x <- sf::st_read(path, layer = layer, quiet = TRUE)
    } else {
      x <- sf::st_read(path, quiet = TRUE)
    }

    if (is.na(sf::st_crs(x))) stop("Input has NA CRS: ", path)
    fix_sf(sf::st_transform(x, target_epsg), epsg = target_epsg)
  }

  # ------------------------------------------------------------------
  # Filters
  # ------------------------------------------------------------------
  build_attr_filter_fun <- function(cols, values, tag = "DATA") {
    function(x) {
      x <- fix_sf(x, epsg = target_epsg)
      cols_present <- cols[cols %in% names(x)]
      if (length(cols_present) == 0) {
        stop(tag, " polygons have none of these columns: ", paste(cols, collapse = ", "))
      }
      col <- cols_present[1]
      keep_idx <- as.character(x[[col]]) %in% as.character(values)
      fix_sf(x[keep_idx %in% TRUE, , drop = FALSE], epsg = target_epsg)
    }
  }

  compose_filters <- function(f1, f2) {
    if (is.null(f1) && is.null(f2)) return(NULL)
    if (is.null(f1)) return(f2)
    if (is.null(f2)) return(f1)
    function(x) f2(f1(x))
  }

  ref_doy_filter_fun <- NULL
  if (isTRUE(apply_ref_doy_filter)) {
    ref_doy_filter_fun <- function(x) {
      x <- fix_sf(x, epsg = target_epsg)
      if (!(ref_doy_col %in% names(x))) stop("Reference polygons have no column: ", ref_doy_col)
      doy <- suppressWarnings(as.integer(x[[ref_doy_col]]))
      keep <- !is.na(doy) & doy >= ref_doy_min & doy <= ref_doy_max
      fix_sf(x[keep %in% TRUE, , drop = FALSE], epsg = target_epsg)
    }
  }

  ref_attr_filter_fun <- NULL
  if (isTRUE(apply_ref_attr_filter)) {
    ref_attr_filter_fun <- build_attr_filter_fun(ref_filter_col, ref_filter_values, tag = "REFERENCE")
  }

  det_filter_fun <- NULL
  if (isTRUE(apply_det_filter)) {
    det_filter_fun <- build_attr_filter_fun(det_filter_col, det_filter_values, tag = "DETECTED")
  }

  ref_filter_fun <- compose_filters(ref_doy_filter_fun, ref_attr_filter_fun)

  # ------------------------------------------------------------------
  # AOI selection helpers
  # ------------------------------------------------------------------
  pick_aoi_paths <- function(year) {
    if (is.null(aoi_mask_dir) || !nzchar(aoi_mask_dir)) return(list(shp = NA_character_, tif = NA_character_))
    year <- as.integer(year)
    if (!is.finite(year)) return(list(shp = NA_character_, tif = NA_character_))

    if (year >= as.integer(aoi_latest_from_year)) {
      list(
        shp = file.path(aoi_mask_dir, aoi_latest_shp),
        tif = file.path(aoi_mask_dir, aoi_latest_tif)
      )
    } else {
      list(
        shp = file.path(aoi_mask_dir, sprintf(aoi_year_shp_template, year)),
        tif = file.path(aoi_mask_dir, sprintf(aoi_year_tif_template, year))
      )
    }
  }

  load_aoi_for_year <- function(year) {
    mp <- pick_aoi_paths(year)
    shp_ok <- is.character(mp$shp) && nzchar(mp$shp) && file.exists(mp$shp)
    tif_ok <- is.character(mp$tif) && nzchar(mp$tif) && file.exists(mp$tif)

    if (shp_ok) {
      aoi <- sf::st_read(mp$shp, quiet = TRUE)
      if (is.na(sf::st_crs(aoi))) stop("AOI shapefile has NA CRS: ", mp$shp)
      aoi <- fix_sf(sf::st_transform(sf::st_make_valid(aoi), target_epsg), epsg = target_epsg)
      aoi_sfc <- normalize_sfc1(sf::st_geometry(aoi), crs = sf::st_crs(aoi))
      return(list(aoi_sf = aoi, aoi_sfc = aoi_sfc, aoi_tif = if (tif_ok) mp$tif else NA_character_))
    }

    if (tif_ok) {
      aoi_r <- terra::rast(mp$tif)
      if (terra::nlyr(aoi_r) > 1) aoi_r <- aoi_r[[1]]
      if (is.na(terra::crs(aoi_r)) || !nzchar(terra::crs(aoi_r))) stop("AOI raster has NA/empty CRS: ", mp$tif)

      aoi_bin <- terra::ifel(aoi_r > 0, 1, NA)
      vv <- terra::as.polygons(aoi_bin, dissolve = TRUE, na.rm = TRUE)
      if (is.null(vv) || nrow(vv) == 0) stop("AOI raster produced 0 polygons: ", mp$tif)

      aoi <- sf::st_as_sf(vv)
      aoi <- fix_sf(sf::st_set_crs(sf::st_make_valid(aoi), target_epsg), epsg = target_epsg)
      aoi_sfc <- normalize_sfc1(sf::st_geometry(aoi), crs = sf::st_crs(aoi))

      return(list(aoi_sf = aoi, aoi_sfc = aoi_sfc, aoi_tif = mp$tif))
    }

    stop("Cannot load AOI for year=", year, ". Provide aoi_shapefile, or enable use_origin_aoi, or set aoi_mask_dir.")
  }

  # ------------------------------------------------------------------
  # Origin-based AOI helpers
  # ------------------------------------------------------------------
  # This block selects complete CCAA/Portugal polygons using the
  # reference fire polygon field `origin`.
  #
  # Explicit origins are matched by dictionary:
  #   Andalucia      -> Andalucía
  #   Cataluña       -> Cataluña/Catalunya
  #   CLM            -> Castilla-La Mancha
  #   Madrid         -> Comunidad de Madrid
  #   Navarra        -> Comunidad Foral de Navarra
  #   Valencia       -> Comunitat Valenciana
  #   Landsat        -> Portugal
  #
  # ------------------------------------------------------------------

  .admin_key <- function(x) {

    x <- as.character(x)
    x[is.na(x)] <- ""

    x <- iconv(
      x,
      from = "",
      to = "ASCII//TRANSLIT"
    )

    x <- tolower(x)
    x <- gsub("&", " y ", x)
    x <- gsub("[^a-z0-9]+", " ", x)
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)

    out <- vapply(
      strsplit(x, " "),
      function(tokens) {

        tokens <- tokens[
          !tokens %in% c(
            "comunidad",
            "comunitat",
            "autonoma",
            "autonoma",
            "region",
            "principado",
            "foral",
            "de",
            "del",
            "la",
            "las",
            "los",
            "el",
            "y"
          )
        ]

        paste(tokens, collapse = "")
      },
      character(1)
    )

    out <- gsub("[^a-z0-9]", "", out)

    out
  }


  .default_origin_to_nameunit <- function() {

    x <- c(
      "Andalucia"     = "Andalucía",
      "Cataluña"      = "Cataluña/Catalunya",
      "CLM"           = "Castilla-La Mancha",
      "Madrid"        = "Comunidad de Madrid",
      "Navarra"       = "Comunidad Foral de Navarra",
      "Valencia"      = "Comunitat Valenciana",
      "Landsat"       = "Portugal"
    )

    names(x) <- .admin_key(names(x))

    x
  }


  .read_borders_once <- local({

    cache <- NULL

    function() {

      if (!is.null(cache)) {
        return(cache)
      }

      if (is.null(borders_shapefile) || !nzchar(borders_shapefile)) {
        stop(
          "Origin-based AOI requires 'borders_shapefile' ",
          "when aoi_shapefile is NULL."
        )
      }

      if (!file.exists(borders_shapefile)) {
        stop("borders_shapefile not found: ", borders_shapefile)
      }

      b <- sf::st_read(borders_shapefile, quiet = TRUE)

      if (is.na(sf::st_crs(b))) {
        stop("borders_shapefile has NA CRS: ", borders_shapefile)
      }

      b <- fix_sf(
        sf::st_transform(sf::st_make_valid(b), target_epsg),
        epsg = target_epsg
      )

      if (!(borders_nameunit_col %in% names(b))) {
        stop(
          "borders_shapefile missing column '", borders_nameunit_col, "'. ",
          "Available fields are: ",
          paste(names(b), collapse = ", ")
        )
      }

      b$NAMEUNIT <- as.character(b[[borders_nameunit_col]])
      b$nameunit_key <- .admin_key(b$NAMEUNIT)

      cache <<- b
      cache
    }
  })


  .resolve_nameunit <- function(target_name, borders_sf) {

    target_key <- .admin_key(target_name)

    direct <- borders_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(nameunit_key == target_key)

    if (nrow(direct) > 0) {
      return(direct$NAMEUNIT[1])
    }

    d <- adist(target_key, borders_sf$nameunit_key)
    best_id <- which.min(d)
    best_dist <- as.numeric(d[best_id])

    if (is.finite(best_dist) && best_dist <= 3) {
      return(borders_sf$NAMEUNIT[best_id])
    }

    NA_character_
  }


  .match_origin_values_to_nameunit <- function(origin_values, borders_sf) {

    origin_values <- sort(unique(as.character(origin_values)))
    origin_values <- origin_values[!is.na(origin_values) & nzchar(origin_values)]

    map_vec <- origin_to_nameunit

    if (is.null(map_vec)) {
      map_vec <- .default_origin_to_nameunit()
    } else {
      names(map_vec) <- .admin_key(names(map_vec))
    }

    out <- lapply(origin_values, function(origin_value) {

      origin_key <- .admin_key(origin_value)

      matched_nameunit <- NA_character_
      match_method <- NA_character_

      # 1) Manual dictionary
      if (origin_key %in% names(map_vec)) {

        target_name <- unname(map_vec[[origin_key]])
        matched_nameunit <- .resolve_nameunit(target_name, borders_sf)
        match_method <- "manual_dictionary"
      }

      # 2) Direct match against NAMEUNIT
      if (is.na(matched_nameunit)) {

        direct <- borders_sf %>%
          sf::st_drop_geometry() %>%
          dplyr::filter(nameunit_key == origin_key)

        if (nrow(direct) > 0) {
          matched_nameunit <- direct$NAMEUNIT[1]
          match_method <- "direct_key_match"
        }
      }

      # 3) Fuzzy fallback for small spelling differences
      if (is.na(matched_nameunit)) {

        d <- adist(origin_key, borders_sf$nameunit_key)
        best_id <- which.min(d)
        best_dist <- as.numeric(d[best_id])

        if (is.finite(best_dist) && best_dist <= 3) {
          matched_nameunit <- borders_sf$NAMEUNIT[best_id]
          match_method <- paste0("fuzzy_adist_", best_dist)
        }
      }

      data.frame(
        origin = origin_value,
        origin_key = origin_key,
        matched_NAMEUNIT = matched_nameunit,
        match_method = match_method,
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(out)
  }


  .read_ref_polys_for_origin_aoi <- function(files, layer = NULL) {

    xs <- lapply(files, function(fp) {

      x <- read_sf_any(fp, layer = layer)

      if (nrow(x) == 0) {
        return(x)
      }

      origin_idx <- which(tolower(names(x)) == tolower(origin_col))

      if (length(origin_idx) == 0) {
        stop(
          "Reference file has no origin column '", origin_col, "': ",
          basename(fp),
          "\nAvailable fields are: ",
          paste(names(x), collapse = ", ")
        )
      }

      origin_real_col <- names(x)[origin_idx[1]]

      x$.__origin_raw <- as.character(x[[origin_real_col]])
      x$.__origin_key <- .admin_key(x$.__origin_raw)
      x$.__source_file <- basename(fp)

      x
    })

    x_all <- rbind_sf_safe(xs, epsg = target_epsg)

    if (nrow(x_all) > 0) {
      x_all$.__origin_row_id <- seq_len(nrow(x_all))
    }

    x_all
  }


  .assign_unmatched_origins_by_largest_overlap <- function(
    fire_sf,
    origin_match,
    borders_sf
  ) {

    if (nrow(fire_sf) == 0) {
      return(list(
        targets = character(0),
        diagnostic = data.frame()
      ))
    }

    unmatched_origins <- origin_match %>%
      dplyr::filter(is.na(matched_NAMEUNIT)) %>%
      dplyr::pull(origin)

    if (length(unmatched_origins) == 0) {
      return(list(
        targets = character(0),
        diagnostic = data.frame()
      ))
    }

    fire_unmatched <- fire_sf %>%
      dplyr::filter(.__origin_raw %in% unmatched_origins)

    fire_unmatched <- fix_sf(fire_unmatched, epsg = target_epsg)

    if (nrow(fire_unmatched) == 0) {
      return(list(
        targets = character(0),
        diagnostic = data.frame()
      ))
    }

    borders_min <- borders_sf[, c("NAMEUNIT")]

    ov <- suppressWarnings(
      sf::st_intersection(
        fire_unmatched[, c(".__origin_row_id", ".__origin_raw")],
        borders_min
      )
    )

    ov <- fix_sf(ov, epsg = target_epsg)

    if (nrow(ov) == 0) {
      return(list(
        targets = character(0),
        diagnostic = data.frame()
      ))
    }

    ov$overlap_ha <- as.numeric(sf::st_area(ov)) / 10000

    best <- ov %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(.__origin_row_id) %>%
      dplyr::arrange(dplyr::desc(overlap_ha), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::rename(
        origin = .__origin_raw,
        matched_NAMEUNIT = NAMEUNIT
      ) %>%
      dplyr::mutate(
        origin_key = .admin_key(origin),
        match_method = "spatial_largest_overlap"
      )

    targets <- sort(unique(best$matched_NAMEUNIT))

    list(
      targets = targets,
      diagnostic = best
    )
  }


  .load_aoi_from_origin_for_year <- function(ref_files_this_year, ref_layer_this_year) {

    borders_sf <- .read_borders_once()

    fire_sf <- .read_ref_polys_for_origin_aoi(
      files = ref_files_this_year,
      layer = ref_layer_this_year
    )

    if (nrow(fire_sf) == 0) {

      if (origin_fallback == "stop") {
        stop("Origin-based AOI failed: reference polygons are empty.")
      }

      geom <- sf::st_union(sf::st_geometry(borders_sf))

      aoi_sf <- fix_sf(
        sf::st_sf(geometry = geom, crs = sf::st_crs(borders_sf)),
        epsg = target_epsg
      )

      aoi_sfc <- normalize_sfc1(
        sf::st_geometry(aoi_sf),
        crs = sf::st_crs(aoi_sf)
      )

      return(list(
        aoi_sf = aoi_sf,
        aoi_sfc = aoi_sfc,
        aoi_tif = NA_character_,
        origin_raw = character(0),
        origin_match = data.frame(),
        origin_spatial_match = data.frame()
      ))
    }

    origin_raw <- sort(unique(fire_sf$.__origin_raw))

    origin_match <- .match_origin_values_to_nameunit(
      origin_values = origin_raw,
      borders_sf = borders_sf
    )

    explicit_targets <- origin_match %>%
      dplyr::filter(!is.na(matched_NAMEUNIT)) %>%
      dplyr::pull(matched_NAMEUNIT) %>%
      unique()

    spatial_match <- .assign_unmatched_origins_by_largest_overlap(
      fire_sf = fire_sf,
      origin_match = origin_match,
      borders_sf = borders_sf
    )

    spatial_targets <- spatial_match$targets

    all_targets <- sort(unique(c(explicit_targets, spatial_targets)))

    unresolved <- origin_match %>%
      dplyr::filter(is.na(matched_NAMEUNIT))

    if (nrow(unresolved) > 0) {

      resolved_spatial_origins <- unique(spatial_match$diagnostic$origin)

      unresolved <- unresolved %>%
        dplyr::filter(!origin %in% resolved_spatial_origins)
    }

    if (nrow(unresolved) > 0) {

      msg(
        "Origin-based AOI unresolved origin values: ",
        paste(unresolved$origin, collapse = "; ")
      )

      if (origin_fallback == "stop") {
        stop(
          "Origin-based AOI failed because some origin values could not be ",
          "matched to NAMEUNIT and could not be resolved spatially: ",
          paste(unresolved$origin, collapse = "; ")
        )
      }
    }

    if (length(all_targets) == 0) {

      if (origin_fallback == "stop") {
        stop("Origin-based AOI failed: no target CCAA/Portugal polygons selected.")
      }

      geom <- sf::st_union(sf::st_geometry(borders_sf))

      aoi_sf <- fix_sf(
        sf::st_sf(geometry = geom, crs = sf::st_crs(borders_sf)),
        epsg = target_epsg
      )

    } else {

      aoi_sf <- borders_sf %>%
        dplyr::filter(NAMEUNIT %in% all_targets)

      aoi_sf <- fix_sf(aoi_sf, epsg = target_epsg)
    }

    if (nrow(aoi_sf) == 0) {

      if (origin_fallback == "stop") {
        stop("Origin-based AOI failed: selected AOI is empty.")
      }

      geom <- sf::st_union(sf::st_geometry(borders_sf))

      aoi_sf <- fix_sf(
        sf::st_sf(geometry = geom, crs = sf::st_crs(borders_sf)),
        epsg = target_epsg
      )
    }

    aoi_sfc <- normalize_sfc1(
      sf::st_geometry(aoi_sf),
      crs = sf::st_crs(aoi_sf)
    )

    msg(
      "Origin-based AOI selected regions: ",
      paste(sort(unique(aoi_sf$NAMEUNIT)), collapse = "; ")
    )

    list(
      aoi_sf = aoi_sf,
      aoi_sfc = aoi_sfc,
      aoi_tif = NA_character_,
      origin_raw = origin_raw,
      origin_match = origin_match,
      origin_spatial_match = spatial_match$diagnostic
    )
  }
  # ------------------------------------------------------------------
  # Burnable picker + loader
  # ------------------------------------------------------------------
  default_corine_year <- function(y) {
    y <- as.integer(y)
    if (y >= 1984 && y <= 1999) "1990"
    else if (y <= 2005) "2000"
    else if (y <= 2011) "2006"
    else if (y <= 2017) "2012"
    else "2018"
  }

  pick_burnable_path <- function(year) {
    if (!is.null(burnable_raster) && nzchar(burnable_raster)) return(burnable_raster)
    if (is.null(burnable_mask_dir) || !nzchar(burnable_mask_dir)) return(NA_character_)

    fun <- corine_year_fun
    if (is.null(fun)) fun <- default_corine_year

    cy <- fun(year)
    file.path(burnable_mask_dir, sprintf(burnable_mask_template, cy))
  }

  load_burnable_for_year <- function(year, aoi_sf) {
    bpath <- pick_burnable_path(year)
    if (!is.character(bpath) || !nzchar(bpath) || !file.exists(bpath)) {
      stop("Burnable raster not found for year=", year, ". Got: ", bpath)
    }

    burn <- terra::rast(bpath)
    if (is.na(terra::crs(burn)) || !nzchar(terra::crs(burn))) stop("burnable raster has NA/empty CRS: ", bpath)

    aoi_v <- terra::vect(fix_sf(aoi_sf, epsg = target_epsg))
    if (!terra::same.crs(burn, aoi_v)) {
      burn <- terra::project(burn, terra::crs(aoi_v), method = "near")
    }

    if (isTRUE(crop_burnable_to_aoi)) {
      burn <- terra::crop(burn, aoi_v, snap = "out")
    }

    if (isTRUE(binary_burnable)) {
      burn[burn < 0.5]  <- NA
      burn[burn >= 0.5] <- 1
    } else {
      if (is.null(burnable_classes)) {
        burn[is.na(burn) | burn <= 0] <- NA
        burn[!is.na(burn)] <- 1
      } else {
        burn <- terra::ifel(burn %in% burnable_classes, 1, NA)
      }
    }
    names(burn) <- "burn"
    burn
  }

  # ------------------------------------------------------------------
  # Clip helpers
  # ------------------------------------------------------------------
  clip_to_aoi <- function(x, aoi_sfc) {
    x <- fix_sf(x, epsg = target_epsg)
    if (nrow(x) == 0) return(x)

    aoi_sfc <- normalize_sfc1(aoi_sfc, crs = sf::st_crs(x))
    if (sf::st_crs(aoi_sfc) != sf::st_crs(x)) {
      aoi_sfc <- sf::st_transform(aoi_sfc, sf::st_crs(x))
      aoi_sfc <- normalize_sfc1(aoi_sfc, crs = sf::st_crs(x))
    }

    x <- sf::st_filter(x, aoi_sfc, .predicate = sf::st_intersects)
    if (nrow(x) == 0) return(x)

    x <- suppressWarnings(sf::st_intersection(x, aoi_sfc))
    fix_sf(x, epsg = target_epsg)
  }

  clip_to_burnable_raster <- function(
    x_sf,
    burn_mask,
    dissolve_by_id = TRUE,
    id_col = ".__id",
    touches = TRUE,
    snap = "out",
    expand_cells = 1L
  ) {
    x_sf <- fix_sf(x_sf, epsg = target_epsg)
    if (nrow(x_sf) == 0) return(x_sf)
    if (is.null(burn_mask)) stop("clip_to_burnable_raster: burn_mask is NULL")

    if (!id_col %in% names(x_sf)) x_sf[[id_col]] <- seq_len(nrow(x_sf))

    rid_col <- ".__rid"
    if (rid_col %in% names(x_sf)) rid_col <- make.unique(c(names(x_sf), rid_col))[length(names(x_sf)) + 1]
    x_sf[[rid_col]] <- seq_len(nrow(x_sf))

    v <- terra::vect(x_sf)

    ext0 <- terra::ext(v)
    if (!is.null(expand_cells) && is.finite(expand_cells) && expand_cells > 0) {
      rxy <- terra::res(burn_mask)
      dx  <- abs(as.numeric(rxy[1])) * expand_cells
      dy  <- abs(as.numeric(rxy[2])) * expand_cells
      ext0[1] <- ext0[1] - dx
      ext0[2] <- ext0[2] + dx
      ext0[3] <- ext0[3] - dy
      ext0[4] <- ext0[4] + dy
    }

    burn_local <- terra::crop(burn_mask, ext0, snap = snap)
    if (is.null(burn_local) || terra::ncell(burn_local) == 0) return(x_sf[0, , drop = FALSE])

    rid_r <- terra::rasterize(v, burn_local, field = rid_col, background = NA, touches = touches)
    rid_r <- rid_r * burn_local

    vv <- terra::as.polygons(rid_r, dissolve = TRUE, na.rm = TRUE)
    if (is.null(vv) || nrow(vv) == 0) return(x_sf[0, , drop = FALSE])

    names(vv) <- rid_col
    out_sf <- sf::st_as_sf(vv)
    out_sf <- fix_sf(out_sf, epsg = target_epsg)
    if (nrow(out_sf) == 0) return(out_sf)

    out_sf[[rid_col]] <- suppressWarnings(as.integer(out_sf[[rid_col]]))

    attrs <- sf::st_drop_geometry(x_sf)
    attrs <- attrs[!duplicated(attrs[[rid_col]]), , drop = FALSE]

    out_sf <- dplyr::left_join(out_sf, attrs, by = rid_col)

    if (isTRUE(dissolve_by_id) && nrow(out_sf) > 0) {
      ids <- out_sf[[id_col]]
      ids_u <- unique(ids)
      first_idx <- match(ids_u, ids)

      # IMPORTANT (robustness): sf::st_union() may return either a single sfg
      # or a length-1 sfc depending on sf version and geometry type. If we keep
      # sfc objects inside a list and then call st_sfc(list_of_sfc), we can end
      # up with nested list geometries that later trigger:
      #   Error in UseMethod("st_bbox"): no applicable method for 'st_bbox' applied to an object of class "list"
      # Force each dissolved geometry to be an sfg before rebuilding the sfc.
      geoms <- lapply(ids_u, function(idv) {
        gg <- sf::st_geometry(out_sf)[ids == idv]
        uu <- suppressWarnings(sf::st_union(gg))
        if (inherits(uu, "sfc")) {
          if (length(uu) == 0) return(sf::st_geometrycollection())
          return(uu[[1]])
        }
        if (!inherits(uu, "sfg")) {
          # Defensive fallback: keep an empty geometry rather than a list
          return(sf::st_geometrycollection())
        }
        uu
      })

      attrs2 <- sf::st_drop_geometry(out_sf)[first_idx, , drop = FALSE]
      attrs2[[id_col]] <- ids_u

      geom_sfc <- do.call(sf::st_sfc, c(geoms, list(crs = sf::st_crs(out_sf))))
      out_sf <- sf::st_sf(attrs2, geometry = geom_sfc)
      out_sf <- fix_sf(out_sf, epsg = target_epsg)
    }

    out_sf
  }

  compute_ratio <- function(before_sf, after_sf, id_col = ".__id") {
    before_sf <- fix_sf(before_sf, epsg = target_epsg)
    after_sf  <- fix_sf(after_sf, epsg = target_epsg)

    if (!id_col %in% names(before_sf)) stop("compute_ratio: before_sf missing ", id_col)
    if (!id_col %in% names(after_sf))  stop("compute_ratio: after_sf missing ", id_col)

    area_before_ha <- as.numeric(sf::st_area(before_sf)) / 10000
    before_id <- before_sf[[id_col]]

    area_after_ha_each <- as.numeric(sf::st_area(after_sf)) / 10000
    after_id <- after_sf[[id_col]]

    after_sum <- tapply(area_after_ha_each, after_id, sum, na.rm = TRUE)
    area_after_ha <- after_sum[match(as.character(before_id), names(after_sum))]
    area_after_ha[is.na(area_after_ha)] <- 0

    burnable_ratio <- ifelse(area_before_ha > 0, area_after_ha / area_before_ha, 0)

    data.frame(
      id = before_id,
      area_before_ha = area_before_ha,
      area_after_ha  = area_after_ha,
      burnable_ratio = burnable_ratio,
      nonburn_ratio  = 1 - burnable_ratio,
      stringsAsFactors = FALSE
    )
  }

  # ------------------------------------------------------------------
  # Rasterize polygons to template grid
  # ------------------------------------------------------------------
  rasterize_polys <- function(polys_sf, template_r, out_tif) {
    if (!isTRUE(save_rasterized)) return(NA_character_)
    if (is.null(template_r)) stop("rasterize_polys: template raster is NULL")

    tr <- template_r
    if (terra::nlyr(tr) > 1) tr <- tr[[1]]

    if (nrow(polys_sf) == 0) {
      r <- tr
      terra::values(r) <- NA
    } else {
      v <- terra::vect(fix_sf(polys_sf, epsg = target_epsg))
      r <- terra::rasterize(v, tr, field = 1, background = NA)
      r <- terra::mask(r, tr)
    }

    terra::writeRaster(
      r,
      filename = out_tif,
      overwrite = TRUE,
      datatype = raster_datatype,
      gdal = raster_gdal
    )
    out_tif
  }

  # ------------------------------------------------------------------
  # Dataset processors
  # ------------------------------------------------------------------
  process_preloaded <- function(tag, x, aoi_sfc, burn, year_out_dir,
                                id_col = ".__id",
                                erase_if_intersects_sfc = NULL,
                                protect_if_intersects_sfc = NULL,
                                template_r_for_rasterize = NULL) {

    x <- fix_sf(x, epsg = target_epsg)
    if (nrow(x) == 0) {
      if (on_empty == "skip") return(NULL)
      stop("0 polygons provided to process_preloaded for: ", tag)
    }

    msg("AOI clip: ", tag)
    x_aoi <- clip_to_aoi(x, aoi_sfc)
    x_aoi <- fix_sf(x_aoi, epsg = target_epsg)

    if (nrow(x_aoi) == 0) {
      if (on_empty == "skip") return(NULL)
      stop("After AOI clip, 0 polygons remain for: ", tag)
    }

    # SAFE erase rule at AOI stage
    n_erased_aoi <- 0L
    if (!is.null(erase_if_intersects_sfc) && length(erase_if_intersects_sfc) > 0) {
      erase_if_intersects_sfc <- normalize_sfc1(erase_if_intersects_sfc, crs = sf::st_crs(x_aoi))
      idx_drop <- lengths(sf::st_intersects(x_aoi, erase_if_intersects_sfc)) > 0

      idx_keep <- rep(FALSE, nrow(x_aoi))
      if (!is.null(protect_if_intersects_sfc) && length(protect_if_intersects_sfc) > 0) {
        protect_if_intersects_sfc <- normalize_sfc1(protect_if_intersects_sfc, crs = sf::st_crs(x_aoi))
        idx_keep <- lengths(sf::st_intersects(x_aoi, protect_if_intersects_sfc)) > 0
      }

      idx_remove <- idx_drop & !idx_keep
      n_erased_aoi <- sum(idx_remove, na.rm = TRUE)

      if (n_erased_aoi > 0) x_aoi <- x_aoi[!idx_remove, , drop = FALSE]
      x_aoi <- fix_sf(x_aoi, epsg = target_epsg)

      if (nrow(x_aoi) == 0) {
        if (on_empty == "skip") return(NULL)
        stop("After SAFE erase (AOI stage), 0 polygons remain for: ", tag)
      }
    }

    if (!id_col %in% names(x_aoi)) x_aoi[[id_col]] <- seq_len(nrow(x_aoi))

    pref <- if (!is.null(output_prefix) && nzchar(output_prefix)) .safe_slug(output_prefix, 40) else ""

    aoi_out2 <- .outfile("polys_aoi", year_out_dir, pref, tag = tag)
    if (isTRUE(save_aoi_only)) {
      .safe_overwrite(aoi_out2)
      st_write_gpkg_arcgis(x_aoi, aoi_out2, layer = out_layer_polys, delete_dsn = TRUE, quiet = TRUE)
    }

    x_burn <- NULL
    ratio_df <- NULL
    burn_out <- NA_character_
    ratio_out <- NA_character_
    rast_out  <- NA_character_
    n_erased_burn <- 0L

    if (isTRUE(do_burnable_clip)) {
      msg("Burnable clip (raster-based): ", tag)
      x_burn <- clip_to_burnable_raster(x_aoi, burn, dissolve_by_id = dissolve_by_id, id_col = id_col)
      x_burn <- fix_sf(x_burn, epsg = target_epsg)

      if (nrow(x_burn) == 0) {
        if (on_empty == "skip") return(NULL)
        stop("After burnable clip, 0 polygons remain for: ", tag)
      }

      # SAFE erase again after burnable clip (defensive)
      if (!is.null(erase_if_intersects_sfc) && length(erase_if_intersects_sfc) > 0) {
        erase_if_intersects_sfc <- normalize_sfc1(erase_if_intersects_sfc, crs = sf::st_crs(x_burn))
        idx_drop2 <- lengths(sf::st_intersects(x_burn, erase_if_intersects_sfc)) > 0

        idx_keep2 <- rep(FALSE, nrow(x_burn))
        if (!is.null(protect_if_intersects_sfc) && length(protect_if_intersects_sfc) > 0) {
          protect_if_intersects_sfc <- normalize_sfc1(protect_if_intersects_sfc, crs = sf::st_crs(x_burn))
          idx_keep2 <- lengths(sf::st_intersects(x_burn, protect_if_intersects_sfc)) > 0
        }

        idx_remove2 <- idx_drop2 & !idx_keep2
        n_erased_burn <- sum(idx_remove2, na.rm = TRUE)

        if (n_erased_burn > 0) x_burn <- x_burn[!idx_remove2, , drop = FALSE]
        x_burn <- fix_sf(x_burn, epsg = target_epsg)

        if (nrow(x_burn) == 0) {
          if (on_empty == "skip") return(NULL)
          stop("After SAFE erase (burnable stage), 0 polygons remain for: ", tag)
        }
      }

      burn_out <- .outfile("polys_burnable", year_out_dir, pref, tag = tag)
      if (isTRUE(save_burnable)) {
        .safe_overwrite(burn_out)
        st_write_gpkg_arcgis(x_burn, burn_out, layer = out_layer_polys, delete_dsn = TRUE, quiet = TRUE)
      }

      if (isTRUE(save_ratio_csv)) {
        ratio_df <- compute_ratio(x_aoi, x_burn, id_col = id_col)
        ratio_out <- .outfile("ratio", year_out_dir, pref, tag = tag)
        .safe_overwrite(ratio_out)
        utils::write.csv(ratio_df, ratio_out, row.names = FALSE)
      }
    }

    if (isTRUE(save_rasterized)) {
      stage_str <- if (isTRUE(do_burnable_clip)) "AOI_BURNABLE" else "AOI"
      out_tif <- .outfile("raster", year_out_dir, pref, tag = tag, stage = stage_str)
      rast_out <- rasterize_polys(
        if (isTRUE(do_burnable_clip)) x_burn else x_aoi,
        template_r_for_rasterize,
        out_tif
      )
    }

    list(
      aoi_path      = if (isTRUE(save_aoi_only)) aoi_out2 else NA_character_,
      burnable_path = if (isTRUE(save_burnable) && isTRUE(do_burnable_clip)) burn_out else NA_character_,
      ratio_path    = if (isTRUE(save_ratio_csv) && isTRUE(do_burnable_clip)) ratio_out else NA_character_,
      raster_path   = if (isTRUE(save_rasterized)) rast_out else NA_character_,
      aoi_sf        = x_aoi,
      burnable_sf   = x_burn,
      ratio_df      = ratio_df,
      n_erased_aoi  = as.integer(n_erased_aoi),
      n_erased_burn = as.integer(n_erased_burn)
    )
  }

  process_one <- function(tag, files, in_layer, aoi_sfc, burn, year_out_dir,
                          filter_fun = NULL, id_col = ".__id",
                          erase_if_intersects_sfc = NULL,
                          protect_if_intersects_sfc = NULL,
                          template_r_for_rasterize = NULL) {

    msg("Reading: ", tag, " | n_files=", length(files))

    xs <- lapply(files, function(fp) {
      read_sf_any(fp, layer = in_layer)
    })

    x <- rbind_sf_safe(xs, epsg = target_epsg)

    if (nrow(x) == 0) {
      if (on_empty == "skip") return(NULL)
      stop("0 polygons after reading for: ", tag)
    }

    if (!id_col %in% names(x)) x[[id_col]] <- seq_len(nrow(x))

    if (!is.null(filter_fun)) {
      x <- filter_fun(x)
      x <- fix_sf(x, epsg = target_epsg)
      if (nrow(x) == 0) {
        if (on_empty == "skip") return(NULL)
        stop("After filtering, 0 polygons remain for: ", tag)
      }
    }

    process_preloaded(
      tag = tag,
      x = x,
      aoi_sfc = aoi_sfc,
      burn = burn,
      year_out_dir = year_out_dir,
      id_col = id_col,
      erase_if_intersects_sfc = erase_if_intersects_sfc,
      protect_if_intersects_sfc = protect_if_intersects_sfc,
      template_r_for_rasterize = template_r_for_rasterize
    )
  }

  process_reference_with_drops <- function(tag, files, in_layer, aoi_sfc, burn, year_out_dir,
                                           filter_fun = NULL, id_col = ".__id",
                                           template_r_for_rasterize = NULL) {

    msg("Reading (reference): ", tag, " | n_files=", length(files))

    xs <- lapply(files, function(fp) {
      read_sf_any(fp, layer = in_layer)
    })

    x_all <- rbind_sf_safe(xs, epsg = target_epsg)

    if (nrow(x_all) == 0) {
      if (on_empty == "skip") return(NULL)
      stop("0 polygons after reading for: ", tag)
    }

    if (!id_col %in% names(x_all)) x_all[[id_col]] <- seq_len(nrow(x_all))

    # Stage A: filters
    x_filt_keep <- x_all
    x_filt_drop <- x_all[0, , drop = FALSE]

    if (!is.null(filter_fun)) {
      x_filt_keep <- filter_fun(x_all)
      x_filt_keep <- fix_sf(x_filt_keep, epsg = target_epsg)

      if (nrow(x_filt_keep) == 0) {
        if (on_empty == "skip") return(NULL)
        stop("After filtering, 0 polygons remain for: ", tag)
      }

      keep_ids <- x_filt_keep[[id_col]]
      x_filt_drop <- x_all[!(x_all[[id_col]] %in% keep_ids), , drop = FALSE]
      x_filt_drop <- fix_sf(x_filt_drop, epsg = target_epsg)
    }

    dropped_by_filter_n <- nrow(x_filt_drop)

    # Stage B: AOI keep
    x_keep_aoi <- clip_to_aoi(x_filt_keep, aoi_sfc)
    x_keep_aoi <- fix_sf(x_keep_aoi, epsg = target_epsg)

    ids_keep_aoi <- unique(x_keep_aoi[[id_col]])
    x_drop_by_aoi <- x_filt_keep[!(x_filt_keep[[id_col]] %in% ids_keep_aoi), , drop = FALSE]
    x_drop_by_aoi <- fix_sf(x_drop_by_aoi, epsg = target_epsg)
    dropped_by_aoi_n <- nrow(x_drop_by_aoi)

    # Stage C: burnable keep (final ref)
    x_keep_final <- x_keep_aoi
    x_drop_by_burn <- x_keep_aoi[0, , drop = FALSE]
    dropped_by_burn_n <- 0L

    if (isTRUE(do_burnable_clip)) {
      x_keep_final <- clip_to_burnable_raster(x_keep_aoi, burn, dissolve_by_id = dissolve_by_id, id_col = id_col)
      x_keep_final <- fix_sf(x_keep_final, epsg = target_epsg)

      ids_keep_burn <- unique(x_keep_final[[id_col]])
      x_drop_by_burn <- x_keep_aoi[!(x_keep_aoi[[id_col]] %in% ids_keep_burn), , drop = FALSE]
      x_drop_by_burn <- fix_sf(x_drop_by_burn, epsg = target_epsg)
      dropped_by_burn_n <- nrow(x_drop_by_burn)
    }

    # Build DROP ANY mask (in AOI space)
    x_filt_drop_stage <- x_filt_drop
    if (nrow(x_filt_drop_stage) > 0) {
      x_filt_drop_stage <- clip_to_aoi(x_filt_drop_stage, aoi_sfc)
      x_filt_drop_stage <- fix_sf(x_filt_drop_stage, epsg = target_epsg)

      if (dropped_ref_stage == "burnable" && isTRUE(do_burnable_clip) && nrow(x_filt_drop_stage) > 0) {
        x_filt_drop_stage <- clip_to_burnable_raster(
          x_filt_drop_stage, burn,
          dissolve_by_id = dissolve_by_id, id_col = id_col
        )
        x_filt_drop_stage <- fix_sf(x_filt_drop_stage, epsg = target_epsg)
      }
    }

    drop_any <- rbind_sf_safe(list(x_filt_drop_stage, x_drop_by_burn), epsg = target_epsg)
    drop_any <- fix_sf(drop_any, epsg = target_epsg)

    dropped_any_n   <- 0L
    dropped_stage_n <- 0L
    if (nrow(drop_any) > 0) {
      dropped_any_n   <- length(unique(drop_any[[id_col]]))
      dropped_stage_n <- nrow(drop_any)
    }

    drop_any_sfc <- NULL
    if (isTRUE(remove_detected_if_intersects_dropped_ref) && nrow(drop_any) > 0) {
      drop_any_sfc <- normalize_sfc1(sf::st_union(sf::st_geometry(drop_any)), crs = sf::st_crs(drop_any))
      if (is.finite(dropped_ref_buffer_m) && dropped_ref_buffer_m > 0) {
        drop_any_sfc <- normalize_sfc1(sf::st_buffer(drop_any_sfc, dropped_ref_buffer_m), crs = sf::st_crs(drop_any_sfc))
      }
    }

    keep_final_sfc <- NULL
    if (nrow(x_keep_final) > 0) {
      keep_final_sfc <- normalize_sfc1(sf::st_union(sf::st_geometry(x_keep_final)), crs = sf::st_crs(x_keep_final))
    }

    keep_res <- process_preloaded(
      tag = tag,
      x = x_filt_keep,
      aoi_sfc = aoi_sfc,
      burn = burn,
      year_out_dir = year_out_dir,
      id_col = id_col,
      erase_if_intersects_sfc = NULL,
      protect_if_intersects_sfc = NULL,
      template_r_for_rasterize = template_r_for_rasterize
    )

    if (is.null(keep_res)) return(NULL)

    keep_res$drop_any_sfc     <- drop_any_sfc
    keep_res$keep_final_sfc   <- keep_final_sfc

    keep_res$dropped_any_n       <- as.integer(dropped_any_n)
    keep_res$dropped_stage_n     <- as.integer(dropped_stage_n)
    keep_res$dropped_by_filter_n <- as.integer(dropped_by_filter_n)
    keep_res$dropped_by_aoi_n    <- as.integer(dropped_by_aoi_n)
    keep_res$dropped_by_burn_n   <- as.integer(dropped_by_burn_n)

    keep_res
  }

  # ------------------------------------------------------------------
  # Decide multi-year mode + expand inputs
  # ------------------------------------------------------------------
  ref_files <- .expand_paths(ref_path)
  det_files <- .expand_paths(det_path)

  auto_multi <- (length(ref_files) != 1L || length(det_files) != 1L ||
                   (length(ref_path) == 1L && dir.exists(ref_path)) ||
                   (length(det_path) == 1L && dir.exists(det_path)))

  if (is.null(multi_year)) multi_year <- auto_multi
  if (!isTRUE(multi_year)) {
    ref_files <- ref_files[1]
    det_files <- det_files[1]
  }

  # ------------------------------------------------------------------
  # S2 handling
  # ------------------------------------------------------------------
  if (isTRUE(disable_s2)) {
    old_s2 <- sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  }

  # ------------------------------------------------------------------
  # Year resolution
  # ------------------------------------------------------------------
  years_ref <- vapply(ref_files, .extract_year_from_path, integer(1), year_regex = year_regex)
  years_det <- vapply(det_files, .extract_year_from_path, integer(1), year_regex = year_regex)

  if (!is.null(year_target) && is.finite(as.integer(year_target))) {
    years_to_process <- as.integer(year_target)
  } else if (isTRUE(multi_year)) {
    yr_ref_u <- sort(unique(years_ref[!is.na(years_ref)]))
    yr_det_u <- sort(unique(years_det[!is.na(years_det)]))

    if (length(yr_ref_u) == 0 || length(yr_det_u) == 0) {
      stop("Multi-year mode requires years in filenames (YYYY). Provide year_target or rename files.")
    }

    years_to_process <- sort(intersect(yr_ref_u, yr_det_u))
    if (length(years_to_process) == 0) {
      stop("No common years between ref and det. Provide year_target or ensure matching YYYY in filenames.")
    }
  } else {
    y1 <- years_ref[1]
    y2 <- years_det[1]
    yy <- if (!is.na(y1)) y1 else y2
    if (is.na(yy)) yy <- NA_integer_
    years_to_process <- yy
  }

  # ------------------------------------------------------------------
  # Run per year
  # ------------------------------------------------------------------
  out_years <- list()

  for (yy in years_to_process) {

    year_str <- if (is.na(yy)) "NA" else as.character(yy)

    year_out_dir <- out_dir
    if (isTRUE(multi_year) && isTRUE(use_year_subdirs) && !is.na(yy)) {
      year_out_dir <- file.path(out_dir, year_str)
      dir.create(year_out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    msg("---- Year: ", year_str, " | out: ", year_out_dir)

    if (isTRUE(multi_year) && !is.na(yy)) {
      ref_sel <- ref_files[years_ref == yy]
      det_sel <- det_files[years_det == yy]
    } else {
      ref_sel <- ref_files
      det_sel <- det_files
    }

    if (length(ref_sel) == 0 || length(det_sel) == 0) {
      if (on_empty == "skip") next
      stop("Missing ref or det files for year=", year_str)
    }

    # ------------------------------------------------------------------
    # AOI selection
    # ------------------------------------------------------------------
    # Rule:
    #   - years < aoi_latest_from_year:
    #       build AOI from fire polygon origin
    #   - years >= aoi_latest_from_year:
    #       use fixed AOI shapefile or AOI mask picker
    # ------------------------------------------------------------------

    use_origin_this_year <- isTRUE(use_origin_aoi) &&
      !is.na(yy) &&
      as.integer(yy) < as.integer(origin_aoi_before_year)

    if (isTRUE(use_origin_this_year)) {

      msg(
        "Using origin-based AOI because year ",
        yy,
        " < ",
        aoi_latest_from_year
      )

      aoi_pack <- .load_aoi_from_origin_for_year(
        ref_files_this_year = ref_sel,
        ref_layer_this_year = ref_layer
      )

    } else {

      msg(
        "Using standard AOI because year ",
        yy,
        " >= ",
        aoi_latest_from_year
      )

      if (!is.null(aoi_shapefile) && nzchar(aoi_shapefile)) {

        if (!file.exists(aoi_shapefile)) {
          stop("AOI shapefile not found: ", aoi_shapefile)
        }

        aoi <- sf::st_read(aoi_shapefile, quiet = TRUE)

        if (is.na(sf::st_crs(aoi))) {
          stop("AOI shapefile has NA CRS: ", aoi_shapefile)
        }

        aoi <- fix_sf(
          sf::st_transform(sf::st_make_valid(aoi), target_epsg),
          epsg = target_epsg
        )

        aoi_sfc <- normalize_sfc1(
          sf::st_geometry(aoi),
          crs = sf::st_crs(aoi)
        )

        aoi_pack <- list(
          aoi_sf = aoi,
          aoi_sfc = aoi_sfc,
          aoi_tif = NA_character_,
          origin_raw = character(0),
          origin_match = data.frame(),
          origin_spatial_match = data.frame()
        )

      } else {

        aoi_pack <- load_aoi_for_year(yy)
      }
    }

    aoi <- aoi_pack$aoi_sf
    aoi_sfc <- normalize_sfc1(aoi_pack$aoi_sfc, crs = sf::st_crs(aoi))
    aoi_tif_path <- aoi_pack$aoi_tif

    pref <- if (!is.null(output_prefix) && nzchar(output_prefix)) .safe_slug(output_prefix, 40) else ""

    # Save AOI (optional)
    aoi_out <- .outfile("aoi_mask", year_out_dir, pref, epsg = target_epsg)
    if (isTRUE(save_aoi_mask)) {
      .safe_overwrite(aoi_out)
      st_write_gpkg_arcgis(
        sf::st_sf(geometry = sf::st_geometry(aoi)),
        aoi_out,
        layer = out_layer_aoi,
        delete_dsn = TRUE,
        quiet = TRUE
      )
    } else {
      aoi_out <- NA_character_
    }

    # Burnable raster ONLY if needed
    burn <- NULL
    needs_burn <- isTRUE(do_burnable_clip) ||
      (isTRUE(save_rasterized) && identical(raster_template, "burnable_raster")) ||
      (isTRUE(remove_detected_if_intersects_dropped_ref) && identical(dropped_ref_stage, "burnable"))

    if (isTRUE(needs_burn)) {
      burn <- load_burnable_for_year(yy, aoi)
    } else {
      msg("Burnable mask not required for this run (do_burnable_clip=FALSE and no burn-based options requested).")
    }

    # Rasterization template
    template_r_for_rasterize <- NULL
    if (isTRUE(save_rasterized)) {
      if (raster_template == "burnable_raster") {
        if (is.null(burn)) stop("save_rasterized=TRUE with raster_template='burnable_raster' requires a burnable raster.")
        template_r_for_rasterize <- burn
      } else {
        if (is.character(aoi_tif_path) && nzchar(aoi_tif_path) && file.exists(aoi_tif_path)) {
          template_r_for_rasterize <- terra::rast(aoi_tif_path)
          if (terra::nlyr(template_r_for_rasterize) > 1) template_r_for_rasterize <- template_r_for_rasterize[[1]]
        } else {
          stop("save_rasterized=TRUE with raster_template='aoi_tif' requires AOI tif to exist. Check aoi_mask_dir.")
        }
      }
    }

    # Tags
    if (is.null(ref_output_tag) || !nzchar(ref_output_tag)) {
      ref_output_tag2 <- paste0("REF_", year_str)
      if (isTRUE(apply_ref_doy_filter)) ref_output_tag2 <- paste0(ref_output_tag2, "_D", ref_doy_min, "_", ref_doy_max)
      if (isTRUE(apply_ref_attr_filter)) ref_output_tag2 <- paste0(ref_output_tag2, "_FILT")
    } else {
      ref_output_tag2 <- ref_output_tag
    }

    if (is.null(det_output_tag) || !nzchar(det_output_tag)) {
      det_output_tag2 <- paste0("DET_", year_str, if (isTRUE(apply_det_filter)) "_FILT" else "_ALL")
    } else {
      det_output_tag2 <- det_output_tag
    }

    ref_tag <- .safe_slug(ref_output_tag2, 90)
    det_tag <- .safe_slug(det_output_tag2, 90)

    # Reference with drops
    ref_res <- process_reference_with_drops(
      tag = ref_tag,
      files = ref_sel,
      in_layer = ref_layer,
      aoi_sfc = aoi_sfc,
      burn = burn,
      year_out_dir = year_out_dir,
      filter_fun = ref_filter_fun,
      id_col = ".__id",
      template_r_for_rasterize = template_r_for_rasterize
    )

    if (is.null(ref_res)) {
      if (on_empty == "skip") next
      stop("Empty reference results for year=", year_str)
    }

    erase_sfc   <- NULL
    protect_sfc <- NULL

    if (isTRUE(remove_detected_if_intersects_dropped_ref)) {
      erase_sfc   <- ref_res$drop_any_sfc
      protect_sfc <- ref_res$keep_final_sfc

      msg(
        "Dropped refs: any_ids=", ref_res$dropped_any_n,
        " | by_filter=", ref_res$dropped_by_filter_n,
        " | by_aoi=", ref_res$dropped_by_aoi_n,
        " | by_burn=", ref_res$dropped_by_burn_n,
        " | drop_geom_n=", ref_res$dropped_stage_n,
        " | SAFE erase: remove = intersects(drop_any) & !intersects(keep_final)"
      )
    }

    # Detected
    det_res <- process_one(
      tag = det_tag,
      files = det_sel,
      in_layer = det_layer,
      aoi_sfc = aoi_sfc,
      burn = burn,
      year_out_dir = year_out_dir,
      filter_fun = det_filter_fun,
      id_col = ".__id",
      erase_if_intersects_sfc   = erase_sfc,
      protect_if_intersects_sfc = protect_sfc,
      template_r_for_rasterize  = template_r_for_rasterize
    )

    if (is.null(det_res)) {
      if (on_empty == "skip") next
      stop("Empty detected results for year=", year_str)
    }

    out_years[[year_str]] <- list(
      out_dir = year_out_dir,
      target_epsg = target_epsg,
      aoi_gpkg = aoi_out,
      reference = list(
        tag = ref_tag,
        aoi_path = ref_res$aoi_path,
        burnable_path = ref_res$burnable_path,
        ratio_path = ref_res$ratio_path,
        raster_path = ref_res$raster_path,
        dropped_any_n = ref_res$dropped_any_n,
        dropped_by_filter_n = ref_res$dropped_by_filter_n,
        dropped_by_aoi_n = ref_res$dropped_by_aoi_n,
        dropped_by_burn_n = ref_res$dropped_by_burn_n,
        drop_geom_n = ref_res$dropped_stage_n
      ),
      detected = list(
        tag = det_tag,
        aoi_path = det_res$aoi_path,
        burnable_path = det_res$burnable_path,
        ratio_path = det_res$ratio_path,
        raster_path = det_res$raster_path,
        erased_by_droppedref_aoi_n = det_res$n_erased_aoi,
        erased_by_droppedref_burn_n = det_res$n_erased_burn
      )
    )
  }

  if (length(out_years) == 1) {
    out <- out_years[[1]]
    out$years <- out_years
    return(out)
  }

  list(
    out_dir = out_dir,
    target_epsg = target_epsg,
    years = out_years
  )
}

utils::globalVariables(c(
  ".data",
  "area_before_ha",
  "area_after_ha",
  "burnable_ratio",
  "nonburn_ratio"
))
