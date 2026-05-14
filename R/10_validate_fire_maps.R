#' Validate fire detections against a reference dataset (pixel, polygon/area, and optional classwise metrics)
#'
#' @description
#' Validates a set of detected (predicted) burned-area polygons against a reference dataset
#' within an AOI mask. The function rasterizes both reference and detections on a common
#' analysis grid (at \code{template_res}) and computes:
#' \itemize{
#'   \item Pixel/grid confusion-matrix metrics (TP/FP/FN/TN and derived scores)
#'   \item Polygon/area overlap summaries (coverage-based detection of reference polygons)
#'   \item Optional class-wise pixel metrics using a categorical raster (e.g., CORINE strata)
#'   \item Optional omission/commission summaries by class (via unmatched polygons)
#'   \item Optional confusion outputs: a pixel-coded confusion raster and/or polygon confusion vectors
#' }
#'
#' @details
#' Detailed behavior is organized in the sections below.
#' @section Pixel (grid) validation (confusion matrix):
#' The AOI mask polygon is rasterized to an analysis grid at \code{template_res} meters.
#' Reference and detected polygons are rasterized on the same grid, producing two binary rasters.
#' Confusion-matrix counts are computed over the grid domain.
#'
#' Main metrics reported:
#' \itemize{
#'   \item \strong{Precision} = TP / (TP + FP)
#'   \item \strong{Recall} = TP / (TP + FN)
#'   \item \strong{F1} = 2 * Precision * Recall / (Precision + Recall)
#'   \item \strong{IoU} = TP / (TP + FP + FN)
#'   \item \strong{ErrorRate} = (FP + FN) / (TP + FP + FN + TN)
#'   \item \strong{Commission} = 1 - Precision
#'   \item \strong{Omission} = 1 - Recall
#' }
#'
#' @section Polygon / area validation:
#' Polygon-level detection is based on coverage of each reference polygon by detections:
#' \eqn{Coverage_i = 100 * Area(ref_i \cap det) / Area(ref_i)}{Coverage_i = 100*A(intersection)/A(ref)}.
#' A reference polygon is detected when \code{Coverage_i >= min_detected_percent}, and completely
#' detected when \code{Coverage_i >= threshold_completely_detected}.
#'
#' @section Dissolve behavior:
#' Reference and detection polygons are clipped to the mask. Optional post-clip dissolve can be applied
#' with \code{dissolve_ref_by} and \code{dissolve_input_by} to recombine fragments by ID.
#'
#' @section Confusion outputs (optional):
#' If \code{write_confusion_pixel_raster=TRUE}, a pixel confusion raster is written:
#' \itemize{
#'   \item 1 = TN
#'   \item 2 = FN
#'   \item 3 = FP
#'   \item 4 = TP
#' }
#' If \code{write_confusion_vectors=TRUE}, a GeoPackage with TP/FP/FN polygon labels is written.
#'
#' @section Class-wise metrics (optional):
#' If \code{class_raster} is provided, confusion-matrix metrics are also computed by class.
#' Optional class-wise omission/commission summaries are written when
#' \code{write_classwise_errors=TRUE}.
#'
#' @section Caching and reproducibility:
#' Reference polygons and rasterized reference masks are cached inside the output folder.
#' Use \code{force_reprocess_ref=TRUE} to force a rebuild.
#'
#' @section ArcGIS-friendly outputs:
#' When \code{arcgis_wkt1=TRUE}, GeoPackages are written with WKT1_ESRI CRS serialization
#' and via a temporary-write/rename strategy to reduce lock issues on Windows/ArcGIS.
#'
#' @param input_shapefile Path(s) to detected polygons (vector), or an \code{sf} object.
#'   If multiple paths are provided, outputs are indexed by \code{InputIndex}.
#' @param ref_shapefile Path to reference polygons, or an \code{sf} object.
#' @param mask_shapefile Path to AOI mask polygon (sf-readable).
#' @param year_target Integer. Target year to filter reference polygons if a year-like field exists
#'   (first match among \code{year/Year/YEAR}).
#' @param validation_dir Character. Output base directory.
#'
#' @param template_res Numeric. Grid resolution (meters) for pixel metrics and for rasterizing classwise errors.
#' @param min_detected_percent Numeric. Reference polygon considered detected if \code{Coverage_i >=} this (\%).
#' @param threshold_completely_detected Numeric. Reference polygon considered completely detected if \code{Coverage_i >=} this (\%).
#'
#' @param min_area_reference_ha Numeric or NULL. Minimum reference polygon area to retain (ha). If NULL, no filter.
#' @param force_reprocess_ref Logical. If TRUE, delete cached reference products and rebuild them.
#' @param force_reprocess_pred Logical. If TRUE, rebuild prediction rasters (per input).
#'
#' @param metrics_type One of \code{c("all","pixel","area")}. Controls which metric blocks to compute.
#' @param dissolve_ref_by Character or NULL. Optional dissolve field for reference polygons \strong{after clipping to the mask}.
#'   Use this to recombine clipped fragments that share the same ID (e.g., multiple parts of the same fire perimeter).
#'   Dissolving aggregates non-geometry attributes within each ID (see \strong{Dissolve behavior} section).
#'   Set to NULL to skip dissolve and preserve clipped fragments as separate rows.
#' @param dissolve_input_by Character or NULL. Optional dissolve field for detection polygons \strong{after clipping to the mask}.
#'   Use this to keep \strong{one feature per original polygon ID} when the mask clip fragments geometries.
#'   Set to NULL to preserve clipped fragments (and avoid attribute aggregation).
#'
#' @param burned_value Numeric/integer. Pixel value representing "burned" in rasterized prediction/reference (default 1).
#' @param pixel_domain One of \code{c("mask","complete_cases")}. Kept for call compatibility; current implementation
#'   evaluates over the mask domain raster (\code{domain_mask}).
#' @param pixel_na_as_zero Logical. Kept for call compatibility; reference/pred are 0/1 inside mask and NA outside.
#' @param strict_binary_check Logical. If TRUE, stop if the prediction raster contains values other than 0/1 (excluding NA).
#'
#' @param class_raster Optional categorical raster for class-wise metrics (path or \code{SpatRaster}).
#' @param class_band Integer. Band index if \code{class_raster} has multiple layers.
#' @param class_values Optional numeric vector. Subset of class codes to keep (others set NA).
#' @param class_lut Optional LUT (data.frame or CSV path). Must define \code{id} and \code{label} columns (or first two columns are used).
#' @param class_domain One of \code{c("class","complete_cases")}. Domain definition for class-wise metrics.
#' @param class_na_as_zero Logical. If TRUE, treat NA in pred/ref as 0 within the classwise domain.
#' @param class_chunk_rows Integer. Chunk size (rows) for class-wise raster streaming (memory/performance control).
#' @param class_min_ref_area_ha Numeric. Minimum reference area per class (ha) to report in classwise metrics.
#' @param class_include_all_lut_ids Logical. If TRUE (and LUT provided), include LUT ids even if absent in raster (filled as 0/NA as appropriate).
#'
#' @param class_shape Deprecated/optional. Vector strata layer (kept for backward compatibility; currently not used).
#' @param class_field Deprecated/optional. Field in \code{class_shape} defining strata (currently not used).
#'
#' @param write_classwise_metrics Logical. If TRUE, write class-wise metrics CSV(s).
#' @param write_classwise_errors Logical. If TRUE, compute omission/commission area by class (CSV outputs).
#'
#' @param input_layer Optional layer name for GeoPackage detections (\code{input_shapefile}).
#'   If NULL and input is a .gpkg, the function auto-picks a preferred layer when available.
#' @param ref_layer Optional layer name for GeoPackage reference (\code{ref_shapefile}).
#'   If NULL and reference is a .gpkg, the function auto-picks a preferred layer when available.
#'
#' @param union_detections_for_area Logical. If TRUE, union detections before area intersection (avoids overlap double-counting).
#' @param arcgis_wkt1 Logical. If TRUE, write GeoPackage outputs with WKT1_ESRI CRS serialization (ArcGIS-friendly).
#'
#' @param output_prefix Character. Short prefix used in filenames.
#' @param output_tag Optional character tag used in folder naming (if \code{tag_as_subdir=TRUE}).
#' @param tag_as_subdir Logical. If TRUE, write outputs into a subfolder \code{prefix_tag} under \code{outputs_subdir}.
#' @param outputs_subdir Character. Subfolder name under \code{validation_dir} for outputs.
#' @param errors_subdir Character. Subfolder name under outputs for error vectors (when enabled).
#' @param write_error_vectors Logical. If TRUE, write omission/commission error vectors (GeoPackage).
#' @param index_width Integer. Width for indexed filenames if multiple inputs are validated.
#'
#' @param digits_metrics Integer. Decimal digits for ratio/percent metrics (e.g., Precision, Recall, F1, IoU, ErrorRate, area recall/precision).
#' @param digits_area Integer. Decimal digits for area fields in hectares.
#' @param digits_res Integer. Decimal digits for resolution-related fields (e.g., TemplateRes_m).
#'
#' @param run_tag Optional character metadata written into output tables.
#' @param ref_tag Optional character metadata written into output tables.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @param write_xlsx Logical. If TRUE, write a summary workbook (requires \pkg{openxlsx}).
#' @param xlsx_path Optional character. Explicit path for the workbook; if NULL, a default name is used in \code{validation_output_dir}.
#' @param xlsx_overwrite Logical. If TRUE, overwrite existing workbook at \code{xlsx_path}.
#'
#' @param write_confusion_vectors Logical. If TRUE, write a GeoPackage with reference + detection polygons labeled TP/FP/FN.
#' @param write_confusion_pixel_raster Logical. If TRUE, write a pixel confusion GeoTIFF coded as 1=TN, 2=FN, 3=FP, 4=TP.
#' @param confusion_subdir Character. Subfolder name under \code{validation_output_dir} for confusion outputs.
#'   If NULL/empty/NA, it defaults to \code{"CONF"}.
#' @param confusion_polys_layer Character. Layer name for the confusion polygon GeoPackage.
#' @param confusion_polys_suffix Character. Filename suffix for confusion polygons (used by the internal naming helper).
#' @param confusion_pixels_suffix Character. Filename suffix for confusion pixels raster (used by the internal naming helper).
#' @param compute_det_overlap_stats Logical. If TRUE, add overlap area and percent (vs reference union) to detection confusion polygons.
#'
#' @param ... Reserved for future use.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{metrics}: \code{data.table} of global pixel metrics (one row per input) or NULL if \code{metrics_type="area"}.
#'   \item \code{polygon_summary}: \code{data.table} of polygon/area summaries (one row per input) or NULL if \code{metrics_type="pixel"}.
#'   \item \code{validation_output_dir}: output folder used.
#'   \item \code{errors_dir}: error vectors folder (created even if empty; depends on options).
#'   \item \code{xlsx_path}: workbook path (if written), otherwise NULL.
#'   \item \code{output_prefix}: sanitized prefix used in filenames.
#'   \item \code{output_tag}: sanitized tag used in folder naming (empty string if none).
#'   \item \code{confusion_polys_paths}: list of GeoPackage paths (per input; NULL entries if not written).
#'   \item \code{confusion_pixels_paths}: list of GeoTIFF paths (per input; NULL entries if not written).
#' }
#'
#' @section Internal helper functions (defined inside the function):
#' \itemize{
#'   \item \code{.safe_slug}, \code{short_id}, \code{round_outputs}, \code{make_valid_sf}
#'   \item \code{transform_to_crs_safe}, \code{dissolve_by_field}, \code{safe_vect}, \code{read_sf_any}
#'   \item \code{sanitize_for_gpkg}, \code{write_gpkg_safe}, \code{out_file}
#'   \item \code{terra_freq_noNA}, \code{classwise_pixel_metrics}
#' }
#'
#' @examples
#' \dontrun{
#' # Single input (GPKG) + classwise CORINE + confusion outputs
#' res <- validate_fire_maps(
#'   input_shapefile = "D:/PROJECT/DET/internal_scored.gpkg",
#'   ref_shapefile   = "D:/PROJECT/REF/effis_ref.gpkg",
#'   mask_shapefile  = "D:/PROJECT/MASK/aoi_3035.shp",
#'   year_target     = 2025,
#'   validation_dir  = "D:/PROJECT/VALIDATION",
#'   input_layer     = "scored",
#'   ref_layer       = "polys",
#'   metrics_type    = "all",
#'   template_res    = 90,
#'   min_detected_percent = 10,
#'   threshold_completely_detected = 50,
#'   min_area_reference_ha = 1,
#'   force_reprocess_ref  = TRUE,
#'   force_reprocess_pred = TRUE,
#'   union_detections_for_area = TRUE,
#'   arcgis_wkt1 = TRUE,
#'   output_prefix = "VAL",
#'   output_tag    = "INT_2025_ALL",
#'   tag_as_subdir = TRUE,
#'   outputs_subdir = "VAL",
#'   errors_subdir  = "ERR",
#'   write_error_vectors = TRUE,
#'   class_raster = "D:/PROJECT/MASK/corine_groups.tif",
#'   class_band   = 1,
#'   class_lut    = "D:/PROJECT/MASK/corine_groups_lut.csv",
#'   write_classwise_metrics = TRUE,
#'   write_classwise_errors  = TRUE,
#'   write_confusion_vectors = TRUE,
#'   write_confusion_pixel_raster = TRUE,
#'   compute_det_overlap_stats = TRUE,
#'   digits_metrics = 3,
#'   digits_area    = 2,
#'   digits_res     = 0,
#'   write_xlsx = TRUE
#' )
#'
#' # Multiple detection inputs (vector of paths); outputs are indexed by InputIndex
#' res_multi <- validate_fire_maps(
#'   input_shapefile = c("D:/PROJECT/DET/runA.gpkg", "D:/PROJECT/DET/runB.gpkg"),
#'   ref_shapefile   = "D:/PROJECT/REF/effis_ref.gpkg",
#'   mask_shapefile  = "D:/PROJECT/MASK/aoi_3035.shp",
#'   year_target     = 2025,
#'   validation_dir  = "D:/PROJECT/VALIDATION",
#'   metrics_type    = "all",
#'   template_res    = 90
#' )
#' }
#'
#' @importFrom sf st_read st_layers st_as_sf st_zm st_make_valid st_is_valid st_is_empty
#' @importFrom sf st_collection_extract st_geometry_type st_geometry st_cast st_crs st_transform
#' @importFrom sf st_filter st_intersects st_intersection st_area st_union st_sfc sf_use_s2
#' @importFrom terra rast vect ext crs res nrow values global rasterize mask freq
#' @importFrom terra compareGeom same.crs project resample ifel writeRaster readStart readStop readValues
#' @importFrom data.table as.data.table data.table setDT set setnames setcolorder setattr rbindlist fread fwrite
#' @importFrom dplyr filter group_by summarise bind_rows
#' @importFrom rlang .data
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils write.table
#' @export
validate_fire_maps <- function(
    input_shapefile,
    ref_shapefile,              # path or sf
    mask_shapefile,             # path (sf)
    year_target,
    validation_dir,

    # --- grid / domain for pixel metrics ---
    template_res = 100,         # meters

    # --- overlap criteria (polygon) ---
    min_detected_percent = 10,
    threshold_completely_detected = 90,

    # --- optional filters ---
    min_area_reference_ha = NULL,
    force_reprocess_ref = FALSE,
    force_reprocess_pred = FALSE,

    metrics_type = c("all", "pixel", "area"),
    dissolve_ref_by = NULL,
    dissolve_input_by = NULL,

    # --- pixel behavior (kept to match your calls) ---
    burned_value = 1,
    pixel_domain = c("mask", "complete_cases"),
    pixel_na_as_zero = TRUE,
    strict_binary_check = FALSE,

    # --- optional classwise strata (raster OR vector) ---
    class_raster = NULL,          # path or SpatRaster
    class_band = 1,
    class_values = NULL,          # optional subset (numeric vector)
    class_lut = NULL,             # data.frame with columns id,label (or path to CSV)
    class_domain = c("class", "complete_cases"),
    class_na_as_zero = TRUE,
    class_chunk_rows = 2048,
    class_min_ref_area_ha = 0,
    class_include_all_lut_ids = FALSE,

    # optional vector strata (kept for backward compatibility)
    class_shape = NULL,
    class_field = NULL,

    write_classwise_metrics = FALSE,
    write_classwise_errors  = FALSE,

    # --- layers for GPKG ---
    input_layer = NULL,
    ref_layer   = NULL,

    # --- avoid double counting (recommended TRUE) ---
    union_detections_for_area = TRUE,

    # --- ArcGIS-friendly outputs ---
    arcgis_wkt1 = TRUE,

    # --- OUTPUT NAMING (short) ---
    output_prefix = "VAL",
    output_tag = NULL,
    tag_as_subdir = TRUE,
    outputs_subdir = "VAL",
    errors_subdir  = "ERR",
    write_error_vectors = TRUE,
    index_width = 2L,

    # --- rounding (optional) ---
    digits_metrics = 3,
    digits_area    = 2,
    digits_res     = 0,

    # --- tags kept ONLY as metadata ---
    run_tag = NULL,
    ref_tag = NULL,
    verbose = TRUE,

    # --- Excel output (optional) ---
    write_xlsx = FALSE,
    xlsx_path = NULL,
    xlsx_overwrite = TRUE,

    # --- CONFUSION OUTPUTS (optional) ---
    write_confusion_vectors = FALSE,        # writes 1 GPKG with ref+det polygons labeled TP/FP/FN (no TN for polygons)
    write_confusion_pixel_raster = FALSE,   # writes a TIFF with TN/FN/FP/TP per pixel
    confusion_subdir = "CONF",
    confusion_polys_layer = "confusion_polys",
    confusion_polys_suffix = "confusion_polys",
    confusion_pixels_suffix = "confusion_pixels",
    compute_det_overlap_stats = FALSE,      # optional extra columns for detections (overlap area/pct vs ref union)

    ...
) {

  metrics_type <- match.arg(metrics_type)
  class_domain <- match.arg(class_domain)
  pixel_domain <- match.arg(pixel_domain)

  suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(dplyr)
    library(data.table)
  })

  msg <- function(...) if (isTRUE(verbose)) message(...)

  # ---- s2 OFF ----
  old_s2 <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  `%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

  year_str <- if (is.null(year_target) || is.na(year_target)) "NA" else as.character(year_target)

  .safe_slug <- function(x, max_len = 40) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    if (nchar(x) > max_len) x <- substr(x, 1, max_len)
    if (!nzchar(x)) x <- "X"
    x
  }

  short_id <- function(x) {
    x <- enc2utf8(as.character(x))
    xi <- utf8ToInt(x)
    if (length(xi) == 0) return("000000")
    h <- sum(((seq_along(xi) * xi) %% 1000003), na.rm = TRUE) %% 16777215
    sprintf("%06x", as.integer(h))
  }

  # -----------------------------
  # Rounding helper (global)
  # -----------------------------
  round_outputs <- function(dt,
                            digits_metrics = 3,
                            digits_area = 2,
                            digits_res = 0) {
    if (is.null(dt) || nrow(dt) == 0) return(dt)
    dt <- data.table::as.data.table(dt)

    # resolution
    for (cc in intersect(c("TemplateRes_m"), names(dt))) {
      if (is.numeric(dt[[cc]])) dt[[cc]] <- round(dt[[cc]], digits_res)
    }

    # areas
    area_cols <- grep("(_ha$|Area_)", names(dt), value = TRUE)
    for (cc in area_cols) {
      if (is.numeric(dt[[cc]])) dt[[cc]] <- round(dt[[cc]], digits_area)
    }

    # metrics
    metric_cols <- intersect(
      c("Precision","Recall","F1","IoU","ErrorRate","Commission","Omission",
        "Recall_Area_percent","Precision_Area_percent","Perc_Detected_Polygons"),
      names(dt)
    )
    for (cc in metric_cols) {
      if (is.numeric(dt[[cc]])) dt[[cc]] <- round(dt[[cc]], digits_metrics)
    }

    dt
  }

  # -----------------------------
  # Helpers: sf validity + CRS
  # -----------------------------
  make_valid_sf <- function(x) {
    x <- sf::st_as_sf(x)
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")

    x <- tryCatch(
      sf::st_make_valid(x),
      error = function(e) suppressWarnings(sf::st_buffer(x, 0))
    )

    x <- x[!sf::st_is_empty(x), , drop = FALSE]

    # Only extract polygons when needed (avoids: "x is already of type POLYGON")
    gt <- sf::st_geometry_type(x, by_geometry = TRUE)
    if (any(!gt %in% c("POLYGON", "MULTIPOLYGON"))) {
      x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON"))
      x <- x[!sf::st_is_empty(x), , drop = FALSE]
    }

    if (nrow(x) > 0) {
      sf::st_geometry(x) <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)
      x <- sf::st_as_sf(x)
    }
    x
  }

  crs_equal_safe <- function(x, crs_target) {
    cx <- sf::st_crs(x)
    if (is.na(cx)) return(NA)
    out <- tryCatch(cx == crs_target, error = function(e) NA)
    if (is.na(out)) NA else isTRUE(out)
  }

  transform_to_crs_safe <- function(x, crs_target, label = "layer") {
    x <- make_valid_sf(x)
    if (nrow(x) == 0) return(x)

    if (is.na(sf::st_crs(x))) stop(label, " has NA CRS. Define/assign CRS before validation.")

    same <- crs_equal_safe(x, crs_target)
    if (is.na(same) || !same) {
      x <- sf::st_transform(x, crs_target)
      x <- make_valid_sf(x)
    }
    x
  }

  dissolve_by_field <- function(x, field) {
    if (is.null(field) || !nzchar(field)) return(x)
    if (!field %in% names(x)) stop("Dissolve field not found: ", field)

    x <- make_valid_sf(x)
    if (nrow(x) == 0) return(x)

    geom_col <- attr(x, "sf_column")
    if (is.null(geom_col) || !nzchar(geom_col)) geom_col <- "geometry"

    # all non-geometry attributes except the dissolve field
    non_geom <- setdiff(names(x), geom_col)
    others   <- setdiff(non_geom, field)

    # aggregation that preserves columns (values become "representative" / concatenated)
    agg_any <- function(v) {
      if (inherits(v, "Date")) {
        vv <- v[!is.na(v)]
        return(if (length(vv)) vv[1] else as.Date(NA))
      }
      if (inherits(v, "POSIXt")) {
        vv <- v[!is.na(v)]
        return(if (length(vv)) vv[1] else as.POSIXct(NA))
      }
      if (is.numeric(v) || is.integer(v)) {
        vv <- v[!is.na(v)]
        return(if (length(vv)) vv[1] else NA_real_)
      }
      # character/factor/other -> concat unique
      vals <- unique(na.omit(as.character(v)))
      if (!length(vals)) return(NA_character_)
      paste(vals, collapse = ";")
    }

    # sf::summarise unions geometries; we keep all attributes via across()
    out <- x |>
      dplyr::group_by(.data[[field]]) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(others), agg_any),
        do_union = TRUE,
        .groups = "drop"
      ) |>
      make_valid_sf()

    out
  }

  safe_vect <- function(sfobj, label = "") {
    tryCatch(
      terra::vect(sfobj),
      error = function(e) {
        sfobj2 <- make_valid_sf(sfobj)
        ok <- sf::st_is_valid(sfobj2)
        if (any(!ok)) {
          warning(sprintf("(%s) Dropping %d invalid geometries.", label, sum(!ok)))
          sfobj2 <- sfobj2[ok, , drop = FALSE]
        }
        terra::vect(sfobj2)
      }
    )
  }

  # -----------------------------
  # Read sf from path or sf
  # -----------------------------
  read_sf_any <- function(x, layer = NULL) {
    if (inherits(x, "sf")) return(make_valid_sf(x))
    if (!is.character(x) || length(x) != 1) stop("Input must be a path (length 1) or an sf object.")
    if (!file.exists(x)) stop("File not found: ", x)

    ext <- tolower(tools::file_ext(x))
    if (ext == "gpkg") {
      lyr_names <- sf::st_layers(x)$name
      lyr <- layer

      if (!is.null(lyr) && nzchar(lyr) && !(lyr %in% lyr_names)) {
        stop("Requested layer '", lyr, "' not found in: ", basename(x),
             "\nAvailable: ", paste(lyr_names, collapse = ", "))
      }

      if (is.null(lyr) || !nzchar(lyr)) {
        pref <- c("polys", "scored", "ref", "internal_flagged")
        hit <- pref[pref %in% lyr_names][1]
        if (is.na(hit) || !nzchar(hit)) hit <- lyr_names[1]
        lyr <- hit
      }

      sf::st_read(x, layer = lyr, quiet = TRUE) |> make_valid_sf()
    } else {
      sf::st_read(x, quiet = TRUE) |> make_valid_sf()
    }
  }

  # -----------------------------
  # GPKG writer (ArcGIS-friendly)
  # -----------------------------
  sanitize_for_gpkg <- function(x,
                                force_geom_name = "geometry",
                                max_field_name_len = 50L,
                                max_fields = 1800L) {
    x <- make_valid_sf(x)

    # ensure single active geometry
    geom_col <- attr(x, "sf_column")
    if (is.null(geom_col) || !nzchar(geom_col) || !(geom_col %in% names(x)) || !inherits(x[[geom_col]], "sfc")) {
      sfc_cols <- names(x)[vapply(x, inherits, logical(1), what = "sfc")]
      if (length(sfc_cols) == 0) stop("No geometry column (sfc) found in sf.")
      geom_col <- sfc_cols[1]
      sf::st_geometry(x) <- geom_col
    }

    # drop extra sfc columns
    sfc_cols <- names(x)[vapply(x, inherits, logical(1), what = "sfc")]
    drop_sfc <- setdiff(sfc_cols, geom_col)
    if (length(drop_sfc) > 0) x[drop_sfc] <- NULL

    # rename geometry to stable name
    if (!identical(geom_col, force_geom_name)) {
      names(x)[names(x) == geom_col] <- force_geom_name
      sf::st_geometry(x) <- force_geom_name
      geom_col <- force_geom_name
    }

    # avoid FID collisions
    nm <- names(x)
    bad_fid <- (tolower(nm) %in% c("fid", "ogr_fid", "ogc_fid")) & (nm != geom_col)
    if (any(bad_fid)) {
      nm[bad_fid] <- paste0(nm[bad_fid], "_src")
      names(x) <- nm
    }

    # flatten list columns (except sfc)
    is_listcol <- vapply(x, is.list, logical(1))
    is_sfc     <- vapply(x, inherits, logical(1), what = "sfc")
    to_fix <- which(is_listcol & !is_sfc)
    if (length(to_fix) > 0) {
      for (i in to_fix) {
        cc <- names(x)[i]
        x[[cc]] <- vapply(x[[cc]], function(z) paste(z, collapse = ","), character(1))
      }
    }

    # convert units columns to numeric
    for (cc in names(x)) {
      if (cc == geom_col) next
      if (inherits(x[[cc]], "units")) x[[cc]] <- as.numeric(x[[cc]])
    }

    make_unique_limited <- function(xn, max_len) {
      out <- character(length(xn))
      seen <- new.env(parent = emptyenv())
      for (i in seq_along(xn)) {
        base <- xn[i]
        if (!nzchar(base)) base <- "x"
        if (nchar(base) > max_len) base <- substr(base, 1, max_len)

        key0 <- base
        if (exists(key0, envir = seen, inherits = FALSE)) {
          k <- get(key0, envir = seen, inherits = FALSE) + 1L
        } else {
          k <- 0L
        }

        cand <- base
        while (cand %in% out[seq_len(i - 1L)]) {
          k <- k + 1L
          suf <- paste0("_", k)
          keep <- max_len - nchar(suf)
          cand <- paste0(substr(base, 1, max(1L, keep)), suf)
        }
        assign(key0, k, envir = seen)
        out[i] <- cand
      }
      out
    }

    nm <- names(x)
    geom_idx <- nm == geom_col

    nm2 <- nm
    nm2[!geom_idx] <- tolower(nm2[!geom_idx])
    nm2[!geom_idx] <- gsub("[^a-z0-9_]+", "_", nm2[!geom_idx])
    nm2[!geom_idx] <- gsub("^_+|_+$", "", nm2[!geom_idx])
    nm2[!geom_idx] <- ifelse(nzchar(nm2[!geom_idx]), nm2[!geom_idx], "x")
    nm2[!geom_idx] <- ifelse(grepl("^[0-9]", nm2[!geom_idx]), paste0("f_", nm2[!geom_idx]), nm2[!geom_idx])

    nm2[!geom_idx] <- substr(nm2[!geom_idx], 1, max_field_name_len)
    nm2[!geom_idx] <- make_unique_limited(nm2[!geom_idx], max_len = max_field_name_len)

    nm2[geom_idx] <- geom_col
    names(x) <- nm2
    sf::st_geometry(x) <- geom_col

    non_geom <- setdiff(names(x), geom_col)
    if (length(non_geom) > max_fields) {
      all_na <- vapply(x[non_geom], function(v) all(is.na(v)), logical(1))
      drop1 <- names(all_na)[all_na]
      keep <- non_geom
      if (length(keep) > max_fields && length(drop1) > 0) keep <- setdiff(keep, drop1)
      if (length(keep) > max_fields) keep <- keep[seq_len(max_fields)]
      x <- x[, c(keep, geom_col), drop = FALSE]
      warning("sanitize_for_gpkg: too many fields; reduced to ", max_fields, " attributes + geometry.")
    }

    x
  }

  write_gpkg_safe <- function(x,
                              path,
                              layer = "polys",
                              fallback_epsg = 3035,
                              max_fields = 10000L) {

    # sanitize (keeps geometry stable; cleans names; flattens list cols)
    x <- sanitize_for_gpkg(x, max_fields = as.integer(max_fields))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

    # ---- FORCE a standard EPSG CRS for ArcGIS compatibility (no coordinate transform) ----
    cx <- sf::st_crs(x)
    epsg_to_write <- NA_integer_

    if (!is.na(cx) && !is.na(cx$epsg)) {
      epsg_to_write <- as.integer(cx$epsg)
    } else if (exists("mask_crs", inherits = TRUE)) {
      mcrs <- get("mask_crs", inherits = TRUE)
      if (!is.na(mcrs) && !is.na(mcrs$epsg)) epsg_to_write <- as.integer(mcrs$epsg)
    }

    if (is.na(epsg_to_write)) epsg_to_write <- as.integer(fallback_epsg)
    sf::st_crs(x) <- sf::st_crs(epsg_to_write)

    # --- if target exists, delete and VERIFY (Windows locks are common)
    if (file.exists(path)) {
      try(unlink(path, force = TRUE), silent = TRUE)
      if (file.exists(path)) {
        stop("Could not overwrite GPKG (file locked?): ", path,
             "\nClose ArcGIS/QGIS and disable Explorer preview, then re-run.")
      }
    }

    # --- write to temp first (avoids partial/locked file issues)
    tmp <- file.path(
      dirname(path),
      paste0(tools::file_path_sans_ext(basename(path)),
             "_tmp_", Sys.getpid(), "_", format(Sys.time(), "%H%M%S"), ".gpkg")
    )
    if (file.exists(tmp)) unlink(tmp, force = TRUE)

    old <- Sys.getenv("OGR_WKT_FORMAT", unset = NA_character_)
    if (isTRUE(arcgis_wkt1)) Sys.setenv(OGR_WKT_FORMAT = "WKT1_ESRI")
    on.exit({
      if (isTRUE(arcgis_wkt1)) {
        if (!is.na(old)) Sys.setenv(OGR_WKT_FORMAT = old) else Sys.unsetenv("OGR_WKT_FORMAT")
      }
    }, add = TRUE)

    ok <- try(
      sf::st_write(x, tmp, driver = "GPKG", layer = layer, delete_dsn = TRUE, quiet = TRUE),
      silent = TRUE
    )
    if (inherits(ok, "try-error")) {
      if (file.exists(tmp)) unlink(tmp, force = TRUE)
      stop("Could not write GPKG temp: ", basename(tmp), "\n", ok)
    }

    ok_mv <- file.rename(tmp, path)
    if (!isTRUE(ok_mv)) {
      stop("Wrote temp GPKG but could not rename to final path.\nTemp: ", tmp, "\nFinal: ", path,
           "\n(Usually a lock/permission issue on the final path.)")
    }

    invisible(TRUE)
  }

  # ---- KEEP ALL ATTRIBUTES + SAFE BIND FOR sf ----
  keep_all_attrs_sf <- function(x,
                                id_col = "._cid",
                                max_fields = 10000L) {
    x <- make_valid_sf(x)
    if (nrow(x) == 0) return(x)

    if (!(id_col %in% names(x))) x[[id_col]] <- seq_len(nrow(x))

    # sanitize for gpkg safety but keep ALL columns (up to max_fields)
    x <- sanitize_for_gpkg(x, max_fields = as.integer(max_fields))
    x
  }

  safe_bind_rows_sf <- function(a,
                                b,
                                max_fields = 10000L) {

    a <- keep_all_attrs_sf(a, id_col = "._cid", max_fields = max_fields)
    b <- keep_all_attrs_sf(b, id_col = "._cid", max_fields = max_fields)

    # sanitize_for_gpkg forces geom name to "geometry"
    gcol <- "geometry"
    if (!(gcol %in% names(a)) || !inherits(a[[gcol]], "sfc")) {
      gcol <- attr(a, "sf_column")
      if (is.null(gcol) || !nzchar(gcol)) gcol <- "geometry"
    }
    if (!(gcol %in% names(b)) || !inherits(b[[gcol]], "sfc")) {
      # last resort: set active geometry from current sf_column
      gb <- attr(b, "sf_column")
      if (!is.null(gb) && nzchar(gb) && gb %in% names(b)) sf::st_geometry(b) <- gb
    }

    # add missing columns
    all_cols <- union(names(a), names(b))

    add_missing <- function(x) {
      miss <- setdiff(all_cols, names(x))
      if (length(miss)) {
        for (m in miss) x[[m]] <- NA
      }
      x <- x[, all_cols, drop = FALSE]
      x
    }
    a <- add_missing(a)
    b <- add_missing(b)

    # 1st attempt
    out_try <- try(dplyr::bind_rows(a, b), silent = TRUE)
    if (!inherits(out_try, "try-error")) return(make_valid_sf(out_try))

    # fallback: coerce problematic mixed-type columns to character
    non_geom <- setdiff(all_cols, gcol)

    is_texty <- function(z) is.character(z) || is.factor(z)
    is_timey <- function(z) inherits(z, "Date") || inherits(z, "POSIXt")

    for (cc in non_geom) {
      if (!cc %in% names(a) || !cc %in% names(b)) next
      if (inherits(a[[cc]], "sfc") || inherits(b[[cc]], "sfc")) next

      ca <- a[[cc]]
      cb <- b[[cc]]

      # If either side is character/factor OR Date/POSIXt, safest is character
      if (is_texty(ca) || is_texty(cb) || is_timey(ca) || is_timey(cb)) {
        a[[cc]] <- as.character(ca)
        b[[cc]] <- as.character(cb)
      }
    }

    out <- dplyr::bind_rows(a, b)
    make_valid_sf(out)
  }

  # -----------------------------
  # OUTPUT NAMES (short)
  # - UPDATED: include output_tag in CSV/XLSX filenames
  # -----------------------------
  prefix_clean <- .safe_slug(output_prefix, 20)
  tag_clean <- if (!is.null(output_tag) && nzchar(output_tag)) .safe_slug(output_tag, 60) else NULL

  idx_str <- function(i) sprintf(paste0("%0", as.integer(index_width), "d"), as.integer(i))

  base_name <- function(suffix, idx = NULL, always_index = FALSE) {
    b <- paste0(prefix_clean, "_", year_str, "_", suffix)
    if (!is.null(idx) && (always_index || idx > 1L)) b <- paste0(b, "_", idx_str(idx))
    b
  }

  # tag_in: only add output_tag to these extensions (keeps GPKG/TIF short)
  out_file <- function(dir, suffix, ext, idx = NULL, always_index = FALSE,
                       tag_in = c("csv", "xlsx")) {
    tag_token <- if (!is.null(tag_clean) && ext %in% tag_in) paste0("_", tag_clean) else ""
    file.path(dir, paste0(base_name(suffix, idx = idx, always_index = always_index),
                          tag_token, ".", ext))
  }

  tag_in_conf_err <- c("csv","xlsx","gpkg","tif")

  validation_output_dir <- file.path(validation_dir, outputs_subdir)
  if (!dir.exists(validation_output_dir)) dir.create(validation_output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(tag_clean) && isTRUE(tag_as_subdir)) {
    validation_output_dir <- file.path(validation_output_dir, paste0(prefix_clean, "_", tag_clean))
    if (!dir.exists(validation_output_dir)) dir.create(validation_output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  errors_dir <- file.path(validation_output_dir, errors_subdir)
  if (isTRUE(write_error_vectors) && !dir.exists(errors_dir)) dir.create(errors_dir, recursive = TRUE, showWarnings = FALSE)
  # -----------------------------
  # CONF outputs dir (optional)
  # -----------------------------
  need_conf <- isTRUE(write_confusion_vectors) || isTRUE(write_confusion_pixel_raster)

  conf_sub <- confusion_subdir
  if (is.null(conf_sub) || length(conf_sub) != 1L || is.na(conf_sub) || !nzchar(conf_sub)) {
    conf_sub <- "CONF"
  }

  conf_dir <- file.path(validation_output_dir, conf_sub)

  if (need_conf && !dir.exists(conf_dir)) {
    dir.create(conf_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # -----------------------------
  # Load mask  (FORCE EPSG if CRS is "custom")
  # -----------------------------
  if (!file.exists(mask_shapefile)) stop("mask_shapefile not found: ", mask_shapefile)

  mask_geom <- sf::st_read(mask_shapefile, quiet = TRUE) |> make_valid_sf()
  mask_crs  <- sf::st_crs(mask_geom)

  if (is.na(mask_crs)) stop("mask_shapefile has NA CRS. Fix CRS before running validation.")

  # If the mask CRS has no EPSG (common in some ESRI WKT exports), ArcGIS may reject GPKG outputs.
  # In your workflow everything is in EPSG:3035, so we safely assign it (no coordinate transform).
  if (is.na(mask_crs$epsg)) {
    sf::st_crs(mask_geom) <- sf::st_crs(3035)
    mask_crs <- sf::st_crs(3035)
  }

  mask_v <- terra::vect(mask_geom)

  ext_m <- terra::ext(mask_v)
  template_domain <- terra::rast(ext = ext_m, res = template_res, crs = mask_crs$wkt)
  domain_mask <- terra::rasterize(mask_v, template_domain, field = 1, background = NA)
  cell_area_ha <- abs(prod(terra::res(domain_mask))) / 10000

  # -----------------------------
  # Load / align class raster (optional)
  # -----------------------------
  class_r <- NULL
  class_lut_df <- NULL

  if (!is.null(class_raster)) {
    if (is.character(class_raster)) {
      if (length(class_raster) != 1) stop("class_raster must be a single path or a SpatRaster.")
      if (!file.exists(class_raster)) stop("class_raster not found: ", class_raster)
    }
    class_r <- if (inherits(class_raster, "SpatRaster")) class_raster else terra::rast(class_raster)
    class_r <- class_r[[as.integer(class_band)]]

    same_geom <- terra::compareGeom(class_r, template_domain, stopOnError = FALSE)
    if (!isTRUE(same_geom)) {
      if (!terra::same.crs(class_r, template_domain)) {
        class_r <- terra::project(class_r, template_domain, method = "near")
      }
      class_r <- terra::resample(class_r, template_domain, method = "near")
    }

    class_r <- terra::mask(class_r, domain_mask)

    if (!is.null(class_values) && length(class_values) > 0) {
      vv <- as.numeric(class_values)
      class_r <- terra::ifel(class_r %in% vv, class_r, NA)
    }

    if (!is.null(class_lut)) {
      if (is.character(class_lut) && length(class_lut) == 1) {
        if (!file.exists(class_lut)) {
          stop("class_lut path not found: ", class_lut,
               "\nCheck the filename/extension (e.g., it should usually end with '.csv').")
        }
        class_lut_df <- data.frame(data.table::fread(class_lut))
      } else {
        class_lut_df <- as.data.frame(class_lut)
      }

      # robust column detection
      nms0 <- names(class_lut_df)
      nms  <- tolower(nms0)

      id_idx <- which(nms %in% c("id", "classvalue", "class_value", "value"))[1]
      lab_idx <- which(nms %in% c("label", "classlabel", "class_label", "name"))[1]

      if (is.na(id_idx) || is.na(lab_idx)) {
        # fallback: use first two columns
        if (ncol(class_lut_df) >= 2) {
          id_idx <- 1
          lab_idx <- 2
          warning("class_lut does not have clear id/label column names; using the first two columns as id,label.")
        } else {
          warning("class_lut must have at least 2 columns (id,label). Ignoring LUT.")
          class_lut_df <- NULL
        }
      }

      if (!is.null(class_lut_df)) {
        class_lut_df <- class_lut_df[, c(id_idx, lab_idx), drop = FALSE]
        names(class_lut_df) <- c("id", "label")
        class_lut_df$id <- suppressWarnings(as.integer(class_lut_df$id))
        class_lut_df$label <- as.character(class_lut_df$label)
        class_lut_df <- class_lut_df[!is.na(class_lut_df$id), , drop = FALSE]
      }
    }
  }

  # -----------------------------
  # Reference cache
  # -----------------------------
  ref_identity <- paste0(
    "REFPATH=", if (is.character(ref_shapefile)) normalizePath(ref_shapefile, winslash = "\\", mustWork = FALSE) else "sf_object",
    "|REFLAYER=", ifelse(is.null(ref_layer), "auto", ref_layer),
    "|MASK=", normalizePath(mask_shapefile, winslash = "\\", mustWork = FALSE),
    "|YEAR=", year_str,
    "|RES=", template_res,
    "|MINDET=", min_detected_percent,
    "|COMP=", threshold_completely_detected,
    "|MINA=", ifelse(is.null(min_area_reference_ha), "none", format(min_area_reference_ha, scientific = FALSE)),
    "|DREF=", ifelse(is.null(dissolve_ref_by), "none", dissolve_ref_by)
  )
  ref_cache_id <- paste0("ref_", short_id(ref_identity))

  ref_vec_cache  <- file.path(validation_output_dir, paste0(ref_cache_id, ".rds"))
  ref_rast_cache <- file.path(validation_output_dir, paste0(ref_cache_id, ".tif"))

  if (isTRUE(force_reprocess_ref)) {
    msg("Force reprocessing reference polygons: deleting cached reference...")
    if (file.exists(ref_vec_cache)) unlink(ref_vec_cache)
    if (file.exists(ref_rast_cache)) unlink(ref_rast_cache)
  }

  # -----------------------------
  # Build/load reference
  # -----------------------------
  if (file.exists(ref_vec_cache) && file.exists(ref_rast_cache)) {
    msg("Loading cached reference (RDS + raster mask)...")
    ref_polygons <- readRDS(ref_vec_cache) |> make_valid_sf()
    ref_mask_r   <- terra::rast(ref_rast_cache)
    ref_polygons <- transform_to_crs_safe(ref_polygons, mask_crs, label = "ref_polygons (cached)")
  } else {
    msg("Processing reference polygons...")
    ref_polygons <- read_sf_any(ref_shapefile, layer = ref_layer)

    # year filter (robust)
    year_col <- intersect(names(ref_polygons), c("year", "Year", "YEAR"))[1]
    if (!is.na(year_col) && nzchar(year_col)) {
      ref_polygons <- dplyr::filter(ref_polygons, .data[[year_col]] == year_target)
      msg("Reference polygons filtered by year (", year_target, "): ", nrow(ref_polygons))
    }
    if (nrow(ref_polygons) == 0) stop("No reference polygons for year_target.")

    ref_polygons <- transform_to_crs_safe(ref_polygons, mask_crs, label = "ref_polygons")

    ref_polygons <- sf::st_filter(ref_polygons, mask_geom, .predicate = sf::st_intersects)
    msg("Reference polygons after mask filter (st_filter): ", nrow(ref_polygons))
    if (nrow(ref_polygons) == 0) stop("No reference polygons inside mask.")

    ref_polygons <- suppressWarnings(sf::st_intersection(ref_polygons, mask_geom)) |> make_valid_sf()
    ref_polygons <- dissolve_by_field(ref_polygons, dissolve_ref_by)

    if (!is.null(min_area_reference_ha)) {
      area_ha <- as.numeric(sf::st_area(ref_polygons)) / 10000
      keep <- area_ha >= min_area_reference_ha
      msg(sprintf("Filtered small reference polygons: %d -> %d (Area >= %.2f ha)",
                  length(keep), sum(keep), min_area_reference_ha))
      ref_polygons <- ref_polygons[keep, , drop = FALSE] |> make_valid_sf()
      if (nrow(ref_polygons) == 0) stop("No reference polygons remain after min_area_reference_ha.")
    }

    ref_v <- safe_vect(ref_polygons, label = "ref")
    ref_mask_r <- terra::rasterize(ref_v, domain_mask, field = 1, background = 0, touches = TRUE)
    ref_mask_r[is.na(domain_mask)] <- NA

    saveRDS(ref_polygons, ref_vec_cache)
    terra::writeRaster(ref_mask_r, ref_rast_cache, overwrite = TRUE,
                       datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255"))
  }

  ref_area_vec_ha <- as.numeric(sf::st_area(ref_polygons)) / 10000

  # -----------------------------
  # Normalize input list
  # -----------------------------
  if (inherits(input_shapefile, "sf")) {
    input_list  <- list(input_shapefile)
    input_names <- "input_sf"
  } else {
    if (!is.character(input_shapefile)) stop("input_shapefile must be path(s) or an sf object.")
    input_list  <- as.list(input_shapefile)
    input_names <- tools::file_path_sans_ext(basename(unlist(input_list)))
  }

  n_inputs <- length(input_list)
  always_index <- (n_inputs > 1L)

  metrics_list <- list()
  polygon_summary_list <- list()

  confusion_polys_paths  <- vector("list", n_inputs)
  confusion_pixels_paths <- vector("list", n_inputs)

  class_metrics_list <- list()
  class_metrics_global_list <- list()
  omission_by_class_list <- list()
  commission_by_class_list <- list()

  terra_freq_noNA <- function(r) {
    fr <- terra::freq(r)
    if (is.null(fr)) return(NULL)
    fr <- as.data.frame(fr)
    if (!("value" %in% names(fr)) && ncol(fr) >= 1) names(fr)[1] <- "value"
    if (!("count" %in% names(fr)) && ncol(fr) >= 2) names(fr)[2] <- "count"
    if ("value" %in% names(fr)) fr <- fr[!is.na(fr$value), , drop = FALSE]
    fr
  }

  # -----------------------------
  # Helper: classwise pixel metrics (chunked)
  # -----------------------------
  classwise_pixel_metrics <- function(pred_r, ref_r, class_r, burned_value = 1) {

    if (!isTRUE(terra::compareGeom(pred_r, ref_r, stopOnError = FALSE))) stop("pred and ref are not aligned.")
    if (!isTRUE(terra::compareGeom(pred_r, class_r, stopOnError = FALSE))) stop("class_raster is not aligned to the validation grid.")

    fr <- terra_freq_noNA(class_r)
    if (is.null(fr) || nrow(fr) == 0) return(NULL)

    ids <- sort(as.integer(fr$value))
    ids <- ids[is.finite(ids)]
    if (!is.null(class_values) && length(class_values) > 0) {
      ids <- ids[ids %in% as.integer(class_values)]
    }
    if (length(ids) == 0) return(NULL)

    K <- length(ids)
    counts_vec <- numeric(K * 4)
    g_TN <- 0; g_FN <- 0; g_FP <- 0; g_TP <- 0

    nr <- terra::nrow(pred_r)
    starts <- seq(1, nr, by = as.integer(class_chunk_rows))

    terra::readStart(pred_r); terra::readStart(ref_r); terra::readStart(class_r)
    on.exit({ terra::readStop(pred_r); terra::readStop(ref_r); terra::readStop(class_r) }, add = TRUE)

    for (r0 in starts) {
      nrs <- min(as.integer(class_chunk_rows), nr - r0 + 1)

      p   <- terra::readValues(pred_r, row = r0, nrows = nrs, mat = FALSE)
      r   <- terra::readValues(ref_r,  row = r0, nrows = nrs, mat = FALSE)
      sid <- terra::readValues(class_r, row = r0, nrows = nrs, mat = FALSE)

      ok <- if (class_domain == "class") {
        !is.na(sid)
      } else {
        !is.na(p) & !is.na(r) & !is.na(sid)
      }
      if (!any(ok)) next

      sid <- sid[ok]
      p   <- p[ok]
      r   <- r[ok]

      if (isTRUE(class_na_as_zero)) {
        p[is.na(p)] <- 0
        r[is.na(r)] <- 0
      } else {
        keep2 <- !is.na(p) & !is.na(r)
        if (!any(keep2)) next
        sid <- sid[keep2]; p <- p[keep2]; r <- r[keep2]
      }

      sid_int <- as.integer(round(sid))
      p01 <- as.integer(p == burned_value)
      r01 <- as.integer(r == burned_value)

      state <- p01 * 2L + r01

      tstate <- tabulate(state + 1L, nbins = 4L)
      g_TN <- g_TN + tstate[1]
      g_FN <- g_FN + tstate[2]
      g_FP <- g_FP + tstate[3]
      g_TP <- g_TP + tstate[4]

      idx <- match(sid_int, ids)
      keep <- !is.na(idx)
      if (!any(keep)) next

      key <- (idx[keep] - 1L) * 4L + state[keep] + 1L
      counts_vec <- counts_vec + tabulate(key, nbins = K * 4L)
    }

    mat <- matrix(counts_vec, nrow = K, ncol = 4, byrow = TRUE)
    TN <- mat[, 1]; FN <- mat[, 2]; FP <- mat[, 3]; TP <- mat[, 4]

    Precision  <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
    Recall     <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
    F1         <- ifelse(!is.na(Precision) & !is.na(Recall) & (Precision + Recall) > 0,
                         2 * Precision * Recall / (Precision + Recall), NA_real_)
    IoU        <- ifelse((TP + FP + FN) > 0, TP / (TP + FP + FN), NA_real_)
    Commission <- ifelse(!is.na(Precision), 1 - Precision, NA_real_)
    Omission   <- ifelse(!is.na(Recall), 1 - Recall, NA_real_)

    Area_Reference_ha    <- (TP + FN) * cell_area_ha
    Area_Detected_ha     <- (TP + FP) * cell_area_ha
    Area_Intersection_ha <- TP * cell_area_ha

    Recall_Area_percent    <- ifelse(Area_Reference_ha > 0, 100 * Area_Intersection_ha / Area_Reference_ha, NA_real_)
    Precision_Area_percent <- ifelse(Area_Detected_ha  > 0, 100 * Area_Intersection_ha / Area_Detected_ha,  NA_real_)

    out <- data.table::data.table(
      ClassValue = ids,
      TP = as.integer(TP), FP = as.integer(FP), FN = as.integer(FN), TN = as.integer(TN),
      Precision = Precision, Recall = Recall, F1 = F1, IoU = IoU,
      Commission = Commission, Omission = Omission,
      Area_Reference_ha = Area_Reference_ha,
      Area_Detected_ha = Area_Detected_ha,
      Area_Intersection_ha = Area_Intersection_ha,
      Recall_Area_percent = Recall_Area_percent,
      Precision_Area_percent = Precision_Area_percent
    )
    data.table::setDT(out)

    if (isTRUE(class_min_ref_area_ha > 0)) {
      out <- out[Area_Reference_ha >= class_min_ref_area_ha]
      data.table::setDT(out)
    }

    # attach labels (and optionally include all LUT ids)
    if (!is.null(class_lut_df)) {

      if (isTRUE(class_include_all_lut_ids)) {

        all_ids <- as.integer(class_lut_df$id)
        all_lab <- as.character(class_lut_df$label)

        if (!is.null(class_values) && length(class_values) > 0) {
          keepi <- all_ids %in% as.integer(class_values)
          all_ids <- all_ids[keepi]
          all_lab <- all_lab[keepi]
        }

        base <- data.table::data.table(ClassValue = all_ids, ClassLabel = all_lab)

        # IMPORTANT: use merge to keep data.table class (avoid shallow-copy warnings later)
        out <- merge(base, out, by = "ClassValue", all.x = TRUE, sort = FALSE)
        data.table::setDT(out)

        # missing classes -> counts/areas = 0, metrics remain NA
        for (cc in c("TP","FP","FN","TN")) {
          ii <- which(is.na(out[[cc]]))
          if (length(ii)) data.table::set(out, i = ii, j = cc, value = 0L)
        }
        for (cc in c("Area_Reference_ha","Area_Detected_ha","Area_Intersection_ha")) {
          ii <- which(is.na(out[[cc]]))
          if (length(ii)) data.table::set(out, i = ii, j = cc, value = 0)
        }

      } else {

        out <- merge(out, data.table::as.data.table(class_lut_df),
                                 by.x = "ClassValue", by.y = "id", all.x = TRUE, sort = FALSE)
        data.table::setDT(out)
        data.table::setnames(out, "label", "ClassLabel", skip_absent = TRUE)
      }
    }

    # ---- global metrics: DO NOT use attr<- on a data.table (causes shallow-copy warning later) ----
    g_precision <- if ((g_TP + g_FP) > 0) g_TP / (g_TP + g_FP) else NA_real_
    g_recall    <- if ((g_TP + g_FN) > 0) g_TP / (g_TP + g_FN) else NA_real_
    g_f1 <- if (!is.na(g_precision) && !is.na(g_recall) && (g_precision + g_recall) > 0) {
      2 * g_precision * g_recall / (g_precision + g_recall)
    } else {
      NA_real_
    }
    g_iou <- if ((g_TP + g_FP + g_FN) > 0) g_TP / (g_TP + g_FP + g_FN) else NA_real_

    global_dt <- data.table::data.table(
      TP = as.integer(g_TP), FP = as.integer(g_FP), FN = as.integer(g_FN), TN = as.integer(g_TN),
      Precision = g_precision,
      Recall    = g_recall,
      F1        = g_f1,
      IoU       = g_iou,
      Area_Reference_ha    = (g_TP + g_FN) * cell_area_ha,
      Area_Detected_ha     = (g_TP + g_FP) * cell_area_ha,
      Area_Intersection_ha = g_TP * cell_area_ha
    )

    data.table::setattr(out, "global", global_dt)

    out
  }
  # -----------------------------
  # Loop inputs
  # -----------------------------
  for (k in seq_along(input_list)) {

    shp <- input_list[[k]]
    input_name_raw <- input_names[k]
    msg("Processing input: ", input_name_raw)

    pred_tif <- out_file(validation_output_dir, "pred", "tif", idx = k, always_index = always_index)
    if (isTRUE(force_reprocess_pred) && file.exists(pred_tif)) unlink(pred_tif)

    det <- if (inherits(shp, "sf")) shp else read_sf_any(shp, layer = input_layer)
    det <- make_valid_sf(det)
    det <- transform_to_crs_safe(det, mask_crs, label = paste0("detections (", input_name_raw, ")"))

    det <- sf::st_filter(det, mask_geom, .predicate = sf::st_intersects)
    if (nrow(det) == 0) {
      warning("Skipping ", input_name_raw, ": no polygons inside mask.")
      next
    }
    det <- suppressWarnings(sf::st_intersection(det, mask_geom)) |> make_valid_sf()
    det <- dissolve_by_field(det, dissolve_input_by)

    if (nrow(det) == 0) {
      warning("Skipping ", input_name_raw, ": empty after dissolve/clip.")
      next
    }

    if (!file.exists(pred_tif)) {
      pred0 <- domain_mask
      pred0[] <- 0
      det_v <- safe_vect(det, label = input_name_raw)

      pred_r <- terra::rasterize(det_v, pred0, field = 1, background = 0, touches = TRUE)
      pred_r[is.na(domain_mask)] <- NA

      terra::writeRaster(pred_r, pred_tif, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255"))
    }

    pred_r <- terra::rast(pred_tif)

    if (isTRUE(strict_binary_check)) {
      vv <- unique(terra::values(pred_r))
      vv <- vv[!is.na(vv)]
      if (length(setdiff(vv, c(0, 1))) > 0) {
        stop("pred raster is not binary (0/1). Set strict_binary_check=FALSE to bypass.")
      }
    }

    # ---- OPTIONAL: PIXEL CONFUSION RASTER (TN/FN/FP/TP) ----
    if (isTRUE(write_confusion_pixel_raster)) {

      # Codes: 1=TN, 2=FN, 3=FP, 4=TP
      # state = (pred==1)*2 + (ref==1)  -> 0..3; then +1 -> 1..4
      conf_r <- terra::ifel(
        is.na(domain_mask),
        NA,
        ((pred_r == burned_value) * 2L + (ref_mask_r == burned_value)) + 1L
      )
      names(conf_r) <- "confusion"

      out_conf_tif <- out_file(conf_dir, confusion_pixels_suffix, "tif",
                               idx = k, always_index = always_index,
                               tag_in = tag_in_conf_err)

      terra::writeRaster(
        conf_r, out_conf_tif, overwrite = TRUE,
        datatype = "INT1U",
        gdal = c("COMPRESS=LZW", "NAflag=255")
      )

      # write LUT once per run (handy for later)
      lut_csv <- file.path(conf_dir, paste0("confusion_LUT", if (!is.null(tag_clean)) paste0("_", tag_clean) else "", ".csv"))
      if (!file.exists(lut_csv)) {
        utils::write.table(
          data.frame(Code = 1:4, Label = c("TN","FN","FP","TP")),
          lut_csv, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
        )
      }

      confusion_pixels_paths[[k]] <- out_conf_tif
    }

    # ---- PIXEL METRICS ----
    if (metrics_type %in% c("all", "pixel")) {

      tp <- as.numeric(terra::global((pred_r == burned_value) & (ref_mask_r == burned_value), "sum", na.rm = TRUE)[[1]][1])
      fp <- as.numeric(terra::global((pred_r == burned_value) & (ref_mask_r != burned_value), "sum", na.rm = TRUE)[[1]][1])
      fn <- as.numeric(terra::global((pred_r != burned_value) & (ref_mask_r == burned_value), "sum", na.rm = TRUE)[[1]][1])
      tn <- as.numeric(terra::global((pred_r != burned_value) & (ref_mask_r != burned_value), "sum", na.rm = TRUE)[[1]][1])

      precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
      recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
      f1        <- ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                          2 * precision * recall / (precision + recall), NA_real_)
      iou       <- ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_)
      error_rate <- ifelse((tp + fp + fn + tn) > 0, (fp + fn) / (tp + fp + fn + tn), NA_real_)

      metrics_list[[k]] <- data.table::data.table(
        Year = year_target,
        InputIndex = k,
        InputName = input_name_raw,
        TemplateRes_m = template_res,
        CellArea_ha = cell_area_ha,
        TP = tp, FP = fp, FN = fn, TN = tn,
        Precision = precision, Recall = recall, F1 = f1, IoU = iou, ErrorRate = error_rate,
        RunTag = run_tag %||% "",
        RefTag = ref_tag %||% ""
      )


      # ---- PIXEL METRICS BY CLASS (optional; chunked) ----
      if (!is.null(class_r) && isTRUE(write_classwise_metrics)) {

        cm <- classwise_pixel_metrics(
          pred_r = pred_r,
          ref_r  = ref_mask_r,
          class_r = class_r,
          burned_value = burned_value
        )

        if (!is.null(cm) && nrow(cm) > 0) {

          # IMPORTANT: grab global attribute BEFORE any coercion
          gcm <- attr(cm, "global")

          cm <- data.table::as.data.table(cm)

          # add metadata WITHOUT :=
          data.table::set(cm, j = "Year",          value = year_target)
          data.table::set(cm, j = "InputIndex",    value = k)
          data.table::set(cm, j = "InputName",     value = input_name_raw)
          data.table::set(cm, j = "TemplateRes_m", value = template_res)
          data.table::set(cm, j = "RunTag",        value = (run_tag %||% ""))
          data.table::set(cm, j = "RefTag",        value = (ref_tag %||% ""))

          front <- c("Year","InputIndex","InputName","TemplateRes_m","ClassValue")
          if ("ClassLabel" %in% names(cm)) front <- c(front, "ClassLabel")
          data.table::setcolorder(cm, c(front, setdiff(names(cm), front)))

          cm <- round_outputs(cm, digits_metrics, digits_area, digits_res)
          out_class_csv <- out_file(validation_output_dir, "metrics_by_class", "csv", idx = k, always_index = always_index)
          data.table::fwrite(cm, out_class_csv)

          if (!is.null(gcm) && nrow(gcm) > 0) {

            gcm <- data.table::as.data.table(gcm)

            data.table::set(gcm, j = "Year",          value = year_target)
            data.table::set(gcm, j = "InputIndex",    value = k)
            data.table::set(gcm, j = "InputName",     value = input_name_raw)
            data.table::set(gcm, j = "TemplateRes_m", value = template_res)
            data.table::set(gcm, j = "RunTag",        value = (run_tag %||% ""))
            data.table::set(gcm, j = "RefTag",        value = (ref_tag %||% ""))

            gcm <- round_outputs(gcm, digits_metrics, digits_area, digits_res)
            out_gcsv <- out_file(validation_output_dir, "metrics_by_class_global", "csv", idx = k, always_index = always_index)
            data.table::fwrite(gcm, out_gcsv)

            class_metrics_list[[k]] <- cm
            class_metrics_global_list[[k]] <- gcm
          } else {
            class_metrics_list[[k]] <- cm
            class_metrics_global_list[[k]] <- NULL
          }
        }
      }
    }

    # ---- AREA / POLYGON METRICS ----
    if (metrics_type %in% c("all", "area")) {

      if (isTRUE(union_detections_for_area)) {
        det_union_sfc <- sf::st_union(sf::st_geometry(det))
        det_union_sfc <- sf::st_sfc(det_union_sfc, crs = mask_crs)
        area_detected_total <- as.numeric(sf::st_area(det_union_sfc)) / 10000
      } else {
        det_union_sfc <- NULL
        area_detected_total <- sum(as.numeric(sf::st_area(det))) / 10000
      }

      ref2 <- ref_polygons
      ref2$._rid <- seq_len(nrow(ref2))

      if (isTRUE(union_detections_for_area)) {
        inter_all <- suppressWarnings(sf::st_intersection(ref2, det_union_sfc))
      } else {
        inter_all <- suppressWarnings(sf::st_intersection(ref2, det))
      }
      inter_all <- make_valid_sf(inter_all)

      int_area_i_ha <- numeric(nrow(ref2))
      int_area_i_ha[] <- 0

      if (nrow(inter_all) > 0) {
        a <- as.numeric(sf::st_area(inter_all)) / 10000
        rid <- inter_all$._rid
        a_by <- tapply(a, rid, sum)
        int_area_i_ha[as.integer(names(a_by))] <- as.numeric(a_by)
      }

      coverage_i <- ifelse(ref_area_vec_ha > 0, (int_area_i_ha / ref_area_vec_ha) * 100, 0)
      detected_i <- coverage_i >= min_detected_percent
      completely_detected <- sum(coverage_i >= threshold_completely_detected)

      # -----------------------------
      # OPTIONAL: CONFUSION POLYGONS (GPKG)
      # - ref polygons: TP (detected_i=TRUE) vs FN
      # - det polygons: TP (any overlap with ref) vs FP
      # (TN is not defined for polygons unless you model an explicit negative set)
      # -----------------------------
      if (isTRUE(write_confusion_vectors)) {

        # detection overlap boolean against ref
        det_overlap_any <- lengths(sf::st_intersects(det, ref_polygons)) > 0

        # ref side
        ref_conf <- ref_polygons
        ref_conf$Source <- "ref"
        ref_conf$Confusion <- ifelse(detected_i, "TP", "FN")
        ref_conf$Coverage_ref_pct <- as.numeric(coverage_i)
        ref_conf$Intersection_ha   <- as.numeric(int_area_i_ha)
        ref_conf$Reference_ha      <- as.numeric(ref_area_vec_ha)

        # det side
        det_conf <- det
        det_conf$Source <- "det"
        det_conf$Confusion <- ifelse(det_overlap_any, "TP", "FP")
        det_conf$Detected_ha <- as.numeric(sf::st_area(det_conf)) / 10000

        # optional overlap stats for detections
        if (isTRUE(compute_det_overlap_stats)) {

          det_conf$._did <- seq_len(nrow(det_conf))
          ref_union <- sf::st_union(sf::st_geometry(ref_polygons))
          ref_union <- sf::st_sfc(ref_union, crs = mask_crs)

          inter_det <- suppressWarnings(sf::st_intersection(det_conf, ref_union))
          inter_det <- make_valid_sf(inter_det)

          ov_ha <- numeric(nrow(det_conf))
          ov_ha[] <- 0

          if (nrow(inter_det) > 0) {
            a_det <- as.numeric(sf::st_area(inter_det)) / 10000
            did   <- inter_det$._did
            a_by_det <- tapply(a_det, did, sum)
            ov_ha[as.integer(names(a_by_det))] <- as.numeric(a_by_det)
          }

          det_conf$Overlap_ref_ha  <- ov_ha
          det_conf$Overlap_ref_pct <- ifelse(det_conf$Detected_ha > 0, 100 * ov_ha / det_conf$Detected_ha, NA_real_)

          det_conf$._did <- NULL
        }

        # add metadata (keep as-is; types will be harmonized safely at bind time)
        ref_conf$Year <- year_target
        ref_conf$InputIndex <- k
        ref_conf$InputName <- input_name_raw
        ref_conf$TemplateRes_m <- template_res

        det_conf$Year <- year_target
        det_conf$InputIndex <- k
        det_conf$InputName <- input_name_raw
        det_conf$TemplateRes_m <- template_res

        # ---- KEEP ALL INPUT VARIABLES + SAFE BIND (no type conflicts) ----
        conf_sf <- safe_bind_rows_sf(ref_conf, det_conf, max_fields = 10000L)

        out_conf_gpkg <- out_file(conf_dir, confusion_polys_suffix, "gpkg",
                                  idx = k, always_index = always_index,
                                  tag_in = tag_in_conf_err)
        write_gpkg_safe(conf_sf, out_conf_gpkg, layer = confusion_polys_layer, max_fields = 10000L)

        confusion_polys_paths[[k]] <- out_conf_gpkg
      }

      n_total <- length(ref_area_vec_ha)
      n_detected <- sum(detected_i)
      n_not_detected <- n_total - n_detected

      area_reference_total <- sum(ref_area_vec_ha)
      area_intersection_total <- sum(int_area_i_ha)

      recall_area <- ifelse(area_reference_total > 0, (area_intersection_total / area_reference_total) * 100, NA_real_)
      precision_area <- ifelse(area_detected_total > 0, (area_intersection_total / area_detected_total) * 100, NA_real_)

      polygon_summary_list[[k]] <- data.table::data.table(
        Year = year_target,
        InputIndex = k,
        InputName = input_name_raw,
        N_Reference_Polygons = n_total,
        N_Completely_Detected = completely_detected,
        N_Detected_Polygons = n_detected,
        N_Not_Detected = n_not_detected,
        Perc_Detected_Polygons = ifelse(n_total > 0, (n_detected / n_total) * 100, NA_real_),
        Area_Reference_ha = area_reference_total,
        Area_Detected_ha = area_detected_total,
        Area_Intersection_ha = area_intersection_total,
        Recall_Area_percent = recall_area,
        Precision_Area_percent = precision_area,
        Detected_Definition = paste0("coverage_ref>= ", min_detected_percent, "%"),
        Completely_Detected_Definition = paste0("coverage_ref>= ", threshold_completely_detected, "%"),
        UnionDetections = isTRUE(union_detections_for_area),
        TemplateRes_m = template_res,
        RunTag = run_tag %||% "",
        RefTag = ref_tag %||% ""
      )

      shrink_error_sf <- function(x, max_fields = 10000L) {
        x <- make_valid_sf(x)
        if (nrow(x) == 0) return(x)

        keep_candidates <- c(
          "year", "YEAR", "Year",
          "origin", "ORIGIN",
          "start_doy", "start_date", "doy",
          "id", "ID", "fireid", "poly_id",
          "class3_id", "merged_flag_value", "flag_internal", "flag_rbr"
        )
        keep <- intersect(names(x), keep_candidates)

        if (!("._eid" %in% names(x))) x$._eid <- seq_len(nrow(x))
        keep <- unique(c("._eid", keep))

        x2 <- x[, unique(c(keep, attr(x, "sf_column"))), drop = FALSE]
        x2 <- sanitize_for_gpkg(x2, max_fields = as.integer(max_fields))
        x2
      }

      if (isTRUE(write_error_vectors) || isTRUE(write_classwise_errors)) {

        ref_not_detected <- ref_polygons[!detected_i, , drop = FALSE]
        det_has_any_overlap <- lengths(sf::st_intersects(det, ref_polygons)) > 0
        det_not_matched <- det[!det_has_any_overlap, , drop = FALSE]

        # ---- helper: keep ALL attributes for error layers ----
        keep_all_error_sf <- function(x,
                                      id_col = "._eid",
                                      max_fields = 10000L) {
          x <- make_valid_sf(x)
          if (nrow(x) == 0) return(x)
          if (!(id_col %in% names(x))) x[[id_col]] <- seq_len(nrow(x))
          sanitize_for_gpkg(x, max_fields = as.integer(max_fields))
        }

        if (isTRUE(write_error_vectors)) {

          if (nrow(ref_not_detected) > 0) {
            out_ref_err <- out_file(errors_dir, "ref_not_detected", "gpkg",
                                    idx = k, always_index = always_index,
                                    tag_in = tag_in_conf_err)
            write_gpkg_safe(
              keep_all_error_sf(ref_not_detected, max_fields = 10000L),
              out_ref_err,
              layer = "ref_not_detected",
              max_fields = 10000L
            )
          }

          if (nrow(det_not_matched) > 0) {
            out_det_err <- out_file(errors_dir, "det_not_matched", "gpkg",
                                    idx = k, always_index = always_index,
                                    tag_in = tag_in_conf_err)
            write_gpkg_safe(
              keep_all_error_sf(det_not_matched, max_fields = 10000L),
              out_det_err,
              layer = "det_not_matched",
              max_fields = 10000L
            )
          }
        }

        if (!is.null(class_r) && isTRUE(write_classwise_errors)) {

          rasterize_err <- function(sfobj) {
            if (is.null(sfobj) || nrow(sfobj) == 0) return(NULL)
            vv <- safe_vect(sfobj, label = "err")
            rr <- terra::rasterize(vv, template_domain, field = 1, background = 0, touches = TRUE)
            rr[is.na(domain_mask)] <- NA
            rr
          }

          area_by_class <- function(err_r) {
            if (is.null(err_r)) return(NULL)
            cls2 <- class_r
            cls2[err_r == 0] <- NA
            fr2 <- terra_freq_noNA(cls2)
            if (is.null(fr2) || nrow(fr2) == 0) return(NULL)
            out <- data.table::data.table(
              ClassValue = as.integer(fr2$value),
              N_cells = as.integer(fr2$count),
              Area_ha = as.numeric(fr2$count) * cell_area_ha
            )
            if (!is.null(class_lut_df)) {
              out <- merge(out, class_lut_df, by.x = "ClassValue", by.y = "id", all.x = TRUE, sort = FALSE)
              data.table::setnames(out, "label", "ClassLabel", skip_absent = TRUE)
            }
            out
          }

          ref_nd_r <- rasterize_err(ref_not_detected)
          det_nm_r <- rasterize_err(det_not_matched)

          om <- area_by_class(ref_nd_r)
          if (!is.null(om) && nrow(om) > 0) {
            data.table::set(om,  j = "Year",       value = year_target)
            data.table::set(om,  j = "InputIndex", value = k)
            data.table::set(om,  j = "InputName",  value = input_name_raw)
            om <- round_outputs(om, digits_metrics, digits_area, digits_res)
            out_om <- out_file(validation_output_dir, "omission_by_class", "csv", idx = k, always_index = always_index)
            data.table::fwrite(om, out_om)
            omission_by_class_list[[k]] <- om
          }

          com <- area_by_class(det_nm_r)
          if (!is.null(com) && nrow(com) > 0) {
            data.table::set(com, j = "Year",       value = year_target)
            data.table::set(com, j = "InputIndex", value = k)
            data.table::set(com, j = "InputName",  value = input_name_raw)
            com <- round_outputs(com, digits_metrics, digits_area, digits_res)
            out_com <- out_file(validation_output_dir, "commission_by_class", "csv", idx = k, always_index = always_index)
            data.table::fwrite(com, out_com)
            commission_by_class_list[[k]] <- com
          }
        }
      }
    }
  }

  # -----------------------------
  # Write summaries
  # -----------------------------
  all_metrics <- NULL
  all_polygon_summary <- NULL

  if (metrics_type %in% c("all", "pixel")) {
    all_metrics <- data.table::rbindlist(metrics_list, fill = TRUE)
    all_metrics <- round_outputs(all_metrics, digits_metrics, digits_area, digits_res)
    out_csv <- out_file(validation_output_dir, "metrics", "csv", idx = NULL)
    data.table::fwrite(all_metrics, out_csv)
  }

  if (metrics_type %in% c("all", "area")) {
    all_polygon_summary <- data.table::rbindlist(polygon_summary_list, fill = TRUE)
    all_polygon_summary <- round_outputs(all_polygon_summary, digits_metrics, digits_area, digits_res)
    out_csv <- out_file(validation_output_dir, "polygon_summary", "csv", idx = NULL)
    data.table::fwrite(all_polygon_summary, out_csv)
  }

  # -----------------------------
  # Optional XLSX export (single workbook)
  # -----------------------------
  xlsx_out <- NULL
  if (isTRUE(write_xlsx) || !is.null(xlsx_path)) {

    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("To write XLSX outputs, install 'openxlsx' (suggested dependency).")
    }

    xlsx_out <- xlsx_path
    if (is.null(xlsx_out) || !nzchar(xlsx_out)) {
      tag_token <- if (!is.null(tag_clean) && nzchar(tag_clean)) paste0("_", tag_clean) else ""
      xlsx_out  <- file.path(
        validation_output_dir,
        paste0(prefix_clean, "_", year_str, "_summary", tag_token, ".xlsx")
      )
    }

    }
    dir.create(dirname(xlsx_out), recursive = TRUE, showWarnings = FALSE)

    wb <- openxlsx::createWorkbook()

    .sheet_name <- function(x) {
      x <- gsub("[\\[\\]\\*\\?/\\\\:]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- gsub("^_|_$", "", x)
      if (!nzchar(x)) x <- "sheet"
      if (nchar(x) > 31) x <- substr(x, 1, 31)
      x
    }

    add_dt_sheet <- function(sheet_name, x, withFilter = TRUE) {
      if (is.null(x)) return(invisible(FALSE))
      x <- as.data.frame(x)
      if (nrow(x) == 0) return(invisible(FALSE))

      sname <- .sheet_name(sheet_name)
      if (sname %in% names(wb)) {
        ii <- 2
        base <- substr(sname, 1, min(28, nchar(sname)))
        while (paste0(base, "_", ii) %in% names(wb)) ii <- ii + 1
        sname <- paste0(base, "_", ii)
      }

      openxlsx::addWorksheet(wb, sname)
      openxlsx::writeDataTable(wb, sname, x, withFilter = withFilter)
      openxlsx::freezePane(wb, sname, firstRow = TRUE, firstCol = FALSE)
      openxlsx::setColWidths(wb, sname, cols = 1:ncol(x), widths = "auto")
      invisible(TRUE)
    }

    add_dt_sheet("metrics", all_metrics, withFilter = TRUE)
    add_dt_sheet("polygon_summary", all_polygon_summary, withFilter = TRUE)

    cm_all  <- if (length(class_metrics_list) > 0) data.table::rbindlist(class_metrics_list, fill = TRUE) else NULL
    cmg_all <- if (length(class_metrics_global_list) > 0) data.table::rbindlist(class_metrics_global_list, fill = TRUE) else NULL
    om_all  <- if (length(omission_by_class_list) > 0) data.table::rbindlist(omission_by_class_list, fill = TRUE) else NULL
    com_all <- if (length(commission_by_class_list) > 0) data.table::rbindlist(commission_by_class_list, fill = TRUE) else NULL

    add_dt_sheet("metrics_by_class", cm_all, withFilter = TRUE)
    add_dt_sheet("metrics_by_class_global", cmg_all, withFilter = TRUE)
    add_dt_sheet("omission_by_class", om_all, withFilter = TRUE)
    add_dt_sheet("commission_by_class", com_all, withFilter = TRUE)

    cfg <- data.frame(
      Parameter = c(
        "year_target", "template_res_m",
        "min_detected_percent", "threshold_completely_detected",
        "min_area_reference_ha", "metrics_type",
        "output_prefix", "output_tag",
        "run_tag", "ref_tag",
        "mask_shapefile",
        "input_shapefile", "ref_shapefile",
        "class_raster", "class_lut",
        "digits_metrics", "digits_area", "digits_res",
        "timestamp"
      ),
      Value = c(
        as.character(year_target), as.character(template_res),
        as.character(min_detected_percent), as.character(threshold_completely_detected),
        ifelse(is.null(min_area_reference_ha), "NULL", as.character(min_area_reference_ha)),
        as.character(metrics_type),
        as.character(output_prefix), as.character(output_tag %||% ""),
        as.character(run_tag %||% ""), as.character(ref_tag %||% ""),
        ifelse(is.character(mask_shapefile), normalizePath(mask_shapefile, winslash = "/", mustWork = FALSE), "sf_object"),
        ifelse(is.character(input_shapefile), paste(input_shapefile, collapse = "; "), "sf_object"),
        ifelse(is.character(ref_shapefile), normalizePath(ref_shapefile, winslash = "/", mustWork = FALSE), "sf_object"),
        if (is.null(class_raster)) "NULL" else if (inherits(class_raster, "SpatRaster")) "SpatRaster" else normalizePath(as.character(class_raster), winslash = "/", mustWork = FALSE),
        if (is.null(class_lut)) "NULL" else if (is.data.frame(class_lut)) "data.frame" else if (is.character(class_lut)) normalizePath(as.character(class_lut), winslash = "/", mustWork = FALSE) else "object",
        as.character(digits_metrics), as.character(digits_area), as.character(digits_res),
        format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ),
      stringsAsFactors = FALSE
    )
    add_dt_sheet("config", cfg, withFilter = FALSE)

    openxlsx::saveWorkbook(wb, xlsx_out, overwrite = isTRUE(xlsx_overwrite))
    msg("Wrote XLSX: ", xlsx_out)


  return(list(
    metrics = if (metrics_type %in% c("all", "pixel")) all_metrics else NULL,
    polygon_summary = if (metrics_type %in% c("all", "area")) all_polygon_summary else NULL,
    validation_output_dir = validation_output_dir,
    errors_dir = errors_dir,
    xlsx_path = xlsx_out,
    output_prefix = prefix_clean,
    output_tag = tag_clean %||% "",
    confusion_polys_paths  = confusion_polys_paths,
    confusion_pixels_paths = confusion_pixels_paths
  ))
}
