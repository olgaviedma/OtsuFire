#' Validate burned-area maps against independent reference fire perimeters
#'
#' @description
#' `validate_fire_maps()` evaluates one or more burned-area prediction
#' layers against an independent set of reference fire perimeters.
#' It is the shared validation entry point in OtsuFire, so it can be
#' used with outputs from the deterministic workflow, the supervised
#' workflow, or any external burned-area map supplied as polygons.
#'
#' The function compares predicted and reference burned areas within a
#' common evaluation domain defined by the study-area mask and the
#' burnable-area raster. Reference polygons are first harmonised to this
#' domain by filtering them to the target year, clipping them to the
#' study area, masking them to burnable land, optionally applying a
#' minimum-area threshold, optionally excluding fires that were not
#' temporally observable, and optionally dissolving polygons by a
#' grouping field. Prediction layers are processed against the same
#' domain so that omission, commission, and agreement are evaluated
#' consistently.
#'
#' The function can compute two complementary types of validation
#' output. Pixel-based metrics quantify agreement between burned and
#' unburned pixels over the burnable domain using a confusion matrix.
#' Polygon- and area-based metrics quantify whether known reference fires
#' were detected, how much of their area was recovered, and how much
#' predicted burned area overlaps the reference data. Optional branches
#' can also summarise omission and commission by an external polygon
#' class, compute pixel-level metrics by raster strata, and write the
#' main outputs to a combined Excel workbook.
#'
#' @param input_shapefile Character vector or `sf` object. One or more
#'   burned-area prediction layers to validate. Each input is processed
#'   independently and contributes one row to the global output tables.
#' @param ref_shapefile Character path or `sf` object. Reference
#'   burned-area polygons.
#' @param mask_shapefile Character path. Study-area boundary polygon used
#'   to clip the evaluation domain.
#' @param burnable_raster Character path. Burnable-domain raster used to
#'   define the validation domain.
#' @param year_target Numeric. Target year used to filter the reference
#'   polygons to the correct fire season.
#' @param validation_dir Character. Root output directory. A
#'   `VALIDATION/` subfolder is created inside it. Main CSV and vector
#'   outputs are currently written directly in that folder.
#' @param binary_burnable Logical. If `TRUE` (default), the burnable
#'   raster is treated as binary and cells at or above
#'   `burnable_threshold` are considered burnable. If `FALSE`,
#'   `burnable_classes` defines which raster codes are burnable.
#' @param burnable_classes Optional numeric vector. Raster values treated
#'   as burnable when `binary_burnable = FALSE`.
#' @param burnable_threshold Numeric scalar in `[0, 1]`. Threshold used
#'   when `binary_burnable = TRUE`. Default `0.5`.
#' @param class_shape Optional character path. Polygon layer used for
#'   polygon-level omission / commission summaries by class. Requires
#'   `class_field`.
#' @param class_field Optional character. Attribute in `class_shape` used
#'   to group polygon-level errors.
#' @param buffer Numeric. Buffer distance in metres applied around
#'   reference polygons before pixel-based comparison. Default `0`.
#' @param threshold_completely_detected Numeric scalar in `[0, 100]`.
#'   Minimum percentage of a reference polygon that must be covered by
#'   detections to count as completely detected. Default `90`.
#' @param threshold_min_detected Numeric scalar in `[0, 100]`. Minimum
#'   percentage of a reference polygon that must be covered by detections
#'   to count as detected at all. This controls `N_Detected_Polygons`.
#'   Default `10`.
#' @param min_area_reference_ha Optional numeric. Minimum reference
#'   polygon area, in hectares after masking to the burnable domain, to
#'   retain in the analysis. Default `NULL`.
#' @param observability_raster Optional `SpatRaster`, raster path, or
#'   `NULL`. When supplied, it is aligned to the validation domain and
#'   used to exclude reference polygons whose latest observable day of
#'   year does not reach the reference fire day. The current
#'   implementation expects a single-layer DOY raster; if a multi-layer
#'   raster is supplied, only the first layer is used.
#' @param ref_end_doy_col Character. Reference attribute storing the fire
#'   end day of year. Used first to derive `obs_ref_doy`. Default
#'   `"end_doy"`.
#' @param ref_start_doy_col Character. Reference attribute storing the
#'   fire start day of year. Used as fallback when `ref_end_doy_col` is
#'   missing for a polygon. Default `"start_doy"`.
#' @param force_reprocess_ref Logical. If `TRUE`, rebuild the cached
#'   reference products even if they already exist.
#' @param force_reprocess_pred Logical. If `TRUE`, rebuild the cached
#'   prediction products even if they already exist.
#' @param metrics_type Character. One of `"all"` (default), `"pixel"`, or
#'   `"area"`.
#' @param dissolve_ref_by Optional character. Field used to dissolve
#'   reference polygons before comparison.
#' @param dissolve_input_by Optional character. Field used to dissolve
#'   prediction polygons before comparison.
#' @param strata_raster Optional `SpatRaster`, raster path, or `NULL`.
#'   When supplied, per-stratum pixel-level confusion matrices and
#'   derived metrics are computed in addition to the global outputs.
#' @param strata_lut Optional lookup table for `strata_raster`. Either a
#'   `data.frame` with columns `id` and `label`, or a CSV / XLSX path
#'   with the same schema.
#' @param chunk_rows Integer. Row chunk size used by the per-stratum
#'   tabulator. Larger values are usually faster but require more memory.
#'   Default `1024L`.
#' @param na_strata_as_zero Logical. How to treat `NA` prediction or
#'   reference pixels inside the strata domain. If `TRUE` (default),
#'   treat them as unburned (`0`) to match the legacy validator. If
#'   `FALSE`, drop those pixels.
#' @param write_excel Logical. If `TRUE`, also write a combined Excel
#'   workbook with the available validation outputs. Requires the
#'   `openxlsx` package. Default `FALSE`.
#' @param excel_filename Optional character. Custom Excel filename.
#'   Defaults to `validation_ALL_<year_target>_res30.xlsx`.
#'
#' @details
#' \strong{Conceptual workflow}
#'
#' `validate_fire_maps()` follows five main steps.
#'
#' First, it builds a common validation domain. The study-area boundary
#' defines the spatial extent, and the burnable raster defines where
#' burned/unburned comparison is meaningful. When
#' `binary_burnable = TRUE`, cells with values greater than or equal to
#' `burnable_threshold` are treated as burnable. When
#' `binary_burnable = FALSE`, the values listed in `burnable_classes`
#' define the burnable domain.
#'
#' Second, it prepares the reference layer. Reference polygons are
#' filtered to `year_target`, clipped to the study area, restricted to
#' the burnable domain and, if requested, filtered by
#' `min_area_reference_ha`. If `observability_raster` is supplied, the
#' function also removes reference polygons whose fire date falls after
#' the latest observable day of year inside the polygon. This avoids
#' penalising a prediction layer for missing fires that could not be
#' observed by the input imagery.
#'
#' Third, it prepares each prediction layer. The object supplied through
#' `input_shapefile` may be a single polygon layer, an `sf` object, or a
#' vector of polygon paths. Each prediction layer is processed
#' independently against the same reference layer and validation domain,
#' allowing several candidate maps, thresholds, or scenarios to be
#' compared in a single call.
#'
#' Fourth, it computes the requested validation metrics. When
#' `metrics_type = "pixel"` or `"all"`, the function rasterises the
#' reference and prediction layers over the burnable domain and computes
#' the pixel-level confusion matrix: true positives, false positives,
#' false negatives, and true negatives. From these counts it derives
#' standard accuracy metrics including precision, recall, F1,
#' intersection over union, specificity, balanced accuracy, and error
#' rate. These outputs are useful for quantifying overall spatial
#' agreement between predicted and reference burned pixels.
#'
#' When `metrics_type = "area"` or `"all"`, the function computes
#' polygon- and area-based validation summaries. These outputs describe
#' how many reference fires were detected, how many were completely
#' detected, how many were missed, and how much reference burned area was
#' recovered by the predictions. A reference polygon is counted as
#' detected when its percentage overlap with the prediction is greater
#' than or equal to `threshold_min_detected`. It is counted as
#' completely detected when its overlap is greater than or equal to
#' `threshold_completely_detected`. This makes the detection rule
#' explicit and avoids treating trivial overlaps as successful detections
#' unless the user chooses that behaviour.
#'
#' Fifth, the function writes the requested outputs into
#' `validation_dir/VALIDATION/`. The main global tables are currently
#' written as `metrics_summary_<year>.csv` and
#' `polygon_summary_<year>.csv`. When observability filtering is used,
#' the summary table is written as
#' `reference_observability_<year>_<tag>.csv`. When stratified metrics
#' are requested, the current file names are
#' `pixel_by_stratum_<year>_<input>.csv`,
#' `stratum_global_<year>_<input>.csv`, and
#' `diagnostics_strata_<year>_<input>.csv`. Error layers and cached
#' intermediate rasters are also written in the same `VALIDATION/`
#' folder using descriptive file names.
#'
#' \strong{Pixel-based validation}
#'
#' Pixel-based validation evaluates map agreement cell by cell over the
#' burnable domain. It is best suited for measuring the spatial match
#' between predicted and reference burned areas, including both omission
#' and commission. The main outputs are:
#' \itemize{
#'   \item `TP`: pixels mapped as burned in both prediction and
#'     reference,
#'   \item `FP`: pixels mapped as burned in the prediction but not in the
#'     reference,
#'   \item `FN`: pixels mapped as burned in the reference but missed by
#'     the prediction,
#'   \item `TN`: pixels mapped as unburned in both prediction and
#'     reference.
#' }
#'
#' The function then derives accuracy metrics such as `Precision`,
#' `Recall`, `F1`, `IoU`, `Specificity`, `BalancedAccuracy`, and
#' `ErrorRate`. These metrics are returned in `metrics` and written to
#' `validation_dir/VALIDATION/metrics_summary_<year>.csv`.
#'
#' \strong{Polygon- and area-based validation}
#'
#' Polygon- and area-based validation evaluates detection at the
#' fire-event level. This branch answers questions such as: How many
#' reference fires were detected? How many were missed? How much
#' reference burned area was recovered? How much predicted area overlaps
#' reference perimeters?
#'
#' The main outputs include `N_Reference_Polygons`,
#' `N_Detected_Polygons`, `N_Completely_Detected`, `N_Not_Detected`,
#' `Area_Reference_ha`, `Area_Detected_ha`, `Area_Intersection_ha`,
#' `Recall_Area_percent`, `Precision_Area_percent`, and coverage
#' summaries. The field `Detected_Definition` records the rule used to
#' count detected polygons, namely
#' `coverage_ref >= threshold_min_detected`.
#'
#' By default, `threshold_min_detected = 10`, meaning that at least 10\%
#' of a reference polygon must overlap the prediction to be counted as
#' detected. This is stricter than the earlier legacy behaviour, which
#' counted any non-zero overlap as detection. To recover that older
#' behaviour, set `threshold_min_detected = 0`.
#'
#' \strong{Temporal observability filter}
#'
#' When `observability_raster` is supplied, `validate_fire_maps()`
#' checks whether each reference fire was observable in the imagery used
#' for mapping. The function compares the fire day of year, derived from
#' `ref_end_doy_col` with fallback to `ref_start_doy_col`, against the
#' maximum observable day of year inside each reference polygon.
#'
#' For each reference polygon, the function records `obs_ref_doy`,
#' `obs_doy_max`, `obs_doy_margin`, `observable_flag`, and
#' `observable_reason`. Non-observable reference polygons are excluded
#' from the validated reference set, and the observability summary table
#' is returned as `reference_observability` and written to
#' `validation_dir/VALIDATION/reference_observability_<year>_<tag>.csv`.
#' This is useful when the prediction map is based on annual or seasonal
#' composites that may not contain valid post-fire observations for all
#' reference events.
#'
#' \strong{Polygon-level omission and commission by class}
#'
#' If `class_shape` and `class_field` are supplied, the function
#' summarises omitted reference polygons and commission polygons by an
#' external polygon classification. This is useful for diagnosing
#' whether errors are concentrated in particular vegetation types,
#' regions, administrative units, or other user-defined classes.
#'
#' This branch is polygon-based. It does not produce a pixel-level
#' confusion matrix by class. For pixel-level stratified validation, use
#' `strata_raster`.
#'
#' \strong{Pixel-level validation by strata}
#'
#' If `strata_raster` is supplied, the function computes a pixel-level
#' confusion matrix separately within each raster stratum. This can be
#' used to evaluate map performance by land-cover class, ecoregion,
#' vegetation group, or any other categorical raster layer.
#'
#' The optional `strata_lut` argument can be used to attach
#' human-readable labels to stratum IDs. Pixels with missing strata are
#' excluded from the stratified aggregation. Within valid strata,
#' missing prediction or reference pixels are treated according to
#' `na_strata_as_zero`: if `TRUE`, they are treated as unburned to match
#' the legacy validator; if `FALSE`, they are dropped from the stratified
#' calculation.
#'
#' \strong{Caching and reprocessing}
#'
#' The function caches intermediate masked reference and prediction
#' layers in `validation_dir/VALIDATION/`. This avoids repeating
#' expensive spatial preprocessing when the same validation setup is run
#' multiple times.
#'
#' Use `force_reprocess_ref = TRUE` when changes to the reference-side
#' settings require rebuilding the reference cache, especially after
#' changing `min_area_reference_ha`, `observability_raster`,
#' `dissolve_ref_by`, or the burnable-domain settings. Use
#' `force_reprocess_pred = TRUE` when the prediction layer or
#' prediction-side settings have changed and the cached prediction
#' products should be rebuilt.
#'
#' @return A named list containing the available validation outputs:
#' \describe{
#'   \item{`metrics`}{Global pixel-based metrics, or `NULL` when
#'     `metrics_type = "area"`.}
#'   \item{`polygon_summary`}{Global polygon- and area-based metrics, or
#'     `NULL` when `metrics_type = "pixel"`.}
#'   \item{`reference_observability`}{Reference-polygon observability
#'     diagnostics, or `NULL` when `observability_raster = NULL`.}
#'   \item{`pixel_by_stratum`}{Per-stratum pixel metrics, or `NULL` when
#'     `strata_raster` is not supplied.}
#'   \item{`stratum_global`}{Aggregate metrics over the stratified
#'     validation domain, or `NULL` when `strata_raster` is not
#'     supplied.}
#'   \item{`diagnostics_strata`}{Diagnostics for the stratified branch,
#'     or `NULL` when `strata_raster` is not supplied.}
#'   \item{`excel_path`}{Path to the combined Excel workbook, or `NULL`
#'     when `write_excel = FALSE`.}
#' }
#'
#' @note
#' This is the canonical v2 validation implementation. Deprecated v1
#' GDAL / Python parameters are no longer part of the public interface.
#'
#' Version 0.2.1 added `Specificity`, `BalancedAccuracy`, and
#' `ErrorRate` to the global and stratified pixel outputs; coverage
#' summaries plus `Detected_Definition` to the polygon outputs; and the
#' `threshold_min_detected`, `dissolve_ref_by`, and `dissolve_input_by`
#' controls. Version 0.5.0 adds the optional temporal observability
#' filter controlled by `observability_raster`.
#'
#' @seealso [run_deterministic_pipeline()] for the deterministic
#'   workflow that commonly produces the prediction layer passed to
#'   `input_shapefile`; [run_oneyear_supervised_pipeline()] for the
#'   supervised workflow whose thresholded `final_map.gpkg` is also a
#'   valid input here. `validate_fire_maps()` performs EXTERNAL,
#'   POST-RUN cartographic validation of PRODUCED maps against
#'   independent reference perimeters (e.g. EFFIS); for the
#'   complementary PRE-RUN validator that audits the supervised
#'   CONFIGURATION + EXECUTION (CRS / overlap / year / cfg
#'   contradictions / cache provenance) BEFORE any map is produced, see
#'   [validate_supervised_execution()].
#'
#' @examples
#' \dontrun{
#' # Minimal call: global pixel + polygon metrics only.
#' validate_fire_maps(
#'   input_shapefile = list.files("shapefiles", pattern = "\\.shp$", full.names = TRUE),
#'   ref_shapefile   = "ref_polygons_2022.shp",
#'   mask_shapefile  = "mask_region.shp",
#'   burnable_raster = "burnable.tif",
#'   year_target     = 2022,
#'   validation_dir  = "validation_results"
#' )
#'
#' # With pixel-level CORINE stratification + combined Excel workbook.
#' # Produces 03_STRATA CSVs plus the optional Excel workbook in
#' # 01_SUMMARY/validation_ALL_2022_res30.xlsx.
#' validate_fire_maps(
#'   input_shapefile = "predicted_2022.gpkg",
#'   ref_shapefile   = "ref_polygons_2022.shp",
#'   mask_shapefile  = "mask_region.shp",
#'   burnable_raster = "burnable.tif",
#'   year_target     = 2022,
#'   validation_dir  = "validation_results",
#'   strata_raster   = "strata_CLC_2018_res30.tif",
#'   strata_lut      = "lut_full_strata8_v1.csv",
#'   chunk_rows      = 1024L,
#'   na_strata_as_zero = TRUE,
#'   write_excel     = TRUE
#' )
#' }
#'
#' @importFrom sf st_read st_write st_crs st_transform st_make_valid st_intersection
#' @importFrom sf st_area st_buffer st_is_empty st_zm st_filter st_drop_geometry
#' @importFrom sf st_is_valid st_join sf_use_s2
#' @importFrom terra rast crop mask project rasterize vect writeRaster global extract res
#' @importFrom data.table data.table fwrite rbindlist
#' @importFrom dplyr filter group_by summarise mutate n
#' @importFrom tools file_path_sans_ext
#' @family modular
#' @export
validate_fire_maps <- function(input_shapefile,
                               ref_shapefile,
                               mask_shapefile,
                               burnable_raster,
                               year_target,
                               validation_dir,
                               binary_burnable = TRUE,
                               burnable_classes = NULL,
                               burnable_threshold = 0.5,
                               class_shape = NULL,
                               class_field = NULL,
                               buffer = 0,
                               threshold_completely_detected = 90,
                               threshold_min_detected = 10,
                               min_area_reference_ha = NULL,
                               observability_raster = NULL,
                               ref_end_doy_col = "end_doy",
                               ref_start_doy_col = "start_doy",
                               force_reprocess_ref = FALSE,
                               force_reprocess_pred = FALSE,
                               metrics_type = c("all", "pixel", "area"),
                               dissolve_ref_by = NULL,
                               dissolve_input_by = NULL,
                               strata_raster = NULL,
                               strata_lut = NULL,
                               chunk_rows = 1024L,
                               na_strata_as_zero = TRUE,
                               write_excel = FALSE,
                               excel_filename = NULL) {
  metrics_type <- match.arg(metrics_type)

  # s2 OFF (avoids wk/loops errors on heavily-invalid geometries)
  old_s2 <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  # ---- helpers ----
  make_valid_sf <- function(x) {
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")
    if (exists("st_make_valid", where = asNamespace("sf"), inherits = FALSE)) {
      x <- sf::st_make_valid(x)
    } else {
      x <- suppressWarnings(sf::st_buffer(x, 0))
    }
    x <- x[!sf::st_is_empty(x), ]
    x
  }

  safe_vect <- function(sfobj, label = "") {
    tryCatch(
      terra::vect(sfobj),
      error = function(e) {
        sfobj2 <- make_valid_sf(sfobj)
        ok <- sf::st_is_valid(sfobj2)
        if (any(!ok)) {
          warning(sprintf("(%s) Dropping %d invalid geometries.", label, sum(!ok)))
          sfobj2 <- sfobj2[ok, ]
        }
        terra::vect(sfobj2)
      }
    )
  }

  read_sf_any <- function(x) {
    if (inherits(x, "sf")) return(x)
    if (!is.character(x) || length(x) != 1) stop("ref_shapefile must be a path (length 1) or an sf.")
    sf::st_read(x, quiet = TRUE)
  }

  dissolve_by_field <- function(x, field) {
    if (is.null(field)) return(x)
    if (!field %in% names(x)) stop("Dissolve field does not exist: ", field)
    x |>
      dplyr::group_by(.data[[field]]) |>
      dplyr::summarise(do_union = TRUE, .groups = "drop")
  }

  # Content-aware cache keys live in R/internal-validate-cache.R as
  # package-internal, side-effect-free helpers (so the cache contract is unit
  # testable without running the full validator):
  #   .vfm_observability_fingerprint(), .vfm_domain_fingerprint(),
  #   .vfm_reference_cache_key(), .vfm_vector_fingerprint().
  # They fold the CONTENT and methodological identity of every input that shapes
  # the reference cache (observability raster content + layer + DOY columns +
  # rule version; burnable raster grid/CRS/content; study-area mask; burnable
  # thresholding; reference vector content; min-area / dissolve options), so a
  # changed input busts the cache automatically. force_reprocess_ref = TRUE
  # still forces a full rebuild regardless of the key.

  align_optional_raster <- function(x, template, mask_v) {
    xr <- if (inherits(x, "SpatRaster")) x else terra::rast(x)
    if (terra::nlyr(xr) > 1L) {
      nm <- tolower(names(xr))
      doy_idx <- which(nm == "doy" | grepl("(^|_)doy($|_)", nm))
      if (!length(doy_idx)) {
        stop(
          "observability_raster has multiple layers but none is named 'doy'. ",
          "Provide a single-layer DOY raster or rename the DOY layer to 'doy'.",
          call. = FALSE
        )
      }
      xr <- xr[[doy_idx[[1]]]]
    }
    if (terra::crs(xr) != terra::crs(template)) {
      xr <- terra::project(xr, terra::crs(template), method = "near")
    }
    if (!terra::compareGeom(xr, template, stopOnError = FALSE)) {
      xr <- terra::resample(xr, template, method = "near")
    }
    xr <- terra::crop(xr, mask_v)
    xr <- terra::mask(xr, mask_v)
    xr[is.na(template)] <- NA
    xr
  }

  build_reference_observability <- function(ref_polygons, ref_mask_r,
                                            observability_raster, domain_mask,
                                            mask_v, cell_area_ha,
                                            ref_end_doy_col,
                                            ref_start_doy_col) {
    has_end_doy <- ref_end_doy_col %in% names(ref_polygons)
    has_start_doy <- ref_start_doy_col %in% names(ref_polygons)
    if (!has_end_doy && !has_start_doy) {
      stop(
        "Reference DOY columns not found. Expected at least one of: ",
        ref_end_doy_col, ", ", ref_start_doy_col,
        call. = FALSE
      )
    }

    obs_r <- align_optional_raster(observability_raster, domain_mask, mask_v)
    ref_v <- safe_vect(ref_polygons, label = "ref_observability")

    obs_extract <- terra::extract(obs_r, ref_v, fun = max, na.rm = TRUE)
    obs_doy_max <- obs_extract[, 2]
    obs_doy_max[!is.finite(obs_doy_max)] <- NA_real_

    ref_area_pix <- terra::extract(ref_mask_r, ref_v, fun = sum, na.rm = TRUE)[, 2]
    ref_area_pix[is.na(ref_area_pix)] <- 0
    ref_area_domain_ha <- ref_area_pix * cell_area_ha

    # Per-fire PARTIAL-observability audit counts (informative only; they do NOT
    # change the whole-fire observability rule, the reference filtering, or any
    # metric). `obs_r` is already aligned to and NA outside the burnable domain,
    # so a finite obs value marks a domain cell that carries an observation.
    # total_pixels      = burnable-domain cells under the fire (= ref_area_pix);
    # observable_pixels = those that ALSO carry a finite observation DOY;
    # non_observable_pixels / observable_fraction follow.
    obs_finite_r <- terra::ifel(is.finite(obs_r), 1, 0)
    obs_finite_r[is.na(domain_mask)] <- NA
    obs_valid_pix <- terra::extract(obs_finite_r, ref_v, fun = sum, na.rm = TRUE)[, 2]
    obs_valid_pix[is.na(obs_valid_pix)] <- 0
    total_pixels <- ref_area_pix
    observable_pixels <- pmin(obs_valid_pix, total_pixels)
    non_observable_pixels <- pmax(total_pixels - observable_pixels, 0)
    observable_fraction <- ifelse(total_pixels > 0,
                                  observable_pixels / total_pixels, NA_real_)

    end_doy <- if (has_end_doy) {
      suppressWarnings(as.numeric(ref_polygons[[ref_end_doy_col]]))
    } else {
      rep(NA_real_, nrow(ref_polygons))
    }
    start_doy <- if (has_start_doy) {
      suppressWarnings(as.numeric(ref_polygons[[ref_start_doy_col]]))
    } else {
      rep(NA_real_, nrow(ref_polygons))
    }
    obs_ref_doy <- ifelse(is.finite(end_doy), end_doy, start_doy)
    obs_ref_doy[!is.finite(obs_ref_doy)] <- NA_real_
    obs_ref_doy_source <- ifelse(
      is.finite(end_doy), ref_end_doy_col,
      ifelse(is.finite(start_doy), ref_start_doy_col, NA_character_)
    )

    obs_doy_margin <- obs_doy_max - obs_ref_doy

    observable_reason <- ifelse(
      is.na(obs_ref_doy), "missing_reference_doy",
      ifelse(
        is.na(obs_doy_max), "no_observability_data",
        ifelse(obs_doy_margin < 0, "obs_before_reference_doy", "observable")
      )
    )
    observable_flag <- observable_reason == "observable"

    out <- data.table::as.data.table(sf::st_drop_geometry(ref_polygons))
    out[, reference_row := seq_len(.N)]
    out[, ref_area_domain_ha := round(ref_area_domain_ha, 4)]
    out[, obs_ref_doy := obs_ref_doy]
    out[, obs_ref_doy_source := obs_ref_doy_source]
    out[, obs_doy_max := obs_doy_max]
    out[, obs_doy_margin := obs_doy_margin]
    out[, observable_flag := observable_flag]
    out[, observable_reason := observable_reason]
    out[, total_pixels := as.integer(round(total_pixels))]
    out[, observable_pixels := as.integer(round(observable_pixels))]
    out[, non_observable_pixels := as.integer(round(non_observable_pixels))]
    out[, observable_fraction := round(observable_fraction, 6)]
    out
  }

  warn_reference_observability <- function(reference_observability) {
    if (is.null(reference_observability) || !nrow(reference_observability)) return(invisible(NULL))
    dropped <- reference_observability[observable_flag == FALSE]
    if (!nrow(dropped)) return(invisible(NULL))
    reason_counts <- dropped[, .N, by = observable_reason][order(-N)]
    reason_text <- paste(sprintf("%s=%d", reason_counts$observable_reason, reason_counts$N), collapse = ", ")
    warning(
      sprintf(
        "Temporal observability filter excluded %d of %d reference polygons (%s).",
        nrow(dropped),
        nrow(reference_observability),
        reason_text
      ),
      call. = FALSE
    )
    invisible(NULL)
  }

  build_not_observable_reference_polygons <- function(ref_polygons,
                                                      reference_observability) {
    if (is.null(reference_observability) || !nrow(reference_observability)) return(NULL)
    if (nrow(ref_polygons) != nrow(reference_observability)) {
      stop(
        "Internal error: reference geometry and observability rows do not align.",
        call. = FALSE
      )
    }

    dropped_idx <- which(reference_observability$observable_flag == FALSE)
    if (!length(dropped_idx)) return(NULL)

    out <- ref_polygons[dropped_idx, , drop = FALSE]
    obs_cols <- c(
      "obs_ref_doy",
      "obs_ref_doy_source",
      "obs_doy_max",
      "obs_doy_margin",
      "observable_flag",
      "observable_reason"
    )
    for (col_nm in obs_cols) {
      out[[col_nm]] <- reference_observability[[col_nm]][dropped_idx]
    }
    out
  }

  # ---- output dir ----
  validation_output_dir <- file.path(validation_dir, "VALIDATION")
  if (!dir.exists(validation_output_dir)) dir.create(validation_output_dir, recursive = TRUE)
  summary_output_dir <- file.path(validation_output_dir, "01_SUMMARY")
  observability_output_dir <- file.path(validation_output_dir, "02_OBSERVABILITY")
  strata_output_dir <- file.path(validation_output_dir, "03_STRATA")
  error_output_dir <- file.path(validation_output_dir, "04_ERROR_LAYERS")
  cache_output_dir <- file.path(validation_output_dir, "_CACHE")
  legacy_output_dir <- file.path(validation_output_dir, "_LEGACY_FLAT_OUTPUTS")
  for (dir_path in c(
    summary_output_dir,
    observability_output_dir,
    strata_output_dir,
    error_output_dir,
    cache_output_dir,
    legacy_output_dir
  )) {
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  }

  validation_output_path <- function(section_dir, filename) {
    file.path(section_dir, filename)
  }

  next_legacy_archive_path <- function(target_dir, filename) {
    ext <- sub("^.*(\\.[^.]+)$", "\\1", filename)
    if (identical(ext, filename)) ext <- ""
    stem <- if (nzchar(ext)) substr(filename, 1L, nchar(filename) - nchar(ext)) else filename
    candidate <- file.path(target_dir, filename)
    idx <- 1L
    while (file.exists(candidate)) {
      candidate <- file.path(target_dir, sprintf("%s__legacy%02d%s", stem, idx, ext))
      idx <- idx + 1L
    }
    candidate
  }

  move_legacy_flat_outputs <- function() {
    legacy_patterns <- c(
      "^metrics_summary_.*\\.csv$",
      "^polygon_summary_.*\\.csv$",
      "^pixel_by_stratum_.*\\.csv$",
      "^stratum_global_.*\\.csv$",
      "^diagnostics_strata_.*\\.csv$",
      "^pixel_metrics_global_.*\\.csv$",
      "^fire_metrics_global_.*\\.csv$",
      "^pixel_metrics_by_stratum_.*\\.csv$",
      "^pixel_metrics_strata_global_.*\\.csv$",
      "^strata_diagnostics_.*\\.csv$",
      "^pred_.*\\.tif$",
      "^predicted_fire_mask_.*\\.tif$",
      "^ref_mask_.*\\.tif$",
      "^reference_fire_mask_.*\\.tif$",
      "^ref_polygons_processed_.*\\.gpkg$",
      "^reference_fires_processed_.*\\.gpkg$",
      "^reference_observability_.*\\.csv$",
      "^reference_fires_observability_.*\\.csv$",
      "^reference_polygons_not_observable_.*\\.gpkg$",
      "^reference_fires_not_observable_.*\\.gpkg$",
      "^ref_polygons_not_detected_.*\\.gpkg$",
      "^reference_fires_omission_.*\\.gpkg$",
      "^input_polygons_not_matched_.*\\.gpkg$",
      "^predicted_fires_commission_.*\\.gpkg$",
      "^omission_by_.*\\.csv$",
      "^commission_by_.*\\.csv$",
      "^reference_fires_omission_by_.*\\.csv$",
      "^predicted_fires_commission_by_.*\\.csv$",
      "^validation_ALL_.*\\.xlsx$"
    )
    root_entries <- list.files(
      validation_output_dir,
      full.names = TRUE,
      recursive = FALSE,
      all.files = FALSE,
      no.. = TRUE
    )
    if (!length(root_entries)) return(invisible(NULL))
    root_files <- root_entries[file.info(root_entries)$isdir == FALSE]
    if (!length(root_files)) return(invisible(NULL))

    base_names <- basename(root_files)
    to_move <- root_files[vapply(
      base_names,
      function(x) any(grepl(paste(legacy_patterns, collapse = "|"), x, perl = TRUE)),
      logical(1)
    )]
    if (!length(to_move)) return(invisible(NULL))

    for (src in to_move) {
      dst <- next_legacy_archive_path(legacy_output_dir, basename(src))
      ok <- file.rename(src, dst)
      if (!ok) {
        ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
        if (ok) unlink(src, force = TRUE)
      }
      if (!ok) {
        warning(
          "Could not move legacy validation output to _LEGACY_FLAT_OUTPUTS/: ",
          basename(src),
          call. = FALSE
        )
      }
    }
    invisible(NULL)
  }

  move_legacy_flat_outputs()

  # ---- load mask & burnable ----
  mask_geom <- sf::st_read(mask_shapefile, quiet = TRUE) |> make_valid_sf()
  # F8b + §N+7.6 (2026-05-29): drop ALL non-geometry attributes from
  # mask_geom. We only need its geometry for clipping ref_polygons and
  # det downstream via st_intersection; mask attributes are never
  # consumed by the validator. Keeping them is unsafe because
  # st_intersection appends them to ref_polygons/det, and the GPKG cache
  # write at line ~921 then fails with
  # "Field count reached: duplicate names present?" whenever an
  # attribute name on mask collides with one on the reference under
  # SQLite's case-insensitive identifier rules (e.g., reference column
  # `id` vs Iberian peninsula mask column `Id`). The previous narrower
  # `!names %in% "FID"` drop covered only one GDAL-generated case
  # (F8b, "Wrong field type for FID"); generalising to geometry-only
  # is the proper fix.
  geom_col <- attr(mask_geom, "sf_column")
  mask_geom <- mask_geom[, geom_col, drop = FALSE]
  burnable  <- terra::rast(burnable_raster)

  if (terra::crs(burnable) != sf::st_crs(mask_geom)$wkt) {
    burnable <- terra::project(burnable, sf::st_crs(mask_geom)$wkt)
  }

  # ---- template domain (mask + burnable) ----
  mask_v <- terra::vect(mask_geom)
  burnable_m <- terra::crop(burnable, mask_v) |> terra::mask(mask_v)

  template_domain <- burnable_m

  if (binary_burnable) {
    template_domain[!is.na(template_domain) & template_domain < burnable_threshold] <- NA
  } else {
    if (is.null(burnable_classes)) stop("binary_burnable=FALSE requires burnable_classes.")
    template_domain[!(template_domain %in% burnable_classes)] <- NA
  }

  domain_mask <- template_domain
  domain_mask[!is.na(domain_mask)] <- 1

  cell_area_ha <- abs(prod(terra::res(domain_mask))) / 10000
  obs_cache_tag <- .vfm_observability_fingerprint(
    observability_raster = observability_raster,
    ref_end_doy_col = ref_end_doy_col,
    ref_start_doy_col = ref_start_doy_col
  )

  # Burnable-domain identity (grid/CRS/extent/content of the burnable raster +
  # study-area mask + burnable thresholding). Shared by the reference and the
  # prediction caches, since both are rasterized onto / masked by this domain.
  dom_cache_tag <- .vfm_domain_fingerprint(
    burnable_raster = burnable_raster,
    mask_shapefile = mask_shapefile,
    binary_burnable = binary_burnable,
    burnable_classes = burnable_classes,
    burnable_threshold = burnable_threshold
  )

  # Content-aware reference cache key. Folds the CONTENT and methodological
  # identity of every input that shapes the cached reference artifact:
  # observability (raster content + layer + DOY columns + rule version),
  # burnable domain (grid/CRS/content + mask + thresholding) and the reference
  # vector content + reference-only options (min area, dissolve field). A change
  # to any of them busts the cache automatically, without requiring
  # force_reprocess_ref = TRUE. `buffer` is intentionally excluded: it is
  # applied downstream to the loaded reference at detection time and never
  # shapes the cached artifact.
  ref_cache_tag <- .vfm_reference_cache_key(
    obs_cache_tag = obs_cache_tag,
    dom_cache_tag = dom_cache_tag,
    ref_shapefile = ref_shapefile,
    min_area_reference_ha = min_area_reference_ha,
    dissolve_ref_by = dissolve_ref_by
  )

  # ---- cache paths ----
  ref_vec_cache  <- validation_output_path(
    cache_output_dir,
    paste0("reference_fires_processed_", year_target, "_", ref_cache_tag, ".gpkg")
  )
  ref_rast_cache <- validation_output_path(
    cache_output_dir,
    paste0("reference_fire_mask_", year_target, "_", ref_cache_tag, ".tif")
  )
  ref_obs_cache  <- if (!is.null(observability_raster)) {
    validation_output_path(
      observability_output_dir,
      paste0("reference_fires_observability_", year_target, "_", ref_cache_tag, ".csv")
    )
  } else {
    NULL
  }
  ref_not_obs_cache <- if (!is.null(observability_raster)) {
    validation_output_path(
      observability_output_dir,
      paste0("reference_fires_not_observable_", year_target, "_", ref_cache_tag, ".gpkg")
    )
  } else {
    NULL
  }

  if (force_reprocess_ref) {
    message("Force reprocessing reference polygons: deleting cached reference...")
    unlink(ref_vec_cache)
    unlink(ref_rast_cache)
    if (!is.null(ref_obs_cache)) unlink(ref_obs_cache)
    if (!is.null(ref_not_obs_cache)) unlink(ref_not_obs_cache)
  }

  # ---- build/load reference (vector + raster mask) ----
  reference_observability <- NULL
  ref_not_observable <- NULL
  observability_audit <- NULL
  observability_applied <- !is.null(observability_raster)
  n_reference_excluded_observability <- 0L
  n_reference_observable <- NA_integer_
  use_cached_reference <- file.exists(ref_vec_cache) &&
    file.exists(ref_rast_cache) &&
    (is.null(ref_obs_cache) || file.exists(ref_obs_cache))

  if (use_cached_reference && !is.null(ref_obs_cache) && file.exists(ref_obs_cache)) {
    reference_observability_meta <- data.table::fread(
      ref_obs_cache,
      select = "observable_flag"
    )
    n_excluded_cached <- sum(reference_observability_meta$observable_flag == FALSE, na.rm = TRUE)
    if (n_excluded_cached > 0L && !file.exists(ref_not_obs_cache)) {
      message(
        "Cached observability diagnostics found without excluded-reference GeoPackage; rebuilding reference cache..."
      )
      use_cached_reference <- FALSE
    }
  }

  if (use_cached_reference) {
    message("Loading cached reference (vector + raster mask)...")
    ref_polygons <- sf::st_read(ref_vec_cache, quiet = TRUE) |> make_valid_sf()
    ref_mask_r   <- terra::rast(ref_rast_cache)
    if (!is.null(ref_obs_cache) && file.exists(ref_obs_cache)) {
      reference_observability <- data.table::fread(ref_obs_cache)
      n_reference_excluded_observability <- sum(!reference_observability$observable_flag, na.rm = TRUE)
      n_reference_observable <- sum(reference_observability$observable_flag, na.rm = TRUE)
      warn_reference_observability(reference_observability)
    }
  } else {
    message("Processing reference polygons...")
    ref_polygons <- read_sf_any(ref_shapefile) |> make_valid_sf()

    if ("year" %in% names(ref_polygons)) {
      ref_polygons <- dplyr::filter(ref_polygons, .data$year == year_target)
      message("Reference polygons filtered by year (", year_target, "): ", nrow(ref_polygons))
    }
    if (nrow(ref_polygons) == 0) stop("No reference polygons for year_target.")

    if (sf::st_crs(ref_polygons) != sf::st_crs(mask_geom)) {
      ref_polygons <- sf::st_transform(ref_polygons, sf::st_crs(mask_geom))
    }

    ref_polygons <- sf::st_filter(ref_polygons, mask_geom, .predicate = sf::st_intersects)
    message("Reference polygons after mask filter: ", nrow(ref_polygons))
    if (nrow(ref_polygons) == 0) stop("No reference polygons inside mask.")
    ref_polygons <- suppressWarnings(sf::st_intersection(ref_polygons, mask_geom)) |> make_valid_sf()
    if (is.null(observability_raster)) {
      ref_polygons <- dissolve_by_field(ref_polygons, dissolve_ref_by) |> make_valid_sf()
    }

    ref_v <- safe_vect(ref_polygons, label = "ref")
    ref_mask_r <- terra::rasterize(ref_v, domain_mask, field = 1, background = 0)
    ref_mask_r[is.na(domain_mask)] <- NA

    if (!is.null(min_area_reference_ha)) {
      ext <- terra::extract(ref_mask_r, ref_v, fun = sum, na.rm = TRUE)
      pix <- ext[, 2]
      pix[is.na(pix)] <- 0
      area_i <- pix * cell_area_ha

      keep <- area_i >= min_area_reference_ha
      message(sprintf(
        "Filtered small reference polygons: %d -> %d (Area >= %.2f ha)",
        length(keep), sum(keep), min_area_reference_ha
      ))

      ref_polygons <- ref_polygons[keep, , drop = FALSE]
      ref_polygons <- make_valid_sf(ref_polygons)
      if (nrow(ref_polygons) == 0) stop("No reference polygons remain after min_area_reference_ha.")

      ref_v <- safe_vect(ref_polygons, label = "ref_filtered")
      ref_mask_r <- terra::rasterize(ref_v, domain_mask, field = 1, background = 0)
      ref_mask_r[is.na(domain_mask)] <- NA
    }

    if (!is.null(observability_raster)) {
      reference_observability <- build_reference_observability(
        ref_polygons = ref_polygons,
        ref_mask_r = ref_mask_r,
        observability_raster = observability_raster,
        domain_mask = domain_mask,
        mask_v = mask_v,
        cell_area_ha = cell_area_ha,
        ref_end_doy_col = ref_end_doy_col,
        ref_start_doy_col = ref_start_doy_col
      )
      ref_not_observable <- build_not_observable_reference_polygons(
        ref_polygons = ref_polygons,
        reference_observability = reference_observability
      )
      n_reference_excluded_observability <- sum(!reference_observability$observable_flag, na.rm = TRUE)
      n_reference_observable <- sum(reference_observability$observable_flag, na.rm = TRUE)
      warn_reference_observability(reference_observability)

      ref_polygons <- ref_polygons[reference_observability$observable_flag, , drop = FALSE]
      ref_polygons <- make_valid_sf(ref_polygons)
      if (nrow(ref_polygons) == 0) {
        stop("No observable reference polygons remain after temporal observability filtering.",
             call. = FALSE)
      }

      ref_v <- safe_vect(ref_polygons, label = "ref_observable")
      ref_mask_r <- terra::rasterize(ref_v, domain_mask, field = 1, background = 0)
      ref_mask_r[is.na(domain_mask)] <- NA
    }

    if (!is.null(observability_raster)) {
      ref_polygons <- dissolve_by_field(ref_polygons, dissolve_ref_by) |> make_valid_sf()
    }
    if (!is.null(dissolve_ref_by) && !is.null(observability_raster)) {
      ref_v <- safe_vect(ref_polygons, label = "ref_dissolved")
      ref_mask_r <- terra::rasterize(ref_v, domain_mask, field = 1, background = 0)
      ref_mask_r[is.na(domain_mask)] <- NA
    }

    sf::st_write(ref_polygons, ref_vec_cache, delete_dsn = TRUE, quiet = TRUE)
    terra::writeRaster(
      ref_mask_r, ref_rast_cache, overwrite = TRUE,
      datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255")
    )
    if (!is.null(ref_obs_cache) && !is.null(reference_observability)) {
      data.table::fwrite(reference_observability, ref_obs_cache)
    }
    if (!is.null(ref_not_obs_cache)) {
      if (is.null(ref_not_observable) || !nrow(ref_not_observable)) {
        if (file.exists(ref_not_obs_cache)) unlink(ref_not_obs_cache)
      } else {
        sf::st_write(ref_not_observable, ref_not_obs_cache, delete_dsn = TRUE, quiet = TRUE)
      }
    }
  }

  if (is.null(reference_observability) && !is.null(observability_raster)) {
    n_reference_observable <- nrow(ref_polygons)
  }

  # ---- partial-observability audit (INFORMATIVE ONLY) ----
  # Summarizes excluded fires and, within KEPT fires, partial spatial
  # observability. Writes three CSVs to 02_OBSERVABILITY and returns the audit
  # on the result. It never alters TP/FP/FN/TN/coverage or the reference filter.
  if (!is.null(reference_observability) && nrow(reference_observability) > 0L) {
    observability_audit <- .vfm_observability_audit(
      reference_observability = reference_observability,
      cell_area_ha = cell_area_ha
    )
    data.table::fwrite(
      observability_audit$summary,
      validation_output_path(observability_output_dir,
        sprintf("OBSERVABILITY_AUDIT_%s.csv", year_target))
    )
    data.table::fwrite(
      observability_audit$excluded,
      validation_output_path(observability_output_dir,
        sprintf("OBSERVABILITY_EXCLUDED_FIRES_%s.csv", year_target))
    )
    data.table::fwrite(
      observability_audit$partial,
      validation_output_path(observability_output_dir,
        sprintf("OBSERVABILITY_PARTIAL_FIRES_%s.csv", year_target))
    )
  }

  # ---- normalize input list ----
  input_list <- NULL
  if (inherits(input_shapefile, "sf")) {
    input_list <- list(input_shapefile)
    input_names <- "input_sf"
  } else {
    if (!is.character(input_shapefile)) stop("input_shapefile must be path(s) or sf.")
    input_list <- as.list(input_shapefile)
    input_names <- tools::file_path_sans_ext(basename(unlist(input_list)))
  }

  metrics_list <- list()
  polygon_summary_list <- list()
  strata_table_list  <- list()
  strata_global_list <- list()
  strata_diag_list   <- list()

  # ---- prepare strata raster + LUT (Block 10c) ----
  strata_r <- NULL
  strata_lut_df <- NULL
  if (!is.null(strata_raster)) {
    strata_r <- if (inherits(strata_raster, "SpatRaster"))
      strata_raster else terra::rast(strata_raster)
    if (terra::nlyr(strata_r) > 1L) strata_r <- strata_r[[1]]
    if (terra::crs(strata_r) != terra::crs(domain_mask)) {
      strata_r <- terra::project(strata_r, terra::crs(domain_mask),
                                  method = "near")
    }
    if (!terra::compareGeom(strata_r, domain_mask, stopOnError = FALSE)) {
      strata_r <- terra::resample(strata_r, domain_mask, method = "near")
    }
    if (!is.null(strata_lut)) {
      if (is.data.frame(strata_lut)) {
        strata_lut_df <- as.data.frame(strata_lut)
      } else if (is.character(strata_lut) && length(strata_lut) == 1L &&
                 file.exists(strata_lut)) {
        ext_lut <- tolower(tools::file_ext(strata_lut))
        strata_lut_df <- if (ext_lut %in% c("xlsx", "xls")) {
          if (!requireNamespace("openxlsx", quietly = TRUE))
            stop("Package 'openxlsx' is required to read XLSX LUT files. ",
                 "Install it with install.packages('openxlsx').", call. = FALSE)
          openxlsx::read.xlsx(strata_lut)
        } else {
          tryCatch(utils::read.csv2(strata_lut, stringsAsFactors = FALSE,
                                     fileEncoding = "UTF-8"),
                   error = function(e)
                     utils::read.csv(strata_lut, stringsAsFactors = FALSE))
        }
      } else {
        stop("'strata_lut' must be a data.frame or a path to a CSV/XLSX.",
             call. = FALSE)
      }
      if (!all(c("id", "label") %in% names(strata_lut_df))) {
        stop("'strata_lut' must have columns 'id' and 'label'.",
             call. = FALSE)
      }
      strata_lut_df$id    <- as.integer(strata_lut_df$id)
      strata_lut_df$label <- as.character(strata_lut_df$label)
    }
  }

  # ---- loop inputs ----
  for (k in seq_along(input_list)) {

    shp <- input_list[[k]]
    input_name <- input_names[k]
    message("Processing input: ", input_name)

    # Content-aware cache key: fold a content/identity fingerprint of the
    # prediction input AND the burnable-domain identity into the filename. A
    # changed input on disk (e.g. keep_only.gpkg rewritten) OR a changed
    # burnable domain (the prediction is rasterized onto / masked by it) busts
    # the cache automatically, even with force_reprocess_pred = FALSE. Unchanged
    # inputs reuse the cache as before.
    pred_in_tag <- paste0(.vfm_vector_fingerprint(shp), "_", dom_cache_tag)
    pred_tif <- validation_output_path(
      cache_output_dir,
      paste0("predicted_fire_mask_", year_target, "_", input_name, "_", pred_in_tag, ".tif")
    )
    if (force_reprocess_pred && file.exists(pred_tif)) {
      message("Force reprocessing prediction raster: deleting ", basename(pred_tif))
      unlink(pred_tif)
    }

    det <- if (inherits(shp, "sf")) shp else sf::st_read(shp, quiet = TRUE)
    det <- det |> make_valid_sf()

    if (sf::st_crs(det) != sf::st_crs(mask_geom)) det <- sf::st_transform(det, sf::st_crs(mask_geom))

    det <- sf::st_filter(det, mask_geom, .predicate = sf::st_intersects)
    if (nrow(det) == 0) {
      warning("Skipping ", input_name, ": no polygons inside mask.")
      next
    }
    det <- suppressWarnings(sf::st_intersection(det, mask_geom)) |> make_valid_sf()

    det <- dissolve_by_field(det, dissolve_input_by) |> make_valid_sf()

    if (!file.exists(pred_tif)) {
      pred0 <- domain_mask
      pred0[] <- 0
      det_v <- safe_vect(det, label = input_name)

      pred_r <- terra::rasterize(det_v, pred0, field = 1, background = 0)
      pred_r[is.na(domain_mask)] <- NA

      terra::writeRaster(
        pred_r, pred_tif, overwrite = TRUE,
        datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255")
      )
    }

    pred_r <- terra::rast(pred_tif)

    # ---- PIXEL METRICS ----
    if (metrics_type %in% c("all", "pixel")) {

      ref_mask_use <- ref_mask_r
      if (buffer > 0) {
        ref_buf <- suppressWarnings(sf::st_buffer(ref_polygons, dist = buffer)) |> make_valid_sf()
        ref_buf_v <- safe_vect(ref_buf, label = "ref_buffer")
        ref_mask_use <- terra::rasterize(ref_buf_v, domain_mask, field = 1, background = 0)
        ref_mask_use[is.na(domain_mask)] <- NA
      }

      tp <- terra::global((pred_r == 1) & (ref_mask_use == 1), "sum", na.rm = TRUE)[[1]]
      fp <- terra::global((pred_r == 1) & (ref_mask_use == 0), "sum", na.rm = TRUE)[[1]]
      fn <- terra::global((pred_r == 0) & (ref_mask_use == 1), "sum", na.rm = TRUE)[[1]]
      tn <- terra::global((pred_r == 0) & (ref_mask_use == 0), "sum", na.rm = TRUE)[[1]]

      precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
      recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
      f1        <- ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                          2 * precision * recall / (precision + recall), NA_real_)
      iou       <- ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_)
      specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
      balanced_accuracy <- mean(c(recall, specificity), na.rm = TRUE)
      error_rate <- ifelse((tp + fp + fn + tn) > 0, (fp + fn) / (tp + fp + fn + tn), NA_real_)

      metrics_list[[input_name]] <- data.table::data.table(
        TP = tp, FP = fp, FN = fn, TN = tn,
        Precision = precision, Recall = recall, F1 = f1, IoU = iou,
        Specificity = specificity, BalancedAccuracy = balanced_accuracy,
        ErrorRate = error_rate,
        Reference_Observability_Filter_Applied = observability_applied,
        N_Reference_Polygons_Observable = n_reference_observable,
        N_Reference_Polygons_Excluded_Observability = n_reference_excluded_observability,
        InputName = input_name, Year = year_target
      )
    }

    # ---- AREA / POLYGON METRICS ----
    if (metrics_type %in% c("all", "area")) {

      ref_v <- safe_vect(ref_polygons, label = "ref_for_extract")
      det_v <- safe_vect(det, label = "det_for_extract")

      ref_pix <- terra::extract(ref_mask_r, ref_v, fun = sum, na.rm = TRUE)[, 2]
      ref_pix[is.na(ref_pix)] <- 0
      ref_area_i <- ref_pix * cell_area_ha

      int_pix <- terra::extract(pred_r, ref_v, fun = sum, na.rm = TRUE)[, 2]
      int_pix[is.na(int_pix)] <- 0
      int_area_i <- int_pix * cell_area_ha

      perc_detected_i <- ifelse(ref_area_i > 0, (int_area_i / ref_area_i) * 100, 0)
      detected_i <- perc_detected_i >= threshold_min_detected

      n_total <- length(ref_area_i)
      n_detected <- sum(detected_i)
      n_not_detected <- n_total - n_detected
      completely_detected <- sum(perc_detected_i >= threshold_completely_detected)

      if (length(perc_detected_i) > 0L) {
        cov_mean   <- round(mean(perc_detected_i, na.rm = TRUE), 2)
        cov_median <- round(stats::median(perc_detected_i, na.rm = TRUE), 2)
        cov_p10    <- round(stats::quantile(perc_detected_i, 0.10,
                                             na.rm = TRUE, names = FALSE), 2)
        cov_p90    <- round(stats::quantile(perc_detected_i, 0.90,
                                             na.rm = TRUE, names = FALSE), 2)
      } else {
        cov_mean <- NA_real_; cov_median <- NA_real_
        cov_p10  <- NA_real_; cov_p90    <- NA_real_
      }
      detected_definition <- sprintf("coverage_ref >= %d%%",
                                      as.integer(threshold_min_detected))

      area_reference_total <- sum(ref_area_i, na.rm = TRUE)
      det_total_pix <- terra::global(pred_r == 1, "sum", na.rm = TRUE)[[1]]
      area_detected_total <- det_total_pix * cell_area_ha
      inter_total_pix <- terra::global((pred_r == 1) & (ref_mask_r == 1), "sum", na.rm = TRUE)[[1]]
      area_intersection_total <- inter_total_pix * cell_area_ha

      recall_area    <- ifelse(area_reference_total > 0, (area_intersection_total / area_reference_total) * 100, NA_real_)
      precision_area <- ifelse(area_detected_total > 0, (area_intersection_total / area_detected_total) * 100, NA_real_)

      polygon_summary_list[[input_name]] <- data.table::data.table(
        InputName = input_name,
        Year = year_target,
        N_Reference_Polygons = n_total,
        N_Completely_Detected = completely_detected,
        N_Detected_Polygons = n_detected,
        N_Not_Detected = n_not_detected,
        Perc_Detected_Polygons = round((n_detected / n_total) * 100, 2),
        Area_Reference_ha = round(area_reference_total, 2),
        Area_Detected_ha = round(area_detected_total, 2),
        Area_Intersection_ha = round(area_intersection_total, 2),
        Recall_Area_percent = round(recall_area, 2),
        Precision_Area_percent = round(precision_area, 2),
        Coverage_mean = cov_mean,
        Coverage_median = cov_median,
        Coverage_p10 = cov_p10,
        Coverage_p90 = cov_p90,
        Reference_Observability_Filter_Applied = observability_applied,
        N_Reference_Polygons_Observable = n_reference_observable,
        N_Reference_Polygons_Excluded_Observability = n_reference_excluded_observability,
        Detected_Definition = detected_definition
      )

      ref_not_detected <- ref_polygons[!detected_i, , drop = FALSE]
      det_ref_pix <- terra::extract(ref_mask_r, det_v, fun = sum, na.rm = TRUE)[, 2]
      det_ref_pix[is.na(det_ref_pix)] <- 0
      det_not_matched <- det[det_ref_pix == 0, , drop = FALSE]

      if (nrow(ref_not_detected) > 0) {
        sf::st_write(
          ref_not_detected,
          validation_output_path(
            error_output_dir,
            paste0("reference_fires_omission_", year_target, "_", input_name, ".gpkg")
          ),
          delete_dsn = TRUE, quiet = TRUE
        )
      }
      if (nrow(det_not_matched) > 0) {
        sf::st_write(
          det_not_matched,
          validation_output_path(
            error_output_dir,
            paste0("predicted_fires_commission_", year_target, "_", input_name, ".gpkg")
          ),
          delete_dsn = TRUE, quiet = TRUE
        )
      }

      if (!is.null(class_shape) && file.exists(class_shape) && !is.null(class_field)) {
        class_map <- sf::st_read(class_shape, quiet = TRUE) |> make_valid_sf()
        if (sf::st_crs(class_map) != sf::st_crs(mask_geom)) {
          class_map <- sf::st_transform(class_map, sf::st_crs(mask_geom))
        }

        if (class_field %in% names(class_map)) {

          if (nrow(ref_not_detected) > 0) {
            ref_nd <- sf::st_join(ref_not_detected, class_map[, class_field, drop = FALSE], left = TRUE)
            # Block 10b: compute area_ha on the sf object directly (works
            # regardless of whether the geometry column is named
            # "geometry", "geom", or anything else); then drop geometry
            # before group_by so the summarise stays scalar.
            ref_nd$area_ha <- as.numeric(sf::st_area(ref_nd)) / 10000
            omission_by_class <- ref_nd |>
              sf::st_drop_geometry() |>
              dplyr::group_by(.data[[class_field]]) |>
              dplyr::summarise(
                N_Omitted = dplyr::n(),
                Area_Omitted_ha = sum(.data$area_ha, na.rm = TRUE),
                .groups = "drop"
              )

            data.table::fwrite(
              omission_by_class,
              validation_output_path(
                error_output_dir,
                paste0("reference_fires_omission_by_", class_field, "_", year_target, "_", input_name, ".csv")
              )
            )
          }

          if (nrow(det_not_matched) > 0) {
            det_nm <- sf::st_join(det_not_matched, class_map[, class_field, drop = FALSE], left = TRUE)
            det_nm$area_ha <- as.numeric(sf::st_area(det_nm)) / 10000
            commission_by_class <- det_nm |>
              sf::st_drop_geometry() |>
              dplyr::group_by(.data[[class_field]]) |>
              dplyr::summarise(
                N_Commission = dplyr::n(),
                Area_Commission_ha = sum(.data$area_ha, na.rm = TRUE),
                .groups = "drop"
              )

            data.table::fwrite(
              commission_by_class,
              validation_output_path(
                error_output_dir,
                paste0("predicted_fires_commission_by_", class_field, "_", year_target, "_", input_name, ".csv")
              )
            )
          }

        } else {
          warning("Field '", class_field, "' not found in class map. Skipping class-wise validation.")
        }
      }
    }

    # ---- PIXEL-LEVEL PER-STRATUM METRICS (Block 10c) ----
    if (!is.null(strata_r)) {
      st_out <- .of_validate_per_stratum_chunks(
        pred_r        = pred_r,
        ref_r         = ref_mask_r,
        strata_r      = strata_r,
        strata_lut_df = strata_lut_df,
        chunk_rows    = chunk_rows,
        na_as_zero    = isTRUE(na_strata_as_zero),
        cell_area_ha  = cell_area_ha,
        input_name    = input_name,
        year_target   = year_target
      )
      strata_table_list[[input_name]]  <- st_out$table
      strata_global_list[[input_name]] <- st_out$global
      strata_diag_list[[input_name]]   <- st_out$diagnostics

      data.table::fwrite(
        st_out$table,
        validation_output_path(
          strata_output_dir,
          paste0("pixel_metrics_by_stratum_", year_target, "_", input_name, ".csv")
        )
      )
      data.table::fwrite(
        st_out$global,
        validation_output_path(
          strata_output_dir,
          paste0("pixel_metrics_strata_global_", year_target, "_", input_name, ".csv")
        )
      )
      data.table::fwrite(
        st_out$diagnostics,
        validation_output_path(
          strata_output_dir,
          paste0("strata_diagnostics_", year_target, "_", input_name, ".csv")
        )
      )
    }
  }

  # ---- write outputs ----
  all_metrics <- NULL
  all_polygon_summary <- NULL

  if (metrics_type %in% c("all", "pixel")) {
    all_metrics <- data.table::rbindlist(metrics_list, fill = TRUE)
    data.table::fwrite(
      all_metrics,
      validation_output_path(
        summary_output_dir,
        paste0("pixel_metrics_global_", year_target, ".csv")
      )
    )
  }

  if (metrics_type %in% c("all", "area")) {
    all_polygon_summary <- data.table::rbindlist(polygon_summary_list, fill = TRUE)
    data.table::fwrite(
      all_polygon_summary,
      validation_output_path(
        summary_output_dir,
        paste0("fire_metrics_global_", year_target, ".csv")
      )
    )
  }

  # ---- combined Excel workbook (Block 10c, optional) ----
  excel_path <- NULL
  if (isTRUE(write_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("Package 'openxlsx' is required for write_excel = TRUE. ",
              "Skipping Excel workbook.", call. = FALSE)
    } else {
      excel_name <- if (!is.null(excel_filename) && nzchar(excel_filename))
        excel_filename else
          sprintf("validation_ALL_%s_res30.xlsx", year_target)
      excel_path <- validation_output_path(summary_output_dir, excel_name)
      wb <- openxlsx::createWorkbook()
      if (!is.null(all_metrics) && nrow(all_metrics) > 0L) {
        openxlsx::addWorksheet(wb, "pixel_global")
        openxlsx::writeDataTable(wb, "pixel_global", as.data.frame(all_metrics),
                                  withFilter = TRUE)
      }
      if (length(strata_table_list) > 0L) {
        openxlsx::addWorksheet(wb, "pixel_by_stratum")
        openxlsx::writeDataTable(
          wb, "pixel_by_stratum",
          as.data.frame(data.table::rbindlist(strata_table_list, fill = TRUE)),
          withFilter = TRUE
        )
      }
      if (length(strata_global_list) > 0L) {
        openxlsx::addWorksheet(wb, "stratum_global")
        openxlsx::writeDataTable(
          wb, "stratum_global",
          as.data.frame(data.table::rbindlist(strata_global_list, fill = TRUE)),
          withFilter = TRUE
        )
      }
      if (length(strata_diag_list) > 0L) {
        openxlsx::addWorksheet(wb, "diagnostics_strata")
        openxlsx::writeDataTable(
          wb, "diagnostics_strata",
          as.data.frame(data.table::rbindlist(strata_diag_list, fill = TRUE)),
          withFilter = TRUE
        )
      }
      if (!is.null(all_polygon_summary) && nrow(all_polygon_summary) > 0L) {
        openxlsx::addWorksheet(wb, "polygon_global")
        openxlsx::writeDataTable(wb, "polygon_global",
                                  as.data.frame(all_polygon_summary),
                                  withFilter = TRUE)
      }
      cache_tag <- sprintf(
        "y%s_res%d_dom%s_na0%s_%s",
        as.character(year_target),
        as.integer(round(terra::res(domain_mask)[1])),
        if (!is.null(strata_r)) "strata" else "burnable",
        isTRUE(na_strata_as_zero),
        obs_cache_tag
      )
      run_info <- data.frame(
        year_target           = year_target,
        validation_dir        = validation_dir,
        binary_burnable       = binary_burnable,
        burnable_threshold    = burnable_threshold,
        buffer                = buffer,
        threshold_completely_detected = threshold_completely_detected,
        threshold_min_detected = threshold_min_detected,
        min_area_reference_ha = if (is.null(min_area_reference_ha)) NA_real_ else min_area_reference_ha,
        observability_used    = observability_applied,
        ref_end_doy_col       = if (is.null(observability_raster)) NA_character_ else ref_end_doy_col,
        ref_start_doy_col     = if (is.null(observability_raster)) NA_character_ else ref_start_doy_col,
        metrics_type          = metrics_type,
        dissolve_ref_by       = if (is.null(dissolve_ref_by))   NA_character_ else dissolve_ref_by,
        dissolve_input_by     = if (is.null(dissolve_input_by)) NA_character_ else dissolve_input_by,
        strata_used           = !is.null(strata_r),
        chunk_rows            = chunk_rows,
        na_strata_as_zero     = isTRUE(na_strata_as_zero),
        cache_tag             = cache_tag,
        produced_at           = format(Sys.time()),
        stringsAsFactors      = FALSE
      )
      openxlsx::addWorksheet(wb, "run_info")
      openxlsx::writeDataTable(wb, "run_info", run_info)
      openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
    }
  }

  list(
    metrics            = if (metrics_type %in% c("all", "pixel")) all_metrics else NULL,
    polygon_summary    = if (metrics_type %in% c("all", "area")) all_polygon_summary else NULL,
    reference_observability = reference_observability,
    observability_audit = observability_audit,
    pixel_by_stratum   = if (length(strata_table_list))  data.table::rbindlist(strata_table_list,  fill = TRUE) else NULL,
    stratum_global     = if (length(strata_global_list)) data.table::rbindlist(strata_global_list, fill = TRUE) else NULL,
    diagnostics_strata = if (length(strata_diag_list))   data.table::rbindlist(strata_diag_list,   fill = TRUE) else NULL,
    excel_path         = excel_path
  )
}

# -----------------------------------------------------------------------------
# Block 10c: chunked per-stratum tabulator. Mirrors the legacy
# by_stratum_metrics() in 00_FUNCTIONS/05_VALIDATION_FUNCTIONS/
# FUNCTION_VALIDATE_FIRE_MAPS.R. Computes TP/FP/FN/TN per stratum over
# the strata-defined evaluation domain, plus the strata-domain global
# aggregate and a diagnostics row.
#
# Inputs:
#   pred_r         : SpatRaster (binary 0/1/NA) on the burnable domain.
#   ref_r          : SpatRaster (binary 0/1/NA) reference mask, aligned.
#   strata_r       : SpatRaster (integer ids; NA outside domain), aligned.
#   strata_lut_df  : data.frame(id, label) or NULL.
#   chunk_rows     : integer chunk size.
#   na_as_zero     : treat NA pred / ref inside strata domain as 0.
#   cell_area_ha   : numeric, hectares per pixel.
#   input_name     : character, used to tag rows.
#   year_target    : numeric/character, used to tag rows.
#
# Returns a list with three data.frames: table, global, diagnostics.
#
#' @keywords internal
#' @noRd
.of_validate_per_stratum_chunks <- function(pred_r, ref_r, strata_r,
                                             strata_lut_df,
                                             chunk_rows, na_as_zero,
                                             cell_area_ha,
                                             input_name, year_target) {
  if (!terra::compareGeom(pred_r, ref_r, stopOnError = FALSE))
    stop("pred and ref rasters are not aligned for per-stratum metrics.",
         call. = FALSE)
  if (!terra::compareGeom(pred_r, strata_r, stopOnError = FALSE))
    stop("strata raster is not aligned with pred/ref.", call. = FALSE)

  if (!is.null(strata_lut_df)) {
    ids <- sort(unique(strata_lut_df$id[is.finite(strata_lut_df$id)]))
  } else {
    f <- terra::freq(strata_r, value = NULL, useNA = "no")
    if (is.null(f) || nrow(f) == 0L)
      stop("Could not infer strata ids from strata_r and no strata_lut supplied.",
           call. = FALSE)
    ids <- sort(as.integer(f$value))
  }
  ids <- ids[ids >= 1L]
  if (length(ids) == 0L)
    stop("No valid strata ids (>= 1).", call. = FALSE)

  K <- length(ids)
  id_to_idx <- stats::setNames(seq_len(K), as.character(ids))
  counts_vec <- numeric(K * 4L)

  g_TN <- 0; g_FN <- 0; g_FP <- 0; g_TP <- 0
  n_eval <- 0L
  n_na_pred_in <- 0L
  n_na_ref_in  <- 0L

  nr <- terra::nrow(pred_r)
  starts <- seq(1L, nr, by = chunk_rows)

  terra::readStart(pred_r); terra::readStart(ref_r); terra::readStart(strata_r)
  on.exit({
    terra::readStop(pred_r); terra::readStop(ref_r); terra::readStop(strata_r)
  }, add = TRUE)

  for (r0 in starts) {
    nrs <- min(chunk_rows, nr - r0 + 1L)
    p   <- terra::readValues(pred_r,   row = r0, nrows = nrs, mat = FALSE)
    rr  <- terra::readValues(ref_r,    row = r0, nrows = nrs, mat = FALSE)
    sid <- terra::readValues(strata_r, row = r0, nrows = nrs, mat = FALSE)

    ok <- !is.na(sid)
    if (!any(ok)) next
    sid <- sid[ok]; p <- p[ok]; rr <- rr[ok]
    n_eval <- n_eval + length(sid)

    if (na_as_zero) {
      n_na_pred_in <- n_na_pred_in + sum(is.na(p))
      n_na_ref_in  <- n_na_ref_in  + sum(is.na(rr))
      p[is.na(p)]   <- 0
      rr[is.na(rr)] <- 0
    } else {
      keep2 <- !is.na(p) & !is.na(rr)
      if (!any(keep2)) next
      sid <- sid[keep2]; p <- p[keep2]; rr <- rr[keep2]
    }

    p01 <- as.integer(p == 1)
    r01 <- as.integer(rr == 1)
    state <- p01 * 2L + r01   # 0=TN, 1=FN, 2=FP, 3=TP

    tstate <- tabulate(state + 1L, nbins = 4L)
    g_TN <- g_TN + tstate[1]
    g_FN <- g_FN + tstate[2]
    g_FP <- g_FP + tstate[3]
    g_TP <- g_TP + tstate[4]

    sid_int <- as.integer(round(sid))
    keep <- sid_int %in% ids
    if (!any(keep)) next
    sid_int <- sid_int[keep]
    state   <- state[keep]

    idx <- unname(id_to_idx[as.character(sid_int)])
    key <- (idx - 1L) * 4L + state + 1L
    counts_vec <- counts_vec + tabulate(key, nbins = K * 4L)
  }

  mat <- matrix(counts_vec, nrow = K, ncol = 4L, byrow = TRUE)
  TN <- mat[, 1]; FN <- mat[, 2]; FP <- mat[, 3]; TP <- mat[, 4]

  Precision  <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  Recall     <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  F1         <- ifelse(!is.na(Precision) & !is.na(Recall) &
                         (Precision + Recall) > 0,
                       2 * Precision * Recall / (Precision + Recall),
                       NA_real_)
  IoU        <- ifelse((TP + FP + FN) > 0, TP / (TP + FP + FN), NA_real_)
  Specificity      <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  BalancedAccuracy <- mapply(function(r, s) mean(c(r, s), na.rm = TRUE),
                              Recall, Specificity)
  Commission <- ifelse(!is.na(Precision), 1 - Precision, NA_real_)
  Omission   <- ifelse(!is.na(Recall),    1 - Recall,    NA_real_)
  ErrorRate  <- ifelse((TP + FP + FN + TN) > 0,
                        (FP + FN) / (TP + FP + FN + TN),
                        NA_real_)

  Area_Reference_ha    <- (TP + FN) * cell_area_ha
  Area_Detected_ha     <- (TP + FP) * cell_area_ha
  Area_Intersection_ha <- TP * cell_area_ha
  Recall_Area_percent    <- ifelse(Area_Reference_ha > 0,
                                    100 * Area_Intersection_ha / Area_Reference_ha,
                                    NA_real_)
  Precision_Area_percent <- ifelse(Area_Detected_ha  > 0,
                                    100 * Area_Intersection_ha / Area_Detected_ha,
                                    NA_real_)

  table_df <- data.frame(
    InputName  = input_name,
    Year       = year_target,
    Stratum_ID = ids,
    TP = as.integer(TP), FP = as.integer(FP),
    FN = as.integer(FN), TN = as.integer(TN),
    Precision = Precision, Recall = Recall, F1 = F1, IoU = IoU,
    Specificity = Specificity, BalancedAccuracy = BalancedAccuracy,
    Commission = Commission, Omission = Omission,
    ErrorRate = ErrorRate,
    Area_Reference_ha = Area_Reference_ha,
    Area_Detected_ha = Area_Detected_ha,
    Area_Intersection_ha = Area_Intersection_ha,
    Recall_Area_percent = Recall_Area_percent,
    Precision_Area_percent = Precision_Area_percent,
    TemplateRes_m = terra::res(pred_r)[1],
    stringsAsFactors = FALSE
  )

  table_df$Stratum_Label <- as.character(table_df$Stratum_ID)
  if (!is.null(strata_lut_df)) {
    lut2 <- strata_lut_df[, c("id", "label"), drop = FALSE]
    table_df <- merge(table_df, lut2, by.x = "Stratum_ID", by.y = "id",
                       all.x = TRUE, sort = FALSE)
    table_df$Stratum_Label <- ifelse(is.na(table_df$label),
                                      table_df$Stratum_Label, table_df$label)
    table_df$label <- NULL
  }

  g_Precision  <- ifelse((g_TP + g_FP) > 0, g_TP / (g_TP + g_FP), NA_real_)
  g_Recall     <- ifelse((g_TP + g_FN) > 0, g_TP / (g_TP + g_FN), NA_real_)
  g_F1         <- ifelse(!is.na(g_Precision) & !is.na(g_Recall) &
                           (g_Precision + g_Recall) > 0,
                         2 * g_Precision * g_Recall / (g_Precision + g_Recall),
                         NA_real_)
  g_IoU        <- ifelse((g_TP + g_FP + g_FN) > 0,
                         g_TP / (g_TP + g_FP + g_FN), NA_real_)
  g_Specificity      <- ifelse((g_TN + g_FP) > 0, g_TN / (g_TN + g_FP), NA_real_)
  g_BalancedAccuracy <- mean(c(g_Recall, g_Specificity), na.rm = TRUE)
  g_Comm       <- ifelse(!is.na(g_Precision), 1 - g_Precision, NA_real_)
  g_Omi        <- ifelse(!is.na(g_Recall),    1 - g_Recall,    NA_real_)
  g_ErrorRate  <- ifelse((g_TP + g_FP + g_FN + g_TN) > 0,
                          (g_FP + g_FN) / (g_TP + g_FP + g_FN + g_TN),
                          NA_real_)
  g_area_ref   <- (g_TP + g_FN) * cell_area_ha
  g_area_det   <- (g_TP + g_FP) * cell_area_ha
  g_area_int   <- g_TP * cell_area_ha

  global_df <- data.frame(
    InputName = input_name,
    Year      = year_target,
    TP = as.integer(g_TP), FP = as.integer(g_FP),
    FN = as.integer(g_FN), TN = as.integer(g_TN),
    Precision = g_Precision, Recall = g_Recall, F1 = g_F1, IoU = g_IoU,
    Specificity = g_Specificity, BalancedAccuracy = g_BalancedAccuracy,
    Commission = g_Comm, Omission = g_Omi,
    ErrorRate = g_ErrorRate,
    Area_Reference_ha    = g_area_ref,
    Area_Detected_ha     = g_area_det,
    Area_Intersection_ha = g_area_int,
    Recall_Area_percent    = ifelse(g_area_ref > 0, 100 * g_area_int / g_area_ref, NA_real_),
    Precision_Area_percent = ifelse(g_area_det > 0, 100 * g_area_int / g_area_det, NA_real_),
    TemplateRes_m = terra::res(pred_r)[1],
    cell_area_ha  = cell_area_ha,
    stringsAsFactors = FALSE
  )

  diag_df <- data.frame(
    InputName             = input_name,
    Year                  = year_target,
    n_eval_pixels         = n_eval,
    na_pred_inside_domain = n_na_pred_in,
    na_ref_inside_domain  = n_na_ref_in,
    frac_na_pred_inside   = if (n_eval > 0) n_na_pred_in / n_eval else NA_real_,
    frac_na_ref_inside    = if (n_eval > 0) n_na_ref_in  / n_eval else NA_real_,
    chunk_rows            = chunk_rows,
    na_as_zero            = isTRUE(na_as_zero),
    stringsAsFactors      = FALSE
  )

  list(table = table_df, global = global_df, diagnostics = diag_df)
}

utils::globalVariables(c(".data", "geometry"))

