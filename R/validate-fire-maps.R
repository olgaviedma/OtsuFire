#' Validate burned-area maps against reference fire perimeters
#'
#' @description
#' Compare one or more burned-area polygon maps with reference fire perimeters
#' within a common evaluation domain.
#'
#' The function supports maps derived from Otsu-guided patches, probabilistic refinement
#' predictions, and external burned-area products. Inputs should contain the
#' polygons classified as burned. For probabilistic refinement outputs, apply the selected
#' score threshold before validation.
#'
#' The evaluation domain is defined by the study-area boundary and
#' burnable-area raster. Reference polygons are filtered to the target year
#' and prepared for comparison within this domain. Optional controls restrict
#' reference size, assess temporal observability, and dissolve polygons into
#' grouped features.
#'
#' Two complementary validation approaches are available:
#' * \strong{Pixel-based validation}: agreement between predicted and
#'   reference burned/unburned cells.
#' * \strong{Polygon- and area-based validation}: reference-fire detection and
#'   spatial coverage.
#'
#' Additional options provide polygon-level error summaries by class,
#' pixel-level metrics by raster strata, and a combined Excel workbook.
#'
#' @param input_shapefile Character vector of polygon-file paths, or an `sf`
#'   object. Burned-area predictions to validate. Each supplied map is
#'   evaluated independently. For a GeoPackage containing multiple layers,
#'   read the intended layer with `sf::st_read()` and pass the resulting
#'   object.
#' @param ref_shapefile Polygon-file path or `sf` object containing reference
#'   fire perimeters.
#' @param mask_shapefile Character scalar. Path to the study-area boundary
#'   polygon.
#' @param burnable_raster Character scalar. Path to the raster defining the
#'   burnable evaluation domain.
#' @param year_target Numeric scalar. Target year used to filter reference
#'   polygons.
#' @param validation_dir Character scalar. Root output directory. Validation
#'   products are written under a `VALIDATION` subdirectory.
#' @param binary_burnable Logical scalar. When `TRUE`, raster values at or
#'   above `burnable_threshold` are burnable. When `FALSE`, burnable values are
#'   selected using `burnable_classes`. Default: `TRUE`.
#' @param burnable_classes Numeric vector of burnable raster codes used when
#'   `binary_burnable = FALSE`. Default: `NULL`.
#' @param burnable_threshold Numeric scalar between `0` and `1`.
#'   Burnable-cell threshold used when `binary_burnable = TRUE`. Default:
#'   `0.5`.
#' @param metrics_type Character scalar. `"all"` computes pixel and
#'   polygon/area metrics; `"pixel"` or `"area"` selects one branch. Default:
#'   `"all"`.
#' @param buffer Numeric scalar. Buffer distance in metres applied around
#'   reference polygons for pixel-based comparison. Default: `0`.
#' @param threshold_completely_detected Numeric scalar between `0` and `100`.
#'   Minimum percentage coverage of a reference polygon required to count it
#'   as completely detected. Default: `90`.
#' @param threshold_min_detected Numeric scalar between `0` and `100`. Minimum
#'   percentage coverage of a reference polygon required to count it as
#'   detected. Default: `10`. Detection thresholds are percentages, whereas
#'   `observability_min_fraction` is a fraction between `0` and `1`.
#' @param min_area_reference_ha Optional numeric scalar. Minimum
#'   reference-polygon area in hectares, measured after restriction to the
#'   burnable domain. Default: `NULL`.
#' @param dissolve_ref_by Optional character scalar. Reference attribute used
#'   to dissolve polygons into grouped features. Default: `NULL`.
#' @param dissolve_input_by Optional character scalar. Prediction attribute
#'   used to dissolve polygons into grouped features. Default: `NULL`.
#' @param observability_raster `terra::SpatRaster`, raster path, or `NULL`.
#'   Day-of-year raster used to assess whether imagery observations reach the
#'   required reference-fire date. A single-layer raster is used directly. A
#'   multi-layer raster must contain a layer named `doy` or with a `_doy` or
#'   `doy_` name component; the first matching layer is selected. Default:
#'   `NULL`, disabling this assessment.
#' @param ref_end_doy_col Character scalar. Reference attribute containing
#'   fire end day of year. Used as the required date unless another supported
#'   date rule applies. Default: `"end_doy"`.
#' @param ref_start_doy_col Character scalar. Reference attribute containing
#'   fire start day of year. Used as a fallback only under `"legacy_any"` when
#'   the end date is unavailable. Default: `"start_doy"`.
#' @param ref_obs_doy_col Optional character scalar. Authoritative reference
#'   attribute containing the day from which a fire can reasonably be expected
#'   to appear in the map. Used only under `"wholefire_fraction"`. When
#'   `NULL`, the required date comes from `ref_end_doy_col`.
#' @param ref_obs_doy_fallback Character scalar. Treatment of missing values
#'   in a supplied `ref_obs_doy_col`: `"none"` retains an undetermined date;
#'   `"end_doy"` allows fallback to `ref_end_doy_col`. Neither option falls
#'   back to the start date in `"wholefire_fraction"` mode. Default: `"none"`.
#' @param observability_undetermined_policy Character scalar. Whether fires
#'   without a valid required date remain in the evaluation: `"evaluate"`
#'   retains and flags them; `"exclude"` removes their reference geometry and
#'   corresponding evaluation territory. Applies only under
#'   `"wholefire_fraction"`. Default: `"evaluate"`.
#' @param observability_mode Character scalar. `"legacy_any"` requires at
#'   least one finite observation DOY at or after the required fire date.
#'   `"wholefire_fraction"` requires a minimum fraction of the fire's
#'   burnable-domain cells to meet that condition. Default: `"legacy_any"`.
#' @param observability_min_fraction Numeric scalar between `0` and `1`.
#'   Minimum temporally observable fraction required under
#'   `"wholefire_fraction"`. Default: `0.75`. Ignored under `"legacy_any"`.
#' @param class_shape Optional character path to a polygon classification
#'   layer used for omission and commission summaries. Requires
#'   `class_field`.
#' @param class_field Optional character scalar. Attribute in `class_shape`
#'   identifying the reporting classes.
#' @param strata_raster `terra::SpatRaster`, raster path, or `NULL`.
#'   Categorical raster defining strata for pixel-level validation.
#' @param strata_lut Optional lookup table containing `id` and `label`
#'   columns. Accepts a `data.frame`, CSV path, or XLSX path.
#' @param chunk_rows Integer scalar. Number of raster rows processed per chunk
#'   by the stratified tabulator. Larger values require more memory. Default:
#'   `1024L`.
#' @param na_strata_as_zero Logical scalar. Within the eligible strata domain,
#'   treat missing prediction/reference values as unburned when `TRUE`, or
#'   omit those cells when `FALSE`. Default: `TRUE`.
#' @param force_reprocess_ref Logical scalar. Force rebuilding of cached
#'   reference products. Default: `FALSE`.
#' @param force_reprocess_pred Logical scalar. Force rebuilding of cached
#'   prediction products. Default: `FALSE`.
#' @param write_excel Logical scalar. Also export available validation tables
#'   to a combined Excel workbook. Requires \pkg{openxlsx}. Default: `FALSE`.
#' @param excel_filename Optional character scalar. Excel filename. Default:
#'   `validation_ALL_<year_target>_res30.xlsx`. The filename itself does not
#'   establish the evaluation resolution.
#'
#' @section Preparing the evaluation:
#' The function prepares a common spatial domain from the study-area boundary
#' and burnable raster.
#'
#' Reference preparation includes target-year filtering, spatial clipping,
#' restriction to burnable land, and any requested area, observability, or
#' dissolve operations. Predictions are prepared for comparison within the
#' same evaluation setup.
#'
#' When several prediction paths are supplied, each map is evaluated
#' independently against the prepared reference.
#'
#' Use consistent domain and reference settings when comparing maps. Changing
#' observability rules, minimum reference area, or grouping can change the
#' evaluated population as well as the resulting metrics.
#'
#' @section Pixel-based metrics:
#' Pixel validation compares rasterised predictions and references within the
#' evaluation domain.
#'
#' | Count | Meaning |
#' |---|---|
#' | `TP` | Burned in both prediction and reference. |
#' | `FP` | Burned in prediction but unburned in reference. |
#' | `FN` | Unburned in prediction but burned in reference. |
#' | `TN` | Unburned in both prediction and reference. |
#'
#' Derived metrics include `Precision`, `Recall`, `F1`, `IoU`,
#' `Specificity`, `BalancedAccuracy`, and `ErrorRate`.
#'
#' These measure agreement with the supplied reference. Where reference
#' coverage is incomplete, apparent commission may include real burns absent
#' from the reference.
#'
#' A non-zero `buffer` changes the reference geometry used for pixel
#' comparison and therefore changes the interpretation of the resulting
#' agreement metrics.
#'
#' @section Polygon- and area-based metrics:
#' Reference-fire detection is based on the percentage of each prepared
#' reference polygon covered by predictions.
#'
#' A reference polygon counts as:
#' * \strong{detected} when coverage meets `threshold_min_detected`;
#' * \strong{completely detected} when coverage meets
#'   `threshold_completely_detected`.
#'
#' Use a minimum-detection threshold no greater than the complete-detection
#' threshold.
#'
#' With the defaults, a fire needs at least 10% coverage to count as detected
#' and at least 90% coverage to count as completely detected.
#'
#' Main fields include:
#'
#' | Field | Meaning |
#' |---|---|
#' | `N_Reference_Polygons` | Number of reference polygons evaluated. |
#' | `N_Detected_Polygons` | Number meeting the minimum detection threshold. |
#' | `N_Completely_Detected` | Number meeting the complete-detection threshold. |
#' | `N_Not_Detected` | Number failing the minimum detection threshold. |
#' | `Area_Reference_ha` | Evaluated reference burned area. |
#' | `Area_Detected_ha` | Predicted burned area included in the area comparison. |
#' | `Area_Intersection_ha` | Area shared by predictions and reference. |
#' | `Recall_Area_percent` | Percentage of reference burned area covered by predictions. |
#' | `Precision_Area_percent` | Percentage of predicted burned area overlapping the reference. |
#' | `Detected_Definition` | Recorded rule used to count detected reference polygons. |
#'
#' The documented rule uses `coverage_ref >= threshold_min_detected`. A
#' threshold of zero should therefore not be interpreted as requiring positive
#' overlap.
#'
#' Polygon counts depend on the reference geometry and any dissolve settings.
#' They represent fire-event counts only when the prepared polygons
#' correspond to individual events.
#'
#' @section Temporal observability:
#' Temporal observability assesses whether the observation dates extend far
#' enough to evaluate a reference fire.
#'
#' Both modes make a whole-fire decision. A retained fire is evaluated
#' throughout its burnable-domain geometry; partially observed portions are
#' not individually removed.
#'
#' \strong{Legacy rule.} With `observability_mode = "legacy_any"`:
#' 1. obtain the required fire date from the end DOY, falling back to the
#'    start DOY;
#' 2. find the maximum finite observation DOY within the fire;
#' 3. retain the fire when that observation reaches or exceeds the required
#'    date.
#'
#' A single qualifying cell is sufficient. Fires lacking a usable reference
#' date or finite observation data are excluded from the reference set.
#'
#' The diagnostics include the maximum observation date, required reference
#' date, and their difference.
#'
#' \strong{Whole-fire fraction rule.} With
#' `observability_mode = "wholefire_fraction"`, the temporally observable
#' fraction is:
#' \preformatted{
#' Number of burnable-domain cells with observation DOY >= required DOY
#' --------------------------------------------------------------------
#'              Total burnable-domain cells in the fire
#' }
#'
#' The fire passes when this fraction meets `observability_min_fraction`.
#'
#' The required date comes from `ref_obs_doy_col` when supplied, or otherwise
#' from `ref_end_doy_col`. Missing values in the authoritative column use the
#' end date only when `ref_obs_doy_fallback = "end_doy"` is explicitly
#' selected.
#'
#' There is no start-date fallback in this mode.
#'
#' Choose an observability raster whose DOY values represent the observation
#' timing relevant to the mapped product. This assessment does not establish
#' that an observation is cloud-free or otherwise suitable unless that
#' information is already represented in the supplied raster.
#'
#' @section Observability status and evaluation membership:
#' Under `"wholefire_fraction"`, observability evidence and inclusion in the
#' evaluation are recorded separately.
#'
#' | Status | Meaning | Evaluation treatment |
#' |---|---|---|
#' | `OUTSIDE_BURNABLE_DOMAIN` | The reference fire contains no burnable-domain cells. | Contributes no evaluated area. |
#' | `UNDETERMINED_OBS_DATE` | Burnable-domain cells exist, but no valid required reference date is available. | Retained with `"evaluate"`; excluded with `"exclude"`. |
#' | `NO_OBSERVABILITY_DATA` | A required date exists, but no usable observation data are available within the fire. | Excluded. |
#' | `NOT_OBSERVABLE` | Observation data exist, but the qualifying fraction is below the required minimum. | Excluded. |
#' | `OBSERVABLE` | The qualifying fraction meets the minimum. | Retained. |
#'
#' `observable_flag` records whether the status is `OBSERVABLE`.
#'
#' `in_evaluation_domain` records whether the fire is included in validation.
#' With the default `"evaluate"` policy, this includes fires whose status is
#' `UNDETERMINED_OBS_DATE`.
#'
#' A missing authoritative date is not evidence that the fire was temporally
#' unobservable. Such fires remain positive reference features under the
#' default policy and are flagged for audit.
#'
#' When a fire is excluded under `"wholefire_fraction"`, its whole geometry is
#' removed from the evaluation domain. Predictions within that excluded
#' territory are not counted as commission, and the territory contributes no
#' confusion-matrix counts.
#'
#' The supplied source files are not modified by this masking.
#'
#' Additional diagnostics include `undetermined_cause` and
#' `observable_flag_legacy_any`, allowing users to inspect missing-date causes
#' and compare observability rules.
#'
#' @section Partial-observability audit:
#' The function reports partial observability within reference fires,
#' including cell counts, observable fractions, and affected areas.
#'
#' This audit is descriptive. It does not introduce pixel-level trimming
#' within retained fires or independently change the metrics.
#'
#' The returned `observability_audit` contains:
#' * `summary`: overall counts and areas;
#' * `excluded`: excluded-fire records;
#' * `partial`: records for partially observable fires.
#'
#' Use these diagnostics to assess how strongly the whole-fire treatment
#' affects the evaluated reference population.
#'
#' @section Polygon classes and raster strata:
#' `class_shape` and `class_field` provide polygon-level omission and
#' commission summaries by an external classification, such as vegetation
#' type or administrative region.
#'
#' `strata_raster` provides pixel-level confusion matrices and metrics within
#' categorical raster strata.
#'
#' These are different reporting units and should not be interpreted
#' interchangeably.
#'
#' Cells without a stratum identifier are omitted from stratified
#' aggregation. Within eligible strata, `na_strata_as_zero` controls the
#' treatment of missing prediction/reference values.
#'
#' The aggregate stratified domain can differ from the global domain where
#' strata are missing. Consult `diagnostics_strata` when comparing global and
#' stratified results.
#'
#' @section Caching and reprocessing:
#' Prepared spatial products are cached under
#' `<validation_dir>/VALIDATION/_CACHE/`.
#'
#' Cache identities incorporate relevant input and processing settings,
#' including the burnable domain, reference preparation, and observability
#' configuration.
#'
#' File fingerprints use file size and modification time, supplemented by an
#' MD5 hash for files below 64 MB. In-memory inputs use structural
#' fingerprints. Consequently, cache identity should not be interpreted as a
#' complete byte-level comparison of every large input.
#'
#' Use `force_reprocess_ref = TRUE` or `force_reprocess_pred = TRUE` when a
#' rebuild is required, particularly after replacing inputs in place.
#'
#' The reference `buffer` is applied downstream of cached reference
#' preparation.
#'
#' @section Outputs:
#' Validation products are stored under `<validation_dir>/VALIDATION/`.
#'
#' Main global tables include:
#'
#' | File | Contents |
#' |---|---|
#' | `metrics_summary_<year>.csv` | Global pixel-based metrics. |
#' | `polygon_summary_<year>.csv` | Polygon- and area-based metrics. |
#'
#' Optional stratified tables include:
#' * `pixel_by_stratum_<year>_<input>.csv`;
#' * `stratum_global_<year>_<input>.csv`;
#' * `diagnostics_strata_<year>_<input>.csv`.
#'
#' Observability outputs include per-fire diagnostics, excluded-fire
#' geometries, and annual audit tables. Summary and audit products include:
#' * `OBSERVABILITY_SUMMARY_<year>.csv`;
#' * `OBSERVABILITY_UNDETERMINED_<year>.csv`;
#' * `OBSERVABILITY_AUDIT_<year>.csv`;
#' * `OBSERVABILITY_EXCLUDED_FIRES_<year>.csv`;
#' * `OBSERVABILITY_PARTIAL_FIRES_<year>.csv`.
#'
#' Additional vector outputs describe spatial disagreement. When requested,
#' an Excel workbook gathers the available validation tables; its location is
#' returned in `excel_path`.
#'
#' @return A named list containing the available outputs.
#'
#' | Field | Contents |
#' |---|---|
#' | `metrics` | Global pixel-based metrics, or `NULL` when `metrics_type = "area"`. |
#' | `polygon_summary` | Global polygon- and area-based metrics, or `NULL` when `metrics_type = "pixel"`. |
#' | `reference_observability` | Per-reference-fire observability diagnostics, including cell counts and observable fractions. `NULL` when no observability raster is supplied. |
#' | `observability_audit` | Descriptive audit containing `summary`, `excluded`, and `partial`. `NULL` when no observability raster is supplied. |
#' | `observability_summary` | Annual counts and areas by observability status, with evaluated area. |
#' | `observability_settings` | Resolved observability settings and domain information, including mode, fraction threshold, required-date column, evaluated and removed areas, and cache tag. |
#' | `pixel_by_stratum` | Per-stratum pixel metrics, or `NULL` when raster stratification is not requested. |
#' | `stratum_global` | Aggregate metrics over the stratified domain, or `NULL` when raster stratification is not requested. |
#' | `diagnostics_strata` | Stratified-domain diagnostics, or `NULL` when raster stratification is not requested. |
#' | `excel_path` | Excel workbook path, or `NULL` when `write_excel = FALSE`. |
#'
#' @seealso [write_thresholded_burned()], [score_supervised_burned_map()],
#'   [run_deterministic_pipeline()], [run_oneyear_supervised_pipeline()],
#'   [validate_supervised_execution()].
#'
#'   [validate_supervised_execution()] checks configuration and inputs before
#'   probabilistic refinement processing. `validate_fire_maps()` evaluates produced maps
#'   against reference perimeters.
#'
#' @examples
#' \dontrun{
#' # Read the specific thresholded prediction layer
#' predicted <- sf::st_read(
#'   "results/thresholded_burned.gpkg",
#'   layer = "thresholded_burned",
#'   quiet = TRUE
#' )
#'
#' # Global pixel and polygon/area validation
#' validation <- validate_fire_maps(
#'   input_shapefile = predicted,
#'   ref_shapefile = "data/reference_fires_2022.gpkg",
#'   mask_shapefile = "data/study_area.gpkg",
#'   burnable_raster = "data/burnable_mask.tif",
#'   year_target = 2022,
#'   validation_dir = "results/validation",
#'   metrics_type = "all",
#'   threshold_min_detected = 10,
#'   threshold_completely_detected = 90
#' )
#'
#' validation$metrics
#' validation$polygon_summary
#'
#' # Whole-fire observability assessment
#' # The reference must contain the specified obs_required_doy field.
#' validation_obs <- validate_fire_maps(
#'   input_shapefile = predicted,
#'   ref_shapefile = "data/reference_fires_2022.gpkg",
#'   mask_shapefile = "data/study_area.gpkg",
#'   burnable_raster = "data/burnable_mask.tif",
#'   year_target = 2022,
#'   validation_dir = "results/validation_observability",
#'   observability_raster = "data/observation_doy_2022.tif",
#'   observability_mode = "wholefire_fraction",
#'   observability_min_fraction = 0.75,
#'   ref_obs_doy_col = "obs_required_doy",
#'   ref_obs_doy_fallback = "none",
#'   observability_undetermined_policy = "evaluate"
#' )
#'
#' validation_obs$observability_summary
#' validation_obs$observability_settings
#'
#' table(
#'   validation_obs$reference_observability$observability_status,
#'   useNA = "ifany"
#' )
#'
#' # Pixel-level land-cover stratification and Excel export
#' validation_strata <- validate_fire_maps(
#'   input_shapefile = predicted,
#'   ref_shapefile = "data/reference_fires_2022.gpkg",
#'   mask_shapefile = "data/study_area.gpkg",
#'   burnable_raster = "data/burnable_mask.tif",
#'   year_target = 2022,
#'   validation_dir = "results/validation_strata",
#'   metrics_type = "all",
#'   strata_raster = "data/land_cover_strata.tif",
#'   strata_lut = "data/strata_labels.csv",
#'   chunk_rows = 1024L,
#'   na_strata_as_zero = TRUE,
#'   write_excel = TRUE,
#'   excel_filename = "validation_2022.xlsx"
#' )
#'
#' validation_strata$pixel_by_stratum
#' validation_strata$diagnostics_strata
#' validation_strata$excel_path
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
                               ref_obs_doy_col = NULL,
                               ref_obs_doy_fallback = c("none", "end_doy"),
                               observability_undetermined_policy = c("evaluate", "exclude"),
                               observability_mode = c("legacy_any", "wholefire_fraction"),
                               observability_min_fraction = 0.75,
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
  observability_mode <- match.arg(observability_mode)
  ref_obs_doy_fallback <- match.arg(ref_obs_doy_fallback)
  observability_undetermined_policy <- match.arg(observability_undetermined_policy)

  # §N+25 (2026-08-19). Whole-fire observability v2.
  #
  # `observability_mode` selects between two whole-fire rules. The unit of
  # decision is the FIRE in both cases; no pixel-level trimming of the
  # reference is introduced by either mode.
  #
  #   "legacy_any"          historical behaviour, preserved bit-for-bit:
  #                         observable <=> max(DOY inside the fire) >= required
  #                         DOY, with the required DOY taken from
  #                         `ref_end_doy_col` and falling back silently to
  #                         `ref_start_doy_col`. Mathematically identical to
  #                         "at least one observable pixel" (verified on
  #                         21,452 reference fires: 0 discrepancies). Excluded
  #                         fires are dropped from the reference but their
  #                         territory stays in the evaluation domain as
  #                         reference = 0, so a detection there scores as FP.
  #
  #   "wholefire_fraction"  observable <=> temporal_observable_fraction >=
  #                         `observability_min_fraction` (default 0.75), where
  #                         the fraction is (# burnable-domain pixels of the
  #                         fire whose composite DOY >= the fire's required
  #                         DOY) / (# burnable-domain pixels of the fire).
  #                         The required DOY is `ref_obs_doy_col` when finite
  #                         (Neves image DOY) and `ref_end_doy_col` otherwise;
  #                         there is NO silent fallback to the start DOY -- a
  #                         fire with no usable end date is flagged
  #                         UNDETERMINED_OBS_DATE. Fires that are not
  #                         observable, undetermined, or without observability
  #                         data have their WHOLE geometry removed from the
  #                         evaluation domain (reference AND prediction AND
  #                         strata are set to NA there) so that they can
  #                         neither create FP nor FN. Everything outside those
  #                         geometries keeps being evaluated normally, so
  #                         commission outside reference fires is unaffected.
  #
  # The mode never rewrites the prediction or the reference on disk: the
  # exclusion happens only on the in-memory rasters used for metrics.
  if (!is.numeric(observability_min_fraction) ||
      length(observability_min_fraction) != 1L ||
      !is.finite(observability_min_fraction) ||
      observability_min_fraction < 0 || observability_min_fraction > 1) {
    stop("observability_min_fraction must be a single number in [0, 1].",
         call. = FALSE)
  }
  if (!is.null(ref_obs_doy_col) &&
      (!is.character(ref_obs_doy_col) || length(ref_obs_doy_col) != 1L)) {
    stop("ref_obs_doy_col must be NULL or a single column name.", call. = FALSE)
  }
  # Only the fraction rule removes territory from the evaluation domain, so
  # only it makes the metric rasters depend on the observability settings.
  observability_masks_domain <-
    identical(observability_mode, "wholefire_fraction") &&
    !is.null(observability_raster)

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
                                            ref_start_doy_col,
                                            ref_obs_doy_col = NULL,
                                            ref_obs_doy_fallback = "none",
                                            observability_undetermined_policy = "evaluate",
                                            observability_mode = "legacy_any",
                                            observability_min_fraction = 0.75) {
    has_end_doy <- ref_end_doy_col %in% names(ref_polygons)
    has_start_doy <- ref_start_doy_col %in% names(ref_polygons)
    has_obs_doy <- !is.null(ref_obs_doy_col) &&
      ref_obs_doy_col %in% names(ref_polygons)
    if (!has_end_doy && !has_start_doy && !has_obs_doy) {
      stop(
        "Reference DOY columns not found. Expected at least one of: ",
        ref_end_doy_col, ", ", ref_start_doy_col,
        if (!is.null(ref_obs_doy_col)) paste0(", ", ref_obs_doy_col) else "",
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

    # Per-fire DATA-COVERAGE counts. `obs_r` is already aligned to and NA
    # outside the burnable domain, so a finite obs value marks a domain cell
    # that carries an observation.
    # total_pixels      = burnable-domain cells under the fire (= ref_area_pix);
    # observable_pixels = those that ALSO carry a finite observation DOY.
    # NOTE ON NAMING (kept for backward compatibility): the historical column
    # `observable_fraction` measures DATA COVERAGE, not temporal adequacy. It
    # is aliased below to the unambiguous `data_coverage_fraction`; the new
    # `temporal_observable_fraction` is a different quantity.
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
    alt_obs_doy <- if (has_obs_doy) {
      suppressWarnings(as.numeric(ref_polygons[[ref_obs_doy_col]]))
    } else {
      rep(NA_real_, nrow(ref_polygons))
    }

    if (identical(observability_mode, "legacy_any")) {
      # Historical cascade: end -> start. Unchanged, including the silent
      # start-DOY fallback, so that the legacy mode reproduces past runs.
      fire_required_doy <- ifelse(is.finite(end_doy), end_doy, start_doy)
      required_doy_source <- ifelse(
        is.finite(end_doy), ref_end_doy_col,
        ifelse(is.finite(start_doy), ref_start_doy_col, NA_character_)
      )
    } else if (has_obs_doy && identical(ref_obs_doy_fallback, "none")) {
      # v2 cascade, AUTHORITATIVE column (the default whenever
      # `ref_obs_doy_col` is supplied). The reference is expected to have
      # resolved the source-specific semantics upstream and to expose a single
      # column with the date from which it is reasonable to require the fire to
      # be mapped. Therefore an NA in that column is taken at face value: it
      # means UNDETERMINED, not "look somewhere else". No other date column is
      # consulted, so a deliberate NA can never be masked by a stale
      # `end_doy`. This is what keeps the methodological decision in the
      # reference instead of in the validator.
      fire_required_doy <- alt_obs_doy
      required_doy_source <- ifelse(is.finite(alt_obs_doy), ref_obs_doy_col,
                                    NA_character_)
    } else {
      # v2 cascade with an EXPLICITLY requested fallback
      # (`ref_obs_doy_fallback = "end_doy"`), or with no observability column at
      # all. Used for reference layers whose observability column is populated
      # only for some sources -- e.g. Reference v2's `DOY`, which exists for the
      # Neves atlas and is NA for EFFIS. NO start-DOY fallback in either case:
      # an unusable end date is reported as UNDETERMINED_OBS_DATE instead of
      # being imputed.
      fire_required_doy <- ifelse(is.finite(alt_obs_doy), alt_obs_doy, end_doy)
      required_doy_source <- ifelse(
        is.finite(alt_obs_doy),
        if (has_obs_doy) ref_obs_doy_col else NA_character_,
        ifelse(is.finite(end_doy), ref_end_doy_col, NA_character_)
      )
    }
    fire_required_doy[!is.finite(fire_required_doy)] <- NA_real_
    required_doy_source[is.na(fire_required_doy)] <- NA_character_

    # Why a fire ended up without a required DOY. Diagnostic only: it never
    # feeds a decision, but it is what tells the reference maintainers whether
    # an UNDETERMINED came from a deliberate NA in the authoritative column or
    # from a missing end date.
    # Each value names exactly what was missing, so none of them says
    # "end_date" when an authoritative observability column was in charge.
    undetermined_cause <- rep(NA_character_, nrow(ref_polygons))
    if (identical(observability_mode, "legacy_any")) {
      # only mode that consults the start DOY at all
      undetermined_cause[is.na(fire_required_doy)] <- "no_end_or_start_doy"
    } else if (has_obs_doy && identical(ref_obs_doy_fallback, "none")) {
      # the reference deliberately has no observability date for this fire
      undetermined_cause[is.na(fire_required_doy)] <-
        sprintf("na_in_authoritative_%s", ref_obs_doy_col)
    } else if (has_obs_doy) {
      # explicit fallback requested and neither column carried a date
      undetermined_cause[is.na(fire_required_doy)] <- "no_obs_doy_and_no_end_doy"
    } else {
      # no observability column supplied at all: the end DOY was the only
      # candidate, and it was missing
      undetermined_cause[is.na(fire_required_doy)] <- "no_end_doy"
    }

    obs_doy_margin <- obs_doy_max - fire_required_doy

    # Temporal observability at pixel level, aggregated back to ONE whole-fire
    # number. The long-form extract gives every burnable-domain cell of every
    # fire exactly once, so the per-fire counts are exact and no rasterization
    # of the required DOY (which would be ambiguous under overlapping fires) is
    # needed. This never trims the reference: the fraction only feeds the
    # whole-fire decision below.
    n_ref <- nrow(ref_polygons)
    obs_long <- terra::extract(obs_r, ref_v)
    long_id <- obs_long[[1]]
    long_val <- obs_long[[2]]
    long_ok <- is.finite(long_val) &
      is.finite(fire_required_doy[long_id]) &
      long_val >= fire_required_doy[long_id]
    n_temporally_observable <- tabulate(long_id[long_ok], nbins = n_ref)
    n_temporally_observable <- pmin(n_temporally_observable, total_pixels)
    temporal_observable_fraction <- ifelse(
      total_pixels > 0, n_temporally_observable / total_pixels, NA_real_
    )

    # A fire with no required DOY has no computable temporal fraction. The
    # predicate above yields 0 for it purely because `is.finite(NA)` is FALSE,
    # and a 0 there would read as "observed nowhere", i.e. as EVIDENCE of
    # non-observability, which is exactly what an undetermined date is not.
    # Report NA instead, so the absence of a date can never be mistaken for the
    # presence of a negative result. `data_coverage_fraction` stays valid: it
    # does not depend on the fire date.
    no_required_doy <- is.na(fire_required_doy)
    n_temporally_observable[no_required_doy] <- NA_integer_
    temporal_observable_fraction[no_required_doy] <- NA_real_

    data_coverage_fraction <- observable_fraction
    no_data_fraction <- ifelse(total_pixels > 0,
                               1 - data_coverage_fraction, NA_real_)
    post_fire_observable_fraction <- temporal_observable_fraction
    pre_fire_or_too_early_fraction <- ifelse(
      total_pixels > 0,
      pmax(data_coverage_fraction - temporal_observable_fraction, 0),
      NA_real_
    )

    if (identical(observability_mode, "legacy_any")) {
      observability_status <- ifelse(
        is.na(fire_required_doy), "UNDETERMINED_OBS_DATE",
        ifelse(
          is.na(obs_doy_max), "NO_OBSERVABILITY_DATA",
          ifelse(obs_doy_margin < 0, "NOT_OBSERVABLE", "OBSERVABLE")
        )
      )
      # Legacy reason strings are preserved verbatim.
      observable_reason <- ifelse(
        is.na(fire_required_doy), "missing_reference_doy",
        ifelse(
          is.na(obs_doy_max), "no_observability_data",
          ifelse(obs_doy_margin < 0, "obs_before_reference_doy", "observable")
        )
      )
      observable_flag <- observable_reason == "observable"
      # legacy never separated evidence from domain membership
      in_evaluation_domain <- observable_flag
    } else {
      # Status cascade, in precedence order. Each state answers a DIFFERENT
      # question, and only NOT_OBSERVABLE is evidence of temporal
      # non-observability:
      #   OUTSIDE_BURNABLE_DOMAIN  the fire has no burnable-domain cell at all.
      #                            A domain fact, not an observability one; it
      #                            contributes nothing whatever the dates say.
      #   UNDETERMINED_OBS_DATE    there IS evaluable surface, but no valid
      #                            authoritative date, so observability cannot
      #                            be verified. NOT evidence of non-observability.
      #   NO_OBSERVABILITY_DATA    there IS evaluable surface and a date, but the
      #                            composite carries no observation anywhere in
      #                            the fire. This IS evidence: nothing was seen.
      #   NOT_OBSERVABLE           date and observations exist, and the measured
      #                            fraction falls below the threshold. Evidence.
      #   OBSERVABLE               fraction reaches the threshold.
      observability_status <- ifelse(
        total_pixels <= 0, "OUTSIDE_BURNABLE_DOMAIN",
        ifelse(
          is.na(fire_required_doy), "UNDETERMINED_OBS_DATE",
          ifelse(
            is.na(obs_doy_max), "NO_OBSERVABILITY_DATA",
            ifelse(temporal_observable_fraction >= observability_min_fraction,
                   "OBSERVABLE", "NOT_OBSERVABLE")
          )
        )
      )
      observable_reason <- tolower(observability_status)
      # EVIDENCE only. This flag no longer decides who is evaluated.
      observable_flag <- observability_status == "OBSERVABLE"
      # DOMAIN MEMBERSHIP. Missing authoritative metadata is not evidence of
      # non-observability, so under the default "evaluate" policy an
      # UNDETERMINED_OBS_DATE fire stays a positive reference feature inside the
      # evaluation domain and contributes TP/FN like any other. "exclude" is
      # kept only so the two can be compared in a sensitivity run.
      in_evaluation_domain <- observability_status == "OBSERVABLE" |
        (identical(observability_undetermined_policy, "evaluate") &
           observability_status == "UNDETERMINED_OBS_DATE")
    }

    # Legacy rule recomputed alongside the active one, always, so that every
    # run carries the diagnostic needed to compare the two rules without a
    # second pass over the rasters.
    legacy_required_doy <- ifelse(is.finite(end_doy), end_doy, start_doy)
    legacy_required_doy[!is.finite(legacy_required_doy)] <- NA_real_
    legacy_margin <- obs_doy_max - legacy_required_doy
    observable_flag_legacy_any <- !is.na(legacy_required_doy) &
      !is.na(obs_doy_max) & legacy_margin >= 0

    out <- data.table::as.data.table(sf::st_drop_geometry(ref_polygons))
    out[, reference_row := seq_len(.N)]
    out[, ref_area_domain_ha := round(ref_area_domain_ha, 4)]
    out[, fire_required_doy := fire_required_doy]
    out[, required_doy_source := required_doy_source]
    # Historical aliases (identical values, kept so downstream readers and the
    # partial-observability audit keep working unchanged).
    out[, obs_ref_doy := fire_required_doy]
    out[, obs_ref_doy_source := required_doy_source]
    out[, obs_doy_max := obs_doy_max]
    out[, obs_doy_margin := obs_doy_margin]
    out[, observability_mode := observability_mode]
    out[, observability_min_fraction := observability_min_fraction]
    out[, observability_undetermined_policy := observability_undetermined_policy]
    out[, observability_status := observability_status]
    out[, undetermined_cause := undetermined_cause]
    out[, observable_flag := observable_flag]
    out[, in_evaluation_domain := in_evaluation_domain]
    out[, observable_reason := observable_reason]
    out[, observable_flag_legacy_any := observable_flag_legacy_any]
    out[, total_pixels := as.integer(round(total_pixels))]
    out[, observable_pixels := as.integer(round(observable_pixels))]
    out[, non_observable_pixels := as.integer(round(non_observable_pixels))]
    out[, observable_fraction := round(observable_fraction, 6)]
    out[, data_coverage_fraction := round(data_coverage_fraction, 6)]
    out[, n_temporally_observable := as.integer(round(n_temporally_observable))]
    out[, temporal_observable_fraction := round(temporal_observable_fraction, 6)]
    out[, no_data_fraction := round(no_data_fraction, 6)]
    out[, pre_fire_or_too_early_fraction := round(pre_fire_or_too_early_fraction, 6)]
    out[, post_fire_observable_fraction := round(post_fire_observable_fraction, 6)]
    out
  }

  warn_reference_observability <- function(reference_observability) {
    if (is.null(reference_observability) || !nrow(reference_observability)) return(invisible(NULL))
    # Report on DOMAIN membership, not on the evidence flag: a fire that is
    # evaluated despite an undetermined date has observable_flag = FALSE and
    # must NOT be announced as excluded.
    dropped <- reference_observability[in_evaluation_domain == FALSE]
    kept_undetermined <- if ("observability_status" %in% names(reference_observability)) {
      sum(reference_observability$in_evaluation_domain &
            reference_observability$observability_status == "UNDETERMINED_OBS_DATE")
    } else 0L
    if (nrow(dropped)) {
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
    }
    if (kept_undetermined > 0L) {
      warning(
        sprintf(
          paste0("%d reference fires have no valid authoritative observability ",
                 "date (UNDETERMINED_OBS_DATE) and are being EVALUATED as ",
                 "positive reference features: missing metadata is not evidence ",
                 "of temporal non-observability. They are flagged for audit in ",
                 "02_OBSERVABILITY/. Use observability_undetermined_policy = ",
                 "\"exclude\" for the sensitivity run."),
          kept_undetermined
        ),
        call. = FALSE
      )
    }
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

    # DOMAIN membership, not evidence: under the "evaluate" policy an
    # UNDETERMINED_OBS_DATE fire has observable_flag = FALSE but stays in.
    dropped_idx <- which(!reference_observability$in_evaluation_domain)
    if (!length(dropped_idx)) return(NULL)

    out <- ref_polygons[dropped_idx, , drop = FALSE]
    obs_cols <- c(
      "obs_ref_doy",
      "obs_ref_doy_source",
      "obs_doy_max",
      "obs_doy_margin",
      "observable_flag",
      "in_evaluation_domain",
      "observable_reason",
      # v2 columns: present since refschema-3, guarded so the helper stays
      # usable with any observability table.
      "fire_required_doy",
      "required_doy_source",
      "observability_status",
      "temporal_observable_fraction",
      "data_coverage_fraction"
    )
    for (col_nm in obs_cols) {
      if (!col_nm %in% names(reference_observability)) next
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
    ref_start_doy_col = ref_start_doy_col,
    observability_mode = observability_mode,
    observability_min_fraction = observability_min_fraction,
    ref_obs_doy_col = ref_obs_doy_col,
    ref_obs_doy_fallback = ref_obs_doy_fallback,
    observability_undetermined_policy = observability_undetermined_policy
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
      select = "in_evaluation_domain"
    )
    n_excluded_cached <- sum(reference_observability_meta$in_evaluation_domain == FALSE, na.rm = TRUE)
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
      n_reference_excluded_observability <- sum(!reference_observability$in_evaluation_domain, na.rm = TRUE)
      n_reference_observable <- sum(reference_observability$in_evaluation_domain, na.rm = TRUE)
      warn_reference_observability(reference_observability)
    }
    # §N+25: the excluded-fire geometries are needed to build the evaluation
    # domain under "wholefire_fraction". On the cache-hit path they were only
    # ever written to disk, never re-read, so read them back here. (The cache
    # is rebuilt when the GeoPackage is missing but exclusions exist -- see the
    # consistency guard above -- so a NULL here means there are none.)
    if (!is.null(ref_not_obs_cache) && file.exists(ref_not_obs_cache)) {
      ref_not_observable <- sf::st_read(ref_not_obs_cache, quiet = TRUE)
      if (!nrow(ref_not_observable)) ref_not_observable <- NULL
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
      # Stable per-fire key so the observability verdict of every ORIGINAL
      # reference fire can be joined back to its detection outcome after the
      # domain subset. It matches `reference_row` in the observability table.
      # Survives the subset and the GPKG round-trip; a dissolve legitimately
      # destroys it, which is how the audit CSV detects that it cannot be built.
      ref_polygons$reference_row_id <- seq_len(nrow(ref_polygons))
      reference_observability <- build_reference_observability(
        ref_polygons = ref_polygons,
        ref_mask_r = ref_mask_r,
        observability_raster = observability_raster,
        domain_mask = domain_mask,
        mask_v = mask_v,
        cell_area_ha = cell_area_ha,
        ref_end_doy_col = ref_end_doy_col,
        ref_start_doy_col = ref_start_doy_col,
        ref_obs_doy_col = ref_obs_doy_col,
        ref_obs_doy_fallback = ref_obs_doy_fallback,
        observability_undetermined_policy = observability_undetermined_policy,
        observability_mode = observability_mode,
        observability_min_fraction = observability_min_fraction
      )
      ref_not_observable <- build_not_observable_reference_polygons(
        ref_polygons = ref_polygons,
        reference_observability = reference_observability
      )
      n_reference_excluded_observability <- sum(!reference_observability$in_evaluation_domain, na.rm = TRUE)
      n_reference_observable <- sum(reference_observability$in_evaluation_domain, na.rm = TRUE)
      warn_reference_observability(reference_observability)

      ref_polygons <- ref_polygons[reference_observability$in_evaluation_domain, , drop = FALSE]
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

  # ---- evaluation domain under "wholefire_fraction" (§N+25) ----
  # Fires that are NOT_OBSERVABLE, UNDETERMINED_OBS_DATE or
  # NO_OBSERVABILITY_DATA have already been removed from `ref_polygons`, so
  # they contribute no reference-positive pixel. Under the legacy rule their
  # territory nonetheless stayed in the domain as reference = 0, which turns a
  # correct detection there into a false positive. Here their WHOLE geometry is
  # additionally removed from the evaluation domain, symmetrically for the
  # reference, the prediction and the strata, so that they can produce neither
  # FP nor FN nor TN. Everything outside those geometries is untouched, so
  # commission outside reference fires keeps being measured exactly as before.
  #
  # The exclusion is applied to in-memory copies only: the cached reference
  # raster, the cached prediction raster and the input layers on disk are never
  # rewritten.
  ref_mask_eval <- ref_mask_r
  strata_r_eval <- NULL
  eval_excluded_v <- NULL
  observability_excluded_area_ha <- 0
  if (observability_masks_domain &&
      !is.null(ref_not_observable) && nrow(ref_not_observable) > 0L) {
    eval_excluded_v <- safe_vect(ref_not_observable, label = "ref_not_observable")
    ref_mask_eval <- terra::mask(ref_mask_r, eval_excluded_v,
                                 inverse = TRUE, updatevalue = NA,
                                 touches = FALSE)
    n_dom_before <- terra::global(!is.na(ref_mask_r), "sum", na.rm = TRUE)[[1]]
    n_dom_after  <- terra::global(!is.na(ref_mask_eval), "sum", na.rm = TRUE)[[1]]
    observability_excluded_area_ha <- (n_dom_before - n_dom_after) * cell_area_ha
  }
  evaluable_area_ha <-
    terra::global(!is.na(ref_mask_eval), "sum", na.rm = TRUE)[[1]] * cell_area_ha

  # `touches = FALSE` is NOT the terra default for a SpatVector mask (terra
  # 1.8.29 removes every cell the polygon TOUCHES, which for the small, ragged
  # fire perimeters here over-excludes by ~19% of their area and would silently
  # drop a fringe of legitimately evaluable territory around each excluded
  # fire). Forcing it to FALSE makes the removed set exactly the set that
  # terra::rasterize()/terra::extract() would mark, i.e. the same cell-centre
  # rule already used to build `ref_mask_r` and `ref_area_domain_ha`, so the
  # exclusion is aligned with the reference rasterization instead of being a
  # slightly dilated version of it.
  apply_observability_domain <- function(r) {
    if (is.null(eval_excluded_v)) return(r)
    terra::mask(r, eval_excluded_v, inverse = TRUE, updatevalue = NA,
                touches = FALSE)
  }

  # ---- annual observability summary ----
  observability_summary <- NULL
  if (!is.null(reference_observability) && nrow(reference_observability) > 0L) {
    ro_sum <- as.data.frame(reference_observability, stringsAsFactors = FALSE)
    st <- if ("observability_status" %in% names(ro_sum)) {
      as.character(ro_sum$observability_status)
    } else {
      ifelse(as.logical(ro_sum$observable_flag), "OBSERVABLE", "NOT_OBSERVABLE")
    }
    ha <- suppressWarnings(as.numeric(ro_sum$ref_area_domain_ha))
    ha[!is.finite(ha)] <- 0
    obs_ok <- st == "OBSERVABLE"
    observability_summary <- data.frame(
      Year = year_target,
      Observability_Mode = observability_mode,
      Observability_Min_Fraction =
        if (identical(observability_mode, "legacy_any")) NA_real_ else
          observability_min_fraction,
      N_Reference_Fires_Total = nrow(ro_sum),
      N_Observable = sum(obs_ok),
      N_Not_Observable = sum(st == "NOT_OBSERVABLE"),
      N_Undetermined_Obs_Date = sum(st == "UNDETERMINED_OBS_DATE"),
      N_No_Observability_Data = sum(st == "NO_OBSERVABILITY_DATA"),
      N_Outside_Burnable_Domain = sum(st == "OUTSIDE_BURNABLE_DOMAIN"),
      # EVALUATED = what actually enters the metrics. Differs from
      # N_Observable whenever undetermined fires are being evaluated.
      Undetermined_Policy = if (identical(observability_mode, "legacy_any"))
        NA_character_ else observability_undetermined_policy,
      N_Evaluated = sum(as.logical(ro_sum$in_evaluation_domain), na.rm = TRUE),
      Area_Evaluated_ha = round(sum(ha[as.logical(ro_sum$in_evaluation_domain)]), 4),
      Area_Reference_Total_ha = round(sum(ha), 4),
      Area_Observable_ha = round(sum(ha[obs_ok]), 4),
      Area_Excluded_ha = round(sum(ha[!as.logical(ro_sum$in_evaluation_domain)]), 4),
      # Per-fire sums, so they double-count any duplicated or overlapping
      # reference geometry; `Domain_Removed_By_Observability_ha` below is the
      # spatial union and is the one that governs the metrics.
      Area_Not_Observable_ha = round(sum(ha[st == "NOT_OBSERVABLE"]), 4),
      Area_Undetermined_Obs_Date_ha =
        round(sum(ha[st == "UNDETERMINED_OBS_DATE"]), 4),
      Area_No_Observability_Data_ha =
        round(sum(ha[st == "NO_OBSERVABILITY_DATA"]), 4),
      Area_Outside_Burnable_Domain_ha =
        round(sum(ha[st == "OUTSIDE_BURNABLE_DOMAIN"]), 4),
      Domain_Removed_By_Observability_ha = round(observability_excluded_area_ha, 4),
      Evaluable_Area_ha = round(evaluable_area_ha, 4),
      Domain_Masked = observability_masks_domain,
      Ref_Obs_Doy_Col = if (is.null(ref_obs_doy_col)) NA_character_ else ref_obs_doy_col,
      Ref_Obs_Doy_Authoritative =
        !is.null(ref_obs_doy_col) && identical(ref_obs_doy_fallback, "none"),
      Undetermined_Causes = if ("undetermined_cause" %in% names(ro_sum)) {
        cz <- ro_sum$undetermined_cause[st == "UNDETERMINED_OBS_DATE"]
        cz <- cz[!is.na(cz)]
        if (!length(cz)) NA_character_ else
          paste(sprintf("%s=%d", names(table(cz)), as.integer(table(cz))),
                collapse = "; ")
      } else NA_character_,
      stringsAsFactors = FALSE
    )
    data.table::fwrite(
      observability_summary,
      validation_output_path(observability_output_dir,
        sprintf("OBSERVABILITY_SUMMARY_%s.csv", year_target))
    )
    # An authoritative observability column that yields UNDETERMINED fires is a
    # deliberate signal from the reference, but it is also exactly what a
    # mis-specified column name or a half-populated column looks like. Say so
    # out loud, with counts, rather than letting a large silent exclusion pass.
    if (!is.null(ref_obs_doy_col) &&
        identical(ref_obs_doy_fallback, "none") &&
        identical(observability_mode, "wholefire_fraction")) {
      n_und_auth <- sum(st == "UNDETERMINED_OBS_DATE")
      if (n_und_auth > 0L) {
        warning(sprintf(
          paste0("Authoritative observability column '%s': %d of %d reference ",
                 "fires (%.1f%%) have no value there and are reported as ",
                 "UNDETERMINED_OBS_DATE (%s). No fallback to '%s' was applied. ",
                 "Pass ref_obs_doy_fallback = \"end_doy\" only if that is a ",
                 "deliberate methodological choice."),
          ref_obs_doy_col, n_und_auth, nrow(ro_sum),
          100 * n_und_auth / nrow(ro_sum),
          if (identical(observability_undetermined_policy, "evaluate"))
            "kept as positive reference features inside the evaluation domain"
          else "excluded from the evaluation domain",
          ref_end_doy_col),
          call. = FALSE)
      }
    }

    # Fires without a valid observability date, reported separately by source.
    und <- st == "UNDETERMINED_OBS_DATE"
    src_col <- intersect(c("origin", "source"), names(ro_sum))
    src <- if (length(src_col)) as.character(ro_sum[[src_col[1]]]) else
      rep(NA_character_, nrow(ro_sum))
    undetermined_by_source <- if (any(und)) {
      agg <- stats::aggregate(
        list(n_fires = rep(1L, sum(und)), area_domain_ha = ha[und]),
        by = list(source = ifelse(is.na(src[und]), "<NA>", src[und])),
        FUN = sum
      )
      data.frame(Year = year_target,
                 observability_status = "UNDETERMINED_OBS_DATE",
                 agg, stringsAsFactors = FALSE)
    } else {
      # Zero-row frame: Year must be length 0 too, otherwise data.frame()
      # refuses to recycle a scalar against empty columns.
      data.frame(Year = year_target[0], observability_status = character(0),
                 source = character(0), n_fires = integer(0),
                 area_domain_ha = numeric(0), stringsAsFactors = FALSE)
    }
    data.table::fwrite(
      undetermined_by_source,
      validation_output_path(observability_output_dir,
        sprintf("OBSERVABILITY_UNDETERMINED_%s.csv", year_target))
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
    # §N+25: the per-stratum tabulator turns NA pred/ref inside the strata
    # domain into 0 when `na_strata_as_zero = TRUE` (the default), which would
    # silently undo the observability exclusion and count the excluded fires as
    # TN. Removing the excluded geometries from the STRATA raster takes those
    # cells out of the strata domain instead, keeping the exclusion effective.
    strata_r_eval <- apply_observability_domain(strata_r)
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
    # §N+25: under "wholefire_fraction" the raster actually used for the metrics
    # depends on the observability configuration (excluded fires are blanked),
    # so the observability fingerprint must be part of the prediction cache key.
    # Under "legacy_any" nothing about the prediction depends on observability,
    # so the key is left byte-identical to previous releases and existing
    # prediction caches keep being reused.
    pred_in_tag <- paste0(.vfm_vector_fingerprint(shp), "_", dom_cache_tag)
    if (observability_masks_domain) {
      pred_in_tag <- paste0(pred_in_tag, "_", obs_cache_tag)
    }
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
      # §N+25: when the observability domain removes territory, the cached
      # raster IS the raster used for the metrics (its key carries the
      # observability fingerprint), so the exclusion is baked in here. The
      # ORIGINAL prediction layer on disk is never modified.
      pred_r <- apply_observability_domain(pred_r)

      terra::writeRaster(
        pred_r, pred_tif, overwrite = TRUE,
        datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255")
      )
    }

    pred_r <- terra::rast(pred_tif)

    # ---- PIXEL METRICS ----
    if (metrics_type %in% c("all", "pixel")) {

      ref_mask_use <- ref_mask_eval
      if (buffer > 0) {
        ref_buf <- suppressWarnings(sf::st_buffer(ref_polygons, dist = buffer)) |> make_valid_sf()
        ref_buf_v <- safe_vect(ref_buf, label = "ref_buffer")
        ref_mask_use <- terra::rasterize(ref_buf_v, domain_mask, field = 1, background = 0)
        ref_mask_use[is.na(domain_mask)] <- NA
        ref_mask_use <- apply_observability_domain(ref_mask_use)
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
        Observability_Mode = observability_mode,
        Observability_Min_Fraction =
          if (identical(observability_mode, "legacy_any")) NA_real_ else
            observability_min_fraction,
        Observability_Domain_Masked = observability_masks_domain,
        Evaluable_Area_ha = round(evaluable_area_ha, 4),
        N_Reference_Polygons_Observable = n_reference_observable,
        N_Reference_Polygons_Excluded_Observability = n_reference_excluded_observability,
        InputName = input_name, Year = year_target
      )
    }

    # ---- AREA / POLYGON METRICS ----
    if (metrics_type %in% c("all", "area")) {

      ref_v <- safe_vect(ref_polygons, label = "ref_for_extract")
      det_v <- safe_vect(det, label = "det_for_extract")

      ref_pix <- terra::extract(ref_mask_eval, ref_v, fun = sum, na.rm = TRUE)[, 2]
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
      inter_total_pix <- terra::global((pred_r == 1) & (ref_mask_eval == 1), "sum", na.rm = TRUE)[[1]]
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
        Observability_Mode = observability_mode,
        Observability_Min_Fraction =
          if (identical(observability_mode, "legacy_any")) NA_real_ else
            observability_min_fraction,
        Observability_Domain_Masked = observability_masks_domain,
        Evaluable_Area_ha = round(evaluable_area_ha, 4),
        N_Reference_Polygons_Observable = n_reference_observable,
        N_Reference_Polygons_Excluded_Observability = n_reference_excluded_observability,
        Detected_Definition = detected_definition
      )

      # ---- per-fire detection audit (§N+26) ----
      # One row per EVALUATED reference fire, carrying its observability verdict
      # next to what the prediction actually did on it. This is what makes it
      # possible to revisit an individual fire -- e.g. every
      # observability_status == "UNDETERMINED_OBS_DATE" -- without recomputing
      # the whole validation. Written only when the fires are still traceable:
      # a dissolve merges several originals into one geometry and destroys the
      # per-fire correspondence, so the join key is gone and the file is skipped.
      if (!is.null(reference_observability) &&
          "reference_row_id" %in% names(ref_polygons)) {
        ro_join <- as.data.frame(reference_observability, stringsAsFactors = FALSE)
        rid <- as.integer(ref_polygons$reference_row_id)
        keep_cols <- intersect(
          c("id", "fire_id", "origin", "source", "source_type", "country",
            "region", "year", "area_ha",
            "observability_status", "undetermined_cause",
            "obs_required_doy", "fire_required_doy", "required_doy_source",
            "observable_flag", "in_evaluation_domain",
            "data_coverage_fraction", "temporal_observable_fraction",
            "obs_doy_max", "total_pixels", "ref_area_domain_ha"),
          names(ro_join))
        det_audit <- data.frame(
          Year = year_target, InputName = input_name,
          reference_row = rid,
          ro_join[rid, keep_cols, drop = FALSE],
          reference_cells = ref_pix,
          reference_area_ha = round(ref_area_i, 4),
          detected_cells = int_pix,
          detected_area_ha = round(int_area_i, 4),
          omitted_cells = pmax(ref_pix - int_pix, 0),
          omitted_area_ha = round(pmax(ref_area_i - int_area_i, 0), 4),
          # per-fire confusion counts restricted to the fire footprint
          TP_cells = int_pix, TP_area_ha = round(int_area_i, 4),
          FN_cells = pmax(ref_pix - int_pix, 0),
          FN_area_ha = round(pmax(ref_area_i - int_area_i, 0), 4),
          fire_recall_percent = round(perc_detected_i, 4),
          detected_flag = detected_i,
          stringsAsFactors = FALSE, row.names = NULL
        )
        data.table::fwrite(
          det_audit,
          validation_output_path(
            observability_output_dir,
            sprintf("REFERENCE_FIRES_DETECTION_%s_%s.csv", year_target, input_name)
          )
        )
        und_audit <- det_audit[
          det_audit$observability_status == "UNDETERMINED_OBS_DATE", , drop = FALSE]
        data.table::fwrite(
          und_audit,
          validation_output_path(
            observability_output_dir,
            sprintf("UNDETERMINED_FIRES_DETECTION_%s_%s.csv", year_target, input_name)
          )
        )
      }

      ref_not_detected <- ref_polygons[!detected_i, , drop = FALSE]
      det_ref_pix <- terra::extract(ref_mask_eval, det_v, fun = sum, na.rm = TRUE)[, 2]
      det_ref_pix[is.na(det_ref_pix)] <- 0
      det_not_matched <- det[det_ref_pix == 0, , drop = FALSE]

      # ---- Block 10c: reconcile size attributes with the written geometry ----
      # The geometry of `ref_polygons` / `det` was clipped by
      # sf::st_intersection(*, mask_geom) (and st_make_valid / optional
      # dissolve), but `area_ha` (and `n_pix`) still carry the PRE-CLIP values
      # inherited from the scored map. Downstream error-layer consumers
      # (EGIF cross-check, commission characterisation, ...) read these
      # attributes as the polygon size, so recompute them from the geometry
      # that is about to be written. `n_pix` is by construction the
      # rasterised-cell count of the polygon (area = n_pix * cell_area_ha),
      # so recomputing it consistently from the corrected area is valid for
      # this post-hoc, vector-clipped layer.
      reconcile_size <- function(g) {
        if (nrow(g) == 0L) return(g)
        a_ha <- suppressWarnings(as.numeric(sf::st_area(g))) / 10000
        if ("area_ha" %in% names(g)) g$area_ha <- a_ha
        if ("n_pix"   %in% names(g)) g$n_pix   <- as.integer(round(a_ha / cell_area_ha))
        g
      }
      ref_not_detected <- reconcile_size(ref_not_detected)
      det_not_matched  <- reconcile_size(det_not_matched)

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
        ref_r         = ref_mask_eval,
        strata_r      = strata_r_eval,
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
    observability_summary = observability_summary,
    observability_settings = list(
      mode = observability_mode,
      min_fraction = if (identical(observability_mode, "legacy_any")) NA_real_
        else observability_min_fraction,
      ref_obs_doy_col = ref_obs_doy_col,
      ref_obs_doy_fallback = ref_obs_doy_fallback,
      ref_obs_doy_authoritative =
        !is.null(ref_obs_doy_col) && identical(ref_obs_doy_fallback, "none"),
      domain_masked = observability_masks_domain,
      evaluable_area_ha = evaluable_area_ha,
      domain_removed_ha = observability_excluded_area_ha,
      cache_tag = obs_cache_tag
    ),
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

