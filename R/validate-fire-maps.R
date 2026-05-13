#' Validate burned-area maps against reference polygons
#'
#' @description
#' Shared, workflow-independent validation utility and the single public
#' validation function in OtsuFire 0.2.x. Computes pixel-based metrics
#' (confusion-matrix statistics over the burnable domain) and area-based
#' metrics (polygon-level detection rates) by comparing one or more
#' burned-area prediction layers against independent reference fire polygons
#' (e.g., EFFIS or a national fire atlas). Reference polygons are clipped to
#' the study area, masked to the burnable domain, optionally filtered by a
#' minimum area threshold, and optionally dissolved by a grouping field
#' before comparison.
#'
#' Optional class-wise and pixel-level stratified outputs are produced
#' when the corresponding inputs are supplied:
#' * **Class-wise polygon breakdown** when both `class_shape` and
#'   `class_field` are set: counts and hectares of un-detected
#'   reference polygons and un-matched detection polygons, grouped by
#'   the polygon class attribute (e.g. ecoregion).
#' * **Pixel-level per-stratum confusion matrices** when
#'   `strata_raster` (and optionally `strata_lut`) are set: TP / FP /
#'   FN / TN per stratum plus derived metrics, computed over the
#'   strata-defined evaluation domain via a chunked tabulator
#'   (matches the legacy `pixel_by_stratum` / `stratum_global` /
#'   `diagnostics_strata` outputs).
#' * **Combined Excel workbook** when `write_excel = TRUE`:
#'   `validation_ALL_<year>_res30.xlsx` with sheets `pixel_global`,
#'   `pixel_by_stratum`, `stratum_global`, `diagnostics_strata`,
#'   `polygon_global`, `run_info`.
#'
#' This function is workflow-independent: it can be called by
#' [run_deterministic_pipeline()] internally, and it can be called directly
#' by users on any compatible burned-area prediction layer, including
#' supervised outputs after thresholding `p_burned` or
#' `p_burned_current_year`. There is no `validate_burned_maps()` wrapper in
#' 0.2.x (see `NEWS.md`).
#'
#' Intermediate masked-reference and masked-detection layers are cached on
#' disk. If `min_area_reference_ha` is changed between runs, set
#' `force_reprocess_ref = TRUE` to force recalculation of the reference
#' cache. Use `force_reprocess_pred = TRUE` to force recalculation of the
#' detection cache.
#'
#' @param input_shapefile Character vector or `sf` object. One or more paths
#'   to shapefiles containing the burned polygons to validate, or a single
#'   in-memory `sf`. Each input is processed independently and produces one
#'   row in the output metrics tables.
#' @param ref_shapefile Character path or `sf` object. Reference burned-area
#'   polygons.
#' @param mask_shapefile Character. Path to the study-area boundary shapefile.
#' @param burnable_raster Character. Path to the burnable raster.
#' @param year_target Numeric. Target year used to filter reference polygons
#'   to the correct fire season.
#' @param validation_dir Character. Root output directory. A `VALIDATION/`
#'   subfolder is created inside it to hold all output files.
#' @param binary_burnable Logical. If `TRUE` (default), treats the burnable
#'   raster as binary: cells at or above `burnable_threshold` are burnable
#'   and cells below are masked as `NA`. If `FALSE`, use `burnable_classes`
#'   to specify which categorical raster values count as burnable.
#' @param burnable_classes Optional numeric vector. Raster values considered
#'   burnable when `binary_burnable = FALSE`.
#' @param burnable_threshold Numeric scalar (0–1). Threshold applied when
#'   `binary_burnable = TRUE`. Default `0.5`.
#' @param class_shape Optional character. Path to a shapefile with landscape
#'   class polygons (e.g., CORINE land-cover classes or ecoregions) used for
#'   class-wise omission/commission breakdowns. Requires `class_field`.
#' @param class_field Optional character. Attribute field name in
#'   `class_shape` to group errors by class.
#' @param buffer Numeric. Buffer distance (metres) applied around reference
#'   polygons before pixel-based comparison. Default `0`.
#' @param threshold_completely_detected Numeric (0–100). Minimum percentage
#'   of a reference polygon area that must be covered by detections to count
#'   it as completely detected. Default `90`.
#' @param threshold_min_detected Numeric scalar (0–100). Minimum percentage
#'   of a reference polygon area that must be covered by detections to
#'   count it as "detected at all". This is the threshold used to compute
#'   `N_Detected_Polygons` (compare with `threshold_completely_detected`,
#'   which is stricter and counts `N_Completely_Detected`). Default `10`.
#' @param min_area_reference_ha Numeric. Minimum area (ha) a reference
#'   polygon must have after masking to the burnable domain to be retained.
#'   Default `NULL` (no filter).
#' @param force_reprocess_ref Logical. Force recalculation of the reference
#'   cache even if a cached version exists. Default `FALSE`.
#' @param force_reprocess_pred Logical. Force recalculation of the detection
#'   cache. Default `FALSE`.
#' @param metrics_type Character. Which metrics to compute. One of `"all"`
#'   (default), `"pixel"`, or `"area"`.
#' @param dissolve_ref_by Optional character. Field used to dissolve
#'   reference polygons before comparison (e.g., `"fire_id"`).
#' @param dissolve_input_by Optional character. Field used to dissolve
#'   detection polygons before comparison (e.g., `"fire_id"`).
#' @param strata_raster Optional. SpatRaster, raster path, or NULL.
#'   When supplied, per-stratum pixel-level confusion matrices
#'   (TP/FP/FN/TN) and derived metrics are computed in addition to the
#'   global pixel + polygon outputs. The raster must be aligned (or
#'   resampleable) to the burnable domain. Strata are integer ids;
#'   pixels with `NA` strata are excluded from per-stratum aggregation
#'   (they still count for the global pixel metrics).
#' @param strata_lut Optional. `data.frame` with columns `id` and
#'   `label`, or a path to a CSV/XLSX with the same schema. Used to
#'   attach human-readable `Stratum_Label` to each per-stratum row.
#' @param chunk_rows Integer. Row chunk size used by the per-stratum
#'   tabulator (memory-vs-speed knob). Default `1024L`.
#' @param na_strata_as_zero Logical. When `TRUE` (default), `NA`
#'   prediction or reference pixels inside the strata domain are
#'   treated as 0 (unburned), matching the legacy chunked validator.
#'   Set `FALSE` to drop them.
#' @param write_excel Logical. When `TRUE`, also writes a combined
#'   Excel workbook (`validation_ALL_<year>_res30.xlsx`) with sheets
#'   `pixel_global`, `pixel_by_stratum`, `stratum_global`,
#'   `diagnostics_strata`, `polygon_global`, `run_info`. Requires the
#'   `openxlsx` package. Default `FALSE`.
#' @param excel_filename Optional character. Override for the Excel
#'   filename. Default builds `validation_ALL_<year_target>_res30.xlsx`
#'   per input.
#'
#' @details
#' ## Pixel-based metrics (`metrics_type` = `"pixel"` or `"all"`)
#' TP / FP / FN / TN over the burnable domain, plus Precision, Recall,
#' F1, IoU, `Specificity = TN / (TN + FP)`,
#' `BalancedAccuracy = mean(c(Recall, Specificity), na.rm = TRUE)`,
#' and `ErrorRate = (FP + FN) / (TP + FP + FN + TN)`. Written to
#' `validation_dir/VALIDATION/metrics_summary_<year>.csv`. Returned in
#' the `metrics` element of the result list.
#'
#' ## Area-based metrics (`metrics_type` = `"area"` or `"all"`)
#' `N_Reference_Polygons`, `N_Completely_Detected`, `N_Detected_Polygons`,
#' `N_Not_Detected`, `Perc_Detected_Polygons`, `Area_Reference_ha`,
#' `Area_Detected_ha`, `Area_Intersection_ha`, `Recall_Area_percent`,
#' `Precision_Area_percent`, plus per-polygon coverage summaries
#' `Coverage_mean`, `Coverage_median`, `Coverage_p10`, `Coverage_p90`
#' (rounded to two decimals; `NA_real_` when no reference polygons
#' remain) and a `Detected_Definition` string documenting the
#' "detected at all" threshold actually applied (`sprintf("coverage_ref
#' >= %d%%", threshold_min_detected)`). Written to
#' `polygon_summary_<year>.csv`. Returned in the `polygon_summary`
#' element of the result list.
#'
#' Note that `N_Detected_Polygons` is governed by
#' `threshold_min_detected` (default `10`). Prior to 0.2.1 the rule
#' was effectively any non-zero overlap (`int_area_i > 0`); the new
#' default is stricter and matches the legacy 1988 validation
#' workbook (`coverage_ref >= 10%`). Pass
#' `threshold_min_detected = 0` to recover the pre-0.2.1 behaviour.
#'
#' ## Class-wise polygon breakdown (`class_shape` + `class_field`)
#' Optional polygon-level summary that groups un-detected reference
#' polygons and un-matched detection polygons by an attribute of an
#' external class layer (e.g. ecoregion shapefile). When both
#' `class_shape` and `class_field` are provided and valid, two CSVs
#' are written to `validation_dir/VALIDATION/`:
#' \itemize{
#'   \item `omission_by_<class_field>_<input>.csv` --
#'     columns `class`, `N_Omitted`, `Area_Omitted_ha`.
#'   \item `commission_by_<class_field>_<input>.csv` --
#'     columns `class`, `N_Commission`, `Area_Commission_ha`.
#' }
#' This branch is **polygon-level**: it counts whole reference /
#' detection polygons and their hectares per class. It does not
#' produce a per-class confusion matrix.
#'
#' ## Pixel-level per-stratum metrics (`strata_raster` + `strata_lut`)
#' Optional pixel-level stratification using an integer-id raster.
#' Mirrors the legacy `validation_ALL_<year>_<scenario>_res30.xlsx`
#' workbook. When `strata_raster` is supplied, three CSVs are
#' written to `validation_dir/VALIDATION/` per input:
#' \itemize{
#'   \item `pixel_by_stratum_<year>_<input>.csv` --
#'     `Stratum_ID`, `Stratum_Label`, `TP`, `FP`, `FN`, `TN`,
#'     `Precision`, `Recall`, `F1`, `IoU`, `Specificity`,
#'     `BalancedAccuracy`, `Commission`, `Omission`, `ErrorRate`,
#'     plus `Area_*_ha` and `*_Area_percent` columns.
#'   \item `stratum_global_<year>_<input>.csv` -- single-row
#'     aggregate over the strata-defined evaluation domain
#'     (counts, derived ratios, areas, `cell_area_ha`).
#'   \item `diagnostics_strata_<year>_<input>.csv` -- evaluation
#'     sample sizes (`n_eval_pixels`), NA counts inside the
#'     strata domain (`na_pred_inside_domain`,
#'     `na_ref_inside_domain`) and their fractions, plus the
#'     `chunk_rows` and `na_as_zero` settings.
#' }
#' The algorithm is a chunked-row tabulator (`chunk_rows`)
#' that streams `pred_r`, `ref_r`, and `strata_r` together; pixels
#' with `NA` strata are excluded from per-stratum aggregation.
#' Inside the strata domain, `NA` predictions / references are
#' treated as 0 (unburned) by default (`na_strata_as_zero = TRUE`,
#' matching the legacy validator). Set `na_strata_as_zero = FALSE`
#' to drop those pixels instead.
#'
#' This branch is **pixel-level**: it produces a full confusion
#' matrix per stratum from the rasterized prediction and reference.
#' It is independent of the polygon-level `class_shape` branch and
#' uses a different data input (a raster, not a shapefile). Both
#' branches can be active in the same call.
#'
#' ## Combined Excel workbook (`write_excel = TRUE`)
#' When `write_excel = TRUE`, an Excel workbook is also written
#' (default name `validation_ALL_<year>_res30.xlsx`, override via
#' `excel_filename`) with the legacy six-sheet layout:
#' `pixel_global`, `pixel_by_stratum`, `stratum_global`,
#' `diagnostics_strata`, `polygon_global`, `run_info`. Sheets that
#' have no data (e.g. `pixel_by_stratum` when `strata_raster` is
#' NULL) are omitted. Requires the `openxlsx` package.
#'
#' The `run_info` sheet captures the call-time arguments needed to
#' reproduce the run: `year_target`, `validation_dir`,
#' `binary_burnable`, `burnable_threshold`, `buffer`,
#' `threshold_completely_detected`, `threshold_min_detected`,
#' `min_area_reference_ha`, `metrics_type`, `dissolve_ref_by`,
#' `dissolve_input_by`, `strata_used`, `chunk_rows`,
#' `na_strata_as_zero`, a derived `cache_tag`
#' (`y<year>_res<m>_dom{strata|burnable}_na0<TRUE|FALSE>`), and a
#' `produced_at` timestamp.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`metrics`}{Pixel-based global confusion `data.table`,
#'     or `NULL` when `metrics_type = "area"`.}
#'   \item{`polygon_summary`}{Area-based polygon-level
#'     `data.table`, or `NULL` when `metrics_type = "pixel"`.}
#'   \item{`pixel_by_stratum`}{Per-stratum confusion + derived
#'     metrics `data.table`, or `NULL` when `strata_raster` is
#'     not supplied.}
#'   \item{`stratum_global`}{Single-row strata-domain aggregate
#'     `data.table`, or `NULL` when `strata_raster` is not
#'     supplied.}
#'   \item{`diagnostics_strata`}{Diagnostics `data.table`
#'     (sample sizes, NA fractions), or `NULL` when
#'     `strata_raster` is not supplied.}
#'   \item{`excel_path`}{Path to the combined Excel workbook,
#'     or `NULL` when `write_excel = FALSE`.}
#' }
#'
#' @note
#' Canonical v2 implementation (post 2026-04-16). The deprecated v1
#' GDAL/Python parameters (`use_gdal`, `python_exe`,
#' `gdal_polygonize_script`) are gone. Includes the F8b GDAL-generated FID
#' drop that was required for GeoPackage writes to succeed after
#' `sf::st_intersection()` propagated the FID field.
#'
#' Block 10c added the optional `strata_raster` / `strata_lut` /
#' `chunk_rows` / `na_strata_as_zero` / `write_excel` /
#' `excel_filename` arguments to restore the legacy
#' pixel-level stratified outputs (`pixel_by_stratum`,
#' `stratum_global`, `diagnostics_strata`) and the combined Excel
#' workbook. All of these arguments are optional and default to
#' the pre-Block-10c behaviour, so existing call sites continue to
#' work unchanged.
#'
#' Version 0.2.1 added `Specificity` / `BalancedAccuracy` /
#' `ErrorRate` to `pixel_global` / `pixel_by_stratum` /
#' `stratum_global`; `Coverage_mean` / `Coverage_median` /
#' `Coverage_p10` / `Coverage_p90` and `Detected_Definition` to
#' `polygon_global`; the new `threshold_min_detected` argument
#' (default `10`); and `dissolve_ref_by`, `dissolve_input_by`,
#' `cache_tag` to `run_info`. **Behaviour change:** the default
#' `threshold_min_detected = 10` is stricter than the pre-0.2.1
#' rule (any non-zero overlap counted as "detected"). To preserve
#' previous behaviour, pass `threshold_min_detected = 0`
#' explicitly.
#'
#' @seealso [run_deterministic_pipeline()] for the deterministic
#'   workflow that produces the prediction layer typically passed to
#'   `input_shapefile`; [run_oneyear_supervised_pipeline()] for the
#'   supervised one-year workflow whose `final_map.gpkg` (after
#'   thresholding `p_burned`) is also a valid input.
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
#' # Produces pixel_by_stratum / stratum_global / diagnostics_strata
#' # CSVs plus validation_ALL_2022_res30.xlsx with the 6 legacy sheets.
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

  # ---- output dir ----
  validation_output_dir <- file.path(validation_dir, "VALIDATION")
  if (!dir.exists(validation_output_dir)) dir.create(validation_output_dir, recursive = TRUE)

  # ---- load mask & burnable ----
  mask_geom <- sf::st_read(mask_shapefile, quiet = TRUE) |> make_valid_sf()
  # F8b: drop GDAL-generated FID column if present; it propagates through
  # st_intersection into ref_polygons and causes a GeoPackage write failure
  # ('Wrong field type for FID').
  mask_geom <- mask_geom[, !(names(mask_geom) %in% "FID"), drop = FALSE]
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

  # ---- cache paths ----
  ref_vec_cache  <- file.path(validation_output_dir, paste0("ref_polygons_processed_", year_target, ".gpkg"))
  ref_rast_cache <- file.path(validation_output_dir, paste0("ref_mask_", year_target, ".tif"))

  if (force_reprocess_ref) {
    message("Force reprocessing reference polygons: deleting cached reference...")
    unlink(ref_vec_cache)
    unlink(ref_rast_cache)
  }

  # ---- build/load reference (vector + raster mask) ----
  if (file.exists(ref_vec_cache) && file.exists(ref_rast_cache)) {
    message("Loading cached reference (vector + raster mask)...")
    ref_polygons <- sf::st_read(ref_vec_cache, quiet = TRUE) |> make_valid_sf()
    ref_mask_r   <- terra::rast(ref_rast_cache)
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

    ref_polygons <- dissolve_by_field(ref_polygons, dissolve_ref_by) |> make_valid_sf()

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

    sf::st_write(ref_polygons, ref_vec_cache, delete_dsn = TRUE, quiet = TRUE)
    terra::writeRaster(
      ref_mask_r, ref_rast_cache, overwrite = TRUE,
      datatype = "INT1U", gdal = c("COMPRESS=LZW", "NAflag=255")
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

    pred_tif <- file.path(validation_output_dir, paste0("pred_", input_name, "_", year_target, ".tif"))
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
        Detected_Definition = detected_definition
      )

      ref_not_detected <- ref_polygons[!detected_i, , drop = FALSE]
      det_ref_pix <- terra::extract(ref_mask_r, det_v, fun = sum, na.rm = TRUE)[, 2]
      det_ref_pix[is.na(det_ref_pix)] <- 0
      det_not_matched <- det[det_ref_pix == 0, , drop = FALSE]

      if (nrow(ref_not_detected) > 0) {
        sf::st_write(
          ref_not_detected,
          file.path(validation_output_dir, paste0("ref_polygons_not_detected_", input_name, ".gpkg")),
          delete_dsn = TRUE, quiet = TRUE
        )
      }
      if (nrow(det_not_matched) > 0) {
        sf::st_write(
          det_not_matched,
          file.path(validation_output_dir, paste0("input_polygons_not_matched_", input_name, ".gpkg")),
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
              file.path(validation_output_dir, paste0("omission_by_", class_field, "_", input_name, ".csv"))
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
              file.path(validation_output_dir, paste0("commission_by_", class_field, "_", input_name, ".csv"))
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
        file.path(validation_output_dir,
                  paste0("pixel_by_stratum_", year_target, "_",
                         input_name, ".csv"))
      )
      data.table::fwrite(
        st_out$global,
        file.path(validation_output_dir,
                  paste0("stratum_global_", year_target, "_",
                         input_name, ".csv"))
      )
      data.table::fwrite(
        st_out$diagnostics,
        file.path(validation_output_dir,
                  paste0("diagnostics_strata_", year_target, "_",
                         input_name, ".csv"))
      )
    }
  }

  # ---- write outputs ----
  all_metrics <- NULL
  all_polygon_summary <- NULL

  if (metrics_type %in% c("all", "pixel")) {
    all_metrics <- data.table::rbindlist(metrics_list, fill = TRUE)
    data.table::fwrite(all_metrics, file.path(validation_output_dir, paste0("metrics_summary_", year_target, ".csv")))
  }

  if (metrics_type %in% c("all", "area")) {
    all_polygon_summary <- data.table::rbindlist(polygon_summary_list, fill = TRUE)
    data.table::fwrite(all_polygon_summary, file.path(validation_output_dir, paste0("polygon_summary_", year_target, ".csv")))
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
      excel_path <- file.path(validation_output_dir, excel_name)
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
        "y%s_res%d_dom%s_na0%s",
        as.character(year_target),
        as.integer(round(terra::res(domain_mask)[1])),
        if (!is.null(strata_r)) "strata" else "burnable",
        isTRUE(na_strata_as_zero)
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
