#' Align a (binary) categorical mask to a template raster, with validation
#'
#' @description
#' Single, reusable internal helper for aligning a categorical / binary mask
#' (e.g. a burnable mask) onto a template raster's grid, CRS, extent and
#' dimensions, and then validating the result. It is the ONE place in the
#' supervised stack where a mask is brought onto a template geometry; every
#' call site that needs this MUST route through here so the behaviour is
#' identical everywhere.
#'
#' @details
#' ## Alignment logic (the core fix)
#' The previous ad-hoc pattern aligned a mask with [terra::resample()] guarded
#' only by [terra::compareGeom()]. When the mask and the template were in
#' **different CRS**, `compareGeom()` reported a mismatch but `resample()` (which
#' assumes a shared CRS and merely re-grids) silently produced an all-`NaN` /
#' zero-burnable result, because there is no spatial overlap in the mask's own
#' coordinate system. This helper fixes that by branching explicitly:
#'
#' \preformatted{
#'   if (!terra::same.crs(mask, template)) {
#'     aligned <- terra::project(mask, template, method = "near")
#'   } else if (!terra::compareGeom(mask, template, stopOnError = FALSE)) {
#'     aligned <- terra::resample(mask, template, method = "near")
#'   } else {
#'     aligned <- mask
#'   }
#' }
#'
#' `method = "near"` is used throughout: a categorical / binary mask must never
#' be continuously interpolated (bilinear would invent fractional class codes).
#' When the mask already matches the template (same CRS and same geometry) it is
#' returned unchanged.
#'
#' ## Binary rules (when `binary = TRUE`)
#' - Only the values `{0, 1, NA}` are accepted. Any other value (e.g. `2` or
#'   `0.5`) is a hard error: the helper does **not** round, clamp or coerce
#'   silently.
#' - `NA` is **not** converted to `0` silently; NA cells stay NA.
#' - A mask with **zero burnable cells** (no cell equal to `1`) is rejected with
#'   a clear error, both before and after alignment.
#'
#' ## Area tolerance
#' Reprojection / re-gridding to a different grid is not area-preserving at the
#' bit level: cells on the boundary are reassigned and the nearest-neighbour
#' rule shifts a thin fringe of cells. We therefore do NOT require the burnable
#' area to be bit-identical before vs after alignment. Instead we record the
#' absolute and relative burnable-area difference and only `stop()` when the
#' RELATIVE difference exceeds `area_tol_rel`.
#'
#' The default tolerance is `0.05` (5 %). Justification: nearest-neighbour
#' resampling/projection perturbs the burnable count by at most a one-cell-wide
#' fringe around the burnable footprint, i.e. on the order of
#' `perimeter / area` of the burnable region; for the compact, large burnable
#' footprints used here that fringe is well under a few percent. 5 % is a
#' deliberately conservative ceiling that catches a genuine alignment failure
#' (e.g. a near-total wipe-out from a silent CRS mismatch, which would show a
#' relative difference approaching 1.0) while tolerating ordinary grid effects.
#' The absolute and relative differences are returned in the report regardless
#' of whether they are within tolerance.
#'
#' @param mask A `SpatRaster` mask to align (single layer). For a binary mask
#'   its values must lie in `{0, 1, NA}`.
#' @param template A `SpatRaster` whose CRS, resolution, extent and dimensions
#'   define the target grid.
#' @param allowed_values Numeric vector of permitted non-NA values. Defaults to
#'   `c(0, 1)`. Ignored when `binary = FALSE` and `allowed_values` is `NULL`.
#' @param binary Logical; when `TRUE` (default) enforce the binary rules above
#'   (`{0,1,NA}` only, at least one burnable cell). When `FALSE` the mask is
#'   treated as a general categorical mask: `allowed_values` (if non-NULL) is
#'   still enforced, but the "burnable" accounting and the zero-burnable error
#'   are skipped.
#' @param area_tol_rel Numeric relative tolerance on the burnable-area change
#'   before vs after alignment. Default `0.05` (5 %). Only meaningful for a
#'   binary mask.
#' @param mask_name Character label used in error messages. Default
#'   `"burnable mask"`.
#'
#' @return A named list with:
#' \describe{
#'   \item{aligned}{The aligned `SpatRaster` (same geometry as `template`).}
#'   \item{report}{A named list: `crs_ok`, `res`, `extent`, `dims`, `overlap`,
#'     `allowed_values_ok`, `n_valid`, `n_burnable`, `area_before`,
#'     `area_after`, `area_abs_diff`, `area_rel_diff`, `within_tol`,
#'     `non_empty`.}
#' }
#'
#' @section Behaviour change:
#' Cross-CRS masks are now **projected** onto the template
#' ([terra::project()], `method = "near"`) instead of being passed straight to
#' [terra::resample()], which silently zeroed all burnable cells whenever the
#' mask CRS differed from the template CRS.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra same.crs compareGeom project resample crs res ext nrow ncol
#'   freq unique values cellSize global ifel as.polygons
.of_align_mask_to_template <- function(mask,
                                       template,
                                       allowed_values = c(0, 1),
                                       binary = TRUE,
                                       area_tol_rel = 0.05,
                                       mask_name = "burnable mask") {

  if (!inherits(mask, "SpatRaster")) {
    stop(sprintf("`%s`: `mask` must be a SpatRaster.", mask_name), call. = FALSE)
  }
  if (!inherits(template, "SpatRaster")) {
    stop(sprintf("`%s`: `template` must be a SpatRaster.", mask_name), call. = FALSE)
  }
  if (!is.numeric(area_tol_rel) || length(area_tol_rel) != 1L ||
      !is.finite(area_tol_rel) || area_tol_rel < 0) {
    stop(sprintf("`%s`: `area_tol_rel` must be a single non-negative number.",
                 mask_name), call. = FALSE)
  }

  # ---- helper: count burnable cells (value == 1) without materialising all values ----
  count_value <- function(r, v) {
    fr <- terra::freq(r, value = v)
    if (is.null(fr) || nrow(fr) == 0L) return(0)
    sum(fr[, "count"])
  }
  # ---- helper: burnable AREA (m^2-equivalent) for cells == 1 ----
  burnable_area <- function(r) {
    cs <- terra::cellSize(r, unit = "m")
    burn01 <- terra::ifel(!is.na(r) & r == 1, 1, NA)
    a <- terra::global(terra::ifel(!is.na(burn01), cs, NA), "sum", na.rm = TRUE)[1, 1]
    if (is.na(a)) 0 else a
  }

  # ---- pre-alignment binary validation (fail before any heavy reproject) ----
  if (isTRUE(binary)) {
    area_before <- burnable_area(mask)
  } else {
    area_before <- NA_real_
  }

  # ---- overlap with template, computed on the ORIGINAL mask BEFORE alignment ----
  # We must test overlap on the input geometry: after resample/project the
  # aligned raster inherits the TEMPLATE extent, so a post-alignment extent test
  # would always "overlap". Project the mask's bounding extent into the template
  # CRS (when CRS differ) and intersect with the template extent.
  if (!terra::same.crs(mask, template)) {
    e_m <- terra::ext(terra::project(
      terra::as.polygons(terra::ext(mask), crs = terra::crs(mask)),
      terra::crs(template)
    ))
  } else {
    e_m <- terra::ext(mask)
  }
  e_t <- terra::ext(template)
  ix_xmin <- max(e_m$xmin, e_t$xmin)
  ix_xmax <- min(e_m$xmax, e_t$xmax)
  ix_ymin <- max(e_m$ymin, e_t$ymin)
  ix_ymax <- min(e_m$ymax, e_t$ymax)
  overlap <- (ix_xmax > ix_xmin) && (ix_ymax > ix_ymin)
  if (!overlap) {
    stop(sprintf(
      "`%s`: mask does not overlap the template extent. The mask and template do not share any spatial footprint (after CRS handling); check the inputs.",
      mask_name), call. = FALSE)
  }

  # =====================================================================
  # CORE ALIGNMENT (the fix): project on CRS mismatch, resample on geom
  # mismatch within same CRS, otherwise return unchanged.
  # =====================================================================
  if (!terra::same.crs(mask, template)) {
    aligned <- terra::project(mask, template, method = "near")
  } else if (!terra::compareGeom(mask, template, stopOnError = FALSE)) {
    aligned <- terra::resample(mask, template, method = "near")
  } else {
    aligned <- mask
  }

  # ---- geometry must match template after alignment ----
  geom_ok <- terra::compareGeom(aligned, template, stopOnError = FALSE)
  if (!geom_ok) {
    stop(sprintf(
      "`%s`: mask geometry still does not match the template after alignment (CRS/res/extent/dims).",
      mask_name), call. = FALSE)
  }

  # ---- allowed-values check (no silent coercion) ----
  # Use terra::unique (EXACT distinct values) rather than terra::freq, which
  # rounds floats and would let a 0.5 slip through binned into {0,1}.
  uvals <- terra::unique(aligned)
  present_vals <- if (is.null(uvals) || nrow(uvals) == 0L) numeric(0) else uvals[, 1]
  present_vals <- present_vals[!is.na(present_vals)]
  if (!is.null(allowed_values)) {
    bad <- setdiff(present_vals, allowed_values)
    allowed_values_ok <- length(bad) == 0L
    if (!allowed_values_ok) {
      stop(sprintf(
        "`%s`: mask contains value(s) outside the allowed set {%s}: {%s}. Values are NOT coerced/rounded silently.",
        mask_name,
        paste(allowed_values, collapse = ", "),
        paste(sort(unique(bad)), collapse = ", ")), call. = FALSE)
    }
  } else {
    allowed_values_ok <- NA
  }

  # ---- valid (non-NA) cell count ----
  n_valid <- {
    fr <- terra::freq(aligned)
    if (is.null(fr) || nrow(fr) == 0L) 0 else
      sum(fr[!is.na(fr[, "value"]), "count"])
  }

  if (isTRUE(binary)) {
    n_burnable <- count_value(aligned, 1)
    if (n_burnable <= 0) {
      stop(sprintf(
        "`%s`: aligned mask has ZERO burnable cells (value == 1). This is the all-NaN/zero-burnable failure the alignment helper guards against; refusing to proceed.",
        mask_name), call. = FALSE)
    }
    area_after   <- burnable_area(aligned)
    area_abs_diff <- abs(area_after - area_before)
    area_rel_diff <- if (area_before > 0) area_abs_diff / area_before else
      if (area_after > 0) 1 else 0
    within_tol <- area_rel_diff <= area_tol_rel
    if (!within_tol) {
      stop(sprintf(
        "`%s`: burnable area changed by %.2f%% during alignment (before = %.3g, after = %.3g, abs diff = %.3g), exceeding the %.1f%% tolerance. This usually indicates a CRS/geometry problem rather than ordinary grid effects.",
        mask_name, 100 * area_rel_diff, area_before, area_after,
        area_abs_diff, 100 * area_tol_rel), call. = FALSE)
    }
  } else {
    n_burnable    <- NA_real_
    area_after    <- NA_real_
    area_abs_diff <- NA_real_
    area_rel_diff <- NA_real_
    within_tol    <- NA
  }

  non_empty <- n_valid > 0
  if (!non_empty) {
    stop(sprintf(
      "`%s`: aligned mask is empty (no valid, non-NA cells).",
      mask_name), call. = FALSE)
  }

  report <- list(
    crs_ok            = terra::same.crs(aligned, template),
    res               = terra::res(aligned),
    extent            = as.vector(terra::ext(aligned)),
    dims              = c(nrow = terra::nrow(aligned), ncol = terra::ncol(aligned)),
    overlap           = overlap,
    allowed_values_ok = allowed_values_ok,
    n_valid           = n_valid,
    n_burnable        = n_burnable,
    area_before       = area_before,
    area_after        = area_after,
    area_abs_diff     = area_abs_diff,
    area_rel_diff     = area_rel_diff,
    within_tol        = within_tol,
    non_empty         = non_empty
  )

  list(aligned = aligned, report = report)
}
