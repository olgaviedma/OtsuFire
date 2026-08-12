# Partial-observability audit for validate_fire_maps().
#
# INFORMATIVE ONLY. This helper summarizes how many reference fires were
# temporally excluded, and — within the fires the contract KEEPS WHOLE — how
# many burnable-domain pixels carry no observation (partial spatial
# observability). It never changes the whole-fire temporal observability rule,
# the reference filtering, or any TP/FP/FN/TN/coverage figure.
#
# Input: the reference_observability table built by build_reference_observability
# (one row per ORIGINAL reference fire), which must carry observable_flag,
# observable_reason, ref_area_domain_ha and the per-fire partial columns
# total_pixels / observable_pixels / non_observable_pixels / observable_fraction.
#
# Returns a list with three data frames:
#   summary  -- one row of audit counts
#   excluded -- the temporally non-observable fires (dropped before metrics)
#   partial  -- the KEPT fires that are only partially observable (kept whole)

.vfm_observability_audit <- function(reference_observability, cell_area_ha,
                                     id_candidates = c("fire_id", "id", "ID",
                                                       "FID", "OBJECTID")) {
  ro <- as.data.frame(reference_observability, stringsAsFactors = FALSE)
  n_original <- nrow(ro)

  id_col <- intersect(id_candidates, names(ro))
  id_col <- if (length(id_col)) id_col[[1]] else NA_character_
  fire_id <- if (!is.na(id_col)) ro[[id_col]] else ro$reference_row

  flag <- as.logical(ro$observable_flag)
  flag[is.na(flag)] <- FALSE
  excl <- !flag

  area <- if ("ref_area_domain_ha" %in% names(ro)) ro$ref_area_domain_ha else
    rep(NA_real_, n_original)

  has_part <- all(c("total_pixels", "observable_pixels",
                    "non_observable_pixels", "observable_fraction") %in% names(ro))

  if (has_part) {
    partial_mask <- flag &
      !is.na(ro$non_observable_pixels) & ro$non_observable_pixels > 0 &
      !is.na(ro$observable_fraction) & ro$observable_fraction < 1
    partial_non_observable_pixels <-
      sum(ro$non_observable_pixels[flag], na.rm = TRUE)
  } else {
    partial_mask <- rep(FALSE, n_original)
    partial_non_observable_pixels <- 0
  }

  summary <- data.frame(
    n_reference_fires_original             = n_original,
    n_reference_fires_observable           = sum(flag),
    n_reference_fires_excluded             = sum(excl),
    n_reference_fires_partially_observable = sum(partial_mask),
    excluded_reference_area_ha             = round(sum(area[excl], na.rm = TRUE), 4),
    partial_non_observable_pixels          = as.integer(round(partial_non_observable_pixels)),
    partial_non_observable_area_ha         = round(partial_non_observable_pixels * cell_area_ha, 4),
    stringsAsFactors = FALSE
  )

  excluded <- data.frame(
    fire_id                     = fire_id[excl],
    reference_row               = ro$reference_row[excl],
    burnable_area_ha            = area[excl],
    total_pixels                = if (has_part) ro$total_pixels[excl] else NA_integer_,
    observable_pixels           = if (has_part) ro$observable_pixels[excl] else NA_integer_,
    observable_fraction         = if (has_part) ro$observable_fraction[excl] else NA_real_,
    exclusion_reason            = ro$observable_reason[excl],
    confirmation_not_in_metrics = rep(
      "excluded before rasterization: 0 TP, 0 FN, 0 evaluated area", sum(excl)),
    stringsAsFactors            = FALSE
  )

  partial <- data.frame(
    fire_id                   = fire_id[partial_mask],
    reference_row             = ro$reference_row[partial_mask],
    burnable_area_ha          = area[partial_mask],
    total_pixels              = if (has_part) ro$total_pixels[partial_mask] else integer(0),
    temporally_valid_pixels   = if (has_part) ro$observable_pixels[partial_mask] else integer(0),
    non_observable_pixels     = if (has_part) ro$non_observable_pixels[partial_mask] else integer(0),
    non_observable_area_ha    = if (has_part)
      round(ro$non_observable_pixels[partial_mask] * cell_area_ha, 4) else numeric(0),
    observable_fraction       = if (has_part) ro$observable_fraction[partial_mask] else numeric(0),
    contract_keeps_whole_fire = if (sum(partial_mask)) TRUE else logical(0),
    stringsAsFactors          = FALSE
  )

  list(summary = summary, excluded = excluded, partial = partial)
}
