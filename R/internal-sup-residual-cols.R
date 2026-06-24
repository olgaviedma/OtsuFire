# Canonical model feature whitelist + legacy residual deny-list.
#
# OtsuFire 0.4.0 architectural refactor (decided by Natalia 2026-05-09):
# the supervised model sees ONLY the columns in `.supervised_feature_cols`
# (plus their auto-generated `<col>_isNA` companions). Everything else
# (administrative metadata, fold assignments, deterministic-stage
# residuals, geometric sampling-biased features, Otsu metadata, OOF
# outputs, labels) is filtered at the LAST mile, inside
# `train_final_model_direct()` and the scoring path.
#
# Upstream stages (`extract_features()`, pool builders, orchestrator)
# preserve all input columns unchanged. The drops that lived inside
# `extract_features()` until 0.3.x have been removed entirely.
#
# History:
#   * 0.3.0 (Agent F): canonical residual deny list defined here as
#     `.deterministic_residual_cols`. Consumed by `extract_features()`
#     (input contract), `train_final_model_direct()` (extra_drop_cols
#     + forbidden detector) and `run_dm_oof_pipeline()` (drop_regex).
#   * 0.3.1 (Agent G): six load-bearing administrative columns
#     (`neg_type`, `cell_id`, `training_selected`, `training_group`,
#     `training_reason`, `intersects_deterministic`) removed from the
#     deny list because they are pipeline-internal metadata, not
#     deterministic-stage residuals.
#   * 0.4.0 (Agent H): the whitelist `.supervised_feature_cols` becomes
#     the single source of truth for what enters the model matrix.
#     `extract_features()` no longer drops anything. The legacy
#     `.deterministic_residual_cols` constant is retained as an audit
#     log of historically known residuals — it is no longer applied as
#     a runtime filter by `extract_features()`, and the matrix-build
#     site uses the whitelist semantics instead.
#   * 0.4.1 (2026-05-09): `hs_any` re-included after Natalia's review
#     decision. Initially removed in 0.4.0 because it is a derived
#     helper synthesised from `hs_used_n > 0`; reverted to maintain
#     historical parity with Phase B and the 30 production runs,
#     and for cleaner reading in feature-importance tables.
#   * 0.6.2 (2026-06-10): `hs_any` REMOVED again from the whitelist
#     (Natalia's final decision). Rationale: it is redundant with
#     `hs_used_n > 0` / `hs_in_poly > 0`, it is NEVER materialised as a
#     base GPKG feature, and it never entered the resolved FULL recipe
#     (the resolved base was already 50 because `hs_any` was never a
#     GPKG column). The whitelist declared 51 but only 50 ever
#     resolved; in NO-HOTSPOT it appeared as a fictitious 13th removed
#     feature. Its synthesis block in `internal-sup-create-matrix.R`
#     was removed. The whitelist is now 50 names; the resolved FULL
#     feature space (50 base + 50 `_isNA` = 100) and its
#     feature-schema fingerprint are UNCHANGED by this removal.
#
# NOT exported.

#' Canonical supervised-model feature whitelist (OtsuFire 0.6.2).
#'
#' The 50 columns the supervised XGBoost model sees, plus their
#' auto-generated `<col>_isNA` missingness companions when the matrix
#' builder synthesises them. Decided by Natalia on 2026-05-09; `hs_any`
#' removed in 0.6.2 (see History above). This list is FIXED — additions
#' or removals are an explicit methodological decision and require a
#' major-version bump and Natalia's approval.
#'
#' Notes:
#'   * `block_id`, `fold_rep1`, `fold_rep2` are NOT features — they are
#'     fold-assignment metadata consumed by the OOF wrapper. They must
#'     survive in the GPKG but must NOT enter the model matrix.
#'   * No ecoregion features (`eco_*`) in this build. Ecoregions were
#'     removed entirely from the supervised phase on 2026-06-05; they
#'     belong to the deterministic delineation stage only.
#'   * Hotspot block (12 columns): `hotspot_available` is a real
#'     structural flag (0 = no MODIS hotspot data for this year, 1 =
#'     available) and is ALWAYS observed (never NA). When the year has
#'     no hotspot data (`hotspot_available == 0`), the features
#'     extractor's `out_nodata()` emits NA (not -9999, fixed in B5,
#'     2026-06-06) for every MEASURED hotspot quantity
#'     (`hs_in_poly`, `hs_in_buffer`, `hs_used_n`, `hs_min_dist_m`,
#'     `hs_frp_sum`, `hs_frp_max`, `hs_conf_mean`, `hs_hiConf_n`) and
#'     for the support-existence flags `hs_support_present` and
#'     `hs_only_buffer_support` (undefined without data). Those NAs are
#'     turned into `<col>_isNA = 1` companions by the matrix builder
#'     and the value is imputed to the column median, so "no hotspots
#'     data available" is encoded as a real signal and a model trained
#'     on hotspot-years follows its learned missing direction when
#'     applied to a hotspot-less year. `hs_no_support_when_available`
#'     stays 0 when `hotspot_available == 0` (its semantics are
#'     conditioned on availability), so its `_isNA` companion is
#'     degenerate and the builder skips it.
#'   * `hs_any` (a binary `hs_used_n > 0` helper) was REMOVED from this
#'     whitelist in 0.6.2 as a redundant, never-materialised helper
#'     (equivalent to `hs_used_n > 0` / `hs_in_poly > 0`), per Natalia's
#'     decision. It was never a base GPKG feature and never entered the
#'     resolved recipe, so the resolved feature space is unchanged.
#'
#' @keywords internal
#' @noRd
.supervised_feature_cols <- c(
  # RBR same-window (5)
  "rbr_valid_frac", "rbr_p10", "rbr_med", "rbr_p90", "rbr_iqr",

  # CORINE land cover (7)
  "cor_open_frac", "cor_wetlands_frac", "cor_water_frac",
  "cor_herbaceous_frac", "cor_urban_frac", "cor_agri_frac",
  "cor_forest_frac",

  # Topography - elevation (7)
  "elev_valid_frac", "elev_p10", "elev_med", "elev_p90",
  "elev_iqr", "elev_mean", "elev_sd",

  # Topography - slope (7)
  "slope_valid_frac", "slope_p10", "slope_med", "slope_p90",
  "slope_iqr", "slope_mean", "slope_sd",

  # DOY (5)
  "doy_valid_frac", "doy_p10", "doy_med", "doy_p90", "doy_iqr",

  # RBR all-window + persistence (7)
  "rbr_aw_valid_frac", "rbr_aw_p10", "rbr_aw_med", "rbr_aw_p90",
  "rbr_aw_iqr", "persist_delta", "persist_ratio",

  # Hotspots (12)
  "hotspot_available", "hs_in_poly", "hs_in_buffer", "hs_used_n",
  "hs_min_dist_m", "hs_frp_sum", "hs_frp_max", "hs_conf_mean",
  "hs_hiConf_n", "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support"
)

#' OPTIONAL shape/size feature block (OtsuFire 0.12.0, OFF by default).
#'
#' The six geometric columns the supervised model MAY additionally see
#' when the experimental `include_shape_features` flag is TRUE. This is
#' a SEPARATE constant from the frozen canonical 50-name
#' `.supervised_feature_cols` (which must NOT grow without a major
#' bump). `area_ha` and `n_pix` ride along as pool-builder columns; the
#' remaining four (`log_area`, `perim_m`, `compactness`, `elongation`)
#' are computed by `shape_features()` only when the flag is ON.
#'
#' SAMPLING-BIAS CAVEAT: the random-background negatives are tiny fixed
#' cells (~0.81 ha squares), so shape/area is partly a sampling
#' artifact, not a physical signal. This block is OFF by default and is
#' intended for experimentation only; judge its effect with EFFIS, not
#' OOF. See NEWS.md 0.12.0.
#'
#' @keywords internal
#' @noRd
.supervised_shape_feature_cols <- c(
  "area_ha", "n_pix", "log_area", "perim_m", "compactness", "elongation"
)

#' THE admissible supervised-feature universe (single source of truth).
#'
#' Returns the canonical 50-name whitelist, optionally extended with the
#' six `.supervised_shape_feature_cols` when `include_shape = TRUE`. This
#' ONE helper is the definition of the admissible feature set everywhere
#' (validation, active-whitelist default, the OOF/FINAL allowed-column
#' invariants). Routing every site through it eliminates the two-filter
#' asymmetry that previously let a column appear "in OOF but not FINAL".
#'
#' @param include_shape Logical. When TRUE the six shape columns are
#'   appended to the canonical list. Default FALSE -> byte-identical to
#'   `.supervised_feature_cols`.
#' @keywords internal
#' @noRd
.supervised_feature_universe <- function(include_shape = FALSE) {
  c(.supervised_feature_cols,
    if (isTRUE(include_shape)) .supervised_shape_feature_cols)
}

# Min/max side length of a (rotated-rectangle) polygon. Package-level
# twin of the deterministic stage's `get_wd_ln_from_mrr()` closure
# (internal-det-polygon-metrics.R), lifted here so the supervised shape
# helper can reuse the SAME MRR length/width math without depending on
# the deterministic function's closure environment. Returns c(wd, ln).
#
# @keywords internal
# @noRd
.of_shape_wd_ln <- function(poly) {
  tryCatch({
    ls <- sf::st_cast(poly, "LINESTRING", warn = FALSE)
    coords <- sf::st_coordinates(ls)[, 1:2, drop = FALSE]
    if (nrow(coords) < 2) return(c(wd = NA_real_, ln = NA_real_))
    if (all(coords[1, ] == coords[nrow(coords), ]))
      coords <- coords[-nrow(coords), , drop = FALSE]
    if (nrow(coords) < 2) return(c(wd = NA_real_, ln = NA_real_))
    prev <- rbind(coords[nrow(coords), , drop = FALSE],
                  coords[-nrow(coords), , drop = FALSE])
    segs <- sqrt(rowSums((coords - prev) ^ 2))
    segs <- segs[is.finite(segs) & segs > 0]
    if (length(segs) == 0) return(c(wd = NA_real_, ln = NA_real_))
    c(wd = min(segs, na.rm = TRUE), ln = max(segs, na.rm = TRUE))
  }, error = function(e) c(wd = NA_real_, ln = NA_real_))
}

# Per-polygon shape/size math for the OPTIONAL shape block (OtsuFire
# 0.12.0). Computes the four shape columns the model MAY use when
# include_shape_features is ON, keyed by `id_col`:
#   * log_area    = log1p(area_ha), area_ha = st_area / 1e4
#   * perim_m     = perimeter via the MULTILINESTRING-length idiom
#                   (utils-polygon-metrics.R:19-21)
#   * compactness = Polsby-Popper 4*pi*A / P^2
#                   (utils-polygon-metrics.R:22-24)
#   * elongation  = min-rotated-rectangle length / width (>= 1), reusing
#                   the deterministic stage's MRR side-length extraction
#                   (.of_shape_wd_ln), with a base sf::st_bbox aspect-
#                   ratio fallback for degenerate geometries.
# `polys` must be a NON-empty sf in a metric CRS (the supervised stack is
# EPSG:3035). area_ha / n_pix already ride along as pool-builder columns,
# so only these four are produced here. Uses base sf only (no lwgeom).
#
# @keywords internal
# @noRd
.of_shape_features <- function(polys, id_col = "fire_uid") {
  stopifnot(inherits(polys, "sf"), nrow(polys) > 0L)
  geom <- sf::st_geometry(polys)

  area_m2 <- as.numeric(sf::st_area(polys))
  area_ha <- area_m2 / 1e4
  log_area <- log1p(area_ha)

  perim_m <- vapply(
    geom,
    function(gi) as.numeric(
      sf::st_length(sf::st_cast(gi, "MULTILINESTRING", warn = FALSE))
    ),
    numeric(1)
  )

  compactness <- ifelse(perim_m > 0,
                        4 * pi * area_m2 / (perim_m^2),
                        NA_real_)

  bbox_ratio <- function(gi) {
    bb <- tryCatch(sf::st_bbox(gi), error = function(e) NULL)
    if (is.null(bb)) return(NA_real_)
    dx <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
    dy <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
    lo <- min(dx, dy); hi <- max(dx, dy)
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0) return(NA_real_)
    hi / lo
  }
  mrr_elong <- function(gi) {
    tryCatch({
      mrr <- sf::st_minimum_rotated_rectangle(gi)
      wl  <- .of_shape_wd_ln(sf::st_geometry(mrr)[[1]])
      wd  <- wl[["wd"]]; ln <- wl[["ln"]]
      if (!is.finite(wd) || !is.finite(ln) || wd <= 0) return(bbox_ratio(gi))
      ln / wd
    }, error = function(e) bbox_ratio(gi))
  }
  elongation <- vapply(geom, mrr_elong, numeric(1))

  out <- data.frame(
    .id_tmp     = as.character(polys[[id_col]]),
    log_area    = as.numeric(log_area),
    perim_m     = as.numeric(perim_m),
    compactness = as.numeric(compactness),
    elongation  = as.numeric(elongation),
    stringsAsFactors = FALSE
  )
  names(out)[1L] <- id_col
  out
}

# Legacy deny list — RETAINED as an audit log of historically known
# deterministic-stage residuals and sampling-biased columns. Under
# 0.4.0 this list is NOT applied as a runtime filter by
# `extract_features()` (which has been refactored to be ADD-only).
# `train_final_model_direct()` no longer consults this list either;
# its filter operates on the whitelist semantics
# (`.supervised_feature_cols`). The constant is kept exported within
# the package namespace because a small number of legacy regression
# tests still reference it; it can be deleted in a future major
# release once those tests are migrated.

#' @keywords internal
#' @noRd
.deterministic_residual_cols <- c(
  # Filter / decisions
  "filter_1", "reason_1", "filter_2", "reason_2", "filter_3", "reason_3",
  "preyear_action", "preyear_reason", "preyear_overlap_frac",
  "class_final", "raw_class", "class_audited",
  "otsu_decision",
  # Geometric features deliberately excluded (sampling-induced size
  # bias documented in NEWS.md 0.3.1).
  "n_pix", "area_ha", "log_area", "perim_m",
  "compactness", "elongation", "n_holes",
  # Spectral residuals from deterministic Otsu/grow + scoring stages
  "median_rbr",
  "p_above_keep_q25", "p_above_keep_ref", "p_above_keep_q05",
  "p_above_keep_q50", "p_above_keep_q75", "p_above_keep_q95",
  "percentile_in_keep",
  "conf_area", "conf_pix",
  "flag_rbr",
  # QA + OOF audit columns
  "qa_changed", "qa_cert", "qa_uncert", "qa_final", "qa_action", "qa_contra",
  "p_oof_mean", "p_oof_med", "p_oof_sd", "n_preds",
  "p_oof_min_long", "p_oof_max_long", "n_oof_long",
  "class_oof", "class_mismatch",
  # Legacy GDAL/polygonize attributes that survive sf reads
  "DN", "N_PIX_PIX", "AREA_HA_HA", "PATCH_ID_I", "N_TOTAL_TO",
  "N_REF_100_", "REF_100_10", "CV_BASE_BA",
  "NUCLEO", "VECINDAD", "DIST_M_M", "W_DIST_DIS",
  "BOOST", "S_PATCH_PA", "T_STAR_STA", "BUF_M_M",
  "THR_CORE_C", "ALPHA", "MB_MIN_MIN", "DPWR",
  "KEEPHI", "DROP_LO_LO", "KEEP_P_P", "DECISION",
  "CVPCT", "WKPCT"
)

# Helper: given a data.frame / sf names vector, return the subset of
# names that are admissible model features under the 0.4.0 whitelist
# semantics (a feature in the active `whitelist`, or its
# auto-generated `<feat>_isNA` companion). Used by
# `train_final_model_direct()`, `run_dm_oof_pipeline()`, and by tests
# that assert the matrix builder's column set is a subset of this.
#
# 0.5.0 (2026-05-09): the `whitelist` argument lets callers pass an
# active feature subset (e.g. `feature_whitelist_override`) while
# defaulting to the canonical `.supervised_feature_cols`.
#
# @keywords internal
# @noRd
.filter_to_supervised_whitelist <- function(nms, whitelist = .supervised_feature_cols) {
  feature_base <- whitelist
  isNA_companions <- paste0(feature_base, "_isNA")
  allowed <- c(feature_base, isNA_companions)
  intersect(nms, allowed)
}
