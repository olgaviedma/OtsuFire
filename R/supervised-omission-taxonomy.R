# =============================================================================
# supervised-omission-taxonomy.R
#
# ONE source of truth for the omission/commission SOURCE taxonomy across the
# whole thesis (deterministic per-preset, supervised unsup-vs-sup figures).
# Packaged from the standalone `omission_cascade_lib.R` (2026-06-24, Phase 1).
#
# CORRECTED, CANDIDATE-FIRST MECHANISM (do not regress):
#   * POPULATION of "omitted" fires = the VALIDATION omission layer
#     (04_ERROR_LAYERS/reference_fires_omission_*.gpkg) of each phase, i.e. the
#     SAME set that defines EFFIS Recall. NEVER an ad-hoc geometric overlap cut,
#     and NEVER pre-filtered by issue_type.
#   * CANDIDATE PRESENCE is decided FIRST. A fire "has a final candidate" if it
#     overlaps ANY deterministic candidate (keep/review/drop). Only when there is
#     NO candidate is the physical RBR/burnability reason consulted. This keeps
#     recoverable cases (review/drop) from being mixed into the "below detection"
#     bucket.
#   * Every NO-CANDIDATE subclass depends ONLY on (burnable, rbr, seed) and is
#     therefore invariant between the unsupervised and supervised models; only
#     the has-candidate buckets may shrink (the model recovering fires).
#
# This 7-class candidate-first taxonomy SUPERSEDES the legacy 5-class taxonomy
# (`ERR_LEVELS`, kept only in the deprecated standalone characterization scripts
# and NOT shipped in this package). The package main path uses only the classes
# defined here.
# =============================================================================

# ---- thresholds -------------------------------------------------------------

#' Omission/commission source-taxonomy thresholds and class levels
#'
#' Fixed thresholds and the class-level/colour vectors that define the
#' candidate-first omission taxonomy (7 classes) and the commission taxonomy
#' (4 classes). They are the values used by
#' \code{\link{of_classify}} and \code{\link{of_classify_commission}}.
#'
#' \describe{
#'   \item{\code{OF_WEAK}}{RBR floor (175); below it the reference burn signal is
#'     considered too weak to be a real omission.}
#'   \item{\code{OF_SLIVER}}{Minimum candidate overlap (ha); below it a fire is
#'     treated as having NO final candidate.}
#'   \item{\code{OF_SEED_THR}}{Per-preset Otsu seed thresholds (named numeric).}
#'   \item{\code{OF_LEVELS}/\code{OF_COLORS}}{7-class omission levels + colours.}
#'   \item{\code{OF_COM_LEVELS}/\code{OF_COM_COLORS}}{4-class commission levels +
#'     colours.}
#' }
#'
#' @format Numeric scalars / named numeric / character / named-character vectors.
#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_WEAK <- 175

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_SLIVER <- 0.01   # ha; below this a candidate overlap is treated as none

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_SEED_THR <- c(very_permissive = 285, permissive = 300, balanced = 310,
                 conservative = 400, very_conservative = 430)

# ---- 7-class candidate-first omission taxonomy ------------------------------
# Stacking / factor order: recoverable & decision buckets on top, then the
# method-limit miss, then the physical/reference floors.

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_LEVELS <- c(
  "Candidate in review (recoverable)",
  "Candidate dropped (rejected as unburned)",
  "Candidate kept (validation-threshold artifact)",
  "No candidate: algorithm failure (RBR>=seed)",
  "No candidate: below detection (RBR<seed)",
  "No candidate: reference too weak (RBR<175)",
  "No candidate: non-burnable (masked)")

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_COLORS <- c(
  "Candidate in review (recoverable)"              = "#E6A700",  # amber  (review)
  "Candidate dropped (rejected as unburned)"       = "#B5524A",  # brick  (drop)
  "Candidate kept (validation-threshold artifact)" = "#2F7D58",  # green  (keep)
  "No candidate: algorithm failure (RBR>=seed)"    = "#D7301F",  # red    (real miss)
  "No candidate: below detection (RBR<seed)"       = "#9970AB",  # purple
  "No candidate: reference too weak (RBR<175)"     = "#969696",  # grey
  "No candidate: non-burnable (masked)"            = "#3182BD")  # blue

# ---- 4-class commission taxonomy --------------------------------------------
# Every commission polygon IS a kept candidate, so the cascade is purely
# physical (burnability + RBR vs WEAK floor and per-preset seed); no
# candidate-first branch applies.

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_COM_LEVELS <- c(
  "Likely real burn (absent from reference)",  # RBR >= seed
  "Borderline signal (RBR < seed)",            # WEAK <= RBR < seed
  "Weak spectral signal",                      # RBR < WEAK
  "Non-burnable")                              # CORINE not burnable

#' @rdname omission_taxonomy_constants
#' @keywords internal
OF_COM_COLORS <- c(
  "Likely real burn (absent from reference)" = "#3A8459",  # green
  "Borderline signal (RBR < seed)"           = "#9970AB",  # purple
  "Weak spectral signal"                     = "#969696",  # grey
  "Non-burnable"                             = "#3182BD")  # blue

# CORINE reclass burnability (codes 1-7 burnable, 8-11 not). Internal detail.
OF_CLC <- data.frame(code = 1:11, burnable = c(rep(TRUE, 7), rep(FALSE, 4)))

# ---- classifiers ------------------------------------------------------------

#' Classify omission sources (candidate-first 7-class taxonomy)
#'
#' Classifies the source of each omitted fire. Candidate presence is
#' evaluated first: when the total candidate overlap is below \code{OF_SLIVER}
#' the fire has no final candidate and is classified by the physical
#' RBR/burnability cascade; otherwise it takes the dominant of its
#' keep/review/drop footprints. Every no-candidate class depends only on
#' \code{(burnable, rbr, seed)} and is invariant between models.
#'
#' @param keep_ha,review_ha,drop_ha Numeric vectors. Per-fire overlap area (ha)
#'   with the hierarchical keep/review/drop candidate footprints (see
#'   \code{\link{of_decompose}}).
#' @param candidate_ha Numeric. Total candidate overlap (\code{keep+review+drop}).
#' @param burnable Logical. Whether the fire footprint is on burnable CORINE.
#' @param rbr Numeric. Per-fire RBR statistic (e.g. median or p90).
#' @param seed Numeric. Per-preset Otsu seed; length-1 is recycled.
#' @return An ordered \code{factor} with levels \code{OF_LEVELS}.
#' @seealso \code{\link{of_classify_commission}}, \code{\link{of_decompose}}
#' @keywords internal
of_classify <- function(keep_ha, review_ha, drop_ha, candidate_ha, burnable, rbr, seed) {
  n <- length(candidate_ha)
  seed <- rep_len(seed, n)
  lab <- character(n)
  for (i in seq_len(n)) {
    if (is.na(candidate_ha[i]) || candidate_ha[i] < OF_SLIVER) {
      # NO final candidate -> physical RBR/burnability cascade
      r <- rbr[i]; b <- burnable[i]
      lab[i] <- if (!is.na(b) && !b) OF_LEVELS[7]
                else if (!is.na(r) && r < OF_WEAK) OF_LEVELS[6]
                else if (!is.na(r) && r < seed[i]) OF_LEVELS[5]
                else OF_LEVELS[4]
    } else {
      dom <- which.max(c(keep_ha[i], review_ha[i], drop_ha[i]))  # 1=keep,2=review,3=drop
      lab[i] <- c(OF_LEVELS[3], OF_LEVELS[1], OF_LEVELS[2])[dom]
    }
  }
  factor(lab, levels = OF_LEVELS)
}

#' Classify commission sources (4-class physical taxonomy)
#'
#' Commission polygons are map detections absent from the reference. Every such
#' polygon is a kept candidate, so the cascade is purely physical: non-burnable,
#' then RBR vs the \code{OF_WEAK} floor and the per-preset seed.
#'
#' @param burnable Logical. Whether the polygon is on burnable CORINE.
#' @param rbr Numeric. Per-polygon RBR statistic.
#' @param seed Numeric. Per-preset Otsu seed; length-1 is recycled.
#' @return An ordered \code{factor} with levels \code{OF_COM_LEVELS}.
#' @seealso \code{\link{of_classify}}
#' @keywords internal
of_classify_commission <- function(burnable, rbr, seed) {
  n <- length(rbr)
  seed <- rep_len(seed, n)
  lab <- character(n)
  for (i in seq_len(n)) {
    b <- burnable[i]; r <- rbr[i]
    lab[i] <- if (!is.na(b) && !b) OF_COM_LEVELS[4]
              else if (!is.na(r) && r < OF_WEAK) OF_COM_LEVELS[3]
              else if (!is.na(r) && r < seed[i]) OF_COM_LEVELS[2]
              else OF_COM_LEVELS[1]
  }
  factor(lab, levels = OF_COM_LEVELS)
}

# ---- geometry / extraction helpers ------------------------------------------

#' Make valid, drop empties and project to ETRS89-LAEA (EPSG:3035)
#' @param x An \code{sf} object.
#' @return The cleaned \code{sf} in EPSG:3035.
#' @keywords internal
of_g3 <- function(x) {
  x <- suppressWarnings(sf::st_make_valid(x))
  x <- x[!sf::st_is_empty(x), ]
  if (is.na(sf::st_crs(x)$epsg) || sf::st_crs(x)$epsg != 3035) x <- sf::st_transform(x, 3035)
  x
}

# numeric coercion of an exact_extract result (internal)
.of_vc <- function(x) { if (is.data.frame(x)) x <- x[[1]]; as.numeric(x) }

#' Per-polygon RBR statistics (median and p90)
#' @param polys An \code{sf} of polygons.
#' @param rbr_rast A \code{SpatRaster} of RBR (first layer used if multi-layer).
#' @return A data.frame with columns \code{median} and \code{p90}.
#' @keywords internal
of_rbr_stats <- function(polys, rbr_rast) {
  rr <- rbr_rast; if (terra::nlyr(rr) > 1) rr <- rr[[1]]
  pp <- sf::st_transform(polys, terra::crs(rr))
  data.frame(
    median = .of_vc(exactextractr::exact_extract(rr, pp, "median", progress = FALSE)),
    p90    = .of_vc(exactextractr::exact_extract(rr, pp, "quantile", quantiles = 0.9, progress = FALSE)))
}

#' Per-polygon burnability from a CORINE-coded raster
#' @param polys An \code{sf} of polygons.
#' @param corine_rast A \code{SpatRaster} of CORINE reclass codes (1-11).
#' @return A logical vector (TRUE = burnable) aligned to \code{polys}.
#' @keywords internal
of_burnable <- function(polys, corine_rast) {
  rc <- corine_rast; if (terra::nlyr(rc) > 1) rc <- rc[[1]]
  cc <- as.integer(round(.of_vc(exactextractr::exact_extract(
        rc, sf::st_transform(polys, terra::crs(rc)), "mode", progress = FALSE))))
  OF_CLC$burnable[match(cc, OF_CLC$code)]
}

#' Decompose omitted fires into hierarchical keep>review>drop candidate areas
#'
#' Builds hierarchical (non-overlapping) keep > review > drop footprints from an
#' \code{internal_decisions} sf, then computes per-omitted-fire overlap area by
#' class. Used when no precomputed diagnosis layer exists (e.g. the supervised
#' omission layer). Toggles s2 off locally and restores it on exit.
#'
#' @param omitted An \code{sf} of omitted reference fires.
#' @param internal An \code{sf} of internal decisions carrying a class field.
#' @param class_field Name of the class column (default \code{"class_final"}).
#' @return \code{omitted} with added \code{keep_ha}, \code{review_ha},
#'   \code{drop_ha}, \code{candidate_ha} columns.
#' @seealso \code{\link{of_classify}}
#' @keywords internal
of_decompose <- function(omitted, internal, class_field = "class_final") {
  prev_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(prev_s2)), add = TRUE)

  internal$.cls <- tolower(trimws(as.character(internal[[class_field]])))
  internal <- internal[internal$.cls %in% c("keep", "review", "drop"), ]
  uni <- function(s) { if (!nrow(s)) return(NULL); g <- suppressWarnings(sf::st_union(sf::st_geometry(s))); sf::st_make_valid(g) }
  dif <- function(a, b) { if (is.null(a)) return(NULL); if (is.null(b)) return(a); sf::st_make_valid(suppressWarnings(sf::st_difference(a, b))) }
  k <- uni(internal[internal$.cls == "keep", ])
  r0 <- uni(internal[internal$.cls == "review", ])
  d0 <- uni(internal[internal$.cls == "drop", ])
  r <- dif(r0, k)
  kr <- if (is.null(k)) r else if (is.null(r)) k else sf::st_union(c(k, r))
  d <- dif(d0, kr)
  # robust per-feature overlap area (st_intersection idx handling varies by version)
  ov <- function(g) { if (is.null(g)) return(rep(0, nrow(omitted)))
    vapply(seq_len(nrow(omitted)), function(i) {
      x <- suppressWarnings(sf::st_intersection(sf::st_geometry(omitted)[i], g))
      if (!length(x)) 0 else sum(as.numeric(sf::st_area(x)) / 1e4)
    }, numeric(1)) }
  omitted$keep_ha <- ov(k); omitted$review_ha <- ov(r); omitted$drop_ha <- ov(d)
  omitted$candidate_ha <- omitted$keep_ha + omitted$review_ha + omitted$drop_ha
  omitted
}
