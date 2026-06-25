#' Optional chain/bridge cleaning of refined deterministic candidates
#'
#' @description
#' \code{apply_chain_cleaning()} (alias \code{apply_chain()}) is an \strong{optional}
#' stage of the deterministic workflow that removes "chaining" artefacts: thin
#' bridges/necks produced by region growing that connect genuinely separate burned
#' cores (and bright-soil false positives) into a single fractal blob. It operates
#' on the \strong{already refined} candidates written by the refine stage
#' (\code{02_REFINE/BA_<year>_REFINE_MERGED_*.gpkg}) and writes a new, separate
#' candidate layer under \code{CHAIN/dechain_<N>m/}, ready to be consumed by the
#' scoring/decision stages WITHOUT re-running grow or refine.
#'
#' @details
#' \strong{Where it sits in the deterministic flow:}
#' \preformatted{
#'   01_GROW -> 02_REFINE -> [apply_chain_cleaning()] -> 03_SCORE_PHASE1 ->
#'   03_SCORE_PHASE2 -> 05_DECISIONS
#' }
#' It is placed AFTER refine (not inside grow) so that the refine \code{merge_overlaps}
#' step cannot re-bridge cut necks, and so existing \code{01_GROW}/\code{02_REFINE}
#' outputs can be reused. Scores and keep/review/drop decisions are therefore computed
#' on the already-cleaned geometry.
#'
#' \strong{Method (selective, do-no-harm).} Only candidates flagged as chained are
#' touched; everything else passes through byte-for-byte. A candidate is chained when
#' its Polsby-Popper inflation \code{1/sqrt(compact_pp)} (with
#' \code{compact_pp = 4*pi*A/P^2}, the package convention) exceeds \code{inflation_thr}
#' (default 22.135, the geometric-sanity-gate threshold), optionally also requiring
#' \code{area_ha >= min_area_ha}. Flagged candidates are rasterised on a
#' \code{pixel_size_m} grid and a morphological \strong{opening} (erosion then
#' dilation with a disk of radius \code{opening_px}) severs necks narrower than the
#' opening diameter; the result is split into connected components and re-polygonised.
#' Untouched candidates are NOT rasterised (no resampling loss).
#'
#' \strong{Distance to pixels.} \code{opening_px = dechain_distance_m / pixel_size_m}.
#' With \code{pixel_size_m = 90}: 90 m = 1 px, 180 m = 2 px, 270 m = 3 px. If the ratio
#' is not an integer a warning is emitted and it is rounded to the nearest integer
#' (\code{round()}); a value below 1 is an error.
#'
#' \strong{Markers.} \code{marker = "connected"} (default, validated) keeps every
#' connected component after opening. \code{marker = "seeds"} additionally anchors
#' pieces to the grow seeds and requires \code{change_index}, \code{vegetation_map} and
#' \code{seed_threshold_by_vegetation}; pieces with no seed are dropped when
#' \code{keep_noseed_fragments = FALSE}. If "seeds" is requested without those inputs
#' the function warns and falls back to "connected".
#'
#' \strong{This function:} is opt-in and does nothing unless you call it; does NOT
#' change the standard pipeline; does NOT re-run grow/refine; does NOT train models;
#' does NOT touch the supervised phase; and never overwrites the original
#' \code{REFINE_MERGED}.
#'
#' \strong{Caveat.} Too aggressive a \code{dechain_distance_m} (or a non-selective
#' threshold) can fragment real fires or erode small cores and raise omission. Always
#' validate the output visually (e.g. in QGIS) before using it downstream.
#'
#' @param refine_path Character. Path to \code{02_REFINE/BA_<year>_REFINE_MERGED_*.gpkg}.
#' @param out_dir Character or NULL. Deterministic run directory under which
#'   \code{CHAIN/dechain_<N>m/} is created. If NULL, inferred as the parent of the
#'   \code{02_REFINE} folder (\code{dirname(dirname(refine_path))}).
#' @param year Integer or NULL. Target year; if NULL it is parsed from the file name.
#' @param dechain_distance_m Numeric. Cleaning intensity in metres (neck width to cut).
#'   Default 180. Typical: 90 / 180 / 270.
#' @param pixel_size_m Numeric. Pixel size used for rasterisation. Default 90.
#' @param inflation_thr Numeric. Only candidates with \code{1/sqrt(compact_pp) >}
#'   this value are de-chained. Default 22.135.
#' @param min_area_ha Numeric or NULL. Optional extra requirement: only de-chain
#'   flagged candidates with \code{area_ha >= min_area_ha}. Default NULL (no minimum).
#' @param marker Character. "connected" (default) or "seeds" (see Details).
#' @param keep_noseed_fragments Logical. With \code{marker = "seeds"}, keep (TRUE) or
#'   drop (FALSE) split pieces that contain no seed. Default TRUE.
#' @param change_index,vegetation_map Character or NULL. Rasters needed only when
#'   \code{marker = "seeds"} (RBR change index and CORINE vegetation classes).
#' @param seed_threshold_by_vegetation Named numeric or NULL. Per-vegetation seed
#'   thresholds (the same used by the grow), needed only for \code{marker = "seeds"}.
#' @param overwrite Logical. Overwrite an existing \code{CHAIN/dechain_<N>m/} output.
#'   Default FALSE.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @return A list with: \code{output_path}, \code{audit_path}, \code{params_path},
#'   \code{n_candidates_before}, \code{n_candidates_after}, \code{n_chained_flagged},
#'   \code{area_before_ha}, \code{area_after_ha}, \code{area_removed_ha},
#'   \code{opening_px}, \code{dechain_distance_m}.
#'
#' @section Continuing to scoring:
#' Feed \code{output_path} to the scoring stage as the burned candidates, e.g.
#' \code{score_burned_patches(burned_candidates = res$output_path, config = cfg, ...)},
#' so \code{03_SCORE_PHASE1}, \code{03_SCORE_PHASE2} and \code{05_DECISIONS} are computed
#' on the cleaned geometry. Grow and refine are not repeated.
#'
#' @examples
#' \dontrun{
#' refine <- file.path("1_DATA/Results/1985/DETERMINISTIC/1985_permissive_omission_push",
#'                     "02_REFINE/BA_1985_REFINE_MERGED_otsu285_d115_seed12_whitebox.gpkg")
#' # cut necks up to 180 m wide on chained blobs only:
#' res <- apply_chain_cleaning(refine, dechain_distance_m = 180)
#' res$output_path   # -> .../CHAIN/dechain_180m/BA_1985_REFINE_MERGED_dechain180m.gpkg
#' # compare a gentler level:
#' res90 <- apply_chain_cleaning(refine, dechain_distance_m = 90)
#' }
#' @export
apply_chain_cleaning <- function(refine_path,
                                 out_dir = NULL,
                                 year = NULL,
                                 dechain_distance_m = 180,
                                 pixel_size_m = 90,
                                 inflation_thr = 22.135,
                                 min_area_ha = NULL,
                                 marker = c("connected", "seeds"),
                                 keep_noseed_fragments = TRUE,
                                 change_index = NULL,
                                 vegetation_map = NULL,
                                 seed_threshold_by_vegetation = NULL,
                                 overwrite = FALSE,
                                 verbose = TRUE) {
  stopifnot(is.character(refine_path), length(refine_path) == 1L, file.exists(refine_path))
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE))
    stop("apply_chain_cleaning() requires the 'sf' and 'terra' packages.")
  marker <- match.arg(marker)
  say <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## --- distance -> pixels (explicit, documented rounding) ---
  if (!is.numeric(pixel_size_m) || pixel_size_m <= 0) stop("pixel_size_m must be a positive number.")
  raw_px <- dechain_distance_m / pixel_size_m
  opening_px <- as.integer(round(raw_px))
  if (abs(raw_px - opening_px) > 1e-9)
    warning(sprintf("dechain_distance_m/pixel_size_m = %.4f is not an integer; rounding to %d px.",
                    raw_px, opening_px))
  if (opening_px < 1L)
    stop(sprintf("opening_px = %d (< 1). dechain_distance_m (%g) must be >= pixel_size_m (%g).",
                 opening_px, dechain_distance_m, pixel_size_m))

  ## --- read refine candidates ---
  ref <- sf::st_read(refine_path, quiet = TRUE)
  ref <- sf::st_make_valid(ref)
  ref <- ref[!sf::st_is_empty(ref), , drop = FALSE]
  if (nrow(ref) == 0L) stop("No candidates in refine_path.")
  if (is.na(sf::st_crs(ref))) stop("refine layer has no CRS.")

  if (is.null(year)) {
    m <- regmatches(basename(refine_path), regexpr("(?<=BA_)[0-9]{4}", basename(refine_path), perl = TRUE))
    year <- if (length(m)) as.integer(m) else NA_integer_
  }
  if (is.null(out_dir)) out_dir <- dirname(dirname(refine_path))  # the run dir (parent of 02_REFINE)
  chain_dir <- file.path(out_dir, "CHAIN", sprintf("dechain_%gm", dechain_distance_m))
  out_path  <- file.path(chain_dir, sprintf("BA_%s_REFINE_MERGED_dechain%gm.gpkg",
                                            ifelse(is.na(year), "NA", year), dechain_distance_m))
  if (file.exists(out_path) && !isTRUE(overwrite))
    stop("Output exists and overwrite = FALSE: ", out_path)
  dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)

  ## --- geometry + inflation (package convention: compact_pp = 4*pi*A/P^2) ---
  A  <- as.numeric(sf::st_area(ref))
  P  <- as.numeric(sf::st_length(sf::st_boundary(sf::st_geometry(ref))))
  compact_pp <- ifelse(P > 0, 4 * pi * A / (P^2), NA_real_)
  inflation  <- 1 / sqrt(compact_pp)
  area_ha    <- A / 1e4
  chained <- is.finite(inflation) & inflation > inflation_thr
  if (!is.null(min_area_ha)) chained <- chained & (area_ha >= min_area_ha)
  n_before <- nrow(ref); n_chained <- sum(chained)
  area_before_ha <- sum(area_ha, na.rm = TRUE)
  say("apply_chain_cleaning: %d candidates | %d flagged as chained (inflation > %.3f) | opening = %d px (%g m)",
      n_before, n_chained, inflation_thr, opening_px, dechain_distance_m)

  ## --- metadata template (carry only run metadata; no stale geometry stats) ---
  meta_cols <- intersect(c("YEAR","RUN","LB","DELTA","MINSEED","INDEX","WORKFLOW","ENGINE","SOURCE_AOI_FILE"),
                         names(ref))
  build <- function(geom, src_tag) {
    s <- sf::st_sf(geometry = sf::st_geometry(geom))
    for (cc in meta_cols) s[[cc]] <- if (cc == "SOURCE_AOI_FILE") src_tag else ref[[cc]][1]
    if ("RUN" %in% meta_cols) s[["RUN"]] <- paste0(ref[["RUN"]][1], sprintf("_dechain%gm", dechain_distance_m))
    s
  }

  ## --- seed mask (only marker = "seeds") ---
  seed_rast <- NULL
  if (identical(marker, "seeds")) {
    if (is.null(change_index) || is.null(vegetation_map) || is.null(seed_threshold_by_vegetation)) {
      warning("marker = 'seeds' needs change_index, vegetation_map and seed_threshold_by_vegetation; ",
              "falling back to marker = 'connected'.")
      marker <- "connected"
    }
  }

  n_necks <- 0L; area_removed_ha <- 0; comp_before <- n_chained; comp_after <- 0L
  pieces <- NULL

  if (n_chained > 0L) {
    ch <- ref[chained, ]
    ## template grid: chained extent + margin, snapped to pixel_size_m
    bb <- sf::st_bbox(ch)
    mg <- (opening_px + 2) * pixel_size_m
    snap <- function(v, down) { k <- if (down) floor(v / pixel_size_m) else ceiling(v / pixel_size_m); k * pixel_size_m }
    tmpl <- terra::rast(xmin = snap(bb["xmin"] - mg, TRUE), xmax = snap(bb["xmax"] + mg, FALSE),
                        ymin = snap(bb["ymin"] - mg, TRUE), ymax = snap(bb["ymax"] + mg, FALSE),
                        resolution = pixel_size_m, crs = sf::st_crs(ref)$wkt)
    mask <- terra::rasterize(terra::vect(ch), tmpl, background = 0)
    mask <- terra::ifel(mask > 0, 1, 0)
    n_burn_before <- as.numeric(terra::global(mask, "sum", na.rm = TRUE)[1, 1])
    ## morphological opening: erosion (focal min) then dilation (focal max) with a disk
    r <- opening_px
    dxy <- seq(-r, r); w <- outer(dxy, dxy, function(i, j) ifelse(i^2 + j^2 <= r^2, 1, NA))
    eroded  <- terra::focal(mask, w = w, fun = "min", na.policy = "all", fillvalue = 0)
    opened  <- terra::focal(eroded, w = w, fun = "max", na.policy = "all", fillvalue = 0)
    opened  <- terra::ifel(opened > 0, 1, NA)
    n_burn_after <- as.numeric(terra::global(terra::ifel(is.na(opened), 0, 1), "sum", na.rm = TRUE)[1, 1])
    area_removed_ha <- (n_burn_before - n_burn_after) * (pixel_size_m^2) / 1e4
    ## connected components -> polygons
    comp <- terra::patches(opened, directions = 8, zeroAsNA = TRUE)
    pv   <- terra::as.polygons(comp, dissolve = TRUE)
    pieces <- sf::st_make_valid(sf::st_as_sf(pv))
    pieces <- pieces[!sf::st_is_empty(pieces), , drop = FALSE]
    comp_after <- nrow(pieces)
    n_necks <- n_chained  # number of chained blobs that were opened
  }

  ## --- combine untouched (unchanged) + split pieces, renumber CLUMP_ID ---
  parts <- list()
  if (any(!chained)) parts[["untouched"]] <- build(ref[!chained, ], "untouched")
  if (!is.null(pieces) && nrow(pieces)) parts[["split"]] <- build(pieces, "dechain_split")
  out <- do.call(rbind, parts)
  out <- sf::st_make_valid(out)
  out <- out[!sf::st_is_empty(out), , drop = FALSE]
  out <- sf::st_cast(out, "MULTIPOLYGON", warn = FALSE)
  out[["CLUMP_ID"]] <- seq_len(nrow(out))
  out <- out[, c("CLUMP_ID", meta_cols)]
  sf::st_crs(out) <- sf::st_crs(ref)  # restore EPSG authority (terra polygonize keeps WKT only)
  n_after <- nrow(out)
  area_after_ha <- sum(as.numeric(sf::st_area(out)), na.rm = TRUE) / 1e4

  ## --- write outputs (never overwrite the original REFINE_MERGED) ---
  for (s in c("", "-wal", "-shm")) if (file.exists(paste0(out_path, s))) unlink(paste0(out_path, s), force = TRUE)
  sf::st_write(out, out_path, quiet = TRUE)

  audit_path  <- file.path(chain_dir, "chain_audit.csv")
  params_path <- file.path(chain_dir, "chain_params.yml")
  audit <- data.frame(
    year = year, n_candidates_before = n_before, n_candidates_after = n_after,
    n_chained_flagged = n_chained, bridges_detected = n_chained,
    components_before = comp_before, components_after = comp_after,
    patches_affected = n_chained, area_before_ha = round(area_before_ha, 2),
    area_after_ha = round(area_after_ha, 2), area_removed_ha = round(area_removed_ha, 2),
    dechain_distance_m = dechain_distance_m, pixel_size_m = pixel_size_m,
    opening_px = opening_px, inflation_thr = inflation_thr,
    min_area_ha = ifelse(is.null(min_area_ha), NA, min_area_ha),
    marker = marker, keep_noseed_fragments = keep_noseed_fragments)
  utils::write.csv(audit, audit_path, row.names = FALSE)
  writeLines(c(
    "# apply_chain_cleaning parameters",
    sprintf("year: %s", year),
    sprintf("refine_path: %s", refine_path),
    sprintf("output_path: %s", out_path),
    sprintf("dechain_distance_m: %g", dechain_distance_m),
    sprintf("pixel_size_m: %g", pixel_size_m),
    sprintf("opening_px: %d", opening_px),
    sprintf("inflation_thr: %g", inflation_thr),
    sprintf("min_area_ha: %s", ifelse(is.null(min_area_ha), "null", as.character(min_area_ha))),
    sprintf("marker: %s", marker),
    sprintf("keep_noseed_fragments: %s", tolower(as.character(keep_noseed_fragments))),
    sprintf("n_candidates_before: %d", n_before),
    sprintf("n_candidates_after: %d", n_after),
    sprintf("n_chained_flagged: %d", n_chained),
    sprintf("area_removed_ha: %.2f", area_removed_ha)
  ), params_path)
  say("apply_chain_cleaning DONE: %d -> %d candidates | removed %.0f ha of necks -> %s",
      n_before, n_after, area_removed_ha, out_path)

  list(output_path = out_path, audit_path = audit_path, params_path = params_path,
       n_candidates_before = n_before, n_candidates_after = n_after,
       n_chained_flagged = n_chained, area_before_ha = area_before_ha,
       area_after_ha = area_after_ha, area_removed_ha = area_removed_ha,
       opening_px = opening_px, dechain_distance_m = dechain_distance_m)
}

#' @rdname apply_chain_cleaning
#' @export
apply_chain <- apply_chain_cleaning
