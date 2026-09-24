#' Remove narrow connections from refined candidate burned patches
#'
#' @description
#' Optionally clean narrow connections between parts of refined candidate
#' burned patches. These connections can arise during region growing and join
#' otherwise separate patches.
#'
#' The function identifies candidates with unusually complex shapes, applies
#' morphological opening to those candidates, and writes a separate layer for
#' subsequent scoring. Candidates not selected for cleaning are retained
#' without rasterisation or morphological processing.
#'
#' Run this function after [detect_burned_patches()] and pass the cleaned
#' output to [score_burned_patches()]. The original refined candidate layer
#' is preserved.
#'
#' `apply_chain()` is an alias for `apply_chain_cleaning()`.
#'
#' @param refine_path Character scalar. Path to the refined candidate
#'   GeoPackage, typically `02_REFINE/BA_<year>_REFINE_MERGED_*.gpkg`. Use the
#'   `refined_patches_path` returned by [detect_burned_patches()].
#' @param out_dir Character scalar or `NULL`. Run directory under which
#'   `CHAIN/dechain_<N>m/` is created. When `NULL`, uses
#'   `dirname(dirname(refine_path))`, assuming the input is inside a
#'   `02_REFINE` folder.
#' @param year Integer scalar or `NULL`. Target year. When `NULL`, the year is
#'   extracted from the input filename.
#' @param dechain_distance_m Numeric scalar. Morphological opening radius in
#'   metres, converted to pixels using `pixel_size_m`. Larger values produce
#'   stronger cleaning. Default: `180`.
#' @param pixel_size_m Numeric scalar. Grid-cell size in metres used to
#'   rasterise selected candidates. Default: `90`.
#' @param inflation_thr Numeric scalar. Shape-complexity threshold used to
#'   select candidates for cleaning. Candidates are selected when
#'   `1 / sqrt(compact_pp)` exceeds this value. Default: `22.135`.
#' @param min_area_ha Numeric scalar or `NULL`. Optional minimum candidate
#'   area, in hectares, for selection. When supplied, candidates must also
#'   satisfy `area_ha >= min_area_ha`. Default: `NULL`, meaning no additional
#'   area requirement.
#' @param marker Character scalar. Component-retention mode: `"connected"` or
#'   `"seeds"`. Default: `"connected"`. See \strong{Component retention}
#'   below.
#' @param keep_noseed_fragments Logical scalar. With `marker = "seeds"`,
#'   whether to retain fragments containing no seed pixels. Default: `TRUE`.
#' @param change_index Character scalar or `NULL`. Path to the RBR
#'   change-index raster used to identify seeds when `marker = "seeds"`.
#' @param vegetation_map Character scalar or `NULL`. Path to the
#'   CORINE-compatible vegetation-class raster used when `marker = "seeds"`.
#' @param seed_threshold_by_vegetation Named numeric vector or `NULL`.
#'   Vegetation-specific seed thresholds matching those used during candidate
#'   generation. Required when `marker = "seeds"`.
#' @param overwrite Logical scalar. Whether existing chain-cleaning outputs may
#'   be replaced. The original refined layer is preserved. Default: `FALSE`.
#' @param verbose Logical scalar. Whether to display progress messages.
#'   Default: `TRUE`.
#'
#' @section Position in the workflow:
#' Chain cleaning is an optional step between refinement and scoring:
#' 1. Generate refined candidates with [detect_burned_patches()].
#' 2. Clean selected candidates with `apply_chain_cleaning()`.
#' 3. Score the cleaned layer with [score_burned_patches()].
#'
#' Cleaning is applied after refinement so that the refinement merging step
#' does not reconnect separated patches. Existing detection outputs can be
#' reused without repeating region growing or refinement.
#'
#' The standard deterministic pipeline does not automatically run this
#' optional step.
#'
#' @section Candidate selection:
#' Candidates are selected using a shape-complexity measure derived from
#' Polsby-Popper compactness:
#' \preformatted{
#' compact_pp = 4 * pi * A / P^2
#' inflation  = 1 / sqrt(compact_pp)
#' }
#' Here, `A` is polygon area and `P` is perimeter, expressed in consistent
#' units.
#'
#' A candidate is selected when `inflation > inflation_thr`. If `min_area_ha`
#' is supplied, the candidate must also satisfy `area_ha >= min_area_ha`.
#'
#' Shape complexity identifies candidates for inspection and cleaning; it
#' does not establish that a narrow connection is an artefact.
#'
#' `min_area_ha` controls which input candidates are processed. It is not a
#' minimum-area filter for the resulting fragments.
#'
#' @section Morphological cleaning:
#' Selected candidates are:
#' 1. rasterised at `pixel_size_m`;
#' 2. processed by morphological opening, consisting of erosion followed by
#'    dilation with a disk-shaped structuring element;
#' 3. separated into connected components;
#' 4. converted back to polygons.
#'
#' Opening can remove narrow connections and small protrusions. It can also
#' remove small components or alter boundaries within selected candidates.
#'
#' Candidates not selected for cleaning are retained without rasterisation or
#' morphological processing.
#'
#' @section Distance and pixel size:
#' The opening radius in pixels is calculated as
#' `opening_px = dechain_distance_m / pixel_size_m`.
#'
#' For a grid-cell size of 90 m:
#'
#' | `dechain_distance_m` | Opening radius in pixels | Nominal opening diameter |
#' |---|---|---|
#' | 90 | 1 | 180 m |
#' | 180 | 2 | 360 m |
#' | 270 | 3 | 540 m |
#'
#' The opening diameter is not an exact cutoff for connection width. Results
#' depend on connection shape, orientation, and rasterisation.
#'
#' If the distance-to-pixel ratio is not an integer, the function issues a
#' warning and rounds it using `round()`. An opening radius below one pixel
#' is rejected.
#'
#' @section Component retention:
#' | Mode | Behaviour |
#' |---|---|
#' | `"connected"` | Retain every connected component remaining after opening. |
#' | `"seeds"` | Use seed support to assess the resulting fragments. Requires `change_index`, `vegetation_map`, and `seed_threshold_by_vegetation`. |
#'
#' With `marker = "seeds"`:
#' * `keep_noseed_fragments = TRUE` retains fragments without seed support;
#' * `keep_noseed_fragments = FALSE` removes fragments without seed support.
#'
#' If seed mode is requested without the required inputs, the function issues
#' a warning and falls back to `"connected"` mode.
#'
#' @section Output location:
#' Cleaned outputs are written under `<out_dir>/CHAIN/dechain_<N>m/`.
#'
#' The original `REFINE_MERGED` input is preserved. To control the output
#' location explicitly, supply `out_dir`, particularly when the input is
#' outside the standard run-folder structure.
#'
#' @section Choosing cleaning settings:
#' Larger opening radii can separate genuine parts of a fire or remove small
#' burned patches. Lowering `inflation_thr` selects more candidates for
#' processing.
#'
#' Inspect the cleaned geometries and compare candidate counts and areas
#' before using the output for scoring. Area removed by cleaning should not
#' automatically be interpreted as corrected false-positive area.
#'
#' @return A named list containing output locations, processing settings, and
#'   before-and-after summaries.
#'
#' | Field | Description |
#' |---|---|
#' | `output_path` | Path to the cleaned candidate layer. |
#' | `audit_path` | Path to the cleaning audit output. |
#' | `params_path` | Path to the saved cleaning parameters. |
#' | `n_candidates_before` | Number of candidates before cleaning. |
#' | `n_candidates_after` | Number of candidates after cleaning. |
#' | `n_chained_flagged` | Number of candidates flagged for cleaning. |
#' | `area_before_ha` | Candidate area before cleaning, in hectares. |
#' | `area_after_ha` | Candidate area after cleaning, in hectares. |
#' | `area_removed_ha` | Reported area removed by cleaning, in hectares. |
#' | `opening_px` | Opening radius used, in pixels. |
#' | `dechain_distance_m` | Cleaning-distance setting, in metres. |
#'
#' Candidate counts may increase when a polygon is split into several
#' components.
#'
#' @seealso [build_burned_mapping_config()], [detect_burned_patches()],
#'   [score_burned_patches()].
#'
#' @examples
#' \dontrun{
#' # Configure an annual RBR workflow
#' config <- build_burned_mapping_config(
#'   change_index = "data/RBR_2022.tif",
#'   vegetation_map = "data/vegetation_classes.tif",
#'   burnable_mask = "data/burnable_mask.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022"
#' )
#'
#' # Generate refined candidates
#' detection <- detect_burned_patches(
#'   config = config,
#'   write_outputs = TRUE
#' )
#'
#' # Clean selected candidates using a two-pixel opening radius
#' cleaned <- apply_chain_cleaning(
#'   refine_path = detection$refined_patches_path,
#'   year = 2022L,
#'   dechain_distance_m = 180,
#'   pixel_size_m = 90
#' )
#'
#' # Inspect changes in candidate counts and area
#' cleaned$n_chained_flagged
#' cleaned$n_candidates_before
#' cleaned$n_candidates_after
#' cleaned$area_removed_ha
#'
#' # Load the cleaned layer for visual inspection
#' cleaned_patches <- sf::st_read(
#'   cleaned$output_path,
#'   quiet = TRUE
#' )
#' plot(sf::st_geometry(cleaned_patches))
#'
#' # Compare a smaller opening radius using the same original input
#' cleaned_90 <- apply_chain_cleaning(
#'   refine_path = detection$refined_patches_path,
#'   year = 2022L,
#'   dechain_distance_m = 90,
#'   pixel_size_m = 90
#' )
#'
#' # After inspection, score the selected cleaned output
#' scored <- score_burned_patches(
#'   burned_candidates = cleaned$output_path,
#'   config = config,
#'   write_outputs = TRUE
#' )
#'
#' table(
#'   scored$internal_decisions$class_final,
#'   useNA = "ifany"
#' )
#' }
#'
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
