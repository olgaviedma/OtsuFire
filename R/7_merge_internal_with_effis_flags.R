#' Merge internal and EFFIS scored polygons and propagate class flags by intersection
#'
#' Combines an internally scored polygon layer (e.g., Stage-2 final scoring)
#' with a scored EFFIS/reference polygon layer by spatial overlap.
#'
#' The function:
#' \itemize{
#'   \item For each internal polygon, finds intersecting EFFIS polygons and assigns
#'         the EFFIS flag from the polygon with the maximum intersection area.
#'   \item Creates an \code{internal_enriched} layer with both internal and matched EFFIS flags,
#'         plus a \code{merged_flag_value} based on \code{assign_rule}.
#'   \item Extracts \code{effis_only} polygons that do not intersect any internal polygon.
#'   \item Writes a GeoPackage with three layers: \code{internal_enriched}, \code{effis_only}, \code{combined}.
#' }
#'
#' @param internal_scored An \code{sf} object with internal scored polygons, or a list with
#'   element \code{scored_sf} containing the \code{sf}.
#' @param effis_scored An \code{sf} object with EFFIS/reference scored polygons, or a list with
#'   element \code{scored_sf} containing the \code{sf}.
#' @param out_dir Character. Output directory where the GeoPackage will be written.
#' @param gpkg_name Character. GeoPackage filename to write inside \code{out_dir}.
#'
#' @param internal_flag_col Character. Column name in \code{internal_scored} containing the internal class flag.
#' @param effis_flag_col Character. Column name in \code{effis_scored} containing the EFFIS/reference class flag.
#'
#' @param assign_rule Character. How to compute \code{merged_flag_value}:
#'   \itemize{
#'     \item \code{"priority"}: choose the best label according to \code{priority_levels}.
#'     \item \code{"prefer_internal"}: internal label preferred (fallback to EFFIS).
#'     \item \code{"prefer_effis"}: EFFIS label preferred (fallback to internal).
#'   }
#' @param priority_levels Character vector defining best-to-worst ordering when \code{assign_rule="priority"}.
#'
#' @param buffer_m Numeric. Optional buffer (meters) applied to EFFIS polygons for intersection tests
#'   (useful to reduce sliver misses). Default 0 (no buffer).
#' @param chunk_size Integer. Chunk size for iterating internal polygons (memory/performance control).
#'
#' @param arcgis_fix Logical. If TRUE, applies GIS safety steps (valid geometries, MULTIPOLYGON casting,
#'   drop Z/M, enforce EPSG) before writing.
#' @param output_epsg Integer or NULL. EPSG code to enforce (default 3035). If NULL, CRS is kept as-is.
#'
#' @param overwrite Logical. If TRUE, overwrites an existing GeoPackage at the same path.
#' @param quiet Logical. Passed to \code{sf::st_write()}.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{gpkg_path}: output GeoPackage path
#'   \item \code{internal_enriched}: internal polygons with matched EFFIS info and merged flag
#'   \item \code{effis_only}: EFFIS polygons not intersecting internal polygons
#'   \item \code{combined}: rbind of \code{internal_enriched} + \code{effis_only} with harmonized schema
#' }
#'
#' @examples
#' \dontrun{
#' # internal_scored <- sf::st_read("phase2_internal_scoring.gpkg", "scored")
#' # effis_scored    <- sf::st_read("effis_scored.gpkg", "scored")
#' res <- merge_internal_with_effis_flags(
#'   internal_scored = internal_scored,
#'   effis_scored    = effis_scored,
#'   out_dir         = "C:/tmp/phase4_merge",
#'   assign_rule     = "priority",
#'   priority_levels = c("good_identified", "ambiguous", "bad_identified"),
#'   buffer_m        = 0
#' )
#' sf::st_layers(res$gpkg_path)
#' }
#'
#' @export
merge_internal_with_effis_flags <- function(
    internal_scored,
    effis_scored,
    out_dir,
    gpkg_name = "internal_plus_effis_scored.gpkg",

    internal_flag_col = "class3_id",
    effis_flag_col    = "class3_id",

    assign_rule = c("priority", "prefer_internal", "prefer_effis"),
    priority_levels = c("good_identified", "ambiguous", "bad_identified"),

    buffer_m = 0,
    chunk_size = 200,

    # ArcGIS / GIS safety
    arcgis_fix = TRUE,
    output_epsg = 3035,

    overwrite = TRUE,
    quiet = TRUE,
    verbose = TRUE
) {
  assign_rule <- match.arg(assign_rule)

  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Could not create out_dir: ", out_dir)

  # allow passing lists (from previous functions)
  if (is.list(internal_scored) && "scored_sf" %in% names(internal_scored)) internal_scored <- internal_scored$scored_sf
  if (is.list(effis_scored)    && "scored_sf" %in% names(effis_scored))    effis_scored    <- effis_scored$scored_sf

  stopifnot(inherits(internal_scored, "sf"), inherits(effis_scored, "sf"))

  if (!internal_flag_col %in% names(internal_scored)) stop("internal_flag_col not found: ", internal_flag_col)
  if (!effis_flag_col %in% names(effis_scored))       stop("effis_flag_col not found: ", effis_flag_col)

  if (nrow(internal_scored) == 0) stop("internal_scored has 0 rows.")
  if (nrow(effis_scored) == 0)    stop("effis_scored has 0 rows.")

  # -------------------------
  # Helpers
  # -------------------------
  drop_empty <- function(x) {
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    x
  }

  polys_only <- function(x) {
    x <- sf::st_make_valid(x)
    x <- drop_empty(x)
    gt <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
    if (any(gt %in% "GEOMETRYCOLLECTION")) {
      x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
      x <- drop_empty(x)
    }
    x <- suppressWarnings(sf::st_cast(x, "MULTIPOLYGON", warn = FALSE))
    x <- drop_empty(x)
    x
  }

  fix_for_gis <- function(x) {
    x <- polys_only(x)
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")

    if (!is.null(output_epsg)) {
      if (is.na(sf::st_crs(x))) {
        sf::st_crs(x) <- output_epsg
      } else {
        x <- sf::st_transform(x, output_epsg)
      }
      sf::st_crs(x) <- sf::st_crs(output_epsg)
    }
    x
  }

  make_gpkg_safe <- function(x, sep = "|") {
    stopifnot(inherits(x, "sf"))
    df <- sf::st_drop_geometry(x)
    is_list <- vapply(df, is.list, logical(1))
    if (any(is_list)) {
      for (col in names(is_list)[is_list]) {
        x[[col]] <- vapply(x[[col]], function(v) {
          if (is.null(v) || length(v) == 0) return(NA_character_)
          if (is.atomic(v)) return(paste(as.character(v), collapse = sep))
          paste(utils::capture.output(str(v)), collapse = " ")
        }, character(1))
      }
    }

    geom_col <- attr(x, "sf_column")
    nm <- names(x)
    idx_attr <- which(nm != geom_col)

    nm2 <- nm
    nm2[idx_attr] <- tolower(nm2[idx_attr])
    nm2[idx_attr] <- gsub("[^a-z0-9_]+", "_", nm2[idx_attr])
    nm2[idx_attr] <- gsub("^_+|_+$", "", nm2[idx_attr])
    nm2[idx_attr] <- ifelse(nm2[idx_attr] == "", paste0("v", idx_attr), nm2[idx_attr])
    nm2[idx_attr] <- make.unique(nm2[idx_attr], sep = "_")

    names(x) <- nm2
    x
  }

  pick_priority <- function(a, b) {
    ra <- match(a, priority_levels); rb <- match(b, priority_levels)
    if (is.na(a) || !nzchar(a)) return(b)
    if (is.na(b) || !nzchar(b)) return(a)
    if (is.na(ra) && is.na(rb)) return(a)
    if (is.na(ra)) return(b)
    if (is.na(rb)) return(a)
    if (ra <= rb) a else b
  }

  harmonize_for_rbind <- function(a, b) {
    geom_a <- attr(a, "sf_column")
    geom_b <- attr(b, "sf_column")
    if (geom_a != geom_b) {
      names(a)[names(a) == geom_a] <- "geom"
      names(b)[names(b) == geom_b] <- "geom"
      attr(a, "sf_column") <- "geom"
      attr(b, "sf_column") <- "geom"
    }

    all_cols <- union(names(a), names(b))
    for (nm in setdiff(all_cols, names(a))) a[[nm]] <- NA
    for (nm in setdiff(all_cols, names(b))) b[[nm]] <- NA

    a <- a[, all_cols, drop = FALSE]
    b <- b[, all_cols, drop = FALSE]

    geom_col <- attr(a, "sf_column")
    for (nm in setdiff(all_cols, geom_col)) {
      ca <- class(a[[nm]])[1]
      cb <- class(b[[nm]])[1]
      if (!identical(ca, cb)) {
        a[[nm]] <- as.character(a[[nm]])
        b[[nm]] <- as.character(b[[nm]])
      }
    }
    list(a = a, b = b)
  }

  # -------------------------
  # Prepare inputs (GIS-safe + CRS)
  # -------------------------
  internal <- fix_for_gis(internal_scored)
  effis    <- fix_for_gis(effis_scored)

  if (sf::st_crs(internal) != sf::st_crs(effis)) {
    effis <- sf::st_transform(effis, sf::st_crs(internal))
  }

  internal$._int_id <- seq_len(nrow(internal))
  effis$._ref_id    <- seq_len(nrow(effis))

  # -------------------------
  # Match EFFIS to internal by max intersection area
  # -------------------------
  ref_geom <- effis
  if (is.numeric(buffer_m) && is.finite(buffer_m) && buffer_m > 0) {
    ref_geom <- sf::st_buffer(ref_geom, dist = buffer_m)
  }

  hits <- sf::st_intersects(internal, ref_geom)

  effis_match_id   <- rep(NA_integer_, nrow(internal))
  effis_match_flag <- rep(NA_character_, nrow(internal))

  for (i0 in seq(1, nrow(internal), by = chunk_size)) {
    i1 <- min(nrow(internal), i0 + chunk_size - 1)

    for (i in i0:i1) {
      cand <- hits[[i]]
      if (length(cand) == 0) next

      best_id <- NA_integer_
      best_area <- -Inf

      for (j in cand) {
        inter <- suppressWarnings(tryCatch(
          sf::st_intersection(internal[i, ], effis[j, ]),
          error = function(e) NULL
        ))
        if (is.null(inter) || nrow(inter) == 0) next

        a <- as.numeric(sf::st_area(inter))
        a <- sum(a[is.finite(a)])
        if (is.finite(a) && a > best_area) {
          best_area <- a
          best_id <- j
        }
      }

      if (is.finite(best_area) && best_area > 0 && is.finite(best_id)) {
        effis_match_id[i]   <- effis$._ref_id[best_id]
        effis_match_flag[i] <- as.character(effis[[effis_flag_col]][best_id])
      }
    }
  }

  internal_enriched <- internal
  internal_enriched$effis_match_id        <- effis_match_id
  internal_enriched$effis_flag_value      <- effis_match_flag
  internal_enriched$internal_flag_value   <- as.character(internal_enriched[[internal_flag_col]])

  merged <- character(nrow(internal_enriched))
  for (i in seq_len(nrow(internal_enriched))) {
    a <- internal_enriched$internal_flag_value[i]
    b <- internal_enriched$effis_flag_value[i]
    if (assign_rule == "prefer_internal") {
      merged[i] <- if (!is.na(a) && nzchar(a)) a else b
    } else if (assign_rule == "prefer_effis") {
      merged[i] <- if (!is.na(b) && nzchar(b)) b else a
    } else {
      merged[i] <- pick_priority(a, b)
    }
  }
  internal_enriched$merged_flag_value <- merged
  internal_enriched$source_layer <- "internal_enriched"

  # -------------------------
  # NEW: copy EFFIS attributes into matched internal polygons
  # (so combined has EFFIS fields filled for internal_enriched rows)
  # -------------------------
  effis_df <- sf::st_drop_geometry(effis)

  # Map internal -> effis row index using the stored match id
  m <- match(internal_enriched$effis_match_id, effis_df$._ref_id)
  ii <- which(!is.na(m))

  make_na_vec <- function(proto, n) {
    if (inherits(proto, "Date")) return(as.Date(rep(NA_character_, n)))
    if (is.integer(proto))       return(rep(NA_integer_, n))
    if (is.numeric(proto))       return(rep(NA_real_, n))
    if (is.logical(proto))       return(rep(NA, n))
    if (is.factor(proto))        return(factor(rep(NA_character_, n), levels = levels(proto)))
    rep(NA_character_, n)
  }

  # Copy all EFFIS columns that are NOT already in internal_enriched
  # (this avoids overwriting internal fields like class3_id)
  cols_to_copy <- setdiff(names(effis_df), c(attr(effis, "sf_column")))  # geometry excluded anyway
  cols_to_copy <- setdiff(cols_to_copy, names(internal_enriched))        # do not overwrite internal columns
  cols_to_copy <- setdiff(cols_to_copy, c("._ref_id"))                   # keep match id separate

  if (length(ii) > 0 && length(cols_to_copy) > 0) {
    for (nm in cols_to_copy) {
      internal_enriched[[nm]] <- make_na_vec(effis_df[[nm]], nrow(internal_enriched))
      internal_enriched[[nm]][ii] <- effis_df[[nm]][m[ii]]
    }
  }

  # -------------------------
  # EFFIS only (no overlap with internal)
  # -------------------------
  hit_int <- lengths(sf::st_intersects(effis, internal)) > 0
  effis_only <- effis[!hit_int, , drop = FALSE]
  effis_only$internal_flag_value <- NA_character_
  effis_only$effis_flag_value    <- as.character(effis_only[[effis_flag_col]])
  effis_only$merged_flag_value   <- effis_only$effis_flag_value
  effis_only$source_layer        <- "effis_only"

  # -------------------------
  # Combined with safe schema
  # -------------------------
  tmp <- harmonize_for_rbind(internal_enriched, effis_only)
  combined <- rbind(tmp$a, tmp$b)

  if (isTRUE(arcgis_fix)) {
    internal_enriched <- fix_for_gis(internal_enriched)
    effis_only        <- fix_for_gis(effis_only)
    combined          <- fix_for_gis(combined)
  }

  internal_enriched_w <- make_gpkg_safe(internal_enriched)
  effis_only_w        <- make_gpkg_safe(effis_only)
  combined_w          <- make_gpkg_safe(combined)

  # -------------------------
  # Write GPKG
  # -------------------------
  gpkg_path <- file.path(out_dir, gpkg_name)
  if (isTRUE(overwrite) && file.exists(gpkg_path)) unlink(gpkg_path)

  sf::st_write(internal_enriched_w, gpkg_path, layer = "internal_enriched",
               driver = "GPKG", quiet = quiet, delete_dsn = TRUE)

  sf::st_write(effis_only_w, gpkg_path, layer = "effis_only",
               driver = "GPKG", quiet = quiet, append = TRUE)

  sf::st_write(combined_w, gpkg_path, layer = "combined",
               driver = "GPKG", quiet = quiet, append = TRUE)

  msg("Wrote GPKG: %s", gpkg_path)

  list(
    gpkg_path = gpkg_path,
    internal_enriched = internal_enriched,
    effis_only = effis_only,
    combined = combined
  )
}
