#' Merge AOI vector outputs into a single vector file
#'
#' Recursively scans an AOI output root folder (e.g., `REFINE_AOI_5pix/`) and
#' merges all matching vector files into a single output file. Designed to be
#' run *after* AOI processing finishes (independent of the refinement function).
#'
#' Supported inputs: ESRI Shapefile (`.shp`), GeoPackage (`.gpkg`), GeoJSON
#' (`.geojson`/`.json`), and FlatGeobuf (`.fgb`). The output can be `.shp` or
#' `.gpkg` (chosen by `out_path` extension).
#'
#' @param base_output_dir Character. Existing directory to scan recursively.
#' @param pattern Character or NULL. Optional regex to filter filenames (matched
#'   against `basename(file)` with `ignore.case = TRUE`). If NULL, only the file
#'   extension filter is applied.
#' @param out_path Character. Output path ending in `.shp` or `.gpkg`.
#' @param dissolve Logical. If TRUE, dissolves all geometries into a single
#'   (multi)polygon feature (attributes are dropped).
#' @param allowed_ext Character vector. Allowed input extensions (lower/upper
#'   case accepted).
#' @param add_source_field Logical. If TRUE, adds a source field indicating the
#'   originating file.
#' @param source_field Character. Name of the source field to add.
#' @param keep_full_path Logical. If TRUE, stores full paths in `source_field`,
#'   otherwise stores `basename(file)`.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return The `out_path` (invisibly). Stops with an error if nothing is found
#'   or if all inputs are empty/invalid.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' merge_aoi_shapefiles(
#'   base_output_dir = "C:/ZENODO/exdata/OTSU_MIN/CORINE_ECOREG_seed_15pix/REFINE_AOI_5pix",
#'   pattern = "^BA_.*\\.(shp|gpkg)$",
#'   out_path = "C:/ZENODO/exdata/OTSU_MIN/CORINE_ECOREG_seed_15pix/BA_2012_REFINE_MERGED.gpkg",
#'   dissolve = FALSE,
#'   verbose = TRUE
#' )
#' }
merge_aoi_shapefiles <- function(base_output_dir,
                                 pattern = NULL,
                                 out_path,
                                 dissolve = FALSE,
                                 allowed_ext = c("shp", "gpkg", "geojson", "json", "fgb"),
                                 add_source_field = TRUE,
                                 source_field = "SOURCE_AOI_FILE",
                                 keep_full_path = FALSE,
                                 verbose = TRUE,

                                 # NEW: robust column dropping (default ON)
                                 strip_ecoregion = TRUE,
                                 ecoregion_field = NULL,
                                 drop_pattern = "(^ens|^enz|ecoreg|ecoregion|ecorreg|(^|_)eco(_|$)|cor_eco|coreco)",
                                 drop_cols = NULL) {

  if (!dir.exists(base_output_dir)) stop("base_output_dir does not exist: ", base_output_dir)
  if (!is.character(out_path) || length(out_path) != 1) stop("out_path must be a single character path.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (!requireNamespace("tools", quietly = TRUE)) stop("Package 'tools' is required.")

  # ----------------------------
  # Helpers
  # ----------------------------
  msg <- function(...) if (isTRUE(verbose)) message(...)

  normalize_path_safe <- function(p) {
    tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE), error = function(e) p)
  }

  delete_shapefile_bundle <- function(path_shp) {
    base <- tools::file_path_sans_ext(path_shp)
    exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx", ".xml")
    for (e in exts) {
      f <- paste0(base, e)
      if (file.exists(f)) try(file.remove(f), silent = TRUE)
    }
    invisible(TRUE)
  }

  choose_driver <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "gpkg") return("GPKG")
    if (ext == "shp")  return("ESRI Shapefile")
    stop("Unsupported out_path extension: '", ext, "'. Use .gpkg or .shp")
  }

  safe_make_valid_polys_local <- function(x) {
    if (is.null(x)) return(NULL)
    x <- sf::st_as_sf(x)
    if (nrow(x) == 0) return(NULL)

    emp <- sf::st_is_empty(sf::st_geometry(x))
    if (any(emp)) x <- x[!emp, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)

    inv <- !sf::st_is_valid(x)
    if (any(inv)) x <- sf::st_make_valid(x)

    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
    if (is.null(x) || nrow(x) == 0) return(NULL)

    sf::st_geometry(x) <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)

    emp2 <- sf::st_is_empty(sf::st_geometry(x))
    if (any(emp2)) x <- x[!emp2, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)

    x
  }

  ensure_geom_col <- function(x, geom_name = "geometry") {
    gcol <- attr(x, "sf_column")
    if (is.null(gcol) || !nzchar(gcol)) gcol <- "geometry"
    if (!identical(gcol, geom_name)) {
      names(x)[names(x) == gcol] <- geom_name
      attr(x, "sf_column") <- geom_name
      sf::st_geometry(x) <- geom_name
    }
    x
  }

  # Robust dropper: pattern + optional explicit ecoregion_field variants + explicit drop_cols
  strip_cols_local <- function(x) {

    if (!isTRUE(strip_ecoregion)) return(x)

    nms <- names(x)
    geom_col <- attr(x, "sf_column")
    if (is.null(geom_col) || !nzchar(geom_col)) geom_col <- "geometry"

    drop <- grep(
      "(^Area_km2|^ens|^enz|ecoreg|ecoregion|ecorreg|eco_|_eco|cor_eco|coreco)",
      nms,
      ignore.case = TRUE,
      value = TRUE
    )

    drop <- setdiff(unique(drop), geom_col)

    if (length(drop) == 0) return(x)

    x[, setdiff(nms, drop), drop = FALSE]
  }

  make_same_cols <- function(x, all_names) {
    miss <- setdiff(all_names, names(x))
    if (length(miss) > 0) {
      for (m in miss) x[[m]] <- NA
    }
    x[, all_names, drop = FALSE]
  }

  # ----------------------------
  # List inputs
  # ----------------------------
  all_files <- list.files(base_output_dir, full.names = TRUE, recursive = TRUE)
  if (!length(all_files)) stop("No files found in base_output_dir.")

  allowed_ext <- tolower(allowed_ext)
  ext <- tolower(tools::file_ext(all_files))
  files <- all_files[ext %in% allowed_ext]

  # If pattern is NULL, and BA_AOI_###.shp exists, auto-prefer those
  if (is.null(pattern) || !nzchar(pattern)) {
    bn <- basename(files)
    ba_only <- files[tolower(tools::file_ext(files)) == "shp" & grepl("^BA_AOI_\\d+\\.shp$", bn)]
    if (length(ba_only) > 0) {
      files <- ba_only
      msg(sprintf("[MERGE] pattern=NULL -> auto-selected %d BA_AOI_*.shp files.", length(files)))
    } else {
      msg(sprintf("[MERGE] pattern=NULL -> using all vector files (%d).", length(files)))
    }
  } else {
    bn <- basename(files)
    files <- files[grepl(pattern, bn, ignore.case = TRUE)]
  }

  if (!length(files)) stop("No vector files found to merge (check base_output_dir/pattern/allowed_ext).")

  # Exclude output if it is inside base_output_dir
  out_norm <- normalize_path_safe(out_path)
  files_norm <- vapply(files, normalize_path_safe, character(1))
  files <- files[files_norm != out_norm]
  if (!length(files)) stop("Only the output file was found; nothing to merge.")

  msg(sprintf("[MERGE] Found %d input vector files.", length(files)))

  # ----------------------------
  # Read + validate + strip cols
  # ----------------------------
  lst <- lapply(files, function(f) {
    x <- tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) NULL)
    x <- safe_make_valid_polys_local(x)
    if (is.null(x) || nrow(x) == 0) return(NULL)

    x <- ensure_geom_col(x, "geometry")

    # Strip eco/intersection columns BEFORE harmonizing
    x <- strip_cols_local(x)

    if (isTRUE(add_source_field)) {
      x[[source_field]] <- if (isTRUE(keep_full_path)) f else basename(f)
    }
    x
  })

  lst <- Filter(function(x) !is.null(x) && inherits(x, "sf") && nrow(x) > 0, lst)
  if (!length(lst)) stop("All candidate vector files were empty/invalid/unreadable.")

  # ----------------------------
  # CRS unify
  # ----------------------------
  crs0 <- sf::st_crs(lst[[1]])
  if (is.na(crs0)) stop("First input has NA CRS; cannot merge safely.")
  lst <- lapply(lst, function(x) sf::st_transform(x, crs0))

  # ----------------------------
  # Harmonize columns
  # ----------------------------
  all_names <- unique(unlist(lapply(lst, names)))
  all_names <- setdiff(all_names, "geometry")
  all_names <- c(all_names, "geometry")

  lst <- lapply(lst, make_same_cols, all_names = all_names)
  merged <- do.call(rbind, lst)
  merged <- sf::st_as_sf(merged)
  sf::st_crs(merged) <- crs0

  # Safety: strip again on merged (in case any slipped through)
  merged <- strip_cols_local(merged)

  # ----------------------------
  # Optional dissolve (geometry-only)
  # ----------------------------
  if (isTRUE(dissolve)) {
    g <- sf::st_union(sf::st_make_valid(sf::st_geometry(merged)))
    g <- sf::st_collection_extract(sf::st_sfc(g, crs = crs0), "POLYGON", warn = FALSE)
    if (length(g) == 0) stop("Dissolve produced no polygon geometry.")
    merged <- sf::st_sf(geometry = sf::st_cast(sf::st_union(g), "MULTIPOLYGON", warn = FALSE))
    sf::st_crs(merged) <- crs0
  }

  # ----------------------------
  # Write output robustly
  # ----------------------------
  driver <- choose_driver(out_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  if (driver == "ESRI Shapefile") {
    delete_shapefile_bundle(out_path)
    sf::st_write(merged, out_path, driver = driver, quiet = !verbose)
  } else {
    if (file.exists(out_path)) try(unlink(out_path, force = TRUE), silent = TRUE)
    sf::st_write(merged, out_path, driver = driver, delete_dsn = TRUE, quiet = !verbose)
  }

  msg(sprintf("[MERGE] Wrote merged file: %s", out_path))
  invisible(out_path)
}
