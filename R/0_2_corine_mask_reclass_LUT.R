#' Crop/mask and reclassify CORINE rasters, with optional burnable masks, LUTs (CSV/XLSX), and vector outputs
#'
#' @description
#' Process one or more CORINE rasters by:
#' \itemize{
#'   \item cropping and masking to an AOI polygon (\code{peninsula_shapefile}),
#'   \item aligning to a template raster grid (\code{template_raster_path}, recommended) or to a target CRS/resolution,
#'   \item optionally reclassifying codes into grouped classes (\code{group_rcl_matrix}),
#'   \item optionally extracting burnable classes and/or producing a binary burnable mask,
#'   \item optionally building a group-level LUT (from CORINE legend) and saving it to CSV/XLSX,
#'   \item optionally polygonizing selected outputs and (optionally) joining LUT fields into the vectors.
#' }
#'
#' @details
#' \strong{Input handling}
#' \itemize{
#'   \item \code{corine_rasters} can be a character vector of paths or a named list of paths.
#'   \item If unnamed, names are derived from filenames (without extension).
#' }
#'
#' \strong{Alignment}
#' \itemize{
#'   \item If \code{template_raster_path} is provided, rasters are projected to the template grid using
#'   nearest-neighbor resampling (\code{method="near"}).
#'   \item Otherwise, rasters are projected to \code{target_crs} at \code{target_res} (also \code{method="near"}).
#' }
#'
#' \strong{Which codes are stored in the CORINE raster?}
#' \itemize{
#'   \item \code{corine_value_col="GRID_CODE"} means the raster stores \code{GRID_CODE} values (commonly 1--44).
#'   \item \code{corine_value_col="CLC_CODE"} means the raster stores \code{CLC_CODE} values (commonly 111--523).
#'   \item \code{group_rcl_matrix[,1]} (old codes) must match the code system used in the raster (\code{corine_value_col}).
#' }
#'
#' \strong{Grouped reclassification and exclusions}
#' \itemize{
#'   \item \code{group_rcl_matrix} must be a 2-column matrix/data.frame: \code{(old_value, new_value)}.
#'   \item \code{exclude_values} can be used to remove any specific \code{old_value} codes before reclassification.
#'   \item If \code{exclude_values} is NULL and \code{exclude_burnt33=TRUE}, a default burnt code is excluded:
#'   \code{33} for \code{GRID_CODE} rasters or \code{334} for \code{CLC_CODE} rasters.
#' }
#'
#' \strong{Burnable outputs}
#' \itemize{
#'   \item To generate burnable outputs, provide \code{burnable_gridcodes} OR \code{burnable_matrix_01}.
#'   \item \code{burnable_matrix_01} must have 2 columns: \code{(class_code, 0/1)}. Codes with 1 are burnable.
#'   \item If both are provided, \code{burnable_matrix_01} takes precedence and is used to derive \code{burnable_gridcodes}.
#'   \item \code{burnable_source} controls whether burnable codes are applied to the \emph{group raster}
#'   (\code{"group"}: group IDs after reclassification) or to the \emph{original} CORINE raster
#'   (\code{"original"}: original \code{GRID_CODE} or \code{CLC_CODE} values).
#'   \item If \code{burnable_source="group"}, \code{group_rcl_matrix} must be provided.
#'   \item \code{write_burnable_classes} saves a raster containing only the selected burnable classes from the chosen source.
#'   \item \code{write_binary_mask} saves a binary raster (1 for burnable, NA otherwise).
#' }
#'
#' \strong{Group LUT (legend-based) and export (CSV/XLSX)}
#' \itemize{
#'   \item The group LUT summarizes which CORINE categories fall into each group ID (from \code{group_rcl_matrix}).
#'   \item \code{corine_grid_table} must be a data.frame or a path to a CSV with columns:
#'   \code{GRID_CODE, CLC_CODE, LABEL1, LABEL2, LABEL3}.
#'   \item CSV legend files using either \code{;} or \code{,} separators are supported (auto-detected).
#'   \item \code{write_group_lut_csv=TRUE} writes \code{lut_min} (and optionally \code{lut_full}) to CSV.
#'   \item \code{write_group_lut_xlsx=TRUE} writes an XLSX with sheets \code{lut_min} and \code{lut_full}
#'   (requires package \pkg{openxlsx}).
#'   \item \code{lut_csv2=TRUE} uses \code{write.csv2()} (semicolon-friendly).
#' }
#'
#' \strong{Vectorization and LUT joins}
#' \itemize{
#'   \item If \code{vectorize=TRUE}, polygons are written for selected rasters in \code{vectorize_which}.
#'   \item \code{vectorize_method="terra"} uses \code{terra::as.polygons()}.
#'   \item \code{vectorize_method="gdal_polygonize"} calls GDAL's \code{gdal_polygonize.py} via Python and
#'   requires \code{python_exe} and \code{gdal_polygonize_script}.
#'   \item If \code{touches=TRUE} and \code{vectorize_method="gdal_polygonize"}, 8-connectivity (\code{-8}) is used.
#'   \item If \code{apply_lut_to_vectors=TRUE}, LUT fields are joined into vectors where applicable
#'   (group vectors always; burnable vectors depend on \code{burnable_source}).
#' }
#'
#' @param corine_rasters Character vector of raster paths or a named list of raster paths.
#' @param peninsula_shapefile Character(1). AOI polygon file path used to crop/mask each CORINE raster.
#' @param template_raster_path Character(1) or NULL. Template raster defining the target grid (recommended).
#' @param target_crs Character(1). Target CRS used only when \code{template_raster_path} is NULL (e.g., "EPSG:3035").
#' @param target_res Numeric(1). Target resolution used only when \code{template_raster_path} is NULL.
#'
#' @param group_rcl_matrix Matrix/data.frame or NULL. Two-column mapping \code{(old_value, new_value)} for grouped classes.
#' @param exclude_burnt33 Logical. If TRUE and \code{exclude_values} is NULL, excludes a default burnt code
#'   (33 for \code{GRID_CODE}, 334 for \code{CLC_CODE}) from \code{group_rcl_matrix} before classifying.
#' @param exclude_values Integer vector or NULL. Optional explicit \code{old_value} codes to exclude from reclassification.
#'
#' @param burnable_gridcodes Integer vector or NULL. Class codes considered burnable (interpreted according to \code{burnable_source}).
#' @param burnable_matrix_01 Matrix/data.frame or NULL. Two columns \code{(class_code, 0/1)}. If provided, overrides
#'   \code{burnable_gridcodes} by selecting \code{class_code} where value==1.
#' @param burnable_source Character(1). Where \code{burnable_gridcodes} apply: \code{"group"} or \code{"original"}.
#' @param corine_value_col Character(1). What code system is stored in the CORINE raster (and used in \code{group_rcl_matrix[,1]}):
#'   \code{"GRID_CODE"} or \code{"CLC_CODE"}.
#'
#' @param output_raster_dir Character(1). Directory to write output rasters.
#' @param output_vector_dir Character(1) or NULL. Directory to write vector outputs (required if \code{vectorize=TRUE}).
#' @param name_prefix Character(1). Prefix for output filenames.
#' @param out_base_names NULL or character. If NULL, base names are derived from input names. If length 1, the same
#'   base name is used for all rasters. If named, names must match \code{names(corine_rasters)}.
#' @param dedup_prefix Logical. If TRUE, avoids repeating \code{name_prefix} when input names already start with it.
#' @param overwrite Logical. If TRUE, existing outputs are removed before writing.
#' @param touches Logical. If TRUE, uses 8-connectivity when vectorizing with GDAL polygonize (\code{-8}).
#' @param compress Character(1). GDAL compression for GeoTIFF (e.g., "LZW").
#'
#' @param write_group_raster Logical. If TRUE, writes grouped reclassification raster (only if \code{group_rcl_matrix} is provided).
#' @param write_burnable_classes Logical. If TRUE, writes burnable classes raster (selected classes only).
#' @param write_binary_mask Logical. If TRUE, writes binary burnable mask (1 for burnable, NA otherwise).
#'
#' @param vectorize Logical. If TRUE, polygonizes selected rasters and writes vectors.
#' @param vectorize_which Character vector. Which rasters to vectorize:
#'   \code{"burnable_classes"}, \code{"binary_mask"}, \code{"group_raster"}.
#' @param vectorize_method Character(1). Polygonization method: \code{"gdal_polygonize"} or \code{"terra"}.
#' @param vector_driver Character(1). OGR driver name (e.g., \code{"ESRI Shapefile"}, \code{"GPKG"}, \code{"GeoJSON"}).
#' @param python_exe Character(1) or NULL. Path to Python executable (required for \code{vectorize_method="gdal_polygonize"}).
#' @param gdal_polygonize_script Character(1) or NULL. Path to \code{gdal_polygonize.py}
#'   (required for \code{vectorize_method="gdal_polygonize"}).
#'
#' @param corine_grid_table NULL, data.frame, or character(1) path to a CSV legend file with columns
#'   \code{GRID_CODE, CLC_CODE, LABEL1, LABEL2, LABEL3}. Separator \code{;} or \code{,} is auto-detected.
#' @param write_group_lut_csv Logical. If TRUE, write group LUT CSV(s) to disk.
#' @param lut_out_dir Character(1) or NULL. Output directory for LUT files (default: \code{output_raster_dir}).
#' @param lut_file_tag Character(1). File stem/tag appended to LUT outputs.
#' @param lut_main_level Character(1). Main label level used for group naming: \code{"LABEL1"}, \code{"LABEL2"}, or \code{"LABEL3"}.
#' @param lut_detail_level Character(1). Detail label level shown in parentheses: \code{"LABEL2"}, \code{"LABEL3"}, or \code{"none"}.
#' @param lut_main_strategy Character(1). How to choose the main label per group: \code{"mode"} or \code{"unique_concat"}.
#' @param lut_main_tie Character(1). Tie-breaking for \code{lut_main_strategy="mode"}: \code{"concat"} or \code{"first"}.
#' @param lut_max_detail Integer(1). Maximum number of detail labels before truncation.
#' @param lut_sep_main Character(1). Separator for main label concatenation.
#' @param lut_sep_detail Character(1). Separator for detail label concatenation.
#' @param lut_prefix Character(1). Prefix used for group tags (e.g., "G" -> "G01", "G02"...).
#' @param lut_pad_id Integer(1). Zero-padding for group IDs in tags.
#' @param lut_sort_main Logical. If TRUE, sort main labels before concatenation.
#' @param lut_sort_detail Logical. If TRUE, sort detail labels before concatenation.
#' @param lut_detail_with_counts Logical. If TRUE, include counts in detail labels (e.g., "Coniferous (n=12)").
#' @param lut_keep_grid_codes Logical. If TRUE, keep \code{GRID_CODE} lists in \code{lut_full}.
#' @param lut_keep_clc_codes Logical. If TRUE, keep \code{CLC_CODE} lists in \code{lut_full}.
#' @param lut_write_full_csv Logical. If TRUE, also writes \code{lut_full} (otherwise only \code{lut_min}).
#' @param write_group_lut_xlsx Logical. If TRUE, write LUT to XLSX with sheets \code{lut_min} and \code{lut_full}
#'   (requires \pkg{openxlsx}).
#' @param lut_xlsx_name Character(1) or NULL. Optional custom XLSX filename; if NULL, uses \code{paste0(prefix,"_",lut_file_tag,".xlsx")}.
#' @param lut_csv2 Logical. If TRUE, write CSV using \code{write.csv2()} (semicolon-friendly). If FALSE, use \code{write.csv()}.
#' @param lut_csv_encoding Character(1). Encoding used for LUT CSV writing (e.g., "UTF-8").
#' @param apply_lut_to_vectors Logical. If TRUE and \code{vectorize=TRUE}, joins LUT fields into polygon outputs where applicable.
#'
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return
#' A named list (one element per input raster). Each element is a list with:
#' \itemize{
#'   \item \code{input}: input raster path,
#'   \item optional raster outputs: \code{group_raster}, \code{burnable_classes_raster}, \code{binary_mask_raster},
#'   \item optional LUT paths (global, same for all elements): \code{group_lut_min_csv}, \code{group_lut_full_csv}, \code{group_lut_xlsx},
#'   \item optional vector outputs: \code{group_vector}, \code{burnable_classes_vector}, \code{binary_mask_vector}.
#' }
#'
#' @examples
#' \dontrun{
#' # -------------------------
#' # Example 1: burnable codes applied to GROUP raster (group IDs)
#' # Also writes group LUT to CSV + XLSX (requires openxlsx for XLSX)
#' # -------------------------
#' my_reclass <- matrix(c(
#'   22, 1,  # agroforestry
#'   26, 2,  # grass
#'   32, 3,  # sparsely vegetated
#'   33, 3,  # burnt areas (often excluded from segmentation)
#'   27, 4,  # shrubs (moors)
#'   28, 4,  # sclerophyllous veg
#'   29, 4,  # transitional wood
#'   23, 5,  # broadleaves
#'   25, 6,  # mixed
#'   24, 7,  # conifers
#'   1,  8,  2,  8,  3,  8,  4,  8,  5,  8,  6,  8,  7,  8,  8,  8,  9,  8, 10, 8, 11, 8,  # 1-11 -> 8 (artificial)
#'   12, 9, 13, 9, 14, 9, 15, 9, 16, 9, 17, 9, 18, 9, 19, 9, 20, 9, 21, 9,  # 12-21 -> 9 (agriculture)
#'   35,10, 36,10, 37,10, 38,10, 39,10, 40,10, 41,10, 42,10, 43,10, 44,10,  # 35-44 -> 10 (water)
#'   30,11, 31,11, 34,11  # bare soils etc. (avoid duplicating 32 here)
#' ), ncol = 2, byrow = TRUE)
#'
#' burnable_01_group <- matrix(c(
#'   1,1, 2,1, 3,1, 4,1, 5,1, 6,1, 7,1, 8,0,
#'   9,0, 10,0, 11,0
#' ), ncol = 2, byrow = TRUE)
#'
#' out <- corine_mask_reclass(
#'   corine_rasters       = "D:/PROJECT/MASKS/CORINE_IBERIA_GRIDCODE.tif",
#'   peninsula_shapefile  = "D:/PROJECT/AOI/iberian_peninsula.shp",
#'   template_raster_path = "D:/PROJECT/TEMPLATE/template_res90m.tif",
#'
#'   corine_value_col     = "GRID_CODE",
#'   group_rcl_matrix     = my_reclass,
#'   exclude_burnt33      = TRUE,  # excludes 33 from the reclass mapping (GRID_CODE convention)
#'
#'   burnable_source      = "group",
#'   burnable_matrix_01   = burnable_01_group,
#'
#'   write_group_raster     = TRUE,
#'   write_burnable_classes = TRUE,
#'   write_binary_mask      = TRUE,
#'
#'   # LUT (legend CSV can be ';' or ',' separated; auto-detected)
#'   corine_grid_table      = "D:/PROJECT/MASKS/clc_legend_gridcode.csv",
#'   write_group_lut_csv    = TRUE,
#'   write_group_lut_xlsx   = TRUE,
#'   lut_out_dir            = "D:/PROJECT/MASKS/LUT",
#'   lut_file_tag           = "strata8_v1",
#'   lut_main_level         = "LABEL2",
#'   lut_detail_level       = "LABEL3",
#'
#'   output_raster_dir    = "D:/PROJECT/MASKS/OUT",
#'   out_base_names       = "corine18",
#'   overwrite            = TRUE,
#'   verbose              = TRUE
#' )
#'
#' # -------------------------
#' # Example 2: burnable codes applied to ORIGINAL raster (original codes)
#' # -------------------------
#' burnable_gridcodes_original <- c(22, 23, 24, 25, 26, 27, 28, 29, 32, 33)
#'
#' out2 <- corine_mask_reclass(
#'   corine_rasters       = "D:/PROJECT/MASKS/CORINE_IBERIA_GRIDCODE.tif",
#'   peninsula_shapefile  = "D:/PROJECT/AOI/iberian_peninsula.shp",
#'   template_raster_path = "D:/PROJECT/TEMPLATE/template_res90m.tif",
#'
#'   corine_value_col     = "GRID_CODE",
#'   burnable_source      = "original",
#'   burnable_gridcodes   = burnable_gridcodes_original,
#'
#'   write_burnable_classes = TRUE,
#'   write_binary_mask      = TRUE,
#'
#'   output_raster_dir    = "D:/PROJECT/MASKS/OUT",
#'   out_base_names       = "corine18_originalBurnable",
#'   overwrite            = TRUE,
#'   verbose              = TRUE
#' )
#' }
#'
#' @importFrom terra rast vect crs project crop mask classify ifel writeRaster as.polygons writeVector nlyr freq
#' @importFrom sf st_read st_crs st_transform st_make_valid st_is_valid st_collection_extract st_as_sf st_write
#' @importFrom tools file_path_sans_ext file_ext
#' @export
utils::globalVariables(c("CLC_CODE", "GRID_CODE"))

  corine_mask_reclass <- function(
    corine_rasters,
    peninsula_shapefile,

    template_raster_path = NULL,
    target_crs = "EPSG:3035",
    target_res = 90,

    group_rcl_matrix = NULL,
    exclude_burnt33 = TRUE,
    exclude_values = NULL,            # NEW: optional explicit values to exclude from reclass

    burnable_gridcodes = NULL,
    burnable_matrix_01 = NULL,

    # NEW: where burnable codes apply
    burnable_source = c("group", "original"),

    # NEW: what codes are stored in the CORINE raster (and in group_rcl_matrix old codes)
    corine_value_col = c("GRID_CODE", "CLC_CODE"),

    output_raster_dir,
    output_vector_dir = NULL,
    name_prefix = "corine",
    out_base_names = NULL,
    dedup_prefix = TRUE,
    overwrite = TRUE,
    touches = TRUE,
    compress = "LZW",

    write_group_raster = TRUE,
    write_burnable_classes = TRUE,
    write_binary_mask = TRUE,

    vectorize = FALSE,
    vectorize_which = c("burnable_classes", "binary_mask", "group_raster"),
    vectorize_method = c("gdal_polygonize", "terra"),
    vector_driver = "ESRI Shapefile",
    python_exe = NULL,
    gdal_polygonize_script = NULL,

    # =========================
    # NEW: LUT functionality (from make_lut_from_gridcode)
    # =========================
    corine_grid_table = NULL,                # data.frame OR path to csv
    write_group_lut_csv = FALSE,             # write LUT CSV(s)
    lut_out_dir = NULL,                      # default: output_raster_dir
    lut_file_tag = "group_lut",              # file stem

    lut_main_level   = c("LABEL1","LABEL2","LABEL3"),
    lut_detail_level = c("LABEL2","LABEL3","none"),
    lut_main_strategy = c("mode", "unique_concat"),
    lut_main_tie      = c("concat", "first"),
    lut_max_detail    = 6,
    lut_sep_main      = " + ",
    lut_sep_detail    = "; ",
    lut_prefix        = "G",
    lut_pad_id        = 2,
    lut_sort_main     = TRUE,
    lut_sort_detail   = TRUE,
    lut_detail_with_counts = FALSE,
    lut_keep_grid_codes = TRUE,
    lut_keep_clc_codes  = TRUE,
    lut_write_full_csv  = TRUE,              # if TRUE writes lut_full, else only lut_min
    write_group_lut_xlsx = FALSE,   # NEW: write LUT to .xlsx (two sheets)
    lut_xlsx_name = NULL,           # NEW: optional custom filename; default uses prefix + lut_file_tag
    lut_csv2 = TRUE,                # NEW: write CSV with ';' (Spain-friendly)
    lut_csv_encoding = "UTF-8" ,      # NEW: encoding for CSV

    apply_lut_to_vectors = TRUE,             # join LUT fields into vectors when vectorizing

    verbose = TRUE
  ) {

    msg <- function(...) if (isTRUE(verbose)) message(...)
    if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
    if (!requireNamespace("tools", quietly = TRUE)) stop("Package 'tools' is required.")
    if (!requireNamespace("utils", quietly = TRUE)) stop("Package 'utils' is required.")

    burnable_source  <- match.arg(burnable_source)
    corine_value_col <- match.arg(corine_value_col)

    stopifnot(file.exists(peninsula_shapefile))
    if (!is.null(template_raster_path)) stopifnot(file.exists(template_raster_path))

    dir.create(output_raster_dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(output_vector_dir)) dir.create(output_vector_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- normalize corine_rasters to named list
    if (is.character(corine_rasters) && length(corine_rasters) >= 1) {
      paths <- corine_rasters
      nm <- tools::file_path_sans_ext(basename(paths))
      corine_rasters <- as.list(paths)
      names(corine_rasters) <- nm
    } else if (is.list(corine_rasters)) {
      if (is.null(names(corine_rasters)) || any(!nzchar(names(corine_rasters)))) {
        nm <- tools::file_path_sans_ext(basename(unlist(corine_rasters)))
        names(corine_rasters) <- nm
      }
    } else {
      stop("corine_rasters must be a character vector of paths or a named list of paths.")
    }

    # ---- clean prefix
    prefix <- if (is.null(name_prefix)) "" else as.character(name_prefix)
    prefix <- gsub("[^A-Za-z0-9_\\-]+", "_", prefix)
    prefix <- gsub("_+", "_", prefix)
    prefix <- sub("^_+|_+$", "", prefix)

    # ---- helper: build output base name without repeated prefix
    make_out_base <- function(nm, prefix, dedup = TRUE) {
      x <- as.character(nm)
      x <- gsub("\\s+", "_", x)
      x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- sub("^_+|_+$", "", x)
      if (!nzchar(x)) x <- "layer"

      if (nzchar(prefix)) {
        if (isTRUE(dedup)) {
          x <- sub(paste0("^", prefix, "_+"), "", x, ignore.case = TRUE)
          x <- sub("^_+|_+$", "", x)
          if (!nzchar(x)) x <- "layer"
        }
        return(paste0(prefix, "_", x))
      } else {
        return(x)
      }
    }

    resolve_out_base <- function(nm, out_base_names) {
      if (is.null(out_base_names)) return(NULL)
      if (is.list(out_base_names)) out_base_names <- unlist(out_base_names)
      if (!is.character(out_base_names)) stop("out_base_names must be NULL or character.")

      if (length(out_base_names) == 1 && nzchar(out_base_names)) return(out_base_names)

      if (!is.null(names(out_base_names)) && nzchar(nm) && (nm %in% names(out_base_names))) {
        val <- out_base_names[[nm]]
        if (length(val) == 1 && nzchar(val)) return(val)
      }

      stop("out_base_names must be length 1 or a named vector/list with names matching corine_rasters.")
    }

    # ---- overwrite helper (handles SHP sidecars)
    safe_overwrite <- function(path) {
      if (!file.exists(path)) return(invisible(TRUE))
      if (!isTRUE(overwrite)) stop("File exists and overwrite=FALSE: ", path)

      ext <- tolower(tools::file_ext(path))
      if (ext == "shp") {
        stem <- sub("\\.shp$", "", path, ignore.case = TRUE)
        sidecars <- c("shp","shx","dbf","prj","cpg","qix","fix","sbn","sbx")
        files <- c(path, file.path(dirname(path), paste0(basename(stem), ".", sidecars)))
        files <- unique(files[file.exists(files)])
        suppressWarnings(file.remove(files))
      } else {
        suppressWarnings(file.remove(path))
        aux <- paste0(path, ".aux.xml")
        if (file.exists(aux)) suppressWarnings(file.remove(aux))
      }
      if (file.exists(path)) stop("Cannot overwrite file (locked?): ", path)
      invisible(TRUE)
    }

    # ---- read peninsula once
    pen_sf0 <- sf::st_read(peninsula_shapefile, quiet = TRUE)
    if (nrow(pen_sf0) == 0) stop("peninsula_shapefile has 0 features: ", peninsula_shapefile)
    if (is.na(sf::st_crs(pen_sf0))) stop("peninsula_shapefile CRS is NA: ", peninsula_shapefile)
    if (any(!sf::st_is_valid(pen_sf0))) pen_sf0 <- sf::st_make_valid(pen_sf0)
    pen_sf0 <- suppressWarnings(sf::st_collection_extract(pen_sf0, "POLYGON", warn = FALSE))
    if (nrow(pen_sf0) == 0) stop("peninsula_shapefile has no POLYGON geometries after extract.")

    # ---- template (optional)
    tmpl <- NULL
    if (!is.null(template_raster_path)) {
      tmpl <- terra::rast(template_raster_path)
      if (is.na(terra::crs(tmpl, proj = TRUE)) || !nzchar(terra::crs(tmpl, proj = TRUE))) {
        stop("Template CRS missing/empty: ", template_raster_path)
      }
    }

    # ---- prepare reclass matrix (2 col)
    prep_rcl_2col <- function(mat, exclude_vals = NULL) {
      if (is.null(mat)) return(NULL)
      mat <- as.matrix(mat)
      if (ncol(mat) != 2) stop("Reclass matrix must have exactly 2 columns (old, new).")
      if (nrow(mat) == 0) stop("Reclass matrix has 0 rows.")
      mat[, 1] <- as.integer(round(mat[, 1]))
      mat[, 2] <- as.integer(round(mat[, 2]))
      if (!is.null(exclude_vals) && length(exclude_vals) > 0) {
        exclude_vals <- as.integer(round(exclude_vals))
        mat <- mat[!(mat[, 1] %in% exclude_vals), , drop = FALSE]
      }
      if (nrow(mat) == 0) stop("Reclass matrix empty after exclusions.")
      if (any(duplicated(mat[, 1]))) mat <- mat[!duplicated(mat[, 1], fromLast = TRUE), , drop = FALSE]
      mat
    }

    # Decide default exclude value if exclude_values is not provided but exclude_burnt33 is TRUE
    if (is.null(exclude_values) && isTRUE(exclude_burnt33)) {
      # Common conventions:
      # - GRID_CODE: burnt areas often coded as 33
      # - CLC_CODE: burnt areas often coded as 334
      exclude_values <- if (corine_value_col == "GRID_CODE") 33L else 334L
    }

    group_rcl_use <- prep_rcl_2col(group_rcl_matrix, exclude_vals = exclude_values)

    # ---- burnable selection
    if (!is.null(burnable_matrix_01)) {
      bm <- as.matrix(burnable_matrix_01)
      if (ncol(bm) != 2) stop("burnable_matrix_01 must have 2 columns: (class, 0/1).")
      bm[, 1] <- as.integer(round(bm[, 1]))
      bm[, 2] <- as.integer(round(bm[, 2]))
      burnable_gridcodes <- bm[bm[, 2] == 1, 1]
    }
    if (!is.null(burnable_gridcodes)) {
      burnable_gridcodes <- sort(unique(as.integer(round(burnable_gridcodes))))
      if (length(burnable_gridcodes) == 0) burnable_gridcodes <- NULL
    }

    # ==========================================================
    # NEW: LUT helper (integrated from make_lut_from_gridcode)
    # ==========================================================
    read_corine_grid_table <- function(x) {
      clean_names <- function(nm) {
        nm <- sub("^\ufeff", "", nm)   # remove UTF-8 BOM if present
        nm <- trimws(nm)
        nm
      }

      required <- c("GRID_CODE","CLC_CODE","LABEL1","LABEL2","LABEL3")

      detect_sep <- function(path) {
        l1 <- readLines(path, n = 1, warn = FALSE)
        if (length(l1) == 0) return(";")
        n_sc <- lengths(regmatches(l1, gregexpr(";", l1, fixed = TRUE)))
        n_cm <- lengths(regmatches(l1, gregexpr(",", l1, fixed = TRUE)))
        n_tb <- lengths(regmatches(l1, gregexpr("\t", l1, fixed = TRUE)))
        # pick the most frequent delimiter
        if (n_tb > n_sc && n_tb > n_cm) return("\t")
        if (n_cm > n_sc) return(",")
        return(";")
      }

      # Attempt to repair a 1-column data.frame that actually contains delimited text
      repair_onecol <- function(df) {
        if (!is.data.frame(df) || ncol(df) != 1) return(df)

        v <- as.character(df[[1]])
        if (length(v) == 0) return(df)

        # detect delimiter from first non-empty row
        idx <- which(nzchar(trimws(v)))[1]
        if (is.na(idx)) return(df)
        row1 <- v[idx]

        delim <- if (grepl(";", row1, fixed = TRUE)) ";" else if (grepl(",", row1, fixed = TRUE)) "," else NULL
        if (is.null(delim)) return(df)

        # derive header fields from current colname (which may have dots from check.names)
        hdr_raw <- clean_names(names(df)[1])
        hdr_fields <- unlist(strsplit(hdr_raw, "[;,.\\t]+"))
        hdr_fields <- hdr_fields[nzchar(hdr_fields)]

        if (length(hdr_fields) < 5 || !all(required %in% hdr_fields)) {
          # fallback to required names (most common legend format)
          hdr_fields <- required
        } else {
          # keep only required order first, then any extras
          hdr_fields <- c(required, setdiff(hdr_fields, required))
        }

        parts <- strsplit(v, delim, fixed = TRUE)

        maxlen <- max(lengths(parts), 0)
        maxlen <- max(maxlen, length(hdr_fields))

        mat <- t(vapply(parts, function(p) {
          p <- trimws(p)
          if (length(p) < maxlen) p <- c(p, rep(NA_character_, maxlen - length(p)))
          p[seq_len(maxlen)]
        }, character(maxlen)))

        out <- as.data.frame(mat, stringsAsFactors = FALSE)
        names(out) <- c(hdr_fields, paste0("X", seq_len(ncol(out) - length(hdr_fields))))
        out
      }

      read_with_sep <- function(path, sep) {
        df <- utils::read.csv(
          path,
          sep = sep,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          fileEncoding = "UTF-8-BOM"
        )
        names(df) <- clean_names(names(df))
        df <- repair_onecol(df)
        names(df) <- clean_names(names(df))
        df
      }

      if (is.null(x)) return(NULL)

      if (is.character(x)) {
        if (!file.exists(x)) stop("corine_grid_table file not found: ", x)

        sep0 <- detect_sep(x)
        df <- read_with_sep(x, sep0)

        # If still not OK, try alternate separator
        if (!(all(required %in% names(df)))) {
          sep1 <- if (sep0 == ";") "," else ";"
          df2 <- read_with_sep(x, sep1)
          if (all(required %in% names(df2))) df <- df2
        }

        return(df)
      }

      if (is.data.frame(x)) {
        df <- x
        names(df) <- clean_names(names(df))
        df <- repair_onecol(df)
        names(df) <- clean_names(names(df))
        return(df)
      }

      stop("corine_grid_table must be NULL, a data.frame, or a path to a CSV.")
    }

    build_group_lut <- function(
    reclass_matrix,
    corine_grid_table,
    corine_value_col,
    main_level, detail_level, main_strategy, main_tie,
    max_detail, sep_main, sep_detail, prefix, pad_id,
    sort_main, sort_detail, detail_with_counts,
    keep_grid_codes, keep_clc_codes
    ) {
      main_level    <- match.arg(main_level,   c("LABEL1","LABEL2","LABEL3"))
      detail_level  <- match.arg(detail_level, c("LABEL2","LABEL3","none"))
      main_strategy <- match.arg(main_strategy, c("mode","unique_concat"))
      main_tie      <- match.arg(main_tie,      c("concat","first"))

      map_df <- as.data.frame(reclass_matrix, stringsAsFactors = FALSE)
      map_df <- map_df[, 1:2, drop = FALSE]
      names(map_df) <- c("KEY", "id")
      map_df$KEY <- suppressWarnings(as.integer(map_df$KEY))
      map_df$id  <- suppressWarnings(as.integer(map_df$id))
      if (anyNA(map_df$KEY) || anyNA(map_df$id)) {
        stop("group_rcl_matrix has non-integer codes after coercion.")
      }

      cor <- corine_grid_table
      needed <- c("GRID_CODE","CLC_CODE","LABEL1","LABEL2","LABEL3")
      miss_cols <- setdiff(needed, names(cor))
      if (length(miss_cols)) stop("corine_grid_table missing columns: ", paste(miss_cols, collapse=", "))

      cor$GRID_CODE <- suppressWarnings(as.integer(cor$GRID_CODE))
      cor$CLC_CODE  <- suppressWarnings(as.integer(cor$CLC_CODE))
      cor$LABEL1    <- as.character(cor$LABEL1)
      cor$LABEL2    <- as.character(cor$LABEL2)
      cor$LABEL3    <- as.character(cor$LABEL3)

      join_col <- corine_value_col
      cor$KEY <- if (join_col == "GRID_CODE") cor$GRID_CODE else cor$CLC_CODE

      x <- merge(map_df, cor, by = "KEY", all.x = TRUE, sort = FALSE)

      if (any(is.na(x$CLC_CODE))) {
        miss <- unique(x$KEY[is.na(x$CLC_CODE)])
        stop("Missing codes in corine_grid_table for: ", paste(miss, collapse = ", "))
      }

      clean_chr <- function(v) {
        v <- as.character(v)
        v <- v[!is.na(v)]
        v <- trimws(v)
        v <- v[nzchar(v)]
        v
      }

      pick_main <- function(v) {
        v <- clean_chr(v)
        if (length(v) == 0) return(NA_character_)

        if (main_strategy == "unique_concat") {
          u <- unique(v)
          if (sort_main) u <- sort(u)
          if (length(u) == 1) return(u)
          return(paste(u, collapse = sep_main))
        }

        tab <- sort(table(v), decreasing = TRUE)
        top_n <- as.integer(tab[1])
        top_vals <- names(tab)[tab == top_n]
        if (sort_main) top_vals <- sort(top_vals)

        if (length(top_vals) == 1) {
          return(top_vals)
        } else {
          if (main_tie == "first") return(top_vals[1])
          return(paste(top_vals, collapse = sep_main))
        }
      }

      pick_detail <- function(v) {
        v <- clean_chr(v)
        if (length(v) == 0) return(NA_character_)

        if (detail_with_counts) {
          tab <- sort(table(v), decreasing = TRUE)
          vals <- names(tab)
          if (sort_detail) vals <- vals[order(-as.integer(tab[vals]), vals)]
          lab <- paste0(vals, " (n=", as.integer(tab[vals]), ")")
        } else {
          vals <- unique(v)
          if (sort_detail) vals <- sort(vals)
          lab <- vals
        }

        if (length(lab) > max_detail) lab <- c(lab[1:max_detail], "...")
        paste(lab, collapse = sep_detail)
      }

      ids <- sort(unique(x$id))
      main_vec   <- character(length(ids))
      detail_vec <- character(length(ids))
      grid_vec   <- character(length(ids))
      clc_vec    <- character(length(ids))

      for (i in seq_along(ids)) {
        id0 <- ids[i]
        xi <- x[x$id == id0, , drop = FALSE]

        main_vec[i] <- pick_main(xi[[main_level]])

        if (detail_level == "none") {
          detail_vec[i] <- NA_character_
        } else {
          detail_vec[i] <- pick_detail(xi[[detail_level]])
        }

        if (keep_grid_codes) {
          grid_vec[i] <- paste(sort(unique(xi$GRID_CODE)), collapse = ",")
        } else {
          grid_vec[i] <- NA_character_
        }

        if (keep_clc_codes) {
          clc_vec[i] <- paste(sort(unique(xi$CLC_CODE)), collapse = ",")
        } else {
          clc_vec[i] <- NA_character_
        }
      }

      id_tag <- paste0(prefix, sprintf(paste0("%0", pad_id, "d"), ids))

      label <- ifelse(
        is.na(detail_vec) | detail_level == "none",
        paste0(id_tag, ": ", main_vec),
        paste0(id_tag, ": ", main_vec, " (", detail_vec, ")")
      )

      lut_full <- data.frame(
        id         = ids,
        id_tag     = id_tag,
        main       = main_vec,
        detail     = detail_vec,
        gridcodes  = if (keep_grid_codes) grid_vec else NA_character_,
        clccodes   = if (keep_clc_codes)  clc_vec  else NA_character_,
        label      = label,
        stringsAsFactors = FALSE
      )

      lut_min <- lut_full[, c("id","label"), drop = FALSE]
      list(lut_min = lut_min, lut_full = lut_full)
    }

    cor_tbl <- read_corine_grid_table(corine_grid_table)

    # Build group LUT once (global), if requested/needed
    need_group_lut <- isTRUE(write_group_lut_csv) || (isTRUE(vectorize) && isTRUE(apply_lut_to_vectors))
    group_lut <- NULL
    group_lut_min_csv <- NULL
    group_lut_full_csv <- NULL
    group_lut_xlsx <- NULL

    if (isTRUE(need_group_lut)) {
      if (is.null(group_rcl_use)) {
        msg("[LUT][INFO] group_rcl_matrix is NULL; group LUT will not be created.")
      } else {
        if (is.null(cor_tbl)) stop("LUT requested but corine_grid_table is NULL. Provide a data.frame or a CSV path.")
        group_lut <- build_group_lut(
          reclass_matrix = group_rcl_use,
          corine_grid_table = cor_tbl,
          corine_value_col = corine_value_col,
          main_level = lut_main_level,
          detail_level = lut_detail_level,
          main_strategy = lut_main_strategy,
          main_tie = lut_main_tie,
          max_detail = lut_max_detail,
          sep_main = lut_sep_main,
          sep_detail = lut_sep_detail,
          prefix = lut_prefix,
          pad_id = lut_pad_id,
          sort_main = lut_sort_main,
          sort_detail = lut_sort_detail,
          detail_with_counts = lut_detail_with_counts,
          keep_grid_codes = lut_keep_grid_codes,
          keep_clc_codes = lut_keep_clc_codes
        )

        if (isTRUE(write_group_lut_csv) || isTRUE(write_group_lut_xlsx)) {

          lut_dir <- if (is.null(lut_out_dir)) output_raster_dir else lut_out_dir
          dir.create(lut_dir, recursive = TRUE, showWarnings = FALSE)

          # Default file names
          group_lut_min_csv  <- file.path(lut_dir, paste0(prefix, "_", lut_file_tag, "_min.csv"))
          group_lut_full_csv <- file.path(lut_dir, paste0(prefix, "_", lut_file_tag, "_full.csv"))

          # CSV (Spain-friendly: ';') if requested
          if (isTRUE(write_group_lut_csv)) {
            safe_overwrite(group_lut_min_csv)
            safe_overwrite(group_lut_full_csv)

            if (isTRUE(lut_csv2)) {
              utils::write.csv2(group_lut$lut_min,  group_lut_min_csv,  row.names = FALSE, fileEncoding = lut_csv_encoding)
              if (isTRUE(lut_write_full_csv)) {
                utils::write.csv2(group_lut$lut_full, group_lut_full_csv, row.names = FALSE, fileEncoding = lut_csv_encoding)
              } else {
                group_lut_full_csv <- NULL
              }
            } else {
              utils::write.csv(group_lut$lut_min,  group_lut_min_csv,  row.names = FALSE, fileEncoding = lut_csv_encoding)
              if (isTRUE(lut_write_full_csv)) {
                utils::write.csv(group_lut$lut_full, group_lut_full_csv, row.names = FALSE, fileEncoding = lut_csv_encoding)
              } else {
                group_lut_full_csv <- NULL
              }
            }

            msg("[LUT] Written: ", group_lut_min_csv)
            if (!is.null(group_lut_full_csv)) msg("[LUT] Written: ", group_lut_full_csv)
          } else {
            # If not writing CSV, keep paths NULL
            group_lut_min_csv  <- NULL
            group_lut_full_csv <- NULL
          }

          # XLSX if requested
          group_lut_xlsx <- NULL
          if (isTRUE(write_group_lut_xlsx)) {
            if (!requireNamespace("openxlsx", quietly = TRUE)) {
              stop("write_group_lut_xlsx=TRUE requires package 'openxlsx'. Install with install.packages('openxlsx').")
            }

            if (is.null(lut_xlsx_name) || !nzchar(lut_xlsx_name)) {
              group_lut_xlsx <- file.path(lut_dir, paste0(prefix, "_", lut_file_tag, ".xlsx"))
            } else {
              group_lut_xlsx <- file.path(lut_dir, lut_xlsx_name)
            }

            safe_overwrite(group_lut_xlsx)

            sheets <- list(lut_min = group_lut$lut_min)
            if (isTRUE(lut_write_full_csv)) {
              sheets$lut_full <- group_lut$lut_full
            } else {
              # Even if you skipped full CSV, you may still want full in XLSX
              sheets$lut_full <- group_lut$lut_full
            }

            openxlsx::write.xlsx(sheets, file = group_lut_xlsx, overwrite = TRUE)
            msg("[LUT] Written: ", group_lut_xlsx)
          }

          # expose paths to outer scope (so you can store in results/out)
          # (make sure these variables exist outside this if-block in your function)
        }
      }
    }

    # Build simple "original code" LUT (optional) for burnable_classes when burnable_source == "original"
    orig_lut <- NULL
    if (!is.null(cor_tbl)) {
      # key depends on corine_value_col
      cor_tbl$KEY <- if (corine_value_col == "GRID_CODE") as.integer(cor_tbl$GRID_CODE) else as.integer(cor_tbl$CLC_CODE)
      # choose a readable label (prefer LABEL3, then LABEL2, then LABEL1)
      lab0 <- cor_tbl$LABEL3
      lab0[is.na(lab0) | !nzchar(trimws(lab0))] <- cor_tbl$LABEL2[is.na(lab0) | !nzchar(trimws(lab0))]
      lab0[is.na(lab0) | !nzchar(trimws(lab0))] <- cor_tbl$LABEL1[is.na(lab0) | !nzchar(trimws(lab0))]
      orig_lut <- unique(data.frame(
        KEY = as.integer(cor_tbl$KEY),
        clc = as.integer(cor_tbl$CLC_CODE),
        grid = as.integer(cor_tbl$GRID_CODE),
        lab = as.character(lab0),
        stringsAsFactors = FALSE
      ))
      orig_lut$label <- paste0("C", orig_lut$KEY, ": ", orig_lut$lab)
    }

    # ---- vector extension based on driver
    vec_ext_from_driver <- function(driver) {
      d <- toupper(driver)
      if (d %in% c("ESRI SHAPEFILE")) return("shp")
      if (d %in% c("GPKG")) return("gpkg")
      if (d %in% c("GEOJSON")) return("geojson")
      # fallback
      "shp"
    }

    # ---- vectorize helper (with optional LUT join)
    vectorize_raster <- function(r_path, out_vec_path, method, driver,
                                 lut_df = NULL, vec_key = "DN", lut_key = NULL) {
      method <- match.arg(method)

      if (method == "terra") {
        r <- terra::rast(r_path)
        p <- terra::as.polygons(r, dissolve = TRUE, values = TRUE, na.rm = TRUE)
        if (is.null(p) || terra::nrow(p) == 0) {
          msg("[VECT] No polygons generated for: ", r_path)
          return(invisible(NULL))
        }

        # normalize the value field to DN
        nm0 <- names(p)
        if (length(nm0) >= 1) names(p)[1] <- vec_key

        p_sf <- sf::st_as_sf(p)

        if (!is.null(lut_df) && !is.null(lut_key)) {
          if (!vec_key %in% names(p_sf)) stop("Vector key column not found: ", vec_key)
          p_sf[[vec_key]] <- suppressWarnings(as.integer(p_sf[[vec_key]]))
          lut_df[[lut_key]] <- suppressWarnings(as.integer(lut_df[[lut_key]]))

          p_sf <- merge(p_sf, lut_df, by.x = vec_key, by.y = lut_key, all.x = TRUE, sort = FALSE)
        }

        safe_overwrite(out_vec_path)
        sf::st_write(p_sf, out_vec_path, driver = driver, quiet = TRUE)
        return(invisible(out_vec_path))
      }

      # gdal_polygonize branch
      if (is.null(python_exe) || !nzchar(python_exe) || !file.exists(python_exe)) {
        stop("vectorize_method='gdal_polygonize' requires python_exe (existing file).")
      }
      if (is.null(gdal_polygonize_script) || !nzchar(gdal_polygonize_script) || !file.exists(gdal_polygonize_script)) {
        stop("vectorize_method='gdal_polygonize' requires gdal_polygonize_script (existing file).")
      }

      # build output first
      safe_overwrite(out_vec_path)
      layer_name <- tools::file_path_sans_ext(basename(out_vec_path))

      args <- c(
        shQuote(gdal_polygonize_script),
        if (isTRUE(touches)) "-8" else NULL,
        shQuote(r_path),
        "-f", shQuote(driver),
        shQuote(out_vec_path),
        shQuote(layer_name),
        vec_key
      )

      out <- suppressWarnings(system2(python_exe, args = args, stdout = TRUE, stderr = TRUE))
      if (!file.exists(out_vec_path)) {
        stop("gdal_polygonize failed.\nLast lines:\n", paste(utils::tail(out, 60), collapse = "\n"))
      }

      if (!is.null(lut_df) && !is.null(lut_key)) {
        # read, join, rewrite
        v <- sf::st_read(out_vec_path, quiet = TRUE)
        if (!vec_key %in% names(v)) stop("Polygonized output has no key field: ", vec_key)

        v[[vec_key]] <- suppressWarnings(as.integer(v[[vec_key]]))
        lut_df[[lut_key]] <- suppressWarnings(as.integer(lut_df[[lut_key]]))

        v2 <- merge(v, lut_df, by.x = vec_key, by.y = lut_key, all.x = TRUE, sort = FALSE)
        safe_overwrite(out_vec_path)
        sf::st_write(v2, out_vec_path, driver = driver, quiet = TRUE)
      }

      invisible(out_vec_path)
    }

    results <- list()

    for (nm in names(corine_rasters)) {

      in_path <- corine_rasters[[nm]]
      if (!file.exists(in_path)) stop("CORINE raster not found: ", in_path)

      forced_base <- resolve_out_base(nm, out_base_names)
      base_out <- if (!is.null(forced_base)) forced_base else make_out_base(nm, prefix, dedup = dedup_prefix)

      msg("========================================")
      msg("[IN] ", nm, " -> ", in_path)
      msg("[OUTBASE] ", base_out)

      cor <- terra::rast(in_path)
      if (terra::nlyr(cor) > 1) {
        msg("[WARN] Input has multiple layers; using the first layer only.")
        cor <- cor[[1]]
      }

      cor_wkt <- terra::crs(cor, proj = TRUE)
      if (is.na(cor_wkt) || !nzchar(cor_wkt)) stop("CORINE raster CRS missing/empty: ", in_path)

      pen_sf <- sf::st_transform(pen_sf0, sf::st_crs(cor_wkt))
      pen_v  <- terra::vect(pen_sf)

      msg("[CLIP] Crop + mask by peninsula...")
      cor2 <- terra::crop(cor, pen_v)
      cor2 <- terra::mask(cor2, pen_v)

      if (!is.null(tmpl)) {
        msg("[ALIGN] Projecting to template grid (near)...")
        cor_aligned <- terra::project(cor2, tmpl, method = "near")
      } else {
        msg("[ALIGN] Projecting to target CRS/res (near)...")
        cor_aligned <- terra::project(cor2, target_crs, method = "near", res = target_res)
      }

      out <- list(
        input = in_path,
        group_lut_min_csv  = group_lut_min_csv,
        group_lut_full_csv = group_lut_full_csv,
        group_lut_xlsx     = group_lut_xlsx
      )

      # ---- build group raster in memory if needed
      need_burnable <- isTRUE(write_burnable_classes) || isTRUE(write_binary_mask)
      need_group_for_burnable <- isTRUE(need_burnable) && burnable_source == "group"

      group_r <- NULL
      group_path <- NULL

      if (!is.null(group_rcl_use) && (isTRUE(write_group_raster) || isTRUE(need_group_for_burnable))) {
        msg("[RECLASS] Building grouped raster...")
        group_r <- terra::classify(cor_aligned, rcl = group_rcl_use, others = NA)

        if (isTRUE(write_group_raster)) {
          group_path <- file.path(output_raster_dir, paste0(base_out, "_reclass.tif"))
          safe_overwrite(group_path)
          terra::writeRaster(
            group_r, group_path, overwrite = TRUE,
            datatype = "INT2U", NAflag = 0,
            gdal = c(paste0("COMPRESS=", compress), "TILED=YES", "BIGTIFF=IF_SAFER")
          )
          out$group_raster <- group_path
        }
      }

      # ---- burnable raster
      burnable_r <- NULL
      burnable_classes_path <- NULL

      if (isTRUE(need_burnable)) {

        if (is.null(burnable_gridcodes)) {
          stop("To compute burnable/binary outputs you must provide burnable_gridcodes or burnable_matrix_01.")
        }

        if (burnable_source == "group") {
          if (is.null(group_rcl_use)) stop("burnable_source='group' requires group_rcl_matrix.")
          if (is.null(group_r)) {
            msg("[RECLASS] Building grouped raster (needed for burnable_source='group')...")
            group_r <- terra::classify(cor_aligned, rcl = group_rcl_use, others = NA)
          }
          msg("[BURNABLE] Selecting burnable classes from GROUP raster...")
          src_r <- group_r
        } else {
          msg("[BURNABLE] Selecting burnable classes from ORIGINAL raster...")
          src_r <- cor_aligned
        }

        # IMPORTANT: avoid %in% with SpatRaster; use classify()
        burnable_r <- terra::classify(src_r, rcl = cbind(burnable_gridcodes, burnable_gridcodes), others = NA)

        if (isTRUE(write_burnable_classes)) {
          burnable_classes_path <- file.path(output_raster_dir, paste0(base_out, "_burnable_classes.tif"))
          safe_overwrite(burnable_classes_path)
          terra::writeRaster(
            burnable_r, burnable_classes_path, overwrite = TRUE,
            datatype = "INT2U", NAflag = 0,
            gdal = c(paste0("COMPRESS=", compress), "TILED=YES", "BIGTIFF=IF_SAFER")
          )
          out$burnable_classes_raster <- burnable_classes_path

          # quick diagnostic
          fv <- terra::freq(burnable_r, digits = 0)
          if (!is.null(fv) && nrow(as.data.frame(fv)) > 0) {
            msg("[BURNABLE] Values present: ", paste(as.data.frame(fv)$value, collapse = ", "))
          } else {
            msg("[BURNABLE][WARN] No values present in burnable raster (all NA).")
          }
        }
      }

      # ---- binary mask
      binary_path <- NULL
      if (isTRUE(write_binary_mask)) {
        if (is.null(burnable_r)) stop("Internal error: burnable_r is NULL but write_binary_mask=TRUE.")
        msg("[BINARY] Creating binary mask (1 for burnable, NA otherwise)...")
        bin <- terra::ifel(!is.na(burnable_r), 1, NA)
        binary_path <- file.path(output_raster_dir, paste0(base_out, "_burnable_mask_binary.tif"))
        safe_overwrite(binary_path)
        terra::writeRaster(
          bin, binary_path, overwrite = TRUE,
          datatype = "INT1U", NAflag = 0,
          gdal = c(paste0("COMPRESS=", compress), "TILED=YES", "BIGTIFF=IF_SAFER")
        )
        out$binary_mask_raster <- binary_path
      }

      # ---- vectorize (optional)
      if (isTRUE(vectorize)) {
        if (is.null(output_vector_dir)) stop("vectorize=TRUE requires output_vector_dir.")
        vmethod <- match.arg(vectorize_method)
        which <- unique(vectorize_which)
        ext_vec <- vec_ext_from_driver(vector_driver)

        # Decide LUT to apply depending on what is being vectorized
        # - group_raster: group LUT (id -> label)
        # - burnable_classes: group LUT if burnable_source == "group"; otherwise original LUT (KEY -> label)
        lut_for_group <- NULL
        lut_key_group <- NULL
        if (isTRUE(apply_lut_to_vectors) && !is.null(group_lut)) {
          lut_for_group <- group_lut$lut_full
          lut_key_group <- "id"
        }

        lut_for_burnable <- NULL
        lut_key_burnable <- NULL
        if (isTRUE(apply_lut_to_vectors)) {
          if (burnable_source == "group" && !is.null(group_lut)) {
            lut_for_burnable <- group_lut$lut_full
            lut_key_burnable <- "id"
          } else if (burnable_source == "original" && !is.null(orig_lut)) {
            # keep only compact fields to avoid shapefile field issues
            lut_for_burnable <- orig_lut[, c("KEY","clc","grid","lab","label"), drop = FALSE]
            lut_key_burnable <- "KEY"
          }
        }

        if ("group_raster" %in% which && !is.null(group_path) && file.exists(group_path)) {
          out_vec <- file.path(output_vector_dir, paste0(base_out, "_group.", ext_vec))
          msg("[VECT] Group raster -> ", out_vec)
          vectorize_raster(
            r_path = group_path,
            out_vec_path = out_vec,
            method = vmethod,
            driver = vector_driver,
            lut_df = lut_for_group,
            vec_key = "DN",
            lut_key = lut_key_group
          )
          out$group_vector <- out_vec
        }

        if ("burnable_classes" %in% which && !is.null(burnable_classes_path) && file.exists(burnable_classes_path)) {
          out_vec <- file.path(output_vector_dir, paste0(base_out, "_burnable_classes.", ext_vec))
          msg("[VECT] Burnable classes -> ", out_vec)
          vectorize_raster(
            r_path = burnable_classes_path,
            out_vec_path = out_vec,
            method = vmethod,
            driver = vector_driver,
            lut_df = lut_for_burnable,
            vec_key = "DN",
            lut_key = lut_key_burnable
          )
          out$burnable_classes_vector <- out_vec
        }

        if ("binary_mask" %in% which && !is.null(binary_path) && file.exists(binary_path)) {
          out_vec <- file.path(output_vector_dir, paste0(base_out, "_burnable_mask_binary.", ext_vec))
          msg("[VECT] Binary mask -> ", out_vec)
          vectorize_raster(
            r_path = binary_path,
            out_vec_path = out_vec,
            method = vmethod,
            driver = vector_driver,
            lut_df = NULL,
            vec_key = "DN",
            lut_key = NULL
          )
          out$binary_mask_vector <- out_vec
        }
      }

      results[[nm]] <- out
      msg("[DONE] ", nm)
    }

    invisible(results)
  }
