#' Flag Burned Polygons Based on Regeneration Overlap
#'
#' @description
#' This function assigns regeneration flags to burned area polygons based on their spatial
#' overlap with post-fire regeneration polygons (e.g., from years 1, 2, or more after fire).
#' A polygon is flagged as `"regenera"` if the total intersected area with regeneration polygons
#' from a specific post-fire year exceeds a user-defined minimum threshold (`min_regen_ratio`)
#' relative to its own area.
#'
#' Filenames of the regeneration layers must contain `"P1"`, `"P2"`, etc., to indicate
#' the number of years after fire. These period labels are automatically extracted and used
#' to name the output columns (e.g., `regen_ratio_P1`, `regenera_flag_P1`).
#'
#' For each detected regeneration year, the function adds:
#' \itemize{
#'   \item `regen_ratio_Px`: proportion of the burned polygon area overlapping with regeneration polygons.
#'   \item `regenera_flag_Px`: `"regenera"` if the threshold is met, `"no_regenera"` otherwise.
#' }
#' A final column `regenera_flag_all` is also created, indicating `"regenera"` only if
#' all selected years in `remove_condition` meet their respective thresholds.
#'
#' If `replace_by_P1 = TRUE`, the function will optionally replace burned polygons in the
#' filtered output by overlapping P1 regeneration polygons that are larger in area. Optionally,
#' a maximum area ratio (`max_area_ratio_p1`) can be defined to restrict replacements to only
#' those P1 polygons not exceeding a given multiple of the burned area.
#'
#' You can also choose to export polygons flagged as `"no_regenera"` using `save_no_regenera`.
#' These can be optionally filtered by a minimum area threshold (`min_area_no_regenera`).
#'
#' If `validate_geometries = TRUE`, geometries are checked and corrected using `sf::st_make_valid()`,
#' which is more robust but significantly slower. Set to `FALSE` by default for performance.
#'
#' @name flag_otsu_regenera
#' @rdname flag_otsu_regenera
#'
#' @param burned_files Character vector with paths to burned area shapefiles (.shp).
#' @param regenera_files Character vector with paths to regeneration shapefiles (.shp).
#'   Filenames must contain `"P1"`, `"P2"`, etc., to indicate the year after fire.
#' @param min_regen_ratio Named numeric vector with minimum area ratio thresholds (between 0 and 1)
#'   for each regeneration year. Names must correspond to period labels (e.g., `c(P1 = 0.05, P2 = 0.10)`).
#'  These values are used:
#' (1) when no classwise_regen_thresholds are provided,
#' (2) or as fallback when a polygon has no match or NA thresholds.
#' @param remove_no_regenera Logical. If `TRUE`, polygons that do not meet all specified
#'   thresholds (see `remove_condition`) will be excluded from the filtered output.
#' @param remove_condition Character vector indicating the post-fire years (e.g., `"P1"`, `"P2"`, `"P3"`)
#'   that must meet their thresholds in `min_regen_ratio` in order for a polygon to be retained.
#'   Used to identify false fire detections due to lack of regeneration.
#' @param output_dir Optional output directory. If `NULL`, output files are saved in the same
#'   folder as the burned shapefiles.
#' @param replace_by_P1 Logical. If `TRUE` (default), filtered burned polygons will be replaced by
#'   overlapping P1 regeneration polygons when the latter are larger in area.
#' @param max_area_ratio_p1 Optional. Maximum allowed ratio between the area of a P1 regeneration polygon and the intersected burned polygon.
#'   If specified, only regenerated polygons with area greater than the burned area and not exceeding \code{max_area_ratio_p1} times that area
#'   will be used to replace burned polygons. If \code{NULL} (default), only the condition \code{area_ha > burn_area_ha} is applied.
#' @param save_no_regenera Character. Options are `"all"` (save all `"no_regenera"` polygons),
#'   `"area_filter"` (save only those with area ? `min_area_no_regenera`), or `"none"` (do not save them).
#' @param min_area_no_regenera Numeric. Minimum area in hectares to retain a `"no_regenera"` polygon when `save_no_regenera = "area_filter"`.
#'   Default is `0.1` ha.
#' @param validate_geometries Logical. If `TRUE`, applies `st_make_valid()` to correct invalid geometries
#'   before processing. Default is `FALSE` for faster performance.
#' @param output_format Character. Output format for saved files: `"shp"` for ESRI Shapefile (default)
#'   or `"geojson"` for GeoJSON.
#'  @param classwise_regen_thresholds Optional `data.frame` with per-class regeneration thresholds.
#' Must include:
#' \itemize{
#'   \item \code{class_field}: field name in the burned polygons (e.g., `"CORINE_CLA"`, `"eco_id"`).
#'   \item \code{class_value}: corresponding class value.
#'   \item One or more columns named \code{P1}, \code{P2}, etc., with numeric thresholds.
#' }
#' If provided, these thresholds will override \code{min_regen_ratio} for polygons matching each class.
#' Polygons without a match will fall back to the values in \code{min_regen_ratio}.


#'
#' @return Saves up to four shapefiles per burned area file:
#' \describe{
#'   \item{Main shapefile (unfiltered)}{All polygons with `regen_ratio_Px`, `regenera_flag_Px`, and `regenera_flag_all`.}
#'   \item{Filtered shapefile}{Only polygons that satisfy all thresholds defined in `remove_condition`. File ends in `_filter.shp`.}
#'   \item{Final shapefile with replacements}{(if `replace_by_P1 = TRUE`) Filtered polygons with some replaced by larger P1 regeneration polygons. File ends in `_new.shp`.}
#'   \item{Shapefile of no_regenera polygons}{(if `save_no_regenera != "none"`) Contains polygons not meeting regeneration thresholds. File ends in `_no_selected.shp`.}
#' }
#'
#' @examples
#' \dontrun{
#' # Apply per-period regeneration thresholds and filter non-regenerating polygons
#' burned_files <- list.files("data/burned", pattern = "\\.shp$", full.names = TRUE)
#' regenera_files <- list.files("data/regenera", pattern = "\\.shp$", full.names = TRUE)
#'
#' flag_otsu_regenera(
#'   burned_files = burned_files,
#'   regenera_files = regenera_files,
#'   min_regen_ratio = c(P1 = 0.05, P2 = 0.10),
#'   remove_no_regenera = TRUE,
#'   remove_condition = c("P1", "P2"),
#'   output_dir = "results/",
#'   replace_by_P1 = TRUE,
#'   max_area_ratio_p1 = 3,
#'   save_no_regenera = "area_filter",
#'   min_area_no_regenera = 100,
#'   validate_geometries = FALSE,
#'   output_format = c("geojson")
#' )
#'
#' # Define class-specific regeneration thresholds
#' class_thresholds <- data.frame(
#'   class_field = "CORINE_CLA",
#'   class_value = c(1, 2),
#'   P1 = c(0.15, 0.15),
#'   P2 = c(0.0, 0.05)
#' )
#'
#' # Run flag_otsu_regenera1 using classwise thresholds and fallback values
#' result <- flag_otsu_regenera(
#'   burned_files = burned_list,
#'   regenera_files = regenera_list,
#'   min_regen_ratio = c(P1 = 0.05, P2 = 0.15),  # fallback for classes not in the table
#'   classwise_regen_thresholds = class_thresholds,
#'   group_field = "CORINE_CLA",
#'   remove_no_regenera = TRUE,
#'   remove_condition = c("P1", "P2"),
#'   replace_by_P1 = TRUE,
#'   max_area_ratio_p1 = 3,
#'   save_no_regenera = TRUE,
#'   min_area_no_regenera = 500,
#'   validate_geometries = FALSE,
#'   output_format = "shp",
#'   output_dir = burned_dir
#' )
#' }
#'
#' @importFrom sf st_read st_write st_crs st_transform st_area st_intersection st_set_geometry st_is_valid st_make_valid st_union st_cast st_as_sf
#' @importFrom dplyr mutate row_number group_by summarise select left_join filter first all_of where
#' @importFrom tidyr replace_na
#' @importFrom tools file_path_sans_ext
#' @importFrom utils glob2rx
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @encoding UTF-8
#' @export


utils::globalVariables(c(
  ".", "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend", "mnbbx_wd",
  "mtext", "na.omit", "ornt_bbox", "p_w_ratio", "par", "png",
  "regen_area_ha", "regen_area_ha_sum", "setNames", "total_burned_area",
  "year", "%>%", ":="
))


flag_otsu_regenera <- function(burned_files,
                                regenera_files,
                                min_regen_ratio = c(P1 = 0.05, P2 = 0.05),
                                remove_no_regenera = FALSE,
                                remove_condition = c("P1", "P2"),
                                output_dir = NULL,
                                replace_by_P1 = TRUE,
                                max_area_ratio_p1 = NULL,
                                save_no_regenera = "none",
                                min_area_no_regenera = 0.1,
                                classwise_regen_thresholds = NULL,
                                 group_field = NULL,
                                validate_geometries = FALSE,
                                output_format = c("shp", "geojson")) {

  # Validar el argumento de formato
  output_format <- match.arg(output_format, choices = c("shp", "geojson"))

  burned_files <- as.character(burned_files)

  # Validaciones iniciales
  if (is.null(names(min_regen_ratio))) {
    stop("min_regen_ratio debe tener nombres (por ejemplo: c(P1 = 0.05, P2 = 0.1))")
  }

  # Detectar los periodos realmente presentes en los archivos de regeneracion
  periods_detected <- unique(regmatches(regenera_files, regexpr("P[0-9]+", regenera_files)))
  periods_detected <- periods_detected[!is.na(periods_detected)]

  # Filtrar remove_condition por los periodos detectados
  remove_condition <- intersect(remove_condition, periods_detected)

  missing_periods <- setdiff(remove_condition, names(min_regen_ratio))
  if (length(missing_periods) > 0) {
    stop("Faltan valores de min_regen_ratio para los siguientes periodos: ",
         paste(missing_periods, collapse = ", "))
  }


  for (burned_file in burned_files) {


    delete_shapefile <- function(shp_path) {
      base <- tools::file_path_sans_ext(shp_path)
      exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbx", ".sbn", ".qix")
      files <- paste0(base, exts)
      files_exist <- files[file.exists(files)]
      if (length(files_exist) > 0) {
        invisible(sapply(files_exist, function(f) tryCatch(file.remove(f), error = function(e) NULL)))
      }
    }


    message("Processing: ", burned_file)

    burned_sf <- sf::st_read(burned_file, quiet = TRUE)
    if (!all(sf::st_is_valid(burned_sf)))
      burned_sf <- sf::st_make_valid(burned_sf)


    # Restaurar 'burned_id' si fue truncado o renombrado
    if (!"burned_id" %in% names(burned_sf)) {
      burned_candidates <- grep("^uid$|burn.*id", names(burned_sf), ignore.case = TRUE, value = TRUE)
      if (length(burned_candidates) == 1) {
        names(burned_sf)[names(burned_sf) == burned_candidates] <- "burned_id"
        message("Restored 'burned_id' column from field: ", burned_candidates)
      } else {
        burned_sf <- dplyr::mutate(burned_sf, burned_id = dplyr::row_number())
        message("Created new 'burned_id' column.")
      }
    }

       burned_sf$burn_area_ha <- as.numeric(sf::st_area(burned_sf)) / 10000


    # Si se define group_field pero no existe en burned_sf
    if (!is.null(group_field) && !(group_field %in% names(burned_sf))) {

      # Buscar columna cuyo nombre contenga "CORINE" (insensible a mayusculas)
      corine_candidates <- grep("CORINE", names(burned_sf), value = TRUE, ignore.case = TRUE)

      if (length(corine_candidates) == 1) {
        message("Field '", group_field, "' not found. Using detected CORINE-like field: ", corine_candidates)
        group_field <- corine_candidates  # sustituir dinamicamente
      } else if (length(corine_candidates) > 1) {
        stop("Multiple columns containing 'CORINE' found: ", paste(corine_candidates, collapse = ", "),
             ". Please specify 'group_field' explicitly.")
      } else {
        stop("The column specified in 'group_field' does not exist in burned_sf, and no 'CORINE' field was found.")
      }
    }

    regen_labels <- c()

    regen_thresh_match <- regmatches(regenera_files, regexpr("thresh[0-9]+", regenera_files))
    regen_thresh <- unique(gsub("thresh", "", regen_thresh_match))
    regen_thresh <- if (length(regen_thresh) == 1)
      regen_thresh else "NA"

    for (regen_file in regenera_files) {
      regenera_sf <- sf::st_read(regen_file, quiet = TRUE)

      if (!sf::st_crs(burned_sf) == sf::st_crs(regenera_sf)) {
        regenera_sf <- sf::st_transform(regenera_sf, sf::st_crs(burned_sf))
      }

      label <- regmatches(regen_file, regexpr("P[0-9]+", regen_file))

      if (length(label) > 0 && !is.na(label)) {

        regen_labels <- union(regen_labels, label)
        ratio_col <- paste0("regen_ratio_", label)

        # Definir group_col primero
        group_col <- if (!is.null(group_field)) group_field else "burned_id"

        # Asegurar que burned_id existe
        if (!"burned_id" %in% names(burned_sf)) {
          burned_candidates <- grep("^uid$|burn.*id", names(burned_sf), ignore.case = TRUE, value = TRUE)
          if (length(burned_candidates) == 1) {
            names(burned_sf)[names(burned_sf) == burned_candidates] <- "burned_id"
            message("Restored 'burned_id' column from field: ", burned_candidates)
          } else {
            burned_sf <- dplyr::mutate(burned_sf, burned_id = seq_len(nrow(burned_sf)))
            message("Generated new 'burned_id' column.")
          }
        }

        # Verificar que group_col existe
        if (!(group_col %in% names(burned_sf))) {
          stop("The column specified in 'group_field' does not exist in burned_sf: ", group_col)
        }

        # Campos a mantener
        fields_to_keep <- intersect(c(group_col, "burned_id", "CORINE_CLA"), names(burned_sf))
        fields_to_keep <- setdiff(fields_to_keep, group_col)  # <- importante para evitar duplicados

        if (length(fields_to_keep) == 0) {
          stop("No valid fields found in 'fields_to_keep'.")
        }

        if (group_col == "burned_id") {
          message("Dissolving burned polygons by 'burned_id'")
          burned_sf <- burned_sf |>
            dplyr::group_by(.data[[group_col]]) |>
            dplyr::summarise(dplyr::across(all_of(fields_to_keep), dplyr::first), .groups = "drop")
        } else {
          message("Using individual polygons (no dissolve). Using '", group_col, "'")
        }

        burned_sf$burn_area_ha <- as.numeric(sf::st_area(burned_sf)) / 10000

        intersec <- suppressWarnings(sf::st_intersection(
          sf::st_make_valid(burned_sf),
          sf::st_make_valid(regenera_sf)
        ))

        if (nrow(intersec) > 0) {
          intersec$int_area_ha <- as.numeric(sf::st_area(intersec)) / 10000

          area_summary <- intersec |>
            sf::st_set_geometry(NULL) |>
            dplyr::group_by(burned_id) |>
            dplyr::summarise(regen_area_ha_sum = sum(int_area_ha, na.rm = TRUE))

          burned_sf <- burned_sf |>
            dplyr::left_join(area_summary, by = "burned_id") |>
            dplyr::mutate(
              regen_area_ha_sum = tidyr::replace_na(regen_area_ha_sum, 0),
              !!ratio_col := pmin(regen_area_ha_sum, burn_area_ha) / burn_area_ha
            ) |>
            dplyr::select(-regen_area_ha_sum)

          message("Added column '", ratio_col, "' with regeneration ratios for ", label)
        } else {
          burned_sf[[ratio_col]] <- 0
          message("No intersections found with regeneration polygons for ", label, ". Set ", ratio_col, " to 0.")
        }
      }
    }

    regen_labels_str <- paste(regen_labels, collapse = "")

    message("Detected regeneration years: ",
            paste(regen_labels, collapse = ", "))

    burned_copy <- burned_sf
    flag_names <- c()
    meets_all_thresholds <- rep(TRUE, nrow(burned_copy))  # inicializar vector conjunto


    if (!is.null(classwise_regen_thresholds)) {
        class_field <- unique(classwise_regen_thresholds$class_field)
        message("class_field '", class_field)

        if (length(class_field) != 1 || !(class_field %in% names(burned_sf))) {
          # Detectar columnas que contengan "CORINE" como fallback
          corine_candidates <- grep("CORINE", names(burned_sf), value = TRUE, ignore.case = TRUE)

          if (length(corine_candidates) == 1) {
            message("class_field '", class_field, "' not found. Using detected CORINE-like field: ", corine_candidates)
            class_field <- corine_candidates
            # Tambien actualizar en classwise_regen_thresholds
            classwise_regen_thresholds$class_field <- class_field
          } else {
            stop("Invalid or missing 'class_field' column in burned polygons and no CORINE-like field found.")
          }
        }

        message("class_field: ", class_field,
                " | burned_copy class: ", class(burned_sf[[class_field]]),
                " | thresholds class: ", class(classwise_regen_thresholds$class_value))
      }


    # Crear flags por periodo y acumular condiciones conjuntas
    for (label in names(min_regen_ratio)) {
      ratio_col <- paste0("regen_ratio_", label)
      flag_col  <- paste0("regenera_flag_", label)


      if (!is.null(classwise_regen_thresholds)) {
        # Determine class field (e.g., CORINE_CLA)
        class_field <- unique(classwise_regen_thresholds$class_field)
        if (length(class_field) != 1 || !(class_field %in% names(burned_copy))) {
          stop("Invalid or missing 'class_field' column in burned polygons.")
        }

        # Coerce both to character to avoid join type mismatch
        # Asegurar tipos
        burned_copy[[class_field]] <- as.character(burned_copy[[class_field]])
        thresholds_df <- classwise_regen_thresholds
        thresholds_df$class_value <- as.character(thresholds_df$class_value)

        # Si la columna del label (P1 o P2) ya existe, eliminala para forzar recreacion limpia
        if (label %in% names(burned_copy)) {
          burned_copy[[label]] <- NULL
        }

        # Crear vector con valor general para todas las filas
        burned_copy[[label]] <- min_regen_ratio[[label]]

        # Extraer solo columnas necesarias del df de thresholds
        thresh_df <- thresholds_df[, c("class_value", label)]
        colnames(thresh_df)[2] <- "special_thresh"

        # Join por clase
        burned_copy <- dplyr::left_join(
          burned_copy,
          thresh_df,
          by = setNames("class_value", class_field)
        )

        # Sobrescribir donde hay valores especiales
        burned_copy[[label]] <- ifelse(!is.na(burned_copy$special_thresh),
                                       burned_copy$special_thresh,
                                       burned_copy[[label]])

        # Eliminar columna auxiliar
        burned_copy$special_thresh <- NULL


        # Inicializa vector de thresholds con valores globales por defecto
        threshold_vec <- rep(min_regen_ratio[[label]], nrow(burned_copy))

        if (!is.null(classwise_regen_thresholds)) {
          # Determina campo de clase
          class_field <- unique(classwise_regen_thresholds$class_field)
          if (length(class_field) != 1 || !(class_field %in% names(burned_copy))) {
            stop("Invalid or missing 'class_field' column in burned polygons.")
          }

          # Forzar tipos compatibles
          burned_copy[[class_field]] <- as.character(burned_copy[[class_field]])
          thresholds_df <- classwise_regen_thresholds
          thresholds_df$class_value <- as.character(thresholds_df$class_value)

          # Join temporal para solo clases definidas
          join_df <- burned_copy[, c("burned_id", class_field)]
          join_df <- dplyr::left_join(join_df, thresholds_df[, c("class_value", label)],
                                      by = setNames("class_value", class_field))

          # Sustituye solo si hay valores definidos para esa clase
          defined_rows <- !is.na(join_df[[label]])
          threshold_vec[defined_rows] <- join_df[[label]][defined_rows]
        }

      } else {
        threshold_vec <- rep(min_regen_ratio[[label]], nrow(burned_copy))
      }


      if (ratio_col %in% names(burned_copy)) {
        # Calcular si cumple el umbral para este periodo
        current_check <- burned_copy[[ratio_col]] >= threshold_vec
        current_check[is.na(current_check)] <- FALSE  # tratar NAs como "no_regenera"

        # Crear flag individual
        burned_copy[[flag_col]] <- ifelse(current_check, "regenera", "no_regenera")
        flag_names <- c(flag_names, flag_col)

        # Actualizar condicion conjunta
        meets_all_thresholds <- meets_all_thresholds & current_check

        message(sprintf("Created '%s' using threshold(s) from '%s' (column: %s)",
                        flag_col,
                        if (!is.null(classwise_regen_thresholds)) class_field else "global",
                        ratio_col))
      } else {
        burned_copy[[flag_col]] <- NA_character_
        meets_all_thresholds <- meets_all_thresholds & FALSE
        warning(sprintf("Column '%s' is missing. Assigned NA to '%s'.",
                        ratio_col, flag_col))
      }
    }

    # Crear flag final conjunta
    burned_copy$regenera_flag_all <- ifelse(meets_all_thresholds, "regenera", "no_regenera")
    message("Create'regenera_flag_all' combining all defined periods.")



    # ---- 1. GUARDAR SHAPEFILE COMPLETO (con todas las banderas) ----

    # Base name of the burned file
    burned_base <- tools::file_path_sans_ext(basename(burned_file))

    # Extract threshold suffix (e.g., "thresh50")
    regen_thresh <- stringr::str_extract(basename(regenera_files[1]), "thresh\\d+")
    regen_thresh_val <- sub("thresh", "", regen_thresh)

    # Detected regeneration periods (e.g., P1P2)
    regen_labels_str <- paste(regen_labels, collapse = "")

    # Build threshold suffix
    rat_suffix <- if (!is.null(classwise_regen_thresholds)) {
      # Include class field and values (e.g., "_classCORINE_CLA_1-2")
      class_field <- unique(classwise_regen_thresholds$class_field)
      class_vals <- unique(classwise_regen_thresholds$class_value)
      class_vals_str <- paste(sort(class_vals), collapse = "-")
      paste0(class_field, "", class_vals_str)
    } else {
      paste(
        paste0("rat", sprintf("%02d", round(100 * min_regen_ratio[names(min_regen_ratio)]))),
        collapse = "_"
      )
    }

    # Optional suffix for max P1 replacement constraint
    area_ratio_suffix <- if (!is.null(max_area_ratio_p1)) {
      paste0("_ratP1_", formatC(max_area_ratio_p1, format = "f", digits = 0))
    } else {
      ""
    }

    # Optional suffix for grouping
    group_suffix <- if (!is.null(group_field)) group_field else ""

    #### SAVING ALL POLYGONS WITH REGENERA DATA ##########

    # Determine extension
    ext <- if (output_format == "geojson") ".geojson" else ".shp"

    # Compose full filename
    full_out_name <- paste0(
      burned_base,
      "_reg_thresh", regen_thresh_val,
      "_", regen_labels_str,
      "_", rat_suffix,
      area_ratio_suffix,
      ext
    )

    # Define output path
    output_dir <- if (!is.null(output_dir)) output_dir else dirname(burned_file)
    full_output_file <- file.path(output_dir, full_out_name)

    # Delete previous output if exists
    if (file.exists(full_output_file)) {
      file.remove(full_output_file)
    }

    message("Writing the complete shapefile (unfiltered): ", full_output_file)
    suppressWarnings(sf::st_write(burned_copy, full_output_file, delete_layer = TRUE, quiet = TRUE))



    # ---- 2. FILTER BY FLAG ----

    # Determine which rows to keep
    keep_rows <- rep(TRUE, nrow(burned_copy))
    if (isTRUE(remove_no_regenera)) {
      if (!"regenera_flag_all" %in% names(burned_copy)) {
        stop("Missing 'regenera_flag_all'. Make sure flags were generated correctly.")
      }
      keep_rows <- burned_copy$regenera_flag_all == "regenera"
    }

    # Drop burned_id (internal) and filter rows
    filtered_output <- burned_copy[keep_rows, ] |>
      dplyr::select(-burned_id)

    # Build base name again (exclude extension)
    base_name <- tools::file_path_sans_ext(full_out_name)

    # Adjust filtered output filename
    filtered_out_name <- paste0(base_name, "_filter", ext)
    filtered_output_file <- file.path(output_dir, filtered_out_name)

    # Remove previous file if it exists
    if (file.exists(filtered_output_file)) {
      file.remove(filtered_output_file)
    }

    # Write filtered output
    message("Saving filtered shapefile: ", filtered_output_file)
    message("Number of polygons after filtering: ", nrow(filtered_output),
            " (removed: ", nrow(burned_copy) - nrow(filtered_output), ")")

    suppressWarnings(sf::st_write(filtered_output, filtered_output_file, delete_layer = TRUE, quiet = TRUE))


    # ---- 2. APLICAR FILTRADO (NO REGENERA) ----

    keep_rows1 <- rep(TRUE, nrow(burned_copy))
    if (isTRUE(remove_no_regenera)) {
      if (!"regenera_flag_all" %in% names(burned_copy)) {
        stop("No se encontro 'regenera_flag_all'. Verifica que hayas generado los flags correctamente.")
      }
      keep_rows1 <- burned_copy$regenera_flag_all == "no_regenera"
    }

    filtered_output1 <- burned_copy[keep_rows1, ] |> dplyr::select(-burned_id)


    # ---- 2. FILTRADO DE POLIGONOS "NO REGENERA" ----

    if (save_no_regenera %in% c("all", "area_filter")) {

      if (!"regenera_flag_all" %in% names(burned_copy)) {
        stop("Missing 'regenera_flag_all'. Make sure regeneration flags were generated.")
      }

      # Select only "no_regenera" polygons
      no_regenera <- burned_copy[burned_copy$regenera_flag_all == "no_regenera", ]

      # Calculate area in ha if not already present
      if (!"burn_area_ha" %in% names(no_regenera)) {
        no_regenera$burn_area_ha <- as.numeric(sf::st_area(no_regenera)) / 10000
      }

      # Apply optional area filter
      if (save_no_regenera == "area_filter") {
        no_regenera <- no_regenera[no_regenera$burn_area_ha >= min_area_no_regenera, ]
      }

     # ---- GUARDADO DE POLIGONOS "NO REGENERA" ----

      # If any rows remain, save them
      if (nrow(no_regenera) > 0) {
        filtered_output1 <- dplyr::select(no_regenera, -burned_id)

        # Determine extension based on format
        ext <- if (output_format == "geojson") ".geojson" else ".shp"


        # Build filename for polygons that did NOT regenerate
        base_name <- tools::file_path_sans_ext(full_out_name)
        filtered_out_name1 <- paste0(base_name, "_noregen", ext)
        filtered_output_file1 <- file.path(output_dir, filtered_out_name1)

        # Remove file if it already exists
        if (file.exists(filtered_output_file1)) {
          file.remove(filtered_output_file1)
        }

        # Write shapefile with non-regenerated polygons
        message("Saving non-regenerated polygons: ", filtered_output_file1)
        message("Number of polygons without regeneration: ", nrow(filtered_output1),
                " (removed: ", nrow(burned_copy) - nrow(filtered_output1), ")")

        suppressWarnings(sf::st_write(filtered_output1, filtered_output_file1, delete_layer = TRUE, quiet = TRUE))
      } else {
        message("No 'no_regenera' polygons met the filter criteria. No file saved.")
      }
    }


  ##### Sustituir poligonos pequenos por P1

    message("Replacing original burned polygons by regenerated polygons in P1: ")

    if (isTRUE(replace_by_P1)) {

      regen_P1_file <- regenera_files[grepl("P1", regenera_files)]
      if (length(regen_P1_file) == 1) {

        regenera_P1 <- sf::st_read(regen_P1_file, quiet = TRUE)

        # Asegurar misma proyeccion
        if (!sf::st_crs(filtered_output) == sf::st_crs(regenera_P1)) {
          regenera_P1 <- sf::st_transform(regenera_P1, sf::st_crs(filtered_output))
        }

        if (isTRUE(validate_geometries)) {
          if (!all(sf::st_is_valid(regenera_P1))) regenera_P1 <- sf::st_make_valid(regenera_P1)
          if (!all(sf::st_is_valid(filtered_output))) filtered_output <- sf::st_make_valid(filtered_output)
        }

        if (!"burn_area_ha" %in% names(filtered_output)) {
          filtered_output$burn_area_ha <- as.numeric(sf::st_area(filtered_output)) / 10000
        }
        filtered_output$fid_final <- seq_len(nrow(filtered_output))

        message("Filtering regenerated polygons in P1: ")
        # Filtrar P1 que realmente se solapan con algun poligono quemado
        regenera_P1 <- sf::st_crop(regenera_P1, sf::st_bbox(filtered_output)) %>%
          dplyr::mutate(
            P1_id = dplyr::row_number(),
            area_ha = as.numeric(sf::st_area(.)) / 10000
          )

        message("Intersecting burned polygons with regenerated ones in P1: ")

        intersec <- suppressWarnings(
          sf::st_intersection(
            filtered_output %>% dplyr::select(fid_final, burn_area_ha),
            regenera_P1 %>% dplyr::select(P1_id, area_ha)
          )
        )

        if (nrow(intersec) > 0) {
          intersec$int_area_ha <- as.numeric(sf::st_area(intersec)) / 10000

          # Asegurar que P1_id y area_ha estan presentes
          if (!all(c("P1_id", "area_ha") %in% names(intersec))) {
            message("Forcing join to recover missing P1_id and area_ha in intersec")
            intersec <- intersec %>%
              dplyr::left_join(
                regenera_P1 %>% sf::st_drop_geometry() %>% dplyr::select(P1_id, area_ha),
                by = "P1_id"
              )
          }

          message("Filtering original burned polygons... ")

          if (!is.null(max_area_ratio_p1)) {
            message("Applying area ratio constraint: only replacing burned polygons where P1 area is greater but not more than ", max_area_ratio_p1, " times the burned area.")
          } else {
            message("No area ratio constraint: replacing all P1 polygons larger than the burned ones.")
          }

          # Extraer pares donde el poligono de regeneracion es mayor
          to_replace <- intersec %>%
            sf::st_set_geometry(NULL) %>%
            dplyr::filter(
              area_ha > burn_area_ha,
              if (!is.null(max_area_ratio_p1)) area_ha / burn_area_ha <= max_area_ratio_p1 else TRUE
            ) %>%
            dplyr::distinct(fid_final, P1_id)


          if (nrow(to_replace) > 0) {

            message("Removing original burned polygons: ")
            # Eliminar poligonos originales
            filtered_output <- filtered_output[!filtered_output$fid_final %in% to_replace$fid_final, ]

            # Insertar nuevos poligonos de P1
            new_geoms <- regenera_P1[to_replace$P1_id, ]
            geom_col <- attr(filtered_output, "sf_column")

            attr_data <- matrix(NA, nrow = nrow(new_geoms), ncol = ncol(filtered_output) - 1) |> as.data.frame()
            names(attr_data) <- names(filtered_output)[names(filtered_output) != geom_col]
            new_rows <- sf::st_sf(attr_data, geometry = sf::st_geometry(new_geoms))

            # Asignar flag de regeneracion si existe
            if ("regenera_flag_P1" %in% names(new_rows)) {
              new_rows$regenera_flag_P1 <- "regenera"
            }

            # Combinar
            filtered_output <- rbind(filtered_output, new_rows)

            # Verificar geometrias vacias o nulas
            check_geom <- sf::st_is_empty(filtered_output) | is.na(sf::st_geometry(filtered_output))

            if (any(check_geom)) {
              message("Se encontraron geometrias vacias o invalidas antes de guardar: ", sum(check_geom))
              filtered_output <- filtered_output[!check_geom, ]
            }

          }
        }

        # LIMPIEZA FINAL DE GEOMETRIAS Y ATRIBUTOS

        message("Cleaning attributes ")
        # 1. Validar y limpiar geometrias vacias
        filtered_output <- sf::st_make_valid(filtered_output)
        filtered_output <- filtered_output[!sf::st_is_empty(filtered_output), ]

        # 2. Extraer solo POLYGONs si vienen de colecciones mixtas
        filtered_output <- sf::st_collection_extract(filtered_output, "POLYGON", warn = FALSE)

        # 3. Separar cada parte del poligono en una fila (esto es lo que elimina el warning)
        filtered_output <- filtered_output |>
          sf::st_cast("POLYGON", group_or_split = TRUE,warn = FALSE)  # Este paso es clave

        # 4. Limpiar atributos

        # Eliminar columnas tipo list
        filtered_output <- filtered_output[, !sapply(filtered_output, is.list), drop = FALSE]

        # Abreviar nombres para Shapefile (max 10 caracteres)
        names(filtered_output) <- make.unique(abbreviate(names(filtered_output), minlength = 10))

        # Convertir factores a caracteres y truncar strings
        filtered_output <- dplyr::mutate(filtered_output, dplyr::across(where(is.factor), as.character))
        filtered_output <- dplyr::mutate(filtered_output, dplyr::across(where(is.character), \(x) substr(x, 1, 254)))

        unique(sf::st_geometry_type(filtered_output))  # Solo deberia ser "POLYGON"


        # --- GUARDAR ---

        # Extract threshold from regeneration file name
        regen_thresh <- stringr::str_extract(basename(regenera_files[1]), "thresh\\d+")

        # Extract base name of the burned file up to before "_metrics" etc.
        burned_base <- basename(burned_file)
        burned_base <- tools::file_path_sans_ext(burned_base)
        short_prefix <- sub("(_metrics.*|_filt.*|_thresh.*)", "", burned_base)

        # Format regeneration ratio string
        ratio_str <- paste0(sprintf("%02d", 100 * min_regen_ratio), collapse = "_")

        # Format area ratio suffix if provided
        if (!is.null(max_area_ratio_p1)) {
          area_ratio_str <- paste0("_ratP1_", formatC(max_area_ratio_p1, format = "f", digits = 0))
        } else {
          area_ratio_str <- ""
        }

        # Final filename including threshold and area ratio
        short_name <- paste0(short_prefix, "_reg_", regen_thresh, "_", regen_labels_str,
                             "_r", ratio_str,
                             area_ratio_str,
                             "_new.shp")

        # Output path
        new_poly_output_file <- file.path(output_dir, short_name)

        delete_shapefile(new_poly_output_file)

        tryCatch({
          suppressWarnings(sf::st_write(filtered_output, new_poly_output_file, delete_layer = TRUE, quiet = TRUE))
          message("Saved new polygons: ", new_poly_output_file)
        }, error = function(e) {
          warning("Failed to write shapefile: ", conditionMessage(e))
        })
      }
    }



  }

}

