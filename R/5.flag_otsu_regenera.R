#' Flag Burned Polygons Based on Regeneration Overlap
#'
#' This function assigns regeneration flags to burned area polygons based on their spatial
#' overlap with post-fire regeneration polygons (e.g., from years 1, 2, or more after fire).
#' A polygon is flagged as `"regenera"` if the total intersected area with regeneration polygons
#' from a specific year exceeds a minimum threshold proportion (`min_regen_ratio`) of the
#' burned polygon's area.
#'
#' Filenames of regeneration layers must contain `"P1"`, `"P2"`, etc., to indicate the year
#' after fire. These are automatically extracted and used to name ratio and flag columns.
#'
#' ## Flags and Output Columns:
#' For each regeneration year detected (e.g., P1, P2, P3), the function adds:
#' - `regen_ratio_P1`, `regen_ratio_P2`, ...: ratio of regeneration area to burned area (ha/ha)
#' - `regenera_flag_P1`, `regenera_flag_P2`, ...: `"regenera"` or `"no_regenera"` depending on threshold
#'
#' If `remove_condition = "year1_and_year2"` is selected, an additional field `regenera_flag_joint` is created,
#' set to `"regenera"` only if:
#' - `regen_ratio_P1` > 0 **and**
#' - `regen_ratio_P2` ? threshold.
#'
#' If `remove_condition = "all_years"` is selected, all years in `all_years_vector` must meet the threshold.
#'
#' @param burned_files Character vector with paths to burned area shapefiles (.shp).
#' @param regenera_files Character vector with paths to regeneration shapefiles (.shp).
#'                       Filenames must contain `"P1"`, `"P2"`, etc., to indicate regeneration year after fire.
#' @param min_regen_ratio Numeric vector of thresholds (between 0 and 1) to define sufficient regeneration.
#'                        Default: `c(0.01, 0.05, 0.20)`.
#' @param remove_no_regenera Logical. If `TRUE`, polygons not flagged as `"regenera"` are excluded from the filtered output.
#' @param remove_condition Condition applied when `remove_no_regenera = TRUE`. One of:
#'   \describe{
#'     \item{`"any_year"`}{Keep polygon if regeneration in **at least one year** meets the threshold.}
#'     \item{`"year1_only"`}{Keep only if regeneration in year 1 (`P1`) meets the threshold.}
#'     \item{`"year2_only"`}{Keep only if regeneration in year 2 (`P2`) meets the threshold.}
#'     \item{`"year1_and_year2"`}{Keep only if `regen_ratio_P1 > 0` **and** `regen_ratio_P2 ? threshold`.}
#'     \item{`"all_years"`}{Keep only if **all** years in `all_years_vector` meet the threshold.}
#'   }
#' @param all_years_vector Character vector indicating the labels of regeneration years to consider for `"all_years"` filtering. Example: `c("P1", "P2", "P3")`.
#'                         Only used when `remove_condition = "all_years"`.
#' @param output_dir Optional output directory. If `NULL`, results are saved in the same location as the burned files.
#'
#' @return Three shapefiles per burned area file and threshold:
#' \describe{
#'   \item{Main shapefile}{With regeneration ratios and flags for each polygon.}
#'   \item{Filtered shapefile}{(optional, if `remove_no_regenera = TRUE`) Only polygons flagged as `"regenera"`. Filename includes `"_filter"`.}
#'   \item{Final shapefile with replacements}{Polygons where burned groups were replaced by P1 polygons if larger. Filename includes `"_new_polys"`.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example using multiple thresholds and custom all_years condition
#' burned_files <- list.files("burned_dir", pattern = "burned_.*\\.shp$", full.names = TRUE)
#' regenera_files <- list.files("burned_dir", pattern = "regenera_.*\\.shp$", full.names = TRUE)
#'
#' flag_otsu_regenera(
#'   burned_files = burned_files,
#'   regenera_files = regenera_files,
#'   min_regen_ratio = c(0.10),
#'   remove_no_regenera = TRUE,
#'   remove_condition = "all_years",
#'   all_years_vector = c("P1", "P2", "P3"),
#'   output_dir = "results"
#' )
#' }
#'
#' @importFrom sf st_read st_write st_crs st_transform st_area st_intersection st_set_geometry st_is_valid st_make_valid st_union st_cast st_as_sf
#' @importFrom dplyr mutate row_number group_by summarise select left_join
#' @importFrom tidyr replace_na
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit setNames
#' @importFrom utils glob2rx
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis legend mtext par
#' @importFrom dplyr filter first
#' @importFrom data.table :=
#' @importFrom magrittr %>%
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
                               min_regen_ratio = c(0.01, 0.05, 0.20),
                               remove_no_regenera = FALSE,
                               remove_condition = c("any_year",
                                                    "year1_only",
                                                    "year2_only",
                                                    "year1_and_year2",
                                                    "all_years"),
                               all_years_vector = NULL,
                               output_dir = NULL) {
  remove_condition <- match.arg(remove_condition)
  burned_files <- as.character(burned_files)

  for (burned_file in burned_files) {
    message("Processing: ", burned_file)

    burned_sf <- sf::st_read(burned_file, quiet = TRUE)
    if (!all(sf::st_is_valid(burned_sf)))
      burned_sf <- sf::st_make_valid(burned_sf)

    needs_dissolve <- any(c("CORINE_CLA", "Ens_name") %in% names(burned_sf))
    if (needs_dissolve) {
      burned_sf <- burned_sf |> sf::st_union() |> sf::st_cast("POLYGON") |> sf::st_as_sf() |> dplyr::mutate(burned_id = dplyr::row_number())
    } else {
      burned_sf <- burned_sf |> dplyr::mutate(burned_id = dplyr::row_number())
    }

    burned_sf$burn_area_ha <- as.numeric(sf::st_area(burned_sf)) / 10000
    regen_labels <- c()

    regen_thresh_match <- regmatches(regenera_files, regexpr("thresh[0-9]+", regenera_files))
    regen_thresh <- unique(gsub("thresh", "", regen_thresh_match))
    regen_thresh <- if (length(regen_thresh) == 1)
      regen_thresh
    else
      "NA"

    for (regen_file in regenera_files) {
      regenera_sf <- sf::st_read(regen_file, quiet = TRUE)
      if (!sf::st_crs(burned_sf) == sf::st_crs(regenera_sf)) {
        regenera_sf <- sf::st_transform(regenera_sf, sf::st_crs(burned_sf))
      }

      label <- regmatches(regen_file, regexpr("P[0-9]+", regen_file))
      if (!is.na(label)) {
        regen_labels <- union(regen_labels, label)
        ratio_col <- paste0("regen_ratio_", label)
        intersec <- suppressWarnings(sf::st_intersection(burned_sf |> dplyr::select(burned_id), regenera_sf))
        if (nrow(intersec) > 0) {
          intersec$int_area_ha <- as.numeric(sf::st_area(intersec)) / 10000
          area_summary <- intersec |> sf::st_set_geometry(NULL) |> dplyr::group_by(burned_id) |> dplyr::summarise(regen_area_ha_sum = sum(int_area_ha, na.rm = TRUE))
          burned_sf <- burned_sf |> dplyr::left_join(area_summary, by = "burned_id") |> dplyr::mutate(
            regen_area_ha_sum = tidyr::replace_na(regen_area_ha_sum, 0),!!ratio_col := regen_area_ha_sum / burn_area_ha
          ) |> dplyr::select(-regen_area_ha_sum)
        }
      }
    }

    regen_labels_str <- paste(regen_labels, collapse = "")

    message("Detected regeneration years: ",
            paste(regen_labels, collapse = ", "))


    for (threshold in min_regen_ratio) {
      burned_copy <- burned_sf
      flag_names <- c()

      # Crear flags para cada ano presente
      for (label in regen_labels) {
        ratio_col <- paste0("regen_ratio_", label)
        flag_col <- paste0("regenera_flag_", label)

        # Asegurarse de que ratio_col existe antes de crear la bandera
        if (ratio_col %in% names(burned_copy)) {
          burned_copy[[flag_col]] <- ifelse(burned_copy[[ratio_col]] >= threshold, "regenera", "no_regenera")
          flag_names <- c(flag_names, flag_col)
        } else {
          warning(paste("Missing column:", ratio_col))
        }
      }

      # Anadir flag conjunto si aplica
      if (remove_condition == "year1_and_year2" &&
          all(c("regen_ratio_P1", "regen_ratio_P2") %in% names(burned_copy))) {
        burned_copy$regenera_flag_joint <- ifelse(
          !is.na(burned_copy$regen_ratio_P1) &
            burned_copy$regen_ratio_P1 > 0 &
            !is.na(burned_copy$regen_ratio_P2) &
            burned_copy$regen_ratio_P2 >= threshold,
          "regenera",
          "no_regenera"
        )
        flag_names <- c(flag_names, "regenera_flag_joint")
      }


      # Anadir flag para "all_years"
      if (remove_condition == "all_years" &&
          !is.null(all_years_vector)) {
        all_flags <- paste0("regen_ratio_", all_years_vector)
        missing_flags <- all_flags[!all_flags %in% names(burned_copy)]

        if (length(missing_flags) > 0) {
          message("Missing columns for 'all_years': ",
                  paste(missing_flags, collapse = ", "))
          stop("Cannot evaluate 'all_years' without required ratio columns.")
        }

        # Comprobar si todas las columnas cumplen el threshold
        flag_matrix <- sapply(all_flags, function(col) {
          val <- burned_copy[[col]]
          val[is.na(val)] <- 0  # tratar NAs como 0 para que no pasen el umbral
          val >= threshold
        })

        if (is.vector(flag_matrix)) {
          keep_all <- flag_matrix
        } else {
          keep_all <- apply(flag_matrix, 1, all)
        }

        burned_copy$regenera_flag_allyears <- ifelse(keep_all, "regenera", "no_regenera")
        flag_names <- c(flag_names, "regenera_flag_allyears")

        message("Resumen flags regenera_flag_allyears:")
        print(table(burned_copy$regenera_flag_allyears, useNA = "always"))
      }


      # Definir sufijo
      suffix_condition <- switch(
        remove_condition,
        any_year = "any",
        year1_only = "P1",
        year2_only = "P2",
        year1_and_year2 = "P1P2",
        all_years = "allyears",
        "filter"
      )

      # Guardar shapefile completo
      full_out_name <- paste0(
        tools::file_path_sans_ext(basename(burned_file)),
        "_reg_thr",
        regen_thresh,
        "_",
        regen_labels_str,
        "_",
        sprintf("rat%02d", round(threshold * 100)),
        ".shp"
      )
      full_output_file <- if (!is.null(output_dir))
        file.path(output_dir, full_out_name)
      else
        file.path(dirname(burned_file), full_out_name)
      suppressWarnings(sf::st_write(
        burned_copy,
        full_output_file,
        delete_layer = TRUE,
        quiet = TRUE
      ))

      message("Resumen flags:")
      print(table(burned_copy$regenera_flag_P1, useNA = "ifany"))
      print(table(burned_copy$regenera_flag_P2, useNA = "ifany"))
      print(table(burned_copy$regenera_flag_allyears, useNA = "ifany"))


      # Filtrar segun condicion
      keep_rows <- rep(TRUE, nrow(burned_copy))
      if (isTRUE(remove_no_regenera)) {
        if (remove_condition == "any_year") {
          keep_rows <- burned_copy$regenera_flag_P1 == "regenera" |
            burned_copy$regenera_flag_P2 == "regenera"
        } else if (remove_condition == "year1_only") {
          keep_rows <- burned_copy$regenera_flag_P1 == "regenera"
        } else if (remove_condition == "year2_only") {
          keep_rows <- burned_copy$regenera_flag_P2 == "regenera"
        } else if (remove_condition == "year1_and_year2") {
          keep_rows <- burned_copy$regenera_flag_joint == "regenera"
        } else if (remove_condition == "all_years") {
          keep_rows <- burned_copy$regenera_flag_allyears == "regenera"
        }
      }

      filtered_output <- burned_copy[keep_rows, ] %>% dplyr::select(-burned_id)

      # Guardar shapefile filtrado
      filtered_out_name <- gsub("\\.shp$",
                                paste0("_", suffix_condition, "_filter.shp"),
                                full_out_name)
      filtered_output_file <- if (!is.null(output_dir)) {
        file.path(output_dir, filtered_out_name)
      } else {
        file.path(dirname(burned_file), filtered_out_name)
      }

      message("Guardando shapefile filtrado: ", filtered_output_file)
      message("N? de poligonos tras el filtrado: ", nrow(filtered_output))
      suppressWarnings(
        sf::st_write(
          filtered_output,
          filtered_output_file,
          delete_layer = TRUE,
          quiet = TRUE
        )
      )

      # Sustituir poligonos pequenos por P1
      regen_P1_file <- regenera_files[grepl("P1", regenera_files)]
      if (length(regen_P1_file) == 1) {
        regenera_P1 <- sf::st_read(regen_P1_file, quiet = TRUE)
        if (!sf::st_crs(filtered_output) == sf::st_crs(regenera_P1)) {
          regenera_P1 <- sf::st_transform(regenera_P1, sf::st_crs(filtered_output))
        }

        regenera_P1$area_ha <- as.numeric(sf::st_area(regenera_P1)) / 10000
        filtered_output$burn_area_ha <- as.numeric(sf::st_area(filtered_output)) / 10000
        filtered_output$fid_final <- seq_len(nrow(filtered_output))

        intersec <- suppressWarnings(
          sf::st_intersection(
            regenera_P1 %>% dplyr::mutate(P1_id = dplyr::row_number()),
            filtered_output %>% dplyr::select(fid_final, burn_area_ha)
          )
        )

        if (nrow(intersec) > 0) {
          intersec$int_area_ha <- as.numeric(sf::st_area(intersec)) / 10000
          summed <- intersec %>%
            sf::st_set_geometry(NULL) %>%
            dplyr::group_by(P1_id) %>%
            dplyr::summarise(
              regen_area_ha = first(area_ha),
              total_burned_area = sum(burn_area_ha, na.rm = TRUE),
              burned_fids = list(fid_final)
            )

          to_replace <- summed %>% dplyr::filter(regen_area_ha > total_burned_area)

          if (nrow(to_replace) > 0) {
            message("Replacing groups of burned polygons by larger regeneration P1 polygons...")
            all_fids_to_remove <- unlist(to_replace$burned_fids)
            filtered_output <- filtered_output[!filtered_output$fid_final %in% all_fids_to_remove, ]

            new_geoms <- regenera_P1[to_replace$P1_id, ]
            geom_colname <- attr(filtered_output, "sf_column")
            attr_data <- matrix(NA,
                                nrow = nrow(new_geoms),
                                ncol = ncol(filtered_output) - 1) |> as.data.frame()
            names(attr_data) <- names(filtered_output)[names(filtered_output) != geom_colname]
            new_rows <- sf::st_sf(attr_data, geometry = sf::st_geometry(new_geoms))

            if ("regenera_flag_P1" %in% names(new_rows)) {
              new_rows$regenera_flag_P1 <- "regenera"
            }

            filtered_output <- rbind(filtered_output, new_rows)
          }
        }
      }

      # Eliminar columnas vacias
      non_empty_cols <- sapply(filtered_output, function(col) {
        if (is.list(col))
          return(TRUE)
        if (is.numeric(col))
          return(!(all(is.na(col)) || all(col == 0, na.rm = TRUE)))
        if (is.character(col))
          return(!(all(is.na(col)) || all(col == "", na.rm = TRUE)))
        TRUE
      })
      filtered_output <- filtered_output[, non_empty_cols, drop = FALSE]

      # Guardar shapefile final con sustitucion
      new_poly_name <- gsub("\\.shp$",
                            paste0("_", suffix_condition, "_new_polys.shp"),
                            full_out_name)
      new_poly_output_file <- if (!is.null(output_dir))
        file.path(output_dir, new_poly_name)
      else
        file.path(dirname(burned_file), new_poly_name)
      suppressWarnings(
        sf::st_write(
          filtered_output,
          new_poly_output_file,
          delete_layer = TRUE,
          quiet = TRUE
        )
      )
      message("Saved new polygons: ", new_poly_output_file)
    }

  }
}
