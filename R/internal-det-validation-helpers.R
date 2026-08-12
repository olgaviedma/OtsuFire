# ============================================================
# VALIDATION SUMMARY HELPERS
# - Reads the per-scenario validation Excel files
# - Extrae hojas globales de pixel y polygon
# - Builds a yearly summary Excel, per scenario
# ============================================================

find_sheet_by_patterns <- function(xlsx_path, patterns) {
  sheets <- openxlsx::getSheetNames(xlsx_path)
  
  for (pat in patterns) {
    hit <- sheets[grepl(pat, sheets, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  
  NULL
}

read_global_sheet <- function(xlsx_path,
                              scenario,
                              year,
                              type = c("pixel", "polygon")) {
  type <- match.arg(type)
  
  patterns <- switch(
    type,
    pixel   = c("^pixel_global$", "^global_pixel$"),
    polygon = c("^polygon_global$", "^global_polygon$")
  )
  
  sheet_name <- find_sheet_by_patterns(xlsx_path, patterns)
  
  if (is.null(sheet_name)) {
    warning(sprintf(
      "No global sheet of type '%s' was found in: %s",
      type, xlsx_path
    ))
    return(NULL)
  }
  
  df <- openxlsx::read.xlsx(xlsx_path, sheet = sheet_name)
  
  if (is.null(df) || nrow(df) == 0) {
    warning(sprintf(
      "Sheet '%s' is empty in: %s",
      sheet_name, xlsx_path
    ))
    return(NULL)
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$year <- year
  df$scenario <- scenario
  df$source_excel <- basename(xlsx_path)
  df$source_sheet <- sheet_name
  
  first_cols <- c("year", "scenario", "source_excel", "source_sheet")
  other_cols <- setdiff(names(df), first_cols)
  df <- df[, c(first_cols, other_cols), drop = FALSE]
  
  df
}

make_yearly_validation_summary <- function(target_year,
                                           scenarios,
                                           data_base,
                                           result_name = "Min_Min",
                                           output_filename = NULL) {
  
  det_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC"
  )
  
  if (!dir.exists(det_dir)) {
    stop("There is no DETERMINISTIC directory for the year: ", det_dir, call. = FALSE)
  }
  
  if (is.null(output_filename)) {
    output_filename <- paste0("validation_summary_", target_year, ".xlsx")
  }
  
  out_file <- file.path(det_dir, output_filename)
  
  pixel_list <- list()
  polygon_list <- list()
  missing_files <- character(0)
  
  for (sc in scenarios) {
    xlsx_path <- file.path(
      det_dir, sc, "06_VALIDATION",
      paste0("validation_ALL_", target_year, "_", sc, "_res30.xlsx")
    )
    
    if (!file.exists(xlsx_path)) {
      missing_files <- c(missing_files, xlsx_path)
      next
    }
    
    px <- read_global_sheet(
      xlsx_path = xlsx_path,
      scenario  = sc,
      year      = target_year,
      type      = "pixel"
    )
    
    pg <- read_global_sheet(
      xlsx_path = xlsx_path,
      scenario  = sc,
      year      = target_year,
      type      = "polygon"
    )
    
    if (!is.null(px)) pixel_list[[sc]] <- px
    if (!is.null(pg)) polygon_list[[sc]] <- pg
  }
  
  pixel_df <- if (length(pixel_list) > 0) {
    dplyr::bind_rows(pixel_list)
  } else {
    data.frame(
      year = target_year,
      scenario = NA_character_,
      source_excel = NA_character_,
      source_sheet = NA_character_,
      note = "No pixel_global results were found",
      stringsAsFactors = FALSE
    )
  }
  
  polygon_df <- if (length(polygon_list) > 0) {
    dplyr::bind_rows(polygon_list)
  } else {
    data.frame(
      year = target_year,
      scenario = NA_character_,
      source_excel = NA_character_,
      source_sheet = NA_character_,
      note = "No polygon_global results were found",
      stringsAsFactors = FALSE
    )
  }
  
  missing_df <- if (length(missing_files) > 0) {
    data.frame(
      year = target_year,
      missing_excel = missing_files,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      year = target_year,
      missing_excel = "Ninguno",
      stringsAsFactors = FALSE
    )
  }
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "pixel_global")
  openxlsx::writeDataTable(wb, "pixel_global", pixel_df)
  
  openxlsx::addWorksheet(wb, "polygon_global")
  openxlsx::writeDataTable(wb, "polygon_global", polygon_df)
  
  openxlsx::addWorksheet(wb, "missing_or_skipped")
  openxlsx::writeDataTable(wb, "missing_or_skipped", missing_df)
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  
  return(out_file)
}
