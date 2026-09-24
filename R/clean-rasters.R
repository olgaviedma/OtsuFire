#' @title Detect and optionally repair raster integrity problems in memory
#' @description
#' Check an in-memory raster for common integrity problems and, if you
#' choose, repair them before running segmentation, mosaicking,
#' thresholding, or feature-extraction workflows.
#'
#' The function inspects an in-memory `SpatRaster` and evaluates whether it is
#' internally consistent according to three validation rules commonly required
#' in OtsuFire workflows:
#' \itemize{
#'   \item the raster declares the expected NoData value,
#'   \item pixel values do not fall below an accepted lower bound,
#'   \item and the raster does not contain non-finite values such as `Inf`,
#'     `-Inf`, or `NaN`.
#' }
#'
#' Depending on the selected `action`, the function can:
#' \itemize{
#'   \item return the raster unchanged,
#'   \item report integrity problems without modifying the raster,
#'   \item automatically clean the raster,
#'   \item or stop execution with a detailed error.
#' }
#'
#' This function is mainly a defensive preprocessing check before
#' downstream operations such as segmentation, compositing,
#' polygonization, or polygon-level feature extraction.
#'
#' For file-based workflows, see [clean_raster_file()].
#' @param r `SpatRaster`. Raster already loaded in memory.
#' @param name Character scalar. Human-readable label used in messages,
#'   warnings, and error outputs.
#' @param expected_nodata Numeric scalar. Expected NoData value that should be
#'   declared in the raster metadata.
#' @param lower_cap Numeric scalar. Minimum accepted pixel value. Pixels below
#'   this threshold are considered invalid when `cap_below = TRUE`.
#' @param cap_below Logical scalar. Whether values below `lower_cap` should be
#'   treated as integrity problems (and capped during cleaning).
#' @param check_finite Logical scalar. Whether non-finite values (`Inf`,
#'   `-Inf`, `NaN`) should be checked and cleaned.
#' @param action Character scalar controlling what to do when a problem is
#'   found. One of `"warn_and_clean"` (warn and clean automatically),
#'   `"report_only"` (report the problem but leave the raster unchanged),
#'   or `"fail"` (stop with a detailed error).
#' @param verbose Logical scalar. If `TRUE`, informative progress messages are
#'   emitted, including confirmation when the raster passes all checks.
#' @details
#' A raster is considered inconsistent when at least one of the following
#' conditions is detected:
#' \enumerate{
#'   \item the declared NoData value differs from `expected_nodata`,
#'   \item the raster minimum falls below `lower_cap` (when
#'     `cap_below = TRUE`),
#'   \item non-finite pixel values are present (when `check_finite = TRUE`).
#' }
#'
#' For source-backed rasters whose physical metadata already declares the
#' expected NoData value in all bands, the non-finite check is skipped. In
#' floating-point rasters (`FLT4S` or `FLT8S`), these values often correspond
#' to correctly encoded NoData pixels rather than true corruption.
#'
#' When cleaning is applied, the function performs the following operations:
#' \enumerate{
#'   \item replaces `Inf`, `-Inf`, and `NaN` values with `NA`,
#'   \item caps pixel values below `lower_cap`,
#'   \item forces the raster NoData metadata to `expected_nodata`.
#' }
#'
#' The cleaning logic mirrors the behaviour used internally by
#' `mosaic_from_tiles()` when writing newly generated mosaics.
#'
#' In large remote-sensing workflows, corrupted raster metadata or invalid
#' pixel values may silently propagate through segmentation, mosaicking,
#' thresholding, or machine-learning stages.
#'
#' Common consequences include:
#' \itemize{
#'   \item failed polygonization,
#'   \item unstable summary statistics,
#'   \item invalid percentile calculations,
#'   \item crashes during model fitting,
#'   \item inconsistent burned-area delineation,
#'   \item or silent propagation of corrupted values into downstream
#'     products.
#' }
#'
#' In practice, this gives you a reproducible way to standardise raster
#' integrity before analysis.
#' @return A `SpatRaster`.
#'
#'   If the raster is already clean, the original raster is returned
#'   unchanged.
#'
#'   If `action = "warn_and_clean"` and integrity problems are detected, the
#'   cleaned raster is returned.
#'
#'   If `action = "report_only"`, the raster is returned unchanged even when
#'   problems are detected.
#'
#'   If `action = "fail"` and integrity problems are detected, the function
#'   aborts with `stop()`.
#' @examples
#' \dontrun{
#' # Example: diagnose and clean a suspicious mosaic before segmentation
#' r <- terra::rast("MinMin_2022_mosaic_res90m.tif")
#'
#' r_clean <- clean_raster_inmem(
#'   r,
#'   name = "2022 RBR mosaic",
#'   action = "warn_and_clean"
#' )
#'
#' # Example: strict validation inside a production workflow
#' clean_raster_inmem(
#'   r,
#'   name = "2022 RBR mosaic",
#'   action = "fail"
#' )
#'
#' # Example: QA audit without modifying the raster
#' clean_raster_inmem(
#'   r,
#'   name = "2022 RBR mosaic",
#'   action = "report_only"
#' )
#' }
#' @keywords internal
clean_raster_inmem <- function(
    r,
    name = "raster",
    expected_nodata = -9999,
    lower_cap = -1000,
    cap_below = TRUE,
    check_finite = TRUE,
    action = c("warn_and_clean", "fail", "report_only"),
    verbose = TRUE
) {
  action <- match.arg(action)

  if (!inherits(r, "SpatRaster")) {
    stop("'r' must be a SpatRaster.")
  }
  if (!is.character(name) || length(name) != 1L) {
    stop("'name' must be a single character string.")
  }
  if (!is.numeric(expected_nodata) || length(expected_nodata) != 1L) {
    stop("'expected_nodata' must be a single numeric value.")
  }
  if (!is.numeric(lower_cap) || length(lower_cap) != 1L) {
    stop("'lower_cap' must be a single numeric value.")
  }
  if (!is.logical(cap_below) || length(cap_below) != 1L) {
    stop("'cap_below' must be TRUE or FALSE.")
  }
  if (!is.logical(check_finite) || length(check_finite) != 1L) {
    stop("'check_finite' must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE or FALSE.")
  }

  current_nodata <- NA_real_
  current_nodata_label <- "not checked"
  nodata_mismatch <- FALSE

  source_paths <- unique(trimws(terra::sources(r)))
  source_paths <- source_paths[
    !is.na(source_paths) & nzchar(source_paths) & file.exists(source_paths)
  ]

  if (length(source_paths) > 0L) {
    current_nodata <- numeric(0)
    current_nodata_raw <- character(0)

    for (source_path in source_paths) {
      source_info <- terra::describe(source_path)
      band_idx <- grep("^Band [0-9]+ ", source_info)

      if (length(band_idx) == 0L) {
        current_nodata <- c(current_nodata, NA_real_)
        current_nodata_raw <- c(current_nodata_raw, "<missing>")
        next
      }

      band_end <- c(band_idx[-1L] - 1L, length(source_info))
      source_nodata <- rep(NA_real_, length(band_idx))
      source_nodata_raw <- rep("<missing>", length(band_idx))

      for (i in seq_along(band_idx)) {
        band_lines <- source_info[band_idx[i]:band_end[i]]
        nodata_line <- grep("NoData Value=", band_lines, value = TRUE)

        if (length(nodata_line) > 0L) {
          source_nodata_raw[i] <- trimws(sub(".*NoData Value=", "", nodata_line[1L]))
          source_nodata[i] <- suppressWarnings(as.numeric(source_nodata_raw[i]))
        }
      }

      current_nodata <- c(current_nodata, source_nodata)
      current_nodata_raw <- c(current_nodata_raw, source_nodata_raw)
    }

    current_nodata_label <- paste(current_nodata_raw, collapse = ", ")
    nodata_mismatch <- any(is.na(current_nodata) | current_nodata != expected_nodata)
  } else if (isTRUE(verbose)) {
    message(sprintf(
      "[clean_raster_inmem] '%s': NoData mismatch check skipped: in-memory raster has no source file.",
      name
    ))
  }

  metadata_nodata_matches_expected <-
    length(source_paths) > 0L &&
    !any(is.na(current_nodata)) &&
    all(current_nodata == expected_nodata)

  raster_min <- NA_real_
  min_below_cap <- FALSE
  if (isTRUE(cap_below)) {
    raster_min <- terra::global(
      terra::ifel(is.finite(r), r, NA),
      "min",
      na.rm = TRUE
    )[1, 1]
    if (is.finite(raster_min) && raster_min < lower_cap) {
      min_below_cap <- TRUE
    }
  }

  n_nonfinite <- 0
  has_nonfinite <- FALSE
  if (isTRUE(check_finite)) {
    if (isTRUE(metadata_nodata_matches_expected)) {
      if (isTRUE(verbose)) {
        message(sprintf(
          "[clean_raster_inmem] '%s': non-finite check skipped: file declares NoData = %s correctly for all bands; in-memory non-finite count would reflect declared NoData pixels, not corruption.",
          name,
          format(expected_nodata)
        ))
      }
    } else {
      n_nonfinite <- terra::global(!is.finite(r), "sum", na.rm = TRUE)[1, 1]
      if (is.finite(n_nonfinite)) {
        has_nonfinite <- n_nonfinite > 0
      }
    }
  }

  failed <- character(0)
  if (nodata_mismatch) {
    failed <- c(
      failed,
      sprintf(
        "NoData mismatch (declared = %s, expected = %s)",
        current_nodata_label,
        format(expected_nodata)
      )
    )
  }
  if (min_below_cap) {
    failed <- c(
      failed,
      sprintf(
        "minimum value %s is below lower_cap = %s",
        format(raster_min), format(lower_cap)
      )
    )
  }
  if (has_nonfinite) {
    failed <- c(
      failed,
      sprintf("found %s non-finite pixel(s) (Inf/-Inf/NaN)", format(n_nonfinite))
    )
  }

  is_dirty <- length(failed) > 0L

  if (!is_dirty) {
    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_inmem] '%s' is clean.", name))
    }
    return(r)
  }

  diag_msg <- sprintf(
    "Raster '%s' is dirty:\n  - %s",
    name,
    paste(failed, collapse = "\n  - ")
  )

  if (action == "fail") {
    stop(diag_msg, call. = FALSE)
  }

  if (action == "report_only") {
    warning(diag_msg, call. = FALSE)
    return(r)
  }

  # action == "warn_and_clean"
  warning(diag_msg, call. = FALSE)

  if (isTRUE(check_finite) && has_nonfinite) {
    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_inmem] '%s': replacing non-finite pixels with NA.", name))
    }
    r <- terra::ifel(is.finite(r), r, NA)
  }

  if (isTRUE(cap_below) && min_below_cap) {
    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_inmem] '%s': capping values below %s.", name, format(lower_cap)))
    }
    r <- terra::ifel(!is.na(r) & r < lower_cap, lower_cap, r)
  }

  terra::NAflag(r) <- expected_nodata
  if (isTRUE(verbose)) {
    message(sprintf(
      "[clean_raster_inmem] '%s': NAflag set in memory to %s; terra may still report NAflag = NaN when reopening FLT4S files from disk.",
      name,
      format(expected_nodata)
    ))
  }

  r
}


#' @title Detect and optionally repair raster integrity problems on disk
#' @description
#' Check a raster file for common integrity problems and, if needed,
#' repair it directly on disk.
#'
#' It loads a raster with `terra::rast()`, checks whether the file is
#' internally consistent, and optionally rewrites it using the standard
#' OtsuFire output settings.
#'
#' Typical problems detected include:
#' \itemize{
#'   \item missing or inconsistent NoData metadata,
#'   \item corrupted or non-finite pixel values (`Inf`, `-Inf`, `NaN`),
#'   \item unexpected extreme negative values below the accepted lower bound.
#' }
#'
#' The goal is to catch problems early and avoid unstable behaviour in
#' larger automated workflows.
#'
#' When the raster is already clean, the function returns the original file
#' path unchanged without modifying disk contents.
#' @param raster_path Character scalar. Path to an existing raster file.
#' @param expected_nodata Numeric scalar. Expected NoData value to enforce in
#'   the raster metadata.
#' @param lower_cap Numeric scalar. Minimum accepted pixel value. Pixels
#'   below this threshold are considered invalid and are capped during
#'   cleaning.
#' @param cap_below Logical scalar. Whether values below `lower_cap` should
#'   be treated as integrity problems.
#' @param check_finite Logical scalar. Whether non-finite values (`Inf`,
#'   `-Inf`, `NaN`) should be checked and cleaned.
#' @param action Character scalar controlling what to do when integrity
#'   problems are found.
#'
#'   Available options:
#'   \itemize{
#'     \item `"ask"`: interactive diagnostic mode. When the raster is dirty,
#'       prompts the user to choose between backup-and-overwrite, overwrite,
#'       skip, or fail.
#'     \item `"backup_and_overwrite"`: safely rewrites the raster after
#'       creating a timestamped backup.
#'     \item `"overwrite"`: destructive overwrite without backup. In
#'       interactive use, the function requires explicit confirmation before
#'       proceeding.
#'     \item `"fail"`: stops execution with a detailed error.
#'   }
#' @param backup_suffix Character scalar. Suffix inserted before the
#'   timestamp when creating backup files.
#'
#'   Resulting filename pattern:
#'   `<stem><backup_suffix>_<YYYYmmdd_HHMMSS><ext>`.
#' @param verbose Logical scalar. If `TRUE`, informative progress messages
#'   are emitted, including confirmation when the raster is already clean.
#' @details
#' A raster is considered inconsistent when at least one of the following
#' conditions is detected:
#' \enumerate{
#'   \item the declared NoData value differs from `expected_nodata`,
#'   \item the raster minimum falls below `lower_cap`,
#'   \item non-finite values (`Inf`, `-Inf`, `NaN`) are present.
#' }
#'
#' When cleaning is applied, the function performs the following
#' operations:
#' \enumerate{
#'   \item replaces non-finite values with `NA`,
#'   \item caps values below `lower_cap`,
#'   \item rewrites raster metadata using the expected NoData value,
#'   \item writes the raster using standard OtsuFire GDAL options.
#' }
#'
#' The raster is always written using:
#' `wopt = list(NAflag = expected_nodata, gdal = c("COMPRESS=LZW",
#' "BIGTIFF=YES", "TILED=YES"))`.
#'
#' This matches the writing strategy used internally by `mosaic_from_tiles()`.
#'
#' Large remote-sensing workflows frequently depend on chained raster
#' operations involving mosaicking, thresholding, polygonization, and
#' machine-learning pipelines.
#'
#' Corrupted raster metadata or invalid pixel values may silently cause:
#' \itemize{
#'   \item failed processing steps,
#'   \item unstable summary statistics,
#'   \item incorrect thresholds,
#'   \item invalid polygons,
#'   \item model instability,
#'   \item or inconsistent burned-area products.
#' }
#'
#' In practice, this gives you a reproducible way to validate and
#' standardise raster integrity before analysis.
#' @return A character scalar containing the raster path.
#'
#'   If the raster is already clean, the original path is returned
#'   unchanged.
#'
#'   If the raster is cleaned and rewritten, the returned path still points
#'   to the cleaned raster location.
#'
#'   If `action = "ask"` and the user chooses `skip`, the original path is
#'   returned unchanged without rewriting the file.
#'
#'   If `action = "fail"` and integrity problems are detected, the function
#'   aborts with `stop()`.
#' @examples
#' \dontrun{
#' # Example: safely clean a raster while preserving a backup
#' clean_raster_file(
#'   raster_path = "MinMin_2022_mosaic_res90m.tif",
#'   action = "backup_and_overwrite"
#' )
#'
#' # Example: strict validation inside a production pipeline
#' clean_raster_file(
#'   raster_path = "MinMin_2022_mosaic_res90m.tif",
#'   action = "fail"
#' )
#'
#' # Example: interactive QA inspection
#' clean_raster_file(
#'   raster_path = "MinMin_2022_mosaic_res90m.tif",
#'   action = "ask"
#' )
#' }
#' @export
clean_raster_file <- function(
    raster_path,
    expected_nodata = -9999,
    lower_cap = -1000,
    cap_below = TRUE,
    check_finite = TRUE,
    action = c("ask", "backup_and_overwrite", "overwrite", "fail"),
    backup_suffix = "_PRECLEAN",
    verbose = TRUE
) {
  action <- match.arg(action)

  if (!is.character(raster_path) || length(raster_path) != 1L) {
    stop("'raster_path' must be a single character string.")
  }
  if (!file.exists(raster_path)) {
    stop("'raster_path' does not exist: ", raster_path)
  }
  if (!is.character(backup_suffix) || length(backup_suffix) != 1L) {
    stop("'backup_suffix' must be a single character string.")
  }

  r <- terra::rast(raster_path)

  diag_failed <- character(0)
  diag_r <- withCallingHandlers(
    clean_raster_inmem(
      r,
      name = basename(raster_path),
      expected_nodata = expected_nodata,
      lower_cap = lower_cap,
      cap_below = cap_below,
      check_finite = check_finite,
      action = "report_only",
      verbose = FALSE
    ),
    warning = function(w) {
      diag_failed <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  is_dirty <- length(diag_failed) > 0L

  if (!is_dirty) {
    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_file] '%s' is clean. Nothing to do.", basename(raster_path)))
    }
    return(raster_path)
  }

  effective_action <- action

  if (action == "ask") {
    message(diag_failed)
    message(sprintf(
      "[clean_raster_file] '%s' is dirty. Choose an action:",
      basename(raster_path)
    ))
    message("  1 = backup_and_overwrite (rename original, write cleaned over the path)")
    message("  2 = overwrite (no backup; destructive)")
    message("  3 = skip (return path unchanged, leave disk untouched)")
    message("  4 = fail (abort with stop)")
    choice <- readline(prompt = "Selection [1-4]: ")
    choice <- trimws(choice)

    if (identical(choice, "1")) {
      effective_action <- "backup_and_overwrite"
    } else if (identical(choice, "2")) {
      effective_action <- "overwrite_confirmed"
    } else if (identical(choice, "3")) {
      if (isTRUE(verbose)) {
        message(sprintf(
          "[clean_raster_file] Skipped '%s'. Disk untouched.",
          basename(raster_path)
        ))
      }
      return(raster_path)
    } else if (identical(choice, "4")) {
      stop(diag_failed, call. = FALSE)
    } else {
      stop(sprintf(
        "Invalid selection '%s'. Expected 1, 2, 3, or 4. Aborting without touching disk.",
        choice
      ), call. = FALSE)
    }
  }

  if (effective_action == "fail") {
    stop(diag_failed, call. = FALSE)
  }

  if (effective_action == "overwrite") {
    message(diag_failed)
    message(sprintf(
      "[clean_raster_file] About to OVERWRITE '%s' WITHOUT BACKUP. This is destructive.",
      raster_path
    ))
    confirm <- readline(prompt = "Type OVERWRITE (uppercase) to proceed: ")
    if (!identical(confirm, "OVERWRITE")) {
      stop("Overwrite not confirmed. Aborting without touching disk.", call. = FALSE)
    }
    effective_action <- "overwrite_confirmed"
  }

  r_clean <- clean_raster_inmem(
    r,
    name = basename(raster_path),
    expected_nodata = expected_nodata,
    lower_cap = lower_cap,
    cap_below = cap_below,
    check_finite = check_finite,
    action = "warn_and_clean",
    verbose = verbose
  )

  wopt <- list(
    NAflag = expected_nodata,
    gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES")
  )

  if (effective_action == "backup_and_overwrite") {
    ext <- tools::file_ext(raster_path)
    stem <- tools::file_path_sans_ext(basename(raster_path))
    dir_path <- dirname(raster_path)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    backup_name <- sprintf("%s%s_%s.%s", stem, backup_suffix, timestamp, ext)
    backup_path <- file.path(dir_path, backup_name)

    if (file.exists(backup_path)) {
      stop("Backup path already exists, refusing to overwrite: ", backup_path,
           call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message(sprintf(
        "[clean_raster_file] Backing up original to: %s",
        backup_path
      ))
    }
    ok <- file.rename(raster_path, backup_path)
    if (!isTRUE(ok)) {
      stop("Failed to rename original file to backup path: ", backup_path,
           call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_file] Writing cleaned raster to: %s", raster_path))
    }
    terra::writeRaster(r_clean, raster_path, overwrite = TRUE, wopt = wopt)
    return(raster_path)
  }

  if (effective_action == "overwrite_confirmed") {
    if (isTRUE(verbose)) {
      message(sprintf("[clean_raster_file] Overwriting (no backup): %s", raster_path))
    }
    terra::writeRaster(r_clean, raster_path, overwrite = TRUE, wopt = wopt)
    return(raster_path)
  }

  stop("Internal error: unhandled effective action '", effective_action, "'.",
       call. = FALSE)
}
