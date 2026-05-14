#' Detect and clean a dirty in-memory raster
#'
#' @title Detect and clean a dirty in-memory raster
#' @description
#' Inspect a `SpatRaster` already loaded in memory and, depending on `action`,
#' either clean it, abort, or just report. A raster is considered "dirty" when
#' at least one of three criteria fails:
#' \enumerate{
#'   \item The declared NoData value does not match `expected_nodata`
#'     (an undeclared value also counts as a mismatch).
#'   \item The global minimum is below `lower_cap` (only checked when
#'     `cap_below = TRUE`).
#'   \item There are non-finite pixel values such as `Inf`, `-Inf`, or
#'     pixel-level `NaN` (only checked when `check_finite = TRUE`).
#'     For source-backed rasters whose physical metadata declares the
#'     expected NoData value in all bands, this criterion is skipped:
#'     in-memory non-finite counts in FLT4S/FLT8S files reflect the
#'     declared NoData pixels rather than corruption.
#' }
#' Cleaning, when applied, mirrors the logic used by `mosaic_from_tiles()`
#' when writing a freshly built mosaic: replace non-finite pixels with `NA`,
#' cap pixels below `lower_cap`, and force the NoData flag to
#' `expected_nodata`.
#'
#' This is an internal helper used by `clean_raster_file()` and by the
#' deterministic and supervised orchestrators of OtsuFire as a safeguard
#' before segmentation and feature extraction. For public use, prefer
#' [clean_raster_file()] which operates on file paths.
#' @param r `SpatRaster`. Raster to inspect.
#' @param name Character scalar. Label used in messages, warnings, and stop
#'   conditions to identify the raster being processed.
#' @param expected_nodata Numeric scalar. Value that should be declared as
#'   NoData on the raster.
#' @param lower_cap Numeric scalar. Lower bound; pixels below this value are
#'   considered dirty (and capped during cleaning when `cap_below = TRUE`).
#' @param cap_below Logical scalar. Whether the minimum-below-`lower_cap`
#'   criterion is evaluated (and, during cleaning, applied).
#' @param check_finite Logical scalar. Whether the non-finite-pixels criterion
#'   is evaluated (and, during cleaning, applied).
#' @param action Character scalar. One of `"warn_and_clean"` (default; emit a
#'   warning and return the cleaned raster), `"fail"` (call `stop()` with a
#'   detailed message), or `"report_only"` (emit a warning but return the
#'   raster unchanged, useful for audits).
#' @param verbose Logical scalar. When `TRUE`, emit informational `message()`
#'   calls (including a confirmation when the raster is clean).
#' @return A `SpatRaster`: the cleaned raster when `action = "warn_and_clean"`
#'   and the input was dirty; the original raster when the input is clean or
#'   when `action = "report_only"`. When `action = "fail"` and the input is
#'   dirty the function aborts with `stop()`.
#' @examples
#' \dontrun{
#' r <- terra::rast("path/to/rbr_summer_2022.tif")
#' r_clean <- clean_raster_inmem(r, name = "rbr_summer_2022")
#'
#' clean_raster_inmem(r, name = "rbr_summer_2022", action = "report_only")
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


#' Detect and clean a dirty raster file on disk
#'
#' @title Detect and clean a dirty raster file on disk
#' @description
#' Function for raster diagnostics and cleanup, useful before running segmentation
#' on a mosaic of uncertain origin. It checks for NoData declared correctly in
#' physical metadata, minimum value above `lower_cap`, and (for float rasters)
#' non-finite pixels.
#'
#' When the raster is clean, the function returns the input path unchanged without
#' touching disk regardless of `action`.
#'
#' The on-disk writer uses
#' `wopt = list(NAflag = expected_nodata, gdal = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"))`,
#' matching `mosaic_from_tiles()`.
#' @param raster_path Character scalar. Path to an existing raster file.
#' @param expected_nodata Numeric scalar. Value that should be declared as
#'   NoData on the raster (also used as `NAflag` when rewriting).
#' @param lower_cap Numeric scalar. Lower bound; pixels below this value are
#'   considered dirty and, during cleaning, capped to this value.
#' @param cap_below Logical scalar. Whether the minimum-below-`lower_cap`
#'   criterion is evaluated (and, during cleaning, applied).
#' @param check_finite Logical scalar. Whether the non-finite-pixels criterion
#'   is evaluated (and, during cleaning, applied).
#' @param action Character scalar. One of:
#'   \itemize{
#'     \item `"ask"` (default): when the raster is dirty, print a summary and
#'       prompt the user via `readline()` to choose between
#'       `1 = backup_and_overwrite`, `2 = overwrite` (no backup),
#'       `3 = skip` (do not touch disk), or `4 = fail`.
#'     \item `"backup_and_overwrite"`: rename the original with
#'       `backup_suffix` plus a timestamp inserted before the extension, then
#'       write the cleaned raster to the original path.
#'     \item `"overwrite"`: destructive overwrite without backup. Requires
#'       explicit interactive confirmation: the user must type `OVERWRITE`
#'       (uppercase) at a `readline()` prompt. Use `"backup_and_overwrite"`
#'       in batch/non-interactive contexts.
#'     \item `"fail"`: call `stop()` with a detailed message.
#'   }
#' @param backup_suffix Character scalar. Suffix inserted before the timestamp
#'   when renaming the original file under `"backup_and_overwrite"` (or the
#'   interactive equivalent). Resulting name pattern:
#'   `<stem><backup_suffix>_<YYYYmmdd_HHMMSS><ext>`.
#' @param verbose Logical scalar. When `TRUE`, emit informational `message()`
#'   calls, including a confirmation when the file is already clean.
#' @return A character scalar with the path to the raster file. The path is
#'   unchanged when the file is clean, when the user picks `skip` under
#'   `"ask"`, or when an overwrite is performed (the final file lives at the
#'   original path).
#' @examples
#' # Synthetic example: a clean raster on disk. Because it is already clean,
#' # the function returns the path without prompting even with action = "ask".
#' tmp_path <- tempfile(fileext = ".tif")
#' r <- terra::rast(
#'   nrows = 10, ncols = 10,
#'   xmin = 0, xmax = 10,
#'   ymin = 0, ymax = 10,
#'   vals = c(rep(NA_real_, 20), seq(-500, 290, length.out = 80))
#' )
#' terra::writeRaster(
#'   r,
#'   tmp_path,
#'   overwrite = TRUE,
#'   wopt = list(
#'     datatype = "FLT4S",
#'     NAflag = -9999,
#'     gdal = c("COMPRESS=LZW", "TILED=YES")
#'   )
#' )
#'
#' clean_raster_file(
#'   raster_path = tmp_path,
#'   expected_nodata = -9999,
#'   lower_cap = -1000,
#'   action = "ask",
#'   verbose = TRUE
#' )
#'
#' \dontrun{
#' # Example with a generic user path
#' clean_raster_file(
#'   raster_path = "path/to/your/raster.tif",
#'   expected_nodata = -9999,
#'   lower_cap = -1000,
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
