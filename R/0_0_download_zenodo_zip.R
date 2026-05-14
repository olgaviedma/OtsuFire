#' Download and extract a Zenodo ZIP dataset to a user-specified directory
#'
#' Downloads a ZIP file (e.g., from Zenodo), validates it by minimum size,
#' extracts it to a temporary folder under \code{output_dir}, and then moves the
#' extracted contents into \code{output_dir}.
#'
#' For large ZIP files (multi-GB), this function preferentially uses the system
#' \code{curl} command (if available) to improve robustness and allow resuming
#' interrupted downloads. If extraction fails (often due to a partial/corrupted ZIP),
#' the function can automatically attempt recovery by resuming and/or re-downloading
#' the ZIP and retrying extraction.
#'
#' @details
#' \strong{Directory safety}
#' \itemize{
#'   \item If \code{output_dir} already contains extracted content and \code{overwrite_data=FALSE},
#'   the function stops to avoid overwriting existing files.
#'   \item Set \code{overwrite_data=TRUE} to remove existing extracted content in \code{output_dir}
#'   before extracting again (the ZIP is kept unless \code{clean_zip=TRUE}).
#' }
#'
#' \strong{ZIP structure}
#' If the ZIP contains a single top-level folder (e.g., \code{DATA/}), the contents of that folder
#' are moved directly under \code{output_dir}. Otherwise, all extracted items are moved under
#' \code{output_dir}.
#'
#' \strong{Automatic recovery}
#' When \code{auto_resume_on_unzip_fail=TRUE} (default), if \code{unzip()} fails:
#' \itemize{
#'   \item First, the function tries to resume the download (if \code{curl} is available).
#'   \item If it still fails, the function deletes the ZIP and downloads it again from scratch
#'   (once), then retries extraction.
#' }
#'
#' @param zenodo_url Character. Direct download URL to a ZIP file (e.g., Zenodo \code{?download=1} link).
#' @param output_dir Character. Target directory where the ZIP and extracted dataset will be stored.
#' @param zip_name Character. Local ZIP filename stored inside \code{output_dir}. Default: \code{"DATA.zip"}.
#' @param min_zip_size_mb Numeric. Minimum acceptable ZIP size (MB) to consider the download valid.
#' @param overwrite_zip Logical. If TRUE, deletes any existing ZIP and downloads again from scratch.
#' @param overwrite_data Logical. If TRUE, removes existing extracted content in \code{output_dir} before extraction.
#' @param clean_zip Logical. If TRUE, deletes the ZIP after a successful extraction.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param use_curl Logical. If TRUE and a system \code{curl} executable is available, use it for downloading.
#' Strongly recommended for large ZIP files.
#' @param curl_retries Integer. Number of \code{curl} retry attempts (only used when \code{use_curl=TRUE}).
#' @param auto_resume_on_unzip_fail Logical. If TRUE, attempts automatic recovery when \code{unzip()} fails
#' (resume and/or re-download).
#' @param max_unzip_retries Integer. Maximum number of extraction attempts (including recovery steps).
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{output_dir}: normalized target directory
#'   \item \code{zip_file}: downloaded ZIP path
#'   \item \code{files}: relative paths of extracted files under \code{output_dir}
#' }
#'
#' @examples
#' \dontrun{
#' out <- download_zenodo_data(
#'   zenodo_url = "https://zenodo.org/records/15322380/files/DATA.zip?download=1",
#'   output_dir = "D:/OLGA/ZENODO/OtsuFire_DATA",
#'   use_curl = TRUE
#' )
#' }
#' @importFrom utils download.file unzip head
#' @export

download_zenodo_data <- function(
    zenodo_url,
    output_dir,
    zip_name = "DATA.zip",
    min_zip_size_mb = 1,
    overwrite_zip = FALSE,
    overwrite_data = FALSE,
    clean_zip = FALSE,
    verbose = TRUE,
    use_curl = TRUE,
    curl_retries = 20,
    auto_resume_on_unzip_fail = TRUE,
    max_unzip_retries = 2
) {
  if (!is.character(zenodo_url) || length(zenodo_url) != 1 || !nzchar(zenodo_url)) {
    stop("'zenodo_url' must be a non-empty character string.")
  }
  if (!is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)) {
    stop("'output_dir' must be a non-empty character string.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  zip_file <- file.path(output_dir, zip_name)
  temp_dir <- file.path(output_dir, "temp_extracted")

  # Optionally clear extracted content (keep ZIP)
  if (isTRUE(overwrite_data)) {
    cur <- list.files(output_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    cur <- cur[basename(cur) != zip_name]
    for (p in cur) unlink(p, recursive = TRUE, force = TRUE)
  }

  curl_bin <- Sys.which("curl")
  have_curl <- isTRUE(use_curl) && nzchar(curl_bin)

  download_zip <- function(resume = TRUE, force_fresh = FALSE) {
    if (isTRUE(force_fresh) && file.exists(zip_file)) unlink(zip_file, force = TRUE)

    if (isTRUE(verbose)) {
      message("[ZENODO] download ", if (resume) "(resume)" else "(fresh)", ":")
      message("  ", zenodo_url)
      message("  -> ", zip_file)
    }

    if (have_curl) {
      args <- c(
        "-L",
        "--fail",
        "--retry", as.character(curl_retries),
        "--retry-delay", "5"
      )
      if (isTRUE(resume)) {
        args <- c(args, "-C", "-")  # resume if partial exists
      }
      args <- c(args, "-o", zip_file, zenodo_url)

      status <- system2(curl_bin, args = args)
      if (!identical(status, 0L)) stop("[ZENODO] curl download failed (exit code ", status, ").")
    } else {
      # Fallback: no resume support, but increase timeout
      old_to <- getOption("timeout")
      on.exit(options(timeout = old_to), add = TRUE)
      options(timeout = max(old_to, 3600))
      utils::download.file(zenodo_url, destfile = zip_file, mode = "wb", quiet = !isTRUE(verbose))
    }

    # basic size check
    fs <- file.info(zip_file)$size
    min_bytes <- as.numeric(min_zip_size_mb) * 1024 * 1024
    if (is.na(fs) || fs < min_bytes) {
      stop("[ZENODO] ZIP too small/empty after download: ", zip_file,
           "\nCurrent size (MB): ", round(fs / (1024 * 1024), 2),
           "\nMin required (MB): ", min_zip_size_mb)
    }
    if (isTRUE(verbose)) message("[ZENODO] ZIP size now: ", round(fs / (1024 * 1024), 2), " MB")
    invisible(TRUE)
  }

  # Download if needed (or overwrite_zip)
  if (isTRUE(overwrite_zip) && file.exists(zip_file)) unlink(zip_file, force = TRUE)
  if (!file.exists(zip_file)) {
    download_zip(resume = FALSE, force_fresh = FALSE)
  } else if (isTRUE(verbose)) {
    message("[ZENODO] ZIP already exists:\n  ", zip_file)
    message("[ZENODO] ZIP size: ", round(file.info(zip_file)$size / (1024 * 1024), 2), " MB")
  }

  # Helper: unzip attempt
  try_extract <- function() {
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE, force = TRUE)
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(verbose)) message("[ZENODO] extracting ZIP...")
    ok <- TRUE

    # capture unzip errors
    res <- tryCatch(
      utils::unzip(zip_file, exdir = temp_dir, overwrite = TRUE),
      error = function(e) { ok <<- FALSE; e }
    )

    extracted <- list.files(temp_dir, recursive = TRUE, all.files = FALSE, no.. = TRUE, full.names = TRUE)
    if (!ok || length(extracted) == 0) {
      # cleanup temp
      unlink(temp_dir, recursive = TRUE, force = TRUE)
      return(list(ok = FALSE, msg = if (inherits(res, "error")) conditionMessage(res) else "no_extracted_files"))
    }

    list(ok = TRUE, extracted = extracted)
  }

  # If unzip fails, auto-resume / redownload
  attempt <- 0L
  extracted_info <- NULL

  repeat {
    attempt <- attempt + 1L
    extracted_info <- try_extract()

    if (isTRUE(extracted_info$ok)) break

    if (!isTRUE(auto_resume_on_unzip_fail) || attempt > max_unzip_retries) {
      stop("[ZENODO] unzip failed: ", extracted_info$msg,
           "\nZIP may be incomplete/corrupt: ", zip_file)
    }

    if (isTRUE(verbose)) {
      message("[ZENODO][WARN] unzip failed (attempt ", attempt, "): ", extracted_info$msg)
      message("[ZENODO] attempting automatic recovery...")
    }

    if (attempt == 1L) {
      # First recovery: resume download (keeps partial file)
      if (have_curl) {
        download_zip(resume = TRUE, force_fresh = FALSE)
      } else {
        # no curl -> best effort is full re-download
        download_zip(resume = FALSE, force_fresh = TRUE)
      }
    } else {
      # Second recovery: delete ZIP and download from scratch
      download_zip(resume = FALSE, force_fresh = TRUE)
    }
  }

  # Move extracted content to output_dir
  move_contents <- function(from_dir, to_dir) {
    items <- list.files(from_dir, full.names = TRUE, recursive = FALSE, all.files = TRUE, no.. = TRUE)
    for (p in items) {
      dest <- file.path(to_dir, basename(p))
      if (file.exists(dest) || dir.exists(dest)) {
        if (!isTRUE(overwrite_data)) stop("[ZENODO] target exists: ", dest, " (set overwrite_data=TRUE).")
        unlink(dest, recursive = TRUE, force = TRUE)
      }
      ok <- file.rename(p, dest)
      if (!isTRUE(ok)) {
        if (dir.exists(p)) {
          ok2 <- file.copy(p, dest, recursive = TRUE)
          if (!isTRUE(ok2)) stop("[ZENODO] failed to move directory: ", p)
          unlink(p, recursive = TRUE, force = TRUE)
        } else {
          ok2 <- file.copy(p, dest, overwrite = TRUE)
          if (!isTRUE(ok2)) stop("[ZENODO] failed to move file: ", p)
          unlink(p, force = TRUE)
        }
      }
    }
  }

  # Decide what to move (single root folder / DATA / flat)
  top_dirs  <- list.dirs(temp_dir, full.names = TRUE, recursive = FALSE)
  top_files <- list.files(temp_dir, full.names = TRUE, recursive = FALSE, all.files = TRUE, no.. = TRUE)

  if (length(top_dirs) == 1 && length(top_files) == 0) {
    move_contents(top_dirs[1], output_dir)
  } else if (dir.exists(file.path(temp_dir, "DATA"))) {
    move_contents(file.path(temp_dir, "DATA"), output_dir)
  } else {
    move_contents(temp_dir, output_dir)
  }

  unlink(temp_dir, recursive = TRUE, force = TRUE)
  if (isTRUE(clean_zip) && file.exists(zip_file)) unlink(zip_file, force = TRUE)

  files <- list.files(output_dir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[files != zip_name]
  if (length(files) == 0) stop("[ZENODO] extraction failed: output_dir contains no extracted files.")

  if (isTRUE(verbose)) {
    message("[ZENODO] done. Files extracted: ", length(files))
    message("[ZENODO] first files:")
    print(utils::head(files))
  }

  invisible(list(
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    zip_file = normalizePath(zip_file, winslash = "/", mustWork = FALSE),
    files = files
  ))
}

