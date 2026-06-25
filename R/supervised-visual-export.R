# =============================================================================
# supervised-visual-export.R  (AUDIT PENDING A)
#
# Public exporter for VISUAL review. Turns a (visually-validated) supervised
# training pool into review-ready layers, split by bucket and with the VISUAL
# column kept editable, so the user can review in QGIS and re-feed the edited
# pool back into apply_visual_validation(). It NEVER trains, scores or rebuilds
# anything: it only reads (via apply_visual_validation), labels and writes.
#
# Output GPKG layers:
#   all_tagged          - the FULL pool, the canonical re-entry layer.
#   keep                - high_confidence_keep (positives).
#   otsu_residual       - otsu drop-only negatives.
#   random_background   - random synthetic background.
#   artifact_hard       - artifact_hard candidates (kept as ONE family; the
#                         per-row decision lives in the visual_action column).
# =============================================================================

# Map pool_source / source (+ artifact_hard_eligible) -> review_bucket.
.of_review_bucket <- function(d) {
  n <- nrow(d)
  elig <- rep(FALSE, n)
  if ("artifact_hard_eligible" %in% names(d)) {
    v <- d[["artifact_hard_eligible"]]
    elig <- if (is.logical(v)) v %in% TRUE else as.character(v) %in% c("TRUE", "true", "1")
  }
  elig[is.na(elig)] <- FALSE
  ps <- if ("pool_source" %in% names(d)) as.character(d[["pool_source"]]) else rep(NA_character_, n)
  sc <- if ("source"      %in% names(d)) as.character(d[["source"]])      else rep(NA_character_, n)
  bucket <- rep(NA_character_, n)
  bucket[ps %in% "high_confidence_keep" | sc %in% "internal_keep_qc"]        <- "keep"
  bucket[ps %in% "otsu"                 | sc %in% "otsu_patch_residual"]      <- "otsu_residual"
  bucket[ps %in% "random"               | sc %in% "random_burnable_background"] <- "random_background"
  bucket[ps %in% "deterministic_drop"   | sc %in% "artifact_hard"]           <- "artifact_hard"
  bucket[elig] <- "artifact_hard"             # eligibility wins (the candidate family)
  bucket[is.na(bucket)] <- "other"
  bucket
}

# Reorder columns so the requested front fields (those present) come first,
# geometry last; sf object in, sf object out.
.of_reorder_front <- function(x, front) {
  gcol <- attr(x, "sf_column")
  nm   <- setdiff(names(x), gcol)
  f    <- front[front %in% nm]
  rest <- setdiff(nm, f)
  x[, c(f, rest, gcol)]
}

#' Export a visually-validated training pool into review-ready layers
#'
#' @description
#' A REVIEW EXPORTER (not a model step). It takes a supervised training pool,
#' resolves the VISUAL decisions through \code{\link{apply_visual_validation}}
#' and writes/returns the pool split by \code{review_bucket} so it is easy to
#' inspect in QGIS. It NEVER trains a model, NEVER scores and NEVER rebuilds the
#' pool; the input object is not modified.
#'
#' @details
#' Typical loop: (1) export with this function, (2) open the GPKG in QGIS and
#' edit the \code{VISUAL} column (\code{0}/\code{1}/blank) on the rows you
#' review, (3) read the edited \code{all_tagged} layer back and pass it to
#' \code{\link{apply_visual_validation}} to get the clean training set. The
#' \code{all_tagged} layer is the canonical RE-ENTRY layer: it carries
#' \code{VISUAL} and \code{artifact_hard_eligible}, so editing \code{VISUAL}
#' there and re-reading preserves the full semantics. The per-bucket layers
#' carry the same columns and the same semantics; they are a convenience for
#' visual review only.
#'
#' \strong{VISUAL semantics} (same as \code{\link{apply_visual_validation}}):
#' \code{VISUAL = 1} confirms; \code{VISUAL = 0} rejects; blank/\code{NA} is
#' unreviewed. For the \code{artifact_hard} family the meaning is special:
#' \itemize{
#'   \item \code{VISUAL = 1}: promote the candidate as a HARD-NEGATIVE (it is a
#'     confirmed false positive, usable as a negative).
#'   \item \code{VISUAL = 0}: treat the candidate as an OMITTED REAL FIRE -- do
#'     NOT use it as a negative; it goes to the omissions set.
#'   \item \code{VISUAL = NA}: not reviewed / not promoted.
#' }
#'
#' @param pool A supervised training pool: a path to a GeoPackage (the
#'   \code{supervised_training_pool} layer is used when present, otherwise the
#'   first layer), an in-memory \code{sf}, or an
#'   \code{otsufire_visual_validation} object (the result of
#'   \code{\link{apply_visual_validation}}).
#' @param out_path Output GeoPackage path. If \code{NULL} (default) nothing is
#'   written and only the in-memory layers + summary are returned.
#' @param overwrite Logical (default \code{FALSE}). If the output exists and
#'   \code{overwrite = FALSE}, the function stops; if \code{TRUE} the file (and
#'   any \code{-wal}/\code{-shm} sidecars) is removed first.
#' @param visual_col,eligible_col,on_invalid Passed to
#'   \code{\link{apply_visual_validation}} when \code{pool} is not already a
#'   validation result (defaults \code{"VISUAL"}, \code{"artifact_hard_eligible"},
#'   \code{"error"}).
#' @param verbose Logical (default \code{TRUE}). Emit progress messages.
#'
#' @return Invisibly, a list with:
#'   \describe{
#'     \item{\code{layers}}{a named list of \code{sf} layers:
#'       \code{all_tagged}, \code{keep}, \code{otsu_residual},
#'       \code{random_background}, \code{artifact_hard}.}
#'     \item{\code{summary}}{a list of count tables: \code{by_bucket},
#'       \code{by_visual}, \code{by_action}, \code{by_training} and
#'       \code{bucket_x_action}.}
#'     \item{\code{path}}{the output path (or \code{NULL}).}
#'     \item{\code{n}}{total row count of \code{all_tagged}.}
#'   }
#'
#' @examples
#' \dontrun{
#' # 1) Export for review:
#' ex <- export_visual_validation_pools(
#'   "1985/.../supervised_training_pool.gpkg",
#'   out_path = "C:/of_tmp/1985_pool_visual.gpkg", overwrite = TRUE)
#' ex$summary$by_bucket
#'
#' # 2) ...edit VISUAL in QGIS on the all_tagged layer...
#'
#' # 3) Re-feed the edited pool:
#' edited <- sf::st_read("C:/of_tmp/1985_pool_visual.gpkg", layer = "all_tagged")
#' vv <- apply_visual_validation(edited)
#' vv$summary
#' }
#' @seealso \code{\link{apply_visual_validation}}
#' @export
export_visual_validation_pools <- function(pool,
                                           out_path     = NULL,
                                           overwrite    = FALSE,
                                           visual_col   = "VISUAL",
                                           eligible_col = "artifact_hard_eligible",
                                           on_invalid   = c("error", "exclude"),
                                           verbose      = TRUE) {
  on_invalid <- match.arg(on_invalid)

  # ---- resolve the annotated pool through apply_visual_validation (no rebuild) ----
  vv <- if (inherits(pool, "otsufire_visual_validation")) {
    pool
  } else {
    apply_visual_validation(pool, visual_col = visual_col,
                            eligible_col = eligible_col, on_invalid = on_invalid)
  }
  tagged <- vv$pool
  if (!inherits(tagged, "sf")) {
    stop("export_visual_validation_pools(): the pool must carry geometry (sf) to export.", call. = FALSE)
  }

  # ---- review_bucket + column ordering (input object never modified) ----
  d <- sf::st_drop_geometry(tagged)
  tagged[["review_bucket"]] <- .of_review_bucket(d)
  front <- c("VISUAL", "review_bucket", "visual_action", "used_for_training",
             "pool_source", "source", "class", "artifact_hard_eligible",
             "fire_uid", "block_id", "poly_id", "source_poly_id", "year")
  tagged <- .of_reorder_front(tagged, front)
  rb <- tagged[["review_bucket"]]

  # ---- split into the per-bucket layers ----
  buckets <- c("keep", "otsu_residual", "random_background", "artifact_hard")
  layers <- list(all_tagged = tagged)
  for (b in buckets) layers[[b]] <- tagged[rb == b, , drop = FALSE]

  # ---- count summaries ----
  td <- sf::st_drop_geometry(tagged)
  tab <- function(x) as.data.frame(table(x, useNA = "ifany"), responseName = "n",
                                   stringsAsFactors = FALSE)
  summary <- list(
    by_bucket       = stats::setNames(tab(rb),                  c("review_bucket",  "n")),
    by_visual       = stats::setNames(tab(td[[visual_col]]),    c("VISUAL",         "n")),
    by_action       = stats::setNames(tab(td[["visual_action"]]),     c("visual_action", "n")),
    by_training     = stats::setNames(tab(td[["used_for_training"]]), c("used_for_training", "n")),
    bucket_x_action = as.data.frame(table(review_bucket = rb,
                                          visual_action = td[["visual_action"]]),
                                    responseName = "n", stringsAsFactors = FALSE)
  )

  # ---- write the GPKG (all_tagged first, then each bucket as a new layer) ----
  if (!is.null(out_path)) {
    if (file.exists(out_path) && !isTRUE(overwrite)) {
      stop("export_visual_validation_pools(): output exists and overwrite=FALSE: ",
           out_path, call. = FALSE)
    }
    for (s in c("", "-wal", "-shm")) if (file.exists(paste0(out_path, s))) unlink(paste0(out_path, s), force = TRUE)
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    write_order <- c("all_tagged", buckets)
    for (lyr in write_order) {
      x <- layers[[lyr]]
      ok <- tryCatch({ sf::st_write(x, out_path, layer = lyr, quiet = TRUE); TRUE },
                     error = function(e) { warning("layer '", lyr, "' not written: ",
                                                   conditionMessage(e), call. = FALSE); FALSE })
      if (isTRUE(verbose)) message(sprintf("  [export] layer %-18s %5d rows  %s",
                                           lyr, nrow(x), if (ok) "written" else "SKIPPED"))
    }
    if (isTRUE(verbose)) message("[export_visual_validation_pools] wrote ", out_path)
  }

  invisible(list(layers = layers, summary = summary, path = out_path, n = nrow(tagged)))
}
