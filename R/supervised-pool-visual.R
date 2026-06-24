# =============================================================================
# supervised-pool-visual.R  (Phase 2A)
#
# Consume a VISUALLY-VALIDATED supervised training pool. Reads an existing
# `supervised_training_pool.gpkg` (the consolidated layer written by the
# training-pool-layer stage, which always carries a blank integer `VISUAL`
# column for human review) and applies the visual decisions to produce a CLEAN
# training pool plus full row-level traceability. It NEVER rebuilds the pool:
# a validated pool passed in is used verbatim as a fixed input.
#
# VISUAL semantics (single source of truth, mirrors assemble_era_pool.R):
#   * NA / blank  -> UNREVIEWED. Positives & background are kept as training
#                    rows; artifact_hard candidates are NOT promoted (off until
#                    explicitly confirmed).
#   * 0           -> visually REJECTED.
#                      - positive / background row  -> DROPPED from training.
#                      - artifact_hard candidate    -> a REAL FIRE among the
#                        drops -> excluded from training and routed to the
#                        omission set (for deterministic-phase review).
#   * 1           -> visually CONFIRMED.
#                      - artifact_hard candidate    -> PROMOTED to a hard-negative
#                        (enters training).
#                      - positive / background row  -> confirmed (kept).
#   * anything else (e.g. 2, "x") -> INVALID (see `on_invalid`).
# =============================================================================

#' Apply visual validation to a supervised training pool
#'
#' Reads/annotates a visually-validated supervised training pool and resolves,
#' per row, whether it enters training, following the \code{VISUAL} column
#' semantics described in this file's header. The pool is used as a FIXED input
#' and is never rebuilt.
#'
#' @param pool A visually-validated pool: either a path to a GeoPackage
#'   (the \code{supervised_training_pool} layer is used when present, otherwise
#'   the first layer), or an in-memory \code{sf} / \code{data.frame}.
#' @param visual_col Name of the visual-review column (default \code{"VISUAL"}).
#' @param eligible_col Name of the logical artifact_hard-candidate column
#'   (default \code{"artifact_hard_eligible"}); missing -> all rows treated as
#'   positives/background (no candidates).
#' @param on_invalid What to do with values outside \{NA, 0, 1\}: \code{"error"}
#'   (default) stops and lists the offending values; \code{"exclude"} flags them
#'   as \code{"invalid_excluded"} and keeps them out of training.
#' @param require_visual If \code{TRUE} (default) a missing \code{visual_col}
#'   is an error (the input is not a visually-validated pool).
#'
#' @return A list with class \code{"otsufire_visual_validation"}:
#'   \describe{
#'     \item{\code{pool}}{the full input with three added traceability columns:
#'       \code{visual_value} (parsed integer), \code{visual_action} (factor),
#'       \code{used_for_training} (logical).}
#'     \item{\code{training}}{the CLEAN subset where \code{used_for_training} is
#'       \code{TRUE} (ready for training; geometry preserved when the input was
#'       \code{sf}).}
#'     \item{\code{omissions}}{candidate rows rejected as real fires
#'       (\code{visual_action == "omitted_real_fire"}).}
#'     \item{\code{dropped}}{positive/background rows visually rejected.}
#'     \item{\code{summary}}{a data.frame of counts per \code{visual_action}.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Use a visually-validated pool as a FIXED training input (no rebuild):
#' vv <- apply_visual_validation(
#'   "1985/.../supervised_training_pool.gpkg")
#' nrow(vv$training)          # clean pool ready for training
#' vv$summary                 # audit: how many kept / dropped / promoted / omitted
#' table(vv$pool$used_for_training)
#' }
#' @export
apply_visual_validation <- function(pool,
                                    visual_col   = "VISUAL",
                                    eligible_col = "artifact_hard_eligible",
                                    on_invalid   = c("error", "exclude"),
                                    require_visual = TRUE) {
  on_invalid <- match.arg(on_invalid)

  # ---- read existing pool (NEVER rebuild) ----
  if (is.character(pool) && length(pool) == 1L) {
    if (!file.exists(pool)) stop("apply_visual_validation(): pool file not found: ", pool, call. = FALSE)
    layers <- sf::st_layers(pool)$name
    lyr <- if ("supervised_training_pool" %in% layers) "supervised_training_pool" else layers[1]
    pool <- sf::st_read(pool, layer = lyr, quiet = TRUE)
  }
  if (!is.data.frame(pool)) {
    stop("apply_visual_validation(): 'pool' must be a path, sf, or data.frame.", call. = FALSE)
  }
  d <- if (inherits(pool, "sf")) sf::st_drop_geometry(pool) else as.data.frame(pool)
  n <- nrow(pool)

  # ---- VISUAL column ----
  if (!visual_col %in% names(d)) {
    if (require_visual) {
      stop("apply_visual_validation(): column '", visual_col, "' not found; this is not a ",
           "visually-validated pool (expected the consolidated supervised_training_pool.gpkg).",
           call. = FALSE)
    }
    chr <- rep(NA_character_, n)
  } else {
    chr <- trimws(as.character(d[[visual_col]]))
    chr[chr == ""] <- NA_character_
  }
  v <- suppressWarnings(as.integer(chr))
  invalid <- !is.na(chr) & (is.na(v) | !(v %in% c(0L, 1L)))
  if (any(invalid)) {
    bad <- sort(unique(chr[invalid]))
    if (on_invalid == "error") {
      stop("apply_visual_validation(): invalid ", visual_col, " value(s) {",
           paste(bad, collapse = ", "), "} in ", sum(invalid),
           " row(s). Valid values are NA/blank (unreviewed), 0 (rejected), 1 (confirmed).",
           call. = FALSE)
    }
    v[invalid] <- NA_integer_   # on_invalid == "exclude": neutralise, flagged below
  }

  # ---- candidate vs positive/background ----
  is_cand <- if (eligible_col %in% names(d)) {
    val <- d[[eligible_col]]
    if (is.logical(val)) val %in% TRUE else as.character(val) %in% c("TRUE", "true", "1")
  } else rep(FALSE, n)
  is_cand[is.na(is_cand)] <- FALSE

  # ---- resolve per-row action + training membership ----
  action <- character(n)
  action[ invalid]                                   <- "invalid_excluded"
  und <- !invalid
  action[und & is.na(v) & !is_cand]                  <- "unreviewed_kept"
  action[und & is.na(v) &  is_cand]                  <- "unreviewed_candidate_off"
  action[und & v == 0L  & !is_cand]                  <- "dropped_visual_reject"
  action[und & v == 0L  &  is_cand]                  <- "omitted_real_fire"
  action[und & v == 1L  &  is_cand]                  <- "promoted_artifact_hard"
  action[und & v == 1L  & !is_cand]                  <- "confirmed"

  used <- action %in% c("unreviewed_kept", "promoted_artifact_hard", "confirmed")

  lvls <- c("confirmed", "unreviewed_kept", "promoted_artifact_hard",
            "unreviewed_candidate_off", "dropped_visual_reject",
            "omitted_real_fire", "invalid_excluded")
  pool[["visual_value"]]       <- v
  pool[["visual_action"]]      <- factor(action, levels = lvls)
  pool[["used_for_training"]]  <- used

  summ <- as.data.frame(table(visual_action = pool[["visual_action"]]),
                        responseName = "n", stringsAsFactors = FALSE)

  structure(list(
    pool      = pool,
    training  = pool[used, , drop = FALSE],
    omissions = pool[action == "omitted_real_fire", , drop = FALSE],
    dropped   = pool[action == "dropped_visual_reject", , drop = FALSE],
    summary   = summ
  ), class = "otsufire_visual_validation")
}
