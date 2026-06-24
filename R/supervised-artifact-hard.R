#' Promote artifact_hard hard-negatives into the training features (Phase 2)
#'
#' @description
#' PHASE 2 (artifact_hard hard-negative mining). ADDITIVE + OFF BY DEFAULT.
#'
#' Augments the supervised training-feature frame with "artifact_hard" negatives
#' mined from the deterministic-DROP universe, so they (a) receive spatial fold
#' assignments, (b) appear in the OOF predictions (and so get a held-out
#' `p_oof_mean`), and (c) are trained on by the FINAL model. The promoted rows
#' are ADDED to `train_features`; `scoring_features` is returned UNCHANGED (the
#' rows STAY in the scoring universe too -- they are never removed).
#'
#' This helper is the single wiring point: it runs AFTER feature extraction and
#' BEFORE the spatial folds are made, so the promoted rows are present when folds
#' are assigned.
#'
#' @details
#' When `config$negative_pool_params$artifact_hard$enabled` is `FALSE` (the
#' default) this function is a STRICT NO-OP: it returns `train_features`
#' byte-identical to the input (same rows, same columns, same order) and an
#' empty audit. This is the OFF-by-default contract: with the feature disabled
#' the training population, folds, OOF, model and scoring are all unchanged.
#'
#' The eligibility rule is the exact, user-fixed rule implemented by
#' \code{.of_select_artifact_hard()} (see the Phase 2 spec, section C). Reliable
#' positives (the rbr_med reference for the runtime quantile threshold) are the
#' training rows with \code{class == "burned"}.
#'
#' @param train_features sf / data.frame. The extracted supervised training
#'   features (the \code{train_features} layer). Must carry \code{class} and the
#'   model feature columns.
#' @param scoring_features sf / data.frame. The deterministic scoring universe
#'   (the \code{scoring_features} layer). Used as the source of artifact_hard
#'   candidates. Returned UNCHANGED.
#' @param config An \code{otsufire_supervised_burned_config} (or any list
#'   carrying \code{negative_pool_params$artifact_hard}). Controls whether the
#'   promotion runs and with which thresholds.
#' @param enable_derived Logical (M3 flag). When \code{TRUE}, the three derived
#'   persistence-shape features are attached to the promoted rows. Default
#'   \code{FALSE} (M1: no extra features).
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{train_features} -- the (possibly augmented) training features.
#'       Unchanged when the feature is disabled.
#'     \item \code{scoring_features} -- always the input \code{scoring_features},
#'       unchanged.
#'     \item \code{artifact_hard} -- the promoted rows (0-row frame when none /
#'       disabled).
#'     \item \code{artifact_uncertain} -- the non-promoted deterministic drops
#'       (audit only; NEVER trained).
#'     \item \code{audit} -- the per-clause selection audit (or an empty frame
#'       when disabled).
#'     \item \code{enabled} -- logical; whether promotion ran.
#'     \item \code{n_promoted} -- integer count of promoted rows.
#'   }
#'
#' @family workflow
#' @export
promote_artifact_hard_negatives <- function(train_features,
                                            scoring_features,
                                            config,
                                            enable_derived = FALSE) {
  if (missing(train_features) || is.null(train_features)) {
    stop("'train_features' is required.", call. = FALSE)
  }
  if (missing(config) || is.null(config)) {
    stop("'config' is required.", call. = FALSE)
  }
  # NOTE: deliberately avoid the package `%||%` here -- it coerces with
  # is.na(x[[1L]]) for length-1 lists, which errors when the single element is a
  # multi-value list (e.g. negative_pool_params = list(artifact_hard = <7>)).
  np <- if (is.null(config$negative_pool_params)) list() else config$negative_pool_params
  ah <- if (is.null(np$artifact_hard)) list(enabled = FALSE) else np$artifact_hard

  empty_audit <- data.frame(clause = character(0), n = integer(0),
                            stringsAsFactors = FALSE)

  # OFF-by-default STRICT NO-OP: return train_features byte-identical.
  if (!isTRUE(ah$enabled)) {
    return(list(
      train_features     = train_features,
      scoring_features   = scoring_features,
      artifact_hard      = train_features[integer(0), , drop = FALSE],
      artifact_uncertain = if (!missing(scoring_features) && !is.null(scoring_features)) {
        scoring_features[integer(0), , drop = FALSE]
      } else NULL,
      audit              = empty_audit,
      enabled            = FALSE,
      n_promoted         = 0L
    ))
  }

  if (missing(scoring_features) || is.null(scoring_features)) {
    stop("'scoring_features' is required when artifact_hard is enabled.",
         call. = FALSE)
  }

  # rbr_med reference pool for the runtime quantile floor. Default "negative":
  # the EXISTING random/otsu negative pool's rbr_med (identifies candidates whose
  # RBR is elevated above the background, which is where the haze sits -- above
  # the negatives but below real fire). "positive": the reliable burned rows.
  tdf <- if (inherits(train_features, "sf")) sf::st_drop_geometry(train_features) else as.data.frame(train_features)
  cls <- if ("class" %in% names(tdf)) as.character(tdf[["class"]]) else character(0)
  rbr_med_train <- if ("rbr_med" %in% names(tdf)) suppressWarnings(as.numeric(tdf[["rbr_med"]])) else rep(NA_real_, length(cls))
  ref_kind <- if (is.null(ah$rbr_med_reference)) "negative" else as.character(ah$rbr_med_reference)
  reference_rbr_med <- if (identical(ref_kind, "positive")) {
    rbr_med_train[!is.na(cls) & cls == "burned"]
  } else {
    rbr_med_train[!is.na(cls) & cls == "unburned"]
  }

  sel <- .of_select_artifact_hard(
    scoring_features  = scoring_features,
    reference_rbr_med = reference_rbr_med,
    params            = ah,
    enable_derived    = enable_derived
  )

  augmented <- train_features
  n_promoted <- nrow(sel$artifact_hard)
  if (n_promoted > 0L) {
    # Bind the promoted rows below the existing training rows. Column union is
    # handled by dplyr::bind_rows (missing columns filled with NA), preserving
    # the schema compatibility .of_select_artifact_hard() guarantees.
    if (inherits(train_features, "sf") || inherits(sel$artifact_hard, "sf")) {
      augmented <- dplyr::bind_rows(
        sf::st_as_sf(train_features),
        sf::st_as_sf(sel$artifact_hard)
      )
    } else {
      augmented <- dplyr::bind_rows(train_features, sel$artifact_hard)
    }
  }

  list(
    train_features     = augmented,
    scoring_features   = scoring_features,
    artifact_hard      = sel$artifact_hard,
    artifact_uncertain = sel$artifact_uncertain,
    audit              = sel$audit,
    enabled            = TRUE,
    n_promoted         = as.integer(n_promoted)
  )
}
