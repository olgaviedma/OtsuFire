#' Promote artifact_hard hard-negatives into the training features
#'
#' @description
#' Mines "artifact_hard" negatives from the deterministic-drop set and adds them
#' to the supervised training features, so the model learns from confusing
#' near-misses. The promoted rows are added to `train_features`; they also stay
#' in the scoring set, so `scoring_features` is returned unchanged.
#'
#' Run it AFTER feature extraction and the spatial folds, and BEFORE the OOF /
#' final-model stages. It selects candidates from `scoring_features` (which only
#' exists post-extraction), so it cannot run before [make_spatial_folds()].
#' Promoted rows receive synthetic per-repeat folds plus a unique
#' `block_id`/`fire_uid` (assigned by the orchestrator after promotion);
#' [make_spatial_folds()] is not re-run. This way the promoted rows still appear
#' in the out-of-fold predictions and are trained on by the final model.
#'
#' This step is additive and off by default. When the `artifact_hard` option is
#' disabled, the function is a strict no-op: it returns `train_features`
#' unchanged (same rows, columns and order) and an empty audit, so the training
#' population, folds, predictions, model and scoring are all unaffected.
#'
#' @details
#' Candidate selection uses a fixed eligibility rule. By default the reference
#' for the runtime RBR quantile floor is the existing negative pool's `rbr_med`
#' (set `rbr_med_reference = "positive"` in the option to use the burned rows
#' instead).
#'
#' @param train_features sf / data.frame. The extracted supervised training
#'   features (the \code{train_features} layer). Must carry \code{class} and the
#'   model feature columns.
#' @param scoring_features sf / data.frame. The deterministic scoring set (the
#'   \code{scoring_features} layer). Used as the source of artifact_hard
#'   candidates. Returned unchanged.
#' @param config An \code{otsufire_supervised_burned_config} (or any list
#'   carrying \code{negative_pool_params$artifact_hard}). Controls whether the
#'   promotion runs and with which thresholds.
#' @param enable_derived Logical. When \code{TRUE}, three derived
#'   persistence-shape features are attached to the promoted rows. Default
#'   \code{FALSE} (no extra features).
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
#'       (audit only; never trained).
#'     \item \code{audit} -- the per-clause selection audit (or an empty frame
#'       when disabled).
#'     \item \code{enabled} -- logical; whether promotion ran.
#'     \item \code{n_promoted} -- integer count of promoted rows.
#'   }
#'
#' @seealso
#' [extract_supervised_features()], [make_spatial_folds()],
#' [build_supervised_training_pools()], [build_supervised_burned_config()]
#'
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L,
#'   negative_pool_params = list(artifact_hard = list(enabled = TRUE))
#' )
#' feats <- extract_supervised_features(folds$train_with_folds, config = cfg)
#' promoted <- promote_artifact_hard_negatives(
#'   feats$train_features, feats$scoring_features, config = cfg
#' )
#' promoted$n_promoted
#' }
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
