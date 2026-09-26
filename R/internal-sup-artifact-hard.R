# =============================================================================
# PHASE 2 (artifact_hard hard-negative mining): selection of "artifact_hard"
# negatives from the deterministic-DROP universe.
#
# ADDITIVE + OFF BY DEFAULT. This file is only ever exercised when the user
# enables `negative_pool_params$artifact_hard$enabled = TRUE`; the exported
# wrapper `promote_artifact_hard_negatives()` is a strict no-op otherwise.
#
# .of_select_artifact_hard() implements the EXACT, user-fixed eligibility rule
# (spec section C) over the scoring universe and partitions the deterministic
# drops into:
#   * artifact_hard       -- promoted hard negatives (stamped as unburned),
#   * artifact_uncertain  -- every OTHER drop (returned for audit, NEVER trained),
#   * audit               -- a one-row-per-clause summary table.
#
# "reliable positives" = the TRAINING rows with class == "burned" (source
# internal_keep_qc); their `rbr_med` distribution supplies the runtime quantile
# threshold (`rbr_med_min_q`). NO EFFIS / external reference is used here
# (contamination overlap is a SEPARATE post-hoc audit, not eligibility).
# =============================================================================

#' Resolve per-row sample weights for the supervised training rows (Phase 2).
#'
#' ADDITIVE + OFF BY DEFAULT. Returns NULL (=> NO weight vector ever set =>
#' byte-identical to today) UNLESS `source_weights` is non-NULL OR an
#' artifact_hard row is present and a target balance is requested.
#'
#' Weight policy (spec E):
#'   * positives (class == "burned")        -> weight 1 (unless source_weights overrides),
#'   * random / otsu negatives              -> weight 1 (unless source_weights overrides),
#'   * artifact_hard negatives              -> per-row weight set so
#'       sum(w[artifact_hard]) == total_weight_ratio * sum(w[positive])
#'       (i.e. total_weight_ratio * n_pos / n_ah),
#'     UNLESS the user pinned an explicit artifact_hard weight in `source_weights`.
#'
#' @param class character vector of class labels (length n).
#' @param source character vector of `source` metadata (length n).
#' @param source_weights NULL or a NAMED numeric (source -> weight) override.
#' @param artifact_hard_source character; source value(s) marking artifact_hard
#'   rows (default "artifact_hard").
#' @param total_weight_ratio numeric scalar: the TOTAL artifact_hard weight as a
#'   fraction of the TOTAL burned (positive) weight (pool-level balance, NOT the
#'   per-polygon weight). Default 0.10. Ignored when an explicit artifact_hard
#'   pin is supplied via `source_weights`.
#' @return list(weights = numeric(n) OR NULL, log = data.frame per-source
#'   n / per_row_weight / total_weight). `weights` is NULL when no weighting
#'   applies (the OFF path).
#' @keywords internal
#' @noRd
.of_resolve_sample_weights <- function(class, source, source_weights = NULL,
                                       artifact_hard_source = "artifact_hard",
                                       total_weight_ratio = 0.10) {
  n <- length(class)
  cls <- as.character(class)
  src <- as.character(source)
  is_pos <- !is.na(cls) & cls == "burned"
  is_ah  <- !is.na(src) & (src %in% artifact_hard_source)

  has_ah <- any(is_ah)
  has_sw <- !is.null(source_weights) && length(source_weights) > 0L

  # OFF path: no artifact_hard rows AND no source_weights override -> NULL.
  if (!has_ah && !has_sw) {
    return(list(weights = NULL, log = .of_sample_weight_log(NULL, src)))
  }

  w <- rep(1.0, n)

  # Source-level overrides take precedence (explicit user pin).
  pinned_ah <- FALSE
  if (has_sw) {
    for (s in names(source_weights)) {
      hit <- !is.na(src) & src == s
      if (any(hit)) w[hit] <- as.numeric(source_weights[[s]])
    }
    pinned_ah <- any(artifact_hard_source %in% names(source_weights))
  }

  # artifact_hard target-balance weight (unless the user explicitly pinned it).
  # sum(w[artifact_hard]) == total_weight_ratio * sum(w[positive]).
  if (has_ah && !pinned_ah) {
    n_ah  <- sum(is_ah)
    sum_w_pos <- sum(w[is_pos])
    per_row <- if (n_ah > 0L) (total_weight_ratio * sum_w_pos) / n_ah else 0
    w[is_ah] <- per_row
  }

  list(weights = w, log = .of_sample_weight_log(w, src))
}

# Build a per-source weight log (n / per-row weight / total weight).
.of_sample_weight_log <- function(w, src) {
  if (is.null(w)) {
    return(data.frame(source = character(0), n = integer(0),
                      per_row_weight = numeric(0), total_weight = numeric(0),
                      stringsAsFactors = FALSE))
  }
  src2 <- ifelse(is.na(src), "<NA>", src)
  us <- sort(unique(src2))
  do.call(rbind, lapply(us, function(s) {
    hit <- src2 == s
    wv <- w[hit]
    data.frame(
      source = s, n = sum(hit),
      per_row_weight = if (length(unique(wv)) == 1L) wv[1L] else NA_real_,
      total_weight = sum(wv), stringsAsFactors = FALSE
    )
  }))
}

# Internal helper: coerce a column to numeric, tolerating comma decimals and
# absent columns (returns all-NA of length n when the column is missing).
.of_ah_num <- function(df, nm, n) {
  if (!nm %in% names(df)) return(rep(NA_real_, n))
  v <- df[[nm]]
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) v <- gsub(",", ".", v, fixed = TRUE)
  suppressWarnings(as.numeric(v))
}

#' Derived artifact_hard features (M3-only; NOT used by M1).
#'
#' The three persistence-shape derived features the spec defines. They are
#' COMPUTED here but only ATTACHED to the promoted rows when `enable_derived =
#' TRUE` (the M3 flag). M1 never sees them, so the M1 feature space / DMatrix is
#' unchanged.
#'
#' @param df data.frame carrying `rbr_med`, `rbr_aw_med`, `rbr_iqr`.
#' @return data.frame with columns persist_loss, persist_abs_delta, rbr_iqr_rel.
#' @keywords internal
#' @noRd
.of_artifact_hard_derived_features <- function(df) {
  n <- nrow(df)
  eps <- 1e-9
  rbr_med    <- .of_ah_num(df, "rbr_med", n)
  rbr_aw_med <- .of_ah_num(df, "rbr_aw_med", n)
  rbr_iqr    <- .of_ah_num(df, "rbr_iqr", n)
  data.frame(
    persist_loss      = (rbr_med - rbr_aw_med) / pmax(abs(rbr_med), eps),
    persist_abs_delta = rbr_med - rbr_aw_med,
    rbr_iqr_rel       = rbr_iqr / pmax(abs(rbr_med), eps),
    stringsAsFactors  = FALSE
  )
}

#' Select artifact_hard / artifact_uncertain rows from the scoring universe.
#'
#' EXACT eligibility (spec C, user-fixed):
#' \preformatted{
#' class_final == "drop" & !is.na(persist_ratio) & !is.na(rbr_med) &
#' rbr_med >= quantile(reference_rbr_med, params$rbr_med_min_q) &
#' persist_ratio <= params$persist_ratio_max &
#' ( reason_1 %in% params$reason_whitelist |
#'   persist_delta <= params$persist_delta_max |
#'   (area_ha >= params$area_ha_min & doy_iqr <= params$doy_iqr_max) )
#' }
#'
#' @param scoring_features sf / data.frame: the Otsu-guided scoring universe.
#' @param reference_rbr_med numeric vector of the REFERENCE pool's rbr_med, used
#'   ONLY to derive the runtime `rbr_med` quantile floor (the
#'   `rbr_med >= quantile(reference_rbr_med, rbr_med_min_q)` gate). The reference
#'   pool is chosen by the caller via `negative_pool_params$artifact_hard$
#'   rbr_med_reference`: `"negative"` (default) -> the existing random/otsu
#'   negative pool's rbr_med (identifies candidates RBR-elevated above the
#'   background); `"positive"` -> the reliable burned positives' rbr_med.
#' @param params the resolved `negative_pool_params$artifact_hard` list.
#' @param enable_derived logical (M3 flag): when TRUE, attach the 3 derived
#'   features to the artifact_hard rows. Default FALSE (M1).
#' @return list(artifact_hard, artifact_uncertain, audit, threshold_rbr_med).
#'   `artifact_hard` rows are stamped class = "unburned", source =
#'   "artifact_hard", neg_type = "artifact_hard_negative"; poly_id / fire_uid and
#'   all feature columns are preserved.
#' @keywords internal
#' @noRd
.of_select_artifact_hard <- function(scoring_features,
                                     reference_rbr_med,
                                     params,
                                     enable_derived = FALSE) {
  df <- if (inherits(scoring_features, "sf")) {
    sf::st_drop_geometry(scoring_features)
  } else {
    as.data.frame(scoring_features)
  }
  n <- nrow(df)

  empty_like <- scoring_features[integer(0), , drop = FALSE]

  # Runtime rbr_med quantile floor from the chosen REFERENCE pool (negatives by
  # default: identifies candidates whose RBR is elevated above the background).
  rp <- suppressWarnings(as.numeric(reference_rbr_med))
  rp <- rp[is.finite(rp)]
  thr_rbr_med <- if (length(rp) > 0L) {
    as.numeric(stats::quantile(rp, probs = params$rbr_med_min_q, names = FALSE,
                               type = 7))
  } else {
    NA_real_
  }

  class_final  <- if ("class_final" %in% names(df)) as.character(df[["class_final"]]) else rep(NA_character_, n)
  persist_ratio <- .of_ah_num(df, "persist_ratio", n)
  rbr_med       <- .of_ah_num(df, "rbr_med", n)
  persist_delta <- .of_ah_num(df, "persist_delta", n)
  area_ha       <- .of_ah_num(df, "area_ha", n)
  doy_iqr       <- .of_ah_num(df, "doy_iqr", n)
  reason_1      <- if ("reason_1" %in% names(df)) as.character(df[["reason_1"]]) else rep(NA_character_, n)

  is_drop <- !is.na(class_final) & class_final == "drop"

  # Clause-by-clause (NA-safe): an NA in any required numeric fails the row.
  c_persist_not_na <- !is.na(persist_ratio)
  c_rbr_not_na     <- !is.na(rbr_med)
  c_rbr_ge_thr     <- !is.na(rbr_med) & is.finite(thr_rbr_med) & (rbr_med >= thr_rbr_med)
  c_persist_le_max <- !is.na(persist_ratio) & (persist_ratio <= params$persist_ratio_max)

  ev_reason  <- !is.na(reason_1) & (reason_1 %in% params$reason_whitelist)
  ev_pdelta  <- !is.na(persist_delta) & (persist_delta <= params$persist_delta_max)
  ev_areadoy <- !is.na(area_ha) & !is.na(doy_iqr) &
    (area_ha >= params$area_ha_min) & (doy_iqr <= params$doy_iqr_max)
  ev_or <- ev_reason | ev_pdelta | ev_areadoy

  eligible <- is_drop & c_persist_not_na & c_rbr_not_na & c_rbr_ge_thr &
    c_persist_le_max & ev_or
  eligible[is.na(eligible)] <- FALSE

  ah_idx  <- which(eligible)
  unc_idx <- which(is_drop & !eligible)

  # Which rule branch produced each eligible row (spec: the consolidated
  # training-pool layer records `artifact_hard_branch`). When more than one
  # evidence clause fires -> "multiple"; otherwise the single firing branch.
  branch_for <- function(i) {
    hits <- c(reason_whitelist = ev_reason[i],
              persist_delta     = ev_pdelta[i],
              large_single_doy  = ev_areadoy[i])
    nm <- names(hits)[which(hits)]
    if (length(nm) > 1L) "multiple" else if (length(nm) == 1L) nm else NA_character_
  }

  artifact_hard <- empty_like
  if (length(ah_idx) > 0L) {
    artifact_hard <- scoring_features[ah_idx, , drop = FALSE]
    # Stamp the label / provenance metadata (schema-compatible with train rows).
    artifact_hard[["class"]]    <- "unburned"
    artifact_hard[["source"]]   <- "artifact_hard"
    artifact_hard[["neg_type"]] <- "artifact_hard_negative"
    # Per-row selection provenance for the consolidated layer (NON-feature
    # columns: dropped by the whitelist filter before the design matrix, so they
    # never reach the model -- they are pure traceability metadata).
    artifact_hard[["artifact_hard_branch"]]        <- vapply(ah_idx, branch_for,
                                                             character(1))
    artifact_hard[["artifact_hard_rbr_threshold"]] <- thr_rbr_med
    if (isTRUE(enable_derived)) {
      derived <- .of_artifact_hard_derived_features(df[ah_idx, , drop = FALSE])
      for (cn in names(derived)) artifact_hard[[cn]] <- derived[[cn]]
    }
  }

  artifact_uncertain <- if (length(unc_idx) > 0L) {
    scoring_features[unc_idx, , drop = FALSE]
  } else {
    empty_like
  }

  audit <- data.frame(
    clause = c("class_final==drop", "persist_ratio not NA", "rbr_med not NA",
               "rbr_med >= q(rbr_med_min_q)", "persist_ratio <= persist_ratio_max",
               "evidence: reason_whitelist", "evidence: persist_delta<=max",
               "evidence: area&doy", "evidence: OR", "ELIGIBLE artifact_hard",
               "artifact_uncertain"),
    n = c(sum(is_drop), sum(is_drop & c_persist_not_na),
          sum(is_drop & c_rbr_not_na), sum(is_drop & c_rbr_ge_thr),
          sum(is_drop & c_persist_le_max), sum(is_drop & ev_reason),
          sum(is_drop & ev_pdelta), sum(is_drop & ev_areadoy),
          sum(is_drop & ev_or), length(ah_idx), length(unc_idx)),
    stringsAsFactors = FALSE
  )

  list(
    artifact_hard      = artifact_hard,
    artifact_uncertain = artifact_uncertain,
    audit              = audit,
    threshold_rbr_med  = thr_rbr_med
  )
}
