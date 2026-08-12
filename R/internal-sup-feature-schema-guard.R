# =============================================================================
# Gate 1E (2026-06-09): RUNTIME feature-schema PARITY GUARD.
#
# Purpose: prevent the OOF (per outer fold), the FINAL refit, and the scoring
# path from silently DIVERGING in feature space on REAL data. This is the
# runtime backstop for the Gate 1D.8 `_isNA` defect (OOF synthesised 51 `_isNA`
# companions while FINAL synthesised none -> the deployed model and the OOF
# diagnostics ran on DIFFERENT feature spaces). Gate 1D.8 unified the three
# paths on ONE shared recipe; this guard ASSERTS, at runtime, that they stay
# unified -- a structural divergence is an ERROR (stop), never a silent warning.
#
# ----------------------------------------------------------------------------
# WHAT THE GUARD COMPARES (the STRICT structural feature-space CONTRACT). These
# MUST be identical across OOF folds, FINAL and scoring:
#   * base_features                (exact names + order)
#   * missing_indicator_features   (exact names + order; the `_isNA` block)
#   * final_feature_order          (exact names + order)
#   * n_base / n_indicators / n_total counts
#   * x_cols / feature_cols        (exact deployed design-matrix column names +
#                                   order, when available)
#   * the TYPE / ENCODING contract (which columns are numeric vs factor + the
#                                   encoding scheme -- NOT the learned level
#                                   VALUES)
#   * the feature-weights POLICY   (the RULE, e.g. "`_isNA` indicators carry the
#                                   canonical weight 1.0" -- NOT a fitted vector)
#   * the schema/recipe CONTRACT VERSION
#       (.OF_SUPERVISED_FEATURE_CONTRACT_VERSION)
#
# WHAT THE GUARD DELIBERATELY DOES *NOT* COMPARE (these legitimately differ
# between an OOF fold and FINAL because they are fit on different training
# sets, and comparing them would produce false-positive aborts):
#   * imputation medians (per-fold / per-final fitted statistic)
#   * learned categorical level VALUES (when variation is legitimate)
#   * scale_pos_weight, best_iteration, seeds, any fold-fitted statistic
#   * wall-clock timestamps (the fingerprint is reproducible / time-independent)
#
# The fingerprint is the SAME deterministic, base-R, wall-clock-FREE style used
# for the fold / pool fingerprints (make_fold_fingerprint /
# otsu_negative_param_fingerprint): a sorted key=value body plus an
# order-stable rolling checksum. No `digest` dependency, no timestamps, no
# data-dependent statistics -> two computations of the SAME contract at
# different wall-clock times produce the IDENTICAL fingerprint.
# =============================================================================

#' Supervised feature-schema / recipe CONTRACT VERSION.
#'
#' Bumped only when the STRUCTURAL feature-space contract changes (the set /
#' order of base features, the `_isNA` indicator scheme, the type/encoding
#' contract, or the feature-weights policy). A change here intentionally
#' invalidates every previously-saved fingerprint, forcing OOF / FINAL / scoring
#' to be recomputed under the new contract. It is part of the hashed payload, so
#' a model trained under a different contract version FAILS the scoring guard.
#'
#' @keywords internal
#' @noRd
.OF_SUPERVISED_FEATURE_CONTRACT_VERSION <- "1.0"

#' The canonical feature-weights POLICY string (the RULE, not a fitted vector).
#'
#' Recorded in (and hashed into) every fingerprint so the guard compares the
#' weight POLICY across OOF / FINAL / scoring -- never the specific learned
#' weight VALUES (which the user may legitimately vary per run via
#' `feature_weights`). The rule is fixed: an `_isNA` indicator carries the
#' canonical weight 1.0 unless the caller names it explicitly; there is no
#' automatic inheritance from the base feature.
#'
#' @keywords internal
#' @noRd
.OF_SUPERVISED_FEATURE_WEIGHTS_POLICY <-
  "isNA_indicators_weight_1_no_base_inheritance"

# Internal: derive the STRUCTURAL feature-space contract payload from a recipe
# (or a column-set list). The payload is the ONLY thing the fingerprint hashes.
# It is intentionally restricted to STRUCTURE: names, order, counts, the
# numeric-vs-factor encoding contract, the weights POLICY string, and the
# contract version. NO medians, NO learned level values, NO seeds, NO
# timestamps. Accepts either:
#   * a fitted recipe from fit_supervised_recipe() / .of_nested_refit_fit()
#     (carries base_features / missing_indicator_features / final_feature_order
#     [+ x_cols / feature_cols when present]); or
#   * the PERSISTED FINAL recipe shape (fields under `$cols`); or
#   * an explicit list(base_features=, missing_indicator_features=,
#     final_feature_order=, x_cols=, feature_cols=).
#
# @keywords internal
# @noRd
.of_feature_schema_payload <- function(recipe_or_cols) {
  x <- recipe_or_cols
  # Unwrap the persisted FINAL recipe (fields live under `$cols`).
  cols <- if (!is.null(x$cols) && is.list(x$cols)) x$cols else x

  pick <- function(nm) {
    v <- cols[[nm]]
    if (is.null(v)) v <- x[[nm]]
    v
  }

  base_features              <- pick("base_features")
  missing_indicator_features <- pick("missing_indicator_features")
  final_feature_order        <- pick("final_feature_order")
  feature_cols               <- pick("feature_cols")
  x_cols                     <- pick("x_cols")

  # Derive any missing partition from final_feature_order (the canonical order).
  if (is.null(final_feature_order)) {
    final_feature_order <- feature_cols %||% x_cols
  }
  if (is.null(missing_indicator_features) && !is.null(final_feature_order)) {
    missing_indicator_features <-
      final_feature_order[grepl("_isNA$", final_feature_order)]
  }
  if (is.null(base_features) && !is.null(final_feature_order)) {
    base_features <- final_feature_order[!grepl("_isNA$", final_feature_order)]
  }

  as_chr <- function(v) if (is.null(v)) character(0) else as.character(v)
  base_features              <- as_chr(base_features)
  missing_indicator_features <- as_chr(missing_indicator_features)
  final_feature_order        <- as_chr(final_feature_order)
  feature_cols               <- as_chr(feature_cols)
  x_cols                     <- as_chr(x_cols)

  # TYPE / ENCODING contract: which columns are factor (one-hot encoded) vs
  # numeric. For the supervised schema since 2026-06-05 there are NO categorical
  # predictors, so every column is numeric -> the contract records the SCHEME and
  # the (possibly empty) set of factor columns by NAME, not their level VALUES.
  factor_cols <- as_chr(pick("factor_cols"))
  if (length(factor_cols) == 0L) {
    fl <- cols[["factor_levels"]]
    if (is.null(fl)) fl <- x[["factor_levels"]]
    if (is.list(fl)) factor_cols <- as_chr(names(fl))
  }
  numeric_cols <- setdiff(final_feature_order, factor_cols)

  list(
    contract_version           = .OF_SUPERVISED_FEATURE_CONTRACT_VERSION,
    base_features              = base_features,
    missing_indicator_features = missing_indicator_features,
    final_feature_order        = final_feature_order,
    feature_cols               = feature_cols,
    x_cols                     = x_cols,
    n_base                     = length(base_features),
    n_indicators               = length(missing_indicator_features),
    n_total                    = length(final_feature_order),
    encoding_scheme            = "sparse_model_matrix_onehot_no_intercept",
    factor_cols                = sort(factor_cols),
    numeric_cols               = numeric_cols,
    weights_policy             = .OF_SUPERVISED_FEATURE_WEIGHTS_POLICY
  )
}

#' Reproducible STRUCTURAL fingerprint of the supervised feature schema.
#'
#' Gate 1E (2026-06-09): the single identity used by the runtime feature-schema
#' parity guard. It hashes ONLY the structural feature-space CONTRACT (see the
#' file header / [.of_feature_schema_payload()]): the base-feature names + order,
#' the `_isNA` indicator names + order, the final feature order, the n_base /
#' n_indicators / n_total counts, the deployed `feature_cols` / `x_cols`, the
#' numeric-vs-factor encoding contract, the feature-weights POLICY string, and
#' the contract version. It DELIBERATELY EXCLUDES imputation medians, learned
#' categorical level values, scale_pos_weight, best_iteration, seeds, any
#' fold-fitted statistic, and any wall-clock timestamp -- so two recipes fit on
#' DIFFERENT training sets but sharing the SAME structural contract (the
#' legitimate OOF-fold vs FINAL case) produce the IDENTICAL fingerprint, while a
#' genuine structural divergence (a reordered / missing / extra column, a
#' dropped `_isNA` indicator) produces a DIFFERENT one.
#'
#' Deterministic, base-R, wall-clock-FREE (mirrors make_fold_fingerprint): a
#' sorted key=value body + an order-stable rolling checksum. No `digest`
#' dependency.
#'
#' @param recipe_or_cols a fitted recipe (fit_supervised_recipe /
#'   .of_nested_refit_fit output), the persisted FINAL recipe (fields under
#'   `$cols`), or an explicit column-set list.
#' @return list(hash = <stable short hash>, payload = <the structured payload it
#'   hashed>, text = <the canonical key=value body>). `hash` is the load-bearing
#'   identity; `payload` lets a manifest show the structure that produced it.
#' @keywords internal
#' @noRd
feature_schema_fingerprint <- function(recipe_or_cols) {
  payload <- .of_feature_schema_payload(recipe_or_cols)

  # Canonical, order-stable key=value body. Vector-valued keys are emitted in
  # their NATURAL order (NOT sorted) because column ORDER is part of the
  # structural contract; scalar keys are sorted by name for stability.
  kv_vec <- function(nm, v) {
    sprintf("%s=[%s]", nm, paste(v, collapse = ","))
  }
  ordered_keys <- c(
    "contract_version",
    "base_features", "missing_indicator_features", "final_feature_order",
    "feature_cols", "x_cols",
    "n_base", "n_indicators", "n_total",
    "encoding_scheme", "factor_cols", "numeric_cols",
    "weights_policy"
  )
  lines <- vapply(ordered_keys, function(nm) {
    v <- payload[[nm]]
    if (length(v) == 1L) sprintf("%s=%s", nm, as.character(v))
    else kv_vec(nm, v)
  }, character(1))

  body <- paste(lines, collapse = "\n")
  bytes <- as.numeric(charToRaw(enc2utf8(body)))
  chk <- 0
  for (b in bytes) chk <- (chk * 31 + b) %% 1000000007
  # A second, independent rolling pass (different multiplier) widens the hash so
  # accidental collisions on the 9-digit space are vanishingly unlikely.
  chk2 <- 0
  for (b in rev(bytes)) chk2 <- (chk2 * 131 + b) %% 1000000007
  hash <- sprintf("fsf1-%09d%09d", as.integer(chk), as.integer(chk2))

  list(hash = hash, payload = payload, text = body)
}

# Internal: assert two STRUCTURAL fingerprints are compatible (identical hash),
# producing a CLEAR, diff-style ERROR (stop) on mismatch. This is the decisive
# assertion shared by OOF (folds vs canonical), FINAL (vs canonical OOF) and
# scoring (produced vs saved). Never a warning: a structural divergence is a
# correctness defect.
#
# @param where character; a short label for the abort message (e.g.
#   "OOF fold rep1/fold2 vs canonical OOF").
# @param expected,produced the two fingerprints (each a feature_schema_fingerprint
#   result, or a bare hash string).
# @keywords internal
# @noRd
.of_assert_schema_fingerprints_equal <- function(where, expected, produced) {
  hx <- if (is.list(expected)) expected$hash else expected
  hp <- if (is.list(produced)) produced$hash else produced
  if (identical(hx, hp)) return(invisible(TRUE))

  # Build a focused diff of the structural payloads to make the abort actionable.
  diff_lines <- character(0)
  if (is.list(expected) && is.list(produced) &&
      !is.null(expected$payload) && !is.null(produced$payload)) {
    pe <- expected$payload; pp <- produced$payload
    keys <- union(names(pe), names(pp))
    for (k in keys) {
      ve <- pe[[k]]; vp <- pp[[k]]
      if (!identical(ve, vp)) {
        fmt <- function(v) {
          if (is.null(v)) return("<absent>")
          if (length(v) > 8L) {
            return(sprintf("[%s, ... (%d items)]",
                           paste(utils::head(v, 8L), collapse = ","), length(v)))
          }
          paste(v, collapse = ",")
        }
        diff_lines <- c(diff_lines,
                        sprintf("  - %s: expected {%s} vs produced {%s}",
                                k, fmt(ve), fmt(vp)))
      }
    }
  }

  stop(sprintf(paste0(
    "Feature-schema parity guard ABORT (%s): the STRUCTURAL feature-space ",
    "contract DIVERGED.\n  expected fingerprint: %s\n  produced fingerprint: %s%s\n",
    "This is the Gate 1D.8 `_isNA` class of defect (OOF / FINAL / scoring must ",
    "share ONE feature space). The guard compares ONLY structure (names, order, ",
    "counts, encoding contract, weights policy, contract version) -- never ",
    "medians / level values / seeds -- so this mismatch is a real divergence."),
    where, hx, hp,
    if (length(diff_lines)) paste0("\n  structural differences:\n",
                                   paste(diff_lines, collapse = "\n")) else ""),
    call. = FALSE)
}

# Internal: establish the CANONICAL OOF fingerprint from a list of per-fold
# fingerprints and ASSERT every fold shares it. A discrepancy is an ERROR (stop).
# Returns the canonical fingerprint (the common one).
#
# @param fold_fps a (named) list of per-fold feature_schema_fingerprint results.
# @keywords internal
# @noRd
.of_oof_canonical_fingerprint <- function(fold_fps) {
  if (length(fold_fps) == 0L) {
    stop(".of_oof_canonical_fingerprint(): no per-fold fingerprints recorded.",
         call. = FALSE)
  }
  canonical <- fold_fps[[1L]]
  nms <- names(fold_fps)
  if (is.null(nms)) nms <- sprintf("fold_%d", seq_along(fold_fps))
  for (i in seq_along(fold_fps)) {
    .of_assert_schema_fingerprints_equal(
      where    = sprintf("OOF %s vs canonical OOF (%s)", nms[[i]], nms[[1L]]),
      expected = canonical,
      produced = fold_fps[[i]]
    )
  }
  canonical
}

# Internal: assemble the run-manifest feature-schema parity record and, for the
# Phase B (nested_refit) path, ABORT the run if the OOF, FINAL and scoring
# structural fingerprints are not compatible. Returns the manifest record
# (a list) regardless; the abort only fires when `phase_b = TRUE` and a
# divergence is found (so the run cannot complete on an incompatible schema).
#
# @param oof_guard the OOF schema_guard (run_oof_xgb / run_dm_oof_pipeline
#   output; NULL on the legacy path / when no OOF ran).
# @param final_fp the FINAL structural fingerprint (feature_schema_fingerprint
#   result, or NULL).
# @param scoring_fp the scoring-matrix structural fingerprint (or NULL).
# @param phase_b logical; TRUE when the run used the nested_refit protocol (the
#   Phase B leakage-free path), which MUST abort on incompatibility.
# @keywords internal
# @noRd
.of_build_schema_parity_manifest <- function(oof_guard = NULL,
                                             final_fp = NULL,
                                             scoring_fp = NULL,
                                             phase_b = FALSE) {
  hash_of <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (is.list(x)) x$hash %||% NA_character_ else as.character(x)
  }
  canonical_oof <- if (!is.null(oof_guard)) oof_guard$canonical else NULL

  per_fold <- NULL
  if (!is.null(oof_guard) && !is.null(oof_guard$per_fold)) {
    per_fold <- vapply(oof_guard$per_fold, hash_of, character(1))
  }

  payload_src <- canonical_oof %||% final_fp %||% scoring_fp
  counts <- if (!is.null(payload_src) && is.list(payload_src) &&
                !is.null(payload_src$payload)) {
    list(n_base = payload_src$payload$n_base,
         n_indicators = payload_src$payload$n_indicators,
         n_total = payload_src$payload$n_total)
  } else list(n_base = NA_integer_, n_indicators = NA_integer_,
              n_total = NA_integer_)

  cv <- if (!is.null(payload_src) && is.list(payload_src) &&
            !is.null(payload_src$payload)) {
    payload_src$payload$contract_version
  } else .OF_SUPERVISED_FEATURE_CONTRACT_VERSION

  # Compatibility of the three legs (only legs that are present participate).
  present <- list(oof = hash_of(canonical_oof),
                  final = hash_of(final_fp),
                  scoring = hash_of(scoring_fp))
  present_hashes <- unlist(present[!vapply(present, is.na, logical(1))])
  compatible <- length(present_hashes) <= 1L ||
    length(unique(present_hashes)) == 1L

  result <- if (compatible) "pass" else "abort"

  manifest <- list(
    contract_version          = cv,
    oof_per_fold_fingerprints = per_fold,
    oof_canonical_fingerprint = hash_of(canonical_oof),
    final_fingerprint         = hash_of(final_fp),
    scoring_fingerprint       = hash_of(scoring_fp),
    n_base                    = counts$n_base,
    n_indicators              = counts$n_indicators,
    n_total                   = counts$n_total,
    phase_b                   = isTRUE(phase_b),
    guard_result              = result
  )

  # Phase B (nested_refit): a structural incompatibility ABORTS the run.
  if (isTRUE(phase_b) && !compatible) {
    stop(sprintf(paste0(
      "Feature-schema parity guard ABORT (Phase B run): the OOF / FINAL / ",
      "scoring STRUCTURAL fingerprints are NOT compatible.\n",
      "  OOF canonical : %s\n  FINAL         : %s\n  scoring       : %s\n",
      "The Phase B (nested_refit) run cannot complete on a divergent feature ",
      "space (the Gate 1D.8 `_isNA` class of defect)."),
      present$oof, present$final, present$scoring), call. = FALSE)
  }

  manifest
}

