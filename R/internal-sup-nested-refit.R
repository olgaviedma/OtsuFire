# =============================================================================
# Shared nested-refit core for the SUPERVISED stage (B1 protocol).
#
# 2026-06-07 (B1 nested_refit): single internal helper shared by BOTH the
# FINAL-model stage (train_final_model_direct) AND each OUTER fold of the OOF
# stage (run_oof_xgb), so the two can never diverge. It fixes the four leakage
# problems found in the prior audit:
#
#   B1a  OOF trained on the FULL negative pool while FINAL trained on the CAPPED
#        pool, so scale_pos_weight differed. -> Caps are now applied UPSTREAM by
#        BOTH callers (FINAL: bucket caps on L_ok; OOF: bucket caps on the outer
#        training negatives) and the SAME helper receives an already-capped
#        training set; spw is computed identically inside.
#   B1b  OOF imputed medians over ALL labelled rows (including the held-out
#        fold). -> Medians here are fit ONLY on the training rows handed in
#        (never on outer_test / held-out data, which the helper never sees).
#   B1c  OOF used the OUTER validation fold for BOTH early stopping AND the
#        reported metric. -> Early stopping uses an INNER validation split carved
#        out of the training set only; the outer test fold is never in any
#        watchlist and is predicted only AFTER the refit.
#   B1d  FINAL picked best_iteration via an internal val split but did NOT refit
#        on all training rows, deploying the tr_idx-only model. -> The helper
#        REFITS a fresh model on ALL of `train_df` at nrounds = best_iteration
#        (no early stopping, no watchlist), and that refit model is what gets
#        deployed / predicted.
#
# The helper NEVER sees or touches any outer_test / held-out data. Imputation
# medians, factor handling, scale_pos_weight, and the inner early-stopping split
# are all fit ONLY on the training rows passed in.
# =============================================================================

# Group split helper, byte-faithful to the FINAL stage's split logic
# (train_final_model_direct lines ~371-396): group split by `block_col` with
# `val_frac`, falling back to a plain row split when a group is missing or the
# resulting split is too small. The caller is responsible for set.seed() before
# calling, so the RNG stream is controlled at the call site (fold_seed).
#
# Returns list(tr_idx, val_idx, mode).
#
# @keywords internal
# @noRd
.of_nested_group_split <- function(n, group_vec, block_col, val_frac) {
  if (!is.null(block_col) && !is.null(group_vec)) {
    g <- as.character(group_vec)
    ug <- unique(g)
    ug <- ug[!is.na(ug)]
    nval_g <- max(1, floor(val_frac * length(ug)))
    val_g <- sample(ug, nval_g)
    val_idx <- which(g %in% val_g)
    tr_idx <- setdiff(seq_len(n), val_idx)
    if (length(val_idx) < 10 || length(tr_idx) < 10) {
      idx <- sample.int(n)
      nval <- max(1, floor(val_frac * n))
      val_idx <- idx[1:nval]
      tr_idx <- idx[(nval + 1):n]
      mode <- "row_split_fallback"
    } else {
      mode <- "group_split"
    }
  } else {
    idx <- sample.int(n)
    nval <- max(1, floor(val_frac * n))
    val_idx <- idx[1:nval]
    tr_idx <- idx[(nval + 1):n]
    mode <- "row_split"
  }
  list(tr_idx = tr_idx, val_idx = val_idx, mode = mode)
}

# Prepare a data.frame of feature columns (factor/character/logical handling)
# WITHOUT imputing numeric NAs. This mirrors the non-numeric coercions performed
# in train_final_model_direct() (lines ~327-339) so the downstream
# sparse.model.matrix() behaves identically, but leaves numeric NAs in place so
# imputation can be fit train-only afterwards.
#
# Returns the coerced data.frame (numeric columns still carry NA where missing,
# and non-finite values coerced to NA).
#
# @keywords internal
# @noRd
.of_nested_coerce_features <- function(X_df, impute_factor_missing) {
  for (nm in names(X_df)) {
    if (is.factor(X_df[[nm]])) X_df[[nm]] <- as.character(X_df[[nm]])
    if (is.character(X_df[[nm]])) {
      X_df[[nm]][is.na(X_df[[nm]])] <- impute_factor_missing
      lv <- unique(X_df[[nm]])
      X_df[[nm]] <- factor(X_df[[nm]], levels = unique(c(lv, "__OTHER__")))
    }
    if (is.factor(X_df[[nm]]) && nlevels(X_df[[nm]]) < 2) {
      lv <- levels(X_df[[nm]])
      levels(X_df[[nm]]) <- unique(c(lv, "__OTHER__"))
    }
    if (is.logical(X_df[[nm]])) X_df[[nm]] <- as.integer(X_df[[nm]])
  }
  X_df
}

# Fit numeric medians on a SUBSET of rows (the training rows) and return them as
# a named list. Non-finite values are treated as NA. DEGENERATE RULE: when a
# numeric feature is all-NA within the fit rows, the median is NA -> we record
# the median as NA (NOT 0, NOT a value fabricated from outside the fit rows) so
# the caller leaves those cells as NA and xgboost treats them as missing. This
# differs from the legacy FINAL path, which substituted 0 for a non-finite
# median; under nested_refit we never fabricate a value from outside the fit
# rows.
#
# @param X_df coerced feature data.frame (output of .of_nested_coerce_features)
# @param fit_rows integer row indices the medians are fit on (train-only)
# @param impute_numeric "median" or "zero"
# @return named list med[[feature]] = numeric scalar OR NA_real_ (degenerate)
# @keywords internal
# @noRd
.of_nested_fit_medians <- function(X_df, fit_rows, impute_numeric) {
  med <- list()
  for (nm in names(X_df)) {
    if (is.numeric(X_df[[nm]]) || is.integer(X_df[[nm]])) {
      v <- X_df[[nm]][fit_rows]
      v[!is.finite(v)] <- NA
      fill <- if (identical(impute_numeric, "zero")) 0 else stats::median(v, na.rm = TRUE)
      # DEGENERATE RULE: all-NA within the fit rows -> keep NA (do not
      # fabricate). xgboost will treat the cell as missing via xgb.DMatrix.
      if (!is.finite(fill)) fill <- NA_real_
      med[[nm]] <- fill
    }
  }
  med
}

# Apply a median recipe to a feature data.frame. For each numeric feature:
#   - coerce non-finite values to NA first,
#   - replace NA with the recipe median, UNLESS the median is NA (degenerate
#     all-NA-in-train column), in which case the cell stays NA so xgboost treats
#     it as missing.
# `_isNA` companions are created per-row BEFORE imputation (row-local), matching
# the FINAL/OOF design-matrix convention but computed on the rows being
# transformed (no cross-row leakage).
#
# @return list(X_df = imputed data.frame). Columns are kept aligned; degenerate
#   columns are NOT dropped (they remain, all-NA, so train/test matrices stay
#   column-aligned).
# @keywords internal
# @noRd
.of_nested_apply_medians <- function(X_df, med) {
  for (nm in names(X_df)) {
    if (is.numeric(X_df[[nm]]) || is.integer(X_df[[nm]])) {
      v <- X_df[[nm]]
      v[!is.finite(v)] <- NA
      fill <- med[[nm]]
      if (!is.null(fill) && is.finite(fill)) {
        v[is.na(v)] <- fill
      }
      # else: degenerate (median NA) -> leave NA -> xgboost missing.
      X_df[[nm]] <- v
    }
  }
  X_df
}

# NA-PRESERVING design-matrix builder. Unlike sparse.model.matrix (which drops
# every row when ANY cell is NA, even with na.action = na.pass), this keeps NA
# cells so they can be handed to xgb.DMatrix as missing values (the
# degenerate-column path). Numeric/integer columns are taken verbatim (NA
# preserved); factor columns are one-hot expanded via model.matrix on the
# complete-case factor (factors never carry numeric NA here because character
# NAs were replaced with the sentinel level upstream). Columns are aligned to
# `ref_cols` when supplied (missing cols added as all-zero, extras dropped,
# reorder).
#
# @return base dense matrix with NA preserved (xgb.DMatrix(missing = NA) treats
#   NA as missing). Returned as a plain matrix so xgboost reads the NAs.
# @keywords internal
# @noRd
.of_nested_build_matrix <- function(X_df, ref_cols = NULL) {
  n <- nrow(X_df)
  parts <- list()
  for (nm in names(X_df)) {
    col <- X_df[[nm]]
    if (is.factor(col)) {
      # One-hot via model.matrix without an intercept; keeps level columns.
      mm <- stats::model.matrix(~ col - 1)
      colnames(mm) <- paste0(nm, sub("^col", "", colnames(mm)))
      # model.matrix may drop rows on NA, but factor cols have no NA here.
      parts[[nm]] <- mm
    } else {
      v <- suppressWarnings(as.numeric(col))
      m <- matrix(v, ncol = 1L, dimnames = list(NULL, nm))
      parts[[nm]] <- m
    }
  }
  M <- if (length(parts) > 0L) do.call(cbind, parts) else matrix(numeric(0), nrow = n, ncol = 0)
  if (is.null(dim(M))) M <- matrix(M, nrow = n)
  if (!is.null(ref_cols)) {
    miss <- setdiff(ref_cols, colnames(M))
    if (length(miss) > 0L) {
      add <- matrix(0, nrow = nrow(M), ncol = length(miss),
                    dimnames = list(NULL, miss))
      M <- cbind(M, add)
    }
    M <- M[, ref_cols, drop = FALSE]
  }
  M
}

#' Shared nested-refit fit for FINAL and each OOF outer fold.
#'
#' Given an ALREADY-CAPPED training data.frame with UN-imputed features, run the
#' leakage-free selection + refit protocol and return the deployable REFIT model
#' plus its recipe (refit medians, feature/x columns) and a full audit record.
#'
#' @param train_df data.frame (geometry dropped) of the capped training rows.
#'   Must carry `label_col` and the model feature columns. May carry
#'   `group_col` / `block_col` (used for the inner split) and bucket-info columns
#'   for the audit.
#' @param feature_cols character vector of feature column names already filtered
#'   to the active whitelist (the model's allowed columns, incl. any `_isNA`).
#' @param label_col character; name of the 0/1-derivable class column (values
#'   compared to "burned").
#' @param group_col character or NULL; unused here directly but recorded.
#' @param block_col character or NULL; grouping column for the inner split.
#' @param val_frac numeric in (0,1); inner validation fraction.
#' @param params_fn function(scale_pos_weight) -> xgboost params list.
#' @param sampling_seed integer; recorded in the audit (the cap RNG seed used
#'   upstream).
#' @param fold_seed integer; seeds BOTH the inner split and both xgb.train calls.
#' @param feature_weights named numeric vector (full per-column weights aligned
#'   to the design matrix) OR NULL.
#' @param nrounds_max integer; max boosting rounds in the selection phase.
#' @param early_stopping_rounds integer; early-stopping patience (selection
#'   phase only).
#' @param impute_numeric "median" (default) or "zero".
#' @param impute_factor_missing sentinel level for character/factor NAs.
#' @param verbose logical.
#'
#' @return list with:
#'   \itemize{
#'     \item `model` — the REFIT xgboost model (deploy this).
#'     \item `best_iteration` — integer from the selection phase.
#'     \item `medians` — the REFIT median recipe (named list; NA = degenerate).
#'     \item `x_cols` — colnames of the refit design matrix.
#'     \item `feature_cols` — the feature columns used.
#'     \item `spw_selection`, `spw_refit` — scale_pos_weight at each phase.
#'     \item `inner_split` — list(tr_idx, val_idx, mode).
#'     \item `audit` — one-row data.frame audit record.
#'   }
#' @keywords internal
#' @noRd
.of_nested_refit_fit <- function(train_df,
                                 feature_cols,
                                 label_col = "class",
                                 group_col = "block_id",
                                 block_col = "block_id",
                                 val_frac = 0.15,
                                 params_fn = .of_canonical_xgb_params,
                                 sampling_seed = 42L,
                                 fold_seed = 42L,
                                 feature_weights = NULL,
                                 nrounds_max = 4000L,
                                 early_stopping_rounds = 80L,
                                 impute_numeric = c("median", "zero"),
                                 impute_factor_missing = "MISSING",
                                 verbose = FALSE) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  stopifnot(requireNamespace("xgboost", quietly = TRUE))
  impute_numeric <- match.arg(impute_numeric)
  if (!is.data.frame(train_df)) stop(".of_nested_refit_fit(): train_df must be a data.frame.")
  if (!label_col %in% names(train_df)) stop(".of_nested_refit_fit(): missing label_col '", label_col, "'.")

  feature_cols <- intersect(feature_cols, names(train_df))
  if (length(feature_cols) == 0L) {
    stop(".of_nested_refit_fit(): no feature columns present in train_df.")
  }

  y_all <- as.integer(train_df[[label_col]] == "burned")
  n_all <- nrow(train_df)
  n_pos <- sum(y_all == 1L)
  n_neg <- sum(y_all == 0L)

  # Coerce non-numeric features once (factor/char/logical) WITHOUT imputing.
  X_df_raw <- .of_nested_coerce_features(train_df[, feature_cols, drop = FALSE],
                                         impute_factor_missing)

  # Build an NA-preserving design matrix (degenerate / missing cells stay NA so
  # xgb.DMatrix(missing = NA) treats them as missing).
  build_mat <- function(X_df_imp, ref_cols = NULL) {
    .of_nested_build_matrix(X_df_imp, ref_cols = ref_cols)
  }

  apply_fw <- function(dmat, x_cols) {
    if (is.null(feature_weights)) return(invisible(NULL))
    fw <- rep(1.0, length(x_cols))
    names(fw) <- x_cols
    in_both <- intersect(names(feature_weights), x_cols)
    if (length(in_both) > 0L) fw[in_both] <- as.numeric(feature_weights[in_both])
    xgboost::setinfo(dmat, "feature_weights", fw)
    invisible(NULL)
  }

  # ---------------------------------------------------------------------------
  # (selection phase)
  # ---------------------------------------------------------------------------
  block_vec <- if (!is.null(block_col) && block_col %in% names(train_df)) {
    train_df[[block_col]]
  } else {
    NULL
  }
  set.seed(fold_seed)
  spl <- .of_nested_group_split(
    n = n_all,
    group_vec = block_vec,
    block_col = if (!is.null(block_vec)) block_col else NULL,
    val_frac = val_frac
  )
  inner_tr <- spl$tr_idx
  inner_val <- spl$val_idx

  # (b) fit medians ONLY on inner_tr; apply to inner_tr AND inner_val.
  med_inner <- .of_nested_fit_medians(X_df_raw, inner_tr, impute_numeric)
  X_inner_imp <- .of_nested_apply_medians(X_df_raw, med_inner)

  # Build inner_tr matrix first to fix the canonical column set, then align
  # inner_val to it (keeps matrices column-aligned even with degenerate cols).
  X_tr_mat  <- build_mat(X_inner_imp[inner_tr, , drop = FALSE], ref_cols = NULL)
  sel_x_cols <- colnames(X_tr_mat)
  X_val_mat <- build_mat(X_inner_imp[inner_val, , drop = FALSE], ref_cols = sel_x_cols)

  y_tr  <- y_all[inner_tr]
  y_val <- y_all[inner_val]

  d_inner_tr  <- xgboost::xgb.DMatrix(X_tr_mat,  label = y_tr,  missing = NA)
  d_inner_val <- xgboost::xgb.DMatrix(X_val_mat, label = y_val, missing = NA)
  apply_fw(d_inner_tr,  sel_x_cols)
  apply_fw(d_inner_val, sel_x_cols)

  # (d) scale_pos_weight on inner_tr only.
  spw_sel <- sum(y_tr == 0L) / max(1, sum(y_tr == 1L))
  params_sel <- params_fn(spw_sel)

  # (e) train with inner_val as the ONLY early-stopping set.
  set.seed(fold_seed)
  m_sel <- xgboost::xgb.train(
    params = params_sel,
    data = d_inner_tr,
    nrounds = nrounds_max,
    watchlist = list(train = d_inner_tr, val = d_inner_val),
    early_stopping_rounds = early_stopping_rounds,
    verbose = if (isTRUE(verbose)) 1 else 0
  )
  best_iteration <- m_sel$best_iteration %||% nrounds_max
  if (!is.finite(best_iteration) || best_iteration < 1L) best_iteration <- 1L
  best_iteration <- as.integer(best_iteration)

  # ---------------------------------------------------------------------------
  # (refit phase)
  # ---------------------------------------------------------------------------
  # (f) refit medians on the FULL train_df; apply to all rows.
  med_refit <- .of_nested_fit_medians(X_df_raw, seq_len(n_all), impute_numeric)
  X_refit_imp <- .of_nested_apply_medians(X_df_raw, med_refit)
  X_refit_mat <- build_mat(X_refit_imp, ref_cols = NULL)
  refit_x_cols <- colnames(X_refit_mat)

  d_refit <- xgboost::xgb.DMatrix(X_refit_mat, label = y_all, missing = NA)
  apply_fw(d_refit, refit_x_cols)

  # (g) scale_pos_weight on the FULL train_df.
  spw_refit <- n_neg / max(1, n_pos)
  params_refit <- params_fn(spw_refit)

  # (h) train a NEW model on ALL rows at nrounds = best_iteration, NO early
  #     stopping, NO watchlist that could leak.
  set.seed(fold_seed)
  m_refit <- xgboost::xgb.train(
    params = params_refit,
    data = d_refit,
    nrounds = best_iteration,
    verbose = if (isTRUE(verbose)) 1 else 0
  )

  # Degenerate columns: features whose REFIT median is NA (all-NA in train_df).
  degenerate_cols <- names(med_refit)[vapply(med_refit, function(z) is.null(z) || is.na(z), logical(1))]

  audit <- data.frame(
    n_train_rows       = n_all,
    n_pos              = n_pos,
    n_neg              = n_neg,
    neg_pos_ratio      = n_neg / max(1, n_pos),
    n_inner_train      = length(inner_tr),
    n_inner_val        = length(inner_val),
    inner_split_mode   = spl$mode,
    spw_selection      = spw_sel,
    spw_refit          = spw_refit,
    best_iteration     = best_iteration,
    fold_seed          = fold_seed,
    sampling_seed      = sampling_seed,
    n_degenerate_cols  = length(degenerate_cols),
    degenerate_cols    = paste(degenerate_cols, collapse = ";"),
    # Structural confirmation: the inner-val set is carved from train_df only,
    # never from any held-out / outer-test data (the helper never receives it).
    val_not_from_heldout = TRUE,
    stringsAsFactors   = FALSE
  )

  list(
    model          = m_refit,
    best_iteration = best_iteration,
    medians        = med_refit,
    x_cols         = refit_x_cols,
    feature_cols   = feature_cols,
    spw_selection  = spw_sel,
    spw_refit      = spw_refit,
    inner_split    = spl,
    degenerate_cols = degenerate_cols,
    audit          = audit
  )
}

# Transform a feature data.frame with a refit recipe (medians + factor handling)
# and return a sparse design matrix aligned to `ref_x_cols`. Used to score the
# outer_test fold with the REFIT medians (no refit on outer_test). Numeric NAs
# left where the median is degenerate (NA) so xgboost treats them as missing.
#
# @keywords internal
# @noRd
.of_nested_transform <- function(df, feature_cols, med, ref_x_cols,
                                 impute_factor_missing = "MISSING") {
  feature_cols <- intersect(feature_cols, names(df))
  X_df_raw <- .of_nested_coerce_features(df[, feature_cols, drop = FALSE],
                                         impute_factor_missing)
  X_imp <- .of_nested_apply_medians(X_df_raw, med)
  # NA-preserving build aligned to the refit column set (degenerate cells stay
  # NA -> xgb.DMatrix(missing = NA) treats them as missing, exactly as the
  # refit model was trained).
  .of_nested_build_matrix(X_imp, ref_cols = ref_x_cols)
}
