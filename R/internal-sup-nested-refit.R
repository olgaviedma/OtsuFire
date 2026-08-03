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
# xgboost 2.0.0 renamed the early-stopping evaluation set from `watchlist` to
# `evals`; the old name still works but warns, and upstream states it will
# become an error. Passing it under the wrong name leaves the evaluation set
# unregistered, after which the early-stopping callback compares an NA score
# and fails. Build the call dynamically so the package behaves identically on
# both API generations.
.of_xgb_train_evals <- function(params, data, nrounds, evals, ...) {
  args <- list(params = params, data = data, nrounds = nrounds, ...)
  evals_arg <- if (utils::packageVersion("xgboost") >= "2.0.0") "evals" else "watchlist"
  args[[evals_arg]] <- evals
  do.call(xgboost::xgb.train, args)
}

# Read the early-stopping iteration off a fitted booster.
#
# Up to xgboost 1.7 the booster was a plain list carrying `best_iteration`.
# From 2.0.0 it is an ALTLIST whose only element is `ptr`, so `model$best_iteration`
# returns NULL and the selected iteration lives in the `early_stop` attribute.
# Reading it through this helper matters: the previous `m$best_iteration %||%
# nrounds_max` fallback would silently resolve to nrounds_max on xgboost >= 2.0,
# i.e. early stopping would be disabled without any error being raised.
.of_xgb_best_iteration <- function(model) {
  bi <- model$best_iteration
  if (is.null(bi)) bi <- attributes(model)$early_stop$best_iteration
  bi
}

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

# =============================================================================
# Gate 1D.8 (2026-06-09): SHARED `_isNA` missingness-indicator synthesis.
#
# THE SINGLE place every supervised preprocessing path (OOF, FINAL, scoring;
# legacy AND nested_refit) creates `<feature>_isNA` companions. Before this gate
# the synthesis lived ONLY inside build_design_matrix_patches() (the OOF design
# matrix), so OOF trained on base + `_isNA` (~102 cols) while FINAL / scoring saw
# only the base whitelist (~51 cols) -- the deployed model and the OOF
# diagnostics ran on DIFFERENT feature spaces. Routing every path through this
# one helper restores OOF == FINAL == scoring parity.
#
# CANONICAL `_isNA` SEMANTICS (code + roxygen on apply_supervised_recipe):
#   * SCOPE: an `_isNA` companion is created for EVERY MEASURED NUMERIC feature
#     (matching the historical B5 / build_design_matrix_patches scope). Logical
#     columns are coerced to integer first and therefore ALSO get a companion
#     (e.g. hotspot_available_isNA, which is all-0 because the flag is never NA).
#     CATEGORICAL (character/factor) features get NO `_isNA` companion (they
#     carry a sentinel level instead) -- this matches the prior behaviour.
#   * VALUE: feature_isNA <- as.integer(is.na(feature)), computed BEFORE
#     imputation, row-local (no cross-row leakage).
#   * DEDUPE / NO SECOND ORDER: a column already ending in `_isNA` never gets a
#     `<x>_isNA_isNA` companion; an `_isNA` column that arrives in the input is
#     REGENERATED from its base feature (if the base is present) and never
#     duplicated.
#   * FEATURE ABSENT / ALL-NA: handled upstream (the base column exists, possibly
#     all-NA) so its companion is simply all-1; nothing special here.
#
# ORDER: base feature columns are kept in their incoming order, then ALL `_isNA`
# companions are appended in base-feature order. This deterministic
# base-block-then-indicator-block layout is identical for every path, which is
# what makes hash(feature_order_OOF) == hash(FINAL) == hash(SCORING).
#
# @param X_df coerced feature data.frame (numeric/integer/factor; logicals
#   already coerced to integer by the caller, or coerce here defensively).
# @return X_df with `<feature>_isNA` companions appended (base block first).
# @keywords internal
# @noRd
.of_synthesize_isna_companions <- function(X_df) {
  base_names <- names(X_df)
  # Numeric/integer base features that are NOT themselves `_isNA` flags.
  isna_targets <- base_names[vapply(base_names, function(nm) {
    if (grepl("_isNA$", nm)) return(FALSE)
    is.numeric(X_df[[nm]]) || is.integer(X_df[[nm]]) || is.logical(X_df[[nm]])
  }, logical(1))]

  companions <- list()
  for (nm in isna_targets) {
    flag_nm <- paste0(nm, "_isNA")
    companions[[flag_nm]] <- as.integer(is.na(X_df[[nm]]))
  }

  # Drop any pre-existing `_isNA` columns whose base feature is present (they are
  # REGENERATED above; never trust / duplicate an input-provided companion). An
  # `_isNA` column whose base feature is ABSENT is left untouched (it is itself a
  # base feature from the recipe's point of view).
  preexisting_isna <- base_names[grepl("_isNA$", base_names)]
  regenerated <- intersect(preexisting_isna, names(companions))
  keep_base <- setdiff(base_names, regenerated)

  out <- X_df[, keep_base, drop = FALSE]
  for (flag_nm in names(companions)) {
    out[[flag_nm]] <- companions[[flag_nm]]
  }
  out
}

#' Fit the SHARED supervised preprocessing recipe on TRAINING data only.
#'
#' Gate 1D.8 (2026-06-09): the ONE recipe object that is the single source of
#' truth for the supervised feature space across OOF, FINAL and scoring. It is
#' fit on TRAINING rows only (no held-out / scoring data) and then APPLIED
#' (never re-fit) to any frame via [apply_supervised_recipe()].
#'
#' CANONICAL PREPROCESSING ORDER (identical everywhere):
#'   original features -> resolve base-feature whitelist -> create
#'   `<feature>_isNA` indicators (BEFORE imputation) -> fit imputation medians /
#'   factor levels on TRAINING -> canonical final column order
#'   (base block, then `_isNA` block).
#'
#' MISSING POLICY: each measured numeric base feature gets an
#' `<feature>_isNA = as.integer(is.na(feature))` indicator and is then imputed to
#' the TRAINING median (`impute_numeric = "median"`; "zero" forces 0). An all-NA
#' training column has an NA median -> the indicator is all-1 and the base is left
#' NA so xgboost treats it as missing (degenerate rule). `_isNA` scope = all
#' measured numeric whitelist features (logicals coerced to integer first, so they
#' get a -- typically all-0 -- companion); categorical features get NO `_isNA`.
#'
#' WEIGHTS POLICY: `_isNA` indicators carry the CANONICAL weight 1.0 unless the
#' caller's `feature_weights` names them explicitly (no automatic inheritance from
#' the base feature). This is applied identically in OOF and FINAL.
#'
#' @param training_df data.frame of TRAINING rows carrying the base feature
#'   columns (`base_features`).
#' @param base_features character; the base whitelist features (no `_isNA`).
#' @param impute_numeric "median" (default) or "zero".
#' @param impute_factor_missing sentinel level for character/factor NAs.
#' @return a recipe list with: `base_features`, `missing_indicator_features`,
#'   `final_feature_order`, `numeric_medians` (named list; NA = degenerate),
#'   `impute_numeric`, `impute_factor_missing`. Apply it with
#'   [apply_supervised_recipe()].
#' @keywords internal
#' @noRd
fit_supervised_recipe <- function(training_df, base_features,
                                  impute_numeric = c("median", "zero"),
                                  impute_factor_missing = "MISSING") {
  impute_numeric <- match.arg(impute_numeric)
  base_features <- intersect(base_features, names(training_df))
  X_raw <- .of_nested_coerce_features(training_df[, base_features, drop = FALSE],
                                      impute_factor_missing, synthesize_isna = TRUE)
  medians <- .of_nested_fit_medians(X_raw, seq_len(nrow(X_raw)), impute_numeric)
  final_order <- names(X_raw)
  list(
    base_features              = base_features,
    missing_indicator_features = final_order[grepl("_isNA$", final_order)],
    final_feature_order        = final_order,
    numeric_medians            = medians,
    impute_numeric             = impute_numeric,
    impute_factor_missing      = impute_factor_missing
  )
}

#' Apply a fitted supervised recipe to ANY frame (no re-fit).
#'
#' Gate 1D.8 (2026-06-09): the deterministic transform half of the shared recipe.
#' Synthesises the SAME `<feature>_isNA` indicators, imputes with the recipe's
#' TRAINING medians (degenerate -> left NA -> xgboost-missing), and returns an
#' NA-preserving design matrix in the recipe's canonical column order. Used
#' identically by OOF (outer-test), FINAL (refit) and scoring.
#'
#' @param df data.frame to transform (any rows; base features may be absent /
#'   all-NA -- the indicator is created and the base follows the missing policy).
#' @param recipe a recipe from [fit_supervised_recipe()] (or the persisted FINAL
#'   recipe, which carries the same fields under `$impute` / `$cols`).
#' @return a base matrix aligned to `recipe$final_feature_order`, NA-preserving.
#' @keywords internal
#' @noRd
apply_supervised_recipe <- function(df, recipe) {
  ifm <- recipe$impute_factor_missing %||% "MISSING"
  meds <- recipe$numeric_medians %||% list()
  ref  <- recipe$final_feature_order
  base <- recipe$base_features %||% setdiff(ref, grep("_isNA$", ref, value = TRUE))
  cols <- intersect(base, names(df))
  X_raw <- .of_nested_coerce_features(df[, cols, drop = FALSE], ifm,
                                      synthesize_isna = TRUE)
  X_imp <- .of_nested_apply_medians(X_raw, meds)
  .of_nested_build_matrix(X_imp, ref_cols = ref)
}

# Prepare a data.frame of feature columns (factor/character/logical handling)
# WITHOUT imputing numeric NAs. This mirrors the non-numeric coercions performed
# in train_final_model_direct() (lines ~327-339) so the downstream
# sparse.model.matrix() behaves identically, but leaves numeric NAs in place so
# imputation can be fit train-only afterwards.
#
# Gate 1D.8 (2026-06-09): after the factor/logical coercions, this now ALSO
# synthesises the shared `<feature>_isNA` companions (.of_synthesize_isna_companions)
# so OOF, FINAL and scoring build the SAME feature space. The `_isNA` flags are
# created BEFORE imputation (the caller fits/applies medians afterwards).
#
# Returns the coerced data.frame (numeric columns still carry NA where missing,
# non-finite values coerced to NA, `_isNA` companions appended).
#
# @param synthesize_isna logical; when TRUE (default) append the shared `_isNA`
#   companions. The scoring path passes TRUE too (its reconciler creates the
#   BASE columns; the companions are synthesised here, identically to training).
# @keywords internal
# @noRd
.of_nested_coerce_features <- function(X_df, impute_factor_missing,
                                       synthesize_isna = TRUE) {
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
  if (isTRUE(synthesize_isna)) {
    X_df <- .of_synthesize_isna_companions(X_df)
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
#' @param sample_weights PHASE 2 (artifact_hard): optional per-ROW weight vector
#'   aligned to the rows of `train_df` (length == nrow(train_df)) OR NULL. When
#'   non-NULL it is applied via `xgboost::setinfo(dmat, "weight", w)` right after
#'   each `xgb.DMatrix(...)` (inner-train, inner-val, and the full refit), and
#'   `scale_pos_weight` is recomputed from the WEIGHTED class sums
#'   (`sum(w[neg]) / sum(w[pos])`) so the positives are not double-weighted. When
#'   NULL, NO weight vector is set and spw uses the raw counts -- byte-identical
#'   to the pre-Phase-2 behaviour.
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
                                 sample_weights = NULL,
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

  # PHASE 2: validate the optional per-row sample weights (aligned to train_df).
  use_sample_weights <- !is.null(sample_weights)
  if (use_sample_weights) {
    sample_weights <- as.numeric(sample_weights)
    if (length(sample_weights) != nrow(train_df)) {
      stop(".of_nested_refit_fit(): sample_weights must align to train_df rows ",
           "(length ", length(sample_weights), " != nrow ", nrow(train_df), ").",
           call. = FALSE)
    }
    if (any(!is.finite(sample_weights)) || any(sample_weights < 0)) {
      stop(".of_nested_refit_fit(): sample_weights must be finite and >= 0.",
           call. = FALSE)
    }
  }

  feature_cols <- intersect(feature_cols, names(train_df))
  if (length(feature_cols) == 0L) {
    stop(".of_nested_refit_fit(): no feature columns present in train_df.")
  }

  # Gate 1D.8: the SHARED recipe synthesises `_isNA` companions itself
  # (.of_nested_coerce_features -> .of_synthesize_isna_companions). Strip any
  # `_isNA` companion that arrives in feature_cols whose BASE feature is also
  # present, so we never carry an input-provided companion into the matrix (it is
  # regenerated). base_features is the recipe's canonical pre-synthesis set.
  base_features <- feature_cols[!(grepl("_isNA$", feature_cols) &
                                  sub("_isNA$", "", feature_cols) %in% feature_cols)]

  y_all <- as.integer(train_df[[label_col]] == "burned")
  n_all <- nrow(train_df)
  n_pos <- sum(y_all == 1L)
  n_neg <- sum(y_all == 0L)

  # Coerce non-numeric features once (factor/char/logical) WITHOUT imputing, then
  # synthesise the shared `_isNA` companions (identical for OOF / FINAL / scoring).
  X_df_raw <- .of_nested_coerce_features(train_df[, base_features, drop = FALSE],
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
  # PHASE 2: per-row weights (inner-train slice). NULL -> no weight set.
  if (use_sample_weights) {
    xgboost::setinfo(d_inner_tr, "weight", sample_weights[inner_tr])
  }
  d_inner_val <- xgboost::xgb.DMatrix(X_val_mat, label = y_val, missing = NA)
  if (use_sample_weights) {
    xgboost::setinfo(d_inner_val, "weight", sample_weights[inner_val])
  }
  apply_fw(d_inner_tr,  sel_x_cols)
  apply_fw(d_inner_val, sel_x_cols)

  # (d) scale_pos_weight on inner_tr only. PHASE 2: when per-row weights are
  # active, compute spw from the WEIGHTED class sums (sum(w[neg])/sum(w[pos]))
  # so the positives are not double-weighted; when NULL, the raw count ratio
  # (byte-identical to today).
  if (use_sample_weights) {
    w_tr <- sample_weights[inner_tr]
    spw_sel <- sum(w_tr[y_tr == 0L]) / max(.Machine$double.eps, sum(w_tr[y_tr == 1L]))
  } else {
    spw_sel <- sum(y_tr == 0L) / max(1, sum(y_tr == 1L))
  }
  params_sel <- params_fn(spw_sel)

  # (e) train with inner_val as the ONLY early-stopping set.
  set.seed(fold_seed)
  m_sel <- .of_xgb_train_evals(
    params = params_sel,
    data = d_inner_tr,
    nrounds = nrounds_max,
    evals = list(train = d_inner_tr, val = d_inner_val),
    early_stopping_rounds = early_stopping_rounds,
    verbose = if (isTRUE(verbose)) 1 else 0
  )
  best_iteration <- .of_xgb_best_iteration(m_sel) %||% nrounds_max
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
  # PHASE 2: per-row weights (full refit). NULL -> no weight set.
  if (use_sample_weights) {
    xgboost::setinfo(d_refit, "weight", sample_weights)
  }
  apply_fw(d_refit, refit_x_cols)

  # (g) scale_pos_weight on the FULL train_df. PHASE 2: weighted sums when
  # per-row weights are active (avoids double weighting); raw counts otherwise.
  if (use_sample_weights) {
    spw_refit <- sum(sample_weights[y_all == 0L]) /
      max(.Machine$double.eps, sum(sample_weights[y_all == 1L]))
  } else {
    spw_refit <- n_neg / max(1, n_pos)
  }
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

  # Gate 1D.8: the recipe's column sets. base_features = the pre-synthesis
  # whitelist features; missing_indicator_features = the synthesised `_isNA`
  # companions; final_feature_order = the canonical post-synthesis feature order
  # (base block then `_isNA` block, BEFORE any factor one-hot expansion).
  # feature_cols (returned + persisted into recipe$cols$feature_cols, which the
  # scoring reconciler aligns to) is the FULL post-synthesis set in canonical
  # order, so OOF/FINAL/scoring all share ONE feature space.
  final_feature_order <- names(X_df_raw)
  missing_indicator_features <- final_feature_order[grepl("_isNA$", final_feature_order)]
  feature_cols_full <- final_feature_order

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
    # Gate 1D.8: feature_cols is now the FULL post-synthesis set (base + `_isNA`),
    # in canonical order -- the single shared feature space. The recipe also
    # records the base / indicator partition separately.
    feature_cols   = feature_cols_full,
    base_features  = base_features,
    missing_indicator_features = missing_indicator_features,
    final_feature_order        = final_feature_order,
    spw_selection  = spw_sel,
    spw_refit      = spw_refit,
    inner_split    = spl,
    degenerate_cols = degenerate_cols,
    audit          = audit
  )
}

# =============================================================================
# Gate 1C.3 (2026-06-08): recipe-driven scoring schema reconciliation.
#
# The H BLOCKER: the supervised SCORING path must work when hotspots = NULL and
# when one or several features are entirely NA -- WITHOUT dropping rows and
# WITHOUT reconstructing the feature schema from the scoring-year data. The
# canonical schema source is the RECIPE saved during the FINAL refit
# (recipe$cols$feature_cols / recipe$cols$x_cols / recipe$impute$numeric_medians
# / recipe$impute$impute_factor_missing), NEVER re-derived from the scoring year.
#
# This helper reconciles a scoring data.frame against the SAVED recipe feature
# schema and returns the reconciled feature data.frame PLUS an explicit,
# REGISTERED reconciliation record (one row per recorded decision; no silent
# corrections). The five registered cases:
#
#   1. EXPECTED FEATURE ABSENT  -> create it as NA_real_ (numeric) so it becomes
#      xgboost-missing; recorded as action = "create_absent_numeric_NA". (Absent
#      columns are created NA -- NOT zero -- so the all-NA-feature / hotspots=NULL
#      path is treated as missing, exactly as the refit model was trained.)
#   2. FEATURE FULLY NA         -> kept as-is (NA preserved). Downstream median
#      imputation (.of_nested_apply_medians) fills it from the recipe median, or
#      leaves it NA when the recipe median is itself degenerate (xgboost-missing).
#      Recorded as action = "feature_fully_NA".
#   3. NEW/UNKNOWN CATEGORICAL LEVEL -> mapped to the recipe sentinel
#      (impute_factor_missing) when the level is unseen, never silently invented.
#      Recorded as action = "remap_unknown_level".
#   4. UNEXPECTED EXTRA COLUMN  -> dropped (the builder aligns to recipe$x_cols so
#      extras cannot shift the matrix), recorded as action = "drop_extra_column".
#   5. INCOMPATIBLE TYPE        -> coerced per the recipe (numeric expected ->
#      as.numeric; non-coercible -> NA + flagged), recorded as
#      action = "coerce_type" or "coerce_type_failed".
#
# @param df scoring feature data.frame (geometry already dropped).
# @param recipe the saved FINAL-refit recipe (must carry $cols$feature_cols and
#   $impute$impute_factor_missing; $cols$x_cols used by the caller).
# @return list(df = reconciled feature data.frame restricted to recipe
#   feature_cols in recipe order, record = data.frame reconciliation log).
# @keywords internal
# @noRd
.of_reconcile_scoring_schema <- function(df, recipe) {
  feature_cols <- recipe$cols$feature_cols
  impute_factor_missing <- recipe$impute$impute_factor_missing %||% "MISSING"
  # The recipe records numeric_medians only for non-degenerate columns; the
  # presence/absence of a name there is provenance, not a coercion target here.

  rec <- list()
  add_rec <- function(column, case, action, detail = NA_character_) {
    rec[[length(rec) + 1L]] <<- data.frame(
      column = column, case = case, action = action,
      detail = detail, stringsAsFactors = FALSE
    )
  }

  present <- names(df)

  # CASE 4: unexpected extra columns (anything not in the recipe feature schema).
  # They are simply not selected below; record each so the drop is not silent.
  extra <- setdiff(present, feature_cols)
  for (nm in extra) add_rec(nm, "extra_column", "drop_extra_column")

  # CASE 1: expected feature absent -> create as NA_real_ (xgboost-missing).
  missing_feat <- setdiff(feature_cols, present)
  for (nm in missing_feat) {
    df[[nm]] <- NA_real_
    add_rec(nm, "expected_absent", "create_absent_numeric_NA")
  }

  # Restrict + reorder to the recipe feature schema (drops the extras).
  df <- df[, feature_cols, drop = FALSE]

  # The supervised feature schema is NUMERIC by construction (the design-matrix
  # builder uses cat_cols = character(0); no categorical predictor survives into
  # the whitelist). A recipe feature therefore carries no factor levels, and any
  # character/factor arriving at scoring is an INCOMPATIBLE TYPE for a numeric
  # feature -> coerce to numeric (CASE 5), never treated as a new categorical
  # level. If a future build introduces a genuine categorical recipe feature
  # (recipe$cols$factor_levels), the unseen-level remap (CASE 3) would apply
  # instead; that branch is retained for forward-compatibility.
  # NOTE: use an explicit NULL test (not `%||%`), because `%||%` treats a
  # length-1 list whose element is multi-valued as a coercion target and would
  # error on `is.na()`; factor_levels is a named list of level vectors.
  recipe_factor_levels <- if (is.null(recipe$cols$factor_levels)) {
    list()
  } else {
    recipe$cols$factor_levels
  }

  for (nm in feature_cols) {
    col <- df[[nm]]
    is_recipe_factor <- nm %in% names(recipe_factor_levels)
    if ((is.character(col) || is.factor(col)) && is_recipe_factor) {
      # CASE 3: genuine categorical recipe feature -> map unseen levels to the
      # recipe sentinel (never invent a level).
      lv <- recipe_factor_levels[[nm]]
      unseen <- setdiff(unique(as.character(col)), c(lv, NA))
      if (length(unseen) > 0L) {
        col_chr <- as.character(col)
        col_chr[col_chr %in% unseen] <- impute_factor_missing
        df[[nm]] <- col_chr
        add_rec(nm, "categorical", "remap_unknown_level",
                detail = paste0("unseen: ", paste(unseen, collapse = ";")))
      }
    } else if (is.character(col) || is.factor(col)) {
      # CASE 5 (numeric schema): a character/factor on a NUMERIC feature is an
      # incompatible type -> coerce to numeric (uncoercible cells -> NA, which
      # the median recipe then imputes / leaves xgboost-missing).
      coerced <- suppressWarnings(as.numeric(as.character(col)))
      df[[nm]] <- coerced
      if (anyNA(coerced) & !all(is.na(coerced))) {
        add_rec(nm, "incompatible_type", "coerce_type",
                "char/factor on numeric feature -> numeric (some cells NA)")
      } else {
        add_rec(nm, "incompatible_type", "coerce_type",
                "char/factor on numeric feature -> numeric")
      }
    } else if (is.logical(col)) {
      df[[nm]] <- as.integer(col)
      add_rec(nm, "incompatible_type", "coerce_type", "logical->integer")
    } else if (!is.numeric(col)) {
      # CASE 5: incompatible type -> coerce per recipe (numeric expected).
      coerced <- suppressWarnings(as.numeric(as.character(col)))
      if (all(is.na(coerced)) && !all(is.na(col))) {
        df[[nm]] <- rep(NA_real_, length(col))
        add_rec(nm, "incompatible_type", "coerce_type_failed",
                "uncoercible to numeric; set NA")
      } else {
        df[[nm]] <- coerced
        add_rec(nm, "incompatible_type", "coerce_type", "->numeric")
      }
    }
    # CASE 2: feature fully NA (after the above) -> record; value left NA so the
    # median recipe imputes it (or leaves it xgboost-missing if degenerate).
    if (all(is.na(df[[nm]]))) {
      add_rec(nm, "feature_fully_NA", "feature_fully_NA")
    }
  }

  record <- if (length(rec) > 0L) {
    do.call(rbind, rec)
  } else {
    data.frame(column = character(0), case = character(0),
               action = character(0), detail = character(0),
               stringsAsFactors = FALSE)
  }
  list(df = df, record = record)
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
