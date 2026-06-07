# =============================================================================
# OOF per-fold trainer.
#
# B1 (2026-06-07): adds the `training_protocol = "nested_refit"` path. For
# "legacy" (default) the per-fold body is byte-identical to the historical
# implementation. For "nested_refit" each OUTER fold:
#   - takes outer_train = idx_tr, outer_test = idx_te (never sub-sampled, never
#     imputed-from);
#   - ASSERTS no fire_uid / block_id overlap between train and test;
#   - applies bucket caps to the outer_train NEGATIVES only (replicating the
#     FINAL bucketing exactly), under fold_seed = seed_base + 1000*r + k
#     (skipped when oof_sampling == "full");
#   - calls the SHARED core .of_nested_refit_fit() on the capped outer_train ->
#     refit model + refit medians (fit inner-train only; inner-val is the only
#     early-stopping set; refit on all outer_train at best_iteration);
#   - transforms the FULL outer_test with the REFIT medians and predicts it;
#   - records a per-fold audit row.
# The 4 cap ratios are REQUIRED formals (no defaults) on the nested path so a
# dropped argument cannot silently revert a bucket to ratio 1.0.
# =============================================================================
run_oof_xgb <- function(
    XL_mat, y, labelled_df,
    fold_cols = c("fold_rep1","fold_rep2"),
    params,
    # Gate 1B (2026-06-07): nrounds_max / early_stop / seed_base are REQUIRED
    # resolved args with NO methodological defaults. Single source of truth =
    # cfg$train_control, threaded down by run_dm_oof_pipeline(). A dropped arg
    # ERRORS in the guard block below.
    nrounds_max,
    early_stop,
    seed_base,
    out_dir = NULL,
    prefix = "oof",
    verbose = 0,
    # 0.5.0: per-feature weights (full vector aligned to colnames(XL_mat),
    # already padded to 1.0 by the caller) applied via xgb.DMatrix
    # `feature_weights` info on every per-fold dtrain / dtest. NULL
    # leaves XGBoost's internal default (uniform 1.0) untouched.
    feature_weights_vector = NULL,
    # 2026-06-06 (partial-write overwrite fix): honor overwrite for the OOF
    # long/agg CSV sidecars. overwrite=TRUE writes exactly as before
    # (byte-identical); overwrite=FALSE skips when the target already exists.
    overwrite = TRUE,
    # B1 (2026-06-07): protocol toggle. "legacy" (default) = historical OOF.
    # "nested_refit" = leakage-free per-fold protocol (see header).
    training_protocol = c("legacy", "nested_refit"),
    # B1: "capped" applies the FINAL bucket caps to the outer-train negatives;
    # "full" uses all outer-train rows. Only consulted on the nested path.
    oof_sampling = c("capped", "full"),
    # B1: prepared (deferred-impute) labelled feature frame + model column set,
    # produced by build_design_matrix_patches(defer_impute=TRUE). REQUIRED on
    # the nested path; ignored on the legacy path.
    prepared_labelled = NULL,
    model_cols = NULL,
    class_col = "class",
    group_col = "block_id",
    # Gate 1B (2026-06-07): val_frac / impute_* are REQUIRED resolved args with
    # NO methodological defaults (single source of truth = cfg$train_control).
    val_frac,
    impute_numeric,
    impute_factor_missing,
    feature_weights = NULL,
    # B1: bucket cap ratios. REQUIRED on the nested path (no defaults) so a
    # dropped argument cannot silently revert a bucket to ratio 1.0.
    contextual_exclusion_to_burned_ratio,
    spectral_hard_negative_to_burned_ratio,
    random_to_burned_ratio,
    otsu_unburned_to_burned_ratio,
    # B1: bucket source / neg_type definitions, mirroring the FINAL stage.
    deterministic_drop_source = c("deterministic_drop_hard"),
    spectral_hard_negative_neg_types = c("spectral_reject_medium"),
    random_background_source = c("random_burnable_background"),
    otsu_unburned_source = c("otsu_patch_residual"),
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep")
) {
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("Instala xgboost")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Instala dplyr")
  # Gate 1B (2026-06-07): required-arg guard. nrounds_max / early_stop /
  # seed_base drive BOTH paths, so they are always required (no silent
  # defaults; single source = cfg$train_control).
  for (.nm in c("nrounds_max", "early_stop", "seed_base")) {
    if (eval(call("missing", as.name(.nm)))) {
      stop("run_oof_xgb(): required resolved arg '", .nm,
           "' is missing (no methodological default; threaded from ",
           "cfg$train_control via run_dm_oof_pipeline()).", call. = FALSE)
    }
  }
  training_protocol <- match.arg(training_protocol)
  oof_sampling <- match.arg(oof_sampling)
  # val_frac / impute_* are consumed ONLY on the nested_refit path (the per-fold
  # leakage-free core); required there, mirroring the cap-ratio pattern below.
  if (identical(training_protocol, "nested_refit")) {
    for (.nm in c("val_frac", "impute_numeric", "impute_factor_missing")) {
      if (eval(call("missing", as.name(.nm)))) {
        stop("run_oof_xgb(nested_refit): required resolved arg '", .nm,
             "' is missing (no methodological default).", call. = FALSE)
      }
    }
  }

  stopifnot("fire_uid" %in% names(labelled_df), "class" %in% names(labelled_df))

  if (identical(training_protocol, "legacy")) {
    # -------------------------------------------------------------------------
    # LEGACY path (byte-identical to the historical implementation).
    # -------------------------------------------------------------------------
    stopifnot(length(y) == nrow(XL_mat))

    if (!is.null(feature_weights_vector)) {
      stopifnot(length(feature_weights_vector) == ncol(XL_mat))
    }

    oof_rows <- list()

    for (r in seq_along(fold_cols)) {
      fold_col <- fold_cols[r]
      stopifnot(fold_col %in% names(labelled_df))

      folds <- as.integer(labelled_df[[fold_col]])
      stopifnot(!anyNA(folds))

      for (k in sort(unique(folds))) {
        idx_te <- which(folds == k)
        idx_tr <- which(folds != k)

        dtrain <- xgboost::xgb.DMatrix(XL_mat[idx_tr, , drop=FALSE], label = y[idx_tr])
        dtest  <- xgboost::xgb.DMatrix(XL_mat[idx_te, , drop=FALSE], label = y[idx_te])

        if (!is.null(feature_weights_vector)) {
          xgboost::setinfo(dtrain, "feature_weights", feature_weights_vector)
          xgboost::setinfo(dtest,  "feature_weights", feature_weights_vector)
        }

        set.seed(seed_base + 1000*r + k)

        m <- xgboost::xgb.train(
          params = params,
          data = dtrain,
          nrounds = nrounds_max,
          watchlist = list(train = dtrain, val = dtest),
          early_stopping_rounds = early_stop,
          verbose = ifelse(verbose > 0, 1, 0)
        )

        p <- predict(m, dtest)

        oof_rows[[length(oof_rows) + 1]] <- data.frame(
          fire_uid = labelled_df$fire_uid[idx_te],
          class    = labelled_df$class[idx_te],
          rep      = r,
          fold     = k,
          p_burned = p,
          best_iter = m$best_iteration,
          stringsAsFactors = FALSE
        )
      }
    }

    oof_long <- do.call(rbind, oof_rows)
    oof_audit <- NULL
  } else {
    # -------------------------------------------------------------------------
    # NESTED_REFIT path (B1).
    # -------------------------------------------------------------------------
    if (is.null(prepared_labelled) || is.null(model_cols)) {
      stop("run_oof_xgb(nested_refit): 'prepared_labelled' and 'model_cols' are ",
           "required (from build_design_matrix_patches(defer_impute=TRUE)).",
           call. = FALSE)
    }
    # Force the 4 cap ratios: required formals -> a dropped arg ERRORS here
    # (cannot silently revert a bucket to ratio 1.0). missing() must be called
    # directly on each formal name.
    if (missing(contextual_exclusion_to_burned_ratio)) {
      stop("run_oof_xgb(nested_refit): required cap 'contextual_exclusion_to_burned_ratio' is missing.", call. = FALSE)
    }
    if (missing(spectral_hard_negative_to_burned_ratio)) {
      stop("run_oof_xgb(nested_refit): required cap 'spectral_hard_negative_to_burned_ratio' is missing.", call. = FALSE)
    }
    if (missing(random_to_burned_ratio)) {
      stop("run_oof_xgb(nested_refit): required cap 'random_to_burned_ratio' is missing.", call. = FALSE)
    }
    if (missing(otsu_unburned_to_burned_ratio)) {
      stop("run_oof_xgb(nested_refit): required cap 'otsu_unburned_to_burned_ratio' is missing.", call. = FALSE)
    }
    stopifnot(nrow(prepared_labelled) == nrow(labelled_df))

    feat_model_cols <- intersect(model_cols, names(prepared_labelled))

    has_source   <- "source"   %in% names(labelled_df)
    has_neg_type <- "neg_type" %in% names(labelled_df)
    has_block    <- group_col  %in% names(labelled_df)

    # Per-fold negative-bucket capping (mirrors FINAL exactly, applied to the
    # OUTER-TRAIN negatives only).
    cap_outer_train <- function(tr_idx, fold_seed) {
      if (identical(oof_sampling, "full")) {
        return(list(keep = tr_idx, audit = NULL))
      }
      cls <- as.character(labelled_df[[class_col]])[tr_idx]
      src <- if (has_source) as.character(labelled_df[["source"]])[tr_idx] else rep(NA_character_, length(tr_idx))
      ngt <- if (has_neg_type) as.character(labelled_df[["neg_type"]])[tr_idx] else rep(NA_character_, length(tr_idx))

      is_burned <- cls == "burned"
      n_burned  <- sum(is_burned)

      pick <- function(mask, ratio) {
        avail <- which(mask)
        if (length(avail) == 0L || !is.finite(ratio)) return(avail)
        target <- ceiling(n_burned * ratio)
        n_take <- min(target, length(avail))
        if (n_take <= 0L) return(integer(0))
        if (n_take >= length(avail)) return(avail)
        sort(sample(avail, n_take))
      }

      ctx_mask <- (!is_burned) & has_source & (src %in% deterministic_drop_source) &
        !(ngt %in% spectral_hard_negative_neg_types)
      shn_mask <- (!is_burned) & has_source & (src %in% deterministic_drop_source) &
        (ngt %in% spectral_hard_negative_neg_types)
      rb_mask  <- (!is_burned) & has_source & (src %in% random_background_source)
      otsu_mask <- (!is_burned) & has_source & (src %in% otsu_unburned_source) &
        !(ngt %in% otsu_unburned_exclude_neg_types)

      set.seed(fold_seed)
      sel_burned <- which(is_burned)
      sel_ctx  <- pick(ctx_mask, contextual_exclusion_to_burned_ratio)
      sel_shn  <- pick(shn_mask, spectral_hard_negative_to_burned_ratio)
      sel_rb   <- pick(rb_mask,  random_to_burned_ratio)
      sel_otsu <- pick(otsu_mask, otsu_unburned_to_burned_ratio)
      # Any negatives not in a known bucket are kept as-is (parity with FINAL,
      # which only excludes via the otsu exclude-list; unknown-source negatives
      # would not be in L_ok at all, but here we keep them to avoid silently
      # dropping rows that FINAL never had a bucket for).
      known_mask <- ctx_mask | shn_mask | rb_mask | otsu_mask | is_burned
      sel_other <- which((!known_mask))

      keep_local <- sort(unique(c(sel_burned, sel_ctx, sel_shn, sel_rb, sel_otsu, sel_other)))
      audit <- list(
        n_burned            = n_burned,
        contextual_available = sum(ctx_mask), contextual_cap = ceiling(n_burned * contextual_exclusion_to_burned_ratio), contextual_selected = length(sel_ctx),
        spectral_available   = sum(shn_mask), spectral_cap = ceiling(n_burned * spectral_hard_negative_to_burned_ratio), spectral_selected = length(sel_shn),
        random_bg_available  = sum(rb_mask),  random_bg_cap = ceiling(n_burned * random_to_burned_ratio), random_bg_selected = length(sel_rb),
        otsu_available       = sum(otsu_mask), otsu_cap = ceiling(n_burned * otsu_unburned_to_burned_ratio), otsu_selected = length(sel_otsu),
        other_kept           = length(sel_other)
      )
      list(keep = tr_idx[keep_local], audit = audit)
    }

    oof_rows <- list()
    audit_rows <- list()

    for (r in seq_along(fold_cols)) {
      fold_col <- fold_cols[r]
      stopifnot(fold_col %in% names(labelled_df))
      folds <- as.integer(labelled_df[[fold_col]])
      stopifnot(!anyNA(folds))

      for (k in sort(unique(folds))) {
        idx_te <- which(folds == k)
        idx_tr <- which(folds != k)
        fold_seed <- as.integer(seed_base + 1000L * r + k)

        # ASSERT no overlap of fire_uid AND block_id between train and test.
        uid_tr <- as.character(labelled_df$fire_uid[idx_tr])
        uid_te <- as.character(labelled_df$fire_uid[idx_te])
        if (length(intersect(uid_tr, uid_te)) > 0L) {
          stop(sprintf("run_oof_xgb(nested_refit): fire_uid overlap between outer ",
                       "train/test (rep %d fold %d).", r, k), call. = FALSE)
        }
        if (has_block) {
          blk_tr <- as.character(labelled_df[[group_col]][idx_tr])
          blk_te <- as.character(labelled_df[[group_col]][idx_te])
          if (length(intersect(blk_tr, blk_te)) > 0L) {
            stop(sprintf("run_oof_xgb(nested_refit): %s overlap between outer ",
                         "train/test (rep %d fold %d).", group_col, r, k),
                 call. = FALSE)
          }
        }

        n_te_before <- length(idx_te)

        # Cap NEGATIVES of the outer-train only (outer-test untouched).
        capped <- cap_outer_train(idx_tr, fold_seed)
        keep_tr <- capped$keep

        # Build the train data.frame the core consumes: model features + class +
        # block col (for the inner split).
        cols_for_core <- unique(c(feat_model_cols, class_col,
                                  if (has_block) group_col else NULL))
        train_df <- prepared_labelled[keep_tr, , drop = FALSE]
        train_df[[class_col]] <- as.character(labelled_df[[class_col]])[keep_tr]
        if (has_block) train_df[[group_col]] <- labelled_df[[group_col]][keep_tr]
        train_df <- train_df[, intersect(cols_for_core, names(train_df)), drop = FALSE]

        fit <- .of_nested_refit_fit(
          train_df              = train_df,
          feature_cols          = feat_model_cols,
          label_col             = class_col,
          group_col             = group_col,
          block_col             = if (has_block) group_col else NULL,
          val_frac              = val_frac,
          params_fn             = .of_canonical_xgb_params,
          sampling_seed         = fold_seed,
          fold_seed             = fold_seed,
          feature_weights       = feature_weights,
          nrounds_max           = nrounds_max,
          early_stopping_rounds = early_stop,
          impute_numeric        = impute_numeric,
          impute_factor_missing = impute_factor_missing,
          verbose               = verbose > 0
        )

        # Transform the FULL outer_test with the REFIT medians (no refit) and
        # predict it.
        te_df <- prepared_labelled[idx_te, feat_model_cols, drop = FALSE]
        X_te <- .of_nested_transform(
          df = te_df, feature_cols = feat_model_cols,
          med = fit$medians, ref_x_cols = fit$x_cols,
          impute_factor_missing = impute_factor_missing
        )
        d_te <- xgboost::xgb.DMatrix(X_te, label = y[idx_te], missing = NA)
        if (!is.null(feature_weights)) {
          fw <- rep(1.0, length(fit$x_cols)); names(fw) <- fit$x_cols
          in_both <- intersect(names(feature_weights), fit$x_cols)
          if (length(in_both) > 0L) fw[in_both] <- as.numeric(feature_weights[in_both])
          xgboost::setinfo(d_te, "feature_weights", fw)
        }
        p <- predict(fit$model, d_te)

        oof_rows[[length(oof_rows) + 1]] <- data.frame(
          fire_uid = labelled_df$fire_uid[idx_te],
          class    = labelled_df$class[idx_te],
          rep      = r,
          fold     = k,
          p_burned = p,
          best_iter = fit$best_iteration,
          stringsAsFactors = FALSE
        )

        # Per-fold audit row.
        a <- fit$audit
        a$stage              <- "OOF"
        a$prefix             <- prefix
        a$rep                <- r
        a$fold               <- k
        a$oof_sampling       <- oof_sampling
        a$n_outer_train_pre  <- length(idx_tr)
        a$n_outer_train_post <- length(keep_tr)
        a$n_outer_test       <- n_te_before
        a$outer_test_capped       <- FALSE
        a$outer_test_used_for_fit <- FALSE
        a$outer_test_rows_unchanged <- (n_te_before == length(idx_te))
        if (!is.null(capped$audit)) {
          ca <- capped$audit
          a$n_burned             <- ca$n_burned
          a$contextual_available <- ca$contextual_available
          a$contextual_cap       <- ca$contextual_cap
          a$contextual_selected  <- ca$contextual_selected
          a$spectral_available   <- ca$spectral_available
          a$spectral_cap         <- ca$spectral_cap
          a$spectral_selected    <- ca$spectral_selected
          a$random_bg_available  <- ca$random_bg_available
          a$random_bg_cap        <- ca$random_bg_cap
          a$random_bg_selected   <- ca$random_bg_selected
          a$otsu_available       <- ca$otsu_available
          a$otsu_cap             <- ca$otsu_cap
          a$otsu_selected        <- ca$otsu_selected
        }
        audit_rows[[length(audit_rows) + 1]] <- a
      }
    }

    oof_long <- do.call(rbind, oof_rows)
    oof_audit <- dplyr::bind_rows(audit_rows)
  }

  oof_agg <- oof_long |>
    dplyr::group_by(fire_uid, class) |>
    dplyr::summarise(
      p_oof_mean = mean(p_burned, na.rm=TRUE),
      p_oof_med  = median(p_burned, na.rm=TRUE),
      p_oof_sd   = sd(p_burned, na.rm=TRUE),
      n_preds    = dplyr::n(),
      .groups="drop"
    )

  stopifnot(nrow(oof_agg) == nrow(labelled_df))
  stopifnot(!anyDuplicated(oof_agg$fire_uid))

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    csv_long <- file.path(out_dir, paste0(prefix, "_oof_long.csv"))
    csv_agg  <- file.path(out_dir, paste0(prefix, "_oof_agg.csv"))
    if (isTRUE(overwrite) || !file.exists(csv_long)) {
      write.csv(oof_long, csv_long, row.names = FALSE)
    }
    if (isTRUE(overwrite) || !file.exists(csv_agg)) {
      write.csv(oof_agg,  csv_agg,  row.names = FALSE)
    }
    # B1: emit the per-fold nested-refit audit CSV (nested path only).
    if (!is.null(oof_audit)) {
      csv_audit <- file.path(out_dir, paste0(prefix, "_oof_nested_refit_audit.csv"))
      if (isTRUE(overwrite) || !file.exists(csv_audit)) {
        write.csv(oof_audit, csv_audit, row.names = FALSE)
      }
    }
  }

  list(oof_long = oof_long, oof_agg = oof_agg, oof_audit = oof_audit)
}
