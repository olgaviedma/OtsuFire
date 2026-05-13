run_oof_xgb <- function(
    XL_mat, y, labelled_df,
    fold_cols = c("fold_rep1","fold_rep2"),
    params,
    nrounds_max = 3000,
    early_stop = 75,
    seed_base = 42,
    out_dir = NULL,
    prefix = "oof",
    verbose = 0,
    # 0.5.0: per-feature weights (full vector aligned to colnames(XL_mat),
    # already padded to 1.0 by the caller) applied via xgb.DMatrix
    # `feature_weights` info on every per-fold dtrain / dtest. NULL
    # leaves XGBoost's internal default (uniform 1.0) untouched.
    feature_weights_vector = NULL
) {
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("Instala xgboost")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Instala dplyr")

  stopifnot(length(y) == nrow(XL_mat))
  stopifnot("fire_uid" %in% names(labelled_df), "class" %in% names(labelled_df))

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
    write.csv(oof_long, csv_long, row.names = FALSE)
    write.csv(oof_agg,  csv_agg,  row.names = FALSE)
  }
  
  list(oof_long = oof_long, oof_agg = oof_agg)
}
