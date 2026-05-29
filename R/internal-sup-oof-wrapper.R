compute_oof_metrics_by_threshold <- function(
    oof_agg,
    class_col = "class",
    prob_col = "p_oof_mean",
    pos_label = "burned",
    thresholds = seq(0.05, 0.95, by = 0.05)
) {
  stopifnot(is.data.frame(oof_agg))
  stopifnot(all(c(class_col, prob_col) %in% names(oof_agg)))

  df <- oof_agg
  df[[class_col]] <- as.character(df[[class_col]])
  df <- df[df[[class_col]] %in% c(pos_label, "unburned"), , drop = FALSE]

  if (!nrow(df)) {
    stop("OOF metrics cannot be computed: no burned/unburned rows found.")
  }

  y_true <- ifelse(df[[class_col]] == pos_label, 1L, 0L)
  p <- suppressWarnings(as.numeric(df[[prob_col]]))
  keep <- is.finite(p)
  y_true <- y_true[keep]
  p <- p[keep]

  if (!length(p)) {
    stop("OOF metrics cannot be computed: no finite probabilities found.")
  }

  metric_rows <- lapply(thresholds, function(thr) {
    y_hat <- ifelse(p >= thr, 1L, 0L)

    tp <- sum(y_true == 1L & y_hat == 1L)
    tn <- sum(y_true == 0L & y_hat == 0L)
    fp <- sum(y_true == 0L & y_hat == 1L)
    fn <- sum(y_true == 1L & y_hat == 0L)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    f1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }
    balanced_accuracy <- mean(c(recall, specificity), na.rm = TRUE)

    data.frame(
      threshold = thr,
      tp = tp,
      fp = fp,
      fn = fn,
      tn = tn,
      accuracy = accuracy,
      precision = precision,
      recall = recall,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      f1 = f1,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(metric_rows)
}

#' Run the diagnostic supervised OOF pipeline on extracted features.
#'
#' OtsuFire 0.5.0: replaces the previous `additional_drop_cols`
#' deny-list mechanism with the symmetric `feature_whitelist_override`
#' / `feature_weights` API used by [train_final_model_direct()]. The
#' OOF stage must train under the same conditions as the final model --
#' otherwise OOF metrics report on a different model.
#'
#' @param feature_whitelist_override Character vector subset of
#'   `.supervised_feature_cols` to use as the active OOF feature
#'   space. NULL (default) means the full canonical whitelist.
#' @param feature_weights Named numeric vector of per-feature weights
#'   forwarded to xgboost via `xgb.DMatrix` `feature_weights` info on
#'   each per-fold dtrain / dtest. NULL leaves the default uniform 1.0.
#' @keywords internal
#' @noRd
run_dm_oof_pipeline <- function(
    labelled,
    burned_like,
    labelled_df,
    params,
    result_dir,
    target_year = 2022,
    fold_cols   = c("fold_rep1", "fold_rep2"),
    nrounds_max = 3000,
    early_stop  = 75,
    seed_base   = 42,
    verbose     = 1,
    out_dir_oof = file.path(result_dir, "05_OOF"),
    prefix      = paste0(target_year, "_patch"),
    id_cols     = c("fire_uid", "class", "source", "poly_id", "block_id", "fold_rep1", "fold_rep2"),
    # 0.4.0 (Agent H): `drop_regex` is DEPRECATED in the OOF stage.
    # The OOF wrapper now applies a whitelist filter against
    # `.supervised_feature_cols` (+ `_isNA` companions) to
    # `labelled` and `burned_like` BEFORE handing them to
    # `build_design_matrix_patches()`. The argument is retained for
    # backward-compatibility; non-default values emit a one-time
    # message and are otherwise ignored.
    drop_regex  = character(0),
    cat_cols    = c("eco_major"),
    hs_n_col    = "hs_used_n",
    hs_conf_col = "hs_conf_mean",
    hs_frp_col  = "hs_frp_max",
    median_from = "labelled",
    save_dir_dm = file.path(result_dir, "04_MATRIX"),
    save_prefix = paste0(target_year, "_patch"),
    overwrite   = TRUE,
    labelled_gpkg = NULL,
    labelled_layer = "train_features",
    # 0.5.0: see Roxygen above. Default NULL preserves canonical
    # behaviour byte-for-byte.
    feature_whitelist_override = NULL,
    feature_weights = NULL,
    ...
) {
  if (!exists("build_design_matrix_patches")) stop("No encuentro build_design_matrix_patches() cargada en el entorno.")
  if (!exists("run_oof_xgb")) stop("No encuentro run_oof_xgb() cargada en el entorno.")

  # 0.5.0 hard removal: `additional_drop_cols` / `extra_drop_cols`
  # were removed in favour of `feature_whitelist_override`. Catch
  # them via `...` to print a clear migration message instead of
  # silently ignoring them.
  .dots <- list(...)
  if ("additional_drop_cols" %in% names(.dots) ||
      "extra_drop_cols" %in% names(.dots)) {
    stop(
      "`additional_drop_cols` / `extra_drop_cols` were removed in ",
      "OtsuFire 0.5.0. Use `feature_whitelist_override` instead. ",
      "Pass a subset of `.supervised_feature_cols`. See NEWS.md and ",
      "the migration note in HANDOFF Section N+20.",
      call. = FALSE
    )
  }
  if (length(.dots) > 0L) {
    stop("Unused arguments passed to run_dm_oof_pipeline(): ",
         paste(names(.dots), collapse = ", "), call. = FALSE)
  }

  dir.create(save_dir_dm, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_oof, recursive = TRUE, showWarnings = FALSE)

  # 0.4.0 deprecation: warn if caller supplies a non-default drop_regex.
  if (length(drop_regex) > 0L) {
    message("run_dm_oof_pipeline(): `drop_regex` is ignored since ",
            "OtsuFire 0.4.0; the OOF column set is controlled by the ",
            "canonical whitelist `.supervised_feature_cols`.")
    drop_regex <- character(0)
  }

  # 0.5.0: validate `feature_whitelist_override` (strict).
  if (!is.null(feature_whitelist_override)) {
    if (!is.character(feature_whitelist_override) ||
        length(feature_whitelist_override) == 0) {
      stop("`feature_whitelist_override` must be a non-empty character vector.",
           call. = FALSE)
    }
    unknown <- setdiff(feature_whitelist_override, .supervised_feature_cols)
    if (length(unknown) > 0) {
      stop(
        "`feature_whitelist_override` contains names not in the canonical ",
        "`.supervised_feature_cols`: ",
        paste(unknown, collapse = ", "), ". The canonical list is fixed; ",
        "this argument can only restrict it, not extend it.",
        call. = FALSE
      )
    }
    feature_whitelist_override <- unique(feature_whitelist_override)
  }
  active_whitelist <- if (is.null(feature_whitelist_override)) {
    .supervised_feature_cols
  } else {
    feature_whitelist_override
  }

  # 0.5.0: validate `feature_weights` (named numeric, finite, >= 0).
  if (!is.null(feature_weights)) {
    if (!is.numeric(feature_weights) || is.null(names(feature_weights))) {
      stop("`feature_weights` must be a named numeric vector.",
           call. = FALSE)
    }
    if (any(is.na(feature_weights)) || any(feature_weights < 0)) {
      stop("`feature_weights` values must be finite and >= 0.",
           call. = FALSE)
    }
  }

  # 0.4.0 architectural refactor + 0.5.0 override: whitelist filter
  # the OOF inputs. Preserve id_cols (fold-assignment metadata,
  # identifiers, class column) and the active sf geometry, plus only
  # those feature columns that belong to the active supervised
  # whitelist (or their `_isNA` companions). This keeps OOF training
  # symmetric with `train_final_model_direct()` -- both consume the
  # same set of model features under the same override.
  apply_whitelist_filter <- function(df, keep_class = FALSE) {
    geom_col <- attr(df, "sf_column")
    if (is.null(geom_col) || !nzchar(geom_col)) geom_col <- "geometry"
    protected <- unique(c(id_cols, geom_col))
    if (isTRUE(keep_class)) protected <- unique(c(protected, "class"))
    allowed_features <- c(
      active_whitelist,
      paste0(active_whitelist, "_isNA")
    )
    keep <- unique(c(
      intersect(names(df), protected),
      intersect(names(df), allowed_features)
    ))
    df[, keep, drop = FALSE]
  }
  labelled    <- apply_whitelist_filter(labelled,    keep_class = TRUE)
  burned_like <- apply_whitelist_filter(burned_like, keep_class = FALSE)

  dm <- build_design_matrix_patches(
    labelled    = labelled,
    burned_like = burned_like,
    id_cols     = id_cols,
    drop_regex  = drop_regex,
    cat_cols    = cat_cols,
    hs_n_col    = hs_n_col,
    hs_conf_col = hs_conf_col,
    hs_frp_col  = hs_frp_col,
    median_from = median_from,
    save_dir    = save_dir_dm,
    save_prefix = save_prefix,
    overwrite   = overwrite
  )

  # 0.5.0: build the per-fold feature_weights vector aligned to the
  # design matrix that `build_design_matrix_patches` actually
  # produced. Names not present after the whitelist filter are
  # silently dropped with a warning.
  feature_weights_applied <- list()
  feature_weights_vector <- NULL
  if (!is.null(feature_weights)) {
    active_cols <- colnames(dm$XL_mat)
    unknown_weighted <- setdiff(names(feature_weights), active_cols)
    if (length(unknown_weighted) > 0) {
      warning(
        "Names in `feature_weights` not present in the active OOF feature ",
        "space (silently dropped): ",
        paste(unknown_weighted, collapse = ", "),
        call. = FALSE
      )
    }
    fw_vector <- rep(1.0, length(active_cols))
    names(fw_vector) <- active_cols
    in_both <- intersect(names(feature_weights), active_cols)
    if (length(in_both) > 0L) {
      fw_vector[in_both] <- as.numeric(feature_weights[in_both])
    }
    feature_weights_vector <- fw_vector
    feature_weights_applied <- list(
      requested = as.list(feature_weights),
      applied = as.list(fw_vector[fw_vector != 1.0]),
      unknown_dropped = unknown_weighted
    )
  }

  oof <- run_oof_xgb(
    XL_mat      = dm$XL_mat,
    y           = dm$y,
    labelled_df = labelled_df,
    fold_cols   = fold_cols,
    params      = params,
    nrounds_max = nrounds_max,
    early_stop  = early_stop,
    seed_base   = seed_base,
    out_dir     = out_dir_oof,
    prefix      = prefix,
    verbose     = verbose,
    feature_weights_vector = feature_weights_vector
  )

  oof_agg_path <- file.path(out_dir_oof, paste0(prefix, "_oof_agg.csv"))
  oof_long_path <- file.path(out_dir_oof, paste0(prefix, "_oof_long.csv"))
  labeled_oof_summary_gpkg <- file.path(out_dir_oof, paste0(prefix, "_labeled_oof_summary.gpkg"))
  oof_metrics_path <- file.path(out_dir_oof, paste0(prefix, "_oof_metrics_by_threshold.csv"))
  oof_metrics_summary_path <- file.path(out_dir_oof, paste0(prefix, "_oof_metrics_summary.txt"))
  oof_best_thresholds_path <- file.path(out_dir_oof, paste0(prefix, "_oof_best_thresholds.csv"))

  oof_metrics <- compute_oof_metrics_by_threshold(oof$oof_agg)
  utils::write.csv(oof_metrics, oof_metrics_path, row.names = FALSE)

  # Bug 9 (0.3.0): when an OOF metric column is degenerate (all NA),
  # `which.max()` returns integer(0), which produces a 0-row subset
  # that crashes the downstream data.frame builder. Guard each pick:
  # if the column is all NA, emit a sentinel one-row data.frame with
  # NAs and warn the caller.
  .pick_best <- function(df, col) {
    v <- df[[col]]
    if (is.null(v) || all(is.na(v))) {
      warning(sprintf(paste0("OOF metric '%s' is degenerate (all NA); writing ",
                             "a sentinel row with NA values."), col),
              call. = FALSE)
      sentinel <- df[NA_integer_, , drop = FALSE]
      if (nrow(sentinel) == 0L) {
        # Build sentinel column-set from the input data.frame schema.
        sentinel <- as.data.frame(
          stats::setNames(
            replicate(ncol(df), NA, simplify = FALSE),
            names(df)
          ),
          stringsAsFactors = FALSE
        )
      } else {
        sentinel <- sentinel[1L, , drop = FALSE]
      }
      return(sentinel)
    }
    df[which.max(v), , drop = FALSE]
  }
  best_acc <- .pick_best(oof_metrics, "accuracy")
  best_bal <- .pick_best(oof_metrics, "balanced_accuracy")
  best_f1  <- .pick_best(oof_metrics, "f1")
  recommended <- best_bal

  best_thresholds <- data.frame(
    prefix = prefix,
    created_at = as.character(Sys.time()),
    n_oof_rows = nrow(oof$oof_agg),
    burned_rows = sum(as.character(oof$oof_agg$class) == "burned", na.rm = TRUE),
    unburned_rows = sum(as.character(oof$oof_agg$class) == "unburned", na.rm = TRUE),
    recommended_metric = "balanced_accuracy",
    recommended_threshold = recommended$threshold,
    recommended_accuracy = recommended$accuracy,
    recommended_precision = recommended$precision,
    recommended_recall = recommended$recall,
    recommended_specificity = recommended$specificity,
    recommended_balanced_accuracy = recommended$balanced_accuracy,
    recommended_f1 = recommended$f1,
    best_accuracy_threshold = best_acc$threshold,
    best_accuracy_value = best_acc$accuracy,
    best_accuracy_precision = best_acc$precision,
    best_accuracy_recall = best_acc$recall,
    best_accuracy_specificity = best_acc$specificity,
    best_accuracy_f1 = best_acc$f1,
    best_balanced_accuracy_threshold = best_bal$threshold,
    best_balanced_accuracy_value = best_bal$balanced_accuracy,
    best_balanced_accuracy_precision = best_bal$precision,
    best_balanced_accuracy_recall = best_bal$recall,
    best_balanced_accuracy_specificity = best_bal$specificity,
    best_balanced_accuracy_f1 = best_bal$f1,
    best_f1_threshold = best_f1$threshold,
    best_f1_value = best_f1$f1,
    best_f1_accuracy = best_f1$accuracy,
    best_f1_precision = best_f1$precision,
    best_f1_recall = best_f1$recall,
    best_f1_specificity = best_f1$specificity,
    best_f1_balanced_accuracy = best_f1$balanced_accuracy,
    stringsAsFactors = FALSE
  )
  utils::write.csv(best_thresholds, oof_best_thresholds_path, row.names = FALSE)

  # 0.5.0: record whether feature_whitelist_override / feature_weights
  # were applied at the OOF stage, for parity with the final-model
  # meta and for reproducibility audits.
  if (is.null(feature_whitelist_override)) {
    oof_whitelist_lines <- c(
      "feature_whitelist_override_applied: FALSE",
      paste0("feature_whitelist_n: ", length(.supervised_feature_cols))
    )
  } else {
    dropped_relative <- setdiff(.supervised_feature_cols,
                                 feature_whitelist_override)
    oof_whitelist_lines <- c(
      "feature_whitelist_override_applied: TRUE",
      paste0("feature_whitelist_override_n: ", length(feature_whitelist_override)),
      paste0("feature_whitelist_dropped_n: ", length(dropped_relative)),
      paste0("feature_whitelist_dropped: ",
             paste(dropped_relative, collapse = ","))
    )
  }
  if (is.null(feature_weights)) {
    oof_weights_lines <- "feature_weights_applied: FALSE"
  } else {
    applied_names <- names(feature_weights_applied$applied %||% list())
    applied_pairs <- if (length(applied_names) > 0L) {
      vapply(applied_names, function(nm) {
        paste0(nm, "=", feature_weights_applied$applied[[nm]])
      }, character(1))
    } else {
      character(0)
    }
    oof_weights_lines <- c(
      "feature_weights_applied: TRUE",
      paste0("feature_weights_n_nondefault: ", length(applied_names)),
      paste0("feature_weights_nondefault: ",
             paste(applied_pairs, collapse = ",")),
      paste0("feature_weights_unknown_dropped: ",
             paste(feature_weights_applied$unknown_dropped %||% character(0),
                   collapse = ","))
    )
  }

  summary_lines <- c(
    paste0("prefix: ", prefix),
    paste0("created_at: ", as.character(Sys.time())),
    paste0("n_oof_rows: ", nrow(oof$oof_agg)),
    paste0("burned_rows: ", sum(as.character(oof$oof_agg$class) == "burned", na.rm = TRUE)),
    paste0("unburned_rows: ", sum(as.character(oof$oof_agg$class) == "unburned", na.rm = TRUE)),
    "",
    "[feature_space_overrides]",
    oof_whitelist_lines,
    oof_weights_lines,
    "",
    "[recommended_threshold]",
    paste0("criterion: balanced_accuracy"),
    paste0("threshold: ", signif(recommended$threshold, 6)),
    paste0("accuracy: ", signif(recommended$accuracy, 6)),
    paste0("precision: ", signif(recommended$precision, 6)),
    paste0("recall: ", signif(recommended$recall, 6)),
    paste0("specificity: ", signif(recommended$specificity, 6)),
    paste0("balanced_accuracy: ", signif(recommended$balanced_accuracy, 6)),
    paste0("f1: ", signif(recommended$f1, 6)),
    "",
    "[best_accuracy]",
    paste0("threshold: ", signif(best_acc$threshold, 6)),
    paste0("accuracy: ", signif(best_acc$accuracy, 6)),
    paste0("precision: ", signif(best_acc$precision, 6)),
    paste0("recall: ", signif(best_acc$recall, 6)),
    paste0("specificity: ", signif(best_acc$specificity, 6)),
    paste0("f1: ", signif(best_acc$f1, 6)),
    "",
    "[best_balanced_accuracy]",
    paste0("threshold: ", signif(best_bal$threshold, 6)),
    paste0("balanced_accuracy: ", signif(best_bal$balanced_accuracy, 6)),
    paste0("precision: ", signif(best_bal$precision, 6)),
    paste0("recall: ", signif(best_bal$recall, 6)),
    paste0("specificity: ", signif(best_bal$specificity, 6)),
    paste0("f1: ", signif(best_bal$f1, 6)),
    "",
    "[best_f1]",
    paste0("threshold: ", signif(best_f1$threshold, 6)),
    paste0("f1: ", signif(best_f1$f1, 6)),
    paste0("accuracy: ", signif(best_f1$accuracy, 6)),
    paste0("precision: ", signif(best_f1$precision, 6)),
    paste0("recall: ", signif(best_f1$recall, 6)),
    paste0("specificity: ", signif(best_f1$specificity, 6)),
    paste0("balanced_accuracy: ", signif(best_f1$balanced_accuracy, 6))
  )
  writeLines(summary_lines, con = oof_metrics_summary_path)

  if (!is.null(labelled_gpkg) && file.exists(labelled_gpkg)) {
    L_sf <- sf::read_sf(labelled_gpkg, layer = labelled_layer, quiet = TRUE)
    if ("fire_uid" %in% names(L_sf)) {
      L_sf$fire_uid <- as.character(L_sf$fire_uid)
      oof$oof_agg$fire_uid <- as.character(oof$oof_agg$fire_uid)
      oof_sf <- dplyr::left_join(L_sf, oof$oof_agg, by = c("fire_uid", "class"))
      if (file.exists(labeled_oof_summary_gpkg)) file.remove(labeled_oof_summary_gpkg)
      sf::st_write(oof_sf, labeled_oof_summary_gpkg, layer = "labeled_oof_summary", quiet = TRUE)
    }
  }

  invisible(list(
    dm = dm,
    oof = oof,
    files = list(
      dm_dir = save_dir_dm,
      oof_dir = out_dir_oof,
      oof_agg_path = oof_agg_path,
      oof_long_path = oof_long_path,
      labeled_oof_summary_gpkg = labeled_oof_summary_gpkg,
      oof_metrics_path = oof_metrics_path,
      oof_metrics_summary_path = oof_metrics_summary_path,
      oof_best_thresholds_path = oof_best_thresholds_path
    )
  ))
}
