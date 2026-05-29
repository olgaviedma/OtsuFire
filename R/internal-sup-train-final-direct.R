
train_final_model_direct <- function(
    qa = NULL,
    labelled_gpkg,
    labelled_layer = "train_features",
    id_col = "fire_uid",
    class_col = "class",
    # --- Deterministic-drop negative groups (replaces old hard_negative_source) ---
    # Source field identifying all deterministic-drop polygons.
    deterministic_drop_source = c("deterministic_drop_hard"),
    # neg_types routed to spectral_hard_negative (genuine spectral boundary cases).
    # All other neg_types from deterministic_drop_source go to contextual_exclusion.
    spectral_hard_negative_neg_types = c("spectral_reject_medium"),
    # Contextual-exclusion cap: spectrally-hot geographic exclusions; capped to
    # limit anti-burned signal (they look hotter than fire). Default 0.25x burned.
    contextual_exclusion_to_burned_ratio = 0.25,
    # Spectral-hard-negative cap: genuine decision-boundary cases; use all available.
    spectral_hard_negative_to_burned_ratio = 1.0,
    # --- Random-background cells (low-RBR cell-scale, easy cold) ---
    random_background_source = c("random_burnable_background"),
    random_to_burned_ratio = 1,
    # --- Otsu current-year unburned patches (patch-scale, easy cold; separate for traceability) ---
    # Kept separate from background_cell so counts, caps, and model diagnostics can
    # be evaluated independently. Both are easy-cold negatives but differ in geometry
    # scale (0.81 ha cells vs ~9.72 ha patches) and derivation.
    otsu_unburned_source = c("otsu_patch_residual"),
    otsu_unburned_to_burned_ratio = 1,
    # Exclude otsu_patch_review (ambiguous; S_PATCH_PA 0.15-0.45) and
    # otsu_patch_keep (S_PATCH_PA 0.45-0.70, substantial burned-pixel coverage).
    # Only otsu_patch_drop (S_PATCH_PA <= 0.15, easy cold) is eligible for training.
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep"),
    sampling_seed = 999,
    group_col = "block_id",
    val_frac = 0.15,
    seed = 999,
    params = NULL,
    nrounds_max = 4000,
    early_stopping_rounds = 80,
    # 0.5.0 (2026-05-09): subset of `.supervised_feature_cols` that
    # the model is allowed to see. NULL (default) means use the full
    # canonical whitelist. Names not in `.supervised_feature_cols`
    # are rejected with a clear error -- the canonical list is fixed
    # and this argument can only restrict it, never extend it.
    feature_whitelist_override = NULL,
    # 0.5.0 (2026-05-09): named numeric vector of per-feature weights
    # forwarded to xgboost via `xgb.DMatrix` `feature_weights` info.
    # Only features the user wants to weight differently from 1.0
    # need to be named; anything else implicitly receives 1.0.
    # Replaces Phase B Exp5a's external `xgboost::setinfo()` block.
    feature_weights = NULL,
    impute_numeric = c("median", "zero"),
    impute_factor_missing = "MISSING",
    out_dir = NULL,
    prefix = "2022_patch_certified_v2",
    overwrite = TRUE,
    verbose = TRUE,
    ...
) {
  impute_numeric <- match.arg(impute_numeric)

  # 0.5.0 hard removal: `extra_drop_cols` / `additional_drop_cols`
  # are no longer accepted. Catch them via `...` so old callers get
  # a clear migration message instead of an "unused argument" error.
  .dots <- list(...)
  if ("extra_drop_cols" %in% names(.dots) ||
      "additional_drop_cols" %in% names(.dots)) {
    stop(
      "`extra_drop_cols` / `additional_drop_cols` were removed in ",
      "OtsuFire 0.5.0. Use `feature_whitelist_override` instead. ",
      "Pass a subset of `.supervised_feature_cols` (e.g., the ",
      "current whitelist minus the columns you want to drop). ",
      "See NEWS.md and the migration note in HANDOFF Section N+20.",
      call. = FALSE
    )
  }
  if (length(.dots) > 0L) {
    stop("Unused arguments passed to train_final_model_direct(): ",
         paste(names(.dots), collapse = ", "), call. = FALSE)
  }

  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

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
        paste(unknown, collapse = ", "), ". Allowed names: see ",
        "OtsuFire:::.supervised_feature_cols. The canonical list is ",
        "fixed; this argument can only restrict it, not extend it.",
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

  if (!file.exists(labelled_gpkg)) stop("No existe labelled_gpkg: ", labelled_gpkg)
  L <- sf::read_sf(labelled_gpkg, layer = labelled_layer, quiet = TRUE)
  if (!id_col %in% names(L)) stop("labelled no tiene id_col=", id_col)
  if (!class_col %in% names(L)) stop("labelled no tiene class_col=", class_col)

  L[[id_col]] <- as.character(L[[id_col]])
  L[[class_col]] <- as.character(L[[class_col]])
  deterministic_drop_source        <- unique(as.character(stats::na.omit(deterministic_drop_source)))
  spectral_hard_negative_neg_types <- unique(as.character(stats::na.omit(spectral_hard_negative_neg_types)))
  random_background_source         <- unique(as.character(stats::na.omit(random_background_source)))
  otsu_unburned_source             <- unique(as.character(stats::na.omit(otsu_unburned_source)))
  otsu_unburned_exclude_neg_types  <- unique(as.character(stats::na.omit(otsu_unburned_exclude_neg_types)))

  burned_pool <- L |>
    dplyr::filter(.data[[class_col]] == "burned")

  # spectral_hard_negative: deterministic drops whose neg_type marks them as
  # genuine spectral boundary cases (e.g. spectral_reject_medium).
  spectral_hard_negative_pool <- L |>
    dplyr::filter(
      .data[[class_col]] == "unburned",
      .data[["source"]] %in% deterministic_drop_source,
      .data[["neg_type"]] %in% spectral_hard_negative_neg_types
    )

  # contextual_exclusion: all other deterministic drops -- geographically or
  # data-quality excluded polygons (e.g. geo_excluded_hot). Capped at a low
  # ratio because these are spectrally hotter than burned, not cold.
  contextual_exclusion_pool <- L |>
    dplyr::filter(
      .data[[class_col]] == "unburned",
      .data[["source"]] %in% deterministic_drop_source,
      !(.data[["neg_type"]] %in% spectral_hard_negative_neg_types)
    )

  # random_background_sampled: low-RBR cell-scale polygons (background_cell neg_type).
  # Spectrally cold and unambiguous. No exclusions needed -- this source never
  # produces otsu_patch_keep neg_types.
  random_background_pool <- L |>
    dplyr::filter(
      .data[[class_col]] == "unburned",
      .data[["source"]] %in% random_background_source
    )

  # otsu_unburned_sampled: Otsu current-year patch-scale unburned objects.
  # Kept separate from background_cell for independent counting, capping, and
  # model diagnostics. Both are easy-cold but differ in geometry scale and source.
  # Only otsu_patch_drop (S_PATCH_PA <= 0.15) is eligible for training.
  # otsu_patch_review (S_PATCH_PA 0.15-0.45, ambiguous) excluded by default.
  # otsu_patch_keep (S_PATCH_PA 0.45-0.70) excluded by default.
  # In deterministic_direct mode this pool is always empty (no otsu_patch_residual
  # rows in that pool) -- no behavioral change for Mode A runs.
  otsu_unburned_pool <- L |>
    dplyr::filter(
      .data[[class_col]] == "unburned",
      .data[["source"]] %in% otsu_unburned_source,
      !(.data[["neg_type"]] %in% otsu_unburned_exclude_neg_types)
    )

  n_burned          <- nrow(burned_pool)
  target_contextual <- ceiling(n_burned * contextual_exclusion_to_burned_ratio)
  target_spectral   <- ceiling(n_burned * spectral_hard_negative_to_burned_ratio)
  target_random_bg  <- ceiling(n_burned * random_to_burned_ratio)
  target_otsu       <- ceiling(n_burned * otsu_unburned_to_burned_ratio)

  set.seed(sampling_seed)
  if (target_contextual > 0L && nrow(contextual_exclusion_pool) > 0L) {
    n_take_ctx <- min(target_contextual, nrow(contextual_exclusion_pool))
    ctx_idx    <- sample.int(nrow(contextual_exclusion_pool), n_take_ctx)
    sampled_contextual <- contextual_exclusion_pool[ctx_idx, , drop = FALSE]
  } else {
    sampled_contextual <- contextual_exclusion_pool[0, , drop = FALSE]
  }

  if (target_spectral > 0L && nrow(spectral_hard_negative_pool) > 0L) {
    n_take_shn <- min(target_spectral, nrow(spectral_hard_negative_pool))
    shn_idx    <- sample.int(nrow(spectral_hard_negative_pool), n_take_shn)
    sampled_spectral <- spectral_hard_negative_pool[shn_idx, , drop = FALSE]
  } else {
    sampled_spectral <- spectral_hard_negative_pool[0, , drop = FALSE]
  }

  if (target_random_bg > 0L && nrow(random_background_pool) > 0L) {
    n_take_rb  <- min(target_random_bg, nrow(random_background_pool))
    rb_idx     <- sample.int(nrow(random_background_pool), n_take_rb)
    sampled_random_bg <- random_background_pool[rb_idx, , drop = FALSE]
  } else {
    sampled_random_bg <- random_background_pool[0, , drop = FALSE]
  }

  if (target_otsu > 0L && nrow(otsu_unburned_pool) > 0L) {
    n_take_otsu <- min(target_otsu, nrow(otsu_unburned_pool))
    otsu_idx    <- sample.int(nrow(otsu_unburned_pool), n_take_otsu)
    sampled_otsu <- otsu_unburned_pool[otsu_idx, , drop = FALSE]
  } else {
    sampled_otsu <- otsu_unburned_pool[0, , drop = FALSE]
  }

  if (nrow(burned_pool)) {
    burned_pool <- dplyr::mutate(
      burned_pool,
      training_selected = 1L,
      training_group = "burned_all",
      training_reason = "all_burned_pool"
    )
  }
  if (nrow(sampled_contextual)) {
    sampled_contextual <- dplyr::mutate(
      sampled_contextual,
      training_selected = 1L,
      training_group = "contextual_exclusion",
      training_reason = "sampled_contextual_exclusion"
    )
  }
  if (nrow(sampled_spectral)) {
    sampled_spectral <- dplyr::mutate(
      sampled_spectral,
      training_selected = 1L,
      training_group = "spectral_hard_negative",
      training_reason = "sampled_spectral_hard_negative"
    )
  }
  if (nrow(sampled_random_bg)) {
    sampled_random_bg <- dplyr::mutate(
      sampled_random_bg,
      training_selected = 1L,
      training_group = "random_background_sampled",
      training_reason = "sampled_random_background"
    )
  }
  if (nrow(sampled_otsu)) {
    sampled_otsu <- dplyr::mutate(
      sampled_otsu,
      training_selected = 1L,
      training_group = "otsu_unburned_sampled",
      training_reason = "sampled_otsu_unburned"
    )
  }

  L_ok <- dplyr::bind_rows(
    burned_pool, sampled_contextual, sampled_spectral, sampled_random_bg, sampled_otsu
  ) |>
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)

  msg("Direct training selection:")
  msg("  burned selected=%d", nrow(burned_pool))
  msg("  contextual_exclusion selected=%d (cap=%.2fx, available=%d)",
      nrow(sampled_contextual), contextual_exclusion_to_burned_ratio, nrow(contextual_exclusion_pool))
  msg("  spectral_hard_negative selected=%d (cap=%.2fx, available=%d)",
      nrow(sampled_spectral), spectral_hard_negative_to_burned_ratio, nrow(spectral_hard_negative_pool))
  msg("  random_background_sampled selected=%d (cap=%.2fx, available=%d)",
      nrow(sampled_random_bg), random_to_burned_ratio, nrow(random_background_pool))
  msg("  otsu_unburned_sampled selected=%d (cap=%.2fx, available=%d)",
      nrow(sampled_otsu), otsu_unburned_to_burned_ratio, nrow(otsu_unburned_pool))
  msg("  training_ok total=%d", nrow(L_ok))

  if (nrow(L_ok) < 20) stop("La seleccion directa dejo muy pocos rows para entrenar.")

  y <- as.integer(L_ok[[class_col]] == "burned")
  L_df <- sf::st_drop_geometry(L_ok)
  geom_col <- attr(L_ok, "sf_column") %||% "geometry"

  # 0.4.0 architectural refactor (Agent H): the supervised model
  # sees ONLY the columns in `.supervised_feature_cols` (plus their
  # `_isNA` companions when present). 0.5.0: when the caller passes
  # `feature_whitelist_override`, the active whitelist becomes that
  # subset (validated above), narrowing the canonical list. Filter
  # the polygon-level frame to the active whitelist intersected with
  # what is actually available.
  feat_cols <- .filter_to_supervised_whitelist(names(L_df),
                                                whitelist = active_whitelist)

  is_listcol <- vapply(L_df[, feat_cols, drop = FALSE], is.list, logical(1))
  if (any(is_listcol)) feat_cols <- feat_cols[!is_listcol]
  if (length(feat_cols) == 0) {
    stop("No supervised feature columns survived whitelist filter. ",
         "Expected at least some of the active whitelist in the ",
         "input GPKG layer '", labelled_layer, "'.")
  }

  # Whitelist assertion: every kept column must be either a member of
  # the active whitelist or a recognised `_isNA` companion. If
  # anything else slipped through, abort loud -- that means
  # `.filter_to_supervised_whitelist()` and the constant fell out of
  # sync.
  .allowed_supervised_cols <- c(
    active_whitelist,
    paste0(active_whitelist, "_isNA")
  )
  forbidden <- setdiff(feat_cols, .allowed_supervised_cols)
  if (length(forbidden) > 0L) {
    stop("Whitelist invariant broken in train_final_model_direct(): ",
         "feat_cols contains columns outside the active whitelist: ",
         paste(forbidden, collapse = ", "))
  }

  X_df <- L_df[, feat_cols, drop = FALSE]
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

  numeric_medians <- list()
  for (nm in names(X_df)) {
    if (is.numeric(X_df[[nm]]) || is.integer(X_df[[nm]])) {
      v <- X_df[[nm]]
      bad_nonfinite <- !is.finite(v)
      if (any(bad_nonfinite, na.rm = TRUE)) {
        v[bad_nonfinite] <- NA
      }
      if (anyNA(v)) {
        fill <- if (impute_numeric == "zero") 0 else stats::median(v, na.rm = TRUE)
        if (!is.finite(fill)) fill <- 0
        v[is.na(v)] <- fill
        X_df[[nm]] <- v
        numeric_medians[[nm]] <- fill
      }
    }
  }

  row_id <- seq_len(nrow(X_df))
  rownames(X_df) <- as.character(row_id)
  X <- Matrix::sparse.model.matrix(~ . - 1, data = X_df)
  if (nrow(X) != length(row_id)) {
    kept <- rownames(X)
    kept <- suppressWarnings(as.integer(kept))
    kept <- kept[!is.na(kept)]
    msg("WARNING: sparse.model.matrix devolvio %d filas (esperaba %d). Realineando...", nrow(X), length(row_id))
    L_ok <- L_ok[kept, ]
    y <- y[kept]
  }

  set.seed(seed)
  n <- nrow(X)
  if (!is.null(group_col) && group_col %in% names(L_ok)) {
    g <- as.character(L_ok[[group_col]])
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
      split_mode <- "row_split_fallback"
    } else {
      split_mode <- "group_split"
    }
  } else {
    idx <- sample.int(n)
    nval <- max(1, floor(val_frac * n))
    val_idx <- idx[1:nval]
    tr_idx <- idx[(nval + 1):n]
    split_mode <- "row_split"
  }

  dtrain <- xgboost::xgb.DMatrix(X[tr_idx, , drop = FALSE], label = y[tr_idx])
  dval <- xgboost::xgb.DMatrix(X[val_idx, , drop = FALSE], label = y[val_idx])

  # 0.5.0: apply per-feature weights via xgboost's `feature_weights`
  # info vector. Names not in the active feature space (after the
  # whitelist filter and override) are silently dropped with a
  # warning. Features the user did not name implicitly receive 1.0.
  feature_weights_applied <- list()
  if (!is.null(feature_weights)) {
    active_cols <- colnames(X)
    unknown_weighted <- setdiff(names(feature_weights), active_cols)
    if (length(unknown_weighted) > 0) {
      warning(
        "Names in `feature_weights` not present in the active feature ",
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
    xgboost::setinfo(dtrain, "feature_weights", fw_vector)
    xgboost::setinfo(dval,   "feature_weights", fw_vector)
    feature_weights_applied <- list(
      requested = as.list(feature_weights),
      applied = as.list(fw_vector[fw_vector != 1.0]),
      unknown_dropped = unknown_weighted
    )
  }

  if (is.null(params)) {
    spw <- sum(y[tr_idx] == 0) / max(1, sum(y[tr_idx] == 1))
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = "logloss",
      eta = 0.05,
      max_depth = 6,
      min_child_weight = 1,
      subsample = 0.8,
      colsample_bytree = 0.8,
      lambda = 1,
      alpha = 0,
      scale_pos_weight = spw
    )
  }

  set.seed(seed)
  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds_max,
    watchlist = list(train = dtrain, val = dval),
    early_stopping_rounds = early_stopping_rounds,
    verbose = if (isTRUE(verbose)) 1 else 0
  )

  files <- list()
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    gpkg_ok <- file.path(out_dir, paste0(prefix, "_training_ok.gpkg"))
    csv_ok <- file.path(out_dir, paste0(prefix, "_training_ok.csv"))
    rds_mod <- file.path(out_dir, paste0(prefix, "_final_model.rds"))
    rds_rec <- file.path(out_dir, paste0(prefix, "_recipe.rds"))
    rds_spl <- file.path(out_dir, paste0(prefix, "_split_idx.rds"))
    txt_meta <- file.path(out_dir, paste0(prefix, "_meta.txt"))
    csv_imp <- file.path(out_dir, paste0(prefix, "_feature_importance.csv"))
    txt_summary <- file.path(out_dir, paste0(prefix, "_model_summary.txt"))

    if (isTRUE(overwrite) && file.exists(gpkg_ok)) file.remove(gpkg_ok)
    sf::st_write(L_ok, gpkg_ok, layer = "training_ok", quiet = TRUE)
    utils::write.csv(sf::st_drop_geometry(L_ok), csv_ok, row.names = FALSE)
    saveRDS(model, rds_mod)
    saveRDS(list(train_idx = tr_idx, val_idx = val_idx, mode = split_mode), rds_spl)

    recipe <- list(
      created_at = as.character(Sys.time()),
      inputs = list(labelled_gpkg = labelled_gpkg, labelled_layer = labelled_layer),
      selection = list(
        selection_mode = "direct_pool_sampling",
        deterministic_drop_source = deterministic_drop_source,
        spectral_hard_negative_neg_types = spectral_hard_negative_neg_types,
        contextual_exclusion_to_burned_ratio = contextual_exclusion_to_burned_ratio,
        spectral_hard_negative_to_burned_ratio = spectral_hard_negative_to_burned_ratio,
        random_background_source = random_background_source,
        random_to_burned_ratio = random_to_burned_ratio,
        otsu_unburned_source = otsu_unburned_source,
        otsu_unburned_to_burned_ratio = otsu_unburned_to_burned_ratio,
        otsu_unburned_exclude_neg_types = otsu_unburned_exclude_neg_types,
        burned_selected = nrow(burned_pool),
        contextual_exclusion_available = nrow(contextual_exclusion_pool),
        contextual_exclusion_selected = nrow(sampled_contextual),
        spectral_hard_negative_available = nrow(spectral_hard_negative_pool),
        spectral_hard_negative_selected = nrow(sampled_spectral),
        random_background_available = nrow(random_background_pool),
        random_background_selected = nrow(sampled_random_bg),
        otsu_unburned_available = nrow(otsu_unburned_pool),
        otsu_unburned_selected = nrow(sampled_otsu)
      ),
      impute = list(
        impute_numeric = impute_numeric,
        numeric_medians = numeric_medians,
        impute_factor_missing = impute_factor_missing
      ),
      cols = list(
        id_col = id_col,
        class_col = class_col,
        group_col = group_col,
        feature_cols = feat_cols,
        x_cols = colnames(X)
      ),
      training = list(
        val_frac = val_frac,
        seed = seed,
        nrounds_max = nrounds_max,
        early_stopping_rounds = early_stopping_rounds,
        best_iteration = model$best_iteration %||% NA_integer_
      ),
      params = params,
      # 0.5.0: record what was applied so reproducibility audits can
      # tell whether a non-default feature space or weight vector was
      # used for this final-model fit.
      feature_whitelist_override = feature_whitelist_override,
      feature_weights = feature_weights_applied
    )
    saveRDS(recipe, rds_rec)

    imp_tbl <- tryCatch(
      xgboost::xgb.importance(model = model, feature_names = recipe$cols$x_cols),
      error = function(e) NULL
    )
    if (is.data.frame(imp_tbl) && nrow(imp_tbl)) {
      utils::write.csv(imp_tbl, csv_imp, row.names = FALSE)
    }

    oof_summary_path <- NULL
    oof_best_thresholds_path <- NULL
    oof_best <- NULL
    if (!is.null(qa) && is.character(qa) && length(qa) == 1L && file.exists(qa)) {
      oof_summary_path <- sub("_oof_agg\\.csv$", "_oof_metrics_summary.txt", qa)
      oof_best_thresholds_path <- sub("_oof_agg\\.csv$", "_oof_best_thresholds.csv", qa)
      if (!identical(oof_best_thresholds_path, qa) && file.exists(oof_best_thresholds_path)) {
        oof_best <- tryCatch(utils::read.csv(oof_best_thresholds_path, stringsAsFactors = FALSE), error = function(e) NULL)
      }
    }

    # 0.5.0: persist whether `feature_whitelist_override` was applied
    # (and which canonical names were dropped relative to the
    # canonical whitelist) and whether `feature_weights` was applied
    # (with the names whose weight differs from 1.0). Required for
    # reproducibility audits.
    if (is.null(feature_whitelist_override)) {
      whitelist_override_lines <- c(
        "feature_whitelist_override_applied: FALSE",
        paste0("feature_whitelist_n: ", length(.supervised_feature_cols))
      )
    } else {
      dropped_relative <- setdiff(.supervised_feature_cols,
                                   feature_whitelist_override)
      whitelist_override_lines <- c(
        "feature_whitelist_override_applied: TRUE",
        paste0("feature_whitelist_override_n: ", length(feature_whitelist_override)),
        paste0("feature_whitelist_dropped_n: ", length(dropped_relative)),
        paste0("feature_whitelist_dropped: ",
               paste(dropped_relative, collapse = ","))
      )
    }
    if (is.null(feature_weights)) {
      feature_weights_lines <- c(
        "feature_weights_applied: FALSE"
      )
    } else {
      applied_names <- names(feature_weights_applied$applied %||% list())
      applied_pairs <- if (length(applied_names) > 0L) {
        vapply(applied_names, function(nm) {
          paste0(nm, "=", feature_weights_applied$applied[[nm]])
        }, character(1))
      } else {
        character(0)
      }
      feature_weights_lines <- c(
        "feature_weights_applied: TRUE",
        paste0("feature_weights_n_nondefault: ", length(applied_names)),
        paste0("feature_weights_nondefault: ",
               paste(applied_pairs, collapse = ",")),
        paste0("feature_weights_unknown_dropped: ",
               paste(feature_weights_applied$unknown_dropped %||% character(0),
                     collapse = ","))
      )
    }

    writeLines(c(
      paste0("prefix: ", prefix),
      paste0("created_at: ", recipe$created_at),
      paste0("labelled_gpkg: ", labelled_gpkg),
      paste0("labelled_layer: ", labelled_layer),
      paste0("selection_mode: direct_pool_sampling"),
      paste0("burned_selected: ", nrow(burned_pool)),
      paste0("contextual_exclusion_available: ", nrow(contextual_exclusion_pool)),
      paste0("contextual_exclusion_selected: ", nrow(sampled_contextual)),
      paste0("contextual_exclusion_to_burned_ratio: ", contextual_exclusion_to_burned_ratio),
      paste0("spectral_hard_negative_available: ", nrow(spectral_hard_negative_pool)),
      paste0("spectral_hard_negative_selected: ", nrow(sampled_spectral)),
      paste0("spectral_hard_negative_to_burned_ratio: ", spectral_hard_negative_to_burned_ratio),
      paste0("random_background_available: ", nrow(random_background_pool)),
      paste0("random_background_selected: ", nrow(sampled_random_bg)),
      paste0("random_to_burned_ratio: ", random_to_burned_ratio),
      paste0("otsu_unburned_available: ", nrow(otsu_unburned_pool)),
      paste0("otsu_unburned_selected: ", nrow(sampled_otsu)),
      paste0("otsu_unburned_to_burned_ratio: ", otsu_unburned_to_burned_ratio),
      paste0("split_mode: ", split_mode),
      paste0("n_training_ok_used: ", nrow(L_ok)),
      paste0("n_feature_cols: ", length(feat_cols)),
      paste0("n_x_cols: ", ncol(X)),
      paste0("best_iteration: ", recipe$training$best_iteration),
      whitelist_override_lines,
      feature_weights_lines
    ), txt_meta)

    summary_lines <- c(
      paste0("prefix: ", prefix),
      paste0("created_at: ", recipe$created_at),
      paste0("labelled_gpkg: ", labelled_gpkg),
      paste0("labelled_layer: ", labelled_layer),
      "",
      "[training_selection]",
      paste0("burned_selected: ", nrow(burned_pool)),
      paste0("contextual_exclusion_available: ", nrow(contextual_exclusion_pool)),
      paste0("contextual_exclusion_selected: ", nrow(sampled_contextual)),
      paste0("contextual_exclusion_to_burned_ratio: ", contextual_exclusion_to_burned_ratio),
      paste0("spectral_hard_negative_available: ", nrow(spectral_hard_negative_pool)),
      paste0("spectral_hard_negative_selected: ", nrow(sampled_spectral)),
      paste0("spectral_hard_negative_to_burned_ratio: ", spectral_hard_negative_to_burned_ratio),
      paste0("random_background_available: ", nrow(random_background_pool)),
      paste0("random_background_selected: ", nrow(sampled_random_bg)),
      paste0("random_to_burned_ratio: ", random_to_burned_ratio),
      paste0("otsu_unburned_available: ", nrow(otsu_unburned_pool)),
      paste0("otsu_unburned_selected: ", nrow(sampled_otsu)),
      paste0("otsu_unburned_to_burned_ratio: ", otsu_unburned_to_burned_ratio),
      paste0("n_training_ok_used: ", nrow(L_ok)),
      paste0("split_mode: ", split_mode),
      paste0("best_iteration: ", recipe$training$best_iteration),
      paste0("n_feature_cols: ", length(feat_cols)),
      paste0("n_x_cols: ", ncol(X)),
      "",
      "[feature_space_overrides]",
      whitelist_override_lines,
      feature_weights_lines
    )

    if (is.data.frame(oof_best) && nrow(oof_best)) {
      summary_lines <- c(
        summary_lines,
        "",
        "[oof_recommended_threshold]",
        paste0("metric: ", oof_best$recommended_metric[1]),
        paste0("threshold: ", oof_best$recommended_threshold[1]),
        paste0("accuracy: ", oof_best$recommended_accuracy[1]),
        paste0("precision: ", oof_best$recommended_precision[1]),
        paste0("recall: ", oof_best$recommended_recall[1]),
        paste0("specificity: ", oof_best$recommended_specificity[1]),
        paste0("balanced_accuracy: ", oof_best$recommended_balanced_accuracy[1]),
        paste0("f1: ", oof_best$recommended_f1[1])
      )
    }

    if (is.data.frame(imp_tbl) && nrow(imp_tbl)) {
      top_imp <- utils::head(imp_tbl, 15)
      summary_lines <- c(
        summary_lines,
        "",
        "[top_feature_importance]",
        vapply(seq_len(nrow(top_imp)), function(i) {
          paste0(
            top_imp$Feature[i],
            " | Gain=", signif(top_imp$Gain[i], 6),
            " | Cover=", signif(top_imp$Cover[i], 6),
            " | Frequency=", signif(top_imp$Frequency[i], 6)
          )
        }, character(1))
      )
    }

    summary_lines <- c(
      summary_lines,
      "",
      "[files]",
      paste0("meta_txt: ", txt_meta),
      paste0("feature_importance_csv: ", csv_imp),
      paste0("oof_summary_txt: ", oof_summary_path %||% NA_character_),
      paste0("oof_best_thresholds_csv: ", oof_best_thresholds_path %||% NA_character_)
    )
    writeLines(summary_lines, txt_summary)

    files <- list(
      training_ok_gpkg = gpkg_ok,
      training_ok_csv = csv_ok,
      model_rds = rds_mod,
      recipe_rds = rds_rec,
      split_rds = rds_spl,
      meta_txt = txt_meta,
      feature_importance_csv = csv_imp,
      model_summary_txt = txt_summary
    )
  }

  invisible(list(
    model = model,
    training_ok_sf = L_ok,
    feature_cols = feat_cols,
    x_cols = colnames(X),
    split = list(train_idx = tr_idx, val_idx = val_idx, mode = split_mode),
    params = params,
    files = files
  ))
}

train_final_model_from_qa <- train_final_model_direct
