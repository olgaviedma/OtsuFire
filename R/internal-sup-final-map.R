
score_burnedlike_and_export_final_map <- function(
    result_dir,
    prefix = "2022_patch_certified_v2",
    year_tag = "2022",
    qa_labelled_gpkg  = file.path(result_dir, "06_QA_LABELS", "2022_patch_labeled_oof_summary.gpkg"),
    qa_labelled_layer = "labeled_oof_summary",
    labelled_features_gpkg  = file.path(result_dir, "03_FEATURES", "features_geometry.gpkg"),
    labelled_features_layer = "train_features",
    model_rds  = file.path(result_dir, "07_FINAL_MODEL_V2", paste0(prefix, "_final_model.rds")),
    recipe_rds = file.path(result_dir, "07_FINAL_MODEL_V2", paste0(prefix, "_recipe.rds")),
    unlabeled_gpkg  = file.path(result_dir, "03_FEATURES", "features_geometry.gpkg"),
    unlabeled_layer = "scoring_features",
    id_col = "fire_uid",
    out_score_dir = file.path(result_dir, "08_SCORED"),
    out_map_dir   = file.path(result_dir, "09_FINAL_MAP"),
    export_burned_like = TRUE,
    burnedlike_basename = paste0(prefix, "_burned_like_scored.gpkg"),
    burnedlike_layer    = "burned_like_scored",
    preyear_overlap_threshold = 0.70,
    hotspot_density_threshold = 0.001,
    temporal_penalty_floor = 0.10,
    overwrite = TRUE,
    verbose = TRUE
) {
  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  dir.create(out_score_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_map_dir, recursive = TRUE, showWarnings = FALSE)

  normalize_geom_name <- function(x, target = "geom") {
    g <- attr(x, "sf_column")
    if (!is.null(g) && g != target && g %in% names(x)) {
      names(x)[names(x) == g] <- target
      attr(x, "sf_column") <- target
    }
    x
  }

  # Bug 2 + AS04 (0.3.0): scoring must mirror the training preprocessing.
  # When `recipe` is supplied, numeric NAs are imputed to the medians
  # recorded at training time (`recipe$impute$numeric_medians`); factor
  # NAs are imputed to `recipe$impute$impute_factor_missing` (default
  # "MISSING"). Columns whose training pass had no NAs simply have no
  # entry in numeric_medians; their NAs are left in place and become
  # implicit-missing under the sparse encoder.
  prep_X_df <- function(df, recipe = NULL) {
    numeric_medians <- if (!is.null(recipe) &&
                           !is.null(recipe$impute$numeric_medians)) {
      recipe$impute$numeric_medians
    } else {
      list()
    }
    impute_factor_missing <- if (!is.null(recipe) &&
                                 !is.null(recipe$impute$impute_factor_missing)) {
      recipe$impute$impute_factor_missing
    } else {
      "MISSING"
    }

    for (nm in names(df)) {
      if (is.factor(df[[nm]])) df[[nm]] <- as.character(df[[nm]])
      if (is.character(df[[nm]])) {
        df[[nm]][is.na(df[[nm]])] <- impute_factor_missing
        lv <- unique(df[[nm]])
        df[[nm]] <- factor(df[[nm]], levels = unique(c(lv, "__OTHER__")))
      }
      if (is.factor(df[[nm]]) && nlevels(df[[nm]]) < 2) {
        lv <- levels(df[[nm]])
        levels(df[[nm]]) <- unique(c(lv, "__OTHER__"))
      }
      if (is.logical(df[[nm]])) df[[nm]] <- as.integer(df[[nm]])
      if (is.numeric(df[[nm]])) {
        v <- df[[nm]]
        v[is.infinite(v)] <- NA_real_
        if (nm %in% names(numeric_medians) && anyNA(v)) {
          v[is.na(v)] <- numeric_medians[[nm]]
        }
        df[[nm]] <- v
      }
    }
    df
  }

  # KB2 (0.3.0): sparse-aware version of the column aligner. Replaces
  # `align_to_x_cols_dense`, which padded missing columns with explicit
  # dense zeros that XGBoost (trained sparse-style) interpreted as real
  # values rather than as missing. Sparse all-zero columns leave the
  # cells unstored, which XGBoost treats as missing — the convention
  # the model was trained under.
  align_to_x_cols <- function(X_new, x_cols) {
    if (!inherits(X_new, "Matrix")) {
      stop("Internal error: align_to_x_cols expects a sparse Matrix; ",
           "got ", paste(class(X_new), collapse = "/"), ".",
           call. = FALSE)
    }
    extra <- setdiff(colnames(X_new), x_cols)
    if (length(extra) > 0) {
      X_new <- X_new[, setdiff(colnames(X_new), extra), drop = FALSE]
    }
    missing_cols <- setdiff(x_cols, colnames(X_new))
    if (length(missing_cols) > 0) {
      Z <- Matrix::Matrix(0, nrow = nrow(X_new),
                          ncol = length(missing_cols),
                          sparse = TRUE)
      colnames(Z) <- missing_cols
      rownames(Z) <- rownames(X_new)
      X_new <- cbind(X_new, Z)
    }
    X_new[, x_cols, drop = FALSE]
  }

  harmonize_cols_for_rbind <- function(A, B) {
    A <- normalize_geom_name(A, "geom")
    B <- normalize_geom_name(B, "geom")

    all_cols <- union(setdiff(names(A), "geom"), setdiff(names(B), "geom"))

    fix_one <- function(x) {
      x <- normalize_geom_name(x, "geom")
      geom <- sf::st_geometry(x)
      x_df <- as.data.frame(sf::st_drop_geometry(x))
      miss <- setdiff(all_cols, names(x_df))
      if (length(miss) > 0) {
        for (m in miss) x_df[[m]] <- NA
      }
      x_df <- x_df[, all_cols, drop = FALSE]
      sf::st_sf(x_df, geom = geom)
    }

    list(A = fix_one(A), B = fix_one(B))
  }

  drop_obsolete_fields <- function(x) {
    obsolete_patterns <- c(
      "^qa_",
      "^p_type$",
      "^score_source$",
      "^p_oof_min_long$",
      "^p_oof_max_long$",
      "^n_oof_long$",
      "^class_oof$",
      "^class_mismatch$"
    )
    obsolete_cols <- names(x)[grepl(paste(obsolete_patterns, collapse = "|"), names(x))]
    keep <- setdiff(names(x), obsolete_cols)
    x[, keep, drop = FALSE]
  }

  make_public_final_map <- function(x) {
    class_input_vec <- if ("class_final" %in% names(x)) {
      as.character(x[["class_final"]])
    } else {
      as.character(x[["class"]])
    }
    x <- dplyr::mutate(x, class_input = class_input_vec)

    public_cols <- c(
      "fire_uid",
      "poly_id",
      "source_poly_id",
      "year",
      "scenario",
      "source_set",
      "source",
      "neg_type",
      "class_input",
      "filter_1",
      "reason_1",
      "filter_2",
      "reason_2",
      "filter_3",
      "reason_3",
      "preyear_action",
      "preyear_reason",
      "preyear_overlap_frac",
      "class_final",
      "area_ha",
      "n_pix",
      "median_rbr",
      "percentile_in_keep",
      "p_above_keep_q25",
      "p_above_keep_ref",
      "hs_used_n",
      "hs_density_ha",
      "hs_frp_max",
      "hs_conf_mean",
      "p_burned",
      "p_burned_model",
      "p_burned_current_year",
      "temporal_penalty",
      "temporal_conflict_flag",
      "p_burned_oof",
      "geom"
    )
    public_cols <- public_cols[public_cols %in% names(x)]
    x[, public_cols, drop = FALSE]
  }

  add_temporal_adjustment <- function(x) {
    x_df <- as.data.frame(sf::st_drop_geometry(x))
    n <- nrow(x_df)

    as_num_local <- function(v) {
      suppressWarnings(as.numeric(gsub(",", ".", as.character(v), fixed = TRUE)))
    }

    area_ha_num <- if ("area_ha" %in% names(x_df)) as_num_local(x_df$area_ha) else rep(NA_real_, n)
    hs_used_num <- if ("hs_used_n" %in% names(x_df)) as_num_local(x_df$hs_used_n) else rep(NA_real_, n)
    overlap_num <- if ("preyear_overlap_frac" %in% names(x_df)) as_num_local(x_df$preyear_overlap_frac) else rep(NA_real_, n)
    p_burned_model <- if ("p_burned" %in% names(x_df)) as_num_local(x_df$p_burned) else rep(NA_real_, n)
    preyear_action_chr <- if ("preyear_action" %in% names(x_df)) as.character(x_df$preyear_action) else rep(NA_character_, n)

    hs_density_ha <- rep(NA_real_, n)
    ok_area <- is.finite(area_ha_num) & area_ha_num > 0 & is.finite(hs_used_num)
    hs_density_ha[ok_area] <- hs_used_num[ok_area] / area_ha_num[ok_area]

    hotspot_available_now <- rep(FALSE, n)
    if ("hotspot_available" %in% names(x_df)) {
      hotspot_available_now <- x_df$hotspot_available %in% c(1, TRUE)
    } else if ("filter_2" %in% names(x_df)) {
      hotspot_available_now <- !is.na(x_df$filter_2) & as.character(x_df$filter_2) != "not_applied"
    }
    hotspot_available_now[is.na(hotspot_available_now)] <- FALSE

    weak_hotspot_density <- rep(NA, n)
    weak_hotspot_density[hotspot_available_now] <-
      !is.finite(hs_density_ha[hotspot_available_now]) |
      hs_density_ha[hotspot_available_now] < hotspot_density_threshold

    preyear_drop_flag <- !is.na(preyear_action_chr) & preyear_action_chr == "drop"
    overlap_proxy <- overlap_num
    overlap_proxy[preyear_drop_flag & !is.finite(overlap_proxy)] <- 1

    temporal_conflict <- (is.finite(overlap_proxy) & overlap_proxy >= preyear_overlap_threshold) | preyear_drop_flag
    weak_current_support <- temporal_conflict & (!hotspot_available_now | weak_hotspot_density %in% TRUE)

    overlap_scaled <- rep(0, n)
    overlap_scaled[temporal_conflict] <-
      pmax(0, pmin(1, (overlap_proxy[temporal_conflict] - preyear_overlap_threshold) /
        max(1e-9, 1 - preyear_overlap_threshold)))

    support_ratio <- rep(1, n)
    if (any(hotspot_available_now)) {
      support_ratio[hotspot_available_now] <-
        pmax(0, pmin(1, hs_density_ha[hotspot_available_now] / hotspot_density_threshold))
      support_ratio[hotspot_available_now & !is.finite(support_ratio)] <- 0
    }
    support_ratio[!hotspot_available_now] <- 0

    temporal_penalty <- rep(1, n)
    temporal_penalty[weak_current_support] <-
      pmax(
        temporal_penalty_floor,
        1 - (1 - temporal_penalty_floor) *
          overlap_scaled[weak_current_support] *
          (1 - support_ratio[weak_current_support])
      )

    p_burned_current_year <- p_burned_model * temporal_penalty

    x$hs_density_ha <- hs_density_ha
    x$temporal_conflict_flag <- temporal_conflict
    x$temporal_penalty <- temporal_penalty
    x$p_burned_model <- p_burned_model
    x$p_burned_current_year <- p_burned_current_year
    x$p_burned <- p_burned_current_year
    x$current_year_public_drop <- weak_current_support

    x
  }

  filter_to_current_year_map <- function(x) {
    if (!"current_year_public_drop" %in% names(x)) return(x)
    keep_idx <- !(x$current_year_public_drop %in% TRUE)
    keep_idx[is.na(keep_idx)] <- TRUE
    x[keep_idx, , drop = FALSE]
  }

  score_with_final_model <- function(sf_obj, model, recipe, id_col) {
    sf_obj <- normalize_geom_name(sf_obj, "geom")
    x_df <- as.data.frame(sf::st_drop_geometry(sf_obj))
    feature_cols <- recipe$cols$feature_cols
    x_cols_train <- recipe$cols$x_cols

    missing_feat <- setdiff(feature_cols, names(x_df))
    if (length(missing_feat) > 0) {
      for (nm in missing_feat) x_df[[nm]] <- NA
      msg("WARNING: se crean %d feature cols ausentes para scoring.", length(missing_feat))
    }

    # KB2 + Bug 2 + AS04 (0.3.0): apply training-time numeric medians
    # BEFORE building the matrix, then use the sparse encoder so the
    # implicit-zero-as-missing convention matches training.
    model_df <- prep_X_df(x_df[, feature_cols, drop = FALSE], recipe = recipe)
    row_id <- seq_len(nrow(model_df))
    model_df$row_id__ <- row_id
    rownames(model_df) <- as.character(row_id)

    form <- stats::as.formula("~ . - 1")
    mf <- stats::model.frame(form, data = model_df,
                              na.action = stats::na.pass)
    X_new <- Matrix::sparse.model.matrix(form, data = mf)

    drop_cols <- grep("^row_id__", colnames(X_new), value = TRUE)
    if (length(drop_cols) > 0) {
      X_new <- X_new[, setdiff(colnames(X_new), drop_cols), drop = FALSE]
    }

    kept_rows <- suppressWarnings(as.integer(rownames(X_new)))
    kept_rows <- kept_rows[!is.na(kept_rows)]
    X_new <- align_to_x_cols(X_new, x_cols_train)

    dmat <- xgboost::xgb.DMatrix(X_new)
    p_kept <- predict(model, dmat)

    p_full <- rep(NA_real_, nrow(sf_obj))
    p_full[kept_rows] <- as.numeric(p_kept)

    sf_obj |>
      dplyr::mutate(
        p_burned = p_full
      )
  }

  if (!file.exists(model_rds)) stop("No existe model_rds: ", model_rds)
  if (!file.exists(recipe_rds)) stop("No existe recipe_rds: ", recipe_rds)
  if (!file.exists(labelled_features_gpkg)) stop("No existe labelled_features_gpkg: ", labelled_features_gpkg)
  if (!file.exists(qa_labelled_gpkg)) stop("No existe qa_labelled_gpkg: ", qa_labelled_gpkg)
  if (!file.exists(unlabeled_gpkg)) stop("No existe unlabeled_gpkg: ", unlabeled_gpkg)

  model <- readRDS(model_rds)
  recipe <- readRDS(recipe_rds)
  S <- sf::read_sf(unlabeled_gpkg, unlabeled_layer, quiet = TRUE)
  S <- normalize_geom_name(S, "geom")
  if (!id_col %in% names(S)) stop("Scoring layer NO tiene id_col=", id_col)
  S[[id_col]] <- as.character(S[[id_col]])

  L_feat <- sf::read_sf(labelled_features_gpkg, labelled_features_layer, quiet = TRUE)
  L_feat <- normalize_geom_name(L_feat, "geom")
  if (!id_col %in% names(L_feat)) stop("train_features NO tiene id_col=", id_col)
  L_feat[[id_col]] <- as.character(L_feat[[id_col]])

  L_qa <- sf::read_sf(qa_labelled_gpkg, qa_labelled_layer, quiet = TRUE) |>
    sf::st_drop_geometry() |>
    as.data.frame()
  if (!id_col %in% names(L_qa)) stop("QA labeled NO tiene id_col=", id_col)
  L_qa[[id_col]] <- as.character(L_qa[[id_col]])

  qa_cols_candidate <- c(id_col, "p_oof_mean", "p_oof_med", "p_oof_sd", "n_preds")
  qa_keep_cols <- intersect(qa_cols_candidate, names(L_qa))
  if (!"p_oof_mean" %in% qa_keep_cols) stop("QA labeled no tiene p_oof_mean (necesario).")
  L_qa2 <- L_qa[, qa_keep_cols, drop = FALSE]
  L_oof <- dplyr::left_join(L_feat, L_qa2, by = id_col) |>
    sf::st_drop_geometry() |>
    as.data.frame()

  oof_join_cols <- intersect(
    c("source_poly_id", "p_oof_mean", "p_oof_med", "p_oof_sd", "n_preds"),
    names(L_oof)
  )
  oof_by_source <- L_oof[, oof_join_cols, drop = FALSE]
  if ("source_poly_id" %in% names(oof_by_source)) {
    oof_by_source$source_poly_id <- as.character(oof_by_source$source_poly_id)
    oof_by_source <- oof_by_source[!duplicated(oof_by_source$source_poly_id), , drop = FALSE]
  }

  if ("source_poly_id" %in% names(S) && "source_poly_id" %in% names(oof_by_source)) {
    S$source_poly_id <- as.character(S$source_poly_id)
    S <- dplyr::left_join(S, oof_by_source, by = "source_poly_id")
  } else {
    S$p_oof_mean <- NA_real_
    S$p_oof_med <- NA_real_
    S$p_oof_sd <- NA_real_
    S$n_preds <- NA_integer_
  }

  final_map_full <- score_with_final_model(S, model, recipe, id_col) |>
    dplyr::mutate(
      p_burned_oof = as.numeric(.data[["p_oof_mean"]]),
      source_set = "deterministic"
    ) |>
    drop_obsolete_fields() |>
    add_temporal_adjustment()
  final_map_full <- normalize_geom_name(final_map_full, "geom")

  if (nrow(final_map_full) != nrow(S)) {
    stop(sprintf(
      "Scored deterministic universe size mismatch: expected %d polygons and got %d.",
      nrow(S), nrow(final_map_full)
    ))
  }

  n_missing_pb <- sum(!is.finite(final_map_full$p_burned))
  if (n_missing_pb > 0) {
    stop(sprintf(
      "Final supervised scoring left %d deterministic polygons without finite p_burned.",
      n_missing_pb
    ))
  }

  final_map <- final_map_full |>
    filter_to_current_year_map() |>
    make_public_final_map()

  scored_gpkg <- file.path(out_score_dir, paste0(prefix, "_deterministic_scored.gpkg"))
  final_gpkg <- file.path(out_map_dir, paste0(prefix, "_final_map.gpkg"))
  counts_csv <- file.path(out_map_dir, paste0(prefix, "_final_map_counts.csv"))

  safe_remove_dataset <- function(path) {
    if (!file.exists(path)) return(invisible(TRUE))
    unlink(path, force = TRUE)
    if (file.exists(path)) {
      stop("Could not overwrite existing dataset: ", path, call. = FALSE)
    }
    invisible(TRUE)
  }

  if (overwrite) safe_remove_dataset(scored_gpkg)
  if (overwrite) safe_remove_dataset(final_gpkg)

  sf::st_write(final_map_full, scored_gpkg, layer = "deterministic_scored", delete_dsn = overwrite, quiet = TRUE)
  sf::st_write(final_map_full, final_gpkg, layer = "deterministic_scored", delete_dsn = overwrite, quiet = TRUE)
  sf::st_write(final_map_full, final_gpkg, layer = "final_map_full", append = TRUE, quiet = TRUE)
  sf::st_write(final_map, final_gpkg, layer = "final_map", append = TRUE, quiet = TRUE)

  tab_counts <- final_map |>
    sf::st_drop_geometry() |>
    dplyr::mutate(has_oof = !is.na(.data[["p_burned_oof"]])) |>
    dplyr::count(source_set, class_input, has_oof, name = "n") |>
    dplyr::arrange(source_set, class_input)
  utils::write.csv(tab_counts, counts_csv, row.names = FALSE)

  burned_like_gpkg <- NULL
  burned_like_counts_csv <- NULL
  burned_like_sf <- NULL

  if (isTRUE(export_burned_like)) {
    burned_like_sf <- final_map_full |>
      dplyr::filter(.data[["source_set"]] == "deterministic", as.character(.data[["class_final"]]) != "keep")

    burned_like_gpkg <- file.path(out_map_dir, burnedlike_basename)
    burned_like_counts_csv <- file.path(
      out_map_dir,
      paste0(tools::file_path_sans_ext(burnedlike_basename), "_counts.csv")
    )

    if (overwrite) safe_remove_dataset(burned_like_gpkg)
    sf::st_write(burned_like_sf, burned_like_gpkg, layer = burnedlike_layer, delete_dsn = overwrite, quiet = TRUE)

    tab_bl <- burned_like_sf |>
      sf::st_drop_geometry() |>
      dplyr::summarise(
        n = dplyr::n(),
        p_burned_min = suppressWarnings(min(.data[["p_burned"]], na.rm = TRUE)),
        p_burned_median = suppressWarnings(stats::median(.data[["p_burned"]], na.rm = TRUE)),
        p_burned_max = suppressWarnings(max(.data[["p_burned"]], na.rm = TRUE))
      )
    utils::write.csv(tab_bl, burned_like_counts_csv, row.names = FALSE)
  }

  msg("Final map rows: total=%d | deterministic=%d",
      nrow(final_map_full),
      sum(sf::st_drop_geometry(final_map_full)$source_set == "deterministic", na.rm = TRUE))
  msg("Current-year public map rows: %d", nrow(final_map))

  invisible(list(
    deterministic_scored = final_map_full,
    final_map_full = final_map_full,
    final_map = final_map,
    burned_like_scored = burned_like_sf,
    files = list(
      scored_gpkg = scored_gpkg,
      final_map_gpkg = final_gpkg,
      counts_csv = counts_csv,
      burned_like_gpkg = burned_like_gpkg,
      burned_like_counts_csv = burned_like_counts_csv
    )
  ))
}
