build_design_matrix_patches <- function(
    labelled, burned_like,
    
    # key columns (kept so that y and the ids can be saved)
    id_col    = "fire_uid",
    class_col = "class",
    pos_lab   = "burned",
    
    # 1) columns that are NOT features
    id_cols = c("fire_uid","class","source","poly_id","block_id","fold_rep1","fold_rep2"),
    
    # 2) pattern filters (regex)
    drop_regex = c("^fold_rep", "^block_id$", "^source$", "^class$", "^fire_uid$", "^poly_id$"),
    
    # 3) categoricals to factor + "(missing)"
    # 2026-06-05: ecoregions removed from supervised phase; no
    # categorical predictors remain by default.
    cat_cols = character(0),
    
    # 4) hotspot logic (set to NA to disable it)
    hs_n_col    = "hs_used_n",
    hs_conf_col = "hs_conf_mean",
    hs_frp_col  = "hs_frp_max",
    
    # 5) imputation
    median_from = c("labelled", "all"),
    
    # 6) debug
    return_prepared_df = FALSE,

    # 6b) B1 (2026-06-07): defer numeric imputation. When TRUE, every step runs
    # EXACTLY as before (hotspot rules, _isNA companions, factor handling,
    # logical->integer) EXCEPT the per-column median imputation and the global
    # sparse-matrix build are SKIPPED. The function returns the prepared
    # (coerced, _isNA-flagged, but NOT median-imputed) labelled / burned_like
    # feature data.frames and the model column set, so the nested_refit OOF path
    # can fit medians per outer fold (B1b leakage fix). Default FALSE keeps the
    # legacy behaviour byte-identical.
    defer_impute = FALSE,

    # 7) SAVE (nuevo)
    save_dir = NULL,                 # if not NULL -> saves the bundle
    save_prefix = "patch_dm",         # prefijo de archivos
    overwrite = TRUE,
    save_matrix_market = FALSE,       # optional: also saves .mtx
    verbose = TRUE
) {
  median_from <- match.arg(median_from)
  
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is missing.")
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("Package 'dplyr' is missing.")
  
  stopifnot(is.data.frame(labelled), is.data.frame(burned_like))
  if (!id_col %in% names(labelled)) stop("labelled has no id_col: ", id_col)
  if (!id_col %in% names(burned_like)) stop("burned_like has no id_col: ", id_col)
  if (!class_col %in% names(labelled)) stop("labelled has no class_col: ", class_col)
  
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---- 0) choose the columns that enter the model ----
  feat_cols <- setdiff(intersect(names(labelled), names(burned_like)), id_cols)
  if (length(feat_cols) < 3) stop("Too few features left after dropping id_cols.")
  
  if (length(drop_regex) > 0) {
    keep <- rep(TRUE, length(feat_cols))
    for (rgx in drop_regex) keep <- keep & !grepl(rgx, feat_cols)
    feat_cols <- feat_cols[keep]
  }
  if (length(feat_cols) < 3) stop("Too few features left after drop_regex. Check id_cols/drop_regex.")
  
  XL_feat  <- labelled[,    feat_cols, drop = FALSE]
  XBL_feat <- burned_like[, feat_cols, drop = FALSE]
  X_all <- dplyr::bind_rows(XL_feat, XBL_feat)
  nL <- nrow(labelled)
  
  # Keep only explicitly requested categorical predictors.
  other_cat_cols <- names(X_all)[vapply(X_all, function(x) is.character(x) || is.factor(x), logical(1))]
  other_cat_cols <- setdiff(other_cat_cols, intersect(cat_cols, names(X_all)))
  if (length(other_cat_cols) > 0) {
    msg(
      "Dropping categorical columns not listed in cat_cols: %s",
      paste(other_cat_cols, collapse = ", ")
    )
    X_all <- X_all[, setdiff(names(X_all), other_cat_cols), drop = FALSE]
  }
  
  # ---- 1) hotspots rules ----
  # 0.6.2 (2026-06-10): the `hs_any` synthesis block was REMOVED here.
  # `hs_any` (= as.integer(hs_used_n > 0)) was a redundant helper, never
  # whitelisted as of 0.6.2 and never a base GPKG feature, so it never
  # reaches the recipe. Removing the synthesis leaves the resolved
  # feature space unchanged (per Natalia's decision).

  # Gate 1D.8 (2026-06-09): the `_isNA` companions for hs_conf / hs_frp are no
  # longer synthesised inline here -- they are produced by the SHARED helper
  # (.of_synthesize_isna_companions) below, identically to the FINAL / scoring
  # paths, so OOF and FINAL share ONE feature space. We still apply the
  # hotspot-specific value rule (when hs_n == 0 a missing conf/frp is a real 0,
  # not "unknown"), and we do it AFTER the shared `_isNA` synthesis captures the
  # ORIGINAL missingness (so the flag reflects the raw NA, not the zeroed value).

  # ---- 2) categoricals ----
  to_factor_with_missing <- function(x) {
    x <- as.factor(x)
    x <- addNA(x)
    lev <- levels(x); lev[is.na(lev)] <- "(missing)"
    levels(x) <- lev
    x
  }
  
  cat_cols_present <- intersect(cat_cols, names(X_all))
  for (nm in cat_cols_present) X_all[[nm]] <- to_factor_with_missing(X_all[[nm]])
  
  factor_levels <- list()
  fac_names <- names(X_all)[vapply(X_all, is.factor, logical(1))]
  for (nm in fac_names) factor_levels[[nm]] <- levels(X_all[[nm]])
  
  # sparse.model.matrix() fails on factors with a single observed level.
  single_level_factors <- fac_names[vapply(
    X_all[fac_names],
    function(x) nlevels(x) < 2,
    logical(1)
  )]
  if (length(single_level_factors) > 0) {
    msg(
      "Dropping single-level factor columns before model matrix: %s",
      paste(single_level_factors, collapse = ", ")
    )
    X_all <- X_all[, setdiff(names(X_all), single_level_factors), drop = FALSE]
    factor_levels <- factor_levels[setdiff(names(factor_levels), single_level_factors)]
    fac_names <- setdiff(fac_names, single_level_factors)
  }
  
  # ---- 3) numerics: NA flag + imputation ----
  logical_cols <- names(X_all)[vapply(X_all, is.logical, logical(1))]
  if (length(logical_cols) > 0) {
    msg(
      "Converting logical columns to integer before model matrix: %s",
      paste(logical_cols, collapse = ", ")
    )
    for (nm in logical_cols) X_all[[nm]] <- as.integer(X_all[[nm]])
  }
  
  # Gate 1D.8 (2026-06-09): SHARED `_isNA` synthesis. Replaces the old inline
  # per-numeric-column loop (and the inline hotspot-block flags above) with the
  # ONE shared helper used by the FINAL / scoring paths, so OOF and FINAL build
  # the IDENTICAL feature space (same `_isNA` companion set, same canonical
  # base-block-then-indicator-block order). Bug-8 dedupe (no `<x>_isNA_isNA`,
  # regenerate any input-provided companion) lives inside the helper.
  X_all <- .of_synthesize_isna_companions(X_all)

  # Hotspot value rule (applied AFTER `_isNA` synthesis so the flag reflects the
  # ORIGINAL missingness): when there were no hotspots (hs_n == 0), a missing
  # conf/frp value is a real 0, not "unknown".
  if (!is.na(hs_n_col) && !is.na(hs_conf_col) &&
      all(c(hs_n_col, hs_conf_col) %in% names(X_all))) {
    X_all[[hs_conf_col]][is.na(X_all[[hs_conf_col]]) & X_all[[hs_n_col]] == 0] <- 0
  }
  if (!is.na(hs_n_col) && !is.na(hs_frp_col) &&
      all(c(hs_n_col, hs_frp_col) %in% names(X_all))) {
    X_all[[hs_frp_col]][is.na(X_all[[hs_frp_col]]) & X_all[[hs_n_col]] == 0] <- 0
  }

  num_names <- names(X_all)[vapply(X_all, is.numeric, logical(1))]
  idx_median <- if (median_from == "labelled") seq_len(nL) else seq_len(nrow(X_all))
  medians_used <- list()

  if (!isTRUE(defer_impute)) {
    for (nm in num_names) {
      med <- stats::median(X_all[[nm]][idx_median], na.rm = TRUE)
      if (is.na(med)) med <- 0
      medians_used[[nm]] <- med
      X_all[[nm]][is.na(X_all[[nm]])] <- med
    }
  }
  # B1: defer_impute leaves numeric NAs in place (per-fold medians are fit later
  # inside the nested-refit core on the outer-train rows only).

  # ---- 3b) B1 deferred-impute early return ----
  # Return the prepared (coerced, _isNA-flagged, hotspot-derived, but NOT
  # median-imputed) labelled / burned_like feature frames + the model column
  # set, so the nested_refit OOF path can fit medians per outer fold. No global
  # matrix is built or saved here (the per-fold core builds its own matrices).
  if (isTRUE(defer_impute)) {
    prepared_labelled    <- X_all[seq_len(nL), , drop = FALSE]
    prepared_burned_like <- X_all[(nL + 1):nrow(X_all), , drop = FALSE]
    y <- ifelse(labelled[[class_col]] == pos_lab, 1L, 0L)
    return(list(
      deferred             = TRUE,
      prepared_labelled    = prepared_labelled,
      prepared_burned_like = prepared_burned_like,
      model_cols           = names(X_all),
      num_cols             = num_names,
      factor_levels        = factor_levels,
      y                    = y,
      id_labelled          = labelled[[id_col]],
      id_burned_like       = burned_like[[id_col]],
      feat_cols_used       = feat_cols,
      n_labelled           = nL
    ))
  }

  # ---- 4) matriz sparse ----
  X_all_mat <- Matrix::sparse.model.matrix(~ . - 1, data = X_all, na.action = stats::na.pass)
  XL_mat  <- X_all_mat[seq_len(nL), , drop = FALSE]
  XBL_mat <- X_all_mat[(nL + 1):nrow(X_all_mat), , drop = FALSE]
  
  prep <- list(
    feat_cols_used  = feat_cols,
    id_col          = id_col,
    class_col       = class_col,
    pos_lab         = pos_lab,
    id_cols_used    = id_cols,
    drop_regex_used = drop_regex,
    cat_cols_used   = cat_cols_present,
    dropped_non_cat_cols = other_cat_cols,
    dropped_single_level_factors = single_level_factors,
    logical_cols_as_integer = logical_cols,
    factor_levels   = factor_levels,
    medians_used    = medians_used,
    median_from     = median_from,
    matrix_colnames = colnames(X_all_mat),
    n_labelled      = nL
  )
  
  # target + ids (kept for saving)
  y <- ifelse(labelled[[class_col]] == pos_lab, 1L, 0L)
  id_labelled    <- labelled[[id_col]]
  id_burned_like <- burned_like[[id_col]]
  
  out <- list(
    XL_mat = XL_mat,
    XBL_mat = XBL_mat,
    y = y,
    id_labelled = id_labelled,
    id_burned_like = id_burned_like,
    prep = prep
  )
  if (return_prepared_df) out$X_all_prepared <- X_all
  
  # ---- 5) SAVE bundle (opcional) ----
  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    
    rds_path <- file.path(save_dir, paste0(save_prefix, "_design_bundle.rds"))
    if (file.exists(rds_path) && !overwrite) stop("Ya existe: ", rds_path)
    
    bundle <- list(
      XL_mat = XL_mat,
      XBL_mat = XBL_mat,
      y = y,
      id_labelled = id_labelled,
      id_burned_like = id_burned_like,
      prep = prep,
      # B6 (2026-06-06): the `created = Sys.time()` field was removed from this
      # checksummed `_design_bundle.rds` artifact so re-runs are
      # byte-reproducible. It was metadata only (nothing reads bundle$meta$
      # created). `R` (R.version.string) is deterministic for a fixed R and is
      # kept for provenance.
      meta = list(
        R = R.version.string
      )
    )
    
    saveRDS(bundle, rds_path)
    msg("Saved bundle: %s", rds_path)
    
    if (isTRUE(save_matrix_market)) {
      Matrix::writeMM(XL_mat,  file.path(save_dir, paste0(save_prefix, "_XL_mat.mtx")))
      Matrix::writeMM(XBL_mat, file.path(save_dir, paste0(save_prefix, "_XBL_mat.mtx")))
      utils::write.csv(data.frame(y = y), file.path(save_dir, paste0(save_prefix, "_y.csv")), row.names = FALSE)
      msg("Saved MatrixMarket + y.csv in: %s", save_dir)
    }
    
    out$saved_bundle_path <- rds_path
  }
  
  out
}
