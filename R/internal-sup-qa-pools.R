# ============================================================
# RELABEL DETERMINISTIC POOLS BEFORE SUPERVISED
# Current active logic:
# - preserve deterministic drop
# - preserve deterministic review
# - demote weak keep -> review
# No promotion of review and no rescue of drop in the active pipeline.
# ============================================================

audit_deterministic_pools <- function(
    internal_sf,
    out_dir = NULL,
    prefix = "pool_qa",
    save_outputs = TRUE,
    verbose = TRUE,
    
    # ------------------------------
    # Columns expected if available
    # ------------------------------
    raw_class_col          = "class_final",
    median_rbr_col         = "median_rbr",
    percentile_col         = "percentile_in_keep",
    pabove_ref_col         = "p_above_keep_ref",
    pabove_q25_col         = "p_above_keep_q25",
    keep_like_col          = "keep_like",
    conf_pix_col           = "conf_pix",
    conf_area_col          = "conf_area",
    npix_col               = "n_pix",
    area_col               = "area_ha",
    
    # ------------------------------
    # Conservative relabelling rules
    # ------------------------------
    # if raw = keep, preserve as keep only if evidence is strong enough
    good_min_percentile    = 0.25,
    good_min_pabove_ref    = 0.50,
    good_min_pabove_q25    = 0.60,
    good_min_area_ha       = 1,
    good_min_npix          = 8,
    
    # force review when confidence fields say low confidence
    force_ambiguous_if_low_conf = TRUE
) {
  
  stopifnot(inherits(internal_sf, "sf"))
  if (nrow(internal_sf) == 0) stop("internal_sf is empty.")
  if (!raw_class_col %in% names(internal_sf)) {
    stop("Missing raw class column: ", raw_class_col)
  }
  
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  x <- internal_sf
  
  # ------------------------------
  # ensure area if missing
  # ------------------------------
  if (!area_col %in% names(x)) {
    x[[area_col]] <- as.numeric(sf::st_area(x)) / 10000
  }
  
  # ------------------------------
  # helper to safely pull numeric/character columns
  # ------------------------------
  get_num <- function(df, col, default = NA_real_) {
    if (!col %in% names(df)) return(rep(default, nrow(df)))
    suppressWarnings(as.numeric(df[[col]]))
  }
  get_chr <- function(df, col, default = NA_character_) {
    if (!col %in% names(df)) return(rep(default, nrow(df)))
    as.character(df[[col]])
  }
  
  raw_class   <- get_chr(x, raw_class_col)
  median_rbr  <- get_num(x, median_rbr_col)
  pct_keep    <- get_num(x, percentile_col)
  p_ref       <- get_num(x, pabove_ref_col)
  p_q25       <- get_num(x, pabove_q25_col)
  keep_like   <- get_num(x, keep_like_col)
  conf_pix    <- get_chr(x, conf_pix_col)
  conf_area   <- get_chr(x, conf_area_col)
  n_pix       <- get_num(x, npix_col)
  area_ha     <- get_num(x, area_col)
  
  if (all(is.na(raw_class))) stop("All raw classes are NA in ", raw_class_col)
  
  # normalize labels
  raw_class <- trimws(raw_class)
  
  # ------------------------------
  # derived logicals
  # ------------------------------
  conf_ok <- rep(TRUE, nrow(x))
  if (isTRUE(force_ambiguous_if_low_conf)) {
    conf_ok <- conf_ok &
      (is.na(conf_pix)  | conf_pix  == "ok") &
      (is.na(conf_area) | conf_area == "ok")
  }
  
  enough_size <- (is.na(area_ha) | area_ha >= good_min_area_ha) &
    (is.na(n_pix) | n_pix >= good_min_npix)
  
  good_signal <- conf_ok &
    enough_size &
    (is.na(pct_keep) | pct_keep >= good_min_percentile) &
    (is.na(p_ref)    | p_ref    >= good_min_pabove_ref) &
    (is.na(p_q25)    | p_q25    >= good_min_pabove_q25)
  
  # ------------------------------
  # audited output label
  # ------------------------------
  class_audited <- raw_class
  qa_reason     <- rep("unchanged", nrow(x))
  
  # ---- keep -> review if weak
  idx_keep <- which(raw_class == "keep")
  if (length(idx_keep) > 0) {
    weak_keep <- idx_keep[!good_signal[idx_keep]]
    if (length(weak_keep) > 0) {
      class_audited[weak_keep] <- "review"
      qa_reason[weak_keep] <- "keep_to_review_weak_signal"
    }
  }
  
  # ---- extra safeguard: if keep_like == 0 and raw keep -> review
  if ("keep_like" %in% names(x)) {
    idx <- which(raw_class == "keep" & !is.na(keep_like) & keep_like == 0)
    if (length(idx) > 0) {
      class_audited[idx] <- "review"
      qa_reason[idx] <- "keep_to_review_keep_like0"
    }
  }
  
  # ------------------------------
  # attach outputs
  # ------------------------------
  x$raw_class     <- raw_class
  x$class_audited <- class_audited
  x$qa_reason     <- qa_reason
  
  # helpful flags
  x$qa_changed <- as.integer(x$raw_class != x$class_audited)
  
  # ------------------------------
  # summaries
  # ------------------------------
  tab_before <- as.data.frame(table(raw_class, useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab_before) <- c("class", "n_before")
  
  tab_after <- as.data.frame(table(class_audited, useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab_after) <- c("class", "n_after")
  
  transitions <- as.data.frame(table(raw_class, class_audited, useNA = "ifany"), stringsAsFactors = FALSE)
  names(transitions) <- c("raw_class", "class_audited", "n")
  transitions <- transitions[transitions$n > 0, , drop = FALSE]
  
  reasons <- as.data.frame(table(qa_reason, useNA = "ifany"), stringsAsFactors = FALSE)
  names(reasons) <- c("qa_reason", "n")
  reasons <- reasons[order(-reasons$n), , drop = FALSE]
  
  summary_df <- merge(tab_before, tab_after, by = "class", all = TRUE)
  summary_df[is.na(summary_df)] <- 0
  
  msg("QA deterministic pools:")
  msg("  raw counts:")
  print(table(raw_class, useNA = "ifany"))
  msg("  audited counts:")
  print(table(class_audited, useNA = "ifany"))
  msg("  changed polygons: %d", sum(x$qa_changed, na.rm = TRUE))
  
  out_paths <- list(
    gpkg = NA_character_,
    summary_csv = NA_character_,
    transitions_csv = NA_character_,
    reasons_csv = NA_character_
  )
  
  if (isTRUE(save_outputs)) {
    if (is.null(out_dir)) stop("out_dir must be provided when save_outputs = TRUE")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    gpkg_path <- file.path(out_dir, paste0(prefix, "_audited_internal.gpkg"))
    summary_csv_path <- file.path(out_dir, paste0(prefix, "_summary.csv"))
    transitions_csv_path <- file.path(out_dir, paste0(prefix, "_transitions.csv"))
    reasons_csv_path <- file.path(out_dir, paste0(prefix, "_reasons.csv"))
    
    if (file.exists(gpkg_path)) file.remove(gpkg_path)
    
    sf::st_write(x, gpkg_path, layer = "audited_internal", quiet = TRUE)
    utils::write.csv(summary_df, summary_csv_path, row.names = FALSE)
    utils::write.csv(transitions, transitions_csv_path, row.names = FALSE)
    utils::write.csv(reasons, reasons_csv_path, row.names = FALSE)
    
    out_paths$gpkg <- gpkg_path
    out_paths$summary_csv <- summary_csv_path
    out_paths$transitions_csv <- transitions_csv_path
    out_paths$reasons_csv <- reasons_csv_path
  }
  
  list(
    audited_sf   = x,
    summary_df   = summary_df,
    transitions  = transitions,
    reasons      = reasons,
    out_paths    = out_paths
  )
}
