# =============================================================================
# PHASE 2 (artifact_hard) -- CONSOLIDATED supervised training-pool layer.
#
# A single, flat, one-row-per-example view of the supervised training pool,
# written alongside (NOT replacing) the existing separate burned / unburned
# pool outputs. It records EXACTLY the rows and weights the training engine
# receives, plus the deterministic origin and the artifact_hard selection
# provenance, so the pools can be inspected and visually validated before any
# further scientific change.
#
# Faithfulness contract (enforced by tests):
#   * used_in_training == TRUE  <=>  the row is in the FINAL capped training set
#     (`training_ok`, the exact frame the refit DMatrix is built from);
#   * sample_weight on the used rows is re-resolved with the SAME
#     .of_resolve_sample_weights() the engine calls, with the SAME inputs
#     (source_weights + total_weight_ratio) -> identical per-source totals;
#   * the layer NEVER carries geometry/list-columns as model features; the model
#     features are plain numeric descriptive columns only.
#
# ADDITIVE + OFF-aware: with artifact_hard disabled the layer still describes the
# baseline pool (positives + random + otsu, all weight 1, no artifact_hard rows);
# it is a NEW output file and does not change the model or any existing artifact.
# =============================================================================

# Canonical column order of the consolidated layer (spec). Separates the three
# distinct artifact_hard states: satisfies-the-rule (artifact_hard_eligible),
# switched-on-in-config (artifact_hard_enabled), effectively-promoted
# (artifact_hard_used); plus whether the final row entered the DMatrix
# (used_in_training). artifact_hard candidates are ALWAYS shown (even with
# enabled = FALSE) so the rule can be inspected without running two pipelines.
.of_training_pool_layer_cols <- function() {
  c("year", "fire_uid", "poly_id", "training_label",
    "original_pool_source", "pool_source",
    "artifact_hard_eligible", "artifact_hard_enabled", "artifact_hard_used",
    "used_in_training", "sample_weight", "source_total_weight",
    "deterministic_class", "deterministic_reason", "fold_id",
    "has_oof_prediction", "p_burned_oof",
    "rbr_med", "rbr_aw_med", "persist_ratio", "persist_delta", "area_ha",
    "doy_iqr", "artifact_hard_branch", "artifact_hard_rbr_threshold",
    "artifact_hard_weight_ratio", "review_priority", "review_tier",
    "label_confidence", "VISUAL")
}

# Weights of the deterministic visual-review priority score (sum to 1; tunable).
# Higher review_priority = more likely a REAL FIRE wrongly flagged as an
# artifact (bright + persistent + large) -> verify first.
.of_review_priority_weights <- function() c(rbr = 0.45, persist = 0.30, area = 0.25)

# Base label-reliability per provenance bucket (tunable). The trusted buckets
# (deterministic keep positives, random/otsu negatives) carry a high flat base;
# artifact_hard candidates get a VARIABLE confidence = 1 - review_priority/100
# (their "unburned" label is least trustworthy exactly when they look fire-like).
.of_label_confidence_base <- function() {
  c(high_confidence_keep = 0.90, random = 0.85, otsu = 0.85)
}

# Resolve the FINAL training membership (the capped frame the engine will train
# on) PURELY from the labelled pool + config, WITHOUT training a model. Mirrors
# the engine's FINAL assembly exactly (same eligibility resolver + capping helper
# + final_sampling_seed), so it can be computed BEFORE OOF/FINAL run. Proven
# equivalent to the engine's `training_ok` by the runner<->wrapper test.
.of_engine_training_ok <- function(train_features, config, id_col = "fire_uid") {
  d <- if (inherits(train_features, "sf")) sf::st_drop_geometry(train_features) else as.data.frame(train_features)
  if (!all(c(id_col, "class") %in% names(d))) {
    stop(".of_engine_training_ok(): train_features needs '", id_col,
         "' and 'class' columns.", call. = FALSE)
  }
  tc <- config$train_control
  ah_on <- isTRUE(config$negative_pool_params$artifact_hard$enabled)
  ids <- as.character(d[[id_col]])
  elig <- .of_resolve_supervised_eligibility(
    id       = ids,
    class    = as.character(d[["class"]]),
    source   = if ("source" %in% names(d)) as.character(d[["source"]]) else rep(NA_character_, nrow(d)),
    neg_type = if ("neg_type" %in% names(d)) as.character(d[["neg_type"]]) else rep(NA_character_, nrow(d)),
    random_background_source = "random_burnable_background",
    otsu_unburned_source = "otsu_patch_residual",
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep"),
    origin_stage = "FINAL",
    artifact_hard_source = if (ah_on) "artifact_hard" else character(0))
  caps <- c(random = tc$caps$random, otsu = tc$caps$otsu)
  if ("artifact_hard" %in% names(elig$negatives_by_bucket)) caps["artifact_hard"] <- Inf
  cp <- .of_cap_negative_buckets(
    positive_idx = elig$positive_idx, negatives_by_bucket = elig$negatives_by_bucket,
    n_burned = length(elig$positive_idx), caps = caps,
    seed = tc$seeds$final_sampling_seed, id = ids, context = "FINAL")
  sel <- cp$selected_indices
  keep <- !duplicated(ids[sel])
  train_features[sel[keep], , drop = FALSE]
}

# Map (class, source) -> pool_source. positives (burned) -> high_confidence_keep
# regardless of source tag; the three negative buckets map by their source tag.
.of_pool_source_of <- function(class, source) {
  cls <- as.character(class); src <- as.character(source)
  out <- rep(NA_character_, length(cls))
  out[!is.na(src) & src == "random_burnable_background"] <- "random"
  out[!is.na(src) & src == "otsu_patch_residual"]        <- "otsu"
  out[!is.na(src) & src == "artifact_hard"]              <- "artifact_hard"
  out[!is.na(cls) & cls == "burned"]                     <- "high_confidence_keep"
  out
}

#' Build the consolidated supervised training-pool layer (Phase 2).
#'
#' @param train_features sf. The AUGMENTED labelled training universe (positives
#'   + random + otsu + any promoted artifact_hard rows), AFTER feature extraction
#'   and fold assignment. One row per available training example.
#' @param scoring_features sf / data.frame, or NULL. The deterministic scoring
#'   universe. Used to compute the artifact_hard ELIGIBLE candidate set ALWAYS
#'   (the selection rule, independent of `enabled`): when `enabled = FALSE` the
#'   eligible candidates are APPENDED to the layer as candidate-only rows
#'   (`artifact_hard_eligible = TRUE`, `artifact_hard_used = FALSE`,
#'   `sample_weight = NA`) so they can be inspected WITHOUT entering training;
#'   when `enabled = TRUE` they are already promoted into `train_features`. NULL
#'   skips candidate computation (eligibility falls back to the promoted rows).
#' @param training_ok sf / data.frame, or NULL. The FINAL capped training set
#'   actually fed to the refit DMatrix (the engine's `training_ok`); defines
#'   `used_in_training` and the per-source weight totals. When NULL (the
#'   PRE-MODEL path) it is derived purely from `train_features` + `config` via
#'   the same eligibility resolver + capping helper + seed the engine uses, so
#'   the consolidated layer can be written BEFORE OOF/FINAL run.
#' @param oof_agg data.frame or CSV path with the per-example OOF aggregate
#'   (id column + `p_oof_mean`), or NULL. NULL (the pre-model path) leaves
#'   `p_burned_oof` NA and `has_oof_prediction` FALSE.
#' @param promo The `promote_artifact_hard_negatives()` result (or NULL when the
#'   feature is OFF). Used only as a fallback source of the per-row branch /
#'   threshold when those columns are not already carried on `train_features`.
#' @param config The resolved supervised config (for source_weights +
#'   artifact_hard$total_weight_ratio, re-resolved identically to the engine).
#' @param target_year integer scalar stamped into the `year` column.
#' @param id_col character; the stable per-example id (default "fire_uid").
#' @return An sf with the canonical consolidated-layer columns + geometry.
#' @keywords internal
#' @noRd
.of_build_supervised_training_pool <- function(train_features,
                                               scoring_features = NULL,
                                               training_ok = NULL,
                                               oof_agg = NULL, promo = NULL,
                                               config, target_year,
                                               id_col = "fire_uid") {
  if (is.null(train_features) || !inherits(train_features, "sf")) {
    train_features <- sf::st_as_sf(train_features)
  }
  ah_cfg     <- config$negative_pool_params$artifact_hard
  ah_enabled <- isTRUE(ah_cfg$enabled)
  twr        <- ah_cfg$total_weight_ratio %||% 0.10
  sw_pin     <- config$negative_pool_params$source_weights

  # ------------------------------------------------------------------------
  # (1) Compute the ELIGIBLE artifact_hard candidate set ALWAYS -- it is the
  # selection RULE, independent of enabled. The rbr_med reference uses the
  # ORIGINAL negatives (random/otsu only, EXCLUDING any already-promoted
  # artifact_hard rows), so the eligible set is IDENTICAL whether enabled or not.
  # ------------------------------------------------------------------------
  elig_polyids   <- character(0)
  branch_by_poly <- character(0)
  thr_rbr_med    <- NA_real_
  cand_orig      <- NULL
  if (!is.null(scoring_features) && nrow(scoring_features) > 0L) {
    tdf  <- sf::st_drop_geometry(train_features)
    cls0 <- if ("class" %in% names(tdf)) as.character(tdf[["class"]]) else character(nrow(tdf))
    src0 <- if ("source" %in% names(tdf)) as.character(tdf[["source"]]) else rep(NA_character_, nrow(tdf))
    rbr0 <- if ("rbr_med" %in% names(tdf)) suppressWarnings(as.numeric(tdf[["rbr_med"]])) else rep(NA_real_, nrow(tdf))
    ref_kind <- if (is.null(ah_cfg$rbr_med_reference)) "negative" else as.character(ah_cfg$rbr_med_reference)
    reference_rbr_med <- if (identical(ref_kind, "positive")) {
      rbr0[!is.na(cls0) & cls0 == "burned"]
    } else {
      rbr0[!is.na(cls0) & cls0 == "unburned" & (is.na(src0) | src0 != "artifact_hard")]
    }
    sel <- .of_select_artifact_hard(
      scoring_features = scoring_features, reference_rbr_med = reference_rbr_med,
      params = ah_cfg, enable_derived = FALSE)
    thr_rbr_med <- sel$threshold_rbr_med
    if (nrow(sel$artifact_hard) > 0L) {
      seldf <- sf::st_drop_geometry(sel$artifact_hard)
      elig_polyids   <- as.character(seldf[["poly_id"]])
      branch_by_poly <- stats::setNames(as.character(seldf[["artifact_hard_branch"]]),
                                        elig_polyids)
      # ORIGINAL (unstamped) scoring rows for the eligible candidates, so a
      # non-promoted candidate keeps its deterministic source/class.
      spoly <- as.character(sf::st_drop_geometry(scoring_features)[["poly_id"]])
      cand_orig <- scoring_features[spoly %in% elig_polyids, , drop = FALSE]
    }
  }

  # ------------------------------------------------------------------------
  # (2) Layer universe = the training pool, PLUS the eligible candidates that are
  # NOT already in the training pool. When enabled, the candidates were promoted
  # INTO train_features (so they are already present); when disabled they live
  # only in the scoring universe and are APPENDED here as candidate-only rows so
  # they ALWAYS appear -- WITHOUT entering training.
  # ------------------------------------------------------------------------
  carry <- c(id_col, "poly_id", "class", "source", "class_final", "reason_1",
             "fold_rep1", "rbr_med", "rbr_aw_med", "persist_ratio",
             "persist_delta", "area_ha", "doy_iqr")
  pick <- function(x) {
    d <- sf::st_drop_geometry(x); nn <- nrow(d)
    cols <- lapply(carry, function(cn) if (cn %in% names(d)) d[[cn]] else rep(NA, nn))
    names(cols) <- carry
    sf::st_sf(as.data.frame(cols, stringsAsFactors = FALSE),
              geometry = sf::st_geometry(x))
  }
  base_sf <- pick(train_features)
  is_cand_only <- rep(FALSE, nrow(base_sf))
  if (!is.null(cand_orig) && nrow(cand_orig) > 0L) {
    base_poly <- as.character(sf::st_drop_geometry(base_sf)[["poly_id"]])
    cand_poly <- as.character(sf::st_drop_geometry(cand_orig)[["poly_id"]])
    add <- which(!(cand_poly %in% base_poly))
    if (length(add) > 0L) {
      cand_sf <- pick(cand_orig[add, , drop = FALSE])
      base_sf <- rbind(base_sf, cand_sf)
      is_cand_only <- c(is_cand_only, rep(TRUE, length(add)))
    }
  }

  # ------------------------------------------------------------------------
  # (3) Derive the FINAL training membership (the capped frame; NULL -> derive
  # purely from train_features + config, no model needed). Candidate-only rows
  # are never in it.
  # ------------------------------------------------------------------------
  if (is.null(training_ok)) {
    training_ok <- .of_engine_training_ok(train_features, config, id_col = id_col)
  }
  tok_uid <- as.character(sf::st_drop_geometry(sf::st_as_sf(training_ok))[[id_col]])

  # ------------------------------------------------------------------------
  # (4) Per-row columns on the assembled universe.
  # ------------------------------------------------------------------------
  df <- sf::st_drop_geometry(base_sf)
  n  <- nrow(df)
  getcol <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(NA, n)
  num <- function(nm) {
    v <- getcol(nm)
    if (is.factor(v)) v <- as.character(v)
    if (is.character(v)) v <- gsub(",", ".", v, fixed = TRUE)
    suppressWarnings(as.numeric(v))
  }
  fire_uid <- as.character(getcol(id_col))
  poly     <- as.character(getcol("poly_id"))
  cls      <- as.character(getcol("class"))
  src      <- as.character(getcol("source"))

  # effective pool_source (after promotion) + original (before promotion).
  pool_source <- .of_pool_source_of(cls, src)
  pool_source[is_cand_only] <- "deterministic_drop"
  original_pool_source <- pool_source
  original_pool_source[!is_cand_only & !is.na(pool_source) &
                         pool_source == "artifact_hard"] <- "deterministic_drop"

  # the three distinct states.
  artifact_hard_eligible <- poly %in% elig_polyids
  artifact_hard_enabled  <- rep(ah_enabled, n)
  artifact_hard_used     <- !is_cand_only & !is.na(pool_source) &
                            pool_source == "artifact_hard"

  # effective training label: the label used in training; NA for candidate-only.
  training_label <- cls
  training_label[is_cand_only] <- NA_character_

  # deterministic origin.
  det_class <- as.character(getcol("class_final"))
  det_class[is.na(det_class) & pool_source == "high_confidence_keep"] <- "keep"
  det_class[is.na(det_class) & (artifact_hard_used | is_cand_only)]   <- "drop"
  det_reason <- as.character(getcol("reason_1"))

  # used_in_training: membership in the FINAL capped frame.
  used_in_training <- fire_uid %in% tok_uid

  # sample_weight: re-resolve with the SAME engine helper/inputs; NA when the row
  # does not enter training (candidate-only OR capped out).
  sample_weight <- rep(NA_real_, n)
  if (length(tok_uid) > 0L) {
    tok <- sf::st_drop_geometry(sf::st_as_sf(training_ok))
    sw_res <- .of_resolve_sample_weights(
      class  = as.character(tok[["class"]]),
      source = if ("source" %in% names(tok)) as.character(tok[["source"]]) else rep(NA_character_, nrow(tok)),
      source_weights = sw_pin, artifact_hard_source = "artifact_hard",
      total_weight_ratio = twr)
    w_used <- sw_res$weights
    if (is.null(w_used)) w_used <- rep(1.0, nrow(tok))  # unweighted => weight 1
    w_by_id <- stats::setNames(w_used, as.character(tok[[id_col]]))
    sample_weight[used_in_training] <- as.numeric(w_by_id[fire_uid[used_in_training]])
  }

  # source_total_weight: sum of USED rows' weight per effective pool_source.
  source_total_weight <- rep(NA_real_, n)
  for (ps in unique(stats::na.omit(pool_source))) {
    in_ps <- !is.na(pool_source) & pool_source == ps
    source_total_weight[in_ps] <- sum(sample_weight[in_ps & used_in_training],
                                      na.rm = TRUE)
  }

  fold_id <- as.character(getcol("fold_rep1"))

  p_burned_oof <- rep(NA_real_, n)
  if (!is.null(oof_agg)) {
    oa <- if (is.character(oof_agg) && length(oof_agg) == 1L && file.exists(oof_agg)) {
      utils::read.csv(oof_agg, stringsAsFactors = FALSE)
    } else as.data.frame(oof_agg)
    if (id_col %in% names(oa) && "p_oof_mean" %in% names(oa)) {
      p_by_id <- stats::setNames(suppressWarnings(as.numeric(oa[["p_oof_mean"]])),
                                 as.character(oa[[id_col]]))
      p_burned_oof <- as.numeric(p_by_id[fire_uid])
    }
  }
  has_oof_prediction <- is.finite(p_burned_oof)

  # artifact_hard provenance for EVERY eligible row (promoted AND candidate-only).
  ah_branch <- if (length(branch_by_poly) > 0L) as.character(branch_by_poly[poly]) else rep(NA_character_, n)
  ah_branch[!artifact_hard_eligible] <- NA_character_
  ah_thr    <- ifelse(artifact_hard_eligible, as.numeric(thr_rbr_med), NA_real_)
  # the configured weight ratio applies to every eligible candidate (it is the
  # value that WOULD be used if/when enabled); NA for non-candidates.
  artifact_hard_weight_ratio <- ifelse(artifact_hard_eligible, as.numeric(twr), NA_real_)

  # ---- review_priority: deterministic visual-review order over candidates ----
  # 0..100; higher = more likely a REAL FIRE wrongly flagged (bright + persistent
  # + large) -> verify first. Computed ONLY for artifact_hard candidates.
  clamp01 <- function(z) pmin(1, pmax(0, z))
  rbr_v     <- num("rbr_med"); persist_v <- num("persist_ratio"); area_v <- num("area_ha")
  rbr_gate  <- as.numeric(thr_rbr_med)
  pos_rbr   <- rbr_v[!is.na(pool_source) & pool_source == "high_confidence_keep"]
  rbr_posref <- stats::median(pos_rbr, na.rm = TRUE)
  persist_max <- ah_cfg$persist_ratio_max %||% 0.35
  area_min    <- ah_cfg$area_ha_min %||% 500
  cand_area   <- area_v[artifact_hard_eligible]
  area_p95    <- if (any(is.finite(cand_area))) {
    stats::quantile(cand_area, 0.95, names = FALSE, na.rm = TRUE)
  } else as.numeric(area_min)
  W <- .of_review_priority_weights()
  review_priority <- rep(NA_real_, n)
  if (any(artifact_hard_eligible) && is.finite(rbr_gate) && is.finite(rbr_posref)) {
    denom_rbr <- max(rbr_posref - rbr_gate, .Machine$double.eps)
    denom_a   <- max(log1p(area_p95) - log1p(as.numeric(area_min)), .Machine$double.eps)
    s_rbr     <- clamp01((rbr_v - rbr_gate) / denom_rbr)
    s_persist <- clamp01(persist_v / persist_max)
    s_area    <- clamp01((log1p(area_v) - log1p(as.numeric(area_min))) / denom_a)
    rp <- 100 * (W[["rbr"]] * s_rbr + W[["persist"]] * s_persist + W[["area"]] * s_area)
    review_priority[artifact_hard_eligible] <- round(as.numeric(rp)[artifact_hard_eligible], 1)
  }
  review_tier <- rep(NA_character_, n)
  hi <- artifact_hard_eligible & is.finite(review_priority)
  review_tier[hi & review_priority >= 66]                        <- "high"
  review_tier[hi & review_priority >= 33 & review_priority < 66] <- "medium"
  review_tier[hi & review_priority < 33]                         <- "low"

  # ---- label_confidence: trust in the (effective/prospective) training label --
  # Trusted buckets get a high flat base; artifact_hard candidates get a VARIABLE
  # confidence = 1 - review_priority/100 (their "unburned" label is least
  # trustworthy exactly when they look like real fire). 0..1.
  base_conf <- .of_label_confidence_base()
  label_confidence <- rep(NA_real_, n)
  for (ps in names(base_conf)) {
    label_confidence[!is.na(pool_source) & pool_source == ps] <- base_conf[[ps]]
  }
  label_confidence[artifact_hard_eligible] <-
    round(clamp01(1 - review_priority[artifact_hard_eligible] / 100), 3)

  out <- data.frame(
    year                        = rep(as.integer(target_year), n),
    fire_uid                    = fire_uid,
    poly_id                     = poly,
    training_label              = training_label,
    original_pool_source        = original_pool_source,
    pool_source                 = pool_source,
    artifact_hard_eligible      = artifact_hard_eligible,
    artifact_hard_enabled       = artifact_hard_enabled,
    artifact_hard_used          = artifact_hard_used,
    used_in_training            = used_in_training,
    sample_weight               = sample_weight,
    source_total_weight         = source_total_weight,
    deterministic_class         = det_class,
    deterministic_reason        = det_reason,
    fold_id                     = fold_id,
    has_oof_prediction          = has_oof_prediction,
    p_burned_oof                = p_burned_oof,
    rbr_med                     = num("rbr_med"),
    rbr_aw_med                  = num("rbr_aw_med"),
    persist_ratio               = num("persist_ratio"),
    persist_delta               = num("persist_delta"),
    area_ha                     = num("area_ha"),
    doy_iqr                     = num("doy_iqr"),
    artifact_hard_branch        = ah_branch,
    artifact_hard_rbr_threshold = ah_thr,
    artifact_hard_weight_ratio  = artifact_hard_weight_ratio,
    review_priority             = review_priority,
    review_tier                 = review_tier,
    label_confidence            = label_confidence,
    # VISUAL: BLANK column for the user to fill during visual review
    # (1 = label confirmed, 0 = not confirmed; NA = not yet reviewed).
    VISUAL                      = rep(NA_integer_, n),
    stringsAsFactors            = FALSE
  )
  out <- out[, .of_training_pool_layer_cols(), drop = FALSE]
  sf::st_sf(out, geometry = sf::st_geometry(base_sf))
}

#' Build + write the consolidated supervised training-pool layer to GPKG.
#'
#' Thin wrapper that builds the layer (\code{.of_build_supervised_training_pool})
#' and writes it as `supervised_training_pool.gpkg` (layer
#' "supervised_training_pool") under `out_dir`, honouring `overwrite`.
#'
#' @return The written file path (invisibly), or NULL when nothing was written.
#' @keywords internal
#' @noRd
.of_write_supervised_training_pool <- function(train_features,
                                               scoring_features = NULL,
                                               training_ok = NULL,
                                               oof_agg = NULL, promo = NULL,
                                               config, target_year, out_dir,
                                               overwrite = TRUE,
                                               id_col = "fire_uid") {
  layer <- .of_build_supervised_training_pool(
    train_features = train_features, scoring_features = scoring_features,
    training_ok = training_ok, oof_agg = oof_agg, promo = promo, config = config,
    target_year = target_year, id_col = id_col)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "supervised_training_pool.gpkg")
  if (!isTRUE(overwrite) && file.exists(path)) return(invisible(path))
  # Robust overwrite: a GPKG left in SQLite WAL mode leaves `-wal` / `-shm`
  # sidecars that keep the dataset "present" even after the .gpkg is unlinked, so
  # a re-run would hit "Dataset already exists". Remove the .gpkg AND every
  # sidecar before writing a clean dataset.
  for (sfx in c("", "-wal", "-shm", "-journal")) {
    f <- paste0(path, sfx)
    if (file.exists(f)) unlink(f, force = TRUE)
  }
  sf::st_write(layer, path, layer = "supervised_training_pool", quiet = TRUE)
  invisible(path)
}
