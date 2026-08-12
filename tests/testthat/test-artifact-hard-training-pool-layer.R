# =============================================================================
# PHASE 2 (artifact_hard): CONSOLIDATED supervised training-pool layer.
#
# The layer ALWAYS computes + shows the artifact_hard candidates (rows that
# satisfy the selection rule), even with artifact_hard$enabled = FALSE, WITHOUT
# letting them enter training. Three states are kept distinct:
#   * artifact_hard_eligible : satisfies the selection rule (rule-only);
#   * artifact_hard_enabled  : the resolved-config switch;
#   * artifact_hard_used     : effectively promoted as artifact_hard;
# plus used_in_training (entered the FINAL DMatrix). All on a small synthetic
# fixture (no rasters / no EFFIS).
# =============================================================================

ns <- asNamespace("OtsuFire")
build_pool   <- ns$.of_build_supervised_training_pool
resolve_w    <- ns$.of_resolve_sample_weights
resolve_elig <- ns$.of_resolve_supervised_eligibility
cap_buckets  <- ns$.of_cap_negative_buckets

FEAT <- c("rbr_med", "elev_med")

# --- training pool (positives + random + otsu; NO artifact_hard) ------------
mk_train_base <- function(seed = 31, n_pos = 24, n_rand = 18, n_otsu = 12) {
  set.seed(seed)
  n <- n_pos + n_rand + n_otsu
  cls <- c(rep("burned", n_pos), rep("unburned", n_rand + n_otsu))
  src <- c(rep("internal_keep_qc", n_pos),
           rep("random_burnable_background", n_rand),
           rep("otsu_patch_residual", n_otsu))
  ngt <- c(rep(NA_character_, n_pos + n_rand), rep("otsu_patch_drop", n_otsu))
  df <- data.frame(
    fire_uid = as.character(seq_len(n)), poly_id = as.character(seq_len(n)),
    class = cls, source = src, neg_type = ngt,
    class_final = c(rep("keep", n_pos), rep(NA_character_, n_rand + n_otsu)),
    reason_1 = NA_character_,
    block_id = paste0("blk_", sample(seq_len(8), n, TRUE)),
    fold_rep1 = sample(1:5, n, TRUE), fold_rep2 = sample(1:5, n, TRUE),
    rbr_med = c(rnorm(n_pos, 0.6, 0.04), rnorm(n_rand, 0.15, 0.02),
                rnorm(n_otsu, 0.2, 0.02)),
    rbr_aw_med = rnorm(n, 0.2, 0.03),
    persist_ratio = rnorm(n, 0.4, 0.05), persist_delta = rep(-10, n),
    area_ha = rep(100, n), doy_iqr = rep(2, n), elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE)
  geom <- sf::st_sfc(lapply(seq_len(n), function(i)
    sf::st_point(c(runif(1, 0, 1e5), runif(1, 0, 1e5)))), crs = 3035)
  sf::st_sf(df, geometry = geom)
}

# --- deterministic scoring universe: 6 eligible drops + 3 non-eligible drops
#     + 2 keep + 1 review. Eligible = drop & rbr high & persist low & pd<=-100. -
mk_scoring <- function(seed = 52) {
  set.seed(seed)
  poly <- as.character(2000:2012); n <- length(poly)              # disjoint ids
  class_final <- c(rep("drop", 9), "keep", "keep", "review", "drop")
  rbr_med <- c(rep(0.45, 6),  0.05, 0.45, 0.45,  0.5, 0.5, 0.5, 0.45)
  persist_ratio <- c(rep(0.15, 6), 0.15, 0.60, 0.15, 0.5, 0.5, 0.5, 0.15)
  persist_delta <- c(rep(-200, 6), -200, -200, -10, -10, -10, -10, -10)
  area_ha <- c(rep(800, 6), 800, 800, 100, 100, 100, 100, 100)
  doy_iqr <- c(rep(0, 6), 0, 0, 5, 0, 0, 0, 5)
  df <- data.frame(
    fire_uid = paste0("SC_", poly), poly_id = poly,
    class = NA_character_, source = NA_character_, neg_type = NA_character_,
    class_final = class_final, reason_1 = NA_character_,
    fold_rep1 = NA_integer_, fold_rep2 = NA_integer_,
    rbr_med = rbr_med, rbr_aw_med = rnorm(n, 0.2, 0.03),
    persist_ratio = persist_ratio, persist_delta = persist_delta,
    area_ha = area_ha, doy_iqr = doy_iqr, elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE)
  geom <- sf::st_sfc(lapply(seq_len(n), function(i)
    sf::st_point(c(runif(1, 0, 1e5), runif(1, 0, 1e5)))), crs = 3035)
  sf::st_sf(df, geometry = geom)
}
# rows 1-6 (poly 2000..2005) are the 6 eligible drops; 13th (poly 2012) is a
# drop that FAILS the evidence clause (pd=-10, area=100, doy=5) -> NOT eligible.
ELIG_POLY <- as.character(2000:2005)

mk_cfg_min <- function(out_dir, enable_ah = FALSE, twr = 0.10, source_weights = NULL) {
  npp <- list(caps = c(random = 1.0, otsu = 1.0),
              artifact_hard = list(enabled = enable_ah, total_weight_ratio = twr,
                persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
                rbr_med_reference = "negative", reason_whitelist = character(0),
                persist_delta_max = -100, area_ha_min = 500, doy_iqr_max = 1))
  if (!is.null(source_weights)) npp$source_weights <- source_weights
  build_supervised_burned_config(
    run_label = "balanced", internal_decisions = file.path(out_dir, "d.gpkg"),
    change_index = file.path(out_dir, "r.tif"), target_year = 1989L,
    output_dir = out_dir, feature_whitelist_override = FEAT,
    negative_pool_params = npp, nrounds_max = 30L, early_stop = 10L)
}

# the engine's exact FINAL membership (capped), mirrored from the helpers.
engine_training_ok <- function(aug, cfg, ah_on) {
  tc <- cfg$train_control; d <- sf::st_drop_geometry(aug)
  elig <- resolve_elig(
    id = as.character(d$fire_uid), class = as.character(d$class),
    source = as.character(d$source), neg_type = as.character(d$neg_type),
    random_background_source = "random_burnable_background",
    otsu_unburned_source = "otsu_patch_residual",
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep"),
    origin_stage = "FINAL",
    artifact_hard_source = if (ah_on) "artifact_hard" else character(0))
  caps <- c(random = tc$caps$random, otsu = tc$caps$otsu)
  if ("artifact_hard" %in% names(elig$negatives_by_bucket)) caps["artifact_hard"] <- Inf
  cp <- cap_buckets(positive_idx = elig$positive_idx,
    negatives_by_bucket = elig$negatives_by_bucket, n_burned = length(elig$positive_idx),
    caps = caps, seed = tc$seeds$final_sampling_seed,
    id = as.character(d$fire_uid), context = "FINAL")
  ids <- as.character(d$fire_uid)[cp$selected_indices]
  aug[cp$selected_indices[!duplicated(ids)], , drop = FALSE]
}

mk_dir <- function(tag) {
  od <- file.path(tempdir(), paste0(tag, as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE); od
}

# ---------------------------------------------------------------------------
# (1) + (9) eligible set is identical OFF vs ON; no candidate disappears OFF.
# ---------------------------------------------------------------------------
test_that("artifact_hard_eligible is the SAME set with enabled FALSE and TRUE", {
  od <- mk_dir("ah_elig_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE, twr = 0.10)
  cfg_on  <- mk_cfg_min(od, enable_ah = TRUE,  twr = 0.10)
  promo <- promote_artifact_hard_negatives(tb, sc, cfg_on)
  aug   <- sf::st_as_sf(promo$train_features)

  Loff <- sf::st_drop_geometry(build_pool(tb,  sc, config = cfg_off, target_year = 1989L))
  Lon  <- sf::st_drop_geometry(build_pool(aug, sc, config = cfg_on,  target_year = 1989L))

  expect_setequal(Loff$poly_id[Loff$artifact_hard_eligible], ELIG_POLY)
  expect_setequal(Lon$poly_id[Lon$artifact_hard_eligible],  ELIG_POLY)
  # OFF: candidates ALL present and NONE entered training.
  expect_equal(sum(Loff$artifact_hard_eligible), length(ELIG_POLY))
  expect_false(any(Loff$artifact_hard_used))
  expect_true(all(!Loff$used_in_training[Loff$artifact_hard_eligible]))
})

# ---------------------------------------------------------------------------
# (3) + (10) OFF: weight_ratio recorded (=0.10) but sample_weight NA for unused
#            candidates; 0.10 also lives in the resolved config.
# ---------------------------------------------------------------------------
test_that("OFF: configured weight_ratio shown, sample_weight NA for candidates", {
  od <- mk_dir("ah_offw_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE, twr = 0.10)
  expect_equal(cfg_off$negative_pool_params$artifact_hard$total_weight_ratio, 0.10)
  L <- sf::st_drop_geometry(build_pool(tb, sc, config = cfg_off, target_year = 1989L))
  cand <- L$artifact_hard_eligible
  expect_true(all(L$artifact_hard_weight_ratio[cand] == 0.10))
  expect_true(all(is.na(L$sample_weight[cand])))
  expect_true(all(!L$artifact_hard_enabled))
  expect_true(all(L$artifact_hard_branch[cand] %in%
                  c("persist_delta", "large_single_doy", "reason_whitelist", "multiple")))
  expect_true(all(is.finite(L$artifact_hard_rbr_threshold[cand])))
  # candidate-only rows keep their deterministic origin.
  expect_true(all(L$original_pool_source[cand] == "deterministic_drop"))
  expect_true(all(L$pool_source[cand] == "deterministic_drop"))
})

# ---------------------------------------------------------------------------
# (2) + (8) OFF: the candidates do NOT perturb the baseline training rows / the
#           training membership computed WITHOUT scoring.
# ---------------------------------------------------------------------------
test_that("OFF: baseline training rows + membership are byte-identical with/without scoring", {
  od <- mk_dir("ah_par_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE, twr = 0.10)
  L_no <- sf::st_drop_geometry(build_pool(tb, scoring_features = NULL,
                                          config = cfg_off, target_year = 1989L))
  L_sc <- sf::st_drop_geometry(build_pool(tb, scoring_features = sc,
                                          config = cfg_off, target_year = 1989L))
  # restrict the scoring-aware layer to the baseline pool rows (drop appended cands)
  L_sc_base <- L_sc[L_sc$fire_uid %in% L_no$fire_uid, ]
  L_sc_base <- L_sc_base[match(L_no$fire_uid, L_sc_base$fire_uid), ]
  for (cl in c("used_in_training", "sample_weight", "pool_source",
               "source_total_weight")) {
    expect_equal(L_sc_base[[cl]], L_no[[cl]], info = cl)
  }
  # the engine membership is unchanged by the presence of candidates.
  tok <- engine_training_ok(tb, cfg_off, ah_on = FALSE)
  expect_setequal(L_sc$fire_uid[L_sc$used_in_training], as.character(tok$fire_uid))
})

# ---------------------------------------------------------------------------
# (4) + (6) + (7) ON: eligibles promoted -> used; membership == DMatrix; weights.
# ---------------------------------------------------------------------------
test_that("ON: eligibles are used; used_in_training == engine frame; weights match", {
  od <- mk_dir("ah_on_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_on <- mk_cfg_min(od, enable_ah = TRUE, twr = 0.10)
  promo <- promote_artifact_hard_negatives(tb, sc, cfg_on)
  aug   <- sf::st_as_sf(promo$train_features)
  tok   <- engine_training_ok(aug, cfg_on, ah_on = TRUE)
  L <- sf::st_drop_geometry(build_pool(aug, sc, training_ok = tok, config = cfg_on,
                                       target_year = 1989L))
  elig <- L$artifact_hard_eligible
  expect_true(all(L$artifact_hard_used[elig]))
  expect_true(all(L$pool_source[elig] == "artifact_hard"))
  expect_true(all(L$training_label[elig] == "unburned"))
  expect_true(all(L$original_pool_source[elig] == "deterministic_drop"))
  # used_in_training matches the engine frame exactly.
  expect_setequal(L$fire_uid[L$used_in_training], as.character(tok$fire_uid))
  # sample_weight matches the engine-resolved weights per id.
  tdf <- sf::st_drop_geometry(tok)
  sw <- resolve_w(class = as.character(tdf$class), source = as.character(tdf$source),
                  source_weights = NULL, artifact_hard_source = "artifact_hard",
                  total_weight_ratio = 0.10)
  w_by_id <- stats::setNames(sw$weights, as.character(tdf$fire_uid))
  used <- L$used_in_training
  expect_equal(L$sample_weight[used], as.numeric(w_by_id[L$fire_uid[used]]))
})

# ---------------------------------------------------------------------------
# (5) + ratio sweep: sum(weight[artifact_hard]) == ratio * sum(weight[burned]);
#     changing the ratio does NOT change the eligible set (req 7).
# ---------------------------------------------------------------------------
test_that("ON: artifact_hard total weight == ratio * burned total, for 0.10/0.25/0.50", {
  od <- mk_dir("ah_ratio_")
  tb <- mk_train_base(); sc <- mk_scoring()
  elig_sets <- list()
  for (twr in c(0.10, 0.25, 0.50)) {
    cfg_on <- mk_cfg_min(od, enable_ah = TRUE, twr = twr)
    promo  <- promote_artifact_hard_negatives(tb, sc, cfg_on)
    aug    <- sf::st_as_sf(promo$train_features)
    tok    <- engine_training_ok(aug, cfg_on, ah_on = TRUE)
    L <- sf::st_drop_geometry(build_pool(aug, sc, training_ok = tok, config = cfg_on,
                                         target_year = 1989L))
    u <- L$used_in_training
    ah_tot  <- sum(L$sample_weight[u & L$pool_source == "artifact_hard"])
    pos_tot <- sum(L$sample_weight[u & L$pool_source == "high_confidence_keep"])
    expect_equal(ah_tot, twr * pos_tot)
    expect_true(all(L$artifact_hard_weight_ratio[L$artifact_hard_eligible] == twr))
    elig_sets[[as.character(twr)]] <- sort(L$poly_id[L$artifact_hard_eligible])
  }
  # the eligible set is identical across ratios (ratio only changes the weight).
  expect_equal(elig_sets[["0.1"]], elig_sets[["0.25"]])
  expect_equal(elig_sets[["0.1"]], elig_sets[["0.5"]])
})

# ---------------------------------------------------------------------------
# Schema hygiene: no geometry / list-columns among the descriptive features.
# ---------------------------------------------------------------------------
test_that("the consolidated layer carries no geometry/list-columns as features", {
  od <- mk_dir("ah_hy_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE)
  layer <- build_pool(tb, sc, config = cfg_off, target_year = 1989L)
  L <- sf::st_drop_geometry(layer)
  is_listcol <- vapply(L, function(z) is.list(z) || inherits(z, "sfc"), logical(1))
  expect_false(any(is_listcol))
  expect_setequal(names(L), ns$.of_training_pool_layer_cols())
})

# ---------------------------------------------------------------------------
# review_priority / review_tier / label_confidence / VISUAL.
# ---------------------------------------------------------------------------
mk_scoring2 <- function() {
  # two eligible drops: a BRIGHT+LARGE+more-persistent one (high priority, low
  # label confidence) and a DIM+SMALL+transient one (low priority, high conf).
  df <- data.frame(
    fire_uid = c("SC_a", "SC_b"), poly_id = c("9001", "9002"),
    class = NA_character_, source = NA_character_, neg_type = NA_character_,
    class_final = c("drop", "drop"), reason_1 = NA_character_,
    fold_rep1 = NA_integer_, fold_rep2 = NA_integer_,
    rbr_med = c(0.55, 0.20), rbr_aw_med = c(0.2, 0.2),
    persist_ratio = c(0.30, 0.15), persist_delta = c(-200, -200),
    area_ha = c(5000, 600), doy_iqr = c(0, 0), elev_med = c(1000, 1000),
    stringsAsFactors = FALSE)
  geom <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 3035)
  sf::st_sf(df, geometry = geom)
}

test_that("review_priority orders fire-like candidates first; label_confidence inverts it", {
  od <- mk_dir("ah_pri_")
  tb <- mk_train_base(); sc <- mk_scoring2()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE, twr = 0.10)
  L <- sf::st_drop_geometry(build_pool(tb, sc, config = cfg_off, target_year = 1989L))
  cand <- L$artifact_hard_eligible
  expect_equal(sum(cand), 2L)
  # priority is in [0,100] for candidates, NA otherwise.
  expect_true(all(L$review_priority[cand] >= 0 & L$review_priority[cand] <= 100))
  expect_true(all(is.na(L$review_priority[!cand])))
  pr <- stats::setNames(L$review_priority[cand], L$poly_id[cand])
  expect_gt(pr[["9001"]], pr[["9002"]])                 # bright+large first
  # label_confidence = 1 - priority/100 for candidates (inverse of fire-risk).
  lc <- stats::setNames(L$label_confidence[cand], L$poly_id[cand])
  expect_equal(unname(lc[["9001"]]), round(1 - pr[["9001"]] / 100, 3))
  expect_lt(lc[["9001"]], lc[["9002"]])                 # fire-like => less trustworthy
  # review_tier bands are consistent with the score.
  expect_true(all((L$review_tier[cand] == "high")   == (L$review_priority[cand] >= 66)))
  expect_true(all((L$review_tier[cand] == "low")    == (L$review_priority[cand] < 33)))
})

test_that("label_confidence uses the per-source base; VISUAL is a blank NA column", {
  od <- mk_dir("ah_conf_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE)
  L <- sf::st_drop_geometry(build_pool(tb, sc, config = cfg_off, target_year = 1989L))
  expect_true(all(L$label_confidence[L$pool_source == "high_confidence_keep"] == 0.90))
  expect_true(all(L$label_confidence[L$pool_source == "random"] == 0.85))
  expect_true(all(L$label_confidence[L$pool_source == "otsu"] == 0.85))
  # VISUAL is present, all-NA, ready for the user to fill (1 = confirm, 0 = not).
  expect_true("VISUAL" %in% names(L))
  expect_true(all(is.na(L$VISUAL)))
})

# ---------------------------------------------------------------------------
# (C) rows that do NOT satisfy the rule keep their normal behaviour.
# ---------------------------------------------------------------------------
test_that("non-eligible rows are not treated as artifact_hard", {
  od <- mk_dir("ah_ne_")
  tb <- mk_train_base(); sc <- mk_scoring()
  cfg_off <- mk_cfg_min(od, enable_ah = FALSE)
  L <- sf::st_drop_geometry(build_pool(tb, sc, config = cfg_off, target_year = 1989L))
  # the baseline pool rows are never eligible / used / weighted as artifact_hard.
  base <- L$pool_source %in% c("high_confidence_keep", "random", "otsu")
  expect_false(any(L$artifact_hard_eligible[base]))
  expect_false(any(L$artifact_hard_used[base]))
  expect_true(all(is.na(L$artifact_hard_weight_ratio[base])))
})
