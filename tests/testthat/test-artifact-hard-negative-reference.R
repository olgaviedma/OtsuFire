# PHASE 2: artifact_hard RBR gate referenced to the EXISTING NEGATIVE pool
# (rbr_med_reference = "negative", the Q90 floor). Self-contained / synthetic --
# no external data files. Verifies: Q90-of-negatives computation, sentinel-like
# capture, trivial-RBR exclusion, explicit NA handling, no-op when disabled, and
# that the scored universe is never reduced.

mk_cfg_ah <- function(reference = "negative", q = 0.90, enabled = TRUE) {
  list(negative_pool_params = list(artifact_hard = list(
    enabled = enabled, persist_ratio_max = 0.35, rbr_med_reference = reference,
    rbr_med_min_q = q, reason_whitelist = character(0), persist_delta_max = -100,
    area_ha_min = 500, doy_iqr_max = 1)))
}

# Training pool: burned positives (rbr_med ~ 440-520) + random/otsu negatives
# (rbr_med ~ 80-230). Negatives' Q90 is the gate floor.
mk_train_pool <- function() {
  data.frame(
    class   = c(rep("burned", 50), rep("unburned", 100)),
    source  = c(rep("internal_keep_qc", 50),
                rep("random_burnable_background", 50), rep("otsu_patch_residual", 50)),
    rbr_med = c(seq(440, 520, length.out = 50),
                seq(80, 180, length.out = 50), seq(120, 230, length.out = 50)),
    stringsAsFactors = FALSE)
}

# Scoring universe: drops (sentinel-like / trivial-RBR / NA), plus keep/review.
mk_scoring_pool <- function() {
  data.frame(
    poly_id       = as.character(1:9),
    class_final   = c("drop","drop","drop","drop","drop","drop","keep","review","drop"),
    rbr_med       = c(310,  320,  150,  NA,   330,  305,  500,  480,  600),
    persist_ratio = c(0.20, 0.10, 0.15, 0.10, NA,   0.30, 0.60, 0.50, 0.10),
    persist_delta = c(-200,-150, -50, -200,-200,-120, -10, -20, -300),
    area_ha       = c(13000,600,  700,  800,  900,  50,   100,  100,  2000),
    doy_iqr       = c(0,    1,    1,    1,    0,    5,    0,    0,    0),
    reason_1      = NA_character_, stringsAsFactors = FALSE)
}

test_that("RBR gate uses the Q90 of the EXISTING negative pool", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  neg_q90 <- as.numeric(stats::quantile(tr$rbr_med[tr$class == "unburned"], 0.90, names = FALSE))
  p <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90))
  # threshold echoed by the selector equals the negatives' Q90 (NOT the positives').
  expect_equal(p$audit$n[p$audit$clause == "rbr_med >= q(rbr_med_min_q)"],
               sum(sc$class_final == "drop" & is.finite(sc$rbr_med) & sc$rbr_med >= neg_q90))
  expect_true(neg_q90 < min(tr$rbr_med[tr$class == "burned"]))  # below fire level
})

test_that("sentinel-like drops are captured; trivial-RBR and NA drops are excluded", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  p <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90))
  ah_ids  <- as.character(sf::st_drop_geometry(p$artifact_hard)$poly_id)
  unc_ids <- as.character(sf::st_drop_geometry(p$artifact_uncertain)$poly_id)
  # eligible: 1,2,6,9 (drop & rbr>=Q90 & persist<=0.35 & (pd<=-100 | area&doy))
  expect_setequal(ah_ids, c("1", "2", "6", "9"))
  expect_equal(p$n_promoted, 4L)
  # trivial RBR (row 3, rbr 150 < Q90) excluded -> uncertain, never promoted.
  expect_false("3" %in% ah_ids); expect_true("3" %in% unc_ids)
  # NA rbr_med (row 4) and NA persist_ratio (row 5) excluded (explicit NA-fail).
  expect_false(any(c("4", "5") %in% ah_ids))
  expect_true(all(c("4", "5") %in% unc_ids))
  # keep/review are never artifact_hard candidates.
  expect_false(any(c("7", "8") %in% c(ah_ids, unc_ids)))
})

test_that("promoted rows are stamped as unburned artifact_hard negatives", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  p <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90))
  ah <- sf::st_drop_geometry(p$artifact_hard)
  expect_true(all(ah$class == "unburned"))
  expect_true(all(ah$source == "artifact_hard"))
  expect_true(all(ah$neg_type == "artifact_hard_negative"))
})

test_that("positive reference is stricter than negative (fewer / no captures here)", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  p_neg <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90))
  p_pos <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("positive", 0.25))
  # positives' Q25 (~450) is above the sentinel-like RBR (~310) -> excludes them.
  expect_gt(p_neg$n_promoted, p_pos$n_promoted)
})

test_that("strict no-op when artifact_hard is disabled", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  p <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90, enabled = FALSE))
  expect_equal(p$n_promoted, 0L)
  expect_identical(p$train_features, tr)   # byte-identical training frame
})

test_that("scored universe is never reduced by promotion", {
  tr <- mk_train_pool(); sc <- mk_scoring_pool()
  p <- promote_artifact_hard_negatives(tr, sc, mk_cfg_ah("negative", 0.90))
  expect_equal(nrow(p$scoring_features), nrow(sc))   # scoring returned unchanged
  expect_gt(nrow(p$train_features), nrow(tr))        # training grew by the promoted rows
  # promotion adds rows, not feature columns -> OOF and FINAL see the same schema.
  expect_true(all(names(tr) %in% names(p$train_features)))
})
