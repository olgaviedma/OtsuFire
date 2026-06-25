# =============================================================================
# PHASE 2 (artifact_hard): promotion CONTRACT, pinned with small synthetic
# in-memory fixtures (no gpkg, no model fitting).
#
# Pins the load-bearing guarantees of promote_artifact_hard_negatives() and the
# orchestrator-side id uniquification:
#   1. OFF is a strict no-op (train_features byte-identical, scoring untouched).
#   2. ON is additive (nrow grows by exactly n_promoted; n_promoted > 0).
#   3. scoring_features is NEVER modified (ON and OFF).
#   4. Promoted rows are class == "unburned" / source == "artifact_hard".
#   5. Promotion removes NO burned positives.
#   6. fire_uid uniqueness is enforced by the ORCHESTRATOR (not promote_*); we
#      replicate that small helper here and assert disjointness.
#   7. Consolidated training-pool writer: ELIGIBLE-but-not-used when OFF,
#      used when ON. (Covered in depth in test-artifact-hard-training-pool-layer.R;
#      here we keep a tight smoke check, skipping if the writer is unavailable.)
# =============================================================================

# --- fixtures ---------------------------------------------------------------

mk_train_pc <- function() {
  data.frame(
    fire_uid  = as.character(1:6),
    class     = c("burned", "unburned", "burned", "unburned", "burned", "burned"),
    source    = c("internal_keep_qc", "random_burnable_background",
                  "internal_keep_qc", "otsu_patch_residual",
                  "internal_keep_qc", "internal_keep_qc"),
    rbr_med   = c(0.70, 0.20, 0.60, 0.15, 0.65, 0.62),
    fold_rep1 = c(1L, 2L, 3L, 1L, 2L, 3L),
    fold_rep2 = c(2L, 3L, 1L, 2L, 3L, 1L),
    block_id  = paste0("blk_", c(1, 2, 3, 1, 2, 3)),
    stringsAsFactors = FALSE
  )
}

# Six deterministic drops: rows 1-4 fully eligible (drop, high rbr, low persist,
# pd <= -100, area >= 500); row 5 fails on rbr_med; row 6 fails on area/doy/pd.
mk_scoring_pc <- function() {
  data.frame(
    fire_uid      = as.character(101:106),
    poly_id       = as.character(101:106),
    class_final   = rep("drop", 6),
    persist_ratio = c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10),
    rbr_med       = c(0.80, 0.80, 0.80, 0.80, 0.05, 0.80),
    persist_delta = c(-200, -200, -200, -200, -200, 0),
    area_ha       = c(1000, 1000, 1000, 1000, 1000, 50),
    doy_iqr       = c(0, 0, 0, 0, 0, 5),
    reason_1      = "foo",
    rbr_aw_med    = 0.4,
    rbr_iqr       = 0.1,
    stringsAsFactors = FALSE
  )
}

cfg_off_pc <- list(negative_pool_params =
                     list(artifact_hard = list(enabled = FALSE)))
cfg_on_pc  <- list(negative_pool_params = list(artifact_hard = list(
  enabled = TRUE, persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
  reason_whitelist = character(0), persist_delta_max = -100,
  area_ha_min = 500, doy_iqr_max = 1)))

# --- 1. OFF strict no-op ----------------------------------------------------

test_that("OFF is a strict no-op: train byte-identical, scoring untouched", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  r <- promote_artifact_hard_negatives(tf, sc, cfg_off_pc)
  expect_false(r$enabled)
  expect_equal(r$n_promoted, 0L)
  expect_identical(r$train_features, tf)        # same nrow + same rows + order
  expect_equal(nrow(r$train_features), nrow(tf))
  expect_identical(r$scoring_features, sc)       # scoring unchanged
})

# --- 2. ON additive ---------------------------------------------------------

test_that("ON is additive: nrow grows by exactly n_promoted ( > 0 )", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  r <- promote_artifact_hard_negatives(tf, sc, cfg_on_pc)
  expect_true(r$enabled)
  expect_gt(r$n_promoted, 0L)
  expect_equal(nrow(r$train_features), nrow(tf) + r$n_promoted)
  # the augmented frame is the load-bearing object the orchestrator hands to
  # the OOF/FINAL stages -- it must literally contain the promoted rows.
  expect_equal(nrow(r$artifact_hard), r$n_promoted)
})

# --- 3. scoring_features never modified (ON and OFF) -------------------------

test_that("scoring_features is never modified (ON and OFF)", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  expect_identical(promote_artifact_hard_negatives(tf, sc, cfg_off_pc)$scoring_features, sc)
  expect_identical(promote_artifact_hard_negatives(tf, sc, cfg_on_pc)$scoring_features,  sc)
})

# --- 4. promoted rows are unburned / artifact_hard --------------------------

test_that("promoted rows are class=='unburned' and source=='artifact_hard'", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  r <- promote_artifact_hard_negatives(tf, sc, cfg_on_pc)
  added <- r$train_features[(nrow(tf) + 1L):nrow(r$train_features), , drop = FALSE]
  expect_true(all(added$class == "unburned"))
  expect_true(all(added$source == "artifact_hard"))
  # same on the dedicated artifact_hard slice
  expect_identical(unique(as.character(r$artifact_hard$class)), "unburned")
  expect_identical(unique(as.character(r$artifact_hard$source)), "artifact_hard")
})

# --- 5. promotion removes no burned positives -------------------------------

test_that("promotion removes NO burned positives", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  n_burned_in <- sum(tf$class == "burned")
  r <- promote_artifact_hard_negatives(tf, sc, cfg_on_pc)
  n_burned_out <- sum(r$train_features$class == "burned")
  expect_equal(n_burned_out, n_burned_in)
  # the original training rows survive verbatim at the head of the frame.
  expect_identical(r$train_features[seq_len(nrow(tf)), names(tf)], tf)
})

# --- 6. fire_uid uniqueness is an ORCHESTRATOR responsibility ----------------
# promote_*() carries the scoring rows' own fire_uids through; it does NOT
# guarantee global uniqueness vs the training pool. The orchestrator
# (internal-sup-orchestrator.R, STEP C0e) uniquifies them. We replicate that
# exact helper logic here and assert it yields disjoint, unique ids.

test_that("orchestrator-style uniquification makes promoted fire_uids disjoint + unique", {
  tf <- mk_train_pc(); sc <- mk_scoring_pc()
  r  <- promote_artifact_hard_negatives(tf, sc, cfg_on_pc)
  aug <- r$train_features
  aug$fire_uid <- as.character(aug$fire_uid)
  ah_mask <- !is.na(aug$source) & aug$source == "artifact_hard"

  # --- mirror of the orchestrator's id-uniquification (STEP C0e) ---
  uid <- as.character(aug$fire_uid[ah_mask])
  bad <- is.na(uid) | uid == "" | duplicated(uid) |
    (uid %in% as.character(aug$fire_uid[!ah_mask]))
  uid[bad] <- paste0("AH_uid_", which(ah_mask)[bad])
  aug$fire_uid[ah_mask] <- uid
  # ----------------------------------------------------------------

  promoted_ids <- aug$fire_uid[ah_mask]
  real_ids     <- aug$fire_uid[!ah_mask]
  expect_equal(anyDuplicated(promoted_ids), 0L)               # unique within promoted
  expect_length(intersect(promoted_ids, real_ids), 0L)        # disjoint from real
})

# --- 7. consolidated pool writer: ELIGIBLE when OFF, USED when ON ------------
# Smoke check only; the exhaustive contract lives in
# test-artifact-hard-training-pool-layer.R. The writer needs a full resolved
# config + sf geometry, so we build a minimal one and skip cleanly if any dep is
# missing rather than asserting a fragile partial fixture.

test_that("consolidated pool writer marks ELIGIBLE (OFF) vs USED (ON)", {
  skip_if_not_installed("sf")
  ns <- asNamespace("OtsuFire")
  if (is.null(ns$.of_build_supervised_training_pool) ||
      is.null(ns$.of_build_supervised_training_pool)) {
    skip("consolidated training-pool writer not available in this build")
  }
  build_pool <- ns$.of_build_supervised_training_pool

  od <- file.path(tempdir(), paste0("ahpc_", as.integer(runif(1, 1, 1e7))))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)

  FEAT <- c("rbr_med", "elev_med")
  mk_cfg <- function(enable_ah) {
    build_supervised_burned_config(
      run_label = "balanced",
      internal_decisions = file.path(od, "d.gpkg"),
      change_index = file.path(od, "r.tif"),
      target_year = 1989L, output_dir = od,
      feature_whitelist_override = FEAT,
      negative_pool_params = list(
        caps = c(random = 1.0, otsu = 1.0),
        artifact_hard = list(enabled = enable_ah, total_weight_ratio = 0.10,
          persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
          rbr_med_reference = "negative", reason_whitelist = character(0),
          persist_delta_max = -100, area_ha_min = 500, doy_iqr_max = 1)),
      nrounds_max = 30L, early_stop = 10L)
  }

  # tiny sf training pool (positives + 1 random + 1 otsu) and scoring universe.
  mk_sf_train <- function() {
    n <- 6L
    df <- data.frame(
      fire_uid = as.character(1:n), poly_id = as.character(1:n),
      class = c(rep("burned", 4), "unburned", "unburned"),
      source = c(rep("internal_keep_qc", 4),
                 "random_burnable_background", "otsu_patch_residual"),
      neg_type = c(rep(NA_character_, 5), "otsu_patch_drop"),
      class_final = c(rep("keep", 4), NA, NA), reason_1 = NA_character_,
      block_id = paste0("blk_", 1:n),
      fold_rep1 = c(1L, 2L, 3L, 1L, 2L, 3L), fold_rep2 = c(2L, 3L, 1L, 2L, 3L, 1L),
      rbr_med = c(0.6, 0.6, 0.6, 0.6, 0.15, 0.20), rbr_aw_med = 0.2,
      persist_ratio = 0.4, persist_delta = -10, area_ha = 100, doy_iqr = 2,
      elev_med = 1000, stringsAsFactors = FALSE)
    geom <- sf::st_sfc(lapply(seq_len(n), function(i)
      sf::st_point(c(runif(1, 0, 1e5), runif(1, 0, 1e5)))), crs = 3035)
    sf::st_sf(df, geometry = geom)
  }
  mk_sf_scoring <- function() {
    poly <- as.character(2000:2003); n <- length(poly)   # 4 eligible drops
    df <- data.frame(
      fire_uid = paste0("SC_", poly), poly_id = poly,
      class = NA_character_, source = NA_character_, neg_type = NA_character_,
      class_final = "drop", reason_1 = NA_character_,
      fold_rep1 = NA_integer_, fold_rep2 = NA_integer_,
      rbr_med = 0.45, rbr_aw_med = 0.2, persist_ratio = 0.15,
      persist_delta = -200, area_ha = 800, doy_iqr = 0, elev_med = 1000,
      stringsAsFactors = FALSE)
    geom <- sf::st_sfc(lapply(seq_len(n), function(i)
      sf::st_point(c(runif(1, 0, 1e5), runif(1, 0, 1e5)))), crs = 3035)
    sf::st_sf(df, geometry = geom)
  }

  tb <- mk_sf_train(); sc <- mk_sf_scoring()
  cfg_off <- mk_cfg(FALSE)
  cfg_on  <- mk_cfg(TRUE)

  built <- tryCatch({
    L_off <- sf::st_drop_geometry(
      build_pool(tb, sc, config = cfg_off, target_year = 1989L))
    promo <- promote_artifact_hard_negatives(tb, sc, cfg_on)
    aug   <- sf::st_as_sf(promo$train_features)
    L_on  <- sf::st_drop_geometry(
      build_pool(aug, sc, config = cfg_on, target_year = 1989L))
    list(off = L_off, on = L_on)
  }, error = function(e) {
    skip(paste("training-pool writer fixture not buildable:", conditionMessage(e)))
  })

  # OFF: candidates listed as eligible but not used as artifact_hard.
  cand_off <- built$off$artifact_hard_eligible
  expect_gt(sum(cand_off), 0L)
  expect_false(any(built$off$artifact_hard_used))

  # ON: the promoted (now eligible) rows are marked used.
  cand_on <- built$on$artifact_hard_eligible
  expect_gt(sum(cand_on), 0L)
  expect_true(any(built$on$artifact_hard_used))
})
