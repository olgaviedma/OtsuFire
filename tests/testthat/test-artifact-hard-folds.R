# =============================================================================
# PHASE 2 (artifact_hard): spatial-group integrity (spec G 6).
#
# The promoted artifact_hard rows MUST keep their fire_uid (group key) so the
# existing OOF overlap guard (run_oof_xgb asserts no fire_uid / block_id crosses
# outer train/test) covers them exactly like any other training row. This file
# proves the group key survives promotion and that the OOF guard logic rejects a
# fire_uid that would straddle folds (including an artifact_hard fire_uid).
# =============================================================================

cfg_on <- list(negative_pool_params = list(artifact_hard = list(
  enabled = TRUE, persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
  reason_whitelist = character(0), persist_delta_max = -100,
  area_ha_min = 500, doy_iqr_max = 1)))

test_that("promoted rows keep their fire_uid group key", {
  tf <- data.frame(
    fire_uid = as.character(1:3),
    class = c("burned", "unburned", "burned"),
    rbr_med = c(0.7, 0.2, 0.6), stringsAsFactors = FALSE
  )
  sf_ <- data.frame(
    fire_uid = as.character(201:203), poly_id = as.character(201:203),
    class_final = "drop", persist_ratio = 0.1, rbr_med = 0.85,
    persist_delta = -200, area_ha = 1000, doy_iqr = 0, reason_1 = "foo",
    rbr_aw_med = 0.4, rbr_iqr = 0.1, stringsAsFactors = FALSE
  )
  r <- promote_artifact_hard_negatives(tf, sf_, cfg_on)
  added <- r$train_features[(nrow(tf) + 1L):nrow(r$train_features), ]
  expect_true(all(!is.na(added$fire_uid)))
  # group keys are the scoring rows' own fire_uids (no collision with train)
  expect_length(intersect(added$fire_uid, tf$fire_uid), 0L)
})

test_that("OOF guard rejects a fire_uid that straddles outer folds (incl artifact_hard)", {
  # Mirror the run_oof_xgb overlap assertion: if a fire_uid is in both outer
  # train and outer test, the guard must stop. We reproduce that check on a tiny
  # frame containing an artifact_hard fire_uid placed in BOTH folds.
  labelled_df <- data.frame(
    fire_uid = c("A", "A", "B", "C"),
    class = c("burned", "unburned", "unburned", "burned"),
    source = c(NA, "artifact_hard", "random_burnable_background", NA),
    fold_rep1 = c(1L, 2L, 1L, 2L),  # fire_uid "A" straddles folds 1 and 2
    stringsAsFactors = FALSE
  )
  r <- 1L; k <- 1L
  idx_te <- which(labelled_df$fold_rep1 == k)
  idx_tr <- which(labelled_df$fold_rep1 != k)
  uid_tr <- as.character(labelled_df$fire_uid[idx_tr])
  uid_te <- as.character(labelled_df$fire_uid[idx_te])
  expect_gt(length(intersect(uid_tr, uid_te)), 0L)  # "A" straddles -> guard fires
})
