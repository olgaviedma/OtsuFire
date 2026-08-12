# =============================================================================
# PHASE 2 (artifact_hard): promote_artifact_hard_negatives() wiring contract.
#
# Covers the OFF-by-default no-op (byte-identical train_features) and the
# enabled augmentation (rows added to train, scoring unchanged), plus the
# guarantee that artifact_uncertain is never folded into training (spec G 2).
# =============================================================================

mk_train <- function() {
  data.frame(
    fire_uid = as.character(1:5),
    class = c("burned", "unburned", "burned", "unburned", "burned"),
    rbr_med = c(0.7, 0.2, 0.6, 0.15, 0.65),
    stringsAsFactors = FALSE
  )
}
mk_scoring <- function() {
  data.frame(
    fire_uid = as.character(101:106), poly_id = as.character(101:106),
    class_final = rep("drop", 6),
    persist_ratio = c(0.1, 0.1, 0.1, 0.9, 0.1, 0.1),
    rbr_med = c(0.8, 0.8, 0.8, 0.8, 0.05, 0.8),
    persist_delta = c(-200, -200, 0, -200, -200, 0),
    area_ha = c(1000, 1000, 1000, 1000, 1000, 50),
    doy_iqr = c(0, 0, 5, 0, 0, 0),
    reason_1 = "foo",
    rbr_aw_med = 0.4, rbr_iqr = 0.1,
    stringsAsFactors = FALSE
  )
}

cfg_off <- list(negative_pool_params = list(artifact_hard = list(enabled = FALSE)))
cfg_on  <- list(negative_pool_params = list(artifact_hard = list(
  enabled = TRUE, persist_ratio_max = 0.35, rbr_med_min_q = 0.25,
  reason_whitelist = character(0), persist_delta_max = -100,
  area_ha_min = 500, doy_iqr_max = 1)))

test_that("disabled: train_features returned byte-identical (strict no-op)", {
  tf <- mk_train(); sf_ <- mk_scoring()
  r <- promote_artifact_hard_negatives(tf, sf_, cfg_off)
  expect_identical(r$train_features, tf)
  expect_identical(r$scoring_features, sf_)
  expect_equal(r$n_promoted, 0L)
  expect_false(r$enabled)
})

test_that("enabled: eligible rows are appended to train, scoring unchanged", {
  tf <- mk_train(); sf_ <- mk_scoring()
  r <- promote_artifact_hard_negatives(tf, sf_, cfg_on)
  expect_true(r$enabled)
  expect_gt(r$n_promoted, 0L)
  expect_equal(nrow(r$train_features), nrow(tf) + r$n_promoted)
  # scoring universe is NOT shrunk (rows stay)
  expect_identical(r$scoring_features, sf_)
  # promoted rows carry the artifact_hard stamp
  added <- r$train_features[(nrow(tf) + 1L):nrow(r$train_features), ]
  expect_true(all(added$class == "unburned"))
  expect_true(all(added$source == "artifact_hard"))
})

test_that("artifact_uncertain rows are NEVER in the augmented training set", {
  tf <- mk_train(); sf_ <- mk_scoring()
  r <- promote_artifact_hard_negatives(tf, sf_, cfg_on)
  unc_ids <- as.character(r$artifact_uncertain$fire_uid)
  expect_length(intersect(unc_ids, as.character(r$train_features$fire_uid)), 0L)
})

test_that("enabled but scoring NULL errors", {
  expect_error(promote_artifact_hard_negatives(mk_train(), NULL, cfg_on),
               "scoring_features.*required")
})
