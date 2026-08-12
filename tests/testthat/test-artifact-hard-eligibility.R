# =============================================================================
# PHASE 2 (artifact_hard): eligibility rule (spec C) + selection contract.
#
# Covers spec G items 1 (each clause; NA-handling; reason whitelist vs
# OR-evidence) and 2 (artifact_uncertain returned but NEVER trained-labelled).
# =============================================================================

ns <- asNamespace("OtsuFire")
select_ah <- ns$.of_select_artifact_hard

base_params <- list(
  enabled            = TRUE,
  persist_ratio_max  = 0.35,
  rbr_med_min_q      = 0.25,
  reason_whitelist   = character(0),
  persist_delta_max  = -100,
  area_ha_min        = 500,
  doy_iqr_max        = 1
)

# A single-row scoring frame factory so each clause can be isolated.
mk <- function(class_final = "drop", persist_ratio = 0.1, rbr_med = 0.9,
               persist_delta = -200, area_ha = 1000, doy_iqr = 0,
               reason_1 = "foo", rbr_aw_med = 0.5, rbr_iqr = 0.1) {
  data.frame(
    fire_uid = "x1", poly_id = "x1",
    class_final = class_final, persist_ratio = persist_ratio, rbr_med = rbr_med,
    persist_delta = persist_delta, area_ha = area_ha, doy_iqr = doy_iqr,
    reason_1 = reason_1, rbr_aw_med = rbr_aw_med, rbr_iqr = rbr_iqr,
    stringsAsFactors = FALSE
  )
}
# reliable positives whose Q25 is ~0.4 (so rbr_med 0.9 passes, 0.1 fails).
rp <- seq(0.2, 0.8, length.out = 50)

test_that("a fully-eligible drop is promoted (OR via persist_delta)", {
  s <- select_ah(mk(), rp, base_params)
  expect_equal(nrow(s$artifact_hard), 1L)
  expect_equal(nrow(s$artifact_uncertain), 0L)
})

test_that("class_final != drop is never eligible", {
  s <- select_ah(mk(class_final = "keep"), rp, base_params)
  expect_equal(nrow(s$artifact_hard), 0L)
  expect_equal(nrow(s$artifact_uncertain), 0L)  # not a drop -> not even uncertain
})

test_that("NA persist_ratio or NA rbr_med fails the row (NA handling)", {
  expect_equal(nrow(select_ah(mk(persist_ratio = NA), rp, base_params)$artifact_hard), 0L)
  expect_equal(nrow(select_ah(mk(rbr_med = NA), rp, base_params)$artifact_hard), 0L)
  # but they ARE drops, so they land in artifact_uncertain
  expect_equal(nrow(select_ah(mk(persist_ratio = NA), rp, base_params)$artifact_uncertain), 1L)
})

test_that("rbr_med below the runtime quantile threshold fails", {
  s <- select_ah(mk(rbr_med = 0.1), rp, base_params)
  expect_equal(nrow(s$artifact_hard), 0L)
  expect_equal(nrow(s$artifact_uncertain), 1L)
})

test_that("persist_ratio above persist_ratio_max fails", {
  s <- select_ah(mk(persist_ratio = 0.9), rp, base_params)
  expect_equal(nrow(s$artifact_hard), 0L)
})

test_that("OR-evidence: reason whitelist alone qualifies", {
  p <- base_params
  p$reason_whitelist <- "artifact_reason"
  # kill the other two OR branches
  s <- select_ah(mk(reason_1 = "artifact_reason", persist_delta = 0,
                    area_ha = 0, doy_iqr = 99), rp, p)
  expect_equal(nrow(s$artifact_hard), 1L)
})

test_that("OR-evidence: area & doy together qualify, but not alone", {
  # area large, doy small -> qualifies (other branches off)
  s_ok <- select_ah(mk(persist_delta = 0, area_ha = 1000, doy_iqr = 0,
                       reason_1 = "foo"), rp, base_params)
  expect_equal(nrow(s_ok$artifact_hard), 1L)
  # area large but doy too large -> fails
  s_bad <- select_ah(mk(persist_delta = 0, area_ha = 1000, doy_iqr = 5,
                        reason_1 = "foo"), rp, base_params)
  expect_equal(nrow(s_bad$artifact_hard), 0L)
})

test_that("no OR-evidence at all -> not eligible (becomes uncertain)", {
  s <- select_ah(mk(persist_delta = 0, area_ha = 0, doy_iqr = 99,
                    reason_1 = "foo"), rp, base_params)
  expect_equal(nrow(s$artifact_hard), 0L)
  expect_equal(nrow(s$artifact_uncertain), 1L)
})

test_that("promoted rows are stamped unburned / artifact_hard / neg_type", {
  s <- select_ah(mk(), rp, base_params)
  expect_identical(unique(as.character(s$artifact_hard$class)), "unburned")
  expect_identical(unique(as.character(s$artifact_hard$source)), "artifact_hard")
  expect_identical(unique(as.character(s$artifact_hard$neg_type)),
                   "artifact_hard_negative")
  # group keys preserved
  expect_true(all(c("poly_id", "fire_uid") %in% names(s$artifact_hard)))
})

test_that("derived features only attached when enable_derived = TRUE (M3 flag)", {
  s1 <- select_ah(mk(), rp, base_params, enable_derived = FALSE)
  s3 <- select_ah(mk(), rp, base_params, enable_derived = TRUE)
  derived <- c("persist_loss", "persist_abs_delta", "rbr_iqr_rel")
  expect_false(any(derived %in% names(s1$artifact_hard)))   # M1 unchanged
  expect_true(all(derived %in% names(s3$artifact_hard)))    # M3 only
})

test_that("every drop is partitioned into exactly hard XOR uncertain", {
  set.seed(11)
  n <- 60
  df <- data.frame(
    fire_uid = as.character(1:n), poly_id = as.character(1:n),
    class_final = c(rep("drop", 40), rep("keep", 20)),
    persist_ratio = runif(n, 0, 0.6), rbr_med = runif(n, 0.05, 0.95),
    persist_delta = rnorm(n, -80, 90), area_ha = runif(n, 50, 1200),
    doy_iqr = sample(0:4, n, TRUE), reason_1 = "foo",
    rbr_aw_med = runif(n, 0.05, 0.95), rbr_iqr = runif(n, 0, 0.4),
    stringsAsFactors = FALSE
  )
  s <- select_ah(df, rp, base_params)
  expect_equal(nrow(s$artifact_hard) + nrow(s$artifact_uncertain), 40L)
  # disjoint poly_ids
  expect_length(intersect(s$artifact_hard$poly_id, s$artifact_uncertain$poly_id), 0L)
})
