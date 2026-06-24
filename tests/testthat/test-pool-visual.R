# Tests for apply_visual_validation() (Phase 2A): consuming a visually-validated
# supervised training pool via the VISUAL column.

mk_pool <- function(VISUAL, eligible, poly_id = seq_along(VISUAL)) {
  data.frame(poly_id = poly_id, artifact_hard_eligible = eligible,
             VISUAL = VISUAL, stringsAsFactors = FALSE)
}

test_that("1. pool without a VISUAL column errors (not a validated pool)", {
  p <- data.frame(poly_id = 1:3, artifact_hard_eligible = c(FALSE, FALSE, TRUE))
  expect_error(apply_visual_validation(p), "not a")
  # ...but require_visual = FALSE proceeds, treating everything as unreviewed
  vv <- apply_visual_validation(p, require_visual = FALSE)
  expect_equal(vv$pool$used_for_training, c(TRUE, TRUE, FALSE))
})

test_that("2. empty VISUAL (all NA) -> unreviewed; positives kept, candidates off", {
  p <- mk_pool(VISUAL = c(NA, NA, NA), eligible = c(FALSE, FALSE, TRUE))
  vv <- apply_visual_validation(p)
  expect_equal(as.character(vv$pool$visual_action),
               c("unreviewed_kept", "unreviewed_kept", "unreviewed_candidate_off"))
  expect_equal(vv$pool$used_for_training, c(TRUE, TRUE, FALSE))
  expect_equal(nrow(vv$omissions), 0L)
  expect_equal(nrow(vv$dropped), 0L)
})

test_that("3. visually accepted rows (VISUAL=1): positive confirmed, candidate promoted", {
  p <- mk_pool(VISUAL = c(1, 1), eligible = c(FALSE, TRUE))
  vv <- apply_visual_validation(p)
  expect_equal(as.character(vv$pool$visual_action), c("confirmed", "promoted_artifact_hard"))
  expect_true(all(vv$pool$used_for_training))
  expect_equal(nrow(vv$training), 2L)
})

test_that("4. visually rejected rows (VISUAL=0): positive dropped, candidate omitted", {
  p <- mk_pool(VISUAL = c(0, 0), eligible = c(FALSE, TRUE))
  vv <- apply_visual_validation(p)
  expect_equal(as.character(vv$pool$visual_action),
               c("dropped_visual_reject", "omitted_real_fire"))
  expect_false(any(vv$pool$used_for_training))
  expect_equal(nrow(vv$dropped), 1L)    # the positive/background
  expect_equal(nrow(vv$omissions), 1L)  # the candidate (real fire among drops)
})

test_that("5. invalid VISUAL values are rejected (error) or excluded", {
  p <- mk_pool(VISUAL = c(1, 2, "x"), eligible = c(FALSE, TRUE, FALSE))
  expect_error(apply_visual_validation(p), "invalid")
  vv <- apply_visual_validation(p, on_invalid = "exclude")
  expect_equal(as.character(vv$pool$visual_action[2:3]),
               c("invalid_excluded", "invalid_excluded"))
  expect_false(any(vv$pool$used_for_training[2:3]))
  expect_true(vv$pool$used_for_training[1])  # the valid confirmed row still trains
})

test_that("6. discarded rows never enter the clean training pool", {
  p <- mk_pool(VISUAL = c(NA, 0, 1, 0, 1), eligible = c(FALSE, FALSE, TRUE, TRUE, FALSE))
  vv <- apply_visual_validation(p)
  expect_true(all(vv$training$used_for_training))
  expect_equal(nrow(vv$training), sum(vv$pool$used_for_training))
  discarded <- vv$pool$poly_id[!vv$pool$used_for_training]
  expect_false(any(discarded %in% vv$training$poly_id))
})

test_that("7. a traceability column (used_for_training) is produced", {
  p <- mk_pool(VISUAL = c(1, 0), eligible = c(FALSE, TRUE))
  vv <- apply_visual_validation(p)
  expect_true(all(c("visual_value", "visual_action", "used_for_training") %in% names(vv$pool)))
  expect_type(vv$pool$used_for_training, "logical")
  expect_s3_class(vv$pool$visual_action, "factor")
})

test_that("8. sf input keeps geometry in the clean training pool", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 3035)
  p <- sf::st_sf(poly_id = 1:2, artifact_hard_eligible = c(FALSE, TRUE),
                 VISUAL = c(1L, 0L), geometry = pts)
  vv <- apply_visual_validation(p)
  expect_s3_class(vv$training, "sf")
  expect_equal(nrow(vv$training), 1L)  # only the confirmed positive survives
})
