# =============================================================================
# PHASE 2 (artifact_hard): BASELINE-REPRODUCTION + additive/off-by-default
# contracts. Covers spec G items 4, 5, 9, 10.
# =============================================================================

ns <- asNamespace("OtsuFire")
resolve_np <- ns$.of_resolve_negative_pool_params
valid_buckets <- ns$.of_valid_negative_buckets
resolve_elig <- ns$.of_resolve_supervised_eligibility
nested_fit <- ns$.of_nested_refit_fit

# The pre-Phase-2 (0.11.0) resolved values for the EXISTING keys. These are the
# hardcoded defaults that must remain byte-identical when artifact_hard is
# disabled and source_weights is NULL.
EXPECTED_0_11_0 <- list(
  random = list(n_cells = 1500L, rbr_quantile = 0.50),
  otsu   = list(candidate_threshold = 0, reference_threshold = 100),
  caps   = c(random = 1.0, otsu = 1.0)
)

test_that("disabled + NULL weights: existing resolved keys == 0.11.0 byte-for-byte", {
  v <- resolve_np(list())$values
  expect_identical(v$random, EXPECTED_0_11_0$random)
  expect_identical(v$otsu,   EXPECTED_0_11_0$otsu)
  expect_identical(v$caps,   EXPECTED_0_11_0$caps)
})

test_that("artifact_hard defaults OFF; source_weights defaults NULL", {
  v <- resolve_np(list())$values
  expect_false(v$artifact_hard$enabled)
  expect_null(v$source_weights)
  # the off-default block has the spec's documented thresholds
  expect_equal(v$artifact_hard$total_weight_ratio, 0.10)
  expect_equal(v$artifact_hard$persist_ratio_max, 0.35)
  expect_identical(v$artifact_hard$rbr_med_reference, "negative")
  expect_equal(v$artifact_hard$rbr_med_min_q, 0.90)
  expect_identical(v$artifact_hard$reason_whitelist, character(0))
  expect_equal(v$artifact_hard$persist_delta_max, -100)
  expect_equal(v$artifact_hard$area_ha_min, 500)
  expect_equal(v$artifact_hard$doy_iqr_max, 1)
})

test_that("unknown artifact_hard key still ERRORS (unknown-key validation kept)", {
  expect_error(resolve_np(list(artifact_hard = list(bogus = 1))),
               "Unknown negative_pool_params\\$artifact_hard")
})

test_that("malformed source_weights ERRORS", {
  expect_error(resolve_np(list(source_weights = c(1, 2))),    # unnamed
               "NAMED numeric")
  expect_error(resolve_np(list(source_weights = c(a = -1))),  # negative
               "finite and")
})

test_that("source_weights / artifact_hard overrides are accepted + recorded", {
  r <- resolve_np(list(
    artifact_hard = list(enabled = TRUE, area_ha_min = 200),
    source_weights = c(artifact_hard = 0.5)
  ))
  expect_true(r$values$artifact_hard$enabled)
  expect_equal(r$values$artifact_hard$area_ha_min, 200)
  expect_equal(r$values$source_weights[["artifact_hard"]], 0.5)
  expect_equal(r$provenance$artifact_hard$enabled, "user")
  expect_equal(r$provenance$source_weights, "user")
})

# ---- bucket enum default unchanged (spec B) -------------------------------
test_that("valid negative buckets default stays the canonical two", {
  expect_identical(valid_buckets(), c("random", "otsu"))
  expect_identical(valid_buckets(include_artifact_hard = FALSE), c("random", "otsu"))
  expect_identical(valid_buckets(include_artifact_hard = TRUE),
                   c("random", "otsu", "artifact_hard"))
})

test_that("eligibility with default artifact_hard_source is unchanged (no 3rd bucket)", {
  e <- resolve_elig(
    id = as.character(1:3),
    class = c("burned", "unburned", "unburned"),
    source = c(NA, "random_burnable_background", "otsu_patch_residual"),
    neg_type = c(NA, NA, NA),
    random_background_source = "random_burnable_background",
    otsu_unburned_source = "otsu_patch_residual",
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep")
  )
  expect_identical(names(e$negatives_by_bucket), c("random", "otsu"))
  expect_false("artifact_hard" %in% e$audit$bucket)
})

test_that("artifact_hard_source rows are carved out of the GATE-6.5 drop error", {
  # An unburned row with a non-random/non-otsu source would normally ERROR; with
  # artifact_hard_source set it resolves to the artifact_hard bucket instead.
  expect_error(
    resolve_elig(
      id = as.character(1:2), class = c("burned", "unburned"),
      source = c(NA, "artifact_hard"), neg_type = c(NA, "artifact_hard_negative"),
      random_background_source = "random_burnable_background",
      otsu_unburned_source = "otsu_patch_residual",
      otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep")
    ),
    "NO valid"
  )
  e <- resolve_elig(
    id = as.character(1:2), class = c("burned", "unburned"),
    source = c(NA, "artifact_hard"), neg_type = c(NA, "artifact_hard_negative"),
    random_background_source = "random_burnable_background",
    otsu_unburned_source = "otsu_patch_residual",
    otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep"),
    artifact_hard_source = "artifact_hard"
  )
  expect_equal(length(e$negatives_by_bucket$artifact_hard), 1L)
})

# ---- DMatrix columns unchanged + reproducibility (spec G 4, 5) -------------
mk_train <- function(n = 80, seed = 7) {
  set.seed(seed)
  data.frame(
    class = c(rep("burned", n / 2), rep("unburned", n / 2)),
    block_id = sample(letters[1:8], n, TRUE),
    rbr_med = c(rnorm(n / 2, 0.6, 0.1), rnorm(n / 2, 0.2, 0.1)),
    elev_med = rnorm(n, 1000, 200),
    stringsAsFactors = FALSE
  )
}
feat <- c("rbr_med", "elev_med")

test_that("promotion/weights do not change the design-matrix columns", {
  td <- mk_train()
  f0 <- nested_fit(td, feature_cols = feat, label_col = "class",
                   block_col = "block_id", val_frac = 0.2,
                   sampling_seed = 1L, fold_seed = 1L,
                   nrounds_max = 25L, early_stopping_rounds = 10L)
  y <- as.integer(td$class == "burned")
  w <- ifelse(y == 1L, 1.0, 0.5)
  fw <- nested_fit(td, feature_cols = feat, label_col = "class",
                   block_col = "block_id", val_frac = 0.2,
                   sampling_seed = 1L, fold_seed = 1L, sample_weights = w,
                   nrounds_max = 25L, early_stopping_rounds = 10L)
  expect_identical(f0$x_cols, fw$x_cols)
  expect_identical(f0$feature_cols, fw$feature_cols)
})

test_that("same seed + config -> identical model bytes (reproducibility)", {
  td <- mk_train()
  a <- nested_fit(td, feature_cols = feat, label_col = "class",
                  block_col = "block_id", val_frac = 0.2,
                  sampling_seed = 3L, fold_seed = 3L,
                  nrounds_max = 25L, early_stopping_rounds = 10L)
  b <- nested_fit(td, feature_cols = feat, label_col = "class",
                  block_col = "block_id", val_frac = 0.2,
                  sampling_seed = 3L, fold_seed = 3L,
                  nrounds_max = 25L, early_stopping_rounds = 10L)
  expect_identical(xgboost::xgb.save.raw(a$model), xgboost::xgb.save.raw(b$model))
})
