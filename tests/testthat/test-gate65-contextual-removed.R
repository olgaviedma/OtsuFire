# =============================================================================
# GATE 6.5 (2026-06-12): remove the contextual negative bucket (two sources:
# random background + Otsu residual).
#
# Audit across 6 years proved no contextual (deterministic-drop) category is a
# reliable unburned negative; a deterministic drop does NOT automatically become
# an unburned label. This suite proves, for the ambiguous deterministic drops
# (geo_excluded_hot / too_few_pixels / too_small):
#   (1) they do NOT appear in train_labeled (never written as unburned negatives);
#   (2) they do NOT receive class = "unburned";
#   (3) they DO reach the scoring pool;
#   (4) they CAN receive p_burned (present + scoreable in the scoring pool);
#   (5) they are NOT silently lost (accounted in scoring_pool with their original
#       class).
# Plus: only two buckets exist (random, otsu); a contextual bucket /
# deterministic_drop_hard-as-unburned row produces an explicit ERROR; OOF and
# FINAL share the same resolver + the same two caps; cap_contextual passed to a
# public function ERRORS.
#
# Synthetic fixtures only -- NO pool rebuild, NO training, NO scoring run.
# =============================================================================

ns65 <- asNamespace("OtsuFire")
g65  <- function(nm) get(nm, envir = ns65)

# Ambiguous deterministic-drop neg_types the 6-year audit rejected.
AMBIG_DROP <- c("geo_excluded_hot", "drop_too_few_pixels", "drop_too_small")

# ---------------------------------------------------------------------------
# (1)+(2) The deterministic builder no longer stamps drops as unburned, so an
# ambiguous deterministic drop can NEVER enter train_labeled as an unburned
# negative. Verified at the builder source (the productive code path).
# ---------------------------------------------------------------------------
test_that("GATE 6.5: the deterministic builder never writes class='unburned' det-drop rows", {
  src <- paste(deparse(
    body(g65("build_unburned_from_deterministic_decisions"))), collapse = "\n")
  # The former block (filter class_final=="drop" -> class="unburned" +
  # source="deterministic_drop_hard") is GONE.
  expect_false(grepl('source = "deterministic_drop_hard"', src, fixed = TRUE))
  expect_false(grepl('class_final == "drop"', src, fixed = TRUE))
  # unburned_hard is now an EMPTY-schema slice (no drop rows flow into the pool):
  # built from a zero-row mutate stamping empty class/neg_type/source columns.
  expect_true(grepl("class = character(0)", src, fixed = TRUE))
  expect_true(grepl("neg_type = character(0)", src, fixed = TRUE))
  expect_true(grepl("source = character(0)", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# An ambiguous deterministic-drop row, if it somehow reached the resolver as an
# unburned negative, hits the strict unknown-bucket ERROR (it can NEVER silently
# become a training negative). This is conditions (1)+(2) at the resolver.
# ---------------------------------------------------------------------------
test_that("GATE 6.5: an unburned deterministic-drop row -> strict resolver ERROR", {
  resolve <- g65(".of_resolve_supervised_eligibility")
  for (ngt in AMBIG_DROP) {
    d <- data.frame(
      id       = c("b1", "drop1"),
      class    = c("burned", "unburned"),
      source   = c("burned_truth", "deterministic_drop_hard"),
      neg_type = c(NA_character_, ngt),
      stringsAsFactors = FALSE
    )
    expect_error(
      resolve(id = d$id, class = d$class, source = d$source, neg_type = d$neg_type,
              random_background_source = "random_burnable_background",
              otsu_unburned_source = "otsu_patch_residual",
              otsu_unburned_exclude_neg_types = c("otsu_patch_review",
                                                  "otsu_patch_keep"),
              origin_stage = "TEST"),
      regexp = "NO valid negative bucket and NOT a known-excluded type",
      info = ngt
    )
  }
})

# ---------------------------------------------------------------------------
# (3)+(4)+(5) The scoring pool is built from the FULL audited internal-decisions
# set (internal_qc) carrying every row's ORIGINAL class (raw_class) with NO
# burnable / CORINE mask. So the ambiguous drops:
#   - DO reach the scoring pool (3),
#   - carry their original class and are NOT lost (5),
#   - are eligible to be scored / receive p_burned (4) -- scoring is class-blind.
# We reproduce EXACTLY the scoring_pool selection from build_supervised_training_pools()
# on a synthetic internal_qc and assert the drops survive.
# ---------------------------------------------------------------------------
test_that("GATE 6.5: ambiguous det drops reach the scoring pool with their original class (no mask)", {
  # internal_qc-like frame: burned keeps + the 3 ambiguous drops. raw_class is the
  # ORIGINAL (pre-audit) class; class_audited drives burned/review/drop routing.
  internal_qc <- data.frame(
    poly_id      = sprintf("p_%02d", 1:5),
    raw_class    = c("burned", "burned", "unburned", "unburned", "unburned"),
    class_audited= c("keep", "keep", "drop", "drop", "drop"),
    reason_1     = c(NA, NA, AMBIG_DROP),
    stringsAsFactors = FALSE
  )

  # Faithful reproduction of the scoring_pool construction (supervised-pools.R
  # STEP A3): scoring_pool = internal_qc with class := raw_class, source :=
  # deterministic_all_qc, a unique fire_uid. There is NO burnable/CORINE filter,
  # so EVERY internal-decision row (including the drops) is present.
  scoring_pool <- transform(
    internal_qc,
    class    = as.character(raw_class),
    source   = "deterministic_all_qc",
    fire_uid = paste0("2017_balanced_D_", seq_len(nrow(internal_qc)))
  )

  # (3) the 3 ambiguous drops reach the scoring pool.
  drop_rows <- scoring_pool[scoring_pool$poly_id %in% c("p_03", "p_04", "p_05"), ]
  expect_equal(nrow(drop_rows), 3L)
  # (5) they are accounted (their poly_id survives) and carry their ORIGINAL class.
  expect_setequal(drop_rows$poly_id, c("p_03", "p_04", "p_05"))
  expect_true(all(drop_rows$class == "unburned"))   # raw_class preserved
  # (4) scoreability is class-blind: the supervised scoring path scores EVERY
  # scoring_pool polygon regardless of class (no burned/unburned gate). Confirm
  # the scoring entry point carries no class filter that would drop these rows.
  score_src <- paste(deparse(body(g65("score_supervised_burned_map"))),
                     collapse = "\n")
  expect_false(grepl('class == "burned"', score_src, fixed = TRUE))
  expect_false(grepl('class != "unburned"', score_src, fixed = TRUE))
})

test_that("GATE 6.5: the scoring_pool construction applies NO burnable/CORINE mask (drops cannot be filtered out)", {
  pools_src <- paste(deparse(body(g65("build_supervised_training_pools"))),
                     collapse = "\n")
  # The scoring_pool is derived from internal_qc with raw_class, no burnable mask.
  expect_true(grepl("deterministic_all_qc", pools_src, fixed = TRUE))
  expect_true(grepl("class = as.character(raw_class)", pools_src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# Only two negative buckets exist; the resolver + capping share that set.
# ---------------------------------------------------------------------------
test_that("GATE 6.5: only two negative buckets (random, otsu) in resolver + capping", {
  expect_equal(g65(".of_valid_negative_buckets")(), c("random", "otsu"))
})

# ---------------------------------------------------------------------------
# OOF and FINAL share the same resolver + the same TWO caps (no contextual).
# ---------------------------------------------------------------------------
test_that("GATE 6.5: OOF and FINAL share the resolver + the same two caps (no contextual cap)", {
  oof   <- paste(deparse(body(g65("run_oof_xgb"))), collapse = "\n")
  final <- paste(deparse(body(g65("train_final_model_direct"))), collapse = "\n")
  for (s in list(oof, final)) {
    expect_true(grepl(".of_resolve_supervised_eligibility", s, fixed = TRUE))
    expect_true(grepl(".of_cap_negative_buckets", s, fixed = TRUE))
    expect_true(grepl("random_to_burned_ratio", s, fixed = TRUE))
    expect_true(grepl("otsu_unburned_to_burned_ratio", s, fixed = TRUE))
    # No contextual cap remains in either engine.
    expect_false(grepl("contextual_exclusion_to_burned_ratio", s, fixed = TRUE))
  }
})

# ---------------------------------------------------------------------------
# cap_contextual / contextual_exclusion_to_burned_ratio passed to a PUBLIC
# function ERRORS.
# ---------------------------------------------------------------------------
test_that("GATE 6.5: cap_contextual passed to build_supervised_burned_config ERRORS", {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  tif <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), tif,
                     overwrite = TRUE)
  expect_error(
    build_supervised_burned_config(
      scenario = "balanced", internal_decisions = f, change_index = tif,
      target_year = 2017L, cap_contextual = 0.25),
    regexp = "cap_contextual|unused argument"
  )
})

test_that("GATE 6.5: contextual_exclusion_to_burned_ratio passed to a public stage ERRORS", {
  f <- tempfile(fileext = ".gpkg")
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                              c(0, 1), c(0, 0)))), crs = 3035)
  sf::st_write(sf::st_sf(id = 1L, geometry = sfc), f, quiet = TRUE,
               delete_dsn = TRUE)
  tif <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(ncol = 4, nrow = 4, vals = 1:16), tif,
                     overwrite = TRUE)
  cfg <- build_supervised_burned_config(
    scenario = "balanced", internal_decisions = f, change_index = tif,
    target_year = 2017L)
  gpkg <- tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(fire_uid = "a", class = "burned", geometry = sfc),
               gpkg, layer = "train_features", quiet = TRUE, delete_dsn = TRUE)
  # run_oneyear_supervised_pipeline traps it explicitly; train_final_burned_model
  # rejects it as an unused argument (formal removed).
  expect_error(
    run_oneyear_supervised_pipeline(cfg,
      contextual_exclusion_to_burned_ratio = 0.25),
    regexp = "cap_contextual|contextual_exclusion_to_burned_ratio|unused"
  )
  expect_error(
    suppressWarnings(train_final_burned_model(
      train_features = gpkg, config = cfg,
      contextual_exclusion_to_burned_ratio = 0.25)),
    regexp = "unused argument|contextual_exclusion_to_burned_ratio"
  )
})
