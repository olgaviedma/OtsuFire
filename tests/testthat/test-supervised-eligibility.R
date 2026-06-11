# =============================================================================
# Supervised ELIGIBILITY resolver contract (2026-06-11 methodological fix).
#
# The supervised training population is defined by EXPLICIT class, NEVER by
# negation. This suite proves `.of_resolve_supervised_eligibility()`:
#   - class=="burned"                       -> positive;
#   - class=="unburned" + a VALID bucket    -> negative tagged with that bucket;
#   - class=="unburned" + otsu exclude type -> EXCLUDED + logged (not error);
#   - class=="unburned" + NO bucket / not excludable -> ERROR (strict);
#   - class NA / unknown                    -> ERROR (strict);
# and that review / keep(ambiguous) / NA / unknown rows can NEVER become
# negatives, regardless of any toggle.
# =============================================================================

ns_elig <- asNamespace("OtsuFire")
resolve_elig <- get(".of_resolve_supervised_eligibility", envir = ns_elig)

# Canonical bucket-definition params (the production defaults).
DET_SRC   <- "deterministic_drop_hard"
SPEC_NGT  <- "spectral_reject_medium"
RAND_SRC  <- "random_burnable_background"
OTSU_SRC  <- "otsu_patch_residual"
OTSU_EXCL <- c("otsu_patch_review", "otsu_patch_keep")

resolve_default <- function(df, origin = "TEST") {
  resolve_elig(
    id = df$id, class = df$class, source = df$source, neg_type = df$neg_type,
    deterministic_drop_source        = DET_SRC,
    spectral_hard_negative_neg_types = SPEC_NGT,
    random_background_source         = RAND_SRC,
    otsu_unburned_source             = OTSU_SRC,
    otsu_unburned_exclude_neg_types  = OTSU_EXCL,
    origin_stage = origin
  )
}

# A row builder mirroring the 2017 balanced DEFAULT structure.
mk_default_pool <- function() {
  rows <- list()
  add <- function(k, cls, src, ngt) for (i in seq_len(k)) {
    rows[[length(rows) + 1L]] <<- data.frame(class = cls, source = src,
                                             neg_type = ngt,
                                             stringsAsFactors = FALSE)
  }
  add(20, "burned",   "burned_truth",               NA_character_)
  add(16, "unburned", DET_SRC,                      "geo_excluded_hot")        # contextual
  add( 8, "unburned", DET_SRC,                      SPEC_NGT)                  # spectral
  add(24, "unburned", RAND_SRC,                     "background_cell")         # random
  add(12, "unburned", OTSU_SRC,                     "otsu_patch_drop")         # otsu
  d <- do.call(rbind, rows)
  d$id <- sprintf("uid_%03d", seq_len(nrow(d)))
  d
}

# ----------------------------------------------------------------------------
test_that("burned -> positive; each unburned bucket -> the right negative bucket", {
  d <- mk_default_pool()
  e <- resolve_default(d)
  expect_equal(length(e$positive_idx), 20L)
  expect_true(all(d$class[e$positive_idx] == "burned"))
  expect_equal(length(e$negatives_by_bucket$contextual), 16L)
  expect_equal(length(e$negatives_by_bucket$spectral),    8L)
  expect_equal(length(e$negatives_by_bucket$random),     24L)
  expect_equal(length(e$negatives_by_bucket$otsu),       12L)
  # No excluded / error in the clean default pool.
  expect_equal(nrow(e$excluded), 0L)
  # Audit table has the canonical rows + actions.
  expect_true(all(c("keep positive", "eligible", "excluded") %in% e$audit$action))
  expect_equal(e$audit$n[e$audit$class == "burned"], 20L)
})

test_that("the four bucket index sets are DISJOINT and cover all eligible negatives", {
  d <- mk_default_pool()
  e <- resolve_default(d)
  all_neg <- unlist(e$negatives_by_bucket, use.names = FALSE)
  expect_equal(anyDuplicated(all_neg), 0L)
  expect_setequal(sort(all_neg), sort(e$negative_idx))
})

test_that("otsu review/keep are EXCLUDED + logged, never negatives, never error", {
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = "unburned", source = OTSU_SRC,
                           neg_type = "otsu_patch_review",
                           id = "rev_1", stringsAsFactors = FALSE))
  d <- rbind(d, data.frame(class = "unburned", source = OTSU_SRC,
                           neg_type = "otsu_patch_keep",
                           id = "keep_1", stringsAsFactors = FALSE))
  e <- resolve_default(d)
  # otsu bucket unchanged (drop only); the 2 excluded rows are logged.
  expect_equal(length(e$negatives_by_bucket$otsu), 12L)
  expect_equal(nrow(e$excluded), 2L)
  expect_setequal(e$excluded$id, c("rev_1", "keep_1"))
  expect_true(all(e$excluded$reason == "otsu_excluded_neg_type"))
  # They are NOT in any negative bucket.
  expect_false(any(c("rev_1", "keep_1") %in% d$id[e$negative_idx]))
})

test_that("class == 'review' is NEVER a negative (ERROR: unknown class)", {
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = "review", source = DET_SRC,
                           neg_type = "ambiguous", id = "rv1",
                           stringsAsFactors = FALSE))
  expect_error(resolve_default(d), regexp = "class NOT in \\{burned, unburned\\}")
})

test_that("class == 'keep' (ambiguous) is NEVER a negative (ERROR: unknown class)", {
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = "keep", source = OTSU_SRC,
                           neg_type = "otsu_patch_keep", id = "kp1",
                           stringsAsFactors = FALSE))
  expect_error(resolve_default(d), regexp = "class NOT in \\{burned, unburned\\}")
})

test_that("class == NA -> ERROR (strict), with context fields", {
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = NA_character_, source = RAND_SRC,
                           neg_type = "background_cell", id = "na1",
                           stringsAsFactors = FALSE))
  expect_error(resolve_default(d), regexp = "class NOT in \\{burned, unburned\\}")
})

test_that("unburned with NO resolvable bucket (and not excludable) -> ERROR (strict)", {
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = "unburned", source = "some_unknown_source",
                           neg_type = "whatever", id = "ub1",
                           stringsAsFactors = FALSE))
  expect_error(resolve_default(d),
               regexp = "NO valid negative bucket and NOT a known-excluded type")
})

test_that("unburned with neg_type NA but a VALID source still buckets (random/otsu)", {
  # random_background and otsu sources are defined by SOURCE; neg_type NA still
  # buckets random; otsu with NA neg_type is not in the exclude list -> otsu.
  d <- data.frame(
    class    = c("burned", "unburned", "unburned"),
    source   = c("burned_truth", RAND_SRC, OTSU_SRC),
    neg_type = c(NA_character_, NA_character_, NA_character_),
    id       = c("b1", "r1", "o1"), stringsAsFactors = FALSE
  )
  e <- resolve_default(d)
  expect_equal(length(e$negatives_by_bucket$random), 1L)
  expect_equal(length(e$negatives_by_bucket$otsu), 1L)
})

test_that("det-drop without spectral neg_type -> contextual; with it -> spectral (upstream classification)", {
  d <- data.frame(
    class    = c("burned", "unburned", "unburned"),
    source   = c("burned_truth", DET_SRC, DET_SRC),
    neg_type = c(NA_character_, "geo_excluded_hot", SPEC_NGT),
    id       = c("b1", "c1", "s1"), stringsAsFactors = FALSE
  )
  e <- resolve_default(d)
  expect_equal(d$id[e$negatives_by_bucket$contextual], "c1")
  expect_equal(d$id[e$negatives_by_bucket$spectral],   "s1")
})

test_that("toggling otsu-exclude membership does NOT turn excluded rows into negatives", {
  # otsu_patch_review is excluded by default. Even if a (hypothetical) different
  # exclude-list were passed, a row whose source has no bucket and is not
  # excluded would ERROR -- it can never silently become a negative. Here we keep
  # the canonical exclude-list: the review row is excluded, not a negative.
  d <- data.frame(
    class    = c("burned", "unburned"),
    source   = c("burned_truth", OTSU_SRC),
    neg_type = c(NA_character_, "otsu_patch_review"),
    id       = c("b1", "rev"), stringsAsFactors = FALSE
  )
  e <- resolve_default(d)
  expect_equal(length(e$negative_idx), 0L)
  expect_equal(nrow(e$excluded), 1L)
})

test_that("det-drop neg_type is NEVER NA upstream (case_when TRUE ~ 'drop_hard' fallback)", {
  # The deterministic builder's case_when has a `TRUE ~ \"drop_hard\"` terminal
  # fallback, so every deterministic_drop_hard row carries a non-NA neg_type ->
  # the unbucketed-det-drop case is VACUOUS by construction. Confirm the fallback
  # literal is present in the builder source.
  src <- paste(deparse(
    body(get("build_unburned_from_deterministic_decisions", envir = ns_elig))),
    collapse = "\n")
  expect_true(grepl('TRUE ~ "drop_hard"', src, fixed = TRUE))
})
