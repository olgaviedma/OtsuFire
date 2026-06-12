# =============================================================================
# Supervised ELIGIBILITY resolver contract (2026-06-11 methodological fix;
# GATE 6.5 2026-06-12: contextual bucket removed).
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
#
# GATE 6.5: the CONTEXTUAL (deterministic-drop) bucket was removed. The valid
# negative buckets are now {random, otsu} ONLY. A deterministic-drop row that
# reaches the resolver as an unburned negative has NO valid bucket and therefore
# hits the strict ERROR (it should never reach here, because det drops are no
# longer written as unburned upstream).
# =============================================================================

ns_elig <- asNamespace("OtsuFire")
resolve_elig <- get(".of_resolve_supervised_eligibility", envir = ns_elig)

# Canonical bucket-definition params (the production defaults).
DET_SRC   <- "deterministic_drop_hard"
RAND_SRC  <- "random_burnable_background"
OTSU_SRC  <- "otsu_patch_residual"
OTSU_EXCL <- c("otsu_patch_review", "otsu_patch_keep")

resolve_default <- function(df, origin = "TEST") {
  resolve_elig(
    id = df$id, class = df$class, source = df$source, neg_type = df$neg_type,
    random_background_source         = RAND_SRC,
    otsu_unburned_source             = OTSU_SRC,
    otsu_unburned_exclude_neg_types  = OTSU_EXCL,
    origin_stage = origin
  )
}

# A row builder mirroring the 2017 balanced DEFAULT structure. GATE 6.5: the
# negative pool is random + otsu ONLY (no deterministic-drop negatives).
mk_default_pool <- function() {
  rows <- list()
  add <- function(k, cls, src, ngt) for (i in seq_len(k)) {
    rows[[length(rows) + 1L]] <<- data.frame(class = cls, source = src,
                                             neg_type = ngt,
                                             stringsAsFactors = FALSE)
  }
  add(20, "burned",   "burned_truth",               NA_character_)
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
  expect_equal(length(e$negatives_by_bucket$random),     24L)
  expect_equal(length(e$negatives_by_bucket$otsu),       12L)
  # GATE 6.5: contextual bucket no longer exists.
  expect_null(e$negatives_by_bucket$contextual)
  # Spectral bucket no longer exists.
  expect_null(e$negatives_by_bucket$spectral)
  expect_setequal(names(e$negatives_by_bucket),
                  c("random", "otsu"))
  # No excluded / error in the clean default pool.
  expect_equal(nrow(e$excluded), 0L)
  # Audit table has the canonical rows + actions.
  expect_true(all(c("keep positive", "eligible", "excluded") %in% e$audit$action))
  expect_equal(e$audit$n[e$audit$class == "burned"], 20L)
  # No contextual / spectral row in the audit table.
  expect_false("contextual" %in% e$audit$bucket)
  expect_false("spectral" %in% e$audit$bucket)
})

test_that("the two bucket index sets are DISJOINT and cover all eligible negatives", {
  d <- mk_default_pool()
  e <- resolve_default(d)
  all_neg <- unlist(e$negatives_by_bucket, use.names = FALSE)
  expect_equal(anyDuplicated(all_neg), 0L)
  expect_setequal(sort(all_neg), sort(e$negative_idx))
})

test_that("ONLY two negative buckets exist (random, otsu)", {
  expect_equal(get(".of_valid_negative_buckets", envir = ns_elig)(),
               c("random", "otsu"))
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

test_that("GATE 6.5: a deterministic_drop_hard unburned row -> strict ERROR (no contextual bucket)", {
  # A deterministic-drop row reaching the resolver as an unburned negative has
  # NO valid bucket (contextual was removed) and is not an otsu-excluded type,
  # so it MUST hit the strict unknown-bucket error. This guards against silently
  # re-creating the contextual bucket.
  d <- mk_default_pool()
  d <- rbind(d, data.frame(class = "unburned", source = DET_SRC,
                           neg_type = "geo_excluded_hot", id = "det1",
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

test_that("toggling otsu-exclude membership does NOT turn excluded rows into negatives", {
  # otsu_patch_review is excluded by default. A row whose source has no bucket
  # and is not excluded would ERROR -- it can never silently become a negative.
  # Here we keep the canonical exclude-list: the review row is excluded, not a
  # negative.
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

test_that("GATE 6.5: the deterministic builder no longer stamps det drops as unburned", {
  # The former block filtered class_final=="drop" and stamped
  # class="unburned" + source="deterministic_drop_hard". GATE 6.5 removed it, so
  # the builder source must no longer turn drops into unburned negatives.
  src <- paste(deparse(
    body(get("build_unburned_from_deterministic_decisions", envir = ns_elig))),
    collapse = "\n")
  expect_false(grepl('source = "deterministic_drop_hard"', src, fixed = TRUE))
  expect_false(grepl('class_final == "drop"', src, fixed = TRUE))
})
