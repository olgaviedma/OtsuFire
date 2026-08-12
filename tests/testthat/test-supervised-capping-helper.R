# =============================================================================
# Shared negative-bucket CAPPING helper contract + OOF<->FINAL equivalence.
#
# `.of_cap_negative_buckets()` is the SINGLE shared capping implementation used
# by BOTH the OOF per-fold trainer and the FINAL pool builder. This suite proves
# the cap semantics, determinism, validation, and that the helper reproduces the
# historical OOF and FINAL selected ids on the default (no-ambiguous) fixture.
#
# GATE 6.5 (2026-06-12): the CONTEXTUAL (deterministic-drop) bucket was removed.
# Valid buckets are now {random, otsu} only.
# =============================================================================

ns_cap <- asNamespace("OtsuFire")
capf       <- get(".of_cap_negative_buckets", envir = ns_cap)
resolve_e  <- get(".of_resolve_supervised_eligibility", envir = ns_cap)
valid_bk   <- get(".of_valid_negative_buckets", envir = ns_cap)()

CANON <- c(random = 1.0, otsu = 1.0)

# Build negatives_by_bucket from a simple spec: list(bucket = n_rows). Positives
# are the first n_burned ids; buckets get disjoint contiguous index blocks.
mk_buckets <- function(n_burned, sizes) {
  pos <- seq_len(n_burned)
  nxt <- n_burned
  byb <- list()
  for (b in valid_bk) {
    k <- sizes[[b]] %||% 0L
    byb[[b]] <- if (k > 0L) (nxt + 1L):(nxt + k) else integer(0)
    nxt <- nxt + k
  }
  total <- nxt
  list(positive_idx = pos, negatives_by_bucket = byb,
       id = sprintf("id_%04d", seq_len(total)))
}
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------------------------------------------------------
test_that("2 buckets present: all-under-cap selects exactly the cap targets", {
  s <- mk_buckets(40, list(random = 100, otsu = 100))
  out <- capf(s$positive_idx, s$negatives_by_bucket, n_burned = 40,
              caps = CANON, seed = 1L, id = s$id)
  a <- out$audit
  # random / otsu cap = ceil(40*1.0)=40.
  expect_equal(a$n_selected[a$bucket == "random"],     40L)
  expect_equal(a$n_selected[a$bucket == "otsu"],       40L)
  # No contextual / spectral bucket exists.
  expect_false("contextual" %in% a$bucket)
  expect_false("spectral" %in% a$bucket)
  # All burned positives kept.
  expect_true(all(seq_len(40) %in% out$selected_indices))
})

test_that("empty bucket -> 0 selected, reason 'empty_bucket'", {
  s <- mk_buckets(10, list(random = 0, otsu = 5))
  out <- capf(s$positive_idx, s$negatives_by_bucket, 10, CANON, seed = 2L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 0L)
  expect_equal(a$reason_when_short[a$bucket == "random"], "empty_bucket")
})

test_that("cap = 0 -> select 0 (reason 'cap=0')", {
  s <- mk_buckets(10, list(random = 20, otsu = 20))
  caps0 <- c(random = 0, otsu = 1)
  out <- capf(s$positive_idx, s$negatives_by_bucket, 10, caps0, seed = 3L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 0L)
  expect_equal(a$reason_when_short[a$bucket == "random"], "cap=0")
})

test_that("cap = Inf -> keep all available (reason 'cap=Inf')", {
  s <- mk_buckets(10, list(random = 20, otsu = 20))
  capsInf <- c(random = Inf, otsu = 1)
  out <- capf(s$positive_idx, s$negatives_by_bucket, 10, capsInf, seed = 4L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 20L)
  expect_equal(a$reason_when_short[a$bucket == "random"], "cap=Inf")
  expect_true(is.na(a$n_cap_max[a$bucket == "random"]))
})

test_that("available < cap_max -> take all (reason 'availability<cap')", {
  s <- mk_buckets(40, list(random = 5, otsu = 5))
  caps_lo <- c(random = 0.25, otsu = 0.25)
  out <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 5L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 5L)  # avail 5 < cap 10
  expect_equal(a$reason_when_short[a$bucket == "random"], "availability<cap")
})

test_that("available > cap_max -> take exactly cap_max (reason 'capped')", {
  s <- mk_buckets(40, list(random = 200, otsu = 5))
  caps_lo <- c(random = 0.25, otsu = 1.0)
  out <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 6L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 10L)  # cap 10 < avail 200
  expect_equal(a$reason_when_short[a$bucket == "random"], "capped")
})

test_that("n_burned == 0 -> all caps select 0, effective_ratio NA", {
  s <- mk_buckets(0, list(random = 10, otsu = 10))
  out <- capf(s$positive_idx, s$negatives_by_bucket, 0, CANON, seed = 7L, id = s$id)
  a <- out$audit
  expect_true(all(a$n_selected == 0L))
  expect_true(all(is.na(a$effective_ratio)))
  expect_true(all(a$reason_when_short == "n_burned=0"))
})

test_that("NULL / NA / non-scalar seed -> ERROR (no silent unseeded draw)", {
  s <- mk_buckets(10, list(random = 20, otsu = 20))
  expect_error(capf(s$positive_idx, s$negatives_by_bucket, 10, CANON,
                    seed = NULL, id = s$id), regexp = "single finite integer")
  expect_error(capf(s$positive_idx, s$negatives_by_bucket, 10, CANON,
                    seed = NA_integer_, id = s$id), regexp = "single finite integer")
  expect_error(capf(s$positive_idx, s$negatives_by_bucket, 10, CANON,
                    seed = c(1L, 2L), id = s$id), regexp = "single finite integer")
})

test_that("unknown bucket -> ERROR (no repair logic inside)", {
  s <- mk_buckets(10, list(random = 5, otsu = 5))
  bad <- s$negatives_by_bucket
  bad[["other"]] <- 100L:101L
  expect_error(capf(s$positive_idx, bad, 10, CANON, seed = 8L, id = s$id),
               regexp = "outside the canonical")
})

test_that("a 'contextual' bucket is now an UNKNOWN bucket -> explicit ERROR (GATE 6.5)", {
  # The contextual negative bucket was removed. Handing the helper a 'contextual'
  # key must hit the strict unknown-bucket error path, NOT be silently accepted.
  s <- mk_buckets(10, list(random = 5, otsu = 5))
  bad <- s$negatives_by_bucket
  bad[["contextual"]] <- 100L:101L
  expect_error(capf(s$positive_idx, bad, 10, CANON, seed = 8L, id = s$id),
               regexp = "outside the canonical")
  # The canonical caps must contain exactly the two buckets.
  expect_setequal(names(CANON), c("random", "otsu"))
})

test_that("duplicate row indices in a bucket are de-duped by row index", {
  s <- mk_buckets(40, list(random = 0, otsu = 0))
  byb <- s$negatives_by_bucket
  # bucket with duplicate row index 41, 41
  byb$random <- c(41L, 41L, 42L)
  id <- sprintf("id_%04d", 1:42)
  out <- capf(s$positive_idx, byb, 40, CANON, seed = 9L, id = id)
  a <- out$audit
  expect_equal(a$n_available[a$bucket == "random"], 2L)  # 41 + 42, deduped
})

test_that("reproducibility: same data+caps+seed -> exactly the same ids", {
  s <- mk_buckets(40, list(random = 200, otsu = 200))
  caps_lo <- c(random = 0.5, otsu = 0.5)
  o1 <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 123L, id = s$id)
  o2 <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 123L, id = s$id)
  expect_identical(o1$selected_ids, o2$selected_ids)
  expect_identical(o1$selected_indices, o2$selected_indices)
})

test_that("different seed may change ids but NOT the allowed per-bucket counts", {
  s <- mk_buckets(40, list(random = 200, otsu = 200))
  caps_lo <- c(random = 0.5, otsu = 0.5)
  o1 <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 1L, id = s$id)
  o2 <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 2L, id = s$id)
  expect_equal(o1$audit$n_selected, o2$audit$n_selected)  # counts invariant
  expect_false(identical(o1$selected_ids, o2$selected_ids))  # ids differ
})

test_that("sampling WITHOUT replacement + deterministic ascending output order", {
  s <- mk_buckets(40, list(random = 200, otsu = 200))
  caps_lo <- c(random = 0.5, otsu = 0.5)
  out <- capf(s$positive_idx, s$negatives_by_bucket, 40, caps_lo, seed = 11L, id = s$id)
  expect_equal(anyDuplicated(out$selected_indices), 0L)        # no row twice
  expect_identical(out$selected_indices, sort(out$selected_indices))  # ascending
})

test_that("singleton bucket (length-1 available, cap >=1) is handled (no sample() gotcha)", {
  # A length-1 available vector with cap >= 1 takes all (availability<cap branch,
  # no draw); with the unsafe base-R sample(x, 1) idiom this would mis-sample.
  s <- mk_buckets(40, list(random = 1, otsu = 0))
  out <- capf(s$positive_idx, s$negatives_by_bucket, 40, CANON, seed = 12L, id = s$id)
  a <- out$audit
  expect_equal(a$n_selected[a$bucket == "random"], 1L)
  expect_true(41L %in% out$selected_indices)
})

test_that("full audit has every documented column + one row per bucket", {
  s <- mk_buckets(40, list(random = 100, otsu = 100))
  out <- capf(s$positive_idx, s$negatives_by_bucket, 40, CANON, seed = 13L,
              id = s$id, context = "ctx_tag")
  a <- out$audit
  expect_setequal(a$bucket, valid_bk)
  expect_true(all(c("bucket", "n_burned", "n_available", "cap_ratio",
                    "n_cap_max", "n_selected", "effective_ratio", "seed",
                    "context", "reason_when_short") %in% names(a)))
  expect_true(all(a$context == "ctx_tag"))
  expect_true(all(a$seed == 13L))
  # effective_ratio = n_selected / n_burned.
  expect_equal(a$effective_ratio, a$n_selected / 40)
})

# ----------------------------------------------------------------------------
# EQUIVALENCE (§8, §12.3): the resolver+helper reproduce the HISTORICAL OOF and
# FINAL selected ids on the default (no-review/keep/unknown/NA) fixture.
# GATE 6.5: random + otsu buckets only.
# ----------------------------------------------------------------------------
RAND_SRC  <- "random_burnable_background"
OTSU_SRC  <- "otsu_patch_residual"
OTSU_EXCL <- c("otsu_patch_review", "otsu_patch_keep")

mk_equiv_pool <- function(seed = 7L) {
  set.seed(seed)
  rows <- list(); add <- function(k, cls, src, ngt) for (i in seq_len(k))
    rows[[length(rows) + 1L]] <<- data.frame(class = cls, source = src,
                                             neg_type = ngt, stringsAsFactors = FALSE)
  # Both buckets OVER-supplied so EVERY bucket samples (no take-all short-circuit)
  # -> historical OOF and FINAL consume the RNG identically and agree.
  add(40, "burned",   "burned_truth", NA_character_)
  add(60, "unburned", RAND_SRC,       "background_cell")
  add(60, "unburned", OTSU_SRC,       "otsu_patch_drop")
  d <- do.call(rbind, rows)
  d <- d[sample.int(nrow(d)), , drop = FALSE]; rownames(d) <- NULL
  d$id <- sprintf("uid_%03d", seq_len(nrow(d)))
  d
}

test_that("NEW resolver+helper == historical OOF == historical FINAL (all-sampling fixture)", {
  d <- mk_equiv_pool()
  caps <- c(random = 0.5, otsu = 0.5); seed <- 4242L

  e <- resolve_e(id = d$id, class = d$class, source = d$source, neg_type = d$neg_type,
    random_background_source = RAND_SRC, otsu_unburned_source = OTSU_SRC,
    otsu_unburned_exclude_neg_types = OTSU_EXCL)
  nb <- length(e$positive_idx)
  new <- capf(e$positive_idx, e$negatives_by_bucket, nb, caps, seed = seed, id = d$id)
  new_ids <- sort(new$selected_ids)

  # Historical OOF reproduction (cap_outer_train pick() logic), random -> otsu.
  old_oof <- local({
    cls <- d$class; src <- d$source; ngt <- d$neg_type
    is_b <- cls == "burned"; n_b <- sum(is_b)
    pick <- function(m, r) { a <- which(m); if (length(a) == 0 || !is.finite(r)) return(a)
      t <- ceiling(n_b * r); nt <- min(t, length(a)); if (nt <= 0) return(integer(0))
      if (nt >= length(a)) return(a); sort(sample(a, nt)) }
    rb  <- (!is_b) & (src %in% RAND_SRC)
    ot  <- (!is_b) & (src %in% OTSU_SRC) & !(ngt %in% OTSU_EXCL)
    set.seed(seed)
    sb <- which(is_b)
    sr <- pick(rb, caps["random"]);      so <- pick(ot, caps["otsu"])
    sort(unique(c(sb, sr, so)))
  })
  old_oof_ids <- sort(d$id[old_oof])

  # Historical FINAL reproduction (per-pool sample.int), random -> otsu.
  old_final_ids <- local({
    bp  <- d[d$class == "burned", ]
    rbp <- d[d$class == "unburned" & d$source %in% RAND_SRC, ]
    otp <- d[d$class == "unburned" & d$source %in% OTSU_SRC & !(d$neg_type %in% OTSU_EXCL), ]
    n_b <- nrow(bp)
    tr <- ceiling(n_b * caps["random"]);     to <- ceiling(n_b * caps["otsu"])
    set.seed(seed)
    sr <- rbp[sample.int(nrow(rbp), min(tr, nrow(rbp))), ]
    so <- otp[sample.int(nrow(otp), min(to, nrow(otp))), ]
    sort(unique(c(bp$id, sr$id, so$id)))
  })

  expect_identical(new_ids, old_oof_ids)
  expect_identical(new_ids, old_final_ids)
  expect_identical(old_oof_ids, old_final_ids)
})

test_that("equivalence: counts + effective ratios match the historical default selection", {
  d <- mk_equiv_pool()
  e <- resolve_e(id = d$id, class = d$class, source = d$source, neg_type = d$neg_type,
    random_background_source = RAND_SRC, otsu_unburned_source = OTSU_SRC,
    otsu_unburned_exclude_neg_types = OTSU_EXCL)
  nb <- length(e$positive_idx)
  out <- capf(e$positive_idx, e$negatives_by_bucket, nb, CANON, seed = 4242L, id = d$id)
  a <- out$audit
  # Both buckets over-supplied (avail 60 >= caps): selected == cap target.
  expect_equal(a$n_selected[a$bucket == "random"],     nb)
  expect_equal(a$n_selected[a$bucket == "otsu"],       nb)
  expect_equal(a$effective_ratio[a$bucket == "random"], 1.0)
})
