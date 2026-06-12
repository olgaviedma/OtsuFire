# =============================================================================
# Integration / symmetry: OOF and FINAL share the SAME eligibility resolver AND
# the SAME capping helper; given the same subset+caps+seed they select the same
# negatives; positives are always kept; scale_pos_weight is computed AFTER
# sampling; outer_test never enters the helper.
# =============================================================================

ns_int <- asNamespace("OtsuFire")
g_int  <- function(nm) get(nm, envir = ns_int)

# ---------------------------------------------------------------------------
# SOURCE-LEVEL: both engines reference BOTH shared helpers; no dead negative
# machinery remains.
# ---------------------------------------------------------------------------
test_that("OOF and FINAL both call the SAME resolver + capping helper", {
  oof   <- paste(deparse(body(g_int("run_oof_xgb"))), collapse = "\n")
  final <- paste(deparse(body(g_int("train_final_model_direct"))), collapse = "\n")
  for (s in list(oof, final)) {
    expect_true(grepl(".of_resolve_supervised_eligibility", s, fixed = TRUE))
    expect_true(grepl(".of_cap_negative_buckets", s, fixed = TRUE))
  }
})

test_that("no productive sel_other / other_kept / !is_burned negative definition remains in supervised engine bodies", {
  # Inspect the DEPARSED bodies of every supervised engine that could carry the
  # old negative-by-negation machinery. This is robust to the test cwd (uses the
  # loaded namespace, not file paths). Comments are stripped by deparse, so any
  # match is PRODUCTIVE code.
  fns <- c("run_oof_xgb", "train_final_model_direct",
           ".of_resolve_supervised_eligibility", ".of_cap_negative_buckets")
  offending <- character(0)
  for (nm in fns) {
    src <- paste(deparse(body(g_int(nm))), collapse = "\n")
    if (grepl("sel_other", src, fixed = TRUE) ||
        grepl("other_kept", src, fixed = TRUE) ||
        grepl("known_mask", src, fixed = TRUE)) {
      offending <- c(offending, nm)
    }
    # `!is_burned` used as a NEGATIVE definition. The resolver legitimately uses
    # `!is_burned & !is_unburned` to identify the ERROR set, which is NOT a
    # negative definition -- allow exactly that idiom.
    if (grepl("!is_burned", src, fixed = TRUE) &&
        !grepl("!is_burned & !is_unburned", src, fixed = TRUE)) {
      offending <- c(offending, paste0(nm, ":!is_burned"))
    }
  }
  expect_equal(offending, character(0),
               info = paste("offending:", paste(offending, collapse = ", ")))
})

# ---------------------------------------------------------------------------
# SYMMETRY: the resolver gives the SAME eligible population to OOF and FINAL for
# the same subset, and the helper gives the same negatives for same subset+seed.
# ---------------------------------------------------------------------------
# GATE 6.5 (2026-06-12): contextual bucket removed; random + otsu only.
RAND_SRC <- "random_burnable_background"; OTSU_SRC <- "otsu_patch_residual"
OTSU_EXCL <- c("otsu_patch_review", "otsu_patch_keep")
CAPS <- c(random = 1.0, otsu = 1.0)

mk_pool <- function(seed = 3L) {
  set.seed(seed)
  rows <- list(); add <- function(k, cls, src, ngt) for (i in seq_len(k))
    rows[[length(rows) + 1L]] <<- data.frame(class = cls, source = src,
                                             neg_type = ngt, stringsAsFactors = FALSE)
  add(30, "burned",   "burned_truth", NA_character_)
  add(45, "unburned", RAND_SRC,       "background_cell")
  add(38, "unburned", OTSU_SRC,       "otsu_patch_drop")
  d <- do.call(rbind, rows)
  d <- d[sample.int(nrow(d)), , drop = FALSE]; rownames(d) <- NULL
  d$id <- sprintf("uid_%03d", seq_len(nrow(d)))
  d
}

resolve_pool <- function(d, origin) {
  g_int(".of_resolve_supervised_eligibility")(
    id = d$id, class = d$class, source = d$source, neg_type = d$neg_type,
    random_background_source = RAND_SRC, otsu_unburned_source = OTSU_SRC,
    otsu_unburned_exclude_neg_types = OTSU_EXCL, origin_stage = origin)
}

test_that("OOF and FINAL resolve the SAME eligible population for the same subset", {
  d <- mk_pool()
  eo <- resolve_pool(d, "OOF")
  ef <- resolve_pool(d, "FINAL")
  expect_identical(eo$positive_idx, ef$positive_idx)
  expect_identical(eo$negatives_by_bucket, ef$negatives_by_bucket)
  expect_identical(eo$negative_idx, ef$negative_idx)
})

test_that("same subset + caps + seed -> EXACTLY the same negatives (OOF vs FINAL helper call)", {
  d <- mk_pool()
  e <- resolve_pool(d, "OOF")
  nb <- length(e$positive_idx)
  capf <- g_int(".of_cap_negative_buckets")
  o1 <- capf(e$positive_idx, e$negatives_by_bucket, nb, CAPS, seed = 555L, id = d$id,
             context = "OOF")
  o2 <- capf(e$positive_idx, e$negatives_by_bucket, nb, CAPS, seed = 555L, id = d$id,
             context = "FINAL")
  expect_identical(o1$selected_ids, o2$selected_ids)
  # Context tag is the only legitimate audit difference.
  expect_identical(o1$audit$n_selected, o2$audit$n_selected)
  expect_false(identical(o1$audit$context, o2$audit$context))
})

test_that("positives are ALWAYS kept regardless of caps", {
  d <- mk_pool()
  e <- resolve_pool(d, "FINAL")
  nb <- length(e$positive_idx)
  capf <- g_int(".of_cap_negative_buckets")
  caps0 <- c(random = 0, otsu = 0)
  out <- capf(e$positive_idx, e$negatives_by_bucket, nb, caps0, seed = 1L, id = d$id)
  expect_true(all(e$positive_idx %in% out$selected_indices))
  # With all caps 0, ONLY positives survive.
  expect_setequal(out$selected_indices, e$positive_idx)
})

test_that("different subsets keep the SAME policy + effective-ratio computation", {
  d1 <- mk_pool(seed = 1L); d2 <- mk_pool(seed = 2L)
  capf <- g_int(".of_cap_negative_buckets")
  e1 <- resolve_pool(d1, "OOF"); e2 <- resolve_pool(d2, "OOF")
  o1 <- capf(e1$positive_idx, e1$negatives_by_bucket, length(e1$positive_idx),
             CAPS, seed = 9L, id = d1$id)
  o2 <- capf(e2$positive_idx, e2$negatives_by_bucket, length(e2$positive_idx),
             CAPS, seed = 9L, id = d2$id)
  # Same effective-ratio rule: effective = n_selected / n_burned in both.
  expect_equal(o1$audit$effective_ratio, o1$audit$n_selected / length(e1$positive_idx))
  expect_equal(o2$audit$effective_ratio, o2$audit$n_selected / length(e2$positive_idx))
})

# ---------------------------------------------------------------------------
# scale_pos_weight computed AFTER sampling, on the data ACTUALLY used.
# ---------------------------------------------------------------------------
test_that("scale_pos_weight is computed inside the shared core AFTER caps (both paths)", {
  # The shared core computes spw on the train rows it RECEIVES, which are the
  # already-capped rows (FINAL: L_ok; OOF: capped outer-train). Confirm the spw
  # rule lives in the shared core, fed the post-cap training frame.
  core <- paste(deparse(body(g_int(".of_nested_refit_fit"))), collapse = "\n")
  expect_true(grepl("scale_pos_weight", core) ||
              grepl("spw", core))
  # FINAL hands the core L_df (derived from the CAPPED L_ok), and the spw used in
  # the deployed params is fit$spw_refit (post-sampling).
  final <- paste(deparse(body(g_int("train_final_model_direct"))), collapse = "\n")
  expect_true(grepl("spw <- fit$spw_refit", final, fixed = TRUE))
  # L_ok is built from the capping helper's selected rows BEFORE the core call.
  expect_true(grepl(".of_cap_negative_buckets", final, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# outer_test never enters the resolver or the helper (OOF).
# ---------------------------------------------------------------------------
test_that("OOF caps only the outer-train; outer_test never enters resolver/helper", {
  oof <- paste(deparse(body(g_int("run_oof_xgb"))), collapse = "\n")
  # cap_outer_train is called on idx_tr (outer-train); idx_te (outer-test) is
  # predicted only AFTER the refit, never capped.
  expect_true(grepl("cap_outer_train(idx_tr", oof, fixed = TRUE))
  expect_false(grepl("cap_outer_train(idx_te", oof, fixed = TRUE))
  # The resolver/helper are invoked inside cap_outer_train (which receives tr_idx).
  expect_true(grepl(".of_resolve_supervised_eligibility", oof, fixed = TRUE))
})
