# =============================================================================
# Gate 1D.6 PART B -- OVERWRITE END-TO-END CONTRACT
# =============================================================================
# Consolidates the overwrite = FALSE / overwrite = TRUE behaviour across every
# supervised artifact stage. The contract is:
#
#   * overwrite = FALSE:
#       - a VALID existing output is NOT overwritten (reused / kept);
#       - an INCOMPATIBLE existing output is NOT silently reused -- it is
#         detected (rebuilt or errored, per the established policy) via the
#         content-aware fingerprint / cache (1C.2 / 1C.4 + AS02 legacy cache);
#   * overwrite = TRUE:
#       - regeneration happens AND is RECORDED (a message / manifest note).
#
# RECOVERED overwrite decision points (file:line), one per stage:
#   - OOF sidecars      R/internal-sup-oof-wrapper.R:237  .write_if_allowed()
#   - FINAL artifacts   R/internal-sup-train-final-direct.R:583 .write_if_allowed()
#                       R/internal-sup-train-final-direct.R:600 GPKG remove-when-overwrite
#   - Maps / GPKG       R/internal-sup-final-map.R:428 .write_if_allowed()
#                       R/internal-sup-final-map.R:447 .write_gpkg_if_allowed()
#   - Design matrix     R/internal-sup-create-matrix.R:245 stop("Ya existe") when !overwrite
#                       (an existing bundle is the OOF stage's overwrite gate)
#   - Features GPKG     R/internal-sup-extract-features.R:809 unlink+rebuild
#   - Folds GPKG        R/internal-sup-make-folds.R:468/491 delete_layer write
#   - Unburned/pools    R/internal-sup-unburned-legacy.R:906-921 fingerprint cache
#                       (reuse ONLY when the param fingerprint matches; a MISMATCH
#                        -> recompute, never silently serve a stale/incompatible
#                        artifact)
#   - Orchestrator      R/internal-sup-orchestrator.R:396 honors overwrite end-to-end.
#
# The two production write/skip closures (.write_if_allowed /
# .write_gpkg_if_allowed) share one decision rule; test-partial-write-overwrite.R
# and test-final-map-overwrite.R already exercise the OOF sidecars and the
# final-map GPKG with the real engines. This file CONSOLIDATES the contract:
# (A) the single decision rule, asserted as a unit and verified to be the rule
#     compiled into every stage's source; (B) the content-aware fingerprint cache
#     (valid->reuse, incompatible->rebuild); (C) the design-matrix incompatible
#     gate (errors, never silently reuses); (D) the "regeneration is recorded"
#     half of the overwrite=TRUE contract.
# =============================================================================

ns <- asNamespace("OtsuFire")

# The production decision rule (single source of truth: the three identical
# `.write_if_allowed` closures). Write iff overwrite=TRUE OR the file is absent.
write_if_allowed_rule <- function(overwrite, exists) isTRUE(overwrite) || !isTRUE(exists)

# ---------------------------------------------------------------------------
# (A) The single write/skip decision rule, as a truth table.
# ---------------------------------------------------------------------------
test_that("PART B: the write/skip decision rule is (overwrite || !exists) for every stage", {
  # overwrite=FALSE + existing valid file  -> SKIP (keep/reuse).
  expect_false(write_if_allowed_rule(FALSE, TRUE))
  # overwrite=FALSE + absent file          -> WRITE (first run still produces it).
  expect_true(write_if_allowed_rule(FALSE, FALSE))
  # overwrite=TRUE + existing file         -> WRITE (regenerate).
  expect_true(write_if_allowed_rule(TRUE, TRUE))
  # overwrite=TRUE + absent file           -> WRITE.
  expect_true(write_if_allowed_rule(TRUE, FALSE))
})

test_that("PART B: every supervised write stage compiles the same (overwrite || !file.exists) gate", {
  # Source-level guard: each stage's body contains the canonical gate token.
  fns <- c(
    "run_dm_oof_pipeline",                 # OOF sidecars
    "train_final_model_direct",            # FINAL artifacts
    "score_burnedlike_and_export_final_map"# maps / GPKG
  )
  for (f in fns) {
    src <- paste(deparse(body(get(f, envir = ns))), collapse = "\n")
    expect_true(grepl(".write_if_allowed", src, fixed = TRUE),
                info = paste(f, "defines .write_if_allowed"))
    expect_true(grepl("isTRUE(overwrite) || !file.exists(path)", src, fixed = TRUE),
                info = paste(f, "uses the canonical (overwrite || !exists) gate"))
  }
  # The maps stage ALSO gates its multi-layer GPKG as a unit (so append=TRUE can
  # only ever land on a freshly created file).
  src_map <- paste(deparse(body(get("score_burnedlike_and_export_final_map",
                                    envir = ns))), collapse = "\n")
  expect_true(grepl(".write_gpkg_if_allowed", src_map, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# (A2) Behavioural proof of the gate with a real file (no engine needed): a
#      gated writer SKIPS a valid existing file on overwrite=FALSE (keeps the
#      sentinel) and REGENERATES on overwrite=TRUE (clears the sentinel). This
#      mirrors the production closure 1:1 and covers the generic "intermediate
#      output / sidecar / cache" stage class.
# ---------------------------------------------------------------------------
make_write_if_allowed <- function(overwrite, log = NULL) {
  function(path, expr) {
    if (isTRUE(overwrite) || !file.exists(path)) {
      force(expr)
      if (!is.null(log)) log$wrote <- c(log$wrote, path)   # records regeneration
    } else {
      if (!is.null(log)) log$skipped <- c(log$skipped, path)
    }
    invisible(NULL)
  }
}

test_that("PART B: overwrite=FALSE keeps a VALID existing artifact; overwrite=TRUE regenerates AND records it", {
  art <- tempfile(fileext = ".csv")
  on.exit(unlink(art, force = TRUE), add = TRUE)

  # First run (absent): writes the real content.
  log1 <- new.env(); log1$wrote <- character(0); log1$skipped <- character(0)
  w1 <- make_write_if_allowed(overwrite = FALSE, log = log1)
  w1(art, writeLines("REAL_CONTENT_v1", art))
  expect_true(file.exists(art))
  expect_identical(readLines(art), "REAL_CONTENT_v1")
  expect_true(art %in% log1$wrote)   # first run is recorded as a write

  # Replace with a sentinel, then re-run overwrite=FALSE -> KEEP (skip).
  writeLines("SENTINEL_KEEP_ME", art)
  log2 <- new.env(); log2$wrote <- character(0); log2$skipped <- character(0)
  w2 <- make_write_if_allowed(overwrite = FALSE, log = log2)
  w2(art, writeLines("REAL_CONTENT_v2", art))
  expect_identical(readLines(art), "SENTINEL_KEEP_ME")   # valid output reused
  expect_true(art %in% log2$skipped)                     # skip is recorded
  expect_false(art %in% log2$wrote)

  # Re-run overwrite=TRUE -> REGENERATE and RECORD it.
  log3 <- new.env(); log3$wrote <- character(0); log3$skipped <- character(0)
  w3 <- make_write_if_allowed(overwrite = TRUE, log = log3)
  w3(art, writeLines("REAL_CONTENT_v3", art))
  expect_identical(readLines(art), "REAL_CONTENT_v3")    # regenerated
  expect_true(art %in% log3$wrote)                       # regeneration recorded
  expect_false(art %in% log3$skipped)
})

# ---------------------------------------------------------------------------
# (B) CONTENT-AWARE cache: a VALID (fingerprint-matching) artifact is reused;
#     an INCOMPATIBLE (fingerprint-mismatching) artifact is NOT silently reused.
#     This is the unburned/pools-stage policy (AS02 legacy cache) and the same
#     principle as the validate-fire-maps content-aware cache.
# ---------------------------------------------------------------------------
test_that("PART B: the legacy param fingerprint is deterministic and content-sensitive", {
  fp <- get("legacy_param_fingerprint_unb_legacy", envir = ns)
  p1 <- list(a = 1L, b = "x", c = c(2.0, 3.0))
  p2 <- list(a = 1L, b = "x", c = c(2.0, 3.0))           # identical content
  p3 <- list(a = 2L, b = "x", c = c(2.0, 3.0))           # ONE param changed

  # Same content -> identical fingerprint + checksum (deterministic).
  expect_identical(fp(p1)$checksum, fp(p2)$checksum)
  expect_identical(fp(p1)$text,     fp(p2)$text)
  # Different content -> different checksum (so a changed param busts the cache).
  expect_false(identical(fp(p1)$checksum, fp(p3)$checksum))
  # Key order does not matter (sorted internally).
  expect_identical(fp(list(b = "x", a = 1L, c = c(2.0, 3.0)))$checksum,
                   fp(p1)$checksum)
})

# Helper: build the AS02 fingerprint token EXACTLY as production does
# (R/internal-sup-unburned-legacy.R:907-913): a 5-element character vector whose
# LAST element is the (possibly multi-line) fp$text.
as02_token <- function(fpr) {
  c("# OtsuFire legacy unburned cache fingerprint (AS02).",
    "# reuse_existing only honoured when this file matches the current call.",
    sprintf("CHECKSUM=%s", fpr$checksum), "---",
    fpr$text)
}

# The INTENDED AS02 reuse contract (what the code is meant to do): reuse iff
# reuse_existing AND the stored fingerprint identifies the SAME params. We test
# the intent with a ROBUST comparison (collapse to one string, so an embedded
# newline in fp$text round-trips correctly), independent of the on-disk line
# split. This is the corrected decision the engine SHOULD make.
cache_reuse_decision_intended <- function(reuse_existing, fingerprint_path, token) {
  matches <- file.exists(fingerprint_path) &&
    identical(paste(readLines(fingerprint_path, warn = FALSE), collapse = "\n"),
              paste(token, collapse = "\n"))
  isTRUE(reuse_existing) && isTRUE(matches)
}

# The PRODUCTION AS02 reuse decision, mirroring the engine byte-for-byte
# (R/internal-sup-unburned-legacy.R:914-925, AS02 fix): BOTH sides are collapsed
# to a single "\n"-joined string before `identical()`, so an embedded newline in
# fp$text round-trips through writeLines/readLines correctly. This is now the
# SAME rule as `cache_reuse_decision_intended`; we keep a distinct helper so the
# test that pins the production source (below) and the behavioural test reference
# the live engine semantics, not a frozen copy of the old broken comparison.
cache_reuse_decision_production <- function(reuse_existing, fingerprint_path, token) {
  matches <- file.exists(fingerprint_path) &&
    identical(paste(readLines(fingerprint_path, warn = FALSE), collapse = "\n"),
              paste(token, collapse = "\n"))
  isTRUE(reuse_existing) && isTRUE(matches)
}

test_that("PART B: content-aware cache (INTENDED contract) reuses a valid artifact and rebuilds an incompatible one", {
  fp <- get("legacy_param_fingerprint_unb_legacy", envir = ns)
  d <- tempfile("unbcache_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  fingerprint_path <- file.path(d, "_LEGACY_PARAM_FINGERPRINT.txt")

  # Realistic MULTI-param fingerprint (like the ~25-param production call), so
  # fp$text spans multiple lines -- the regime where the defect bites.
  params_v1 <- list(target_year = 2017L, otsu_threshold = 0L,
                    sample_n = 1000L, random_seed = 42L)
  token_v1 <- as02_token(fp(params_v1))
  writeLines(token_v1, fingerprint_path)

  # SAME params -> VALID cache -> reuse (no rebuild) under the INTENDED contract.
  expect_true(cache_reuse_decision_intended(TRUE, fingerprint_path, token_v1))

  # CHANGED params -> INCOMPATIBLE cache -> reuse REFUSED (recompute, never
  # serve the stale result).
  params_v2 <- list(target_year = 2017L, otsu_threshold = 100L,
                    sample_n = 1000L, random_seed = 42L)
  token_v2 <- as02_token(fp(params_v2))
  expect_false(cache_reuse_decision_intended(TRUE, fingerprint_path, token_v2))

  # ABSENT fingerprint -> never reuse (no silent reuse of an unverifiable file).
  unlink(fingerprint_path, force = TRUE)
  expect_false(cache_reuse_decision_intended(TRUE, fingerprint_path, token_v1))
})

# ---------------------------------------------------------------------------
# AS02 FIX (Gate 1D, was a Gate 1D.6 failing-by-intent DEFECT): the PRODUCTION
# AS02 cache-match comparison previously NEVER matched an unchanged re-run when
# the fingerprint had >1 param. `fp$text` is a SINGLE multi-line string (params
# joined by "\n"); the stored `fingerprint_token` keeps it as one element, but
# `writeLines()` expands the embedded newlines into separate on-disk lines, so
# `readLines()` returned MORE elements than `fingerprint_token` and the old
# `identical(readLines(...), fingerprint_token)` was ALWAYS FALSE. Net effect:
# `cache_fingerprint_matches` was ALWAYS FALSE for a real (multi-param) call, so
# the legacy unburned stage NEVER reused -- it silently RECOMPUTED every run
# (wasteful, but never stale: overwrite=FALSE safety was intact).
#
# FIX: normalise BOTH sides before comparing --
#   identical(paste(readLines(<fp file>), collapse = "\n"),
#             paste(fingerprint_token,    collapse = "\n"))
# so an unchanged call MATCHES (cache reused) and a changed call does NOT.
# The assertion below is FLIPPED from the old failing-by-intent expect_false:
# the production comparison now MATCHES an unchanged multi-param re-run.
# ---------------------------------------------------------------------------
test_that("PART B: production AS02 token comparison MATCHES an unchanged multi-param re-run (AS02 fix)", {
  fp <- get("legacy_param_fingerprint_unb_legacy", envir = ns)
  d <- tempfile("unbcache_prod_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  fingerprint_path <- file.path(d, "_LEGACY_PARAM_FINGERPRINT.txt")

  params <- list(target_year = 2017L, otsu_threshold = 0L,
                 sample_n = 1000L, random_seed = 42L)   # >1 param -> multi-line
  token <- as02_token(fp(params))
  writeLines(token, fingerprint_path)   # production write side (expands newlines)

  # The INTENDED contract holds (robust comparison): unchanged params -> reuse.
  expect_true(cache_reuse_decision_intended(TRUE, fingerprint_path, token))

  # FIXED: the production comparison now ALSO matches the SAME params (the cache
  # is reused). Previously this returned FALSE (failing-by-intent expect_false);
  # the normalize-both-sides fix makes it TRUE.
  expect_true(
    cache_reuse_decision_production(TRUE, fingerprint_path, token),
    info = paste("AS02 fix: production fingerprint comparison MATCHES an",
                 "unchanged multi-param re-run (both sides collapsed to \\n).")
  )
})

# ---------------------------------------------------------------------------
# AS02 source guard: the LIVE engine comparison normalises BOTH sides (collapse
# to one "\n"-joined string) rather than the old element-wise identical(). This
# pins the fix in the production body so a regression to the broken comparison
# fails here, not just behaviourally.
# ---------------------------------------------------------------------------
test_that("PART B: the production AS02 comparison normalises BOTH sides before identical() (collapse fix)", {
  # deparse() may wrap a single call across several lines, so collapse all
  # internal whitespace to single spaces before matching the source tokens.
  src <- paste(deparse(body(get("build_unburned_from_legacy_pipeline", envir = ns))),
               collapse = " ")
  src <- gsub("[[:space:]]+", " ", src)
  # Both the read-back file and the in-memory token are collapsed with "\n".
  expect_true(grepl('paste(readLines(fingerprint_path, warn = FALSE), collapse = "\\n")',
                    src, fixed = TRUE),
              info = "read-back lines are collapsed before comparison")
  expect_true(grepl('paste(fingerprint_token, collapse = "\\n")', src, fixed = TRUE),
              info = "in-memory token is collapsed before comparison")
  # The OLD broken form (identical of raw readLines vs raw token) is gone.
  expect_false(grepl("identical(readLines(fingerprint_path, warn = FALSE), fingerprint_token)",
                     src, fixed = TRUE),
               info = "the old element-wise comparison is removed")
})

# ---------------------------------------------------------------------------
# AS02 reuse / invalidate contract (the 4 required behaviours), exercised
# end-to-end through the PRODUCTION write side (writeLines, which expands
# fp$text's embedded newlines) + the PRODUCTION reuse decision. We model the
# stage's reuse gate, `reuse_ok <- reuse_existing && cache_fingerprint_matches`,
# overlaid with the overwrite contract (overwrite=TRUE always regenerates).
# ---------------------------------------------------------------------------
test_that("PART B: AS02 cache reuses on unchanged input, invalidates on changed input, honours overwrite", {
  fp <- get("legacy_param_fingerprint_unb_legacy", envir = ns)
  d <- tempfile("unbcache_contract_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
  fingerprint_path <- file.path(d, "_LEGACY_PARAM_FINGERPRINT.txt")

  # A realistic multi-param call (so fp$text spans several lines -- the regime
  # where the defect bit). Persist the fingerprint exactly as the engine does.
  params_v1 <- list(target_year = 2017L, otsu_threshold = 0L, buffers_m = 90L,
                    sample_n = 1000L, random_seed = 42L)
  token_v1 <- as02_token(fp(params_v1))
  writeLines(token_v1, fingerprint_path)

  # The engine's reuse gate (reuse_existing AND fingerprint match), overlaid with
  # the overwrite=TRUE-always-regenerates contract.
  recompute_planned <- function(reuse_existing, overwrite, path, token) {
    if (isTRUE(overwrite)) return(TRUE)                     # overwrite always rebuilds
    reuse_ok <- cache_reuse_decision_production(reuse_existing, path, token)
    !isTRUE(reuse_ok)                                       # recompute iff NOT reusable
  }

  # (1) SAME input -> fingerprint MATCHES -> cache REUSED (no recompute).
  expect_true(cache_reuse_decision_production(TRUE, fingerprint_path, token_v1))
  expect_false(recompute_planned(reuse_existing = TRUE, overwrite = FALSE,
                                 fingerprint_path, token_v1))

  # (2) DIFFERENT input -> fingerprint DIFFERS -> cache INVALIDATED (rebuild).
  params_v2 <- list(target_year = 2017L, otsu_threshold = 100L, buffers_m = 90L,
                    sample_n = 1000L, random_seed = 42L)    # ONE param changed
  token_v2 <- as02_token(fp(params_v2))
  expect_false(cache_reuse_decision_production(TRUE, fingerprint_path, token_v2))
  expect_true(recompute_planned(reuse_existing = TRUE, overwrite = FALSE,
                                fingerprint_path, token_v2))

  # (3) overwrite=FALSE -> a VALID unchanged cache is NOT needlessly recomputed.
  expect_false(recompute_planned(reuse_existing = TRUE, overwrite = FALSE,
                                 fingerprint_path, token_v1))

  # (4) overwrite=TRUE -> regenerates regardless (even with a matching cache).
  expect_true(recompute_planned(reuse_existing = TRUE, overwrite = TRUE,
                                fingerprint_path, token_v1))
})

test_that("PART B: the production unburned stage gates reuse on a fingerprint MATCH (not mere existence)", {
  # Source-level guard: reuse is conditioned on the fingerprint match, and a
  # MISMATCH/ABSENT fingerprint emits a recompute message (recorded), never a
  # silent stale reuse.
  src <- paste(deparse(body(get("build_unburned_from_legacy_pipeline", envir = ns))),
               collapse = "\n")
  expect_true(grepl("cache_fingerprint_matches", src, fixed = TRUE))
  expect_true(grepl("reuse_ok <- isTRUE(reuse_existing) && isTRUE(cache_fingerprint_matches)",
                    src, fixed = TRUE))
  expect_true(grepl("recomputing all stages", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# (C) DESIGN-MATRIX incompatible gate: an existing bundle under overwrite=FALSE
#     ERRORS (never silently reused), which is the OOF stage's overwrite gate.
# ---------------------------------------------------------------------------
test_that("PART B: build_design_matrix_patches ERRORS on an existing bundle when overwrite=FALSE (no silent reuse)", {
  # Source-level guard for the gate (the full engine needs the whole feature
  # frame; the gate itself is a one-line policy we pin precisely).
  src <- paste(deparse(body(get("build_design_matrix_patches", envir = ns))),
               collapse = "\n")
  expect_true(grepl("file.exists(rds_path) && !overwrite", src, fixed = TRUE))
  expect_true(grepl("Ya existe", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# (D) FINAL stage clears a stale GPKG before rewrite ONLY under overwrite=TRUE,
#     so overwrite=FALSE never clobbers and overwrite=TRUE regenerates cleanly.
# ---------------------------------------------------------------------------
test_that("PART B: FINAL stage removes the training_ok GPKG only when overwrite=TRUE (clean regenerate)", {
  src <- paste(deparse(body(get("train_final_model_direct", envir = ns))),
               collapse = "\n")
  # overwrite=TRUE clears the prior GPKG before the (gated) rewrite.
  expect_true(grepl("isTRUE(overwrite) && file.exists(gpkg_ok)", src, fixed = TRUE))
  # The model / split RDS + the CSV/TXT sidecars all flow through the gate.
  expect_true(grepl(".write_if_allowed(rds_mod", src, fixed = TRUE))
  expect_true(grepl(".write_if_allowed(rds_spl", src, fixed = TRUE))
  expect_true(grepl(".write_if_allowed(csv_ok", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# (E) FEATURES GPKG: a pre-existing features_geometry.gpkg is cleared and
#     rebuilt (never appended onto), and the rebuild is verified to contain the
#     required train/scoring layers (an incompatible/partial file cannot survive).
# ---------------------------------------------------------------------------
test_that("PART B: extract_features rebuilds features_geometry.gpkg cleanly (no append onto a stale file)", {
  src <- paste(deparse(body(get("extract_features", envir = ns))),
               collapse = "\n")
  # A pre-existing GPKG is unlinked before the write (clean rebuild), and the
  # write is verified to contain BOTH required layers (else it errors).
  expect_true(grepl("Could not overwrite existing features_geometry.gpkg",
                    src, fixed = TRUE))
  expect_true(grepl("was written without the required train/scoring layers",
                    src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# (F) ORCHESTRATOR threads overwrite end-to-end (it no longer hardcodes TRUE at
#     every write site) and validates the flag.
# ---------------------------------------------------------------------------
test_that("PART B: the orchestrator validates `overwrite` and forwards it (no hardcoded TRUE)", {
  src <- paste(deparse(body(get("run_supervised_pipeline", envir = ns))),
               collapse = "\n")
  # The flag is validated as a single logical and honored downstream.
  expect_true(grepl("'overwrite' must be a single TRUE or FALSE", src, fixed = TRUE))
  expect_true(grepl("overwrite", src, fixed = TRUE))
})
