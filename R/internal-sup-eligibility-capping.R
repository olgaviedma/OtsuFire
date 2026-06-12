# =============================================================================
# Shared supervised ELIGIBILITY resolver + shared negative-bucket CAPPING helper.
#
# Methodological fix (2026-06-11): the supervised training population is defined
# by EXPLICIT class, NEVER by negation. A row is a:
#   - POSITIVE   iff class == "burned";
#   - NEGATIVE   iff class == "unburned" AND it maps to exactly ONE of the two
#                valid negative buckets (random, otsu);
#   - EXCLUDED   iff class == "unburned" but its (source, neg_type) is a KNOWN
#                excluded type (the otsu exclude-list, e.g. otsu_patch_review /
#                otsu_patch_keep) -> logged, never trained, never an error;
#   - ERROR      in every other case (unburned with no resolvable bucket and not
#                a known-excluded type; class NA; class not in {burned,unburned}).
#
# GATE 6.5 (2026-06-12): the CONTEXTUAL negative bucket was REMOVED entirely. An
# audit across 6 years proved no contextual category is a reliable negative
# (geo_excluded_hot / too_few_pixels / too_small carry no non-burn evidence;
# drop_outside_burnable, the only conceptually-valid negative, has zero cases).
# A deterministic drop does NOT automatically become an unburned label. The
# final negative architecture is TWO sources only: random background + Otsu
# residual. Deterministic-drop rows are therefore NEVER written as unburned
# negatives upstream, so an unburned row routed to "contextual" (or carrying
# source deterministic_drop_hard) now hits the strict unknown-bucket ERROR.
#
# The former OOF predicate `is_negative <- !is_burned` (plus its `sel_other` /
# `other_kept` "5th category") is METHODOLOGICALLY WRONG and is removed: review /
# keep(ambiguous) / NA / unknown rows can NEVER silently become negatives.
#
# Two PURE helpers, shared by BOTH the OOF per-fold trainer (run_oof_xgb) and the
# FINAL pool builder (train_final_model_direct):
#   .of_resolve_supervised_eligibility()  -- eligibility + bucket assignment,
#                                            UPSTREAM of capping.
#   .of_cap_negative_buckets()            -- capping ONLY, the single shared
#                                            implementation.
# PURE means: no globals / options / paths / session state read inside; every
# input arrives via an explicit argument; no file/RNG side effects except the
# single, explicitly-seeded sampling inside the capping helper.
# =============================================================================

# Canonical valid negative buckets, in canonical order. Used by BOTH helpers so
# the ordering of the seeded draws (and therefore the selected ids) is one fixed
# convention: random -> otsu. GATE 6.5 (2026-06-12): the "contextual" bucket was
# removed (deterministic drops are no longer training negatives).
.of_valid_negative_buckets <- function() {
  c("random", "otsu")
}

# Resolve the supervised training ELIGIBILITY of every row by EXPLICIT class +
# its (source, neg_type) metadata. Pure: no globals / options / paths read.
#
# @param id        character/atomic vector of row ids (id_col values), length n.
# @param class     character vector of class labels, length n.
# @param source    character vector of `source` metadata, length n (NA allowed).
# @param neg_type  character vector of `neg_type` metadata, length n (NA allowed).
# @param random_background_source         sources defining the RANDOM bucket.
# @param otsu_unburned_source             sources defining the OTSU bucket.
# @param otsu_unburned_exclude_neg_types  neg_types (within the otsu source) that
#        are KNOWN-EXCLUDED (otsu_patch_review / otsu_patch_keep): excluded +
#        logged, never trained, never an error.
# @param origin_stage character tag ("OOF" / "FINAL") woven into error messages.
#
# @return list(
#   positive_idx        : integer row indices of class=="burned" positives,
#   negative_idx        : integer row indices of eligible negatives,
#   negative_bucket     : character bucket tag aligned to negative_idx,
#   negatives_by_bucket : named list(bucket -> integer row indices),
#   excluded            : data.frame(id, class, source, neg_type, reason),
#   audit               : data.frame(class, bucket, n, action)
# )
# The helper does NOT infer or repair classes/buckets. It classifies by explicit
# metadata and errors if an unbucketed unburned negative ever appears (GATE 6.5:
# a deterministic-drop row reaching here as an unburned negative is exactly such
# an error, because det drops are no longer written as unburned upstream).
#
# @keywords internal
# @noRd
.of_resolve_supervised_eligibility <- function(
    id, class, source, neg_type,
    random_background_source,
    otsu_unburned_source,
    otsu_unburned_exclude_neg_types,
    origin_stage = "supervised"
) {
  n <- length(class)
  if (length(id) != n || length(source) != n || length(neg_type) != n) {
    stop(".of_resolve_supervised_eligibility(): id / class / source / neg_type ",
         "must be the SAME length (n=", n, ").", call. = FALSE)
  }
  cls <- as.character(class)
  src <- as.character(source)
  ngt <- as.character(neg_type)
  idv <- as.character(id)

  valid_buckets <- .of_valid_negative_buckets()

  # --- bucket predicates (mirror the FINAL pool filters exactly) -------------
  in_rand  <- !is.na(src) & (src %in% random_background_source)
  in_otsu  <- !is.na(src) & (src %in% otsu_unburned_source)
  otsu_excl <- !is.na(ngt) & (ngt %in% otsu_unburned_exclude_neg_types)

  bucket_of <- function(i) {
    if (in_rand[i])               return("random")
    if (in_otsu[i] && !otsu_excl[i]) return("otsu")
    NA_character_
  }

  is_burned   <- !is.na(cls) & cls == "burned"
  is_unburned <- !is.na(cls) & cls == "unburned"

  positive_idx <- which(is_burned)

  negative_idx    <- integer(0)
  negative_bucket <- character(0)
  excl_rows       <- list()

  n_excluded_otsu_review <- 0L

  for (i in which(is_unburned)) {
    b <- bucket_of(i)
    if (!is.na(b)) {
      negative_idx    <- c(negative_idx, i)
      negative_bucket <- c(negative_bucket, b)
      next
    }
    # No valid bucket. Either a KNOWN-excluded type (log) or an ERROR (strict).
    if (in_otsu[i] && otsu_excl[i]) {
      excl_rows[[length(excl_rows) + 1L]] <- data.frame(
        id = idv[i], class = cls[i], source = src[i], neg_type = ngt[i],
        reason = "otsu_excluded_neg_type", stringsAsFactors = FALSE
      )
      n_excluded_otsu_review <- n_excluded_otsu_review + 1L
      next
    }
    stop(sprintf(paste0(
      "[%s] .of_resolve_supervised_eligibility(): unburned row with NO valid ",
      "negative bucket and NOT a known-excluded type (id=%s, class=%s, ",
      "source=%s, neg_type=%s). An unburned negative MUST resolve to exactly ",
      "one of {%s} or be an otsu-excluded type. This is a metadata defect to ",
      "be fixed UPSTREAM at pool-build time, never inferred here."),
      origin_stage, idv[i], cls[i],
      ifelse(is.na(src[i]), "NA", src[i]),
      ifelse(is.na(ngt[i]), "NA", ngt[i]),
      paste(valid_buckets, collapse = ", ")), call. = FALSE)
  }

  # Anything not burned and not unburned (NA / unknown class) is a hard ERROR:
  # it can NEVER silently become a negative.
  bad <- which(!is_burned & !is_unburned)
  if (length(bad) > 0L) {
    i <- bad[1L]
    stop(sprintf(paste0(
      "[%s] .of_resolve_supervised_eligibility(): row with class NOT in ",
      "{burned, unburned} (id=%s, class=%s, source=%s, neg_type=%s). class is ",
      "either NA or an unknown label; it can never be a training row. %d such ",
      "row(s) found."),
      origin_stage, idv[i],
      ifelse(is.na(cls[i]), "NA", cls[i]),
      ifelse(is.na(src[i]), "NA", src[i]),
      ifelse(is.na(ngt[i]), "NA", ngt[i]),
      length(bad)), call. = FALSE)
  }

  excluded <- if (length(excl_rows)) {
    do.call(rbind, excl_rows)
  } else {
    data.frame(id = character(0), class = character(0), source = character(0),
               neg_type = character(0), reason = character(0),
               stringsAsFactors = FALSE)
  }

  negatives_by_bucket <- stats::setNames(
    lapply(valid_buckets, function(b) sort(negative_idx[negative_bucket == b])),
    valid_buckets
  )

  # --- audit table -----------------------------------------------------------
  n_burned_pos <- length(positive_idx)
  audit <- rbind(
    data.frame(class = "burned", bucket = "NA/n.a.", n = n_burned_pos,
               action = "keep positive", stringsAsFactors = FALSE),
    data.frame(class = "unburned", bucket = "random",
               n = length(negatives_by_bucket[["random"]]),
               action = "eligible", stringsAsFactors = FALSE),
    data.frame(class = "unburned", bucket = "otsu",
               n = length(negatives_by_bucket[["otsu"]]),
               action = "eligible", stringsAsFactors = FALSE),
    data.frame(class = "unburned", bucket = "otsu-excluded",
               n = n_excluded_otsu_review,
               action = "excluded", stringsAsFactors = FALSE)
  )

  list(
    positive_idx        = positive_idx,
    negative_idx        = negative_idx,
    negative_bucket     = negative_bucket,
    negatives_by_bucket = negatives_by_bucket,
    excluded            = excluded,
    audit               = audit
  )
}

# Apply the negative-bucket caps to ALREADY-RESOLVED positives + negatives. Pure
# except for the single explicitly-seeded draw. Receives positives and negatives
# grouped by VALID bucket; FAILS if handed any bucket outside the canonical two.
# Performs NO eligibility/repair logic.
#
# @param positive_idx       integer row indices of burned positives (all kept).
# @param negatives_by_bucket named list(bucket -> integer row indices), buckets
#        MUST be a subset of the canonical two.
# @param n_burned           number of burned positives (drives every cap).
# @param caps               named numeric(random, otsu)
#        cap ratios. ceiling(n_burned * cap) is the per-bucket max.
# @param seed               integer RNG seed (REQUIRED; NULL/invalid -> error;
#        no silent unseeded draw).
# @param id                 character/atomic vector of row ids, indexable by the
#        row indices above (for dedup + the selected_ids return).
# @param context            character tag woven into the audit (e.g. fold key).
#
# Cap semantics (canonical):
#   n_cap_max  <- ceiling(n_burned * cap_ratio)
#   n_selected <- min(n_available, n_cap_max)
#   cap == 0            -> select 0                       (reason "cap=0")
#   cap Inf/non-finite  -> keep all available             (reason "cap=Inf")
#   n_available == 0    -> empty                          (reason "empty_bucket")
#   n_available < n_cap_max -> take all                   (reason "availability<cap")
#   n_burned == 0       -> all caps 0, effective_ratio NA (reason "n_burned=0")
# All burned positives are kept; only negatives are subsampled; sampling WITHOUT
# replacement; deterministic ascending-row-index output order; same data+caps+
# seed -> exactly the same ids.
#
# @return list(selected_indices, selected_ids, audit) where audit has one row per
#   bucket: bucket, n_burned, n_available, cap_ratio, n_cap_max, n_selected,
#   effective_ratio, seed, context, reason_when_short.
#
# @keywords internal
# @noRd
.of_cap_negative_buckets <- function(
    positive_idx, negatives_by_bucket, n_burned, caps, seed, id,
    context = NA_character_
) {
  valid_buckets <- .of_valid_negative_buckets()

  # --- strict validation: NO repair logic inside -----------------------------
  if (is.null(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop(".of_cap_negative_buckets(): a single finite integer `seed` is ",
         "REQUIRED (no silent unseeded draw).", call. = FALSE)
  }
  if (!is.list(negatives_by_bucket)) {
    stop(".of_cap_negative_buckets(): `negatives_by_bucket` must be a named ",
         "list keyed by the canonical buckets.", call. = FALSE)
  }
  bnames <- names(negatives_by_bucket)
  if (is.null(bnames) || any(!nzchar(bnames))) {
    stop(".of_cap_negative_buckets(): `negatives_by_bucket` must be NAMED by ",
         "bucket.", call. = FALSE)
  }
  bad_bucket <- setdiff(bnames, valid_buckets)
  if (length(bad_bucket) > 0L) {
    stop(".of_cap_negative_buckets(): received bucket(s) outside the canonical ",
         "two {", paste(valid_buckets, collapse = ", "), "}: ",
         paste(bad_bucket, collapse = ", "),
         ". The caller must resolve eligibility FIRST; this helper performs no ",
         "repair.", call. = FALSE)
  }
  if (is.null(caps) || is.null(names(caps)) ||
      !all(valid_buckets %in% names(caps))) {
    stop(".of_cap_negative_buckets(): `caps` must be a named numeric with all ",
         "two buckets {", paste(valid_buckets, collapse = ", "), "}.",
         call. = FALSE)
  }

  cap_one <- function(avail, n_cap_max) {
    if (length(avail) == 0L) return(list(sel = integer(0), reason = "empty_bucket"))
    if (n_cap_max <= 0L)      return(list(sel = integer(0), reason = "cap=0"))
    if (n_cap_max >= length(avail)) {
      return(list(sel = avail, reason = "availability<cap"))
    }
    # Draw WITHOUT replacement over the available rows in ascending row-index
    # order, via sample.int (the singleton-safe primitive: sample(avail, k) is
    # avail[sample.int(length(avail), k)], so this preserves historical ids and
    # avoids the base-R sample() length-1 gotcha).
    sel <- avail[sample.int(length(avail), n_cap_max)]
    list(sel = sort(sel), reason = "capped")
  }

  set.seed(as.integer(seed))

  selected_neg <- integer(0)
  audit_rows   <- list()

  for (b in valid_buckets) {
    avail <- negatives_by_bucket[[b]]
    if (is.null(avail)) avail <- integer(0)
    avail <- sort(unique(as.integer(avail)))   # dedup by row index
    cap_ratio <- as.numeric(caps[[b]])

    if (n_burned <= 0L) {
      n_cap_max <- 0L
    } else if (!is.finite(cap_ratio)) {
      n_cap_max <- length(avail)               # Inf disables the cap
    } else {
      n_cap_max <- as.integer(ceiling(n_burned * cap_ratio))
    }

    if (length(avail) == 0L) {
      sel <- integer(0); reason <- "empty_bucket"
    } else if (n_burned <= 0L) {
      sel <- integer(0); reason <- "n_burned=0"
    } else if (!is.finite(cap_ratio)) {
      sel <- avail;       reason <- "cap=Inf"
    } else {
      cc <- cap_one(avail, n_cap_max)
      sel <- cc$sel; reason <- cc$reason
    }

    selected_neg <- c(selected_neg, sel)
    n_sel <- length(sel)
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      bucket          = b,
      n_burned        = n_burned,
      n_available     = length(avail),
      cap_ratio       = cap_ratio,
      n_cap_max       = if (is.finite(cap_ratio) && n_burned > 0L) n_cap_max else NA_integer_,
      n_selected      = n_sel,
      effective_ratio = if (n_burned > 0L) n_sel / n_burned else NA_real_,
      seed            = as.integer(seed),
      context         = as.character(context),
      reason_when_short = reason,
      stringsAsFactors = FALSE
    )
  }

  # All burned positives kept; only negatives were subsampled. Deterministic
  # ascending output order.
  selected_indices <- sort(unique(c(as.integer(positive_idx), selected_neg)))

  # A negative cannot appear twice (different buckets never share a row index;
  # assert it so a future bucket-overlap regression fails loud).
  if (anyDuplicated(selected_neg)) {
    stop(".of_cap_negative_buckets(): a negative row index was selected by more ",
         "than one bucket. Buckets must be disjoint.", call. = FALSE)
  }

  selected_ids <- as.character(id)[selected_indices]
  if (anyDuplicated(selected_ids)) {
    stop(".of_cap_negative_buckets(): duplicate ids in the selected set after ",
         "capping; id_col must be unique over the selected rows.", call. = FALSE)
  }

  list(
    selected_indices = selected_indices,
    selected_ids     = selected_ids,
    audit            = do.call(rbind, audit_rows)
  )
}
