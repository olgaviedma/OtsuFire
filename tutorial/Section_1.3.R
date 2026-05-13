
# =====================================================================
# Section 1.3: Reviews rescued by the supervised stage
# =====================================================================
#
# Reviews are the only class where the supervised model can change
# the deterministic outcome at the polygon level. Counting how many
# reviews are "rescued" (i.e., classified as burned with high
# probability) is the operational measure of supervised value.

review_p <- by_class[["review"]]

cat("=== Reviews rescued by p_burned_model threshold ===\n")
for (thr in c(0.5, 0.7, 0.9, 0.95, 0.99)) {
  n_rescued <- sum(review_p > thr, na.rm = TRUE)
  pct <- 100 * n_rescued / length(review_p)
  cat(sprintf("  p > %.2f : %4d / %d reviews rescued (%.1f%%)\n",
              thr, n_rescued, length(review_p), pct))
}

# Expected output:

#  p > 0.50 :   30 / 1114 reviews rescued (2.7%)
#  p > 0.70 :   30 / 1114 reviews rescued (2.7%)
#  p > 0.90 :   30 / 1114 reviews rescued (2.7%)
#  p > 0.95 :   30 / 1114 reviews rescued (2.7%)
#  p > 0.99 :   27 / 1114 reviews rescued (2.4%)

# Distribution of reviews across the [0,1] probability range
breaks <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
review_bins <- cut(review_p, breaks = breaks, include.lowest = TRUE,
                   right = TRUE)
print(table(review_bins))

# Expected output:
#
#   review_bins
#     [0,0.1] (0.1,0.3] (0.3,0.5] (0.5,0.7] (0.7,0.9]   (0.9,1]
#        1084         0         0         0         0        30
#
# Reading the output
# ---------------------------------------------------------------------
#
# The review distribution is strictly bimodal. 1084 reviews (97.3%)
# fall below 0.1, and 30 reviews (2.7%) fall above 0.9. The interval
# [0.1, 0.9] is empty — the supervised model commits to one side or
# the other, with no ambiguous polygons in between.
#
# This bimodality has two practical consequences:
#
# 1. The choice of public-map threshold is robust. Any threshold
#    between 0.1 and 0.9 yields the same rescued set (30 polygons).
#    The pipeline is therefore not sensitive to threshold tuning,
#    which is desirable for reproducibility.
#
# 2. The supervised stage replaces 1114 unsupervised "uncertain"
#    labels with 1084 confirmed unburned + 30 confirmed burned.
#    Without the supervised stage, this resolution would not be
#    available.


# =====================================================================
# Section 1.3.1: Profile of the rescued reviews
# =====================================================================
#
# The 30 rescued reviews are the operational output of the
# supervised stage. Inspecting their characteristics is essential
# both for scientific interpretation and for the manuscript.

rescued_idx <- fm$class_final == "review" & fm$p_burned_model > 0.5
fm_rescued  <- fm[rescued_idx, ]

cat(sprintf("=== %d reviews rescued (p_burned_model > 0.5) ===\n\n",
            nrow(fm_rescued)))

# Hotspot evidence
cat("--- hotspot_available ---\n")
print(table(fm_rescued$hotspot_available, useNA = "always"))

cat("\n--- hs_in_poly stats ---\n")
print(summary(fm_rescued$hs_in_poly))

# Spectral signature
cat("\n--- rbr_med stats ---\n")
print(summary(fm_rescued$rbr_med))

cat("\n--- rbr_p90 stats ---\n")
print(summary(fm_rescued$rbr_p90))

# Patch size
cat("\n--- area_ha stats ---\n")
print(summary(fm_rescued$area_ha))

# Deterministic provenance — why was each rescued polygon initially
# flagged as "review"?
cat("\n--- preyear_action ---\n")
print(table(fm_rescued$preyear_action, useNA = "always"))

cat("\n--- preyear_reason ---\n")
print(table(fm_rescued$preyear_reason, useNA = "always"))

# Expected output (2005 / balanced):
#
#   hotspot_available: 1 (all 30)
#   hs_in_poly:        median 1, max 6 (every rescued polygon has
#                      at least one MODIS hotspot inside)
#   rbr_med:           median 389, max 544 (high spectral change)
#   rbr_p90:           median 508, max 628
#   area_ha:           median 49, max 590, mean 81 (substantial size)
#
#   preyear_action:    drop=1, keep=12, not_applied=14, review=3
#   preyear_reason:    no_overlap_previous_year=12,
#                      not_applied=14,
#                      overlap_previous_year_removed=3,
#                      previous_year_conflict=1
#
# Reading the output
# ---------------------------------------------------------------------
#
# The rescued reviews share a consistent profile:
#
#   - All 30 are MODIS-supported. Every rescued polygon contains at
#     least one validated hotspot, and on average they contain
#     between one and two.
#   - All 30 carry strong spectral evidence of fire. Their RBR
#     percentiles cluster well above the values typical of unburned
#     vegetation.
#   - All 30 are operationally relevant patches. Median area is
#     ~49 ha (and the 30 polygons together represent ~2400 ha of
#     burned area that would otherwise be lost to the deterministic
#     "review" bucket).
#
# The deterministic provenance reveals why these polygons could not
# be resolved by the deterministic filters alone: the previous-year
# audit step (`preyear`) introduced ambiguity for most of them
# (12 had a previous-year keep that did not spatially overlap, 14
# had no previous-year audit applied at all, 3 had previous-year
# data that was removed by sanitisation, and 1 was in active
# previous-year conflict). In all such cases the deterministic
# stage refused to commit, and labelled the polygon "review" so
# the supervised stage could resolve it.
#
# The supervised model resolves this uncertainty using information
# the deterministic filters do not access: the joint pattern of
# hotspot evidence and full RBR percentile distribution. The result
# is 30 high-confidence rescues that the burned-area product would
# otherwise miss.
#
# A reviewer-facing question that this analysis raises: are the 30
# rescued polygons actually burned in the EFFIS reference? This
# question is addressed in Chapter 5 (validation against external
# references); here we only confirm that the rescues are
# methodologically defensible.