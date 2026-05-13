# =====================================================================
# OtsuFire Tutorial — Chapter 1: The Standard Workflow
# Section 1.2: Probability distributions by deterministic class
# =====================================================================
#
# The supervised model assigns a `p_burned_model` probability to every
# polygon in the scoring universe. The first scientific diagnostic is
# to inspect how this probability distributes across the deterministic
# classes (drop / keep / review).
#
# Use `final_map_full` for analysis throughout this chapter — it
# carries the complete column schema and the full polygon universe.

library(sf)

fm <- sf::read_sf(final_map_path, layer = "final_map_full")

# Distribution by class
by_class <- split(fm$p_burned_model, fm$class_final)

cat("=== p_burned_model distribution by class_final ===\n\n")
for (cls in names(by_class)) {
  cat(sprintf("--- %s (n=%d) ---\n", cls, length(by_class[[cls]])))
  print(summary(by_class[[cls]]))
  cat(sprintf("  sd  = %.5f\n", sd(by_class[[cls]], na.rm = TRUE)))
  cat(sprintf("  p10 = %.5f\n", quantile(by_class[[cls]], 0.10, na.rm = TRUE)))
  cat(sprintf("  p90 = %.5f\n\n", quantile(by_class[[cls]], 0.90, na.rm = TRUE)))
}

#  output (2005 / balanced):
#
#   --- drop (n=247) ---
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.0005711 0.0018709 0.0022149 0.0034024 0.0032391 0.0242665 
#   sd  = 0.00380
#   p10 = 0.00116
#   p90 = 0.00603

#  --- keep (n=544) ---
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9878  0.9987  0.9991  0.9988  0.9993  0.9993 
#  sd  = 0.00117
#  p10 = 0.99830
#  p90 = 0.99931
#
#   --- review (n=1114) ---
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.0005849 0.0050525 0.0104789 0.0408152 0.0224186 0.9992946 
#   sd  = 0.15939
#   p10 = 0.00287
#   p90 = 0.03849
#
# Reading the output
# ---------------------------------------------------------------------
#
# Drop (n=247). The model is highly confident these are NOT burned.
# All drops have p_burned_model < 0.02. The distribution is compact
# (sd = 0.003), with no overlap with the keep distribution. This
# confirms the deterministic stage was right to discard them, and
# the supervised model adds no surprise here.
#
# Keep (n=544). The model is highly confident these ARE burned. All
# keeps have p_burned_model > 0.987 (sd = 0.001). The two clusters
# (drop and keep) are perfectly separated — there is no probability
# overlap whatsoever. This is the signature of a correctly fitted
# model that learned a clean decision boundary on the labelled
# training set.
#
# Review (n=1114). This is where the supervised stage adds value.
# Reviews are the unsupervised "I don't know" class — polygons
# that the unsupervised filters could not classify with confidence.
# The supervised model resolves this uncertainty into a granular
# probability:
#   - the median review has p_burned_model = 0.006 (very unlikely
#     burned)
#   - p90 = 0.027, but max = 0.999, indicating a small subset is
#     classified as burned with high confidence
#   - the standard deviation is 0.16 — two orders of magnitude
#     larger than drop or keep, reflecting genuine bimodality in
#     this group
#
# The next sub-section quantifies how many reviews the supervised
# stage rescues at different probability thresholds.

