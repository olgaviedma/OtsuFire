# =====================================================================
# OtsuFire Tutorial — Chapter 1: The Standard Workflow
# Section 1.4: Feature importance — what drives the supervised model
# =====================================================================
#
# Once the supervised model is trained, the next diagnostic is to
# inspect *which* features the model relies on. XGBoost reports three
# importance metrics per feature:
#
#   - Gain:      improvement in the objective when a feature is used
#                in a split. The most informative metric for tree
#                models.
#   - Cover:     relative number of observations affected by splits
#                on this feature.
#   - Frequency: how often the feature is selected for splits.
#
# This section focuses on Gain. The package writes a ready-made
# `<prefix>_feature_importance.csv` to `07_FINAL_MODEL_V2/`, which we
# read directly rather than extracting from the model object.

library(sf)


# ---------------------------------------------------------------------
# 1.4.1 Loading the importance table
# ---------------------------------------------------------------------

prefix <- "2005_balanced_patch_certified"
model_dir <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results",
  "2005/Min_Min/SUPERVISED/balanced/07_FINAL_MODEL_V2"
)

fi <- read.csv(
  file.path(model_dir, paste0(prefix, "_feature_importance.csv")),
  stringsAsFactors = FALSE
)

cat("=== Global dimensions ===\n")
cat("Total features with Gain > 0:", nrow(fi), "\n")
cat("Total Gain (must sum to 1):", sum(fi$Gain), "\n")
cat("Cumulative Gain top 1:", sum(fi$Gain[1]), "\n")
cat("Cumulative Gain top 3:", sum(fi$Gain[1:3]), "\n")
cat("Cumulative Gain top 10:", sum(fi$Gain[1:10]), "\n")

# Expected output (2005 / balanced, OtsuFire 0.5.0):
#
#   Total features with Gain > 0: 18
#   Total Gain: 1
#   Cumulative Gain top 1: 0.7311
#   Cumulative Gain top 3: 0.9834
#   Cumulative Gain top 10: 0.9984
#
# Reading the output
# ---------------------------------------------------------------------
#
# Of the 51 features in the canonical whitelist
# (`OtsuFire:::.supervised_feature_cols`), only 18 produced any split
# during training. The remaining 33 were technically available to the
# tree builder but never selected — XGBoost found nothing useful to
# split on with them at any tree depth. This is informative: when a
# feature is dominant and well-placed, the rest of the feature space
# can become structurally redundant.
#
# A single feature accounts for 73% of the Gain. Three features
# account for 98%. The model is overwhelmingly concentrated.


# ---------------------------------------------------------------------
# 1.4.2 Top features
# ---------------------------------------------------------------------

cat("\n=== Top 15 features (Gain) ===\n")
print(head(fi[, c("Feature", "Gain", "Cover", "Frequency")], 15))

# Expected output (2005 / balanced, OtsuFire 0.5.0):
#
#             Feature         Gain        Cover   Frequency
#   1      hs_in_poly 0.7311236650 0.5305013791 0.333333333
#   2   hs_min_dist_m 0.2011433719 0.1260969291 0.047619048
#   3    hs_in_buffer 0.0510948242 0.0224831863 0.007326007
#   4         rbr_p90 0.0070389237 0.0539272053 0.117216117
#   5      rbr_aw_p90 0.0022225446 0.0073297511 0.076923077
#   6         rbr_p10 0.0016257487 0.0781858085 0.069597070
#   7         rbr_med 0.0013954300 0.0404584486 0.040293040
#   8    hs_conf_mean 0.0012028755 0.0029785795 0.051282051
#   9   cor_agri_frac 0.0011059278 0.0120414480 0.014652015
#   10     rbr_aw_iqr 0.0004503976 0.0474609124 0.058608059
#   11        elev_sd 0.0003626257 0.0015052862 0.043956044
#   12     hs_frp_sum 0.0003473229 0.0005533592 0.003663004
#   13        rbr_iqr 0.0002452919 0.0532769482 0.043956044
#   14 cor_forest_frac 0.0001988562 0.0022593311 0.040293040
#   15     rbr_aw_med 0.0001563405 0.0119981621 0.018315018
#
# Reading the output
# ---------------------------------------------------------------------
#
# The top three features are all hotspot-derived. Together they
# account for 98.3% of the Gain:
#
#   - `hs_in_poly` (73.1%): number of MODIS hotspots inside the
#     polygon. The dominant signal.
#   - `hs_min_dist_m` (20.1%): distance to the nearest hotspot.
#     Probably refining the binary "has hotspot" decision into a
#     graded "how close" decision near the boundary.
#   - `hs_in_buffer` (5.1%): hotspots in the spatial buffer around
#     the polygon.
#
# The first non-hotspot feature is `rbr_p90` at rank 4 with Gain
# 0.7%. Spectral information (RBR) carries less than 2% of the
# total Gain. Land cover (CORINE), topography, and DOY between them
# carry less than 0.2%.


# ---------------------------------------------------------------------
# 1.4.3 Aggregated Gain by feature family
# ---------------------------------------------------------------------
#
# It is more informative to aggregate the Gain across the six
# feature families defined in the whitelist:
#
#   G1_RBR_SW    — RBR same-window statistics (5 features)
#   G2_RBR_AW    — RBR all-window + persistence (7 features)
#   G3_CORINE    — Land cover fractions (7 features)
#   G4_TOPO      — Elevation + slope (14 features)
#   G5_DOY       — Day-of-year statistics (5 features)
#   G6_HOTSPOTS  — Hotspot evidence and distance (13 features)

classify_family <- function(feat) {
  is_isna <- grepl("_isNA$", feat)
  base <- sub("_isNA$", "", feat)
  
  fam <- if (base %in% c("rbr_valid_frac", "rbr_p10", "rbr_med",
                         "rbr_p90", "rbr_iqr")) {
    "G1_RBR_SW"
  } else if (base %in% c("rbr_aw_valid_frac", "rbr_aw_p10",
                         "rbr_aw_med", "rbr_aw_p90", "rbr_aw_iqr",
                         "persist_delta", "persist_ratio")) {
    "G2_RBR_AW"
  } else if (grepl("^cor_", base)) {
    "G3_CORINE"
  } else if (grepl("^(elev_|slope_)", base)) {
    "G4_TOPO"
  } else if (grepl("^doy_", base)) {
    "G5_DOY"
  } else if (grepl("^(hs_|hotspot_)", base)) {
    "G6_HOTSPOTS"
  } else {
    "OTHER"
  }
  list(family = fam, is_isna = is_isna)
}

fi$family  <- sapply(fi$Feature, function(x) classify_family(x)$family)
fi$is_isna <- sapply(fi$Feature, function(x) classify_family(x)$is_isna)

cat("\n=== Gain aggregated by family ===\n")
fam_gain <- aggregate(Gain ~ family, data = fi, FUN = sum)
fam_gain$pct <- 100 * fam_gain$Gain / sum(fam_gain$Gain)
fam_gain$n_features <- aggregate(Feature ~ family,
                                 data = fi, FUN = length)$Feature
fam_gain <- fam_gain[order(-fam_gain$Gain), ]
print(fam_gain)

# Expected output (2005 / balanced, OtsuFire 0.5.0):
#
#         family         Gain         pct n_features
#   6 G6_HOTSPOTS 0.9849809701 98.49809701          6
#   1   G1_RBR_SW 0.0103053943  1.03053943          4
#   2   G2_RBR_AW 0.0029227722  0.29227722          4
#   3   G3_CORINE 0.0013047840  0.13047840          2
#   4     G4_TOPO 0.0003626257  0.03626257          1
#   5      G5_DOY 0.0001234536  0.01234536          1
#
# Reading the output
# ---------------------------------------------------------------------
#
# The HOTSPOTS family carries 98.5% of the Gain. Even though only 6
# of the 13 hotspot features were used in splits, those 6 dominate
# the model entirely. The remaining 7 hotspot features (which the
# tree builder did not select) are nonetheless retained in the
# whitelist; the 0.5.0 architecture deliberately presents the full
# canonical feature space to the model and lets the tree builder
# decide.
#
# The RBR families together carry 1.3%, with same-window slightly
# more important than all-window. CORINE, topography and DOY are
# operationally negligible at the model level — though they remain
# methodologically necessary as control variables.
#
# This concentration is a structural property of the model under the
# 2005 / balanced configuration: when hotspot evidence is available
# and reliable, the model has very little incentive to use anything
# else.


# ---------------------------------------------------------------------
# 1.4.4 Sanity check: missing-value indicators carry no signal
# ---------------------------------------------------------------------
#
# The matrix builder synthesises a `<feature>_isNA` companion for
# every original feature, encoding whether the value was imputed.
# These flags are admissible features for the model. If any of them
# carried Gain, it would mean the *absence* of a measurement is
# itself predictive — a potential leakage pathway worth auditing.

cat("\n=== _isNA companions with Gain > 0 ===\n")
isna_only <- fi[fi$is_isna == TRUE & fi$Gain > 0, ]
if (nrow(isna_only) > 0) {
  print(isna_only[order(-isna_only$Gain),
                  c("Feature", "Gain", "family")])
} else {
  cat("None. Missing-value indicators contribute no Gain.\n")
}

# Expected output (2005 / balanced, OtsuFire 0.5.0):
#
#   None. Missing-value indicators contribute no Gain.
#
# Reading the output
# ---------------------------------------------------------------------
#
# No `_isNA` flag was ever selected as a split. This is a positive
# sanity check: the model does not rely on the pattern of missingness
# itself, which would be a leakage signature. Missingness is treated
# strictly as imputed numerical noise, not as an information
# channel.


# ---------------------------------------------------------------------
# 1.4.5 Methodological implications
# ---------------------------------------------------------------------
#
# The Baseline configuration concentrates 98.5% of the explanatory
# power on six hotspot-derived features. Three observations follow:
#
#   1. The supervised model under the canonical configuration is
#      essentially a hotspot detector with mild spectral refinement.
#      Its decisions are dominated by the binary signal "has at
#      least one MODIS hotspot inside the polygon", graded by the
#      distance to the nearest hotspot.
#
#   2. Spectral, topographic and land-cover features carry
#      negligible weight in the trained Baseline model — but this
#      does not mean they carry no information. It means the tree
#      builder, given full access to hotspots, finds it more
#      efficient to split on hotspots first. Whether spectral
#      features could replace hotspots if hotspots were unavailable
#      is a separate question, addressed in the experiments of
#      Chapter 3.
#
#   3. The 33 features in the canonical whitelist that received zero
#      Gain in this run are not necessarily redundant in absolute
#      terms — only in this configuration, with this training set,
#      and given the dominance of hotspot evidence. Removing them
#      preemptively would limit the model's adaptability to other
#      years or scenarios where hotspot evidence may be sparser or
#      absent.
#
# Reproducibility note. The `feature_whitelist_override_applied`
# field in `<prefix>_meta.txt` confirms that no override was active
# in this Baseline run; the model saw all 51 canonical features.
# Likewise, `feature_weights_applied: FALSE` confirms uniform
# unit weights. Subsequent experiments (Chapter 3) modify these
# settings and the same diagnostic recipe applies to their
# importance tables.