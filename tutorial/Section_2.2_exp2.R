# =====================================================================
# OtsuFire Tutorial — Chapter 3: Custom feature-engineering experiments
# Experiment 2 (E2): Removing hotspots + amplifying temporal persistence
# =====================================================================
#
# Experiment E1 (Section 3.1) demonstrated that removing the hotspot
# block collapses the supervised model into a generic spectral-change
# detector that over-rescues review polygons by a factor of ~31. The
# closing question of E1 was whether features representing temporal
# persistence — burned scars typically remain spectrally darkened
# for 2-3 years, while many non-fire spectral changes (agricultural
# cycles, ephemeral hydric events) reverse within a season — could
# recover the lost specificity if amplified explicitly.
#
# Experiment E2 tests this hypothesis by re-running E1 with the same
# hotspot-free whitelist, but adding a 5x weight to the two
# persistence features in the canonical feature space:
#   - `persist_delta`  (RBR same-window minus RBR all-window)
#   - `persist_ratio`  (RBR all-window / RBR same-window)
#
# Both belong to the RBR all-window family (G2) and capture, in
# different parameterisations, whether the spectral change observed
# in the target year is anomalous relative to the surrounding years
# or whether the polygon already showed reduced reflectance before
# the target year (more consistent with persistent, fire-related
# scarring) or recovered quickly after (more consistent with
# transient phenomena).


# ---------------------------------------------------------------------
# 3.2.1 Empirical pre-screening: do the persistence features
#       discriminate fire from non-fire candidates?
# ---------------------------------------------------------------------
#
# Before launching the experiment, we examine whether `persist_delta`
# and `persist_ratio` actually differ between the three groups
# defined by the E1 / Baseline comparison:
#
#   - 16 polygons rescued by both Baseline and E1 (real fires;
#     sufficient signal even without hotspots).
#   - 14 polygons rescued by Baseline only (real fires that need
#     the hotspot to be recovered; lost in E1).
#   - 912 polygons rescued only by E1 (suspected false positives;
#     spectral change without hotspot evidence).
#
# Pre-screening output (2005 / balanced):
#
#   persist_delta median (n, IQR):
#     16 solid:           -247   (n=16,  IQR -302 to -186)
#     14 lost by E1:      -261   (n=14,  IQR -270 to -249)
#     912 new in E1:      -227   (n=912, IQR -297 to -156, 3 NA)
#
#   persist_ratio median:
#     16 solid:            0.348
#     14 lost by E1:       0.325
#     912 new in E1:       0.450
#
# Reading
# ---------------------------------------------------------------------
#
# The three group medians for `persist_delta` cluster within a
# narrow band (-227 to -261) and their interquartile ranges overlap
# almost entirely. `persist_ratio` shows a slightly larger separation
# (0.45 for suspected false positives versus 0.32-0.35 for real
# fires), but again with substantial overlap.
#
# This pre-screening already suggests that persistence carries
# limited discriminatory power between real fires and the suspected
# false positives in this year. The persistence amplification could
# nonetheless help if XGBoost finds non-trivial, non-mean-based
# splits that exploit subtle interactions. Experiment E2 tests
# this empirically.


# ---------------------------------------------------------------------
# 3.2.2 Setup and execution
# ---------------------------------------------------------------------

# Whitelist: identical to E1 (51 canonical minus 13 hotspot)
canonical <- OtsuFire:::.supervised_feature_cols
HOTSPOT_BLOCK <- c(
  "hotspot_available",
  "hs_in_poly", "hs_in_buffer", "hs_used_n", "hs_min_dist_m",
  "hs_frp_sum", "hs_frp_max", "hs_conf_mean", "hs_hiConf_n",
  "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support",
  "hs_any"
)
E2_whitelist <- setdiff(canonical, HOTSPOT_BLOCK)

# Weights: 5x amplification on the two persistence features
E2_weights <- c(persist_delta = 5.0, persist_ratio = 5.0)

# (Pipeline call omitted; same as E1 but with feature_weights = E2_weights.)


# ---------------------------------------------------------------------
# 3.2.3 OOF performance: indistinguishable from E1
# ---------------------------------------------------------------------
#
# E2 reproduces E1's behaviour at every observable level. The model
# converged in 953 trees (vs 949 in E1, vs 344 in Baseline). The
# meta confirms both arguments were applied:
#
#   feature_whitelist_override_applied: TRUE
#   feature_whitelist_override_n: 38
#   feature_weights_applied: TRUE
#   feature_weights_n_nondefault: 2
#   feature_weights_nondefault: persist_delta=5,persist_ratio=5
#
# The drop / keep / review distributions and the rescued review
# count are essentially identical to E1:
#
#   Drops mediana:   0.0007 (E1: 0.0006)
#   Keeps mediana:   0.9998 (E1: 0.9998)
#   Reviews mediana: 0.9971 (E1: 0.9971)
#
#   Reviews rescatadas:
#     Baseline (with hotspots):           30
#     E1 (no hotspots):                   928
#     E2 (no hotspots + persist 5x):      929
#
# A single polygon difference in 1114 reviews. E2 = E1 operationally.


# ---------------------------------------------------------------------
# 3.2.4 Feature importance: the boosted features remain marginal
# ---------------------------------------------------------------------
#
# Top 10 features in E2 (Gain):
#
#                Feature        Gain
#   1            rbr_p90 0.577954646
#   2         rbr_aw_p90 0.128869714
#   3            rbr_med 0.099019705
#   4      cor_agri_frac 0.093019892
#   5            rbr_p10 0.046370544
#   6    cor_forest_frac 0.022915432
#   7            elev_sd 0.007752607
#   8         rbr_aw_iqr 0.004599378
#   9           elev_iqr 0.002740531
#   10           rbr_iqr 0.002559999
#
# The persistence features (persist_ratio, persist_delta) sit at
# rank 13 and 24 with Gain 0.0022 and 0.0002 respectively.
#
# Did the 5x weight have any effect? It did, at the algorithmic
# level: persist_ratio's Gain rose from 0.0005 in E1 to 0.0022 in E2
# (a 4.4x increase), and persist_delta moved from 0.0004 to 0.0002.
# XGBoost did try to use the boosted features more often during
# tree construction. But their absolute contribution remains
# negligible (combined under 0.3% of total Gain), and the model's
# top-level structure remains identical to E1: spectral change
# magnitude (rbr_p90) refined by land cover context (cor_agri_frac,
# cor_forest_frac).
#
# The interpretation is unambiguous. XGBoost's `feature_weights`
# parameter biases the algorithm to *consider* a feature more often
# during split selection, but it cannot manufacture discriminatory
# information that the feature does not carry. When the underlying
# distributions overlap (as the pre-screening showed for
# persistence), no amount of upweighting recovers separation.


# ---------------------------------------------------------------------
# 3.2.5 Methodological implication
# ---------------------------------------------------------------------
#
# Combining E1 and E2, two complementary strategies for compensating
# the loss of hotspots have now been ruled out empirically:
#
#   - Threshold calibration (E1 sweep): does not separate the
#     overlapping probability distributions.
#   - Explicit amplification of temporal persistence (E2): does
#     not provide additional information beyond what RBR
#     already captures.
#
# In this feature space and for the 2005 / balanced year, hotspots
# provide a non-recoverable specificity layer. Neither probabilistic
# thresholding nor within-year feature engineering reconstitutes
# the discrimination capacity that hotspots add to spectral
# magnitude.
#
# This negative result is itself informative for the manuscript:
# it bounds the operational scope of the supervised stage. For
# years with hotspot data (2000 onwards), the Baseline configuration
# provides the calibrated solution. For pre-MODIS years (1985-1999)
# where hotspot data are unavailable, two strategies remain on the
# table:
#
#   - S3: post-hoc filters on geometric and land-cover features
#     (e.g., area_ha < threshold, cor_agri_frac < threshold).
#     Tested in Section 3.3.
#
#   - S6: cross-year transfer learning. Train the supervised model
#     on hotspot-aware years (2000-2025) and apply it to pre-MODIS
#     years with hotspot features set to NA (i.e., letting XGBoost's
#     missing-value handling do the work). This is a LOYO-framework
#     experiment and is addressed in Chapter 4.
#
# The empirical demonstration that within-year feature engineering
# alone cannot recover Baseline-equivalent specificity is the
# central finding of this experiment. It motivates the cross-year
# transfer-learning approach as the principled next step rather
# than as an alternative explored among many.