# =====================================================================
# OtsuFire Tutorial — Chapter 3: Custom feature-engineering experiments
# Experiment 1 (E1): Removing the hotspot block
# =====================================================================
#
# Chapter 1 established that the canonical Baseline supervised model
# is overwhelmingly hotspot-driven (98.5% of feature Gain in the
# HOTSPOTS family — Section 1.4). The natural follow-up question is
# whether the model can classify burned vs unburned without access
# to hotspot features at all.
#
# This question is not academic. Natalia's PhD covers 1985-2025, but
# MODIS hotspots only exist from 2000 onwards (TERRA from late 1999,
# AQUA from 2002). Fifteen years (1985-1999) cannot use the
# hotspot-dependent Baseline configuration. Whether a hotspot-free
# model is operationally viable for those years is therefore a
# central methodological concern of the manuscript.
#
# Experiment E1 removes all 13 features of the HOTSPOTS family (G6)
# from the model's feature space and otherwise keeps everything
# identical to Baseline (same pools, folds, training-test partition,
# hyperparameters, ratios). The OtsuFire 0.5.0 API supports this
# directly via `feature_whitelist_override`.

library(sf)


# ---------------------------------------------------------------------
# 3.1.1 Setup and execution
# ---------------------------------------------------------------------

PKG_ROOT <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild"
suppressMessages(pkgload::load_all(PKG_ROOT, quiet = TRUE))

# Define the override: canonical 51 features minus the 13 of G6
canonical <- OtsuFire:::.supervised_feature_cols
HOTSPOT_BLOCK <- c(
  "hotspot_available",
  "hs_in_poly", "hs_in_buffer", "hs_used_n", "hs_min_dist_m",
  "hs_frp_sum", "hs_frp_max", "hs_conf_mean", "hs_hiConf_n",
  "hs_support_present", "hs_no_support_when_available",
  "hs_only_buffer_support",
  "hs_any"
)
E1_whitelist <- setdiff(canonical, HOTSPOT_BLOCK)

# Run E1 with reuse_upstream = TRUE so the only difference vs
# Baseline is the active feature space.
# (Code for the full pipeline call omitted here; see the runner
# script. The key argument is `feature_whitelist_override =
# E1_whitelist`.)


# ---------------------------------------------------------------------
# 3.1.2 OOF performance: the model still classifies, but differently
# ---------------------------------------------------------------------
#
# E1 was launched on 2005 / balanced with `reuse_upstream = TRUE`.
# Inspecting the OOF metrics summary written by the package:

# (Reading from the snapshot at
#  _PHASE_B_RESULTS/2005_balanced/E1_no_hotspots/05_OOF/.)

# Expected output (2005 / balanced, OtsuFire 0.5.0):
#
#   feature_whitelist_override_applied: TRUE
#   feature_whitelist_override_n: 38
#   feature_whitelist_dropped_n: 13
#   [recommended_threshold]
#   criterion: balanced_accuracy
#   threshold: 0.15            # vs 0.85 in Baseline
#   accuracy: 0.998526
#   precision: 0.959364
#   recall: 0.998162
#   specificity: 0.998539
#   balanced_accuracy: 0.99835
#   f1: 0.978378
#
# Reading the output
# ---------------------------------------------------------------------
#
# The model's OOF performance remains very high. F1 dropped from
# 0.999 (Baseline) to 0.978 (E1) — a relative decrease of about 2%.
# Recall (sensitivity) is essentially preserved at 99.8%. The model
# is still capable of discriminating burned from unburned in the
# training universe.
#
# What changed substantially is the optimal decision threshold:
# 0.15 in E1 versus 0.85 in Baseline. This shift signals that the
# probability distributions of the two classes are no longer cleanly
# separated. In Baseline the gap between the lowest burned p_oof
# (0.879) and the highest unburned p_oof (0.867) was wide. In E1,
# unburned polygons spill upward into the 0.1-0.5 range, requiring
# a lower decision boundary.


# ---------------------------------------------------------------------
# 3.1.3 Feature importance: RBR and CORINE take over
# ---------------------------------------------------------------------
#
# Without hotspots, the tree builder redistributes its splits across
# the spectral and contextual features. The top 10 features and the
# Gain by family for E1:

# Expected output (top 10):
#
#                Feature         Gain
#   1            rbr_p90 0.5341822673
#   2         rbr_aw_p90 0.2211570897
#   3      cor_agri_frac 0.0884757662
#   4            rbr_med 0.0632891384
#   5            rbr_p10 0.0281703628
#   6    cor_forest_frac 0.0279894340
#   7            rbr_iqr 0.0065567639
#   8            elev_sd 0.0059253010
#   9           elev_iqr 0.0042309215
#   10        rbr_aw_iqr 0.0041522928
#
# Gain by family (E1):
#
#         family         Gain         pct
#   1   G1_RBR_SW   0.6321         63.2%
#   2   G2_RBR_AW   0.2484         24.8%
#   3   G3_CORINE   0.1182         11.8%
#   4     G4_TOPO   0.0102          1.0%
#   5      G5_DOY   0.0001          ~0%
#   6 G6_HOTSPOTS   0.0000          0%
#
# The Gain has redistributed sharply away from hotspots:
#
#   - RBR same-window now carries 63% of the Gain (from 1% in
#     Baseline). `rbr_p90` alone accounts for over half.
#   - RBR all-window contributes 25%, with `rbr_aw_p90` as the
#     second-most-important feature.
#   - CORINE land cover (agricultural and forest fractions) jumps
#     to 12%.
#   - Topography contributes 1%, DOY is negligible.
#
# 29 of the 38 available features are used in at least one split
# (versus 18 of 50 in Baseline). The model relies on a broader
# feature combination, but the dominant signal is now spectral
# change magnitude (`rbr_p90`) refined by land cover context.
#
# More telling is the tree count: best_iteration = 949 in E1
# versus 344 in Baseline. The model needed almost three times as
# many trees to converge. This is the empirical signature of a
# genuinely harder learning task — the signal is real but more
# distributed.


# ---------------------------------------------------------------------
# 3.1.4 Reviews rescued: the critical operational difference
# ---------------------------------------------------------------------
#
# OOF metrics describe how the model classifies its own training
# universe. The operational question is what happens when the
# trained model is applied to the scoring set — the polygons the
# model will encounter in production. In OtsuFire, the most
# revealing comparison is the rescued reviews: polygons the
# deterministic stage classified as ambiguous and the supervised
# stage promotes to burned.

fm_baseline <- sf::read_sf(
  file.path("...", "_PHASE_B_RESULTS/2005_balanced/Baseline",
            "09_FINAL_MAP/2005_balanced_patch_certified_final_map.gpkg"),
  layer = "final_map_full"
)
fm_e1 <- sf::read_sf(
  file.path("...", "_PHASE_B_RESULTS/2005_balanced/E1_no_hotspots",
            "09_FINAL_MAP/2005_balanced_patch_certified_final_map.gpkg"),
  layer = "final_map_full"
)

reviews_b  <- fm_baseline[fm_baseline$class_final == "review", ]
reviews_e1 <- fm_e1[fm_e1$class_final == "review", ]

rescued_b  <- reviews_b$poly_id[reviews_b$p_burned_model > 0.5]
rescued_e1 <- reviews_e1$poly_id[reviews_e1$p_burned_model > 0.5]

# Expected output (2005 / balanced):
#
#   Rescatados Baseline:  30
#   Rescatados E1:       928
#
#   In both (consensus):   16
#   Only Baseline (lost):  14
#   Only E1 (new):        912
#
# Reading the output
# ---------------------------------------------------------------------
#
# Removing hotspots transforms the rescue behaviour qualitatively,
# not just quantitatively:
#
#   - 16 polygons are rescued by both models. These are the
#     "consensus rescues": large patches with extreme RBR (median
#     rbr_med ~400, max rbr_med 718) and small hotspot evidence
#     (median hs_in_poly = 1). For these cases, the spectral
#     signal alone is strong enough to drive a confident burned
#     classification, with or without the hotspot.
#
#   - 14 polygons are rescued by Baseline only. E1 drops their
#     probability from p_burned ≈ 0.997 (Baseline median) to
#     p_burned ≈ 0.008 (E1 median; max 0.155). All 14 had hotspots
#     in Baseline (hs_in_poly: 9 had 1, 4 had 2, 1 had 6) but
#     moderate RBR (median 404). These are operationally relevant
#     fires that the spectral signal alone cannot recover. Without
#     hotspot evidence to confirm "this elevated RBR is fire", the
#     model treats them as unburned.
#
#   - 912 polygons are rescued by E1 only — none of them had a
#     hotspot inside in the Baseline data. RBR characteristics:
#     median rbr_med 402, max 717. Area characteristics: median
#     64 ha, max 1357 ha. These are polygons with strong spectral
#     change but no fire-specific confirmation. Almost all are
#     likely false positives — spectral changes that are not fire:
#     deforestation, intensive agriculture cycles, hydric damage,
#     extreme phenology, etc.
#
# The 912-fold inflation of rescued reviews is the central finding
# of this experiment.


# ---------------------------------------------------------------------
# 3.1.5 Threshold tuning does not solve the problem
# ---------------------------------------------------------------------
#
# A natural reaction is to ask whether a stricter probability
# threshold would filter out the false positives. The sweep below
# tests thresholds from 0.5 to 0.995.

# Expected output:
#
#    Thresh |  n_E1 | n_BL30 |   common |  lost |   new
#   --------+-------+--------+----------+-------+------
#   p>0.500 |   928 |     30 |       16 |    14 |   912
#   p>0.900 |   830 |     30 |       15 |    15 |   815
#   p>0.950 |   790 |     30 |       13 |    17 |   777
#   p>0.990 |   669 |     30 |        8 |    22 |   661
#   p>0.995 |   608 |     30 |        8 |    22 |   600
#
#   Fraction of E1 rescues that have hotspot in Baseline data:
#     p > 0.50:  928 rescatados |  16 con hotspot (1.7%)
#     p > 0.90:  830 rescatados |  15 con hotspot (1.8%)
#     p > 0.95:  790 rescatados |  13 con hotspot (1.6%)
#     p > 0.99:  669 rescatados |   8 con hotspot (1.2%)
#
# Reading the output
# ---------------------------------------------------------------------
#
# Raising the threshold from 0.5 to 0.99 reduces E1 rescues from
# 928 to 669, a modest 28% decrease. Crucially, the fraction of
# rescues that coincide with hotspot evidence does not improve —
# it actually decreases slightly (1.7% → 1.2%). The threshold sweep
# discards both real fires and false positives in the same
# proportion: there is no probability range in which the model's
# rescues become predominantly hotspot-supported.
#
# At p > 0.99, the model has already lost 22 of the 30 Baseline
# rescues while still proposing 661 new (likely-false-positive)
# rescues that have no hotspot support.
#
# Methodological conclusion. Threshold calibration is insufficient
# to recover Baseline-quality specificity from the hotspot-free
# E1 model. The two probability distributions (real fires and
# spectral-change false positives) overlap substantially in E1's
# probability output. No probabilistic decision boundary separates
# them cleanly.


# ---------------------------------------------------------------------
# 3.1.6 Methodological implication and the question for E2
# ---------------------------------------------------------------------
#
# Hotspots and RBR are not redundant. They carry information of
# different kinds:
#
#   - RBR captures spectral-change magnitude. Strong RBR change
#     indicates that something happened in the polygon between
#     the reference and target periods. It does not, by itself,
#     identify the type of change.
#
#   - Hotspots capture thermal anomaly attribution. Their presence
#     specifies that the spectral change observed is due to active
#     combustion at the time of acquisition, distinguishing fire
#     from other change types (deforestation, agriculture, hydric
#     events).
#
# In the Baseline model, hotspot features dominate the Gain
# precisely because they provide this specificity layer that RBR
# cannot. Removing them collapses the model to a generic
# spectral-change detector, which is sensitive but not specific
# to fire.
#
# This raises a clear question for the next experiment: can the
# model's lost specificity be recovered without hotspots, by
# leveraging features that distinguish fire from other forms of
# spectral change through complementary mechanisms?
#
# A natural candidate is temporal persistence. Burned scars
# typically persist 2-3 years before vegetation recovery; many
# non-fire spectral changes (agricultural cycles, ephemeral
# hydric events) reverse within one year. The OtsuFire feature
# space carries two persistence features in the RBR all-window
# block (G2): `persist_delta` and `persist_ratio`. In the E1
# model, these features were available but received negligible
# weight (Gain together < 0.01). Whether explicit amplification
# of persistence — via the `feature_weights` mechanism — can
# restore Baseline-comparable specificity is the question
# addressed by Experiment E2 (Section 3.2).