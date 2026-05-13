# =====================================================================
# OtsuFire Tutorial — Chapter 3: Custom feature-engineering experiments
# Experiment 4 (E4): Removing hotspots + amplifying CORINE composition
# =====================================================================
#
# Experiments E1 (no hotspots) and E2 (no hotspots + persistence
# boosted 5x) demonstrated that within-year feature engineering
# cannot recover Baseline-equivalent specificity through probability
# tuning or temporal persistence amplification. Section 3.3
# (Experiment E3) then showed that simple post-hoc CORINE filters
# applied to E1 output reduce — but do not eliminate — the
# false-positive inflation: even the best recall-preserving rule
# leaves ~9x as many rescues as Baseline.
#
# A natural follow-up question is whether amplifying CORINE inside
# the model (rather than as a post-hoc filter) would do better.
# The empirical pre-screening of Section 3.3 had revealed that
# CORINE features show a genuinely separating signal between real
# fires (mosaic agro-forestal landscapes) and suspected false
# positives (pure forest landscapes). Unlike persistence — which
# did not separate the groups — CORINE composition does.
#
# Experiment E4 tests this: same hotspot-free whitelist as E1 / E2,
# but with all 7 CORINE features amplified by a factor of 5 via
# `feature_weights`. The hypothesis: forced to consider land-cover
# composition more, the model should learn splits that exploit
# the agro-forestal vs pure-forest distinction empirically
# observed in the pre-screening.


# ---------------------------------------------------------------------
# 3.4.1 Setup and execution
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
E4_whitelist <- setdiff(canonical, HOTSPOT_BLOCK)

# Weights: 5x amplification on all 7 CORINE features
E4_weights <- c(
  cor_agri_frac       = 5.0,
  cor_forest_frac     = 5.0,
  cor_open_frac       = 5.0,
  cor_herbaceous_frac = 5.0,
  cor_urban_frac      = 5.0,
  cor_water_frac      = 5.0,
  cor_wetlands_frac   = 5.0
)

# (Pipeline call omitted; same as E1 but with feature_weights = E4_weights.)


# ---------------------------------------------------------------------
# 3.4.2 OOF performance and behaviour: indistinguishable from E1
# ---------------------------------------------------------------------
#
# E4 again reproduces E1's high-level behaviour despite the
# feature-weight intervention. Tree count grew further:
# best_iteration = 1473 (versus 949 in E1, 953 in E2, 344 in
# Baseline). The model converged after far more iterations,
# suggesting the algorithm searched harder for splits but did
# not find a qualitatively different solution.
#
# Drop / keep / review distributions and rescued reviews are
# essentially identical to E1 / E2:
#
#   Drops mediana:   0.0005 (E1: 0.0006, E2: 0.0007)
#   Keeps mediana:   0.9998 (E1: 0.9998, E2: 0.9998)
#   Reviews mediana: 0.9975 (E1: 0.9971, E2: 0.9971)
#
#   Reviews rescatadas:
#     Baseline (with hotspots):           30
#     E1 (no hotspots):                   928
#     E2 (no hotspots + persistence 5x):  929
#     E4 (no hotspots + CORINE 5x):       931
#
# The three hotspot-free configurations differ by at most three
# polygons across 1114 reviews. From an operational standpoint,
# E4 = E2 = E1.


# ---------------------------------------------------------------------
# 3.4.3 Feature importance: CORINE rises modestly, RBR remains dominant
# ---------------------------------------------------------------------
#
# Top 10 features in E4:
#
#                Feature         Gain
#   1            rbr_p90 0.4725689395
#   2            rbr_med 0.1738024818
#   3         rbr_aw_p90 0.1428638752
#   4      cor_agri_frac 0.1240745634
#   5            rbr_p10 0.0479678757
#   6    cor_forest_frac 0.0108309801
#   7            rbr_iqr 0.0066798299
#   8         rbr_aw_iqr 0.0050157613
#   9           elev_iqr 0.0042020000
#   10           elev_sd 0.0041718374
#
# Comparing CORINE features between E1 and E4:
#
#                  Feature      Gain_E1      Gain_E4    ratio
#   1        cor_agri_frac     0.0884758    0.1240746     1.40
#   2      cor_forest_frac     0.0279894    0.0108310     0.39
#   3  cor_herbaceous_frac     0.0017586    0.0027205     1.55
#   4       cor_open_frac     0.0000650    0.0001063     1.64
#   5      cor_water_frac     0.0007910    0.0006087     0.77
#   6      cor_urban_frac        ---       0.0000165      ---
#
# Reading
# ---------------------------------------------------------------------
#
# The 5x weight produced two distinct effects within the CORINE
# block:
#
#   - `cor_agri_frac` increased from Gain 0.088 to 0.124 (+40%),
#     consistent with the algorithm using it more often.
#     `cor_herbaceous_frac` and `cor_open_frac` also rose,
#     though their absolute Gain remains tiny.
#
#   - `cor_forest_frac` decreased from 0.028 to 0.011 (-61%).
#     The algorithm, asked to prioritise CORINE more, decided
#     that forest fraction was less informative than agricultural
#     fraction — consistent with the pre-screening observation
#     that the suspected false positives are characterised by
#     "not having agriculture", not by "having forest".
#
# Critically, CORINE remains a secondary signal. Total CORINE Gain
# rose from ~12% in E1 to ~14% in E4. Total RBR Gain (same-window +
# all-window) remains around 78% of the model. Even with explicit
# 5x amplification, the algorithm cannot rebalance the model away
# from RBR-dominated splits because the underlying discriminating
# power of RBR exceeds that of CORINE within the same training
# universe.


# ---------------------------------------------------------------------
# 3.4.4 The fundamental limit identified by E4
# ---------------------------------------------------------------------
#
# E4 closes the within-year feature-engineering investigation with
# a complementary finding to those of E2 and E3:
#
#   - Persistence (E2) does not separate fires from non-fire
#     spectral changes because the underlying distributions
#     overlap: persistence values cluster similarly in both
#     groups regardless of weight.
#
#   - CORINE (E4) does separate fires from non-fire spectral
#     changes (the pre-screening confirmed this empirically),
#     but XGBoost cannot exploit this separation against the
#     stronger RBR signal even when explicitly upweighted.
#     Amplifying CORINE 5x increases its consideration during
#     tree construction but does not displace the RBR-dominated
#     decision structure.
#
# Together, E2 and E4 illustrate two distinct limits:
#
#   - When a feature's underlying signal does not separate the
#     groups (E2 case), no amount of upweighting helps.
#
#   - When a feature's underlying signal does separate the groups
#     (E4 case) but a competing feature has stronger separation,
#     upweighting alone is insufficient — the algorithm always
#     prefers the stronger split.
#
# This is consistent with how `feature_weights` is intended to
# work in XGBoost: it biases sampling probabilities at split-time,
# not the algorithm's split-quality decisions. If a competing
# feature offers a higher gain on a given node, it will still
# be selected.
#
# The combined message of Chapter 3 is therefore that within
# the feature space available in OtsuFire — RBR (immediate +
# delayed), CORINE land cover, topography, DOY — there is no
# in-model substitute for the hotspot signal. The specificity
# that hotspots provide for distinguishing fire from spectrally
# similar non-fire change is, within this feature universe,
# non-recoverable.
#
# This finding has direct manuscript implications. For pre-MODIS
# years (1985-1999, before active-fire hotspot data are available),
# operational use of the supervised stage requires either
# external supplementary information (validation against
# independent fire perimeters, additional remote-sensing
# products) or a cross-year transfer-learning approach that
# leverages the patterns learned in hotspot-aware years to
# inform predictions in hotspot-absent years. The latter
# strategy is the subject of Chapter 4.