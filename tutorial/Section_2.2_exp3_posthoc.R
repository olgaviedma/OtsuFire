# =====================================================================
# OtsuFire Tutorial — Chapter 3: Custom feature-engineering experiments
# Experiment 3 (E3): Post-hoc CORINE filters on the hotspot-free model
# =====================================================================
#
# Experiments E1 (no hotspots) and E2 (no hotspots + persistence
# amplified) demonstrated that within-year feature engineering
# cannot recover Baseline-equivalent specificity when hotspot
# evidence is unavailable. E1 produced 928 rescued reviews where
# Baseline produced 30, and E2 reproduced E1 essentially without
# difference.
#
# Experiment E3 takes a different approach. Instead of modifying
# the model itself, it asks whether a simple post-hoc rule applied
# to E1's output can filter the false-positive rescues using
# variables that distinguish fire-prone landscapes from
# spectral-change-prone landscapes that are not fire.
#
# E3 is not a model retraining. It is an analysis of the predictions
# already produced by E1, applying landscape-based filters and
# evaluating their precision and recall against the Baseline-rescued
# set treated as ground truth.


# ---------------------------------------------------------------------
# 3.3.1 Empirical pre-screening: do the CORINE features separate
#       fires from suspected false positives?
# ---------------------------------------------------------------------
#
# The 30 polygons rescued by Baseline (validated by hotspot evidence)
# are treated as the reference set of "real fires" for the purposes
# of this analysis. The 912 polygons rescued only by E1 (without
# hotspot inside any of them) are treated as "suspected false
# positives".
#
# Two CORINE features show clearly different distributions between
# the two groups:
#
#   cor_agri_frac:
#     real fires    (n=30):   median 0.50, mean 0.54
#     suspect (n=912):        median 0.02, mean 0.09
#
#   cor_forest_frac:
#     real fires    (n=30):   median 0.21, mean 0.29
#     suspect (n=912):        median 0.91, mean 0.77
#
# Reading
# ---------------------------------------------------------------------
#
# The two groups occupy genuinely different landscape niches:
#
#   - Real fires sit predominantly in agro-forestal mosaics
#     (median ~50% agricultural, ~21% forest). This pattern is
#     consistent with the typical Iberian fire regime: fires
#     are most often anthropogenic in origin, ignited at the
#     wildland-urban or wildland-agricultural interface, and
#     they spread through patchy mosaic landscapes.
#
#   - Suspected false positives sit predominantly in pure forest
#     (median 91% forest, 2% agricultural). These are polygons
#     where the spectral change observed by RBR is real but is
#     not a fire — likely candidates include wind damage,
#     drought stress, insect outbreaks, recent industrial
#     logging, or extreme phenology.
#
# This pattern suggests CORINE-based rules can act as a
# specificity layer post hoc, replacing some of what hotspots
# provided in Baseline.
#
# A methodological caveat. The "ground truth" used here is not
# external (no EFFIS validation). It is the set of polygons that
# Baseline rescued, which themselves are presumed-fires because
# they have hotspot evidence. The analysis is therefore internally
# consistent but not externally validated. EFFIS validation
# remains future work.


# ---------------------------------------------------------------------
# 3.3.2 Four candidate filtering rules
# ---------------------------------------------------------------------
#
# Four rules were evaluated, each combining cor_agri_frac and
# cor_forest_frac in different ways:

# R1: cor_agri_frac > 0.20
# R2: cor_forest_frac < 0.70
# R3: (cor_agri_frac > 0.20) | (cor_forest_frac < 0.50)
# R4: (cor_agri_frac > 0.20) & (cor_forest_frac < 0.50)

# Expected output (2005 / balanced):
#
#                                   Real recovered    Suspects passed
#                                       (n / 30)         (n / 912)
#   R1: agri > 0.20                       22                 155
#   R2: forest < 0.70                     27                 246
#   R3: agri > 0.20 | forest < 0.50       27                 266
#   R4: agri > 0.20 & forest < 0.50       18                  33
#
# Reading
# ---------------------------------------------------------------------
#
# None of the four rules reproduces Baseline-equivalent specificity.
# The trade-offs:
#
#   - R2 and R3 achieve 90% recall (27/30 fires recovered) but
#     still let through 246-266 suspect rescues. Even with these
#     filters, the operational rescue count would be ~270 versus
#     Baseline's 30 — an inflation factor of ~9x.
#
#   - R4 (the strict AND) filters 96% of suspects but drops 12 of
#     the 30 real fires (recall = 60%). This is the closest to
#     Baseline-like restraint but loses too many real fires.
#
#   - R1 (agriculture only) is intermediate but pays a recall
#     cost of 27% for moderate gain in specificity.
#
# For pre-MODIS application, R3 may be defensible if the goal is
# conservative recall, but the user must accept that the
# burned-area area estimate will be ~9x larger than what Baseline
# would produce in equivalent post-MODIS years. The choice between
# R3 and R4 depends on the operational priority: R3 favours
# recall (don't miss real fires), R4 favours precision (don't
# over-attribute fire to non-fire spectral changes).


# ---------------------------------------------------------------------
# 3.3.3 Methodological implication: hotspot specificity is bounded
# ---------------------------------------------------------------------
#
# Combining the three experiments now closes the within-year
# feature-engineering investigation:
#
#   E1 — Probability threshold tuning: no probability cut separates
#        real fires from spectral-change false positives.
#
#   E2 — Explicit persistence amplification: persistence features
#        do not carry independent discriminating information beyond
#        what RBR already captures.
#
#   E3 — Post-hoc CORINE filters: better than E1/E2 but still
#        insufficient. The best recall-preserving rule (R3) leaves
#        ~9x as many rescues as Baseline; the best precision-
#        preserving rule (R4) loses 40% of real fires.
#
# Three internal mechanisms have been ruled out as full substitutes
# for the hotspot signal. The features available within a single
# year (RBR, persistence, CORINE, topography, DOY) cannot, in any
# combination tested, reconstitute the specificity that hotspot
# evidence provides for the supervised stage.
#
# The remaining strategy on the table is cross-year transfer:
# train on years where hotspot evidence exists, predict on years
# where it does not (with hotspot features set to NA, letting
# XGBoost's missing-value handling decide). This requires the
# leave-one-year-out (LOYO) framework rather than the one-year
# one, and is the subject of Chapter 4. The within-year experiments
# of Chapter 3 are therefore complete, with their conclusion being
# negative but principled: they establish the lower bound of what
# pre-MODIS years can achieve with single-year feature engineering
# alone, and motivate the cross-year transfer approach as the
# necessary next step.


# ---------------------------------------------------------------------
# 3.3.4 Note on the absence of EFFIS validation
# ---------------------------------------------------------------------
#
# The conclusion that the 912 E1-only rescues are "false positives"
# rests on the assumption that the 30 Baseline rescues (with
# hotspot evidence) constitute the reference set of real fires.
# This assumption is reasonable but not externally validated.
#
# A possibility worth considering for the manuscript: a non-trivial
# fraction of the 912 E1-only rescues might be real fires that
# happened to evade MODIS hotspot detection (small fires, fires
# under cloud cover, fires occurring between MODIS overpasses).
# In that case, E1 would not be over-rescuing as much as it
# appears — it would be recovering fires that Baseline missed
# because the hotspot signal was absent for non-meteorological
# reasons.
#
# Resolving this ambiguity requires external validation against
# EFFIS or another independent fire perimeter source. Until such
# validation is performed, the precision figures reported in this
# section are upper-bound estimates of the false-positive rate.
# The relative comparison between Baseline and E1/E2/E3, however,
# remains valid: regardless of how many of the 912 are real fires,
# the contrast in operational behaviour between hotspot-aware
# and hotspot-free configurations is the central methodological
# finding of this chapter.