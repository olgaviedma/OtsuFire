# =============================================================================
# OtsuFire - PACKAGE DESCRIPTION
# =============================================================================
#
# Contents
# --------
#   1. Introduction
#   2. OtsuFire
#       2.1 unsupervised phase
#       2.2 Supervised phase
#   3. Inputs
#
# =============================================================================


# =============================================================================
# 1. INTRODUCTION
# =============================================================================
#
# This tutorial introduces the OtsuFire framework through the
# accompanying R package of the same name. The package provides a
# reference implementation of the framework, and the tutorial uses it
# as a practical entry point: each conceptual step of the workflow is
# explained alongside the package function that realises it, so that
# users learn the method and the tool together.
#
# OtsuFire targets annual burned-area mapping from medium-resolution
# optical imagery in settings where externally compiled training
# labels are sparse, inconsistent or unavailable. Its design priority
# is interpretability. Every candidate burned patch produced by the
# workflow carries the chain of decisions that led to its
# classification, which allows users to trace, audit and reproduce
# each mapping outcome. The framework is therefore positioned within
# the family of intrinsically interpretable GeoAI approaches, in
# which interpretability is built into the mapping architecture
# rather than added afterwards through post-hoc explanations.
#
# The following sections describe the framework conceptually
# (Section 2) and summarise the input information it relies on
# (Section 3). The remainder of the tutorial then opens up the
# package step by step, illustrating how each component of the
# framework is configured and executed in practice.


# =============================================================================
# 2. OtsuFire
# =============================================================================
#
# OtsuFire is a hybrid, two-stage framework for annual burned-area
# mapping. It is built on the premise that burned-area cartography
# benefits from combining an unsupervised stage with a transparent rule-based core with a
# supervised stage layer, rather than from either approach in
# isolation. The rule-based core delineates a candidate-patch
# universe from spectral-change information using Otsu thresholding
# and explicit spatial rules. The supervised stage, when applied,
# refines this candidate universe by assigning a burned probability
# to each patch using a supervised classifier trained on labels
# derived from the rule-based core itself.
#
# The two stages are hierarchical rather than interchangeable. The
# unsupervised phase produces the candidate universe and the
# auditable evidence supporting each polygon-level decision. The
# supervised phase operates on that same universe, re-ranking
# patches according to their probability of being burned, but never
# redefining the set of candidates. As a result, the supervised
# layer remains dependent on, and constrained by, the spatial
# decisions of the unsupervised stage. This separation allows the
# framework to combine the traceability of explicit rules with the
# discriminative power of supervised learning, without inheriting
# the full dependence of conventional supervised pipelines on
# externally compiled training data.
#
# The framework was developed and calibrated for the Iberian
# Peninsula using the Relativised Burn Ratio (RBR) as the primary
# spectral-change index, but its architecture is index-agnostic and
# transferable. Application to other regions, sensors or
# spectral-change metrics requires recalibration of vegetation stratum-specific
# thresholds, growth parameters and lower bounds, but does not
# require changes to the structure of the workflow.


# -----------------------------------------------------------------------------
# 2.1 Unsupervised stage
# -----------------------------------------------------------------------------
#
# The unsupervised stage transforms an annual spectral-change
# composite into an auditable layer of candidate burned patches,
# each assigned to one of three classes that summarise the strength
# of the evidence supporting it.
#
# The stage begins with stratified Otsu thresholding. The RBR
# distribution is partitioned into strata defined by the
# intersection of vegetation classes and ecoregions, so that
# thresholds reflect the local background reflectance, land cover
# and fire effects rather than a global histogram. For each
# stratum, Otsu's method selects the threshold that maximises
# between-class variance in the change-index distribution, and this
# stratum-specific threshold is then constrained by global and
# class-specific lower bounds to prevent over-permissive detection
# in spectrally noisy strata.
#
# Pixels above the resulting seed threshold are treated as
# high-confidence burned seeds and are subsequently expanded into
# contiguous patches through region growing under a more permissive
# class-specific growth threshold. Patches are retained only when
# they contain a minimum number of seed pixels, which ensures that
# every mapped patch is anchored to a robust spectral core and
# limits the influence of isolated false positives. The result is
# an initial raster of seed-grown candidate patches.
#
# Each candidate patch is then evaluated by a sequence of
# rule-based filters. The first filter assesses internal support,
# including seed evidence, minimum patch plausibility, dominance of
# burnable land-cover classes and the absence of excessive missing
# or non-burnable cover. The second filter evaluates external
# corroboration using active-fire hotspots when reliable hotspot
# detections are available for the target year; when they are not,
# this filter is recorded as not applied, so
# that pre-hotspot years are not penalised for the absence of an
# input that does not exist for them. The third filter assesses
# spectral consistency by comparing each candidate patch against a
# burned-like reference distribution derived from high-confidence
# patches. Together, these filters integrate into a canonical
# unsupervised decision layer that assigns every patch to one of
# three classes: keep, indicating high-confidence burned candidates
# with strong internal or external support; drop, indicating likely
# unburned areas or false positives; and review, indicating
# ambiguous cases often associated with haze, shadows,
# drought-stressed herbaceous systems, land-water adjacency effects
# or other persistent spectral confounders.
#
# All filter outcomes, support metrics and stratification metadata
# are stored as polygon-level attributes. The full decision path of
# every filter is logged at the patch level, so that the reasons
# behind each keep, drop or review assignment can be reconstructed
# afterwards. This is what makes the sule based unsupervised phase
# auditable: the candidate-patch layer is not only a map but also a
# record of the spatial reasoning that produced it.


# -----------------------------------------------------------------------------
# 2.2 Supervised phase
# -----------------------------------------------------------------------------
#
# The supervised phase refines the candidate-patch universe
# produced by the unsupervised phase by assigning a burned
# probability to each patch. Its purpose is to reduce residual
# omission and commission errors in spectrally ambiguous settings
# without redefining the candidate universe and without depending
# on externally compiled training labels.
#
# Training pools are derived directly from the unsupervised
# output. High-confidence keep patches form the burned training
# pool, while the unburned training pool is assembled from
# complementary negative evidence: drop patches rejected by the
# rule-based filters, random background samples drawn from
# burnable areas subject to exclusion buffers around keep patches
# and previously burned areas, and lower-confidence Otsu-derived
# candidates. The contribution of each source is capped relative
# to the burned pool to prevent class imbalance from biasing the
# classifier. Review patches are withheld from training, since
# they represent the cases for which unsupervised evidence is
# ambiguous, and are scored only at inference time.
#
# Patch-level predictors are aggregated from pixel-level
# information within each candidate patch. They include spectral
# summaries derived from the immediate and, where applicable,
# delayed change-index composites; contextual descriptors
# capturing vegetation composition, ecoregion membership,
# elevation and slope; and active-fire variables describing the
# presence, intensity and spatial configuration of hotspot
# detections when these are available. Day-of-year metrics
# associated with the selected change-index observations are
# included as temporal predictors.
#
# The supervised classifier is implemented as a gradient-boosted
# tree model (XGBoost) and trained under spatial block
# cross-validation to mitigate spatial autocorrelation between
# training and validation splits. The framework supports three
# complementary modelling variants. The one-year variant trains
# and evaluates within the same calendar year and provides the
# operational within-year baseline. The out-of-fold variant uses
# spatial cross-validation within a single year to produce
# predictions for which each patch has been scored by a model that
# did not see it during training. The leave-one-year-out variant
# trains on the remaining years and predicts on the held-out year,
# providing the strict cross-year transfer estimate.
#
# Because the unsupervised phase has already filtered out most
# non-burned candidates, the supervised stage operates on a
# pre-curated universe in which the majority of patches are
# already plausible burned candidates. Predicted probabilities
# therefore tend to concentrate close to one, with the operational
# decision threshold also lying closer to one than to the
# conventional 0.5 cut-off. The supervised stage in OtsuFire is
# accordingly best understood as a re-ranker over a pre-curated
# candidate set rather than as an independent burned-area
# classifier.


# =============================================================================
# 3. INPUTS
# =============================================================================
#
# OtsuFire derives patch-level information from six conceptual
# input groups that jointly support the unsupervised filters and
# the supervised predictors. These groups are summarised below;
# their formal definition, sources and detailed variable lists are
# provided in the Supplementary Material of the OtsuFire paper.
#
# Immediate change-index composites are derived from
# medium-resolution optical imagery (Landsat for earlier years
# and Sentinel-2 from 2018 onwards), restricted to the summer
# fire season (1 June to 30 September). For each target year,
# pre- and post-fire composites are obtained by selecting the
# minimum NBR observation within predefined temporal windows, and
# RBR is computed from these temporally matched composites as the
# primary measure of spectral change. The framework is
# index-agnostic; RBR is used as the calibration reference for
# the implementation described here.
#
# Delayed change-index composites are derived from autumn-winter
# observations following the target fire season. They capture
# persistence of the spectral change beyond the immediate
# post-fire period and help to separate true burned areas from
# transient signals such as crop harvesting or seasonal drying.
#
# Land-cover information is taken from CORINE Land Cover
# products, reclassified into eleven ecologically relevant
# classes. These classes underpin three workflow components: the
# burnable mask that constrains both candidate generation and
# validation, the vegetation strata used in combination with
# ecoregions to stratify Otsu thresholding, and the categorical
# vegetation predictors used in the supervised phase. Ecoregions
# are derived from the WWF Terrestrial Ecoregions of the World
# and are intersected with the vegetation classes to define
# ecologically homogeneous units for stratification.
#
# Topographic information comprises elevation and slope derived
# from a 30-metre digital elevation model. Topography is used as
# a set of continuous predictors in the supervised phase, where it
# helps account for systematic variations in illumination,
# moisture availability and vegetation structure, but it is not
# used to stratify the unsupervised phase.
#
# Active-fire information is obtained from satellite hotspot
# detections, available reliably from 1995 onwards. Hotspot
# variables describe temporal availability, spatial presence,
# intensity statistics and contextual signals at the patch level,
# and play two complementary roles in the workflow. In the
# unsupervised phase, they contribute to the external
# corroboration filter, providing additional evidence for
# candidate patches intersected by, or located near, retained
# detections. In the supervised phase, they enter as patch-level
# predictors alongside spectral, contextual and temporal
# variables. Because hotspot availability varies across the
# historical period, a temporal-availability flag is used to
# distinguish years with reliable detections from periods
# predating the products, allowing the workflow to remain
# applicable to the full historical time series.
#
# Day-of-year metrics summarise the acquisition dates of the
# observations selected for the immediate and delayed
# change-index composites. They are included as temporal
# predictors in the supervised phase.
#
# Reference burned-area perimeters from official national,
# regional and pan-European sources (EFFIS) are used exclusively for
# validation. They are treated as independent external benchmarks
# rather than as training labels, and are filtered to the summer
# fire season and to burnable land-cover classes to ensure
# comparability with the mapped outputs.
