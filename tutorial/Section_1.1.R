# =====================================================================
# OtsuFire Tutorial — Chapter 1: The Standard Workflow
# Section 1.1: Understanding the supervised pipeline output
# =====================================================================
#
# After `run_oneyear_supervised_pipeline()` finishes, the canonical
# output of the supervised stage lives in `09_FINAL_MAP/`. The most
# important file is `<year>_<scenario>_patch_certified_final_map.gpkg`,
# a GeoPackage that carries three layers — each representing a
# different view of the same scoring universe.
#
# This section explains what each layer contains and when to use it.
#
# Setup

# ---------------------------------------------------------------------
# Convention used in this tutorial
# ---------------------------------------------------------------------
#
# All analytical examples in subsequent sections use `final_map_full`.
# This layer carries the complete column schema (51 model features
# plus deterministic provenance, OOF outputs, and temporal scoring
# metadata) and the full polygon universe (no public-drop filter
# applied). The `final_map` layer is the user-facing product and is
# only relevant when distributing the burned-area output; it is
# discussed in Section 1.X.

# ---------------------------------------------------------------------

library(sf)

# The path produced by run_oneyear_supervised_pipeline() is returned
# in the result object as `result$final_map_gpkg`. For this tutorial
# we point directly to the canonical Baseline output for 2005 / balanced.
final_map_path <- file.path(
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results",
  "2005/Min_Min/SUPERVISED/balanced/09_FINAL_MAP",
  "2005_balanced_patch_certified_final_map.gpkg"
)

# Inspect the layers in the GPKG
sf::st_layers(final_map_path)
#> Driver: GPKG
#> Available layers:
#>             layer_name geometry_type features fields                      crs_name
#> 1 deterministic_scored Multi Polygon     1905     90 ETRS89-extended / LAEA Europe
#> 2       final_map_full Multi Polygon     1905     90 ETRS89-extended / LAEA Europe
#> 3            final_map Multi Polygon     1888     34 ETRS89-extended / LAEA Europe

# The three layers
# ---------------------------------------------------------------------
#
# `deterministic_scored` (1905 × 90).
#   The complete set of polygons the deterministic stage handed to the
#   supervised stage for scoring, augmented with the probability columns
#   (`p_burned`, `p_burned_oof`, `p_burned_model`, `p_burned_current_year`)
#   produced by the supervised model. Each polygon carries its full
#   feature set (the 51 model inputs) plus deterministic provenance
#   columns (`filter_*`, `reason_*`, `class_audited`, `qa_*`).
#   Use this layer when comparing deterministic vs supervised decisions,
#   or when auditing why a particular polygon received a given score.
#
# `final_map_full` (1905 × 90).
#   The same universe and the same column schema as
#   `deterministic_scored`. In the present implementation, the supervised
#   stage does not reassign `class_final` — it only attaches
#   probabilities. The two layers therefore hold byte-equivalent
#   information. The duplication is a forward-compatibility convention:
#   future supervised behaviour that overrides deterministic decisions
#   will diverge between the two layers, with `final_map_full` carrying
#   the supervised-adjusted view.
#   Use this layer for the same purposes as `deterministic_scored`;
#   prefer this name when you want to make the supervised provenance
#   explicit.
#
# `final_map` (1888 × 34).
#   The public-facing burned-area product. This layer is a subset of
#   `final_map_full` along two axes: (i) only the columns relevant for
#   end users are retained — class labels, key probabilities, and
#   deterministic metadata — and (ii) polygons flagged as
#   `current_year_public_drop = TRUE` are excluded.
#   In this run, 17 review-class polygons were excluded by the
#   public-drop filter. All 17 share two characteristics: a
#   `temporal_conflict_flag` (a current-year inconsistency between the
#   patch and its hotspot context) and a low `p_burned_model` (median
#   0.0038), meaning the supervised model already rates them as
#   non-burned. The public layer applies this combined safety filter
#   so that the published burned-area product never carries
#   temporally-inconsistent low-probability candidates.
#   Use this layer for downstream analysis, validation against
#   external references, or when distributing the burned-area product.