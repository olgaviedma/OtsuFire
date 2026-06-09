# OtsuFire (development version)

## Corrected defect: OOF / FINAL `_isNA` divergence

A feature-recipe defect was found and corrected during Gate 1D.8, before the
smoke run and before the definitive run:

* Previously, the OOF diagnostics generated the 51 `<feature>_isNA`
  missing-indicator companions, but the FINAL refit and the scoring path did
  NOT use them. OOF therefore evaluated a 102-column feature space (51 base +
  51 indicators) while the deployed FINAL model and scoring ran on a different
  (51-column) space.
* As a result, the prior OOF diagnostics did **not** exactly represent the
  deployed model. The prior FINAL was not necessarily invalid, but it was a
  **different model** than the one OOF evaluated.
* The correction was made **before the smoke and before the definitive 2017
  run**. OOF, FINAL and scoring now share the **same recipe** (one shared
  builder; 51 base + 51 `_isNA` = 102 columns), enforced at runtime by the
  feature-schema parity guard below.
* **Definitive results require RETRAINING both models** (OOF and FINAL) under
  the corrected common recipe.
* Any comparison baseline produced under the current code must be named the
  **"legacy training protocol under the corrected common feature recipe"** — it
  is **not** an exact reproduction of the old 51-feature FINAL, because the
  feature recipe (and therefore the feature space) changed when the `_isNA`
  indicators were unified across the three paths.

## Methodology documentation refresh (Gate 1E.2/1E.3)

The narrative docs (METHODS_BOOK chapters, the workflow-overview vignette) were
updated to match the current code and to remove obsolete content:

* Documented internal validation (OOF, against deterministic-decision labels)
  vs external validation (`validate_fire_maps()` against EFFIS — omission,
  commission, F1, IoU); the two are computed against different references and
  are not comparable.
* Clarified that `p_burned` / `p_burned_current_year` is a MODEL SCORE under the
  training distribution, not a calibrated probability.
* Documented the burnable-only negative pool, the burnable-domain B4 random
  background, the content-aware `neg_pool_fingerprint`, and that CORINE ×
  ecoregion stratification is used only in the deterministic delineation stage
  (not the supervised negative pool; ecoregions are not a supervised feature).
* Documented `training_protocol` (legacy vs nested_refit) as a training-procedure
  choice (not a feature-space change), `cfg$model_params` / `cfg$train_control`
  as the single source of truth, the caps, and `scale_pos_weight`.
* Documented the shared feature recipe / `_isNA` indicators and the runtime
  feature-schema parity guard.
* Removed obsolete content: the rejected/removed supervised burned-like
  registry, `deterministic_direct` as a current path, ecoregions in the
  supervised pool, the standalone-script "bridge" framing (the pipeline is
  package-internal), and references to the removed `negative_pool_policy` knob.

## Runtime feature-schema parity guard

Gate 1E (2026-06-09) adds a RUNTIME feature-schema parity guard that prevents the
out-of-fold (OOF), FINAL-refit, and scoring stages from silently diverging in
feature space on real data. It is the runtime backstop for the Gate 1D.8 `_isNA`
defect (OOF synthesised 51 `_isNA` companions while FINAL synthesised none, so the
deployed model and the OOF diagnostics ran on different feature spaces).

* New internal helper `feature_schema_fingerprint()` computes a reproducible,
  wall-clock-FREE structural fingerprint (base-R rolling checksum, no `digest`
  dependency, mirroring the fold/pool fingerprinters). It hashes ONLY the
  STRUCTURAL feature-space contract: the base-feature names + order, the `_isNA`
  indicator names + order, the `final_feature_order`, the `feature_cols` /
  `x_cols`, the `n_base` / `n_indicators` / `n_total` counts, the
  numeric-vs-factor encoding contract, the feature-weights POLICY (the rule
  "`_isNA` indicators carry the canonical weight 1.0", not a fitted vector), and
  the contract version `.OF_SUPERVISED_FEATURE_CONTRACT_VERSION` (`"1.0"`). It
  DELIBERATELY EXCLUDES imputation medians, learned categorical level values,
  `scale_pos_weight`, `best_iteration`, seeds, any fold-fitted statistic, and any
  timestamp — these legitimately differ between an OOF fold and FINAL because they
  are fit on different training sets. The guard checks the structural CONTRACT,
  not the equality of fitted recipes.

* OOF (nested_refit path): each outer fold's refit recipe is fingerprinted; the
  guard ASSERTS all folds satisfy the SAME structural contract (identical
  fingerprint) and establishes the canonical OOF fingerprint. A discrepancy is an
  ERROR (stop), not a warning. A per-fold fingerprint sidecar CSV is written.

* FINAL: the refit recipe records its structural fingerprint and, BEFORE training
  or saving the FINAL model, ASSERTS it equals the canonical OOF contract (when
  the run carried one). Mismatch is an ERROR. The fingerprint is saved together
  with the model/recipe (`recipe$schema_fingerprint`).

* Scoring: before predicting, the scoring matrix produced from the saved FINAL
  recipe is checked against `final_feature_order`, the column count, the column
  names + order, and the SAVED structural fingerprint. Any discrepancy ABORTS
  (stop) with a clear message — the schema expected by the FINAL model must equal
  the schema produced during scoring.

* Run manifest: the package writes the OOF per-fold fingerprints, the canonical
  OOF fingerprint, the FINAL fingerprint, the scoring-matrix fingerprint, the
  `n_base` / `n_indicators` / `n_total` counts, the guard result (pass/abort), and
  the contract version into a `_feature_schema_parity.txt` summary and the run
  return value. For the Phase B (nested_refit) path the RUN ABORTS if the OOF,
  FINAL and scoring structural fingerprints are not compatible. The guard logic
  lives in the package, so the runner inherits it.
