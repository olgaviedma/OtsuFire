# OtsuFire (development version)

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
