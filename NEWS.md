# OtsuFire 0.10.1 (2026-06-14)

## Observability gate — `validate_fire_maps()` only (no methodological contract change)

This is a patch release. It touches `validate_fire_maps()`, its observability
preprocessing, cache robustness, audit/reporting, documentation and tests only.
The deterministic phase, the supervised model, scoring and the `p_burned`
threshold are unchanged, and the 2017 reference metrics are bit-for-bit
identical (TP=450746, FP=52060, FN=241280, TN=37704260, Precision=0.8964611,
Recall=0.6513426, F1=0.7544927, IoU=0.6057714, N_observable=663, N_excluded=14,
detected=431).

* **Frozen temporal observability semantics.** Observability in
  `validate_fire_maps()` is defined temporally at the reference-fire level: a
  reference fire with no valid observable DOY, or whose date is after the last
  observation, is excluded from the reference **before** rasterization, so it
  contributes 0 TP / 0 FN / 0 evaluated area and never enters omission, recall,
  F1, IoU or per-fire detection. It is **not** a pixel-level cloud / valid-data
  mask; partially observed fires are kept whole and no per-pixel denominator is
  introduced.
* **Content-aware reference/observability cache.** The reference cache key now
  folds the CONTENT and methodological identity of every input that shapes the
  cached artifact (observability raster content + band + DOY columns + rule
  version; burnable raster grid / CRS / extent / content; study-area mask;
  burnable thresholding; reference content; min-area / dissolve options). The
  filename tag is `obs-<hash>_dom-<hash>_ref-<hash>`; the prediction cache also
  folds in the burnable-domain identity. A changed input busts the cache
  automatically; `force_reprocess_ref = TRUE` still forces a full rebuild.
  `buffer` is intentionally excluded (it is applied downstream, not cached).
  New internal helpers in `R/internal-validate-cache.R`.
* **Partial-observability audit (informative only).** `validate_fire_maps()`
  now returns `observability_audit` and writes `OBSERVABILITY_AUDIT_<year>.csv`,
  `OBSERVABILITY_EXCLUDED_FIRES_<year>.csv` and
  `OBSERVABILITY_PARTIAL_FIRES_<year>.csv` to `02_OBSERVABILITY`, plus per-fire
  `total_pixels` / `observable_pixels` / `non_observable_pixels` /
  `observable_fraction` columns on `reference_observability`. These outputs
  never change TP/FP/FN/TN or coverage. New helper in
  `R/internal-validate-observability-audit.R`.
* **GPKG overwrite hotfix (deterministic, I/O only).** The three stage-2
  scoring writers propagate `delete_dsn = isTRUE(overwrite)` for a robust GPKG
  replace on Windows; geometries, attributes, counts and deterministic results
  are unchanged.
* **Tests.** Added the 10 mandatory cache-invalidation cases
  (`test-validate-cache-key.R`), the frozen-semantics + audit tests
  (`test-observability-audit.R`) and the GPKG overwrite tests
  (`test-deterministic-gpkg-overwrite-deldsn.R`).

Diagnostic note (not a code change): the 2017 omission is **not** explained by
observability. Of the omitted observable fires, 232 of 236 have no
deterministic candidate (98.3% of omitted fires, 94.3% of omitted area). The
next gate is the deterministic candidate-generation ceiling, addressed
separately.

# OtsuFire 0.10.0 (2026-06-12)

## BREAKING CHANGE — final two-source negative architecture (random + Otsu residual)

* **The supervised negative pool now has exactly TWO sources: random background
  and Otsu residual.** GATE 6 removed the contextual and spectral negative
  buckets entirely; a deterministic drop is no longer a training negative (those
  ambiguous drops stay scoreable but unlabelled). An audit across six years
  showed no contextual category is a reliable negative and the spectral bucket's
  supply is effectively zero in modern data years.
* **`cap_contextual` / `cap_spectral` removed** (G6.1, G6.5), together with
  `contextual_exclusion_to_burned_ratio` / `spectral_to_burned_ratio` and the
  contextual/spectral negative generators. Passing any of them now errors.
* **Otsu residual generator renamed** `legacy_*` → `otsu_negative_*`
  ("Otsu residual negative generator", `R/internal-sup-otsu-negative.R`) across
  the full public + internal surface (G6.6). **Old caches are invalidated; there
  is no migration** — the negative-pool fingerprint changes, so a stale cache is
  rebuilt rather than silently reused.
* **`legacy_sample_n` removed** (G6.2): the Otsu generation-side pre-thinning is
  gone; `cap_otsu` (`otsu_unburned_to_burned_ratio`) is now the SOLE, seeded Otsu
  selector. A new **Otsu > random spatial dedup** runs before capping: a location
  that is both a random background cell and an Otsu patch is kept once, with Otsu
  given priority (`otsu_random_dedup_audit`, folded into the fingerprint).
* **Dead Otsu review/keep parameters removed** (G6.4).
* **Typed public negative-pool API** (G6.7): `negative_pool_params = list(random
  = list(n_cells, rbr_quantile), otsu = list(candidate_threshold,
  reference_threshold), caps = c(random, otsu))` and `runtime_options =
  list(reuse_existing, write_outputs, verbose)`. **Caps are single-source**
  inside `negative_pool_params$caps` (no top-level `cap_*` args; passing them
  errors) and are mirrored into `cfg$train_control$caps`. The negative-pool
  random-background seed `random_seed` lives in the canonical seeds block
  `cfg$train_control$seeds$random_seed`. The runtime toggles are technical only
  and are EXCLUDED from the methodological fingerprint; unknown keys in either
  typed block error. The remaining low-level Otsu-residual knobs (mode
  `burnable_only`, min pixels, buffers, core threshold, boosts, distance power,
  keep/drop confidence, exclusion buffer, min area, `use_drop`,
  `drop_max_s_patch`, `allow_empty_otsu_pool`) are FIXED internal validated
  defaults, not user-settable.
* **Documentation aligned to the final two-source architecture** (G6.8): roxygen,
  `man/*.Rd`, the vignette, the METHODS_BOOK supervised + one-year chapters, the
  `inst/scripts` (incl. the manual and long-run runners) and the in-tree READMEs
  now describe only random background + Otsu residual; every trace of the old
  contextual/spectral/legacy/four-caps API was removed from the public docs.

# OtsuFire 0.9.0 (2026-06-11)

## Methodological fix — supervised eligibility defined by EXPLICIT class

* **The supervised training population is now defined by EXPLICIT class, never by
  negation.** A row is a POSITIVE iff `class == "burned"`; a NEGATIVE iff
  `class == "unburned"` AND it resolves to exactly one of the four valid negative
  buckets (`contextual`, `spectral`, `random`, `otsu`). Otsu review / keep
  patches are EXCLUDED and logged (never trained, never an error). Any other
  case is a strict ERROR: an `unburned` row with `NA` / unknown / "other" /
  empty bucket (and not a known-excluded type), or a row whose `class` is `NA`
  or not in `{burned, unburned}` — each error names the id, class, source,
  neg_type and origin stage (OOF / FINAL). **Review / keep / `NA` / unknown /
  "other" rows can therefore never silently become negatives.** The old OOF
  predicate `is_negative <- !is_burned` and its `sel_other` / `other_kept` 5th
  "other" negative category have been removed.

## Structural — single shared internal capping helper for OOF and FINAL

* The OOF per-fold trainer and the FINAL pool builder now route their negatives
  through ONE shared eligibility resolver and ONE shared capping implementation
  (two PURE internal `@noRd` helpers): `.of_resolve_supervised_eligibility()`
  (eligibility + bucket assignment + audit table + excluded-ids log, upstream of
  capping) and `.of_cap_negative_buckets()` (the single canonical cap:
  `cap = ceiling(n_burned * ratio)`; take-all when available <= cap; all
  positives kept; negatives sampled without replacement in a deterministic
  bucket order; per-bucket audit). This deduplicates the two prior parallel
  capping blocks. The ONLY differences between OOF and FINAL are the retained
  outer-test fold (OOF) and the seed (OOF per-fold fold seed; FINAL
  `cfg$train_control$seeds`). **Bucket definitions, caps, features, thresholds,
  and model params are unchanged.**

## Behaviour-preserving for the default 2017 balanced run

* For the default 2017 balanced run (`use_review = FALSE` / `use_keep = FALSE`),
  the problematic review/keep/unknown/NA set is empty, so counts, policy,
  effective ratios and audit are identical and the OOF selected ids are
  preserved exactly. **Transparency note:** historically OOF and FINAL already
  produced different negative IDs on take-all buckets (FINAL churned the RNG via
  `sample.int` even when taking all; OOF short-circuited). The unified helper
  adopts the deterministic OOF short-circuit, so OOF ids are preserved exactly
  while a future re-run's FINAL negative IDs (not counts) on take-all buckets may
  differ — equivalent to a seed reshuffle, not a methodological change.

# OtsuFire 0.8.0 (2026-06-11)

## BREAKING CHANGE — single capped OOF negative-sampling policy

* **`oof_sampling` removed from the public API.** It is no longer an argument of
  `build_supervised_burned_config()`, `run_oof_diagnostics()`,
  `run_oneyear_supervised_pipeline()` or `validate_supervised_execution()`, and
  the OOF "full" negative-sampling path was deleted. **Passing `oof_sampling` now
  errors** (R's "unused argument" on the plain functions; an explicit guard on
  the `...`-bearing entry point `run_oneyear_supervised_pipeline()`). OOF now
  ALWAYS uses the same capped negative-sampling policy as the FINAL model,
  applied independently within each training fold. The policy survives
  internally only as a fixed, non-settable provenance constant
  (`cfg$train_control$oof_sampling = "capped"`). New per-fold capping audit
  fields record the policy applied to each training fold. **Bucket definitions,
  caps, features, thresholds, and model params are unchanged.**
* In-package scripts and docs updated accordingly: the long-run runners no
  longer build a "capped vs full" OOF comparison — there is one capped OOF
  negative-sampling policy, identical to the FINAL model's, applied per fold.

# OtsuFire 0.7.0 (2026-06-11)

## BREAKING CHANGE — single supervised training procedure

* **`training_protocol` removed from the public API.** It is no longer a formal
  argument of `build_supervised_burned_config()`, `run_oof_diagnostics()`,
  `train_final_burned_model()`, `run_oneyear_supervised_pipeline()` or
  `validate_supervised_execution()`. **Passing `training_protocol` now errors**
  (R's "unused argument" on the plain functions; an explicit guard on the two
  `...`-bearing entry points).
* **The `legacy` training path was deleted.** OtsuFire now always uses ONE
  training procedure: the number of boosting rounds is selected with an inner
  validation split and early stopping, then the recipe and model are refitted on
  all available training observations before prediction. In OOF, each outer fold
  selects `best_iteration` on an inner validation split and refits on all
  outer-train rows before predicting the untouched outer-test fold (the
  outer-test fold never enters imputation, feature selection, levels, weights,
  round selection or refit). OOF and FINAL share one leakage-free core.
* The procedure survives internally only as a fixed, non-settable provenance
  constant (`cfg$train_control$training_protocol = "nested_refit"`); the saved
  recipe records `training_method = "inner_early_stopping_full_refit"`. The
  user-facing `print()` shows `training = early-stopping selection + full-data
  refit`.
* In-package scripts and docs updated accordingly: the protocol-pair example
  scripts (`03_supervised_protocol_legacy.R` /
  `04_supervised_protocol_nested_refit.R`) were removed; the long-run runners no
  longer build "Config A legacy vs B nested" — there is one training procedure.

# OtsuFire 0.6.2

Cleanup + documentation + tooling release. **No methodology, caps, thresholds,
recipe, model, or prediction behaviour changed.** The RESOLVED supervised feature
space (50 base + 50 `_isNA` = 100) and the `feature_schema_fingerprint` are
UNCHANGED, so this release is **NON-BREAKING**: no model, recipe, or scoring
output changes.

## CLEANUP / WHITELIST

* **`hs_any` removed from the canonical supervised whitelist.**
  `.supervised_feature_cols` now lists **50** names (was 51). `hs_any` was a
  redundant derived helper (= `hs_used_n > 0`): it was never materialised as a
  base GPKG feature and never entered the recipes, so the RESOLVED feature space
  (50 base + 50 `_isNA` = 100) and the `feature_schema_fingerprint` are
  **UNCHANGED**. The `hs_any` synthesis block was removed. This is a pure
  cleanup — **no model / recipe / prediction change**.

## DOCS / TESTS

* **Final-map output-layer contract documented and regression-tested.**
  `deterministic_scored` / `final_map_full` carries ALL scored candidates
  (1274 for 2017); `final_map` is the current-year public/thresholded map
  (1273). The 1-row difference is the intentional `current_year_public_drop`
  temporal filter (the current-year public split), **not** a silent loss. The
  contract is now asserted by a dedicated regression test.

## SCRIPTS

* **Canonical versioned long-run runners** under `inst/scripts/long_run/`
  (explicit cfg routes, no junctions; `validate_shared_inputs` fail-fast).
* **Function-by-function manual script**
  `inst/scripts/manual/SUPERVISED_2017_FUNCTION_BY_FUNCTION.R` (FULL /
  NO_HOTSPOT profiles via explicit feature config).
* **Low-memory EFFIS modes** (`STANDARD_90M` default; `NATIVE_HIGH_MEMORY`
  guarded).

# OtsuFire 0.6.1

Patch release. A pure geometry-sanitization **bug fix**: **no methodology,
training parameters, caps, recipe, or configuration changed**. This release only
makes supervised polygon sanitization robust against empty geometries. The defect
was caught by the real-2017 smoke run against the v0.6.0 snapshot (reference HEAD
`bf612a7` / `GATE1_FINAL_v0.6.0_20260609_141533`).

## BUG FIX (geometry sanitization)

* **Supervised polygon sanitize is now empty-geometry-safe.**
  `st_collection_extract(..., "POLYGON")` is applied ONLY to real
  `GEOMETRYCOLLECTION` rows, never to a whole layer already made of
  `POLYGON`/`MULTIPOLYGON`. Previously, under sf 1.0-20, a single empty geometry
  in `internal_decisions` zeroed the entire layer and aborted pool construction
  ("internal_sf is empty"). The fix drops pre-existing empties before
  `make_valid`, repairs, drops post-`make_valid` empties, then extracts polygonal
  parts only from `GEOMETRYCOLLECTION`s (discarding and counting collections with
  no polygonal component), preserving attributes / IDs / order; it errors
  informatively only if nothing usable remains. For 2017: 1278 -> 1274 (4 empties
  dropped).

## VALIDATION

* **`validate_supervised_execution()` gains a non-destructive sanitize dry-run on
  `internal_decisions`.** Empties are reported as `PASS_WITH_SANITIZATION` with
  explicit counts (e.g. 1278 -> 1274); it is blocking only if sanitization would
  leave the layer completely empty.

# OtsuFire 0.6.0

Supervised-closure milestone (Gate 1B–1E). This release unifies the supervised
feature recipe across the out-of-fold (OOF), FINAL-refit, and scoring stages,
restricts the random background to the burnable domain, removes several
abandoned supervised code paths, and adds the public
`validate_supervised_execution()` entry point plus a runtime feature-schema
parity guard. For a 0.x package this is shipped as a MINOR bump (0.5.0 ->
0.6.0); 0.x minors may carry breaking changes, and the items below are breaking.

## BREAKING CHANGES

* **Shared `_isNA` feature recipe (OOF / FINAL / scoring now share ONE recipe).**
  Previously the OOF diagnostics built the 51 `<feature>_isNA` missing-indicator
  companions (a 102-column space: 51 base + 51 indicators) while the FINAL refit
  and the scoring path did NOT, so they ran on a different 51-column space. The
  three stages now build features through a single shared recipe builder; the
  FINAL model and the scoring matrix therefore use the **102-column** space
  (51 base + 51 `_isNA`). **Models trained under the old recipe must be
  retrained** — the deployed feature space changed.
* **B4 `random_background` restricted to the burnable domain.** The random
  background percentile selection and sampling now operate only inside the
  burnable mask (with CRS alignment) and contribute to a content-aware
  `neg_pool_fingerprint`. The **negative pool changed**, so **prior negative-pool
  caches are invalidated** and will be rebuilt.
* **Removed abandoned supervised paths and knobs:**
  * `deterministic_direct` (the pipeline is always the `all_sources` negative
    pool; the user-selectable path no longer exists).
  * the supervised burned-like **registry** (API argument, orchestrator code,
    and cfg field).
  * the `negative_pool_policy` parameter (the policy is implicitly and always
    `all_sources`).
  * the supervised **ecoregion Otsu branch** (`burnable_only` is canonical; the
    deterministic CORINE × ecoregion stratification in the delineation stage is
    preserved — ecoregions are not a supervised feature).
* **Mask / CRS now projected on mismatch.** The burnable-mask alignment helper
  reprojects on a CRS mismatch instead of silently zeroing cross-CRS masks
  (which previously produced empty masks without error).

## DEPRECATIONS

* **Function-level methodological / training-control parameter shims** on the
  supervised entry points (the four `*_to_burned_ratio` caps,
  `feature_whitelist_override`, `feature_weights`, `sampling_seed`, `seed`,
  `val_frac`, `group_col`, `nrounds_max`, `early_stopping_rounds`, `impute_*`,
  `training_protocol`) are DEPRECATED. The canonical single source of truth is
  `build_supervised_burned_config()` (`cfg$model_params` / `cfg$train_control`).
  A non-`NULL` override of a canonical-default field now emits a deprecation
  warning of class `"otsufire_deprecated_param"`; an override that conflicts with
  an EXPLICIT builder value is an error. Removal plan: deprecated now (warn) ->
  scheduled for removal in **0.7.0**, consistent with the roxygen note.

## NEW

* **`validate_supervised_execution()` is now a public export.** Workflow-level
  validator with a structured report, running on the single shared supervised
  engine.
* **`cfg$model_params` / `cfg$train_control` are the single source of truth**
  for supervised parameters, with **per-field provenance** (default vs
  user-set), so partial overrides (e.g. a partial xgboost override) are tracked
  field by field.
* **Runtime feature-schema parity guard.** Gate 1E (2026-06-09) adds a runtime
  guard that prevents the OOF, FINAL-refit, and scoring stages from silently
  diverging in feature space on real data — the runtime backstop for the `_isNA`
  defect. A new internal helper `feature_schema_fingerprint()` computes a
  reproducible, wall-clock-FREE structural fingerprint (base-R rolling checksum,
  no `digest` dependency) that hashes ONLY the structural feature-space contract:
  base-feature names + order, the `_isNA` indicator names + order,
  `final_feature_order`, `feature_cols` / `x_cols`, the `n_base` /
  `n_indicators` / `n_total` counts, the numeric-vs-factor encoding contract, the
  feature-weights POLICY (the rule "`_isNA` indicators carry the canonical weight
  1.0"), and the contract version `.OF_SUPERVISED_FEATURE_CONTRACT_VERSION`
  (`"1.0"`). It DELIBERATELY EXCLUDES imputation medians, learned categorical
  levels, `scale_pos_weight`, `best_iteration`, seeds, any fold-fitted statistic,
  and any timestamp (these legitimately differ between an OOF fold and FINAL). In
  the nested_refit (Phase B) path each outer fold's refit recipe is fingerprinted
  and the guard asserts all folds satisfy the SAME contract (mismatch = stop);
  FINAL asserts its fingerprint equals the canonical OOF contract before training
  or saving and persists it as `recipe$schema_fingerprint`; scoring re-checks
  `final_feature_order`, the column count / names / order, and the saved
  fingerprint before predicting (mismatch = stop). The fingerprints, counts,
  guard result, and contract version are written to a
  `_feature_schema_parity.txt` summary and the run return value; the Phase B run
  ABORTS if OOF, FINAL, and scoring fingerprints are not compatible.
* **Three-level reproducibility contract** (A exact identity / B numeric /
  C cartographic) governing what must be bit-identical versus numerically or
  cartographically equivalent across runs.
* **NA-preserving, recipe-driven scoring.** Scoring no longer drops rows;
  `hotspots = NULL` and all-NA features are supported end to end (missing
  `hs_*` are treated as missing, not as a sentinel value).
* **Semantic / spatial fail-fast validation** of supervised inputs and output
  routes, with cfg isolation guarantees (no hidden fallbacks or convention
  reconstruction).

## FIXES

* **AS02 cache fingerprint.** Fixed the legacy-pool cache fingerprint multiline
  comparison so honest cache reuse works (the cache filenames previously encoded
  only a few parameters; the full parameter fingerprint is now persisted and
  compared).
* **The 9 pre-existing warnings eliminated** (tibble row-names, GDAL shapefile
  warnings, etc.), behaviour-preserving.
* **Fingerprints exclude timestamps.** The fold, negative-pool, and run
  fingerprints are wall-clock-FREE, so reruns reproduce identical fingerprints.
* **`overwrite` default documented as `FALSE`** and honored end to end, including
  partial writes.

## Note on the corrected defect and baseline naming

The `_isNA` divergence above was found and corrected during Gate 1D.8, **before**
the smoke run and **before** the definitive run. As a result the prior OOF
diagnostics did **not** exactly represent the deployed model: OOF evaluated the
102-column space while the deployed FINAL/scoring ran on a 51-column space — the
prior FINAL was a **different model** than the one OOF evaluated. Definitive
results require RETRAINING both models (OOF and FINAL) under the corrected common
recipe. Any comparison baseline produced under the current code must be named the
**"legacy training protocol under the corrected common feature recipe"** — it is
**not** an exact reproduction of the old 51-feature FINAL, because the feature
recipe (and therefore the feature space) changed when the `_isNA` indicators were
unified across the three paths.
