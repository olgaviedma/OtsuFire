# LONG 2017 BALANCED — FULL vs NO-HOTSPOT comparison (OtsuFire 0.6.1)

**Run root:** `C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/00_LONG_RUNS/LONG_2017_BALANCED_20260609_215851`
**Compiled:** 2026-06-10 (READ-ONLY on run outputs; no model/score/map/threshold/whitelist touched).
**Snapshot:** GATE1_FINAL_v0.6.1_20260609_173925 — tarball SHA256 `63dc3be25a4733ab4012883503d396153bbface88826cc59c3299a7c86fbc694` (declared==computed).
**Isolated lib:** `00_SMOKE/SMOKE_2017_v0.6.1_20260609_175335/Rlib` — OtsuFire 0.6.1 asserted. sf 1.0.20 / terra 1.8.29 / xgboost 1.7.9.1 / data.table 1.17.0 (user lib).

> **One-line take:** Both profiles PASS. Removing the 12 thermal/hotspot features in 2017 costs **F1 −0.6 pts (0.739→0.733)** and **IoU −0.7 pts (0.586→0.579)**; the loss is concentrated in **commission/precision (5.08%→7.50%)**, with **omission essentially unchanged (39.51%→39.29%)**. The marginal EFFIS value of the thermal features is **small**. This does **NOT** establish pre-MODIS validity — the balanced OOF labels were partly built with hotspot support, so the comparison is not a clean pre-thermal test.

---

## 1. FULL status — PASS
`FULL_SUCCESS.marker`: `FULL_OK 2026-06-10 09:51:33`. Deployed FINAL = **Config B nested_refit** (capped). **best_iteration = 124**. feature_schema_fingerprint **`f78f1841bd6b6c7da4d418e94dc1fb78`**. scored_rows = 1274. Configs trained: A (legacy OOF+FINAL), B (nested_refit OOF capped + FINAL, DEPLOYED), C (nested_refit OOF full; reuses B FINAL).

## 2. NO-HOTSPOT status — PASS
`NO_HOTSPOT_SUCCESS.marker`: `NO_HOTSPOT_OK 2026-06-10 10:03:50`. FINAL = nested_refit (capped), `use_hotspots=FALSE`. **best_iteration = 166** (FINAL nested_refit audit). fingerprint **`a65b2ecdbb4e27fdb627c8c758168066`**. n_base=38, removed=12, scored_rows=1274. OOF nested capped + FINAL + OOF-full diagnostic (reuses FINAL). No legacy config.

## 3. Timings (MASTER_LOG + RESUME_LOG)
| Stage | Wall time |
|---|---|
| SHARED upstream (pools+folds+features, built ONCE) | 22:16:35 → 23:09:54 ≈ **52.3 min** |
| FULL resume (A/B/C retrain ~3 min + SCORE+MAP+EFFIS) | 09:39:04 → 09:51:33 ≈ **12.5 min** |
| NO-HOTSPOT (train 38-feat + OOF-full + SCORE+MAP+EFFIS) | 09:51:37 → 10:03:50 ≈ **12.2 min** |
| **Total compute** (SHARED + both tails) | **≈ 77 min** |

History: original master run reached SHARED OK then **FAILED at FULL SCORE+MAP** (runner path-wiring bug — SHARED `03_FEATURES/features_geometry.gpkg` not visible at the profile's ENGINE_ROUTES; NOT a package/model/methodology fault). Fixed by junctioning SHARED routes into each profile and resuming only the cheap tail (no retraining of the 52-min SHARED part). See `_ORCHESTRATION/FAILURE_DIAGNOSIS.md`.

## 4. Resources (env_snapshot.txt)
- R 4.3.1 (ucrt); OtsuFire 0.6.1 (isolated lib).
- xgboost **nthread = 6** (FIXED both profiles); OMP_NUM_THREADS=6.
- free disk C: 1027.5 GB; free physical RAM 49.5 GB.
- nrounds_max=4000; early_stop=80; seed=42; caps contextual=0.25 / spectral=2.0 / random=1.0 / otsu=1.0.
- Tarball SHA256 verified (declared==computed). 11 resolved 2017 input SHA256 hashes recorded (internal_decisions, change_index/MinMin, delayed_change_index/Autumn, hotspots, topo, corine raster/strata/lut, burnable masks, peninsula mask, EFFIS reference).

## 5. Feature profiles
- **FULL:** 50 base + 50 `_isNA` = 100 model cols (whitelist size 51 incl. `hs_any` which is dropped/synthesized → 50 resolved base). fingerprint `f78f1841…`.
- **NO-HOTSPOT:** 38 base + 38 `_isNA` = 76 model cols. fingerprint `a65b2ec…`. whitelist override n=38 (provenance=user); 13 dropped at training (12 hotspot-lineage **+ hs_any**).
- **12 removed hotspot-lineage features:** `hotspot_available, hs_conf_mean, hs_frp_max, hs_frp_sum, hs_hiConf_n, hs_in_buffer, hs_in_poly, hs_min_dist_m, hs_no_support_when_available, hs_only_buffer_support, hs_support_present, hs_used_n`.
- The 38 retained: rbr_samewindow (5), rbr_allwindow (5), persistence (2), doy_spectral (5), corine (7), topo elevation (7), topo slope (7). Full lineage in `SHARED/feature_lineage.csv` and `_ORCHESTRATION/FEATURE_LINEAGE.md`.

## 6. OOF metrics (per config) — see `FULL_vs_NOHOTSPOT_OOF.csv`
Computed from each config's own `*_oof.rds` `oof_agg` (532 burned / 14,750 unburned = 15,282 labeled rows).

| Config | AUC ROC | PR-AUC | best-F1 thr | best F1 | F1@0.50 | P@0.50 | R@0.50 | F1@0.25 |
|---|---|---|---|---|---|---|---|---|
| FULL A (legacy) | 1.000 | 0.9986 | 0.93 | 0.9972 | 0.9944 | 0.9907 | 0.9981 | 0.9925 |
| FULL B (nested capped, **DEPLOYED**) | 1.000 | 0.9997 | 0.85 | 0.9962 | 0.9907 | 0.9833 | 0.9981 | 0.9898 |
| FULL C (nested full) | 1.000 | 0.9990 | 0.84 | 0.9962 | 0.9953 | 0.9907 | 1.0000 | 0.9925 |
| NOHS (nested capped, **DEPLOYED**) | 1.000 | 0.9988 | 0.77 | 0.9925 | 0.9897 | 0.9815 | 0.9981 | 0.9907 |
| NOHS (nested full diagnostic) | 1.000 | 0.9990 | 0.96 | 0.9944 | 0.9907 | 0.9815 | 1.0000 | 0.9898 |

> **CAVEAT (critical):** OOF here is **saturated** (AUC=1.0, F1≈0.99) and is **NOT independent skill**. The balanced keep/review OOF labels were partly built with hotspot support (`intersects_hotspot -> keep`), so the OOF over-credits the thermal features and inflates separability for both profiles. Use **EFFIS** (item 7), not OOF, as the operational comparison. The per-profile on-disk `05_OOF/*metrics_summary.txt` reflects only the LAST config written to that path (C for FULL); the table above is config-specific because it is recomputed from each RDS.

## 7. EFFIS metrics (validate_fire_maps, both profiles)
Reference observability filter applied: 652 observable / 11 excluded of 663.

| | FULL | NO-HOTSPOT |
|---|---|---|
| Precision | 0.9492 | 0.9250 |
| Recall | 0.6049 | 0.6071 |
| F1 | 0.7389 | 0.7330 |
| IoU | 0.5859 | 0.5786 |
| Omission % (=1−R) | 39.51 | 39.29 |
| Commission % (=1−P) | 5.08 | 7.50 |
| TP / FP / FN | 418607 / 22403 / 273419 | 420117 / 34084 / 271909 |
| Balanced accuracy | 0.8022 | 0.8031 |

## 8. FULL vs NO-HOTSPOT comparison table — see `FULL_vs_NOHOTSPOT_EFFIS.csv`
| Metric | FULL | NO-HOTSPOT | Δ abs | Δ rel % |
|---|---|---|---|---|
| Precision | 0.9492 | 0.9250 | −0.0242 | −2.55 |
| Recall | 0.6049 | 0.6071 | +0.0022 | +0.36 |
| F1 | 0.7389 | 0.7330 | −0.0059 | −0.79 |
| IoU | 0.5859 | 0.5786 | −0.0073 | −1.25 |
| Omission % | 39.51 | 39.29 | −0.22 | −0.55 |
| Commission % | 5.08 | 7.50 | +2.42 | +47.7 |
| TP (area px) | 418607 | 420117 | +1510 | +0.36 |
| FP (commission area px) | 22403 | 34084 | +11681 | +52.1 |
| FN (omission area px) | 273419 | 271909 | −1510 | −0.55 |
| Polygons kept @0.5 | 753 | 930 | +177 | +23.5 |

**Read:** the entire cost of dropping thermal is **extra commission** (FP +52%): the no-hotspot model keeps +177 more polygons at thr 0.5, slightly raising recall/TP but adding far more false-positive area. F1/IoU move only fractionally.

## 9. Score (p_burned) distribution — see `P_BURNED_DIST_FULL_vs_NOHOTSPOT.csv` (+ `_AREA`)
Both maps: 1274 polygons, joined by `fire_uid`.

| | n | min | P05 | P10 | P25 | median | mean | P90 | max | sd | %>0.25 | %>0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FULL | 1274 | 0.0017 | 0.0056 | 0.0096 | 0.119 | **0.928** | 0.600 | 0.996 | 0.996 | 0.435 | 60.6 | 59.1 |
| NO-HOTSPOT | 1274 | 0.0019 | 0.0070 | 0.044 | **0.455** | **0.992** | 0.754 | 0.998 | 0.998 | 0.379 | 80.0 | 73.0 |

7-bin polygon counts `[0,.05) [.05,.10) [.10,.25) [.25,.50) [.50,.75) [.75,.90) [.90,1]`:
- FULL: 203, 71, 228, 19, 69, 20, **664**
- NOHS: 138, 34, 83, 89, 33, 19, **878**

**Characterisation:** the NO-HOTSPOT distribution shifts **up and gets less bimodal**. FULL is sharply bimodal (mass at the bottom bins + a big top bin, sparse middle). NOHS thins the bottom (203→138 in [0,.05)) and the lower-middle ([.10,.25): 228→83) and piles probability into the very top bin (664→878) and the [.25,.50) shoulder (19→89). The median climbs 0.928→0.992 and P25 jumps 0.119→0.455. Without the decisive thermal "no-support" signal that pushed marginal polygons low, the model becomes **less decisive at the low end and more permissive overall** — hence the 753→930 jump at thr 0.5.

## 10. Classification changes FULL↔NO-HOTSPOT @0.5 — see `CLASSIFICATION_CHANGES.csv` (+ `_byclass`)
Cross-tab (1274 polygons, by `fire_uid`):

| | NOHS burned | NOHS not | total |
|---|---|---|---|
| **FULL burned** | 643 | 110 | 753 |
| **FULL not** | 287 | 234 | 521 |
| total | 930 | 344 | 1274 |

- **Gained (FULL not → NOHS burned): 287 polys, 31,643 ha.** By deterministic class: **review 283 (31,493 ha)**, drop 3 (130 ha), keep 1 (19 ha).
- **Lost (FULL burned → NOHS not): 110 polys, 21,803 ha.** By class: keep 53 (14,628 ha), review 57 (7,174 ha).
- **Net +177** burned = +287 − 110. The **+177 extra NO-HOTSPOT burned polygons are overwhelmingly `review` class** (the deterministically ambiguous middle); the no-hotspot model resolves the review backlog toward "burned". The 110 it loses include 53 confident `keep` polygons that FULL kept above 0.5 — these are real precision losses.
- Burned@0.5 by class: FULL = keep 610 + review 143; NOHS = keep 558 + review 369 + drop 3.

## 11. Feature importance — see `NO_HOTSPOT_feature_importance.csv`, `IMPORTANCE_FULL_vs_NOHOTSPOT.csv`
**FULL FINAL (xgb.importance, Gain):** `hs_in_poly` 0.456 + `hs_min_dist_m` 0.353 **dominate (≈81% of Gain)**, then `rbr_med` 0.159, `rbr_p90` 0.027. The thermal pair carries the FULL model.

**NO-HOTSPOT FINAL (Gain):** signal collapses onto spectral RBR: `rbr_med` **0.704** + `rbr_p90` **0.212** (≈92% combined), then `rbr_iqr` 0.040, `rbr_aw_med` 0.022, `cor_forest_frac` 0.008, `cor_agri_frac` 0.003, `elev_sd`/`slope_sd` ≈0.002. **Without hotspots the burned signal is carried by same-window RBR (median + p90), with minor contributions from all-window RBR, CORINE forest/agri fraction, and topographic SD.** This is methodologically the expected fallback and explains why recall holds (RBR detects the burn scar) while precision drops (no thermal confirmation to suppress spectral false positives).

## 12. Warnings / errors (logs)
Both `*.log` and `*_resume_run.log` scanned. **No unexpected warnings or errors.**
- Expected: `1278 -> 1274` sanitize in SHARED (`[sanitize_polygons] n_input=1278 -> n_output=1274 | empty_before=4`) — 4 empty geometries dropped once at feature build; both profiles score 1274.
- Expected: `Temporal observability filter excluded 11 of 663 reference polygons` (EFFIS, both) — by design (obs_before_reference_doy=6, no_observability_data=5).
- Expected/benign: `package 'sf'/'terra' was built under R version 4.3.3` (build-version note); histogram `1% sample … NA` notes; `st_centroid assumes attributes constant`.
- Operational note (`WARNINGS.txt`, both profiles): the on-disk negative-pool cache fingerprint (CHECKSUM 985390489) did not match the current cfg (350250422) → **pool rebuilt this run** (sidecar overwritten). Expected because the cfg changed vs an older cached run; not an error.
- The original-run FULL failure (route-wiring) is documented and was resolved by the resume; not a model/methodology issue.

## 13. Generated files inventory (key outputs)
**FULL:** `SUMMARY.md`, `EFFIS_metrics.csv`, `FULL_SUCCESS.marker`, `FULL_feature_list.txt`, `WARNINGS.txt`, `full.log`, `full_resume_run.log`, `cfgs_full.rds`, `effis.rds`; models `A_model.rds`, `B_model.rds` (DEPLOYED FINAL), `A/B/C_oof.rds`, `B_scored.rds`; `CONFIG_B/07_FINAL_MODEL/*final_model.rds + *feature_importance.csv + *recipe.rds + *meta.txt + *nested_refit_audit.csv`; `CONFIG_B/08_SCORED/*deterministic_scored.gpkg`; `CONFIG_B/09_FINAL_MAP/*final_map.gpkg` (layers: deterministic_scored, final_map_full[1274], final_map[1273]) + `*burned_like_scored.gpkg` + counts CSVs; `ENGINE_ROUTES/.../05_OOF/*` (oof_agg/long/metrics/best_thresholds/audit).
**NO-HOTSPOT:** `SUMMARY.md`, `EFFIS_metrics.csv`, `NO_HOTSPOT_SUCCESS.marker`, `NO_HOTSPOT_feature_list.txt`, `WARNINGS.txt`, `no_hotspot.log`, `no_hotspot_resume_run.log`, `cfgs_no_hotspot.rds`, `effis.rds`; `B_model.rds` (DEPLOYED FINAL), `B_oof.rds`, `C_oof.rds`, `B_scored.rds`; `CONFIG_B/07_FINAL_MODEL/*` (model+importance+recipe+meta+audit); `CONFIG_B/08_SCORED` + `09_FINAL_MAP` GPKGs + counts; `ENGINE_ROUTES/.../05_OOF/*`.
**SHARED:** `feats.rds`, `folds.rds`, `pools.rds`, `cfg_shared_upstream.rds`, `feature_lineage.csv`, `shared_manifest.csv`, `SHARED_SUCCESS.marker`, `shared_WARNINGS.txt`, logs, `ENGINE_ROUTES/.../03_FEATURES/features_geometry.gpkg`.
**_ORCHESTRATION:** `env_snapshot.txt`, `feature_lineage.csv`/`FEATURE_LINEAGE.md`, `FULL_whitelist_50_base.txt`, `NO_HOTSPOT_whitelist_38_base.txt`, `MASTER_LOG.txt`, `RESUME_LOG.txt`, `RESUME_COMPLETE.marker`, `FAILURE_DIAGNOSIS.md`, runner scripts (00_common/05/10/20/30/99 .R, MASTER_RUN.bat, RESUME_RUN.bat).
**COMPARISON (new, this deliverable):** this report; `FULL_vs_NOHOTSPOT_EFFIS.csv`, `FULL_vs_NOHOTSPOT_OOF.csv`, `P_BURNED_DIST_FULL_vs_NOHOTSPOT.csv` (+`_AREA`), `CLASSIFICATION_CHANGES.csv` (+`_byclass`), `NO_HOTSPOT_feature_importance.csv`, `IMPORTANCE_FULL_vs_NOHOTSPOT.csv`; figures `P_BURNED_DIST_*.png`, `EFFIS_*.png`, `IMPORTANCE_*.png`, `CLASSIFICATION_CHANGES.png`; analysis scripts `_analysis.R`, `_figures.R`, cache `_cache.rds`.

## 14. Interpretation
The marginal EFFIS value of the thermal/hotspot features in 2017 is **small**: F1 −0.6 pts (0.739→0.733), IoU −0.7 pts (0.586→0.579). The effect is **concentrated in commission/precision** (commission 5.08%→7.50%, precision −2.4 pts) and **NOT in recall** (omission ~unchanged, 39.51→39.29; recall +0.2 pts). Mechanistically: thermal features were the FULL model's top-2 by Gain and acted as a **false-positive suppressor**; removing them shifts the score distribution up, resolves the deterministic `review` backlog toward "burned" (+283 review polygons → burned), and trades a little extra true-positive area for substantially more false-positive area. The loss classification at removing thermal in 2017 = **SMALL on F1/IoU, none on recall, MODERATE on commission**.

> **Do NOT conclude pre-MODIS validity.** The OOF labels (balanced keep/review) were hotspot-informed (`intersects_hotspot -> keep`), and EFFIS itself is an active-fire-derived reference. This experiment shows the *model* tolerates thermal removal with small EFFIS cost in a hotspot-rich year; it does **not** test transfer to a genuinely pre-MODIS year with non-hotspot-informed labels.

## 15. Decisions requiring Natalia's review
a. **`hs_any` cleanup** — `hs_any` is in the 51-entry whitelist but dropped at training (the 13th drop in NO-HOTSPOT). Decide: remove it from the canonical whitelist, or synthesize it as a real feature. (Whitelist NOT touched here.)
b. **Operating threshold** — deployed at 0.50. OOF "recommended" (balanced-accuracy) lands ~0.45 and best-F1 ~0.77–0.93; the task's expected ~0.25 over-detects on this saturated OOF. Recommend a proper EFFIS-based threshold sweep rather than trusting the hotspot-informed OOF threshold.
c. **Pre-MODIS transferability follow-up** — the no-hotspot model recovers recall but at +commission; a true pre-MODIS test needs non-hotspot-informed labels and a reference that is not active-fire-derived.
d. **Upstream deterministic candidate-generation limitation** — the ~39% omission is dominated by no-candidate omissions (≈223 fires with no deterministic candidate); thermal removal does not touch this — it is an upstream candidate-generation ceiling, not a supervised-model ceiling.
e. **Versioning** — the long-run runner scripts (`_ORCHESTRATION/*.R/.bat`) and these COMPARISON reports live in `00_LONG_RUNS`, outside the package git repo. Decide whether to version them.
f. **EFFIS-lineage caveat for publication** — hotspot is both a model feature AND (indirectly) the EFFIS active-fire reference lineage; this circularity must be stated explicitly in any publication claim about thermal feature value.
