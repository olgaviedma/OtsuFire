# LONG 2017 BALANCED — FAILURE DIAGNOSIS (FULL profile)

**Status:** FULL failed at SCORE+MAP. NO-HOTSPOT was NOT launched (failure policy honored). SHARED + FULL trained models are intact. NO patch applied, NO relaunch performed — awaiting Natalia's authorization.

## 1. Stage
FULL profile, **SCORE + MAP** step (deploy B nested_refit FINAL), AFTER all three configs trained successfully:
- CONFIG A (legacy OOF + FINAL): done, FINAL best_iteration=124.
- CONFIG B (nested_refit OOF capped + FINAL): done, FINAL best_iteration=124 (the deployed FINAL).
- CONFIG C (nested_refit OOF full, reuses B FINAL): OOF done.
Timeline: MASTER START 22:16:35 -> SHARED OK 23:09:54 (~52 min) -> FULL train ~3 min -> FAILURE at SCORE+MAP 23:12:50.

## 2. Traceback
```
[23:11:51] [full] SCORE + MAP (deploy B nested_refit FINAL) ...
[Gate 1D.2] feature-schema check PASS: recipe schema (100 features) covered by scoring inputs.
Error in score_burnedlike_and_export_final_map(result_dir = result_dir, ...):
  No existe labelled_features_gpkg:
  .../LONG_2017_BALANCED_FULL/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/03_FEATURES/features_geometry.gpkg
Calls: ... score_supervised_burned_map -> score_burnedlike_and_export_final_map
```

## 3. Root cause — RUNNER SCRIPT path-wiring bug (NOT package / model / methodology)
- Features were built ONCE under **SHARED**: `SHARED/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/03_FEATURES/features_geometry.gpkg` (exists, valid; pools+folds also in SHARED).
- The FULL profile's SCORE step reconstructs the features path under the **FULL profile's** ENGINE_ROUTES (its result_dir) -> `.../LONG_2017_BALANCED_FULL/ENGINE_ROUTES/.../03_FEATURES/` which contains only `05_OOF/` (no features gpkg).
- The installed package 0.6.1 behaves CORRECTLY: score_burnedlike_and_export_final_map expects the labelled features GPKG at `<result_dir>/03_FEATURES/`. The bug is purely in the orchestration runner's "shared-upstream-once + per-profile-routes" wiring: the SHARED 01_POOLS/02_FOLDS/03_FEATURES were not made visible at each profile's scoring route. No package code, model, feature space, caps, threshold, seeds, or methodology is involved. Training (which received the features explicitly) succeeded; only the scoring step's path-reconstruction failed.

## 4. cfg (resolved, saved cfgs_full.rds)
training_protocol=nested_refit (B deployed); caps contextual=0.25 / spectral=2.0 / random=1.0 / otsu=1.0; FULL 50 base + 50 _isNA = 100 cols; seeds 42; nthread=6; nrounds_max=4000 / early_stop=80.

## 5. Version
OtsuFire 0.6.1 (isolated lib, verified). Snapshot GATE1_FINAL_v0.6.1_20260609_173925; tarball SHA256 63dc3be25a4733ab4012883503d396153bbface88826cc59c3299a7c86fbc694.

## 6. Hashes
All 11 resolved 2017 inputs hashed in `_ORCHESTRATION/env_snapshot.txt`; SHARED pools/folds/features built from them (SHARED_SUCCESS.marker present).

## 7. Last valid output (INTACT)
- SHARED upstream: pools, folds, features (SHARED/ENGINE_ROUTES/...), SHARED_SUCCESS.marker — INTACT and reusable.
- FULL trained models: `A_model.rds`, `B_model.rds` (deployed FINAL, best_iter=124), `A_oof.rds`, `B_oof.rds`, `C_oof.rds` — INTACT. CONFIG_A/B/C dirs present.
- Not produced: FULL scoring, final map, EFFIS; entire NO-HOTSPOT profile.

## 8. Fix proposal (RUNNER SCRIPT ONLY — no package/model/methodology change)
Before the SCORE step of each profile, make the SHARED engine routes visible at the profile's ENGINE_ROUTES path:
- copy (or create a directory junction `mklink /J`) of `SHARED/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/{01_POOLS,02_FOLDS,03_FEATURES}` into `<profile>/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/`; OR pass the explicit SHARED `features_geometry.gpkg` path to the scoring call if the API accepts it.
Then RELAUNCH only the cheap tail:
- FULL: SCORE + MAP + EFFIS reusing `B_model.rds` (NO retraining; ~10-15 min);
- then NO-HOTSPOT (full profile, ~1-3 h).
SHARED (the 52-min part) is reused untouched. This is an orchestration-path fix; it does not alter the experiment's substance.

**HELD: no patch, no relaunch until Natalia authorizes. A/B/C models + SHARED preserved.**
