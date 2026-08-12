# RUN MANIFEST — LONG 2017 BALANCED (FULL vs NO-HOTSPOT)

Provenance for the versioned audit reports in this folder, and references to the
**heavy run artifacts** that are intentionally **NOT copied** into the package
repo (kept by absolute path + SHA256 so they can be located and integrity-checked
without bloating the repo).

## Run identity

- Run root: `C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/00_LONG_RUNS/LONG_2017_BALANCED_20260609_215851`
- Frozen code: OtsuFire **0.6.1** (isolated lib), snapshot
  `GATE1_FINAL_v0.6.1_20260609_173925`, tarball SHA256
  `63dc3be25a4733ab4012883503d396153bbface88826cc59c3299a7c86fbc694`.
- Caps: contextual 0.25 / spectral 2.0 / random 1.0 / otsu 1.0. Seed 42.
  nrounds_max 4000 / early_stop 80 / nthread 6.
- FULL: 50 base + 50 _isNA = 100 features. NO-HOTSPOT: 38 base + 38 _isNA = 76
  (the 12 hotspot-lineage features removed via `feature_whitelist_override`).

## Versioned reports (copied here — small)

| File | Source |
|------|--------|
| `LONG_2017_FULL_vs_NOHOTSPOT_REPORT.md` | `<run>/COMPARISON/` |
| `ROUTE_FIX_NOTE.md` | `<run>/_ORCHESTRATION/` (G2.3 route-failure diagnosis + fix) |
| `FAILURE_DIAGNOSIS.md` | `<run>/_ORCHESTRATION/` (original SCORE+MAP failure) |
| `FEATURE_LINEAGE.md` | `<run>/_ORCHESTRATION/` |
| `feature_lineage.csv` | `<run>/_ORCHESTRATION/` |
| `FULL_vs_NOHOTSPOT_OOF.csv` | `<run>/COMPARISON/` |
| `FULL_vs_NOHOTSPOT_EFFIS.csv` | `<run>/COMPARISON/` |
| `IMPORTANCE_FULL_vs_NOHOTSPOT.csv` | `<run>/COMPARISON/` |
| `P_BURNED_DIST_FULL_vs_NOHOTSPOT.csv` | `<run>/COMPARISON/` |
| `P_BURNED_DIST_AREA_FULL_vs_NOHOTSPOT.csv` | `<run>/COMPARISON/` |
| `CLASSIFICATION_CHANGES.csv` | `<run>/COMPARISON/` |
| `CLASSIFICATION_CHANGES_byclass.csv` | `<run>/COMPARISON/` |
| `validate_FULL.csv` | `<run>/_ORCHESTRATION/` |
| `validate_NO_HOTSPOT.csv` | `<run>/_ORCHESTRATION/` |

(The audit reports are ALSO versioned at
`C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/00_AUDITS/AUDIT_SCRIPT_SUPERVISADO/LONG_2017/`.)

## Heavy artifacts — referenced by path + SHA256 (NOT copied)

| Artifact | bytes | SHA256 |
|----------|------:|--------|
| `<run>/SHARED/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/03_FEATURES/features_geometry.gpkg` | 47906816 | `6515039625c7c2794bb43d4c809ac8207ab53f2723092e3e30de1ac162e45b18` |
| `<run>/SHARED/pools.rds` | 27789575 | `9473c845a75748a8980d26f06778f27292283fc026b2327e8cd2f5cbb873ec96` |
| `<run>/SHARED/folds.rds` | 7992006 | `5ebe677b812a8b420f796148f43cc4f49e38df157c4c6128532557729cef1c22` |
| `<run>/SHARED/feats.rds` | 16793081 | `bf022470b19ebc8f73c7a28a9b780f33464ca0761b16f05505809e1816c48f9f` |
| `<run>/LONG_2017_BALANCED_FULL/B_model.rds` (deployed FULL FINAL) | 15040443 | `f2905a63fe3a443828a6c359e8b07fdc2fc7bb4b34a72c5d18832b952ed1ce63` |
| `<run>/LONG_2017_BALANCED_NO_HOTSPOT/B_model.rds` (deployed NO-HOTSPOT FINAL) | 15052707 | `51bf86268b7e000e5eb56344a47ae08fe429c73bc66b187434c3d6bc4067395d` |

`<run>` = the run root above. The full per-input SHA256 set is in
`<run>/_ORCHESTRATION/env_snapshot.txt`; per-profile shared-input provenance is
in each profile's `SHARED_INPUTS_PROVENANCE.csv`.
