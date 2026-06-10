# LONG 2017 BALANCED — ROUTE PROPAGATION FIX (G2.3)

**Scope:** runner-script only (orchestration). NO OtsuFire package change, NO model /
methodology / feature-space / caps / seeds / threshold change. Frozen code remains
installed OtsuFire 0.6.1 (isolated lib, version-pinned). Nothing heavy was re-run.

## 1. Original failure
The paired LONG run built the upstream products (pools, folds, **features GPKG**)
ONCE under `SHARED`:

```
SHARED/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/03_FEATURES/features_geometry.gpkg
```

Each profile (FULL, NO-HOTSPOT) then built its OWN cfg with a DIFFERENT
`output_dir` (its own `ENGINE_ROUTES`). The SCORE step
(`score_supervised_burned_map` -> `score_burnedlike_and_export_final_map`,
package `R/internal-sup-final-map.R`) resolves the **labelled-features** GPKG from
`config$output_routes$features_geometry_gpkg` when `labelled_features = NULL`. That
route reconstructs to the **profile's** `ENGINE_ROUTES`, NOT to `SHARED`. The
features GPKG was never written there, so FULL failed at SCORE+MAP:

```
Error: No existe labelled_features_gpkg:
  .../LONG_2017_BALANCED_FULL/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/03_FEATURES/features_geometry.gpkg
```

The package behaves correctly; the bug was purely in the runner's
"shared-upstream-once + per-profile-routes" wiring (the SHARED upstream was not
made visible at each profile's scoring route).

## 2. Workaround used by the resume (NOT canonical)
The resume created **8 Windows directory junctions** (`mklink /J`), 4 per profile:

```
<profile>/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/{01_POOLS,02_FOLDS,03_FEATURES,_LEGACY_UNBURNED}
  -> SHARED/ENGINE_ROUTES/2017/Min_Min/SUPERVISED/balanced/{...}
```

so the convention path resolved through the reparse point. Junctions are an
out-of-band filesystem trick: the canonical path must NOT depend on them.

## 3. Canonical fix (explicit cfg routing + fail-fast validation)
**API surface (no package change needed):** `score_supervised_burned_map()` already
accepts an explicit `labelled_features` argument (package
`R/supervised-score-map.R:163`); when `NULL` it falls back to
`config$output_routes$features_geometry_gpkg`
(`R/supervised-score-map.R:288-289`). The scoring universe (`scoring_features`) and
`oof_summary` are likewise explicit. So option **(ii)** applies: the runner passes
the SHARED absolute path explicitly. (Option (i) — overriding `cfg$output_routes`
at build time — is NOT supported: `output_routes` is derived from
`output_dir/target_year/run_name/scenario` and not overridable. Option (iii), a
package change, was NOT required.)

The runners now:

1. **Define the SHARED upstream products as explicit absolute paths** in
   `00_common.R` (`SHARED_FEATURES_GPKG`, `SHARED_FOLDS_GPKG`,
   `SHARED_POOLS_GPKG`, the RDS handoffs, and a `shared_inputs_map()`), derived
   from the SHARED `out_base` — NOT reconstructed from each profile's route.
2. **Pass `labelled_features = SHARED_FEATURES_GPKG` explicitly** to
   `score_supervised_burned_map()` in BOTH `20_full.R` and `30_no_hotspot.R`, so
   the score step reads the SHARED features purely from cfg-level routing (no
   convention fallback, no junction).
3. **Fail-fast validate + SHA256 every shared input BEFORE any heavy compute**
   (`validate_shared_inputs()` in `00_common.R`, called at the top of each profile
   after loading SHARED). It aborts with a clear error if any path is unresolvable
   or unhashable, and writes a per-profile `SHARED_INPUTS_PROVENANCE.csv`.

Profile-specific outputs (models, recipes, scoring, maps, EFFIS) still go to each
profile's own folder.

## 4. Junctions removed
All **8** directory junctions were removed (only the reparse points; the SHARED
targets were left intact — verified `features_geometry.gpkg` size 47,906,816 bytes
before and after). The canonical path no longer depends on them.

## 5. Dry-verify (config-level, no heavy compute)
`98_dry_verify_routes.R` (reusable) builds both profiles' deployed cfg and prints
Test-Path results. WITH junctions removed:

- `cfg$output_routes$features_geometry_gpkg` (the OLD convention path) -> **FALSE**
  for both profiles (the original failure path is gone);
- the EFFECTIVE explicit `labelled_features = SHARED_FEATURES_GPKG` -> **TRUE** for
  both profiles (the SHARED file the score step will actually read);
- `SHARED_FOLDS_GPKG` / `SHARED_POOLS_GPKG` -> **TRUE**;
- no upstream subdirs remain under either profile's ENGINE_ROUTES.

VERDICT: both profiles' scoring resolves the SHARED features purely from explicit
cfg routing, without junctions. No training / scoring / feature extraction was run.
Console captured in `dry_verify_routes_console.txt`.

## 6. Files changed (runner working copies; versioned later in G2)
- `00_common.R` — SHARED_* explicit route constants + `shared_inputs_map()` +
  `validate_shared_inputs()` (existence + SHA256 fail-fast).
- `20_full.R` — fail-fast shared-input validation after SHARED load; explicit
  `labelled_features = SHARED_FEATURES_GPKG` at the score call.
- `30_no_hotspot.R` — same two edits.
- `98_dry_verify_routes.R` — new reusable config-level dry-verify.
- `ROUTE_FIX_NOTE.md` — this note.

**No package commit** (script-only fix). To re-run the cheap tail, rerun the
profiles' SCORE+MAP+EFFIS; SHARED (the ~52-min upstream) is reused untouched.
