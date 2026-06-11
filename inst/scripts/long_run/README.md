# OtsuFire — canonical LONG-RUN paired-experiment runners (`inst/scripts/long_run/`)

These are the **CANONICAL, VERSIONED** runner templates for the paired
**FULL vs NO-HOTSPOT** long-run supervised experiment. They live inside the
package repository so the official long-run recipe is under version control.

The matching **working copies** live under a run's `_ORCHESTRATION/` folder
(e.g. `00_LONG_RUNS/LONG_<YEAR>_BALANCED_<STAMP>/_ORCHESTRATION/`). Those copies
carry a banner pointing back here; **this folder is the source of truth**.

## Architecture: SHARED -> FULL -> NO-HOTSPOT

```
                 05_env_snapshot.R   (provenance; SHA256 of every input)
                          |
                 10_shared_upstream.R
                 candidates -> labels -> pools -> spatial folds ->
                 extract_supervised_features (FULL feature set incl hotspots)
                 written ONCE under SHARED/ENGINE_ROUTES, FULL 50-base feature GPKG
                          |
        +-----------------+------------------+
        |                                    |
   20_full.R                          30_no_hotspot.R   (gated on FULL_SUCCESS)
   DEPLOYED OOF (capped) + FINAL      feature_whitelist_override = 38 non-hotspot
   score + map + EFFIS                use_hotspots = FALSE
                                      DEPLOYED OOF (capped) + FINAL
                                      score + map + EFFIS
```

The **SHARED upstream is built once** and **reused by both profiles**. The two
profiles differ **only by explicit configuration**:

- **FULL** = canonical 50-base-feature whitelist, hotspots ON.
- **NO-HOTSPOT** = `feature_whitelist_override` = the 38 non-hotspot base
  features (50 minus the 12 hotspot-lineage features) + `use_hotspots = FALSE`.

Nothing upstream is regenerated for NO-HOTSPOT; it selects a feature subset from
the SAME SHARED feature GPKG.

## Explicit-cfg-routes contract (no junctions) — the G2.3 fix

Because the SHARED upstream is built once under `SHARED/ENGINE_ROUTES` and each
profile builds its OWN cfg with a different `output_dir`, the score step's
**convention** route `cfg$output_routes$features_geometry_gpkg` reconstructs to
the *profile's* `ENGINE_ROUTES` — where no features GPKG exists. The original
long run failed there (`No existe labelled_features_gpkg: ...`).

The **canonical fix** uses NO directory junctions:

1. `00_common.R` exposes the SHARED upstream products as **explicit absolute
   routes** (`SHARED_FEATURES_GPKG`, `SHARED_FOLDS_GPKG`, `SHARED_POOLS_GPKG`,
   the RDS handoffs, and `shared_inputs_map()`).
2. `20_full.R` / `30_no_hotspot.R` pass
   **`labelled_features = SHARED_FEATURES_GPKG`** explicitly to
   `score_supervised_burned_map()`, so the score step never falls back to the
   convention route.
3. `validate_shared_inputs()` **fail-fast validates + SHA256-hashes every shared
   input BEFORE any heavy compute**, writing a per-profile
   `SHARED_INPUTS_PROVENANCE.csv`.

`98_dry_verify_routes.R` proves this at config level (no heavy compute): the
convention path is ABSENT, the explicit `labelled_features` resolves to the real
SHARED file, and no junctions remain. See `ROUTE_FIX_NOTE.md`.

## EFFIS (memory-safe)

`run_long_effis()` in `00_common.R` is the dedicated EFFIS job: it sets
`terraOptions(todisk = TRUE, modest memfrac, dedicated tempdir)`, thresholds the
final map at the operating threshold, runs `validate_fire_maps()` on the 3035
90 m composite-aligned burnable mask (memory-safe), then cleans up the terra
tempdir. This is the long-run counterpart of the manual script's
`EFFIS_MODE = "STANDARD_90M"`; both call the **same** `validate_fire_maps()` with
the **same** metric definitions.

## Files

| File | Purpose |
|------|---------|
| `00_common.R` | Sourced by every stage: CONFIG block (placeholders), version pin, paths/caps/seeds, feature lineage, explicit SHARED routes + `validate_shared_inputs()`, cfg builders, `run_long_effis()`. |
| `05_env_snapshot.R` | Provenance snapshot (disk/mem, version pin, tarball + input SHA256) before any heavy step. |
| `10_shared_upstream.R` | Builds the SHARED upstream ONCE (pools/folds/features, FULL feature set). |
| `20_full.R` | FULL profile: DEPLOYED OOF (capped) + FINAL, deploy FINAL, score + map + EFFIS. OOF uses the same capped policy as FINAL. |
| `30_no_hotspot.R` | NO-HOTSPOT profile (38-feature override), gated on FULL success. |
| `98_dry_verify_routes.R` | Config-level dry verify of the explicit-route resolution (no heavy compute, no junctions). |
| `MASTER_RUN.bat` | Unattended orchestrator: env -> SHARED -> FULL -> (gate) -> NO-HOTSPOT. |
| `RESUME_RUN.bat` | Skip the SHARED rebuild; re-run the cheap tail (no junctions). |
| `MANIFEST.md` / `MANIFEST.csv` | Content hash per canonical file. |
| `reports/` | Versioned audit reports (route-fix note, failure diagnosis, run manifest). Heavy artifacts are referenced by absolute path + SHA256, NOT copied. |

## How to run

1. Copy this folder into your run's `_ORCHESTRATION/`, OR edit in place.
2. In `00_common.R`, fill every `<PATH_TO_...>` / `<...>` placeholder in the
   CONFIG block (LONG_ROOT, RLIB, PIN_VERSION, snapshot id/dir/sha, data base,
   tool paths). Set `target_year`/`scenario` to your run.
3. In each stage script and `.bat`, set `<PATH_TO_ORCHESTRATION>` /
   `<PATH_TO_Rscript.exe>` / `<PATH_TO_LONG_RUN_ROOT>` / `<SNAPSHOT_ID>`.
4. (Optional) `Rscript 98_dry_verify_routes.R` once SHARED exists, to confirm
   route resolution.
5. Launch `MASTER_RUN.bat` (or `RESUME_RUN.bat` to re-run only the tail).

## Frozen-code contract

The runners load the **INSTALLED** OtsuFire from an **isolated library** (`RLIB`)
and **pin the version** (`.assert_version_pin()`); they NEVER call
`pkgload::load_all()`. The profiles differ only by explicit config — never by
package edits. On any failure the policy is: STOP, write `FAILURE_REPORT.txt`,
do not patch the package, do not relaunch with changes.
