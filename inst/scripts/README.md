# OtsuFire — canonical supervised usage scripts (`inst/scripts/`)

These are the **CANONICAL, VERSIONED templates** for the official supervised
burned-area mapping workflow. They live **inside the package repository** so the
official usage recipes are no longer only-outside-version-control.

Target: **OtsuFire >= 0.5.0**.

## What "canonical" means here

`inst/scripts/` is the **single source of truth**. Any other copy of these
scripts — in particular the external working copies under
`C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/00_USAGE/02_SUPERVISED_USAGE/`
— is a **WORKING COPY**, not the source. The sync direction is one-way:

```
inst/scripts/ (canonical)  ->  00_USAGE/02_SUPERVISED_USAGE/ (working copies)
```

Each external working-copy script carries a banner at its top pointing back here.
`MANIFEST.csv` records a content hash of every canonical script — the
`git hash-object` blob SHA-1 (stable in-repo identity, survives CRLF
normalization) plus the LF-content md5 — so that divergence between canonical and
working copies can be detected. Prefer the blob SHA-1 for divergence checks:
`git hash-object inst/scripts/<file>`.

## Design rules every script follows

1. **Public API only** + `build_supervised_burned_config()` as the
   **single source of truth** for all methodological params (the 4 negative-pool
   caps, `training_protocol`, `oof_sampling`, seeds). No deprecated function-level
   parameter shims are used (those emit `otsufire_deprecated_param` warnings).
2. **Portable.** A clearly-marked `CONFIG` block at the top holds `<PATH_TO_...>`
   placeholders and a `paths <- list(...)` / `tool_paths <- list(...)` block with
   obvious `# TODO` markers. No project absolute data paths are hardcoded — the
   user edits the placeholders.
3. **Safe to source for a syntax check.** The real pipeline lives in `main()` and
   is guarded behind `RUN <- FALSE`. Sourcing the file with `RUN = FALSE` does
   **not** load the package or run any compute; it just prints a reminder. Set
   `RUN <- TRUE` after editing the CONFIG block to actually run.
4. **`p_burned` is a MODEL SCORE, not a calibrated probability** (per the package
   docs). Use it for ranking / thresholding; pick the operating threshold from
   external (EFFIS) validation, not from the score scale.
5. `validate_supervised_execution(cfg, strict = TRUE)` is the **pre-run** check
   (structural / schema / inputs). It is distinct from `validate_fire_maps()`,
   which is the **post-run** accuracy validation against an external reference.

## The scripts

| # | File | Purpose (one line) |
|---|------|--------------------|
| 01 | `01_supervised_oneyear_legacy.R` | Corrected common-recipe one-year pipeline (legacy training protocol), full 6-stage modular chain. |
| 02 | `02_supervised_phaseB_config.R` | Phase B config: `cap_spectral = 2.0` via the builder + `nested_refit`, with the abort-on-cap-mismatch guard and the inherited runtime feature-schema parity guard. |
| 03 | `03_supervised_protocol_legacy.R` | Minimal `training_protocol = "legacy"` example. |
| 04 | `04_supervised_protocol_nested_refit.R` | Minimal `training_protocol = "nested_refit"` example (controlled pair with 03). |
| 05 | `05_supervised_effis_validation.R` | EFFIS validation via `validate_fire_maps()` on a thresholded map (distinct from `validate_supervised_execution()`). |
| 06 | `06_supervised_no_hotspots_reduced_historical.R` | No-hotspots / reduced-historical (pre-MODIS) run: `hotspots = NULL`, `use_hotspots = FALSE`. |

## How to run one

1. Open the script and edit the `CONFIG` block: replace every `<PATH_TO_...>`
   placeholder and set the run-level settings (`target_year`, `scenario`, ...).
2. Set `RUN <- TRUE`.
3. `source("inst/scripts/0X_....R")` (or `Rscript inst/scripts/0X_....R` after
   setting `RUN <- TRUE`).

The pipeline loads the development tree with `pkgload::load_all()` (NOT
`library(OtsuFire)`, which may pick up an older installed copy that lacks the
modular API), then runs the 6 exported stages, reading every methodological knob
from the cfg.

## Known placeholders the user MUST fill

Every `<PATH_TO_...>` in the CONFIG block is required:

- `pkg_root` — the OtsuFire_v02_rebuild development tree.
- `internal_decisions`, `change_index` — the two **required** run inputs.
- `hotspots` — required post-2000; `NULL` for pre-MODIS years (script 06).
- `reference_burned`, `burnable_raster`, `validation_mask`, `strata_*` — only for
  the EFFIS validation script (05); `strata_*` may be `NULL` to skip strata.
- `data_base`, `composite_base`, `output_dir` — roots forwarded to the legacy
  all_sources Otsu builder; `output_dir` is the **root** (the package appends
  `<year>/<run_name>/SUPERVISED/<scenario>`).
- `tool_paths` (`python_exe`, `gdal_polygonize_script`, `gdalwarp_path`,
  `ogr2ogr_exe`) — the external GDAL/Python tools the all_sources legacy Otsu
  builder requires.
