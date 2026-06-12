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
   caps, seeds). OtsuFire always uses one training procedure
   (inner-early-stopping selection + full-data refit) — there is no
   training-protocol choice, and OOF always uses the same capped negative-sampling
   policy as the FINAL model (applied independently within each training fold).
   No deprecated function-level parameter shims are
   used (those emit `otsufire_deprecated_param` warnings).
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

## Training eligibility + the single shared capping helper

The supervised model uses **two explicit negative sources: random background and
Otsu residuals.** There are no other negative buckets (the contextual + spectral
buckets were removed in GATE 6). The two per-bucket caps live in the typed PUBLIC
`negative_pool_params$caps = c(random = 1.0, otsu = 1.0)` (the single source of
truth), mirrored into `cfg$train_control$caps`.

Training eligibility is defined by **explicit class, never by negation**:

- `class == "burned"` → **positive** (always kept);
- `class == "unburned"` **and** the row resolves to exactly one of the two
  valid negative buckets (`random`, `otsu`) →
  **negative**, eligible for capping;
- otsu review / keep rows (the otsu exclude-list) → **excluded + logged**, never
  trained;
- anything else (`class` `NA` / unknown, or an unburned row with no valid
  bucket) → **hard error** (strict), surfaced with id / class / source /
  neg_type / origin-stage.

Review / keep / `NA` / unknown rows can **never** silently become negatives. The
old "negatives = everything that is not burned" rule (the OOF `!is_burned`
predicate with its `sel_other` / `other_kept` 5th category) has been removed.

The OOF folds and the FINAL model share **one** internal eligibility resolver
and **one** capping helper (`cap = ceiling(n_burned * ratio)`, take-all when
available ≤ cap, all positives kept, negatives sampled without replacement in a
deterministic bucket order). The only differences between OOF and FINAL are the
retained outer-test fold (OOF) and the RNG seed. This is behaviour-preserving
for the default 2017 balanced run.

## The scripts

| # | File | Purpose (one line) |
|---|------|--------------------|
| 01 | `01_supervised_oneyear_otsu_negative.R` | Common-recipe one-year pipeline, full 6-stage modular chain. |
| 02 | `02_supervised_phaseB_config.R` | Phase B config: two-source caps (`random` + `otsu`) via the typed `negative_pool_params` builder block, with the abort-on-cap-mismatch guard and the inherited runtime feature-schema parity guard. |
| 03 | `05_supervised_effis_validation.R` | EFFIS validation via `validate_fire_maps()` on a thresholded map (distinct from `validate_supervised_execution()`). |
| 04 | `06_supervised_no_hotspots_reduced_historical.R` | No-hotspots / reduced-historical (pre-MODIS) run: `hotspots = NULL`, `use_hotspots = FALSE`. |

> The single OtsuFire training procedure (inner-validation early stopping to
> select the number of boosting rounds, then a full-data refit) is fixed; there
> is no `training_protocol` argument. The former protocol-pair example scripts
> (`03_supervised_protocol_legacy.R` / `04_supervised_protocol_nested_refit.R`)
> were removed because they demonstrated a choice that no longer exists.

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
- `data_base`, `composite_base`, `output_dir` — roots forwarded to the
  all_sources Otsu residual negative builder; `output_dir` is the **root** (the package appends
  `<year>/<run_name>/SUPERVISED/<scenario>`).
- `tool_paths` (`python_exe`, `gdal_polygonize_script`, `gdalwarp_path`,
  `ogr2ogr_exe`) — the external GDAL/Python tools the all_sources Otsu residual
  negative builder requires.
