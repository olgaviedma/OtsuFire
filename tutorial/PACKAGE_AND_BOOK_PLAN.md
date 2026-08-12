# Package and Quarto Book Plan

> **HISTORICAL PLANNING DOCUMENT (superseded).** This file captured an early
> packaging/book plan and does **not** describe the current package. Several
> items here are obsolete: the proposed API names (`run_supervised_pipeline()`,
> `run_master_pipeline()`, `run_loyo_pipeline()`, `build_burned_like_registry()`)
> are not the shipped API; the supervised burned-like registry and the
> `registry.R` module were removed (2026-06-05); ecoregions are not a supervised
> feature; and the supervised stage is a fully package-internal runtime, not a
> bridge to standalone scripts. For the authoritative API see the package
> `NAMESPACE`, the function reference, the METHODS_BOOK chapters, and `NEWS.md`.
> This file is retained only for historical traceability.

## Recommendation

Use the existing `OtsuFire` package as the starting scaffold, but treat the current
deterministic + supervised pipeline as a new cleaned layer to integrate on top of it.

This avoids rebuilding package metadata, tests, and documentation structure from zero,
while keeping the updated workflow separate from legacy code during refactoring.

## Immediate Goal

Freeze a clean, package-ready API around the current active workflow:

- deterministic yearly mapping
- supervised yearly mapping
- consistency checks
- LOYO training/evaluation
- historical application using burned-like registry

## Proposed Public API

High-level exported functions:

- `run_deterministic_pipeline()`
- `run_supervised_pipeline()`
- `run_master_pipeline()`
- `run_consistency_check()`
- `run_loyo_pipeline()` (next phase)
- `build_burned_like_registry()` (next phase)

Internal helpers should stay unexported.

## Proposed Package Modules

Inside `R/`:

- `paths-config.R`
  - path resolution
  - project config objects
- `deterministic-run.R`
  - high-level deterministic wrapper
- `deterministic-filters.R`
  - canonical filters and integration
- `deterministic-rbr.R`
  - RBR keep-pool logic and scoring
- `supervised-run.R`
  - high-level supervised wrapper
- `supervised-unburned.R`
  - unburned strategies
- `supervised-features.R`
  - feature extraction
- `supervised-model.R`
  - OOF, final model, export
- `registry.R`
  - burned-like registry helpers
- `consistency-checks.R`
  - deterministic/supervised consistency checks
- `utils-sf.R`
  - geometry safety helpers
- `utils-raster.R`
  - raster alignment helpers

## What Must Be Cleaned Before Packaging

- no script should auto-run on `source()`
- no reliance on global variables like `DO_*`, `years_to_run`, `scenario_name`
- no hardcoded year defaults inside package functions
- no absolute paths inside reusable functions
- all functions should take explicit arguments or a config object
- outputs should be controlled by function arguments, not script-side globals

## Workflows Outside the Package

Keep user-facing orchestration scripts outside exported package code:

- yearly deterministic run
- yearly supervised run
- multiyear batch run
- LOYO run

These can live in:

- `inst/workflows/`
or
- a separate `workflows/` folder next to the package

## Quarto Book Structure

Recommended book chapters:

1. `index.qmd`
   - project overview
   - deterministic + supervised concept

2. `01-setup.qmd`
   - required inputs
   - folder structure
   - package installation

3. `02-one-year-oof.qmd`
   - run one year
   - inspect OOF
   - inspect `p_burned`
   - inspect final outputs

4. `03-loyo.qmd`
   - LOYO design
   - threshold selection
   - multiyear evaluation

5. `04-historical-application.qmd`
   - burned-like registry
   - applying recent-year model logic to older years

6. `05-troubleshooting.qmd`
   - common failures
   - missing layers
   - raster/grid mismatches

## Recommended Order of Work

1. confirm `2025 balanced` outputs
2. confirm `2022 balanced` outputs
3. freeze current active logic
4. inventory public vs internal functions
5. move active functions into package structure
6. make workflow scripts call package functions
7. start Quarto book

## Best Starting Point

Start from `OtsuFire`, but do not merge everything immediately.

First create a clean branch of the package structure that contains:

- current active deterministic code
- current active supervised code
- current active registry/consistency code

Then, once stable, merge or port pieces from the older package as needed.
