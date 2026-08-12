# OtsuFire 2.0.0 — step-by-step tutorial

A short, runnable tour of the whole workflow, one stage per file. Written
against the **2.0.0 API** (43 exported functions, config + run pattern).

All paths in these files are placeholders (`/path/to/...`). Edit
`00_SETUP.R` once, and the other four files pick everything up from there.

## The files

| File | Stage | What it produces |
|---|---|---|
| `00_SETUP.R` | — | Loads the package, declares every path. Runs no analysis. |
| `01_MOSAIC.R` | Mosaic | One change-index raster for the year. |
| `02_DETERMINISTIC.R` | Deterministic | Candidate patches labelled keep / review / drop. |
| `03_SUPERVISED.R` | Supervised | A trained model and `p_burned` for every patch. |
| `04_VALIDATION.R` | Validation | A burned-area map and metrics against EFFIS. |

Run them in order. Inside each file, run block by block rather than sourcing
the whole thing: the point is to look at what comes out of each step.

## How the stages connect

```
tiles ──01──> change_index.tif
                   │
                   ├──02──> internal_decisions.gpkg   (keep / review / drop)
                   │              │
                   └──03──────────┴──> final_map.gpkg (p_burned per patch)
                                              │
                                        04──> burned map + metrics
```

Stage 2 produces the labels stage 3 learns from. Stage 4 is the only stage
that compares anything against an independent reference.

## Four things that are easy to get wrong

**`p_burned` is a score, not a probability.** 0.8 does not mean "80 % chance
of burned". Do not read it as calibrated, and do not compare raw scores
across years without calibrating first.

**Pick the threshold per year.** The best cut varies year to year. `04` picks
it from the data by maximising pixel F1 against the reference. Using 0.5
because it is round costs real accuracy.

**Out-of-fold metrics are not ground truth.** They measure how well the model
reproduces stage 2's *rules*. Only stage 4's external metrics answer "does
this match real fires?". The two are not comparable — a high OOF AUC with a
mediocre external F1 is perfectly consistent.

**Folds are spatial for a reason.** Fires cluster in space. Random folds put a
patch in training and its neighbour in test, and the metrics come out far too
optimistic. Block folds prevent that.

## Requirements

- R, plus the package and its dependencies (`sf`, `terra`, `xgboost`, …).
- **xgboost 1.7.x.** The package also supports 2.x/3.x, but a model saved by
  one generation cannot be read by the other. Do not upgrade mid-project.
- Input data: yearly severity composites, CORINE land cover and a burnable
  mask, a study-area polygon, ecoregions, terrain, and an external reference
  (EFFIS) for validation. Hotspots are optional and only exist from 2000.

## Where to go deeper

- `?build_burned_mapping_config` — every deterministic parameter, including
  the per-vegetation overrides.
- `?build_supervised_burned_config` — the negative pool, caps, seeds and
  training control.
- `?validate_fire_maps` — every validation option, including temporal
  observability.
- `vignette("workflow-overview", package = "OtsuFire")` — the same workflow
  at API level, with the design rationale.
