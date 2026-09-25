<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/OtsuFire_logo.png" width="400"/>
</p>

<h1 align="center">OtsuFire: a self-labelling framework for wall-to-wall burned-area mapping across sensors and decades</h1>

**Authors:** Natalia Quintero, Olga Viedma, Hammadi Achour, Jose Manuel Moreno

Reproducible multi-year burned-area mapping from change-index rasters. OtsuFire
takes annual RBR/dNBR composites and produces per-patch burned-area maps through
four stages: a **mosaic** stage, an **Otsu-guided segmentation** stage, a
one-year **probabilistic refinement** stage, and a **workflow-independent
validation** utility.

The stage names follow the accompanying paper. In the code, the Otsu-guided
segmentation stage is run with `run_deterministic_pipeline()` (configured with
`build_burned_mapping_config()`) and writes to `DETERMINISTIC/`; the
probabilistic refinement stage, which trains a supervised XGBoost model, is run
with `run_oneyear_supervised_pipeline()` (configured with
`build_supervised_burned_config()`) and writes to `SUPERVISED/`.

The package itself is sensor-agnostic and only requires a change-index raster.
The Landsat (through 2016) and Sentinel-2 (from 2017) imagery refers to the
reference workflow the package was developed and validated with (the Iberian
Peninsula), not to a requirement of the package.

---

## ⚠️ Version 2.0.0 is a complete new release

**2.0.0 is a complete new release, not an update of the 0.1.x line.** The
package has been rebuilt from scratch: its public API is 30 exported functions,
and runs are driven by configuration objects rather than long argument lists.

Code written against 0.1.x will not run unchanged. 0.1.4 is currently the
release on CRAN, for anyone who needs it.

The probabilistic refinement stage — training pools, spatial folds, out-of-fold diagnostics
and gradient-boosted scoring — is entirely new in 2.0.0 and has no counterpart
in 0.1.x.

---

## Installation

```r
# Development version (2.0.0):
# install.packages("remotes")
remotes::install_github("olgaviedma/OtsuFire", dependencies = TRUE)

# CRAN currently ships 0.1.4, NOT this version:
install.packages("OtsuFire")
```

## System requirements

- **GDAL** (>= 3.0) and **PROJ** (>= 6.0).
- **Python** (>= 3.8) with GDAL bindings — required only by the supervised
  workflow's negative-pool construction, which shells out to
  `gdal_polygonize.py` and `ogr2ogr`.
- **WhiteboxTools**, through the R package `whitebox` — required by the
  Otsu-guided segmentation stage, which uses it to label seed-supported candidate
  components. Its location can be given with `options$whitebox_exe` in
  `build_burned_mapping_config()`.

Install GDAL into a Python environment with conda and point OtsuFire at the
executables:

```r
python_exe            <- "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"
```

---

## Workflow

### 1. Mosaic stage

`change_index_mosaic()` builds an annual change-index raster (RBR or dNBR) from
raster tiles — one GeoTIFF per year, which becomes the `change_index` input of
both downstream configurations.

```r
library(OtsuFire)

mosaic <- change_index_mosaic(
  folder_path    = "/path/to/yearly/tiles",
  mask           = "/path/to/study_area.gpkg",
  year           = 2020L,
  raster_pattern = "RBR_2020_*.tif",
  output_dir     = "/path/to/mosaics"
)
mosaic_path <- mosaic$mosaic_path
```

### 2. Otsu-guided segmentation stage

Otsu thresholding per land-cover and ecoregion unit, growth of seed-supported
candidate patches, rule-based filtering, and a decision layer
(`internal_decisions.gpkg`) that classifies each patch as `keep`, `review` or
`drop`.

```r
det_cfg <- build_burned_mapping_config(
  change_index   = mosaic_path,
  vegetation_map = "/path/to/vegetation_classes.tif",
  burnable_mask  = "/path/to/burnable_mask.tif",
  target_year    = 2020L,
  output_dir     = "/path/to/results",
  run_name       = "rbr_2020_custom",
  detect_params = list(
    seed_threshold           = 310,
    growth_delta             = 90,
    minimum_growth_threshold = 240,
    minimum_seed_pixels      = 30L
  ),
  refine_params  = list(aoi_buffer_m = 5000, minimum_detected_area_m2 = 10),
  scoring_params = list(minimum_remaining_area_m2 = 20000, reference_buffer_m = 90)
)

det_run <- run_deterministic_pipeline(det_cfg)
det_run$result_paths$internal_decisions
```

Detection, refinement and scoring settings are passed through
`detect_params`, `refine_params` and `scoring_params`. Use the
`*_by_vegetation` keys when seed or growth behaviour must differ by land-cover
class — see `?build_burned_mapping_config`.

For step-by-step access the stage decomposes into `detect_burned_patches()` and
`score_burned_patches()`.

### 3. One-year probabilistic refinement stage

Trains a gradient-boosted model on labels derived from the deterministic
decisions, runs out-of-fold diagnostics, and writes a `final_map.gpkg` whose
`p_burned` field is the per-patch burned model score.

```r
sup_cfg <- build_supervised_burned_config(
  run_label          = "balanced",
  internal_decisions = det_run$result_paths$internal_decisions,
  change_index       = mosaic_path,
  hotspots           = "/path/to/hotspots_2020.gpkg",
  target_year        = 2020L,
  output_dir         = "/path/to/results",
  options = list(
    data_base      = "/path/to/data",
    composite_base = "/path/to/composites"
  )
)

sup_run <- run_oneyear_supervised_pipeline(sup_cfg, run_consistency = TRUE)
```

Two points that matter methodologically:

- **Three differentiated unburned pools.** Burned labels are the keep patches.
  Unburned labels come from a *burnable background* pool (weak change-index
  response, outside a 500-m exclusion zone around all candidate patches), a
  *moderately burned-like* pool (drop patches with `S_PATCH_PA <= 0.15`) and an
  optional *strongly burned-like* pool (drop patches with a stronger
  burned-like response; off by default, and optionally filtered by visual
  validation with `apply_visual_validation()`). The first two are capped
  relative to the burned count; all eligible strongly burned-like patches are
  kept and balanced by weight. Review patches are withheld from training and
  scored afterwards. Training uses explicit labels only: review, keep, `NA` or
  unknown rows never become negatives.
- **`p_burned` is a model score, not a calibrated probability.** Do not treat it
  as one without an explicit calibration step.

The sub-stages are individually available: `build_supervised_training_pools()`,
`make_spatial_folds()`, `extract_supervised_features()`, `run_oof_diagnostics()`,
`train_final_burned_model()` and `score_supervised_burned_map()`.

### 4. Validation

There are **two distinct validations, against two different references, and they
are not comparable.**

**Internal (OOF).** Out-of-fold diagnostics measure how well the model
reproduces the *Otsu-guided segmentation decision* labels using spatial block folds. OOF
AUC / PR-AUC are **not ground truth** — the labels are decisions, not field
observations.

**External.** `validate_fire_maps()` is workflow-independent: threshold
`p_burned` into a burned-area layer and compare it against an independent
reference such as EFFIS. It reports omission, commission, F1 and IoU.

```r
fm <- sf::st_read(sup_run$final_map_gpkg, layer = "final_map_full", quiet = TRUE)
burned <- fm[is.finite(fm$p_burned_model) & fm$p_burned_model >= 0.5, ]
sf::st_write(burned, "/path/to/thresholded_burned.gpkg", quiet = TRUE)

val <- validate_fire_maps(
  input_shapefile = "/path/to/thresholded_burned.gpkg",
  ref_shapefile   = "/path/to/reference_fires_2020.gpkg",
  mask_shapefile  = "/path/to/validation_mask_2020.gpkg",
  burnable_raster = "/path/to/burnable_mask.tif",
  year_target     = 2020L,
  validation_dir  = "/path/to/validation",
  observability_raster = "/path/to/observation_doy_2020.tif",
  observability_mode   = "wholefire_fraction",
  ref_obs_doy_col      = "obs_required_doy"
)
```

When an `observability_raster` (day-of-year) is supplied, observability is
assessed **temporally at the reference-fire level**: with
`observability_mode = "wholefire_fraction"` a reference fire is evaluated only
if at least 75 % of its burnable pixels were observed on or after its required
date, and a fire that fails is removed from both the reference and the
evaluation domain. This is *not* a pixel-level cloud mask — fires are kept or
removed whole.

`check_supervised_consistency()` compares the Otsu-guided segmentation and probabilistic refinement
outputs of the same run and writes per-issue diagnostic artefacts.

---

## Outputs

A one-year probabilistic refinement run writes under
`<output_dir>/<target_year>/<run_name>/SUPERVISED/<run_label>/`:

| Directory | Contents |
|---|---|
| `01_POOLS/` | burned + unburned training pools |
| `02_FOLDS/` | spatial folds for OOF |
| `03_FEATURES/` | full feature matrix |
| `04_MATRIX/` | design-matrix snapshot |
| `05_OOF/` | out-of-fold diagnostics |
| `07_FINAL_MODEL_V2/` | trained model |
| `09_FINAL_MAP/` | final per-patch `p_burned` map |
| `11_CONSISTENCY_CHECKS/` | consistency artefacts, when enabled |

The Otsu-guided segmentation stage writes to the parallel `DETERMINISTIC/<run_name>/` tree.

---

## Function reference

**Mosaic and rasters** — `change_index_mosaic()`, `mosaic_from_tiles()`,
`clean_raster_file()`

**Otsu-guided segmentation stage** — `build_burned_mapping_config()`,
`run_deterministic_pipeline()`, `detect_burned_patches()`,
`score_burned_patches()`, `apply_chain_cleaning()` (also exported under the
short alias `apply_chain()`)

**Probabilistic refinement stage** — `build_supervised_burned_config()`,
`run_oneyear_supervised_pipeline()`, `build_supervised_training_pools()`,
`make_spatial_folds()`, `extract_supervised_features()`, `run_oof_diagnostics()`,
`train_final_burned_model()`, `score_supervised_burned_map()`,
`write_thresholded_burned()`, `validate_supervised_execution()`

**Training-pool curation** — `diagnose_training_pools()`,
`export_visual_validation_pools()`, `apply_visual_validation()`,
`collect_keep_reference_samples()`, `build_keep_pool_from_samples()`,
`assemble_era_training_pool()`, `promote_artifact_hard_negatives()`,
`subsample_random_background()`, `add_pool_shape_features()`

**Validation and consistency** — `validate_fire_maps()`,
`check_supervised_consistency()`

## Documentation

```r
vignette("OtsuFire", package = "OtsuFire")           # API tour, toy data
vignette("workflow-overview", package = "OtsuFire")  # four-layer workflow
```

The workflow vignette is not evaluated: the Otsu-guided segmentation and probabilistic refinement stages
operate on full-scale annual mosaics (Iberian Peninsula, 90 m) that are too
large to ship with the package. The `OtsuFire` vignette covers the API surface
on synthetic data with no external downloads.

---

## Acknowledgements

We gratefully acknowledge funding from project INFORICAM (PID2020-119402RB-I00),
funded by the Spanish MCIN/AEI/10.13039/501100011033 and by the "European Union
NextGenerationEU/PRTR". Carlos Silva was supported by NASA's Carbon Monitoring
System funding (CMS, grant 22-CMS22-0015).

## Reporting issues

Please report any issue regarding the OtsuFire package to Dr. Olga Viedma
([olga.viedma@uclm.es](mailto:olga.viedma@uclm.es)).

## Citing OtsuFire

Quintero, N.; Viedma, O.; Achour, H.; and Moreno, J. M. (2026). *OtsuFire: a
self-labelling framework for wall-to-wall burned-area mapping across sensors and
decades.* R package version 2.0.0.
<https://github.com/olgaviedma/OtsuFire>

OtsuFire 0.1.4 is currently the release on CRAN:
<https://CRAN.R-project.org/package=OtsuFire>

## Disclaimer

**OtsuFire comes with no guarantee, expressed or implied, and the authors hold
no responsibility for its use or the reliability of its outputs.**
