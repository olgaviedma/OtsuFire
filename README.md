<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/OtsuFire_logo.png" width="400"/>
</p>

<h1 align="center">OtsuFire: Otsu-Based Burned-Area Mapping</h1>

**Authors:** Natalia Quintero, Olga Viedma, Hammadi Achour, Jose Manuel Moreno

Reproducible multi-year burned-area mapping from change-index rasters. OtsuFire
takes annual RBR/dNBR composites and produces per-patch burned-area maps through
four stages: a **mosaic** stage, a **deterministic** Otsu/grow segmentation
workflow, a **supervised** one-year probabilistic workflow, and a
**workflow-independent validation** utility.

Imagery is Landsat (through 2016) and Sentinel-2 (from 2017); the package itself
is sensor-agnostic and only requires a change-index raster.

---

## ⚠️ Version 2.0.0 is a complete rewrite

**2.0.0 is a complete new release, not an update of the 0.1.x line.** The
package has been rebuilt from scratch: its public API is 43 exports (36
functions and 7 constants), and runs are driven by configuration objects rather
than long argument lists.

Code written against 0.1.x will not run unchanged. The 0.1.4 release remains
available on CRAN for anyone who needs it.

The supervised stage — training pools, spatial folds, out-of-fold diagnostics
and gradient-boosted scoring — is entirely new in 2.0.0 and has no counterpart
in 0.1.x.

---

## Installation

```r
# Development version (2.0.0):
# install.packages("remotes")
remotes::install_github("olgaviedma/OtsuFire", dependencies = TRUE)

# CRAN currently ships the 0.1.x line, NOT this version:
install.packages("OtsuFire")
```

## System requirements

- **GDAL** (>= 3.0) and **PROJ** (>= 6.0).
- **Python** (>= 3.8) with GDAL bindings — required only by the supervised
  workflow's negative-pool construction, which shells out to
  `gdal_polygonize.py` and `ogr2ogr`.
- **WhiteboxTools** — optional, only used when `grow_engine = "whitebox"` in the
  deterministic stage.

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
per-tile inputs — one GeoTIFF per year, which becomes the `change_index`
argument to both downstream configs.

```r
library(OtsuFire)

mosaic_path <- change_index_mosaic(
  input_dir   = "/path/to/yearly/tiles",
  output_dir  = "/path/to/Composites_90m/Min_Min",
  target_year = 2020L,
  index_type  = "RBR"
)
```

### 2. Deterministic stage

Otsu segmentation per CORINE × ecoregion unit, growth of seed-supported
components, multi-stage filtering, and an `internal_decisions.gpkg` output.

```r
det_cfg <- build_burned_mapping_config(
  change_index   = mosaic_path,
  vegetation_map = "/path/to/Corine_Masks/CLC_2018_peninsula.tif",
  burnable_mask  = "/path/to/Corine_Masks/burnable_mask_binary.tif",
  target_year    = 2020L,
  output_dir     = "/path/to/Results",
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
det_run$internal_decisions_gpkg
```

There are no closed scenario presets: ecological controls are passed explicitly
through `detect_params`, `refine_params` and `scoring_params`. Use the
`*_by_vegetation` keys when seed or growth behaviour must differ by land-cover
class — see `?build_burned_mapping_config`.

For step-by-step access the stage decomposes into `detect_burned_patches()` and
`score_burned_patches()`.

### 3. Supervised one-year stage

Trains a gradient-boosted model on the deterministic decisions plus a
burnable-only unburned negative pool, runs out-of-fold diagnostics, and emits a
`final_map.gpkg` whose `p_burned` field is the per-patch burned model score.

```r
sup_cfg <- build_supervised_burned_config(
  run_label          = "balanced",
  internal_decisions = det_run$internal_decisions_gpkg,
  change_index       = mosaic_path,
  hotspots           = "/path/to/Hotspots/hotspots_iberia_2020.geojson",
  target_year        = 2020L,
  output_dir         = "/path/to/Results",
  options = list(
    data_base      = "/path/to/Data",
    composite_base = "/path/to/Imagery/Composites_90m",
    result_name    = "Min_Min"
  )
)

sup_run <- run_oneyear_supervised_pipeline(sup_cfg, run_consistency = TRUE)
```

Two points that matter methodologically:

- **Negatives come from exactly two sources** — random burnable background and
  Otsu current-year residuals — each with its own cap relative to the burned
  count. Training eligibility is defined by *explicit class, never by negation*:
  review/keep/`NA`/unknown rows can never silently become negatives, and an
  unresolvable row is a hard error.
- **`p_burned` is a model score, not a calibrated probability.** Do not treat it
  as one without an explicit calibration step.

The sub-stages are individually available: `build_supervised_training_pools()`,
`make_spatial_folds()`, `extract_supervised_features()`, `run_oof_diagnostics()`,
`train_final_burned_model()` and `score_supervised_burned_map()`.

### 4. Validation

There are **two distinct validations, against two different references, and they
are not comparable.**

**Internal (OOF).** Out-of-fold diagnostics measure how well the model
reproduces the *deterministic decision* labels using spatial block folds. OOF
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
  ref_shapefile   = "/path/to/Validation_fires_burneable_verano_FINAL/Effis_CA_2020_maskKeep_summer.gpkg",
  mask_shapefile  = "/path/to/Mask_StudyArea_FINAL/mask_3035_2020_final.shp",
  burnable_raster = "/path/to/burnable_mask_binary.tif",
  year_target     = 2020L,
  validation_dir  = "/path/to/VALIDATION",
  observability_raster = "/path/to/MinMin_2020_mosaic_res90m.tif",  # "doy" band
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

`check_supervised_consistency()` compares the deterministic and supervised
outputs of the same run and writes per-issue diagnostic artefacts.

---

## Outputs

A supervised one-year run writes under
`<output_dir>/<target_year>/<run_name>/SUPERVISED/<scenario>/`:

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

The deterministic stage writes to the parallel `DETERMINISTIC/<run_name>/` tree.

---

## Function reference

**Mosaic and rasters** — `change_index_mosaic()`, `mosaic_from_tiles()`,
`clean_raster_file()`

**Deterministic stage** — `build_burned_mapping_config()`,
`run_deterministic_pipeline()`, `detect_burned_patches()`,
`score_burned_patches()`, `apply_chain_cleaning()`

**Supervised stage** — `build_supervised_burned_config()`,
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

**Omission / commission diagnostics** — `of_classify()`,
`of_classify_commission()`, `of_decompose()`, `of_burnable()`, `of_rbr_stats()`,
`of_g3()`, plus the `OF_*` taxonomy constants

## Documentation

```r
vignette("OtsuFire", package = "OtsuFire")           # API tour, toy data
vignette("workflow-overview", package = "OtsuFire")  # four-layer workflow
```

The workflow vignette is not evaluated: the deterministic and supervised stages
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

Quintero, N.; Viedma, O.; Achour, H.; and Moreno, J. M. (2026). *OtsuFire:
Otsu-Based Burned-Area Mapping with Rule-Based Unsupervised Filters and
Supervised Classification.* R package version 2.0.0.
<https://github.com/olgaviedma/OtsuFire>

The 0.1.x line remains available on CRAN:
<https://CRAN.R-project.org/package=OtsuFire>

## Disclaimer

**OtsuFire comes with no guarantee, expressed or implied, and the authors hold
no responsibility for its use or the reliability of its outputs.**
