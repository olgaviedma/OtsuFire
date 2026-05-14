# OtsuFire: Fire Mapping Assessment Toolkit

<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/OtsuFire_logo.png" width="400"/>
</p>

<h1 align="center">OtsuFire: Fire Scars and Severity Mapping Toolkit</h1>

**Authors:** Natalia Quintero, Olga Viedma, Hammadi Achour, Jose Manuel Moreno

Automated tools to map fire scars using Landsat-based RBR/dNBR composites and Otsu thresholding. Includes mosaic preparation, segmentation, refinement, scoring, merging, and validation workflows.

## Quick Workflow Summary
1. Download example data from Zenodo.
2. Mosaic, reproject, mask and resample raster inputs.
3. Optionally reclassify CORINE and create burnable masks.
4. Run Otsu segmentation with class-based seed/growth thresholds.
5. Refine segmentation by AOIs and merge outputs.
6. Score internal and reference burned polygons.
7. Compute polygon metrics and apply geometric filters.
8. (Optional) Merge internal and reference scored polygons.
9. (Optional) Run histogram/tests diagnostics.
10. Prepare datasets and validate fire maps.
11. Compute validation statistics and separability diagnostics.
---

## System Requirements

### Python & GDAL

Several functions in `OtsuFire` call the GDAL Python script `gdal_polygonize.py`. To use these functions, ensure:

- **Python** is installed (e.g., via Anaconda or Miniconda).
- **GDAL** is installed in your Python environment:
  - Via **Conda**: `conda install -c conda-forge gdal`
  - Or via **pip**: `pip install gdal` (less recommended)
- You know the full paths to:
  - Your Python executable (`python.exe`)
  - The `gdal_polygonize.py` script (typically in `Scripts/` on Windows)

 **Example configuration:**

```r
python_exe <- "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

```

###GDAL Installation Guide

# Open Anaconda Prompt or CMD
conda install -c conda-forge gdal

# Check version
gdalinfo --version

# Confirm availability
where gdalinfo

### Optional: Create a Dedicated GDAL Environment

conda create --name gdal_env -c conda-forge gdal
conda activate gdal_env
gdalinfo --version

### Validate GDAL in Python

python

import osgeo.gdal
print(osgeo.gdal.__version__)

### Additional Notes

- All input rasters must be projected and masked consistently.
- Functions support large-area mosaics, tiling, and region-specific thresholding.
- Always ensure gdal_polygonize.py is accessible and executable.

<small> Note: All input files are automatically downloaded and extracted from Zenodo (https://doi.org/10.5281/zenodo.15322380) into a local ZENODO/exdata directory. Paths in this document assume that directory structure.

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("OtsuFire")

# The development version:
#install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/OtsuFire", dependencies = TRUE)


```

## Libraries
```r

# Load the OtsuFire package and dependencies
library(OtsuFire)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(glue)
library(purrr)
library(tools)
library(data.table)
library(histogram)
library(OtsuSeg)
library(stats)
library(whitebox)

```

## Download and Prepare Example Data from Zenodo

This package uses synthetic raster and shapefile data for testing and demonstration purposes. The data is archived on Zenodo and will be downloaded and extracted into a local folder (`inst/extdata/ZENODO`) when running the following code.

The following script:

- Downloads the dataset from Zenodo (only if not already downloaded)
- Validates the download
- Extracts the content
- Places the extracted files in a consistent folder structure

## 0.0. Import data from ZENODO
```r

library(utils)
out <- download_zenodo_data(
  zenodo_url = "https://zenodo.org/records/15322380/files/DATA.zip?download=1",
  output_dir = "ZENODO/exdata",
  overwrite_zip = FALSE,     # NO borra el parcial
  overwrite_data = FALSE,
  use_curl = TRUE
)

```

## 0.1. Mosaic and Resample Raster Data
```r

rasters_dir <- "ZENODO/exdata"
mask_file <-"ZENODO/exdata/iberian_peninsula_proj_final.shp"


out_mosaic <- mask_mosaic_raster(
  folder_path    = "ZENODO/exdata",
  mask_shapefile = "ZENODO/exdata/iberian_peninsula_proj_final.shp",
  apply_mask     = TRUE,
  year           = 2012,
  tif_glob       = "MinMin_*.tif",
  exclude_glob   = c("*_mosaic*.tif", "*_proj.tif", "*_res*m.tif"),
  crs_target      = 3035,     # <- ETRS89_LAEA
  proj_res_m      = NULL,
  out_res_m       = 90,      
  band_names      = c("rbr", "doy"),
  gdalwarp_path   = "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe",
  return_stage   = "resampled",
  keep_intermediates = "none"
)
out_mosaic

```

## 0.2. OPTIONAL: Reclassify and mask CORINE raster
```r

# --- Your reclass matrix ---
my_reclass <- matrix(c(
   22, 1,  # agroforestry
   26, 2,  # grass
   32, 3,  # sparsely vegetated
   33, 3,  # burnt areas (often excluded from segmentation)
   27, 4,  # shrubs (moors)
   28, 4,  # sclerophyllous veg
   29, 4,  # transitional wood
   23, 5,  # broadleaves
   25, 6,  # mixed
   24, 7,  # conifers
   1,  8,  2,  8,  3,  8,  4,  8,  5,  8,  6,  8,  7,  8,  8,  8,  9,  8, 10, 8, 11, 8,  # 1-11 -> 8 (artificial)
   12, 9, 13, 9, 14, 9, 15, 9, 16, 9, 17, 9, 18, 9, 19, 9, 20, 9, 21, 9,  # 12-21 -> 9 (agriculture)
   35,10, 36,10, 37,10, 38,10, 39,10, 40,10, 41,10, 42,10, 43,10, 44,10,  # 35-44 -> 10 (water)
   30,11, 31,11, 34,11  # bare soils etc. (avoid duplicating 32 here)
 ), ncol = 2, byrow = TRUE)

burnable_01_group <- matrix(c(
   1,1, 2,1, 3,1, 4,1, 5,1, 6,1, 7,1, 8,0,
   9,0, 10,0, 11,0
 ), ncol = 2, byrow = TRUE)


cor_groups <- list( agroforestry = c(22), grasslands = c(26), open_areas = c(32, 33), shrublands = c(27,28, 29), broadleaves = c(23), mixed = c(25), conifers = c(24), bare_soils = c(30,31,34), artificial = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), agriculture = c(12, 13, 14, 15, 16, 17, 19, 20, 21), water = c(35, 36, 37, 38, 39,40, 41, 42, 43, 44) )

base_dir <- "ZENODO/exdata"

out <- corine_mask_reclass(
  corine_rasters       = file.path(base_dir, "U2018_CLC2018_V2020_20u1_IBERIAN_ETRS89.tif"),
  peninsula_shapefile  = file.path(base_dir,"iberian_peninsula_proj_final.shp"),
  template_raster_path = file.path(base_dir,"MinMin_2025_mosaic_res90m.tif"),

  corine_value_col   = "GRID_CODE",
  group_rcl_matrix   = my_reclass,
  exclude_burnt33    = FALSE,

  burnable_matrix_01 = burnable_01_group,
  burnable_source    = "group",

  output_raster_dir  = base_dir,
  out_base_names     = "corine18",

  write_group_raster     = TRUE,
  write_burnable_classes = TRUE,
  write_binary_mask      = TRUE,

  # LUT
  corine_grid_table      = file.path(base_dir,"clc_legend_gridcode.csv"),
  lut_main_level   = c("LABEL2"),
  lut_detail_level = c("LABEL3"),
  lut_main_strategy = c("unique_concat"),
  write_group_lut_csv    = TRUE,
  write_group_lut_xlsx   = TRUE,
  lut_out_dir            = file.path(base_dir,"LUT1"),
  lut_file_tag           = "strata8_v1",

  overwrite = TRUE,
  verbose   = TRUE
)

```

## 1. Apply Otsu Thresholding to RBR Raster
```r


base_out <- "ZENODO/exdata/2025/OTSU_MIN"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

base_dir <- "ZENODO/exdata"

raster_path <- file.path(base_dir, "MinMin_2025_mosaic_res90m.tif")
corine_raster_path <- file.path(base_dir,"corine18_burnable_classes.tif")
ecoregion_shapefile_path <- file.path(base_dir,"ecoregiones_olson.shp")
peninsula_shapefile <- file.path(base_dir,"iberian_peninsula_proj_final.shp")
  
 #  22, 1,  # agroforestry
 #  26, 2,  # grass
 #  32, 3,  # sparsely vegetated
  # 33, 3,  # burnt areas (often excluded from segmentation)
 #  27, 4,  # shrubs (moors)
 #  28, 4,  # sclerophyllous veg
 #  29, 4,  # transitional wood
 #  23, 5,  # broadleaves
 #  25, 6,  # mixed
 #  24, 7  # conifers)

  # 1) SEEDS BY CLASS (LB_CLASS)
  #otsu_min_by_class <- list(corine = list( "2" = 300, "3" = 300, "4" = 300,  "6" = 300,"7" = 300,"8" = 300))
  #otsu_min_by_class <- list(corine = list( "1" = 200,  "2" = 100, "3" = 200, "4" = 350, "5" = 350,"6" = 350,"7" = 350))
  otsu_min_by_class <- list(corine = list( "1" = 300,  "2" = 350, "3" = 350, "4" = 350, "5" = 350,"6" = 350,"7" = 350))
  
  # 2) DELTA BY CLASS
  #grow_delta_by_class <- list(corine = list("2" = 100, "3" = 100, "4" = 100,"6" = 100, "7" = 100,"8" = 100))
  #grow_delta_by_class <- list(corine = list( "1" = 100, "2" = 50, "3" = 100, "4" = 50,"5" = 50, "6" = 50,"7" = 50) )
  grow_delta_by_class <- list(corine = list( "1" = 100, "2" = 50, "3" = 100, "4" = 100,"5" = 50, "6" = 50,"7" = 50) )
  
  # 3) GROW_FLOOR_USED BY CLASS)
  #min_grow_threshold_by_class <- list(corine = c( "2" = 200, "3" = 200,"4" = 200, "6" = 200, "7" = 200, "8" = 200 ) )
 # min_grow_threshold_by_class <- list(corine = c( "1" = 100, "2" = 50,"3" = 100, "4" = 300,"5" = 300,"6" = 300, "7" = 300 ))
  min_grow_threshold_by_class <- list(corine = c( "1" = 200, "2" = 300,"3" = 250, "4" = 250,"5" = 300,"6" = 300, "7" = 300 ))
  
  out_growth <- file.path(base_out, "CORINE_ECOREG_seed_30pix")
  dir.create(out_growth, showWarnings = FALSE, recursive = TRUE)
  
  
  res <- process_otsu_rasters_grow(
    raster_path = raster_path,
    output_dir  = out_growth,
    year = 2025,
    
    peninsula_shapefile = peninsula_shapefile,
    corine_raster_path = corine_raster_path,
    reclassify_corine  = FALSE,
    reclass_matrix     = NULL,
    output_corine_raster_dir = NULL,
    corine_classes = NULL,
    
    ecoregion_shapefile_path = ecoregion_shapefile_path,
    ecoregion_field = "EnZ_name",
    segment_by_intersection = TRUE,
    ecoregion_touches = TRUE,
    
    otsu_value_range = c(0, 1500),
    trim_percentiles = list(min = 0.05, max = 0.99),
    
    
    # LB_GLOBAL
    otsu_thresholds = c(300),
    
    # Seeds (LB_CLASS)
    otsu_min_by_class = otsu_min_by_class,
    
    # Growing
    grow_delta = 100,  # default for not specified classes
    grow_delta_by_class = grow_delta_by_class,
    min_grow_threshold_value = 200, #200
    min_grow_threshold_by_class = min_grow_threshold_by_class,
    
    mask_edge_n_pixels = 1,
    min_seed_pixels_per_component = 30L,
    
    grow_engine = "whitebox",
    wbt_diag = FALSE,       # more conservative
    crop_to_units = TRUE,
    save_debug_rasters = TRUE
  )
  

```

## 2 & 3. Otsu segmentation refinement by AOIs and merging
```r

base_dir <- "ZENODO/exdata"
raster_path <- file.path(base_dir, "MinMin_2025_mosaic_res90m.tif")
corine_raster_path <- file.path(base_dir,"corine18_reclass.tif")
ecoregion_shapefile_path <- file.path(base_dir,"ecoregiones_olson.shp")
peninsula_shapefile <- file.path(base_dir,"iberian_peninsula_proj_final.shp")
ecoregion_field <- "EnZ_name"


base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
out_growth  <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

detected_polys <- file.path(out_growth,"BA_2025_OTSUGROW_CORI_ECOREG_otsu300_d100_seed30_whitebox.shp")

otsu_min_by_class <- list(corine = list( "4" = 350, "5" = 350, "9" = 350))
grow_delta_by_class <- list(corine = list(  "4" = 50, "5" = 100,"9" = 100))
min_grow_threshold_by_class <- list(corine = c( "4" = 300,"5" = 250,"9" = 250))


process_fun <- process_otsu_rasters_grow

res_ref <- segmentation_refinement(
  n_workers = 4,
  raster_path = raster_path,
  base_output_dir = file.path(out_growth, "REFINE_AOI_30_30_classes"),
  year = 2025,

  # AOI source
  detected_polys = detected_polys,
  reference_polys = NULL,
  binary_final_raster = NULL,

  # CRS / AOI
  target_crs = "EPSG:3035",
  aoi_buffer_m = 5000,
  omission_buffer_m = 0,
  top_n_if_no_ref = Inf,
  min_detected_area_m2 = 10,
  clip_aois_to_raster_extent = TRUE,
  merge_overlaps = TRUE,
  merge_overlaps_buffer_m = 0,
  verbose = TRUE,

  # Aux inputs (global; wrapper will crop/align per AOI)
  corine_raster_path = corine_raster_path,
  corine_classes = NULL,
  ecoregion_shapefile_path = ecoregion_shapefile_path,
  ecoregion_field = ecoregion_field,

  # Main processor
  process_fun = process_fun,

  # Output hygiene
  write_aoi_per_folder_shp = TRUE,
  standardize_ba_output = TRUE,
  keep_original_ba = FALSE,     # <- avoids duplicates: original output is moved/renamed to BA_AOI_###.shp
  stop_on_unused_args = TRUE,   # <- keep TRUE until everything runs cleanly

  # Forwarded args to process_fun
  rescue_args = list(
    trim_percentiles = list(min = 0.05, max = 0.99),
    otsu_value_range = c(0, 1500),
    otsu_thresholds = 300,

    otsu_min_by_class = otsu_min_by_class,
    grow_delta = 100,                 # default for classes not listed
    min_grow_threshold_value = 200,
    grow_delta_by_class = grow_delta_by_class,

    mask_edge_n_pixels = 1L,
    min_seed_pixels_per_component = 30L,
    grow_engine = "whitebox",
    wbt_diag = FALSE,

    ecoregion_touches = TRUE,
    segment_by_intersection = TRUE,
    crop_to_units = TRUE,

    save_debug_rasters = TRUE,
    tile = FALSE,

    # used by the wrapper as project_res for AOI projected rasters
    resolution = 90
  )
)

final_merged <- merge_aoi_shapefiles(
    base_output_dir = file.path(out_growth, "REFINE_AOI_30_30_classes"),
    pattern = "^BA_AOI_\\d+\\.shp$",
    out_path = file.path(out_growth, "BA_2025_REFINE_MERGED_300_100_seed30_30_classes.gpkg"),
    dissolve = FALSE,
    verbose = TRUE
  )


```


## 4. Burned area flagging I: Internal scoring by spatial issues
```r

  burnable_corine_path <- "ZENODO/exdata/corine18_burnable_mask_binary.tif"
  base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
  out_growth  <- file.path(base_out, "CORINE_ECOREG_seed_30pix")
  
  polys_stage1_path <- file.path(out_growth, "BA_2025_OTSUGROW_CORI_ECOREG_otsu300_d100_seed30_whitebox.shp")
  polys_stage2_path <- file.path(out_growth, "BA_2025_REFINE_MERGED_300_100_seed30_30_classes.gpkg")
  
  base_dir <- "ZENODO/exdata"
  rbr_path <- file.path(base_dir, "MinMin_2025_mosaic_res90m.tif")
  
  # If you have reference burned polygons as file, use path. If already loaded (ref_effis), skip this.
  ref_effis_path <- file.path(base_dir,"Effis_CA_2025.shp")
  ref_effis_pre_path <- file.path(base_dir,"Effis_CA_2024.shp")
  
  # Load
    ref_effis <- sf::st_read(ref_effis_path, quiet = TRUE)
    ref_effis <- sf::st_make_valid(ref_effis)

    ref_effis_pre <- sf::st_read(ref_effis_pre_path, quiet = TRUE)
    ref_effis_pre <- sf::st_make_valid(ref_effis_pre)
  
    burnable_corine <- terra::rast(burnable_corine_path)
  
  polys_stage1 <- sf::st_read(polys_stage1_path, quiet = TRUE)
  polys_stage2 <- sf::st_read(polys_stage2_path, quiet = TRUE)
  
  polys_stage1 <- sf::st_make_valid(polys_stage1)
  polys_stage2 <- sf::st_make_valid(polys_stage2)
  
  rbr_rast <- terra::rast(rbr_path)
  
  # Ensure single-layer RBR (if multiband, pick the right band)
  if (terra::nlyr(rbr_rast) > 1) {
    names(rbr_rast) <- make.names(names(rbr_rast), unique = TRUE)
    # Choose the first layer by default; change if needed:
    rbr_rast <- rbr_rast[[1]]
  }
  
  
source(file.path("R", "4_scoring_burned_area_wrapper.R"))
  
# Output dir
out_score <- file.path(out_growth, "phase1_scoring_internal_and_ref_classes_30")
dir.create(out_score, recursive = TRUE, showWarnings = FALSE)

year <- 2025

res <- scoring_burned_area_stage2(
  # -------------------------
  # Internal scoring
  # -------------------------
  polys_stage2     = polys_stage2,
  polys_stage1     = polys_stage1,
  burnable_corine  = burnable_corine,
  burnable_classes = NULL,          # NULL if burnable_corine is 0/1
  rbr_rast         = rbr_rast,

  corine_na_scope    = "all",
  max_corine_na_frac = 0.65,

  support_buffer_m = 0,

  # -------------------------
  # Pre-year erase (applies to internal and reference if score_ref_polys=TRUE)
  # -------------------------
  preyear_polys       = ref_effis_pre,
  erase_mask_buffer_m = 90,
  erase_post_shave_m  = 90,
  erase_min_area_m2   = 20000,

  # -------------------------
  # Reference polys (EFFIS) -> enables:
  #   (a) scoring_reference_burned_area()  [if score_ref_polys=TRUE]
  #   (b) scoring_validation_burned_area() [if ref_polys is not NULL]
  # -------------------------
  ref_polys          = ref_effis,
  score_ref_polys    = TRUE,   # <- KEY: computes comparable scoring on ref_effis
  ref_buffer_m       = 90,
  internal_keep_value = "keep",

  # Validation: uses "clean" geometries (post-erase) by default
  validate_use_clean      = TRUE,
  validate_ref_use_clean  = TRUE,
  validate_ref_use_kept   = TRUE,  # use ref_kept* (filtered by burnable/NA) instead of ref_flagged*

  # -------------------------
  # Saving
  # -------------------------
  save_outputs = TRUE,
  out_dir      = out_score,
  prefix       = paste0("BA_", year),
  driver       = "GPKG",
  overwrite    = TRUE,
  quiet        = TRUE,
  save_splits  = FALSE,

  suppress_sf_warnings = TRUE
)

# 
res$internal$out_paths$container    # -> .../BA_2025_phase1.gpkg
res$reference$out_paths$container   # -> .../BA_2025_reference.gpkg
res$validation$out_paths$container  # -> .../BA_2025_validation.gpkg

```



## 5. Burned area flagging II: internal and reference scoring by RBR percentiles
```r

base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring  <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

gpkg_phase1 <- file.path(out_scoring, "phase1_scoring_internal_and_ref_classes_30/BA_2025_validation.gpkg")
sf::st_layers(gpkg_phase1)
stage2_keep_review <- sf::st_read(gpkg_phase1, layer = "internal_flagged", quiet = TRUE)

base_dir <- "ZENODO/exdata"
rbr_path <- file.path(base_dir, "MinMin_2025_mosaic_res90m.tif")
rbr_rast <- terra::rast(rbr_path)
rbr_rast <- rbr_rast[[1]]


  res_internal <- score_rbr_keep_classes(
  polys_sf  = stage2_keep_review,
  rbr_rast  = rbr_rast,
  out_dir   = file.path(out_growth, "phase2_internal_scoring_classes_30"),
  prefix    = "internal_2025_p10_p5",

  # keep pool
  keep_pool = NULL,
  use_keep_common_pool = TRUE,
  keep_common_col   = "flag_internal_effis",
  keep_common_value = "keep_common",
  keep_fallback_col = "flag_internal",
  keep_fallback_value = "keep",

  # sampling
  max_keep_samples = 5e6,

  # ref quantile & decision thresholds
  above_keep_prob  = 0.05, #minimum quantile
  promote_p_above_ref = 0.50, #in percentage
  
  # decision thresholds
  promote_percentile  = 0.10,

  # quality thresholds
  min_area_ha = 10,
  min_pix = NULL,
  

  # class3
  do_class3 = TRUE,
  out_col = "class3_id",
  good_percentile = 0.05,
  bad_percentile  = 0.01,
  require_conf_ok = TRUE,
  
  fill_na_corine = TRUE,
  fill_na_corine_value = 0,

  # performance
  chunk_size = 500,
  verbose = TRUE,

  # saving + arcgis
  save_outputs = TRUE,
  driver = "GPKG",
  overwrite = TRUE,
  quiet = TRUE,
  save_splits = FALSE,
  arcgis_fix = TRUE,
  output_epsg = 3035
)
  
ref_effis_flagged <- sf::st_read(gpkg_phase1, layer = "ref_flagged", quiet = TRUE)
keep_pool <- res_internal$keep_pool
  
res_external <- score_rbr_keep_classes(
    polys_sf  = ref_effis_flagged, #ref_effis,   # your EFFIS sf to score
    rbr_rast  = rbr_rast,
    out_dir   = file.path(out_growth, "phase2_external_scoring_classes_30"),
    prefix    = "effis_2025_p10_p5",
    
    keep_pool = keep_pool,     # << reuse
    
    do_class3 = TRUE,
    out_col = "class3_id",
    good_percentile = 0.05,
    bad_percentile  = 0.01,
    require_conf_ok = TRUE,
    
    save_outputs = TRUE,
    driver = "GPKG",
    overwrite = TRUE,
    quiet = TRUE
  )
  
  effis_scored <- res_external$scored_sf
  


```

## 6. Calculate Polygon Metrics and Apply Geometric Filters
```r

library(sf)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

phase2_int_dir <- file.path(out_scoring, "phase2_internal_scoring_classes_30")
phase2_ext_dir <- file.path(out_scoring, "phase2_external_scoring_classes_30")

gpkg_phase2_int <- file.path(phase2_int_dir, "internal_2025_p10_p5_scored.gpkg")
gpkg_phase2_ext <- file.path(phase2_ext_dir, "effis_2025_p10_p5_scored.gpkg")

# Layer name in your scored gpkg
internal_layer <- "scored"
external_layer <- "scored"

# ------------------------------------------------------------
# 1) INTERNAL metrics
# ------------------------------------------------------------
internal_scored <- sf::st_read(gpkg_phase2_int, layer = internal_layer, quiet = TRUE)

# Create a single-layer gpkg (optional but OK)
metrics_in_int <- file.path(phase2_int_dir, "internal_scored_only.gpkg")
if (file.exists(metrics_in_int)) unlink(metrics_in_int, force = TRUE)
sf::st_write(internal_scored, metrics_in_int, layer = internal_layer, delete_dsn = TRUE, quiet = TRUE)

res_metrics_internal <- calculate_polygon_metrics(
  shapefile_paths = metrics_in_int,
  input_layer     = internal_layer,   # IMPORTANT for GPKG

  output_dir    = file.path(phase2_int_dir, "poly_metrics"),
  output_format = "gpkg",
  gpkg_layer    = "scored_metrics",

  dissolve     = TRUE,
  filter_logic = NULL,

  compute_area      = TRUE,
  compute_bbox      = TRUE,
  compute_perim     = TRUE,
  compute_mrr       = TRUE,
  compute_compact   = TRUE,
  compute_rectfill  = TRUE,

  keep_legacy = FALSE,

  # Keep ALL original attributes + metrics
  # If calculate_polygon_metrics2 expects a PATH, use metrics_in_int instead of internal_scored
  overlay_polygons_path = internal_scored,
  columns_to_keep       = NULL
)

# ------------------------------------------------------------
# 2) EXTERNAL metrics
# ------------------------------------------------------------
external_scored <- sf::st_read(gpkg_phase2_ext, layer = external_layer, quiet = TRUE)

metrics_in_ext <- file.path(phase2_ext_dir, "external_scored_only.gpkg")
if (file.exists(metrics_in_ext)) unlink(metrics_in_ext, force = TRUE)
sf::st_write(external_scored, metrics_in_ext, layer = external_layer, delete_dsn = TRUE, quiet = TRUE)

res_metrics_external <- calculate_polygon_metrics(
  shapefile_paths = metrics_in_ext,
  input_layer     = external_layer,   # IMPORTANT for GPKG

  output_dir    = file.path(phase2_ext_dir, "poly_metrics"),
  output_format = "gpkg",
  gpkg_layer    = "scored_metrics",

  dissolve     = TRUE,
  filter_logic = NULL,

  compute_area      = TRUE,
  compute_bbox      = TRUE,
  compute_perim     = TRUE,
  compute_mrr       = TRUE,
  compute_compact   = TRUE,
  compute_rectfill  = TRUE,

  keep_legacy = FALSE,

  # If calculate_polygon_metrics2 expects a PATH, use metrics_in_ext instead of external_scored
  overlay_polygons_path = external_scored,
  columns_to_keep       = NULL
)
# Optional: quick access to written outputs
# res_metrics_internal[[metrics_in_int]]$metrics
# res_metrics_external[[metrics_in_ext]]$metrics


```





## 7. OPTIONAL: Burned area flagging III: merging internal and reference fire polys
```r

base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring  <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

gpkg_phase2_int <- file.path(out_scoring, "phase2_internal_scoring_classes_30/internal_2025_p10_p5_scored.gpkg")
sf::st_layers(gpkg_phase2_int)
internal_scored <- sf::st_read(gpkg_phase2_int, layer = "scored", quiet = TRUE)

gpkg_phase2_ref <- file.path(out_scoring, "phase2_external_scoring_classes_30/effis_2025_p10_p5_scored.gpkg")
sf::st_layers(gpkg_phase2_ref)
effis_scored <- sf::st_read(gpkg_phase2_ref, layer = "scored", quiet = TRUE)

out_dir_merge <- file.path(out_scoring, "phase3_merge_internal_external")
dir.create(out_dir_merge, recursive = TRUE, showWarnings = FALSE)
  
  res <- merge_internal_with_effis_flags(
    internal_scored = internal_scored,
    effis_scored    = effis_scored,
    out_dir         = out_dir_merge,
    gpkg_name       = "merged_internal_reference_scored.gpkg",
    
    internal_flag_col = "class3_id",
    effis_flag_col    = "class3_id",
    
    assign_rule = "priority",
    priority_levels = c("good_identified", "ambiguous", "bad_identified"),
    
    buffer_m = 0,
    arcgis_fix = TRUE,
    output_epsg = 3035,
    
    overwrite = TRUE,
    quiet = TRUE,
    verbose = TRUE
  )

  
```

## 8. OPTIONAL: Histograms and statistics of Burned area flagging: internal and reference fire polys
```r

base_out    <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring  <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

gpkg_phase2_int <- file.path(out_scoring, "phase2_internal_scoring_classes_30/internal_2025_p10_p5_scored.gpkg")
sf::st_layers(gpkg_phase2_int)
internal_scored <- sf::st_read(gpkg_phase2_int, layer = "scored", quiet = TRUE)

gpkg_phase2_ref <- file.path(out_scoring, "phase2_external_scoring_classes_30/effis_2025_p10_p5_scored.gpkg")
sf::st_layers(gpkg_phase2_ref)
effis_scored <- sf::st_read(gpkg_phase2_ref, layer = "scored", quiet = TRUE)


rbr_path <- "ZENODO/exdata/MinMin_2025_mosaic_res90m.tif"
rbr_rast <- terra::rast(rbr_path)
rbr_rast <- rbr_rast[[1]]

out_dir_merge <- file.path(out_scoring, "phase4_histograms_tests")
dir.create(out_dir_merge, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "8_histograms_tests.R"))
  
 # 1) median_rbr existss in both layers
res <- histograms_internal_external(
  internal_sf = internal_scored,
  ref_sf      = effis_scored,
  rbr_rast    = rbr_rast,
  out_dir     = out_dir_merge,
  value_col   = "median_rbr",
  verbose     = TRUE
)

out_dir_merge1 <- file.path(out_scoring, "phase4a_histograms_tests")
dir.create(out_dir_merge1, recursive = TRUE, showWarnings = FALSE)

# 2) If RBFR median doesnÃƒâ€šÃ‚Â´t exist or is empty-> use raster "rbr_rast""
res <- phase4_diagnostics_internal_external(
  internal_sf = internal_scored,
  ref_sf      = effis_scored,
  rbr_rast    = rbr_rast,          # terra::SpatRaster
  out_dir     = out_dir_merge1,
  value_cols  = c("median_rbr"), #c("median_rbr", "p_above_keep_ref", "area_ha"),
  do_histograms = TRUE,
  do_tests      = TRUE,
  min_n_hist    = 10,
  min_n_test    = 20
)

out_dir_merge2 <- file.path(out_scoring, "phase4b_histograms_tests")
dir.create(out_dir_merge2, recursive = TRUE, showWarnings = FALSE)

# 3) With class-based comparisons (class3_id) and pairs
res <- phase4_diagnostics_internal_external(
  internal_sf = internal_scored,
  ref_sf      = effis_scored,
  rbr_rast    = rbr_rast,
  out_dir     = out_dir_merge2,
  value_cols  = c("median_rbr"),

  # flags
  internal_flag_col = "flag_internal",  # or NULL if flag_internal / flag_rbr exists
  internal_keep_value   = "keep",
  internal_review_value = "review",
  ref_flag_col = "flag_ref",
  ref_keep_value   = "keep_ref",
  ref_review_value = "review_ref",

  # keep_common preferred via cross-flag
  use_internal_effis_flag = TRUE,
  internal_effis_col = "flag_internal_effis",
  keep_common_value = "keep_common",
  review_int_no_effis_value = "review_int_no_effis",
  make_extra_no_effis = TRUE,

  # classes
  compare_class3 = TRUE,
  class_col_internal = "class3_id",
  class_col_ref      = "class3_id",
  hist_class_pairwise = TRUE,
  hist_class_pairwise_mode ="all_pairs", # or "cross_dataset"
  hist_class_pairwise_max_pairs = 200,

  test_class_pairwise = "same_class_internal_vs_ref",  # or "all"
  drop_na_class = TRUE
)



  
```

## 9 & 10. Validate Detected Burned Areas with Reference Data: Preparation & Validation
```r

################################
## PREPARING DATA FOR VALIDATION
## A) Fixed AOI, fixed burnable raster, year 2025
################################

base_dir    <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring <- file.path(base_dir, "CORINE_ECOREG_seed_30pix")

validation_dir <- file.path(out_scoring, "PREP_VALIDATION")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

ref_dir <- file.path(
  out_scoring,
  "phase2_external_scoring_classes_30",
  "poly_metrics"
)

det_dir <- file.path(
  out_scoring,
  "phase2_internal_scoring_classes_30",
  "poly_metrics"
)

ref_path <- file.path(
  ref_dir,
  "external_scored_only_orig_metrics.gpkg"
)

det_path <- file.path(
  det_dir,
  "internal_scored_only_orig_metrics.gpkg"
)

aoi_shapefile <- "ZENODO/exdata/iberian_peninsula_proj_final.shp"

burnable_raster <- "ZENODO/exdata/corine18_burnable_mask_binary.tif"

res_paths <- prepare_fire_polys_for_validation(
  ref_path = ref_path,
  det_path = det_path,
  out_dir  = validation_dir,
  
  # --- CRS / year context
  target_epsg = 3035,
  year_target = 2025,
  multi_year = FALSE,
  use_year_subdirs = FALSE,
  
  # --- Layers
  ref_layer = "scored_metrics",
  det_layer = "scored_metrics",
  
  # --- Fixed AOI
  aoi_shapefile = aoi_shapefile,
  use_origin_aoi = FALSE,
  aoi_mask_dir = NULL,
  
  # --- Fixed burnable raster
  do_burnable_clip = TRUE,
  burnable_raster = burnable_raster,
  burnable_mask_dir = NULL,
  dissolve_by_id = TRUE,
  crop_burnable_to_aoi = TRUE,
  binary_burnable = TRUE,
  burnable_classes = NULL,
  
  # --- Reference filters: summer DOY
  apply_ref_doy_filter = TRUE,
  ref_doy_col = "start_doy",
  ref_doy_min = 152,
  ref_doy_max = 274,
  
  # --- Optional attribute filters
  apply_ref_attr_filter = FALSE,
  ref_filter_col = c("class3_id"),
  ref_filter_values = c("good_identified"),
  
  apply_det_filter = FALSE,
  det_filter_col = c("class3_id"),
  det_filter_values = c("good_identified", "ambiguous"),
  
  # --- SAFE erase:
  # Remove detected polygons that intersect reference polygons dropped
  # by filtering or clipping, unless they also intersect final kept refs.
  remove_detected_if_intersects_dropped_ref = TRUE,
  dropped_ref_stage = "aoi",
  dropped_ref_buffer_m = 0,
  
  # --- Output names
  ref_output_tag = "REF_2025_D152_274_burneable_class",
  det_output_tag = "DET_2025_burneable_class_all_polys",
  output_prefix = "PREP_",
  
  # --- Outputs
  save_aoi_mask = FALSE,
  save_aoi_only = FALSE,
  save_burnable = TRUE,
  save_ratio_csv = FALSE,
  save_rasterized = FALSE,
  
  out_scheme = "default",
  out_patterns = NULL,
  
  overwrite = TRUE,
  disable_s2 = TRUE,
  on_empty = "stop",
  verbose = TRUE
)

res_paths



################################
## PREPARING DATA FOR VALIDATION
## B) Pre-2000 origin-based AOI
################################

ref_dir <- "ZENODO/exdata/Validation_fires_burneable_verano"

det_dir <- "ZENODO/exdata/Detected_fires_1984_1999"

validation_dir <- "ZENODO/exdata/PREP_VALIDATION_1984_1999"

dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

borders_shp <- "ZENODO/exdata/iberian_peninsula_comunidades.shp"

burnable_raster <- "ZENODO/exdata/corine90_burnable_mask_binary.tif"

origin_to_nameunit <- c(
  "Andalucia" = "AndalucÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â­a",
  "CataluÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â±a"  = "CataluÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â±a/Catalunya",
  "CLM"       = "Castilla-La Mancha",
  "Madrid"    = "Comunidad de Madrid",
  "Navarra"   = "Comunidad Foral de Navarra",
  "Valencia"  = "Comunitat Valenciana",
  "Landsat"   = "Portugal"
)

res_paths <- prepare_fire_polys_for_validation(
  ref_path = ref_dir,
  det_path = det_dir,
  out_dir  = validation_dir,
  
  # --- Multi-year mode
  multi_year = TRUE,
  use_year_subdirs = TRUE,
  year_regex = "\\d{4}",
  
  # --- CRS
  target_epsg = 3035,
  
  # --- Layers
  ref_layer = NULL,
  det_layer = "scored",
  
  # --- Origin-based AOI only before 2000
  aoi_shapefile = NULL,
  use_origin_aoi = TRUE,
  origin_aoi_before_year = 2000L,
  
  borders_shapefile = borders_shp,
  borders_nameunit_col = "NAMEUNIT",
  origin_col = "origin",
  origin_to_nameunit = origin_to_nameunit,
  resolve_unmatched_origin_by_overlap = TRUE,
  origin_fallback = "stop",
  
  # --- No standard AOI mask picker for pre-2000 in this run
  aoi_mask_dir = NULL,
  
  # --- Fixed burnable raster
  burnable_raster = burnable_raster,
  burnable_mask_dir = NULL,
  do_burnable_clip = TRUE,
  dissolve_by_id = TRUE,
  crop_burnable_to_aoi = TRUE,
  binary_burnable = TRUE,
  burnable_classes = NULL,
  
  # --- Reference summer DOY filter
  apply_ref_doy_filter = TRUE,
  ref_doy_col = "start_doy",
  ref_doy_min = 152,
  ref_doy_max = 274,
  
  apply_ref_attr_filter = FALSE,
  apply_det_filter = FALSE,
  
  remove_detected_if_intersects_dropped_ref = TRUE,
  dropped_ref_stage = "aoi",
  dropped_ref_buffer_m = 0,
  
  output_prefix = "PREP_",
  
  save_aoi_mask = TRUE,
  save_aoi_only = FALSE,
  save_burnable = TRUE,
  save_ratio_csv = FALSE,
  save_rasterized = FALSE,
  
  out_scheme = "default",
  overwrite = TRUE,
  disable_s2 = TRUE,
  on_empty = "stop",
  verbose = TRUE
)

res_paths

#############################
### VALIDATION #####
#############################

# --- your paths (adapt) ---
base_out      <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring   <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

#validation_dir <- file.path(out_scoring, "VALIDATION_BURNEABLE_LULC_GOOD_AMB")
validation_dir <- file.path(out_scoring, "VALIDATION_BURNEABLE_LULC_ALL_POLYS")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

prep_dir <- file.path(out_scoring, "PREP_VALIDATION")

# detected polygons (prediction)
map3035_path <- file.path(prep_dir, "PREP_DET_2025_burneable_class_all_polys__AOI_BURNABLE.gpkg")

# reference polygons (EFFIS, etc.)
ref3035_path <- file.path(prep_dir, "PREP_REF_2025_D152_274_burneable_class__AOI_BURNABLE.gpkg")

# mask AOI
mask_3035_path <- "ZENODO/exdata/iberian_peninsula_proj_final.shp"

# class raster + LUT  
class_raster_path <- "ZENODO/exdata/corine18_burnable_classes.tif"   #( ALL classes: corine18_reclass.tif; burneable: corine18_burnable_classes.tif)
class_lut_path    <- "ZENODO/exdata/LUT/corine_strata8_v1_min.csv"

res_val <- validate_fire_maps(
  input_shapefile = map3035_path,
  ref_shapefile   = ref3035_path,
  mask_shapefile  = mask_3035_path,
  year_target     = 2025,
  validation_dir  = validation_dir,

  # --- layers (GPKG): the best is to leave it as NULL for auto-detection ---
  input_layer = NULL,
  ref_layer   = NULL,

  # --- main controls ---
  metrics_type = "all",
  template_res = 90,
  min_detected_percent = 5,
  threshold_completely_detected = 50,
  min_area_reference_ha = 1,
  force_reprocess_ref  = TRUE,
  force_reprocess_pred = TRUE,
  union_detections_for_area = TRUE,
  arcgis_wkt1 = TRUE,

  # --- output naming ---
  output_prefix = "VAL",
  output_tag    = "2025_burneable_LULC_all_polys",
  tag_as_subdir = TRUE,
  outputs_subdir = "VAL",
  errors_subdir  = "ERR",
  write_error_vectors = TRUE,

  # --- classwise (CORINE strata) ---
  class_raster = class_raster_path,
  class_band   = 1,
  class_values = NULL,
  class_lut    = class_lut_path,
  class_domain = "class",
  class_na_as_zero = TRUE,
  class_chunk_rows = 2048,
  class_min_ref_area_ha = 0,
  class_include_all_lut_ids = TRUE,
  write_classwise_metrics = TRUE,
  write_classwise_errors  = TRUE,

  # --- NEW: confusion outputs ---
  write_confusion_vectors = TRUE,
  write_confusion_pixel_raster = TRUE,
  compute_det_overlap_stats = TRUE,   # optional
  confusion_subdir = NULL,          # optional
  confusion_polys_layer = "confusion_polys",  # optional
  
  dissolve_input_by = NULL, #"ID",   # <-- your real field
  dissolve_ref_by   = NULL,   # or "ID_ref" if you also want to rebuild the reference

  # --- rounding ---
  digits_metrics = 3,
  digits_area    = 2,
  digits_res     = 0,

  # --- XLSX summary ---
  write_xlsx = TRUE,
  xlsx_path = NULL,
  xlsx_overwrite = TRUE,

  verbose = TRUE
)

```

## 11. Calculate Validation statistics (separability)
```r

base_out      <- "ZENODO/exdata/2025/OTSU_MIN"
out_scoring   <- file.path(base_out, "CORINE_ECOREG_seed_30pix")

## VALIDATE ONLY GOOD & AMBIGUOUS FIRE POLYS
#validation_dir <- file.path(out_scoring, "VALIDATION_BURNEABLE_LULC_GOOD_AMB/VAL/VAL_2025_burneable_LULC_good_amb_polys/CONF")

## VALIDATE ALL FIRE POLYS
validation_dir <- file.path(out_scoring, "VALIDATION_BURNEABLE_LULC_ALL_POLYS/VAL/VAL_2025_burneable_LULC_all_polys/CONF")
polys <- file.path(validation_dir,"VAL_2025_confusion_polys_2025_burneable_LULC_all_polys.gpkg")

rbr_fire <- "ZENODO/exdata/MinMin_2025_mosaic_res90m.tif"

res <- validation_statistics(
  polys = polys,
  gpkg_layer = "confusion_polys",
  x_var = "confusion",
  y_mode = "attribute", #"raster" if you have multiple rasters
  y_var  = c("median_rbr", "area_ha", "elong_mrr"),
  y_name = c("RBR", "Area (ha)", "Elongation"),
   log_area_var = "area_ha",
  log_area = TRUE,            # apply log10 to Area_ha
  facet_scales = "fixed", 
  panel_vars = NULL,              # <- without panels. With panels by c("source"),
  out_dir = validation_dir,
  file_prefix = "by_confusion"
)
  


```
## Acknowledgements

We gratefully acknowledge funding from project INFORICAM (PID2020-119402RB-I00), funded by the Spanish MCIN/AEI/ 10.13039/501100011033 and by the "European Union NextGenerationEU/PRTR". Carlos Silva was supported by the NASA's Carbon Monitoring System funding (CMS, grant 22-CMS22-0015).

## Reporting Issues

Please report any issue regarding the OtsuFire package to Dr. Olga Viedma ([olga.viedma@uclm.es](mailto:olga.viedma@uclm.es)).

## Citing OtsuFire

Quintero, N.; Viedma, O.; Achour, H.; and Moreno, J. M. (2025). OtsuFire: Fire Scars and Severity Mapping Using 'Otsu' Thresholding. Version 0.1.4, accessed on June 24, 2025. Available at: <https://cran.r-project.org/web/packages/OtsuFire/index.html>.

## Disclaimer

**OtsuFire package comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**










