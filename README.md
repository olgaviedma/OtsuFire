<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/OtsuFire_logo.png" width="400"/>
</p>

<h1 align="center">OtsuFire: Fire Mapping and Regeneration Assessment Toolkit</h1>
================

``` r
# Carga el paquete desde el código fuente del proyecto
devtools::load_all(".")
```

**Authors:** Natalia Quintero, Olga Viedma, Hammadi Achour, and Jose
Manuel Moreno

Automated tools to map fire scars using Landsat-based RBR/dNBR
composites and Otsu thresholding. Supports large-area mosaics,
polygonization, filtering, regeneration assessment, and validation
workflows.

------------------------------------------------------------------------

## 🔧 System Requirements

### 🐍 Python & GDAL

Several functions in `OtsuFire` call the GDAL Python script
`gdal_polygonize.py`. To use these functions, ensure:

- **Python** is installed (e.g., via Anaconda or Miniconda).
- **GDAL** is installed in your Python environment:
  - Via **Conda**: `conda install -c conda-forge gdal`
  - Or via **pip**: `pip install gdal` (less recommended)
- You know the full paths to:
  - Your Python executable (`python.exe`)
  - The `gdal_polygonize.py` script (typically in `Scripts/` on Windows)

🔧 **Example configuration:**

``` r
python_exe <- "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"
```

### 🛠️ GDAL Installation Guide

### Open Anaconda Prompt or CMD

conda install -c conda-forge gdal

### Check version

gdalinfo –version

### Confirm availability

where gdalinfo

### Optional: Create a Dedicated GDAL Environment

conda create –name gdal_env -c conda-forge gdal conda activate gdal_env
gdalinfo –version

### Validate GDAL in Python

python

import osgeo.gdal print(osgeo.gdal.\_\_version\_\_)

### ✅ Additional Notes

- All input rasters must be projected and masked consistently.
- Functions support large-area mosaics, tiling, and region-specific
  thresholding.
- Always ensure gdal_polygonize.py is accessible and executable.

## Getting Started

## Installation

``` r
#The CRAN version:
install.packages("OtsuFire")

# The development version:
#install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/OtsuFire", dependencies = TRUE)

# Install 'otsuSeg' dependency (if not already installed)
if (!requireNamespace("otsuSeg", quietly = TRUE)) {
  remotes::install_github("olgaviedma/otsuSeg")
}


```

## Libraries

``` r
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
library(otsuSeg)
library(stats)
```

## 📥 Download and Prepare Example Data from Zenodo

This package includes examples that require raster and vector data
stored externally on Zenodo. These files are not included in the CRAN
version of the package to avoid exceeding size limits.

The following code downloads and extracts a ZIP archive into a local
folder called `ZENODO/exdata`, located in the root of your working
directory. This folder is ignored by R during package build and check
processes (via `.Rbuildignore`).

Once downloaded and extracted, the data can be used in examples
throughout this README or in your local scripts.

The following script:

- Downloads the dataset from Zenodo (only if not already downloaded)
- Validates the download
- Extracts the content
- Places the extracted files in a consistent folder structure

## 0. Import data from ZENODO

``` r
library(curl)

# Define Zenodo download URL and local paths
zenodo_url <- "https://zenodo.org/records/15332851/files/DATA.zip?download=1"
output_dir <- "ZENODO/exdata"
zip_file <- file.path(output_dir, "DATA.zip")

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Set up extended timeout using curl handle
h <- curl::new_handle()
curl::handle_setopt(h, timeout = 900)  # Timeout in seconds (30 minutes)

# Attempt download if file does not already exist
if (!file.exists(zip_file)) {
  message("Downloading DATA.zip from Zenodo...")
  tryCatch({
    curl::curl_download(url = zenodo_url, destfile = zip_file, handle = h, quiet = FALSE)
    message("Download complete.")
  }, error = function(e) {
    message("Download failed: ", conditionMessage(e))
    if (file.exists(zip_file)) file.remove(zip_file)
  })
} else {
  message("DATA.zip already exists. Skipping download.")
}

# unzip if download was successful
if (file.exists(zip_file)) {
  unzip(zip_file, exdir = output_dir)
  message("Unzipped contents to ", output_dir)
}

if (file.exists(zip_file)) {
  file.remove(zip_file)
  message("DATA.zip has been removed.")
} else {
  message("DATA.zip does not exist.")
}
```

## 0. Create a burnable mask from CORINE maps

``` r
folder_path <- "ZENODO/exdata"
# Define input CORINE rasters (year-labeled)
 corine_rasters <- list(
   "corine06" = "ZENODO/exdata/U2012_CLC2006_V2020_20u1.tif")

 # Path to shapefile of the Iberian Peninsula
 peninsula_shp <- "ZENODO/exdata/iberian_peninsula_proj_final.shp"

 # Output folders
 output_raster_dir <- folder_path
 output_vector_dir <- folder_path

 # Paths to required GDAL utilities
 gdalwarp_path <- "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe"
 python_exe <- "C:/ProgramData/anaconda3/python.exe"
 gdal_polygonize <- "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

 # Run the function
 generate_burnable_mask(
   corine_rasters = corine_rasters,
   peninsula_shapefile = peninsula_shp,
   output_raster_dir = output_raster_dir,
   output_vector_dir = output_vector_dir,
   gdalwarp_path = gdalwarp_path,
   python_exe = python_exe,
   gdal_polygonize_script = gdal_polygonize,
   reproject = TRUE,
   res = 90,
   to_wgs84 = TRUE,
   vectorize = TRUE
 )
 
 ## PLOTTING

# Path to the reclassified raster
r <- terra::rast("ZENODO/exdata/burneable_classes_def1_corine06_ETRS89.tif")

# Definir clases y colores
corine_classes <- data.frame(
  value = c(22, 23, 24, 25, 26, 27, 28, 29, 32, 33),
  label = c(
    "Agro-forestry areas",
    "Broad-leaved forest",
    "Coniferous forest",
    "Mixed forest",
    "Natural grasslands",
    "Moors and heathland",
    "Sclerophyllous vegetation",
    "Transitional woodland-shrub",
    "Sparsely vegetated areas",
    "Burnt areas"
  )
)

colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
  "#e5c494", "#b3b3b3", "#a65628", "#ff6666", "#999999"
)

par(mar = c(1, 1, 3, 10))  # bottom, left, top, right

# Plot raster without legend
plot(r,
     col = colors,
     breaks = corine_classes$value - 0.5,
     main = "CORINE Burnable Classes (Reclassified)", axes = FALSE, legend = FALSE)

# Add the legend
legend("bottomright",
       legend = corine_classes$label,fill = colors, cex = 0.7,
       bty = "n", inset = c(-0.05, 0.04), xpd = TRUE)            
```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig1_CORINE_burneable_classes.png" width="500"/>
</p>

## 1. Mosaic and Resample Landsat composites (raster data). Optionally, mask by CORINE map and study area

``` r
## if masked with burneable CORINE classes from GEE
folder_path     <- "ZENODO/exdata"
gdalwarp_path   <- "C:/ProgramData/anaconda3/Library/bin/gdalwarp.exe"
gdal_calc_path   <- "C:/ProgramData/anaconda3/Scripts/gdal_calc.py"
gdal_merge_path <- "C:/ProgramData/anaconda3/Scripts/gdal_merge.py"

mosaic_tiles_max_clean(
  folder_path = "ZENODO/exdata",
  year = 2012,
  raster_pattern = "IBERIAN_MinMin_all_year_2012_*.tif",
  mask_path = file.path(folder_path, "iberian_peninsula_proj_final.shp"),
  crs_target = "EPSG:3035",
  res_target = 90,
  nodata_value = -9999
)


## if Landsat Composites are NOT masked with the burneable CORINE classes from GEE. Masked them!

folder_path <- "ZENODO/exdata"
mosaic_path <- file.path(folder_path, "IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif")
mask_raster_path <- file.path(folder_path, "burneable_mask_binary_corine06_ETRS89.tif")
shapefile_clip <- file.path(folder_path, "iberian_peninsula_proj_final.shp")

mask_to_mosaic(
  mosaic_path      = mosaic_path,
  mask_raster_path = mask_raster_path,
  shapefile_clip   = shapefile_clip
)

##PLOTTING
r <- rast(mosaic_path)
summary(r)
rbr <- r[[1]]
doy <- r[[2]]

# layout configuration
par(mfrow = c(1, 2), mar = c(4, 4, 4, 5)) 

# Reclasificación
rbr_classes <- terra::classify(rbr, rbind(
  c(0,   100,   1),   # Low
  c(100, 300,   2),   # Medium
  c(300, 500,   3),   # High
  c(500, 9000,  4)    # Extreme
))

# Convertir a factor
rbr_classes <- terra::as.factor(rbr_classes)

# Establecer tabla de niveles (data.frame con columnas 'value' y 'label')
lvls <- data.frame(value = 1:4, label = c("Low", "Medium", "High", "Extreme"))
levels(rbr_classes) <- lvls


# Colores y layout
cols <- c("darkgreen", "orange", "red", "brown")
par(mfrow = c(1, 2), mar = c(4, 4, 4, 5))
plot(rbr_classes, col = cols, main = "RBR Composite (90m)")
plot(doy, main = "DOY Composite (90m)", col = hcl.colors(100, "Viridis"), maxcell = 500000)
```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig2_mosaicking_masking1.png" width="900"/>
</p>


## 2. Fire mapping: Apply Otsu Thresholding to RBR Mosaic Raster

``` r
folder_path <- "ZENODO/exdata"
raster_path = file.path(folder_path, "IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif")
corine_raster_path = file.path(folder_path,"burneable_classes_def1_corine06_ETRS89.tif")

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

process_otsu_rasters(
  raster_path = raster_path,
  output_dir = folder_path,
  year= 2012,
  otsu_thresholds = c(0),
  use_original = FALSE,
  python_exe = python_exe,
  gdal_polygonize_script = gdal_polygonize_script,
  tile = TRUE,
  index_type = NULL,
  corine_raster_path = corine_raster_path,
  ecoregion_shapefile_path = NULL,
  min_otsu_threshold_value = 100
)
```

## 3. Calculate Polygon Metrics and Apply Geometric Filters

``` r
folder_path <- "ZENODO/exdata"
burned_files <- list.files(folder_path, pattern = "burned_areas_2012_otsu_*.shp$", full.names = TRUE)

calculate_polygon_metrics(
  shapefile_paths = burned_files,
  output_dir = folder_path,
  area_min_ha = 10,
  bbox_h_min = NULL, #630
  mnbbx_wd_min = NULL, #820
  p_w_ratio_min = NULL, #4.49
  h_w_ratio_min = NULL, #0.25
)
```

## 4. Detect Regeneration Areas using Otsu Thresholding on Negative RBR

``` r
folder_path <- "ZENODO/exdata"

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

process_otsu_regenera(
  rbr_post = list(
    P1 = file.path(folder_path,"IBERIAN_MinMin_2013_all_year_mosaic_masked_res90m.tif"),
    P2 = file.path(folder_path,"IBERIAN_MinMin_2014_all_year_mosaic_masked_res90m")
  ),
  output_dir = folder_path,
  fire_year = 2012,
  rbr_date = file.path(folder_path,"IBERIAN_MinMin_2012_all_year_mosaic_masked_res90m.tif"),
  regen_year = c(1, 2),
  use_fixed_threshold = TRUE,
  fixed_threshold_value = -100,
  trim_percentiles = NULL,
  python_exe = "path/to/python.exe",
  gdal_polygonize_script = "path/to/gdal_polygonize.py"
)
```

## 5. Flag Burned Polygons as Regenerating or Not

``` r
folder_path <- "ZENODO/exdata"

flag_otsu_regenera(
  burned_files = list.files(folder_path, pattern = "*_filt_area.shp$", full.names = TRUE),
  regenera_files = list.files(folder_path, pattern = "*_thresh100.shp$", full.names = TRUE),
  min_regen_ratio = 0.10,
  remove_no_regenera = TRUE,
  remove_condition = "year1_and_year2",
  all_years_vector = NULL, #c("P1", "P2")
  output_dir = folder_path
)
```

## 6. Extract DOY Statistics from Raster for Each Polygon

``` r
folder_path <- "ZENODO/exdata"

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

raster_path = file.path(folder_path,"IBERIAN_MinMin_2012_all_year_mosaic_masked_res90m.tif")

calculate_doy_flags(
  raster = terra::rast(raster_path),
  doy_band = 2,
  polygons_sf = list.files(folder_path, pattern = "*_rat*_filter.shp$", full.names = TRUE),
  output_dir = folder_path,
  year = 2012,
  doy_thresholds = c(10),
  polygonize = TRUE,
  stats = "mode",
  keep_all_polygons = TRUE,
  calc_percentiles = TRUE,
  percentiles = c(0.05, 0.25, 0.95),
  python_exe = python_exe,
  gdal_polygonize_script = gdal_polygonize_script
)
```

## 7. Validate Detected Burned Areas with Reference Data

``` r
folder_path <- "ZENODO/exdata"

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

polygons_sf <- list.files(
  burned_dir,
  pattern = glob2rx("burned_areas_2012_otsu_*_thr100_P1P2_rat*.shp"),
  full.names = TRUE
)

ref_shapefile = file.path(folder_path,"fires_2012.shp")
mask_shapefile = file.path(folder_path,"PORTUGAL_SHAPE.shp")
burnable_raster = file.path(folder_path,"burneable_mask_binary_corine06_ETRS89.tif")

validate_fire_maps(
  input_shapefile = polygons_sf,
  ref_shapefile = ref_shapefile,
  mask_shapefile = mask_shapefile,
  burnable_raster = burnable_raster,
  year_target = 2012,
  validation_dir = folder_path,
  binary_burnable = TRUE,
  burnable_classes = NULL,
  min_area_reference_ha = 2,
  buffer = 0,
  threshold_completely_detected = 90,
  force_reprocess_ref = TRUE,
  use_gdal = TRUE,
  python_exe = python_exe,
  gdal_polygonize_script = gdal_polygonize_script
)
```
