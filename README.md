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

<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/flow_chart_OtsuFire.png" width="700"/>
</p>

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

#### Open Anaconda Prompt or CMD

conda install -c conda-forge gdal

#### Check version

gdalinfo –version

#### Confirm availability

where gdalinfo

#### Optional: Create a Dedicated GDAL Environment

conda create –name gdal_env -c conda-forge gdal conda activate gdal_env
gdalinfo –version

#### Validate GDAL in Python

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

# Install 'OtsuSeg' dependency (if not already installed)
if (!requireNamespace("OtsuSeg", quietly = TRUE)) {
  remotes::install_github("olgaviedma/OtsuSeg")
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
library(OtsuSeg)
library(stats)
library(tidyr)
library(readr)
library(PMCMRplus)
library(multcompView)

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

mosaic_reproject_resample(
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

<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig3_RBR_2012_burned_corine_ge0.png" width="900"/>
</p>

## 3. Calculate Polygon Metrics and Apply Geometric Filters

``` r
folder_path <- "ZENODO/exdata"
burned_files <- list.files(folder_path, pattern = "burned_areas_2012_otsu_corine_.*\\.shp$", full.names = TRUE)

calculate_polygon_metrics(
  shapefile_paths = burned_files,
  output_dir = folder_path,
  area_min_ha = 10,
  bbox_h_min = 630, #NULL
  mnbbx_wd_min = 820, #NULL
  p_w_ratio_min = 4.99, #NULL
  h_w_ratio_min = 0.35 #NULL
)

```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig4_RBR_2012_burned_corine_ge0_metrics.png" width="900"/>
</p>

## 4. Detect Regeneration Areas using Otsu Thresholding on Negative RBR

``` r
## RENAME THE LANDSAT RBR-DOY MOSAIC FOR A CORRECT USE OF THE FUCNTION
file.rename(
  from = file.path(folder_path, "IBERIAN_MinMin_all_year_2012_mosaic_res90m.tif"),
  to   = file.path(folder_path, "IBERIAN_MinMin_2012_all_year_mosaic_masked_res90m.tif")
)


folder_path <- "ZENODO/exdata"

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

process_otsu_regenera(
    rbr_post = list(
        P1 = list.files(folder_path, pattern = "IBERIAN_MinMin_2013.*\\.tif$", full.names = TRUE),
        P2 = list.files(folder_path, pattern = "IBERIAN_MinMin_2014.*\\.tif$", full.names = TRUE)
    ),
    output_dir = folder_path,
    fire_year = 2012,
    rbr_date = list.files(
        folder_path,
        pattern = "IBERIAN_MinMin_2012.*\\.tif$",
        full.names = TRUE
    ),
    regen_year = c(1, 2),
    use_fixed_threshold = TRUE,
    fixed_threshold_value = -100,
    trim_percentiles = NULL,
    bind_all = FALSE,
    python_exe = python_exe,
    gdal_polygonize_script = gdal_polygonize_script
)

```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig5_burned_areas_corine_ge0_filter_all_regeneraP1P2.png" width="900"/>
</p>

## 5. Flag Burned Polygons as Regenerating or Not

``` r
folder_path <- "ZENODO/exdata"

flag_otsu_regenera(
  burned_files = list.files(folder_path, pattern = "*_filt_all.shp$", full.names = TRUE),
  regenera_files = list.files(folder_path, pattern = "burned_areas_2012_regenera.*\\.shp$", full.names = TRUE),
  min_regen_ratio = c(0.10, 0.20, 0.30),
  min_regen_ratio_P1 = 0.10,
  remove_no_regenera = TRUE,
  remove_condition = "year1_and_year2",
  all_years_vector = NULL, #  c("P1", "P2") 
  output_dir = folder_path
)

```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig6_cleaning_burned_polygons_regeneraY2.png" width="900"/>
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig7_rebuilding_burned_polygons_regeneraY1.png" width="700"/>
</p>


## 6. Extract DOY Statistics from Raster for Each Polygon

``` r
folder_path <- "ZENODO/exdata"

python_exe = "C:/ProgramData/anaconda3/python.exe"
gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py"

raster_path <- list.files(
  folder_path,
  pattern = "IBERIAN_MinMin_2012_.*\\.tif$",
  full.names = TRUE
)

calculate_doy_flags(
  raster = terra::rast(raster_path),
  doy_band = 2,
  polygons_sf <- list.files(
  folder_path,
  pattern = "burned_areas_2012_.*_rat30_P1P2_new_polys\\.shp$", #"burned_areas_2012_.*_rat[0-9]+_P1P2_new_polys\\.shp$"
  full.names = TRUE),
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
  folder_path,
  pattern = glob2rx("burned_areas_2012_otsu_*_thr100_P1P2_rat30_*.shp"),
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
  validation_dir = "folder_path",
  binary_burnable = TRUE,
  burnable_classes = NULL,
  buffer = 0,
  threshold_completely_detected = 75,
  min_area_reference_ha = 2,
  use_gdal = TRUE,
  python_exe = "C:/ProgramData/anaconda3/python.exe",
  gdal_polygonize_script = "C:/ProgramData/anaconda3/Scripts/gdal_polygonize.py",
  force_reprocess_ref = TRUE,
  metrics_type = "area"  # Options: "all", "pixel", or "area"
)

```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig8_cleaning_Portugal_burned_polygons_regeneraY12.png" width="900"/>
</p>
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig9_accuray_area_polygons_metrics.PNG" width="900"/>
</p>

## 8. Plotting Burned polygons as Regenerating or Not based on RBR values and CORINE classes

``` r
# INPUTS
folder_path <- "ZENODO/exdata"
corine_raster <- file.path(folder_path, "burneable_classes_def1_corine06_ETRS89.tif")

polygons_sf <- list.files(
  folder_path,
  pattern = glob2rx("burned_areas_2012_otsu_corine_ge0_metrics_filt_all_reg_thr100_P1P2_rat30.shp"),
  full.names = TRUE
)

rbr_years <- list(
  fire = file.path(folder_path, "IBERIAN_MinMin_2012_all_year_mosaic_masked_res90m.tif"),
  p1 = list.files(folder_path, pattern = "IBERIAN_MinMin_2013.*\\.tif$", full.names = TRUE),
  p2 = list.files(folder_path, pattern = "IBERIAN_MinMin_2014.*\\.tif$", full.names = TRUE)
)

corine_ras <- terra::rast(corine_raster)

for (burned_shp in polygons_sf) {
  message("Processing shapefile: ", burned_shp)
  
  burned_sf <- sf::st_read(burned_shp, quiet = TRUE)
  burned_vect <- terra::vect(burned_sf)
  burned_vect$ID <- seq_len(nrow(burned_vect))
  
  result_list <- list()
  
  for (yr in names(rbr_years)) {
    message("  Extracting RBR for year: ", yr)
    rbr_ras <- terra::rast(rbr_years[[yr]])[[1]]
    rbr_aligned <- terra::project(rbr_ras, corine_ras, method = "bilinear")
    
    vals_rbr <- terra::extract(rbr_aligned, burned_vect, fun = mean, na.rm = TRUE)
    vals_corine <- terra::extract(corine_ras, burned_vect)
    names(vals_corine)[2] <- "CORINE"
    
    corine_df <- as.data.frame(vals_corine) |>
      dplyr::group_by(ID) |>
      dplyr::summarise(CORINE_class = {
        ux <- na.omit(unique(CORINE))
        if (length(ux) == 0) NA_integer_ else ux[which.max(tabulate(match(CORINE, ux)))]
      })
    
    rbr_summary <- dplyr::tibble(
      ID = vals_rbr$ID,
      mean_RBR = vals_rbr[[2]],
      RBR_year = yr
    ) |>
      dplyr::left_join(corine_df, by = "ID")
    
    result_list[[yr]] <- rbr_summary
  }
  
  final_df <- dplyr::bind_rows(result_list)
  
  rbr_wide <- final_df |>
    dplyr::select(ID, RBR_year, mean_RBR) |>
    tidyr::pivot_wider(names_from = RBR_year, values_from = mean_RBR)
  
  corine_class <- final_df |>
    dplyr::select(ID, CORINE_class) |>
    dplyr::distinct()
  
  df <- burned_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(ID = seq_len(n())) |>
    dplyr::left_join(rbr_wide, by = "ID") |>
    dplyr::left_join(corine_class, by = "ID") |>
    dplyr::mutate(
      regen_P1 = ifelse(rgnr__P1 == "regenera", "Yes", "No"),
      regen_P2 = ifelse(rgnr__P2 == "regenera", "Yes", "No")
    ) |>
    tidyr::pivot_longer(cols = c(fire, p1, p2), names_to = "RBR_year", values_to = "RBR") |>
    tidyr::pivot_longer(cols = c(regen_P1, regen_P2), names_to = "regen_phase", values_to = "Regenerated") |>
    dplyr::filter(
      (RBR_year == "p1" & regen_phase == "regen_P1") |
      (RBR_year == "p2" & regen_phase == "regen_P2")
    ) |>
    dplyr::filter(!is.na(RBR))
    
    df <- df %>%
  dplyr::mutate(rgnr_f= as.factor(rgnr_f_))


# Levels
df_valid <- df %>%
  dplyr::filter(!is.na(RBR), !is.na(CORINE_class), !is.na(RBR_year)) %>%
  dplyr::mutate(RBR_year = factor(RBR_year, levels = c("fire", "p1", "p2")))


group_plots <- list()

for (regen_status in unique(df_valid$rgnr_f)) {
  df_group <- df_valid %>% filter(rgnr_f == regen_status)

  # Summarize and letters for CORINE_class
  all_letters <- list()
  summary_data <- df_group %>%
    group_by(CORINE_class, RBR_year) %>%
    summarise(
      median = median(RBR, na.rm = TRUE),
      sd = sd(RBR, na.rm = TRUE),
      .groups = "drop"
    )

  for (class_val in unique(df_group$CORINE_class)) {
    df_sub <- df_group %>% filter(CORINE_class == class_val)

    if (length(unique(df_sub$RBR_year)) < 2) next

    df_sub <- df_sub %>%
      filter(!is.na(RBR_year), !is.na(RBR)) %>%
      mutate(RBR_year = factor(RBR_year, levels = c("fire", "p1", "p2")))

    kw <- try(
      kwAllPairsDunnTest(RBR ~ RBR_year, data = df_sub, p.adjust.method = "bonferroni"),
      silent = TRUE
    )

    if (inherits(kw, "try-error")) next

    pvals <- kw$p.value
    valid <- which(!is.na(as.vector(pvals)))
    if (length(valid) == 0) next

    pvals_vector <- as.vector(pvals)[valid]
    comp_names <- apply(expand.grid(rownames(pvals), colnames(pvals)), 1, function(x) paste(sort(x), collapse = "-"))[valid]
    names(pvals_vector) <- comp_names

    sig_letters <- multcompLetters(p.adjust(pvals_vector, method = "bonferroni"))
    letter_df <- data.frame(
      CORINE_class = class_val,
      RBR_year = names(sig_letters$Letters),
      letter = sig_letters$Letters
    )

    all_letters[[as.character(class_val)]] <- letter_df
  }

  # Join letters with data
  letter_all <- bind_rows(all_letters)

  plot_df <- summary_data %>%
    left_join(letter_all, by = c("CORINE_class", "RBR_year")) %>%
    filter(!is.na(median))

  # Plot panel
  p <- ggplot(plot_df, aes(x = RBR_year, y = median, fill = RBR_year)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0.2) +
    geom_text(aes(y = median + 1.5 * sd, label = letter), size = 4) +
    facet_wrap(~ CORINE_class, ncol = 5) +
    labs(
      title = paste("RBR median ± SD - Group:", regen_status),
      x = "Year",
      y = "RBR median"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")

  # Save
  out_file <- file.path(folder_path, paste0(
    tools::file_path_sans_ext(basename(burned_shp)),
    "_RBR_panel_CORINE_", regen_status, ".png"
  ))

  ggsave(out_file, p, width = 12, height = 6)
  message("Guardado: ", out_file)
}
}

```
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig10_burned_areas_2012_otsu_corine_ge0_metrics_filt_all_reg_thr100_P1P2_rat30_RBR_regenera.png" width="900"/>
</p>
<p align="center">
  <img src="https://raw.githubusercontent.com/olgaviedma/OtsuFire/main/README/Fig11_burned_areas_2012_otsu_corine_ge0_metrics_filt_all_reg_thr100_P1P2_rat30_RBR_no_regenera.png" width="900"/>
</p>

