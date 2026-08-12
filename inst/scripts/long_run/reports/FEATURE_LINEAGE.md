# Feature lineage: FULL vs NO-HOTSPOT (2017 balanced)

Derived from FROZEN code (OtsuFire 0.6.1 SOURCE at run time): R/internal-sup-residual-cols.R
(.supervised_feature_cols, 50-name whitelist in 0.6.2 after hs_any was removed),
R/internal-sup-extract-features.R (hotspot_features + doy_post=change-index band2).
NOTE: this run was produced under 0.6.1 where the whitelist still nominally listed
hs_any (#51, a GPKG-absent helper synthesised from hs_used_n>0). In 0.6.2 hs_any was
REMOVED from the whitelist and its create-matrix synthesis deleted (redundant, never
materialised); the resolved feature space is UNCHANGED (still 50 base + 50 _isNA).

FULL base = 50 (canonical resolved whitelist; hs_any was always naturally absent).
NO-HOTSPOT base = 38 = 50 minus the 12 hotspot-lineage features.

| Feature | Family | FULL | NO-HOTSPOT | Source |
|---|---|---|---|---|
| rbr_valid_frac | rbr_samewindow | in | in | change_index RBR summer band 1 (spectral) |
| rbr_p10 | rbr_samewindow | in | in | change_index RBR summer band 1 (spectral) |
| rbr_med | rbr_samewindow | in | in | change_index RBR summer band 1 (spectral) |
| rbr_p90 | rbr_samewindow | in | in | change_index RBR summer band 1 (spectral) |
| rbr_iqr | rbr_samewindow | in | in | change_index RBR summer band 1 (spectral) |
| cor_open_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_wetlands_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_water_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_herbaceous_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_urban_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_agri_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| cor_forest_frac | corine | in | in | CORINE land-cover fractions (zonal) |
| elev_valid_frac | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_p10 | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_med | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_p90 | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_iqr | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_mean | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| elev_sd | topo_elevation | in | in | DEM elevation zonal stats (topography) |
| slope_valid_frac | topo_slope | in | in | Slope zonal stats (topography) |
| slope_p10 | topo_slope | in | in | Slope zonal stats (topography) |
| slope_med | topo_slope | in | in | Slope zonal stats (topography) |
| slope_p90 | topo_slope | in | in | Slope zonal stats (topography) |
| slope_iqr | topo_slope | in | in | Slope zonal stats (topography) |
| slope_mean | topo_slope | in | in | Slope zonal stats (topography) |
| slope_sd | topo_slope | in | in | Slope zonal stats (topography) |
| doy_valid_frac | doy_spectral | in | in | change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS) |
| doy_p10 | doy_spectral | in | in | change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS) |
| doy_med | doy_spectral | in | in | change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS) |
| doy_p90 | doy_spectral | in | in | change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS) |
| doy_iqr | doy_spectral | in | in | change_index DOY band 2 (doy_post; SPECTRAL/temporal; STAYS) |
| rbr_aw_valid_frac | rbr_allwindow | in | in | delayed_change_index RBR all-window composite (spectral) |
| rbr_aw_p10 | rbr_allwindow | in | in | delayed_change_index RBR all-window composite (spectral) |
| rbr_aw_med | rbr_allwindow | in | in | delayed_change_index RBR all-window composite (spectral) |
| rbr_aw_p90 | rbr_allwindow | in | in | delayed_change_index RBR all-window composite (spectral) |
| rbr_aw_iqr | rbr_allwindow | in | in | delayed_change_index RBR all-window composite (spectral) |
| persist_delta | persistence | in | in | RBR same-window vs all-window persistence (spectral) |
| persist_ratio | persistence | in | in | RBR same-window vs all-window persistence (spectral) |
| hotspot_available | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_in_poly | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_in_buffer | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_used_n | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_min_dist_m | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_frp_sum | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_frp_max | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_conf_mean | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_hiConf_n | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_support_present | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_no_support_when_available | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |
| hs_only_buffer_support | hotspot | in | out | MODIS/active-fire hotspots (extract-features hotspot_features); REMOVED in NO-HOTSPOT |

FULL base=50  NO-HOTSPOT base=38  removed=12 (all hotspot family).
Each base feature gets an auto _isNA companion -> FULL=100 model cols, NO-HOTSPOT=76.
hs_any was REMOVED from the whitelist in 0.6.2 (redundant with hs_used_n>0; never a GPKG column, never synthesized into the resolved recipe); the NO-HOTSPOT removed-count is exactly 12 genuine hotspot features (no phantom 13th).
