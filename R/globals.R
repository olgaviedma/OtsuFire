# Block 7: centralised globals + missing importFrom for the migrated
# engine code. The migrated engine uses unquoted column names with dplyr
# verbs, raw stats/utils functions, and a few helpers that R CMD check
# cannot resolve from static analysis. Declaring them here keeps the
# engine code byte-identical while satisfying CRAN-style checks.

#' @importFrom stats aggregate median predict runif sd weighted.mean IQR
#' @importFrom utils capture.output str tail write.csv installed.packages
#' @importFrom dplyr across arrange count desc distinct full_join
#'   if_else rename row_number transmute slice_sample case_when
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
NULL

# Variables used by dplyr/data.table NSE inside the migrated engine.
# These do not exist as R bindings — they are column names referenced
# inside the data.frame argument to the verb. Declaring them silences
# "no visible binding for global variable" NOTEs without changing
# scientific behaviour.
utils::globalVariables(c(
  ".", ".year", ".data",
  "a", "across", "aggregate", "arrange",
  "block_id",
  "class_audited", "class_final", "class_input", "count",
  "desc", "distinct",
  "eco_major", "elapsed_sec",
  "fire_uid", "fold", "fold_rep1",
  "get_fire_mapping_registry_paths",
  "has_oof",
  "hotspot_available",
  "hs_conf_mean", "hs_frp_max", "hs_frp_sum", "hs_hiConf_n",
  "hs_in_buffer", "hs_in_poly", "hs_min_dist_m", "hs_used_n",
  "j",
  "n_burned_units", "n_pos", "n_pos_blocks", "n_total",
  "p_burned", "poly_id",
  "raw_class", "rbr_aw_med", "rbr_med", "rename", "row_number",
  "source_set",
  "tibble", "transmute", "y",
  "geometry",
  # carryover from internal-sup-process-otsu-unburned.R top-level
  # globalVariables() that was inlined in roxygen — moved here.
  "ID", "P1_id", "abline", "area_ha", "area_m2", "axis",
  "burn_area_ha", "burned_id", "dev.off", "doy", "fid_final",
  "first", "glob2rx", "h_w_ratio", "int_area_ha", "legend",
  "mnbbx_wd", "mtext", "na.omit", "ornt_bbox", "p_w_ratio",
  "par", "png", "regen_area_ha", "regen_area_ha_sum", "setNames",
  "total_burned_area", "year", "%>%", ":=",
  "ecoregion_classes_internal", "unit_id2", "CORINE_CLASS",
  "CORINE_YEAR", "ECO_CLASS",
  # carryover from internal-det-segmentation-refinement.R
  "AREA_M2", "REFINE_REASON", "AOI_ID", "AOI_GROUP",
  # carryover from internal-det-process-otsu-grow.R
  "CLUMP_ID", "ECO_ID_INTERNAL", "UNIT_ID"
))
