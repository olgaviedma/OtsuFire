#' Configure the probabilistic refinement workflow
#'
#' @description
#' Create a configuration object for the probabilistic refinement workflow. The
#' function collects and validates the decision layer for Otsu-guided
#' patches, raster inputs, training-pool settings, feature controls, model
#' parameters, and output locations.
#'
#' Pass the resulting configuration to the individual probabilistic refinement workflow
#' functions or to [run_oneyear_supervised_pipeline()].
#'
#' This function prepares the configuration only. It does not build training
#' pools, extract features, train a model, or generate predictions.
#'
#' @param internal_decisions `sf` polygon layer or GeoPackage path. Required
#'   decision layer for Otsu-guided patches from the target year and
#'   scenario. Patches classified as `"keep"` provide the initial burned
#'   labels.
#' @param change_index `terra::SpatRaster` or raster path. Required annual
#'   immediate change-index input. The expected layout is band 1: summer RBR;
#'   band 2: post-fire day of year. Used for negative-pool construction and
#'   immediate-response features.
#' @param delayed_change_index `terra::SpatRaster`, raster path, or `NULL`.
#'   Delayed autumn-winter change-index input used for delayed-response and
#'   persistence features. If omitted, a conventional path may be resolved
#'   from `options$composite_base`; otherwise it remains `NULL` and features
#'   requiring this input are skipped.
#' @param immediate_change_index_path Alternative argument name for
#'   `change_index`. Supply one name for this input. Supplying both is
#'   accepted only when their values are identical.
#' @param delayed_change_index_path Alternative argument name for
#'   `delayed_change_index`. Supply one name for this input. Supplying both is
#'   accepted only when their values are identical.
#' @param hotspots `sf` point layer, file path, or `NULL`. Optional
#'   target-year active-fire observations used to build hotspot features. When
#'   `NULL`, those features are not built.
#' @param reference_burned_map `sf` polygon layer, raster, file path, or
#'   `NULL`. Optional external burned-area reference for subsequent map
#'   validation with [validate_fire_maps()]. It does not supply the supervised
#'   training labels.
#' @param peninsula_shapefile `sf` object, `terra::SpatVector`, file path, or
#'   `NULL`. Optional study-area boundary. Used by the negative-pool builder
#'   under the `"corine"` Otsu mode; typically unused with the fixed default
#'   mode, `"burnable_only"`.
#' @param topo `terra::SpatRaster`, raster path, or `NULL`. Two-band topography
#'   raster containing elevation in band 1 and slope in band 2. Used for
#'   topographic features.
#' @param corine_raster `terra::SpatRaster`, raster path, or `NULL`. CORINE
#'   land-cover raster used for land-cover features.
#' @param burnable_mask `terra::SpatRaster`, raster path, or `NULL`. Binary
#'   burnable-area mask used to construct background and Otsu-based negative
#'   pools.
#' @param run_label Non-empty character scalar used for output naming and
#'   routing. Has no methodological effect. Default: `"balanced"`.
#' @param target_year Integer scalar between `1900` and `2100`. Required. Used
#'   for output naming, temporal filtering, and conventional CORINE-epoch
#'   selection.
#' @param output_dir Character scalar. Root output directory. Default:
#'   `tempdir()`. Specify a persistent directory to retain results beyond the
#'   temporary session.
#' @param run_name Non-empty character scalar identifying the run. Default:
#'   `"supervised_burned_map"`.
#' @param flat_output_routes Logical scalar. Use the standard nested output
#'   layout when `FALSE`, or place stage folders directly under `output_dir`
#'   when `TRUE`. Default: `FALSE`.
#' @param min_burned_pool_n Non-negative integer. Minimum number of
#'   burned-label polygons required after quality checks. Training stops if
#'   fewer remain. Default: `5L`.
#' @param nrounds_max Maximum number of XGBoost boosting rounds. Default when
#'   `NULL`: `4000`.
#' @param early_stop Early-stopping patience, in boosting rounds. Default when
#'   `NULL`: `80`.
#' @param oof_seed_base Random seed for out-of-fold training. Default when
#'   `NULL`: `42`.
#' @param final_sampling_seed Random seed for final negative-pool sampling.
#'   Default when `NULL`: `42`.
#' @param final_seed Random seed for the final train/validation split and
#'   model training. Default when `NULL`: `42`.
#' @param random_seed Random seed for background-negative sampling. Affects
#'   which background examples are selected. Default when `NULL`: `42`.
#' @param val_frac Inner validation fraction used for early stopping. Must be
#'   between `0` and `1`, excluding the endpoints. Default when `NULL`:
#'   `0.15`.
#' @param group_col Column used to keep grouped observations together in the
#'   inner train/validation split. Default when `NULL`: `"block_id"`.
#' @param impute_numeric Numeric imputation rule: `"median"` or `"zero"`.
#'   Default when `NULL`: `"median"`.
#' @param impute_factor_missing Category used to represent missing factor
#'   values. Default when `NULL`: `"MISSING"`.
#' @param negative_pool_params Named list controlling negative-pool sampling,
#'   caps, and optional strongly burned-like negatives. See
#'   \strong{Negative-pool settings}. Default: `list()`.
#' @param runtime_options Named list containing `reuse_existing`,
#'   `write_outputs`, and `verbose`, each a logical scalar. All resolve to
#'   `TRUE` by default. Unknown names raise an error.
#' @param feature_whitelist_override Character vector restricting the model to
#'   a subset of the package feature set, or `NULL` to use the default
#'   permitted set. Feature availability also depends on supplied inputs and
#'   `include_shape_features`.
#' @param feature_weights Named numeric vector of per-feature weights, or
#'   `NULL` for equal weights. Applied consistently during out-of-fold and
#'   final training.
#' @param include_shape_features Logical scalar. Whether to compute and allow
#'   the optional shape and size features. Default: `FALSE`. See
#'   \strong{Shape and size features}.
#' @param model_params Named list of partial XGBoost parameter overrides, or
#'   `NULL` for package defaults. For example, `list(eta = 0.03)` changes only
#'   `eta`. Setting `scale_pos_weight` here is rejected because it is
#'   calculated during training.
#' @param options Named list of advanced paths and settings. See
#'   \strong{Advanced options}. Default: `list()`.
#'
#' @details
#' The two change-index argument pairs resolve to `config$inputs$change_index`
#' and `config$inputs$delayed_change_index`, respectively. Examples below use
#' `change_index` and `delayed_change_index`.
#'
#' An input that is optional when constructing the configuration may still be
#' required by a downstream stage or feature block.
#'
#' @section Workflow stages:
#' The configuration is used by the following stages:
#'
#' | Stage | Function |
#' |---|---|
#' | Build labelled and scoring pools | [build_supervised_training_pools()] |
#' | Assign spatial cross-validation folds | [make_spatial_folds()] |
#' | Extract polygon features | [extract_supervised_features()] |
#' | Generate out-of-fold diagnostics | [run_oof_diagnostics()] |
#' | Train the final model | [train_final_burned_model()] |
#' | Score candidate polygons | [score_supervised_burned_map()] |
#' | Validate the resulting map against an external reference | [validate_fire_maps()] |
#'
#' [run_oneyear_supervised_pipeline()] runs the first six stages from one
#' configuration.
#'
#' @section Training labels and pools:
#' The positive pool consists of Otsu-guided patches classified as `"keep"`,
#' represented as burned training labels. These are automatically derived
#' labels and are not equivalent to independent ground truth.
#'
#' Negative labels come from up to three pools:
#'
#' | Pool | Internal name | Description |
#' |---|---|---|
#' | Burnable background | `random` | Cells with weak immediate change-index responses, sampled within the burnable domain and outside a 500 m exclusion zone around Otsu-guided candidate patches. |
#' | Moderately burned-like | `otsu` | Otsu-guided patches rejected by the rule-based filters, with `S_PATCH_PA <= 0.15`. |
#' | Strongly burned-like | `artifact_hard` | Selected Otsu-guided patches classified as `"drop"`, with strong burned-like responses and additional evidence of potential artefacts. Disabled by default. |
#'
#' Background selection uses values at or below the configured
#' `rbr_quantile`, calculated within the sampling domain.
#'
#' Review-class candidates are withheld from training and retained for
#' subsequent scoring. Training requires explicit burned or unburned labels.
#' A candidate is not treated as unburned simply because it is absent from
#' the burned pool.
#'
#' When visual validation is applied through [apply_visual_validation()],
#' only confirmed candidates with `VISUAL = 1` are admitted through that
#' visual-review step.
#'
#' @section Negative-pool settings:
#' `negative_pool_params` accepts the following blocks. Unknown parameter
#' names raise an error.
#'
#' \strong{Background and Otsu settings}
#'
#' | Setting | Description | Default |
#' |---|---|---|
#' | `random$n_cells` | Positive integer specifying the background sampling size. | `1500L` |
#' | `random$rbr_quantile` | RBR quantile used to select background cells. Must be between `0` and `1`. | `0.50` |
#' | `otsu$candidate_threshold` | Otsu-residual severity threshold for candidates. | `0` |
#' | `otsu$reference_threshold` | Otsu-residual severity threshold for the reference. | `100` |
#' | `caps` | Named numeric vector with exactly `random` and `otsu`. Each value specifies the maximum pool size as a multiple of the burned-label count. Values must be non-negative; `Inf` disables a cap. | `c(random = 1, otsu = 1)` |
#'
#' For example:
#' \preformatted{
#' negative_pool_params = list(
#'   random = list(
#'     n_cells = 1500L,
#'     rbr_quantile = 0.50
#'   ),
#'   otsu = list(
#'     candidate_threshold = 0,
#'     reference_threshold = 100
#'   ),
#'   caps = c(
#'     random = 1.0,
#'     otsu = 1.0
#'   )
#' )
#' }
#'
#' The same capping policy is used during out-of-fold and final training,
#' applied to the relevant training subset.
#'
#' \strong{Strongly burned-like negatives}
#'
#' The optional `artifact_hard` block controls selection and weighting of
#' additional hard negatives.
#'
#' | Setting | Description | Default |
#' |---|---|---|
#' | `enabled` | Enable promotion of eligible candidates to unburned training labels. | `FALSE` |
#' | `total_weight_ratio` | Total hard-negative weight relative to total burned-pool weight. Must be non-negative. | `0.10` |
#' | `rbr_med_reference` | Pool used to calculate the RBR quantile threshold: `"negative"` or `"positive"`. | `"negative"` |
#' | `rbr_med_min_q` | Quantile used to establish the RBR floor. Must be between `0` and `1`. | `0.90` |
#' | `persist_ratio_max` | Persistence-ratio threshold used in eligibility rules. | `0.35` |
#' | `persist_delta_max` | Persistence-difference threshold used in eligibility rules. | `-100` |
#' | `area_ha_min` | Area threshold, in hectares, used in eligibility rules. | `500` |
#' | `doy_iqr_max` | Day-of-year interquartile-range threshold used in eligibility rules. | `1` |
#' | `reason_whitelist` | Classification-reason values for Otsu-guided patches used by the reason-based eligibility branch. | `character(0)` |
#'
#' See [promote_artifact_hard_negatives()] for how the eligibility conditions
#' are combined.
#'
#' All eligible hard negatives are retained under this policy; their
#' contribution is controlled by weight rather than by the `random` and
#' `otsu` caps.
#'
#' Under the pool-ratio weighting rule:
#' \preformatted{
#' total hard-negative weight =
#'   total_weight_ratio x total burned-pool weight
#' }
#'
#' The weight is distributed across the hard-negative examples. Thus,
#' `total_weight_ratio = 0.10` does not assign a weight of `0.10` to each
#' polygon.
#'
#' The default is a starting value, not a universally suitable balance.
#' Review the selected examples and evaluate the effect against independent
#' reference data.
#'
#' When `enabled = FALSE`, eligible candidates may still appear in audit
#' outputs, but they are not promoted to training through the hard-negative
#' mechanism.
#'
#' @section Model fitting and evaluation:
#' \strong{Out-of-fold and final training.} Both training stages use an inner
#' validation split to select the number of boosting rounds through early
#' stopping, followed by a refit on the available training set.
#'
#' For out-of-fold diagnostics, this procedure is repeated within each
#' spatial training fold. The held-out outer fold is reserved for prediction
#' and evaluation. Preprocessing, feature selection, and training must use
#' the corresponding training partition.
#'
#' For final training, the same procedure is applied to the eligible labelled
#' dataset, subject to sampling and caps.
#'
#' \strong{Feature preparation.} The feature recipe adds missingness
#' indicators, imputes numeric and categorical values, and fixes the
#' predictor-column order. The trained model and recipe are saved together so
#' that scoring uses the same transformations.
#'
#' \strong{Hotspot features.} To omit hotspot features, leave
#' `hotspots = NULL`. A feature whitelist can also explicitly exclude hotspot
#' predictors.
#'
#' Supplying hotspots makes their feature block available. It does not
#' require using those features if a whitelist excludes them.
#'
#' \strong{Shape and size features.} With `include_shape_features = TRUE`,
#' the following features are computed and made eligible for model use:
#' \preformatted{
#' c(
#'   "area_ha",
#'   "n_pix",
#'   "log_area",
#'   "perim_m",
#'   "compactness",
#'   "elongation"
#' )
#' }
#'
#' These features can encode differences in how training examples were
#' sampled. For example, small fixed background cells can be distinguished
#' from burned polygons by size and shape alone.
#'
#' Treat this block as experimental. Evaluate its effect using independent
#' map references as well as out-of-fold diagnostics.
#'
#' \strong{Scores and validation.} The scoring stage returns `p_burned`, a
#' model score that is not necessarily a calibrated probability. Select a
#' binary-map threshold with this distinction in mind.
#'
#' Out-of-fold diagnostics measure agreement with the internal training
#' labels. External validation with [validate_fire_maps()] assesses the
#' thresholded map against a separate reference, such as EFFIS.
#'
#' @section Consolidated training-pool output:
#' The workflow provides a consolidated GeoPackage named
#' `supervised_training_pool.gpkg`, with layer name `supervised_training_pool`,
#' for inspecting available examples, candidate hard negatives, and training
#' participation.
#'
#' Presence in this layer does not imply that a row was used for training.
#'
#' | Field or group | Meaning |
#' |---|---|
#' | `year`, `fire_uid`, `poly_id` | Year and example identifiers. |
#' | `training_label` | Effective training label; may be `NA` for candidates not assigned a training label. |
#' | `original_pool_source` | Provenance before hard-negative promotion. |
#' | `pool_source` | Effective source: `high_confidence_keep`, `random`, `otsu`, `artifact_hard`, or `deterministic_drop`. |
#' | `artifact_hard_eligible` | Whether the candidate satisfies the hard-negative selection rule. |
#' | `artifact_hard_enabled` | Whether hard-negative promotion is enabled in the configuration. |
#' | `artifact_hard_used` | Whether the candidate was effectively promoted through that mechanism. |
#' | `used_in_training` | Whether the row entered the final training set after caps and filters. |
#' | `sample_weight`, `source_total_weight` | Effective example weight and source-level weight information. |
#' | `deterministic_class`, `deterministic_reason` | Original classification and decision reason for the Otsu-guided patch. |
#' | `fold_id`, `has_oof_prediction`, `p_burned_oof` | Cross-validation assignment and prediction information. |
#' | `rbr_med`, `rbr_aw_med`, `persist_ratio`, `persist_delta`, `area_ha`, `doy_iqr` | Descriptive features used to inspect candidates. |
#' | `artifact_hard_branch` | Eligibility branch: `persist_delta`, `large_single_doy`, `reason_whitelist`, or `multiple`. |
#' | `artifact_hard_rbr_threshold` | RBR threshold used in candidate selection. |
#' | `artifact_hard_weight_ratio` | Configured pool-weight ratio, including for eligible candidates when promotion is disabled. |
#' | `review_priority`, `review_tier` | Visual-review priority score from 0 to 100 and corresponding tier. Higher priority indicates a greater concern that a proposed negative may be a real fire. |
#' | `label_confidence` | Heuristic label-confidence indicator. For hard-negative candidates, derived as `1 - review_priority / 100`. |
#' | `VISUAL` | Initially blank review field: `1` means the proposed label is confirmed; `0` means it is not confirmed. |
#'
#' `review_priority` and `label_confidence` are review aids, not calibrated
#' probabilities. For an unburned candidate, `VISUAL = 1` confirms the
#' proposed unburned label.
#'
#' Use the reported workflow outputs to locate the consolidated layer.
#'
#' @section Output locations and optional input paths:
#' With `flat_output_routes = FALSE`, outputs use
#' `<output_dir>/<target_year>/<run_name>/SUPERVISED/<run_label>/`.
#'
#' With `flat_output_routes = TRUE`, stage folders are placed directly under
#' `output_dir`.
#'
#' Some omitted inputs are resolved using conventional paths when the
#' corresponding base directory is supplied:
#'
#' | Input | Base option | Conventional relative path |
#' |---|---|---|
#' | `delayed_change_index` | `composite_base` | `Autumn/mean_mean_<target_year>_mosaic.tif` |
#' | `peninsula_shapefile` | `data_base` | `Borders/Iberian_peninsula.shp` |
#' | `topo` | `data_base` | `Topography/elevation_slope.tif` |
#' | `corine_raster` | `data_base` | `Corine_Masks/CLC_<corine_year>_peninsula.tif` |
#' | `burnable_mask` | `data_base` | `Corine_Masks/burneable_mask_binary_corine_<corine_year>_ETRS89.tif` |
#'
#' These names are package conventions. Supply explicit paths when your files
#' use different names or cover another study area. The conventional mask
#' filename retains the spelling `burneable`.
#'
#' The CORINE epoch is resolved from `target_year`. If the required base
#' option is absent, the optional input remains `NULL`.
#'
#' @section Advanced options:
#' Common entries in `options` include:
#'
#' | Option | Purpose |
#' |---|---|
#' | `data_base` | Root directory for conventional supporting-input paths. |
#' | `composite_base` | Root directory for conventional delayed-index paths. |
#' | `result_name` | Identifier used in convention-based paths. |
#' | `unburned_base_dir` | Shared root for negative-pool products. |
#' | `python_exe` | Path to the Python executable. |
#' | `gdal_polygonize_script` | Path to the GDAL polygonisation script. |
#' | `gdalwarp_path` | Path to `gdalwarp`. |
#' | `ogr2ogr_exe` | Path to `ogr2ogr`. |
#' | `tool_paths` | Nested list containing external-tool paths. |
#' | `currentyear_preyear_overlap_thr` | Temporal-adjustment overlap threshold. Default: `0.70`. |
#' | `currentyear_hotspot_density_thr` | Temporal-adjustment hotspot-density threshold. Default: `0.001`. |
#' | `currentyear_temporal_penalty_floor` | Temporal-adjustment penalty floor. Default: `0.10`. |
#'
#' When `unburned_base_dir` is supplied, shared outputs include
#' `<unburned_base_dir>/UNBURNED/<year>_unburned.gpkg` and
#' `<unburned_base_dir>/_OTSU_NEGATIVE/`.
#'
#' Reuse negative pools only when their source inputs and selection settings
#' are compatible. Matching the year alone does not establish compatibility.
#'
#' Use `negative_pool_params`, `random_seed`, and `runtime_options` for their
#' documented controls. Other low-level negative-pool settings remain fixed
#' internally.
#'
#' \strong{Compatibility notes}
#' * Negative-pool caps must be supplied through `negative_pool_params$caps`.
#'   Former top-level arguments such as `cap_random`, `cap_otsu`,
#'   `cap_contextual`, and `cap_spectral` are not supported.
#' * Legacy `options$unb_*` and `options$otsu_negative_*` settings are no
#'   longer read.
#' * The negative-pool policy is fixed to `"all_sources"`. The compatibility
#'   setting `options$negative_pool_policy = "all_sources"` has no effect;
#'   other values are rejected.
#' * Former engine-location options `engine_root`, `supervised_engine_root`,
#'   and `scripts_root` are no longer used.
#'
#' @return An S3 object of class `otsufire_supervised_burned_config`.
#'
#' | Field | Contents |
#' |---|---|
#' | `scenario`, `target_year`, `run_name` | Run identifiers. |
#' | `inputs` | Resolved spatial and labelled inputs. |
#' | `output_dir`, `output_routes` | Output root and stage locations. |
#' | `model_params` | Resolved XGBoost parameters. |
#' | `train_control` | Resolved training, sampling, seed, imputation, and feature settings. |
#' | `negative_pool_params` | Public negative-pool settings and caps. |
#' | `negative_pool_runtime` | Runtime controls for negative-pool processing. |
#' | `negative_pool_internal` | Fixed internal negative-pool settings. |
#' | `options`, `tool_paths` | Advanced settings and external-tool locations. |
#' | `resolved_params_provenance` | Information on how parameter values were resolved. |
#' | `negative_pool_policy` | Informational value, always `"all_sources"`. |
#' | `min_burned_pool_n` | Minimum required burned-pool size. |
#'
#' Use the constructor arguments to configure the workflow rather than
#' editing internal fields directly.
#'
#' @seealso [build_supervised_training_pools()], [make_spatial_folds()],
#'   [extract_supervised_features()], [run_oof_diagnostics()],
#'   [train_final_burned_model()], [score_supervised_burned_map()],
#'   [run_oneyear_supervised_pipeline()], [promote_artifact_hard_negatives()],
#'   [apply_visual_validation()], [validate_supervised_execution()],
#'   [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' # Configure a probabilistic refinement run with explicit input paths
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced",
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   reference_burned_map = "data/reference_burned_2022.gpkg",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022",
#'   negative_pool_params = list(
#'     caps = c(random = 1.0, otsu = 1.0)
#'   ),
#'   runtime_options = list(
#'     reuse_existing = TRUE,
#'     write_outputs = TRUE,
#'     verbose = TRUE
#'   )
#' )
#'
#' # Inspect the resolved settings
#' cfg$train_control
#' cfg$negative_pool_params
#' cfg$output_routes
#'
#' # Build the training pools
#' pools <- build_supervised_training_pools(cfg)
#'
#' # Configure a run without hotspot predictors
#' cfg_no_hotspots <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_1994.gpkg",
#'   change_index = "data/RBR_1994.tif",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_1994.tif",
#'   burnable_mask = "data/burnable_mask_1994.tif",
#'   hotspots = NULL,
#'   target_year = 1994L,
#'   output_dir = "results",
#'   run_name = "RBR_1994_no_hotspots",
#'   feature_whitelist_override = c(
#'     "rbr_med", "rbr_p90", "elev_med", "slope_med"
#'   )
#' )
#'
#' # Enable strongly burned-like negatives with an explicit pool weight
#' cfg_hard_negatives <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022_hard_negatives",
#'   negative_pool_params = list(
#'     artifact_hard = list(
#'       enabled = TRUE,
#'       total_weight_ratio = 0.10,
#'       rbr_med_reference = "negative",
#'       rbr_med_min_q = 0.90,
#'       persist_ratio_max = 0.35
#'     )
#'   )
#' )
#'
#' # Run the probabilistic refinement workflow
#' result <- run_oneyear_supervised_pipeline(cfg_hard_negatives)
#'
#' # Inspect the consolidated layer at its reported output location
#' pool <- sf::st_read(
#'   "path/to/supervised_training_pool.gpkg",
#'   layer = "supervised_training_pool",
#'   quiet = TRUE
#' )
#'
#' # Inspect eligible hard-negative candidates
#' eligible <- pool[
#'   which(pool$artifact_hard_eligible %in% TRUE),
#' ]
#' table(eligible$artifact_hard_branch, useNA = "ifany")
#'
#' # Inspect the examples used for final training
#' trained <- pool[
#'   which(pool$used_in_training %in% TRUE),
#' ]
#' table(trained$pool_source, useNA = "ifany")
#' }
#'
#' @family workflow
#' @export
build_supervised_burned_config <- function(
    run_label = "balanced",
    internal_decisions,
    change_index,
    delayed_change_index = NULL,
    # Gate 1B PIECE 3 (2026-06-07): canonical immediate/delayed change-index
    # ROUTES. `immediate_change_index_path` is the canonical name for the
    # required immediate change index; it is an alias of `change_index` (the
    # historical name, kept for back-compat AND shared with the deterministic
    # phase). `delayed_change_index_path` is the canonical name for the optional
    # delayed (autumn-winter) change index; alias of `delayed_change_index`.
    # Supplying BOTH forms of a route is allowed only if they are identical;
    # a genuine conflict ERRORS (never silently ignored). Internally both
    # resolve onto the single cfg$inputs$change_index / $delayed_change_index
    # fields, so there is ONE source of truth per route.
    immediate_change_index_path = NULL,
    delayed_change_index_path   = NULL,
    hotspots = NULL,
    reference_burned_map = NULL,
    target_year,
    output_dir = tempdir(),
    run_name = "supervised_burned_map",
    flat_output_routes = FALSE,
    min_burned_pool_n = 5L,
    # §N+25 (2026-06-05): supervised-RUN inputs wired to the config. All
    # default NULL; when NULL they fall back to the EXACT convention path
    # previously built by the orchestrator / unburned helpers, so passing
    # nothing reproduces the historical run byte-for-byte. The
    # VALIDATION-only inputs (strata / study-area mask / EFFIS reference)
    # are NOT part of the supervised run; they belong to the separate
    # validate_fire_maps() call, exactly like the deterministic phase.
    peninsula_shapefile = NULL,
    topo = NULL,
    corine_raster = NULL,
    burnable_mask = NULL,
    # ---- Gate 1B (2026-06-07): resolved methodological / training params ----
    # The cfg is the SINGLE SOURCE OF TRUTH for these. Every argument defaults
    # to NULL = "use the canonical default" (the canonical defaults live ONLY
    # in .of_canonical_train_control() / .of_canonical_model_params(), nowhere
    # else). A non-NULL value OVERRIDES the canonical default and is validated +
    # stored in cfg$train_control / cfg$model_params, which is what every
    # downstream stage then reads.
    nrounds_max                = NULL,
    early_stop                 = NULL,
    oof_seed_base              = NULL,
    final_sampling_seed        = NULL,
    final_seed                 = NULL,
    # GATE 6.7 (2026-06-12): `random_seed` is the negative-pool random-background
    # RNG seed. It lives in the CANONICAL SEEDS block (cfg$train_control$seeds)
    # alongside the three training seeds, NOT in runtime_options: it changes WHICH
    # negatives are selected, so it is methodological + reproducibility-affecting
    # and is folded into the methodological fingerprint. NULL = canonical 42.
    random_seed                = NULL,
    val_frac                   = NULL,
    group_col                  = NULL,
    impute_numeric             = NULL,
    impute_factor_missing      = NULL,
    # GATE 6.7 (2026-06-12): the typed PUBLIC negative-pool block. This is the
    # SINGLE source of truth for the two negative-bucket caps and the small set
    # of user-settable negative-pool methodological knobs. The former top-level
    # `cap_random` / `cap_otsu` builder args were REMOVED: passing them now hits
    # the unused-argument guard. See .of_resolve_negative_pool_params().
    negative_pool_params       = list(),
    # GATE 6.7 (2026-06-12): TECHNICAL runtime options (reuse_existing /
    # write_outputs / verbose). These do NOT alter the methodological fingerprint.
    runtime_options            = list(),
    feature_whitelist_override = NULL,
    feature_weights            = NULL,
    # OPTIONAL shape/size feature block (OtsuFire 0.12.0, OFF by default).
    # When TRUE it makes feature extraction COMPUTE the six shape columns
    # AND makes them ELIGIBLE model features (the active feature universe
    # becomes the canonical 50 + the 6 shape names). FALSE (default)
    # leaves the resolved feature universe and every model byte-identical
    # to 0.11.0. Experimental; see the SAMPLING-BIAS caveat in NEWS.md.
    include_shape_features     = FALSE,
    model_params               = NULL,
    options = list()
) {
  if (!is.character(run_label) || length(run_label) != 1L || is.na(run_label) ||
      !nzchar(run_label)) {
    stop("'run_label' must be a single non-empty character string.", call. = FALSE)
  }

  if (missing(target_year) || is.null(target_year)) {
    stop("'target_year' is required.", call. = FALSE)
  }
  target_year <- suppressWarnings(as.integer(target_year)[1L])
  if (is.na(target_year) || target_year < 1900L || target_year > 2100L) {
    stop("'target_year' must be an integer in [1900, 2100].", call. = FALSE)
  }

  if (missing(internal_decisions) || is.null(internal_decisions)) {
    stop("'internal_decisions' is required.", call. = FALSE)
  }

  # ---- Gate 1B PIECE 3: resolve the canonical change-index ROUTES -----------
  # Each route has a canonical name (immediate/delayed_change_index_path) and a
  # historical alias (change_index / delayed_change_index). Resolve each to a
  # single value: if only one form is given, use it; if BOTH are given they must
  # be identical (a genuine conflict ERRORS - never silently ignored). The
  # resolved value flows into the single cfg$inputs field for that route.
  .resolve_route <- function(canonical, alias, canonical_nm, alias_nm,
                             alias_missing) {
    has_canon <- !is.null(canonical)
    has_alias <- !alias_missing && !is.null(alias)
    if (has_canon && has_alias) {
      if (!identical(canonical, alias)) {
        stop("Conflicting change-index routes: '", canonical_nm,
             "' and its alias '", alias_nm, "' were both supplied with ",
             "different values. Pass only one (they name the same input).",
             call. = FALSE)
      }
      return(canonical)
    }
    if (has_canon) return(canonical)
    if (has_alias) return(alias)
    NULL
  }
  # immediate: canonical = immediate_change_index_path, alias = change_index (req)
  change_index <- .resolve_route(
    canonical = immediate_change_index_path,
    alias     = if (missing(change_index)) NULL else change_index,
    canonical_nm = "immediate_change_index_path", alias_nm = "change_index",
    alias_missing = missing(change_index)
  )
  # delayed: canonical = delayed_change_index_path, alias = delayed_change_index
  delayed_change_index <- .resolve_route(
    canonical = delayed_change_index_path,
    alias     = delayed_change_index,
    canonical_nm = "delayed_change_index_path", alias_nm = "delayed_change_index",
    alias_missing = FALSE
  )

  if (is.null(change_index)) {
    stop("'change_index' (a.k.a. 'immediate_change_index_path') is required.",
         call. = FALSE)
  }

  if (!is.character(run_name) || length(run_name) != 1L || !nzchar(run_name)) {
    stop("'run_name' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("'output_dir' must be a single non-empty character string.", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  min_burned_pool_n <- suppressWarnings(as.integer(min_burned_pool_n)[1L])
  if (is.na(min_burned_pool_n) || min_burned_pool_n < 0L) {
    stop("'min_burned_pool_n' must be a non-negative integer.", call. = FALSE)
  }

  if (!is.list(options)) {
    stop("'options' must be a list().", call. = FALSE)
  }
  if (length(options) > 0L && (is.null(names(options)) || any(!nzchar(names(options))))) {
    stop("'options' must be a named list (all elements must have names).", call. = FALSE)
  }

  # ---- §N+25: resolve the supervised-RUN input paths --------------------
  # Centralized convention construction (`.of_supervised_convention_paths`)
  # so the cfg, the orchestrator and the unburned helpers agree on EXACTLY
  # the same default paths. The convention requires `data_base` (and, for
  # delayed_change_index, `composite_base` + `result_name`), which live in
  # `options`. When `data_base` is absent (e.g. an ahead-of-batch config or a
  # unit test that does not exercise the run), the convention path cannot be
  # built; in that case the user-supplied value is normalized as-is and NULL
  # stays NULL, deferring resolution to runtime exactly as before.
  conv <- .of_supervised_convention_paths(
    data_base      = options$data_base,
    composite_base = options$composite_base,
    result_name    = options$result_name %||% "Min_Min",
    target_year    = target_year
  )

  delayed_change_index <- delayed_change_index %||% conv$delayed_change_index
  peninsula_shapefile  <- peninsula_shapefile  %||% conv$peninsula_shapefile
  topo                 <- topo                 %||% conv$topo
  corine_raster        <- corine_raster        %||% conv$corine_raster
  burnable_mask        <- burnable_mask        %||% conv$burnable_mask

  inputs <- list(
    internal_decisions   = .of_normalize_input_spec(internal_decisions,
                                                     "internal_decisions",
                                                     allow_null = FALSE),
    change_index         = .of_normalize_input_spec(change_index,
                                                     "change_index",
                                                     allow_null = FALSE),
    delayed_change_index = .of_normalize_input_spec(delayed_change_index,
                                                     "delayed_change_index",
                                                     allow_null = TRUE),
    hotspots             = .of_normalize_input_spec(hotspots, "hotspots",
                                                     allow_null = TRUE),
    reference_burned_map = .of_normalize_input_spec(reference_burned_map,
                                                     "reference_burned_map",
                                                     allow_null = TRUE),
    # §N+25 RUN inputs wired to cfg (convention defaults applied above).
    peninsula_shapefile  = .of_normalize_input_spec(peninsula_shapefile,
                                                     "peninsula_shapefile",
                                                     allow_null = TRUE),
    topo                 = .of_normalize_input_spec(topo, "topo",
                                                     allow_null = TRUE),
    corine_raster        = .of_normalize_input_spec(corine_raster,
                                                     "corine_raster",
                                                     allow_null = TRUE),
    burnable_mask        = .of_normalize_input_spec(burnable_mask,
                                                     "burnable_mask",
                                                     allow_null = TRUE)
  )

  # ---- §N+25: fail-fast validation of the newly-wired RUN inputs --------
  # When a wired RUN input resolves to an on-disk PATH spec, assert the file
  # exists NOW (config-construction is <1s) instead of 5-10 min into the run.
  # This only fires for path specs; in-memory objects and NULL (deferred) are
  # skipped and re-validated downstream. The convention default for a
  # `data_base`-less config is NULL, so such configs are not blocked here.
  # `internal_decisions` / `change_index` / `hotspots` keep their historical
  # config-time behaviour (existence enforced downstream, not here). The
  # VALIDATION-only inputs (strata / mask / EFFIS reference) are intentionally
  # absent here: they belong to the separate validate_fire_maps() call.
  # Gate 1B PIECE 4 (2026-06-08): the supervised ecoregion / corine_ecoregion
  # Otsu modes (and their `ecoregion_shapefile` cfg input) were removed. The
  # canonical negative-pool Otsu mode is `burnable_only`; `"corine"` is also
  # supported. CORINE×ecoregion stratification lives only in the deterministic
  # stage, which is configured separately.
  for (.nm in c("delayed_change_index", "peninsula_shapefile", "topo",
                "corine_raster", "burnable_mask")) {
    .of_check_input_file(inputs[[.nm]], .nm)
  }

  # BUG 3 Phase 1a (2026-06-05): the supervised engine is fully in-package, so
  # the engine-root / scripts-root discovery here was inert metadata that only
  # fed the now-removed external PROJECT_PATHS.R fallback. Removed the
  # `engine_root`, `supervised_engine_root` and `scripts_root` config fields,
  # their getwd()-walk finders (`.of_find_supervised_engine_root()`,
  # `.of_find_supervised_scripts_root()`) and the `.of_find_engine_root()`
  # call from the supervised path. `.of_find_engine_root()` itself is retained
  # in deterministic-config.R for the deterministic stage.

  # External tool paths (AS01, OtsuFire 0.3.0). Mirrors the deterministic
  # stage convention from `deterministic-config.R`. All NULL by default;
  # the Otsu residual negative pipeline (`all_sources` mode) errors loudly if
  # the path it needs is missing. Callers may pass either flat options
  # (e.g. options$python_exe) or a nested options$tool_paths list.
  tool_paths <- list(
    python_exe             = options$python_exe             %||%
                              options$tool_paths$python_exe,
    gdal_polygonize_script = options$gdal_polygonize_script %||%
                              options$tool_paths$gdal_polygonize_script,
    gdalwarp_path          = options$gdalwarp_path          %||%
                              options$tool_paths$gdalwarp_path,
    ogr2ogr_exe            = options$ogr2ogr_exe            %||%
                              options$tool_paths$ogr2ogr_exe
  )

  # Negative-pool policy (§N+26, 2026-06-05): FULLY REMOVED as a user knob.
  # The supervised pipeline is ALWAYS `all_sources` (deterministic drops +
  # random burnable background + Otsu current-year patches). The resolver
  # `.of_resolve_negative_pool_policy()`, the deprecated
  # `otsu_unburned_generation` alias, and the `deterministic_direct` operational
  # mode are gone. HONEST GUARD: a config that carries a stale
  # `negative_pool_policy` / `otsu_unburned_generation` option errors clearly
  # rather than being silently ignored. The historical canonical value
  # `"all_sources"` is still accepted as a no-op so existing scripts keep
  # working. `neg_pool` is then pinned to the fixed constant below.
  for (.stale_key in c("negative_pool_policy", "negative_pool_strategy",
                       "otsu_unburned_generation")) {
    if (!is.null(options[[.stale_key]]) &&
        !identical(options[[.stale_key]], "all_sources")) {
      stop(
        "`negative_pool_policy`/`otsu_unburned_generation` was removed in ",
        "2026-06; the supervised pipeline is always 'all_sources'. ",
        "Remove this option (you passed `", .stale_key, " = ",
        as.character(options[[.stale_key]])[1L], "`).",
        call. = FALSE
      )
    }
  }
  neg_pool <- "all_sources"

  output_routes <- .of_build_supervised_output_routes(
    output_dir = output_dir, target_year = target_year,
    run_name = run_name, scenario = run_label,
    flat = flat_output_routes
  )

  # ---- Gate 1B (2026-06-07): resolve cfg$model_params + cfg$train_control ----
  # The cfg becomes the SINGLE SOURCE OF TRUTH for all supervised
  # methodological / training-control parameters. Canonical defaults come from
  # .of_canonical_model_params() / .of_canonical_train_control() (the ONLY
  # places they are defined); user-supplied builder args override them and are
  # validated here. Every downstream stage reads cfg$model_params /
  # cfg$train_control, so the 5-layer default duplication is eliminated.
  # Capture whether model_params was explicitly supplied BEFORE it is folded
  # onto the canonical block (provenance "user" vs "default").
  .model_params_user_set <- !is.null(model_params)
  .model_params_requested <- model_params
  model_params <- .of_resolve_supervised_model_params(model_params)
  # Gate 1B (2026-06-07): PER-FIELD provenance for cfg$model_params. The builder
  # accepts a PARTIAL override list and folds it onto the canonical block (see
  # .of_resolve_supervised_model_params), so provenance MUST be recorded per xgb
  # field exactly like train_control: each field carries its canonical value, the
  # requested value (NA when the field was not in the partial override), the
  # resolved value, and provenance ("user" when that specific field was supplied,
  # "default" otherwise). scale_pos_weight is intentionally absent (site-specific,
  # rejected by the resolver, never a cfg$model_params field).
  mp_provenance <- .of_model_params_provenance(.model_params_requested)

  # ---- GATE 6.7 (2026-06-12): typed negative-pool + runtime resolution -------
  # negative_pool_params (PUBLIC methodological block) is the SINGLE source of
  # truth for the two caps and the user-settable negative-pool knobs. Its caps
  # flow into cfg$train_control$caps (the one place every stage reads caps from),
  # so there is NO second caps source. runtime_options (reuse_existing /
  # write_outputs / verbose) is technical and is EXCLUDED from the methodological
  # fingerprint. The remaining negative-pool knobs are fixed internal defaults
  # (.of_supervised_negative_pool_internal_defaults()), NOT public.
  np_res      <- .of_resolve_negative_pool_params(negative_pool_params)
  rt_res      <- .of_resolve_runtime_options(runtime_options)
  np_internal <- .of_supervised_negative_pool_internal_defaults()

  train_control <- .of_resolve_supervised_train_control(
    nrounds_max                = nrounds_max,
    early_stop                 = early_stop,
    oof_seed_base              = oof_seed_base,
    final_sampling_seed        = final_sampling_seed,
    final_seed                 = final_seed,
    random_seed                = random_seed,
    val_frac                   = val_frac,
    group_col                  = group_col,
    impute_numeric             = impute_numeric,
    impute_factor_missing      = impute_factor_missing,
    # GATE 6.7: caps come from the typed negative_pool_params block, the SINGLE
    # caps source. The train_control resolver no longer accepts cap_* args.
    caps                       = np_res$values$caps,
    feature_whitelist_override = feature_whitelist_override,
    feature_weights            = feature_weights,
    include_shape_features     = include_shape_features
  )

  # ---- Gate 1B / Precision 1 (2026-06-07): provenance of each methodological
  # train_control field. Tagged "user" when the caller explicitly supplied the
  # corresponding builder argument (non-NULL), "default" when the canonical
  # default was used. The public-boundary shim resolver
  # (.of_resolve_methodological_shim) consults this to (a) ERROR when a
  # function-level override conflicts with an EXPLICIT user value, and (b) emit a
  # deprecation warning when a function-level override is applied over a
  # canonical default. The provenance record is exposed on the cfg AND written
  # into the run-level summary the pipeline emits.
  tc_provenance <- list(
    nrounds_max           = if (is.null(nrounds_max))           "default" else "user",
    early_stop            = if (is.null(early_stop))            "default" else "user",
    oof_seed_base         = if (is.null(oof_seed_base))         "default" else "user",
    final_sampling_seed   = if (is.null(final_sampling_seed))   "default" else "user",
    final_seed            = if (is.null(final_seed))            "default" else "user",
    # GATE 6.7: random_seed provenance lives in the canonical-seeds provenance
    # alongside the three training seeds.
    random_seed           = if (is.null(random_seed))           "default" else "user",
    val_frac              = if (is.null(val_frac))              "default" else "user",
    group_col             = if (is.null(group_col))             "default" else "user",
    impute_numeric        = if (is.null(impute_numeric))        "default" else "user",
    impute_factor_missing = if (is.null(impute_factor_missing)) "default" else "user",
    # GATE 6.7: caps provenance now sourced from negative_pool_params$caps. The
    # per-bucket provenance is in resolved_params_provenance$negative_pool below;
    # a single block-level flag is kept here for the shim conflict-error logic.
    cap_random            = np_res$provenance$caps$random %||% "default",
    cap_otsu              = np_res$provenance$caps$otsu   %||% "default",
    feature_whitelist_override =
      if (is.null(feature_whitelist_override)) "default" else "user",
    feature_weights       = if (is.null(feature_weights))       "default" else "user",
    # OPTIONAL shape/size block (OtsuFire 0.12.0). "user" only when the
    # caller flipped it ON; the FALSE default is "default" so the OFF
    # baseline provenance is unchanged.
    include_shape_features = if (isTRUE(include_shape_features)) "user" else "default",
    model_params          = if (.model_params_user_set)         "user"    else "default"
  )

  # §N+27 RESOLVED (2026-06-05): the supervised burned-like registry was an
  # abandoned research line. The `burned_like_registry_path` parameter, its
  # default-path resolution, and the cfg field are removed. The deterministic
  # phase keeps its own registry READER (build_rbr_keep_pool_from_registry),
  # which receives an explicit `registry_path` from its caller; there is no
  # implicit default-path resolution helper anywhere in the package
  # (R/utils-registry.R and its .resolve_registry_path() were deleted in
  # 1E.6 as dead code with no live caller).

  cfg <- list(
    scenario                 = run_label,
    target_year              = target_year,
    inputs                   = inputs,
    run_name                 = run_name,
    output_dir               = output_dir,
    output_routes            = output_routes,
    # §N+26: informational constant only. No longer user-settable; always
    # "all_sources". Kept on the cfg so saveRDS(cfg) still records the policy
    # and downstream introspection / print() stays self-documenting.
    negative_pool_policy     = neg_pool,
    min_burned_pool_n        = min_burned_pool_n,
    # Gate 1B (2026-06-07): resolved methodological / training params. THE
    # single source of truth read by every supervised stage.
    model_params             = model_params,
    train_control            = train_control,
    # GATE 6.7 (2026-06-12): typed PUBLIC negative-pool block (random / otsu
    # knobs + caps). Caps are mirrored into train_control$caps (single source);
    # this block surfaces the full public negative-pool surface for print() /
    # manifests / introspection. The FIXED internal defaults (not user-settable)
    # live in $negative_pool_internal. The TECHNICAL runtime toggles live in
    # $negative_pool_runtime and are EXCLUDED from the methodological fingerprint.
    negative_pool_params     = np_res$values,
    negative_pool_runtime    = rt_res$values,
    negative_pool_internal   = np_internal,
    # Gate 1B / Precision 1: per-field provenance ("default" canonical vs "user"
    # explicitly set in the builder). The public-boundary shim resolver enriches
    # this into a full requested/resolved record at run time; the static builder
    # record below is the authoritative "was this field user-set?" source the
    # conflict-error logic relies on.
    #   - $train_control: a flat per-field "user"/"default" map (consumed by the
    #     .of_resolve_methodological_shim path at the public boundary). It also
    #     carries a single block-level $model_params flag for back-compat.
    #   - $model_params: a PER-FIELD record (canonical / requested / resolved /
    #     provenance per xgb field) because the builder folds a PARTIAL override
    #     onto the canonical block. This mirrors the train_control shim record and
    #     is what a downstream manifest renders as the model_params provenance
    #     table.
    resolved_params_provenance = list(train_control  = tc_provenance,
                                      model_params   = mp_provenance,
                                      # GATE 6.7: per-parameter provenance of the
                                      # typed negative_pool_params + runtime_options
                                      # (user vs default), mirroring the existing
                                      # train_control / model_params records.
                                      negative_pool  = np_res$provenance,
                                      runtime_options = rt_res$provenance),
    tool_paths               = tool_paths,
    options                  = options
  )
  class(cfg) <- c("otsufire_supervised_burned_config", "list")
  cfg
}

#' @export
print.otsufire_supervised_burned_config <- function(x, ...) {
  cat("<otsufire_supervised_burned_config>\n")
  cat("  scenario        :", x$scenario, "\n")
  cat("  target_year     :", x$target_year, "\n")
  cat("  run_name        :", x$run_name, "\n")
  cat("  output_dir      :", x$output_dir, "\n")
  cat("  negative_pool   :", x$negative_pool_policy, "\n")
  cat("  min_burned_pool :", x$min_burned_pool_n, "\n")
  if (!is.null(x$train_control)) {
    tc <- x$train_control
    cat("  train_control   : training=early-stopping selection + full-data refit",
        " nrounds=", tc$nrounds_max, " early_stop=", tc$early_stop,
        " val_frac=", tc$val_frac, "\n", sep = "")
    sd <- tc$seeds %||% list()
    cat("    seeds         : oof=", sd$oof_seed_base,
        " final_sampling=", sd$final_sampling_seed,
        " final=", sd$final_seed,
        " random=", sd$random_seed, "\n", sep = "")
  }
  # GATE 6.7: typed PUBLIC negative-pool block + caps + technical runtime.
  if (!is.null(x$negative_pool_params)) {
    np <- x$negative_pool_params
    cat("  negative_pool   :\n")
    cat("    random        : n_cells=", np$random$n_cells,
        " rbr_quantile=", np$random$rbr_quantile, "\n", sep = "")
    cat("    otsu          : candidate_threshold=", np$otsu$candidate_threshold,
        " reference_threshold=", np$otsu$reference_threshold, "\n", sep = "")
    cat("    caps          : random=", np$caps[["random"]],
        " otsu=", np$caps[["otsu"]], "\n", sep = "")
  }
  if (!is.null(x$negative_pool_runtime)) {
    rt <- x$negative_pool_runtime
    cat("  runtime_options : reuse_existing=", rt$reuse_existing,
        " write_outputs=", rt$write_outputs,
        " verbose=", rt$verbose, "\n", sep = "")
  }
  cat("  inputs          :\n")
  for (nm in names(x$inputs)) {
    v <- x$inputs[[nm]]
    cat(sprintf("    - %-22s: %s\n", nm,
                if (is.null(v)) "<NULL>" else
                  sprintf("%s [%s]",
                          if (!is.null(v$path)) v$path else "<in-memory>",
                          v$type)))
  }
  invisible(x)
}

# --- internal helpers --------------------------------------------------

#' Resolve + validate cfg$model_params (Gate 1B).
#'
#' Canonical default = .of_canonical_model_params() (the canonical xgb block
#' minus scale_pos_weight). A user override (named list) is merged on top of
#' the canonical block and validated. scale_pos_weight is rejected here: it is
#' site-specific (computed at train time), never a cfg field.
#'
#' @keywords internal
#' @noRd
.of_resolve_supervised_model_params <- function(model_params) {
  base <- .of_canonical_model_params()
  if (is.null(model_params)) return(base)
  if (!is.list(model_params) || is.null(names(model_params)) ||
      any(!nzchar(names(model_params)))) {
    stop("'model_params' must be NULL or a named list of XGBoost params.",
         call. = FALSE)
  }
  if ("scale_pos_weight" %in% names(model_params)) {
    stop("'model_params' must not set 'scale_pos_weight': it is computed at ",
         "train time from the training labels, not a cfg field.", call. = FALSE)
  }
  for (nm in names(model_params)) base[[nm]] <- model_params[[nm]]
  base
}

#' Per-field provenance for cfg$model_params (Gate 1B).
#'
#' The builder folds a PARTIAL `model_params` override onto the canonical xgb
#' block (`.of_canonical_model_params()`), so each xgb field's provenance is
#' tracked INDIVIDUALLY, exactly like `cfg$train_control`. For every canonical
#' field this returns a one-row-per-field record with:
#'   - `canonical`: the canonical default value (formatted one-line);
#'   - `requested`: the value the caller supplied for THIS field, or `NA` when
#'     the field was absent from the partial override;
#'   - `resolved`: the value actually stored in `cfg$model_params` (requested
#'     when supplied, else canonical);
#'   - `provenance`: `"user"` when this specific field was supplied, else
#'     `"default"`.
#' `scale_pos_weight` is intentionally NOT a field here: it is site-specific
#' (computed at train time) and is rejected by
#' `.of_resolve_supervised_model_params()`, so it is never a `cfg$model_params`
#' field. A user override containing only fields NOT in the canonical block would
#' already have failed downstream merge expectations; such extra fields are
#' surfaced here as additional rows with `canonical = NA` so the record stays a
#' faithful audit of what was requested.
#'
#' @param requested The raw `model_params` argument as passed to the builder
#'   (NULL = no override, else a named list of the partial fields), BEFORE it is
#'   folded onto the canonical block.
#' @return A named list keyed by xgb field; each element is a list with
#'   `canonical`, `requested`, `resolved`, `provenance`. The all-default case
#'   marks every field `"default"` with `requested = NA`.
#' @keywords internal
#' @noRd
.of_model_params_provenance <- function(requested) {
  canon <- .of_canonical_model_params()
  req   <- if (is.null(requested)) list() else requested
  # Union of canonical fields and any (atypical) extra requested fields, in a
  # stable order: canonical fields first (canonical block order), then extras.
  fields <- c(names(canon), setdiff(names(req), names(canon)))
  out <- list()
  for (nm in fields) {
    in_req      <- nm %in% names(req)
    canon_val   <- if (nm %in% names(canon)) canon[[nm]] else NULL
    resolved_val <- if (in_req) req[[nm]] else canon_val
    out[[nm]] <- list(
      canonical  = if (is.null(canon_val)) NA_character_ else .of_shim_fmt(canon_val),
      requested  = if (in_req) .of_shim_fmt(req[[nm]]) else NA_character_,
      resolved   = if (is.null(resolved_val)) NA_character_ else .of_shim_fmt(resolved_val),
      provenance = if (in_req) "user" else "default"
    )
  }
  out
}

#' Resolve + validate cfg$train_control (Gate 1B).
#'
#' Each argument is NULL (= canonical default from .of_canonical_train_control())
#' or a user override. Validates numeric ranges, match.arg-style enums, and
#' scalar types, then returns the resolved list stored as cfg$train_control.
#'
#' @keywords internal
#' @noRd
.of_resolve_supervised_train_control <- function(
    nrounds_max, early_stop, oof_seed_base, final_sampling_seed, final_seed,
    random_seed,
    val_frac, group_col, impute_numeric, impute_factor_missing,
    caps,
    feature_whitelist_override, feature_weights,
    include_shape_features = FALSE) {

  tc <- .of_canonical_train_control()

  # --- positive-integer / positive-numeric controls --------------------------
  chk_pos_int <- function(v, nm) {
    vi <- suppressWarnings(as.integer(v)[1L])
    if (is.na(vi) || vi <= 0L) {
      stop(sprintf("'%s' must be a single positive integer.", nm),
           call. = FALSE)
    }
    vi
  }
  if (!is.null(nrounds_max)) tc$nrounds_max <- chk_pos_int(nrounds_max, "nrounds_max")
  if (!is.null(early_stop))  tc$early_stop  <- chk_pos_int(early_stop, "early_stop")

  if (!is.null(oof_seed_base))
    tc$seeds$oof_seed_base <- chk_pos_int(oof_seed_base, "oof_seed_base")
  if (!is.null(final_sampling_seed))
    tc$seeds$final_sampling_seed <- chk_pos_int(final_sampling_seed,
                                                "final_sampling_seed")
  if (!is.null(final_seed))
    tc$seeds$final_seed <- chk_pos_int(final_seed, "final_seed")
  # GATE 6.7 (2026-06-12): random_seed (negative-pool random-background seed)
  # lives in the canonical seeds block alongside the three training seeds.
  if (!is.null(random_seed))
    tc$seeds$random_seed <- chk_pos_int(random_seed, "random_seed")

  # --- val_frac in (0, 1) ----------------------------------------------------
  if (!is.null(val_frac)) {
    if (!is.numeric(val_frac) || length(val_frac) != 1L || is.na(val_frac) ||
        val_frac <= 0 || val_frac >= 1) {
      stop("'val_frac' must be a single number in (0, 1).", call. = FALSE)
    }
    tc$val_frac <- as.numeric(val_frac)
  }

  # --- group_col: single non-empty string or NULL ----------------------------
  if (!is.null(group_col)) {
    if (!is.character(group_col) || length(group_col) != 1L ||
        !nzchar(group_col)) {
      stop("'group_col' must be a single non-empty character string or NULL.",
           call. = FALSE)
    }
    tc$group_col <- group_col
  }

  # --- imputation rules ------------------------------------------------------
  if (!is.null(impute_numeric)) {
    tc$impute_numeric <- match.arg(impute_numeric, c("median", "zero"))
  }
  if (!is.null(impute_factor_missing)) {
    if (!is.character(impute_factor_missing) ||
        length(impute_factor_missing) != 1L ||
        is.na(impute_factor_missing) || !nzchar(impute_factor_missing)) {
      stop("'impute_factor_missing' must be a single non-empty character ",
           "string.", call. = FALSE)
    }
    tc$impute_factor_missing <- impute_factor_missing
  }

  # --- negative-bucket caps (GATE 6.7): the resolved caps are the SINGLE source
  # of truth, supplied by .of_resolve_negative_pool_params() from the typed
  # negative_pool_params$caps block. The train_control resolver no longer parses
  # cap_* args; it just stores the already-validated named numeric here so every
  # downstream stage reads tc$caps.
  if (!is.null(caps)) {
    if (!is.numeric(caps) || is.null(names(caps)) ||
        !setequal(names(caps), c("random", "otsu"))) {
      stop("internal: train_control 'caps' must be a named numeric over exactly ",
           "{random, otsu}.", call. = FALSE)
    }
    tc$caps$random <- as.numeric(caps[["random"]])
    tc$caps$otsu   <- as.numeric(caps[["otsu"]])
  }

  # --- feature whitelist / weights (pass-through, light validation) ----------
  if (!is.null(feature_whitelist_override)) {
    if (!is.character(feature_whitelist_override) ||
        length(feature_whitelist_override) < 1L) {
      stop("'feature_whitelist_override' must be NULL or a non-empty ",
           "character vector.", call. = FALSE)
    }
    tc$feature_whitelist_override <- feature_whitelist_override
  }
  if (!is.null(feature_weights)) {
    if (!is.numeric(feature_weights) || is.null(names(feature_weights)) ||
        any(!nzchar(names(feature_weights)))) {
      stop("'feature_weights' must be NULL or a NAMED numeric vector.",
           call. = FALSE)
    }
    tc$feature_weights <- feature_weights
  }

  # --- OPTIONAL shape/size block (OtsuFire 0.12.0, OFF by default) ------------
  # Single boolean stored on the resolved train_control. THE single source of
  # truth read by extract_supervised_features() (compute the columns) and by
  # the OOF + FINAL engines (make them eligible model features via
  # .supervised_feature_universe()). FALSE -> baseline unchanged.
  if (!is.logical(include_shape_features) ||
      length(include_shape_features) != 1L || is.na(include_shape_features)) {
    stop("'include_shape_features' must be a single TRUE or FALSE.",
         call. = FALSE)
  }
  tc$include_shape_features <- isTRUE(include_shape_features)

  # OtsuFire OOF always uses the capped negative-sampling policy (the SAME one
  # the FINAL model uses), applied independently within each training fold.
  # tc$oof_sampling is a FIXED traceability constant ("capped") from
  # .of_canonical_train_control(); it is not user-settable.

  tc
}

# =============================================================================
# GATE 6.7 (2026-06-12): typed negative-pool config + runtime options.
#
# Mirrors the DETERMINISTIC builder's restrained surface (typed list-blocks,
# .of_check_unknown_keys whitelist with unknown-keys-error, visible defaults,
# per-parameter provenance). The PUBLIC negative-pool surface is intentionally
# small: only the methodological knobs that change WHICH negatives are selected
# (random n_cells / rbr_quantile, otsu candidate/reference thresholds) plus the
# two caps. Everything else is a FIXED internal default (NOT public).
# =============================================================================

#' Visible defaults for the PUBLIC typed negative_pool_params block (GATE 6.7).
#'
#' @return Named list with `random` (n_cells, rbr_quantile), `otsu`
#'   (candidate_threshold, reference_threshold) and `caps` (a named numeric over
#'   \{random, otsu\}).
#' @keywords internal
#' @noRd
.of_negative_pool_params_defaults <- function() {
  list(
    random = list(
      n_cells      = 1500L,
      rbr_quantile = 0.50
    ),
    otsu = list(
      candidate_threshold = 0,
      reference_threshold = 100
    ),
    caps = c(random = 1.0, otsu = 1.0),
    # PHASE 2 (artifact_hard hard-negative mining): ADDITIVE + OFF BY DEFAULT.
    # With enabled = FALSE this sub-block selects NO rows and the resolver MUST
    # NOT alter caps / enum / anything observable -- the resolved value stays
    # byte-identical to the pre-Phase-2 default (asserted by a baseline test).
    artifact_hard = list(
      enabled            = FALSE,           # master switch (OFF by default)
      # total_weight_ratio: the TOTAL effective weight of the artifact_hard pool
      # divided by the TOTAL weight of the burned (positive) pool, i.e.
      # sum(w[artifact_hard]) == total_weight_ratio * sum(w[burned]). This is the
      # POOL-LEVEL balance knob, NOT the per-polygon weight: each artifact_hard
      # row gets total_weight_ratio * sum_w_pos / n_artifact_hard. Configurable;
      # the 0.10 default is a conservative starting point, NOT a universally
      # recommended value (the weight x year sweep showed the right weight is
      # year-dependent). An explicit `source_weights` pin overrides this.
      total_weight_ratio = 0.10,
      persist_ratio_max  = 0.35,
      rbr_med_reference  = "negative",      # "negative" (existing neg pool) | "positive"
      rbr_med_min_q      = 0.90,            # quantile of the reference pool's rbr_med (runtime)
      reason_whitelist   = character(0),    # explicit artifact reason_1 values (user fills)
      persist_delta_max  = -100,
      area_ha_min        = 500,
      doy_iqr_max        = 1
    ),
    # PHASE 2: optional per-row source weighting. NULL (default) => NO weight
    # vector is ever set => byte-identical to today.
    source_weights = NULL
  )
}

#' Visible defaults for the technical runtime_options block (GATE 6.7).
#'
#' These are TECHNICAL toggles only; they MUST NOT enter the methodological
#' fingerprint.
#' @keywords internal
#' @noRd
.of_runtime_options_defaults <- function() {
  list(
    reuse_existing = TRUE,
    write_outputs  = TRUE,
    verbose        = TRUE
  )
}

#' FIXED internal negative-pool defaults (GATE 6.7) -- NOT public, NOT settable.
#'
#' The lower-level negative-pool / Otsu-residual engine knobs that are pinned to
#' validated operational constants (mirroring the Otsu-guided engine's fixed
#' settings). They MAY enter the methodological fingerprint (they affect the
#' pool) but are NOT user-settable and NOT part of the public block. The clean
#' public names in `negative_pool_params` map onto the historical
#' `unb_*` / `otsu_negative_*` bindings; these are the rest.
#' @keywords internal
#' @noRd
.of_supervised_negative_pool_internal_defaults <- function() {
  list(
    # random background
    random_exclusion_buffer_m = 500,
    random_patch_size_cells   = 3,
    # otsu residual negative
    otsu_mode                 = "burnable_only",
    otsu_min_threshold_value  = 0,
    otsu_min_pixels           = 8,
    otsu_buffers_m            = 90,
    otsu_core_thr             = 0.60,
    otsu_alpha_boost          = 0.25,
    otsu_min_base_boost       = 0.35,
    otsu_dist_power           = 1,
    otsu_keep_hi              = 0.45,
    otsu_drop_lo              = 0.15,
    otsu_excl_buffer_m        = 0,
    otsu_min_area_ha          = 0,
    otsu_use_drop             = TRUE,
    otsu_drop_max_s_patch     = 0.15,
    allow_empty_otsu_pool     = FALSE
  )
}

#' Resolve + validate the typed PUBLIC negative_pool_params block (GATE 6.7).
#'
#' Whitelisted per sub-block (unknown key anywhere -> ERROR). Validates types
#' (n_cells integer > 0; rbr_quantile in [0,1]; thresholds numeric; caps a named
#' numeric over exactly {random, otsu}, each >= 0). Returns the resolved values
#' AND a per-parameter provenance record ("user" vs "default").
#'
#' @param negative_pool_params The user `negative_pool_params` list (may be
#'   empty / NULL).
#' @return list(values, provenance).
#' @keywords internal
#' @noRd
.of_resolve_negative_pool_params <- function(negative_pool_params) {
  defs <- .of_negative_pool_params_defaults()
  np   <- if (is.null(negative_pool_params)) list() else negative_pool_params
  .of_check_named_list(np, "negative_pool_params")
  .of_check_unknown_keys(np, names(defs), "negative_pool_params")

  # --- random sub-block ------------------------------------------------------
  rnd <- if (is.null(np$random)) list() else np$random
  .of_check_named_list(rnd, "negative_pool_params$random")
  .of_check_unknown_keys(rnd, names(defs$random), "negative_pool_params$random")
  n_cells <- defs$random$n_cells
  if (!is.null(rnd$n_cells)) {
    ni <- suppressWarnings(as.integer(rnd$n_cells)[1L])
    if (is.na(ni) || ni <= 0L) {
      stop("'negative_pool_params$random$n_cells' must be a single positive ",
           "integer.", call. = FALSE)
    }
    n_cells <- ni
  }
  rbr_q <- defs$random$rbr_quantile
  if (!is.null(rnd$rbr_quantile)) {
    q <- rnd$rbr_quantile
    if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 0 || q > 1) {
      stop("'negative_pool_params$random$rbr_quantile' must be a single number ",
           "in [0, 1].", call. = FALSE)
    }
    rbr_q <- as.numeric(q)
  }

  # --- otsu sub-block --------------------------------------------------------
  ots <- if (is.null(np$otsu)) list() else np$otsu
  .of_check_named_list(ots, "negative_pool_params$otsu")
  .of_check_unknown_keys(ots, names(defs$otsu), "negative_pool_params$otsu")
  chk_num1 <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf("'%s' must be a single finite numeric value.", nm),
           call. = FALSE)
    }
    as.numeric(v)
  }
  cand_thr <- defs$otsu$candidate_threshold
  if (!is.null(ots$candidate_threshold)) {
    cand_thr <- chk_num1(ots$candidate_threshold,
                         "negative_pool_params$otsu$candidate_threshold")
  }
  ref_thr <- defs$otsu$reference_threshold
  if (!is.null(ots$reference_threshold)) {
    ref_thr <- chk_num1(ots$reference_threshold,
                        "negative_pool_params$otsu$reference_threshold")
  }

  # --- caps (single source of truth) -----------------------------------------
  caps <- defs$caps
  caps_user <- list(random = FALSE, otsu = FALSE)
  if (!is.null(np$caps)) {
    uc <- np$caps
    if (!is.numeric(uc) || is.null(names(uc)) ||
        !setequal(names(uc), c("random", "otsu"))) {
      stop("'negative_pool_params$caps' must be a NAMED numeric over exactly ",
           "{random, otsu}.", call. = FALSE)
    }
    chk_cap <- function(v, nm) {
      if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
        stop(sprintf("'%s' must be a single number >= 0 (Inf disables the cap).",
                     nm), call. = FALSE)
      }
      as.numeric(v)
    }
    caps[["random"]] <- chk_cap(uc[["random"]], "negative_pool_params$caps['random']")
    caps[["otsu"]]   <- chk_cap(uc[["otsu"]],   "negative_pool_params$caps['otsu']")
    caps_user$random <- TRUE
    caps_user$otsu   <- TRUE
  }

  # --- PHASE 2: artifact_hard sub-block (ADDITIVE, OFF BY DEFAULT) ------------
  # When the user supplies nothing this resolves to the OFF default block and
  # selects NO rows. Unknown keys ERROR (mirrors every other sub-block). The
  # resolved scalar types are validated so a malformed override fails fast, but
  # NONE of this alters the caps / enum / training rows when enabled = FALSE.
  ah_def <- defs$artifact_hard
  ah_in  <- if (is.null(np$artifact_hard)) list() else np$artifact_hard
  .of_check_named_list(ah_in, "negative_pool_params$artifact_hard")
  .of_check_unknown_keys(ah_in, names(ah_def), "negative_pool_params$artifact_hard")
  chk_lgl1 <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be a single TRUE/FALSE.", nm), call. = FALSE)
    }
    as.logical(v)
  }
  chk_num1_ah <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf("'%s' must be a single finite numeric value.", nm),
           call. = FALSE)
    }
    as.numeric(v)
  }
  ah_enabled <- ah_def$enabled
  if (!is.null(ah_in$enabled)) ah_enabled <- chk_lgl1(ah_in$enabled,
    "negative_pool_params$artifact_hard$enabled")
  # total_weight_ratio: pool-level balance (artifact_hard total weight / burned
  # total weight). Single finite numeric >= 0. Not used when enabled = FALSE.
  ah_total_weight_ratio <- ah_def$total_weight_ratio
  if (!is.null(ah_in$total_weight_ratio)) {
    twr <- ah_in$total_weight_ratio
    if (!is.numeric(twr) || length(twr) != 1L || !is.finite(twr) || twr < 0) {
      stop("'negative_pool_params$artifact_hard$total_weight_ratio' must be a ",
           "single finite numeric >= 0 (artifact_hard total weight / burned ",
           "total weight).", call. = FALSE)
    }
    ah_total_weight_ratio <- as.numeric(twr)
  }
  ah_persist_ratio_max <- ah_def$persist_ratio_max
  if (!is.null(ah_in$persist_ratio_max)) ah_persist_ratio_max <- chk_num1_ah(
    ah_in$persist_ratio_max, "negative_pool_params$artifact_hard$persist_ratio_max")
  ah_rbr_med_min_q <- ah_def$rbr_med_min_q
  if (!is.null(ah_in$rbr_med_min_q)) {
    q <- ah_in$rbr_med_min_q
    if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 0 || q > 1) {
      stop("'negative_pool_params$artifact_hard$rbr_med_min_q' must be a single ",
           "number in [0, 1].", call. = FALSE)
    }
    ah_rbr_med_min_q <- as.numeric(q)
  }
  ah_rbr_med_reference <- ah_def$rbr_med_reference
  if (!is.null(ah_in$rbr_med_reference)) {
    rr <- ah_in$rbr_med_reference
    if (!is.character(rr) || length(rr) != 1L || !(rr %in% c("negative", "positive"))) {
      stop("'negative_pool_params$artifact_hard$rbr_med_reference' must be ",
           "\"negative\" or \"positive\".", call. = FALSE)
    }
    ah_rbr_med_reference <- rr
  }
  ah_reason_whitelist <- ah_def$reason_whitelist
  if (!is.null(ah_in$reason_whitelist)) {
    rw <- ah_in$reason_whitelist
    if (!is.character(rw)) {
      stop("'negative_pool_params$artifact_hard$reason_whitelist' must be a ",
           "character vector (possibly empty).", call. = FALSE)
    }
    ah_reason_whitelist <- as.character(rw)
  }
  ah_persist_delta_max <- ah_def$persist_delta_max
  if (!is.null(ah_in$persist_delta_max)) ah_persist_delta_max <- chk_num1_ah(
    ah_in$persist_delta_max, "negative_pool_params$artifact_hard$persist_delta_max")
  ah_area_ha_min <- ah_def$area_ha_min
  if (!is.null(ah_in$area_ha_min)) ah_area_ha_min <- chk_num1_ah(
    ah_in$area_ha_min, "negative_pool_params$artifact_hard$area_ha_min")
  ah_doy_iqr_max <- ah_def$doy_iqr_max
  if (!is.null(ah_in$doy_iqr_max)) ah_doy_iqr_max <- chk_num1_ah(
    ah_in$doy_iqr_max, "negative_pool_params$artifact_hard$doy_iqr_max")
  artifact_hard <- list(
    enabled            = ah_enabled,
    total_weight_ratio = ah_total_weight_ratio,
    persist_ratio_max  = ah_persist_ratio_max,
    rbr_med_reference  = ah_rbr_med_reference,
    rbr_med_min_q      = ah_rbr_med_min_q,
    reason_whitelist   = ah_reason_whitelist,
    persist_delta_max  = ah_persist_delta_max,
    area_ha_min        = ah_area_ha_min,
    doy_iqr_max        = ah_doy_iqr_max
  )

  # --- PHASE 2: source_weights (optional per-row weighting; NULL = OFF) -------
  source_weights <- NULL
  if (!is.null(np$source_weights)) {
    sw <- np$source_weights
    if (!is.numeric(sw) || is.null(names(sw)) || any(!nzchar(names(sw)))) {
      stop("'negative_pool_params$source_weights' must be a NAMED numeric ",
           "vector (source -> weight) or NULL.", call. = FALSE)
    }
    if (any(is.na(sw)) || any(sw < 0)) {
      stop("'negative_pool_params$source_weights' values must be finite and ",
           ">= 0.", call. = FALSE)
    }
    source_weights <- sw
  }

  values <- list(
    random = list(n_cells = n_cells, rbr_quantile = rbr_q),
    otsu   = list(candidate_threshold = cand_thr, reference_threshold = ref_thr),
    caps   = caps,
    # PHASE 2 additive keys. With the defaults (enabled = FALSE, weights NULL)
    # these are inert: no row is selected, no weight vector is built.
    artifact_hard  = artifact_hard,
    source_weights = source_weights
  )
  provenance <- list(
    random = list(
      n_cells      = if (is.null(rnd$n_cells))      "default" else "user",
      rbr_quantile = if (is.null(rnd$rbr_quantile)) "default" else "user"
    ),
    otsu = list(
      candidate_threshold = if (is.null(ots$candidate_threshold)) "default" else "user",
      reference_threshold = if (is.null(ots$reference_threshold)) "default" else "user"
    ),
    caps = list(
      random = if (caps_user$random) "user" else "default",
      otsu   = if (caps_user$otsu)   "user" else "default"
    ),
    artifact_hard = list(
      enabled            = if (is.null(ah_in$enabled))            "default" else "user",
      total_weight_ratio = if (is.null(ah_in$total_weight_ratio)) "default" else "user",
      persist_ratio_max = if (is.null(ah_in$persist_ratio_max)) "default" else "user",
      rbr_med_reference = if (is.null(ah_in$rbr_med_reference)) "default" else "user",
      rbr_med_min_q     = if (is.null(ah_in$rbr_med_min_q))     "default" else "user",
      reason_whitelist  = if (is.null(ah_in$reason_whitelist))  "default" else "user",
      persist_delta_max = if (is.null(ah_in$persist_delta_max)) "default" else "user",
      area_ha_min       = if (is.null(ah_in$area_ha_min))       "default" else "user",
      doy_iqr_max       = if (is.null(ah_in$doy_iqr_max))       "default" else "user"
    ),
    source_weights = if (is.null(np$source_weights)) "default" else "user"
  )
  list(values = values, provenance = provenance)
}

#' Resolve + validate the technical runtime_options block (GATE 6.7).
#'
#' Whitelisted (reuse_existing / write_outputs / verbose only; unknown key ->
#' ERROR). Each must be a single TRUE/FALSE. Runtime changes MUST NOT alter the
#' methodological fingerprint, so these are kept in their own cfg block and
#' excluded from every fingerprint.
#' @keywords internal
#' @noRd
.of_resolve_runtime_options <- function(runtime_options) {
  defs <- .of_runtime_options_defaults()
  ro   <- if (is.null(runtime_options)) list() else runtime_options
  .of_check_named_list(ro, "runtime_options")
  .of_check_unknown_keys(ro, names(defs), "runtime_options")
  chk_flag <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'runtime_options$%s' must be a single TRUE/FALSE.", nm),
           call. = FALSE)
    }
    v
  }
  values <- defs
  prov   <- list()
  for (nm in names(defs)) {
    if (!is.null(ro[[nm]])) {
      values[[nm]] <- chk_flag(ro[[nm]], nm)
      prov[[nm]]   <- "user"
    } else {
      prov[[nm]]   <- "default"
    }
  }
  list(values = values, provenance = prov)
}

#' Centralized supervised-RUN input-path conventions (§N+25)
#'
#' Single source of truth for the by-convention default paths of the
#' supervised-RUN inputs. `build_supervised_burned_config()` calls this to
#' fill the cfg defaults; the orchestrator and the unburned helpers fall back
#' to the SAME constructions when the cfg does not carry an explicit path, so
#' the cfg-default and the runtime-fallback always agree (byte-identical).
#'
#' Returns a named list whose entries are character paths, or NULL when
#' `data_base` (resp. `composite_base`) is not available to build them. The
#' CORINE epoch is resolved with `.of_corine_year_for()` (the same mapping the
#' orchestrator's `get_corine_year()` uses).
#'
#' @keywords internal
#' @noRd
.of_supervised_convention_paths <- function(data_base, composite_base,
                                            result_name, target_year) {
  has_db <- is.character(data_base) && length(data_base) == 1L &&
            nzchar(data_base)
  has_cb <- is.character(composite_base) && length(composite_base) == 1L &&
            nzchar(composite_base)
  corine_year <- .of_corine_year_for(target_year)

  list(
    # AW (delayed) change-index mosaic lives under composite_base/Autumn.
    delayed_change_index = if (has_cb) {
      file.path(composite_base, "Autumn",
                paste0("mean_mean_", target_year, "_mosaic.tif"))
    } else NULL,
    peninsula_shapefile = if (has_db) {
      file.path(data_base, "Borders", "Iberian_peninsula.shp")
    } else NULL,
    topo = if (has_db) {
      file.path(data_base, "Topography", "elevation_slope.tif")
    } else NULL,
    corine_raster = if (has_db) {
      file.path(data_base, "Corine_Masks",
                paste0("CLC_", corine_year, "_peninsula.tif"))
    } else NULL,
    burnable_mask = if (has_db) {
      file.path(data_base, "Corine_Masks",
                paste0("burneable_mask_binary_corine_", corine_year,
                       "_ETRS89.tif"))
    } else NULL
  )
}

#' Resolve a cfg$inputs on-disk path (Gate 1B input inventory).
#'
#' Single shared accessor for the supervised-input subsystem: returns the
#' normalized on-disk PATH of `config$inputs[[name]]` when that input is a
#' `type = "path"` spec, else `NULL` (the input is absent, in-memory, or the
#' config itself is NULL). Every probabilistic refinement stage that needs an input's file
#' path (orchestrator, pool builder, feature extractor) consumes the input
#' THROUGH this accessor, so `cfg$inputs` is the single source of truth and no
#' stage reconstructs an input path by filename convention behind the cfg's
#' back. A `NULL` return lets the caller apply its documented optional-input
#' behaviour (skip / convention-fallback for the standalone-script path).
#'
#' @param config probabilistic refinement config S3 object, or `NULL`.
#' @param name character scalar input name (a key of `cfg$inputs`).
#' @return character path or `NULL`.
#'
#' @keywords internal
#' @noRd
.of_sup_input_path <- function(config, name) {
  if (is.null(config) || is.null(config$inputs)) return(NULL)
  sp <- config$inputs[[name]]
  if (is.null(sp)) return(NULL)
  if (!is.null(sp$path) && length(sp$path) == 1L && !is.na(sp$path) &&
      nzchar(sp$path) && identical(sp$type, "path")) {
    sp$path
  } else {
    NULL
  }
}

#' CORINE epoch for a target year (§N+25 convention helper).
#'
#' Mirrors the orchestrator's `get_corine_year()` and the unburned helpers'
#' `get_corine_year_unb*()` so every site resolves the same CORINE raster /
#' burnable-mask filename. Returns a character year string.
#'
#' @keywords internal
#' @noRd
.of_corine_year_for <- function(y) {
  y <- suppressWarnings(as.integer(y)[1L])
  if (!is.na(y) && y >= 1984L && y <= 1999L) "1990"
  else if (!is.na(y) && y <= 2005L) "2000"
  else if (!is.na(y) && y <= 2011L) "2006"
  else if (!is.na(y) && y <= 2017L) "2012"
  else "2018"
}

# BUG 3 Phase 1a (2026-06-05): removed the dead supervised getwd()-walk
# finders `.of_find_supervised_engine_root()` and
# `.of_find_supervised_scripts_root()`. They only populated inert
# `supervised_engine_root` / `scripts_root` config metadata that fed the
# now-removed external PROJECT_PATHS.R fallback. The supervised engine is
# in-package; nothing reads these roots anymore.

# §N+26 (2026-06-05): `.of_resolve_negative_pool_policy()` REMOVED. The
# supervised pipeline no longer resolves a negative-pool policy — it is always
# `all_sources`. The user-facing knob, the `otsu_unburned_generation` alias and
# the `deterministic_direct` mode are gone; a stale option now triggers the
# honest guard in `build_supervised_burned_config()` instead of a silent
# redirect. (This is unrelated to the D4a empty-Otsu-pool robustness fallback
# in internal-sup-otsu-negative.R, which is retained.)

#' @keywords internal
#' @noRd
.of_build_supervised_output_routes <- function(output_dir, target_year,
                                                run_name, scenario,
                                                flat = FALSE) {
  # `flat = TRUE` hangs every product DIRECTLY off output_dir (output_dir/01_POOLS,
  # output_dir/03_FEATURES, ...) instead of the canonical deep
  # <output_dir>/<year>/<run_name>/SUPERVISED/<run_label> tree. Use it when
  # output_dir already IS the per-run scenario folder, so the layout is not
  # duplicated. Default FALSE = the canonical nested layout (unchanged).
  base <- if (isTRUE(flat)) {
    output_dir
  } else {
    file.path(output_dir, as.character(target_year), run_name,
              "SUPERVISED", scenario)
  }
  prefix_oof <- sprintf("%d_%s_patch", target_year, scenario)
  prefix     <- sprintf("%d_%s_patch_certified", target_year, scenario)
  # Bug 1 (0.3.0): folder names match the operational runtime
  # (`07_FINAL_MODEL_V2`, `08_SCORED`, `09_FINAL_MAP`,
  # `11_CONSISTENCY_CHECKS`). Previously the config advertised
  # `06_FINAL_MODEL`/`07_SCORED`/`08_FINAL_MAP`/`09_CONSISTENCY`,
  # which did not exist on disk — public callers reading
  # `config$output_routes` directly were getting nonexistent paths.
  list(
    base                = base,
    pools_dir           = file.path(base, "01_POOLS"),
    folds_dir           = file.path(base, "02_FOLDS"),
    features_dir        = file.path(base, "03_FEATURES"),
    matrix_dir          = file.path(base, "04_MATRIX"),
    oof_dir             = file.path(base, "05_OOF"),
    final_model_dir     = file.path(base, "07_FINAL_MODEL_V2"),
    scored_dir          = file.path(base, "08_SCORED"),
    final_map_dir       = file.path(base, "09_FINAL_MAP"),
    consistency_dir     = file.path(base, "11_CONSISTENCY_CHECKS"),
    logs_dir            = file.path(base, "99_LOGS"),
    pools_gpkg          = file.path(base, "01_POOLS",
                                     sprintf("%d_%s_pools.gpkg",
                                             target_year, scenario)),
    train_with_folds_gpkg = file.path(base, "02_FOLDS",
                                       sprintf("%d_train_with_folds_5000m.gpkg",
                                               target_year)),
    features_geometry_gpkg = file.path(base, "03_FEATURES",
                                        "features_geometry.gpkg"),
    oof_agg_csv         = file.path(base, "05_OOF",
                                     paste0(prefix_oof, "_oof_agg.csv")),
    final_model_rds     = file.path(base, "07_FINAL_MODEL_V2",
                                     paste0(prefix, "_final_model.rds")),
    final_map_gpkg      = file.path(base, "09_FINAL_MAP",
                                     paste0(prefix, "_final_map.gpkg")),
    burned_like_gpkg    = file.path(base, "09_FINAL_MAP",
                                     paste0(prefix, "_burned_like_scored.gpkg")),
    timing_csv          = file.path(base, "99_LOGS",
                                     paste0(target_year, "_timing_steps.csv"))
  )
}
