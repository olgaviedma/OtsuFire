#' @title Build a configuration object for the supervised burned-area workflow
#'
#' @description
#' Create, validate, and normalize the configuration used by the supervised
#' OtsuFire workflow.
#'
#' Think of this as the setup step for the supervised (machine-learning)
#' pipeline. It gathers:
#' \itemize{
#'   \item the spatial and labelled inputs (the deterministic decision layer,
#'     change-index rasters, topography, land cover, optional hotspots and an
#'     external reference),
#'   \item the negative-pool and training/sampling controls,
#'   \item the feature configuration (whitelist and per-feature weights),
#'   \item the XGBoost hyper-parameters,
#'   \item the output folder structure,
#'   \item and the technical runtime options
#' }
#' into one validated object of class
#' `otsufire_supervised_burned_config` that you then pass to the supervised
#' stages:
#' \code{\link[=build_supervised_training_pools]{build_supervised_training_pools()}},
#' \code{\link[=make_spatial_folds]{make_spatial_folds()}},
#' \code{\link[=extract_supervised_features]{extract_supervised_features()}},
#' \code{\link[=run_oof_diagnostics]{run_oof_diagnostics()}},
#' \code{\link[=train_final_burned_model]{train_final_burned_model()}},
#' \code{\link[=score_supervised_burned_map]{score_supervised_burned_map()}},
#' or the single-call orchestrator
#' \code{\link[=run_oneyear_supervised_pipeline]{run_oneyear_supervised_pipeline()}}.
#'
#' It does not build pools, train models or score polygons by itself. Its job
#' is to prepare the workflow cleanly and fail early when required inputs or
#' settings are inconsistent.
#'
#' The arguments are documented below in four conceptual groups, even though
#' the function signature exposes each argument individually:
#' \enumerate{
#'   \item spatial and labelled inputs,
#'   \item training and sampling controls,
#'   \item feature and model controls,
#'   \item runtime and advanced options.
#' }
#'
#' \strong{Optional inputs and convention paths}
#'
#' Several raster and vector inputs are optional. When such an input is left
#' `NULL` and the advanced option `data_base` (and, for the delayed change
#' index, `composite_base`) is supplied, the builder fills in a default path
#' by naming convention so that a complete run can be configured with a small
#' set of arguments. When `data_base` is not available (for example an
#' ahead-of-batch configuration or a unit test), the optional input simply
#' stays `NULL` and is resolved later by the stage that needs it. See the
#' per-argument notes for the exact `NULL` behaviour.
#'
#' @param scenario Character scalar. Methodological preset that names the run
#'   and its output sub-folders. One of `"balanced"` (default), `"original"`,
#'   `"lax"`, `"restrictive"`. Consumed throughout the workflow for naming and
#'   routing.
#'
#' @param internal_decisions sf POLYGON layer or a GPKG path. The deterministic
#'   decision layer (the labelled output of the deterministic stage) produced
#'   for the same year and scenario. \strong{Required.} Provides the positive
#'   (keep-class) labels and seeds the labelled training pool. Consumed by
#'   \code{build_supervised_training_pools()}.
#'
#' @param change_index SpatRaster or raster path. The main annual change-index
#'   raster (band 1 = summer RBR, band 2 = post-fire day-of-year).
#'   \strong{Required.} Consumed end-to-end from `cfg$inputs$change_index` by
#'   the pool builder (legacy-Otsu severity raster) and the feature extractor
#'   (the same-window RBR features). `immediate_change_index_path` is the
#'   canonical alias for this argument (see below).
#'
#' @param delayed_change_index SpatRaster or raster path. Optional delayed
#'   (autumn-winter) change-index raster used for persistence-style features.
#'   When supplied it is consumed as the `rbr_aw` feature layer by
#'   \code{extract_supervised_features()}. When `NULL` and `composite_base` is
#'   available it defaults to the convention path
#'   `<composite_base>/Autumn/mean_mean_<target_year>_mosaic.tif`; otherwise it
#'   stays `NULL` and the all-window RBR features are skipped.
#'   `delayed_change_index_path` is the canonical alias for this argument.
#'
#' @param immediate_change_index_path,delayed_change_index_path Canonical names
#'   for the two change-index inputs. `immediate_change_index_path` is the
#'   canonical name of the required immediate change index (deprecated-friendly
#'   alias of `change_index`); `delayed_change_index_path` is the canonical name
#'   of the optional delayed change index (alias of `delayed_change_index`).
#'   Supply either the canonical name or its alias for a given input; supplying
#'   both is allowed only if they are identical, otherwise the call errors with
#'   a clear conflict message. Both resolve onto the single
#'   `cfg$inputs$change_index` / `cfg$inputs$delayed_change_index` field.
#'   Deprecated aliases retained for backward compatibility. Prefer
#'   `change_index` / `delayed_change_index`.
#'
#' @param hotspots sf POINT layer or path. Optional hotspot (active-fire) layer
#'   for the target year. When present it powers the hotspot feature block in
#'   \code{extract_supervised_features()}; when `NULL` the hotspot features are
#'   not built (see the no-hotspot profile in \strong{Details}). Optional;
#'   stays `NULL` when not supplied.
#'
#' @param reference_burned_map sf POLYGON / raster / path. Optional external
#'   burned-area reference (for example an EFFIS perimeter layer) used by
#'   \code{\link[=validate_fire_maps]{validate_fire_maps()}} on the thresholded
#'   supervised output. Not used during the run itself; optional; stays `NULL`
#'   when not supplied.
#'
#' @param peninsula_shapefile sf / SpatVector / path. Optional study-area border
#'   polygon (e.g. the Iberian peninsula). Consumed by the negative-pool
#'   (legacy-Otsu) builder only when `options$legacy_otsu_mode = "corine"`, to
#'   crop CORINE to the study area. When `NULL` and `data_base` is available it
#'   defaults to `<data_base>/Borders/Iberian_peninsula.shp`; otherwise it stays
#'   `NULL`.
#'
#' @param topo SpatRaster / path. Optional two-band topography raster
#'   (band 1 = DEM, band 2 = slope). Consumed by
#'   \code{extract_supervised_features()} for the elevation and slope feature
#'   blocks. When `NULL` and `data_base` is available it defaults to
#'   `<data_base>/Topography/elevation_slope.tif`; otherwise it stays `NULL`.
#'   The DEM and slope live in a single file by convention, so they are exposed
#'   as one `topo` input rather than two separate paths.
#'
#' @param corine_raster SpatRaster / path. Optional CORINE land-cover raster for
#'   the target year's CORINE epoch. Consumed by the feature extractor for the
#'   land-cover feature block. When `NULL` and `data_base` is available it
#'   defaults to `<data_base>/Corine_Masks/CLC_<corine_year>_peninsula.tif`
#'   (the epoch is chosen automatically from `target_year`); otherwise it stays
#'   `NULL`.
#'
#' @param burnable_mask SpatRaster / path. Optional binary burnable-area mask
#'   for the target year's CORINE epoch. Consumed by the negative-pool builders
#'   (random-background and Otsu buckets). When `NULL` and `data_base` is
#'   available it defaults to
#'   `<data_base>/Corine_Masks/burneable_mask_binary_corine_<corine_year>_ETRS89.tif`;
#'   otherwise it stays `NULL`.
#'
#' @param target_year Integer scalar in `[1900, 2100]`. The target year, used
#'   to name and route outputs, to filter temporal layers, and to resolve the
#'   CORINE epoch of the convention paths. \strong{Required.}
#'
#' @param output_dir Character scalar. Root output directory. The builder
#'   creates a `<output_dir>/<year>/<run_name>/SUPERVISED/<scenario>/...`
#'   layout. Defaults to `tempdir()`.
#'
#' @param run_name Character scalar. Stable identifier for the supervised run,
#'   used to name output folders and files. Single non-empty string. Defaults
#'   to `"supervised_burned_map"`.
#'
#' @param min_burned_pool_n Integer `>= 0`. Sparse-year guard: the pipeline
#'   aborts before training when fewer than this many keep-class (burned-label)
#'   polygons remain after QA. Consumed by the pool builder / orchestrator.
#'   Default `5L`.
#'
#' @param nrounds_max,early_stop Integer or `NULL`. XGBoost maximum number of
#'   boosting rounds and early-stopping patience, used at both the OOF and FINAL
#'   training stages. `NULL` uses the canonical defaults (`4000` / `80`). Stored
#'   in `cfg$train_control`.
#'
#' @param oof_seed_base,final_sampling_seed,final_seed Integer or `NULL`. The
#'   three distinct RNG seeds: out-of-fold block cross-validation; FINAL
#'   negative-pool sampling; and FINAL train/validation split plus XGBoost
#'   training. Keep them set for full reproducibility. `NULL` uses `42` for
#'   each. Stored in `cfg$train_control$seeds`.
#'
#' @param val_frac Numeric in `(0, 1)` or `NULL`. Inner validation fraction used
#'   by the inner train/validation split that selects the number of boosting
#'   rounds (early stopping) before the full-data refit, at both the OOF and
#'   FINAL stages. `NULL` uses `0.15`. Stored in `cfg$train_control`.
#'
#' @param group_col Character or `NULL`. Grouping column for the grouped
#'   train/validation split (keeps spatially grouped observations together).
#'   `NULL` uses `"block_id"`. Stored in `cfg$train_control`.
#'
#' @param impute_numeric,impute_factor_missing Character or `NULL`. The feature
#'   recipe's numeric imputation rule (`"median"` or `"zero"`) and the sentinel
#'   level used for missing factor values. `NULL` uses `"median"` / `"MISSING"`.
#'   Stored in `cfg$train_control` and applied identically at OOF and FINAL.
#'
#' @param cap_contextual,cap_spectral,cap_random,cap_otsu Numeric `>= 0`
#'   (`Inf` disables the cap) or `NULL`. The four negative-bucket caps, each
#'   expressed as a multiple of the number of burned labels (`n_burned`). They
#'   control the size of the contextual, spectral, random-background and Otsu
#'   negative buckets respectively (see \strong{Details}). `NULL` uses the
#'   general package defaults `0.25` / `1.0` / `1.0` / `1.0`. Stored in
#'   `cfg$train_control$caps`.
#'
#' @param feature_whitelist_override,feature_weights Optional feature controls,
#'   or `NULL`. `feature_whitelist_override` is a character subset of the
#'   package feature set that restricts which features the model may use (for
#'   example, dropping the hotspot block to define a no-hotspot profile).
#'   `feature_weights` is a named numeric vector of per-feature weights. Both
#'   are stored in `cfg$train_control` and consumed identically by the OOF and
#'   FINAL stages. `NULL` uses the full feature set with equal weights.
#'
#' @details
#' OOF uses the same capped negative-sampling policy as the final model, applied
#' independently within each training fold (inner early-stopping selection +
#' full-data refit for both the OOF folds and the FINAL model).
#'
#' Training eligibility is defined by EXPLICIT class, never by negation: only
#' explicit burned rows (positives) and explicit unburned rows that resolve to a
#' valid negative bucket (contextual, spectral, random, otsu) enter training;
#' review / keep / `NA` / unknown rows never become negatives. OOF and FINAL
#' share one internal eligibility resolver and one capping helper.
#'
#' @param model_params Named list or `NULL`. Optional \emph{partial} override of
#'   the canonical XGBoost hyper-parameter block stored in `cfg$model_params`
#'   (the canonical block without `scale_pos_weight`, which is computed at train
#'   time from the actual training labels). Supplying it merges your fields onto
#'   the canonical block field by field; for example `list(eta = 0.03)` changes
#'   only `eta` and leaves `max_depth`, `subsample`, ... at their canonical
#'   values. Setting `scale_pos_weight` here is rejected. Per-field provenance
#'   is recorded on `cfg$resolved_params_provenance$model_params`. `NULL` uses
#'   the canonical block unchanged.
#'
#' @param options Named list of technical and advanced options (not the core
#'   methodological controls, which have dedicated arguments above). All keys
#'   are optional. See the \strong{Options} section in \strong{Details} for the
#'   recognized keys, their defaults, and which ones affect reproducibility.
#'
#' @details
#' \subsection{Workflow structure}{
#'   The supervised workflow is organised as a sequence of stages, each driven
#'   by the configuration object this builder produces:
#'   \enumerate{
#'     \item \strong{Training pools} --- build the labelled and scoring pools,
#'       including the negative (unburned) pool;
#'       \code{build_supervised_training_pools()}.
#'     \item \strong{Spatial folds} --- assign spatial cross-validation folds;
#'       \code{make_spatial_folds()}.
#'     \item \strong{Feature extraction} --- compute the per-polygon feature
#'       table for the labelled and scoring pools;
#'       \code{extract_supervised_features()}.
#'     \item \strong{OOF diagnostics} --- out-of-fold cross-validated
#'       diagnostics of the model family; \code{run_oof_diagnostics()}.
#'     \item \strong{Final model} --- final feature selection plus refit on all
#'       labelled data; \code{train_final_burned_model()}.
#'     \item \strong{Scoring} --- score the candidate polygons with the final
#'       model; \code{score_supervised_burned_map()}.
#'     \item \strong{External map validation} --- optionally compare the
#'       thresholded map against an external reference;
#'       \code{validate_fire_maps()}.
#'   }
#'   \code{run_oneyear_supervised_pipeline()} runs stages 1--6 in one call from
#'   the same configuration object.
#' }
#'
#' \subsection{Negative pools}{
#'   The negative (unburned) training pool is assembled from four buckets:
#'   \itemize{
#'     \item \strong{contextual} --- deterministic-decision drops and
#'       context-derived negatives;
#'     \item \strong{spectral} --- spectrally selected negatives;
#'     \item \strong{random} --- random burnable-background cells;
#'     \item \strong{otsu} --- current-year Otsu-derived unburned patches.
#'   }
#'   Each bucket is capped by `cap_contextual` / `cap_spectral` / `cap_random` /
#'   `cap_otsu`, expressed as a multiple of the number of burned labels; an
#'   `Inf` cap disables capping for that bucket. The general package defaults
#'   are `0.25` / `1.0` / `1.0` / `1.0`. Specific experiments may use other
#'   values; for example a Phase B configuration raises the spectral cap to
#'   `2.0`.
#' }
#'
#' \subsection{OOF and FINAL}{
#'   OtsuFire uses a single supervised training procedure: the number of
#'   boosting rounds is selected using an inner validation split and early
#'   stopping, and the recipe and model are then refitted on all available
#'   training observations before prediction. There is no protocol choice.
#'
#'   In the out-of-fold (OOF) stage this procedure is applied per spatial fold:
#'   each outer fold selects `best_iteration` on an inner validation split
#'   drawn from the outer-train rows, then the recipe and model are refit on
#'   ALL outer-train rows at that `best_iteration` before predicting the
#'   held-out outer-test fold. The outer-test fold never enters imputation,
#'   feature selection, factor levels, weights, round selection or the refit,
#'   so the cross-validated diagnostics are leakage-free.
#'
#'   The FINAL stage applies the same procedure once to the whole labelled set:
#'   it selects the feature set, picks `best_iteration` on an inner validation
#'   split, and refits the recipe and model on all labelled observations. The
#'   OOF and FINAL stages share the same feature contract and the same recipe,
#'   so the diagnostics describe the same model family that is ultimately fit.
#' }
#'
#' \subsection{Feature recipe}{
#'   The feature recipe builds the base features, adds `_isNA` missingness
#'   indicators, imputes missing values (numeric rule from `impute_numeric`,
#'   factor sentinel from `impute_factor_missing`), and enforces a fixed final
#'   column order. Missing columns are handled gracefully so the same recipe
#'   applies to both the labelled and scoring tables. The trained model and its
#'   recipe are saved together so scoring reproduces the exact feature
#'   construction.
#' }
#'
#' \subsection{Full and no-hotspot profiles}{
#'   The "full" and "no-hotspot" profiles are not separate arguments. A
#'   no-hotspot profile is defined by leaving `hotspots = NULL` (and disabling
#'   the hotspot feature block at extraction time) and/or by restricting the
#'   active features with `feature_whitelist_override`. A full profile simply
#'   supplies hotspots and uses the complete feature set.
#' }
#'
#' \subsection{Scores}{
#'   The scoring stage attaches a per-polygon score `p_burned`. This is a model
#'   score, not necessarily a calibrated probability; thresholds for producing a
#'   binary map should be chosen with that in mind.
#' }
#'
#' \subsection{Validation}{
#'   Two distinct notions of validation apply. OOF diagnostics evaluate the
#'   model against the internal labels (the deterministic decision layer).
#'   \code{validate_fire_maps()} compares the thresholded supervised map against
#'   an external reference such as EFFIS, supplied via `reference_burned_map`.
#' }
#'
#' \subsection{Options}{
#'   The `options` list carries technical and advanced settings. Keys fall into
#'   three groups.
#'
#'   \strong{Common advanced options} (typically needed to wire a real run):
#'   \itemize{
#'     \item `data_base`, `composite_base` --- roots used to build the
#'       convention default paths for the optional inputs (see
#'       \strong{Description}).
#'     \item `result_name` --- run label used in convention paths.
#'     \item external tool paths (`python_exe`, `gdal_polygonize_script`,
#'       `gdalwarp_path`, `ogr2ogr_exe`), or a nested `tool_paths` list with the
#'       same names; surfaced on `cfg$tool_paths`.
#'     \item `legacy_otsu_mode` --- negative-pool Otsu mode; `"burnable_only"`
#'       (default behaviour) or `"corine"` (crops to `peninsula_shapefile`).
#'     \item `unb_verbose` --- verbose logging of the negative-pool builders.
#'   }
#'
#'   \strong{Reproducibility-sensitive options} (set them to make a run fully
#'   reproducible from the configuration alone): the negative-pool seeds
#'   `unb_random_seed` (default `42`) and `legacy_random_seed` (default `42`);
#'   the Otsu confidence thresholds `legacy_keep_hi` (`0.45`) and
#'   `legacy_drop_lo` (`0.15`); and the current-year temporal-adjustment
#'   thresholds `currentyear_preyear_overlap_thr` (`0.70`),
#'   `currentyear_hotspot_density_thr` (`0.001`) and
#'   `currentyear_temporal_penalty_floor` (`0.10`). The remaining
#'   negative-pool and legacy-Otsu knobs (`unb_excl_buffer_m`,
#'   `unb_n_random_cells`, `unb_random_rbr_q`, `unb_random_patch_size_cells`,
#'   `legacy_min_pixels`, `legacy_buffers_m`, `legacy_core_thr`, and related)
#'   default to their historical operational values, so omitting them
#'   reproduces the standard run.
#'
#'   \strong{Deprecated / unsupported options}: the negative-pool policy is
#'   fixed to `"all_sources"` and is no longer user-settable. Passing
#'   `negative_pool_policy = "all_sources"` is accepted as a no-op; any other
#'   value errors. The former `engine_root` / `supervised_engine_root` /
#'   `scripts_root` keys are no longer used (the engine is in-package). See
#'   \code{NEWS.md} for the full history.
#' }
#'
#' @return An S3 object of class `otsufire_supervised_burned_config`.
#'
#'   Stable public fields include:
#'   \itemize{
#'     \item `scenario`
#'     \item `target_year`
#'     \item `run_name`
#'     \item `inputs`
#'     \item `output_dir`
#'     \item `output_routes`
#'     \item `model_params`
#'     \item `train_control`
#'     \item `options`
#'     \item `tool_paths`
#'     \item `resolved_params_provenance`
#'     \item `negative_pool_policy` (informational constant, always
#'       `"all_sources"`)
#'     \item `min_burned_pool_n`
#'   }
#'
#'   `cfg$model_params` (the XGBoost block, without `scale_pos_weight`) and
#'   `cfg$train_control` (caps, seeds, rounds, validation fraction, grouping,
#'   imputation rules, feature whitelist/weights) are the
#'   resolved methodological parameters read by every downstream stage.
#'   Internal implementation details beyond these stable public fields are not a
#'   stable API and should not be relied upon by downstream user code.
#'
#' @examples
#' \dontrun{
#' ## Full configuration: balanced scenario, all inputs.
#' cfg <- build_supervised_burned_config(
#'   scenario             = "balanced",
#'   internal_decisions   = "2017/internal_decisions.gpkg",
#'   change_index         = "MinMin_2017_mosaic_res90m.tif",
#'   delayed_change_index = "Autumn/mean_mean_2017_mosaic.tif",
#'   hotspots             = "hotspots_2017.gpkg",
#'   topo                 = "Topography/elevation_slope.tif",
#'   corine_raster        = "Corine_Masks/CLC_2012_peninsula.tif",
#'   burnable_mask        = "Corine_Masks/burneable_mask_binary_corine_2012_ETRS89.tif",
#'   reference_burned_map = "EFFIS_2017.gpkg",
#'   target_year          = 2017,
#'   output_dir           = "results/",
#'   run_name             = "balanced_2017",
#'   cap_contextual       = 0.25,
#'   cap_spectral         = 2.0,
#'   cap_random           = 1.0,
#'   cap_otsu             = 1.0
#' )
#'
#' pools <- build_supervised_training_pools(cfg)
#'
#' ## No-hotspot profile: omit hotspots and drop the hotspot feature block via
#' ## a whitelist override (see inst/scripts/ for the canonical versioned
#' ## examples).
#' cfg_nohs <- build_supervised_burned_config(
#'   scenario                   = "balanced",
#'   internal_decisions         = "1994/internal_decisions.gpkg",
#'   change_index               = "MinMin_1994_mosaic_res90m.tif",
#'   hotspots                   = NULL,
#'   target_year                = 1994,
#'   output_dir                 = "results/",
#'   run_name                   = "no_hotspot_1994",
#'   feature_whitelist_override = c("rbr_med", "rbr_p90", "elev_med", "slope_med")
#' )
#' }
#'
#' @seealso
#' Supervised workflow, in order:
#' \itemize{
#'   \item \code{\link[=build_supervised_training_pools]{build_supervised_training_pools()}}
#'   \item \code{\link[=make_spatial_folds]{make_spatial_folds()}}
#'   \item \code{\link[=extract_supervised_features]{extract_supervised_features()}}
#'   \item \code{\link[=run_oof_diagnostics]{run_oof_diagnostics()}}
#'   \item \code{\link[=train_final_burned_model]{train_final_burned_model()}}
#'   \item \code{\link[=score_supervised_burned_map]{score_supervised_burned_map()}}
#'   \item \code{\link[=run_oneyear_supervised_pipeline]{run_oneyear_supervised_pipeline()}}
#'   \item \code{\link[=validate_supervised_execution]{validate_supervised_execution()}}
#'   \item \code{\link[=validate_fire_maps]{validate_fire_maps()}}
#' }
#'
#' @family workflow
#' @export
build_supervised_burned_config <- function(
    scenario = c("balanced", "original", "lax", "restrictive"),
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
    val_frac                   = NULL,
    group_col                  = NULL,
    impute_numeric             = NULL,
    impute_factor_missing      = NULL,
    cap_contextual             = NULL,
    cap_spectral               = NULL,
    cap_random                 = NULL,
    cap_otsu                   = NULL,
    feature_whitelist_override = NULL,
    feature_weights            = NULL,
    model_params               = NULL,
    options = list()
) {
  scenario <- match.arg(scenario)

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
  # the legacy unburned pipeline (`all_sources` mode) errors loudly if
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
    run_name = run_name, scenario = scenario
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
  train_control <- .of_resolve_supervised_train_control(
    nrounds_max                = nrounds_max,
    early_stop                 = early_stop,
    oof_seed_base              = oof_seed_base,
    final_sampling_seed        = final_sampling_seed,
    final_seed                 = final_seed,
    val_frac                   = val_frac,
    group_col                  = group_col,
    impute_numeric             = impute_numeric,
    impute_factor_missing      = impute_factor_missing,
    cap_contextual             = cap_contextual,
    cap_spectral               = cap_spectral,
    cap_random                 = cap_random,
    cap_otsu                   = cap_otsu,
    feature_whitelist_override = feature_whitelist_override,
    feature_weights            = feature_weights
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
    val_frac              = if (is.null(val_frac))              "default" else "user",
    group_col             = if (is.null(group_col))             "default" else "user",
    impute_numeric        = if (is.null(impute_numeric))        "default" else "user",
    impute_factor_missing = if (is.null(impute_factor_missing)) "default" else "user",
    cap_contextual        = if (is.null(cap_contextual))        "default" else "user",
    cap_spectral          = if (is.null(cap_spectral))          "default" else "user",
    cap_random            = if (is.null(cap_random))            "default" else "user",
    cap_otsu              = if (is.null(cap_otsu))              "default" else "user",
    feature_whitelist_override =
      if (is.null(feature_whitelist_override)) "default" else "user",
    feature_weights       = if (is.null(feature_weights))       "default" else "user",
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
    scenario                 = scenario,
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
    resolved_params_provenance = list(train_control = tc_provenance,
                                      model_params  = mp_provenance),
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
    cat("    caps          : contextual=", tc$caps$contextual,
        " spectral=", tc$caps$spectral, " random=", tc$caps$random,
        " otsu=", tc$caps$otsu, "\n", sep = "")
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
    val_frac, group_col, impute_numeric, impute_factor_missing,
    cap_contextual, cap_spectral, cap_random, cap_otsu,
    feature_whitelist_override, feature_weights) {

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

  # --- negative-bucket caps: numeric >= 0 (Inf allowed = disable) -------------
  chk_cap <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop(sprintf("'%s' must be a single number >= 0 (Inf disables the cap).",
                   nm), call. = FALSE)
    }
    as.numeric(v)
  }
  if (!is.null(cap_contextual)) tc$caps$contextual <- chk_cap(cap_contextual, "cap_contextual")
  if (!is.null(cap_spectral))   tc$caps$spectral   <- chk_cap(cap_spectral, "cap_spectral")
  if (!is.null(cap_random))     tc$caps$random     <- chk_cap(cap_random, "cap_random")
  if (!is.null(cap_otsu))       tc$caps$otsu       <- chk_cap(cap_otsu, "cap_otsu")

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

  # OtsuFire OOF always uses the capped negative-sampling policy (the SAME one
  # the FINAL model uses), applied independently within each training fold.
  # tc$oof_sampling is a FIXED traceability constant ("capped") from
  # .of_canonical_train_control(); it is not user-settable.

  tc
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
#' config itself is NULL). Every supervised stage that needs an input's file
#' path (orchestrator, pool builder, feature extractor) consumes the input
#' THROUGH this accessor, so `cfg$inputs` is the single source of truth and no
#' stage reconstructs an input path by filename convention behind the cfg's
#' back. A `NULL` return lets the caller apply its documented optional-input
#' behaviour (skip / convention-fallback for the standalone-script path).
#'
#' @param config supervised config S3 object, or `NULL`.
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
# in internal-sup-unburned-legacy.R, which is retained.)

#' @keywords internal
#' @noRd
.of_build_supervised_output_routes <- function(output_dir, target_year,
                                                run_name, scenario) {
  base <- file.path(output_dir, as.character(target_year), run_name,
                    "SUPERVISED", scenario)
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
