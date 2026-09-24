#' Extract predictor features for training and scoring polygons
#'
#' @description
#' Extract polygon-level predictor features for the labelled training pool
#' and the Otsu-guided patches prepared for supervised scoring.
#'
#' The function loads and aligns the supporting rasters, computes features for
#' both polygon sets, and optionally includes hotspot and shape features.
#'
#' Run this stage after [make_spatial_folds()] and before out-of-fold
#' diagnostics and final model training. The resulting feature layers are
#' used by the subsequent supervised stages.
#'
#' @param train_with_folds `sf` polygon object or GeoPackage path containing
#'   the labelled training polygons and fold assignments. Must contain
#'   `fire_uid` and `class`. When a path is supplied, the function reads the
#'   `train_with_folds` layer.
#' @param scoring_pool `sf` polygon object or GeoPackage path containing the
#'   Otsu-guided patches prepared for scoring. When a path is supplied, the
#'   function reads the `scoring_pool` layer.
#' @param config Required object of class `otsufire_supervised_burned_config`,
#'   created with [build_supervised_burned_config()]. Supplies raster inputs,
#'   the change-index template, optional supporting data, feature settings,
#'   and output locations.
#' @param use_hotspots Logical scalar. When `FALSE`, hotspot features are
#'   disabled. When `TRUE`, they are enabled only if the loaded hotspot layer
#'   contains observations and has a usable CRS. Default: `FALSE`.
#' @param use_shape Logical scalar. Set to `TRUE` to enable shape-feature
#'   extraction. With `FALSE`, shape features may still be enabled through
#'   `config$train_control$include_shape_features` or
#'   `config$options$use_shape`. Default: `FALSE`.
#' @param out_dir Character scalar or `NULL`. Output directory. When `NULL`,
#'   uses `config$output_routes$features_dir`, normally the run's
#'   `03_FEATURES` folder.
#' @param write_outputs Logical scalar. Whether to write the feature
#'   GeoPackage and RDS files. Default: `TRUE`.
#' @param .aligned_rasters Internal argument. Optional named list of already
#'   aligned rasters supplied by the pipeline. Leave as `NULL` for standalone
#'   use.
#' @param .hotspots_sf Internal argument. Optional preloaded hotspot layer
#'   supplied by the pipeline. Leave as `NULL` for standalone use.
#' @param .cor_groups Internal argument. Optional CORINE group lookup supplied
#'   by the pipeline. When `NULL`, uses the package lookup.
#'
#' @section Feature-extraction workflow:
#' The function:
#' 1. Loads the configured raster inputs and aligns them to the main
#'    change-index template.
#' 2. Loads hotspot observations when available.
#' 3. Runs `extract_features()` for the training and scoring polygons.
#' 4. Adds shape features when enabled.
#' 5. Writes the feature layers and RDS files when requested.
#' 6. Returns the training and scoring feature objects.
#'
#' Feature extraction does not train a model or select the final predictor
#' subset. Model-level feature selection and weighting are controlled through
#' the supervised configuration and applied during training.
#'
#' @section Raster inputs:
#' The supporting raster layers include:
#'
#' | Layer | Purpose |
#' |---|---|
#' | Summer RBR | Immediate spectral-response features. |
#' | Post-fire day of year | Timing features. |
#' | Autumn-winter RBR | Delayed-response and persistence features. |
#' | Elevation | Topographic features. |
#' | Slope | Topographic features. |
#' | CORINE land cover | Land-cover features. |
#'
#' The main change-index raster defines the alignment template.
#'
#' Inputs are resolved from the configuration, including conventional paths
#' where applicable. Feature availability depends on the supplied inputs; see
#' [build_supervised_burned_config()] for optional-input behaviour.
#'
#' When called directly, the function loads and aligns its supporting data.
#' The full pipeline can supply preloaded objects through the internal
#' arguments to avoid repeating those operations.
#'
#' @section Hotspot features:
#' For direct calls, hotspot extraction requires both:
#' * `use_hotspots = TRUE`;
#' * a loaded hotspot layer containing observations and a usable CRS.
#'
#' Supplying `hotspots` in the configuration alone does not enable this
#' feature block in a direct call. The default `use_hotspots = FALSE` forces
#' it off.
#'
#' Extracted hotspot features are available for modelling only if they are
#' also permitted by the active model feature whitelist.
#'
#' @section Shape and size features:
#' When enabled, the function computes and joins four geometric features to
#' both feature layers: `log_area`, `perim_m`, `compactness` and
#' `elongation`.
#'
#' The associated `area_ha` and `n_pix` fields are carried forward from the
#' pool-building stage.
#'
#' Shape extraction is enabled by `use_shape = TRUE` or by the corresponding
#' configuration settings. To keep this block disabled, leave
#' `use_shape = FALSE` and ensure that the configuration settings are also
#' disabled.
#'
#' For consistent extraction and model eligibility, configure shape features
#' through `include_shape_features = TRUE` in
#' [build_supervised_burned_config()].
#'
#' Shape and size can reflect how training examples were sampled. For
#' example, small background cells may differ systematically from burned
#' polygons. Evaluate this feature block against independent map references
#' as well as out-of-fold diagnostics.
#'
#' @section Existing feature outputs:
#' An existing `features_geometry.gpkg` is rebuilt when:
#' * it is older than the pool or fold input files;
#' * either the `train_features` or `scoring_features` layer is missing;
#' * its layer row counts differ from the corresponding input polygon counts.
#'
#' Otherwise, it is reused.
#'
#' These checks do not establish that every raster input, geometry, or
#' feature setting is unchanged. When changing supporting rasters, hotspot
#' settings, or shape-feature settings, use a separate output directory or
#' ensure that outdated feature products are regenerated.
#'
#' @section Output files:
#' When `write_outputs = TRUE`, the function writes the following files under
#' `out_dir`:
#'
#' | File | Contents |
#' |---|---|
#' | `features_geometry.gpkg` | Spatial feature layers named `train_features` and `scoring_features`. |
#' | `train_features.rds` | Training feature object. |
#' | `scoring_features.rds` | Scoring feature object. |
#'
#' The default output directory is the configured `03_FEATURES` folder.
#'
#' @return A named list containing feature objects and output paths.
#'
#' | Field | Contents |
#' |---|---|
#' | `train_features` | Labelled training features as an `sf` object. |
#' | `scoring_features` | Features for the candidate scoring pool as an `sf` object. |
#' | `features_geometry_gpkg` | Path to the feature GeoPackage when written. |
#' | `train_features_rds` | Path to the training-feature RDS file when written. |
#' | `scoring_features_rds` | Path to the scoring-feature RDS file when written. |
#'
#' Use the returned feature objects directly when working in memory. With
#' `write_outputs = FALSE`, do not assume that the returned output locations
#' contain newly written files.
#'
#' @seealso [build_supervised_burned_config()],
#'   [build_supervised_training_pools()], [make_spatial_folds()],
#'   `extract_features()`, [run_oof_diagnostics()],
#'   [train_final_burned_model()], [score_supervised_burned_map()],
#'   [validate_supervised_execution()], [run_oneyear_supervised_pipeline()].
#'
#' @examples
#' \dontrun{
#' # Configure the supervised workflow
#' config <- build_supervised_burned_config(
#'   internal_decisions = "data/internal_decisions_2022.gpkg",
#'   change_index = "data/RBR_2022.tif",
#'   delayed_change_index = "data/RBR_delayed_2022.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   topo = "data/elevation_slope.tif",
#'   corine_raster = "data/land_cover_2022.tif",
#'   burnable_mask = "data/burnable_mask_2022.tif",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022",
#'   include_shape_features = FALSE
#' )
#'
#' # Build training and scoring pools
#' pools <- build_supervised_training_pools(
#'   config = config,
#'   write_outputs = TRUE
#' )
#'
#' # Assign spatial folds
#' folds <- make_spatial_folds(
#'   train_labelled = pools$pools_gpkg,
#'   config = config
#' )
#'
#' # Review the partition before continuing
#' if (!isTRUE(folds$selected$ok)) {
#'   stop("Review the fallback partition before continuing.")
#' }
#'
#' # Extract features, explicitly enabling hotspot features
#' features <- extract_supervised_features(
#'   train_with_folds = folds$train_with_folds_gpkg,
#'   scoring_pool = pools$pools_gpkg,
#'   config = config,
#'   use_hotspots = TRUE,
#'   write_outputs = TRUE
#' )
#'
#' # Inspect the training feature columns
#' names(features$train_features)
#'
#' # Inspect the scoring feature table
#' head(sf::st_drop_geometry(features$scoring_features))
#'
#' # Locate the written feature GeoPackage
#' features$features_geometry_gpkg
#'
#' # Extract a separate version without hotspot features
#' # A separate directory avoids reusing the previous feature file.
#' features_no_hotspots <- extract_supervised_features(
#'   train_with_folds = folds$train_with_folds,
#'   scoring_pool = pools$scoring_pool,
#'   config = config,
#'   use_hotspots = FALSE,
#'   out_dir = file.path(
#'     config$output_routes$features_dir,
#'     "no_hotspots"
#'   ),
#'   write_outputs = TRUE
#' )
#' }
#'
#' @family workflow
#' @export
extract_supervised_features <- function(train_with_folds, scoring_pool, config,
                                        use_hotspots = FALSE,
                                        use_shape = FALSE,
                                        out_dir = NULL,
                                        write_outputs = TRUE,
                                        .aligned_rasters = NULL,
                                        .hotspots_sf = NULL,
                                        .cor_groups = NULL) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(config) || is.null(config) ||
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (missing(train_with_folds) || is.null(train_with_folds)) {
    stop("'train_with_folds' is required.", call. = FALSE)
  }
  if (!(inherits(train_with_folds, c("sf", "data.frame", "SpatVector")) ||
        (is.character(train_with_folds) && length(train_with_folds) == 1L))) {
    stop("'train_with_folds' must be sf / data.frame, SpatVector, or a GPKG ",
         "path.", call. = FALSE)
  }
  if (missing(scoring_pool) || is.null(scoring_pool)) {
    stop("'scoring_pool' is required.", call. = FALSE)
  }
  if (!(inherits(scoring_pool, c("sf", "data.frame", "SpatVector")) ||
        (is.character(scoring_pool) && length(scoring_pool) == 1L))) {
    stop("'scoring_pool' must be sf / data.frame, SpatVector, or a GPKG path.",
         call. = FALSE)
  }
  if (!is.logical(use_hotspots) || length(use_hotspots) != 1L ||
      is.na(use_hotspots)) {
    stop("'use_hotspots' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(use_shape) || length(use_shape) != 1L ||
      is.na(use_shape)) {
    stop("'use_shape' must be TRUE or FALSE.", call. = FALSE)
  }
  # OPTIONAL shape/size block (OtsuFire 0.12.0). When the explicit
  # use_shape arg is left at its FALSE default, drive it from the config:
  # the resolved train_control flag include_shape_features (canonical
  # single source of truth) OR the convenience options$use_shape toggle.
  # An explicit use_shape = TRUE always wins. OFF everywhere ->
  # byte-identical extraction.
  if (!isTRUE(use_shape)) {
    use_shape <- isTRUE(config$train_control$include_shape_features) ||
      isTRUE(config$options$use_shape)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L ||
      is.na(write_outputs)) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  .of_check_input_file(config$inputs$change_index, "change_index")

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve out_dir + target_year/scenario + canonical output path from
  #    config. out_dir == config$output_routes$features_dir (the 03_FEATURES
  #    folder); features_gpkg is always <out_dir>/features_geometry.gpkg, which
  #    is what the orchestrator (and the OOF / final-model stages downstream)
  #    consume.
  # ---------------------------------------------------------------------------
  if (is.null(out_dir)) out_dir <- config$output_routes$features_dir
  if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L ||
      !nzchar(out_dir)) {
    stop("Could not resolve 'out_dir' (03_FEATURES) from config.",
         call. = FALSE)
  }
  target_year <- config$target_year
  scenario    <- config$scenario
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  features_gpkg <- file.path(out_dir, "features_geometry.gpkg")

  # ---------------------------------------------------------------------------
  # 2) Safe-IO kit. Built here via make_sf_safety_kit() so the function does NOT
  #    depend on the orchestrator's closure environment. We use only its
  #    safe_read_gpkg / msg helpers (the same closures STEP B3 used).
  # ---------------------------------------------------------------------------
  sfkit <- make_sf_safety_kit(
    dirs        = list(`99_LOGS_EMPTY` = file.path(out_dir, "99_LOGS_EMPTY")),
    result_dir  = config$output_routes$base,
    target_year = target_year
  )
  msg            <- sfkit$msg
  safe_read_gpkg <- sfkit$safe_read_gpkg

  # ---------------------------------------------------------------------------
  # 3) Resolve train_with_folds / scoring_pool to on-disk GPKG paths. The
  #    staleness check + the engine read both go through safe_read_gpkg() on a
  #    GPKG path, so when an in-memory object is supplied we materialise it to a
  #    temp GPKG under the canonical layer name first.
  # ---------------------------------------------------------------------------
  resolve_layer_gpkg <- function(x, layer, what) {
    if (is.character(x) && length(x) == 1L) {
      if (!file.exists(x)) {
        stop(sprintf("'%s' path does not exist: %s", what, x), call. = FALSE)
      }
      return(x)
    }
    xs <- if (inherits(x, "SpatVector")) {
      sf::st_as_sf(x)
    } else if (inherits(x, "sf")) {
      x
    } else {
      sf::st_as_sf(as.data.frame(x))
    }
    p <- tempfile(fileext = ".gpkg")
    sf::st_write(xs, p, layer = layer, quiet = TRUE, delete_dsn = TRUE)
    p
  }
  train_with_folds_gpkg <- resolve_layer_gpkg(train_with_folds,
                                              "train_with_folds",
                                              "train_with_folds")
  gpkg_pools_out        <- resolve_layer_gpkg(scoring_pool, "scoring_pool",
                                              "scoring_pool")

  # ---------------------------------------------------------------------------
  # 4) Load + align the raster support stack, hotspots, and CORINE groups.
  #    Self-contained by default: resolve INPUT paths from config$options
  #    (data_base / composite_base / result_name) exactly as the orchestrator
  #    SETUP block does, and run the SAME align_to_template() alignment. When
  #    the orchestrator delegates it passes its already-aligned objects via the
  #    dot-prefixed args -> no double load, byte-identical engine call.
  # ---------------------------------------------------------------------------
  cor_groups_use <- .cor_groups %||% cor_groups

  if (is.null(.aligned_rasters)) {
    data_base      <- config$options$data_base
    composite_base <- config$options$composite_base
    result_name    <- config$options$result_name %||% "Min_Min"
    if (is.null(data_base) || !nzchar(data_base)) {
      stop("extract_supervised_features() needs config$options$data_base to ",
           "locate the raster support stack (topography / CORINE / autumn ",
           "composite). Provide it via build_supervised_burned_config(",
           "options = list(data_base = ..., composite_base = ...)).",
           call. = FALSE)
    }
    if (is.null(composite_base) || !nzchar(composite_base)) {
      stop("extract_supervised_features() needs config$options$composite_base ",
           "to locate the change-index / autumn composites.", call. = FALSE)
    }

    # Gate 1B: single shared cfg-input accessor (.of_sup_input_path,
    # supervised-config.R). Convention fallback only when the cfg carries no
    # explicit path (single source of truth).
    .cfg_input_path <- function(name) .of_sup_input_path(config, name)

    corine_year         <- sfkit$get_corine_year(target_year)
    # Gate 1B PIECE 2: feature-support paths consumed from cfg$inputs;
    # convention fallback only when the cfg carries no explicit path.
    topo_path           <- .cfg_input_path("topo") %||%
      file.path(data_base, "Topography", "elevation_slope.tif")
    corine_raster_path  <- .cfg_input_path("corine_raster") %||%
      file.path(data_base, "Corine_Masks",
                paste0("CLC_", corine_year, "_peninsula.tif"))
    # change_index: the REQUIRED, validated cfg$inputs$change_index field,
    # CONSUMED from cfg (not reconstructed by filename convention). Convention
    # is the fallback only for an in-memory / pathless change_index spec.
    one_year_tif <- .cfg_input_path("change_index") %||% {
      cand <- file.path(
        composite_base, result_name,
        paste0("MinMin_", target_year, "_mosaic_res90m.tif")
      )
      if (!file.exists(cand)) {
        cand <- file.path(
          composite_base, "Min_Min",
          paste0("MinMin_", target_year, "_mosaic_res90m.tif")
        )
      }
      cand
    }
    # delayed_change_index: optional cfg$inputs field; convention fallback.
    rbr_aw_tif <- .cfg_input_path("delayed_change_index") %||% file.path(
      composite_base, "Autumn",
      paste0("mean_mean_", target_year, "_mosaic.tif")
    )
    stopifnot(file.exists(one_year_tif))
    stopifnot(file.exists(corine_raster_path))
    stopifnot(file.exists(topo_path))
    stopifnot(file.exists(rbr_aw_tif))

    topo       <- terra::rast(topo_path)
    rbr_stack  <- terra::rast(one_year_tif)   # band1 = RBR summer, band2 = DOY
    rbr_summer <- rbr_stack[[1]]
    doy_post   <- rbr_stack[[2]]

    template_r <- rbr_summer

    doy_post <- align_to_template(
      r = doy_post, template = template_r, method = "near",
      name = "doy_post", verbose = TRUE
    )
    dem_r <- align_to_template(
      r = topo[[1]], template = template_r, method = "bilinear",
      name = "dem", verbose = TRUE
    )
    slope_r <- align_to_template(
      r = topo[[2]], template = template_r, method = "bilinear",
      name = "slope", verbose = TRUE
    )
    corine_r <- align_to_template(
      r = terra::rast(corine_raster_path), template = template_r,
      method = "near", name = "corine_r", verbose = TRUE
    )
    rbr_aw <- align_to_template(
      r = terra::rast(rbr_aw_tif)[[1]], template = template_r,
      method = "bilinear", name = "rbr_aw", verbose = TRUE
    )
  } else {
    ar <- .aligned_rasters
    rbr_summer <- ar$rbr_summer
    doy_post   <- ar$doy_post
    rbr_aw     <- ar$rbr_aw
    dem_r      <- ar$dem_r
    slope_r    <- ar$slope_r
    corine_r   <- ar$corine_r
  }

  if (is.null(.hotspots_sf)) {
    data_base     <- config$options$data_base
    # hotspots: OPTIONAL cfg$inputs$hotspots (allow_null = TRUE so pre-MODIS
    # years run with hotspots = NULL). CONSUMED from cfg; convention path is the
    # fallback only when the cfg carries no explicit hotspots path AND a
    # data_base is available. NULL -> treated as an absent layer below.
    hotspots_path <- .of_sup_input_path(config, "hotspots") %||%
      (if (!is.null(data_base) && nzchar(data_base))
         file.path(data_base, "Hotspots",
                   paste0("hotspots_iberia_", target_year, ".geojson"))
       else NULL)
    hotspots_sf_base <- if (!is.null(hotspots_path) && file.exists(hotspots_path)) {
      sf::read_sf(hotspots_path)
    } else {
      msg("WARNING: hotspots not found: %s (use_hotspots is disabled)",
          hotspots_path %||% "<none: cfg$inputs$hotspots = NULL>")
      sf::st_sf(
        frp = numeric(),
        confidence = numeric(),
        year = integer(),
        month = integer(),
        geometry = sf::st_sfc(crs = sf::NA_crs_)
      )[0, ]
    }
  } else {
    hotspots_sf_base <- .hotspots_sf
  }

  # ---------------------------------------------------------------------------
  # 5) Staleness-rebuild: decide reuse vs rebuild. Verbatim MOVE of the
  #    orchestrator STEP B3 logic (features older than pools/folds, missing
  #    layers, or row-count mismatch -> rebuild).
  # ---------------------------------------------------------------------------
  feature_layers <- function(dsn) {
    if (!file.exists(dsn)) return(character(0))
    tryCatch(as.character(sf::st_layers(dsn)$name),
             error = function(e) character(0))
  }
  feature_layer_n <- function(dsn, layer) {
    if (!file.exists(dsn)) return(NA_integer_)
    if (!layer %in% feature_layers(dsn)) return(NA_integer_)
    out <- tryCatch(
      nrow(sf::read_sf(dsn, layer = layer, quiet = TRUE)),
      error = function(e) NA_integer_
    )
    as.integer(out)
  }
  safe_remove_features_gpkg <- function(path) {
    if (!file.exists(path)) return(invisible(TRUE))
    unlink(path, force = TRUE)
    if (file.exists(path)) file.remove(path)
    if (file.exists(path)) {
      stop("No se pudo borrar features_geometry.gpkg antes de reconstruir: ",
           path, call. = FALSE)
    }
    invisible(TRUE)
  }

  stopifnot(!is.null(train_with_folds_gpkg),
            file.exists(train_with_folds_gpkg))
  stopifnot(file.exists(gpkg_pools_out))

  should_rebuild_features <- FALSE
  rebuild_reasons <- character(0)

  if (file.exists(features_gpkg)) {
    feat_time  <- file.info(features_gpkg)$mtime
    pools_time <- file.info(gpkg_pools_out)$mtime
    folds_time <- file.info(train_with_folds_gpkg)$mtime
    feat_layers <- feature_layers(features_gpkg)

    if (isTRUE(feat_time < pools_time)) {
      should_rebuild_features <- TRUE
      rebuild_reasons <- c(rebuild_reasons,
                           "features_geometry.gpkg is older than 01_POOLS")
    }
    if (isTRUE(feat_time < folds_time)) {
      should_rebuild_features <- TRUE
      rebuild_reasons <- c(rebuild_reasons,
                           "features_geometry.gpkg is older than 02_FOLDS")
    }

    n_train_expected <- tryCatch(
      nrow(safe_read_gpkg(train_with_folds_gpkg, "train_with_folds",
                          "train_with_folds_count_check")),
      error = function(e) NA_integer_
    )
    n_unl_expected <- tryCatch(
      nrow(safe_read_gpkg(gpkg_pools_out, "scoring_pool",
                          "scoring_pool_count_check")),
      error = function(e) NA_integer_
    )
    n_train_feat <- feature_layer_n(features_gpkg, "train_features")
    n_unl_feat   <- feature_layer_n(features_gpkg, "scoring_features")

    if (!all(c("train_features", "scoring_features") %in% feat_layers)) {
      should_rebuild_features <- TRUE
      rebuild_reasons <- c(
        rebuild_reasons,
        sprintf(
          "features_geometry.gpkg does not contain every required layer (current layers: %s)",
          paste(feat_layers, collapse = ", ")
        )
      )
    }

    if (!is.na(n_train_expected) && !is.na(n_train_feat) &&
        n_train_expected != n_train_feat) {
      should_rebuild_features <- TRUE
      rebuild_reasons <- c(rebuild_reasons,
        sprintf("train_features has %d rows and train_with_folds %d",
                n_train_feat, n_train_expected))
    }
    if (!is.na(n_unl_expected) && !is.na(n_unl_feat) &&
        n_unl_expected != n_unl_feat) {
      should_rebuild_features <- TRUE
      rebuild_reasons <- c(rebuild_reasons,
        sprintf("scoring_features has %d rows and scoring_pool %d",
                n_unl_feat, n_unl_expected))
    }

    if (isTRUE(should_rebuild_features)) {
      msg("STEP B3 - features_geometry.gpkg exists but is stale. Rebuilding it.")
      for (rr in unique(rebuild_reasons)) msg("  * %s", rr)
      safe_remove_features_gpkg(features_gpkg)
    } else {
      msg("STEP B3 - features_geometry.gpkg ya existe, SKIP: %s", features_gpkg)
    }
  }

  # ---------------------------------------------------------------------------
  # 6) Build features when the canonical gpkg is absent (verbatim MOVE of the
  #    orchestrator STEP B3.1 - B3.3 block). When it already exists and is fresh
  #    we reuse it (feat stays NULL and objects are read back below).
  # ---------------------------------------------------------------------------
  feat <- NULL

  if (!file.exists(features_gpkg)) {
    msg("STEP B3.1 - Read train_with_folds + scoring_pool")
    train_folds <- safe_read_gpkg(train_with_folds_gpkg, "train_with_folds",
                                  "train_with_folds")
    scoring_candidates <- safe_read_gpkg(gpkg_pools_out, "scoring_pool",
                                         "scoring_pool")

    hotspots_sf <- hotspots_sf_base
    if (nrow(hotspots_sf) > 0) {
      if (is.na(sf::st_crs(hotspots_sf))) {
        msg("WARNING: hotspots_sf has NA CRS; disabling hotspots.")
        hotspots_sf <- hotspots_sf[0, ]
      } else if (sf::st_crs(hotspots_sf) != sf::st_crs(train_folds)) {
        hotspots_sf <- sf::st_transform(hotspots_sf, sf::st_crs(train_folds))
      }
    }

    # Runtime hotspot flag: the orchestrator's use_hotspots_flag was
    # nrow(hotspots_sf) > 0. The public use_hotspots arg can force it OFF; the
    # orchestrator delegation passes its runtime flag so behaviour is unchanged.
    use_hotspots_flag <- isTRUE(use_hotspots) && nrow(hotspots_sf) > 0

    msg("STEP B3.3 - extract_features (train + unlabeled)")

    msg("STEP B3.2b - Check raster geometry before extract_features")
    raster_list <- list(
      rbr_summer = rbr_summer,
      doy_post   = doy_post,
      rbr_aw     = rbr_aw,
      dem_r      = dem_r,
      slope_r    = slope_r,
      corine_r   = corine_r
    )
    for (nm in names(raster_list)) {
      rr <- raster_list[[nm]]
      ok <- terra::compareGeom(rr, rbr_summer, stopOnError = FALSE)
      msg("Raster check: %s | ok=%s | nrow=%s ncol=%s | res=(%s,%s) | origin=(%s,%s)",
          nm, ok, nrow(rr), ncol(rr),
          terra::res(rr)[1], terra::res(rr)[2],
          terra::origin(rr)[1], terra::origin(rr)[2])
      if (!isTRUE(ok)) {
        stop(sprintf("Raster geometry mismatch before extract_features: %s does not match rbr_summer", nm))
      }
    }

    feat <- extract_features(
      train_folds = train_folds,
      unlabeled   = scoring_candidates,

      id_col    = "fire_uid",
      class_col = "class",
      pos_lab   = "burned",
      neg_lab   = "unburned",

      build_features = TRUE,

      rbr_summer = rbr_summer,
      doy_post   = doy_post,
      rbr_aw     = rbr_aw,
      nbr_pre    = NULL,
      nbr_post   = NULL,
      dnbr       = NULL,
      dem        = dem_r,
      slope      = slope_r,
      corine_r   = corine_r,

      hotspots   = hotspots_sf,
      frp_col    = "frp",
      conf_col   = "confidence",
      year_col   = "year",
      month_col  = "month",
      hs_start_month = 6,
      hs_end_month   = 10,
      hi_conf_thr    = 0.8,

      year_target = target_year,
      hs_year_min_available = 2000,
      hs_buffer_m = 1000,
      hs_missing_value = -9999,
      hs_use_season_filter = TRUE,

      cor_groups = cor_groups_use,

      use_doy      = TRUE,
      use_aw       = TRUE,
      use_nbr      = FALSE,
      use_hotspots = use_hotspots_flag,
      # OPTIONAL shape/size block (OtsuFire 0.12.0). OFF by default ->
      # no shape join in the engine -> byte-identical features.
      use_shape    = use_shape,

      max_cells_in_memory      = NULL,
      return_features          = TRUE,
      return_features_geometry = TRUE,
      # write_outputs gates on-disk materialisation. When TRUE (default, the
      # orchestrator path) save_features_dir == out_dir and the gpkg/rds are
      # written exactly as before -> byte-identical. When FALSE nothing is
      # written and only the in-memory feature objects are returned.
      save_features_dir        = if (isTRUE(write_outputs)) out_dir else NULL,
      save_features_format     = "rds",
      save_features_gpkg       = isTRUE(write_outputs),
      train_features_basename  = "train_features",
      scoring_features_basename = "scoring_features",
      train_features_layer     = "train_features",
      scoring_features_layer   = "scoring_features",
      verbose                  = TRUE
    )

    if (isTRUE(write_outputs)) {
      stopifnot(file.exists(features_gpkg))
      msg("DONE - features saved: %s", features_gpkg)
    }
  }

  # ---------------------------------------------------------------------------
  # 7) Assemble the documented return (objects + paths). When the gpkg was
  #    (re)built this run, the engine return carries the feature sf objects;
  #    when it was reused we read them back from the canonical layers so the
  #    return shape is identical either way.
  # ---------------------------------------------------------------------------
  train_features_rds   <- file.path(out_dir, "train_features.rds")
  scoring_features_rds <- file.path(out_dir, "scoring_features.rds")

  if (!is.null(feat)) {
    train_features   <- feat$train_feat
    scoring_features <- feat$unl_feat
  } else {
    train_features <- tryCatch(
      sf::read_sf(features_gpkg, layer = "train_features", quiet = TRUE),
      error = function(e) NULL
    )
    scoring_features <- tryCatch(
      sf::read_sf(features_gpkg, layer = "scoring_features", quiet = TRUE),
      error = function(e) NULL
    )
  }

  list(
    train_features         = train_features,
    scoring_features       = scoring_features,
    features_geometry_gpkg = features_gpkg,
    train_features_rds     = train_features_rds,
    scoring_features_rds   = scoring_features_rds
  )
}
