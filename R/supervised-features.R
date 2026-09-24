#' Extract predictor features for labelled and scoring polygons
#'
#' @description
#' Computes the per-polygon predictor features the burned-area model is trained
#' and scored on, for both the labelled training polygons and the deterministic
#' scoring universe. This is the feature-extraction stage of the supervised
#' pipeline; its `features_geometry.gpkg` output then feeds the OOF and
#' final-model stages.
#'
#' Run it after the folds stage and before [run_oof_diagnostics()]. Pass the
#' train-with-folds layer, the scoring pool, and the configuration from
#' [build_supervised_burned_config()]. Calling this function directly produces
#' the same `03_FEATURES` outputs as the equivalent stage of the full pipeline
#' [run_oneyear_supervised_pipeline()].
#'
#' @section What it does:
#' \enumerate{
#'   \item Loads and aligns the supervised raster support stack (RBR summer,
#'     post-fire DOY, autumn-winter RBR, DEM, slope, CORINE) to the
#'     change-index master template.
#'   \item Optionally loads the target-year hotspot layer.
#'   \item Runs the patch-level feature engine on the
#'     labelled and scoring polygons with the pipeline's standard settings.
#'   \item Optionally appends shape/size features (see `use_shape`).
#'   \item Writes `features_geometry.gpkg` (carrying the `train_features` and
#'     `scoring_features` layers) plus `train_features.rds` and
#'     `scoring_features.rds`, and returns the in-memory objects.
#' }
#'
#' @details
#' \strong{Raster loading.} The function resolves the raster input paths from
#' `config$options$data_base` / `composite_base` / `result_name` and aligns
#' them to the change-index template itself.
#'
#' \strong{Staleness rebuild.} When `features_geometry.gpkg` already exists it
#' is rebuilt when it is older than the pools / folds inputs, when it is missing
#' the `train_features` / `scoring_features` layers, or when its layer row
#' counts do not match the train-with-folds / scoring-pool row counts; otherwise
#' it is reused.
#'
#' \strong{Hotspots.} The engine's hotspot block is enabled only when the loaded
#' hotspot layer has rows with a usable CRS. Setting `use_hotspots = FALSE`
#' (the default) forces the block off regardless of the loaded layer.
#'
#' @param train_with_folds sf POLYGON layer OR a single GPKG path. The labelled
#'   training layer with fold columns produced by the folds stage (the
#'   `train_with_folds` object / the `train_with_folds` layer of
#'   `02_FOLDS/<year>_train_with_folds_<bs>m.gpkg`). Must carry `fire_uid` and
#'   `class`. When a path is supplied the `train_with_folds` layer is read.
#' @param scoring_pool sf POLYGON layer OR a single GPKG path. The deterministic
#'   scoring universe produced by the pools stage (the `scoring_pool` object /
#'   the `scoring_pool` layer of `01_POOLS/<year>_<run_label>_pools.gpkg`). When
#'   a path is supplied the `scoring_pool` layer is read.
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Used to derive `out_dir` (the
#'   `03_FEATURES` folder), the raster INPUT paths
#'   (`config$options$data_base` / `composite_base` / `result_name`), the
#'   change-index master template, and `target_year` / `scenario`.
#' @param use_hotspots Logical. When `FALSE` (default) the hotspot feature block
#'   is forced off. When `TRUE` the block is enabled if and only if the loaded
#'   hotspot layer has rows with a usable CRS. Default `FALSE`.
#' @param use_shape Logical. Optional shape/size feature block, off by default
#'   (`FALSE`). When `TRUE` the engine computes four geometric columns
#'   (`log_area`, `perim_m`, `compactness`, `elongation`; `area_ha` / `n_pix`
#'   already ride along as pool-builder columns) and left-joins them onto both
#'   the train and scoring feature layers. When left at its `FALSE` default it
#'   is driven from the config (`cfg$train_control$include_shape_features` or
#'   the convenience `options$use_shape`); an explicit `TRUE` always wins. With
#'   the block off no join runs and the extraction output is unchanged. See the
#'   sampling-bias note in [build_supervised_burned_config()].
#' @param out_dir Character or `NULL`. Output folder for the `03_FEATURES`
#'   outputs. Defaults to `config$output_routes$features_dir`.
#' @param write_outputs Logical. Whether to write the `features_geometry.gpkg`
#'   and `*_features.rds` outputs. Default `TRUE`.
#' @param .aligned_rasters,.hotspots_sf,.cor_groups For internal use by
#'   [run_oneyear_supervised_pipeline()], which passes rasters and layers it has
#'   already loaded. Leave them as `NULL` (the default).
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `train_features` — the labelled training features (sf), the
#'       `train_features` layer contents.
#'     \item `scoring_features` — the scoring-universe features (sf), the
#'       `scoring_features` layer contents.
#'     \item `features_geometry_gpkg` — written path of `features_geometry.gpkg`.
#'     \item `train_features_rds`, `scoring_features_rds` — written RDS paths.
#'   }
#'
#' @seealso
#' [build_supervised_burned_config()], [make_spatial_folds()],
#' [run_oof_diagnostics()], [train_final_burned_model()],
#' [score_supervised_burned_map()], [validate_supervised_execution()],
#' [run_oneyear_supervised_pipeline()]
#'
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   run_label = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L,
#'   options = list(data_base = "D:/FIRE", composite_base = "D:/FIRE/Composites")
#' )
#' feats <- extract_supervised_features(
#'   train_with_folds = "02_FOLDS/2017_train_with_folds_2000m.gpkg",
#'   scoring_pool     = "01_POOLS/2017_balanced_pools.gpkg",
#'   config = cfg
#' )
#' feats$features_geometry_gpkg
#' }
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
