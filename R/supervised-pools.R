#' Build burned and negative training pools for supervised learning
#'
#' @description
#' Supervised-pipeline stage A (training-pool construction). This is the
#' standalone, exported implementation of the pools stage that the one-year
#' orchestrator [run_oneyear_supervised_pipeline()] delegates to, so calling it
#' directly produces byte-identical `01_POOLS` outputs to a full run.
#'
#' Unlike the other modular stages (which wrap a single engine call), this is a
#' SEQUENCE (orchestrator steps A1-A6): it reads the deterministic
#' `internal_decisions` layer, runs the conservative QA relabelling
#' ([audit_deterministic_pools()]), assembles the burned / review / scoring
#' pools, builds the scenario `all_sources` negative pool from BOTH unburned
#' builders ([build_unburned_from_deterministic_decisions()] = deterministic
#' drops + random burnable background, and
#' [build_unburned_from_legacy_pipeline()] = Otsu current-year residual
#' patches), merges the labelled training set (`burned + unburned`), and writes
#' the canonical `01_POOLS/<year>_<scenario>_pools.gpkg` (layers `burned_pool`,
#' `unburned_pool`, `review_pool`, `scoring_pool`, `train_labeled`, plus the
#' auxiliary unburned layers) and the QA CSV/GPKG audit files.
#'
#' @details
#' **Self-contained config resolution.** Every value the inlined orchestrator
#' STEP A consumed as an injected bare-name local is resolved here straight
#' from `config`: the ~30 `UNB_*` / `UNB_LEGACY_*` negative-pool parameters from
#' `config$options$*` (with the SAME defaults the engine bindings use, so
#' default behaviour is byte-identical), the four external tool paths
#' (`python_exe`, `gdal_polygonize_script`, `gdalwarp_path`, `ogr2ogr_exe`) from
#' `config$tool_paths$*`, the INPUT imagery roots
#' (`config$options$data_base` / `composite_base` / `result_name`), and the
#' `01_POOLS` output folder + canonical `pools.gpkg` path from
#' `config$output_routes`. The deterministic decisions GPKG is the
#' user-supplied path (`config$inputs$internal_decisions`), exactly as the
#' dispatcher resolves it.
#'
#' **RNG determinism.** The random burnable-background sampling
#' (`UNB_RANDOM_SEED = 42`) and the legacy Otsu sampling
#' (`UNB_LEGACY_RANDOM_SEED = 42`) are forwarded verbatim and the two builders
#' are invoked in the SAME order (deterministic-decisions first, legacy second)
#' so the RNG stream — and therefore the sampled negatives — are identical to a
#' full orchestrator run.
#'
#' @param config Required `otsufire_supervised_burned_config` (from
#'   [build_supervised_burned_config()]). Source of every input path, output
#'   route, `UNB_*` parameter and tool path used by the stage.
#' @param deterministic_decisions Optional sf, SpatVector, or single GPKG path
#'   to the deterministic decisions (layer `internal_decisions`). When `NULL`
#'   (default) the canonical `config$inputs$internal_decisions` path is used
#'   (the single source of truth already threaded into the config).
#' @param write_outputs Logical. Whether to write the pools GPKG + QA files.
#'   Default `TRUE`.
#' @param overwrite Logical. Forwarded to the deterministic-decisions unburned
#'   builder (`overwrite_output`). Default `FALSE`.
#'
#' @return A named list with both the objects and the written paths:
#'   \itemize{
#'     \item `burned_pool`, `unburned_pool`, `review_pool`, `scoring_pool` —
#'       the four pool sf objects.
#'     \item `train_labeled` — the merged labelled training set
#'       (`burned + unburned`), the `train_labeled` layer contents consumed by
#'       [make_spatial_folds()].
#'     \item `pools_gpkg` — written path of
#'       `01_POOLS/<year>_<scenario>_pools.gpkg` (consumed by the folds and
#'       feature stages).
#'     \item `audited_internal_gpkg`, `qa_summary_csv`, `qa_transitions_csv`,
#'       `qa_reasons_csv` — the QA audit file paths.
#'     \item `unburned_hard`, `unburned_random`, `unburned_final_raw`,
#'       `exclusion_buffer` — the auxiliary unburned sub-pools.
#'     \item `b4_audit` — Gate 1C.2 audit record for the random burnable
#'       background: cell accounting at each stage (valid index, burnable,
#'       per-rule exclusions), the percentile domain + value, eligible cells,
#'       finally selected observations, and an explicit confirmation that NO
#'       selected observation falls outside the burnable domain.
#'     \item `neg_pool_fingerprint` — deterministic `list(text, checksum)`
#'       identity of the negative pool. It folds in the B4 burnable-domain
#'       decision, the aligned-mask hash, the percentile config + value, the
#'       RNG seeds, the exclusions and the relevant inputs, so any cache/pool
#'       whose fingerprint predates the burnable-domain restriction is
#'       INVALIDATED rather than silently reused. Also written as a sidecar
#'       `01_POOLS/<year>_<scenario>_neg_pool_fingerprint.txt`.
#'   }
#'
#' @section Burnable-domain restriction (Gate 1C.2, INTENTIONAL pool change):
#' The B4 random burnable-background bucket now computes its change-index
#' percentile AND draws its sample ONLY within the (aligned) burnable domain.
#' Previously the percentile was taken over the WHOLE change-index raster, so a
#' run on this code produces a DIFFERENT negative pool than older runs. This is
#' a deliberate methodological correction; the `neg_pool_fingerprint` changes
#' accordingly so stale pools cannot be reused. The legacy single-fit and the
#' nested_refit consumers both operate on this one upstream-built pool.
#'
#' @family workflow
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- build_supervised_burned_config(
#'   scenario = "balanced", internal_decisions = "decisions.gpkg",
#'   change_index = "rbr.tif", target_year = 2017L,
#'   options = list(data_base = "D:/FIRE", composite_base = "D:/FIRE/Composites")
#' )
#' pools <- build_supervised_training_pools(cfg)
#' pools$pools_gpkg
#' }
build_supervised_training_pools <- function(config,
                                            deterministic_decisions = NULL,
                                            write_outputs = TRUE,
                                            overwrite = FALSE) {
  # ---------------------------------------------------------------------------
  # 0) Validation
  # ---------------------------------------------------------------------------
  if (missing(config) || is.null(config) ||
      !inherits(config, "otsufire_supervised_burned_config")) {
    stop("'config' must be created by build_supervised_burned_config().",
         call. = FALSE)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L ||
      is.na(write_outputs)) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(deterministic_decisions)) {
    if (!(inherits(deterministic_decisions, c("sf", "SpatVector")) ||
          (is.character(deterministic_decisions) &&
           length(deterministic_decisions) == 1L &&
           file.exists(deterministic_decisions)))) {
      stop("'deterministic_decisions' must be sf, SpatVector, or an existing ",
           "path.", call. = FALSE)
    }
  }
  .of_check_input_file(config$inputs$internal_decisions, "internal_decisions")

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }

  # ---------------------------------------------------------------------------
  # 1) Resolve year / scenario / output routes from config. result_dir is the
  #    SUPERVISED output base (== config$output_routes$base, the dispatcher's
  #    supervised_output_base / the orchestrator's result_dir). dirs mirrors the
  #    orchestrator subfolders that STEP A touches.
  # ---------------------------------------------------------------------------
  target_year <- config$target_year
  scenario    <- config$scenario

  result_dir <- config$output_routes$base
  if (is.null(result_dir) || !is.character(result_dir) ||
      length(result_dir) != 1L || !nzchar(result_dir)) {
    stop("Could not resolve output base (config$output_routes$base).",
         call. = FALSE)
  }
  pools_dir  <- config$output_routes$pools_dir %||%
    file.path(result_dir, "01_POOLS")
  pools_gpkg <- config$output_routes$pools_gpkg %||%
    file.path(pools_dir, sprintf("%d_%s_pools.gpkg", target_year, scenario))

  dirs <- list(
    `01_POOLS`      = pools_dir,
    `99_LOGS_EMPTY` = file.path(result_dir, "99_LOGS_EMPTY")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  # ---------------------------------------------------------------------------
  # 2) Resolve INPUT imagery roots from config$options (the INPUT authority,
  #    exactly as the orchestrator / extract_supervised_features do). The
  #    unburned builders need these to locate composites / masks.
  # ---------------------------------------------------------------------------
  data_base      <- config$options$data_base
  composite_base <- config$options$composite_base
  result_name    <- config$options$result_name %||% "Min_Min"
  if (is.null(data_base) || !nzchar(data_base)) {
    stop("build_supervised_training_pools() needs config$options$data_base to ",
         "locate the deterministic / composite inputs for the negative-pool ",
         "builders. Provide it via build_supervised_burned_config(",
         "options = list(data_base = ..., composite_base = ...)).",
         call. = FALSE)
  }
  if (is.null(composite_base) || !nzchar(composite_base)) {
    stop("build_supervised_training_pools() needs config$options$composite_base ",
         "to locate the change-index / autumn composites.", call. = FALSE)
  }

  # §N+25 (2026-06-05): the burnable mask / CORINE raster / peninsula border
  # are wired RUN inputs. Read the on-disk path from config$inputs when the
  # cfg carries one (the convention default is filled at config-build time);
  # pass NULL otherwise so each unburned builder falls back to EXACTLY its
  # historical convention path (byte-identical when nothing is supplied).
  # Gate 1B: single shared cfg-input accessor (.of_sup_input_path,
  # supervised-config.R) — same resolver the orchestrator + feature extractor use.
  .cfg_input_path <- function(name) .of_sup_input_path(config, name)
  cfg_burnable_mask_path  <- .cfg_input_path("burnable_mask")
  cfg_corine_raster_path  <- .cfg_input_path("corine_raster")
  cfg_peninsula_shapefile <- .cfg_input_path("peninsula_shapefile")

  # one_year_tif: the MAIN change-index raster, used here as the severity raster
  # for the legacy Otsu builder. It is the REQUIRED, validated cfg$inputs$change_index
  # field; CONSUMED from cfg$inputs (single source of truth), NOT reconstructed
  # by the MinMin_<year>_mosaic_res90m.tif filename convention. The convention is
  # the fallback ONLY when the cfg carries no on-disk change_index path (e.g. an
  # in-memory SpatRaster spec) so behaviour stays byte-identical otherwise.
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

  # Unburned output roots, replicating the orchestrator ROOT PATHS block.
  deterministic_dir <- file.path(
    data_base, "Results", target_year, result_name, "DETERMINISTIC", scenario
  )
  legacy_unb_root_dir <- file.path(result_dir, "_LEGACY_UNBURNED")
  unb_out_gpkg <- file.path(
    deterministic_dir, "UNBURNED",
    sprintf("%d_%s_unburned.gpkg", target_year, scenario)
  )
  dir.create(dirname(unb_out_gpkg), recursive = TRUE, showWarnings = FALSE)
  dir.create(legacy_unb_root_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # 3) Resolve the deterministic decisions GPKG path. Prefer the explicit
  #    `deterministic_decisions` arg; else the canonical config path (the same
  #    $path the dispatcher extracts). In-memory sf/SpatVector inputs are
  #    materialised to a temp GPKG under the canonical `internal_decisions`
  #    layer so the read + both builders consume an on-disk path identically.
  # ---------------------------------------------------------------------------
  internal_decisions_path <- if (!is.null(deterministic_decisions)) {
    if (is.character(deterministic_decisions)) {
      deterministic_decisions
    } else {
      xs <- if (inherits(deterministic_decisions, "SpatVector")) {
        sf::st_as_sf(deterministic_decisions)
      } else {
        deterministic_decisions
      }
      p <- tempfile(fileext = ".gpkg")
      sf::st_write(xs, p, layer = "internal_decisions", quiet = TRUE,
                   delete_dsn = TRUE)
      p
    }
  } else {
    config$inputs$internal_decisions$path
  }
  if (is.null(internal_decisions_path) ||
      !is.character(internal_decisions_path) ||
      length(internal_decisions_path) != 1L ||
      !nzchar(internal_decisions_path)) {
    stop("Could not resolve an on-disk 'internal_decisions' GPKG path from ",
         "'deterministic_decisions' / config$inputs$internal_decisions.",
         call. = FALSE)
  }
  if (!file.exists(internal_decisions_path)) {
    stop(sprintf("internal_decisions GPKG does not exist: %s",
                 internal_decisions_path), call. = FALSE)
  }
  internal_decisions_gpkg <- internal_decisions_path

  # ---------------------------------------------------------------------------
  # 4) Resolve all UNB_* / UNB_LEGACY_* parameters + tool paths from config,
  #    with the SAME defaults the engine bindings use (byte-identical default
  #    behaviour). See .of_supervised_engine_bindings() for the canonical list.
  # ---------------------------------------------------------------------------
  min_burned_pool_n <- config$min_burned_pool_n %||% 5L

  UNB_VERBOSE       <- config$options$unb_verbose %||% TRUE

  # deterministic-decisions builder (det drops + random burnable background)
  UNB_EXCL_BUFFER_M        <- config$options$unb_excl_buffer_m %||% 500
  UNB_N_RANDOM_CELLS       <- config$options$unb_n_random_cells %||% 1500
  UNB_RANDOM_RBR_Q         <- config$options$unb_random_rbr_q %||% 0.50
  UNB_RANDOM_SEED          <- config$options$unb_random_seed %||% 42
  UNB_RANDOM_PATCH_SIZE_CELLS <-
    config$options$unb_random_patch_size_cells %||% 3

  # legacy Otsu builder
  UNB_LEGACY_CODE_DIR      <- config$options$legacy_code_dir
  UNB_LEGACY_OTSU_MODE     <- config$options$legacy_otsu_mode %||% "burnable_only"
  UNB_LEGACY_OTSU_THRESHOLD <- config$options$legacy_otsu_threshold %||% 0
  UNB_LEGACY_REFERENCE_OTSU_THRESHOLD <-
    config$options$legacy_reference_otsu_threshold %||% 100
  UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE <-
    config$options$legacy_min_otsu_threshold_value %||% 0
  UNB_LEGACY_MIN_PIXELS    <- config$options$legacy_min_pixels %||% 8
  UNB_LEGACY_BUFFERS_M     <- config$options$legacy_buffers_m %||% 90
  UNB_LEGACY_CORE_THR      <- config$options$legacy_core_thr %||% 0.60
  UNB_LEGACY_ALPHA_BOOST   <- config$options$legacy_alpha_boost %||% 0.25
  UNB_LEGACY_MIN_BASE_BOOST <- config$options$legacy_min_base_boost %||% 0.35
  UNB_LEGACY_DIST_POWER    <- config$options$legacy_dist_power %||% 1
  UNB_LEGACY_KEEP_HI       <- config$options$legacy_keep_hi %||% 0.45
  UNB_LEGACY_DROP_LO       <- config$options$legacy_drop_lo %||% 0.15
  UNB_LEGACY_EXCL_BUFFER_M <- config$options$legacy_excl_buffer_m %||% 0
  UNB_LEGACY_MIN_AREA_HA   <- config$options$legacy_min_area_ha %||% 0
  UNB_LEGACY_SAMPLE_N      <- config$options$legacy_sample_n %||% 2000
  UNB_LEGACY_SAMPLE_PROPS  <- config$options$legacy_sample_props %||%
    c(drop = 0.70, review = 0.25, keep = 0.05)
  UNB_LEGACY_REUSE_EXISTING <- config$options$legacy_reuse_existing %||% TRUE
  UNB_LEGACY_WRITE_OUTPUT   <- config$options$legacy_write_output %||% TRUE
  UNB_LEGACY_USE_DROP      <- config$options$legacy_use_drop   %||% TRUE
  UNB_LEGACY_USE_REVIEW    <- config$options$legacy_use_review %||% FALSE
  UNB_LEGACY_USE_KEEP      <- config$options$legacy_use_keep   %||% FALSE
  UNB_LEGACY_RANDOM_SEED   <- config$options$legacy_random_seed %||% 42L
  UNB_LEGACY_DROP_MAX_S_PATCH   <- config$options$legacy_drop_max_s_patch   %||% 0.15
  UNB_LEGACY_REVIEW_MAX_S_PATCH <- config$options$legacy_review_max_s_patch %||% 0.45
  UNB_LEGACY_KEEP_MAX_S_PATCH   <- config$options$legacy_keep_max_s_patch   %||% 0.70
  # D4a (2026-06-05): default FALSE -> empty Otsu legacy pool ERRORS instead of
  # silently degrading to deterministic_direct. Opt back in via
  # config$options$allow_empty_otsu_pool = TRUE.
  UNB_ALLOW_EMPTY_OTSU_POOL     <- config$options$allow_empty_otsu_pool %||% FALSE

  # external tool paths
  python_exe             <- config$tool_paths$python_exe
  gdal_polygonize_script <- config$tool_paths$gdal_polygonize_script
  gdalwarp_path          <- config$tool_paths$gdalwarp_path
  ogr2ogr_exe            <- config$tool_paths$ogr2ogr_exe

  # ---------------------------------------------------------------------------
  # 5) Safe-IO kit. Built here via make_sf_safety_kit() so the function does NOT
  #    depend on the orchestrator's closure environment. STEP A used these exact
  #    closures.
  # ---------------------------------------------------------------------------
  sfkit <- make_sf_safety_kit(
    dirs        = dirs,
    result_dir  = result_dir,
    target_year = target_year
  )
  msg                  <- sfkit$msg
  sanitize_polygons    <- sfkit$sanitize_polygons
  ensure_area_ha       <- sfkit$ensure_area_ha
  drop_empty_sf        <- sfkit$drop_empty
  check_sf             <- sfkit$check_sf
  check_sf_if_nonempty <- sfkit$check_sf_if_nonempty
  to_crs_safe          <- sfkit$to_crs_safe
  safe_write_gpkg      <- sfkit$safe_write_gpkg
  time_step            <- sfkit$time_step

  mutate     <- dplyr::mutate
  filter     <- dplyr::filter
  bind_rows  <- dplyr::bind_rows
  row_number <- dplyr::row_number

  gpkg_pools_out <- pools_gpkg

  # ===========================================================================
  # A1. Read internal_decisions from deterministic pipeline (scenario)
  # ===========================================================================
  msg("STEP A1 - Read internal_decisions")
  internal_sf <- time_step("A1 Read internal_decisions", {
    msg("Using internal_decisions_gpkg: %s", internal_decisions_gpkg)

    x <- sf::read_sf(internal_decisions_gpkg, layer = "internal_decisions") |>
      mutate(
        class = as.character(class_final)
      ) |>
      sanitize_polygons() |>
      ensure_area_ha()

    x <- drop_empty_sf(x, tag = "internal_sf", dump_dir = dirs$`99_LOGS_EMPTY`)
    check_sf(x, "internal_sf")
    x
  })

  crs_master <- sf::st_crs(internal_sf)
  stopifnot(!is.na(crs_master))

  # ===========================================================================
  # A2. Build clean internal set
  # ===========================================================================
  msg("STEP A2 - Build internal_clean")
  internal_clean <- time_step("A2 Build internal_clean", {
    x <- internal_sf |>
      mutate(source = "internal") |>
      sanitize_polygons() |>
      ensure_area_ha()

    x <- drop_empty_sf(x, tag = "internal_clean", dump_dir = dirs$`99_LOGS_EMPTY`)
    check_sf(x, "internal_clean")
    x
  })

  # ===========================================================================
  # A3. QA burned/review pools from deterministic decisions
  # ===========================================================================
  msg("STEP A3 - QA deterministic pools and build burned_pool + review_pool")
  poolsA <- time_step("A3 QA deterministic pools + build supervised pools", {
    qa_det <- audit_deterministic_pools(
      internal_sf = internal_clean,
      out_dir = dirs$`01_POOLS`,
      prefix = sprintf("%d_%s_deterministic_pool_qa", target_year, scenario),
      save_outputs = TRUE,
      verbose = TRUE,
      raw_class_col = "class"
    )

    internal_qc <- qa_det$audited_sf |>
      mutate(
        raw_class = as.character(raw_class),
        class = as.character(class_audited)
      ) |>
      sanitize_polygons() |>
      ensure_area_ha()

    burned_pool <- internal_qc |>
      filter(class == "keep") |>
      mutate(source = "internal_keep_qc") |>
      sanitize_polygons() |>
      ensure_area_ha() |>
      mutate(
        class = "burned",
        fire_uid = paste0(target_year, "_", scenario, "_B_", row_number())
      )

    burned_pool <- drop_empty_sf(
      burned_pool,
      tag = "burned_pool",
      dump_dir = dirs$`99_LOGS_EMPTY`
    )
    check_sf_if_nonempty(burned_pool, "burned_pool")

    review_pool <- internal_qc |>
      filter(class == "review") |>
      mutate(source = "internal_review_qc") |>
      sanitize_polygons() |>
      ensure_area_ha() |>
      mutate(
        class = "review",
        fire_uid = paste0(target_year, "_", scenario, "_R_", row_number())
      )

    review_pool <- drop_empty_sf(
      review_pool,
      tag = "review_pool",
      dump_dir = dirs$`99_LOGS_EMPTY`
    )
    check_sf_if_nonempty(review_pool, "review_pool")

    scoring_pool <- internal_qc |>
      mutate(
        class = as.character(raw_class),
        source = "deterministic_all_qc",
        fire_uid = paste0(target_year, "_", scenario, "_D_", row_number())
      ) |>
      sanitize_polygons() |>
      ensure_area_ha()

    scoring_pool <- drop_empty_sf(
      scoring_pool,
      tag = "scoring_pool",
      dump_dir = dirs$`99_LOGS_EMPTY`
    )
    check_sf_if_nonempty(scoring_pool, "scoring_pool")

    list(
      internal_qc = internal_qc,
      burned_pool = burned_pool,
      review_pool = review_pool,
      scoring_pool = scoring_pool,
      qa_det = qa_det
    )
  })

  internal_qc  <- poolsA$internal_qc
  burned_pool  <- poolsA$burned_pool
  review_pool  <- poolsA$review_pool
  scoring_pool <- poolsA$scoring_pool

  # ===========================================================================
  # A3b. Pool sanity checks
  # ===========================================================================
  msg("STEP A3b - Sanity checks on pools")
  time_step("A3b Pool sanity checks", {

    n_burned <- nrow(burned_pool)
    n_rev    <- nrow(review_pool)

    msg("Pools: burned=%d | review=%d", n_burned, n_rev)

    if (n_burned < min_burned_pool_n) {
      stop(sprintf(
        "Supervised run aborted: insufficient burned pool (n=%d, required>=%d). Year likely has too few detected fires for supervised modeling.",
        n_burned, min_burned_pool_n
      ))
    }

    if (n_rev == 0) {
      stop("Supervised run aborted: review pool is empty.")
    }

    NULL
  })

  # ===========================================================================
  # A4. Generate unburned datasets (all_sources policy):
  #   build_unburned_from_deterministic_decisions() -> det drops + random bg
  #   build_unburned_from_legacy_pipeline()         -> Otsu current-year patches
  # RNG: forwarded UNB_RANDOM_SEED / UNB_LEGACY_RANDOM_SEED; det builder runs
  # BEFORE the legacy builder, identical to the orchestrator, so the random
  # stream is preserved.
  # ===========================================================================
  msg("STEP A4 - Generate unburned datasets [all_sources]")
  unb <- time_step("A4 Generate unburned datasets", {
    # ---- Part 1: deterministic drops + random burnable background ----
    res_det <- build_unburned_from_deterministic_decisions(
      target_year             = target_year,
      scenario_name           = scenario,
      data_base               = data_base,
      result_name             = result_name,
      composite_base          = composite_base,
      exclude_buffer_m        = UNB_EXCL_BUFFER_M,
      n_random_cells          = UNB_N_RANDOM_CELLS,
      random_rbr_q            = UNB_RANDOM_RBR_Q,
      random_seed             = UNB_RANDOM_SEED,
      random_patch_size_cells = UNB_RANDOM_PATCH_SIZE_CELLS,
      overwrite_output        = overwrite,
      out_gpkg                = unb_out_gpkg,
      internal_decisions_path = internal_decisions_path,
      # §N+25: wired RUN input (NULL -> convention inside the builder).
      burnable_mask_path      = cfg_burnable_mask_path,
      # Gate 1B PIECE 3: immediate change-index route from cfg (== one_year_tif).
      severity_raster_path    = one_year_tif,
      verbose                 = UNB_VERBOSE
    )
    out_path <- res_det$out_gpkg
    stopifnot(file.exists(out_path))

    # Gate 1C.2: B4 audit + aligned-burnable-mask hash, carried up so the
    # negative-pool fingerprint (below) can fold in the ACTUAL burnable domain
    # used for the random background (not just the mask path).
    b4_audit           <- res_det$b4_audit
    burnable_mask_hash <- res_det$burnable_mask_hash

    unburned_hard   <- to_crs_safe(res_det$unburned_hard, crs_master) |>
      mutate(source = "deterministic_drop_hard")
    unburned_random <- to_crs_safe(res_det$unburned_random, crs_master) |>
      mutate(source = "random_burnable_background")
    det_final_raw   <- to_crs_safe(res_det$unburned_final, crs_master) |>
      mutate(source = as.character(source))
    exclusion_buffer <- to_crs_safe(res_det$exclusion_buffer, crs_master)

    # ---- Part 2: Otsu current-year unburned patches ----
    res_otsu <- build_unburned_from_legacy_pipeline(
      target_year              = target_year,
      scenario_name            = scenario,
      data_base                = data_base,
      result_name              = result_name,
      composite_base           = composite_base,
      severity_raster_path     = one_year_tif,
      internal_decisions_path  = internal_decisions_path,
      # §N+25: wired RUN inputs (NULL -> convention inside the builder).
      burnable_mask_path       = cfg_burnable_mask_path,
      corine_raster_path       = cfg_corine_raster_path,
      peninsula_shapefile      = cfg_peninsula_shapefile,
      legacy_code_dir          = UNB_LEGACY_CODE_DIR,
      python_exe               = python_exe,
      gdal_polygonize_script   = gdal_polygonize_script,
      gdalwarp_path            = gdalwarp_path,
      ogr2ogr_exe              = ogr2ogr_exe,
      otsu_mode                = UNB_LEGACY_OTSU_MODE,
      otsu_threshold           = UNB_LEGACY_OTSU_THRESHOLD,
      reference_otsu_threshold = UNB_LEGACY_REFERENCE_OTSU_THRESHOLD,
      min_otsu_threshold_value = UNB_LEGACY_MIN_OTSU_THRESHOLD_VALUE,
      min_pixels               = UNB_LEGACY_MIN_PIXELS,
      buffers_m                = UNB_LEGACY_BUFFERS_M,
      core_thr                 = UNB_LEGACY_CORE_THR,
      alpha_boost              = UNB_LEGACY_ALPHA_BOOST,
      min_base_boost           = UNB_LEGACY_MIN_BASE_BOOST,
      dist_power               = UNB_LEGACY_DIST_POWER,
      keep_hi                  = UNB_LEGACY_KEEP_HI,
      drop_lo                  = UNB_LEGACY_DROP_LO,
      use_drop                 = UNB_LEGACY_USE_DROP,
      use_review               = UNB_LEGACY_USE_REVIEW,
      use_keep                 = UNB_LEGACY_USE_KEEP,
      drop_max_s_patch         = UNB_LEGACY_DROP_MAX_S_PATCH,
      review_max_s_patch       = UNB_LEGACY_REVIEW_MAX_S_PATCH,
      keep_max_s_patch         = UNB_LEGACY_KEEP_MAX_S_PATCH,
      sample_n                 = UNB_LEGACY_SAMPLE_N,
      sample_props             = UNB_LEGACY_SAMPLE_PROPS,
      exclude_buffer_m         = UNB_LEGACY_EXCL_BUFFER_M,
      min_area_ha              = UNB_LEGACY_MIN_AREA_HA,
      random_seed              = UNB_LEGACY_RANDOM_SEED,
      reuse_existing           = UNB_LEGACY_REUSE_EXISTING,
      write_unburned           = UNB_LEGACY_WRITE_OUTPUT,
      out_root_dir             = legacy_unb_root_dir,
      allow_empty_otsu_pool    = UNB_ALLOW_EMPTY_OTSU_POOL,
      verbose                  = UNB_VERBOSE
    )
    otsu_pool    <- to_crs_safe(res_otsu$unburned$legacy_unburned_pool, crs_master)
    otsu_sampled <- to_crs_safe(res_otsu$unburned$legacy_unburned_sampled, crs_master)
    # AS08 (0.5.0): the former defensive fallback
    #   if (nrow(otsu_sampled) == 0L && nrow(otsu_pool) > 0L) otsu_sampled <- otsu_pool
    # was DEAD: sample_stratified_legacy_unburned() returns its input `x`
    # unchanged when target_n <= 0 or nrow(x) == 0, and otherwise returns a
    # non-empty subset of x. So `otsu_sampled` is empty only when `otsu_pool`
    # (the same `combined` pool) is also empty -> the branch condition
    # (sampled empty AND pool non-empty) can never hold. Replaced with the
    # invariant assertion the audit recommended; never fires on any run, so
    # behaviour is byte-identical.
    stopifnot(nrow(otsu_sampled) > 0L || nrow(otsu_pool) == 0L)

    otsu_raw <- otsu_sampled |>
      mutate(
        source = "otsu_patch_residual",
        neg_type = dplyr::case_when(
          as.character(.data$legacy_decision) == "drop"   ~ "otsu_patch_drop",
          as.character(.data$legacy_decision) == "review" ~ "otsu_patch_review",
          as.character(.data$legacy_decision) == "keep"   ~ "otsu_patch_keep",
          TRUE ~ as.character(.data$neg_type)
        )
      )

    # ---- Part 3: combine all sources ----
    # Geometry-safe binding: drop geometry, bind data frames, re-attach.
    .bind_sf_safe <- function(a, b) {
      df_a <- sf::st_drop_geometry(a)
      df_b <- sf::st_drop_geometry(b)
      geom <- c(sf::st_geometry(a), sf::st_geometry(b))
      combined <- dplyr::bind_rows(df_a, df_b)
      combined[["geometry"]] <- geom
      sf::st_as_sf(combined, sf_column_name = "geometry",
                   crs = sf::st_crs(a))
    }
    unburned_final_raw <- .bind_sf_safe(det_final_raw, otsu_raw)

    check_sf(unburned_final_raw, "unburned_final_raw")
    check_sf_if_nonempty(unburned_hard, "unburned_hard")
    check_sf_if_nonempty(unburned_random, "unburned_random")
    check_sf_if_nonempty(exclusion_buffer, "exclusion_buffer")

    unburned_pool <- unburned_final_raw |>
      sanitize_polygons() |>
      ensure_area_ha() |>
      mutate(
        class = "unburned",
        fire_uid = paste0(target_year, "_", scenario, "_U_", row_number())
      )

    unburned_pool <- drop_empty_sf(
      unburned_pool,
      tag = "unburned_pool",
      dump_dir = dirs$`99_LOGS_EMPTY`
    )
    check_sf(unburned_pool, "unburned_pool")

    msg("UNB out path: %s", out_path %||% NA_character_)
    msg("Unburned final count: %d", nrow(unburned_pool))

    list(
      wrapper_out_gpkg   = out_path,
      unburned_hard      = unburned_hard,
      unburned_random    = unburned_random,
      unburned_final_raw = unburned_final_raw,
      exclusion_buffer   = exclusion_buffer,
      unburned_pool      = unburned_pool,
      b4_audit           = b4_audit,
      burnable_mask_hash = burnable_mask_hash
    )
  })

  unburned_hard      <- unb$unburned_hard
  unburned_random    <- unb$unburned_random
  unburned_final_raw <- unb$unburned_final_raw
  exclusion_buffer   <- unb$exclusion_buffer
  unburned_pool      <- unb$unburned_pool
  b4_audit           <- unb$b4_audit
  burnable_mask_hash <- unb$burnable_mask_hash

  # ===========================================================================
  # A4b. Negative-pool fingerprint (Gate 1C.2). Deterministic identity of the
  # negative pool sufficient to INVALIDATE any cached/stale pool whose B4 domain
  # decision, burnable mask, percentile/config, seed, exclusions or relevant
  # inputs predate this run. Built via the same base-R fingerprinter the legacy
  # cache uses (legacy_param_fingerprint_unb_legacy). It folds in AT LEAST:
  #   - the B4 domain decision (burnable-restricted) + the aligned-mask hash,
  #   - the percentile config (random_rbr_q) + its computed value,
  #   - the RNG seeds (det + legacy),
  #   - the exclusions (exclude_buffer_m, patch size, n_random_cells),
  #   - the relevant inputs (change index, burnable mask path, decisions path).
  # Because the pool is built ONCE here and SHARED downstream, both the legacy
  # (single-fit) and nested_refit consumers operate on exactly this pool /
  # fingerprint — there is no second, independently-built pool.
  # ===========================================================================
  neg_pool_fingerprint <- legacy_param_fingerprint_unb_legacy(list(
    neg_pool_policy          = "all_sources",
    b4_domain_decision       = b4_audit$domain_decision %||% "burnable_restricted",
    b4_burnable_mask_hash    = burnable_mask_hash %||% NA_character_,
    b4_burnable_mask_path    = cfg_burnable_mask_path %||% "",
    b4_random_rbr_q          = UNB_RANDOM_RBR_Q,
    b4_percentile_value      = b4_audit$percentile_value %||% NA_real_,
    b4_random_seed           = UNB_RANDOM_SEED,
    b4_n_random_cells        = UNB_N_RANDOM_CELLS,
    b4_random_patch_size     = UNB_RANDOM_PATCH_SIZE_CELLS,
    b4_exclude_buffer_m      = UNB_EXCL_BUFFER_M,
    b4_n_percentile_domain   = b4_audit$n_percentile_domain %||% NA_integer_,
    b4_n_eligible_cells      = b4_audit$n_eligible_cells %||% NA_integer_,
    legacy_random_seed       = UNB_LEGACY_RANDOM_SEED,
    legacy_sample_n          = UNB_LEGACY_SAMPLE_N,
    legacy_otsu_mode         = UNB_LEGACY_OTSU_MODE,
    target_year              = target_year,
    scenario                 = scenario,
    change_index             = one_year_tif,
    internal_decisions_path  = internal_decisions_path
  ))
  msg("Negative-pool fingerprint (Gate 1C.2): CHECKSUM=%s",
      neg_pool_fingerprint$checksum)

  # ===========================================================================
  # A5. Final labelled training set = burned + unburned
  # ===========================================================================
  msg("STEP A5 - Build train_labeled = burned + unburned")
  train_labeled <- time_step("A5 Build train_labeled", {
    x <- bind_rows(
      burned_pool,
      unburned_pool
    ) |>
      sanitize_polygons() |>
      ensure_area_ha()

    x <- drop_empty_sf(x, tag = "train_labeled", dump_dir = dirs$`99_LOGS_EMPTY`)
    check_sf(x, "train_labeled")
    stopifnot(!anyDuplicated(x$fire_uid))

    msg("train_labeled counts: burned=%d | unburned=%d",
        sum(as.character(x$class) == "burned"),
        sum(as.character(x$class) == "unburned"))
    x
  })

  # ===========================================================================
  # A6. Write pools + auxiliary unburned datasets
  # ===========================================================================
  if (isTRUE(write_outputs)) {
    msg("STEP A6 - Write pools to 01_POOLS")
    time_step("A6 Write pools GPKG", {
      safe_write_gpkg(burned_pool,       gpkg_pools_out, layer = "burned_pool",       tag = "burned_pool_write")
      safe_write_gpkg(unburned_pool,     gpkg_pools_out, layer = "unburned_pool",     tag = "unburned_pool_write")
      safe_write_gpkg(review_pool,       gpkg_pools_out, layer = "review_pool",       tag = "review_pool_write")
      safe_write_gpkg(scoring_pool,      gpkg_pools_out, layer = "scoring_pool",      tag = "scoring_pool_write")
      safe_write_gpkg(train_labeled,     gpkg_pools_out, layer = "train_labeled",     tag = "train_labeled_write")

      if (nrow(unburned_hard) > 0) {
        safe_write_gpkg(unburned_hard, gpkg_pools_out, layer = "unburned_hard", tag = "unburned_hard_write")
      }
      if (nrow(unburned_random) > 0) {
        safe_write_gpkg(unburned_random, gpkg_pools_out, layer = "unburned_random", tag = "unburned_random_write")
      }
      if (nrow(unburned_final_raw) > 0) {
        safe_write_gpkg(unburned_final_raw, gpkg_pools_out, layer = "unburned_final_raw", tag = "unburned_final_raw_write")
      }
      if (nrow(exclusion_buffer) > 0) {
        safe_write_gpkg(exclusion_buffer, gpkg_pools_out, layer = "exclusion_buffer", tag = "exclusion_buffer_write")
      }

      NULL
    })

    # Gate 1C.2: persist the negative-pool fingerprint sidecar next to the pools
    # GPKG. It is the authoritative identity of THIS (burnable-domain-restricted)
    # negative pool. Any cache/consumer that compares against a fingerprint
    # predating this decision will MISMATCH and must rebuild rather than silently
    # reuse a burnable-domain-naive pool. Overwritten every run with the pool.
    fp_path <- file.path(
      pools_dir,
      sprintf("%d_%s_neg_pool_fingerprint.txt", target_year, scenario)
    )
    fp_token <- c(
      "# OtsuFire negative-pool fingerprint (Gate 1C.2).",
      "# B4 random background percentile+sampling restricted to burnable domain.",
      "# A cache whose fingerprint differs from this MUST rebuild (no silent reuse).",
      sprintf("CHECKSUM=%s", neg_pool_fingerprint$checksum),
      "---",
      neg_pool_fingerprint$text
    )
    tryCatch(
      writeLines(fp_token, fp_path),
      error = function(e)
        warning("Could not write negative-pool fingerprint sidecar (",
                conditionMessage(e), ").", call. = FALSE)
    )

    msg("DONE - pools saved: %s", gpkg_pools_out)
  }

  # ---------------------------------------------------------------------------
  # 6) Assemble the documented return (objects + paths). QA audit file paths
  #    mirror the audit_deterministic_pools(save_outputs=TRUE) naming under the
  #    01_POOLS folder with the same prefix.
  # ---------------------------------------------------------------------------
  qa_prefix <- sprintf("%d_%s_deterministic_pool_qa", target_year, scenario)

  list(
    burned_pool           = burned_pool,
    unburned_pool         = unburned_pool,
    review_pool           = review_pool,
    scoring_pool          = scoring_pool,
    train_labeled         = train_labeled,
    pools_gpkg            = gpkg_pools_out,
    unburned_hard         = unburned_hard,
    unburned_random       = unburned_random,
    unburned_final_raw    = unburned_final_raw,
    exclusion_buffer      = exclusion_buffer,
    # Gate 1C.2: B4 audit record + negative-pool fingerprint (burnable-domain
    # restricted). The fingerprint INVALIDATES any cached/stale pool whose B4
    # domain / mask / percentile / seed / exclusions / inputs differ.
    b4_audit              = b4_audit,
    neg_pool_fingerprint  = neg_pool_fingerprint,
    audited_internal_gpkg = file.path(
      dirs$`01_POOLS`, paste0(qa_prefix, "_audited_internal.gpkg")
    ),
    qa_summary_csv        = file.path(
      dirs$`01_POOLS`, paste0(qa_prefix, "_summary.csv")
    ),
    qa_transitions_csv    = file.path(
      dirs$`01_POOLS`, paste0(qa_prefix, "_transitions.csv")
    ),
    qa_reasons_csv        = file.path(
      dirs$`01_POOLS`, paste0(qa_prefix, "_reasons.csv")
    )
  )
}
