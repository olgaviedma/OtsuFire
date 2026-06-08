# ============================================================================
# Gate 1C.4: semantic + spatial fail-fast validation for the supervised
# one-year pipeline, plus the cfg-isolation fingerprint helper.
#
# This is the SECOND validation layer. Gate 1B already covers, at the cfg
# boundary and at the top of the orchestrator:
#   - REQUIRED inputs present / not NULL (build_supervised_burned_config);
#   - on-disk PATH specs exist (.of_check_input_file, the orchestrator's
#     stopifnot(file.exists(...)) block);
#   - methodological params consumed only from cfg (single source of truth);
#   - no hidden convention fallbacks on REQUIRED routes (C3 fail-fast).
#
# Gate 1C.4 ADDS the SEMANTIC + SPATIAL layer on top, run ONCE, EARLY and
# CHEAPLY (metadata / geometry-level checks; rasters are opened lazily via
# terra::rast(), which reads only the header, never the full grid). It hooks
# AFTER the cfg is fully resolved but BEFORE any heavy compute (pools / folds /
# features / model). See validate_supervised_execution() below for the 10
# checks and run_supervised_pipeline() for the single invocation site.
# ============================================================================

#' Semantic + spatial fail-fast validation of a resolved supervised run
#'
#' @description
#' Single canonical validator for the one-year supervised pipeline. It runs
#' ONCE, immediately after the configuration is resolved and BEFORE any heavy
#' compute (training-pool construction, spatial folds, feature extraction,
#' OOF / final model). It FAILS FAST with a clear, actionable message — naming
#' the offending input and the relevant cfg field — for any of the following
#' SEMANTIC / SPATIAL problems, so a misconfigured run aborts in <1s instead of
#' 5-10 min into the chain after producing partial artefacts.
#'
#' @details
#' The 10 checks (each `stop()`s with `call. = FALSE` on failure):
#' \enumerate{
#'   \item \strong{CRS absent or invalid} on any spatial input that resolves to
#'     an on-disk path (change_index / delayed_change_index / topo /
#'     corine_raster / burnable_mask rasters; internal_decisions / hotspots
#'     vectors). A raster/vector whose CRS is `NA`/empty is rejected.
#'   \item \strong{No spatial overlap} between the change-index template and
#'     each other spatial input (mask / template / decisions), tested on bounding
#'     geometries reprojected into the template CRS — reuses the Gate 1C.1
#'     extent-intersection logic.
#'   \item \strong{Empty raster} — a raster input with zero rows or columns.
#'   \item \strong{Mask with ZERO burnable cells} — the burnable mask aligned to
#'     the change-index template via the Gate 1C.1 helper
#'     [.of_align_mask_to_template()] (which itself rejects a zero-burnable
#'     mask). Skipped only when no burnable mask is resolvable.
#'   \item \strong{Wrong year} — `cfg$target_year` inconsistent with the inputs:
#'     a `<name>_<year>` token embedded in a resolved input PATH (change_index,
#'     hotspots, ...) that disagrees with `target_year`; OR a hotspots layer that
#'     carries a `year` column with NO row for `target_year`; OR an
#'     internal_decisions layer carrying a `year`/`fire_year` column whose values
#'     never equal `target_year`.
#'   \item \strong{Required columns absent} — the internal_decisions layer is
#'     missing its label column (`class_final`, the column the pools stage reads).
#'   \item \strong{Incompatible types} — the internal_decisions label column
#'     cannot serve its role (e.g. a numeric/list column where a
#'     character/factor class label is required).
#'   \item \strong{Incompatible feature schema} — when a model/recipe is involved
#'     (the `model`/`recipe` arguments, e.g. a reuse-upstream scoring run), the
#'     scoring feature columns do not cover the schema the recipe/model expects.
#'   \item \strong{Contradictory configuration} — mutually exclusive / inconsistent
#'     cfg settings (e.g. `reuse_upstream = TRUE` with no upstream artefacts;
#'     `training_protocol = "legacy"` combined with `oof_sampling = "full"`, which
#'     is only meaningful under `nested_refit`; a feature_whitelist_override that
#'     is not a subset of the canonical whitelist).
#'   \item \strong{Cache belonging to ANOTHER cfg} — a persisted negative-pool
#'     fingerprint sidecar (Gate 1C.2) whose checksum does not match the
#'     fingerprint recomputed from the CURRENT cfg. A stale/foreign cache is NOT
#'     silently reused: under `reuse_upstream = TRUE` this is a hard error
#'     (the run would consume a foreign pool); otherwise it is a warning (the
#'     pool will be rebuilt this run and the sidecar overwritten).
#' }
#'
#' @section Cheapness contract:
#' Rasters are opened with [terra::rast()] (header-only; no cell read) except
#' check (4), which must align the burnable mask to detect a zero-burnable /
#' all-NaN result — that is the cheapest faithful way to catch the silent
#' CRS-mismatch wipe-out the Gate 1C.1 helper guards against, and it is the only
#' input read at cell level. Vectors are read header-first; the decisions/hotspot
#' attribute checks read the attribute table without geometry where the driver
#' allows. The validator performs NO write and creates NO artefact, so a failure
#' leaves the output tree untouched (asserted by the fail-fast tests).
#'
#' @param config A resolved `otsufire_supervised_burned_config`.
#' @param target_year Integer. The resolved target year (defaults to
#'   `config$target_year`).
#' @param reuse_upstream Logical. The resolved reuse-upstream toggle. Drives the
#'   contradictory-config and cache-belonging checks.
#' @param feature_whitelist_override Character or `NULL`. The RESOLVED whitelist
#'   override (cfg value with any boundary override folded in). Validated as a
#'   subset of the canonical whitelist.
#' @param training_protocol,oof_sampling Character. The resolved protocol toggles.
#' @param model,recipe Optional fitted model / recipe (e.g. a reuse-upstream
#'   scoring run). When supplied, check (8) asserts the scoring feature schema is
#'   compatible.
#' @param scoring_feature_names Character or `NULL`. The available scoring-feature
#'   column names, for check (8). `NULL` skips the schema check.
#' @param data_base,composite_base,result_name Character or `NULL`. The resolved
#'   convention roots, so the validator resolves the SAME input paths the
#'   orchestrator will consume (cfg path %||% convention).
#'
#' @return Invisibly, a named list `report` summarising each check's outcome
#'   (`"ok"`, `"skipped:<reason>"`). The function's contract is its side effect:
#'   it either returns cleanly or `stop()`s.
#'
#' @section cfg isolation:
#' The validator is PURE w.r.t. session state: it reads ONLY `config` (and the
#' convention roots derived from it) and the on-disk inputs that cfg points at.
#' It consults NO global/option/tempdir state and writes nothing, so two
#' different cfgs validated back-to-back (or interleaved) cannot contaminate
#' each other. The cache-belonging check (10) is what makes a cfg refuse a pool
#' built for a DIFFERENT cfg: the fingerprint is recomputed from the current cfg
#' via [.of_cfg_neg_pool_fingerprint_preview()] and compared to the persisted
#' sidecar; timestamps NEVER enter the fingerprint (see
#' [.of_cfg_run_fingerprint()]).
#'
#' @keywords internal
#' @noRd
validate_supervised_execution <- function(config,
                                          target_year = config$target_year,
                                          reuse_upstream = FALSE,
                                          feature_whitelist_override = NULL,
                                          training_protocol = "legacy",
                                          oof_sampling = "capped",
                                          model = NULL,
                                          recipe = NULL,
                                          scoring_feature_names = NULL,
                                          data_base = NULL,
                                          composite_base = NULL,
                                          result_name = "Min_Min") {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop("validate_supervised_execution(): 'config' must be a resolved ",
         "otsufire_supervised_burned_config.", call. = FALSE)
  }
  report <- list()
  mark <- function(name, val) { report[[name]] <<- val; invisible(NULL) }

  ty <- suppressWarnings(as.integer(target_year)[1L])
  if (is.na(ty)) {
    stop("validate_supervised_execution(): could not resolve a target_year.",
         call. = FALSE)
  }

  # --- resolve the input paths the orchestrator will actually consume --------
  # cfg path %||% convention. We resolve ONLY paths here (in-memory specs and
  # NULL are tolerated and skipped at the spatial level).
  conv <- .of_supervised_convention_paths(
    data_base = data_base, composite_base = composite_base,
    result_name = result_name %||% "Min_Min", target_year = ty)

  ci_path     <- .of_sup_input_path(config, "change_index")
  dci_path    <- .of_sup_input_path(config, "delayed_change_index") %||%
                   conv$delayed_change_index
  topo_path   <- .of_sup_input_path(config, "topo")          %||% conv$topo
  corine_path <- .of_sup_input_path(config, "corine_raster") %||% conv$corine_raster
  mask_path   <- .of_sup_input_path(config, "burnable_mask") %||% conv$burnable_mask
  id_path     <- .of_sup_input_path(config, "internal_decisions")
  hs_path     <- .of_sup_input_path(config, "hotspots")

  raster_inputs <- list(
    change_index         = ci_path,
    delayed_change_index = dci_path,
    topo                 = topo_path,
    corine_raster        = corine_path,
    burnable_mask        = mask_path
  )
  raster_inputs <- raster_inputs[!vapply(raster_inputs, is.null, logical(1))]

  # ===========================================================================
  # (1)+(3) CRS valid + non-empty raster, per on-disk raster input.
  # ===========================================================================
  open_rast <- function(p, nm) {
    if (!file.exists(p)) {
      stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                          "raster path does not exist: %s"), nm, p),
           call. = FALSE)
    }
    r <- tryCatch(terra::rast(p), error = function(e)
      stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                          "cannot open raster '%s' (%s)."),
                   nm, p, conditionMessage(e)),
           call. = FALSE))
    r
  }
  template_r <- NULL
  for (nm in names(raster_inputs)) {
    r <- open_rast(raster_inputs[[nm]], nm)
    # (3) empty raster
    if (terra::nrow(r) <= 0L || terra::ncol(r) <= 0L) {
      stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                          "empty raster (0 rows/cols). cfg field ",
                          "cfg$inputs$%s."), nm, nm),
           call. = FALSE)
    }
    # (1) CRS absent / invalid
    cr <- tryCatch(terra::crs(r, describe = FALSE), error = function(e) "")
    if (is.na(cr) || !nzchar(cr)) {
      stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                          "raster has an absent/empty CRS. Every spatial input ",
                          "must carry a valid CRS. cfg field cfg$inputs$%s."),
                   nm, nm), call. = FALSE)
    }
    if (identical(nm, "change_index")) template_r <- r
  }
  mark("crs_rasters", if (length(raster_inputs)) "ok" else "skipped:no-raster-paths")
  mark("empty_raster", if (length(raster_inputs)) "ok" else "skipped:no-raster-paths")

  # ===========================================================================
  # (1) CRS for the vector inputs (internal_decisions, hotspots).
  # ===========================================================================
  decisions_crs <- NULL
  decisions_bbox <- NULL
  if (!is.null(id_path)) {
    if (!file.exists(id_path)) {
      stop(sprintf("validate_supervised_execution() [input='internal_decisions']: ",
                   "GPKG does not exist: %s. cfg field cfg$inputs$internal_decisions.",
                   id_path), call. = FALSE)
    }
    id_geom <- tryCatch(
      sf::st_read(id_path, quiet = TRUE, layer = "internal_decisions"),
      error = function(e)
        stop("validate_supervised_execution() [input='internal_decisions']: ",
             "cannot read layer 'internal_decisions' (", conditionMessage(e),
             ").", call. = FALSE))
    decisions_crs <- sf::st_crs(id_geom)
    decisions_bbox <- tryCatch(sf::st_bbox(id_geom), error = function(e) NULL)
    if (is.null(decisions_crs) || is.na(decisions_crs)) {
      stop("validate_supervised_execution() [input='internal_decisions']: layer ",
           "has an absent/invalid CRS. cfg field cfg$inputs$internal_decisions.",
           call. = FALSE)
    }
  }
  mark("crs_vectors", if (!is.null(id_path)) "ok" else "skipped:in-memory-or-null")

  # ===========================================================================
  # (2) Spatial overlap: each other spatial input vs the change-index template.
  #     Bounding-geometry intersection in the template CRS (Gate 1C.1 logic).
  # ===========================================================================
  if (!is.null(template_r)) {
    e_t <- terra::ext(template_r)
    crs_t <- terra::crs(template_r)
    overlap_ext <- function(e_in, crs_in, nm) {
      if (!terra::same.crs(crs_in, crs_t)) {
        e_in <- terra::ext(terra::project(
          terra::as.polygons(e_in, crs = crs_in), crs_t))
      }
      ok <- (min(e_in$xmax, e_t$xmax) > max(e_in$xmin, e_t$xmin)) &&
            (min(e_in$ymax, e_t$ymax) > max(e_in$ymin, e_t$ymin))
      if (!ok) {
        stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                            "does NOT spatially overlap the change-index ",
                            "template. The inputs do not share a footprint ",
                            "(after CRS handling); check the cfg paths."), nm),
             call. = FALSE)
      }
    }
    for (nm in setdiff(names(raster_inputs), "change_index")) {
      r <- terra::rast(raster_inputs[[nm]])
      overlap_ext(terra::ext(r), terra::crs(r), nm)
    }
    if (!is.null(id_path) && !is.null(decisions_bbox)) {
      e_in <- terra::ext(decisions_bbox[["xmin"]], decisions_bbox[["xmax"]],
                         decisions_bbox[["ymin"]], decisions_bbox[["ymax"]])
      overlap_ext(e_in, terra::crs(decisions_crs$wkt), "internal_decisions")
    }
    mark("overlap", "ok")
  } else {
    mark("overlap", "skipped:no-template")
  }

  # ===========================================================================
  # (4) Mask with ZERO burnable cells — align to template (Gate 1C.1 helper),
  #     which itself rejects a zero-burnable / empty / all-NaN aligned mask.
  # ===========================================================================
  if (!is.null(mask_path) && !is.null(template_r) && file.exists(mask_path)) {
    mr <- terra::rast(mask_path)
    # The Gate 1C.1 helper fails fast on: zero burnable cells, empty mask, no
    # overlap, CRS-mismatch wipe-out, out-of-{0,1} values.
    invisible(.of_align_mask_to_template(
      mask = mr, template = template_r, allowed_values = c(0, 1),
      binary = TRUE, mask_name = "burnable_mask (cfg$inputs$burnable_mask)"))
    mark("zero_burnable_mask", "ok")
  } else {
    mark("zero_burnable_mask", "skipped:no-mask-or-template")
  }

  # ===========================================================================
  # (5) Wrong year: filename token / hotspot rows / decisions year column.
  # ===========================================================================
  token_year <- function(p) {
    if (is.null(p)) return(NULL)
    m <- regmatches(basename(p), gregexpr("(?<![0-9])(19|20)[0-9]{2}(?![0-9])",
                                          basename(p), perl = TRUE))[[1]]
    if (length(m) == 0L) return(NULL)
    suppressWarnings(as.integer(m))
  }
  for (nm in names(raster_inputs)) {
    yrs <- token_year(raster_inputs[[nm]])
    # corine_raster / burnable_mask filenames encode the CORINE EPOCH year
    # (e.g. CLC_2012_*), NOT the target year — skip those two.
    if (nm %in% c("corine_raster", "burnable_mask")) next
    if (!is.null(yrs) && length(yrs) == 1L && !(ty %in% yrs)) {
      stop(sprintf(paste0("validate_supervised_execution() [input='%s']: path ",
                          "encodes year %d but cfg$target_year = %d. ",
                          "Wrong-year input; check the cfg path for this ",
                          "route."), nm, yrs[1L], ty),
           call. = FALSE)
    }
  }
  if (!is.null(hs_path) && file.exists(hs_path)) {
    yrs <- token_year(hs_path)
    if (!is.null(yrs) && length(yrs) == 1L && !(ty %in% yrs)) {
      stop(sprintf(paste0("validate_supervised_execution() [input='hotspots']: ",
                          "path encodes year %d but cfg$target_year = %d."),
                   yrs[1L], ty),
           call. = FALSE)
    }
    hs_attr <- tryCatch(sf::st_drop_geometry(sf::st_read(hs_path, quiet = TRUE)),
                        error = function(e) NULL)
    if (!is.null(hs_attr) && "year" %in% names(hs_attr) && nrow(hs_attr) > 0L) {
      hy <- suppressWarnings(as.integer(hs_attr[["year"]]))
      if (all(is.na(hy)) || !(ty %in% hy)) {
        stop(sprintf(paste0("validate_supervised_execution() ",
                            "[input='hotspots']: layer has a 'year' column but ",
                            "NO row for cfg$target_year = %d (years present: ",
                            "%s). Wrong-year hotspots."), ty,
                     paste(sort(unique(hy[!is.na(hy)])), collapse = ", ")),
             call. = FALSE)
      }
    }
    mark("wrong_year", "ok")
  } else {
    mark("wrong_year", "ok")
  }

  # ===========================================================================
  # (6)+(7) Required columns present + correct type in internal_decisions.
  #     The pools stage reads `class_final` (-> class label). Read attributes
  #     only (one row sample), no geometry materialisation needed.
  # ===========================================================================
  if (!is.null(id_path) && file.exists(id_path)) {
    id_attr <- tryCatch(
      sf::st_drop_geometry(sf::st_read(id_path, quiet = TRUE,
                                       layer = "internal_decisions")),
      error = function(e)
        stop("validate_supervised_execution() [input='internal_decisions']: ",
             "cannot read layer 'internal_decisions' (", conditionMessage(e),
             ").", call. = FALSE))
    # (6) required column
    if (!("class_final" %in% names(id_attr))) {
      stop("validate_supervised_execution() [input='internal_decisions']: ",
           "required column 'class_final' is absent (columns: ",
           paste(names(id_attr), collapse = ", "), "). The supervised pools ",
           "stage reads class_final as the decision label.", call. = FALSE)
    }
    # (7) incompatible type for the class label
    cf <- id_attr[["class_final"]]
    if (!(is.character(cf) || is.factor(cf))) {
      stop(sprintf(paste0("validate_supervised_execution() ",
                          "[input='internal_decisions']: column 'class_final' ",
                          "has incompatible type '%s'; it must be character or ",
                          "factor (the class label is coerced via ",
                          "as.character())."), class(cf)[1L]), call. = FALSE)
    }
    # (5b) decisions year column, when present, must include target_year.
    yr_col <- intersect(c("year", "fire_year"), names(id_attr))
    if (length(yr_col) >= 1L && nrow(id_attr) > 0L) {
      dy <- suppressWarnings(as.integer(id_attr[[yr_col[1L]]]))
      if (!all(is.na(dy)) && !(ty %in% dy)) {
        stop(sprintf(paste0("validate_supervised_execution() ",
                            "[input='internal_decisions']: layer has a '%s' ",
                            "column but NO row for cfg$target_year = %d (years ",
                            "present: %s). Wrong-year decisions."),
                     yr_col[1L], ty,
                     paste(sort(unique(dy[!is.na(dy)])), collapse = ", ")),
             call. = FALSE)
      }
    }
    mark("required_columns", "ok")
    mark("incompatible_types", "ok")
  } else {
    mark("required_columns", "skipped:in-memory-or-null")
    mark("incompatible_types", "skipped:in-memory-or-null")
  }

  # ===========================================================================
  # (8) Incompatible feature schema (only when a model/recipe is involved).
  # ===========================================================================
  if ((!is.null(model) || !is.null(recipe)) && !is.null(scoring_feature_names)) {
    expected <- .of_model_expected_features(model, recipe)
    if (length(expected)) {
      missing_cols <- setdiff(expected, scoring_feature_names)
      if (length(missing_cols)) {
        stop(sprintf(paste0("validate_supervised_execution(): scoring inputs ",
                            "are missing %d feature(s) the model/recipe ",
                            "expects: %s. Incompatible feature schema."),
                     length(missing_cols),
                     paste(utils::head(missing_cols, 10L), collapse = ", ")),
             call. = FALSE)
      }
    }
    mark("feature_schema", "ok")
  } else {
    mark("feature_schema", "skipped:no-model")
  }

  # ===========================================================================
  # (9) Contradictory configuration.
  # ===========================================================================
  training_protocol <- training_protocol %||% "legacy"
  oof_sampling      <- oof_sampling      %||% "capped"
  if (identical(training_protocol, "legacy") &&
      identical(oof_sampling, "full")) {
    stop("validate_supervised_execution(): contradictory configuration — ",
         "oof_sampling = 'full' is only meaningful under training_protocol = ",
         "'nested_refit', but training_protocol = 'legacy'. Set oof_sampling = ",
         "'capped' or switch to nested_refit.", call. = FALSE)
  }
  if (!is.null(feature_whitelist_override)) {
    canon <- tryCatch(get(".supervised_feature_cols",
                          envir = asNamespace("OtsuFire")),
                      error = function(e) NULL)
    if (!is.null(canon)) {
      extra <- setdiff(feature_whitelist_override, canon)
      if (length(extra)) {
        stop(sprintf(paste0("validate_supervised_execution(): contradictory ",
                            "configuration — feature_whitelist_override ",
                            "contains name(s) not in the canonical whitelist: ",
                            "%s. The override can only RESTRICT the canonical ",
                            "51-feature list."),
                     paste(utils::head(extra, 10L), collapse = ", ")),
             call. = FALSE)
      }
    }
  }
  if (isTRUE(reuse_upstream)) {
    routes <- config$output_routes
    need <- c(routes$pools_gpkg, routes$features_geometry_gpkg)
    have_folds <- length(list.files(routes$folds_dir,
                          pattern = sprintf("^%d_train_with_folds_\\d+m\\.gpkg$",
                                            ty), full.names = TRUE)) > 0L
    missing_up <- need[!file.exists(need)]
    if (length(missing_up) || !have_folds) {
      stop("validate_supervised_execution(): contradictory configuration — ",
           "reuse_upstream = TRUE but the required upstream artefacts are ",
           "absent (run a baseline first). Missing: ",
           paste(c(missing_up, if (!have_folds) "train_with_folds gpkg"),
                 collapse = "; "), ".", call. = FALSE)
    }
  }
  mark("contradictory_config", "ok")

  # ===========================================================================
  # (10) Cache belonging to ANOTHER cfg: compare the persisted negative-pool
  #      fingerprint sidecar (Gate 1C.2) to the fingerprint recomputed from the
  #      CURRENT cfg. A mismatch means the on-disk pool was built for a DIFFERENT
  #      cfg — never reuse it silently.
  # ===========================================================================
  fp_path <- file.path(
    config$output_routes$pools_dir,
    sprintf("%d_%s_neg_pool_fingerprint.txt", ty, config$scenario))
  if (file.exists(fp_path)) {
    persisted <- .of_read_neg_pool_checksum(fp_path)
    current   <- .of_cfg_neg_pool_fingerprint_preview(config, ty)
    if (!is.null(persisted) && !is.null(current) &&
        !identical(persisted, current$checksum)) {
      msg <- sprintf(
        paste0("validate_supervised_execution(): the on-disk negative-pool ",
               "cache at %s carries fingerprint CHECKSUM=%s, which does NOT ",
               "match the fingerprint recomputed from the CURRENT cfg ",
               "(CHECKSUM=%s). The cached pool belongs to a DIFFERENT cfg."),
        fp_path, persisted, current$checksum)
      if (isTRUE(reuse_upstream)) {
        stop(msg, " reuse_upstream = TRUE would consume this foreign pool; ",
             "refusing. Rebuild the pool for this cfg.", call. = FALSE)
      } else {
        warning(msg, " It will be REBUILT this run (the sidecar is overwritten).",
                call. = FALSE)
      }
    }
    mark("foreign_cache", "ok")
  } else {
    mark("foreign_cache", "skipped:no-cache")
  }

  invisible(list(report = report))
}

#' Expected feature names of a fitted model / recipe (check 8 helper).
#'
#' Best-effort extraction of the feature schema a model/recipe expects, used to
#' assert the scoring inputs cover it. Returns `character(0)` when the schema
#' cannot be introspected (the check then passes vacuously rather than guessing).
#'
#' @keywords internal
#' @noRd
.of_model_expected_features <- function(model = NULL, recipe = NULL) {
  out <- character(0)
  if (!is.null(recipe)) {
    fn <- recipe$feature_names %||% recipe$features %||%
          recipe$expected_features
    if (is.character(fn)) out <- fn
  }
  if (!length(out) && !is.null(model)) {
    fn <- tryCatch({
      if (inherits(model, "xgb.Booster")) {
        model$feature_names
      } else {
        model$feature_names %||% NULL
      }
    }, error = function(e) NULL)
    if (is.character(fn)) out <- fn
  }
  out[!is.na(out) & nzchar(out)]
}

#' Read the CHECKSUM line from a Gate 1C.2 negative-pool fingerprint sidecar.
#'
#' @keywords internal
#' @noRd
.of_read_neg_pool_checksum <- function(fp_path) {
  ln <- tryCatch(readLines(fp_path, warn = FALSE), error = function(e) NULL)
  if (is.null(ln)) return(NULL)
  hit <- grep("^CHECKSUM=", ln, value = TRUE)
  if (!length(hit)) return(NULL)
  sub("^CHECKSUM=", "", hit[1L])
}

#' Preview the negative-pool fingerprint from the cfg ALONE (check 10 helper).
#'
#' Recomputes the cfg-determined SUBSET of the Gate 1C.2 negative-pool
#' fingerprint — the fields that are knowable from the cfg BEFORE the pool is
#' built (policy, seeds, percentile config, exclusions, year/scenario, and the
#' resolved input paths). This is INTENTIONALLY a deterministic function of the
#' cfg only; it folds in NO timestamp and NO runtime-measured quantity, so two
#' equivalent cfgs validated at different wall-clock times yield the IDENTICAL
#' checksum, while a cfg field change (year, seed, decisions path, ...) changes
#' it. The data-dependent fields the full Gate 1C.2 fingerprint adds at build
#' time (B4 percentile value, eligible-cell counts, mask hash) are deliberately
#' OMITTED here: the goal is to detect a cache built for a DIFFERENT cfg, which a
#' divergence in any cfg-determined field already proves.
#'
#' @keywords internal
#' @noRd
.of_cfg_neg_pool_fingerprint_preview <- function(config, target_year) {
  o <- config$options %||% list()
  ci_path <- .of_sup_input_path(config, "change_index") %||% ""
  id_path <- .of_sup_input_path(config, "internal_decisions") %||% ""
  bm_path <- .of_sup_input_path(config, "burnable_mask") %||% ""
  legacy_param_fingerprint_unb_legacy(list(
    neg_pool_policy         = "all_sources",
    b4_domain_decision      = "burnable_restricted",
    b4_burnable_mask_path   = bm_path,
    b4_random_rbr_q         = o$unb_random_rbr_q %||% 0.50,
    b4_random_seed          = o$unb_random_seed %||% 42,
    b4_n_random_cells       = o$unb_n_random_cells %||% 1500,
    b4_random_patch_size    = o$unb_random_patch_size_cells %||% 3,
    b4_exclude_buffer_m     = o$unb_excl_buffer_m %||% 500,
    legacy_random_seed      = o$legacy_random_seed %||% 42L,
    legacy_sample_n         = o$legacy_sample_n %||% 2000,
    legacy_otsu_mode        = o$legacy_otsu_mode %||% "burnable_only",
    target_year             = as.integer(target_year),
    scenario                = config$scenario,
    change_index            = ci_path,
    internal_decisions_path = id_path
  ))
}

#' Reproducible cfg-run fingerprint (cfg isolation; NO timestamps).
#'
#' @description
#' Deterministic identity of a supervised RUN, derived from the cfg fields that
#' DETERMINE its outputs, with the EXPLICIT guarantee that it folds in NO
#' timestamp and NO wall-clock / session state. Two equivalent cfgs fingerprinted
#' at different times produce the IDENTICAL checksum; any consequential cfg field
#' change (target_year, scenario, run_name, output_dir, the resolved input
#' paths, the methodological caps/seeds/protocol, the model_params) produces a
#' DIFFERENT checksum.
#'
#' This is the cfg-isolation anchor: a manifest / cache / output keyed on this
#' fingerprint is tied to the SPECIFIC cfg, so changing a cfg path actually
#' changes the fingerprint (and thus the consumed input), and a fingerprint
#' built for cfg A never matches cfg B. It deliberately reuses the SAME base-R
#' fingerprinter as the fold (Gate 1C.1) and pool (Gate 1C.2) fingerprints, and
#' obeys the same no-timestamps rule.
#'
#' @param config A resolved `otsufire_supervised_burned_config`.
#' @return `list(text, checksum)` — the canonical fingerprint token.
#'
#' @keywords internal
#' @noRd
.of_cfg_run_fingerprint <- function(config) {
  if (!inherits(config, "otsufire_supervised_burned_config")) {
    stop(".of_cfg_run_fingerprint(): 'config' must be a resolved ",
         "otsufire_supervised_burned_config.", call. = FALSE)
  }
  tc <- config$train_control %||% list()
  inp <- function(nm) .of_sup_input_path(config, nm) %||% ""
  caps <- tc$caps %||% list()
  seeds <- tc$seeds %||% list()
  legacy_param_fingerprint_unb_legacy(list(
    scenario                = config$scenario,
    target_year             = as.integer(config$target_year),
    run_name                = config$run_name %||% "",
    output_dir              = config$output_dir %||% "",
    negative_pool_policy    = config$negative_pool_policy %||% "all_sources",
    min_burned_pool_n       = as.integer(config$min_burned_pool_n %||% 5L),
    in_internal_decisions   = inp("internal_decisions"),
    in_change_index         = inp("change_index"),
    in_delayed_change_index = inp("delayed_change_index"),
    in_topo                 = inp("topo"),
    in_corine_raster        = inp("corine_raster"),
    in_burnable_mask        = inp("burnable_mask"),
    in_hotspots             = inp("hotspots"),
    cap_contextual          = caps$contextual %||% NA_real_,
    cap_spectral            = caps$spectral %||% NA_real_,
    cap_random              = caps$random %||% NA_real_,
    cap_otsu                = caps$otsu %||% NA_real_,
    nrounds_max             = tc$nrounds_max %||% NA_integer_,
    early_stop              = tc$early_stop %||% NA_integer_,
    oof_seed_base           = seeds$oof_seed_base %||% NA_integer_,
    final_sampling_seed     = seeds$final_sampling_seed %||% NA_integer_,
    final_seed              = seeds$final_seed %||% NA_integer_,
    val_frac                = tc$val_frac %||% NA_real_,
    group_col               = tc$group_col %||% "",
    impute_numeric          = tc$impute_numeric %||% "",
    impute_factor_missing   = tc$impute_factor_missing %||% "",
    training_protocol       = tc$training_protocol %||% "legacy",
    oof_sampling            = tc$oof_sampling %||% "capped",
    model_params            = config$model_params %||% list()
  ))
}
