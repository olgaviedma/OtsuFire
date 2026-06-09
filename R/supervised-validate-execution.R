# ============================================================================
# Gate 1C.4 / 1D.1: semantic + spatial validation for the supervised
# one-year pipeline, plus the cfg-isolation fingerprint helper.
#
# Architecture (Gate 1D.1 — ONE implementation of the 10 checks):
#   .of_validate_supervised_execution_engine()  -- internal ENGINE. Runs all 10
#        checks and RETURNS a STRUCTURED report (one record per check). It does
#        NOT stop() on a failing check; it records the outcome.
#   validate_supervised_execution()             -- PUBLIC wrapper (@export).
#        Calls the engine; in strict = TRUE raises a single aggregated error if
#        any BLOCKING check FAILed, otherwise returns the report.
#   run_supervised_pipeline() (orchestrator)    -- calls the public wrapper in
#        STRICT mode at the same early hook (after cfg resolve, before heavy
#        compute), preserving the Gate 1C.4 fail-fast behaviour. NO second copy
#        of the checks.
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
# features / model).
# ============================================================================

# ----------------------------------------------------------------------------
# Gate 1D.3 — year-validation classification + evidence hierarchy (check 5).
# ----------------------------------------------------------------------------

#' Classify a supervised input as YEAR-SPECIFIC or ATEMPORAL (check-5 lookup).
#'
#' @description
#' EXPLICIT lookup (not ad hoc) driving the robust year check. YEAR-SPECIFIC
#' inputs MUST agree with the resolved cfg target year; ATEMPORAL inputs MUST
#' NOT be required to carry the run year — they are SKIPPED by the year check.
#'
#' YEAR-SPECIFIC: `internal_decisions`, `change_index` (immediate change index),
#'   `delayed_change_index` (delayed change index), `hotspots`,
#'   `reference_burned_map` (the run's reference burned map, when supplied).
#' ATEMPORAL: `topo`, `peninsula_shapefile`, `burnable_mask`, and the CORINE
#'   products (`corine_raster`) — CORINE encodes its OWN epoch/year (e.g. 2012),
#'   not the run target year, so a CORINE filename token must never be mistaken
#'   for the run year.
#'
#' @param name character scalar input name (a key of `cfg$inputs`).
#' @return `"year_specific"`, `"atemporal"`, or `"unknown"` (treated as
#'   atemporal/skipped by the caller — never failed for lacking the run year).
#' @keywords internal
#' @noRd
.of_vse_year_class <- function(name) {
  year_specific <- c("internal_decisions", "change_index",
                     "delayed_change_index", "hotspots",
                     "reference_burned_map")
  atemporal     <- c("topo", "peninsula_shapefile", "burnable_mask",
                     "corine_raster")
  if (name %in% year_specific) "year_specific"
  else if (name %in% atemporal) "atemporal"
  else "unknown"
}

#' Extract an unambiguous standalone 4-digit year token from a filename.
#'
#' A "year-like" token is a standalone `(19|20)\\d\\d` not glued to other digits
#' (so `res90m`, `_res90m_`, a CRS code, or a resolution number is NOT matched),
#' constrained to a plausible run-year magnitude window. Returns the UNIQUE
#' year tokens found (so an ambiguous filename carrying two different years
#' yields >1 element and is treated as non-determinable by the caller).
#'
#' @keywords internal
#' @noRd
.of_vse_filename_year <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(integer(0))
  bn <- basename(path)
  m <- regmatches(bn, gregexpr("(?<![0-9])(19|20)[0-9]{2}(?![0-9])",
                               bn, perl = TRUE))[[1]]
  if (length(m) == 0L) return(integer(0))
  y <- suppressWarnings(as.integer(m))
  y <- y[!is.na(y) & y >= 1900L & y <= 2100L]
  sort(unique(y))
}

#' Resolve a YEAR-SPECIFIC input's year via the Gate 1D.3 evidence hierarchy.
#'
#' @description
#' For ONE year-specific input, attempt to determine its year, IN ORDER:
#'   1. a `year` / `fire_year` (case-insensitive) COLUMN for vector inputs;
#'   2. explicit layer/raster METADATA that records a year (best-effort);
#'   3. an UNAMBIGUOUS standalone 4-digit FILENAME token.
#' The FIRST level that yields year(s) decides the outcome:
#'   - the resolved year(s) MATCH the expected year -> `status = "match"`;
#'   - they DIFFER -> `status = "mismatch"` (blocking FAIL in the caller);
#'   - if NO level yields a year -> `status = "not_verifiable"` (NOT a false
#'     PASS; the caller records NOT_VERIFIABLE, non-blocking by itself).
#' Honesty rule: absence of evidence is recorded as NOT_VERIFIABLE, NEVER as a
#' PASS — the function does not fabricate certainty.
#'
#' NOT IMPLEMENTED (future extension): a fourth hierarchy level — a
#' cfg-registered DECLARED year per input (`spec$declared_year`) — was scoped
#' (Gate 1D.3) but never wired: the config builder
#' [build_supervised_burned_config()] does NOT populate a `declared_year` on any
#' `cfg$inputs[[*]]` spec, so the branch could never fire. It is REMOVED from
#' the active flow (1E.5, 2026-06-09) rather than left as a dormant branch that
#' would imply a capability the package does not have. If a per-input declared
#' year is ever needed, (a) have the builder record `spec$declared_year`, and
#' (b) re-add a level-4 branch here that consults it; until then the hierarchy
#' is column -> metadata -> filename only.
#'
#' @param name character input name.
#' @param path resolved on-disk path (or `NULL`).
#' @param expected_year integer expected (resolved cfg target) year.
#' @return list(status, years, evidence, checked) where `evidence` names the
#'   hierarchy level used (or, for not_verifiable, the levels checked).
#' @keywords internal
#' @noRd
.of_vse_resolve_input_year <- function(name, path, expected_year) {
  ty <- suppressWarnings(as.integer(expected_year)[1L])
  decide <- function(years, evidence) {
    years <- sort(unique(years[!is.na(years)]))
    if (!length(years)) return(NULL)
    list(status = if (ty %in% years) "match" else "mismatch",
         years = years, evidence = evidence)
  }
  is_vector <- name %in% c("internal_decisions", "hotspots",
                           "reference_burned_map")

  # --- level 1: a year / fire_year column (vector inputs only) --------------
  if (is_vector && !is.null(path) && file.exists(path)) {
    attrs <- tryCatch(
      sf::st_drop_geometry(sf::st_read(path, quiet = TRUE)),
      error = function(e) NULL)
    if (!is.null(attrs) && nrow(attrs) > 0L) {
      yr_col <- names(attrs)[tolower(names(attrs)) %in% c("year", "fire_year")]
      if (length(yr_col) >= 1L) {
        yv <- suppressWarnings(as.integer(attrs[[yr_col[1L]]]))
        d <- decide(yv, sprintf("column:%s", yr_col[1L]))
        if (!is.null(d)) return(c(d, list(checked = "column")))
      }
    }
  }

  # --- level 2: explicit raster/layer metadata (best-effort) ----------------
  # terra/GDAL may expose a year via per-layer names or metadata tags. We treat
  # a metadata-derived year as evidence only when it is unambiguous.
  if (!is_vector && !is.null(path) && file.exists(path)) {
    meta_years <- tryCatch({
      r <- terra::rast(path)
      toks <- c(terra::names(r),
                tryCatch(terra::time(r), error = function(e) NULL),
                unlist(tryCatch(terra::metags(r), error = function(e) NULL)))
      toks <- toks[!is.na(toks)]
      yy <- integer(0)
      for (tk in as.character(toks)) {
        yy <- c(yy, .of_vse_filename_year(tk))
      }
      sort(unique(yy))
    }, error = function(e) integer(0))
    if (length(meta_years) == 1L) {
      d <- decide(meta_years, "metadata")
      if (!is.null(d)) return(c(d, list(checked = "metadata")))
    }
  }

  # --- level 3: unambiguous filename token ----------------------------------
  if (!is.null(path)) {
    ftok <- .of_vse_filename_year(path)
    if (length(ftok) == 1L) {
      d <- decide(ftok, "filename")
      if (!is.null(d)) return(c(d, list(checked = "filename")))
    }
  }

  # NOTE (1E.5): a level-4 "cfg-registered declared year" (`spec$declared_year`)
  # is NOT IMPLEMENTED — the builder never records one, so it is intentionally
  # absent from this hierarchy (see the function header). The hierarchy is
  # column -> metadata -> filename only.

  list(status = "not_verifiable", years = integer(0),
       evidence = "no year via column/metadata/filename",
       checked = "column,metadata,filename")
}

# ----------------------------------------------------------------------------
# Report-record constructor. ONE row of the structured report.
# ----------------------------------------------------------------------------
.of_vse_record <- function(check, status, message = NA_character_,
                           evidence = NA_character_, severity = "blocking",
                           verifiable = TRUE) {
  data.frame(
    check      = as.character(check),
    status     = as.character(status),
    message    = as.character(message),
    evidence   = as.character(evidence),
    severity   = as.character(severity),
    verifiable = as.logical(verifiable),
    stringsAsFactors = FALSE
  )
}

#' Internal ENGINE: run the 10 semantic/spatial checks, return a report.
#'
#' @description
#' The SINGLE implementation of the 10 supervised pre-run checks. It NEVER
#' stop()s on a failing check (the foreign-cache check may still `warning()`
#' in the non-blocking branch, preserving the Gate 1C.4 behaviour); each check
#' appends one record to a structured report. Both the public wrapper
#' [validate_supervised_execution()] and the orchestrator consume THIS engine —
#' there is no duplicated check logic.
#'
#' Each record carries: `check`, `status` (PASS / FAIL / NOT_VERIFIABLE /
#' SKIPPED), `message`, `evidence`, `severity` (blocking / warning / info) and
#' `verifiable` (logical). A check that genuinely cannot be evaluated given the
#' inputs is recorded as NOT_VERIFIABLE (verifiable = FALSE) rather than a false
#' PASS; a check that is legitimately inapplicable (e.g. no model supplied for
#' the feature-schema check) is recorded as SKIPPED.
#'
#' @inheritParams validate_supervised_execution
#' @return A `data.frame` with one row per check (the columns above).
#' @keywords internal
#' @noRd
.of_validate_supervised_execution_engine <- function(
    config,
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
    stop(".of_validate_supervised_execution_engine(): 'config' must be a ",
         "resolved otsufire_supervised_burned_config.", call. = FALSE)
  }

  records <- list()
  add <- function(rec) { records[[length(records) + 1L]] <<- rec; invisible(NULL) }

  # A check runner: evaluate `expr` (which on a violation calls stop() exactly
  # like the original Gate 1C.4 logic). If it stop()s, record FAIL with the
  # error message as `message`. If it returns cleanly, the block is expected to
  # have already added its PASS/SKIPPED/NOT_VERIFIABLE record(s) itself.
  run_check <- function(check, severity, expr) {
    tryCatch(
      force(expr),
      error = function(e) {
        add(.of_vse_record(
          check = check, status = "FAIL", message = conditionMessage(e),
          evidence = NA_character_, severity = severity, verifiable = TRUE))
      })
    invisible(NULL)
  }

  ty <- suppressWarnings(as.integer(target_year)[1L])
  if (is.na(ty)) {
    stop(".of_validate_supervised_execution_engine(): could not resolve a ",
         "target_year.", call. = FALSE)
  }

  # --- resolve the input paths the orchestrator will actually consume --------
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
  if (length(raster_inputs)) {
    run_check("crs_rasters", "blocking", {
      run_check("empty_raster", "blocking", {
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
                                "raster has an absent/empty CRS. Every spatial ",
                                "input must carry a valid CRS. cfg field ",
                                "cfg$inputs$%s."), nm, nm), call. = FALSE)
          }
          if (identical(nm, "change_index")) template_r <- r
        }
        add(.of_vse_record("empty_raster", "PASS",
                           "No raster input is empty (>0 rows/cols).",
                           evidence = paste(names(raster_inputs), collapse = ", "),
                           severity = "blocking"))
      })
      add(.of_vse_record("crs_rasters", "PASS",
                         "Every on-disk raster input carries a valid CRS.",
                         evidence = paste(names(raster_inputs), collapse = ", "),
                         severity = "blocking"))
    })
  } else {
    add(.of_vse_record("crs_rasters", "NOT_VERIFIABLE",
                       "No on-disk raster path resolved (in-memory or NULL).",
                       evidence = "no raster paths", severity = "blocking",
                       verifiable = FALSE))
    add(.of_vse_record("empty_raster", "NOT_VERIFIABLE",
                       "No on-disk raster path resolved (in-memory or NULL).",
                       evidence = "no raster paths", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (1) CRS for the vector inputs (internal_decisions, hotspots).
  # ===========================================================================
  decisions_crs <- NULL
  decisions_bbox <- NULL
  if (!is.null(id_path)) {
    run_check("crs_vectors", "blocking", {
      if (!file.exists(id_path)) {
        stop(sprintf(paste0("validate_supervised_execution() ",
                     "[input='internal_decisions']: GPKG does not exist: %s. ",
                     "cfg field cfg$inputs$internal_decisions."),
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
        stop("validate_supervised_execution() [input='internal_decisions']: ",
             "layer has an absent/invalid CRS. cfg field ",
             "cfg$inputs$internal_decisions.", call. = FALSE)
      }
      add(.of_vse_record("crs_vectors", "PASS",
                         "internal_decisions carries a valid CRS.",
                         evidence = sprintf("CRS=%s",
                           tryCatch(decisions_crs$input, error = function(e) NA)),
                         severity = "blocking"))
    })
  } else {
    add(.of_vse_record("crs_vectors", "NOT_VERIFIABLE",
                       "internal_decisions is in-memory or NULL; CRS not on disk.",
                       evidence = "in-memory or NULL", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (2) Spatial overlap: each other spatial input vs the change-index template.
  # ===========================================================================
  if (!is.null(template_r)) {
    run_check("overlap", "blocking", {
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
      add(.of_vse_record("overlap", "PASS",
                         "All spatial inputs overlap the change-index template.",
                         evidence = "change_index template footprint",
                         severity = "blocking"))
    })
  } else {
    add(.of_vse_record("overlap", "NOT_VERIFIABLE",
                       "No change-index template raster resolved; overlap not evaluable.",
                       evidence = "no template", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (4) Mask with ZERO burnable cells — align to template (Gate 1C.1 helper).
  # ===========================================================================
  if (!is.null(mask_path) && !is.null(template_r) && file.exists(mask_path)) {
    run_check("zero_burnable_mask", "blocking", {
      mr <- terra::rast(mask_path)
      invisible(.of_align_mask_to_template(
        mask = mr, template = template_r, allowed_values = c(0, 1),
        binary = TRUE, mask_name = "burnable_mask (cfg$inputs$burnable_mask)"))
      add(.of_vse_record("zero_burnable_mask", "PASS",
                         "Burnable mask aligns to the template with >0 burnable cells.",
                         evidence = basename(mask_path), severity = "blocking"))
    })
  } else {
    add(.of_vse_record("zero_burnable_mask", "NOT_VERIFIABLE",
                       "No burnable mask or no template resolvable.",
                       evidence = "no mask or template", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (5) Wrong year — ROBUST, evidence-hierarchy year validation (Gate 1D.3).
  #
  # CONTRACT (see .of_vse_year_class() / .of_vse_resolve_input_year() below):
  #   * INPUT CLASSIFICATION is explicit, not ad hoc. YEAR-SPECIFIC inputs MUST
  #     agree with the resolved cfg target year; ATEMPORAL inputs (topo,
  #     peninsula, burnable_mask, CORINE products — CORINE encodes its OWN epoch,
  #     not the run target) are SKIPPED (never failed for lacking the run year).
  #   * For each YEAR-SPECIFIC input we resolve a year via a strict HIERARCHY:
  #       (1) a `year` / `fire_year` column (vector inputs);
  #       (2) explicit layer/raster metadata (if present);
  #       (3) an UNAMBIGUOUS standalone 4-digit filename token of the right
  #           magnitude (CORINE epoch / resolution numbers are excluded by
  #           construction — those inputs are ATEMPORAL and never reach here).
  #     (A level-4 cfg-registered declared year was scoped but is NOT
  #      IMPLEMENTED — the builder records no `spec$declared_year`; it was
  #      removed from the active flow in 1E.5. See .of_vse_resolve_input_year().)
  #   * MATCH    -> PASS, evidence names the hierarchy level used.
  #   * MISMATCH -> FAIL (blocking), naming input + found year + expected year.
  #   * NO EVIDENCE at any level -> NOT_VERIFIABLE (verifiable = FALSE), NOT a
  #     false PASS and NOT blocking by itself; the message lists the levels
  #     checked and why none yielded a year.
  # The engine collapses the per-input outcomes into ONE `wrong_year` record:
  # FAIL if any input mismatches (the first mismatch raises, recorded blocking);
  # otherwise PASS/NOT_VERIFIABLE summarising the per-input evidence.
  # ===========================================================================
  year_targets <- list(
    change_index         = raster_inputs[["change_index"]],
    delayed_change_index = raster_inputs[["delayed_change_index"]],
    hotspots             = hs_path,
    internal_decisions   = id_path,
    reference_burned_map = .of_sup_input_path(config, "reference_burned_map")
  )
  year_targets <- year_targets[
    vapply(names(year_targets),
           function(nm) .of_vse_year_class(nm) == "year_specific" &&
                        !is.null(year_targets[[nm]]),
           logical(1))]

  run_check("wrong_year", "blocking", {
    rows <- list()
    for (nm in names(year_targets)) {
      yr <- .of_vse_resolve_input_year(
        nm, year_targets[[nm]], expected_year = ty)
      if (yr$status == "mismatch") {
        stop(sprintf(paste0("validate_supervised_execution() [input='%s']: ",
                            "resolved year %s (evidence: %s) does NOT match the ",
                            "expected cfg$target_year = %d. Wrong-year input; ",
                            "check the cfg route for this input."),
                     nm,
                     paste(yr$years, collapse = "/"), yr$evidence, ty),
             call. = FALSE)
      }
      rows[[nm]] <- yr
    }
    matched <- names(rows)[vapply(rows, function(r) r$status == "match",
                                  logical(1))]
    unverif <- names(rows)[vapply(rows, function(r) r$status == "not_verifiable",
                                  logical(1))]
    if (length(matched) > 0L && length(unverif) == 0L) {
      add(.of_vse_record(
        "wrong_year", "PASS",
        sprintf("Every year-specific input agrees with cfg$target_year = %d.",
                ty),
        evidence = paste(sprintf("%s[%s]", matched,
                                 vapply(rows[matched],
                                        function(r) r$evidence, character(1))),
                         collapse = "; "),
        severity = "blocking"))
    } else if (length(matched) > 0L) {
      add(.of_vse_record(
        "wrong_year", "PASS",
        sprintf(paste0("All year-resolvable inputs agree with cfg$target_year ",
                       "= %d; %d input(s) carried no verifiable year evidence ",
                       "(%s) — recorded as not-verifiable, not a failure."),
                ty, length(unverif), paste(unverif, collapse = ", ")),
        evidence = paste(c(
          sprintf("%s[%s]", matched,
                  vapply(rows[matched], function(r) r$evidence, character(1))),
          sprintf("%s[NOT_VERIFIABLE]", unverif)), collapse = "; "),
        severity = "blocking"))
    } else if (length(unverif) > 0L) {
      add(.of_vse_record(
        "wrong_year", "NOT_VERIFIABLE",
        sprintf(paste0("No year-specific input carried verifiable year ",
                       "evidence (checked column -> metadata -> filename ",
                       "token). Cannot confirm agreement with ",
                       "cfg$target_year = %d; recorded NOT_VERIFIABLE rather ",
                       "than a false PASS. Inputs: %s."),
                ty, paste(unverif, collapse = ", ")),
        evidence = paste(sprintf("%s[%s]", unverif,
                                 vapply(rows[unverif],
                                        function(r) r$evidence, character(1))),
                         collapse = "; "),
        severity = "warning", verifiable = FALSE))
    } else {
      add(.of_vse_record(
        "wrong_year", "SKIPPED",
        "No year-specific input resolved (all inputs absent or atemporal).",
        evidence = "no year-specific inputs", severity = "info",
        verifiable = FALSE))
    }
  })

  # ===========================================================================
  # (6)+(7) Required columns present + correct type in internal_decisions.
  # ===========================================================================
  if (!is.null(id_path) && file.exists(id_path)) {
    run_check("required_columns", "blocking", {
      run_check("incompatible_types", "blocking", {
        id_attr <- tryCatch(
          sf::st_drop_geometry(sf::st_read(id_path, quiet = TRUE,
                                           layer = "internal_decisions")),
          error = function(e)
            stop("validate_supervised_execution() ",
                 "[input='internal_decisions']: cannot read layer ",
                 "'internal_decisions' (", conditionMessage(e), ").",
                 call. = FALSE))
        # (6) required column
        if (!("class_final" %in% names(id_attr))) {
          stop("validate_supervised_execution() [input='internal_decisions']: ",
               "required column 'class_final' is absent (columns: ",
               paste(names(id_attr), collapse = ", "), "). The supervised ",
               "pools stage reads class_final as the decision label.",
               call. = FALSE)
        }
        add(.of_vse_record("required_columns", "PASS",
                           "internal_decisions has the required 'class_final' column.",
                           evidence = paste(names(id_attr), collapse = ", "),
                           severity = "blocking"))
        # (7) incompatible type for the class label
        cf <- id_attr[["class_final"]]
        if (!(is.character(cf) || is.factor(cf))) {
          stop(sprintf(paste0("validate_supervised_execution() ",
                              "[input='internal_decisions']: column ",
                              "'class_final' has incompatible type '%s'; it ",
                              "must be character or factor (the class label is ",
                              "coerced via as.character())."), class(cf)[1L]),
               call. = FALSE)
        }
        # NOTE: the internal_decisions year-column check moved to the unified,
        # evidence-hierarchy check 5 (wrong_year) above (Gate 1D.3); it is no
        # longer duplicated here.
        add(.of_vse_record("incompatible_types", "PASS",
                           "'class_final' is character/factor (valid class label type).",
                           evidence = sprintf("class=%s", class(cf)[1L]),
                           severity = "blocking"))
      })
    })
  } else {
    add(.of_vse_record("required_columns", "NOT_VERIFIABLE",
                       "internal_decisions is in-memory or NULL; columns not on disk.",
                       evidence = "in-memory or NULL", severity = "blocking",
                       verifiable = FALSE))
    add(.of_vse_record("incompatible_types", "NOT_VERIFIABLE",
                       "internal_decisions is in-memory or NULL; type not on disk.",
                       evidence = "in-memory or NULL", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (8) Incompatible feature schema (only when a model/recipe is involved).
  #
  # Gate 1D.2: the saved FINAL-refit RECIPE is the CANONICAL schema source. When
  # a recipe carrying `$cols$feature_cols` is supplied, the check is RECIPE-
  # DRIVEN and aligned with the Gate 1C.3 5-case reconciliation policy: a
  # recipe feature that is MISSING / EXTRA / a TYPE change / a new LEVEL /
  # all-NA is RECOVERABLE by .of_reconcile_scoring_schema() and therefore PASSES
  # (the reconciler creates/drops/coerces/remaps and RECORDS the decision); the
  # check only FAILS on the genuinely-INCOMPATIBLE case the reconciler cannot
  # recover (a feature that coerces to all-NA -> "coerce_type_failed"). The
  # schema is taken FROM the recipe, never re-derived from the scoring year.
  #
  # When only a legacy/alternative schema (recipe$feature_names, or a model's
  # feature_names) is available — i.e. no `$cols$feature_cols` to drive the
  # reconciler — the check keeps the strict "scoring must cover the schema"
  # semantics (any missing expected name FAILs).
  # ===========================================================================
  if ((!is.null(model) || !is.null(recipe)) && !is.null(scoring_feature_names)) {
    run_check("feature_schema", "blocking", {
      expected <- .of_model_expected_features(model, recipe)
      recipe_feature_cols <- tryCatch(recipe$cols$feature_cols,
                                       error = function(e) NULL)
      recipe_driven <- is.character(recipe_feature_cols) &&
        length(recipe_feature_cols) > 0L
      if (recipe_driven) {
        # Dry-run the SAME reconciler the scoring path uses, on a 1-row frame
        # carrying exactly the scoring feature columns (values are placeholders;
        # the reconciliation DECISIONS depend on column presence/type/levels,
        # which we model from the names). This proves recipe -> reconciliation
        # -> exact recipe schema, and classifies each expected feature.
        probe <- as.data.frame(
          stats::setNames(
            rep(list(NA_real_), length(scoring_feature_names)),
            scoring_feature_names),
          stringsAsFactors = FALSE)
        if (nrow(probe) == 0L) probe <- probe[1, , drop = FALSE]
        rec <- .of_reconcile_scoring_schema(probe, recipe)
        # After reconciliation the frame MUST carry exactly the recipe feature
        # schema, in recipe order — the structural proof of check 8 on the main
        # path.
        if (!identical(names(rec$df), recipe_feature_cols)) {
          stop(sprintf(paste0("validate_supervised_execution(): post-",
                              "reconciliation feature schema does NOT match the ",
                              "recipe (expected %d cols in recipe order; got a ",
                              "different set/order). Incompatible feature schema."),
                       length(recipe_feature_cols)), call. = FALSE)
        }
        unrecoverable <- rec$record[rec$record$action == "coerce_type_failed", ,
                                    drop = FALSE]
        if (nrow(unrecoverable) > 0L) {
          stop(sprintf(paste0("validate_supervised_execution(): %d scoring ",
                              "feature(s) are INCOMPATIBLE with the recipe and ",
                              "cannot be reconciled (uncoercible to numeric): ",
                              "%s. Incompatible feature schema."),
                       nrow(unrecoverable),
                       paste(utils::head(unrecoverable$column, 10L),
                             collapse = ", ")),
               call. = FALSE)
        }
        created <- sum(rec$record$action == "create_absent_numeric_NA")
        dropped <- sum(rec$record$action == "drop_extra_column")
        add(.of_vse_record("feature_schema", "PASS",
                           sprintf(paste0("Scoring inputs reconcile to the recipe ",
                                          "schema (%d expected features; %d ",
                                          "created-absent, %d extra dropped)."),
                                   length(recipe_feature_cols), created, dropped),
                           evidence = sprintf("%d recipe features (recipe-driven)",
                                              length(recipe_feature_cols)),
                           severity = "blocking"))
      } else if (length(expected)) {
        missing_cols <- setdiff(expected, scoring_feature_names)
        if (length(missing_cols)) {
          stop(sprintf(paste0("validate_supervised_execution(): scoring inputs ",
                              "are missing %d feature(s) the model/recipe ",
                              "expects: %s. Incompatible feature schema."),
                       length(missing_cols),
                       paste(utils::head(missing_cols, 10L), collapse = ", ")),
               call. = FALSE)
        }
        add(.of_vse_record("feature_schema", "PASS",
                           "Scoring features cover the model/recipe schema.",
                           evidence = sprintf("%d expected features",
                                              length(expected)),
                           severity = "blocking"))
      } else {
        add(.of_vse_record("feature_schema", "NOT_VERIFIABLE",
                           "Model/recipe schema could not be introspected.",
                           evidence = "no expected features extractable",
                           severity = "blocking", verifiable = FALSE))
      }
    })
  } else {
    add(.of_vse_record("feature_schema", "SKIPPED",
                       "No model/recipe + scoring feature names supplied.",
                       evidence = "no model", severity = "blocking",
                       verifiable = FALSE))
  }

  # ===========================================================================
  # (9) Contradictory configuration.
  # ===========================================================================
  training_protocol <- training_protocol %||% "legacy"
  oof_sampling      <- oof_sampling      %||% "capped"
  run_check("contradictory_config", "blocking", {
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
    add(.of_vse_record("contradictory_config", "PASS",
                       "No contradictory cfg setting detected.",
                       evidence = sprintf("protocol=%s; oof=%s; reuse=%s",
                                          training_protocol, oof_sampling,
                                          isTRUE(reuse_upstream)),
                       severity = "blocking"))
  })

  # ===========================================================================
  # (10) Cache belonging to ANOTHER cfg.
  # ===========================================================================
  fp_path <- file.path(
    config$output_routes$pools_dir,
    sprintf("%d_%s_neg_pool_fingerprint.txt", ty, config$scenario))
  if (file.exists(fp_path)) {
    run_check("foreign_cache",
              if (isTRUE(reuse_upstream)) "blocking" else "warning", {
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
          add(.of_vse_record("foreign_cache", "FAIL",
                             paste0(msg, " (non-blocking: pool will be rebuilt)."),
                             evidence = basename(fp_path), severity = "warning"))
        }
      } else {
        add(.of_vse_record("foreign_cache", "PASS",
                           "Negative-pool cache fingerprint matches the current cfg.",
                           evidence = basename(fp_path), severity = "info"))
      }
    })
  } else {
    add(.of_vse_record("foreign_cache", "SKIPPED",
                       "No persisted negative-pool fingerprint cache on disk.",
                       evidence = "no cache", severity = "info",
                       verifiable = FALSE))
  }

  report <- do.call(rbind, records)
  rownames(report) <- NULL
  report
}

#' Semantic + spatial pre-run validation of a resolved supervised configuration
#'
#' @description
#' `validate_supervised_execution()` is the canonical PRE-RUN validator for the
#' one-year supervised pipeline. It accepts an ALREADY-RESOLVED configuration
#' (the object returned by [build_supervised_burned_config()]) and runs ONCE,
#' cheaply, BEFORE any heavy compute (training-pool construction, spatial folds,
#' feature extraction, OOF / final model). It inspects only metadata and
#' geometry (rasters are opened header-only via [terra::rast()]; the only
#' cell-level read is the burnable-mask alignment), so a misconfigured run is
#' caught in well under a second instead of 5-10 min into the chain after
#' partial artefacts are produced.
#'
#' Unlike the previous internal-only contract (which `stop()`ed silently from
#' inside the orchestrator), this function returns a STRUCTURED REPORT — one row
#' per check — so a caller can INSPECT exactly which checks passed, failed, or
#' could not be evaluated. The orchestrator calls this same function in
#' `strict = TRUE` mode at the same early hook, so the report engine and the
#' fail-fast behaviour share ONE implementation of the checks.
#'
#' @details
#' The 10 checks:
#' \enumerate{
#'   \item \strong{CRS absent or invalid} on any spatial input that resolves to
#'     an on-disk path (rasters: change_index / delayed_change_index / topo /
#'     corine_raster / burnable_mask; vectors: internal_decisions / hotspots).
#'   \item \strong{No spatial overlap} between the change-index template and
#'     each other spatial input, tested on bounding geometries reprojected into
#'     the template CRS (Gate 1C.1 extent-intersection logic).
#'   \item \strong{Empty raster} — a raster input with zero rows or columns.
#'   \item \strong{Mask with ZERO burnable cells} — the burnable mask aligned to
#'     the change-index template via [.of_align_mask_to_template()].
#'   \item \strong{Wrong year} (Gate 1D.3, robust evidence hierarchy). Inputs
#'     are CLASSIFIED explicitly (\code{.of_vse_year_class()}): YEAR-SPECIFIC
#'     inputs (internal_decisions, change_index, delayed_change_index, hotspots,
#'     reference_burned_map) MUST agree with the resolved `config$target_year`;
#'     ATEMPORAL inputs (topo, peninsula_shapefile, burnable_mask, and the CORINE
#'     products — CORINE encodes its OWN epoch, not the run year) are SKIPPED and
#'     never failed for lacking the run year. For each year-specific input the
#'     year is resolved via a strict HIERARCHY
#'     (\code{.of_vse_resolve_input_year()}): (1) a `year`/`fire_year` column;
#'     (2) explicit layer/raster metadata; (3) an UNAMBIGUOUS standalone 4-digit
#'     filename token (CORINE epochs / resolution numbers are excluded by
#'     construction). (A fourth level — a cfg-registered declared year
#'     (`spec$declared_year`) — was scoped but is \strong{NOT IMPLEMENTED}: the
#'     builder records no such field, so it was removed from the active flow in
#'     1E.5 and is reserved as a future extension.)
#'     A determinable year that DIFFERS from the expected year FAILs (blocking,
#'     naming input + found + expected); a MATCH PASSes with evidence of the
#'     level used; and an input with NO verifiable year at ANY level is recorded
#'     NOT_VERIFIABLE (verifiable = FALSE) — never a false PASS, and not blocking
#'     by itself.
#'   \item \strong{Required columns absent} — internal_decisions missing its
#'     label column `class_final`.
#'   \item \strong{Incompatible types} — the `class_final` label column is not
#'     character/factor.
#'   \item \strong{Incompatible feature schema} — when a `model`/`recipe` is
#'     supplied with `scoring_feature_names`, the scoring inputs are checked
#'     against the model/recipe schema. Gate 1D.2 makes this RECIPE-DRIVEN and
#'     ENFORCED on the main scoring path: the canonical schema is read from
#'     `recipe$cols$feature_cols` (never re-derived from the scoring year), and
#'     the check is aligned with the Gate 1C.3 5-case reconciliation policy — a
#'     missing / extra / type-changed / unknown-level / all-NA feature is
#'     RECOVERABLE and PASSES (the reconciler heals + records it), while a
#'     genuinely incompatible feature (uncoercible) FAILs (blocking). When only a
#'     legacy schema (`recipe$feature_names` / `model$feature_names`) is
#'     available the strict coverage semantics apply (any uncovered name FAILs).
#'     The check is SKIPPED only when no model/recipe + names are supplied, and
#'     NOT_VERIFIABLE when no schema can be introspected — never a false PASS.
#'     [score_supervised_burned_map()] runs this check (`strict = TRUE`) before
#'     building the scoring matrix, so it is no longer best-effort on the main
#'     path.
#'   \item \strong{Contradictory configuration} — mutually exclusive cfg
#'     settings (e.g. `reuse_upstream = TRUE` with no upstream artefacts;
#'     `training_protocol = "legacy"` with `oof_sampling = "full"`; a
#'     `feature_whitelist_override` that is not a subset of the canonical list).
#'   \item \strong{Cache belonging to ANOTHER cfg} — a persisted negative-pool
#'     fingerprint sidecar (Gate 1C.2) whose checksum disagrees with the
#'     fingerprint recomputed from the CURRENT cfg. Under `reuse_upstream = TRUE`
#'     this is a BLOCKING failure; otherwise it is a WARNING (the pool is rebuilt
#'     this run).
#' }
#'
#' @section Structured report:
#' The return value is a `data.frame` with one row per check and the columns:
#' \describe{
#'   \item{`check`}{Check id (e.g. `"crs_rasters"`, `"overlap"`).}
#'   \item{`status`}{One of `"PASS"`, `"FAIL"`, `"NOT_VERIFIABLE"` (the check
#'     could not be evaluated given the inputs — e.g. an in-memory input with no
#'     on-disk path), or `"SKIPPED"` (the check is legitimately inapplicable —
#'     e.g. no model supplied for the feature-schema check).}
#'   \item{`message`}{Human-readable outcome / failure message.}
#'   \item{`evidence`}{What was inspected (path / CRS / year / column list).}
#'   \item{`severity`}{`"blocking"`, `"warning"`, or `"info"`.}
#'   \item{`verifiable`}{Logical — could this check actually be evaluated for the
#'     given inputs. A NOT_VERIFIABLE / SKIPPED row has `verifiable = FALSE`,
#'     making it impossible to mistake an unevaluable check for a true PASS.}
#' }
#'
#' @section Strict mode:
#' With `strict = TRUE` the function raises a single aggregated error if ANY
#' row with `severity == "blocking"` has `status == "FAIL"`, listing every
#' failing check and its message; otherwise it returns the report invisibly.
#' With `strict = FALSE` (default) it ALWAYS returns the report (visibly) and
#' never errors on a check FAIL — the caller decides what to do. (A non-blocking
#' foreign-cache mismatch still emits a `warning()` in both modes, preserving the
#' Gate 1C.4 behaviour.) The orchestrator invokes `strict = TRUE` so a
#' misconfigured run aborts before heavy compute.
#'
#' @section Distinction from [validate_fire_maps()]:
#' These two validators are complementary and do NOT overlap:
#' \itemize{
#'   \item `validate_supervised_execution()` is a \strong{pre-run} validator: it
#'     checks the supervised CONFIGURATION and its INPUTS (CRS, overlap, year,
#'     decision columns, cfg contradictions, cache provenance) BEFORE any map is
#'     produced. It reads only metadata/geometry and writes nothing.
#'   \item [validate_fire_maps()] is a \strong{post-run} cartographic validator:
#'     it scores PRODUCED burned-area maps against independent reference fire
#'     perimeters (e.g. EFFIS), computing pixel/polygon agreement, omission and
#'     commission. It operates on OUTPUT layers, not on the cfg.
#' }
#'
#' @param config A resolved `otsufire_supervised_burned_config`.
#' @param strict Logical. If `TRUE`, raise a clear aggregated error when any
#'   blocking-severity check FAILs; if `FALSE` (default), return the report
#'   without erroring on a check FAIL.
#' @param target_year Integer. The resolved target year (defaults to
#'   `config$target_year`).
#' @param reuse_upstream Logical. The resolved reuse-upstream toggle. Drives the
#'   contradictory-config and cache-belonging checks.
#' @param feature_whitelist_override Character or `NULL`. The RESOLVED whitelist
#'   override, validated as a subset of the canonical whitelist.
#' @param training_protocol,oof_sampling Character. The resolved protocol toggles.
#' @param model,recipe Optional fitted model / recipe. When supplied, check (8)
#'   asserts the scoring feature schema is compatible.
#' @param scoring_feature_names Character or `NULL`. The available scoring-feature
#'   column names, for check (8). `NULL` skips the schema check.
#' @param data_base,composite_base,result_name Character or `NULL`. The resolved
#'   convention roots, so the validator resolves the SAME input paths the
#'   orchestrator will consume (cfg path %||% convention).
#'
#' @return A `data.frame` — the structured report (see \strong{Structured
#'   report}). Returned invisibly when `strict = TRUE` and all blocking checks
#'   pass; returned visibly when `strict = FALSE`.
#'
#' @seealso [validate_fire_maps()] for post-run cartographic validation of
#'   produced maps against reference perimeters (incl. EFFIS);
#'   [build_supervised_burned_config()] for the configuration this function
#'   validates.
#'
#' @export
validate_supervised_execution <- function(config,
                                          strict = FALSE,
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
  report <- .of_validate_supervised_execution_engine(
    config                     = config,
    target_year                = target_year,
    reuse_upstream             = reuse_upstream,
    feature_whitelist_override = feature_whitelist_override,
    training_protocol          = training_protocol,
    oof_sampling               = oof_sampling,
    model                      = model,
    recipe                     = recipe,
    scoring_feature_names      = scoring_feature_names,
    data_base                  = data_base,
    composite_base             = composite_base,
    result_name                = result_name)

  blocking_fail <- report$severity == "blocking" & report$status == "FAIL"
  if (isTRUE(strict) && any(blocking_fail)) {
    fr <- report[blocking_fail, , drop = FALSE]
    bullet <- paste(sprintf("  - [%s] %s", fr$check, fr$message),
                    collapse = "\n")
    stop(sprintf(paste0("validate_supervised_execution(): %d blocking ",
                        "check(s) FAILED before heavy compute:\n%s"),
                 nrow(fr), bullet), call. = FALSE)
  }
  if (isTRUE(strict)) invisible(report) else report
}

#' Expected feature names of a fitted model / recipe (check 8 helper).
#'
#' Extraction of the feature schema a model/recipe expects, used to assert the
#' scoring inputs cover it. The CANONICAL source is the saved FINAL-refit
#' RECIPE: `recipe$cols$feature_cols` is the exact pre-design-matrix feature set
#' (incl. any `_isNA` companions) the model was trained on, and is what the
#' scoring path reconciles against via [.of_reconcile_scoring_schema()] (Gate
#' 1C.3). This MUST therefore be the schema check 8 enforces on the main path
#' (Gate 1D.2): `recipe of refit model -> final scoring -> exact same schema`.
#'
#' Resolution order (recipe is authoritative; never reconstructed from the
#' scoring year):
#' \enumerate{
#'   \item `recipe$cols$feature_cols` — the canonical refit feature schema.
#'   \item legacy/alternative recipe shapes (`feature_names` / `features` /
#'     `expected_features`) — forward/back-compat only.
#'   \item the fitted xgboost `model$feature_names` — the EXPANDED design-matrix
#'     columns (== `recipe$cols$x_cols`); used only when no recipe schema is
#'     available, since a scoring frame carries pre-expansion feature columns.
#' }
#' Returns `character(0)` only when NONE of these can be introspected (the check
#' then records NOT_VERIFIABLE rather than guessing).
#'
#' @keywords internal
#' @noRd
.of_model_expected_features <- function(model = NULL, recipe = NULL) {
  out <- character(0)
  if (!is.null(recipe)) {
    # CANONICAL: the refit recipe's pre-design-matrix feature columns.
    fc <- tryCatch(recipe$cols$feature_cols, error = function(e) NULL)
    if (is.character(fc) && length(fc)) {
      out <- fc
    } else {
      # Forward/back-compat alternative recipe shapes.
      fn <- recipe$feature_names %||% recipe$features %||%
            recipe$expected_features
      if (is.character(fn)) out <- fn
    }
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
