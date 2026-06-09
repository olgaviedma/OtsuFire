# Block 7: top-level library() calls removed. Packages resolved via
# DESCRIPTION Imports.
# AS07 (0.5.0): the former top-level `sf::sf_use_s2(FALSE)` was removed. It was
# DEAD in installed-package (library(OtsuFire)) mode -- top-level expressions in
# R/*.R run at build time, not at load time, so the planar setting only took
# effect under pkgload::load_all(). The setting is now scoped (with on.exit
# restore) inside the legacy-stack functions that actually need planar geometry:
# process_otsu_rasters_(), build_unburned_from_legacy_pipeline(),
# build_unburned_from_legacy_decisions(), build_unburned_from_legacy_patches()
# (run_scenarios() already did so at internal-sup-unburned-bloques.R).

sanitize_unb_legacy <- function(x) {
  x <- sf::st_make_valid(x)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  x
}

remove_vector_sidecars_unb_legacy <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "shp") {
    stem <- tools::file_path_sans_ext(path)
    side <- paste0(stem, c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx", ".shp.xml"))
    side <- side[file.exists(side)]
    if (length(side)) file.remove(side)
  } else if (ext == "gpkg" && file.exists(path)) {
    file.remove(path)
  }
  invisible(path)
}

unique_vector_path_unb_legacy <- function(path) {
  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  candidate <- path
  i <- 1L
  while (file.exists(candidate)) {
    candidate <- sprintf("%s_%02d.%s", stem, i, ext)
    i <- i + 1L
  }
  candidate
}

make_shapefile_safe_unb_legacy <- function(x) {
  nm <- names(x)
  gcol <- attr(x, "sf_column")
  idxg <- which(nm == gcol)
  attrs <- if (length(idxg)) nm[-idxg] else nm
  attrs <- make.unique(toupper(substr(attrs, 1, 10)), sep = "_")
  if (length(idxg)) nm[-idxg] <- attrs else nm <- attrs
  names(x) <- nm
  x
}

ensure_area_ha_unb_legacy <- function(x) {
  x$area_ha <- as.numeric(sf::st_area(x)) / 1e4
  x
}

read_vector_unb_legacy <- function(path, layer = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "gpkg") {
    sf::st_read(path, layer = layer %||% "internal_decisions", quiet = TRUE)
  } else {
    sf::st_read(path, quiet = TRUE)
  }
}

# BUG 3 Phase 1a (2026-06-05): removed the dead PROJECT_PATHS.R discovery
# helpers `find_project_paths_file_unb_legacy()` and
# `get_default_fire_mapping_paths_unb_legacy()`. They only fed the
# getwd()-walk + sys.source fallback in build_unburned_from_legacy_pipeline(),
# which is now a fail-fast stop() because the dispatcher always supplies
# data_base / result_name / composite_base explicitly.

# C3 (Gate 1B piece 5, 2026-06-08): the first-existing-path candidate search
# helper `resolve_first_existing_path_unb_legacy()` was REMOVED. It only fed the
# severity-mosaic filesystem glob in build_unburned_from_legacy_pipeline(), which
# now fails fast on a missing explicit `severity_raster_path` instead of guessing
# across MinMin/DOY convention candidates.

# AS02 (0.5.0): build a deterministic, human-readable fingerprint of all
# decision-affecting legacy parameters. The cached intermediate outputs
# (OTSU raster, patches, coverage, decisions) only encode a handful of params
# in their filenames; this fingerprint covers the rest so `reuse_existing`
# cannot silently serve stale results when a non-filename param changed.
# Base R only (no `digest` in Imports): sorted key=value text plus a small
# rolling checksum so the stored token is compact and order-stable.
legacy_param_fingerprint_unb_legacy <- function(params) {
  flat <- vapply(params, function(v) {
    if (is.null(v)) return("NULL")
    v <- unlist(v, use.names = TRUE)
    if (!is.null(names(v)) && any(nzchar(names(v)))) {
      v <- v[order(names(v))]
      paste(sprintf("%s=%s", names(v), format(v, trim = TRUE, scientific = FALSE)),
            collapse = ",")
    } else {
      paste(format(v, trim = TRUE, scientific = FALSE), collapse = ",")
    }
  }, character(1))
  kv <- sprintf("%s=%s", names(params), flat)
  kv <- kv[order(names(params))]
  body <- paste(kv, collapse = "\n")
  # Simple deterministic checksum over the UTF-8 bytes (rolling mod prime).
  # Use double arithmetic (not 32-bit integer) to avoid overflow: every
  # intermediate `chk * 31 + b` stays < 31 * (1e9+7) + 255 ~ 3.1e10, which is
  # exactly representable as a double, so the modulo result is deterministic.
  bytes <- as.numeric(charToRaw(enc2utf8(body)))
  chk <- 0
  for (b in bytes) chk <- (chk * 31 + b) %% 1000000007
  list(text = body, checksum = sprintf("%09d", as.integer(chk)))
}

get_corine_year_unb_legacy <- function(y) {
  if (y >= 1984 && y <= 1999) "1990"
  else if (y <= 2005) "2000"
  else if (y <= 2011) "2006"
  else if (y <= 2017) "2012"
  else "2018"
}

make_corine_reclass_matrix_unb_legacy <- function() {
  matrix(c(
    1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,
    12,2,13,2,14,2,15,2,16,2,17,2,18,2,19,2,20,2,21,2,
    22,3,23,4,24,5,25,6,
    26,7,31,7,32,7,
    27,8,28,8,29,8,
    30,9,33,10,34,11,
    35,12,36,12,37,12,38,12,39,12,40,12,41,12,42,12,43,12,44,12
  ), ncol = 2, byrow = TRUE)
}

source_legacy_unburned_helpers <- function(
  legacy_code_dir = NULL,
  data_base = NULL,
  result_name = "Min_Min",
  target_year = NULL,
  scenario_name = NULL
) {
  # Block 5: the four required helpers (process_otsu_rasters_,
  # polygonize_Otsu, coverage_by_patch_raster, run_scenarios) are now
  # package-internal \u2014 see R/internal-sup-process-otsu-unburned.R and
  # R/internal-sup-unburned-bloques.R. The package namespace already
  # exposes them via lexical scoping, so the legacy runtime sys.source()
  # of external 2_SCRIPTS/para_unburned files is no longer needed and
  # has been removed entirely.
  #
  # We still verify that the four functions resolve in the calling env
  # so the error message stays informative if the package is broken.
  req <- c(
    "process_otsu_rasters_",
    "polygonize_Otsu",
    "coverage_by_patch_raster",
    "run_scenarios"
  )
  ns <- asNamespace("OtsuFire")
  miss_fun <- req[!vapply(req, function(nm)
    exists(nm, envir = ns, mode = "function", inherits = FALSE),
    logical(1))]
  if (length(miss_fun)) {
    stop(
      "Legacy unburned helpers not available in package namespace: ",
      paste(miss_fun, collapse = ", "),
      ". Block 5 expected these to be migrated into R/internal-sup-*.R.",
      call. = FALSE
    )
  }

  invisible(character(0))
}

sanitize_decision_pool_unb_legacy <- function(x, internal, exclude_buffer_m = 0) {
  if (sf::st_crs(x) != sf::st_crs(internal)) {
    x <- sf::st_transform(x, sf::st_crs(internal))
  }

  exclude <- internal
  if (is.finite(exclude_buffer_m) && exclude_buffer_m > 0) {
    exclude <- sf::st_buffer(exclude, exclude_buffer_m)
    exclude <- sanitize_unb_legacy(exclude)
  }

  hits <- lengths(sf::st_intersects(x, exclude)) > 0
  x |>
    dplyr::mutate(intersects_deterministic = hits) |>
    dplyr::filter(!.data$intersects_deterministic)
}

sample_stratified_legacy_unburned <- function(
  x,
  target_n = NULL,
  props = c(drop = 0.70, review = 0.25, keep = 0.05),
  random_seed = 42
) {
  if (!nrow(x)) return(x)
  if (is.null(target_n) || !is.finite(target_n) || target_n <= 0) return(x)

  target_n <- min(as.integer(target_n), nrow(x))
  props <- props[names(props) %in% unique(x$legacy_decision)]
  props <- props[is.finite(props) & props > 0]
  if (!length(props)) {
    set.seed(random_seed)
    return(dplyr::slice_sample(x, n = target_n))
  }
  props <- props / sum(props)

  split_idx <- split(seq_len(nrow(x)), x$legacy_decision)
  wanted <- setNames(rep(0L, length(split_idx)), names(split_idx))
  for (nm in names(props)) {
    if (nm %in% names(split_idx)) {
      wanted[nm] <- floor(target_n * props[[nm]])
    }
  }

  picked <- integer()
  set.seed(random_seed)
  for (nm in names(split_idx)) {
    idx <- split_idx[[nm]]
    n_take <- min(length(idx), wanted[[nm]])
    if (n_take > 0) {
      picked <- c(picked, sample(idx, n_take))
    }
  }

  remaining <- setdiff(seq_len(nrow(x)), picked)
  n_left <- target_n - length(picked)
  if (n_left > 0 && length(remaining) > 0) {
    picked <- c(picked, sample(remaining, min(n_left, length(remaining))))
  }

  x[sort(unique(picked)), , drop = FALSE]
}

build_unburned_from_legacy_decisions <- function(
  legacy_patches_path,
  internal_decisions_path,
  out_gpkg = NULL,
  internal_layer = "internal_decisions",
  use_drop = TRUE,
  # AS12 (0.3.0): defaults for review/keep flipped to FALSE so the
  # legacy pool only carries patches that the training stage actually
  # consumes. Review and keep rows are excluded by
  # `train_final_model_direct()`'s
  # `otsu_unburned_exclude_neg_types = c("otsu_patch_review", "otsu_patch_keep")`,
  # so any sampling budget spent on them was wasted.
  use_review = FALSE,
  use_keep = FALSE,
  drop_max_s_patch = 0.15,
  review_max_s_patch = 0.45,
  keep_max_s_patch = 0.70,
  exclude_buffer_m = 0,
  min_area_ha = 0,
  sample_n = NULL,
  sample_props = c(drop = 0.70, review = 0.25, keep = 0.05),
  random_seed = 42,
  # D4a (2026-06-05): when the Otsu legacy pool is empty after sanitisation,
  # all_sources silently degraded to deterministic_direct semantics (a
  # methodologically different negative pool). That silent degradation is now
  # an ERROR by default. Set TRUE to opt back into the historical warn+degrade
  # behaviour (still writes the `_LEGACY_POOL_EMPTY.txt` audit file).
  allow_empty_otsu_pool = FALSE,
  verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # AS07 (0.5.0): planar geometry for the sf ops below (st_make_valid via
  # sanitize_unb_legacy, st_intersects / st_buffer). Scoped + restored on exit
  # so direct callers (incl. the test suite) get the same planar behaviour the
  # dead top-level `sf_use_s2(FALSE)` only provided under pkgload::load_all().
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  stopifnot(file.exists(legacy_patches_path))
  stopifnot(file.exists(internal_decisions_path))

  legacy <- read_vector_unb_legacy(legacy_patches_path) |>
    sanitize_unb_legacy()

  internal <- read_vector_unb_legacy(internal_decisions_path, layer = internal_layer) |>
    sanitize_unb_legacy()

  required_cols <- c("DECISION", "S_PATCH_PA")
  miss_cols <- setdiff(required_cols, names(legacy))
  if (length(miss_cols)) {
    stop(
      "Legacy patches layer is missing required columns: ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # AS09 (0.5.0): consolidate the fragile area-column trail to a SINGLE
  # canonical column. The decision shapefile arrives with ESRI-truncated
  # area duplicates (`AREA_HA`, and a `make.unique` collision `AREA_HA_1`)
  # carried over from polygonize_Otsu() -> run_scenarios() shapefile
  # round-trips. `ensure_area_ha_unb_legacy()` then (re)computes the
  # authoritative lowercase `area_ha` directly from geometry, so those
  # truncated copies are redundant. Drop them explicitly here so only
  # `area_ha` survives. Byte-identical: the truncated copies never drove a
  # decision (the `min_area_ha` filter, summary, and downstream training all
  # use the geometry-derived `area_ha`; the deny lists already removed the
  # uppercase copies before XGBoost), so removing them earlier changes no
  # produced value -- it only removes dead duplicate columns.
  area_dups <- intersect(c("AREA_HA", "AREA_HA_1", "area_ha"), names(legacy))
  if (length(area_dups)) {
    legacy <- legacy[, setdiff(names(legacy), area_dups), drop = FALSE]
  }

  legacy <- ensure_area_ha_unb_legacy(legacy)
  if (is.finite(min_area_ha) && min_area_ha > 0) {
    legacy <- legacy |>
      dplyr::filter(.data$area_ha >= min_area_ha)
  }

  pieces <- list()

  if (isTRUE(use_drop)) {
    pieces$drop <- legacy |>
      dplyr::filter(.data$DECISION == "drop") |>
      dplyr::filter(is.finite(.data$S_PATCH_PA), .data$S_PATCH_PA <= drop_max_s_patch)
  }
  if (isTRUE(use_review)) {
    pieces$review <- legacy |>
      dplyr::filter(.data$DECISION == "review") |>
      dplyr::filter(is.finite(.data$S_PATCH_PA), .data$S_PATCH_PA <= review_max_s_patch)
  }
  if (isTRUE(use_keep)) {
    pieces$keep <- legacy |>
      dplyr::filter(.data$DECISION == "keep") |>
      dplyr::filter(is.finite(.data$S_PATCH_PA), .data$S_PATCH_PA <= keep_max_s_patch)
  }

  pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0]
  if (!length(pieces)) {
    stop("No legacy candidates remain after decision/score filters.", call. = FALSE)
  }

  combined <- dplyr::bind_rows(
    lapply(names(pieces), function(nm) {
      pieces[[nm]] |>
        dplyr::mutate(
          legacy_decision = nm,
          class = "unburned",
          source = "otsu_patch_residual",
          neg_type = paste0("otsu_patch_", nm)
        )
    })
  )

  combined <- sanitize_decision_pool_unb_legacy(
    x = combined,
    internal = internal,
    exclude_buffer_m = exclude_buffer_m
  ) |>
    ensure_area_ha_unb_legacy()

  # AS06 (0.3.0) / D4a (2026-06-05): if every legacy patch was excluded by the
  # deterministic exclusion buffer, the all_sources mode would silently degrade
  # to deterministic_direct semantics (a methodologically different negative
  # pool: no Otsu current-year patches). That silent degradation is now an
  # ERROR by default -- the user gets a clear message instead of an unflagged
  # change of negative-class composition. Set `allow_empty_otsu_pool = TRUE` to
  # opt back into the historical warn+degrade behaviour, which still writes the
  # `_LEGACY_POOL_EMPTY.txt` audit file so post-hoc analyses can flag affected
  # years. Normal (non-degenerate) years are byte-identical either way.
  if (nrow(combined) == 0L) {
    if (!isTRUE(allow_empty_otsu_pool)) {
      stop(
        paste0(
          "Otsu legacy pool is EMPTY after sanitisation (exclude_buffer_m = ",
          as.character(exclude_buffer_m), "): every legacy patch was removed ",
          "by the deterministic exclusion buffer. all_sources would silently ",
          "degrade to deterministic_direct semantics (no Otsu current-year ",
          "patches in the negative pool), which changes the model's negative ",
          "class. Aborting. To proceed anyway with the historical ",
          "warn+degrade behaviour, set allow_empty_otsu_pool = TRUE (config ",
          "option `options$allow_empty_otsu_pool`). ",
          "Legacy patches path: ", legacy_patches_path, ". ",
          "Internal decisions path: ", internal_decisions_path, "."
        ),
        call. = FALSE
      )
    }
    warning(
      paste0("Otsu legacy pool empty after sanitisation; all_sources ",
             "degrading to deterministic_direct semantics for this ",
             "year/scenario (allow_empty_otsu_pool = TRUE)."),
      call. = FALSE
    )
    if (!is.null(out_gpkg) && nzchar(out_gpkg)) {
      out_dir_audit <- dirname(out_gpkg)
      tryCatch({
        dir.create(out_dir_audit, recursive = TRUE, showWarnings = FALSE)
        writeLines(
          c(
            "OtsuFire 0.3.0 -- Otsu legacy pool empty after sanitisation (AS06).",
            sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            sprintf("Legacy patches path: %s", legacy_patches_path),
            sprintf("Internal decisions path: %s", internal_decisions_path),
            sprintf("exclude_buffer_m = %s", as.character(exclude_buffer_m)),
            "",
            "All_sources mode degraded to deterministic_direct ",
            "semantics for this year/scenario (allow_empty_otsu_pool = TRUE)."
          ),
          file.path(out_dir_audit, "_LEGACY_POOL_EMPTY.txt")
        )
      }, error = function(e) {
        warning("Could not write _LEGACY_POOL_EMPTY.txt: ",
                conditionMessage(e), call. = FALSE)
      })
    }
  }

  sampled <- sample_stratified_legacy_unburned(
    x = combined,
    target_n = sample_n,
    props = sample_props,
    random_seed = random_seed
  )

  summary_tbl <- combined |>
    sf::st_drop_geometry() |>
    dplyr::count(.data$legacy_decision, name = "n_unburned_pool") |>
    dplyr::arrange(.data$legacy_decision)

  summary_meta <- data.frame(
    legacy_path = legacy_patches_path,
    internal_path = internal_decisions_path,
    drop_max_s_patch = drop_max_s_patch,
    review_max_s_patch = review_max_s_patch,
    keep_max_s_patch = keep_max_s_patch,
    exclude_buffer_m = exclude_buffer_m,
    min_area_ha = min_area_ha,
    n_unburned_pool = nrow(combined),
    n_unburned_sampled = nrow(sampled),
    stringsAsFactors = FALSE
  )

  if (!is.null(out_gpkg) && nzchar(out_gpkg)) {
    dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
    write_ok <- tryCatch({
      remove_vector_sidecars_unb_legacy(out_gpkg)
      sf::st_write(combined, out_gpkg, layer = "legacy_unburned_pool", delete_layer = TRUE, quiet = TRUE)
      sf::st_write(sampled, out_gpkg, layer = "legacy_unburned_sampled", delete_layer = TRUE, quiet = TRUE)
      TRUE
    }, error = function(e) {
      message("GPKG write failed, falling back to shapefiles: ", conditionMessage(e))
      FALSE
    })

    if (!isTRUE(write_ok)) {
      stem <- file.path(dirname(out_gpkg), tools::file_path_sans_ext(basename(out_gpkg)))
      pool_shp <- paste0(stem, "_pool.shp")
      sampled_shp <- paste0(stem, "_sampled.shp")
      remove_vector_sidecars_unb_legacy(pool_shp)
      remove_vector_sidecars_unb_legacy(sampled_shp)
      if (file.exists(pool_shp)) pool_shp <- unique_vector_path_unb_legacy(pool_shp)
      if (file.exists(sampled_shp)) sampled_shp <- unique_vector_path_unb_legacy(sampled_shp)
      combined_shp <- make_shapefile_safe_unb_legacy(combined)
      sampled_shp_sf <- make_shapefile_safe_unb_legacy(sampled)
      sf::st_write(combined_shp, pool_shp, quiet = TRUE)
      tryCatch(
        sf::st_write(sampled_shp_sf, sampled_shp, quiet = TRUE),
        error = function(e) message("Sampled shapefile write failed; continuing with pool only: ", conditionMessage(e))
      )
      out_gpkg <- pool_shp
    }

    utils::write.csv(
      summary_tbl,
      file.path(dirname(out_gpkg), paste0(tools::file_path_sans_ext(basename(out_gpkg)), "_by_decision.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      summary_meta,
      file.path(dirname(out_gpkg), paste0(tools::file_path_sans_ext(basename(out_gpkg)), "_summary.csv")),
      row.names = FALSE
    )
  }

  msg("Legacy residual unburned pool: %d", nrow(combined))
  msg("Legacy residual sampled subset: %d", nrow(sampled))

  invisible(list(
    summary = summary_tbl,
    meta = summary_meta,
    legacy_unburned_pool = combined,
    legacy_unburned_sampled = sampled,
    out_gpkg = out_gpkg
  ))
}

build_unburned_from_legacy_patches <- function(
  legacy_patches_path,
  internal_decisions_path,
  out_gpkg = NULL,
  internal_layer = "internal_decisions",
  legacy_decisions = c("drop"),
  max_s_patch = 0.15,
  exclude_buffer_m = 0,
  min_area_ha = 0,
  sample_n = NULL,
  random_seed = 42,
  verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # AS07 (0.5.0): planar geometry for the sf ops below (st_make_valid via
  # sanitize_unb_legacy, st_intersects / st_buffer). Scoped + restored on exit
  # so direct callers (incl. the test suite) get the same planar behaviour the
  # dead top-level `sf_use_s2(FALSE)` only provided under pkgload::load_all().
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  stopifnot(file.exists(legacy_patches_path))
  stopifnot(file.exists(internal_decisions_path))

  legacy <- read_vector_unb_legacy(legacy_patches_path) |>
    sanitize_unb_legacy()

  internal <- read_vector_unb_legacy(internal_decisions_path, layer = internal_layer) |>
    sanitize_unb_legacy()

  required_cols <- c("DECISION", "S_PATCH_PA")
  miss_cols <- setdiff(required_cols, names(legacy))
  if (length(miss_cols)) {
    stop(
      "Legacy patches layer is missing required columns: ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  legacy_decisions <- unique(as.character(legacy_decisions))
  legacy <- legacy |>
    dplyr::filter(.data$DECISION %in% legacy_decisions)

  if (!is.null(max_s_patch) && "S_PATCH_PA" %in% names(legacy)) {
    legacy <- legacy |>
      dplyr::filter(is.finite(.data$S_PATCH_PA), .data$S_PATCH_PA <= max_s_patch)
  }

  if (sf::st_crs(legacy) != sf::st_crs(internal)) {
    legacy <- sf::st_transform(legacy, sf::st_crs(internal))
  }

  legacy <- ensure_area_ha_unb_legacy(legacy)
  if (is.finite(min_area_ha) && min_area_ha > 0) {
    legacy <- legacy |>
      dplyr::filter(.data$area_ha >= min_area_ha)
  }

  exclude <- internal
  if (is.finite(exclude_buffer_m) && exclude_buffer_m > 0) {
    exclude <- sf::st_buffer(exclude, exclude_buffer_m)
    exclude <- sanitize_unb_legacy(exclude)
  }

  hits <- lengths(sf::st_intersects(legacy, exclude)) > 0

  legacy_unburned <- legacy |>
    dplyr::mutate(
      intersects_deterministic = hits,
      class = "unburned",
      source = "otsu_patch_residual",
      neg_type = paste0("otsu_patch_", tolower(.data$DECISION))
    ) |>
    dplyr::filter(!.data$intersects_deterministic)

  legacy_unburned <- ensure_area_ha_unb_legacy(legacy_unburned)

  sampled <- legacy_unburned
  if (is.numeric(sample_n) && length(sample_n) == 1L && is.finite(sample_n) && sample_n > 0) {
    sample_n <- min(as.integer(sample_n), nrow(legacy_unburned))
    set.seed(random_seed)
    sampled <- legacy_unburned |>
      dplyr::slice_sample(n = sample_n)
  }

  summary_tbl <- data.frame(
    legacy_path = legacy_patches_path,
    internal_path = internal_decisions_path,
    legacy_decisions = paste(sort(unique(legacy_decisions)), collapse = ","),
    max_s_patch = if (is.null(max_s_patch)) NA_real_ else as.numeric(max_s_patch),
    exclude_buffer_m = as.numeric(exclude_buffer_m),
    min_area_ha = as.numeric(min_area_ha),
    n_legacy_filtered = nrow(legacy),
    n_intersecting_deterministic = sum(hits),
    n_unburned_pool = nrow(legacy_unburned),
    n_unburned_sampled = nrow(sampled),
    stringsAsFactors = FALSE
  )

  if (!is.null(out_gpkg) && nzchar(out_gpkg)) {
    dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)

    sf::st_write(legacy_unburned, out_gpkg, layer = "legacy_unburned_pool", delete_layer = TRUE, quiet = TRUE)
    sf::st_write(sampled, out_gpkg, layer = "legacy_unburned_sampled", delete_layer = TRUE, quiet = TRUE)
    sf::st_write(
      internal |>
        dplyr::select(dplyr::any_of(c("poly_id", "source_poly_id", "class_final")), geometry),
      out_gpkg,
      layer = "deterministic_reference",
      delete_layer = TRUE,
      quiet = TRUE
    )

    utils::write.csv(
      summary_tbl,
      file.path(dirname(out_gpkg), paste0(tools::file_path_sans_ext(basename(out_gpkg)), "_summary.csv")),
      row.names = FALSE
    )
  }

  msg("Legacy filtered candidates: %d", nrow(legacy))
  msg("Legacy candidates intersecting deterministic: %d", sum(hits))
  msg("Legacy residual unburned pool: %d", nrow(legacy_unburned))
  if (nrow(sampled) != nrow(legacy_unburned)) {
    msg("Legacy residual sampled subset: %d", nrow(sampled))
  }

  invisible(list(
    summary = summary_tbl,
    legacy_filtered = legacy,
    legacy_unburned_pool = legacy_unburned,
    legacy_unburned_sampled = sampled,
    out_gpkg = out_gpkg
  ))
}

build_unburned_from_legacy_pipeline <- function(
  target_year,
  scenario_name,
  data_base = NULL,
  result_name = NULL,
  composite_base = NULL,
  severity_raster_path = NULL,
  internal_decisions_path = NULL,
  legacy_code_dir = NULL,
  otsu_mode = c("burnable_only", "corine"),
  # AS01 (0.3.0): NULL by default. Plumbed from config$tool_paths via the
  # dispatcher. The downstream `process_otsu_rasters_()` call asserts
  # non-NULL. Provide explicitly when calling this function directly.
  python_exe = NULL,
  gdal_polygonize_script = NULL,
  gdalwarp_path = NULL,
  ogr2ogr_exe = NULL,
  otsu_threshold = 0,
  reference_otsu_threshold = 100,
  min_otsu_threshold_value = 0,
  min_pixels = 8,
  target_epsg = 3035,
  buffers_m = 90,
  core_thr = 0.60,
  alpha_boost = 0.25,
  min_base_boost = 0.35,
  dist_power = 1,
  keep_hi = 0.45,
  drop_lo = 0.15,
  dist_mode = "centroid",
  near_mode = "centroid",
  use_drop = TRUE,
  use_review = TRUE,
  use_keep = TRUE,
  drop_max_s_patch = 0.15,
  review_max_s_patch = 0.45,
  keep_max_s_patch = 0.70,
  sample_n = 2000,
  sample_props = c(drop = 0.70, review = 0.25, keep = 0.05),
  exclude_buffer_m = 0,
  min_area_ha = 0,
  random_seed = 42,
  reuse_existing = TRUE,
  write_unburned = TRUE,
  out_root_dir = NULL,
  # §N+25 (2026-06-05): burnable mask, CORINE raster and peninsula border are
  # now wired RUN inputs. The caller (supervised-pools.R) threads
  # config$inputs$burnable_mask / corine_raster / peninsula_shapefile here.
  # When NULL each falls back to EXACTLY the historical convention path below,
  # so passing nothing is byte-identical.
  burnable_mask_path = NULL,
  corine_raster_path = NULL,
  peninsula_shapefile = NULL,
  # D4a (2026-06-05): forwarded to build_unburned_from_legacy_decisions().
  # FALSE (default) errors when the Otsu legacy pool is empty after
  # sanitisation instead of silently degrading to deterministic_direct
  # semantics; TRUE keeps the historical warn+degrade behaviour.
  allow_empty_otsu_pool = FALSE,
  verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  otsu_mode <- match.arg(otsu_mode)

  # AS07 (0.5.0): planar geometry for the whole legacy unburned chain
  # (process_otsu_rasters_, polygonize_Otsu, coverage_by_patch_raster,
  # run_scenarios, and the build_unburned_from_legacy_decisions sf ops:
  # st_make_valid / st_intersects / st_buffer). The historical top-level
  # `sf::sf_use_s2(FALSE)` was dead in installed-package mode; scope it here
  # (restored on exit) so behaviour matches the pkgload::load_all path in both
  # modes without leaking the setting into the caller's session.
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  # BUG 3 Phase 1a (2026-06-05): the production dispatcher ALWAYS supplies
  # data_base / result_name / composite_base via the engine bindings, so the
  # former PROJECT_PATHS.R / getwd()-walk fallback was dead. Fail fast instead
  # of silently resolving paths from the working directory.
  if (is.null(data_base) || is.null(result_name) || is.null(composite_base)) {
    stop(
      "build_unburned_from_legacy_pipeline() requires 'data_base', ",
      "'result_name' and 'composite_base' to be supplied explicitly. ",
      "(The legacy PROJECT_PATHS.R / getwd() fallback was removed.)",
      call. = FALSE
    )
  }

  if (!dir.exists(data_base)) {
    stop("Could not resolve 'data_base': directory does not exist: ",
         data_base, call. = FALSE)
  }

  helper_files <- source_legacy_unburned_helpers(
    legacy_code_dir = legacy_code_dir,
    data_base = data_base,
    result_name = result_name,
    target_year = target_year,
    scenario_name = scenario_name
  )

  corine_year <- get_corine_year_unb_legacy(target_year)
  reclass_matrix <- make_corine_reclass_matrix_unb_legacy()
  # C3 (Gate 1B piece 5, 2026-06-08): the severity / change-index mosaic is a
  # REQUIRED supervised input. It MUST be supplied explicitly — in production
  # the pools stage threads cfg$inputs$change_index here as `severity_raster_path`.
  # The former behaviour silently GLOBBED the filesystem across six convention
  # candidates (MinMin/DOY under composite_base or data_base/Imagery) and used
  # the first that happened to exist, which could pick a stale or wrong-product
  # raster behind the cfg's back. That candidate search is removed: if no
  # explicit path is given we FAIL FAST naming the cfg field to set, BEFORE any
  # heavy raster compute.
  if (is.null(severity_raster_path) || !nzchar(severity_raster_path)) {
    stop(
      "build_unburned_from_legacy_pipeline() requires 'severity_raster_path' ",
      "(the change-index / severity mosaic). In a supervised run this is the ",
      "cfg$inputs$change_index route, threaded by the pools stage; set it via ",
      "build_supervised_burned_config(change_index = ...). There is no ",
      "convention-based filesystem search for this raster.",
      call. = FALSE
    )
  }
  one_year_tif <- severity_raster_path
  one_year_tif_source <- "explicit (from config)"
  if (is.null(burnable_mask_path) || !nzchar(burnable_mask_path)) {
    burnable_mask_path <- file.path(
      data_base, "Corine_Masks",
      paste0("burneable_mask_binary_corine_", corine_year, "_ETRS89.tif")
    )
  }
  if (is.null(corine_raster_path) || !nzchar(corine_raster_path)) {
    corine_raster_path <- file.path(
      data_base, "Corine_Masks",
      paste0("CLC_", corine_year, "_peninsula.tif")
    )
  }
  if (is.null(peninsula_shapefile) || !nzchar(peninsula_shapefile)) {
    peninsula_shapefile <- file.path(data_base, "Borders", "Iberian_peninsula.shp")
  }
  # The deterministic decisions GPKG is the user-supplied path, period. No
  # convention reconstruction. The caller (orchestrator) always threads the
  # user path down; require it here.
  if (is.null(internal_decisions_path) || !nzchar(internal_decisions_path)) {
    stop(
      "build_unburned_from_legacy_pipeline() requires ",
      "'internal_decisions_path' (the deterministic decisions .gpkg). ",
      "There is no convention-based fallback.",
      call. = FALSE
    )
  }

  if (!file.exists(one_year_tif)) {
    stop(
      "Severity / change-index mosaic for year ", target_year,
      " does not exist on disk: ", one_year_tif, ". ",
      "It was supplied explicitly via 'severity_raster_path' ",
      "(cfg$inputs$change_index). Check that path.",
      call. = FALSE
    )
  }
  one_year_tif <- normalizePath(one_year_tif, winslash = "/", mustWork = TRUE)
  msg("Legacy unburned severity raster [%s]: %s", one_year_tif_source, one_year_tif)
  stopifnot(file.exists(burnable_mask_path))
  if (identical(otsu_mode, "corine")) {
    stopifnot(file.exists(corine_raster_path))
    stopifnot(file.exists(peninsula_shapefile))
  }
  stopifnot(file.exists(internal_decisions_path))

  # C3 (Gate 1B piece 5, 2026-06-08): the legacy-unburned output ROOT must be
  # resolved from the cfg output routes, never invented from a filename
  # convention. In a supervised run the pools stage derives
  # `legacy_unb_root_dir` from config$output_routes$base and threads it here as
  # `out_root_dir`. The former fallback silently fabricated a
  # data_base/Results/.../SUPERVISED/.../_LEGACY_UNBURNED directory (a different
  # output tree from the cfg-configured one), so cached stages and the unburned
  # pool could be written/read under a path the caller never asked for. Require
  # `out_root_dir`; fail fast BEFORE any heavy compute if it is absent.
  if (is.null(out_root_dir) || !is.character(out_root_dir) ||
      length(out_root_dir) != 1L || !nzchar(out_root_dir)) {
    stop(
      "build_unburned_from_legacy_pipeline() requires 'out_root_dir' (the ",
      "legacy-unburned output root). In a supervised run this is derived from ",
      "config$output_routes$base by the pools stage. There is no ",
      "convention-based output directory fallback; pass 'out_root_dir' ",
      "explicitly or supply config$output_routes.",
      call. = FALSE
    )
  }
  dir.create(out_root_dir, recursive = TRUE, showWarnings = FALSE)

  dirs <- list(
    otsu = file.path(out_root_dir, "01_OTSU"),
    patches = file.path(out_root_dir, "02_PATCHES"),
    coverage = file.path(out_root_dir, "03_COVERAGE"),
    decisions = file.path(out_root_dir, "04_PATCH_DECISIONS"),
    unburned = file.path(out_root_dir, "05_UNBURNED")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  otsu_file_stem <- switch(
    otsu_mode,
    burnable_only = sprintf("BA_%s_ge%d", target_year, otsu_threshold),
    corine = sprintf("BA_%s_otsu_CORI_ge%d", target_year, otsu_threshold)
  )
  ref_file_stem <- switch(
    otsu_mode,
    burnable_only = sprintf("BA_%s_ge%d", target_year, reference_otsu_threshold),
    corine = sprintf("BA_%s_otsu_CORI_ge%d", target_year, reference_otsu_threshold)
  )
  otsu_raster_path <- file.path(dirs$otsu, sprintf("%s_binary.tif", otsu_file_stem))
  ref_raster_path <- file.path(dirs$otsu, sprintf("%s_binary.tif", ref_file_stem))
  patch_path <- file.path(
    dirs$patches,
    sprintf("%s_patches.shp", otsu_file_stem)
  )
  coverage_path <- file.path(
    dirs$coverage,
    sprintf("%s_ref%d_binary_cov.shp", otsu_file_stem, reference_otsu_threshold)
  )
  decision_slug <- sprintf("cv100_b%dm_t%02d", as.integer(buffers_m), round(core_thr * 100))
  decision_path <- file.path(
    dirs$decisions,
    sprintf("sc_%s_patches.shp", decision_slug)
  )
  unburned_out_gpkg <- file.path(
    dirs$unburned,
    sprintf("%d_%s_legacy_unburned.gpkg", target_year, scenario_name)
  )

  # AS02 (0.5.0): honest reuse. The cached stage filenames only encode a few
  # of the ~25 governing parameters (target_year, otsu_mode, otsu_threshold,
  # reference_otsu_threshold, buffers_m, core_thr). Changing any OTHER param
  # while keeping reuse_existing = TRUE previously served STALE cached outputs
  # silently. We now write a parameter fingerprint next to the cache and only
  # reuse when it MATCHES; on mismatch (or when no fingerprint exists, e.g. a
  # pre-AS02 cache) we recompute rather than reuse. For a normal production
  # re-run with unchanged params the fingerprint matches and the reuse path is
  # byte-identical to before.
  fp <- legacy_param_fingerprint_unb_legacy(list(
    target_year              = target_year,
    scenario_name            = scenario_name,
    otsu_mode                = otsu_mode,
    otsu_threshold           = otsu_threshold,
    reference_otsu_threshold = reference_otsu_threshold,
    min_otsu_threshold_value = min_otsu_threshold_value,
    min_pixels               = min_pixels,
    target_epsg              = target_epsg,
    buffers_m                = buffers_m,
    core_thr                 = core_thr,
    alpha_boost              = alpha_boost,
    min_base_boost           = min_base_boost,
    dist_power               = dist_power,
    keep_hi                  = keep_hi,
    drop_lo                  = drop_lo,
    dist_mode                = dist_mode,
    near_mode                = near_mode,
    use_drop                 = use_drop,
    use_review               = use_review,
    use_keep                 = use_keep,
    drop_max_s_patch         = drop_max_s_patch,
    review_max_s_patch       = review_max_s_patch,
    keep_max_s_patch         = keep_max_s_patch,
    sample_n                 = sample_n,
    sample_props             = sample_props,
    exclude_buffer_m         = exclude_buffer_m,
    min_area_ha              = min_area_ha,
    random_seed              = random_seed,
    one_year_tif             = one_year_tif,
    burnable_mask_path       = burnable_mask_path,
    corine_raster_path       = corine_raster_path,
    peninsula_shapefile      = peninsula_shapefile,
    internal_decisions_path  = internal_decisions_path
  ))
  fingerprint_path <- file.path(out_root_dir, "_LEGACY_PARAM_FINGERPRINT.txt")
  fingerprint_token <- c(
    "# OtsuFire legacy unburned cache fingerprint (AS02).",
    "# reuse_existing only honoured when this file matches the current call.",
    sprintf("CHECKSUM=%s", fp$checksum),
    "---",
    fp$text
  )
  # AS02 fix: `fingerprint_token` is a character vector whose LAST element
  # (`fp$text`) is itself a multi-line string with embedded "\n". `writeLines`
  # expands those embedded newlines into separate on-disk lines, so a naive
  # `identical(readLines(...), fingerprint_token)` always returned FALSE (the
  # round-trip changes the element count) and the cache was NEVER reused.
  # Normalise BOTH sides to a single "\n"-joined string before comparing so an
  # unchanged call MATCHES and a changed call does NOT.
  cache_fingerprint_matches <- file.exists(fingerprint_path) &&
    identical(
      paste(readLines(fingerprint_path, warn = FALSE), collapse = "\n"),
      paste(fingerprint_token, collapse = "\n")
    )
  reuse_ok <- isTRUE(reuse_existing) && isTRUE(cache_fingerprint_matches)
  if (isTRUE(reuse_existing) && !isTRUE(cache_fingerprint_matches)) {
    msg(paste0("Legacy cache fingerprint %s; recomputing all stages instead ",
               "of reusing (AS02: params changed or no fingerprint present)."),
        if (file.exists(fingerprint_path)) "MISMATCH" else "ABSENT")
  }

  if (!isTRUE(reuse_ok) || !file.exists(otsu_raster_path) || !file.exists(ref_raster_path)) {
    msg("STEP 1 - Legacy OTSU raster [%s] | candidate ge%d | reference ge%d", otsu_mode, otsu_threshold, reference_otsu_threshold)
    process_otsu_rasters_(
      raster_path = one_year_tif,
      output_dir = dirs$otsu,
      year = target_year,
      otsu_thresholds = sort(unique(c(otsu_threshold, reference_otsu_threshold))),
      min_otsu_threshold_value = min_otsu_threshold_value,
      use_original = FALSE,
      trim_percentiles = NULL,
      corine_raster_path = if (identical(otsu_mode, "corine")) corine_raster_path else NULL,
      peninsula_shapefile = if (identical(otsu_mode, "corine")) peninsula_shapefile else NULL,
      reclassify_corine = identical(otsu_mode, "corine"),
      reclass_matrix = if (identical(otsu_mode, "corine")) reclass_matrix else NULL,
      corine_classes = NULL,
      output_corine_raster_dir = if (identical(otsu_mode, "corine")) file.path(dirs$otsu, "output_corine") else NULL,
      reproject = TRUE,
      resolution = 90,
      python_exe = python_exe,
      gdal_polygonize_script = gdal_polygonize_script,
      gdalwarp_path = gdalwarp_path,
      output_format = "geojson",
      vectorize = FALSE,
      write_raster = TRUE,
      burnable_mask = terra::rast(burnable_mask_path)
    )
  } else {
    msg("STEP 1 - Reusing legacy OTSU rasters [%s] | candidate ge%d | reference ge%d", otsu_mode, otsu_threshold, reference_otsu_threshold)
  }
  stopifnot(file.exists(otsu_raster_path))
  stopifnot(file.exists(ref_raster_path))

  if (!isTRUE(reuse_ok) || !file.exists(patch_path)) {
    msg("STEP 2 - Polygonize legacy OTSU raster")
    patches_sf <- polygonize_Otsu(
      burn_raster = otsu_raster_path,
      python_exe = python_exe,
      gdal_polygonize_script = gdal_polygonize_script,
      tile = TRUE,
      n_rows = 2,
      n_cols = 3,
      tile_overlap = 1000,
      dissolve_tiles = FALSE,
      min_pixels = min_pixels,
      out_path = NULL,
      ogr2ogr_exe = ogr2ogr_exe
    )
    patches_sf <- sanitize_unb_legacy(patches_sf)
    remove_vector_sidecars_unb_legacy(patch_path)
    sf::st_write(patches_sf, patch_path, quiet = TRUE)
  } else {
    msg("STEP 2 - Reusing polygonized patches")
  }
  stopifnot(file.exists(patch_path))

  if (!isTRUE(reuse_ok) || !file.exists(coverage_path)) {
    msg("STEP 3 - Coverage by patch")
    coverage_by_patch_raster(
      patches = patch_path,
      ref_raster = list(terra::rast(ref_raster_path)),
      ref_names = c("ref_100"),
      out_path = coverage_path
    )
  } else {
    msg("STEP 3 - Reusing patch coverage")
  }
  stopifnot(file.exists(coverage_path))

  if (!isTRUE(reuse_ok) || !file.exists(decision_path)) {
    msg("STEP 4 - Legacy patch decisions")
    run_scenarios(
      patches_path = coverage_path,
      out_dir = dirs$decisions,
      cv_fields = "ref_100",
      buffers_m = buffers_m,
      core_thr = core_thr,
      alpha_boost = alpha_boost,
      min_base_boost = min_base_boost,
      dist_power = dist_power,
      keep_hi = keep_hi,
      drop_lo = drop_lo,
      target_epsg = target_epsg,
      dist_mode = dist_mode,
      near_mode = near_mode,
      write_audit = FALSE,
      write_all = TRUE,
      write_keep = TRUE,
      write_review = TRUE,
      keep_all_attrs = TRUE,
      keep_aux = TRUE,
      driver = "ESRI Shapefile"
    )
  } else {
    msg("STEP 4 - Reusing patch decisions")
  }
  stopifnot(file.exists(decision_path))

  # AS02 (0.5.0): persist the parameter fingerprint next to the cache. On the
  # next call, reuse is honoured only if this token still matches (see the
  # `reuse_ok` gate above). Written unconditionally so a previously absent /
  # mismatched fingerprint is refreshed after the stages were (re)computed.
  if (!isTRUE(cache_fingerprint_matches)) {
    tryCatch(
      writeLines(fingerprint_token, fingerprint_path),
      error = function(e)
        warning("Could not write legacy cache fingerprint (",
                conditionMessage(e), ")", call. = FALSE)
    )
  }

  msg("STEP 5 - Build unburned from legacy patch decisions")
  res_unb <- build_unburned_from_legacy_decisions(
    legacy_patches_path = decision_path,
    internal_decisions_path = internal_decisions_path,
    out_gpkg = if (isTRUE(write_unburned)) unburned_out_gpkg else NULL,
    use_drop = use_drop,
    use_review = use_review,
    use_keep = use_keep,
    drop_max_s_patch = drop_max_s_patch,
    review_max_s_patch = review_max_s_patch,
    keep_max_s_patch = keep_max_s_patch,
    exclude_buffer_m = exclude_buffer_m,
    min_area_ha = min_area_ha,
    sample_n = sample_n,
    sample_props = sample_props,
    random_seed = random_seed,
    allow_empty_otsu_pool = allow_empty_otsu_pool,
    verbose = verbose
  )

  invisible(list(
    target_year = target_year,
    scenario_name = scenario_name,
    out_root_dir = out_root_dir,
    otsu_raster_path = otsu_raster_path,
    ref_raster_path = ref_raster_path,
    helper_files = helper_files,
    patch_path = patch_path,
    coverage_path = coverage_path,
    decision_path = decision_path,
    unburned_out_gpkg = unburned_out_gpkg,
    unburned = res_unb
  ))
}
