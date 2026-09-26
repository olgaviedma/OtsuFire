# Content-aware cache keys for validate_fire_maps().
#
# The reference-side cache of validate_fire_maps() stores the masked, dissolved,
# observability-filtered reference polygons, the rasterized reference mask, and
# the observability diagnostics. Re-using that cache is only safe when the
# *content and methodological identity* of every input that shaped it is
# unchanged. Earlier keys folded only the observability raster BASENAME, the
# layer names, the DOY columns and a fingerprint of the reference vector. A
# changed burnable mask, study-area mask, resolution/CRS/grid, observability
# rule, or burnable threshold could silently serve a stale cache.
#
# These helpers compute a content-aware key. They are package-internal (not
# exported) and side-effect free so the cache contract can be unit tested
# without running the full (heavy) validator. `force_reprocess_ref = TRUE` in
# validate_fire_maps() still deletes the cache and forces a full rebuild
# regardless of the key.
#
# Identity proxy: for FILE inputs we hash basename + size + mtime + (for files
# at or below `max_md5_bytes`) the full md5 of the bytes; rasters additionally
# fold a header fingerprint (dims, resolution, extent, CRS, layer names) so a
# changed grid/CRS/band busts the key even when md5 is skipped for a very large
# file. For IN-MEMORY objects (SpatRaster / sf with no source on disk) we fold a
# structural fingerprint only; callers needing a guaranteed rebuild after an
# in-place mutation should pass force_reprocess_ref = TRUE.

# Bump this string whenever the temporal observability RULE itself changes
# (predicate for a valid observability pixel, the obs_doy_max / ref_doy
# comparison, or the observable_flag definition). It is folded into the
# observability fingerprint so a rule change invalidates every reference cache.
# obsrule-2 (0.11.1, 2026-08-19): the observability rule became selectable.
# The FAMILY of rules changed (a second whole-fire rule and a second required-
# DOY cascade were added, and the cached table gained v2 columns), so the
# family version is bumped once here; the per-run choice of mode / threshold /
# observability DOY column is folded separately in
# .vfm_observability_fingerprint() below.
.VFM_OBS_RULE_VERSION <- "obsrule-3:wholefire{legacy_any|fraction}+undetpolicy"

# Bump this whenever the SCHEMA of the cached reference artifact changes (e.g.
# new columns on the cached reference_observability table). It is folded into the
# reference cache key so a schema change rebuilds the cache instead of serving a
# structurally incompatible one. v2 adds the per-fire partial-observability
# audit columns (total_pixels / observable_pixels / non_observable_pixels /
# observable_fraction). v3 adds the whole-fire v2 columns
# (fire_required_doy / required_doy_source / observability_status /
# observability_mode / observability_min_fraction / observable_flag_legacy_any /
# data_coverage_fraction / n_temporally_observable /
# temporal_observable_fraction / no_data_fraction /
# pre_fire_or_too_early_fraction / post_fire_observable_fraction). v4 renames
# the observability state UNDETERMINED_END_DATE -> UNDETERMINED_OBS_DATE and
# adds `undetermined_cause`; bumping it guarantees that no table written with
# the 2.1.x vocabulary can ever be read back by this version, so there is no
# code path on which the two vocabularies could mix.
.VFM_REF_CACHE_SCHEMA <- "refschema-5"

.vfm_cache_safe_tag <- function(x) {
  x <- paste(x, collapse = "_")
  x <- gsub("[^A-Za-z0-9]+", "-", x)
  x <- gsub("(^-+|-+$)", "", x)
  if (!nzchar(x)) "na" else substr(x, 1L, 80L)
}

# Deterministic 8-hex polynomial hash (no external deps, stable across
# platforms/sessions). Identical algorithm to the historical in-function hash.
.vfm_short_cache_hash <- function(...) {
  x <- paste(unlist(list(...), use.names = FALSE), collapse = "|")
  ints <- utf8ToInt(enc2utf8(x))
  if (!length(ints)) return("00000000")
  mod <- 2147483647
  h <- 0
  for (ii in ints) {
    h <- (h * 131 + ii) %% mod
  }
  hx <- toupper(as.character(as.hexmode(h)))
  hx <- paste0(strrep("0", 8L), hx)
  substr(hx, nchar(hx) - 7L, nchar(hx))
}

# Content token for a single file on disk: basename + size + mtime + (for
# small-enough files) full md5 of the bytes. Returns a stable string.
.vfm_file_content_token <- function(path, max_md5_bytes = 64L * 1024L * 1024L) {
  fi <- file.info(path)
  size <- fi$size
  mt <- tryCatch(format(fi$mtime, "%Y%m%d%H%M%S", tz = "UTC"),
                 error = function(e) "namt")
  md5 <- if (isTRUE(is.finite(size)) && size <= max_md5_bytes) {
    tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  paste(
    basename(path),
    format(size, scientific = FALSE),
    mt,
    if (is.na(md5)) "nomd5" else md5,
    sep = "|"
  )
}

# Header fingerprint of a raster: dims, resolution, extent, CRS and layer names.
# Read lazily (no cell values loaded), so cheap even for huge rasters.
.vfm_raster_header <- function(r) {
  tryCatch(
    paste(
      "dim", paste(c(terra::nrow(r), terra::ncol(r), terra::nlyr(r)),
                   collapse = "x"),
      "res", paste(format(terra::res(r), trim = TRUE, scientific = FALSE),
                   collapse = "x"),
      "ext", paste(format(as.vector(terra::ext(r)), trim = TRUE,
                          scientific = FALSE), collapse = ","),
      "crs", terra::crs(r, describe = FALSE),
      "names", paste(names(r), collapse = "+"),
      sep = "|"
    ),
    error = function(e) "rheader-na"
  )
}

# Content/identity token for a raster input (path or SpatRaster).
.vfm_raster_token <- function(x) {
  if (is.null(x)) return("rnone")
  if (inherits(x, "SpatRaster")) {
    src <- tryCatch(terra::sources(x), error = function(e) character())
    src <- src[nzchar(src)]
    if (length(src) && file.exists(src[[1]])) {
      return(paste("rsrc", .vfm_file_content_token(src[[1]]),
                   .vfm_raster_header(x), sep = "|"))
    }
    return(paste("rmem", .vfm_raster_header(x), sep = "|"))
  }
  if (is.character(x) && length(x) == 1L && nzchar(x)) {
    np <- normalizePath(x, winslash = "/", mustWork = FALSE)
    if (file.exists(np)) {
      hdr <- tryCatch(.vfm_raster_header(terra::rast(np)),
                      error = function(e) "rheader-na")
      return(paste("rpath", .vfm_file_content_token(np), hdr, sep = "|"))
    }
    return(paste("rpath-missing", basename(np), sep = "|"))
  }
  paste("rclass", class(x)[1], sep = "|")
}

# Content/identity token for a vector input (path or sf). For shapefiles the
# .dbf (attributes), .prj (CRS) and .shx sidecars are folded in too.
.vfm_vector_token <- function(x) {
  if (is.null(x)) return("vnone")
  if (is.character(x) && length(x) == 1L && nzchar(x)) {
    np <- normalizePath(x, winslash = "/", mustWork = FALSE)
    parts <- c("vpath", basename(np))
    files <- np
    if (identical(tolower(tools::file_ext(np)), "shp")) {
      stem <- sub("\\.shp$", "", np, ignore.case = TRUE)
      files <- c(np, paste0(stem, c(".dbf", ".prj", ".shx")))
    }
    for (f in files) if (file.exists(f)) {
      parts <- c(parts, .vfm_file_content_token(f))
    }
    return(paste(parts, collapse = "|"))
  }
  if (inherits(x, "sf")) {
    bb <- tryCatch(sf::st_bbox(x), error = function(e) NULL)
    bb_str <- if (is.null(bb)) "nobbox" else
      paste(format(as.numeric(bb), trim = TRUE, nsmall = 0), collapse = ",")
    crs_str <- tryCatch(
      format(sf::st_crs(x)$wkt %||% sf::st_crs(x)$input %||% "nocrs"),
      error = function(e) "nocrs"
    )
    geom_str <- tryCatch(
      paste(as.character(sf::st_geometry_type(x, by_geometry = FALSE)),
            collapse = "+"),
      error = function(e) "nogeom"
    )
    return(paste("vsf", format(nrow(x), scientific = FALSE),
                 bb_str, crs_str, geom_str, sep = "|"))
  }
  paste("vclass", class(x)[1], sep = "|")
}

# Public-style "in-XXXX" fingerprint of a vector input (used for the prediction
# cache key, preserving the historical tag shape).
.vfm_vector_fingerprint <- function(x) {
  if (is.null(x)) return("in-none")
  paste0("in-", .vfm_short_cache_hash(.vfm_cache_safe_tag(.vfm_vector_token(x))))
}

# "obs-XXXX": observability raster content + layer + DOY columns + rule version
# + the per-run rule configuration (mode, threshold, observability DOY column).
#
# `observability_min_fraction` is folded ONLY under the fraction mode: under
# "legacy_any" the threshold is inert, and folding it there would invalidate
# caches for a parameter that cannot change the result. `ref_obs_doy_col` is
# folded whenever it is set, because it changes the required DOY cascade.
.vfm_observability_fingerprint <- function(observability_raster,
                                           ref_end_doy_col,
                                           ref_start_doy_col,
                                           observability_mode = "legacy_any",
                                           observability_min_fraction = 0.75,
                                           ref_obs_doy_col = NULL,
                                           ref_obs_doy_fallback = "none",
                                           observability_undetermined_policy = "evaluate") {
  if (is.null(observability_raster)) return("obs-none")
  layer_id <- tryCatch({
    if (inherits(observability_raster, "SpatRaster")) {
      paste(names(observability_raster), collapse = "+")
    } else if (is.character(observability_raster) &&
               length(observability_raster) == 1L) {
      paste(names(terra::rast(observability_raster)), collapse = "+")
    } else {
      class(observability_raster)[1]
    }
  }, error = function(e) "layer-na")
  frac_tag <- if (identical(observability_mode, "legacy_any")) {
    "minfrac-inert"
  } else {
    paste0("minfrac", format(observability_min_fraction, scientific = FALSE))
  }
  paste0("obs-", .vfm_short_cache_hash(
    .vfm_raster_token(observability_raster),
    .vfm_cache_safe_tag(layer_id),
    .vfm_cache_safe_tag(ref_end_doy_col),
    .vfm_cache_safe_tag(ref_start_doy_col),
    paste0("mode-", .vfm_cache_safe_tag(observability_mode)),
    frac_tag,
    paste0("obsdoycol-", if (is.null(ref_obs_doy_col)) "none" else
      .vfm_cache_safe_tag(ref_obs_doy_col)),
    # The authoritative-vs-fallback policy changes which fires end up
    # UNDETERMINED, so it must bust the reference cache. Inert when no
    # observability column is supplied.
    paste0("obsdoyfb-", if (is.null(ref_obs_doy_col)) "inert" else
      .vfm_cache_safe_tag(ref_obs_doy_fallback)),
    # The undetermined policy decides which fires stay in the evaluation
    # domain, so it reshapes the cached reference artifact. Inert under legacy.
    paste0("undet-", if (identical(observability_mode, "legacy_any")) "inert" else
      .vfm_cache_safe_tag(observability_undetermined_policy)),
    .VFM_OBS_RULE_VERSION
  ))
}

# "dom-XXXX": the burnable domain identity (burnable raster grid/CRS/content +
# study-area mask + burnable thresholding settings). Shared by the reference
# and prediction caches because both are rasterized onto / masked by this
# domain.
.vfm_domain_fingerprint <- function(burnable_raster,
                                    mask_shapefile,
                                    binary_burnable,
                                    burnable_classes,
                                    burnable_threshold) {
  paste0("dom-", .vfm_short_cache_hash(
    .vfm_raster_token(burnable_raster),
    .vfm_vector_token(mask_shapefile),
    paste0("bin", as.integer(isTRUE(binary_burnable))),
    paste0("bcls", if (is.null(burnable_classes)) "none" else
      paste(burnable_classes, collapse = "-")),
    paste0("bthr", format(burnable_threshold, scientific = FALSE))
  ))
}

# Full reference cache tag: obs identity + domain identity + reference vector
# content + reference-only options (min area, dissolve field). `buffer` is NOT
# folded in: it is applied downstream to the loaded reference at detection time
# and never shapes the cached artifact, so it must not bust the reference cache.
.vfm_reference_cache_key <- function(obs_cache_tag,
                                     dom_cache_tag,
                                     ref_shapefile,
                                     min_area_reference_ha,
                                     dissolve_ref_by) {
  opts <- paste0(
    "minha", if (is.null(min_area_reference_ha)) "none" else
      format(min_area_reference_ha, scientific = FALSE),
    "|dslv", if (is.null(dissolve_ref_by)) "none" else
      paste(dissolve_ref_by, collapse = "-"),
    "|", .VFM_REF_CACHE_SCHEMA
  )
  paste0(
    obs_cache_tag, "_", dom_cache_tag, "_ref-",
    .vfm_short_cache_hash(.vfm_vector_token(ref_shapefile), opts)
  )
}
