# ==============================================================================
# make_blockcv_folds()
# Create spatial blocks (grid) + assign k-folds by block (repeated) for OOF / block-CV
#
# Mini “for dummies”
# - fold: one partition (1..k) used as VALIDATION; the rest is TRAINING
# - CV: repeat train/validation rotating the folds (k rounds)
# - OOF: a row's prediction made ONLY while that row sat in validation (unseen during training)
# - block-CV: folds assigned by spatial BLOCKS to avoid spatial leakage
# ==============================================================================

make_block_folds <- function(
    # --- Input: provide either an sf object OR a gpkg path+layer ---
  train_labelled_sf     = NULL,   # sf (polygons)
  train_labelled        = NULL,   # path to .gpkg
  train_labelled_layer  = NULL,   # layer name inside gpkg
  
  # --- Columns / labels ---
  class_col   = "class",
  pos_lab     = "burned",
  neg_lab     = "unburned",
  fire_id_col = "fire_uid",  # used if split_unit="fire"
  
  # --- What unit should be kept together in CV? ---
  # "fire"   = each fire_uid belongs to exactly one fold (recommended for transferability)
  # "polygon"= each polygon belongs to one fold (OK if you truly have per-polygon independent units)
  split_unit = c("fire", "polygon"),
  
  # --- Auto-tuning candidates (tried in this order: stricter -> looser) ---
  block_sizes_m = c(5000, 3000, 2000),
  k_candidates  = c(5, 4, 3),
  
  # --- Acceptance thresholds (per fold, repeat 1) ---
  # For split_unit="fire": burned units = burned fires
  # For split_unit="polygon": burned units = burned polygons
  min_burned_units_per_fold = 3,
  min_pos_blocks_per_fold   = 10,
  
  # --- CV repeats ---
  n_repeats = 2,
  seed_base = 42,
  balance_unburned_by = c("blocks", "polys"),
  
  # --- CRS handling (block size is meters) ---
  lonlat_action = c("warn", "stop", "transform"),
  crs_projected = 3035,  # used only if lonlat_action="transform"
  
  # --- Outputs ---
  out_dir       = NULL,
  target_year   = NULL,
  out_prefix    = NULL,
  write_blocks_gpkg         = TRUE,
  write_folds_csv           = TRUE,
  write_train_with_folds_gpkg = TRUE, # writes polygons + fold_rep* to GPKG
  folds_csv_id_col          = NULL,    # default: fire_id_col if split_unit="fire", else "poly_id"
  
  # --- Verbosity ---
  verbose = TRUE
) {
  split_unit <- match.arg(split_unit)
  lonlat_action <- match.arg(lonlat_action)
  balance_unburned_by <- match.arg(balance_unburned_by)
  
  # Block 7: library() calls removed; deps in DESCRIPTION Imports.
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---------------------------------------------------------------------------
  # 1) Read input
  # ---------------------------------------------------------------------------
  if (is.null(train_labelled_sf)) {
    if (is.null(train_labelled) || is.null(train_labelled_layer)) {
      stop("Provide either `train_labelled_sf` OR (`train_labelled` + `train_labelled_layer`).")
    }
    train_labelled_sf <- sf::st_read(train_labelled, layer = train_labelled_layer, quiet = TRUE)
  }
  if (!inherits(train_labelled_sf, "sf")) stop("train_labelled_sf must be an sf object.")
  
  # Basic label checks
  if (!(class_col %in% names(train_labelled_sf))) stop(sprintf("Missing column `%s`.", class_col))
  labs <- unique(train_labelled_sf[[class_col]])
  if (!all(c(pos_lab, neg_lab) %in% labs)) {
    stop(sprintf("`%s` must include labels '%s' and '%s'. Found: %s",
                 class_col, pos_lab, neg_lab, paste(labs, collapse = ", ")))
  }
  
  # Ensure polygon unique id exists (needed if split_unit="polygon" or for safe joins)
  if (!("poly_id" %in% names(train_labelled_sf))) {
    train_labelled_sf$poly_id <- seq_len(nrow(train_labelled_sf))
  }
  
  # If split by fire: need fire_id_col
  if (split_unit == "fire" && !(fire_id_col %in% names(train_labelled_sf))) {
    stop(sprintf("split_unit='fire' requires fire_id_col='%s' present in data.", fire_id_col))
  }
  
  # Decide what ID column to write to the folds CSV
  if (is.null(folds_csv_id_col)) {
    folds_csv_id_col <- if (split_unit == "fire") fire_id_col else "poly_id"
  }
  if (!(folds_csv_id_col %in% names(train_labelled_sf))) {
    # For fire-level, folds_csv_id_col is in polygons table, so ok.
    # If user passed something else, stop.
    stop(sprintf("folds_csv_id_col='%s' not found in train_labelled_sf.", folds_csv_id_col))
  }
  
  # ---------------------------------------------------------------------------
  # 2) CRS sanity (block size is meters)
  # ---------------------------------------------------------------------------
  if (sf::st_is_longlat(train_labelled_sf)) {
    if (lonlat_action == "stop") {
      stop("CRS is lon/lat (degrees). Transform to projected CRS in meters (e.g. EPSG:3035).")
    } else if (lonlat_action == "transform") {
      .msg("CRS is lon/lat -> transforming to EPSG:%s for block-CV", crs_projected)
      train_labelled_sf <- sf::st_transform(train_labelled_sf, crs_projected)
    } else {
      warning("CRS is lon/lat (degrees). block_size_m is meters -> results may be wrong. Consider lonlat_action='transform'.")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3) Helpers: blocks + assignment + fold allocation
  # ---------------------------------------------------------------------------
  
  # --- build block grid on bbox (padded) ---
  .make_blocks <- function(sfobj, block_size_m) {
    bb <- sf::st_bbox(sfobj)
    bb["xmin"] <- bb["xmin"] - block_size_m
    bb["ymin"] <- bb["ymin"] - block_size_m
    bb["xmax"] <- bb["xmax"] + block_size_m
    bb["ymax"] <- bb["ymax"] + block_size_m
    
    bb_sfc <- sf::st_as_sfc(bb)
    sf::st_crs(bb_sfc) <- sf::st_crs(sfobj)
    
    grid <- sf::st_make_grid(bb_sfc, cellsize = block_size_m, what = "polygons", square = TRUE)
    blocks <- sf::st_sf(block_id = seq_along(grid), geometry = grid)
    sf::st_crs(blocks) <- sf::st_crs(sfobj)
    blocks
  }
  
  # --- assign block_id to sf rows using point_on_surface on geometry only (no warning) ---
  .assign_block_id <- function(sfobj, blocks) {
    geom_valid <- sf::st_make_valid(sf::st_geometry(sfobj))
    pts <- sf::st_point_on_surface(geom_valid)                 # sfc
    pts <- sf::st_as_sf(pts)                                   # sf with only geometry
    sf::st_crs(pts) <- sf::st_crs(sfobj)
    
    idx <- sf::st_intersects(pts, blocks)
    block_id <- vapply(idx, function(x) if (length(x)) blocks$block_id[x[1]] else NA_integer_, integer(1))
    
    na <- which(is.na(block_id))
    if (length(na) > 0) {
      nearest <- sf::st_nearest_feature(pts[na, , drop = FALSE], blocks)
      block_id[na] <- blocks$block_id[nearest]
    }
    block_id
  }
  
  # --- two-stage fold assignment by block_id ---
  .assign_folds_by_block_two_stage <- function(df_units, k_folds, n_repeats, seed, balance_unburned_by) {
    stopifnot(all(c("block_id", "y") %in% names(df_units)))  # y: 1=pos, 0=neg
    
    blk <- df_units %>%
      group_by(block_id) %>%
      summarise(
        n_total = n(),
        n_pos   = sum(y == 1),
        n_neg   = sum(y == 0),
        .groups = "drop"
      )
    
    blk_pos <- blk %>% filter(n_pos > 0)
    blk_neg <- blk %>% filter(n_pos == 0)
    
    one_repeat <- function(rseed) {
      set.seed(rseed)
      
      fold_pos_blocks <- integer(k_folds)
      fold_blocks_tot <- integer(k_folds)
      fold_units_tot  <- integer(k_folds)
      
      pick_fold <- function(cand) cand[sample.int(length(cand), 1)]
      
      # Stage 1: distribute POS blocks
      blk_pos2 <- blk_pos %>%
        mutate(j = runif(n())) %>%
        arrange(desc(n_pos + 1e-6*j), desc(n_total), j)
      
      fold_for_pos <- integer(nrow(blk_pos2))
      
      for (i in seq_len(nrow(blk_pos2))) {
        cand <- which(fold_pos_blocks == min(fold_pos_blocks))
        if (length(cand) > 1) cand <- cand[fold_units_tot[cand] == min(fold_units_tot[cand])]
        f <- pick_fold(cand)
        
        fold_for_pos[i] <- f
        fold_pos_blocks[f] <- fold_pos_blocks[f] + 1L
        fold_blocks_tot[f] <- fold_blocks_tot[f] + 1L
        fold_units_tot[f]  <- fold_units_tot[f]  + blk_pos2$n_total[i]
      }
      
      map_pos <- blk_pos2 %>% transmute(block_id, fold = fold_for_pos)
      
      # Stage 2: distribute NEG-only blocks
      blk_neg2 <- blk_neg %>%
        mutate(j = runif(n())) %>%
        arrange(desc(n_total + 1e-6*j), j)
      
      fold_for_neg <- integer(nrow(blk_neg2))
      
      for (i in seq_len(nrow(blk_neg2))) {
        if (balance_unburned_by == "blocks") {
          cand <- which(fold_blocks_tot == min(fold_blocks_tot))
          if (length(cand) > 1) cand <- cand[fold_units_tot[cand] == min(fold_units_tot[cand])]
        } else {
          cand <- which(fold_units_tot == min(fold_units_tot))
          if (length(cand) > 1) cand <- cand[fold_blocks_tot[cand] == min(fold_blocks_tot[cand])]
        }
        
        f <- pick_fold(cand)
        
        fold_for_neg[i] <- f
        fold_blocks_tot[f] <- fold_blocks_tot[f] + 1L
        fold_units_tot[f]  <- fold_units_tot[f]  + blk_neg2$n_total[i]
      }
      
      map_neg <- blk_neg2 %>% transmute(block_id, fold = fold_for_neg)
      
      bind_rows(map_pos, map_neg)
    }
    
    out <- df_units
    for (r in seq_len(n_repeats)) {
      map_r <- one_repeat(seed + 10000*r)
      col   <- paste0("fold_rep", r)
      out <- out %>% left_join(rename(map_r, !!col := fold), by = "block_id")
      stopifnot(!anyNA(out[[col]]))
    }
    
    out
  }
  
  # ---------------------------------------------------------------------------
  # 4) Build “units” table: either fires (recommended) or polygons
  # ---------------------------------------------------------------------------
  .build_units <- function(train_sf, split_unit, fire_id_col, class_col, pos_lab, neg_lab) {
    if (split_unit == "polygon") {
      units_sf <- train_sf %>%
        select(poly_id, !!class_col) %>%
        mutate(y = ifelse(.data[[class_col]] == pos_lab, 1L, 0L))
      list(units_sf = units_sf, unit_id_col = "poly_id")
    } else {
      # Fire-level units: one geometry per fire_uid (FAST: st_combine, not union)
      # Bug 6 (0.3.0): resolve the active sf geometry column name (which
      # may be `geom` after a GPKG round-trip) and rename it to
      # `geometry` before grouping, instead of hard-coding `geometry`
      # in the summarise call. Downstream code keeps the stable
      # `geometry` contract.
      .geom_col <- attr(train_sf, "sf_column") %||% "geometry"
      if (!identical(.geom_col, "geometry") &&
          .geom_col %in% names(train_sf)) {
        names(train_sf)[names(train_sf) == .geom_col] <- "geometry"
        attr(train_sf, "sf_column") <- "geometry"
      }
      units_sf <- train_sf %>%
        group_by(.data[[fire_id_col]]) %>%
        summarise(
          # label a fire as positive if it has ANY burned polygon
          y = as.integer(any(.data[[class_col]] == pos_lab)),
          .groups = "drop",
          geometry = sf::st_combine(geometry)
        ) %>%
        sf::st_as_sf()

      names(units_sf)[names(units_sf) == fire_id_col] <- "fire_uid_tmp"
      list(units_sf = units_sf, unit_id_col = "fire_uid_tmp")
    }
  }
  
  units_pack <- .build_units(train_labelled_sf, split_unit, fire_id_col, class_col, pos_lab, neg_lab)
  units_sf   <- units_pack$units_sf
  unit_id_col <- units_pack$unit_id_col
  
  n_pos_units <- sum(units_sf$y == 1)
  if (n_pos_units == 0) stop("No positive units found (no burned). Block-CV cannot be created.")
  
  .msg("Split unit: %s | positive units=%d | total units=%d", split_unit, n_pos_units, nrow(units_sf))
  
  # ---------------------------------------------------------------------------
  # 5) AUTO-TUNE: try (k, block_size) combos and pick the best transferable config
  # ---------------------------------------------------------------------------
  
  .evaluate_config <- function(k, bs) {
    blocks <- .make_blocks(units_sf, bs)
    units_sf2 <- units_sf
    units_sf2$block_id <- .assign_block_id(units_sf2, blocks)
    
    df_units <- sf::st_drop_geometry(units_sf2) %>%
      select(all_of(unit_id_col), y, block_id)
    
    df_folds <- .assign_folds_by_block_two_stage(
      df_units = df_units,
      k_folds = k,
      n_repeats = n_repeats,
      seed = seed_base,
      balance_unburned_by = balance_unburned_by
    )
    
    # Metrics on repeat 1
    burned_units_by_fold <- df_folds %>%
      filter(y == 1) %>%
      distinct(.data[[unit_id_col]], fold_rep1) %>%
      count(fold_rep1, name = "n_burned_units")
    
    # Ensure all folds appear (if some fold has 0 positives, it will be missing)
    all_folds <- tibble(fold_rep1 = 1:k)
    burned_units_by_fold <- all_folds %>%
      left_join(burned_units_by_fold, by = "fold_rep1") %>%
      mutate(n_burned_units = ifelse(is.na(n_burned_units), 0L, n_burned_units))
    
    pos_blocks_by_fold <- df_folds %>%
      filter(y == 1) %>%
      distinct(block_id, fold_rep1) %>%
      count(fold_rep1, name = "n_pos_blocks")
    
    pos_blocks_by_fold <- all_folds %>%
      left_join(pos_blocks_by_fold, by = "fold_rep1") %>%
      mutate(n_pos_blocks = ifelse(is.na(n_pos_blocks), 0L, n_pos_blocks))
    
    min_burned_units <- min(burned_units_by_fold$n_burned_units)
    min_pos_blocks   <- min(pos_blocks_by_fold$n_pos_blocks)
    
    ok <- (min_burned_units >= min_burned_units_per_fold) &&
      (min_pos_blocks   >= min_pos_blocks_per_fold)
    
    # Score: prefer stricter configs IF ok; otherwise maximize minima
    score <- if (ok) {
      # prefer larger blocks + larger k
      (bs/1000) * 10 + k
    } else {
      # still keep the “best available” if none pass
      min_burned_units * 1000 + min_pos_blocks
    }
    
    list(
      ok = ok,
      score = score,
      k = k,
      block_size_m = bs,
      min_burned_units = min_burned_units,
      min_pos_blocks   = min_pos_blocks,
      burned_units_by_fold = burned_units_by_fold,
      pos_blocks_by_fold   = pos_blocks_by_fold,
      blocks = blocks,
      df_folds_units = df_folds  # folds at UNIT level
    )
  }
  
  # Bug 4 (0.3.0): track best-OK and best-non-OK configurations
  # separately, and only fall back to a non-OK candidate when no OK
  # candidate exists. On fallback emit a loud warning() and write a
  # `_FOLD_FALLBACK.txt` audit file in the run directory.
  best_ok <- NULL
  best_nonok <- NULL
  for (k in k_candidates) {
    for (bs in block_sizes_m) {
      .msg("Trying config: block_size_m=%d | k_folds=%d ...", bs, k)
      ev <- .evaluate_config(k, bs)
      .msg("  -> min burned units/fold=%d | min pos blocks/fold=%d | ok=%s",
           ev$min_burned_units, ev$min_pos_blocks, ev$ok)

      if (isTRUE(ev$ok)) {
        if (is.null(best_ok) || ev$score > best_ok$score) best_ok <- ev
      } else {
        if (is.null(best_nonok) || ev$score > best_nonok$score) best_nonok <- ev
      }
    }
  }

  if (!is.null(best_ok)) {
    best <- best_ok
    .msg("Selected (OK): block_size_m=%d | k_folds=%d | min burned units/fold=%d | min pos blocks/fold=%d",
         best$block_size_m, best$k, best$min_burned_units, best$min_pos_blocks)
  } else {
    if (is.null(best_nonok)) {
      stop("make_block_folds(): no fold-tuning candidate could be evaluated.",
           call. = FALSE)
    }
    best <- best_nonok
    fallback_msg <- sprintf(
      paste0("make_block_folds(): no candidate (k, block_size_m) passed the ",
             "acceptance gate (min_burned_units_per_fold=%d, ",
             "min_pos_blocks_per_fold=%d). Falling back to best non-OK ",
             "candidate: block_size_m=%d | k_folds=%d (min burned units/fold=%d; ",
             "min pos blocks/fold=%d)."),
      min_burned_units_per_fold, min_pos_blocks_per_fold,
      best$block_size_m, best$k, best$min_burned_units, best$min_pos_blocks
    )
    warning(fallback_msg, call. = FALSE)
    if (!is.null(out_dir) && nzchar(out_dir)) {
      tryCatch({
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines(
          c(
            "OtsuFire 0.3.0 fold-tuning soft fallback (Bug 4)",
            sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            sprintf("min_burned_units_per_fold = %d", min_burned_units_per_fold),
            sprintf("min_pos_blocks_per_fold   = %d", min_pos_blocks_per_fold),
            "",
            sprintf("Selected configuration (BEST NON-OK):"),
            sprintf("  block_size_m         = %d", best$block_size_m),
            sprintf("  k_folds              = %d", best$k),
            sprintf("  min burned units/fold = %d", best$min_burned_units),
            sprintf("  min pos blocks/fold   = %d", best$min_pos_blocks),
            "",
            "All candidates were rejected by the acceptance gate.",
            "Re-run with stricter input data or relaxed thresholds."
          ),
          file.path(out_dir, "_FOLD_FALLBACK.txt")
        )
      }, error = function(e) {
        warning("Could not write _FOLD_FALLBACK.txt: ",
                conditionMessage(e), call. = FALSE)
      })
    }
  }
  
  # ---------------------------------------------------------------------------
  # 6) Propagate UNIT-level folds back to POLYGON table
  # ---------------------------------------------------------------------------
  units_folds <- best$df_folds_units
  
  if (split_unit == "polygon") {
    # Join by poly_id
    fold_cols <- paste0("fold_rep", seq_len(n_repeats))
    train_with_folds <- train_labelled_sf %>%
      left_join(units_folds %>% select(poly_id, all_of(fold_cols), block_id),
                by = "poly_id")
  } else {
    # Join by fire_uid_tmp -> back to original fire_id_col
    fold_cols <- paste0("fold_rep", seq_len(n_repeats))
    map_fire <- units_folds %>%
      rename(!!fire_id_col := all_of(unit_id_col)) %>%
      select(all_of(fire_id_col), block_id, all_of(fold_cols))
    
    train_with_folds <- train_labelled_sf %>%
      left_join(map_fire, by = fire_id_col)
  }
  
  # sanity
  stopifnot(!anyNA(train_with_folds[[paste0("fold_rep", 1)]]))
  
  # ---------------------------------------------------------------------------
  # 7) Write outputs
  # ---------------------------------------------------------------------------
  saved <- list(blocks_gpkg = NULL, folds_csv = NULL, train_gpkg = NULL)
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (is.null(out_prefix)) {
      out_prefix <- if (!is.null(target_year)) sprintf("%d", target_year) else "blockcv"
    }
    
    # 7.1 blocks gpkg
    if (isTRUE(write_blocks_gpkg)) {
      saved$blocks_gpkg <- file.path(out_dir, sprintf("%s_blocks_%dm.gpkg", out_prefix, best$block_size_m))
      sf::st_write(best$blocks, saved$blocks_gpkg, layer = "blocks", delete_layer = TRUE, quiet = TRUE)
      .msg("Saved blocks GPKG: %s", saved$blocks_gpkg)
    }
    
    # 7.2 folds csv (POLYGON level, because that’s what you’ll use downstream)
    if (isTRUE(write_folds_csv)) {
      fold_cols <- paste0("fold_rep", seq_len(n_repeats))
      
      folds_df <- sf::st_drop_geometry(train_with_folds) %>%
        select(all_of(folds_csv_id_col), block_id, all_of(class_col), all_of(fold_cols))
      
      saved$folds_csv <- file.path(
        out_dir,
        sprintf("%s_folds_%dm_k%d_r%d_unit-%s.csv",
                out_prefix, best$block_size_m, best$k, n_repeats, split_unit)
      )
      write.csv(folds_df, saved$folds_csv, row.names = FALSE)
      .msg("Saved folds CSV: %s", saved$folds_csv)
    }
    
    # 7.3 optional gpkg with polygons + folds
    if (isTRUE(write_train_with_folds_gpkg)) {
      saved$train_gpkg <- file.path(out_dir, sprintf("%s_train_with_folds_%dm.gpkg", out_prefix, best$block_size_m))
      sf::st_write(train_with_folds, saved$train_gpkg, layer = "train_with_folds", delete_layer = TRUE, quiet = TRUE)
      .msg("Saved train_with_folds GPKG: %s", saved$train_gpkg)
    }
  }
  
  # ---------------------------------------------------------------------------
  # 8) Return
  # ---------------------------------------------------------------------------
  # Item 10 (Gate 1B, 2026-06-07): REPRODUCIBLE fold fingerprint.
  # `make_fold_fingerprint()` builds a deterministic identity of the fold
  # configuration + the resulting fold assignments. It DELIBERATELY EXCLUDES
  # any wall-clock time, so two equivalent fold runs (same inputs, same seed,
  # same config) produce IDENTICAL fingerprints. The wall-clock timestamp is
  # kept as a SEPARATE `created_at` metadata field that is NOT hashed into the
  # fingerprint (it documents WHEN the folds were built, never WHICH folds).
  fold_params <- list(
    split_unit = split_unit,
    fire_id_col = fire_id_col,
    class_col = class_col,
    pos_lab = pos_lab,
    neg_lab = neg_lab,
    block_sizes_m = block_sizes_m,
    k_candidates = k_candidates,
    n_repeats = n_repeats,
    seed_base = seed_base,
    balance_unburned_by = balance_unburned_by,
    lonlat_action = lonlat_action,
    crs_projected = crs_projected
  )
  fold_fingerprint <- make_fold_fingerprint(
    params           = fold_params,
    selected         = list(
      block_size_m = best$block_size_m,
      k_folds      = best$k,
      ok           = best$ok
    ),
    train_with_folds = train_with_folds,
    fold_cols        = paste0("fold_rep", seq_len(n_repeats)),
    id_col           = folds_csv_id_col
  )

  list(
    train_with_folds = train_with_folds,  # polygons + fold_rep*
    blocks           = best$blocks,       # chosen grid
    selected         = list(
      block_size_m = best$block_size_m,
      k_folds      = best$k,
      ok           = best$ok,
      min_burned_units_per_fold = best$min_burned_units,
      min_pos_blocks_per_fold   = best$min_pos_blocks
    ),
    diagnostics       = list(
      burned_units_by_fold = best$burned_units_by_fold,
      pos_blocks_by_fold   = best$pos_blocks_by_fold
    ),
    saved = saved,
    params = fold_params,
    # Item 10: the REPRODUCIBLE identity of this fold set (no wall-clock).
    fingerprint = fold_fingerprint,
    # Item 10: SEPARATE metadata. `created_at` is wall-clock provenance only;
    # it is NOT part of `fingerprint` and never participates in any
    # reproducibility / equivalence comparison.
    created_at = Sys.time()
  )
}

# Item 10 (Gate 1B, 2026-06-07): deterministic, wall-clock-FREE fingerprint of
# a fold set. Two equivalent runs (same inputs, same seed, same config) yield an
# IDENTICAL fingerprint; the wall-clock `created_at` is recorded separately by
# the caller and is intentionally NOT fed in here. Base R only (no `digest` in
# Imports): sorted key=value text plus a small order-stable rolling checksum,
# mirroring `otsu_negative_param_fingerprint()`.
make_fold_fingerprint <- function(params, selected, train_with_folds,
                                  fold_cols, id_col = NULL) {
  flat1 <- function(v) {
    if (is.null(v)) return("NULL")
    v <- unlist(v, use.names = TRUE)
    if (!is.null(names(v)) && any(nzchar(names(v)))) {
      v <- v[order(names(v))]
      paste(sprintf("%s=%s", names(v),
                    format(v, trim = TRUE, scientific = FALSE)),
            collapse = ",")
    } else {
      paste(format(v, trim = TRUE, scientific = FALSE), collapse = ",")
    }
  }

  # Config + selected layer.
  cfg_kv <- vapply(params, flat1, character(1))
  cfg_kv <- sprintf("%s=%s", names(params), cfg_kv)
  cfg_kv <- cfg_kv[order(names(params))]
  sel_kv <- vapply(selected, flat1, character(1))
  sel_kv <- sprintf("selected.%s=%s", names(selected), sel_kv)
  sel_kv <- sel_kv[order(names(selected))]

  # Fold-assignment layer: the actual per-unit fold labels, ordered by a stable
  # id (when available) so row order cannot perturb the fingerprint. This is
  # what makes the fingerprint identify the ACTUAL folds, not merely the config.
  assign_kv <- character(0)
  df <- tryCatch(sf::st_drop_geometry(train_with_folds),
                 error = function(e) as.data.frame(train_with_folds))
  fold_cols <- intersect(fold_cols, names(df))
  if (length(fold_cols) > 0L) {
    ord <- if (!is.null(id_col) && id_col %in% names(df)) {
      order(as.character(df[[id_col]]))
    } else {
      seq_len(nrow(df))
    }
    for (fc in sort(fold_cols)) {
      vals <- as.integer(df[[fc]])[ord]
      assign_kv <- c(assign_kv,
                     sprintf("assign.%s=%s", fc,
                             paste(vals, collapse = ",")))
    }
  }

  body <- paste(c(cfg_kv, sel_kv, assign_kv), collapse = "\n")
  bytes <- as.numeric(charToRaw(enc2utf8(body)))
  chk <- 0
  for (b in bytes) chk <- (chk * 31 + b) %% 1000000007
  list(text = body, checksum = sprintf("%09d", as.integer(chk)))
}
