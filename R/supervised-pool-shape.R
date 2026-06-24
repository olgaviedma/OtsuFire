# =============================================================================
# supervised-pool-shape.R  (Phase 2B)
#
# Champion pool-build helpers (packaged from build_champion_pool.R):
#   * add_pool_shape_features()       -> the 6 shape/size columns + random
#                                        background DECORRELATION.
#   * subsample_random_background()   -> SUBR05 = random reduced to
#                                        random_to_burned_ratio x n_positives.
# Neither trains a model nor rebuilds the geometric pool; they only post-process
# an already-assembled training pool.
# =============================================================================

# The 6 champion shape/size columns (read by training, never recomputed there).
.of_champion_shape_cols <- function() {
  c("area_ha", "n_pix", "log_area", "perim_m", "compactness", "elongation")
}

#' Add the 6 champion shape/size features (+ decorrelate the random background)
#'
#' Mirrors \code{build_champion_pool.R}: \code{area_ha} from geometry,
#' \code{n_pix} from the existing column, and \code{log_area}/\code{perim_m}/
#' \code{compactness}/\code{elongation} from \code{.of_shape_features()}. With
#' \code{decorrelate_random = TRUE} each \code{random_burnable_background} row's
#' FULL shape tuple is replaced by a bootstrap copy of a randomly drawn
#' burned-class row's tuple (so synthetic squares carry no geometric giveaway,
#' while keep/otsu/hard-neg keep their real shape). Geometry is preserved.
#'
#' @param pool An \code{sf} training pool with \code{id_col}, \code{n_pix},
#'   \code{source} and \code{class} columns.
#' @param id_col Join id (default \code{"fire_uid"}).
#' @param decorrelate_random Logical (default \code{TRUE}).
#' @param seed Integer seed for the decorrelation bootstrap (default \code{1985L}).
#' @param random_source Source tag of the random background bucket.
#' @param positive_class Class tag of positives (default \code{"burned"}).
#' @return \code{pool} with the 6 shape columns added (geometry preserved).
#' @export
add_pool_shape_features <- function(pool, id_col = "fire_uid",
                                    decorrelate_random = TRUE, seed = 1985L,
                                    random_source = "random_burnable_background",
                                    positive_class = "burned") {
  if (!inherits(pool, "sf")) stop("add_pool_shape_features(): 'pool' must be an sf with geometry.", call. = FALSE)
  SHAPE  <- .of_champion_shape_cols()
  COMPUT <- c("log_area", "perim_m", "compactness", "elongation")
  if (is.na(sf::st_crs(pool)) || is.na(sf::st_crs(pool)$epsg) || sf::st_crs(pool)$epsg != 3035) {
    pool <- sf::st_transform(pool, 3035)
  }
  geom <- sf::st_geometry(pool)
  d <- sf::st_drop_geometry(pool)                       # plain data.frame from here
  if (!id_col %in% names(d)) stop("add_pool_shape_features(): pool lacks id column '", id_col, "'.", call. = FALSE)
  if (!"n_pix" %in% names(d)) stop("add_pool_shape_features(): pool lacks 'n_pix'.", call. = FALSE)
  d$area_ha <- as.numeric(sf::st_area(pool)) / 1e4
  d$n_pix   <- suppressWarnings(as.numeric(d$n_pix))

  shp <- .of_shape_features(pool, id_col = id_col)       # 4 geometry-derived cols
  idx <- match(d[[id_col]], shp[[id_col]])
  for (col in COMPUT) d[[col]] <- as.numeric(shp[[col]][idx])

  if (isTRUE(decorrelate_random)) {
    is_rand <- as.character(d$source) == random_source
    is_burn <- as.character(d$class)  == positive_class
    nb <- sum(is_burn); nr <- sum(is_rand)
    if (nb == 0L) stop("add_pool_shape_features(): no burned rows to sample shape from.", call. = FALSE)
    if (nr > 0L) {
      set.seed(seed)
      donor <- sample(which(is_burn), size = nr, replace = TRUE)
      don <- d[donor, SHAPE, drop = FALSE]               # snapshot donor tuples first
      for (col in SHAPE) d[[col]][is_rand] <- as.numeric(don[[col]])
    }
  }
  for (col in SHAPE) d[[col]] <- as.numeric(d[[col]])     # enforce plain numeric
  sf::st_sf(d, geometry = geom)
}

#' Subsample the random background bucket (SUBR05)
#'
#' SUBR05 = the random background is reduced to
#' \code{ceiling(random_to_burned_ratio * n_positives)} rows. Verified against
#' the champion pool on disk: random 51693 -> 714 = 0.5 x 1428 positives (NOT
#' 0.5 of the random bucket). Only \code{random_source} rows are touched;
#' positives, Otsu residual and artifact_hard are left intact. Reproducible
#' under \code{seed}.
#'
#' The default \code{seed = 42L} is the project's canonical seed and reproduces
#' the champion \code{premodis_training_pool_CHAMPION_RAND05.gpkg} EXACTLY
#' (714/714 random rows, verified by per-row \code{fire_uid} identity). NOTE:
#' this is a DIFFERENT seed from the shape-decorrelation step in
#' \code{\link{add_pool_shape_features}} (seed \code{1985L}); the two operations
#' are independent and must keep their own seeds.
#'
#' @param pool An \code{sf} / \code{data.frame} training pool with \code{source}
#'   and \code{class} columns.
#' @param random_to_burned_ratio Numeric (default \code{0.5}). Target random
#'   count as a multiple of the positive (burned) count.
#' @param seed Integer seed for the random selection (default \code{42L}, the
#'   canonical champion-pool seed).
#' @param random_source Source tag of the random background bucket.
#' @param positive_class Class tag of positives (default \code{"burned"}).
#' @return \code{pool} with the random bucket subsampled (other rows unchanged).
#' @export
subsample_random_background <- function(pool, random_to_burned_ratio = 0.5, seed = 42L,
                                        random_source = "random_burnable_background",
                                        positive_class = "burned") {
  if (!is.numeric(random_to_burned_ratio) || random_to_burned_ratio < 0) {
    stop("subsample_random_background(): 'random_to_burned_ratio' must be >= 0.", call. = FALSE)
  }
  d <- if (inherits(pool, "sf")) sf::st_drop_geometry(pool) else as.data.frame(pool)
  is_rand <- which(as.character(d$source) == random_source)
  n_rand  <- length(is_rand)
  n_burn  <- sum(as.character(d$class) == positive_class, na.rm = TRUE)
  target  <- as.integer(ceiling(n_burn * random_to_burned_ratio))
  if (n_rand <= target) return(pool)                    # nothing to subsample
  set.seed(seed)
  keep_rand <- sample(is_rand, size = target)
  drop_idx  <- setdiff(is_rand, keep_rand)
  pool[-drop_idx, , drop = FALSE]
}
