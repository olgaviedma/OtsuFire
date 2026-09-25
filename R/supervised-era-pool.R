# =============================================================================
# supervised-era-pool.R  (Phase 2B)
#
# Assemble a single TRAINABLE cross-year ERA training pool from per-year,
# visually-confirmed pools (packaged from assemble_era_pool.R). The VISUAL marks
# are consumed THROUGH apply_visual_validation() (the official Phase 2A entry):
#   * positives / background = train_features minus the VISUAL=0 (dropped) rows;
#   * hard-negs   = the VISUAL=1 artifact_hard CANDIDATES, pulled with full
#                   features from scoring_features and stamped class="unburned",
#                   source="artifact_hard", neg_type="artifact_hard_negative";
#   * omissions   = the VISUAL=0 CANDIDATES (real fires among the drops), routed
#                   to a separate layer with omit_reason="VISUAL0_real_fire".
# fire_uid and block_id are year-prefixed so they never collide across years;
# hard-neg folds are assigned reproducibly (seed 20260615).
# =============================================================================

# sf-safe row bind across parts with DIFFERENT column sets (fills missing with
# NA, aligns column order, preserves geometry). Returns NULL if all parts empty.
.of_rbind_sf <- function(parts) {
  parts <- Filter(function(z) !is.null(z) && nrow(z) > 0L, parts)
  if (!length(parts)) return(NULL)
  gcol <- attr(parts[[1]], "sf_column")
  allcols <- unique(unlist(lapply(parts, function(z) setdiff(names(z), attr(z, "sf_column")))))
  parts2 <- lapply(parts, function(z) {
    for (m in setdiff(allcols, names(z))) z[[m]] <- NA
    z[, c(allcols, attr(z, "sf_column")), drop = FALSE]
  })
  out <- do.call(rbind, parts2)
  if (!is.na(gcol)) attr(out, "sf_column") <- gcol
  out
}

# Assemble ONE year. tf/sc are sf (train_features / scoring_features); cl is the
# consolidated pool (path | sf | data.frame). Returns list(part, omissions).
.of_assemble_era_year <- function(year, train_features, scoring_features,
                                  consolidated_pool, hard_neg_fold_seed = 20260615L,
                                  visual_col = "VISUAL", on_invalid = "error") {
  y <- as.character(year)
  read_sf_if_path <- function(x, layer) if (is.character(x)) sf::st_read(x, layer = layer, quiet = TRUE) else x
  tf <- read_sf_if_path(train_features,   "train_features")
  sc <- read_sf_if_path(scoring_features, "scoring_features")

  # ---- VISUAL decisions THROUGH the official Phase 2A entry ----
  vv  <- apply_visual_validation(consolidated_pool, visual_col = visual_col, on_invalid = on_invalid)
  cld <- if (inherits(vv$pool, "sf")) sf::st_drop_geometry(vv$pool) else as.data.frame(vv$pool)
  act <- as.character(cld$visual_action)
  drop_uids  <- as.character(cld$fire_uid[act == "dropped_visual_reject"])
  conf_polys <- as.character(cld$poly_id[act == "promoted_artifact_hard"])
  omit_polys <- as.character(cld$poly_id[act == "omitted_real_fire"])

  # ---- positives / background (train_features minus VISUAL=0 drops) ----
  tfd <- sf::st_drop_geometry(tf)
  keep_row <- !(as.character(tfd$fire_uid) %in% drop_uids)
  positives  <- tf[as.character(tfd$class) == "burned"   & keep_row, , drop = FALSE]
  background <- tf[as.character(tfd$class) == "unburned" & keep_row, , drop = FALSE]

  # ---- confirmed hard negatives + omissions from scoring_features ----
  scd  <- sf::st_drop_geometry(sc)
  hard <- sc[as.character(scd$poly_id) %in% conf_polys, , drop = FALSE]
  if (nrow(hard) > 0L) {
    hard$class <- "unburned"; hard$source <- "artifact_hard"; hard$neg_type <- "artifact_hard_negative"
  }
  omit <- sc[as.character(scd$poly_id) %in% omit_polys, , drop = FALSE]

  # ---- folds: year-prefix blocks (leakage-safe); hard-negs get fresh unique
  # blocks + a reproducible distributed CV fold (mirrors the engine promotion).
  if (nrow(positives)  > 0L) positives$block_id  <- paste0(y, "_", as.character(positives$block_id))
  if (nrow(background) > 0L) background$block_id <- paste0(y, "_", as.character(background$block_id))
  if (nrow(hard) > 0L) {
    hard$block_id <- paste0(y, "_AH_", seq_len(nrow(hard)))
    kv <- sort(unique(c(as.integer(positives$fold_rep1), as.integer(background$fold_rep1))))
    kv <- kv[is.finite(kv)]; if (!length(kv)) kv <- 1:5
    set.seed(hard_neg_fold_seed)
    hard$fold_rep1 <- sample(kv[(seq_len(nrow(hard)) - 1L) %% length(kv) + 1L])
    hard$fold_rep2 <- sample(kv[(seq_len(nrow(hard)) - 1L) %% length(kv) + 1L])
  }

  # ---- year-prefix fire_uid + provenance ----
  pref <- function(z) { if (nrow(z) == 0L) return(z)
    z$fire_uid <- paste0(y, "_", as.character(z$fire_uid)); z$source_year <- as.integer(y); z }
  part <- .of_rbind_sf(list(pref(positives), pref(background), pref(hard)))
  omits <- NULL
  if (nrow(omit) > 0L) { o <- pref(omit); o$omit_reason <- "VISUAL0_real_fire"; omits <- o }
  list(part = part, omissions = omits)
}

#' Assemble a cross-year training pool from visually-validated per-year pools
#'
#' @description
#' Merges several years of visually-validated supervised pools into one
#' trainable cross-year pool. Use it when you want to train a single model
#' across multiple fire years instead of one model per year.
#'
#' For each year it keeps the positives and background from
#' \code{train_features} (minus the rows a reviewer marked \code{VISUAL = 0}),
#' promotes the \code{VISUAL = 1} artifact_hard candidates from
#' \code{scoring_features} as hard negatives, and routes the \code{VISUAL = 0}
#' candidates (real fires found among the drops) to a separate omissions layer.
#' \code{fire_uid} and \code{block_id} are year-prefixed so they never collide
#' across years, and hard-negative folds are assigned reproducibly. The visual
#' decisions are read through \code{\link{apply_visual_validation}}.
#'
#' @param years A list of per-year inputs. Each element is a list with
#'   \code{year} and either (a) \code{base_dir} (the standard
#'   \code{03_FEATURES/features_geometry.gpkg} +
#'   \code{01_POOLS/supervised_training_pool.gpkg} are read from it), or
#'   (b) explicit \code{train_features}, \code{scoring_features} and
#'   \code{consolidated_pool} (each a path, \code{sf}, or \code{data.frame}).
#' @param hard_neg_fold_seed Integer seed for hard-negative CV-fold assignment
#'   (default \code{20260615L}).
#' @param visual_col Visual-review column name (default \code{"VISUAL"}).
#' @param on_invalid Passed to \code{\link{apply_visual_validation}}
#'   (\code{"error"} or \code{"exclude"}).
#' @return A list with \code{pool} (the merged \code{sf} training pool),
#'   \code{omissions} (\code{sf} of \code{VISUAL = 0} real fires, or
#'   \code{NULL}) and \code{summary} (per-year row counts).
#'
#' @seealso
#' \code{\link{apply_visual_validation}},
#' \code{\link{export_visual_validation_pools}},
#' \code{\link{build_supervised_training_pools}}
#'
#' @examples
#' \dontrun{
#' era <- assemble_era_training_pool(list(
#'   list(year = 2020, base_dir = "results/2020/SUPERVISED/balanced"),
#'   list(year = 2021, base_dir = "results/2021/SUPERVISED/balanced")))
#' nrow(era$pool); nrow(era$omissions)
#' }
#' @export
assemble_era_training_pool <- function(years, hard_neg_fold_seed = 20260615L,
                                       visual_col = "VISUAL", on_invalid = "error") {
  if (!is.list(years) || !length(years)) stop("assemble_era_training_pool(): 'years' must be a non-empty list.", call. = FALSE)
  resolve <- function(yr) {
    if (!is.null(yr$base_dir)) {
      fg <- file.path(yr$base_dir, "03_FEATURES", "features_geometry.gpkg")
      cg <- file.path(yr$base_dir, "01_POOLS", "supervised_training_pool.gpkg")
      list(year = yr$year, train_features = fg, scoring_features = fg, consolidated_pool = cg)
    } else {
      list(year = yr$year, train_features = yr$train_features,
           scoring_features = yr$scoring_features, consolidated_pool = yr$consolidated_pool)
    }
  }
  parts <- list(); omits <- list(); summ <- list()
  for (yr in years) {
    r <- resolve(yr)
    res <- .of_assemble_era_year(r$year, r$train_features, r$scoring_features,
                                 r$consolidated_pool, hard_neg_fold_seed = hard_neg_fold_seed,
                                 visual_col = visual_col, on_invalid = on_invalid)
    parts[[as.character(r$year)]] <- res$part
    if (!is.null(res$omissions)) omits[[as.character(r$year)]] <- res$omissions
    pd <- if (is.null(res$part)) NULL else sf::st_drop_geometry(res$part)
    summ[[length(summ) + 1L]] <- data.frame(
      year       = as.integer(r$year),
      positives  = if (is.null(pd)) 0L else sum(as.character(pd$class) == "burned", na.rm = TRUE),
      rows       = if (is.null(pd)) 0L else nrow(pd),
      hard_neg   = if (is.null(pd)) 0L else sum(as.character(pd$source) == "artifact_hard", na.rm = TRUE),
      omissions  = if (is.null(res$omissions)) 0L else nrow(res$omissions))
  }
  list(pool      = .of_rbind_sf(parts),
       omissions = .of_rbind_sf(omits),
       summary   = do.call(rbind, summ))
}
