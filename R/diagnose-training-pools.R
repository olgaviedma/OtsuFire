# =============================================================================
# diagnose-training-pools.R
#
# Pre-fit DIAGNOSTIC for the supervised training pool. Statistically assesses
# whether the training-pool CATEGORIES (e.g. pool_source / training_label) are
# coherent and separable across the numeric features that the supervised model
# will later see. It is a *diagnostic only*: it NEVER trains or calls the final
# supervised model (no train_supervised_model / xgboost). The optional "probe
# classifier" is an explicitly-labelled separability probe (LDA), not the
# production model, and its predictions never feed anything downstream.
#
# Modules:
#   A. Composition & quality        .pool_composition()
#   B. Univariate separability      .pool_univariate()
#   C. Z-score fingerprint          .pool_zscore()
#   D. PCA (diagnostic only)        .pool_pca()
#   E. PERMANOVA (optional)         .pool_permanova()
#   F. RS separability (optional)   .pool_rs_separability()
#   G. Probe classifier (optional)  .pool_probe()
#   H. Correlation / redundancy     .pool_correlation()
#   + console summary               .pool_console_summary()
# =============================================================================

# ggplot2 aes() uses non-standard evaluation; these are column names referenced
# inside the plotting data.frames, not real R bindings. Declaring them silences
# "no visible binding for global variable" NOTEs without changing behaviour.
utils::globalVariables(c("group", "value", "feature", "median_z", "PC1", "PC2",
                         "f1", "f2", "r"))

# ---- small internal utilities ----------------------------------------------

.pool_msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(sprintf(...))
}

# case-insensitive membership against exclude list
.pool_in_excl <- function(nm, excl) {
  tolower(nm) %in% tolower(excl)
}

# Coerce any supported input to a plain data.frame WITHOUT geometry, never
# mutating the caller's object.
.pool_as_df <- function(pools) {
  if (inherits(pools, "sf")) {
    pools <- sf::st_drop_geometry(pools)
  }
  if (inherits(pools, "data.table")) {
    pools <- as.data.frame(pools)
  }
  if (inherits(pools, "tbl_df") || inherits(pools, "tbl")) {
    pools <- as.data.frame(pools)
  }
  as.data.frame(pools, stringsAsFactors = FALSE)
}

# ---- A. composition & quality -----------------------------------------------

.pool_composition <- function(df, group_col, label_col, visual_col, used_col,
                              features) {
  comp <- as.data.frame(table(group = df[[group_col]], useNA = "ifany"),
                        stringsAsFactors = FALSE)
  names(comp) <- c("group", "n")
  comp <- comp[order(-comp$n), , drop = FALSE]

  cross <- list()
  add_cross <- function(nm, col) {
    if (!is.null(col) && col %in% names(df)) {
      tt <- as.data.frame(table(group = df[[group_col]], level = df[[col]],
                                useNA = "ifany"), stringsAsFactors = FALSE)
      names(tt) <- c("group", col, "n")
      tt <- tt[tt$n > 0, , drop = FALSE]
      cross[[col]] <<- tt
    }
  }
  add_cross("label", label_col)
  add_cross("visual", visual_col)
  add_cross("used", used_col)

  # missingness per feature
  miss <- data.frame(
    feature = features,
    n               = nrow(df),
    n_missing       = vapply(features, function(f) sum(is.na(df[[f]])), integer(1)),
    stringsAsFactors = FALSE
  )
  miss$frac_missing <- miss$n_missing / pmax(miss$n, 1L)
  rownames(miss) <- NULL

  list(composition = comp, cross_tabs = cross, missingness = miss)
}

# high-correlation pairs (informational)
.pool_highcor <- function(df, features, thr = 0.9) {
  if (length(features) < 2) return(NULL)
  m <- as.matrix(df[, features, drop = FALSE])
  cc <- suppressWarnings(stats::cor(m, use = "pairwise.complete.obs"))
  out <- data.frame(feature1 = character(0), feature2 = character(0),
                    cor = numeric(0), stringsAsFactors = FALSE)
  for (i in seq_along(features)) {
    for (j in seq_len(i - 1L)) {
      r <- cc[i, j]
      if (!is.na(r) && abs(r) >= thr) {
        out <- rbind(out, data.frame(feature1 = features[i],
                                     feature2 = features[j],
                                     cor = r, stringsAsFactors = FALSE))
      }
    }
  }
  if (nrow(out) == 0) return(NULL)
  out[order(-abs(out$cor)), , drop = FALSE]
}

# ---- H. correlation / redundancy --------------------------------------------
# Full correlation matrix + high-|r| pairs + greedy redundancy clusters + a
# suggested NON-REDUNDANT reduced feature set. Features whose |r| >= thr are put
# in the same redundancy cluster (connected components of the |r| >= thr graph);
# the representative of each cluster is the feature with the HIGHEST importance
# (epsilon^2 from the univariate ranking), falling back to feature order when
# importance is missing/tied. The reduced set is all cluster representatives plus
# all singletons. Dependency-light (base stats::cor); never hard-fails.
.pool_correlation <- function(df, features, family_fun = .pool_feature_family,
                              importance = NULL, method = "spearman",
                              thr = 0.9) {
  tryCatch({
    features <- intersect(features, names(df))
    if (length(features) < 2) return(NULL)
    m  <- as.matrix(df[, features, drop = FALSE])
    cc <- suppressWarnings(stats::cor(m, use = "pairwise.complete.obs",
                                      method = method))
    if (is.null(cc) || all(is.na(cc))) return(NULL)
    np <- length(features)

    # family lookup (feature -> family)
    fam    <- family_fun(features)
    fam_of <- stats::setNames(fam$family, fam$feature)

    # importance lookup (feature -> epsilon^2; higher = better). Accept either a
    # named numeric vector or a feature_ranking-style data.frame.
    imp_of <- stats::setNames(rep(NA_real_, np), features)
    if (!is.null(importance)) {
      if (is.data.frame(importance) && "feature" %in% names(importance)) {
        ic <- intersect(c("epsilon_squared", "importance", "value"),
                        names(importance))[1]
        if (!is.na(ic)) {
          keep <- importance$feature %in% features
          imp_of[importance$feature[keep]] <- importance[[ic]][keep]
        }
      } else if (!is.null(names(importance))) {
        keep <- names(importance) %in% features
        imp_of[names(importance)[keep]] <- as.numeric(importance[keep])
      }
    }

    # high-correlation pairs (|r| >= thr), family-tagged
    pr <- data.frame(feature1 = character(0), feature2 = character(0),
                     cor = numeric(0), family1 = character(0),
                     family2 = character(0), stringsAsFactors = FALSE)
    for (i in seq_len(np)) {
      for (j in seq_len(i - 1L)) {
        r <- cc[i, j]
        if (!is.na(r) && abs(r) >= thr) {
          pr <- rbind(pr, data.frame(
            feature1 = features[i], feature2 = features[j], cor = r,
            family1 = unname(fam_of[features[i]]),
            family2 = unname(fam_of[features[j]]),
            stringsAsFactors = FALSE))
        }
      }
    }
    if (nrow(pr)) pr <- pr[order(-abs(pr$cor)), , drop = FALSE]
    rownames(pr) <- NULL

    # greedy clusters: union-find over the |r| >= thr graph (connected comps)
    parent <- seq_len(np)
    findr  <- function(i) { while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
    unite  <- function(a, b) { ra <- findr(a); rb <- findr(b)
      if (ra != rb) parent[ra] <<- rb }
    if (nrow(pr)) {
      for (k in seq_len(nrow(pr))) {
        unite(match(pr$feature1[k], features), match(pr$feature2[k], features))
      }
    }
    roots      <- vapply(seq_len(np), findr, integer(1))
    cluster_id <- match(roots, sort(unique(roots)))   # compact 1..K

    # representative per cluster = highest importance; NA -> -Inf; ties -> order
    rep_of <- character(np)
    for (cl in sort(unique(cluster_id))) {
      idx  <- which(cluster_id == cl)
      impv <- imp_of[features[idx]]; impv[is.na(impv)] <- -Inf
      best <- if (all(is.infinite(impv))) idx[1] else idx[which.max(impv)]
      rep_of[idx] <- features[best]
    }

    clusters <- data.frame(
      feature           = features,
      family            = unname(fam_of[features]),
      cluster           = cluster_id,
      representative    = rep_of,
      importance        = unname(imp_of[features]),
      is_representative = features == rep_of,
      dropped           = features != rep_of,
      stringsAsFactors  = FALSE)
    clusters <- clusters[order(clusters$cluster, !clusters$is_representative,
                               clusters$feature), , drop = FALSE]
    rownames(clusters) <- NULL

    # reduced set = representatives + singletons (in original feature order)
    reduced_set <- features[features %in% rep_of]

    # per-cluster representative table
    ucl <- sort(unique(cluster_id))
    representative <- data.frame(
      cluster        = ucl,
      representative = vapply(ucl, function(cl) rep_of[which(cluster_id == cl)[1]],
                              character(1)),
      n_members      = vapply(ucl, function(cl) sum(cluster_id == cl), integer(1)),
      stringsAsFactors = FALSE)
    representative$family     <- unname(fam_of[representative$representative])
    representative$importance <- unname(imp_of[representative$representative])

    # dropped features mapped to their representative (+ |r| to it)
    drop_idx <- which(features != rep_of)
    dropped  <- data.frame(
      feature        = features[drop_idx],
      representative = rep_of[drop_idx],
      cluster        = cluster_id[drop_idx],
      cor            = vapply(drop_idx, function(i) cc[i, match(rep_of[i], features)],
                              numeric(1)),
      family         = unname(fam_of[features[drop_idx]]),
      stringsAsFactors = FALSE)
    if (nrow(dropped)) dropped <- dropped[order(-abs(dropped$cor)), , drop = FALSE]
    rownames(dropped) <- NULL

    list(method = method, threshold = thr, matrix = cc, high_pairs = pr,
         clusters = clusters, representative = representative,
         reduced_set = reduced_set, dropped = dropped,
         n_features = np, n_reduced = length(reduced_set))
  }, error = function(e) NULL)
}

# ---- B. univariate separability ---------------------------------------------

# Manual effect sizes for Kruskal-Wallis, correctly labelled:
#   epsilon^2 = H / (n - 1)                    (rank epsilon-squared)
#   eta^2[H]  = (H - k + 1) / (n - k)          (rank eta-squared)
.kw_effects_manual <- function(H, n, k) {
  eps <- if (n > 1) H / (n - 1) else NA_real_
  eta <- if (n > k) (H - k + 1) / (n - k) else NA_real_
  eta <- max(eta, 0, na.rm = TRUE)
  list(epsilon_squared = eps, eta_squared = eta)
}

.pool_univariate <- function(df, group_col, features, run_pairwise,
                             have_effectsize) {
  g <- factor(df[[group_col]])
  res <- vector("list", length(features))
  for (i in seq_along(features)) {
    f <- features[i]
    x <- df[[f]]
    ok <- !is.na(x) & !is.na(g)
    xx <- x[ok]; gg <- droplevels(g[ok])
    k  <- nlevels(gg)
    n  <- length(xx)
    if (k < 2 || n < 3) {
      res[[i]] <- data.frame(feature = f, statistic = NA_real_, df = NA_real_,
                             p = NA_real_, epsilon_squared = NA_real_,
                             eta_squared = NA_real_, n = n, k = k,
                             stringsAsFactors = FALSE)
      next
    }
    kw  <- suppressWarnings(stats::kruskal.test(xx, gg))
    H   <- unname(kw$statistic)
    eff <- if (isTRUE(have_effectsize)) {
      e1 <- tryCatch(as.data.frame(effectsize::rank_epsilon_squared(xx, groups = gg)),
                     error = function(e) NULL)
      e2 <- tryCatch(as.data.frame(effectsize::rank_eta_squared(xx, groups = gg)),
                     error = function(e) NULL)
      eps <- if (!is.null(e1)) e1[[grep("epsilon", names(e1), ignore.case = TRUE)[1]]] else NA_real_
      eta <- if (!is.null(e2)) e2[[grep("eta",     names(e2), ignore.case = TRUE)[1]]] else NA_real_
      if (is.na(eps) || is.na(eta)) {
        man <- .kw_effects_manual(H, n, k)
        if (is.na(eps)) eps <- man$epsilon_squared
        if (is.na(eta)) eta <- man$eta_squared
      }
      list(epsilon_squared = eps, eta_squared = eta)
    } else {
      .kw_effects_manual(H, n, k)
    }
    res[[i]] <- data.frame(feature = f, statistic = H, df = unname(kw$parameter),
                           p = unname(kw$p.value),
                           epsilon_squared = eff$epsilon_squared,
                           eta_squared = eff$eta_squared,
                           n = n, k = k, stringsAsFactors = FALSE)
  }
  kruskal <- do.call(rbind, res)
  kruskal$p_adj <- stats::p.adjust(kruskal$p, method = "BH")
  kruskal <- kruskal[order(-kruskal$epsilon_squared,
                           kruskal$p_adj,
                           na.last = TRUE), , drop = FALSE]
  rownames(kruskal) <- NULL

  pairwise <- NULL
  if (isTRUE(run_pairwise)) {
    pairwise <- .pool_pairwise(df, group_col, features)
  }
  list(kruskal = kruskal, pairwise = pairwise)
}

# Vectorised-ish pairwise Wilcoxon across all group pairs x features.
.pool_pairwise <- function(df, group_col, features) {
  g <- factor(df[[group_col]])
  levs <- levels(g)
  if (length(levs) < 2) return(NULL)
  combs <- utils::combn(levs, 2)
  rows <- vector("list", ncol(combs) * length(features))
  idx <- 1L
  for (ci in seq_len(ncol(combs))) {
    g1 <- combs[1, ci]; g2 <- combs[2, ci]
    i1 <- which(g == g1); i2 <- which(g == g2)
    for (f in features) {
      x1 <- df[[f]][i1]; x2 <- df[[f]][i2]
      x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
      if (length(x1) < 1 || length(x2) < 1) {
        stat <- NA_real_; p <- NA_real_; mdiff <- NA_real_
      } else {
        wt <- suppressWarnings(stats::wilcox.test(x1, x2))
        stat <- unname(wt$statistic); p <- unname(wt$p.value)
        mdiff <- stats::median(x1) - stats::median(x2)
      }
      rows[[idx]] <- data.frame(group1 = g1, group2 = g2, feature = f,
                                statistic = stat, p = p,
                                median_diff = mdiff,
                                direction = if (is.na(mdiff)) NA_character_
                                            else if (mdiff > 0) paste0(g1, ">", g2)
                                            else if (mdiff < 0) paste0(g1, "<", g2)
                                            else "equal",
                                stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  pw <- do.call(rbind, rows)
  pw$p_adj <- stats::p.adjust(pw$p, method = "BH")
  pw <- pw[, c("group1", "group2", "feature", "statistic", "p", "p_adj",
               "median_diff", "direction")]
  pw <- pw[order(pw$p_adj, na.last = TRUE), , drop = FALSE]
  rownames(pw) <- NULL
  pw
}

# ---- C. z-score fingerprint -------------------------------------------------

.pool_feature_family <- function(features) {
  fam <- rep("other", length(features))
  lf <- tolower(features)
  fam[grepl("rbr_aw", lf)]                   <- "rbr_aw"
  fam[grepl("rbr", lf) & fam == "other"]     <- "rbr"
  fam[grepl("persist", lf)]                  <- "persist"
  fam[grepl("corine|cor_|^cor$", lf)]        <- "corine"
  fam[grepl("elev", lf)]                     <- "elev"
  fam[grepl("slope", lf)]                    <- "slope"
  fam[grepl("doy", lf)]                      <- "doy"
  fam[grepl("hotspot", lf)]                  <- "hotspot"
  fam[grepl("shape|inflation|polsby|compact", lf)] <- "shape"
  fam[grepl("area", lf)]                     <- "area"
  data.frame(feature = features, family = fam, stringsAsFactors = FALSE)
}

.pool_zscore <- function(df, group_col, features) {
  z <- scale(as.matrix(df[, features, drop = FALSE]))
  zdf <- as.data.frame(z)
  zdf[[group_col]] <- df[[group_col]]
  groups <- sort(unique(df[[group_col]]))
  out <- data.frame(group = character(0), feature = character(0),
                    median_z = numeric(0), stringsAsFactors = FALSE)
  for (grp in groups) {
    sub <- zdf[zdf[[group_col]] == grp & !is.na(zdf[[group_col]]), features, drop = FALSE]
    med <- vapply(sub, function(col) stats::median(col, na.rm = TRUE), numeric(1))
    out <- rbind(out, data.frame(group = as.character(grp), feature = features,
                                 median_z = unname(med), stringsAsFactors = FALSE))
  }
  fam <- .pool_feature_family(features)
  out <- merge(out, fam, by = "feature", all.x = TRUE, sort = FALSE)
  out <- out[, c("group", "feature", "family", "median_z")]
  rownames(out) <- NULL
  out
}

# ---- D. PCA (diagnostic only) -----------------------------------------------

# Imputation strategy: MEDIAN imputation (per feature). Documented in roxygen.
.pool_pca <- function(df, group_col, features, seed) {
  if (length(features) < 2) return(NULL)
  X <- as.matrix(df[, features, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    col <- X[, j]
    if (anyNA(col)) {
      m <- stats::median(col, na.rm = TRUE)
      col[is.na(col)] <- m
      X[, j] <- col
    }
  }
  keep <- apply(X, 2, function(c) {
    v <- stats::var(c, na.rm = TRUE)
    !is.na(v) && v > 0
  })
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2) return(NULL)

  pc <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  sdev <- pc$sdev
  expl <- sdev^2 / sum(sdev^2)
  npc  <- length(sdev)

  scores <- as.data.frame(pc$x)
  scores$group <- df[[group_col]]

  loadings <- as.data.frame(pc$rotation)
  loadings$feature <- rownames(pc$rotation)
  rownames(loadings) <- NULL

  explained <- data.frame(
    pc = paste0("PC", seq_len(npc)),
    sdev = sdev,
    prop_var = expl,
    cum_var = cumsum(expl),
    stringsAsFactors = FALSE
  )

  top_pc1 <- loadings[order(-abs(loadings$PC1)), c("feature", "PC1")]
  top_pc2 <- if (npc >= 2) loadings[order(-abs(loadings$PC2)), c("feature", "PC2")] else NULL

  # centroids in PC1-PC2 and in full scaled space
  cent_dist <- .pool_centroid_dist(scores, X, df[[group_col]],
                                    pc_dims = min(2L, npc))

  list(
    scores = scores,
    loadings = loadings,
    explained = explained,
    top_pc1 = utils::head(top_pc1, 10),
    top_pc2 = if (!is.null(top_pc2)) utils::head(top_pc2, 10) else NULL,
    centroid_distances = cent_dist,
    imputation = "median"
  )
}

.pool_centroid_dist <- function(scores, X_scaled, group, pc_dims = 2L) {
  groups <- sort(unique(group[!is.na(group)]))
  if (length(groups) < 2) return(NULL)
  pc_cols <- paste0("PC", seq_len(pc_dims))
  Xs <- scale(X_scaled)
  cent_pc <- lapply(groups, function(g) {
    colMeans(scores[scores$group == g & !is.na(scores$group), pc_cols, drop = FALSE],
             na.rm = TRUE)
  })
  cent_full <- lapply(groups, function(g) {
    colMeans(Xs[group == g & !is.na(group), , drop = FALSE], na.rm = TRUE)
  })
  names(cent_pc) <- names(cent_full) <- groups
  cb <- utils::combn(groups, 2)
  out <- data.frame(group1 = character(0), group2 = character(0),
                    dist_pc12 = numeric(0), dist_full = numeric(0),
                    stringsAsFactors = FALSE)
  for (ci in seq_len(ncol(cb))) {
    g1 <- cb[1, ci]; g2 <- cb[2, ci]
    d_pc <- sqrt(sum((cent_pc[[g1]] - cent_pc[[g2]])^2))
    d_fl <- sqrt(sum((cent_full[[g1]] - cent_full[[g2]])^2))
    out <- rbind(out, data.frame(group1 = g1, group2 = g2,
                                 dist_pc12 = d_pc, dist_full = d_fl,
                                 stringsAsFactors = FALSE))
  }
  out <- out[order(-out$dist_full), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- E. PERMANOVA (optional) ------------------------------------------------

.pool_permanova <- function(df, group_col, features, seed, verbose) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    .pool_msg(verbose, "PERMANOVA skipped: package 'vegan' not installed.")
    return(list(available = FALSE, reason = "vegan not installed",
                adonis2 = NULL, dispersion = NULL))
  }
  X <- as.matrix(df[, features, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    col <- X[, j]; if (anyNA(col)) col[is.na(col)] <- stats::median(col, na.rm = TRUE)
    X[, j] <- col
  }
  Xs <- scale(X)
  g  <- factor(df[[group_col]])
  d  <- stats::dist(Xs)
  set.seed(seed)
  ad <- tryCatch(vegan::adonis2(d ~ g), error = function(e) NULL)
  bd <- tryCatch({
    bb <- vegan::betadisper(d, g)
    vegan::permutest(bb)
  }, error = function(e) NULL)
  list(available = TRUE, reason = NA_character_,
       adonis2 = if (!is.null(ad)) as.data.frame(ad) else NULL,
       dispersion = bd)
}

# ---- F. RS separability (optional) ------------------------------------------

# Bhattacharyya, Jeffries-Matusita, transformed divergence, Mahalanobis between
# group pairs, computed from class means/covariances (Gaussian assumption).
.pair_rs_metrics <- function(m1, m2, S1, S2) {
  dm <- m1 - m2
  Sm <- (S1 + S2) / 2
  Sm_inv <- tryCatch(solve(Sm), error = function(e) MASS::ginv(Sm))
  term1 <- as.numeric(0.125 * t(dm) %*% Sm_inv %*% dm)
  d1 <- det(S1); d2 <- det(S2); dm_ <- det(Sm)
  term2 <- if (d1 > 0 && d2 > 0 && dm_ > 0) 0.5 * log(dm_ / sqrt(d1 * d2)) else NA_real_
  bhatt <- term1 + ifelse(is.na(term2), 0, term2)
  jm <- 2 * (1 - exp(-bhatt))
  # transformed divergence
  S1i <- tryCatch(solve(S1), error = function(e) MASS::ginv(S1))
  S2i <- tryCatch(solve(S2), error = function(e) MASS::ginv(S2))
  div <- 0.5 * sum(diag((S1 - S2) %*% (S2i - S1i))) +
         0.5 * as.numeric(t(dm) %*% (S1i + S2i) %*% dm)
  td <- 2 * (1 - exp(-div / 8))
  # Mahalanobis (pooled)
  maha <- sqrt(as.numeric(t(dm) %*% Sm_inv %*% dm))
  c(bhattacharyya = bhatt, jeffries_matusita = jm,
    transformed_divergence = td, mahalanobis = maha)
}

.pool_rs_separability <- function(df, group_col, features, verbose) {
  if (length(features) < 2) {
    return(list(available = FALSE, reason = "need >= 2 features",
                pairwise = NULL))
  }
  has_pkg <- requireNamespace("spatialEco", quietly = TRUE) ||
             requireNamespace("varSel", quietly = TRUE)
  # We can always compute manually with MASS::ginv as fallback; flag which.
  g <- factor(df[[group_col]])
  levs <- levels(g)
  if (length(levs) < 2) {
    return(list(available = FALSE, reason = "need >= 2 groups", pairwise = NULL))
  }
  X <- df[, features, drop = FALSE]
  stats_by_grp <- lapply(levs, function(l) {
    sub <- as.matrix(X[g == l, , drop = FALSE])
    sub <- sub[stats::complete.cases(sub), , drop = FALSE]
    if (nrow(sub) <= length(features)) return(NULL)
    list(m = colMeans(sub), S = stats::cov(sub), n = nrow(sub))
  })
  names(stats_by_grp) <- levs
  cb <- utils::combn(levs, 2)
  rows <- list(); idx <- 1L
  for (ci in seq_len(ncol(cb))) {
    g1 <- cb[1, ci]; g2 <- cb[2, ci]
    s1 <- stats_by_grp[[g1]]; s2 <- stats_by_grp[[g2]]
    if (is.null(s1) || is.null(s2)) next
    mm <- tryCatch(.pair_rs_metrics(s1$m, s2$m, s1$S, s2$S),
                   error = function(e) c(bhattacharyya = NA, jeffries_matusita = NA,
                                         transformed_divergence = NA, mahalanobis = NA))
    rows[[idx]] <- data.frame(group1 = g1, group2 = g2,
                              bhattacharyya = mm[["bhattacharyya"]],
                              jeffries_matusita = mm[["jeffries_matusita"]],
                              transformed_divergence = mm[["transformed_divergence"]],
                              mahalanobis = mm[["mahalanobis"]],
                              stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  if (length(rows) == 0) {
    return(list(available = FALSE, reason = "insufficient per-group sample size",
                pairwise = NULL))
  }
  pw <- do.call(rbind, rows)
  pw <- pw[order(-pw$jeffries_matusita), , drop = FALSE]
  rownames(pw) <- NULL
  engine <- if (has_pkg) "manual (Gaussian); spatialEco/varSel available" else
            "manual (Gaussian); spatialEco/varSel not installed"
  list(available = TRUE, reason = NA_character_, engine = engine, pairwise = pw)
}

# ---- G. probe classifier (optional) -----------------------------------------
# Lightweight LDA separability probe with grouped (leave-fire-out / leave-year-
# out) CV. EXPLICITLY NOT the final supervised model.

.pool_probe <- function(df, group_col, features, fold_key, seed, verbose) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    .pool_msg(verbose, "Probe skipped: 'MASS' (LDA) not installed.")
    return(list(available = FALSE, reason = "MASS not installed",
                method = NA, cv = "none", balanced_accuracy = NA_real_,
                macro_f1 = NA_real_, confusion = NULL))
  }
  d <- df[, c(group_col, features), drop = FALSE]
  d[[group_col]] <- factor(d[[group_col]])
  cc <- stats::complete.cases(d)
  d <- d[cc, , drop = FALSE]
  key <- if (!is.null(fold_key)) fold_key[cc] else NULL
  if (nlevels(droplevels(d[[group_col]])) < 2 || nrow(d) < 20) {
    return(list(available = FALSE, reason = "too few complete cases/groups",
                method = "LDA", cv = "none", balanced_accuracy = NA_real_,
                macro_f1 = NA_real_, confusion = NULL))
  }
  d[[group_col]] <- droplevels(d[[group_col]])
  lvls <- levels(d[[group_col]])

  # build folds: grouped if key available, else random 5-fold
  set.seed(seed)
  if (!is.null(key) && length(unique(key[!is.na(key)])) >= 2) {
    ukey <- unique(key)
    fold_of <- stats::setNames(sample(rep_len(seq_len(min(5L, length(ukey))),
                                              length(ukey))), as.character(ukey))
    folds <- fold_of[as.character(key)]
    cv_type <- "grouped (leave-fire/year-out)"
  } else {
    folds <- sample(rep_len(seq_len(5L), nrow(d)))
    cv_type <- "random 5-fold"
  }
  nf <- length(unique(folds))
  preds <- rep(NA_character_, nrow(d))
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- folds == k
    if (length(unique(d[[group_col]][tr])) < 2) next
    fit <- tryCatch(MASS::lda(stats::as.formula(paste(group_col, "~ .")),
                              data = d[tr, , drop = FALSE]),
                    error = function(e) NULL)
    if (is.null(fit)) next
    pr <- tryCatch(as.character(stats::predict(fit, d[te, , drop = FALSE])$class),
                   error = function(e) NULL)
    if (!is.null(pr)) preds[te] <- pr
  }
  obs <- as.character(d[[group_col]])
  ok <- !is.na(preds)
  obs <- obs[ok]; preds <- preds[ok]
  cm <- table(observed = factor(obs, lvls), predicted = factor(preds, lvls))
  # per-class recall -> balanced accuracy; macro F1
  recalls <- numeric(0); f1s <- numeric(0)
  for (l in lvls) {
    tp <- cm[l, l]
    fn <- sum(cm[l, ]) - tp
    fp <- sum(cm[, l]) - tp
    rec <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    f1 <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
            2 * prec * rec / (prec + rec) else 0
    recalls <- c(recalls, rec); f1s <- c(f1s, f1)
  }
  list(available = TRUE, reason = NA_character_, method = "LDA probe (NOT final model)",
       cv = cv_type, n = length(obs),
       balanced_accuracy = mean(recalls, na.rm = TRUE),
       macro_f1 = mean(f1s, na.rm = TRUE),
       confusion = as.data.frame.matrix(cm))
}

# ---- console summary --------------------------------------------------------

# burned / hard-negative flexible detection
.pool_detect_burned <- function(levs) {
  levs[grepl("burned|keep|positive|high_confidence", levs, ignore.case = TRUE)]
}
.pool_detect_hardneg <- function(levs) {
  levs[grepl("hard|artifact|artifact_hard|hard_negative", levs, ignore.case = TRUE)]
}

.pool_console_summary <- function(out, top_n_features, print_pairwise_top_n) {
  L <- character(0)
  add <- function(...) L[[length(L) + 1L]] <<- sprintf(...)

  add("=== diagnose_training_pools(): pool separability diagnostic ===")
  add("Group column: %s   |   n rows: %d   |   features: %d",
      out$settings$group_col, out$settings$n_rows, length(out$feature_info$features_used))
  add("")

  add("-- A. Composition --")
  comp <- out$summary$composition
  for (i in seq_len(nrow(comp))) {
    add("   %-40s n = %d", as.character(comp$group[i]), comp$n[i])
  }
  if (length(out$summary$dropped_small_groups)) {
    add("   (dropped from tests, n < %d: %s)",
        out$settings$min_n_per_group,
        paste(out$summary$dropped_small_groups, collapse = ", "))
  }
  add("")

  add("-- B. Top separating features (by epsilon^2, BH-FDR) --")
  if (!is.null(out$kruskal)) {
    kk <- utils::head(out$kruskal, top_n_features)
    for (i in seq_len(nrow(kk))) {
      add("   %-16s eps^2 = %.3f  eta^2 = %.3f  p_adj = %.3g",
          kk$feature[i], kk$epsilon_squared[i], kk$eta_squared[i], kk$p_adj[i])
    }
  }
  add("")

  # burned vs hard-neg contrast
  add("-- Burned vs hard-negative contrast --")
  levs <- as.character(out$summary$composition$group)
  bcat <- .pool_detect_burned(levs); hcat <- .pool_detect_hardneg(levs)
  if (length(bcat) == 0 || length(hcat) == 0) {
    add("   Contrast unavailable: burned (%s) or hard-neg (%s) category not present; skipped.",
        if (length(bcat)) paste(bcat, collapse = "/") else "none",
        if (length(hcat)) paste(hcat, collapse = "/") else "none")
  } else if (!is.null(out$pairwise)) {
    sel <- out$pairwise[
      (out$pairwise$group1 %in% bcat & out$pairwise$group2 %in% hcat) |
      (out$pairwise$group1 %in% hcat & out$pairwise$group2 %in% bcat), , drop = FALSE]
    sel <- utils::head(sel[order(sel$p_adj), , drop = FALSE], print_pairwise_top_n)
    if (nrow(sel) == 0) add("   No pairwise rows for this contrast.")
    for (i in seq_len(nrow(sel))) {
      add("   %s vs %s : %-14s p_adj = %.3g  (%s)",
          sel$group1[i], sel$group2[i], sel$feature[i], sel$p_adj[i], sel$direction[i])
    }
  } else {
    add("   (pairwise not run)")
  }
  add("")

  add("-- D. PCA --")
  if (!is.null(out$pca)) {
    ev <- out$pca$explained
    p1 <- ev$prop_var[1] * 100
    p2 <- if (nrow(ev) >= 2) ev$prop_var[2] * 100 else 0
    add("   PC1 = %.1f%%   PC2 = %.1f%%   sum = %.1f%%   (imputation: %s)",
        p1, p2, p1 + p2, out$pca$imputation)
    drv <- utils::head(out$pca$top_pc1, 3)
    add("   PC1 drivers: %s", paste(sprintf("%s(%.2f)", drv$feature, drv$PC1), collapse = ", "))
    cd <- out$pca$centroid_distances
    if (!is.null(cd)) {
      top <- utils::head(cd, 3)
      for (i in seq_len(nrow(top))) {
        add("   largest separation: %s vs %s  full-space d = %.2f (PC1-2 d = %.2f)",
            top$group1[i], top$group2[i], top$dist_full[i], top$dist_pc12[i])
      }
      if (length(bcat) && length(hcat)) {
        bh <- cd[(cd$group1 %in% bcat & cd$group2 %in% hcat) |
                 (cd$group1 %in% hcat & cd$group2 %in% bcat), , drop = FALSE]
        if (nrow(bh)) {
          sep <- if (bh$dist_full[1] > stats::median(cd$dist_full)) "remains separable" else "weakly separated"
          add("   Burned vs hard-neg in full scaled space: %s (d = %.2f).",
              sep, bh$dist_full[1])
        }
      }
    }
  } else {
    add("   PCA not run / not enough features.")
  }
  add("")

  add("-- H. Correlation / redundancy --")
  if (!is.null(out$correlation)) {
    cr <- out$correlation
    add("   method = %s   threshold |r| >= %.2f", cr$method, cr$threshold)
    add("   high-correlation pairs: %d   redundant features dropped: %d",
        nrow(cr$high_pairs), nrow(cr$dropped))
    add("   suggested non-redundant set: %d -> %d features",
        cr$n_features, cr$n_reduced)
  } else {
    add("   (correlation module not run / < 2 features)")
  }
  add("")

  add("-- Interpretation --")
  topo_feats <- c("elev", "slope")
  if (!is.null(out$kruskal)) {
    top3 <- utils::head(out$kruskal$feature, 3)
    fam3 <- tolower(top3)
    if (any(grepl(paste(topo_feats, collapse = "|"), fam3))) {
      add("   * Topographic features dominate separation -> possible SAMPLING BIAS")
      add("     (groups may differ by terrain, not by burn signal). Inspect pool geography.")
    } else {
      add("   * Separation is driven by burn-signal features (rbr/persist), as desired.")
    }
  }
  if (!is.null(out$pca) && nrow(out$pca$explained) >= 2 &&
      sum(out$pca$explained$prop_var[1:2]) < 0.4) {
    add("   * Low PC1+PC2 variance: structure is high-dimensional; single axes won't tell the story.")
  }
  add("   * This is a PRE-FIT diagnostic only; no supervised model was trained.")

  paste(unlist(L), collapse = "\n")
}

# ---- plotting (guarded) -----------------------------------------------------

.pool_subsample <- function(df, group_col, max_per_group, seed) {
  set.seed(seed)
  groups <- unique(df[[group_col]])
  keep <- unlist(lapply(groups, function(g) {
    idx <- which(df[[group_col]] == g)
    if (length(idx) > max_per_group) sample(idx, max_per_group) else idx
  }))
  df[sort(keep), , drop = FALSE]
}

.pool_colors <- function(groups) {
  if (exists("OF_COLORS") && !is.null(names(OF_COLORS)) &&
      all(as.character(groups) %in% names(OF_COLORS))) {
    return(OF_COLORS[as.character(groups)])
  }
  stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"),
                  as.character(groups))
}

.pool_make_plots <- function(out, df, group_col, features, fig_dir, prefix,
                             seed, verbose) {
  paths <- list()
  has_gg <- requireNamespace("ggplot2", quietly = TRUE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  sub <- .pool_subsample(df, group_col, 8000L, seed)
  cols <- .pool_colors(sort(unique(df[[group_col]])))

  # 1. top-feature boxplots
  top_feats <- if (!is.null(out$kruskal))
    utils::head(stats::na.omit(out$kruskal$feature), min(6L, length(features))) else
    utils::head(features, 6L)
  p1 <- file.path(fig_dir, paste0(prefix, "_top_feature_boxplots.png"))
  ok1 <- tryCatch({
    if (has_gg) {
      long <- do.call(rbind, lapply(top_feats, function(f)
        data.frame(group = sub[[group_col]], feature = f, value = sub[[f]],
                   stringsAsFactors = FALSE)))
      g <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = value, fill = group)) +
        ggplot2::geom_boxplot(outlier.size = 0.3) +
        ggplot2::facet_wrap(~feature, scales = "free_y") +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                       legend.position = "none") +
        ggplot2::labs(title = "Top separating features by pool group",
                      x = NULL, y = NULL)
      ggplot2::ggsave(p1, g, width = 9, height = 6, dpi = 120)
    } else {
      grDevices::png(p1, width = 1100, height = 750)
      op <- graphics::par(mfrow = c(2, 3), mar = c(7, 4, 2, 1))
      on.exit(graphics::par(op), add = TRUE)
      for (f in top_feats)
        graphics::boxplot(sub[[f]] ~ sub[[group_col]], main = f, las = 2,
                          col = cols, xlab = "", ylab = "")
      grDevices::dev.off()
    }
    TRUE
  }, error = function(e) { .pool_msg(verbose, "boxplot failed: %s", conditionMessage(e)); FALSE })
  if (ok1) paths$top_feature_boxplots <- p1

  # 2. z-score fingerprint heatmap
  if (!is.null(out$zscore_fingerprint)) {
    p2 <- file.path(fig_dir, paste0(prefix, "_zscore_fingerprint_heatmap.png"))
    ok2 <- tryCatch({
      zf <- out$zscore_fingerprint
      if (has_gg) {
        g <- ggplot2::ggplot(zf, ggplot2::aes(x = feature, y = group, fill = median_z)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)) +
          ggplot2::labs(title = "Median z-score fingerprint", x = NULL, y = NULL)
        ggplot2::ggsave(p2, g, width = 9, height = 5, dpi = 120)
      } else {
        m <- stats::reshape(zf[, c("group", "feature", "median_z")],
                            idvar = "group", timevar = "feature", direction = "wide")
        rn <- m$group; m$group <- NULL
        mm <- as.matrix(m); rownames(mm) <- rn
        grDevices::png(p2, width = 1000, height = 600)
        graphics::image(t(mm), axes = FALSE,
                        col = grDevices::hcl.colors(20, "Blue-Red"))
        grDevices::dev.off()
      }
      TRUE
    }, error = function(e) { .pool_msg(verbose, "heatmap failed: %s", conditionMessage(e)); FALSE })
    if (ok2) paths$zscore_fingerprint_heatmap <- p2
  }

  # 3 & 4 PCA plots
  if (!is.null(out$pca)) {
    sc <- out$pca$scores
    sc_sub <- .pool_subsample(sc, "group", 8000L, seed)
    p3 <- file.path(fig_dir, paste0(prefix, "_pca_scatter.png"))
    ok3 <- tryCatch({
      if (has_gg) {
        g <- ggplot2::ggplot(sc_sub, ggplot2::aes(x = PC1, y = PC2, color = group)) +
          ggplot2::geom_point(size = 0.5, alpha = 0.5) +
          ggplot2::scale_color_manual(values = .pool_colors(sort(unique(sc$group)))) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::labs(title = "PCA scores by pool group")
        ggplot2::ggsave(p3, g, width = 7, height = 6, dpi = 120)
      } else {
        gcol <- .pool_colors(sort(unique(sc$group)))
        grDevices::png(p3, width = 800, height = 700)
        graphics::plot(sc_sub$PC1, sc_sub$PC2, col = gcol[as.character(sc_sub$group)],
                       pch = 19, cex = 0.4, xlab = "PC1", ylab = "PC2",
                       main = "PCA scores by pool group")
        graphics::legend("topright", legend = names(gcol), col = gcol, pch = 19, cex = 0.7)
        grDevices::dev.off()
      }
      TRUE
    }, error = function(e) { .pool_msg(verbose, "pca scatter failed: %s", conditionMessage(e)); FALSE })
    if (ok3) paths$pca_scatter <- p3

    p4 <- file.path(fig_dir, paste0(prefix, "_pca_loadings.png"))
    ok4 <- tryCatch({
      ld <- out$pca$loadings
      if (has_gg && all(c("PC1", "PC2") %in% names(ld))) {
        g <- ggplot2::ggplot(ld, ggplot2::aes(x = PC1, y = PC2)) +
          ggplot2::geom_segment(ggplot2::aes(xend = PC1, yend = PC2, x = 0, y = 0),
                                arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))) +
          ggplot2::geom_text(ggplot2::aes(label = feature), size = 3, vjust = -0.4) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::labs(title = "PCA loadings (PC1 vs PC2)")
        ggplot2::ggsave(p4, g, width = 7, height = 6, dpi = 120)
        TRUE
      } else FALSE
    }, error = function(e) FALSE)
    if (ok4) paths$pca_loadings <- p4
  }

  # 5. correlation heatmap (ggplot2 only; skipped gracefully if ggplot2 absent)
  if (has_gg && !is.null(out$correlation) && !is.null(out$correlation$matrix)) {
    p5 <- file.path(fig_dir, paste0(prefix, "_correlation_heatmap.png"))
    ok5 <- tryCatch({
      cm <- out$correlation$matrix
      fe <- rownames(cm)
      long <- data.frame(
        f1 = factor(rep(fe, times = length(fe)), levels = fe),
        f2 = factor(rep(fe, each  = length(fe)), levels = rev(fe)),
        r  = as.vector(cm), stringsAsFactors = FALSE)
      g <- ggplot2::ggplot(long, ggplot2::aes(x = f1, y = f2, fill = r)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", limits = c(-1, 1)) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = sprintf("Feature correlation (%s)",
                                      out$correlation$method),
                      x = NULL, y = NULL, fill = "r")
      ggplot2::ggsave(p5, g, width = 8, height = 7, dpi = 120)
      TRUE
    }, error = function(e) {
      .pool_msg(verbose, "corr heatmap failed: %s", conditionMessage(e)); FALSE })
    if (ok5) paths$correlation_heatmap <- p5
  }

  paths
}

# =============================================================================
# MAIN EXPORTED FUNCTION
# =============================================================================

#' Diagnose supervised training-pool separability (pre-fit)
#'
#' Statistically assesses whether the categories of a supervised training pool
#' (e.g. \code{pool_source} or \code{training_label}) are coherent and separable
#' across the numeric features the supervised model will later use, \emph{before}
#' any model is fit. It is a \strong{diagnostic only}: it never trains or calls
#' the final supervised model (no \code{xgboost}).
#' The optional probe classifier (module G) is an explicitly-labelled separability
#' probe (LDA), not the production model, and its predictions are never used
#' downstream.
#'
#' @details
#' Modules:
#' \describe{
#'   \item{A. Composition & quality}{counts per group; cross-tabs with label /
#'     visual / used columns when supplied; per-feature missingness; zero-variance
#'     report; optional high-correlation pairs.}
#'   \item{B. Univariate separability}{per-feature Kruskal-Wallis across the group
#'     column with BH-FDR. Effect sizes use \code{effectsize::rank_epsilon_squared}
#'     (epsilon^2 = H/(n-1)) and \code{effectsize::rank_eta_squared}
#'     (eta^2\[H\] = (H-k+1)/(n-k)); if \pkg{effectsize} is unavailable the same two
#'     quantities are computed manually and reported under the SAME column names
#'     (\code{epsilon_squared}, \code{eta_squared}). Features ranked primarily by
#'     epsilon^2. If \code{run_pairwise}, pairwise Wilcoxon between categories with
#'     BH-FDR.}
#'   \item{C. Z-score fingerprint}{global-scaled median z per group x feature
#'     (+ feature family grouping), optional heatmap.}
#'   \item{D. PCA (diagnostic only)}{scaled features with \strong{median
#'     imputation} of missing values; returns scores,
#'     loadings, explained variance, top PC1/PC2 loadings, and centroid distances
#'     between groups in PC1-PC2 and in the full scaled space. Never feeds the
#'     final model.}
#'   \item{E. PERMANOVA (optional)}{\code{vegan::adonis2} + dispersion via
#'     \code{vegan::betadisper}; skipped gracefully if \pkg{vegan} absent.}
#'   \item{F. RS separability (optional)}{pairwise Bhattacharyya, Jeffries-Matusita,
#'     transformed divergence and Mahalanobis (Gaussian assumption); uses
#'     \pkg{spatialEco}/\pkg{varSel} if present, otherwise a built-in computation.}
#'   \item{G. Probe classifier (optional)}{lightweight LDA separability probe with
#'     grouped (leave-fire/year-out) CV when a fire/year key is present; returns
#'     balanced accuracy, macro-F1 and a confusion matrix. Explicitly NOT the
#'     final model.}
#'   \item{H. Correlation / redundancy (optional)}{a feature correlation matrix
#'     (\code{corr_method}, default Spearman -- robust to the heavy-tailed
#'     RBR/area distributions), the high-correlation pairs
#'     (\eqn{|r| \ge} \code{corr_threshold}, family-tagged), greedy redundancy
#'     clusters (connected components of the high-correlation graph), a
#'     representative per cluster chosen as the feature with the HIGHEST
#'     importance (epsilon^2 from module B; ties / missing importance fall back to
#'     feature order), and a suggested NON-REDUNDANT \code{reduced_set} = all
#'     cluster representatives plus all singletons. Dependency-light (base
#'     \code{stats::cor}); never hard-fails. When \code{output_dir} is set it
#'     writes \code{*_15_correlation_high_pairs.csv},
#'     \code{*_16_correlation_clusters.csv} and
#'     \code{*_17_correlation_reduced_set.csv}; an optional correlation heatmap
#'     PNG is added when \code{make_plots} and \pkg{ggplot2} are available.}
#' }
#'
#' @param pools A \code{data.frame}, \code{data.table}, \code{tibble} or \code{sf}
#'   object holding the training pool. \code{sf} geometry is dropped internally;
#'   the caller's object is never mutated.
#' @param group_col Column defining the categories to test (default
#'   \code{"pool_source"}).
#' @param label_col,visual_col,used_col Optional columns for cross-tabulation
#'   (and \code{used_col} for the \code{only_used} filter).
#' @param feature_cols Optional explicit numeric feature columns. If \code{NULL},
#'   numeric columns are auto-detected and \code{exclude_cols} (matched
#'   case-insensitively), zero-variance and grouping/label/ID columns removed.
#' @param exclude_cols Columns to exclude from auto-detected features
#'   (case-insensitive).
#' @param only_used If \code{TRUE} and \code{used_col} given, restrict to rows
#'   where that column is truthy.
#' @param visual_filter Optional values; if given with \code{visual_col}, keep
#'   only rows whose visual value is in this set.
#' @param min_n_per_group Groups with fewer rows are dropped from tests but
#'   reported.
#' @param top_n_features Number of top features shown in the console summary.
#' @param print_pairwise_top_n Number of pairwise contrasts shown in the summary.
#' @param run_pairwise,run_pca,run_permanova,run_rs_separability,run_probe_classifier
#'   Module toggles.
#' @param run_correlation If \code{TRUE} (default), compute the correlation /
#'   redundancy module (module H) whenever there are at least two numeric
#'   features; degrades gracefully (returns \code{NULL}) on degenerate input.
#' @param corr_method Correlation method passed to \code{stats::cor} for module H
#'   (default \code{"spearman"}; \code{"pearson"}/\code{"kendall"} also accepted).
#' @param corr_threshold Absolute correlation threshold (default \code{0.9}) at or
#'   above which two features are treated as redundant (placed in the same
#'   cluster) in module H.
#' @param balance_n Optional per-group cap; groups are downsampled to this size
#'   (seeded) before stats.
#' @param n_repeats Reserved repeat count for stochastic modules (currently 1).
#' @param output_dir If given, CSVs / figures / console summary are written here.
#' @param prefix Filename prefix for outputs.
#' @param make_plots Whether to attempt figures (guarded by \pkg{ggplot2}; falls
#'   back to base graphics; never hard-fails).
#' @param seed RNG seed (set once at top, passed to stochastic parts).
#' @param verbose Print progress messages.
#'
#' @return An object of class \code{otsufire_pool_diagnostics}: a list with
#'   \code{summary}, \code{settings}, \code{feature_info}, \code{missingness},
#'   \code{zero_variance}, \code{kruskal}, \code{pairwise}, \code{feature_ranking},
#'   \code{correlation}, \code{zscore_fingerprint}, \code{pca}, \code{permanova},
#'   \code{rs_separability}, \code{probe_classifier}, \code{plot_paths} and
#'   \code{output_dir}. Modules that are skipped are \code{NULL}. The
#'   \code{correlation} component (module H) is a list with \code{matrix},
#'   \code{high_pairs}, \code{clusters}, \code{representative}, \code{reduced_set}
#'   (the suggested non-redundant feature set) and \code{dropped}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' df <- data.frame(
#'   pool_source    = sample(c("otsu", "random", "high_confidence_keep"), n, TRUE),
#'   training_label = sample(c("burned", "unburned"), n, TRUE),
#'   used_in_training = sample(c(TRUE, FALSE), n, TRUE),
#'   VISUAL         = sample(c(0, 1), n, TRUE),
#'   rbr_med        = rnorm(n, 200, 60),
#'   rbr_aw_med     = rnorm(n, 180, 55),
#'   persist_ratio  = runif(n),
#'   persist_delta  = rnorm(n),
#'   area_ha        = rexp(n, 0.1),
#'   doy_iqr        = rpois(n, 20)
#' )
#' diagnose_training_pools(df, group_col = "pool_source")
#' diagnose_training_pools(df, group_col = "training_label")
#' diagnose_training_pools(df, only_used = TRUE, used_col = "used_in_training")
#' diagnose_training_pools(df, visual_col = "VISUAL", visual_filter = c(0, 1))
#' diagnose_training_pools(df, output_dir = tempfile("pooldiag_"))
#' }
#'
#' @export
diagnose_training_pools <- function(
  pools,
  group_col       = "pool_source",
  label_col       = NULL,
  visual_col      = NULL,
  used_col        = NULL,
  feature_cols    = NULL,
  exclude_cols    = c(
    "geometry", "geom",
    "fire_id", "patch_id", "candidate_id", "sample_id",
    "year", "date", "doy",
    "label", "training_label", "target",
    "pool_source", "pool_bucket", "source",
    "VISUAL", "visual", "used", "used_in_training",
    "sample_weight"
  ),
  only_used            = FALSE,
  visual_filter        = NULL,
  min_n_per_group      = 10,
  top_n_features       = 12,
  print_pairwise_top_n = 5,
  run_pairwise         = TRUE,
  run_pca              = TRUE,
  run_permanova        = FALSE,
  run_rs_separability  = FALSE,
  run_probe_classifier = FALSE,
  run_correlation      = TRUE,
  corr_method          = "spearman",
  corr_threshold       = 0.9,
  balance_n            = NULL,
  n_repeats            = 1,
  output_dir           = NULL,
  prefix               = "training_pools",
  make_plots           = TRUE,
  seed                 = 123,
  verbose              = TRUE
) {
  t0 <- Sys.time()
  set.seed(seed)

  # --- coerce input (never mutate caller) ---
  df <- .pool_as_df(pools)
  if (!group_col %in% names(df)) stop("group_col '", group_col, "' not found in 'pools'.")

  # --- filters ---
  if (isTRUE(only_used)) {
    if (is.null(used_col) || !used_col %in% names(df)) {
      .pool_msg(verbose, "only_used=TRUE but used_col missing/absent; ignored.")
    } else {
      u <- df[[used_col]]
      keep <- if (is.logical(u)) u %in% TRUE else as.character(u) %in% c("TRUE", "1", "yes", "Y")
      df <- df[which(keep), , drop = FALSE]
    }
  }
  if (!is.null(visual_filter) && !is.null(visual_col) && visual_col %in% names(df)) {
    df <- df[df[[visual_col]] %in% visual_filter, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No rows left after filtering.")

  # --- feature detection ---
  if (is.null(feature_cols)) {
    is_num <- vapply(df, is.numeric, logical(1))
    cand <- names(df)[is_num]
    excl_all <- unique(c(exclude_cols, group_col, label_col, visual_col, used_col))
    cand <- cand[!.pool_in_excl(cand, excl_all)]
    excluded_by_rule <- setdiff(names(df)[is_num], cand)
  } else {
    cand <- feature_cols[feature_cols %in% names(df)]
    cand <- cand[vapply(df[cand], is.numeric, logical(1))]
    excluded_by_rule <- setdiff(feature_cols, cand)
  }
  if (length(cand) == 0) stop("No numeric feature columns available after exclusions.")

  # zero-variance removal (report)
  variances <- vapply(cand, function(f) stats::var(df[[f]], na.rm = TRUE), numeric(1))
  zv <- cand[is.na(variances) | variances == 0]
  features <- setdiff(cand, zv)
  zero_variance <- data.frame(feature = zv, stringsAsFactors = FALSE)
  if (length(features) == 0) stop("All candidate features are zero-variance.")

  # --- drop small groups (report) ---
  gtab <- table(df[[group_col]], useNA = "no")
  small_groups <- names(gtab)[gtab < min_n_per_group]
  if (length(small_groups)) {
    .pool_msg(verbose, "Dropping %d small group(s) from tests (n < %d): %s",
              length(small_groups), min_n_per_group, paste(small_groups, collapse = ", "))
    df_test <- df[!df[[group_col]] %in% small_groups, , drop = FALSE]
  } else {
    df_test <- df
  }
  df_test[[group_col]] <- as.character(df_test[[group_col]])
  df_test <- df_test[!is.na(df_test[[group_col]]), , drop = FALSE]

  # --- optional balancing (seeded) ---
  if (!is.null(balance_n)) {
    df_test <- .pool_subsample(df_test, group_col, balance_n, seed)
  }

  have_effectsize <- requireNamespace("effectsize", quietly = TRUE)

  # === A. composition & quality ===
  comp <- .pool_composition(df_test, group_col, label_col, visual_col, used_col, features)
  highcor <- .pool_highcor(df_test, features)

  # === B. univariate ===
  uni <- .pool_univariate(df_test, group_col, features, run_pairwise, have_effectsize)
  feature_ranking <- uni$kruskal[, c("feature", "statistic", "epsilon_squared",
                                     "eta_squared", "p", "p_adj")]
  feature_ranking$rank <- seq_len(nrow(feature_ranking))

  # === H. correlation / redundancy (importance from feature_ranking) ===
  correlation <- NULL
  if (isTRUE(run_correlation) && length(features) >= 2) {
    imp_vec <- stats::setNames(feature_ranking$epsilon_squared,
                               feature_ranking$feature)
    correlation <- .pool_correlation(df_test, features,
                                     family_fun = .pool_feature_family,
                                     importance = imp_vec,
                                     method = corr_method, thr = corr_threshold)
  }

  # === C. z-score fingerprint ===
  zfp <- .pool_zscore(df_test, group_col, features)

  # === D. PCA ===
  pca <- if (isTRUE(run_pca)) .pool_pca(df_test, group_col, features, seed) else NULL

  # === E. PERMANOVA ===
  permanova <- if (isTRUE(run_permanova))
    .pool_permanova(df_test, group_col, features, seed, verbose) else NULL

  # === F. RS separability ===
  rs <- if (isTRUE(run_rs_separability))
    .pool_rs_separability(df_test, group_col, features, verbose) else NULL

  # === G. probe ===
  fold_key <- NULL
  for (kc in c("fire_id", "fire_uid", "year")) {
    if (kc %in% names(df_test)) { fold_key <- as.character(df_test[[kc]]); break }
  }
  probe <- if (isTRUE(run_probe_classifier))
    .pool_probe(df_test, group_col, features, fold_key, seed, verbose) else NULL

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out <- list(
    summary = list(
      composition         = comp$composition,
      cross_tabs          = comp$cross_tabs,
      high_correlations   = highcor,
      dropped_small_groups = small_groups,
      n_rows_tested       = nrow(df_test),
      elapsed_sec         = elapsed
    ),
    settings = list(
      group_col = group_col, label_col = label_col, visual_col = visual_col,
      used_col = used_col, only_used = only_used, visual_filter = visual_filter,
      min_n_per_group = min_n_per_group, balance_n = balance_n,
      n_repeats = n_repeats, seed = seed, n_rows = nrow(df),
      effectsize_available = have_effectsize, prefix = prefix
    ),
    feature_info = list(
      features_used   = features,
      excluded_columns = excluded_by_rule,
      zero_variance_removed = zv,
      effect_size_source = if (have_effectsize) "effectsize package" else
        "manual (epsilon^2 = H/(n-1); eta^2 = (H-k+1)/(n-k))"
    ),
    missingness        = comp$missingness,
    zero_variance      = zero_variance,
    kruskal            = uni$kruskal,
    pairwise           = uni$pairwise,
    feature_ranking    = feature_ranking,
    correlation        = correlation,
    zscore_fingerprint = zfp,
    pca                = pca,
    permanova          = permanova,
    rs_separability    = rs,
    probe_classifier   = probe,
    plot_paths         = list(),
    output_dir         = output_dir
  )
  class(out) <- c("otsufire_pool_diagnostics", class(out))

  console <- .pool_console_summary(out, top_n_features, print_pairwise_top_n)
  attr(out, "console_summary") <- console
  if (isTRUE(verbose)) cat(console, "\n")

  # --- plots ---
  if (isTRUE(make_plots) && !is.null(output_dir)) {
    fig_dir <- file.path(output_dir, "figures")
    out$plot_paths <- tryCatch(
      .pool_make_plots(out, df_test, group_col, features, fig_dir, prefix, seed, verbose),
      error = function(e) { .pool_msg(verbose, "plotting failed: %s", conditionMessage(e)); list() })
  }

  # --- file outputs ---
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    wr <- function(obj, fn) {
      if (is.null(obj) || (is.data.frame(obj) && nrow(obj) == 0)) return(invisible())
      utils::write.csv(obj, file.path(output_dir, fn), row.names = FALSE)
    }
    writeLines(console, file.path(output_dir, "00_console_summary.txt"))
    wr(out$summary$composition,        paste0(prefix, "_01_pool_composition.csv"))
    wr(out$missingness,                paste0(prefix, "_02_feature_missingness.csv"))
    wr(out$zero_variance,              paste0(prefix, "_03_zero_variance_features.csv"))
    wr(out$kruskal,                    paste0(prefix, "_04_kruskal_results.csv"))
    wr(out$pairwise,                   paste0(prefix, "_05_pairwise_wilcoxon.csv"))
    wr(out$feature_ranking,            paste0(prefix, "_06_feature_ranking.csv"))
    wr(out$zscore_fingerprint,         paste0(prefix, "_07_zscore_fingerprint.csv"))
    if (!is.null(pca)) {
      wr(pca$scores,                   paste0(prefix, "_08_pca_scores.csv"))
      wr(pca$loadings,                 paste0(prefix, "_09_pca_loadings.csv"))
      wr(pca$explained,                paste0(prefix, "_10_pca_explained_variance.csv"))
      wr(pca$centroid_distances,       paste0(prefix, "_11_centroid_distances.csv"))
    }
    if (!is.null(permanova) && isTRUE(permanova$available) && !is.null(permanova$adonis2))
      wr(permanova$adonis2,            paste0(prefix, "_12_permanova_results.csv"))
    if (!is.null(rs) && isTRUE(rs$available))
      wr(rs$pairwise,                  paste0(prefix, "_13_rs_separability.csv"))
    if (!is.null(probe) && isTRUE(probe$available)) {
      pr <- data.frame(method = probe$method, cv = probe$cv, n = probe$n,
                       balanced_accuracy = probe$balanced_accuracy,
                       macro_f1 = probe$macro_f1, stringsAsFactors = FALSE)
      wr(pr,                           paste0(prefix, "_14_probe_classifier_results.csv"))
    }
    if (!is.null(correlation)) {
      wr(correlation$high_pairs,       paste0(prefix, "_15_correlation_high_pairs.csv"))
      wr(correlation$clusters,         paste0(prefix, "_16_correlation_clusters.csv"))
      reduced_df <- data.frame(
        feature    = correlation$reduced_set,
        family     = .pool_feature_family(correlation$reduced_set)$family,
        stringsAsFactors = FALSE)
      wr(reduced_df,                   paste0(prefix, "_17_correlation_reduced_set.csv"))
    }
  }

  out
}

# Print method for pool diagnostics (registered S3 method; like the package's
# other print methods it has no help page of its own).
#' @exportS3Method print otsufire_pool_diagnostics
print.otsufire_pool_diagnostics <- function(x, ...) {
  cs <- attr(x, "console_summary")
  if (is.null(cs)) cs <- .pool_console_summary(x, 12L, 5L)
  cat(cs, "\n")
  invisible(x)
}
