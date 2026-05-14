#' Compare a numeric variable across polygon categories with Kruskal-Wallis, Dunn letters, and plots
#'
#' @description
#' Given polygons (sf or path), this function builds a long data.frame where:
#' - X is a categorical column (e.g., Confusion, flag_ref, flag_internal)
#' - Y is a numeric variable computed from:
#'   - an existing attribute column (y_mode="attribute")
#'   - polygon area in hectares (y_mode="area_ha")
#'   - extraction from one or more rasters (y_mode="raster") (e.g., RBR fire/p1/p2)
#'   - percent overlap with an overlay layer (y_mode="overlay_percent") (e.g., % regenerated)
#'
#' It runs Kruskal-Wallis per panel (optional), then Dunn post-hoc (PMCMRplus) with letters
#' (multcompView), and produces a ggplot (boxplot + optional jitter) with letters.
#'
#' @param polys sf object OR path to a vector file (.gpkg/.shp/.geojson).
#' @param gpkg_layer Optional layer name if polys is a .gpkg. If NULL, auto-picks a layer.
#' @param x_var Character. Name of the categorical column to use on the X axis.
#' @param y_mode One of c("attribute","area_ha","raster","overlay_percent").
#' @param y_var Character. If y_mode="attribute", the name of a numeric column in polys.
#' @param y_name Character. Output column name for the numeric variable (default "Y").
#'
#' @param raster_list Named list of raster paths (or SpatRaster). Only used if y_mode="raster".
#'   Example: list(fire="...2012.tif", p1="...2013.tif", p2="...2014.tif")
#' @param raster_stat One of c("mean","median"). Statistic per polygon for raster extraction.
#' @param raster_template Optional SpatRaster used as template for alignment (res/crs).
#'   If NULL and raster_list provided, the first raster becomes the template.
#' @param raster_panel_name Character. Name of the panel column created for raster_list names (default "Period").
#'
#' @param overlay sf object OR path. Only used if y_mode="overlay_percent".
#' @param overlay_union Logical. If TRUE, union overlay geometries before intersection (avoids double counting).
#' @param overlay_percent_base One of c("polys","overlay"). If "polys": percent of each polygon covered by overlay.
#'
#' @param panel_vars Character vector of columns used to define separate panels/tests (optional).
#'   Examples: c("Source","Period") or c("RBR_year","Source"). If NULL, single global test.
#'
#' @param x_levels Optional explicit order for x_var factor levels.
#' @param alpha Numeric. Significance level for deciding whether to compute Dunn letters (default 0.05).
#' @param p_adjust Character. P adjustment method for Dunn p-values (default "bonferroni").
#' @param show_points Logical. Add jittered points over boxplots (default TRUE).
#'
#' @param facet_type One of c("wrap","grid"). If grid and panel_vars has 2 vars, uses facet_grid(row~col).
#' @param out_dir Optional directory to save outputs (plot + CSV summaries). If NULL, nothing is saved.
#' @param file_prefix Output prefix for saved files (default "compare_by_group").
#' @param width,height,dpi Plot saving parameters.
#'
#' @return A list with:
#' - data: long data.frame (all original attributes + new columns)
#' - kw: Kruskal-Wallis summary per panel
#' - letters: letter table per panel (if computed)
#' - plot: ggplot object
#' - paths: saved file paths (if out_dir provided)
#'
#' @import sf
#' @importFrom stats kruskal.test
#' @importFrom PMCMRplus kwAllPairsDunnTest
#' @importFrom multcompView multcompLetters
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_text labs theme_bw theme facet_wrap ggsave
#' @importFrom dplyr mutate filter select left_join group_by summarise distinct across all_of
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom data.table as.data.table setDT set setcolorder fwrite rbindlist
#' @export

validation_statistics <- function(
    polys,
    gpkg_layer = NULL,
    x_var,
    y_mode = c("attribute","area_ha","raster","overlay_percent"),
    y_var = NULL,
    y_name = "Y",

    raster_list = NULL,
    raster_stat = c("mean","median"),
    raster_template = NULL,
    raster_panel_name = "Period",

    overlay = NULL,
    overlay_union = TRUE,
    overlay_percent_base = c("polys","overlay"),

    panel_vars = NULL,

    x_levels = NULL,
    alpha = 0.05,
    p_adjust = "bonferroni",
    show_points = TRUE,

    facet_type = c("wrap","grid"),
    facet_scales = c("fixed","free_y"),

    out_dir = NULL,
    file_prefix = "compare_by_group",
    width = 12,
    height = 7,
    dpi = 300,

    # --- attribute transforms (only used when y_mode == "attribute") ---
    log_area = FALSE,
    log_area_var = "Area_ha",
    log_area_base = 10,
    log_eps = 1e-9,
    y_transforms = NULL
) {
  # ----------------------------
  # args
  # ----------------------------
  y_mode <- match.arg(y_mode)
  raster_stat <- match.arg(raster_stat)
  overlay_percent_base <- match.arg(overlay_percent_base)
  facet_type <- match.arg(facet_type)
  facet_scales <- match.arg(facet_scales)
  facet_scales_gg <- if (facet_scales == "free_y") "free_y" else "fixed"

  if (!is.character(x_var) || length(x_var) != 1L || !nzchar(trimws(x_var))) {
    stop("x_var must be a non-empty character scalar.")
  }
  x_var <- trimws(x_var)

  # ----------------------------
  # helpers
  # ----------------------------
  .safe_slug <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    if (!nzchar(x)) "out" else x
  }

  read_sf_any <- function(path, layer = NULL) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "gpkg") {
      lyr_names <- sf::st_layers(path)$name
      lyr <- layer

      if (!is.null(lyr) && nzchar(lyr) && !(lyr %in% lyr_names)) {
        stop("Requested layer '", lyr, "' not found. Available: ", paste(lyr_names, collapse = ", "))
      }
      if (is.null(lyr) || !nzchar(lyr)) {
        pref <- c("confusion_polys","polys","scored","ref","internal_flagged")
        hit <- pref[pref %in% lyr_names][1]
        if (is.na(hit) || !nzchar(hit)) hit <- lyr_names[1]
        lyr <- hit
      }
      sf::st_read(path, layer = lyr, quiet = TRUE)
    } else {
      sf::st_read(path, quiet = TRUE)
    }
  }

  make_valid_sf <- function(x) {
    x <- sf::st_as_sf(x)
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")
    x <- tryCatch(sf::st_make_valid(x), error = function(e) suppressWarnings(sf::st_buffer(x, 0)))
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    gt <- sf::st_geometry_type(x, by_geometry = TRUE)
    if (any(!gt %in% c("POLYGON","MULTIPOLYGON"))) {
      x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON"))
      x <- x[!sf::st_is_empty(x), , drop = FALSE]
    }
    if (nrow(x) > 0) {
      sf::st_geometry(x) <- sf::st_cast(sf::st_geometry(x), "MULTIPOLYGON", warn = FALSE)
      x <- sf::st_as_sf(x)
    }
    x
  }

  normalize_x_factor <- function(df, x_var, x_levels) {
    xx <- as.character(df[[x_var]])
    xx <- trimws(xx)
    keep <- !is.na(xx) & nzchar(xx)
    df <- df[keep, , drop = FALSE]
    xx <- xx[keep]

    if (!is.null(x_levels) && length(x_levels) > 0) {
      x_levels <- as.character(x_levels)
      x_levels <- trimws(x_levels)
      x_levels <- x_levels[nzchar(x_levels)]
      x_levels <- unique(x_levels)
      # keep any extra observed levels to avoid dropping data
      extra <- setdiff(unique(xx), x_levels)
      levs <- c(x_levels, extra)
      df[[x_var]] <- factor(xx, levels = levs)
    } else {
      df[[x_var]] <- factor(xx)
    }
    df
  }

  normalize_panel_vars <- function(df, panel_vars) {
    if (is.null(panel_vars) || !length(panel_vars)) return(NULL)

    panel_vars <- as.character(panel_vars)
    panel_vars <- trimws(panel_vars)
    panel_vars <- panel_vars[nzchar(panel_vars)]
    panel_vars <- unique(panel_vars)
    if (!length(panel_vars)) return(NULL)

    miss <- setdiff(panel_vars, names(df))
    if (length(miss) > 0) {
      warning("panel_vars not found in data (ignored): ", paste(miss, collapse = ", "))
      panel_vars <- intersect(panel_vars, names(df))
    }
    if (!length(panel_vars)) return(NULL)

    for (pv in panel_vars) {
      df[[pv]] <- as.factor(df[[pv]])
      df[[pv]] <- droplevels(df[[pv]])
    }
    panel_vars
  }

  get_panel_key <- function(df, panel_vars) {
    if (is.null(panel_vars) || !length(panel_vars)) {
      rep("ALL", nrow(df))
    } else {
      interaction(df[, panel_vars, drop = FALSE], drop = TRUE, sep = " | ")
    }
  }

  # letters deps (optional)
  has_pmcmr <- requireNamespace("PMCMRplus", quietly = TRUE)
  has_mcl   <- requireNamespace("multcompView", quietly = TRUE)
  do_letters <- isTRUE(has_pmcmr) && isTRUE(has_mcl)
  warned_letters <- FALSE

  dunn_letters <- function(df_sub, group_col, value_col, alpha, p_adjust) {
    lev <- levels(droplevels(df_sub[[group_col]]))
    if (length(lev) < 2) {
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    kw0 <- try(stats::kruskal.test(df_sub[[value_col]] ~ df_sub[[group_col]]), silent = TRUE)
    if (inherits(kw0, "try-error") || is.na(kw0$p.value) || kw0$p.value >= alpha) {
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    if (!do_letters) {
      if (!warned_letters) {
        warning("Letters require packages 'PMCMRplus' and 'multcompView'. Returning 'a' for all groups.")
        warned_letters <<- TRUE
      }
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    kw <- try(
      PMCMRplus::kwAllPairsDunnTest(
        as.formula(paste(value_col, "~", group_col)),
        data = df_sub,
        p.adjust.method = p_adjust
      ),
      silent = TRUE
    )
    if (inherits(kw, "try-error") || is.null(kw$p.value)) {
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    pvals <- kw$p.value
    rr <- rownames(pvals); cc <- colnames(pvals)
    if (is.null(rr) || is.null(cc)) {
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    grid_names <- apply(expand.grid(rr, cc), 1, function(x) paste(sort(x), collapse = "-"))
    pv <- as.vector(pvals)
    keep <- which(!is.na(pv))
    if (!length(keep)) {
      return(data.frame(level = lev, letter = "a", stringsAsFactors = FALSE))
    }

    pv2 <- pv[keep]
    names(pv2) <- grid_names[keep]

    sig <- multcompView::multcompLetters(pv2)
    out <- data.frame(level = names(sig$Letters), letter = sig$Letters, stringsAsFactors = FALSE)
    out <- merge(data.frame(level = lev, stringsAsFactors = FALSE), out, by = "level", all.x = TRUE, sort = FALSE)
    out$letter[is.na(out$letter)] <- "a"
    out
  }

  apply_transform <- function(x, tr, eps = 1e-9) {
    x <- suppressWarnings(as.numeric(x))
    if (tr == "none") return(x)
    if (tr %in% c("log10","ln")) {
      x[x <= 0] <- NA_real_
      if (tr == "log10") return(log10(pmax(x, eps)))
      return(log(pmax(x, eps)))
    }
    if (tr == "log1p") {
      x[x < 0] <- NA_real_
      return(log1p(pmax(x, 0)))
    }
    if (tr == "sqrt") {
      x[x < 0] <- NA_real_
      return(sqrt(x))
    }
    x
  }

  normalize_attribute_targets <- function(df_long, y_var, y_name, file_prefix,
                                          log_area, log_area_var, log_area_base, y_transforms) {
    # y_var as vector
    y_vars <- as.character(y_var)
    y_vars <- trimws(y_vars)
    y_vars <- y_vars[nzchar(y_vars)]
    y_vars <- unique(y_vars)
    if (!length(y_vars)) stop("y_var is empty after normalization.")

    miss_y <- setdiff(y_vars, names(df_long))
    if (length(miss_y) > 0) stop("These y_var columns are missing in data: ", paste(miss_y, collapse = ", "))

    # y_name: length 1, same length as y_vars, or named by y_var
    y_names <- y_name
    if (is.null(y_names) || !length(y_names)) y_names <- "Y"

    if (!is.null(names(y_names)) && any(nzchar(names(y_names)))) {
      y_names <- as.character(y_names)
      y_names <- y_names[y_vars]
    } else {
      y_names <- as.character(y_names)
      if (length(y_names) == 1L && length(y_vars) > 1L) {
        y_names <- rep(y_names, length(y_vars))
      } else if (length(y_names) != length(y_vars)) {
        stop("If y_var is a vector, y_name must be length 1, length(y_var), or a named vector keyed by y_var.")
      }
    }
    if (any(!nzchar(y_names))) stop("Some y_name entries are empty after normalization.")

    # file_prefix: scalar, aligned vector, or named by y_var
    fp <- file_prefix
    if (is.null(fp) || !length(fp) || !nzchar(as.character(fp)[1])) fp <- "plot"

    if (!is.null(names(fp)) && any(nzchar(names(fp)))) {
      fp <- as.character(fp)
      file_prefixes <- fp[y_vars]
      for (i in seq_along(y_vars)) {
        if (is.na(file_prefixes[i]) || !nzchar(file_prefixes[i])) {
          file_prefixes[i] <- paste0(as.character(fp[1]), "_", y_vars[i])
        }
      }
    } else {
      fp <- as.character(fp)
      if (length(fp) == 1L) {
        file_prefixes <- paste0(fp, "_", y_vars)
      } else if (length(fp) == length(y_vars)) {
        file_prefixes <- fp
      } else {
        stop("file_prefix must be length 1, length(y_var), or a named vector keyed by y_var.")
      }
    }

    # transform map
    allowed_tr <- c("none","log10","ln","log1p","sqrt")
    tr_map <- setNames(rep("none", length(y_vars)), y_vars)

    # convenience: log_area
    if (isTRUE(log_area) && (log_area_var %in% y_vars)) {
      tr_map[log_area_var] <- if (isTRUE(log_area_base == 10)) "log10" else "ln"
    }

    # explicit overrides
    if (!is.null(y_transforms) && length(y_transforms) > 0) {
      yt <- y_transforms
      if (!is.null(names(yt)) && any(nzchar(names(yt)))) {
        yt <- as.character(yt)
        hit <- intersect(names(yt), y_vars)
        tr_map[hit] <- yt[hit]
      } else {
        yt <- as.character(yt)
        if (length(yt) == 1L && length(y_vars) > 1L) yt <- rep(yt, length(y_vars))
        if (length(yt) != length(y_vars)) {
          stop("y_transforms must be NULL, length 1, length(y_var), or a named vector keyed by y_var.")
        }
        tr_map[] <- yt
      }
    }

    bad_tr <- setdiff(unique(tr_map), allowed_tr)
    if (length(bad_tr) > 0) {
      stop("Invalid y_transforms value(s): ", paste(bad_tr, collapse = ", "),
           "\nAllowed: ", paste(allowed_tr, collapse = ", "))
    }

    list(y_vars = y_vars, y_names = y_names, file_prefixes = file_prefixes, tr_map = tr_map)
  }

  compute_stats_and_plot <- function(df_in, y_col, x_var, panel_vars,
                                     facet_type, facet_scales_gg,
                                     alpha, p_adjust, show_points,
                                     width, height, dpi,
                                     out_dir, file_prefix, title_y) {

    # panel key
    pk <- get_panel_key(df_in, panel_vars)
    panels <- levels(factor(pk))

    # KW + letters per panel
    kw_list <- list()
    let_list <- list()

    for (pp in panels) {
      sub <- df_in[pk == pp, , drop = FALSE]
      # need >=2 groups
      ng <- as.data.frame(table(sub[[x_var]]))
      ng <- ng[ng$Freq > 0, , drop = FALSE]
      if (nrow(ng) < 2) next

      kw0 <- try(stats::kruskal.test(sub[[y_col]] ~ sub[[x_var]]), silent = TRUE)
      p_kw <- if (inherits(kw0, "try-error")) NA_real_ else as.numeric(kw0$p.value)

      kw_row <- data.frame(panel = pp, p_kw = p_kw, n = nrow(sub), stringsAsFactors = FALSE)
      if (!is.null(panel_vars) && length(panel_vars)) {
        # store panel vars values
        key_vals <- strsplit(pp, " \\| ", fixed = FALSE)[[1]]
        # when interaction() is used, pp is already a label; safest: rebuild from first row of sub
        for (pv in panel_vars) kw_row[[pv]] <- as.character(sub[[pv]][1])
      }
      kw_list[[length(kw_list) + 1]] <- kw_row

      let <- dunn_letters(sub, x_var, y_col, alpha, p_adjust)
      let[[x_var]] <- let$level
      let$level <- NULL
      let$panel <- pp
      if (!is.null(panel_vars) && length(panel_vars)) {
        for (pv in panel_vars) let[[pv]] <- as.character(sub[[pv]][1])
      }
      let_list[[length(let_list) + 1]] <- let
    }

    kw_df <- if (length(kw_list)) do.call(rbind, kw_list) else NULL
    letters_df <- if (length(let_list)) do.call(rbind, let_list) else NULL

    # letters y-position per panel + group
    if (!is.null(letters_df) && nrow(letters_df) > 0) {
      rng <- range(df_in[[y_col]], na.rm = TRUE)
      off <- 0.04 * diff(rng)
      if (!is.finite(off) || off == 0) off <- 1

      if (is.null(panel_vars) || !length(panel_vars)) {
        ymax <- aggregate(df_in[[y_col]], by = list(x = df_in[[x_var]]), FUN = max, na.rm = TRUE)
        names(ymax)[2] <- "ymax"
        letters_df <- merge(letters_df, ymax, by.x = x_var, by.y = "x", all.x = TRUE, sort = FALSE)
        letters_df$ypos <- letters_df$ymax + off
      } else {
        # compute max per panel vars + x
        key_df <- df_in[, c(panel_vars, x_var, y_col), drop = FALSE]
        key_df$y <- key_df[[y_col]]
        ymax <- aggregate(key_df$y, by = key_df[, c(panel_vars, x_var), drop = FALSE], FUN = max, na.rm = TRUE)
        names(ymax)[ncol(ymax)] <- "ymax"
        letters_df <- merge(letters_df, ymax, by = c(panel_vars, x_var), all.x = TRUE, sort = FALSE)
        letters_df$ypos <- letters_df$ymax + off
      }
    }

    # common ylim (only when fixed)
    coord_ylim <- NULL
    if (facet_scales_gg == "fixed") {
      yv <- df_in[[y_col]]
      yv <- yv[is.finite(yv)]
      if (length(yv)) {
        y_min <- min(yv, na.rm = TRUE)
        y_max <- max(yv, na.rm = TRUE)
        pad <- if ((y_max - y_min) > 0) 0.10 * (y_max - y_min) else 1
        # ensure letters not clipped (if present)
        if (!is.null(letters_df) && nrow(letters_df) > 0 && any(is.finite(letters_df$ypos))) {
          y_max2 <- max(c(y_max, letters_df$ypos), na.rm = TRUE)
          coord_ylim <- c(y_min, y_max2 + pad)
        } else {
          coord_ylim <- c(y_min, y_max + pad)
        }
      }
    }

    # plot
    p <- ggplot2::ggplot(df_in, ggplot2::aes(x = .data[[x_var]], y = .data[[y_col]], fill = .data[[x_var]])) +
      ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.25) +
      ggplot2::labs(
        title = paste0(title_y, " by ", x_var, " (Kruskal-Wallis + Dunn letters)"),
        x = x_var,
        y = title_y
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(legend.position = "none")

    if (!is.null(coord_ylim)) {
      p <- p + ggplot2::coord_cartesian(ylim = coord_ylim)
    }

    if (isTRUE(show_points)) {
      p <- p + ggplot2::geom_jitter(width = 0.15, alpha = 0.15, size = 0.6)
    }

    # facets
    if (!is.null(panel_vars) && length(panel_vars) > 0) {
      if (facet_type == "grid" && length(panel_vars) == 2) {
        p <- p + ggplot2::facet_grid(
          stats::as.formula(paste(panel_vars[1], "~", panel_vars[2])),
          scales = facet_scales_gg
        )
      } else if (length(panel_vars) == 1) {
        p <- p + ggplot2::facet_wrap(
          stats::as.formula(paste("~", panel_vars[1])),
          scales = facet_scales_gg
        )
      } else {
        df_in$.panel <- apply(df_in[, panel_vars, drop = FALSE], 1,
                              function(z) paste(paste(panel_vars, z, sep = "="), collapse = " | "))
        p <- p %+% df_in + ggplot2::facet_wrap(~ .panel, scales = facet_scales_gg)

        if (!is.null(letters_df) && nrow(letters_df) > 0) {
          letters_df$.panel <- apply(letters_df[, panel_vars, drop = FALSE], 1,
                                     function(z) paste(paste(panel_vars, z, sep = "="), collapse = " | "))
        }
      }
    }

    # letters
    if (!is.null(letters_df) && nrow(letters_df) > 0) {
      p <- p + ggplot2::geom_text(
        data = letters_df,
        ggplot2::aes(x = .data[[x_var]], y = ypos, label = letter),
        inherit.aes = FALSE,
        size = 4
      )
    }

    # save (optional)
    paths <- list(plot_png = NULL, kw_csv = NULL, letters_csv = NULL)
    if (!is.null(out_dir) && nzchar(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      paths$plot_png <- file.path(out_dir, paste0(.safe_slug(file_prefix), ".png"))
      ggplot2::ggsave(paths$plot_png, p, width = width, height = height, dpi = dpi)

      if (!is.null(kw_df) && nrow(kw_df) > 0) {
        paths$kw_csv <- file.path(out_dir, paste0(.safe_slug(file_prefix), "_kw.csv"))
        utils::write.table(kw_df, paths$kw_csv, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }

      if (!is.null(letters_df) && nrow(letters_df) > 0) {
        paths$letters_csv <- file.path(out_dir, paste0(.safe_slug(file_prefix), "_letters.csv"))
        utils::write.table(letters_df, paths$letters_csv, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }

    list(data = df_in, kw = kw_df, letters = letters_df, plot = p, paths = paths)
  }

  # ----------------------------
  # read polygons
  # ----------------------------
  sf_polys <- if (inherits(polys, "sf")) polys else read_sf_any(polys, layer = gpkg_layer)
  sf_polys <- make_valid_sf(sf_polys)

  if (!(x_var %in% names(sf_polys))) stop("x_var not found in polygons: ", x_var)

  sf_polys$.__pid <- seq_len(nrow(sf_polys))
  base_df <- sf::st_drop_geometry(sf_polys)

  # ----------------------------
  # build df_long (shared)
  # ----------------------------
  df_long <- NULL

  if (y_mode == "attribute") {
    if (is.null(y_var) || !length(y_var)) stop("y_var must be provided for y_mode='attribute'.")
    # NOTE: y_var can be vector, handled later; here we just keep all attributes

    df_long <- base_df

  } else if (y_mode == "area_ha") {
    area_ha <- as.numeric(sf::st_area(sf_polys)) / 10000
    df_long <- base_df
    df_long[[y_name]] <- area_ha

  } else if (y_mode == "raster") {
    if (!requireNamespace("terra", quietly = TRUE)) stop("y_mode='raster' requires package 'terra'.")
    if (is.null(raster_list) || length(raster_list) == 0) stop("raster_list must be a named list for y_mode='raster'.")
    if (is.null(names(raster_list)) || any(!nzchar(names(raster_list)))) stop("raster_list must be a *named* list.")

    # template
    if (is.null(raster_template)) {
      first <- raster_list[[1]]
      raster_template <- if (inherits(first, "SpatRaster")) first else terra::rast(first)
      raster_template <- raster_template[[1]]
    } else {
      raster_template <- if (inherits(raster_template, "SpatRaster")) raster_template else terra::rast(raster_template)
      raster_template <- raster_template[[1]]
    }

    # align CRS
    rcrs <- terra::crs(raster_template, proj = TRUE)
    if (is.na(sf::st_crs(sf_polys))) stop("Polygons have NA CRS.")
    if (!is.na(rcrs) && sf::st_crs(sf_polys)$wkt != rcrs) {
      sf_polys <- sf::st_transform(sf_polys, rcrs)
      base_df <- sf::st_drop_geometry(sf_polys)
    }

    vpolys <- terra::vect(sf_polys)
    fun_extract <- if (raster_stat == "mean") mean else median

    res_list <- list()
    for (nm in names(raster_list)) {
      rr <- raster_list[[nm]]
      rras <- if (inherits(rr, "SpatRaster")) rr else terra::rast(rr)
      if (terra::nlyr(rras) > 1) rras <- rras[[1]]

      if (!terra::same.crs(rras, raster_template)) {
        rras <- terra::project(rras, raster_template, method = "bilinear")
      }
      if (!isTRUE(terra::compareGeom(rras, raster_template, stopOnError = FALSE))) {
        rras <- terra::resample(rras, raster_template, method = "bilinear")
      }

      vals <- terra::extract(rras, vpolys, fun = fun_extract, na.rm = TRUE)
      vals <- as.data.frame(vals)
      names(vals)[1] <- ".__pid"
      names(vals)[2] <- y_name
      vals[[raster_panel_name]] <- nm
      res_list[[nm]] <- vals
    }

    y_df <- do.call(rbind, res_list)
    df_long <- base_df
    df_long <- merge(df_long, y_df, by = ".__pid", all.x = TRUE, sort = FALSE)

  } else if (y_mode == "overlay_percent") {
    if (is.null(overlay)) stop("overlay must be provided for y_mode='overlay_percent'.")
    ov_sf <- if (inherits(overlay, "sf")) overlay else read_sf_any(overlay, layer = NULL)
    ov_sf <- make_valid_sf(ov_sf)

    if (is.na(sf::st_crs(ov_sf))) stop("Overlay has NA CRS.")
    if (sf::st_crs(ov_sf) != sf::st_crs(sf_polys)) {
      ov_sf <- sf::st_transform(ov_sf, sf::st_crs(sf_polys))
      ov_sf <- make_valid_sf(ov_sf)
    }

    poly_area <- as.numeric(sf::st_area(sf_polys))
    if (isTRUE(overlay_union)) {
      ov_geom <- sf::st_union(sf::st_geometry(ov_sf))
      ov_geom <- sf::st_sfc(ov_geom, crs = sf::st_crs(sf_polys))
      overlay_base_area <- as.numeric(sf::st_area(ov_geom))
    } else {
      ov_geom <- ov_sf
      overlay_base_area <- sum(as.numeric(sf::st_area(ov_sf)))
    }

    tmp <- sf_polys[, c(".__pid")]
    inter <- suppressWarnings(sf::st_intersection(tmp, ov_geom))
    inter <- make_valid_sf(inter)

    inter_area <- rep(0, nrow(sf_polys))
    if (nrow(inter) > 0) {
      a <- as.numeric(sf::st_area(inter))
      pid <- inter$.__pid
      a_by <- tapply(a, pid, sum)
      inter_area[as.integer(names(a_by))] <- as.numeric(a_by)
    }

    if (overlay_percent_base == "polys") {
      pct <- ifelse(poly_area > 0, 100 * inter_area / poly_area, NA_real_)
    } else {
      pct <- ifelse(overlay_base_area > 0, 100 * inter_area / overlay_base_area, NA_real_)
    }

    df_long <- base_df
    df_long[[y_name]] <- as.numeric(pct)
  }

  # ----------------------------
  # normalize x + panel vars
  # ----------------------------
  df_long <- normalize_x_factor(df_long, x_var, x_levels)
  panel_vars <- normalize_panel_vars(df_long, panel_vars)

  # ----------------------------
  # run (single or multi output)
  # ----------------------------
  results <- list()

  if (y_mode == "attribute") {
    targets <- normalize_attribute_targets(
      df_long = df_long,
      y_var = y_var,
      y_name = y_name,
      file_prefix = file_prefix,
      log_area = log_area,
      log_area_var = log_area_var,
      log_area_base = log_area_base,
      y_transforms = y_transforms
    )

    y_vars <- targets$y_vars
    y_names <- targets$y_names
    file_prefixes <- targets$file_prefixes
    tr_map <- targets$tr_map

    for (i in seq_along(y_vars)) {
      yv <- y_vars[i]
      yn <- y_names[i]
      tr <- tr_map[[yv]]

      dfp <- df_long
      # build numeric Y column used for plotting/stats
      dfp[[yn]] <- apply_transform(dfp[[yv]], tr = tr, eps = log_eps)
      dfp <- dfp[!is.na(dfp[[yn]]), , drop = FALSE]
      if (nrow(dfp) == 0) {
        warning("Skipping y_var='", yv, "' because all values are NA after transform/coercion.")
        next
      }

      # title label
      title_y <- if (tr == "none") yn else paste0(yn, " (", tr, ")")

      # compute + plot
      res_i <- compute_stats_and_plot(
        df_in = dfp,
        y_col = yn,
        x_var = x_var,
        panel_vars = panel_vars,
        facet_type = facet_type,
        facet_scales_gg = facet_scales_gg,
        alpha = alpha,
        p_adjust = p_adjust,
        show_points = show_points,
        width = width,
        height = height,
        dpi = dpi,
        out_dir = out_dir,
        file_prefix = file_prefixes[i],
        title_y = title_y
      )

      results[[yv]] <- c(res_i, list(y_var = yv, y_name = yn, y_transform = tr))
    }

  } else {
    # single output modes
    if (!is.numeric(df_long[[y_name]]) && !is.integer(df_long[[y_name]])) {
      df_long[[y_name]] <- suppressWarnings(as.numeric(df_long[[y_name]]))
    }
    df_long <- df_long[!is.na(df_long[[y_name]]), , drop = FALSE]
    if (nrow(df_long) == 0) stop("No valid rows after filtering NA y values.")

    title_y <- y_name
    res_i <- compute_stats_and_plot(
      df_in = df_long,
      y_col = y_name,
      x_var = x_var,
      panel_vars = panel_vars,
      facet_type = facet_type,
      facet_scales_gg = facet_scales_gg,
      alpha = alpha,
      p_adjust = p_adjust,
      show_points = show_points,
      width = width,
      height = height,
      dpi = dpi,
      out_dir = out_dir,
      file_prefix = file_prefix,
      title_y = title_y
    )
    results[[y_name]] <- c(res_i, list(y_var = NULL, y_name = y_name, y_transform = "none"))
  }

  invisible(list(
    y_mode = y_mode,
    x_var = x_var,
    panel_vars = panel_vars,
    results = results
  ))
}
