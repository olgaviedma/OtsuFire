# Block 6: deterministic-vs-supervised consistency engine, migrated from
# 2_SCRIPTS/00_USAGE/03_ANALYSIS/CHECK_DETERMINISTIC_SUPERVISED_CONSISTENCY.R
# (~448 lines). The migrated form drops the script-level path-discovery
# preamble (find_project_paths_file, source(PROJECT_PATHS.R), top-level
# data_base/result_name globals, run_consistency_batch loop) and exposes
# only the per-run checker. Inputs are passed in explicitly so the engine
# is independent of working directory and external scripts_root.
#
# All four legacy helper functions plus the per-run checker live in this
# file. None of them is exported from the package; they are reached via
# the internal dispatcher and check_supervised_consistency() wrapper.

#' @keywords internal
#' @noRd
.of_consistency_safe_read_layer <- function(dsn, layer) {
  if (is.null(dsn) || is.na(dsn) || !nzchar(dsn) || !file.exists(dsn))
    return(NULL)
  tryCatch(
    sf::read_sf(dsn, layer = layer, quiet = TRUE),
    error = function(e) NULL
  )
}

#' @keywords internal
#' @noRd
.of_consistency_ensure_join_id <- function(x) {
  if (is.null(x)) return(x)
  # Bug 7 (0.3.0): also record which underlying column the join_id
  # was derived from so callers can assert that compared layers
  # resolve to the same namespace before computing row-difference
  # checks. The attribute is read by `.of_run_consistency_check()`.
  if ("source_poly_id" %in% names(x)) {
    x$join_id <- as.character(x$source_poly_id)
    attr(x, "join_id_source") <- "source_poly_id"
  } else if ("poly_id" %in% names(x)) {
    x$join_id <- as.character(x$poly_id)
    attr(x, "join_id_source") <- "poly_id"
  } else if ("fire_uid" %in% names(x)) {
    x$join_id <- as.character(x$fire_uid)
    attr(x, "join_id_source") <- "fire_uid"
  } else {
    stop("No suitable join id (source_poly_id / poly_id / fire_uid) ",
         "found in object.", call. = FALSE)
  }
  x
}

#' @keywords internal
#' @noRd
.of_consistency_count_duplicates <- function(x) {
  if (is.null(x) || !"join_id" %in% names(x)) return(NA_integer_)
  as.integer(sum(duplicated(x$join_id)))
}

#' @keywords internal
#' @noRd
.of_consistency_compute_no_hs_support <- function(x) {
  no_hs_support <- rep(FALSE, nrow(x))
  if ("hs_no_support_when_available" %in% names(x)) {
    no_hs_support <- x$hs_no_support_when_available == 1
    no_hs_support[is.na(no_hs_support)] <- FALSE
  } else if (all(c("hotspot_available", "hs_used_n") %in% names(x))) {
    no_hs_support <- x$hotspot_available == 1 & x$hs_used_n == 0
    no_hs_support[is.na(no_hs_support)] <- FALSE
  }
  no_hs_support
}

#' @keywords internal
#' @noRd
.of_consistency_append_issue <- function(issue_list, sf_obj, issue_type,
                                         issue_detail) {
  if (is.null(sf_obj) || !inherits(sf_obj, "sf") || !nrow(sf_obj))
    return(issue_list)
  sf_obj$issue_type <- issue_type
  sf_obj$issue_detail <- issue_detail
  issue_list[[length(issue_list) + 1L]] <- sf_obj
  issue_list
}

#' @keywords internal
#' @noRd
.of_consistency_calc_pburned_stats <- function(df, group_cols = NULL) {
  if (is.null(group_cols) || !length(group_cols)) {
    return(dplyr::summarise(df,
      n               = dplyr::n(),
      p_burned_min    = min(.data$p_burned, na.rm = TRUE),
      p_burned_q01    = as.numeric(stats::quantile(.data$p_burned, 0.01,
                                                    na.rm = TRUE)),
      p_burned_q05    = as.numeric(stats::quantile(.data$p_burned, 0.05,
                                                    na.rm = TRUE)),
      p_burned_q25    = as.numeric(stats::quantile(.data$p_burned, 0.25,
                                                    na.rm = TRUE)),
      p_burned_median = stats::median(.data$p_burned, na.rm = TRUE),
      p_burned_mean   = mean(.data$p_burned, na.rm = TRUE),
      p_burned_q75    = as.numeric(stats::quantile(.data$p_burned, 0.75,
                                                    na.rm = TRUE)),
      p_burned_q95    = as.numeric(stats::quantile(.data$p_burned, 0.95,
                                                    na.rm = TRUE)),
      p_burned_q99    = as.numeric(stats::quantile(.data$p_burned, 0.99,
                                                    na.rm = TRUE)),
      p_burned_max    = max(.data$p_burned, na.rm = TRUE),
      p_burned_sd     = stats::sd(.data$p_burned, na.rm = TRUE),
      n_lt_0_10       = sum(.data$p_burned < 0.10, na.rm = TRUE),
      n_lt_0_25       = sum(.data$p_burned < 0.25, na.rm = TRUE),
      n_lt_0_50       = sum(.data$p_burned < 0.50, na.rm = TRUE),
      n_gt_0_75       = sum(.data$p_burned > 0.75, na.rm = TRUE),
      n_gt_0_90       = sum(.data$p_burned > 0.90, na.rm = TRUE),
      n_gt_0_99       = sum(.data$p_burned > 0.99, na.rm = TRUE)
    ))
  }
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n               = dplyr::n(),
      p_burned_min    = min(.data$p_burned, na.rm = TRUE),
      p_burned_q05    = as.numeric(stats::quantile(.data$p_burned, 0.05,
                                                    na.rm = TRUE)),
      p_burned_q25    = as.numeric(stats::quantile(.data$p_burned, 0.25,
                                                    na.rm = TRUE)),
      p_burned_median = stats::median(.data$p_burned, na.rm = TRUE),
      p_burned_mean   = mean(.data$p_burned, na.rm = TRUE),
      p_burned_q75    = as.numeric(stats::quantile(.data$p_burned, 0.75,
                                                    na.rm = TRUE)),
      p_burned_q95    = as.numeric(stats::quantile(.data$p_burned, 0.95,
                                                    na.rm = TRUE)),
      p_burned_max    = max(.data$p_burned, na.rm = TRUE),
      p_burned_sd     = stats::sd(.data$p_burned, na.rm = TRUE),
      n_lt_0_50       = sum(.data$p_burned < 0.50, na.rm = TRUE),
      n_gt_0_90       = sum(.data$p_burned > 0.90, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Internal: per-run Otsu-guided segmentation-vs-probabilistic refinement consistency check.
#'
#' All inputs are explicit; the function does not consult getwd(),
#' PROJECT_PATHS.R, or any 00_USAGE/ script. Mirrors the per-run logic
#' of `run_consistency_check()` in the legacy
#' CHECK_DETERMINISTIC_SUPERVISED_CONSISTENCY.R.
#'
#' @param target_year Integer year.
#' @param scenario Scenario name (e.g. "balanced").
#' @param internal_decisions_gpkg Path to the deterministic
#'   `internal_decisions.gpkg` (must contain layer `internal_decisions`).
#' @param final_map_gpkg Path to the probabilistic refinement final-map GPKG
#'   (must contain layers `final_map_full` and `final_map`).
#' @param out_dir Directory where consistency artefacts are written;
#'   created if missing.
#' @param prefix_base Suffix used to construct the per-run prefix
#'   (default `"patch_certified"` — the legacy convention).
#' @return One-row `data.frame` summarising the consistency check;
#'   side effect: writes `*_consistency_issues.gpkg`,
#'   `*_consistency_summary.csv`, `*_consistency_summary.txt`,
#'   `*_pburned_overall_summary.csv`, optionally
#'   `*_pburned_by_class_input.csv` and `*_pburned_by_source_set.csv`.
#'
#' @keywords internal
#' @noRd
.of_run_consistency_check <- function(target_year, scenario,
                                       internal_decisions_gpkg,
                                       final_map_gpkg,
                                       out_dir,
                                       prefix_base = "patch_certified") {
  stopifnot(file.exists(internal_decisions_gpkg))
  stopifnot(file.exists(final_map_gpkg))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- sprintf("%d_%s_%s", target_year, scenario, prefix_base)

  issues_gpkg          <- file.path(out_dir,
                                    sprintf("%d_%s_consistency_issues.gpkg",
                                            target_year, scenario))
  summary_csv          <- file.path(out_dir,
                                    sprintf("%d_%s_consistency_summary.csv",
                                            target_year, scenario))
  summary_txt          <- file.path(out_dir,
                                    sprintf("%d_%s_consistency_summary.txt",
                                            target_year, scenario))
  pburned_overall_csv  <- file.path(out_dir,
                                    sprintf("%d_%s_pburned_overall_summary.csv",
                                            target_year, scenario))
  pburned_by_class_csv <- file.path(out_dir,
                                    sprintf("%d_%s_pburned_by_class_input.csv",
                                            target_year, scenario))
  pburned_by_source_csv <- file.path(out_dir,
                                    sprintf("%d_%s_pburned_by_source_set.csv",
                                            target_year, scenario))

  D <- .of_consistency_ensure_join_id(.of_consistency_safe_read_layer(
    internal_decisions_gpkg, "internal_decisions"))
  F_full <- .of_consistency_ensure_join_id(.of_consistency_safe_read_layer(
    final_map_gpkg, "final_map_full"))
  F_public <- .of_consistency_ensure_join_id(.of_consistency_safe_read_layer(
    final_map_gpkg, "final_map"))

  if (is.null(D)) stop("Could not read internal_decisions.", call. = FALSE)
  if (is.null(F_full)) stop("Could not read final_map_full.", call. = FALSE)
  if (is.null(F_public)) stop("Could not read final_map.", call. = FALSE)

  # Bug 7 (0.3.0): assert all three layers resolve their join_id from
  # the same underlying column (source_poly_id / poly_id / fire_uid).
  # If they don't, setdiff()-based comparisons produce meaningless
  # row-correspondence reports.
  .join_sources <- c(
    internal_decisions = attr(D, "join_id_source"),
    final_map_full     = attr(F_full, "join_id_source"),
    final_map          = attr(F_public, "join_id_source")
  )
  if (length(unique(.join_sources)) != 1L) {
    stop(sprintf(paste0("Consistency check aborted: layers resolved their join ",
                        "key from different namespaces (%s). Provide a shared ",
                        "join column across all compared layers."),
                 paste(sprintf("%s=%s", names(.join_sources), .join_sources),
                       collapse = "; ")),
         call. = FALSE)
  }

  det_ids    <- unique(D$join_id)
  full_ids   <- unique(F_full$join_id)
  public_ids <- unique(F_public$join_id)

  missing_in_full_ids <- setdiff(det_ids, full_ids)
  extra_in_full_ids   <- setdiff(full_ids, det_ids)

  missing_in_full <- D      |> dplyr::filter(.data$join_id %in% missing_in_full_ids)
  extra_in_full   <- F_full |> dplyr::filter(.data$join_id %in% extra_in_full_ids)

  invalid_pb <- F_full |> dplyr::filter(!is.finite(.data$p_burned))

  pb_ok <- F_full |>
    sf::st_drop_geometry() |>
    dplyr::filter(is.finite(.data$p_burned))

  pburned_overall <- .of_consistency_calc_pburned_stats(pb_ok)
  pburned_by_class <- if ("class_input" %in% names(pb_ok)) {
    .of_consistency_calc_pburned_stats(pb_ok, "class_input")
  } else if ("class_final" %in% names(pb_ok)) {
    .of_consistency_calc_pburned_stats(pb_ok, "class_final")
  } else {
    data.frame()
  }
  pburned_by_source <- if ("source_set" %in% names(pb_ok)) {
    .of_consistency_calc_pburned_stats(pb_ok, "source_set")
  } else {
    data.frame()
  }

  degenerate_high <- nrow(pburned_overall) == 1 &&
    isTRUE(pburned_overall$p_burned_min > 0.95)
  degenerate_low  <- nrow(pburned_overall) == 1 &&
    isTRUE(pburned_overall$p_burned_max < 0.05)
  degenerate_flat <- nrow(pburned_overall) == 1 &&
    isTRUE(pburned_overall$p_burned_sd < 0.01)

  class_mismatch <- D |>
    sf::st_drop_geometry() |>
    dplyr::select("join_id", class_final_det = "class_final") |>
    dplyr::left_join(
      F_full |>
        sf::st_drop_geometry() |>
        dplyr::select("join_id", class_final_sup = "class_final"),
      by = "join_id"
    ) |>
    dplyr::filter(!is.na(.data$class_final_sup),
                  .data$class_final_det != .data$class_final_sup)

  class_mismatch_sf <- F_full |>
    dplyr::filter(.data$join_id %in% class_mismatch$join_id) |>
    dplyr::left_join(class_mismatch, by = "join_id")

  filter_signature_mismatch <- D |>
    sf::st_drop_geometry() |>
    dplyr::select("join_id",
                  filter_1_det = "filter_1",
                  filter_2_det = "filter_2",
                  filter_3_det = "filter_3") |>
    dplyr::left_join(
      F_full |>
        sf::st_drop_geometry() |>
        dplyr::select("join_id",
                      filter_1_sup = "filter_1",
                      filter_2_sup = "filter_2",
                      filter_3_sup = "filter_3"),
      by = "join_id"
    ) |>
    dplyr::filter(
      !is.na(.data$filter_1_sup),
      .data$filter_1_det != .data$filter_1_sup |
        .data$filter_2_det != .data$filter_2_sup |
        .data$filter_3_det != .data$filter_3_sup
    )

  filter_signature_mismatch_sf <- F_full |>
    dplyr::filter(.data$join_id %in% filter_signature_mismatch$join_id) |>
    dplyr::left_join(filter_signature_mismatch, by = "join_id")

  no_hs_support <- .of_consistency_compute_no_hs_support(F_full)
  preyear_overlap <- if ("preyear_overlap_frac" %in% names(F_full))
    F_full$preyear_overlap_frac else rep(NA_real_, nrow(F_full))
  expected_drop_public <- is.finite(preyear_overlap) &
    preyear_overlap >= 0.7 & no_hs_support
  expected_public_ids  <- F_full$join_id[!expected_drop_public]
  expected_drop_ids    <- F_full$join_id[expected_drop_public]

  unexpected_missing_public_ids <- setdiff(expected_public_ids, public_ids)
  unexpected_kept_public_ids    <- intersect(expected_drop_ids, public_ids)

  unexpected_missing_public <- F_full   |>
    dplyr::filter(.data$join_id %in% unexpected_missing_public_ids)
  unexpected_kept_public    <- F_public |>
    dplyr::filter(.data$join_id %in% unexpected_kept_public_ids)

  empty_geom_det    <- sum(sf::st_is_empty(D))
  empty_geom_full   <- sum(sf::st_is_empty(F_full))
  empty_geom_public <- sum(sf::st_is_empty(F_public))

  summary_row <- data.frame(
    year = target_year,
    scenario = scenario,
    internal_decisions_rows = nrow(D),
    final_map_full_rows = nrow(F_full),
    final_map_rows = nrow(F_public),
    duplicate_join_id_det    = .of_consistency_count_duplicates(D),
    duplicate_join_id_full   = .of_consistency_count_duplicates(F_full),
    duplicate_join_id_public = .of_consistency_count_duplicates(F_public),
    missing_in_final_map_full = length(missing_in_full_ids),
    extra_in_final_map_full   = length(extra_in_full_ids),
    invalid_p_burned_full     = nrow(invalid_pb),
    class_final_mismatch      = nrow(class_mismatch),
    filter_signature_mismatch = nrow(filter_signature_mismatch),
    expected_public_drop      = sum(expected_drop_public, na.rm = TRUE),
    unexpected_missing_in_public = length(unexpected_missing_public_ids),
    unexpected_kept_in_public    = length(unexpected_kept_public_ids),
    p_burned_min    = if (nrow(pburned_overall)) pburned_overall$p_burned_min    else NA_real_,
    p_burned_q05    = if (nrow(pburned_overall)) pburned_overall$p_burned_q05    else NA_real_,
    p_burned_median = if (nrow(pburned_overall)) pburned_overall$p_burned_median else NA_real_,
    p_burned_q95    = if (nrow(pburned_overall)) pburned_overall$p_burned_q95    else NA_real_,
    p_burned_max    = if (nrow(pburned_overall)) pburned_overall$p_burned_max    else NA_real_,
    p_burned_sd     = if (nrow(pburned_overall)) pburned_overall$p_burned_sd     else NA_real_,
    degenerate_high_scores = degenerate_high,
    degenerate_low_scores  = degenerate_low,
    degenerate_flat_scores = degenerate_flat,
    empty_geom_det    = empty_geom_det,
    empty_geom_full   = empty_geom_full,
    empty_geom_public = empty_geom_public,
    summary_csv          = summary_csv,
    summary_txt          = summary_txt,
    issues_gpkg          = issues_gpkg,
    pburned_overall_csv  = pburned_overall_csv,
    pburned_by_class_csv = pburned_by_class_csv,
    pburned_by_source_csv = pburned_by_source_csv,
    stringsAsFactors = FALSE
  )

  issue_layers <- list()
  issue_layers <- .of_consistency_append_issue(issue_layers, missing_in_full,
    "missing_in_final_map_full",
    "Present in deterministic, absent in final_map_full")
  issue_layers <- .of_consistency_append_issue(issue_layers, extra_in_full,
    "extra_in_final_map_full",
    "Present in final_map_full, absent in deterministic")
  issue_layers <- .of_consistency_append_issue(issue_layers, invalid_pb,
    "invalid_p_burned",
    "Missing or non-finite p_burned in final_map_full")
  issue_layers <- .of_consistency_append_issue(issue_layers, class_mismatch_sf,
    "class_final_mismatch",
    "class_final differs between deterministic and final_map_full")
  issue_layers <- .of_consistency_append_issue(issue_layers,
    filter_signature_mismatch_sf,
    "filter_signature_mismatch",
    "filter_1/filter_2/filter_3 differ between deterministic and final_map_full")
  issue_layers <- .of_consistency_append_issue(issue_layers,
    unexpected_missing_public,
    "unexpected_missing_in_public",
    "Expected in final_map but missing")
  issue_layers <- .of_consistency_append_issue(issue_layers,
    unexpected_kept_public,
    "unexpected_kept_in_public",
    "Should have been removed from final_map by temporal overlap rule")

  if (file.exists(issues_gpkg)) file.remove(issues_gpkg)
  if (length(issue_layers) > 0) {
    for (i in seq_along(issue_layers)) {
      layer_name <- sprintf("issue_%02d", i)
      sf::st_write(issue_layers[[i]], issues_gpkg, layer = layer_name,
                   append = i > 1L, quiet = TRUE)
    }
  }

  utils::write.csv(summary_row,     summary_csv,         row.names = FALSE)
  utils::write.csv(pburned_overall, pburned_overall_csv, row.names = FALSE)
  if (nrow(pburned_by_class) > 0)
    utils::write.csv(pburned_by_class,  pburned_by_class_csv,  row.names = FALSE)
  if (nrow(pburned_by_source) > 0)
    utils::write.csv(pburned_by_source, pburned_by_source_csv, row.names = FALSE)

  counts_det <- D |>
    sf::st_drop_geometry() |>
    dplyr::count(.data$class_final, name = "n_det")
  counts_full <- F_full |>
    sf::st_drop_geometry() |>
    dplyr::count(.data$class_final, name = "n_final_map_full")
  counts_public <- F_public |>
    sf::st_drop_geometry() |>
    dplyr::count(.data$class_final, name = "n_final_map")
  count_table <- counts_det |>
    dplyr::full_join(counts_full,   by = "class_final") |>
    dplyr::full_join(counts_public, by = "class_final")

  lines_out <- c(
    sprintf("Consistency check | year=%s | scenario=%s", target_year, scenario),
    "",
    sprintf("internal_decisions rows: %d", nrow(D)),
    sprintf("final_map_full rows: %d", nrow(F_full)),
    sprintf("final_map rows: %d", nrow(F_public)),
    "",
    sprintf("missing_in_final_map_full: %d", length(missing_in_full_ids)),
    sprintf("extra_in_final_map_full: %d", length(extra_in_full_ids)),
    sprintf("invalid_p_burned_full: %d", nrow(invalid_pb)),
    sprintf("class_final_mismatch: %d", nrow(class_mismatch)),
    sprintf("filter_signature_mismatch: %d", nrow(filter_signature_mismatch)),
    sprintf("expected_public_drop: %d", sum(expected_drop_public, na.rm = TRUE)),
    sprintf("unexpected_missing_in_public: %d", length(unexpected_missing_public_ids)),
    sprintf("unexpected_kept_in_public: %d", length(unexpected_kept_public_ids)),
    "",
    "p_burned overall:",
    sprintf("min=%.6f | q05=%.6f | median=%.6f | q95=%.6f | max=%.6f | sd=%.6f",
            pburned_overall$p_burned_min,
            pburned_overall$p_burned_q05,
            pburned_overall$p_burned_median,
            pburned_overall$p_burned_q95,
            pburned_overall$p_burned_max,
            pburned_overall$p_burned_sd),
    sprintf("n<0.10=%d | n<0.25=%d | n<0.50=%d | n>0.90=%d | n>0.99=%d",
            pburned_overall$n_lt_0_10,
            pburned_overall$n_lt_0_25,
            pburned_overall$n_lt_0_50,
            pburned_overall$n_gt_0_90,
            pburned_overall$n_gt_0_99),
    sprintf("degenerate_high_scores=%s | degenerate_low_scores=%s | degenerate_flat_scores=%s",
            degenerate_high, degenerate_low, degenerate_flat),
    "",
    "Counts by class_final:"
  )
  lines_out <- c(lines_out,
                 utils::capture.output(print(as.data.frame(count_table),
                                             row.names = FALSE)))
  if (nrow(pburned_by_class) > 0) {
    lines_out <- c(lines_out, "", "p_burned by class_input/class_final:",
                   utils::capture.output(print(as.data.frame(pburned_by_class),
                                               row.names = FALSE)))
  }
  if (nrow(pburned_by_source) > 0) {
    lines_out <- c(lines_out, "", "p_burned by source_set:",
                   utils::capture.output(print(as.data.frame(pburned_by_source),
                                               row.names = FALSE)))
  }
  if (length(issue_layers) > 0) {
    lines_out <- c(lines_out, "",
                   sprintf("Issue layers written to: %s", issues_gpkg))
  }
  writeLines(lines_out, summary_txt)

  summary_row
}
