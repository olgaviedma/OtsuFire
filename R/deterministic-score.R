#' Score candidate burned patches and assign final decisions
#'
#' @description
#' Evaluate candidate burned patches using the deterministic OtsuFire scoring
#' rules and assign each patch a final decision: `"keep"`, `"review"`, or
#' `"drop"`.
#'
#' The function combines seed support, burnable-land context, optional
#' active-fire hotspot evidence, temporal assessment, and spectral support.
#' Filter outcomes, supporting metrics, and decision reasons are retained in
#' the output layer.
#'
#' Use this function after [detect_burned_patches()]. Validation against an
#' external burned-area map is performed separately with
#' [validate_fire_maps()].
#'
#' @param burned_candidates Candidate polygons produced by the detection and
#'   refinement workflow. Accepts an `sf` polygon object, a polygon
#'   `terra::SpatVector`, or a path to a readable vector file. The
#'   `refined_patches_path` returned by [detect_burned_patches()] can be
#'   supplied directly.
#' @param config An object of class `otsufire_burned_mapping_config`, created
#'   with [build_burned_mapping_config()]. Contains the supporting inputs,
#'   scoring settings, temporal-cleanup settings, and output locations.
#' @param keep_pool Optional spectral-support reference object compatible with
#'   [score_rbr_keep_classes()]. Explicit pools must contain reference samples
#'   from years other than the target year. When `NULL`, the function
#'   constructs a local reference from the early high-confidence candidate
#'   core for the target year. Default: `NULL`.
#' @param write_outputs Logical scalar. Whether to write scoring products to
#'   disk. When `FALSE`, returned path fields may be `NA`. Default: `TRUE`.
#' @param overwrite Logical scalar. Whether existing output files may be
#'   replaced. Default: `FALSE`.
#'
#' @section Scoring workflow:
#' Candidates are evaluated through three filters and a separate temporal
#' assessment.
#'
#' | Component | Evidence assessed | Output fields |
#' |---|---|---|
#' | Filter 1 | Seed support, burnable-land context, vegetation availability, and non-burnable contamination. | `filter_1`, `reason_1` |
#' | Filter 2 | Corroboration from active-fire hotspots, when available. | `filter_2`, `reason_2` |
#' | Temporal assessment | Overlap and conflicts with previous-year burned areas. | `preyear_overlap_frac`, `preyear_action`, `preyear_reason` |
#' | Filter 3 | Spectral support relative to a high-confidence reference distribution, together with minimum area and pixel support. | `filter_3`, `reason_3` |
#'
#' Final decisions are stored in `class_final`. Individual filter outcomes
#' represent intermediate assessments.
#'
#' @section Filter 1 (seed support and burnable-land context):
#' Filter 1 assesses internal seed support and the plausibility of the
#' candidate's land-cover context.
#'
#' Seed-supported candidates remain on the keep pathway. Candidates without
#' seed support are assessed against the burnable domain:
#' * burnable but unsupported candidates are assigned to review;
#' * non-burnable candidates are assigned to drop;
#' * candidates with excessive missing vegetation information may also be
#'   dropped.
#'
#' This filter is independent of the spectral-support reference.
#'
#' @section Filter 2 (hotspot corroboration):
#' When hotspot observations are available for the target year:
#' * hotspot-supported candidates receive `filter_2 = "keep"`;
#' * candidates without hotspot support receive `filter_2 = "review"`.
#'
#' The absence of hotspot support does not, by itself, produce a drop
#' decision in this filter.
#'
#' When hotspot data are unavailable, Filter 2 is recorded as not applied.
#'
#' @section Temporal assessment:
#' When a previous-year burned-area map is supplied, overlap is assessed to
#' identify potential residual scars, repeated detections, and temporal
#' conflicts.
#'
#' Overlap fractions, temporal actions, and their reasons are retained
#' separately from the numbered filters.
#'
#' @section Filter 3 (spectral support):
#' Filter 3 compares candidate spectral responses with a high-confidence
#' reference distribution.
#'
#' The following metrics are recorded:
#'
#' | Field | Description |
#' |---|---|
#' | `percentile_in_keep` | Percentile position within the spectral-support reference distribution. |
#' | `p_above_keep_ref` | Fraction of pixels exceeding the reference support threshold. |
#' | `p_above_keep_q25` | Fraction of pixels exceeding the lower reference support threshold. |
#' | `area_ha` | Patch area in hectares. |
#' | `n_pix` | Number of valid change-index pixels. |
#' | `conf_area` | Flag indicating whether minimum area requirements are met. |
#' | `conf_pix` | Flag indicating whether minimum pixel-support requirements are met. |
#'
#' Candidates receive `filter_3 = "keep"` when they meet the percentile and
#' exceedance criteria and satisfy the minimum area and pixel-support
#' requirements.
#'
#' Borderline candidates generally receive `filter_3 = "review"`. Drop
#' outcomes are mainly associated with insufficient geometric support, such
#' as `conf_area = "too_small"` or `conf_pix = "too_few_pix"`.
#'
#' @section Spectral-support references:
#' Two reference modes are supported.
#'
#' \strong{Local reference from the target year.} Use `keep_pool = NULL`, the
#' default.
#'
#' The function constructs a reference from the early, seed-supported
#' candidate core: polygons entering the intermediate keep/review layer with
#' `flag_internal = "keep"`.
#'
#' This reference is constructed before final confidence assignment. It is
#' not selected from the final `class_final = "keep"` output and does not use
#' an external burned-area registry or validation map.
#'
#' The reference provides the thresholds and empirical distribution used by
#' Filter 3, including `qref_keep`, `q25_keep`, and `keep_medians`.
#'
#' \strong{Explicit reference from other years.} Supply a reference object
#' prepared from trusted decision layers for years other than
#' `config$target_year`.
#'
#' Use [collect_keep_reference_samples()] and [build_keep_pool_from_samples()]
#' to prepare reference samples and summaries. The supplied object must match
#' the structure expected by [score_rbr_keep_classes()].
#'
#' When an explicit pool is supplied, Filter 3 uses it directly instead of
#' constructing the local reference.
#'
#' The standard helper workflow selects reference samples meeting all of the
#' following conditions:
#' * `class_final == "keep"`;
#' * `filter_1 == "keep"`;
#' * `filter_2 == "keep"`;
#' * `filter_3 == "keep"`;
#' * `conf_area == "ok"`;
#' * `conf_pix == "ok"`;
#' * compliance with the minimum-area requirement;
#' * no temporal drop when `preyear_action` is available.
#'
#' These selection rules differ from those used for the local reference.
#'
#' Explicit pools containing the target year are not supported and are
#' rejected when helper-generated metadata identify that overlap. Ensure that
#' the target year is excluded when preparing an explicit pool.
#'
#' Legacy objects of the form `list(keep_pool = summary)` are also accepted.
#'
#' @section Final decisions:
#' Filter outcomes are integrated into the `class_final` field.
#'
#' | Value | Interpretation |
#' |---|---|
#' | `"keep"` | High-confidence burned candidate under the deterministic rules. |
#' | `"review"` | Ambiguous candidate or candidate with incomplete or conflicting support. |
#' | `"drop"` | Candidate rejected by the deterministic rules. |
#'
#' If exactly one active filter assigns drop and all remaining active filters
#' assign keep, the final decision is softened to review.
#'
#' These classes are rule-based decisions, not estimated probabilities of
#' burning. Assess map accuracy separately using suitable reference data.
#'
#' @return A named list containing scoring results, output paths, and
#'   diagnostics.
#'
#' | Field | Description |
#' |---|---|
#' | `internal_decisions` | Main candidate decision layer, containing filter outcomes, support metrics, temporal diagnostics, decision reasons, and `class_final`. |
#' | `internal_decisions_path` | Path to the written candidate decision layer. |
#' | `reference_decisions` | Reference-side deterministic decision layer, when available. |
#' | `reference_decisions_path` | Path to the written reference-side decision layer, when available. |
#' | `keep_pool_summary` | Summary of the spectral-support reference used during scoring. |
#' | `scoring_diagnostics` | Additional diagnostics and scoring metadata. |
#'
#' Path fields may be `NA` when the corresponding products are not written.
#' The reference-side decision outputs do not constitute external map
#' validation.
#'
#' @seealso [build_burned_mapping_config()], [detect_burned_patches()],
#'   [run_deterministic_pipeline()], [collect_keep_reference_samples()],
#'   [build_keep_pool_from_samples()], [score_rbr_keep_classes()],
#'   [validate_fire_maps()].
#'
#' @examples
#' \dontrun{
#' # Configure an annual RBR workflow
#' config <- build_burned_mapping_config(
#'   change_index = "data/RBR_2022.tif",
#'   vegetation_map = "data/vegetation_classes.tif",
#'   burnable_mask = "data/burnable_mask.tif",
#'   hotspots = "data/hotspots_2022.gpkg",
#'   previous_year_burned = "data/burned_2021.gpkg",
#'   target_year = 2022L,
#'   output_dir = "results",
#'   run_name = "RBR_2022"
#' )
#'
#' # Detect and refine candidate patches
#' detection <- detect_burned_patches(
#'   config = config,
#'   write_outputs = TRUE
#' )
#'
#' # Score candidates using the local target-year reference
#' scored <- score_burned_patches(
#'   burned_candidates = detection$refined_patches_path,
#'   config = config,
#'   keep_pool = NULL,
#'   write_outputs = TRUE,
#'   overwrite = FALSE
#' )
#'
#' # Summarise final decisions
#' table(
#'   scored$internal_decisions$class_final,
#'   useNA = "ifany"
#' )
#'
#' # Inspect diagnostics and the decision-layer path
#' scored$scoring_diagnostics
#' scored$internal_decisions_path
#'
#' # Use a previously prepared reference from other years
#' # The saved object must have the required keep-pool structure
#' # and must exclude the target year, 2022.
#' historical_pool <- readRDS("data/keep_pool_other_years.rds")
#'
#' scored_external <- score_burned_patches(
#'   burned_candidates = detection$refined_patches_path,
#'   config = config,
#'   keep_pool = historical_pool,
#'   write_outputs = FALSE
#' )
#'
#' table(
#'   scored_external$internal_decisions$class_final,
#'   useNA = "ifany"
#' )
#' }
#'
#' @export
score_burned_patches <- function(burned_candidates, config, keep_pool = NULL,
                                 write_outputs = TRUE, overwrite = FALSE) {

  if (!inherits(config, "otsufire_burned_mapping_config")) {
    stop("'config' must be created by build_burned_mapping_config().",
         call. = FALSE)
  }
  if (!is.logical(write_outputs) || length(write_outputs) != 1L) {
    stop("'write_outputs' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  if (missing(burned_candidates) || is.null(burned_candidates)) {
    stop("'burned_candidates' is required.", call. = FALSE)
  }

  # Validate config content (change_index) before touching burned_candidates
  # on disk - if the config itself is not usable, we fail fast with a
  # config-level error instead of a misleading path-existence error.
  if (is.null(config$inputs$change_index)) {
    stop(
      "score_burned_patches(): config$inputs$change_index is NULL.\n",
      "Rebuild the config with an explicit 'change_index' argument.",
      call. = FALSE
    )
  }
  .of_check_input_file(config$inputs$change_index, "change_index")

  if (is.character(burned_candidates) && length(burned_candidates) == 1L) {
    if (!file.exists(burned_candidates)) {
      stop("'burned_candidates' path does not exist: ", burned_candidates,
           call. = FALSE)
    }
  } else if (!inherits(burned_candidates, c("sf", "SpatVector"))) {
    stop(
      "'burned_candidates' must be an sf POLYGON, SpatVector, or a readable path.",
      call. = FALSE
    )
  }

  keep_pool <- .of_normalize_keep_pool_arg(keep_pool)
  keep_pool <- .of_validate_keep_pool_target_year(
    keep_pool,
    config$target_year
  )

  res <- .of_run_scoring(burned_candidates, config, keep_pool = keep_pool,
                         write_outputs = isTRUE(write_outputs),
                         overwrite = isTRUE(overwrite))

  list(
    internal_decisions       = res$internal_decisions,
    internal_decisions_path  = res$internal_decisions_path %||% NA_character_,
    reference_decisions      = res$reference_decisions,
    reference_decisions_path = res$reference_decisions_path %||% NA_character_,
    keep_pool_summary        = res$keep_pool_summary,
    scoring_diagnostics      = res$scoring_diagnostics
  )
}
