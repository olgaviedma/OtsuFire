#' Assign deterministic filters and final decisions to candidate burned patches
#'
#' @description
#' Public modular scoring function. Receives candidate burned patches
#' produced by
#' \code{\link[=detect_burned_patches]{detect_burned_patches()}} and
#' applies the canonical deterministic decision workflow.
#'
#' The scoring workflow evaluates candidate patches through a sequential
#' rule-based system combining: (i) internal seed support and
#' burnable-context plausibility (Filter 1), (ii) optional active-fire
#' hotspot corroboration when hotspot data are available (Filter 2),
#' (iii) temporal-conflict assessment against `previous_year_burned`, and
#' (iv) spectral-support evaluation relative to a high-confidence
#' keep-like reference distribution (Filter 3). The outputs of these
#' components are then integrated into a final deterministic confidence
#' decision assigning each patch to `keep`, `review`, or `drop`.
#'
#' Spectral support is evaluated using patch-level change-index metrics
#' derived from a high-confidence keep-like reference pool. This reference
#' is supplied explicitly through `keep_pool` or, when omitted, is
#' automatically constructed from high-confidence current-year candidate
#' patches through the deterministic scoring local fallback.
#' Deterministic scoring does not read or depend on any implicit external
#' burned-like registry.
#'
#' The function records the full decision path of every candidate patch,
#' including filter outcomes, supporting metrics, temporal-overlap
#' diagnostics, and final decisions, producing an auditable burned-area
#' decision layer.
#'
#' Contract source: `DETERMINISTIC_PUBLIC_FUNCTION_CONTRACTS.csv`,
#' `DETERMINISTIC_OUTPUTS_FINAL.csv`, `DETERMINISTIC_DATA_MODEL_FINAL.csv`.
#'
#' Validation of the decision layer against an external burned reference
#' is NOT part of this function; users should subsequently call
#' \code{\link[=validate_fire_maps]{validate_fire_maps()}} on the
#' resulting `internal_decisions` layer (or on supervised thresholded
#' outputs).
#'
#' @details
#' \strong{Rule-based deterministic scoring workflow}
#'
#' Candidate burned patches are evaluated through a sequential filter
#' system designed to distinguish high-confidence burned patches from
#' ambiguous or spurious detections while preserving full decision
#' traceability.
#'
#' \code{Filter 1} evaluates internal support using burned-seed
#' intersection, burnable-context plausibility, and screening for
#' excessive non-burnable or missing vegetation information. Patches with
#' strong internal support remain on the keep pathway, whereas spectrally
#' ambiguous but burnable candidates are downgraded to review and
#' contextually implausible patches are dropped.
#'
#' \code{Filter 2} evaluates external corroboration using active-fire
#' hotspots when hotspot data are available for the target year.
#' Candidate patches intersecting retained hotspots receive stronger
#' fire-support evidence, whereas unsupported patches remain under review
#' unless already rejected by previous filters. When hotspot products are
#' unavailable, this filter is recorded as not applied.
#'
#' The temporal consistency assessment evaluates overlap with
#' previous-year burned areas. It is a distinct step, separate from
#' \code{Filter 1}, \code{Filter 2}, and \code{Filter 3}: candidate
#' patches conflicting with recent burns may be retained, downgraded, or
#' removed depending on overlap magnitude and remaining valid area after
#' cleanup.
#'
#' \code{Filter 3} evaluates spectral support relative to a
#' high-confidence keep-like reference distribution. The scoring workflow
#' computes metrics describing spectral consistency, spatial coherence of
#' burned signal, and minimum geometric support. Candidate patches with
#' strong support remain keep, borderline cases are reviewed, and
#' spectrally inconsistent patches are dropped.
#'
#' The outputs of all filters and the temporal-consistency assessment are
#' then integrated into a final deterministic confidence decision and
#' associated reasoning fields.
#'
#' \strong{Spectral-support metrics}
#'
#' \code{Filter 3} evaluates patch-level spectral consistency using
#' metrics derived from the keep-support distribution, including:
#' \itemize{
#'   \item \code{percentile_in_keep}: percentile position of the
#'     candidate patch within the burned-support distribution;
#'   \item \code{p_above_keep_ref}: fraction of patch pixels exceeding
#'     the burned reference threshold;
#'   \item \code{p_above_keep_q25}: fraction of patch pixels exceeding a
#'     lower support threshold;
#'   \item \code{area_ha}: patch area in hectares;
#'   \item \code{n_pix}: number of valid change-index pixels.
#' }
#'
#' These metrics are used jointly to determine whether candidate patches
#' show sufficient burned-like spectral support and spatial consistency
#' to remain within the high-confidence burned pool.
#'
#' @param burned_candidates sf POLYGON, SpatVector, or path to a readable
#'   vector file. Candidate burned patches generated by
#'   \code{\link[=detect_burned_patches]{detect_burned_patches()}} after
#'   seed-and-grow segmentation.
#' @param config An `otsufire_burned_mapping_config` object returned by
#'   \code{\link[=build_burned_mapping_config]{build_burned_mapping_config()}}.
#' @param keep_pool Optional explicit keep-like reference pool used for
#'   spectral-support scoring in `Filter 3`. This should be a pre-built
#'   keep-pool summary object matching the contract consumed by
#'   \code{\link[=score_rbr_keep_classes]{score_rbr_keep_classes()}},
#'   for example the object returned by
#'   \code{\link[=build_keep_pool_from_samples]{build_keep_pool_from_samples()}}.
#'
#'   When supplied, the scoring workflow uses this reference directly.
#'   For backward compatibility, a legacy wrapper of the form
#'   \code{list(keep_pool = <summary>)} is also accepted and unwrapped
#'   internally.
#'
#'   When `NULL`, the scoring engine automatically constructs a local
#'   same-year reference distribution from high-confidence current-year
#'   candidate patches (the deterministic scoring local fallback).
#'
#'   Deterministic scoring does not read or depend on external
#'   burned-like registries: support resolution is strictly explicit
#'   `keep_pool` or local fallback.
#' @param write_outputs Logical scalar. Whether scoring products should be
#'   written to disk. Default `TRUE`.
#' @param overwrite Logical scalar. Whether existing scoring outputs may
#'   be replaced. Default `FALSE`.
#'
#' @return A named list with fields:
#' \describe{
#'   \item{`internal_decisions`}{Canonical deterministic decision layer
#'     with full filter outcomes, support metrics, temporal-overlap
#'     diagnostics, and final unsupervised decisions.}
#'   \item{`internal_decisions_path`}{Path to the written internal
#'     decision layer.}
#'   \item{`reference_decisions`}{Reference-side deterministic decision
#'     layer, when available.}
#'   \item{`reference_decisions_path`}{Path to the written reference
#'     decision layer, when available.}
#'   \item{`keep_pool_summary`}{Summary object describing the keep-like
#'     spectral-support reference distribution used during scoring.}
#'   \item{`scoring_diagnostics`}{Auxiliary diagnostics and scoring
#'     metadata generated during deterministic filtering.}
#' }
#'
#' @seealso
#' \code{\link[=build_burned_mapping_config]{build_burned_mapping_config()}},
#' \code{\link[=detect_burned_patches]{detect_burned_patches()}},
#' \code{\link[=run_deterministic_pipeline]{run_deterministic_pipeline()}},
#' \code{\link[=validate_fire_maps]{validate_fire_maps()}}
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
  # on disk — if the config itself is not usable, we fail fast with a
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
