# =============================================================================
# Gate 1B / Precision 1 (2026-06-07): deprecated function-level methodological
# parameter SHIMS.
#
# The CANONICAL way to set a supervised methodological / training-control
# parameter is build_supervised_burned_config() -> cfg$train_control /
# cfg$model_params (the single source of truth, Gate 1B keystone). The
# matching function-level arguments on the PUBLIC functions
# (run_oneyear_supervised_pipeline, run_oof_diagnostics,
# train_final_burned_model) are DEPRECATED COMPATIBILITY SHIMS, kept only so
# pre-Gate-1B caller scripts keep working.
#
# This file holds the ONE resolver every public boundary uses to fold a
# function-level override into the resolved cfg value. It guarantees the eight
# Precision-1 acceptance conditions:
#   - resolution happens ONLY at the public boundary (this function is called
#     from the public wrappers; the dispatcher / orchestrator / engines consume
#     the already-resolved scalar);
#   - a non-NULL override over a CANONICAL DEFAULT emits a deprecation warning
#     with the stable class "otsufire_deprecated_param";
#   - a non-NULL override that conflicts with an EXPLICIT user value
#     (provenance "user") AND differs from it -> stop() (two explicit sources);
#   - each application is recorded as a provenance row (cfg value, requested
#     override, resolved value, provenance label) for the run-level artifact.
# =============================================================================

#' Resolve a single deprecated methodological shim at the public boundary.
#'
#' @param override The function-level argument value (NULL = not supplied).
#' @param cfg_value The value already resolved on the cfg
#'   (cfg$train_control / cfg$model_params).
#' @param param Character. The user-facing builder field name (the provenance
#'   key, e.g. "cap_spectral", "nrounds_max", "val_frac"), used in messages.
#' @param arg_name Character. The function-level argument name (e.g.
#'   "spectral_hard_negative_to_burned_ratio"), used in messages. Defaults to
#'   `param`.
#' @param cfg_provenance Character. "user" or "default" — the builder-time
#'   provenance of the cfg field (from cfg$resolved_params_provenance).
#' @param record An environment used to accumulate provenance rows; each
#'   applied/observed field appends one entry. May be NULL (no recording).
#'
#' @return The resolved value (the override when non-NULL, else the cfg value).
#'
#' @details
#' Equality for the conflict check uses `isTRUE(all.equal(...))` so numeric
#' near-equality (e.g. 2.0 vs 2L) does not spuriously error. A non-NULL override
#' that is *equal* to an explicit user value is accepted silently (no conflict,
#' no deprecation warning — the two explicit sources agree).
#'
#' @keywords internal
#' @noRd
.of_resolve_methodological_shim <- function(override, cfg_value, param,
                                            arg_name = param,
                                            cfg_provenance = "default",
                                            record = NULL) {
  prov <- if (is.null(override)) {
    # No override: the cfg value flows through with its builder provenance.
    if (identical(cfg_provenance, "user")) "cfg-user" else "cfg-default"
  } else {
    same <- isTRUE(all.equal(override, cfg_value))
    if (identical(cfg_provenance, "user") && !same) {
      # Two EXPLICIT, incompatible sources -> hard error (never silently pick).
      stop(
        sprintf(
          paste0(
            "Conflicting supervised parameter '%s': the builder explicitly set ",
            "it (build_supervised_burned_config(%s=...) -> %s) AND a ",
            "function-level override '%s=%s' was passed, and the two DIFFER. ",
            "These are two explicit, incompatible sources. Set the parameter in ",
            "ONE place: prefer build_supervised_burned_config(%s=...) (the ",
            "canonical path) and drop the deprecated function-level override."),
          param, param, .of_shim_fmt(cfg_value), arg_name,
          .of_shim_fmt(override), param),
        call. = FALSE
      )
    }
    if (!same) {
      # Override applied over a CANONICAL DEFAULT -> deprecation warning.
      warning(
        structure(
          class = c("otsufire_deprecated_param", "warning", "condition"),
          list(
            message = sprintf(
              paste0(
                "The function-level argument '%s' is DEPRECATED. It overrides ",
                "the canonical default for '%s' (%s -> %s). Set this parameter ",
                "in build_supervised_burned_config(%s=...) instead; the ",
                "function-level shim will be removed in a future version."),
              arg_name, param, .of_shim_fmt(cfg_value),
              .of_shim_fmt(override), param),
            call = NULL
          )
        )
      )
      "function-override"
    } else {
      # Override equals the cfg value -> no conflict, no deprecation noise.
      if (identical(cfg_provenance, "user")) "cfg-user" else "function-override-noop"
    }
  }

  resolved <- override %||% cfg_value

  if (!is.null(record) && is.environment(record)) {
    record$rows[[length(record$rows) + 1L]] <- list(
      param          = param,
      arg_name       = arg_name,
      cfg_value      = .of_shim_fmt(cfg_value),
      cfg_provenance = cfg_provenance,
      requested      = if (is.null(override)) NA_character_ else .of_shim_fmt(override),
      resolved       = .of_shim_fmt(resolved),
      provenance     = prov
    )
  }

  resolved
}

#' Compact one-line formatter for a parameter value (for messages/records).
#' @keywords internal
#' @noRd
.of_shim_fmt <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.numeric(x) && length(x) == 1L) return(format(x, trim = TRUE))
  if (is.character(x) && length(x) == 1L) return(x)
  if (length(x) <= 6L) return(paste(format(x, trim = TRUE), collapse = ","))
  paste0("<", class(x)[1L], "[", length(x), "]>")
}

#' Build a fresh provenance-record accumulator environment.
#' @keywords internal
#' @noRd
.of_new_shim_record <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e
}

# =============================================================================
# Precision 2 (2026-06-07): PACKAGE-LEVEL spectral-cap parity guard.
#
# The Phase B experiments train the supervised model with a spectral
# hard-negative cap of 2.0 (the general package default is 1.0, set in
# .of_canonical_train_control()). The methodological requirement is that BOTH
# the OOF diagnostics stage and the FINAL model train under the SAME spectral
# cap (otherwise OOF metrics describe a different model family than the deployed
# one). The B1_PHASE2 runner has a cap-mismatch guard, but that lives in the
# RUNNER. This guard lives in the PACKAGE so the supervised pipeline itself
# fails fast on a silent cap reversion, regardless of which runner drives it.
# =============================================================================

#' Assert the spectral cap reaching OOF and FINAL agree with the resolved cfg.
#'
#' Aborts (stop) when the spectral hard-negative cap value about to be handed to
#' the OOF stage differs from the one about to be handed to the FINAL stage, or
#' when either differs from the cap resolved on the cfg. This detects a silent
#' reversion (e.g. a residual default quietly returning the cap to 1.0) on the
#' capped training path before any model is fit.
#'
#' @param oof_spectral Numeric. Spectral cap forwarded to the OOF stage.
#' @param final_spectral Numeric. Spectral cap forwarded to the FINAL stage.
#' @param cfg_spectral Numeric. Spectral cap resolved on
#'   `cfg$train_control$caps$spectral`.
#' @param tol Numeric tolerance for the floating-point comparison.
#'
#' @return Invisibly the agreed spectral cap value (errors otherwise).
#' @keywords internal
#' @noRd
.of_assert_spectral_cap_parity <- function(oof_spectral, final_spectral,
                                           cfg_spectral, tol = 1e-9) {
  bad <- function(v) is.null(v) || !is.numeric(v) || length(v) != 1L || is.na(v)
  if (bad(oof_spectral) || bad(final_spectral) || bad(cfg_spectral)) {
    stop(".of_assert_spectral_cap_parity(): spectral caps must each be a single ",
         "finite numeric (OOF / FINAL / cfg).", call. = FALSE)
  }
  oof_final_ok <- abs(oof_spectral - final_spectral) <= tol
  oof_cfg_ok   <- abs(oof_spectral - cfg_spectral)   <= tol
  final_cfg_ok <- abs(final_spectral - cfg_spectral) <= tol
  if (!(oof_final_ok && oof_cfg_ok && final_cfg_ok)) {
    stop(sprintf(
      paste0(
        "Spectral-cap parity guard FAILED (silent reversion detected): the ",
        "spectral hard-negative cap reaching OOF (%s), FINAL (%s) and the ",
        "resolved cfg (cfg$train_control$caps$spectral=%s) are NOT all equal. ",
        "OOF and FINAL on the capped training path MUST receive the SAME ",
        "spectral cap as the cfg. Rebuild the cfg with ",
        "build_supervised_burned_config(cap_spectral=...) and do not let a ",
        "function-level override or residual default revert it."),
      format(oof_spectral, trim = TRUE), format(final_spectral, trim = TRUE),
      format(cfg_spectral, trim = TRUE)),
      call. = FALSE)
  }
  invisible(oof_spectral)
}

#' Convert an accumulated shim record into a tidy data.frame for artifacts.
#' @keywords internal
#' @noRd
.of_shim_record_to_df <- function(record) {
  if (is.null(record) || length(record$rows) == 0L) {
    return(data.frame(
      param = character(0), arg_name = character(0), cfg_value = character(0),
      cfg_provenance = character(0), requested = character(0),
      resolved = character(0), provenance = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(record$rows, function(r) {
    data.frame(r, stringsAsFactors = FALSE)
  }))
}
