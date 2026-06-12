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
#'   key, e.g. "cap_random", "nrounds_max", "val_frac"), used in messages.
#' @param arg_name Character. The function-level argument name (e.g.
#'   "random_to_burned_ratio"), used in messages. Defaults to
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

#' Render the per-field cfg$model_params provenance as a tidy data.frame.
#'
#' Gate 1B (2026-06-07): `cfg$resolved_params_provenance$model_params` is a named
#' list keyed by xgb field, each element a list with `canonical`, `requested`,
#' `resolved`, `provenance` (built by `.of_model_params_provenance()` in the
#' builder, since the builder folds a PARTIAL override onto the canonical block).
#' This flattens it to one row per field for the run-level manifest, mirroring
#' `.of_shim_record_to_df()`. Renders the exact provenance table a downstream
#' manifest documents (param / canonical / requested / resolved / provenance).
#'
#' @param mp_prov The `cfg$resolved_params_provenance$model_params` list (or
#'   NULL).
#' @return A data.frame with columns param, canonical, requested, resolved,
#'   provenance (zero rows when `mp_prov` is NULL/empty).
#' @keywords internal
#' @noRd
.of_model_params_provenance_to_df <- function(mp_prov) {
  empty <- data.frame(
    param = character(0), canonical = character(0), requested = character(0),
    resolved = character(0), provenance = character(0), stringsAsFactors = FALSE
  )
  if (is.null(mp_prov) || length(mp_prov) == 0L) return(empty)
  rows <- lapply(names(mp_prov), function(nm) {
    r <- mp_prov[[nm]]
    data.frame(
      param      = nm,
      canonical  = r$canonical  %||% NA_character_,
      requested  = r$requested  %||% NA_character_,
      resolved   = r$resolved   %||% NA_character_,
      provenance = r$provenance %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
