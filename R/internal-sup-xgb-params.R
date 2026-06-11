# =============================================================================
# Canonical XGBoost hyperparameters for the supervised burned-area model.
#
# FRENTE 1 (2026-06-05): SINGLE SOURCE OF TRUTH for the XGBoost parameter
# block shared by BOTH the OOF diagnostics stage (run_oof_diagnostics ->
# run_dm_oof_pipeline -> run_oof_xgb) and the FINAL-model stage
# (train_final_burned_model -> train_final_model_direct). Before this, each
# stage built its own inline `params` list when `params = NULL`, and the two
# lists had DIVERGED (OOF: eta 0.06 / subsample 0.85 / no lambda+alpha;
# FINAL: eta 0.05 / depth 6 / min_child_weight 1 / colsample 0.8 / single
# "logloss" eval_metric). Natalia signed off on UNIFYING them to one
# canonical set so the OOF diagnostics report on the SAME model family the
# FINAL stage actually fits. This is RESULT-AFFECTING and intentional: the
# FINAL model changes vs the historical baseline; it is first evaluated in
# the 2017 run. See HANDOFF section "FRENTE 1 RESOLVED".
#
# Both inline-`params` sites call THIS builder when `params = NULL`, so the
# two can never silently diverge again. `scale_pos_weight` is the ONLY
# field that is NOT fixed here: each caller computes it from its OWN training
# labels (OOF: per OOF training set; FINAL: per its train split) and passes
# it in. Everything else is fixed.
# =============================================================================

#' Canonical XGBoost params for the supervised burned model (single source of truth).
#'
#' @param scale_pos_weight Numeric. Class-imbalance weight computed by the
#'   caller from its OWN training labels, i.e.
#'   `sum(y == 0) / max(1, sum(y == 1))`. Not defaulted here on purpose: it
#'   must reflect the actual training set at each fit site.
#' @return A named list of XGBoost parameters, identical for OOF and FINAL
#'   except for `scale_pos_weight`.
#' @keywords internal
#' @noRd
.of_canonical_xgb_params <- function(scale_pos_weight) {
  if (missing(scale_pos_weight) || !is.numeric(scale_pos_weight) ||
      length(scale_pos_weight) != 1L || !is.finite(scale_pos_weight)) {
    stop(".of_canonical_xgb_params(): 'scale_pos_weight' must be a single ",
         "finite numeric computed from the training labels.", call. = FALSE)
  }
  .of_canonical_xgb_block(scale_pos_weight = scale_pos_weight)
}

#' Canonical XGBoost block WITHOUT scale_pos_weight (Gate 1B single source).
#'
#' Gate 1B (2026-06-07): the methodological xgb hyperparameters live in exactly
#' ONE place. `build_supervised_burned_config()` calls this (with the sentinel
#' weight) and strips `scale_pos_weight` to populate `cfg$model_params`, so the
#' cfg block and the engine block are sourced from the same definition and can
#' never silently diverge.
#'
#' @param scale_pos_weight Numeric scalar merged into the returned block.
#' @return Named list of XGBoost parameters.
#' @keywords internal
#' @noRd
.of_canonical_xgb_block <- function(scale_pos_weight) {
  list(
    booster          = "gbtree",
    objective        = "binary:logistic",
    # FRENTE 1: eval_metric is a VECTOR with logloss FIRST so early stopping
    # is driven by logloss (which does not saturate); aucpr is kept second so
    # it stays visible in the per-iteration eval_log only. The FINAL stage
    # previously used a single "logloss"; the OOF stage already used this
    # vector. NOTE (handoff): an earlier aucpr-FIRST ordering made AUCPR drive
    # early stopping and collapsed best_iter to 1 (AUCPR saturates instantly on
    # this imbalanced set); logloss-first restores a calibrated best_iter ~137.
    eval_metric      = c("logloss", "aucpr"),
    eta              = 0.05,
    max_depth        = 5,
    min_child_weight = 5,
    subsample        = 0.8,
    colsample_bytree = 0.75,
    gamma            = 0,
    lambda           = 1,
    alpha            = 0,
    scale_pos_weight = scale_pos_weight
  )
}

# =============================================================================
# Gate 1B (2026-06-07): SINGLE SOURCE OF TRUTH for the supervised
# methodological / training-control parameters.
#
# Before Gate 1B these defaults were duplicated across 5 layers (the public
# run_oneyear_supervised_pipeline, the internal dispatch / orchestrator, the
# OOF / FINAL public wrappers, and the internal-pure engines). They now live
# in EXACTLY two builders below, which `build_supervised_burned_config()`
# resolves into `cfg$model_params` and `cfg$train_control`. Every downstream
# consumer reads those resolved cfg sections; the internal-pure engines take
# the resolved values as REQUIRED args with no methodological defaults.
# =============================================================================

#' Canonical methodological XGBoost block (Gate 1B), WITHOUT scale_pos_weight.
#'
#' This is the block stored in `cfg$model_params`. It is the canonical xgb
#' block (`.of_canonical_xgb_block()`) with the site-specific
#' `scale_pos_weight` removed: spw is computed at train time from the actual
#' training labels at each fit site and is therefore NOT a cfg field.
#'
#' @return Named list of XGBoost params minus `scale_pos_weight`.
#' @keywords internal
#' @noRd
.of_canonical_model_params <- function() {
  blk <- .of_canonical_xgb_block(scale_pos_weight = 1)
  blk[["scale_pos_weight"]] <- NULL
  blk
}

#' Canonical supervised training-control defaults (Gate 1B single source).
#'
#' The ONE place the supervised training controls are defined. Resolved into
#' `cfg$train_control` by `build_supervised_burned_config()`. The three seeds
#' are kept distinct on purpose (they have always been 42 each, but address
#' different RNG streams): `oof_seed_base` (block-CV OOF), `final_sampling_seed`
#' (FINAL negative-pool sampling) and `final_seed` (FINAL train/val split +
#' xgboost training).
#'
#' @return Named list with: nrounds_max, early_stop, seeds (list of
#'   oof_seed_base / final_sampling_seed / final_seed), val_frac, group_col,
#'   impute_numeric, impute_factor_missing, caps (list of contextual / random
#'   / otsu), feature_whitelist_override, feature_weights,
#'   training_protocol (fixed internal constant "nested_refit"), oof_sampling.
#' @keywords internal
#' @noRd
.of_canonical_train_control <- function() {
  list(
    nrounds_max = 4000L,
    early_stop  = 80L,
    seeds = list(
      oof_seed_base       = 42L,
      final_sampling_seed = 42L,
      final_seed          = 42L
    ),
    val_frac    = 0.15,
    group_col   = "block_id",
    impute_numeric        = "median",
    impute_factor_missing = "MISSING",
    caps = list(
      contextual = 0.25,
      random     = 1.0,
      otsu       = 1.0
    ),
    feature_whitelist_override = NULL,
    feature_weights            = NULL,
    # Fixed internal constant (traceability / provenance only; NOT user-settable
    # and NOT an argument). OtsuFire always uses inner-early-stopping selection +
    # full-data refit. Kept as a constant so manifests / run fingerprints stay
    # self-documenting and stable.
    training_protocol = "nested_refit",
    oof_sampling      = "capped"
  )
}
