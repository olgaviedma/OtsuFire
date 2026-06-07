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
