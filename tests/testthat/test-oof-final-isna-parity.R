# ============================================================================
# Gate 1D.7 (2026-06-08): OOF <-> FINAL `_isNA` COMPANION PARITY, END-TO-END
# FROM ONE `extract_supervised_features` OUTPUT.
#
# Gate 1D.4 proved OOF and FINAL share the resolved recipe/params and call the
# SAME leakage-free core. It explicitly DEFERRED one boundary check (see the
# NOTE at test-oof-final-symmetry.R around line 343): the `_isNA` companion
# column parity. 1D.4 asserted only the BASE feature set (whitelist members,
# `_isNA` stripped) was identical, because the two captured `feature_cols`
# vectors differed in their `_isNA` companions FOR THAT FIXTURE -- the OOF chain
# ran build_design_matrix_patches() (which SYNTHESISES an `_isNA` flag per
# numeric feature) while the FINAL engine read `_isNA` companions straight from
# its input GPKG, which carried none. 1D.4 attributed this to "how each path is
# fed here, NOT a recipe asymmetry", reasoning that in a real run BOTH stages
# consume the SAME features GPKG produced once by extract_supervised_features(),
# so the `_isNA` set would be identical too.
#
# THIS suite RE-ASSERTS that deferred boundary directly: it drives BOTH the OOF
# design-matrix build and the FINAL refit-core from ONE upstream feature frame
# (the analogue of a single extract_supervised_features() GPKG feeding both
# stages) and asks whether the deployed `_isNA` companion SET + the full design
# matrix x_cols are IDENTICAL across the two stages.
#
# RESULT (recorded for the director): they are NOT identical. The `_isNA`
# companion synthesis is ASYMMETRIC:
#   * OOF  -> build_design_matrix_patches() (R/internal-sup-create-matrix.R:156-165)
#            SYNTHESISES `<feat>_isNA` for EVERY numeric feature; these become the
#            OOF design-matrix model_cols the per-fold core trains on.
#   * FINAL -> feat_cols = .filter_to_supervised_whitelist(names(L_df))
#            (R/internal-sup-train-final-direct.R:332) picks up ONLY the `_isNA`
#            columns PHYSICALLY PRESENT in the input frame. The shared nested
#            core (.of_nested_build_matrix, R/internal-sup-nested-refit.R:173)
#            does NOT synthesise any `_isNA`; it consumes feature_cols verbatim.
#   * extract_supervised_features() does NOT write `_isNA` companions into the
#            features GPKG (the synthesis is documented as a downstream
#            matrix-builder step, R/internal-sup-extract-features.R:431-437), so
#            FINAL receives ZERO `_isNA` columns from the same upstream frame.
#
# NET: fed the SAME extract_supervised_features() output, OOF trains a design
# matrix with 51 `_isNA` companions while FINAL trains one with 0. The deployed
# feature spaces are NOT identical -- a real `_isNA` parity defect at the
# OOF/FINAL boundary, exactly the one 1D.4 flagged as not re-asserted.
#
# Per the Gate 1D.7 scope limit (contract-test/consolidation; do NOT redesign):
# this suite (a) DOCUMENTS the intended parity contract with an executable check
# that currently FAILS-BY-INTENT (xfail), recording the exact divergence, and
# (b) PROVES the base-feature parity that DOES hold. The fix (make FINAL
# synthesise `_isNA` symmetrically, or make extract_supervised_features() write
# the companions once for both stages) is a methodological change left to the
# director.
# ============================================================================

ns_ip <- asNamespace("OtsuFire")
g_ip  <- function(nm) get(nm, envir = ns_ip)

# One upstream feature frame == the analogue of a single
# extract_supervised_features() GPKG feeding BOTH stages. Carries the canonical
# whitelist features (some with injected NA so `_isNA` companions are
# meaningful) + the admin / fold-assignment columns both stages preserve. NO
# `_isNA` columns are pre-written (extract_supervised_features() does not write
# them), so each stage's synthesis behaviour is what's under test.
make_isna_upstream_df <- function(n = 60L, seed = 7L) {
  set.seed(seed)
  wl <- g_ip(".supervised_feature_cols")
  feat <- as.data.frame(
    matrix(stats::runif(n * length(wl)), nrow = n,
           dimnames = list(NULL, wl)))
  # Inject NA into a spread of features (RBR / topo / hotspot) so the synthesised
  # companions are non-degenerate and the asymmetry is visible.
  for (t in c("rbr_med", "elev_med", "slope_med", "hs_conf_mean", "hs_frp_sum")) {
    feat[[t]][sample.int(n, 12L)] <- NA
  }
  admin <- data.frame(
    fire_uid  = sprintf("uid_%03d", seq_len(n)),
    class     = rep(c("burned", "unburned"), length.out = n),
    source    = rep(c("burned_truth", "random_burnable_background",
                      "deterministic_drop_hard", "otsu_patch_residual"),
                    length.out = n),
    neg_type  = rep(c(NA_character_, "background_cell",
                      "spectral_reject_medium", "otsu_patch_drop"),
                    length.out = n),
    poly_id   = sprintf("p_%03d", seq_len(n)),
    block_id  = rep(seq_len(12), length.out = n),
    fold_rep1 = rep(c(1L, 2L, 3L), length.out = n),
    fold_rep2 = rep(c(2L, 3L, 1L), length.out = n),
    stringsAsFactors = FALSE)
  cbind(admin, feat)
}

# Resolve the OOF deployed design-matrix model_cols (what the OOF per-fold core
# trains on) from one upstream frame, replicating the OOF wrapper's whitelist
# filter + deferred design-matrix build.
oof_model_cols_from <- function(df) {
  wl <- g_ip(".supervised_feature_cols")
  allowed <- c(wl, paste0(wl, "_isNA"))
  id_cols <- c("fire_uid", "poly_id", "block_id", "fold_rep1", "fold_rep2",
               "class", "source", "neg_type")
  apply_wl <- function(d) {
    keep <- unique(c(intersect(names(d), id_cols),
                     intersect(names(d), allowed)))
    d[, keep, drop = FALSE]
  }
  labelled    <- apply_wl(df)
  burned_like <- labelled[labelled$class == "unburned", , drop = FALSE]
  dm <- g_ip("build_design_matrix_patches")(
    labelled = labelled, burned_like = burned_like, id_cols = id_cols,
    defer_impute = TRUE, save_dir = NULL, verbose = FALSE)
  dm$model_cols
}

# Resolve the FINAL deployed design-matrix x_cols (what the FINAL refit model is
# trained on + the recipe$cols$x_cols scoring aligns to) from the SAME upstream
# frame, replicating train_final_model_direct()'s feat_cols derivation + the
# shared nested-refit core.
final_x_cols_from <- function(df) {
  wl <- g_ip(".supervised_feature_cols")
  feat_cols <- g_ip(".filter_to_supervised_whitelist")(names(df), whitelist = wl)
  params_fn <- function(spw) {
    p <- g_ip(".of_canonical_xgb_params")(spw); p$nthread <- 1L; p
  }
  fit <- suppressMessages(suppressWarnings(g_ip(".of_nested_refit_fit")(
    train_df = df, feature_cols = feat_cols, label_col = "class",
    group_col = "block_id", block_col = "block_id", val_frac = 0.2,
    params_fn = params_fn, sampling_seed = 42L, fold_seed = 42L,
    nrounds_max = 8L, early_stopping_rounds = 4L,
    impute_numeric = "median", impute_factor_missing = "MISSING",
    verbose = FALSE)))
  fit$x_cols
}

isna_of <- function(v) sort(grep("_isNA$", v, value = TRUE))
base_of <- function(v) v[!grepl("_isNA$", v)]

# ---------------------------------------------------------------------------
# (1) BASE-FEATURE parity DOES hold (the load-bearing 1D.4 guarantee): stripped
#     of `_isNA`, OOF and FINAL train on the SAME canonical whitelist features
#     in the SAME canonical order. This is asserted here end-to-end from ONE
#     upstream frame (not just from the captured args of 1D.4).
# ---------------------------------------------------------------------------
test_that("base feature set (whitelist, _isNA stripped) is identical OOF vs FINAL from one upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  wl <- g_ip(".supervised_feature_cols")

  df <- make_isna_upstream_df()
  oof_cols   <- oof_model_cols_from(df)
  final_cols <- final_x_cols_from(df)

  oof_base   <- base_of(oof_cols)
  final_base <- base_of(final_cols)

  # Same base features (set) ...
  expect_setequal(oof_base, final_base)
  # ... and both follow the canonical whitelist order.
  expect_identical(oof_base,   intersect(wl, oof_base))
  expect_identical(final_base, intersect(wl, final_base))
  # The base sets are byte-identical (same names, same order) across stages.
  expect_identical(oof_base, final_base)
})

# ---------------------------------------------------------------------------
# (2) THE DEFERRED 1D.4 BOUNDARY, re-asserted: the `_isNA` companion SET + the
#     FULL deployed x_cols must be IDENTICAL across OOF and FINAL when both are
#     fed the SAME upstream extract_supervised_features() frame.
#
#     This is the intended contract. It currently FAILS (documented defect):
#     OOF synthesises 51 `_isNA` companions; FINAL synthesises none. We mark it
#     xfail and assert the EXACT divergence so the contract is navigable and the
#     regression is locked until the director fixes the synthesis asymmetry.
# ---------------------------------------------------------------------------
test_that("[xfail/DEFECT] _isNA companion set + x_cols should be IDENTICAL OOF vs FINAL from one upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")

  df <- make_isna_upstream_df()
  oof_cols   <- oof_model_cols_from(df)
  final_cols <- final_x_cols_from(df)

  oof_isna   <- isna_of(oof_cols)
  final_isna <- isna_of(final_cols)

  # --- The INTENDED contract (what a correct, symmetric pipeline guarantees) ---
  # When the engine is fixed, flip these two `expect_false`/divergence asserts to
  # the commented `expect_identical`/`expect_setequal` and delete the xfail note.
  #   expect_setequal(oof_isna, final_isna)          # _isNA companion SET parity
  #   expect_identical(oof_cols, final_cols)          # full x_cols parity (order)

  # --- The DEFECT, asserted explicitly (recorded + reported to the director) ---
  # OOF synthesises a companion for every numeric feature; FINAL synthesises none
  # because extract_supervised_features() writes no `_isNA` columns and the
  # nested core never fabricates them.
  expect_true(length(oof_isna)  > 0L,
              info = "OOF synthesises _isNA companions (build_design_matrix_patches)")
  expect_identical(length(final_isna), 0L,
                   info = paste("DEFECT: FINAL synthesises NO _isNA companions",
                                "from the same upstream frame"))
  # The companion sets are NOT identical -> the deployed feature spaces differ.
  expect_false(setequal(oof_isna, final_isna),
               info = paste("DEFECT: OOF/FINAL _isNA companion SET differs when",
                            "fed the same extract_supervised_features output."))
  # The full design-matrix x_cols are NOT identical either (neither set nor order)
  expect_false(setequal(oof_cols, final_cols),
               info = "DEFECT: OOF/FINAL deployed x_cols differ (set).")
  expect_false(identical(oof_cols, final_cols),
               info = "DEFECT: OOF/FINAL deployed x_cols differ (order).")

  # Pin the EXACT divergence so the fix is unambiguous: every OOF `_isNA` is one
  # FINAL lacks; FINAL has no `_isNA` FINAL-only.
  expect_setequal(setdiff(final_isna, oof_isna), character(0))
  expect_true(length(setdiff(oof_isna, final_isna)) == length(oof_isna))
})

# ---------------------------------------------------------------------------
# (3) DETERMINISM of the boundary: re-running each stage's column derivation on
#     the SAME upstream frame is byte-identical (so the parity comparison above
#     is stable, and a future fix can be asserted exactly).
# ---------------------------------------------------------------------------
test_that("each stage's deployed columns are deterministic for a fixed upstream frame", {
  skip_if_not_installed("sf"); skip_if_not_installed("xgboost")
  skip_if_not_installed("Matrix")
  df <- make_isna_upstream_df()
  expect_identical(oof_model_cols_from(df), oof_model_cols_from(df))
  expect_identical(final_x_cols_from(df),   final_x_cols_from(df))
})
