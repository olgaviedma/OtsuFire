# Bug 8 (OtsuFire 0.3.0): build_design_matrix_patches must NOT generate
# second-order `<x>_isNA_isNA` flag columns.

# --- T16 ------------------------------------------------------------
test_that("T16: prep$matrix_colnames does not contain *_isNA_isNA", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("xgboost")

  fn <- get("build_design_matrix_patches", envir = asNamespace("OtsuFire"))

  # Tiny synthetic labelled + burned-like frames.
  set.seed(1)
  n <- 8L
  L_df <- data.frame(
    fire_uid = sprintf("L%02d", seq_len(n)),
    poly_id  = sprintf("p%02d", seq_len(n)),
    class    = rep(c("burned", "unburned"), each = n / 2),
    block_id = rep(c(1L, 2L), each = n / 2),
    fold_rep1 = rep(c(1L, 2L), n / 2),
    rbr_med  = c(10, 12, 14, NA, 11, NA, 13, 15),
    elev_mean = c(100, 200, 150, 250, 110, 220, NA, 270),
    rbr_med_isNA = as.integer(c(0, 0, 0, 1, 0, 1, 0, 0)),
    stringsAsFactors = FALSE
  )
  BL_df <- data.frame(
    fire_uid = sprintf("BL%02d", seq_len(n)),
    poly_id  = sprintf("q%02d", seq_len(n)),
    class    = rep("unburned", n),
    block_id = NA_integer_,
    fold_rep1 = NA_integer_,
    rbr_med  = c(10, 12, 14, NA, 11, NA, 13, 15),
    elev_mean = c(100, 200, 150, 250, NA, 220, NA, 270),
    rbr_med_isNA = as.integer(c(0, 0, 0, 1, 0, 1, 0, 0)),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(fn(
    labelled = L_df,
    burned_like = BL_df,
    id_cols = c("fire_uid", "class", "poly_id", "block_id", "fold_rep1"),
    drop_regex = c("^fold_rep", "^block_id$", "^poly_id$"),
    cat_cols = character(0),
    hs_n_col = NA_character_,
    hs_conf_col = NA_character_,
    hs_frp_col = NA_character_,
    median_from = "labelled",
    verbose = FALSE
  ))
  cn <- out$prep$matrix_colnames
  if (is.null(cn)) cn <- colnames(out$XL_mat)
  expect_false(any(grepl("_isNA_isNA$", cn)),
               info = paste("found _isNA_isNA columns:",
                            paste(cn[grepl("_isNA_isNA$", cn)], collapse = ", ")))
})
