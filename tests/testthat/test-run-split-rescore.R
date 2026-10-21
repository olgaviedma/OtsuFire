# run_split_rescore(): parameter validation, a synthetic end-to-end case (feature extraction and
# prediction stubbed) and, when the real 2006 data exist, parity with the supervised split flow.

mk_rescore_inputs <- function() {
  fx <- mk_doy_fixture()
  cand <- fx$cand
  cand$poly_id <- c(11L, 12L)
  cand$p_burned_model <- c(0.9, 0.2)
  dummy <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(nrows = 2, ncols = 2, vals = 1:4), dummy)
  list(fx = fx, cand = cand, dummy = dummy,
       recipe = list(cols = list(x_cols = "a"), training = list(seed = 1L), params = list(eta = 0.05)))
}

call_rescore <- function(inp, ..., out_dir = tempfile("rescore")) {
  args <- list(year = 2006L, scored_candidates = inp$cand, change_index = inp$fx$ci,
               delayed_change_index = inp$dummy, topo = inp$dummy, corine_raster = inp$dummy,
               model = list(stub = TRUE), recipe = inp$recipe, threshold = 0.5, out_dir = out_dir, verbose = FALSE)
  extra <- list(...); args[names(extra)] <- extra   # (modifyList would merge lists such as 'recipe')
  do.call(run_split_rescore, args)
}

test_that("run_split_rescore validates its parameters", {
  inp <- mk_rescore_inputs()
  expect_error(call_rescore(inp, threshold = 1.5), "threshold")
  expect_error(call_rescore(inp, threshold = NA_real_), "threshold")
  expect_error(call_rescore(inp, year = 2006.5), "year")
  expect_error(call_rescore(inp, doy_gap_days = -1), "doy_gap_days")
  expect_error(call_rescore(inp, min_group_share = 1), "min_group_share")
  expect_error(call_rescore(inp, window_margin_m = -5), "window_margin_m")
  expect_error(call_rescore(inp, artifact_hard_weight_ratio = "a"), "artifact_hard_weight_ratio")
  expect_error(call_rescore(inp, overwrite = NA), "overwrite")
  expect_error(call_rescore(inp, topo = "no_such_file.tif"), "topo")
  expect_error(call_rescore(inp, recipe = list(a = 1)), "recipe")
  expect_error(call_rescore(inp, scored_candidates = 5), "scored_candidates")
  bad <- inp$cand; bad$p_burned_model <- NULL
  expect_error(call_rescore(inp, scored_candidates = bad), "p_burned_model")
})

test_that("run_split_rescore splits, rescoring only the new sub-patches, and writes recipe and audit", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  inp <- mk_rescore_inputs()
  # stubs: every sub-patch gets p = 0.8 when its DOY median is late (>= 200), 0.1 otherwise
  stub_extract <- function(g, ...) { g$fire_uid <- g$sub_id; g$doy_med_feat <- g$doy_med; g }
  stub_predict <- function(df, recipe, model) ifelse(df$doy_med_feat >= 200, 0.8, 0.1)
  testthat::local_mocked_bindings(.of_rescore_extract = stub_extract, .of_rescore_predict = stub_predict)
  od <- tempfile("rescore")
  res <- call_rescore(inp, out_dir = od, artifact_hard_weight_ratio = 0.5)

  expect_equal(res$split$n_mixed, 1L)
  expect_equal(res$split$n_subpatches, 2L)
  m <- sf::st_read(res$final_map_path, quiet = TRUE)
  expect_equal(nrow(m), 3L)                                  # 1 homogeneous + 2 sub-patches
  expect_true(12L %in% m$poly_id)                            # homogeneous candidate untouched
  expect_equal(m$p_burned_model[m$poly_id == 12L], 0.2)      # its stored score is not recomputed
  expect_false(11L %in% m$poly_id)                           # split parent replaced
  kept <- sf::st_read(res$thresholded_path, quiet = TRUE)
  expect_equal(nrow(kept), 1L)                               # only the late sub-patch passes 0.5
  expect_equal(res$n_kept, 1L)
  a <- utils::read.csv(file.path(od, "run_split_rescore_audit.csv"))
  expect_equal(a$n_split_parents, 1L)
  expect_equal(a$n_kept_original, 1L)
  expect_equal(a$n_kept_split, 1L)
  rec <- readLines(file.path(od, "RECIPE.yml"))
  expect_true(any(grepl("doy_gap_days: 30", rec, fixed = TRUE)))
  expect_true(any(grepl("split_by_doy: true", rec, ignore.case = TRUE)))
  expect_true(any(grepl("artifact_hard_weight_ratio: 0.5", rec, fixed = TRUE)))
  expect_true(any(grepl("git_commit", rec, fixed = TRUE)))
  expect_true(file.exists(file.path(od, "sessionInfo.txt")))

  # no overwrite by default; allowed on request
  expect_error(call_rescore(inp, out_dir = od), "overwrite")
  expect_silent(call_rescore(inp, out_dir = od, overwrite = TRUE))
})

test_that("run_split_rescore with no mixed candidate leaves the map unchanged", {
  skip_if_not_installed("sf"); skip_if_not_installed("terra")
  inp <- mk_rescore_inputs()
  inp$cand <- inp$cand[2, ]                                  # only the homogeneous candidate
  testthat::local_mocked_bindings(.of_rescore_extract = function(...) stop("must not be called"),
                                  .of_rescore_predict = function(...) stop("must not be called"))
  res <- call_rescore(inp)
  expect_equal(res$split$n_mixed, 0L)
  expect_equal(nrow(sf::st_read(res$final_map_path, quiet = TRUE)), 1L)
  expect_true(is.na(res$subpatches_path))
})

test_that("recipe serialiser nests lists and vectors", {
  y <- .of_to_yaml(list(a = 1, b = list(c = "x", d = c(1, 2)), e = NULL))
  expect_true("a: 1" %in% y)
  expect_true("b:" %in% y)
  expect_true("  c: \"x\"" %in% y)
  expect_true("  d: [1, 2]" %in% y)
})

test_that("real 2006 case: same sub-patch scores as the supervised split flow (skipped without data)", {
  skip_on_cran()
  b <- "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/"
  d <- file.path(b, "1_DATA")
  mdir <- file.path(d, "Results/1989/SUPERVISED/FINAL_MAPS/POOL_GLOBAL_SUPPL/FINAL/1989_very_permissive_keepPool_2017_2022_2025_patch_certified_")
  cand <- list.files(file.path(d, "Results/2006/SUPERVISED/FINAL_MAPS/POOL_GLOBAL_SUPPL/SCORE"), "deterministic_scored[.]gpkg$", full.names = TRUE)
  ref <- file.path(b, "00_HANDOFF/agent_runs/2026-10-21/validacion_paquete/2006/subpatches_G30_pkg_scored.gpkg")
  paths <- c(file.path(d, "Imagery/Composites_90m/Min_Min/MinMin_2006_mosaic_res90m.tif"),
             file.path(d, "Imagery/Composites_90m/Autumn/mean_mean_2006_mosaic.tif"),
             file.path(d, "Topography/elevation_slope.tif"),
             file.path(d, "Corine_Masks/CLC_2006_peninsula.tif"), paste0(mdir, "final_model.rds"), paste0(mdir, "recipe.rds"), ref)
  skip_if_not(length(cand) == 1L && all(file.exists(paths)), "real 2006 data not available")
  s <- sf::st_read(cand, quiet = TRUE)
  small <- s[as.character(s$poly_id) %in% c("605", "1294"), ]
  res <- run_split_rescore(year = 2006L, scored_candidates = small, change_index = paths[1],
                           delayed_change_index = paths[2], topo = paths[3], corine_raster = paths[4],
                           model = paths[5], recipe = paths[6], threshold = 0.51,
                           out_dir = tempfile("real2006"), verbose = FALSE)
  new <- sf::st_read(res$subpatches_path, quiet = TRUE)
  old <- sf::st_read(ref, quiet = TRUE)
  old <- old[as.character(old$parent_id) %in% c("605", "1294"), ]
  expect_setequal(new$sub_id, old$sub_id)
  dp <- abs(new$p_burned[match(old$sub_id, new$sub_id)] - old$p_burned)
  # topography dispersion of very small sub-patches moves slightly with the cropped window
  expect_lt(stats::median(dp), 0.005)
  expect_gt(mean(dp < 0.01), 0.9)
  expect_lt(max(dp), 0.1)
})
