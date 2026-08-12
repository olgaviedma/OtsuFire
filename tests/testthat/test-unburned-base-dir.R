# config$options$unburned_base_dir (0.10.2): the negative pool only depends on
# the YEAR, so a caller may pin the random/deterministic unburned GPKG and the
# Otsu-negative working root to a shared, scenario-independent base. When the
# option is absent the historical scenario-scoped paths are returned unchanged.

ub_resolve <- get(".of_resolve_unburned_paths", envir = asNamespace("OtsuFire"))

test_that("default (no unburned_base_dir): negative paths are scenario-scoped", {
  cfg <- list(options = list(data_base = "D:/base"))
  p <- ub_resolve(cfg, 1985L, "Min_Min", "balanced", "R:/run")
  expect_equal(
    p$unb_out_gpkg,
    file.path("D:/base", "Results", 1985L, "Min_Min", "DETERMINISTIC",
              "balanced", "UNBURNED", "1985_balanced_unburned.gpkg")
  )
  expect_equal(p$otsu_negative_root_dir, file.path("R:/run", "_OTSU_NEGATIVE"))
})

test_that("unburned_base_dir set: year-scoped paths, shared across scenarios", {
  base <- "D:/base/Results/1985/Min_Min/DETERMINISTIC/_UNBURNED"
  cfg <- list(options = list(data_base = "D:/base", unburned_base_dir = base))

  p_bal <- ub_resolve(cfg, 1985L, "Min_Min", "balanced", "R:/run_bal")
  p_kp  <- ub_resolve(cfg, 1985L, "Min_Min",
                      "balanced_keepPool_2017_2022_2025", "R:/run_kp")

  expect_equal(p_bal$unb_out_gpkg,
               file.path(base, "UNBURNED", "1985_unburned.gpkg"))
  expect_equal(p_bal$otsu_negative_root_dir,
               file.path(base, "_OTSU_NEGATIVE"))

  # Different scenarios / result_dirs resolve to the SAME shared negative paths.
  expect_identical(p_bal, p_kp)

  # The shared paths carry NO scenario string.
  expect_false(grepl("balanced", p_bal$unb_out_gpkg, fixed = TRUE))
  expect_false(grepl("keepPool", p_bal$unb_out_gpkg, fixed = TRUE))
  expect_false(grepl("balanced", p_bal$otsu_negative_root_dir, fixed = TRUE))
})

test_that("empty or NULL unburned_base_dir falls back to the default", {
  default <- ub_resolve(list(options = list(data_base = "D:/base")),
                        1985L, "Min_Min", "balanced", "R:/run")
  cfg_empty <- list(options = list(data_base = "D:/base", unburned_base_dir = ""))
  cfg_null  <- list(options = list(data_base = "D:/base", unburned_base_dir = NULL))
  expect_identical(ub_resolve(cfg_empty, 1985L, "Min_Min", "balanced", "R:/run"),
                   default)
  expect_identical(ub_resolve(cfg_null, 1985L, "Min_Min", "balanced", "R:/run"),
                   default)
})
