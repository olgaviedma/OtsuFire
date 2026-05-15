# Tests for build_burned_mapping_config(). Lightweight — no engine runs.

mk_tmp_tif <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = 1:64)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mk_tmp_mask <- function() {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(ncol = 8, nrow = 8, vals = rep(c(0L, 1L), 32))
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("config builder returns the expected S3 class and required fields", {
  ci <- mk_tmp_tif()
  bm <- mk_tmp_mask()

  cfg <- build_burned_mapping_config(
    change_index = ci,
    burnable_mask = bm,
    target_year = 2025,
    output_dir = tempdir()
  )

  expect_s3_class(cfg, "otsufire_burned_mapping_config")
  expect_false("scenario" %in% names(cfg))
  expect_identical(cfg$target_year, 2025L)
  expect_identical(cfg$run_name, "deterministic_burned_map")
  expect_type(cfg$inputs, "list")
  expect_named(cfg$inputs,
               c("change_index", "vegetation_map", "burnable_mask",
                 "hotspots", "previous_year_burned", "reference_burned_map"))
  expect_identical(cfg$deterministic_seed, 12345L)
  expect_type(cfg$output_routes, "list")
  expect_true(all(c("grow_vector", "refine_merged_gpkg",
                    "internal_decisions", "timing_csv") %in%
                  names(cfg$output_routes)))
})

test_that("scenario has been fully removed from the public surface", {
  expect_false("scenario" %in% names(formals(build_burned_mapping_config)))

  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    output_dir = tempdir())
  expect_false("scenario" %in% names(cfg))

  # Passing scenario= is rejected as an unused argument.
  expect_error(
    build_burned_mapping_config(
      scenario = "anything", change_index = ci, burnable_mask = bm,
      target_year = 2025L),
    regexp = "unused argument"
  )
})

test_that("run_name is the single visible experiment label and drives routes", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()

  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2012L,
    run_name = "my_custom_run_2025", output_dir = tempdir()
  )
  expect_identical(cfg$run_name, "my_custom_run_2025")
  # Single run_name segment: <year>/DETERMINISTIC/<run_name> (no duplicate).
  expect_match(cfg$output_routes$base,
               "2012/DETERMINISTIC/my_custom_run_2025",
               fixed = TRUE)
  expect_false(grepl("my_custom_run_2025/DETERMINISTIC/my_custom_run_2025",
                      cfg$output_routes$base, fixed = TRUE))
  expect_match(cfg$output_routes$validation_workbook,
               "my_custom_run_2025_res30.xlsx", fixed = TRUE)
})

test_that("registry path uses toupper(run_name)", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    run_name = "expA", output_dir = tempdir()
  )
  expect_match(cfg$registry_path,
               "/EXPA/RBR_TRAINING_REGISTRY\\.gpkg$")
})

test_that("empty run_name fails", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      run_name = ""),
    regexp = "run_name"
  )
})

test_that("run_name containing / or backslash fails", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      run_name = "a/b"),
    regexp = "run_name"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      run_name = "a\\b"),
    regexp = "run_name"
  )
})

test_that("config builder rejects missing target_year and invalid year values", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(change_index = ci, burnable_mask = bm),
    regexp = "target_year"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 1500L
    ),
    regexp = "target_year"
  )
})

test_that("config builder rejects missing burnable_mask", {
  ci <- mk_tmp_tif()
  expect_error(
    build_burned_mapping_config(change_index = ci, target_year = 2025L),
    regexp = "burnable_mask"
  )
})

test_that("config builder rejects unnamed options list entries", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      options = list(1, 2, 3)
    ),
    regexp = "named list"
  )
})

test_that("options$change_index_validation block is structurally validated", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()

  # Accepted: a named list using only the recognized keys.
  expect_s3_class(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      options = list(change_index_validation = list(lower_cap = -250))),
    "otsufire_burned_mapping_config")

  # Rejected: not a list.
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      options = list(change_index_validation = "nope")),
    regexp = "change_index_validation"
  )
  # Rejected: unknown sub-key.
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      options = list(change_index_validation = list(bogus = 1))),
    regexp = "Unknown options\\$change_index_validation"
  )
})

test_that("detect/refine/scoring params must be named lists", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      detect_params = list(310, 90)),
    regexp = "named list"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      refine_params = list(5000)),
    regexp = "named list"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      scoring_params = list(90)),
    regexp = "named list"
  )
})

test_that("unknown keys in the public param blocks fail", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      detect_params = list(seed_threshold = 300, nope = 1)),
    regexp = "Unknown detect_params"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      refine_params = list(aoi_buffer_m = 1000, bogus = 1)),
    regexp = "Unknown refine_params"
  )
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      scoring_params = list(reference_buffer_m = 30, ghost = 1)),
    regexp = "Unknown scoring_params"
  )
})

test_that("resolved blocks are present in the returned config (defaults)", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L)

  expect_type(cfg$detect_params, "list")
  expect_type(cfg$refine_params, "list")
  expect_type(cfg$scoring_params, "list")
  expect_type(cfg$rescue_params, "list")

  # detect defaults (internal-name shape, == legacy balanced numbers).
  expect_identical(cfg$detect_params$otsu_thresholds, 310)
  expect_identical(cfg$detect_params$grow_delta, 90)
  expect_identical(cfg$detect_params$min_grow_threshold_value, 240)
  expect_identical(cfg$detect_params$min_seed_pixels_per_component, 30L)
  expect_named(cfg$detect_params$otsu_min_by_class, as.character(1:11))
  expect_identical(unname(cfg$detect_params$otsu_min_by_class["2"]), 420)
  expect_identical(unname(cfg$detect_params$grow_delta_by_class["4"]), 100)
  expect_identical(
    unname(cfg$detect_params$min_grow_threshold_by_class["10"]), 560)

  # refine / scoring defaults (internal-name shape).
  expect_identical(cfg$refine_params$aoi_buffer_m, 5000)
  expect_identical(cfg$refine_params$min_detected_area_m2, 10)
  expect_true(cfg$refine_params$merge_overlaps)
  expect_identical(cfg$scoring_params$support_buffer_m, 0)
  expect_identical(cfg$scoring_params$erase_mask_buffer_m, 90)
  expect_identical(cfg$scoring_params$erase_min_area_m2, 20000)
  expect_identical(cfg$scoring_params$ref_buffer_m, 90)

  # rescue derived from detect by-class via the legacy relaxation logic:
  #   class 1 (relax): seed -20, delta +10, floor -20
  expect_identical(unname(cfg$rescue_params$otsu_min_by_class["1"]), 290)
  expect_identical(unname(cfg$rescue_params$grow_delta_by_class["1"]), 90)
  expect_identical(
    unname(cfg$rescue_params$min_grow_threshold_by_class["1"]), 210)
  #   class 3: seed -10, delta +5, floor -10
  expect_identical(unname(cfg$rescue_params$otsu_min_by_class["3"]), 370)
  #   class 2 (unchanged)
  expect_identical(unname(cfg$rescue_params$otsu_min_by_class["2"]), 420)
})

test_that("global detect values broadcast to all 11 vegetation classes", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    detect_params = list(seed_threshold = 333, growth_delta = 44,
                         minimum_growth_threshold = 222))

  expect_identical(cfg$detect_params$otsu_thresholds, 333)
  expect_true(all(cfg$detect_params$otsu_min_by_class == 333))
  expect_true(all(cfg$detect_params$grow_delta_by_class == 44))
  expect_true(all(cfg$detect_params$min_grow_threshold_by_class == 222))
  expect_named(cfg$detect_params$otsu_min_by_class, as.character(1:11))
})

test_that("partial by-vegetation vectors are completed correctly", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()

  # Case A: partial vector + matching global -> missing classes use global.
  cfgA <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    detect_params = list(
      seed_threshold = 500,
      seed_threshold_by_vegetation = c("1" = 100, "5" = 200)))
  v <- cfgA$detect_params$otsu_min_by_class
  expect_identical(unname(v["1"]), 100)
  expect_identical(unname(v["5"]), 200)
  expect_identical(unname(v["2"]), 500)   # global broadcast for missing
  expect_identical(unname(v["11"]), 500)

  # Case B: partial vector, no global -> missing classes use legacy default.
  cfgB <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm, target_year = 2025L,
    detect_params = list(
      growth_delta_by_vegetation = c("4" = 999)))
  d <- cfgB$detect_params$grow_delta_by_class
  expect_identical(unname(d["4"]), 999)
  expect_identical(unname(d["2"]), 20)    # legacy default for class 2
  expect_identical(unname(d["10"]), 5)    # legacy default for class 10

  # Unnamed / non-numeric by-vegetation vector is rejected.
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      detect_params = list(seed_threshold_by_vegetation = c(100, 200))),
    regexp = "named numeric vector"
  )
  # Unknown vegetation class name is rejected.
  expect_error(
    build_burned_mapping_config(
      change_index = ci, burnable_mask = bm, target_year = 2025L,
      detect_params = list(seed_threshold_by_vegetation = c("99" = 1))),
    regexp = "unknown vegetation class"
  )
})

test_that("config builder accepts NULL change_index at build time", {
  bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = NULL, burnable_mask = bm, target_year = 2025L
  )
  expect_s3_class(cfg, "otsufire_burned_mapping_config")
  expect_null(cfg$inputs$change_index)
})

test_that("output routes respect run_name and year", {
  ci <- mk_tmp_tif(); bm <- mk_tmp_mask()
  cfg <- build_burned_mapping_config(
    change_index = ci, burnable_mask = bm,
    target_year = 1985L, run_name = "smoke_x",
    output_dir = tempdir()
  )
  expect_match(cfg$output_routes$base,
               "1985/DETERMINISTIC/smoke_x", fixed = FALSE)
  expect_match(cfg$output_routes$internal_decisions,
               "internal_decisions\\.gpkg$")
})
