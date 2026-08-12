# Phase B (block 3): 15 implemented, 0 stubbed.
# §N+27 (2026-06-05): update_burned_like_registry removed (abandoned
# supervised burned-like registry research line). Count dropped 16 -> 15.

implemented <- c(
  "change_index_mosaic", "validate_fire_maps",
  "build_burned_mapping_config", "detect_burned_patches",
  "score_burned_patches", "run_deterministic_pipeline",
  "build_supervised_burned_config", "build_supervised_training_pools",
  "make_spatial_folds", "extract_supervised_features",
  "run_oof_diagnostics", "train_final_burned_model",
  "score_supervised_burned_map",
  "check_supervised_consistency", "run_oneyear_supervised_pipeline"
)

test_that("Phase B-3: all 15 public names exist in the namespace", {
  for (nm in implemented) {
    expect_true(
      exists(nm, where = asNamespace("OtsuFire"), inherits = FALSE),
      info = paste("missing:", nm)
    )
    expect_true(
      is.function(get(nm, envir = asNamespace("OtsuFire"), inherits = FALSE)),
      info = paste("not a function:", nm)
    )
  }
})

test_that("Phase B-3: implemented functions do not return NotYetImplemented", {
  for (nm in implemented) {
    err <- tryCatch(do.call(nm, list()), error = function(e) conditionMessage(e))
    expect_false(
      grepl("not implemented yet", err, ignore.case = TRUE),
      info = paste("stub still active for implemented function:", nm)
    )
  }
})

test_that("validate_burned_maps is intentionally NOT a public export", {
  ns_path <- file.path(
    "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild",
    "NAMESPACE"
  )
  if (file.exists(ns_path)) {
    ns_txt <- readLines(ns_path, warn = FALSE)
    expect_false(any(grepl("^\\s*export\\(\\s*validate_burned_maps\\s*\\)\\s*$",
                           ns_txt)),
                 info = "validate_burned_maps must not be reintroduced in 0.2.x")
  }
})

test_that("run_supervised_pipeline is intentionally NOT a public export", {
  # The legacy multi-year orchestrator was migrated as an internal helper in
  # Block 4B (lives in the package namespace, not exported).
  ns_path <- file.path(
    "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild",
    "NAMESPACE"
  )
  if (file.exists(ns_path)) {
    ns_txt <- readLines(ns_path, warn = FALSE)
    expect_false(any(grepl("^\\s*export\\(\\s*run_supervised_pipeline\\s*\\)\\s*$",
                           ns_txt)),
                 info = "run_supervised_pipeline must not be a public export")
  }
})

test_that("%||% is internal and not a 0.2.x public export", {
  ns <- asNamespace("OtsuFire")
  expect_true(exists("%||%", envir = ns, inherits = FALSE))
  ns_path <- file.path(
    "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/2_SCRIPTS/OtsuFire_v02_rebuild",
    "NAMESPACE"
  )
  if (file.exists(ns_path)) {
    ns_txt <- readLines(ns_path, warn = FALSE)
    expect_false(any(grepl("^\\s*export\\(\\s*[\"`]?%\\|\\|%[\"`]?\\s*\\)\\s*$",
                           ns_txt)))
  }
})
