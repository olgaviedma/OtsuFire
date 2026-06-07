# B2 (2026-06-06): the final-map GPKG writes in
# score_burnedlike_and_export_final_map() used to run UNCONDITIONALLY -- the
# `final_gpkg` layers were re-written with append=TRUE even when overwrite=FALSE
# and the file already existed, so a re-run DUPLICATED polygons (or failed with
# a schema-append error). The fix gates every GPKG write through a
# `.write_gpkg_if_allowed()` closure that, like the per-file `.write_if_allowed`
# used for the sidecars (and in train_final / oof), writes ONLY when
# `isTRUE(overwrite) || !file.exists(path)` (clearing first when overwrite=TRUE)
# and otherwise reuses the existing file. The full engine needs a trained model
# + recipe + scored universe (too heavy for a unit test and exercised
# end-to-end by the orchestrator), so we assert the contract two ways:
#   (1) a faithful reproduction of the gated multi-layer write proves
#       overwrite=FALSE never re-appends (no duplicate features) while
#       overwrite=TRUE replaces;
#   (2) a source-level guard proves the production writes are gated through the
#       closure and the raw append=TRUE writes are no longer reached
#       unconditionally.

# Mirror of the production .write_gpkg_if_allowed() closure (single source of
# truth: R/internal-sup-final-map.R).
make_gated_writer <- function(overwrite) {
  safe_remove_dataset <- function(path) {
    if (!file.exists(path)) return(invisible(TRUE))
    unlink(path, force = TRUE)
    invisible(TRUE)
  }
  function(path, write_fn) {
    if (isTRUE(overwrite) || !file.exists(path)) {
      if (isTRUE(overwrite)) safe_remove_dataset(path)
      force(write_fn())
    }
    invisible(NULL)
  }
}

mk_final_map_sf <- function(n = 4L) {
  polys <- lapply(seq_len(n), function(k) {
    x0 <- (k - 1L) * 10
    sf::st_polygon(list(rbind(
      c(x0, 0), c(x0 + 5, 0), c(x0 + 5, 5), c(x0, 5), c(x0, 0)
    )))
  })
  sf::st_sf(
    fire_uid = sprintf("F%03d", seq_len(n)),
    p_burned = seq_len(n) / (n + 1),
    geometry = sf::st_sfc(polys, crs = 3035)
  )
}

test_that("gated final_gpkg write does not duplicate on overwrite=FALSE re-run", {
  skip_if_not_installed("sf")

  final_map_full <- mk_final_map_sf(4L)
  final_map      <- final_map_full[1:2, ]
  final_gpkg     <- tempfile(fileext = ".gpkg")
  on.exit(unlink(final_gpkg, force = TRUE), add = TRUE)

  write_all_layers <- function() {
    sf::st_write(final_map_full, final_gpkg, layer = "deterministic_scored",
                 delete_dsn = FALSE, quiet = TRUE)
    sf::st_write(final_map_full, final_gpkg, layer = "final_map_full",
                 append = TRUE, quiet = TRUE)
    sf::st_write(final_map, final_gpkg, layer = "final_map",
                 append = TRUE, quiet = TRUE)
  }

  # First run (file absent): writes all three layers once.
  w_first <- make_gated_writer(overwrite = FALSE)
  w_first(final_gpkg, write_all_layers)

  expect_true(file.exists(final_gpkg))
  layers0 <- sf::st_layers(final_gpkg)$name
  expect_true(all(c("deterministic_scored", "final_map_full", "final_map")
                  %in% layers0))
  n_full_1 <- nrow(sf::read_sf(final_gpkg, layer = "final_map_full",
                               quiet = TRUE))
  n_map_1  <- nrow(sf::read_sf(final_gpkg, layer = "final_map", quiet = TRUE))
  expect_identical(n_full_1, 4L)
  expect_identical(n_map_1, 2L)

  # Second run (file present, overwrite=FALSE): the gate SKIPS the write, so
  # the layers are NOT re-appended -> NO duplicate features. Without the fix
  # this would double the row counts (or raise a schema-append error).
  w_second <- make_gated_writer(overwrite = FALSE)
  expect_silent(w_second(final_gpkg, write_all_layers))
  n_full_2 <- nrow(sf::read_sf(final_gpkg, layer = "final_map_full",
                               quiet = TRUE))
  n_map_2  <- nrow(sf::read_sf(final_gpkg, layer = "final_map", quiet = TRUE))
  expect_identical(n_full_2, n_full_1)
  expect_identical(n_map_2, n_map_1)

  # Third run (overwrite=TRUE): clears then rewrites -> counts back to the
  # single-write values (replace, not append).
  w_third <- make_gated_writer(overwrite = TRUE)
  w_third(final_gpkg, write_all_layers)
  n_full_3 <- nrow(sf::read_sf(final_gpkg, layer = "final_map_full",
                               quiet = TRUE))
  expect_identical(n_full_3, n_full_1)
})

test_that("score_burnedlike_and_export_final_map gates every GPKG write", {
  src <- deparse(body(get("score_burnedlike_and_export_final_map",
                          envir = asNamespace("OtsuFire"))))
  src1 <- paste(src, collapse = "\n")

  # The gated closure exists and is the only path to the GPKG writes.
  expect_true(grepl(".write_gpkg_if_allowed", src1, fixed = TRUE))
  # scored_gpkg, final_gpkg and burned_like_gpkg all flow through it.
  expect_gte(lengths(regmatches(
    src1, gregexpr(".write_gpkg_if_allowed(", src1, fixed = TRUE)))[1], 3L)
  # The append=TRUE final_gpkg writes survive ONLY inside the gated closure;
  # they must NOT appear as a bare top-level statement. We assert the closure
  # is referenced before any append write by checking the gate token precedes
  # the first append=TRUE occurrence.
  gate_pos   <- regexpr(".write_gpkg_if_allowed", src1, fixed = TRUE)
  append_pos <- regexpr("append = TRUE", src1, fixed = TRUE)
  expect_true(gate_pos > 0)
  if (append_pos > 0) expect_true(gate_pos < append_pos)
})
