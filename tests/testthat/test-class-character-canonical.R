# Bug 5 (OtsuFire 0.3.0): `class` is character end-to-end.
#
# The orchestrator no longer coerces `class` to factor at A4 / A5 / C0d.
# A GPKG round-trip + downstream join should therefore see character
# on both sides.

# --- T14 ------------------------------------------------------------
test_that("T14: GPKG round-trip preserves character `class`", {
  skip_if_not_installed("sf")

  sfc <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))),
    sf::st_polygon(list(rbind(c(2,0), c(3,0), c(3,1), c(2,1), c(2,0)))),
    crs = 3035
  )
  x <- sf::st_sf(
    fire_uid = c("a", "b"),
    class = c("burned", "unburned"),
    geometry = sfc
  )
  expect_type(x$class, "character")

  f <- tempfile(fileext = ".gpkg")
  sf::st_write(x, f, layer = "test_layer", quiet = TRUE,
               delete_dsn = TRUE)
  y <- sf::read_sf(f, layer = "test_layer")

  expect_type(y$class, "character")

  # Simulated downstream join: both sides character → straight match.
  joined <- merge(
    sf::st_drop_geometry(y),
    data.frame(fire_uid = c("a", "b"),
               class = c("burned", "unburned"),
               extra = c("X", "Y"),
               stringsAsFactors = FALSE),
    by = c("fire_uid", "class")
  )
  expect_equal(nrow(joined), 2L)
  expect_setequal(joined$extra, c("X", "Y"))
})

test_that("orchestrator no longer coerces character columns to factor", {
  # The pre-0.3.0 orchestrator had multiple
  # `mutate(across(where(is.character), as.factor))` calls. None must
  # remain. This guard reads the package source on disk, which only
  # exists during `testthat::test_dir()` from the source tree; under
  # R CMD check the installed copy has no R/ directory, so skip there.
  testthat::skip_on_cran()
  src_path <- file.path(testthat::test_path("..", ".."),
                        "R", "internal-sup-orchestrator.R")
  testthat::skip_if_not(file.exists(src_path),
    "source file not available outside the package source tree")
  src <- readLines(src_path)
  hits <- grep("mutate\\(across\\(where\\(is\\.character\\), as\\.factor\\)\\)",
               src)
  expect_length(hits, 0L)
})
