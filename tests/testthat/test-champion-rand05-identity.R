# Integration test (Phase 2B): row-for-row identity of SUBR05 against the real
# champion ERA pools on disk. This is the regression guard for the SUBR05 seed
# (canonical 42L). It is SKIPPED whenever the real GeoPackages are not present,
# so the unit-test suite never depends on local machine paths.
#
# Reference artifacts (provenance: _ERA_POOLS):
#   premodis_training_pool_CHAMPION.gpkg        -> pre-SUBR05 (random = 51693)
#   premodis_training_pool_CHAMPION_RAND05.gpkg -> SUBR05      (random = 714)
# Override the directory holding them with env var OF_ERA_POOLS_DIR if needed.

era_dir <- Sys.getenv(
  "OF_ERA_POOLS_DIR",
  "C:/00_NATALIA_DOCTORADO/00_FIRE_MAPPING/1_DATA/Results/_ERA_POOLS")
champ_pool  <- file.path(era_dir, "premodis_training_pool_CHAMPION.gpkg")
rand05_pool <- file.path(era_dir, "premodis_training_pool_CHAMPION_RAND05.gpkg")

test_that("SUBR05 seed=42L reproduces the champion RAND05 pool row-for-row", {
  skip_if_not(file.exists(champ_pool),  "champion CHAMPION pool not on disk")
  skip_if_not(file.exists(rand05_pool), "champion RAND05 pool not on disk")
  skip_if_not_installed("sf")

  ch <- sf::st_read(champ_pool,  layer = "train_features", quiet = TRUE)
  r5 <- sf::st_drop_geometry(sf::st_read(rand05_pool, layer = "train_features", quiet = TRUE))

  ref_rand <- sort(as.character(r5$fire_uid[r5$source == "random_burnable_background"]))
  expect_equal(length(ref_rand), 714L)          # documented champion count

  ss  <- subsample_random_background(ch, random_to_burned_ratio = 0.5, seed = 42L)
  ssd <- sf::st_drop_geometry(ss)
  got_rand <- sort(as.character(ssd$fire_uid[ssd$source == "random_burnable_background"]))

  # (1) random subset is EXACTLY the champion's 714 (per-row fire_uid identity)
  expect_identical(got_rand, ref_rand)

  # (2) whole-pool identity (all sources), so nothing else was perturbed
  expect_identical(sort(as.character(ssd$fire_uid)), sort(as.character(r5$fire_uid)))
  expect_equal(nrow(ss), nrow(r5))

  # (3) positives / Otsu residual / artifact_hard are untouched vs CHAMPION
  chd <- sf::st_drop_geometry(ch)
  for (src in c("internal_keep_qc", "otsu_patch_residual", "artifact_hard")) {
    expect_equal(sum(ssd$source == src), sum(chd$source == src),
                 info = paste("bucket untouched:", src))
  }
  expect_equal(sum(ssd$class == "burned"), sum(chd$class == "burned"))
})

test_that("SUBR05 default seed (no seed arg) also reproduces champion RAND05", {
  skip_if_not(file.exists(champ_pool),  "champion CHAMPION pool not on disk")
  skip_if_not(file.exists(rand05_pool), "champion RAND05 pool not on disk")
  skip_if_not_installed("sf")

  ch <- sf::st_read(champ_pool,  layer = "train_features", quiet = TRUE)
  r5 <- sf::st_drop_geometry(sf::st_read(rand05_pool, layer = "train_features", quiet = TRUE))
  ref_rand <- sort(as.character(r5$fire_uid[r5$source == "random_burnable_background"]))

  ssd <- sf::st_drop_geometry(subsample_random_background(ch, random_to_burned_ratio = 0.5))
  got_rand <- sort(as.character(ssd$fire_uid[ssd$source == "random_burnable_background"]))
  expect_identical(got_rand, ref_rand)
})
