make_rect_eco <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(matrix(
    c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ),
    ncol = 2,
    byrow = TRUE
  )))
}

make_extract_features_fixture <- function() {
  crs_epsg <- 3035

  train_folds <- sf::st_sf(
    fire_uid = c("t1", "t2", "t3"),
    class = c("burned", "unburned", "burned"),
    block_id = c(1L, 1L, 2L),
    fold_rep1 = c(1L, 1L, 2L),
    geometry = sf::st_sfc(
      make_rect_eco(0, 0, 1, 1),
      make_rect_eco(3, 0, 4, 1),
      make_rect_eco(1, 3, 2, 4),
      crs = crs_epsg
    )
  )

  unlabeled <- sf::st_sf(
    fire_uid = c("u1", "u2", "u3"),
    geometry = sf::st_sfc(
      make_rect_eco(0, 1, 1, 2),
      make_rect_eco(3, 1, 4, 2),
      make_rect_eco(4, 3, 5, 4),
      crs = crs_epsg
    )
  )

  template <- terra::rast(
    nrows = 5, ncols = 5,
    xmin = 0, xmax = 5, ymin = 0, ymax = 5,
    crs = "EPSG:3035"
  )

  rbr_summer <- template
  terra::values(rbr_summer) <- seq_len(terra::ncell(rbr_summer))

  dem <- template
  terra::values(dem) <- seq(100, 124)

  slope <- template
  terra::values(slope) <- seq(5, 29)

  corine_r <- template
  terra::values(corine_r) <- rep(c(1, 1, 1, 2, 2), times = 5)

  ecoregions <- sf::st_sf(
    eco_id = c("west", "east"),
    geometry = sf::st_sfc(
      make_rect_eco(0, 0, 2.5, 5),
      make_rect_eco(2.5, 0, 5, 5),
      crs = crs_epsg
    )
  )

  list(
    train_folds = train_folds,
    unlabeled = unlabeled,
    rbr_summer = rbr_summer,
    dem = dem,
    slope = slope,
    corine_r = corine_r,
    cor_groups = list(natural = c(1), developed = c(2)),
    ecoregions = ecoregions
  )
}

run_extract_features_fixture <- function(use_ecoregions, ecoregions = NULL) {
  fx <- make_extract_features_fixture()
  extract_features_internal <- get("extract_features", envir = asNamespace("OtsuFire"))

  extract_features_internal(
    train_folds = fx$train_folds,
    unlabeled = fx$unlabeled,
    build_features = TRUE,
    rbr_summer = fx$rbr_summer,
    doy_post = NULL,
    rbr_aw = NULL,
    nbr_pre = NULL,
    nbr_post = NULL,
    dnbr = NULL,
    dem = fx$dem,
    slope = fx$slope,
    corine_r = fx$corine_r,
    ecoregions = ecoregions,
    cor_groups = fx$cor_groups,
    use_doy = FALSE,
    use_aw = FALSE,
    use_nbr = FALSE,
    use_hotspots = FALSE,
    use_ecoregions = use_ecoregions,
    return_features = TRUE,
    return_features_geometry = FALSE,
    save_features_gpkg = FALSE,
    verbose = FALSE
  )
}

test_that("extract_features skips eco_major when use_ecoregions is FALSE", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  out <- run_extract_features_fixture(
    use_ecoregions = FALSE,
    ecoregions = NULL
  )

  expect_false("eco_major" %in% names(out$train_feat))
  expect_false("eco_major" %in% names(out$unl_feat))
})

test_that("extract_features adds eco_major when use_ecoregions is TRUE", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  fx <- make_extract_features_fixture()
  out <- run_extract_features_fixture(
    use_ecoregions = TRUE,
    ecoregions = fx$ecoregions
  )

  expect_true("eco_major" %in% names(out$train_feat))
  expect_equal(
    as.character(out$train_feat$eco_major),
    c("west", "east", "west")
  )
  expect_equal(
    as.character(out$unl_feat$eco_major),
    c("west", "east", "east")
  )
})

test_that("extract_features errors when use_ecoregions is TRUE but ecoregions is NULL", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  expect_error(
    run_extract_features_fixture(
      use_ecoregions = TRUE,
      ecoregions = NULL
    ),
    regexp = "use_ecoregions=TRUE requires ecoregions"
  )
})
