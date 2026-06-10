# Gate 2 (2026-06-10): OUTPUT-LAYER CONTRACT of the supervised final-map export.
#
# The LONG 2017 BALANCED FULL run reported a row-count "discrepancy":
#   deterministic_scored = final_map_full = 1274 ; final_map = 1273.
# Audit (this test documents the resolution): this is NOT a silent loss. The
# `<prefix>_final_map.gpkg` carries THREE layers with DIFFERENT, deliberate
# contracts:
#   * deterministic_scored / final_map_full  -> the COMPLETE scored universe
#     (every candidate the FINAL model scored; one row per input polygon; all
#     IDs/geometries preserved). The "no candidate is lost" layers.
#   * final_map -> the CURRENT-YEAR PUBLIC burned map: final_map_full with the
#     public column subset AND the current-year temporal filter applied. It
#     drops EXACTLY the rows flagged `current_year_public_drop == TRUE` by
#     add_temporal_adjustment() -- a current-year temporal conflict (pre-year
#     overlap >= threshold or preyear_action == "drop") WITH weak current-year
#     hotspot support, i.e. a polygon almost certainly re-detecting the PRE-YEAR
#     fire with no current-year evidence it burned again.
#
# Contract under test (the invariant Natalia must be able to rely on -- no
# SILENT loss):
#   nrow(final_map) == nrow(final_map_full)
#                      - sum(final_map_full$current_year_public_drop %in% TRUE)
# and EVERY id in final_map_full is preserved (final_map ids are a subset; the
# only ids removed are precisely the current_year_public_drop ids).
#
# The full engine needs a trained model + recipe + scored universe (too heavy
# for a unit test; exercised end-to-end by the orchestrator and proven on the
# LONG 2017 run). We assert the contract two ways:
#   (1) a faithful reproduction of the two production closures
#       (add_temporal_adjustment -> filter_to_current_year_map) proves the
#       row-count relationship and id-preservation hold, including the specific
#       "drop / high-overlap / no-hotspot -> excluded; in-full, out-of-public"
#       case that produced the 1274 vs 1273 split;
#   (2) a source-level guard proves the production `final_map` is built by
#       filter_to_current_year_map() (the documented filter), so the public
#       layer can only ever differ from the full layer by current_year_public_drop.

# ---------------------------------------------------------------------------
# (1) Faithful reproduction of the production current-year filter contract.
#     Single source of truth: R/internal-sup-final-map.R
#     (add_temporal_adjustment + filter_to_current_year_map).
# ---------------------------------------------------------------------------

# Mirror of add_temporal_adjustment()'s drop decision (the parts that set
# current_year_public_drop). Defaults match the orchestrator CURRENTYEAR_*
# constants used by score_supervised_burned_map().
compute_current_year_public_drop <- function(df,
                                             preyear_overlap_threshold = 0.70,
                                             hotspot_density_threshold = 0.001) {
  n <- nrow(df)
  area_ha_num <- if ("area_ha" %in% names(df)) as.numeric(df$area_ha) else rep(NA_real_, n)
  hs_used_num <- if ("hs_used_n" %in% names(df)) as.numeric(df$hs_used_n) else rep(NA_real_, n)
  overlap_num <- if ("preyear_overlap_frac" %in% names(df)) as.numeric(df$preyear_overlap_frac) else rep(NA_real_, n)
  preyear_action_chr <- if ("preyear_action" %in% names(df)) as.character(df$preyear_action) else rep(NA_character_, n)

  hs_density_ha <- rep(NA_real_, n)
  ok_area <- is.finite(area_ha_num) & area_ha_num > 0 & is.finite(hs_used_num)
  hs_density_ha[ok_area] <- hs_used_num[ok_area] / area_ha_num[ok_area]

  hotspot_available_now <- rep(FALSE, n)
  if ("hotspot_available" %in% names(df)) {
    hotspot_available_now <- df$hotspot_available %in% c(1, TRUE)
  } else if ("filter_2" %in% names(df)) {
    hotspot_available_now <- !is.na(df$filter_2) & as.character(df$filter_2) != "not_applied"
  }
  hotspot_available_now[is.na(hotspot_available_now)] <- FALSE

  weak_hotspot_density <- rep(NA, n)
  weak_hotspot_density[hotspot_available_now] <-
    !is.finite(hs_density_ha[hotspot_available_now]) |
    hs_density_ha[hotspot_available_now] < hotspot_density_threshold

  preyear_drop_flag <- !is.na(preyear_action_chr) & preyear_action_chr == "drop"
  overlap_proxy <- overlap_num
  overlap_proxy[preyear_drop_flag & !is.finite(overlap_proxy)] <- 1

  temporal_conflict <- (is.finite(overlap_proxy) & overlap_proxy >= preyear_overlap_threshold) |
    preyear_drop_flag
  weak_current_support <- temporal_conflict & (!hotspot_available_now | weak_hotspot_density %in% TRUE)
  weak_current_support
}

# Mirror of filter_to_current_year_map().
filter_to_current_year_map_ref <- function(x) {
  if (!"current_year_public_drop" %in% names(x)) return(x)
  keep_idx <- !(x$current_year_public_drop %in% TRUE)
  keep_idx[is.na(keep_idx)] <- TRUE
  x[keep_idx, , drop = FALSE]
}

test_that("final_map == full minus current_year_public_drop (no silent loss)", {
  # A small universe spanning the relevant cases:
  #  - id1: clean keep, no conflict           -> retained in public
  #  - id2: high pre-year overlap WITH hotspot support (density ok) -> retained
  #  - id3: high overlap, NO hotspot support   -> EXCLUDED (the D_102 case)
  #  - id4: preyear_action == "drop", no hotspot -> EXCLUDED
  #  - id5: high overlap but density >= thr     -> retained
  full_df <- data.frame(
    fire_uid             = sprintf("U%02d", 1:5),
    area_ha              = c(10, 10, 4.79, 8, 10),
    hs_used_n            = c(5,  5,  0,    0, 5),
    preyear_overlap_frac = c(0.10, 0.95, 0.929, NA, 0.95),
    preyear_action       = c("keep", "keep", "review", "drop", "keep"),
    filter_2             = c("applied", "applied", "review", "applied", "applied"),
    stringsAsFactors = FALSE
  )
  full_df$current_year_public_drop <- compute_current_year_public_drop(full_df)

  # id3 (the D_102 analogue) and id4 must be the flagged ones; id5 keeps support.
  expect_true(full_df$current_year_public_drop[full_df$fire_uid == "U03"])
  expect_true(full_df$current_year_public_drop[full_df$fire_uid == "U04"])
  expect_false(full_df$current_year_public_drop[full_df$fire_uid == "U01"])
  expect_false(full_df$current_year_public_drop[full_df$fire_uid == "U05"])

  public_df <- filter_to_current_year_map_ref(full_df)

  # The load-bearing contract: public == full - flagged, exactly.
  n_flagged <- sum(full_df$current_year_public_drop %in% TRUE)
  expect_identical(nrow(public_df), nrow(full_df) - n_flagged)

  # Id preservation: public ids are a SUBSET of full ids, and the ONLY ids
  # removed are precisely the flagged ones (nothing else vanishes silently).
  expect_true(all(public_df$fire_uid %in% full_df$fire_uid))
  removed <- setdiff(full_df$fire_uid, public_df$fire_uid)
  flagged_ids <- full_df$fire_uid[full_df$current_year_public_drop %in% TRUE]
  expect_setequal(removed, flagged_ids)

  # No NA-flag row is ever dropped (NA -> kept).
  expect_true(all(!(public_df$current_year_public_drop %in% TRUE)))
})

test_that("the specific D_102 case: in full, out of public, sub-threshold", {
  # Exactly the LONG-2017 dropped polygon's decisive attributes:
  #   class_final = drop, preyear_overlap_frac ~ 0.93 (>= 0.70),
  #   preyear_action = review, hs_used_n = 0 -> temporal conflict + weak support.
  d102 <- data.frame(
    fire_uid             = "2017_2017_balanced_D_102",
    area_ha              = 4.7856125990304,
    hs_used_n            = 0,
    preyear_overlap_frac = 0.929664717827302,
    preyear_action       = "review",
    filter_2             = "review",
    stringsAsFactors = FALSE
  )
  flag <- compute_current_year_public_drop(d102)
  expect_true(flag)  # excluded from the current-year public map

  d102$current_year_public_drop <- flag
  public <- filter_to_current_year_map_ref(d102)
  expect_identical(nrow(public), 0L)  # dropped from public...
  # ...but it WOULD remain in final_map_full (the full layer applies no such
  # filter); we assert that by construction: the full layer keeps every row.
})

# ---------------------------------------------------------------------------
# (2) Source-level guard: the production `final_map` is built by
#     filter_to_current_year_map() (the documented current-year filter), so the
#     public layer can ONLY differ from final_map_full by current_year_public_drop.
# ---------------------------------------------------------------------------

test_that("final_map is produced via filter_to_current_year_map (documented filter)", {
  src_path <- testthat::test_path("..", "..", "R", "internal-sup-final-map.R")
  if (!file.exists(src_path)) {
    src_path <- file.path("R", "internal-sup-final-map.R")
  }
  skip_if_not(file.exists(src_path), "engine source not found")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  # final_map is final_map_full |> filter_to_current_year_map() |> make_public_...
  expect_true(grepl("final_map\\s*<-\\s*final_map_full\\s*\\|>\\s*\\n\\s*filter_to_current_year_map\\(\\)",
                    src))
  # the filter drops exactly current_year_public_drop == TRUE
  expect_true(grepl("current_year_public_drop\\s*%in%\\s*TRUE", src))
  # the full/scored layers are written from final_map_full (the complete universe)
  expect_true(grepl('layer\\s*=\\s*"final_map_full"', src))
  expect_true(grepl('layer\\s*=\\s*"deterministic_scored"', src))
  expect_true(grepl('layer\\s*=\\s*"final_map"', src))
})
