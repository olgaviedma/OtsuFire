# Tests for the packaged candidate-first omission/commission taxonomy
# (R/supervised-omission-taxonomy.R, Phase 1).

test_that("of_classify returns the 7 candidate-first classes correctly", {
  seed <- 310
  # has-candidate -> dominant of keep/review/drop
  expect_equal(as.character(of_classify(5, 0, 0, 5, TRUE, 320, seed)), OF_LEVELS[3]) # keep
  expect_equal(as.character(of_classify(0, 5, 0, 5, TRUE, 320, seed)), OF_LEVELS[1]) # review
  expect_equal(as.character(of_classify(0, 0, 5, 5, TRUE, 320, seed)), OF_LEVELS[2]) # drop
  # no-candidate -> physical RBR/burnability cascade
  expect_equal(as.character(of_classify(0, 0, 0, 0, TRUE,  320, seed)), OF_LEVELS[4]) # algo failure RBR>=seed
  expect_equal(as.character(of_classify(0, 0, 0, 0, TRUE,  200, seed)), OF_LEVELS[5]) # below detection
  expect_equal(as.character(of_classify(0, 0, 0, 0, TRUE,  100, seed)), OF_LEVELS[6]) # reference too weak
  expect_equal(as.character(of_classify(0, 0, 0, 0, FALSE, 320, seed)), OF_LEVELS[7]) # non-burnable
  # all results are factors with the canonical 7 levels
  expect_identical(levels(of_classify(0, 0, 0, 0, TRUE, 320, seed)), OF_LEVELS)
})

test_that("final-candidate absence is evaluated BEFORE any RBR rule (candidate-first)", {
  seed <- 310
  # review candidate present but RBR=100 (< WEAK): RBR-first would call it
  # 'reference too weak'; candidate-first MUST call it 'Candidate in review'.
  expect_equal(as.character(of_classify(0, 5, 0, 5, TRUE, 100, seed)), OF_LEVELS[1])
  # drop candidate present, RBR<seed: 'Candidate dropped', NOT 'below detection'.
  expect_equal(as.character(of_classify(0, 0, 5, 5, TRUE, 200, seed)), OF_LEVELS[2])
  # a present candidate dominates even non-burnable / very low RBR (candidate-first).
  expect_equal(as.character(of_classify(5, 0, 0, 5, FALSE, 50, seed)), OF_LEVELS[3])
  # sliver boundary: below OF_SLIVER -> no candidate; at/above -> candidate.
  expect_equal(as.character(of_classify(0, 0.005, 0, 0.005, TRUE, 320, seed)), OF_LEVELS[4])
  expect_equal(as.character(of_classify(0, 0.02,  0, 0.02,  TRUE, 320, seed)), OF_LEVELS[1])
})

test_that("no-candidate classes are invariant unsup vs sup on the same EFFIS population", {
  seed     <- rep(310, 5)
  burnable <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
  rbr      <- c(320,  200,  320,  100,  320)   # fires 3/4/5 are no-candidate physics

  # UNSUP: fires 1 (review) & 2 (drop) have candidates; 3,4,5 have none.
  keep_u <- c(0, 0, 0, 0, 0); rev_u <- c(5, 0, 0, 0, 0); drop_u <- c(0, 5, 0, 0, 0)
  lu <- of_classify(keep_u, rev_u, drop_u, keep_u + rev_u + drop_u, burnable, rbr, seed)

  # SUP: model recovered fire 1 (review -> keep); no-candidate fires 3,4,5 unchanged.
  keep_s <- c(5, 0, 0, 0, 0); rev_s <- c(0, 0, 0, 0, 0); drop_s <- c(0, 5, 0, 0, 0)
  ls <- of_classify(keep_s, rev_s, drop_s, keep_s + rev_s + drop_s, burnable, rbr, seed)

  # the no-candidate subset is IDENTICAL between models...
  expect_equal(as.character(lu[3:5]), as.character(ls[3:5]))
  # ...and equals the pure physical cascade.
  expect_equal(as.character(lu[3:5]), c(OF_LEVELS[4], OF_LEVELS[6], OF_LEVELS[7]))
  # only the recoverable has-candidate bucket changed (review -> keep).
  expect_equal(as.character(lu[1]), OF_LEVELS[1])
  expect_equal(as.character(ls[1]), OF_LEVELS[3])
})

test_that("the legacy 5-class taxonomy is not shipped; the 7-class is canonical", {
  expect_length(OF_LEVELS, 7L)
  expect_length(OF_COM_LEVELS, 4L)
  legacy <- c("Below detection threshold", "Detection failure (algorithm)",
              "Decision failure (recoverable)", "Reference error: weak/no burn signal",
              "Reference error: non-burnable")
  expect_false(any(legacy %in% OF_LEVELS))
  exports <- getNamespaceExports("OtsuFire")
  expect_false(any(c("ERR_LEVELS", "ERR_COLORS") %in% exports))  # no legacy taxonomy exported
  expect_true("of_classify" %in% exports)                       # the new classifier IS exported
})

test_that("of_classify_commission returns the 4 physical classes", {
  seed <- 310
  expect_equal(as.character(of_classify_commission(TRUE,  320, seed)), OF_COM_LEVELS[1]) # likely real
  expect_equal(as.character(of_classify_commission(TRUE,  200, seed)), OF_COM_LEVELS[2]) # borderline
  expect_equal(as.character(of_classify_commission(TRUE,  100, seed)), OF_COM_LEVELS[3]) # weak signal
  expect_equal(as.character(of_classify_commission(FALSE, 320, seed)), OF_COM_LEVELS[4]) # non-burnable
})

test_that("a scalar seed is recycled across fires", {
  res <- of_classify(c(0, 0), c(0, 0), c(0, 0), c(0, 0),
                     c(TRUE, TRUE), c(320, 200), seed = 310)
  expect_equal(as.character(res), c(OF_LEVELS[4], OF_LEVELS[5]))
})
