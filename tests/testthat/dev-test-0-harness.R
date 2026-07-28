test_that("identity: a fixture equals itself (NaN-as-equal)", {
  fixture <- LoadFixture("Smplx_E3_block_3sp")   # has NA cells in row 1
  expect_true(anyNA(fixture$Predictions))         # confirm NA present
  cmp <- CompareEDM(fixture, fixture, digits = 6)
  expect_true(cmp$equal)
})

test_that("sub-threshold perturbation still compares equal at 6 digits", {
  fixture <- LoadFixture("Smplx_E3_block_3sp")
  perturbed <- fixture; perturbed$Predictions[2] <- perturbed$Predictions[2] + 1e-9
  expect_true(CompareEDM(perturbed, fixture, digits = 6)$equal)
})

test_that("over-threshold perturbation is detected", {
  fixture <- LoadFixture("Smplx_E3_block_3sp")
  perturbed <- fixture; perturbed$Predictions[2] <- perturbed$Predictions[2] + 1e-3
  cmp <- CompareEDM(perturbed, fixture, digits = 6)
  expect_false(cmp$equal)
  expect_equal(nrow(cmp$diffs), 1L)
})

test_that("NaN mismatch is detected", {
  fixture <- LoadFixture("Smplx_E3_block_3sp")
  perturbed <- fixture; perturbed$Predictions[2] <- NA   # was finite
  expect_false(CompareEDM(perturbed, fixture, digits = 6)$equal)
})

test_that("column-name mismatch fails; reconcile fixes theta/Theta", {
  reference <- LoadFixture("PredictNonlinear")   # columns: theta, rho
  expect_identical(colnames(reference), c("theta", "rho"))
  result <- reference; colnames(result)[1] <- "Theta"   # rEDM emits capitalized
  expect_false(CompareEDM(result, reference, digits = 6)$equal)
  expect_true(CompareEDM(result, reference, digits = 6,
                         reconcile = c(theta = "Theta"))$equal)
})

test_that("CCM fixture round-2 identity holds", {
  fixture <- LoadFixture("CCM_anch_sst")
  expect_true(CompareEDM(fixture, fixture, digits = 2)$equal)
})

test_that("row slicing restricts the comparison", {
  fixture <- LoadFixture("Smplx_E3_block_3sp")
  perturbed <- fixture; perturbed$Predictions[1] <- 999   # differs outside slice
  expect_false(CompareEDM(perturbed, fixture, digits = 6)$equal)
  expect_true(CompareEDM(perturbed, fixture, digits = 6,
                         rows = 2:nrow(fixture))$equal)
})

test_that("ParseIndexPairs expands 1-based (start end) pairs", {
  expect_equal(ParseIndexPairs("1 100"), 1:100)
  expect_equal(ParseIndexPairs("1 100 201 300"), c(1:100, 201:300))
  expect_equal(ParseIndexPairs(c(1, 5)), 1:5)
  expect_error(ParseIndexPairs("1 2 3"))     # odd count
  expect_error(ParseIndexPairs("0 5"))       # <1 start
})

test_that("FlattenToString preserves rEDM delimiter convention", {
  expect_equal(FlattenToString(c(1, 2, 3)), "1 2 3")
  expect_equal(FlattenToString(c("a", "b")), "a b")
})
