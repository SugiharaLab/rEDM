# data.frame validation and base-graphics plotting.

test_that("isValidDataFrame and ValidateDataFrame classify inputs", {
  cir <- LoadSampleData("circle.csv")
  expect_true(isValidDataFrame(cir))
  expect_false(isValidDataFrame(data.frame(t = 1:3, v = c("a","b","c"),
                                           stringsAsFactors = FALSE)))
  expect_false(isValidDataFrame(list(a = 1)))
  expect_true(ValidateDataFrame(cir, "x", "y"))
  expect_false(ValidateDataFrame(cir, "nope", "y"))
  expect_false(ValidateDataFrame(cir, "x", "missing"))
})

test_that("the API rejects invalid dataFrame / columns / target", {
  cir <- LoadSampleData("circle.csv")
  expect_error(Simplex(cir, "nope", "x", c(1,100), c(101,195), E = 2))
  expect_error(SMap(cir, "x", "gone", c(1,100), c(110,160), E = 2, theta = 3))
  expect_error(CCM(cir, "x", "absent", E = 2, libSizes = c(20,40), sample = 2))
})

test_that("showPlot renders without error", {
  cir <- LoadSampleData("circle.csv")
  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  # Coefficient labels use the partial-derivative symbol; base graphics may warn
  # in non-UTF-8 locales, so assert absence of errors rather than silence.
  suppressWarnings({
    expect_error(
      Simplex(cir, "x", "x", c(1,100), c(101,195), E = 2, backend = "brute", showPlot = TRUE), NA)
    expect_error(
      SMap(cir, "x", "x", c(1,100), c(110,160), E = 2, theta = 3, backend = "brute", showPlot = TRUE), NA)
    d <- CCM(cir, "x", "y", E = 2, Tp = 1, tau = -1, libSizes = c(20,60,120),
             sample = 3, seed = 1, numProcess = 1, backend = "brute", showPlot = TRUE)
  })
  expect_true(is.data.frame(d))
})
