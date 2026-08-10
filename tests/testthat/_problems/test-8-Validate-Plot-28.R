# Extracted from test-8-Validate-Plot.R:28

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "rEDM", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
cir <- LoadSampleData("circle.csv")
pdf(tempfile())
on.exit(dev.off(), add = TRUE)
suppressWarnings({
    expect_error(
      Simplex(cir, "x", "x", c(1,100), c(101,195), E = 2, backend = "brute", showPlot = TRUE), NA)
    expect_error(
      SMap(cir, "x", "x", c(1,100), c(110,160), E = 2, theta = 3, backend = "brute", showPlot = TRUE), NA)
    d <- CCM(cir, "x", "y", E = 2, Tp = 1, tau = -1, libSizes = c(20,60,120),
             sample = 3, seed = 1, numProcess = 1, backend = "brute", showPlot = TRUE)
  })
