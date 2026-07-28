# Diagnostic sweep validation against pyEDM fixtures (RANN production backend).
# Column names follow the rEDM convention; the swept-parameter column and rho
# are compared to the pyEDM fixture (lower-case theta / 'Exclusion radius').

BK <- "RANN"

SweepMatch <- function(df, fixture, digits = 6) {
  ref <- LoadFixture(fixture)
  isTRUE(all.equal(round(df[[1]], digits), round(ref[[1]], digits))) &&
  isTRUE(all.equal(round(df$rho, digits), round(ref$rho, digits)))
}

test_that("EmbedDimension reproduces the Lorenz5D fixture", {
  lz <- LoadSampleData("LorenzData1000.csv")
  df <- EmbedDimension(lz, "V1", "V1", c(1,500), c(501,800), maxE=12, Tp=15,
                       tau=-5, exclusionRadius=20, numProcess=1, backend=BK)
  expect_equal(names(df), c("E", "rho"))
  expect_true(SweepMatch(df, "EmbedDim"))
})

test_that("PredictInterval reproduces the block_3sp fixture", {
  b3 <- LoadSampleData("block_3sp.csv")
  df <- PredictInterval(b3, "x_t", "x_t", c(1,150), c(151,198), maxTp=15, E=3,
                        tau=-1, numProcess=1, backend=BK)
  expect_equal(names(df), c("Tp", "rho"))
  expect_true(SweepMatch(df, "PredictInterval"))
})

test_that("PredictExclusionRadius reproduces the SumFlow fixture", {
  s12 <- LoadSampleData("S12CD-S333-SumFlow_1980-2005.csv")
  df <- PredictExclusionRadius(s12, "S12.C.D.S333", "S12.C.D.S333", c(1,900),
                               c(1,900), E=7, tau=-3, Tp=0, numProcess=1, backend=BK)
  expect_equal(names(df), c("ExclusionRadius", "rho"))
  expect_true(SweepMatch(df, "PredictExclusionRadius"))
})

test_that("PredictNonlinear reproduces the TentMap fixture", {
  tm <- LoadSampleData("TentMapNoise.csv")
  df <- PredictNonlinear(tm, "TentMap", "TentMap", c(1,500), c(501,800), E=4,
                         Tp=1, tau=-1,
                         theta=c(0.01,0.1,0.3,0.5,0.75,1,1.5,2,3,4,5,6,7,8,9,10,15,20),
                         numProcess=1, backend=BK)
  expect_equal(names(df), c("Theta", "rho"))
  expect_true(SweepMatch(df, "PredictNonlinear"))
})

test_that("parallel and sequential sweeps agree", {
  b3 <- LoadSampleData("block_3sp.csv")
  a <- PredictInterval(b3, "x_t", "x_t", c(1,150), c(151,198), maxTp=15, E=3,
                       tau=-1, numProcess=1, backend=BK)
  b <- PredictInterval(b3, "x_t", "x_t", c(1,150), c(151,198), maxTp=15, E=3,
                       tau=-1, numProcess=4, backend=BK)
  expect_equal(a, b)
})
