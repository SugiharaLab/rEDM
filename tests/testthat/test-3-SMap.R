# SMap validation against the pyEDM SMap_* fixtures (RANN production backend).

BK <- "RANN"

SMapDiffs <- function(res, fixture, a, b, digits = 6) {
  ref <- LoadFixture(fixture)
  df  <- res$predictions
  b   <- min(b, nrow(df), nrow(ref))
  r1 <- round(df$Predictions[(a + 1):b], digits)
  r2 <- round(ref$Predictions[(a + 1):b], digits)
  sum(abs(r1 - r2) > 10^-digits | (is.na(r1) != is.na(r2)), na.rm = TRUE)
}

test_that("circle E4 (theta 3) reproduces fixture", {
  cir <- LoadSampleData("circle.csv")
  s <- SMap(cir, "x", "x", c(1,100), c(110,160), E=4, Tp=1, tau=-1,
            theta=3, backend=BK)
  expect_named(s, c("predictions","coefficients","singularValues"))
  expect_equal(SMapDiffs(s, "SMap_circle_E4", 1, 50), 0)
})

test_that("embedded E2 reproduces predictions and coefficient means", {
  cir <- LoadSampleData("circle.csv")
  s <- SMap(cir, c("x","y"), "x", c(1,200), c(1,200), E=2, Tp=1, tau=-1,
            embedded=TRUE, theta=3, backend=BK)
  expect_equal(SMapDiffs(s, "SMap_circle_E2_embd", 1, 195), 0)
  # coefficient columns are named with the partial-derivative convention
  cxx <- round(mean(s$coefficients[["\u2202x/\u2202x"]], na.rm = TRUE), 5)
  cxy <- round(mean(s$coefficients[["\u2202x/\u2202y"]], na.rm = TRUE), 5)
  expect_equal(cxx, 0.99801)
  expect_equal(cxy, 0.06311)
})

test_that("NaN in target vector reproduces fixture", {
  c2 <- LoadSampleData("circle.csv"); c2[c(6,7,13),2] <- NA; c2[c(11,12,18),3] <- NA
  s <- SMap(c2, "x", "y", c(1,50), c(1,50), E=2, Tp=1, tau=-1, theta=3, backend=BK)
  expect_equal(SMapDiffs(s, "SMap_nan", 1, 50), 0)
})

test_that("noTime index reproduces fixture", {
  cnt <- LoadSampleData("circle_noTime.csv")
  s <- SMap(cnt, "x", "y", c(1,100), c(101,150), E=2, theta=3, noTime=TRUE, backend=BK)
  expect_equal(SMapDiffs(s, "SMap_noTime", 1, 50), 0)
})

test_that("SMapSolve reproduces the OLS solution and returns E+1 singular values", {
  set.seed(1)
  X  <- cbind(1, matrix(rnorm(30), 10, 3))
  bb <- as.numeric(X %*% c(2, -1, 0.5, 3))
  sol <- SMapSolve(X, bb)
  expect_lt(max(abs(sol$C - unname(coef(lm(bb ~ X - 1))))), 1e-9)
  expect_equal(length(sol$SV), 4L)
})
