# Extracted from test-2-Simplex.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "rEDM", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
BK <- "RANN"
SimplexDiffs <- function(df, fixture, a, b, digits = 6, exclude = integer(0)) {
  ref  <- LoadFixture(fixture)
  b    <- min(b, nrow(df), nrow(ref))
  rows <- setdiff((a + 1):b, exclude)
  r1 <- round(df$Predictions[rows], digits)
  r2 <- round(ref$Predictions[rows], digits)
  sum(abs(r1 - r2) > 10^-digits | (is.na(r1) != is.na(r2)), na.rm = TRUE)
}

# test -------------------------------------------------------------------------
e <- Embed(data.frame(x = 1:5), E = 3, tau = -1, columns = "x")
expect_s3_class(e, "data.frame")
expect_equal(names(e), c("x(t-0)", "x(t-1)", "x(t-2)"))
expect_equal(e[[1]], as.numeric(1:5))
expect_equal(names(Embed(data.frame(x = 1:8), E = 3, tau = -2, columns = "x")),
               c("x(t-0)", "x(t-2)", "x(t-4)"))
