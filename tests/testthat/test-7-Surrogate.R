# SurrogateData validation (rEDM signature; statistical/structural checks).

test_that("ebisuzaki preserves the power spectrum and variance", {
  x <- as.numeric(LoadSampleData("circle.csv")$x)
  set.seed(1)
  m <- SurrogateData(x, method = "ebisuzaki", num_surr = 20)
  expect_equal(dim(m), c(length(x), 20L))
  expect_lt(abs(mean(apply(m, 2, sd)) / sd(x) - 1), 0.02)
  expect_gt(cor(sort(Mod(fft(x))[-1]), sort(Mod(fft(m[,1]))[-1])), 0.999)
})

test_that("random_shuffle produces permutations", {
  x <- as.numeric(LoadSampleData("circle.csv")$x)
  set.seed(2)
  m <- SurrogateData(x, method = "random_shuffle", num_surr = 5)
  expect_true(all(apply(m, 2, function(c) identical(sort(c), sort(x)))))
})

test_that("seasonal and odd-length inputs are finite and correctly shaped", {
  y <- as.numeric(LoadSampleData("circle.csv")$y)
  set.seed(3)
  s <- SurrogateData(y, method = "seasonal", T_period = 12, num_surr = 5, alpha = 0.1)
  expect_equal(dim(s), c(length(y), 5L))
  expect_true(all(is.finite(s)))
  o <- SurrogateData(y[1:199], method = "ebisuzaki", num_surr = 3)
  expect_false(anyNA(o))
  expect_equal(nrow(o), 199L)
})

test_that("invalid method and non-finite input error", {
  x <- as.numeric(LoadSampleData("circle.csv")$x)
  expect_error(SurrogateData(x, method = "bogus"))
  expect_error(SurrogateData(c(1, 2, NA, 4), method = "ebisuzaki"))
})
