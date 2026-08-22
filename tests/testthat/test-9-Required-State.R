# Required E (pyEDM parity) and includeState / $internal.

test_that("E is required for Simplex, SMap, PredictInterval, PredictNonlinear", {
  cir <- LoadSampleData("circle.csv")
  expect_error(Simplex(cir, "x", "x", c(1,100), c(101,195), backend = "brute"),
               "E is required")
  expect_error(SMap(cir, "x", "x", c(1,100), c(110,160), theta = 3, backend = "brute"),
               "E is required")
  expect_error(PredictInterval(cir, "x", "x", c(1,100), c(101,195),
                               backend = "brute", numProcess = 1), "E is required")
  expect_error(PredictNonlinear(cir, "x", "x", c(1,100), c(101,195),
                                backend = "brute", numProcess = 1), "E is required")
})

test_that("embedded = TRUE infers E and needs no E argument", {
  cir <- LoadSampleData("circle.csv")
  s <- SMap(cir, c("x","y"), "y", c(1,100), c(101,200), theta = 3,
            embedded = TRUE, backend = "brute")
  expect_equal(ncol(s$coefficients), 4L)   # C0 + 2 partials + time
})

test_that("includeState returns the internal engine state", {
  cir <- LoadSampleData("circle.csv")
  s <- Simplex(cir, "x", "x", c(1,100), c(101,195), E = 2,
               backend = "brute", includeState = TRUE)
  expect_true(all(c("predictions", "internal") %in% names(s)))
  st <- s$internal
  expect_setequal(names(st), c("knn_neighbors","knn_distances","lib_i",
                               "pred_i","targetVec","embedding"))
  expect_identical(dim(st$knn_neighbors), dim(st$knn_distances))
  expect_equal(nrow(st$knn_neighbors), length(st$pred_i))

  sm <- SMap(cir, "x", "x", c(1,100), c(110,160), E = 2, theta = 3,
             backend = "brute", includeState = TRUE)
  expect_true(all(c("predictions","coefficients","singularValues","internal")
                  %in% names(sm)))
})

test_that("includeState composes with parameterList; plain call is inert", {
  cir <- LoadSampleData("circle.csv")
  b <- Simplex(cir, "x", "x", c(1,100), c(101,195), E = 2, backend = "brute",
               parameterList = TRUE, includeState = TRUE)
  expect_true(all(c("predictions","parameters","internal") %in% names(b)))
  expect_false(any(c("includeState","pathIn","showPlot") %in% names(b$parameters)))
  expect_true("E" %in% names(b$parameters))
  expect_s3_class(Simplex(cir, "x", "x", c(1,100), c(101,195), E = 2,
                          backend = "brute"), "data.frame")
})
