# Multiview validation against pyEDM fixtures (RANN production backend).
# Built on the deterministic Simplex path, so predictions and the ranked view
# table reproduce the fixtures exactly.

BK <- "RANN"

test_that("Multiview reproduces averaged predictions and view ranking", {
  d <- LoadSampleData("block_3sp.csv")
  M <- Multiview(d, "x_t y_t z_t", "x_t", lib = c(1,100), pred = c(101,198),
                 D = 0, E = 3, Tp = 1, tau = -1, multiview = 0, trainLib = FALSE,
                 numProcess = 1, backend = BK)
  expect_named(M, c("Predictions", "View"))

  refp <- LoadFixture("Multiview_pred")
  expect_equal(round(M$Predictions$Predictions, 4), round(refp$Predictions, 4))

  refc <- LoadFixture("Multiview_combos")
  expect_equal(nrow(M$View), 9L)
  expect_equal(round(M$View$rho,  4), round(refc$rho,  4))
  expect_equal(round(M$View$MAE,  4), round(refc$MAE,  4))
  expect_equal(round(M$View$RMSE, 4), round(refc$RMSE, 4))
})

test_that("Multiview parallel and sequential agree", {
  d <- LoadSampleData("block_3sp.csv")
  a <- Multiview(d, "x_t y_t z_t", "x_t", lib=c(1,100), pred=c(101,198),
                 D=0, E=3, Tp=1, tau=-1, trainLib=FALSE, numProcess=1, backend=BK)
  b <- Multiview(d, "x_t y_t z_t", "x_t", lib=c(1,100), pred=c(101,198),
                 D=0, E=3, Tp=1, tau=-1, trainLib=FALSE, numProcess=4, backend=BK)
  expect_equal(a$Predictions$Predictions, b$Predictions$Predictions)
  expect_equal(a$View$rho, b$View$rho)
})
