# End-to-end Simplex validation against the 10 pyEDM Smplx_* fixtures.
# Uses the RANN production backend. Neighbours/exclusion validated in Phase 1.

BK <- "RANN"

# Count Predictions mismatches over the pyEDM slice [a:b] (R rows (a+1):b).
SimplexDiffs <- function(df, fixture, a, b, digits = 6, exclude = integer(0)) {
  ref  <- LoadFixture(fixture)
  b    <- min(b, nrow(df), nrow(ref))
  rows <- setdiff((a + 1):b, exclude)
  r1   <- round(df$Predictions[rows], digits)
  r2   <- round(ref$Predictions[rows], digits)
  sum(abs(r1 - r2) > 10^-digits | (is.na(r1) != is.na(r2)), na.rm = TRUE)
}

test_that("E3 block_3sp reproduces fixture exactly", {
  d <- LoadSampleData("block_3sp.csv")
  df <- Simplex(d, "x_t", "x_t", c(1,100), c(101,195), E=3, Tp=1, backend=BK)
  expect_equal(names(df), c("time","Observations","Predictions","Pred_Variance"))
  expect_equal(SimplexDiffs(df, "Smplx_E3_block_3sp", 1, 95), 0)
})

test_that("embedded input reproduces fixture", {
  d <- LoadSampleData("block_3sp.csv")
  df <- Simplex(d, "x_t y_t z_t", "x_t", c(1,99), c(100,198), E=3, embedded=TRUE, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_E3_embd_block_3sp", 1, 98), 0)
})

test_that("negative Tp reproduces fixture", {
  d <- LoadSampleData("block_3sp.csv")
  df <- Simplex(d, "x_t", "y_t", c(1,100), c(50,80), E=3, Tp=-2, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_negTp_block_3sp", 1, 98), 0)
})

test_that("validLib reproduces fixture", {
  cir <- LoadSampleData("circle.csv"); vlib <- cir$x > 0.5 | cir$x < -0.5
  df <- Simplex(cir, "x", "x", c(1,200), c(1,200), E=2, Tp=1, validLib=vlib, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_validLib", 1, 195), 0)
})

test_that("disjoint library reproduces fixture", {
  cir <- LoadSampleData("circle.csv")
  df <- Simplex(cir, "x", "x", c(1,40,50,130), c(80,170), E=2, Tp=1, tau=-3, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_disjointLib", 1, 195), 0)
})

test_that("disjoint pred with NaN reproduces fixture", {
  lz <- LoadSampleData("LorenzData1000.csv"); lz[c(9,51,502), c(2,3)] <- NA
  df <- Simplex(lz, "V1", "V2", c(1,50,101,200,251,500),
                c(1,10,151,155,551,555,881,885,991,1000), E=5, Tp=2, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_disjointPred_nan", 1, 195, digits=5), 0)
})

test_that("exclusionRadius reproduces fixture", {
  cir <- LoadSampleData("circle.csv")
  df <- Simplex(cir, "x", "y", c(1,100), c(21,81), E=2, Tp=1, exclusionRadius=5, backend=BK)
  expect_equal(SimplexDiffs(df, "Smplx_exclRadius", 1, 60), 0)
})

test_that("NaN propagation reproduces fixtures (both directions)", {
  c2 <- LoadSampleData("circle.csv"); c2[c(6,7,13),2] <- NA; c2[c(11,12,18),3] <- NA
  df1 <- Simplex(c2, "x", "y", c(1,100), c(1,95),  E=2, Tp=1, backend=BK)
  df2 <- Simplex(c2, "y", "x", c(1,200), c(1,195), E=2, Tp=1, backend=BK)
  expect_equal(SimplexDiffs(df1, "Smplx_nan",  1, 90),  0)
  expect_equal(SimplexDiffs(df2, "Smplx_nan2", 1, 190), 0)
})

test_that("DateTime index reproduces fixture (deterministic tie-break)", {
  s12 <- LoadSampleData("S12CD-S333-SumFlow_1980-2005.csv")
  df  <- Simplex(s12, "S12.C.D.S333", "S12.C.D.S333", c(1, 800), c(801, 1001),
                 E = 3, Tp = 1, backend = BK)
  # Formerly excluded rows 184-185 (a large exact-distance tie); now resolved
  # deterministically and equal to the regenerated reference fixture.
  expect_equal(SimplexDiffs(df, "Smplx_DateTime", 1, 200), 0)
})

test_that("MakeBlock builds the delay matrix", {
  mb <- MakeBlock(data.frame(x = 1:5), E = 3, tau = -1, columns = "x")
  expect_true(is.matrix(mb))
  expect_equal(mb[3, ], c(3, 2, 1))
  expect_true(is.na(mb[1, 2]) && is.na(mb[2, 3]))
})

test_that("Embed returns a named data.frame (pyEDM convention)", {
  e <- Embed(data.frame(x = 1:5), E = 3, tau = -1, columns = "x")
  expect_s3_class(e, "data.frame")
  expect_equal(names(e), c("x(t-0)", "x(t-1)", "x(t-2)"))
  expect_equal(e[[1]], as.numeric(1:5))
  # lag scales with |tau|
  expect_equal(names(Embed(data.frame(x = 1:8), E = 3, tau = -2, columns = "x")),
               c("x(t-0)", "x(t-2)", "x(t-4)"))
  # includeTime prepends the time column with its original name
  et <- Embed(data.frame(Date = 101:105, x = 1:5), E = 2, tau = -1,
              columns = "x", includeTime = TRUE)
  expect_equal(names(et), c("Date", "x(t-0)", "x(t-1)"))
  expect_equal(et[[1]], 101:105)
  # dataFile input matches in-memory input (round trip)
  d  <- data.frame(Time = 1:6, x = c(2, 4, 6, 8, 10, 12))
  tf <- tempfile(fileext = ".csv"); write.csv(d, tf, row.names = FALSE)
  expect_equal(
    Embed(dataFile = basename(tf), pathIn = dirname(tf), columns = "x", E = 3),
    Embed(d, E = 3, tau = -1, columns = "x"))
})

test_that("ComputeError matches pyEDM reductions", {
  ce <- ComputeError(as.numeric(1:10), as.numeric(1:10))
  expect_equal(ce$rho, 1); expect_equal(ce$RMSE, 0)
  expect_true(is.na(ComputeError(1:3, 1:3)$rho))     # < 5 pairs
})

test_that("neighbour ties (E=1 PSF salmon) reproduce fixture: RANN backend", {
  d  <- LoadSampleData("Smplx_NeighborTies.csv")
  df <- Simplex(d, "1137", "1135", c(1, 96), c(1, 96),
                E = 1, Tp = 0, tau = -1, exclusionRadius = 1,
                noTime = TRUE, backend = BK)

  expect_equal(names(df),
               c("Time", "Observations", "Predictions", "Pred_Variance"))
  # Full Predictions column, deterministic tie-break; no excluded rows.
  expect_equal(SimplexDiffs(df, "Smplx_NeighborTies", 0, 96, digits = 5), 0)
})

test_that("neighbour ties: brute reproduces fixture and equals RANN", {
  d <- LoadSampleData("Smplx_NeighborTies.csv")

  rann  <- Simplex(d, "1137", "1135", c(1, 96), c(1, 96),
                   E = 1, Tp = 0, tau = -1, exclusionRadius = 1,
                   noTime = TRUE, backend = "RANN")
  brute <- Simplex(d, "1137", "1135", c(1, 96), c(1, 96),
                   E = 1, Tp = 0, tau = -1, exclusionRadius = 1,
                   noTime = TRUE, backend = "brute")

  # brute reproduces the pyEDM-generated fixture ...
  expect_equal(SimplexDiffs(brute, "Smplx_NeighborTies", 0, 96, digits = 5), 0)

  # ... and both backends agree (tie-break is backend-independent by design).
  expect_equal(round(rann$Predictions, 5), round(brute$Predictions, 5))
})
