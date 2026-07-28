# CCM validation (RANN production backend).
#
# pyEDM subsamples with NumPy PCG64/SeedSequence, which R cannot bit-match, so
# validation is tiered: (1) EXACT at full library (RNG-independent) against an
# independent Simplex leave-one-out cross-map; (2) DISTRIBUTIONAL Monte-Carlo
# bands vs the pyEDM fixtures; (3) reproducibility / dispatch invariance;
# (4) structure (multivariate space-named columns).

BK <- "RANN"

CcmBand <- function(df, fixture, tol) {
  ref <- LoadFixture(fixture)
  max(abs(df[[2]] - ref[[2]]), abs(df[[3]] - ref[[3]])) < tol
}

test_that("full-library CCM equals an independent Simplex LOO cross-map", {
  cir <- LoadSampleData("circle.csv"); N <- nrow(cir)
  ccm <- CCM(cir, "x", "y", E=2, Tp=1, tau=-1, libSizes=c(N), sample=1, seed=1,
             numProcess=1, backend=BK)
  sx  <- Simplex(cir, "x", "y", c(1,N), c(1,N), E=2, Tp=1, tau=-1, backend=BK)
  expect_equal(round(ccm[[2]], 6),
               round(ComputeError(sx$Observations, sx$Predictions)$rho, 6))
})

test_that("fixed seed is reproducible and dispatch-invariant", {
  cir <- LoadSampleData("circle.csv")
  a <- CCM(cir,"x","y",E=2,Tp=1,tau=-1,libSizes=c(20,60,120),sample=20,seed=42,numProcess=1,backend=BK)
  b <- CCM(cir,"x","y",E=2,Tp=1,tau=-1,libSizes=c(20,60,120),sample=20,seed=42,numProcess=1,backend=BK)
  p <- CCM(cir,"x","y",E=2,Tp=1,tau=-1,libSizes=c(20,60,120),sample=20,seed=42,numProcess=4,backend=BK)
  expect_equal(a, b)
  expect_equal(a, p)
})

test_that("CCM means fall within Monte-Carlo bands of the pyEDM fixtures", {
  sa  <- LoadSampleData("sardine_anchovy_sst.csv")
  cir <- LoadSampleData("circle.csv")
  lz  <- LoadSampleData("LorenzData1000.csv")
  c2  <- LoadSampleData("circle.csv"); c2[c(6,7,13),2] <- NA; c2[c(11,12,18),3] <- NA

  expect_true(CcmBand(
    CCM(sa,"anchovy","np_sst",E=3,Tp=0,tau=-1,libSizes=c(10,20,30,40,50,60,70,75),
        sample=100,seed=777,numProcess=1,backend=BK), "CCM_anch_sst_valid", 0.06))
  expect_true(CcmBand(
    CCM(lz,"V3 V5","V1",E=5,Tp=10,tau=-5,libSizes=c(20,200,500,950),
        sample=30,seed=777,numProcess=1,backend=BK), "CCM_Lorenz5D_MV_valid", 0.05))
  expect_true(CcmBand(
    CCM(c2,"x","y",E=2,Tp=5,tau=-1,libSizes=c(10,190,10),
        sample=20,seed=777,numProcess=1,backend=BK), "CCM_nan_valid", 0.08))
  expect_true(CcmBand(
    CCM(cir,"x","y",E=2,Tp=-5,tau=-1,libSizes=c(20,200,50),
        sample=10,seed=777,numProcess=1,backend=BK), "CCM_NegativeTp", 0.08))
  expect_true(CcmBand(
    CCM(lz,c("V1","V2","V3"),"V5",E=0,Tp=0,tau=-1,embedded=TRUE,libSizes=c(50,950,100),
        sample=20,seed=777,numProcess=1,backend=BK), "CCM_embedded", 0.06))
  expect_true(CcmBand(
    CCM(cir,"x","y",E=2,Tp=3,tau=-1,exclusionRadius=5,libSizes=c(20,200,30),
        sample=10,seed=777,numProcess=1,backend=BK), "CCM_exclusionRadius", 0.08))
  expect_true(CcmBand(
    CCM(cir,"x","y",E=2,Tp=0,tau=3,libSizes=c(20,200,30),
        sample=10,seed=777,numProcess=1,backend=BK), "CCM_positiveTau", 0.08))
})

test_that("multivariate column names with spaces produce correct output columns", {
  sp <- LoadSampleData("columnNameSpace.csv")
  df <- CCM(sp, c("Var 1","Var3","Var 5 1"), c("Var 2","Var 4 A"),
            E=5, Tp=0, tau=-1, libSizes=c(20,50,90), sample=1, seed=777,
            numProcess=1, backend=BK)
  expect_equal(names(df), c("LibSize", "Var 1:Var 2", "Var 2:Var 1"))
  expect_true(all(is.finite(df[[2]])) && all(is.finite(df[[3]])))
})

test_that("includeData adds sample-variance columns", {
  cir <- LoadSampleData("circle.csv")
  df <- CCM(cir,"x","y",E=2,Tp=1,tau=-1,libSizes=c(20,60),sample=15,seed=7,
            includeData=TRUE,numProcess=1,backend=BK)
  expect_true(all(c("x:y_var","y:x_var") %in% names(df)))
})
