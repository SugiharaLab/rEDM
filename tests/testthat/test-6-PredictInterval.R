# NOTE: Numerical tests are performed in cppEDM unit tests

context("Predict Interval test")

data( TentMap )

test_that("PredictInterval works", {
    df <- PredictInterval( dataFrame = TentMap,
                           lib = "1 100", pred = "201 500", E = 2,
                           columns = "TentMap", target = "TentMap",
                           showPlot = FALSE )
    expect_s3_class(df, "data.frame")
    expect_true("Tp"  %in% names(df))
    expect_true("rho" %in% names(df))
    expect_equal( dim(df), c(10,2) )
})

test_that("PredictInterval errors", {
    expect_error( PredictInterval() )
    expect_error( PredictInterval( dataFrame = TentMap,
                                   lib = "1 100", pred = "201 500", E = 2,
                                   columns = "", target = "TentMap",
                                   showPlot = FALSE ) )
})
