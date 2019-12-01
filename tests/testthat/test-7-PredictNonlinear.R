# NOTE: Numerical tests are performed in cppEDM unit tests

context("Predict Nonlinear test")

data( TentMapNoise )

test_that("PredictNonlinear works", {
    df <- PredictNonlinear( dataFrame = TentMapNoise,
                            E = 2, lib = "1 100", pred = "201 500",
                            columns = "TentMap", target = "TentMap",
                            showPlot = FALSE )
    expect_s3_class(df, "data.frame")
    expect_true("Theta" %in% names(df))
    expect_true("rho"   %in% names(df))
    expect_equal( dim(df), c(15,2) )
})

test_that("PredictNonlinear errors", {
    expect_error( PredictNonlinear() )
    expect_error( PredictNonlinear( dataFrame = TentMapNoise,
                                    E = 2, lib = "1 100", pred = "201 500",
                                    columns = "", target = "TentMap",
                                    showPlot = FALSE ) )
})
