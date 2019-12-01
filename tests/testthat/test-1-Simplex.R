# NOTE: Numerical tests are performed in cppEDM unit tests

context("Simplex test")

data( block_3sp )

test_that("Simplex embedded works", {
    S.df <- Simplex( dataFrame = block_3sp,
                     lib = "1 99", pred = "100 195",
                     E = 3, embedded = TRUE, showPlot = FALSE,
                     columns = "x_t y_t z_t", target = "x_t" )
    expect_s3_class(S.df, "data.frame")
    expect_true("time"         %in% names(S.df))
    expect_true("Observations" %in% names(S.df))
    expect_true("Predictions"  %in% names(S.df))
    expect_equal( dim(S.df), c(97,3) )
    Err <- ComputeError( S.df $ Observations, S.df $ Predictions )
    expect_true("MAE"  %in% names(Err))
    expect_true("rho"  %in% names(Err))
    expect_true("RMSE" %in% names(Err))
})

test_that("Simplex embedding works", {
    S.df <- Simplex( dataFrame = block_3sp,
                     lib = "1 99", pred = "100 195",
                     E = 3, embedded = FALSE, showPlot = FALSE,
                     columns = "x_t", target = "x_t" )
    expect_s3_class(S.df, "data.frame")
    expect_true("time"         %in% names(S.df))
    expect_true("Observations" %in% names(S.df))
    expect_true("Predictions"  %in% names(S.df))
    expect_equal( dim(S.df), c(97,3) )
})

test_that("Simplex errors", {
    expect_error( Simplex() )
    expect_error( Simplex( dataFrame = block_3sp ) )
    expect_error( Simplex( dataFrame = block_3sp,
                           lib = "1 99", pred = "100 195",
                           E = 3, columns = "x_t y_t z_t", target = "None" ) )
    expect_error( Simplex( dataFrame = block_3sp,
                           lib = "1 99", pred = "100 195",
                           E = 3, columns = "None", target = "x_t" ) )
    expect_error( Simplex( dataFrame = block_3sp,
                           lib = "1 99", pred = "100 200",
                           E = 3, columns = "x_t y_t z_t", target = "x_t" ) )
})
