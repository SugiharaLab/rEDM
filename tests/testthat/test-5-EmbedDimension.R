# NOTE: Numerical tests are performed in cppEDM unit tests

context("Embed Dimension test")

data( TentMap )

test_that("EmbedDimension works", {
    df <- EmbedDimension( dataFrame = TentMap, lib = "1 100", pred = "201 500",
                          columns = "TentMap", target = "TentMap",
                          showPlot = FALSE )
    expect_s3_class(df, "data.frame")
    expect_true("E"   %in% names(df))
    expect_true("rho" %in% names(df))
    expect_equal( dim(df), c(10,2) )
})

test_that("EmbedDimension errors", {
    expect_error( EmbedDimension() )
    expect_error( EmbedDimension( dataFrame = TentMap,
                                  lib = "1 100", pred = "201 500",
                                  columns = "TentMap", target = "None",
                                  showPlot = FALSE ) )
})
