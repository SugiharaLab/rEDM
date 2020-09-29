# NOTE: Numerical tests are performed in cppEDM unit tests

context("CCM test")

data( sardine_anchovy_sst )

test_that("CCM works", {
    C.df = CCM( dataFrame = sardine_anchovy_sst,
                E = 3, Tp = 0, columns = "anchovy", target = "np_sst",
                libSizes = "10 70 10", sample = 100 )
    expect_s3_class(C.df, "data.frame")
    expect_true("LibSize"         %in% names(C.df))
    expect_true("anchovy:np_sst"  %in% names(C.df))
    expect_true("np_sst:anchovy"  %in% names(C.df))
    expect_equal( dim(C.df), c(7,3) )
})

test_that("CCM errors", {
    expect_error( CCM() )
    expect_error( CCM( dataFrame = sardine_anchovy_sst,
                       E = 3, Tp = 0, columns = "anchovy", target = "",
                       libSizes = "10 70 10", sample = 100 ) )
    expect_error( CCM( dataFrame = sardine_anchovy_sst,
                       E = 3, Tp = 0, columns = "", target = "np_sst",
                       libSizes = "10 70 10", sample = 100 ) )
    expect_error( CCM( dataFrame = sardine_anchovy_sst,
                       E = 3, Tp = 0, columns = "X", target = "np_sst",
                       libSizes = "10 70 10", sample = 100 ) )
    expect_error( CCM( dataFrame = sardine_anchovy_sst,
                       E = 3, Tp = 0, columns = "anchovy", target = "np_sst",
                       libSizes = "10 70 80", sample = 100 ) )
})
