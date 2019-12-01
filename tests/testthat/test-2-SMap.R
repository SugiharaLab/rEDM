# NOTE: Numerical tests are performed in cppEDM unit tests

context("SMap test")

data( circle )

test_that("SMap works", {
    S.List = SMap( dataFrame = circle,
                   lib = "1 100", pred = "110 190", theta = 4, E = 2,
                   embedded = TRUE, columns = "x y", target = "x" )
    expect_type(S.List, "list")
    expect_true("predictions"  %in% names(S.List))
    expect_true("coefficients" %in% names(S.List))
    expect_equal( dim(S.List $ predictions  ), c(82,3) )
    expect_equal( dim(S.List $ coefficients ), c(82,4) )
})

test_that("SMap errors", {
    expect_error( SMap() )
    expect_error( SMap( dataFrame = circle,
                        lib = "1 100", pred = "110 190", theta = 4, E = 2,
                        embedded = TRUE, columns = "x y", target = "None" ) )
    expect_error( SMap( dataFrame = circle,
                        lib = "1 100", pred = "110 190", theta = 4, E = 2,
                        embedded = TRUE, columns = "None", target = "x" ) )
    expect_error( SMap( dataFrame = circle,
                        lib = "1 100", pred = "110 201", theta = 4, E = 2,
                        embedded = TRUE, columns = "x y", target = "x" ) )
})
