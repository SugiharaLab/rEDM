# NOTE: Numerical tests are performed in cppEDM unit tests

context("Multiview test")

data( block_3sp )

test_that("Multiview works", {
  M.List = Multiview( dataFrame = block_3sp,
                      lib = "1 99", pred = "105 190",
                      E = 3, columns = "x_t y_t z_t", target = "x_t" )
    
  expect_type(M.List, "list")
  expect_true("View"        %in% names(M.List))
  expect_true("Predictions" %in% names(M.List))
  expect_equal( dim(M.List $ View), c(9,9) )
  expect_equal( dim(M.List $ Predictions), c(87,3) )
})

test_that("Multiview errors", {
  expect_error( Multiview() )
  expect_error( Multiview( dataFrame = block_3sp,
                           lib = "1 99", pred = "105 190",
                           E = 3, columns = "x_t y_t z_t", target = "None" ) )
  expect_error( Multiview( dataFrame = block_3sp,
                           lib = "1 99", pred = "105 190",
                           E = 3, columns = "None", target = "x_t" ) )
  expect_error( Multiview( dataFrame = block_3sp,
                           lib = "1 99", pred = "105 201",
                           E = 3, columns = "None", target = "x_t" ) )
})
