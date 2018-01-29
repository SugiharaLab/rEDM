context("Check helper functions")

test_that("make_block produces desired output", {
    df <- data.frame(x = c(1, 4, 5, 8, 7, 8, 4, 2, 5, 2, 5, 7 ),
                     y = c(5, 7, 3, 9, 3, 2, 5, 1, 0, 8, 4, 6 ))
    lib <- matrix(c(1, 4, 5, 12), ncol = 2, byrow = TRUE)
    
    lag_one_test <- data.frame(
        time   = c( 1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11, 12),
        x      = c( 1, 4, 5, 8,  7, 8, 4, 2, 5,  2,  5,  7),
        x_1    = c(NA, 1, 4, 5, NA, 7, 8, 4, 2,  5,  2,  5),
        y      = c( 5, 7, 3, 9,  3, 2, 5, 1, 0,  8,  4,  6),
        y_1    = c(NA, 5, 7, 3, NA, 3, 2, 5, 1,  0,  8,  4)
    )
    lag_one <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 1)
    testthat::expect_equal(lag_one, lag_one_test)
    
    lag_two_test <- data.frame(
        time   = c( 1,  2, 3, 4,  5,  6, 7, 8, 9, 10, 11, 12),
        x      = c( 1,  4, 5, 8,  7,  8, 4, 2, 5,  2,  5,  7),
        x_2    = c(NA, NA, 1, 4, NA, NA, 7, 8, 4,  2,  5,  2),
        y      = c( 5,  7, 3, 9,  3,  2, 5, 1, 0,  8,  4,  6),
        y_2    = c(NA, NA, 5, 7, NA, NA, 3, 2, 5,  1,  0,  8)
    )
    lag_two <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 2)
    testthat::expect_equal(lag_two, lag_two_test)
    
    lag_three_test <- data.frame(
        time   = c( 1,  2,  3, 4,  5,  6,  7, 8, 9, 10, 11, 12),
        x      = c( 1,  4,  5, 8,  7,  8,  4, 2, 5,  2,  5,  7),
        x_3    = c(NA, NA, NA, 1, NA, NA, NA, 7, 8,  4,  2,  5),
        y      = c( 5,  7,  3, 9,  3,  2,  5, 1, 0,  8,  4,  6),
        y_3    = c(NA, NA, NA, 5, NA, NA, NA, 3, 2,  5,  1,  0)
    )
    lag_three <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 3)
    testthat::expect_equal(lag_three, lag_three_test)
})
