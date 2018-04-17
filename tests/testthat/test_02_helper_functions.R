context("Check helper functions")

test_that("rEDM_warning filters warnings", {
    expect_warning(rEDM_warning("test ABC123"))
    expect_warning(rEDM_warning("test ABC123"), "^test ABC123$")
    expect_warning(rEDM_warning("test ABC123", silent = TRUE), NA)
})

test_that("check_params_against_lib produces desired output", {
    lib <- matrix(c(1, 5), ncol = 2)
    expect_true(check_params_against_lib(3, 1, 1, lib))
    expect_warning(value <- check_params_against_lib(5, 1, 1, lib))
    expect_false(value)
    expect_true(check_params_against_lib(3, 2, 0, lib))
    expect_warning(value <- check_params_against_lib(3, 2, 1, lib))
    expect_false(value)
    expect_true(check_params_against_lib(1, 1, 4, lib))
    expect_true(check_params_against_lib(1, 1, -4, lib))
    expect_warning(value <- check_params_against_lib(1, 1, 5, lib))
    expect_false(value)
    expect_warning(value <- check_params_against_lib(1, 1, -5, lib))
    expect_false(value)
})

test_that("convert_to_column_indices produces desired output", {
    df <- data.frame(a = 1:10, b = 2:11, c = rep(0, 10))
    expect_equal(convert_to_column_indices(1, df), 1)
    expect_equal(convert_to_column_indices(1:3, df), 1:3)
    expect_warning(cols <- convert_to_column_indices(c(1, 2, 4), df))
    expect_equal(cols, 1:2)
    expect_warning(convert_to_column_indices(1:5, df, silent = TRUE), NA)
    expect_error(convert_to_column_indices(6:7, df, silent = TRUE))
    
    expect_equal(convert_to_column_indices("a", df), 1)
    expect_equal(convert_to_column_indices(c("a", "b", "c"), df), 1:3)
    expect_warning(cols <- convert_to_column_indices(c("a", "b", "d"), df))
    expect_equal(cols, 1:2)
    expect_warning(convert_to_column_indices(letters[1:5], df, silent = TRUE), 
                   NA)
    expect_error(convert_to_column_indices(letters[6:7], df, silent = TRUE))
})

test_that("coerce_lib produces desired output", {
    expect_equal(coerce_lib(c(1, 20)), matrix(c(1, 20), ncol = 2))
    lib <- c(20, 10)
    expect_warning(coerce_lib(lib), "the lib argument")
    pred <- c(50, 5)
    expect_warning(coerce_lib(pred), "the pred argument")
})

test_that("make_block produces desired output", {
    out_actual <- data.frame(
        time = 1:100, 
        col1 = 1:100, 
        col1_1 = c(NA, 1:99), 
        col1_2 = c(NA, NA, 1:98)
    )
    expect_error(out <- make_block(1:100), NA)
    expect_equal(out, out_actual)
    
    df <- data.frame(x = c(1, 4, 5, 8, 7, 8, 4, 2, 5, 2, 5, 7 ),
                     y = c(5, 7, 3, 9, 3, 2, 5, 1, 0, 8, 4, 6 ))
    lib <- matrix(c(1, 4, 5, 12), ncol = 2, byrow = TRUE)
    
    lag_one_actual <- data.frame(
        time   = c( 1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11, 12),
        x      = c( 1, 4, 5, 8,  7, 8, 4, 2, 5,  2,  5,  7),
        x_1    = c(NA, 1, 4, 5, NA, 7, 8, 4, 2,  5,  2,  5),
        y      = c( 5, 7, 3, 9,  3, 2, 5, 1, 0,  8,  4,  6),
        y_1    = c(NA, 5, 7, 3, NA, 3, 2, 5, 1,  0,  8,  4)
    )
    lag_one <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 1)
    expect_equal(lag_one, lag_one_actual)
    
    lag_one_short <- make_block(df, max_lag = 2, t = NULL, lib = c(1, 4), tau = 1)
    expect_equal(lag_one_short, lag_one_actual[1:4, ])
    lag_one_long <- make_block(df, max_lag = 2, t = NULL, tau = 1, 
                               lib =　matrix(c(1, 5, 4, 8), ncol = 2), 
                               restrict_to_lib = FALSE)
    expect_equal(lag_one_long, lag_one_actual)
    
    lag_two_actual <- data.frame(
        time   = c( 1,  2, 3, 4,  5,  6, 7, 8, 9, 10, 11, 12),
        x      = c( 1,  4, 5, 8,  7,  8, 4, 2, 5,  2,  5,  7),
        x_2    = c(NA, NA, 1, 4, NA, NA, 7, 8, 4,  2,  5,  2),
        y      = c( 5,  7, 3, 9,  3,  2, 5, 1, 0,  8,  4,  6),
        y_2    = c(NA, NA, 5, 7, NA, NA, 3, 2, 5,  1,  0,  8)
    )
    lag_two <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 2)
    expect_equal(lag_two, lag_two_actual)
    
    lag_three_actual <- data.frame(
        time   = c( 1,  2,  3, 4,  5,  6,  7, 8, 9, 10, 11, 12),
        x      = c( 1,  4,  5, 8,  7,  8,  4, 2, 5,  2,  5,  7),
        x_3    = c(NA, NA, NA, 1, NA, NA, NA, 7, 8,  4,  2,  5),
        y      = c( 5,  7,  3, 9,  3,  2,  5, 1, 0,  8,  4,  6),
        y_3    = c(NA, NA, NA, 5, NA, NA, NA, 3, 2,  5,  1,  0)
    )
    lag_three <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 3)
    expect_equal(lag_three, lag_three_actual)
})

test_that("make_surrogate_shuffle works", {
    set.seed(42)
    expect_error(dat <- make_surrogate_shuffle(1:100, 15), NA)
    expect_equal(nrow(dat), 100)
    expect_equal(ncol(dat), 15)
    expect_equal(colSums(dat), rep.int(5050, 15))
    set.seed(42)
    expect_error(dat2 <- make_surrogate_data(1:100, "random_shuffle", 15), NA)
    expect_equal(dat, dat2)
})

test_that("make_surrogate_ebisuzaki works", {
    set.seed(42)
    expect_error(dat <- make_surrogate_ebisuzaki(1:100, 15), NA)
    expect_equal(nrow(dat), 100)
    expect_equal(ncol(dat), 15)
    set.seed(42)
    expect_error(dat2 <- make_surrogate_data(1:100, "ebisuzaki", 15), NA)
    expect_equal(dat, dat2)
})

test_that("make_surrogate_seasonal works", {
    set.seed(42)
    expect_error(dat <- make_surrogate_seasonal(1:100, 15, T_period = 4), NA)
    expect_equal(nrow(dat), 100)
    expect_equal(ncol(dat), 15)
    set.seed(42)
    expect_error(dat2 <- make_surrogate_data(1:100, "seasonal", 15, T_period = 4), NA)
    expect_equal(dat, dat2)
    set.seed(42)
    expect_error(dat3 <- make_surrogate_data(1:100, "seasonal", 15, T_period = 5), NA)
    expect_true(identical(dat, dat2))
    expect_false(identical(dat, dat3))
})

test_that("make_surrogate_twin works", {
    ts <- rnorm(100) + sin(1:100 * pi / 6)
    set.seed(42)
    expect_error(dat <- make_surrogate_twin(ts, 15, T_period = 12), NA)
    expect_equal(nrow(dat), 100)
    expect_equal(ncol(dat), 15)
    set.seed(42)
    expect_error(dat2 <- make_surrogate_data(ts, "twin", 15, T_period = 12), NA)
    expect_equal(dat, dat2)
    set.seed(42)
    expect_error(dat3 <- make_surrogate_data(ts, "twin", 15, T_period = 13))
})

