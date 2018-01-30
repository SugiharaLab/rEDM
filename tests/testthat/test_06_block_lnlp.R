context("Check block_lnlp function")

data("two_species_model")
block <- two_species_model[1:200, ]

test_that("block_lnlp works", {
    expect_error(output <- block_lnlp(block, columns = c("x", "y"),
                                      first_column_time = TRUE, 
                                      silent = TRUE), 
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("embedding" %in% names(output))
    expect_true("tp" %in% names(output))
    expect_true("nn" %in% names(output))
    expect_true("num_pred" %in% names(output))
    expect_true("rho" %in% names(output))
    expect_true("mae" %in% names(output))
    expect_true("rmse" %in% names(output))
    expect_equal(NROW(output), 1)
    expect_equal(digest::digest(output), 
                 "b37f7bbfa0551ffff92b1c1e6fda3e28")
})

test_that("block_lnlp model_output works", {
    expect_warning(output <- block_lnlp(block, columns = c("x", "y"),
                                        first_column_time = TRUE, 
                                        stats_only = FALSE))
    expect_s3_class(output, "data.frame")
    expect_true("model_output" %in% names(output))
    expect_true(is.list(output$model_output))
    expect_error(model_output <- output$model_output[[1]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    expect_equal(dim(model_output), c(199, 4))
    expect_equal(digest::digest(output), 
                 "1dc9e215ad9c665f0606e1079242b939")
})

test_that("block_lnlp smap_coefficients works", {
    expect_warning(output <- block_lnlp(block, columns = c("x", "y"),
                                        first_column_time = TRUE, 
                                        method = "s-map", theta = 1, 
                                        save_smap_coefficients = TRUE))
    expect_s3_class(output, "data.frame")
    expect_true("smap_coefficients" %in% names(output))
    expect_true(is.list(output$smap_coefficients))
    expect_error(smap_coefficients <- output$smap_coefficients[[1]], NA)
    expect_s3_class(smap_coefficients, "data.frame")
    expect_true("c_1" %in% names(smap_coefficients))
    expect_true("c_2" %in% names(smap_coefficients))
    expect_true("c_0" %in% names(smap_coefficients))
    expect_equal(dim(smap_coefficients), c(199, 3))
    expect_equal(digest::digest(output), 
                 "91c5a2164e204dd65f600584c8348fce")
})

test_that("block_lnlp error checking works", {
    df <- data.frame(a = 1:5, b = 0:4)
    expect_warning(block_lnlp(df))
    expect_warning(block_lnlp(df, columns = 1:3))
    expect_warning(block_lnlp(df, columns = list(1, 1:3)))
    expect_error(block_lnlp(df, columns = list(0, 4:5), silent = TRUE))
    expect_error(block_lnlp(df, tp = 5, silent = TRUE))
    expect_error(block_lnlp(df, tp = -5, silent = TRUE))
    expect_error(block_lnlp(df, target_column = 0, silent = TRUE))
})