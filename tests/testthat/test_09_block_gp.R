context("Check block_gp function")

data("two_species_model")
block <- two_species_model[1:200, ]

test_that("block_gp works", {
    expect_error(output <- block_gp(block, columns = c("x", "y"), 
                                    phi = 0.5, 
                                    v_e = seq(from = -300, to = 300, by = 50),
                                    eta = 7,
                                    fit_params = FALSE,
                                    first_column_time = TRUE, silent = TRUE),
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("embedding" %in% names(output))
    expect_true("tp" %in% names(output))
    expect_true("phi" %in% names(output))
    expect_true("v_e" %in% names(output))
    expect_true("eta" %in% names(output))
    expect_true("fit_params" %in% names(output))
    expect_true("num_pred" %in% names(output))
    expect_true("rho" %in% names(output))
    expect_true("mae" %in% names(output))
    expect_true("rmse" %in% names(output))
    expect_equal(NROW(output), 13)
    expect_equal(digest::digest(round(output$rho, 4)),
                 "00056d3065ac8649775ff8ee949385a5")
})

test_that("block_gp model_output works", {
    expect_error(output <- block_gp(block, columns = c("x", "y"), 
                                    fit_params = TRUE,
                                    phi = 0.5, v_e = -170, eta = 7, 
                                    first_column_time = TRUE, 
                                    stats_only = FALSE, 
                                    silent = TRUE),
                 NA)
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
    expect_equal(digest::digest(round(model_output, 4)),
                 "5faf8205dd2871e5c4cb9c84bc429651")
})

test_that("block_gp covariance matrix works", {
    expect_error(output <- block_gp(block, columns = c("x", "y"), 
                                    phi = 0.5, v_e = -170, eta = 7, 
                                    first_column_time = TRUE, 
                                    save_covariance_matrix = TRUE, 
                                    silent = TRUE),
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("model_output" %in% names(output))
    expect_true(is.list(output$model_output))
    expect_true("covariance_matrix" %in% names(output))
    expect_true(is.list(output$covariance_matrix))
    expect_error(covariance_matrix <- output$covariance_matrix[[1]], NA)
    expect_is(covariance_matrix, "matrix")
    expect_equal(dim(covariance_matrix), c(199, 199))
    expect_equal(digest::digest(round(covariance_matrix, 4)),
                 "df5a7a72abacf6f4ecb9791928177888")
})