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
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "3f457044ffc65c1b7a6e54f47ebbcc97")
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
    expect_equal(dim(model_output), c(200, 4))
    expect_equal(digest::digest(round(model_output, 4)), 
                 "56a6aac91d858ff5a81c619ef0bf92f6")
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
    expect_equal(dim(smap_coefficients), c(200, 3))
    expect_equal(digest::digest(round(smap_coefficients, 4)), 
                 "82a3b6164cfcfe8d69d98689b95d04c7")
})

test_that("block_lnlp smap_coefficient_covariances works", {
    expect_warning(output <- block_lnlp(block, columns = c("x", "y"),
                                        first_column_time = TRUE, 
                                        method = "s-map", theta = 1, 
                                        save_smap_coefficients = TRUE))
    expect_s3_class(output, "data.frame")
    expect_true("smap_coefficient_covariances" %in% names(output))
    expect_true(is.list(output$smap_coefficient_covariances))
    expect_error(smap_coeff_covariances <- output$smap_coefficient_covariances[[1]], NA)
    expect_true(is.list(smap_coeff_covariances))
    expect_equal(length(smap_coeff_covariances), 200)
    expect_null(smap_coeff_covariances[[200]])
    expect_equal(sapply(smap_coeff_covariances[1:199], dim), 
                 matrix(3, nrow = 2, ncol = 199))
    expect_error(covariance_mat <- do.call(rbind, smap_coeff_covariances[1:199]), NA)
    expect_equal(digest::digest(round(covariance_mat, 4)), 
                 "440dcb1311b66571f942575c77aa0eb8")
})

test_that("block_lnlp works on multivariate time series", {
    expect_warning(output <- block_lnlp(EuStockMarkets, columns = c("DAX", "SMI"),
                                        target_column = "CAC", 
                                        method = "s-map", theta = 1, 
                                        stats_only = FALSE))
    model_output <- round(output$model_output[[1]], 4)
    expect_equal(digest::digest(model_output), "6e2cff4cc75b5251a40b8955b75b321d")
    
    output <- output[, !(names(output) %in% "model_output")]
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "708342ad3fd0a0250345f637b5cc67ec")
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
    expect_error(block_lnlp(df[, 1], first_column_time = TRUE))
    expect_error(block_lnlp(sunspot.year, first_column_time = TRUE))
})